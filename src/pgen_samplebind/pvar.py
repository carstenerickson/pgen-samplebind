""".pvar parsing, biallelic-SNP filter, multi-allelic startup check, chromosome normalization.

Per LLD §3.4. Pandas-driven for the speed and memory wins documented in the pgenlib-verify-report.
"""

from __future__ import annotations

import io
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO

import pandas as pd

from .errors import InvariantViolation, IOFailure


def is_zst_path(path: Path) -> bool:
    """True if path ends in `.zst` (a zstd-compressed file).

    Used to gate decompression. Plink2 v2.0.0-a.6+ and reference panels
    distributed via HuggingFace / Dataverse routinely ship `.pvar.zst`
    rather than uncompressed `.pvar` (cuts ~10x file size for typical
    HGDP+1kGP panels).
    """
    return path.suffix == ".zst"


@contextmanager
def open_pvar_text(path: Path) -> Iterator[IO[str]]:
    """Open a `.pvar` or `.pvar.zst` as a text-mode line iterator.

    Streams zstandard decompression so we don't materialize the
    uncompressed file in memory. Used by `count_raw_variants` and
    `_find_header_line`; `read_pvar` uses pandas's built-in compression
    handling instead, since it pipes directly into the C parser.

    Raises:
        IOFailure: path unreadable.
    """
    try:
        if is_zst_path(path):
            import zstandard

            with open(path, "rb") as raw:
                stream = zstandard.ZstdDecompressor().stream_reader(raw)
                yield io.TextIOWrapper(stream, encoding="utf-8")
        else:
            with open(path, encoding="utf-8") as f:
                yield f
    except OSError as e:
        raise IOFailure(f"cannot read {path}: {e}") from e


_VALID_NUCLEOTIDES = frozenset({"A", "C", "G", "T"})


# Letter → int mapping shared by scalar and vectorized normalizers.
_CHROM_LETTER_MAP: dict[str, int] = {"X": 23, "Y": 24, "MT": 26, "M": 26}


def normalize_chrom(s: str | int) -> int:
    """Chromosome string → int per HLD §Variant alignment.

    1-22 numeric, X/chrX/23 → 23, Y/chrY/24 → 24, MT/chrMT/chrM/26 → 26.
    `chr` prefix stripped if present.

    Scalar variant. For hot-path bulk normalization (read_pvar at 1240k
    scale, ASCII-EIGENSTRAT loader), call `normalize_chrom_series` against
    the Series directly — the pandas `.map(normalize_chrom)` path is a
    per-cell Python apply and is ~5-10x slower at scale.

    Raises:
        InvariantViolation: unparseable string.
    """
    if isinstance(s, int):
        if 1 <= s <= 26:
            return s
        raise InvariantViolation(f"chromosome integer out of range 1-26: {s}")
    raw = str(s).strip()
    if raw.lower().startswith("chr"):
        raw = raw[3:]
    raw_upper = raw.upper()
    if raw.isdigit():
        n = int(raw)
        if 1 <= n <= 26:
            return n
        raise InvariantViolation(f"chromosome integer out of range 1-26: {raw}")
    if raw_upper in _CHROM_LETTER_MAP:
        return _CHROM_LETTER_MAP[raw_upper]
    raise InvariantViolation(f"unparseable chromosome string: {s!r}")


def normalize_chrom_series(s: pd.Series) -> pd.Series:
    """Vectorized chromosome normalization for the read_pvar / ASCII-EIGENSTRAT
    hot paths.

    Same semantics as `normalize_chrom` but in a single C-level pandas pass:
    - Strip optional `chr` prefix (case-insensitive).
    - Numeric values in [1, 26] pass through.
    - X/Y/MT/M (case-insensitive) map to 23 / 24 / 26 / 26.
    - Anything else raises `InvariantViolation` naming the first offender.

    Returns an int8 Series. ~5-10x faster than `.map(normalize_chrom)` at
    1240k scale because the per-cell Python overhead is removed.
    """
    raw = s.astype(str).str.strip()
    # Case-insensitive `chr` prefix strip (handles chr/CHR/Chr/cHr).
    no_prefix = raw.str.replace(r"(?i)^chr", "", regex=True, n=1)
    upper = no_prefix.str.upper()

    numeric_mask = upper.str.fullmatch(r"\d+", na=False)
    letter_mask = upper.isin(_CHROM_LETTER_MAP.keys())
    unknown_mask = ~(numeric_mask | letter_mask)
    if unknown_mask.any():
        first_bad = s.loc[unknown_mask].iloc[0]
        raise InvariantViolation(f"unparseable chromosome string: {first_bad!r}")

    out = pd.Series(0, index=s.index, dtype="int64")
    if numeric_mask.any():
        out.loc[numeric_mask] = upper.loc[numeric_mask].astype("int64")
    if letter_mask.any():
        out.loc[letter_mask] = upper.loc[letter_mask].map(_CHROM_LETTER_MAP).astype("int64")

    out_of_range = (out < 1) | (out > 26)
    if out_of_range.any():
        first_idx = out.loc[out_of_range].index[0]
        raw_val = s.loc[first_idx]
        raise InvariantViolation(
            f"chromosome integer out of range 1-26: {int(out.loc[first_idx])} (from {raw_val!r})"
        )

    return out.astype("int8")


def count_raw_variants(path: Path) -> int:
    """Pre-filter row count of a .pvar (data lines, excluding header lines).

    Header lines are any line starting with `#` (both `##`-prefixed metadata
    and the single-`#` column header). Accepts `.pvar` or `.pvar.zst`.
    """
    n = 0
    with open_pvar_text(path) as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n


def _find_header_line(path: Path) -> int:
    """Locate the line index (0-based) of the `#CHROM` header in a .pvar.

    Plink2 .pvar files may have ``##``-prefixed metadata lines before the
    column header (which starts with a single `#`). Accepts `.pvar` or
    `.pvar.zst`.
    """
    with open_pvar_text(path) as f:
        for i, line in enumerate(f):
            if line.startswith("#CHROM"):
                return i
    raise IOFailure(f"no #CHROM header found in {path}")


def read_pvar(path: Path) -> pd.DataFrame:
    """Parse a .pvar via pandas. Apply the biallelic-SNP filter vectorized.

    Returns a DataFrame with columns: chrom (int8), pos (int64), id (object),
    ref (object, len 1, ACGT), alt (object, len 1, ACGT), cm (float64).
    Failed-filter rows are dropped. The `cm` column carries genetic-position
    in centiMorgans (preserved from input .pvar's CM column when present;
    0.0 otherwise). cM is required end-to-end so that downstream tools using
    Morgan-spaced jackknife blocks (e.g., AT2 `extract_f2 blgsize=0.05`)
    partition variants the same way as if they'd consumed the original input.

    Raises:
        IOFailure: file unreadable or unparseable.
    """
    header_idx = _find_header_line(path)
    read_csv_kwargs: dict[str, object] = dict(
        sep="\t",
        skiprows=header_idx,
        header=0,
        dtype=str,  # parse everything as string first; we'll cast
        na_filter=False,
    )
    if is_zst_path(path):
        read_csv_kwargs["compression"] = "zstd"
    try:
        df = pd.read_csv(path, **read_csv_kwargs)  # type: ignore[call-overload]
    except (OSError, pd.errors.ParserError) as e:
        raise IOFailure(f"cannot parse {path}: {e}") from e

    # Strip the leading '#' from the first column (e.g., '#CHROM' → 'CHROM')
    df.columns = [c.lstrip("#") for c in df.columns]

    required = {"CHROM", "POS", "ID", "REF", "ALT"}
    missing = required - set(df.columns)
    if missing:
        raise IOFailure(f"{path} missing required columns: {sorted(missing)}")

    # Biallelic-SNP filter (vectorized): single-character ACGT REF and ALT.
    is_biallelic_snp = (
        df["REF"].str.len().eq(1)
        & df["ALT"].str.len().eq(1)
        & df["REF"].isin(_VALID_NUCLEOTIDES)
        & df["ALT"].isin(_VALID_NUCLEOTIDES)
    )
    df = df.loc[is_biallelic_snp].copy()

    # Normalize chrom to int8; pos to int64; uppercase REF/ALT.
    df["CHROM"] = normalize_chrom_series(df["CHROM"])
    df["POS"] = df["POS"].astype("int64")
    df["REF"] = df["REF"].str.upper()
    df["ALT"] = df["ALT"].str.upper()

    # CM column: optional in plink2 .pvar (emitted by --eigfile / --bfile
    # converters when source had cM, omitted otherwise). Default to 0.0.
    if "CM" in df.columns:
        df["CM"] = pd.to_numeric(df["CM"], errors="coerce").fillna(0.0).astype("float64")
    else:
        df["CM"] = 0.0

    return df.rename(
        columns={
            "CHROM": "chrom",
            "POS": "pos",
            "ID": "id",
            "REF": "ref",
            "ALT": "alt",
            "CM": "cm",
        }
    ).reset_index(drop=True)


def validate_unique_keys(df: pd.DataFrame, key: str) -> None:
    """Assert no duplicate keys in canonical input. Per LLD §3.4 / HLD §Variant alignment.

    Called only for input[0] before alignment. Catches the corner where input[0] was
    itself produced by a buggy merge or malformed source file.

    Vectorized via pandas `duplicated`: a single C-level scan instead of the
    n_variants-step Python set/loop the v0.1 implementation used. ~5-10x
    faster at 1240k scale; constant per-row Python overhead removed.

    Raises:
        InvariantViolation: duplicate found; message names the first duplicate key.
    """
    if key == "chr_pos":
        subset = ["chrom", "pos"]
        col_label = "(chrom, pos)"
    elif key == "id":
        subset = ["id"]
        col_label = "id"
    else:
        raise InvariantViolation(f"unknown variant_key: {key!r}")

    dup_mask = df.duplicated(subset=subset, keep="first")
    if not dup_mask.any():
        return

    first_dup_row = df.loc[dup_mask].iloc[0]
    if key == "chr_pos":
        first_dup_key: tuple[int, int] | str = (
            int(first_dup_row["chrom"]),
            int(first_dup_row["pos"]),
        )
    else:
        first_dup_key = str(first_dup_row["id"])
    raise InvariantViolation(
        f"duplicate canonical {col_label} key: {first_dup_key!r}. Input[0] must have unique "
        f"variant keys (HLD §Variant alignment)."
    )


def check_max_alleles(pgen_path: Path) -> None:
    """Open .pgen via pgenlib.PvarReader and assert get_max_allele_ct() == 2.

    Per HLD §Non-SNP handling: multi-allelic input causes uncatchable C-layer
    SIGSEGV in PgenReader.read_range, so we must reject at startup.

    Raises:
        InvariantViolation: max_allele_ct > 2.
    """
    import pgenlib  # lazy import; pgenlib's load is non-trivial

    # Plink2 may emit either `.pvar` or `.pvar.zst` next to the .pgen; pgenlib's
    # PvarReader handles both transparently (libzstd-linked at the C layer).
    base = pgen_path.with_suffix("")
    pvar_path = Path(str(base) + ".pvar")
    if not pvar_path.exists():
        zst_path = Path(str(base) + ".pvar.zst")
        if zst_path.exists():
            pvar_path = zst_path
    try:
        pv = pgenlib.PvarReader(str(pvar_path).encode())
    except Exception as e:
        raise IOFailure(f"cannot open {pvar_path} for max-allele check: {e}") from e
    try:
        max_alleles = pv.get_max_allele_ct()
    finally:
        # pgenlib readers don't always have an explicit close; rely on GC.
        pass
    if max_alleles > 2:
        raise InvariantViolation(
            f"{pgen_path} contains multi-allelic variants (max_allele_ct={max_alleles}); "
            f"pgen-samplebind requires biallelic input. Preprocess with: "
            f"plink2 --pfile {pgen_path.with_suffix('')} --max-alleles 2 --snps-only "
            f"--make-pgen --out preprocessed"
        )
