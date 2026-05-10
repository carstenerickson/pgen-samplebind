""".pvar parsing, biallelic-SNP filter, multi-allelic startup check, chromosome normalization.

Per LLD §3.4. Pandas-driven for the speed and memory wins documented in the pgenlib-verify-report.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .errors import InvariantViolation, IOFailure

_VALID_NUCLEOTIDES = frozenset({"A", "C", "G", "T"})


def normalize_chrom(s: str | int) -> int:
    """Chromosome string → int per HLD §Variant alignment.

    1-22 numeric, X/chrX/23 → 23, Y/chrY/24 → 24, MT/chrMT/chrM/26 → 26.
    `chr` prefix stripped if present.

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
    if raw_upper == "X":
        return 23
    if raw_upper == "Y":
        return 24
    if raw_upper in {"MT", "M"}:
        return 26
    raise InvariantViolation(f"unparseable chromosome string: {s!r}")


def count_raw_variants(path: Path) -> int:
    """Pre-filter row count of a .pvar (data lines, excluding header lines).

    Header lines are any line starting with `#` (both `##`-prefixed metadata
    and the single-`#` column header).
    """
    n = 0
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if not line.startswith("#"):
                    n += 1
    except OSError as e:
        raise IOFailure(f"cannot read {path}: {e}") from e
    return n


def _find_header_line(path: Path) -> int:
    """Locate the line index (0-based) of the `#CHROM` header in a .pvar.

    Plink2 .pvar files may have ``##``-prefixed metadata lines before the
    column header (which starts with a single `#`).
    """
    try:
        with open(path, encoding="utf-8") as f:
            for i, line in enumerate(f):
                if line.startswith("#CHROM"):
                    return i
    except OSError as e:
        raise IOFailure(f"cannot read {path}: {e}") from e
    raise IOFailure(f"no #CHROM header found in {path}")


def read_pvar(path: Path) -> pd.DataFrame:
    """Parse a .pvar via pandas. Apply the biallelic-SNP filter vectorized.

    Returns a DataFrame with columns: chrom (int8), pos (int64), id (object),
    ref (object, len 1, ACGT), alt (object, len 1, ACGT). Failed-filter rows
    are dropped. Pre-filter count via count_raw_variants() if needed.

    Raises:
        IOFailure: file unreadable or unparseable.
    """
    header_idx = _find_header_line(path)
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            skiprows=header_idx,
            header=0,
            dtype=str,  # parse everything as string first; we'll cast
            na_filter=False,
        )
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
    df["CHROM"] = df["CHROM"].map(normalize_chrom).astype("int8")
    df["POS"] = df["POS"].astype("int64")
    df["REF"] = df["REF"].str.upper()
    df["ALT"] = df["ALT"].str.upper()

    return df.rename(
        columns={"CHROM": "chrom", "POS": "pos", "ID": "id", "REF": "ref", "ALT": "alt"}
    ).reset_index(drop=True)


def validate_unique_keys(df: pd.DataFrame, key: str) -> None:
    """Assert no duplicate keys in canonical input. Per LLD §3.4 / HLD §Variant alignment.

    Called only for input[0] before alignment. Catches the corner where input[0] was
    itself produced by a buggy merge or malformed source file.

    Raises:
        InvariantViolation: duplicate found; message names the first duplicate key.
    """
    if key == "chr_pos":
        keys = list(zip(df["chrom"].tolist(), df["pos"].tolist(), strict=True))
        col_label = "(chrom, pos)"
    elif key == "id":
        keys = df["id"].tolist()
        col_label = "id"
    else:
        raise InvariantViolation(f"unknown variant_key: {key!r}")

    seen: set[object] = set()
    for k in keys:
        if k in seen:
            raise InvariantViolation(
                f"duplicate canonical {col_label} key: {k!r}. Input[0] must have unique "
                f"variant keys (HLD §Variant alignment)."
            )
        seen.add(k)


def check_max_alleles(pgen_path: Path) -> None:
    """Open .pgen via pgenlib.PvarReader and assert get_max_allele_ct() == 2.

    Per HLD §Non-SNP handling: multi-allelic input causes uncatchable C-layer
    SIGSEGV in PgenReader.read_range, so we must reject at startup.

    Raises:
        InvariantViolation: max_allele_ct > 2.
    """
    import pgenlib  # lazy import; pgenlib's load is non-trivial

    pvar_path = pgen_path.with_suffix(".pvar")
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
