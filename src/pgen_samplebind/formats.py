"""Format detection, EIGENSTRAT/BFILE pre-conversion, tempdir lifecycle.

Per LLD §3.3.

EIGENSTRAT inputs in two flavors:
- PACKEDANCESTRYMAP (binary, `GENO `/`TGENO ` header): converted to PFILE via
  plink2 `--eigfile --make-pgen` shell-out.
- ASCII per-line (one digit per sample-variant cell, no header): parsed
  natively in `_convert_ascii_eigenstrat` since plink2 doesn't read it.

Both routes write a PFILE triplet into a per-invocation
`tempfile.TemporaryDirectory`; cleanup on success OR failure happens via
the TemporaryDirectory context exit.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

from .errors import IOFailure, UsageError
from .pvar import normalize_chrom_series
from .types import InputDescriptor, InputFormat

# Tempdir prefix used for both EIGENSTRAT and BFILE pre-conversion.
# Tests grep $TMPDIR for this prefix to verify cleanup-on-failure.
_TEMPDIR_PREFIX = "pgen-samplebind-"

_KNOWN_SUFFIXES = frozenset(
    {".pgen", ".pvar", ".psam", ".bed", ".bim", ".fam", ".geno", ".snp", ".ind"}
)


def strip_known_suffix(prefix: Path) -> Path:
    """If prefix has a recognized PFILE/BFILE/EIGENSTRAT extension, strip it.

    Lets users pass either `data` or `data.pgen` and detect identically.
    Also reused by `pseudohaploid.read_sidecar` to locate the
    `<base>.pseudohaploid.json` sidecar against either prefix form.
    """
    if prefix.suffix in _KNOWN_SUFFIXES:
        return prefix.with_suffix("")
    return prefix


def detect_format(prefix: Path) -> InputFormat:
    """Infer format from companion-file presence per HLD §Format detection.

    Tries PFILE first (.pgen/.pvar/.psam triplet), then BFILE (.bed/.bim/.fam),
    then EIGENSTRAT (.geno/.snp/.ind).

    Raises:
        UsageError: prefix doesn't resolve to any supported triplet, or a
            partial triplet is found (e.g., .pgen exists but .pvar is missing).
    """
    base = strip_known_suffix(prefix)
    base_str = str(base)

    pgen = Path(base_str + ".pgen")
    pvar = Path(base_str + ".pvar")
    psam = Path(base_str + ".psam")
    if pgen.exists():
        if pvar.exists() and psam.exists():
            return InputFormat.PFILE
        missing = [str(p) for p in (pvar, psam) if not p.exists()]
        raise UsageError(f"PFILE detected ({pgen}) but missing companion file(s): {missing}")

    bed = Path(base_str + ".bed")
    bim = Path(base_str + ".bim")
    fam = Path(base_str + ".fam")
    if bed.exists():
        if bim.exists() and fam.exists():
            return InputFormat.BFILE
        missing = [str(p) for p in (bim, fam) if not p.exists()]
        raise UsageError(f"BFILE detected ({bed}) but missing companion file(s): {missing}")

    geno = Path(base_str + ".geno")
    snp = Path(base_str + ".snp")
    ind = Path(base_str + ".ind")
    if geno.exists():
        if snp.exists() and ind.exists():
            return InputFormat.EIGENSTRAT
        missing = [str(p) for p in (snp, ind) if not p.exists()]
        raise UsageError(f"EIGENSTRAT detected ({geno}) but missing companion file(s): {missing}")

    raise UsageError(
        f"no PFILE/BFILE/EIGENSTRAT triplet found for prefix {prefix} "
        f"(looked for .pgen, .bed, .geno + companions)"
    )


def check_plink2_available() -> str | None:
    """Return plink2 version string if on PATH, None otherwise.

    Used at startup if any input is EIGENSTRAT/BFILE, and to stamp
    RunSummary.plink2_version.
    """
    plink2 = shutil.which("plink2")
    if plink2 is None:
        return None
    try:
        # plink2 --version prints version on first line of stdout.
        result = subprocess.run(
            [plink2, "--version"],
            shell=False,
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip().splitlines()[0] if result.stdout else None


def _format_chrom_range(include_chrom: tuple[int, ...]) -> str:
    """Convert a tuple of normalized chrom ints to plink2's --chr syntax.

    (1, 2, ..., 22)            → "1-22"
    (1, 2, ..., 22, 23, 24, 26) → "1-22,X,Y,MT"
    Falls back to comma-listing on non-contiguous ranges.
    """
    if not include_chrom:
        return ""
    int_to_named = {23: "X", 24: "Y", 25: "XY", 26: "MT"}
    sorted_chroms = sorted(set(include_chrom))
    nums = [c for c in sorted_chroms if c <= 22]
    others = [c for c in sorted_chroms if c > 22]
    parts: list[str] = []
    if nums:
        if nums == list(range(nums[0], nums[-1] + 1)):
            parts.append(f"{nums[0]}-{nums[-1]}")
        else:
            parts.extend(str(c) for c in nums)
    parts.extend(int_to_named.get(c, str(c)) for c in others)
    return ",".join(parts)


def _run_plink2_convert(
    fmt: InputFormat,
    in_prefix: Path,
    out_prefix: Path,
    include_chrom: tuple[int, ...],
) -> None:
    """Shell out to plink2 to convert BFILE/EIGENSTRAT → PFILE at out_prefix.

    Per LLD §3.3 subprocess hardening pin: shell=False, list args,
    check=False (we capture stderr ourselves), capture_output=True.

    Raises:
        IOFailure: plink2 missing from PATH; plink2 returned non-zero
            (last 20 lines of stderr surfaced in message).
    """
    plink2 = shutil.which("plink2")
    if plink2 is None:
        raise IOFailure(
            f"plink2 not found on PATH; required for {fmt.value} input. "
            f"Install plink2 v2.x (the HLD-verified version is v2.0.0-a.7.1) "
            f"or use PFILE input only."
        )

    flag = "--eigfile" if fmt is InputFormat.EIGENSTRAT else "--bfile"
    chrom_arg = _format_chrom_range(include_chrom)
    cmd: list[str] = [
        plink2,
        flag,
        str(in_prefix),
        "--make-pgen",
        "--out",
        str(out_prefix),
        "--allow-extra-chr",
    ]
    if chrom_arg:
        cmd.extend(["--chr", chrom_arg])

    try:
        result = subprocess.run(
            cmd,
            shell=False,
            check=False,
            capture_output=True,
            text=True,
        )
    except (subprocess.SubprocessError, OSError) as e:
        raise IOFailure(f"plink2 subprocess failed to launch: {e}") from e

    if result.returncode != 0:
        stderr_tail = "\n".join(result.stderr.splitlines()[-20:])
        if not stderr_tail.strip():
            stderr_tail = "(no stderr; check stdout)\n" + "\n".join(
                result.stdout.splitlines()[-20:]
            )
        raise IOFailure(
            f"plink2 subprocess failed converting {fmt.value} input "
            f"{in_prefix} (exit {result.returncode}):\n"
            f"command: {' '.join(cmd)}\n"
            f"stderr (last 20 lines):\n{stderr_tail}"
        )


def _is_packedancestrymap(geno_path: Path) -> bool:
    """Sniff first 6 bytes of an EIGENSTRAT .geno: PACKEDANCESTRYMAP files
    start with `GENO ` (binary header) or `TGENO ` (transposed binary).
    Anything else is treated as ASCII per-line.
    """
    try:
        with open(geno_path, "rb") as f:
            head = f.read(6)
    except OSError:
        return False
    return head.startswith(b"GENO ") or head.startswith(b"TGENO ")


def _convert_ascii_eigenstrat(
    in_prefix: Path,
    out_prefix: Path,
    include_chrom: tuple[int, ...],
) -> None:
    """Native parser for ASCII per-line EIGENSTRAT → PFILE at out_prefix.

    Format spec (no plink2 dependency):
    - `.ind`: one line per sample. Whitespace-delimited 3 cols: IID sex pop.
    - `.snp`: one line per variant. Whitespace-delimited 6 cols:
        rsID chrom cM pos REF ALT.
    - `.geno`: one line per variant. Each line = N chars (one per sample),
        values 0/1/2 = count of REF allele, 9 = missing. Trailing newline.

    Conversion to plink/pgenlib convention:
    - genotype = 2 - eig_value for {0,1,2}; missing → -9.
    - Apply the autosome filter (`include_chrom`) at parse time.
    - Preserve the cM column in the output `.pvar`.

    Memory: reads the full `.geno` as a (n_variants x (n_samples+1)) uint8
    matrix via `np.fromfile`. For 1240k-scale (1.2M variants x ~1000 samples
    ≈ 1.3 GB) this fits comfortably in 8GB+ environments. Streamed write via
    PgenWriter keeps RSS dominated by the parsed matrix, not by buffered
    output.

    Raises:
        IOFailure: file unreadable; .geno size doesn't match
            (n_variants_raw x (n_samples + 1)) bytes; pgenlib write failure.
    """
    import pgenlib

    geno_path = Path(str(in_prefix) + ".geno")
    snp_path = Path(str(in_prefix) + ".snp")
    ind_path = Path(str(in_prefix) + ".ind")

    # 1. .ind — sample table.
    try:
        ind_df = pd.read_csv(
            ind_path,
            sep=r"\s+",
            header=None,
            names=["IID", "SEX", "PHENO1"],
            dtype=str,
            engine="python",
        )
    except (OSError, pd.errors.ParserError) as e:
        raise IOFailure(f"cannot parse ASCII EIGENSTRAT .ind {ind_path}: {e}") from e
    n_samples = len(ind_df)

    # Write .psam with FID=0 (orchestrator's add_fid_from_pop sets FID=POP later).
    # Sex mapping: M→1, F→2, U/anything→0 (plink2 spec).
    psam_df = pd.DataFrame(
        {
            "#FID": ["0"] * n_samples,
            "IID": ind_df["IID"].values,
            "SEX": ind_df["SEX"].map({"M": "1", "F": "2"}).fillna("0").values,
            "PHENO1": ind_df["PHENO1"].values,
        }
    )
    out_psam = Path(str(out_prefix) + ".psam")
    try:
        psam_df.to_csv(out_psam, sep="\t", index=False, lineterminator="\n")
    except OSError as e:
        raise IOFailure(f"cannot write {out_psam}: {e}") from e

    # 2. .snp — variant table.
    try:
        snp_df = pd.read_csv(
            snp_path,
            sep=r"\s+",
            header=None,
            names=["ID", "CHROM_RAW", "CM", "POS", "REF", "ALT"],
            dtype={
                "ID": str,
                "CHROM_RAW": str,
                "CM": float,
                "POS": int,
                "REF": str,
                "ALT": str,
            },
            engine="python",
        )
    except (OSError, pd.errors.ParserError) as e:
        raise IOFailure(f"cannot parse ASCII EIGENSTRAT .snp {snp_path}: {e}") from e
    n_variants_raw = len(snp_df)

    # Normalize chromosome ints, then apply autosome filter via include_chrom.
    snp_df["CHROM"] = normalize_chrom_series(snp_df["CHROM_RAW"])
    keep_mask = snp_df["CHROM"].isin(include_chrom).to_numpy()
    snp_kept = snp_df.loc[keep_mask].copy()
    n_variants = len(snp_kept)

    out_pvar = Path(str(out_prefix) + ".pvar")
    pvar_df = pd.DataFrame(
        {
            "#CHROM": snp_kept["CHROM"].astype(int),
            "POS": snp_kept["POS"],
            "ID": snp_kept["ID"],
            "REF": snp_kept["REF"].str.upper(),
            "ALT": snp_kept["ALT"].str.upper(),
            "CM": snp_kept["CM"],
        }
    )
    try:
        pvar_df.to_csv(out_pvar, sep="\t", index=False, lineterminator="\n")
    except OSError as e:
        raise IOFailure(f"cannot write {out_pvar}: {e}") from e

    # 3. .geno — fixed-width per-line genotype matrix.
    try:
        raw = np.fromfile(geno_path, dtype=np.uint8)
    except OSError as e:
        raise IOFailure(f"cannot read ASCII EIGENSTRAT .geno {geno_path}: {e}") from e

    expected = n_variants_raw * (n_samples + 1)
    if raw.size != expected:
        raise IOFailure(
            f"ASCII EIGENSTRAT .geno size mismatch at {geno_path}: expected "
            f"{n_variants_raw} variants x ({n_samples} samples + 1 newline) = "
            f"{expected} bytes, got {raw.size} bytes. Possible causes: line-count "
            f"disagreement with .snp ({snp_path}), trailing newlines, CRLF line "
            f"endings, or PACKEDANCESTRYMAP binary file mis-detected as ASCII."
        )

    # Drop the trailing-newline column and convert ASCII '0'/'1'/'2'/'9' digits.
    matrix_chars = raw.reshape(n_variants_raw, n_samples + 1)[:, :n_samples]
    matrix_int = matrix_chars.astype(np.int16) - ord("0")  # int16 to hold 9 cleanly
    # eig 0/1/2 (count of REF) → plink 0/1/2 (count of ALT) via 2 - x; 9 → -9.
    matrix_int = np.where(matrix_int == 9, -9, 2 - matrix_int).astype(np.int8)
    matrix_kept = matrix_int[keep_mask]

    # 4. Stream variants to PFILE via PgenWriter.
    out_pgen = Path(str(out_prefix) + ".pgen")
    try:
        writer = pgenlib.PgenWriter(
            str(out_pgen).encode(),
            sample_ct=n_samples,
            variant_ct=n_variants,
            nonref_flags=False,
        )
    except Exception as e:
        raise IOFailure(f"cannot open PgenWriter for {out_pgen}: {e}") from e
    try:
        # PgenWriter.append_biallelic_batch takes (block_size, n_samples) int8.
        # Stream in 1024-variant blocks to amortize the call overhead.
        block = 1024
        for start in range(0, n_variants, block):
            end = min(start + block, n_variants)
            chunk = np.ascontiguousarray(matrix_kept[start:end], dtype=np.int8)
            writer.append_biallelic_batch(chunk)
    finally:
        writer.close()


@contextmanager
def prepared_input(
    prefix: Path,
    is_target: bool = False,
    include_chrom: tuple[int, ...] = tuple(range(1, 23)),
) -> Iterator[InputDescriptor]:
    """Resolve a user prefix to an InputDescriptor with valid PFILE paths.

    PFILE: pure resolution. BFILE/EIGENSTRAT: shells out to
    `plink2 --bfile/--eigfile <prefix> --make-pgen` in a per-invocation
    `tempfile.TemporaryDirectory`. Tempdir is cleaned on context exit
    (success OR failure) — `TemporaryDirectory.__exit__` runs even on
    uncaught exceptions, so partial-conversion artifacts never leak.

    Per HLD §EIGENSTRAT a7.x quirks: plink2 v2.0.0-a.7.x's
    `--eigfile --make-pgen` preserves the population label as `PHENO1`
    (not stripped), emits no FID column, and doesn't need the
    .ind-re-read awk dance. The orchestrator's existing detect_population_column
    + rename_to_pop + add_fid_from_pop flow handles `PHENO1 → POP` and
    `FID = POP` post-conversion (no special path needed here).

    Per HLD §Format detection: `--chr 1-22 --allow-extra-chr` is the
    default chrom filter (autosomes); override via `include_chrom` for
    workflows that need sex chromosomes.

    Raises:
        UsageError: format unrecognized.
        IOFailure: plink2 missing from PATH while BFILE/EIGENSTRAT supplied;
            plink2 returned non-zero (last 20 lines of stderr surfaced).
    """
    fmt = detect_format(prefix)
    base = strip_known_suffix(prefix)
    base_str = str(base)

    if fmt is InputFormat.PFILE:
        desc = InputDescriptor(
            path=prefix,
            pgen_path=Path(base_str + ".pgen"),
            pvar_path=Path(base_str + ".pvar"),
            psam_path=Path(base_str + ".psam"),
            fmt=fmt,
            is_target=is_target,
            eigfile_tempdir=None,
        )
        yield desc
        return

    # BFILE and EIGENSTRAT: convert to PFILE in a per-invocation tempdir.
    with tempfile.TemporaryDirectory(prefix=_TEMPDIR_PREFIX) as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        out_prefix = tmpdir / base.name

        # EIGENSTRAT splits into PACKEDANCESTRYMAP (binary; plink2 reads it)
        # vs ASCII per-line (plink2 doesn't; native parser handles it).
        if fmt is InputFormat.EIGENSTRAT and not _is_packedancestrymap(Path(base_str + ".geno")):
            _convert_ascii_eigenstrat(base, out_prefix, include_chrom)
        else:
            _run_plink2_convert(fmt, base, out_prefix, include_chrom)

        desc = InputDescriptor(
            path=prefix,
            pgen_path=Path(str(out_prefix) + ".pgen"),
            pvar_path=Path(str(out_prefix) + ".pvar"),
            psam_path=Path(str(out_prefix) + ".psam"),
            fmt=fmt,
            is_target=is_target,
            eigfile_tempdir=tmpdir,
        )
        yield desc
        # tmpdir auto-removed on context exit (TemporaryDirectory)
