"""Format detection, EIGENSTRAT/BFILE pre-conversion, tempdir lifecycle.

Per LLD §3.3. Day 6 lights up BFILE and EIGENSTRAT input via a plink2
shell-out into a per-invocation `tempfile.TemporaryDirectory`. Cleanup
on success OR failure happens via the TemporaryDirectory context exit.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from .errors import IOFailure, UsageError
from .types import InputDescriptor, InputFormat

# Tempdir prefix used for both EIGENSTRAT and BFILE pre-conversion.
# Tests grep $TMPDIR for this prefix to verify cleanup-on-failure.
_TEMPDIR_PREFIX = "pgen-samplebind-"

_KNOWN_SUFFIXES = frozenset(
    {".pgen", ".pvar", ".psam", ".bed", ".bim", ".fam", ".geno", ".snp", ".ind"}
)


def _strip_known_suffix(prefix: Path) -> Path:
    """If prefix has a recognized PFILE/BFILE/EIGENSTRAT extension, strip it.

    Lets users pass either `data` or `data.pgen` and detect identically.
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
    base = _strip_known_suffix(prefix)
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
    base = _strip_known_suffix(prefix)
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
