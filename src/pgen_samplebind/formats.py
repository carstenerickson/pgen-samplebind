"""Format detection, EIGENSTRAT/BFILE pre-conversion, tempdir lifecycle.

Per LLD §3.3. Day 1 implements PFILE detection only; BFILE and EIGENSTRAT
shell out to plink2 (deferred to later days per HLD project plan Day 6).
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from .errors import UsageError
from .types import InputDescriptor, InputFormat

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
    (success OR failure).

    Day 1 status: PFILE works end-to-end. BFILE and EIGENSTRAT raise
    NotImplementedError (deferred to project Day 6 per HLD).

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

    # BFILE and EIGENSTRAT both need plink2 conversion to PFILE.
    raise NotImplementedError(
        f"{fmt.value} input is deferred to project Day 6 per HLD §Project plan. "
        f"For Day 1, only PFILE input is supported. To work around, convert manually: "
        f"plink2 --{fmt.value} {base} --make-pgen --out converted"
    )
