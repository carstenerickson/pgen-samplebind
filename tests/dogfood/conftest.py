"""Dogfood-test fixtures and tool-availability detection.

The dogfood test exercises pgen-samplebind on a real AADR-derived panel
(see README.md). Three test tiers, gated by tool availability:

  Tier 1 (default CI): pgen-samplebind merge + schema/count checks.
                       No external tools required.
  Tier 2 (plink2):    PFILE -> BED conversion + .bim sanity.
  Tier 3 (R + AT2):   full qpAdm shootout + md5-compare against
                       the vendored reference.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

FIXTURE_DIR = Path(__file__).parent / "fixtures"


def _has_plink2() -> bool:
    return shutil.which("plink2") is not None


def _has_admixtools() -> bool:
    """True if `R` is on PATH and `library(admixtools)` succeeds.

    Cheap-but-loud probe: runs an R one-liner that exits 0 only when the
    package loads cleanly. Cached per process via module-level memoization.
    """
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False
    try:
        result = subprocess.run(
            [rscript, "-e", "suppressPackageStartupMessages(library(admixtools))"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, OSError):
        return False


HAS_PLINK2 = _has_plink2()
HAS_ADMIXTOOLS = _has_admixtools()


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    """Path to the vendored fixture directory."""
    assert FIXTURE_DIR.is_dir(), f"fixture dir missing: {FIXTURE_DIR}"
    return FIXTURE_DIR


@pytest.fixture(scope="session")
def panel_prefix(fixture_dir: Path) -> Path:
    return fixture_dir / "panel_v66_subset"


@pytest.fixture(scope="session")
def brit_prefix(fixture_dir: Path) -> Path:
    return fixture_dir / "brit_subset_subset"


@pytest.fixture(scope="session")
def target_prefix(fixture_dir: Path) -> Path:
    return fixture_dir / "target_individual"


@pytest.fixture(scope="session")
def qpadm_reference_path(fixture_dir: Path) -> Path:
    return fixture_dir / "qpadm_reference.tsv"
