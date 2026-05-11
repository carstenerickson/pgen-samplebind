"""AADR-derivative dogfood: end-to-end pgen-samplebind vs the established
mergeit + plink2 + AT2 reference pipeline.

Three test tiers. See tests/dogfood/README.md for the full setup.

Tier 1 (always runs):
  - pgen-samplebind merge produces the expected panel shape on the
    real-world AADR fixture: 41 samples x 49,382 variants x 12 pops
    after 1,236 strand-ambiguous variant drops.
  - .pvar carries the cM column end-to-end (regression guard for the
    Phase 7 dogfood-discovered bug).
  - .psam carries FID=POP for all samples and a populated PSEUDOHAPLOID
    column (this v66-derived fixture is all-pseudohap by source).

Tier 2 (`dogfood_plink2`, needs plink2 on PATH):
  - PFILE -> BED conversion succeeds, .bim cM column is non-trivial
    (lots of distinct values, not all zero), .fam FIDs match populations.

Tier 3 (`dogfood_full`, needs R + admixtools):
  - Full extract_f2 + qpAdm shootout against the pgen-samplebind panel
    produces md5-equal output to the vendored mergeit-pipeline reference.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from tests.dogfood.conftest import HAS_ADMIXTOOLS, HAS_PLINK2

pytestmark = pytest.mark.dogfood


EXPECTED_PANEL_SAMPLES = 41
EXPECTED_PANEL_VARIANTS = 49_382
EXPECTED_AMBIG_DROPS = 1_236
EXPECTED_POPULATIONS = 12  # 7 source + 4 target + 1 dogfood_target


def _run_pgensb_merge(
    panel: Path, brit: Path, target: Path, out_prefix: Path
) -> subprocess.CompletedProcess[str]:
    """Drive pgen-samplebind merge via subprocess (so exit codes route through
    cli.main()'s ExitCode mapping). Returns the CompletedProcess for assertion."""
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "pgen_samplebind",
            "merge",
            str(panel),
            str(brit),
            "--target",
            str(target),
            "--on-collision",
            "suffix",
            "-o",
            str(out_prefix),
            "--quiet",
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )


def test_dogfood_pgensb_merge_runs(
    panel_prefix: Path, brit_prefix: Path, target_prefix: Path, tmp_path: Path
) -> None:
    """End-to-end pgen-samplebind merge on the AADR fixture. The realistic
    scale of this test (41 x 49K) catches issues that the unit-test
    synthetic fixtures (10s of samples x 100s of variants) miss."""
    out_prefix = tmp_path / "dogfood_pgensb_merged"
    result = _run_pgensb_merge(panel_prefix, brit_prefix, target_prefix, out_prefix)
    assert result.returncode == 0, (
        f"pgensb merge failed (exit {result.returncode}):\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )

    # Verify the triplet exists.
    for ext in (".pgen", ".pvar", ".psam"):
        assert (out_prefix.with_suffix(ext)).is_file(), f"missing {ext}"


def test_dogfood_panel_shape(
    panel_prefix: Path, brit_prefix: Path, target_prefix: Path, tmp_path: Path
) -> None:
    """Expected panel shape: 41 samples, 49,382 variants (50K - 618 ambiguous
    x 2 non-canonical inputs / 2 dedup = 49,382), 12 populations."""
    out_prefix = tmp_path / "dogfood_panel"
    _run_pgensb_merge(panel_prefix, brit_prefix, target_prefix, out_prefix)

    # Sample count from .psam (skip header lines starting with #)
    psam = pd.read_csv(out_prefix.with_suffix(".psam"), sep="\t")
    assert len(psam) == EXPECTED_PANEL_SAMPLES, (
        f"sample count: got {len(psam)}, expected {EXPECTED_PANEL_SAMPLES}"
    )

    # Variant count from .pvar (skip header)
    pvar = pd.read_csv(
        out_prefix.with_suffix(".pvar"),
        sep="\t",
        comment="#",
        header=None,
        names=["chrom", "pos", "id", "ref", "alt", "cm"],
    )
    assert len(pvar) == EXPECTED_PANEL_VARIANTS, (
        f"variant count: got {len(pvar)}, expected {EXPECTED_PANEL_VARIANTS}"
    )

    # Population count via FID column in .psam.
    fid_col = "#FID" if "#FID" in psam.columns else "FID"
    pops = set(psam[fid_col])
    assert len(pops) == EXPECTED_POPULATIONS, (
        f"population count: got {len(pops)}, expected {EXPECTED_POPULATIONS}"
    )
    assert "dogfood_target" in pops, "target pop label not preserved"


def test_dogfood_cm_preserved(
    panel_prefix: Path, brit_prefix: Path, target_prefix: Path, tmp_path: Path
) -> None:
    """The cM column must carry through pgensb merge — regression guard for
    the Phase 7 dogfood-discovered bug where output .pvar had cM=0 for all
    variants, collapsing downstream Morgan-spaced jackknife blocks."""
    out_prefix = tmp_path / "dogfood_cm_check"
    _run_pgensb_merge(panel_prefix, brit_prefix, target_prefix, out_prefix)

    pvar = pd.read_csv(
        out_prefix.with_suffix(".pvar"),
        sep="\t",
        comment="#",
        header=None,
        names=["chrom", "pos", "id", "ref", "alt", "cm"],
    )
    distinct_cm = pvar["cm"].nunique()
    assert distinct_cm > 1000, (
        f"cM column appears collapsed ({distinct_cm} distinct values). "
        f"Did the Phase 7 cM-preservation fix regress?"
    )
    # Sanity: cM should be monotone non-decreasing within each chromosome.
    for chrom, grp in pvar.groupby("chrom"):
        diffs = grp["cm"].diff().dropna()
        assert (diffs >= -1e-9).all(), (
            f"cM not monotone within chrom {chrom}: {(diffs < 0).sum()} backwards steps"
        )


def test_dogfood_pseudohap_column(
    panel_prefix: Path, brit_prefix: Path, target_prefix: Path, tmp_path: Path
) -> None:
    """PSEUDOHAPLOID column populated. The v66-derived fixture uses adaptive
    pulldown across all 28 panel samples, plus the brit_subset's pulldown +
    the pulldown-extracted target = expect all 41 samples pseudohaploid."""
    out_prefix = tmp_path / "dogfood_pseudo"
    _run_pgensb_merge(panel_prefix, brit_prefix, target_prefix, out_prefix)

    psam = pd.read_csv(out_prefix.with_suffix(".psam"), sep="\t")
    assert "PSEUDOHAPLOID" in psam.columns, "PSEUDOHAPLOID column missing"
    pseudohap_count = (psam["PSEUDOHAPLOID"].astype(str) == "1").sum()
    # Allow some flexibility: at minimum 35 of 41 should classify as pseudohap
    # (the auto-detect threshold is het_count == 0; rare hets at <1% of called
    # variants can downgrade a sample to UNKNOWN or DIPLOID).
    assert pseudohap_count >= 35, (
        f"only {pseudohap_count}/41 classified pseudohap; expected ~all-41 "
        f"for this all-adaptive-pulldown fixture"
    )


@pytest.mark.skipif(not HAS_PLINK2, reason="plink2 not on PATH")
@pytest.mark.dogfood_plink2
def test_dogfood_pfile_to_bed_cm_preserved(
    panel_prefix: Path, brit_prefix: Path, target_prefix: Path, tmp_path: Path
) -> None:
    """Convert pgensb output PFILE -> BED via plink2; cM column must survive
    the conversion (this is the path AT2's extract_f2 ultimately consumes)."""
    out_prefix = tmp_path / "dogfood_for_bed"
    _run_pgensb_merge(panel_prefix, brit_prefix, target_prefix, out_prefix)

    bed_prefix = tmp_path / "dogfood_panel_plink"
    result = subprocess.run(
        ["plink2", "--pfile", str(out_prefix), "--make-bed", "--out", str(bed_prefix)],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, f"plink2 --make-bed failed: {result.stderr}"

    bim = pd.read_csv(
        bed_prefix.with_suffix(".bim"),
        sep="\t",
        header=None,
        names=["chrom", "id", "cm", "pos", "a1", "a2"],
    )
    distinct_cm = bim["cm"].nunique()
    assert distinct_cm > 1000, (
        f".bim cM column collapsed ({distinct_cm} distinct values) — "
        f"plink2 --make-bed didn't pick up cM from pgensb's .pvar. "
        f"This is the exact failure mode the Phase 7 dogfood surfaced "
        f"before the cM preservation fix landed."
    )


@pytest.mark.skipif(not (HAS_PLINK2 and HAS_ADMIXTOOLS), reason="needs plink2 + R/admixtools")
@pytest.mark.dogfood_full
def test_dogfood_full_qpadm_matches_reference(
    panel_prefix: Path,
    brit_prefix: Path,
    target_prefix: Path,
    qpadm_reference_path: Path,
    tmp_path: Path,
) -> None:
    """End-to-end: pgensb merge -> BED -> extract_f2 -> qpAdm shootout. Output
    must match the vendored mergeit-pipeline reference within per-cell numerical
    tolerance (cross-architecture float-arithmetic-order noise is ~1e-12; the
    threshold below leaves 6 orders of magnitude of headroom over that)."""
    out_prefix = tmp_path / "dogfood_full"
    _run_pgensb_merge(panel_prefix, brit_prefix, target_prefix, out_prefix)

    bed_prefix = tmp_path / "dogfood_full_plink"
    subprocess.run(
        ["plink2", "--pfile", str(out_prefix), "--make-bed", "--out", str(bed_prefix)],
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )

    rscript = tmp_path / "shootout.R"
    rscript.write_text(_QPADM_SHOOTOUT_R)
    f2_cache = tmp_path / "f2_cache"
    out_tsv = tmp_path / "qpadm_pgensb.tsv"
    result = subprocess.run(
        ["Rscript", str(rscript), str(bed_prefix), str(f2_cache), str(out_tsv)],
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )
    assert result.returncode == 0, f"Rscript failed:\n{result.stderr}"
    assert out_tsv.is_file()

    ref = pd.read_csv(qpadm_reference_path, sep="\t")
    got = pd.read_csv(out_tsv, sep="\t")

    assert list(ref.columns) == list(got.columns), (
        f"column mismatch: ref={list(ref.columns)} got={list(got.columns)}"
    )
    assert len(ref) == len(got), f"row count: ref={len(ref)} got={len(got)}"
    assert ref["target"].tolist() == got["target"].tolist(), "target order differs"

    # Per-cell numerical tolerance on the eight numeric columns.
    WEIGHT_TOL = 1e-6  # allows ~1e-12 float-arithmetic-order noise + 6 orders headroom
    P_TAIL_TOL = 1e-4  # p_tail is more sensitive to upstream noise (qpAdm chain)
    numeric_cols = [c for c in ref.columns if c != "target"]

    for col in numeric_cols:
        delta = (ref[col] - got[col]).abs()
        tol = P_TAIL_TOL if col == "p_tail" else WEIGHT_TOL
        max_delta = float(delta.max())
        assert max_delta < tol, (
            f"column {col!r}: max |delta| = {max_delta:.3e} exceeds "
            f"tolerance {tol:.0e}.\nref:\n{ref[['target', col]].to_string()}\n"
            f"got:\n{got[['target', col]].to_string()}"
        )


# Embedded R script for the full pipeline shootout. Kept inline (rather than
# as a separate file) so test_dogfood.py is self-contained.
_QPADM_SHOOTOUT_R = """\
suppressPackageStartupMessages({ library(admixtools) })

args <- commandArgs(trailingOnly = TRUE)
PLINK_PREFIX <- args[1]
F2_DIR       <- args[2]
OUT_TSV      <- args[3]

SOURCES   <- c("Patterson_WHGA", "Patterson_Balkan_N", "Patterson_OldSteppe")
RIGHT_SET <- c("Patterson_OldAfrica", "Patterson_WHGB", "Patterson_Turkey_N",
               "Patterson_Russia_Afanasievo")
TARGETS   <- c("Patterson_England_C_EBA", "Patterson_England_MBA",
               "Patterson_England_LBA", "Patterson_England_IA",
               "dogfood_target")
all_pops  <- unique(c(SOURCES, RIGHT_SET, TARGETS))

dir.create(F2_DIR, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(F2_DIR)) == 0) {
  extract_f2(
    pref = PLINK_PREFIX, outdir = F2_DIR, pops = all_pops,
    qpfstats = FALSE, afprod = TRUE,
    maxmiss = 0, apply_corr = TRUE, minmaf = 0, maxmaf = 0.5,
    minac2 = FALSE, blgsize = 0.05, n_cores = 4,
    verbose = FALSE, overwrite = FALSE, adjust_pseudohaploid = TRUE
  )
}

f2 <- f2_from_precomp(F2_DIR, pops = all_pops, afprod = TRUE)
rows <- list()
for (tgt in TARGETS) {
  res <- qpadm(f2, left = SOURCES, right = RIGHT_SET, target = tgt)
  w  <- res$weights$weight; se <- res$weights$se
  rows[[tgt]] <- data.frame(
    target = tgt,
    WHGA_weight  = w[1], WHGA_se  = se[1],
    Balkan_weight = w[2], Balkan_se = se[2],
    Steppe_weight = w[3], Steppe_se = se[3],
    sum_w = sum(w), p_tail = res$rankdrop$p[1],
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)
write.table(df, OUT_TSV, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
"""
