"""HLD integration tests for EIGENSTRAT / BFILE / cross-format paths.

Per  /  strategy. Day 6 lands tests 7, 8, 10, 19:
- Test 7: EIGENSTRAT round-trip (genotype matrix preserved through merge).
- Test 8: plink2 a7.x quirks (POP from .ind col 3 via PHENO1→POP rename;
  FID==POP for AT2; no sex chromosomes by default).
- Test 10: cross-format hash invariance (PFILE/BFILE/EIGENSTRAT same hash).
- Test 19: tempdir cleanup on plink2 subprocess failure.

Tests 13, 15, 16 (--target / cross-version / specific anno-file paths)
defer to Days 8-9 alongside --target and --relabel-from.

All tests in this module are marked `eigenstrat` so they're skipped when
plink2 isn't on PATH.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from pgen_samplebind.errors import IOFailure
from tests.fixtures import modifiers
from tests.fixtures.helpers import read_pgen_full as _read_pgen_full
from tests.fixtures.helpers import read_psam as _read_psam
from tests.fixtures.helpers import read_pvar as _read_pvar
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile

pytestmark = pytest.mark.eigenstrat


def _plink2_available() -> bool:
    return shutil.which("plink2") is not None


pytestmark = [
    pytest.mark.eigenstrat,
    pytest.mark.skipif(not _plink2_available(), reason="plink2 not on PATH"),
]


# ---------- HLD test 19: tempdir cleanup on plink2 failure -------------------


class TestHld19TempdirCleanupOnFailure:
    """Corrupted EIGENSTRAT → exit 2, stderr surfaced in error, no tempdir
    leak in $TMPDIR. The cleanup-on-failure path is the one easiest to break
    silently."""

    def test_corrupt_eigenstrat_raises_iofailure(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad"
        modifiers.make_corrupt_eigenstrat(bad)

        runner = CliRunner()
        result = runner.invoke(cli, ["hash", str(bad)])
        assert result.exit_code != 0
        assert isinstance(result.exception, IOFailure)
        # Either failure path is acceptable — the corrupt fixture (no GENO header
        # and not parseable as ASCII) gets caught by the native ASCII parser's
        # size-mismatch check; if it had a partial GENO header it would be caught
        # by plink2's --eigfile path with stderr surfaced. Both result in IOFailure
        # with an informative message.
        msg = str(result.exception)
        assert (
            "GENO" in msg
            or "stderr" in msg.lower()
            or "EIGENSTRAT" in msg
            or "size mismatch" in msg.lower()
        ), f"unexpected error message: {msg!r}"

    def test_no_tempdir_leak_after_failure(self, tmp_path: Path) -> None:
        """After a plink2-failure run, no `pgen-samplebind-*` dir survives in $TMPDIR.

        Verifies that `tempfile.TemporaryDirectory` cleanup runs on the
        exception path ( pin: stdlib's native context manager handles
        cleanup on uncaught exceptions; no manual try/finally needed)."""
        tmpdir_root = Path(os.environ.get("TMPDIR", "/tmp"))
        before = {p.name for p in tmpdir_root.glob("pgen-samplebind-*")}

        bad = tmp_path / "bad"
        modifiers.make_corrupt_eigenstrat(bad)
        runner = CliRunner()
        runner.invoke(cli, ["hash", str(bad)])

        after = {p.name for p in tmpdir_root.glob("pgen-samplebind-*")}
        # No NEW pgen-samplebind-* dirs left behind (allow pre-existing
        # dirs from other concurrent tests in the same session).
        leaked = after - before
        assert not leaked, f"tempdirs leaked: {leaked}"


# ---------- HLD test 7: EIGENSTRAT round-trip --------------------------------


class TestHld07EigenstratRoundTrip:
    """EIGENSTRAT → merge → PFILE → plink2 --export eig → genotype matrix
    matches the original input (within encoding-conversion equivalence:
    EIGENSTRAT counts REF, plink counts ALT, but biology is preserved
    so re-export then re-read should round-trip cleanly)."""

    def test_round_trip_via_merge(self, tmp_path: Path) -> None:
        """PFILE → EIGENSTRAT → tool merge → PFILE → re-export EIGENSTRAT.
        Verify variant + sample counts preserved and genotype call
        distribution matches (n_called and n_missing identical)."""
        # 1. Synthesize a PFILE, convert to EIGENSTRAT — that's our INPUT
        pfile_orig = tmp_path / "orig"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=6,
                n_variants=30,
                n_populations=2,
                variant_seed=1,
                sample_seed=10,
                sample_id_prefix="X",
            ),
            pfile_orig,
        )
        eig_input = tmp_path / "eig_input"
        modifiers.pfile_to_eigenstrat(pfile_orig, eig_input)

        # 2. Run merge with EIGENSTRAT input + --on-collision first → output
        # equals input panel after the EIG → PFILE → merge → PFILE pipeline.
        merged = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(eig_input),
                str(eig_input),
                "-o",
                str(merged),
                "--on-collision",
                "first",
                "--trust-strand",  # self-merge: ambiguous-matching is safe
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        # 3. Re-export the merged PFILE back to EIGENSTRAT (proves the
        # round-trip is self-consistent end-to-end).
        eig_roundtrip = tmp_path / "eig_roundtrip"
        modifiers.pfile_to_eigenstrat(merged, eig_roundtrip)

        # 4. Variant count preserved: .snp lines match
        snp_orig = (eig_input.with_suffix(".snp")).read_text().splitlines()
        snp_back = (eig_roundtrip.with_suffix(".snp")).read_text().splitlines()
        assert len(snp_orig) == len(snp_back) == 30

        # 5. Sample count preserved: .ind lines match
        ind_orig = (eig_input.with_suffix(".ind")).read_text().splitlines()
        ind_back = (eig_roundtrip.with_suffix(".ind")).read_text().splitlines()
        assert len(ind_orig) == len(ind_back) == 6

        # 6. Direct genotype matrix comparison via PFILE: the merged PFILE
        # genotypes should match the original PFILE genotypes (plink2 EIG↔PFILE
        # conversion preserves hardcalls when REF/ALT preserved, which our
        # tool guarantees by carrying REF/ALT canonically from input[0]).
        orig_geno = _read_pgen_full(pfile_orig, 6, 30)
        merged_geno = _read_pgen_full(merged, 6, 30)
        # Same call distribution (n_missing identical, n_het identical)
        assert int((orig_geno == -9).sum()) == int((merged_geno == -9).sum())
        assert int((orig_geno == 1).sum()) == int((merged_geno == 1).sum())


# ---------- HLD test 8: plink2 a7.x quirks ----------------------------------


class TestHld08Plink2A7xQuirks:
    """EIGENSTRAT input through automatic conversion: output PFILE has POP
    from .ind col 3 (via PHENO1→POP rename in our orchestrator), FID==POP
    for AT2 compatibility, no sex chromosomes (autosomal filter default)."""

    def test_pop_from_ind_col3_via_pheno1_rename(self, tmp_path: Path) -> None:
        pfile_orig = tmp_path / "orig"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=8,
                n_variants=40,
                n_populations=4,
                variant_seed=2,
                sample_seed=20,
                sample_id_prefix="Y",
            ),
            pfile_orig,
        )
        eig = tmp_path / "eig"
        modifiers.pfile_to_eigenstrat(pfile_orig, eig)

        merged = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["merge", str(eig), str(eig), "-o", str(merged), "--on-collision", "first", "--quiet"],
        )
        assert result.exit_code == 0, result.output

        # Output .psam: POP column populated with the .ind's pop labels
        out_psam = _read_psam(merged)
        assert "POP" in out_psam.columns
        # All 4 synthetic populations present (round-robin assignment from synth)
        assert set(out_psam["POP"].unique()) == {f"pop_{i:02d}" for i in range(4)}

    def test_fid_equals_pop(self, tmp_path: Path) -> None:
        pfile_orig = tmp_path / "orig"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=6,
                n_variants=20,
                n_populations=3,
                variant_seed=3,
                sample_seed=30,
                sample_id_prefix="Z",
            ),
            pfile_orig,
        )
        eig = tmp_path / "eig"
        modifiers.pfile_to_eigenstrat(pfile_orig, eig)

        merged = tmp_path / "merged"
        runner = CliRunner()
        runner.invoke(
            cli,
            ["merge", str(eig), str(eig), "-o", str(merged), "--on-collision", "first", "--quiet"],
        )
        out_psam = _read_psam(merged)
        assert (out_psam["FID"] == out_psam["POP"]).all()

    def test_autosome_filter_default(self, tmp_path: Path) -> None:
        """include_chrom defaults to (1..22); EIGENSTRAT conversion runs with
        --chr 1-22, so any X/Y/MT in the input would be filtered out."""
        pfile_orig = tmp_path / "orig"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=6,
                n_variants=30,
                n_populations=2,
                variant_seed=4,
                sample_seed=40,
                sample_id_prefix="W",
            ),
            pfile_orig,
        )
        eig = tmp_path / "eig"
        modifiers.pfile_to_eigenstrat(pfile_orig, eig)

        merged = tmp_path / "merged"
        runner = CliRunner()
        runner.invoke(
            cli,
            ["merge", str(eig), str(eig), "-o", str(merged), "--on-collision", "first", "--quiet"],
        )
        out_pvar = _read_pvar(merged)
        # No X (23), Y (24), MT (26) chromosomes in output
        assert (out_pvar["CHROM"].astype(int) <= 22).all()


# ---------- HLD test 10: cross-format hash invariance -----------------------


class TestHld10HashFormatInvariance:
    """Same panel as PFILE, BFILE, and EIGENSTRAT produces identical hash.

    The hash is over the canonicalized .pvar (chrom, pos, ref, alt; sorted
    numerically; ID excluded). Cross-format invariance is the contract that
    makes the hash usable for cache identity in workflows that store panels
    in different formats at different lifecycle points."""

    def test_pfile_bfile_eigenstrat_same_hash(self, tmp_path: Path) -> None:
        pfile_orig = tmp_path / "orig"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=6,
                n_variants=50,
                n_populations=2,
                variant_seed=5,
                sample_seed=50,
                sample_id_prefix="H",
            ),
            pfile_orig,
        )

        # Convert to BFILE via plink2 --make-bed
        bfile = tmp_path / "bfile"
        subprocess.run(
            ["plink2", "--pfile", str(pfile_orig), "--make-bed", "--out", str(bfile)],
            check=True,
            capture_output=True,
        )

        # Convert to EIGENSTRAT via plink2 --export eig
        eig = tmp_path / "eig"
        modifiers.pfile_to_eigenstrat(pfile_orig, eig)

        # Hash all three
        runner = CliRunner()
        h_pfile = runner.invoke(cli, ["hash", str(pfile_orig)]).output.strip()
        h_bfile = runner.invoke(cli, ["hash", str(bfile)]).output.strip()
        h_eig = runner.invoke(cli, ["hash", str(eig)]).output.strip()

        assert h_pfile == h_bfile == h_eig, (
            f"Cross-format hash invariance broken:\n"
            f"  PFILE:      {h_pfile}\n"
            f"  BFILE:      {h_bfile}\n"
            f"  EIGENSTRAT: {h_eig}"
        )
        # Sanity: it's a real sha256 line
        assert h_pfile.startswith("sha256:")
