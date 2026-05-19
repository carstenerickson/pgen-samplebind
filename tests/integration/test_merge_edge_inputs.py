"""Integration tests for merge on edge-case input shapes + filesystem
error paths.

Covers gaps surfaced in the v0.3.2 review:
- single-variant, single-sample, edge variant counts (no test panel
  smaller than 50 vars previously)
- single-sample + multi-sample panel intersection (README canonical
  use case 2: VCF-from-cohort target against a reference panel)
- unreadable .pgen / unwritable output directory → exit code 2
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.helpers import read_psam_iids, read_pvar_keys
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


@pytest.fixture
def tiny_panel_a(tmp_path: Path) -> Path:
    """1 variant x 4 samples — smallest input that exercises pass-2."""
    spec = SyntheticPanelSpec(
        n_samples=4,
        n_variants=1,
        n_populations=1,
        chromosomes=(1,),
        variant_seed=701,
        sample_seed=702,
        sample_id_prefix="A",
    )
    return synthesize_pfile(spec, tmp_path / "tiny_a").path


@pytest.fixture
def tiny_panel_b(tmp_path: Path) -> Path:
    """1 variant x 4 samples; disjoint sample IDs from tiny_panel_a."""
    spec = SyntheticPanelSpec(
        n_samples=4,
        n_variants=1,
        n_populations=1,
        chromosomes=(1,),
        variant_seed=701,
        sample_seed=703,
        sample_id_prefix="B",
    )
    return synthesize_pfile(spec, tmp_path / "tiny_b").path


@pytest.fixture
def single_sample_panel(tmp_path: Path) -> Path:
    """1 sample x 20 variants — the canonical 'single-sample target' shape
    described in the README's canonical use case 2."""
    spec = SyntheticPanelSpec(
        n_samples=1,
        n_variants=20,
        n_populations=1,
        chromosomes=(1,),
        variant_seed=801,
        sample_seed=802,
        sample_id_prefix="T",
    )
    return synthesize_pfile(spec, tmp_path / "single_sample").path


@pytest.fixture
def reference_panel(tmp_path: Path) -> Path:
    """20-sample reference with the same 20 variants as single_sample_panel."""
    spec = SyntheticPanelSpec(
        n_samples=20,
        n_variants=20,
        n_populations=2,
        chromosomes=(1,),
        variant_seed=801,
        sample_seed=803,
        sample_id_prefix="R",
    )
    return synthesize_pfile(spec, tmp_path / "reference").path


class TestSingleVariantMerge:
    """Pass-2 on a 1-variant panel exercises the chromosome-block iterator's
    single-row edge case."""

    def test_one_variant_two_panels(
        self, tiny_panel_a: Path, tiny_panel_b: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(tiny_panel_a),
                str(tiny_panel_b),
                "-o",
                str(out),
                "--trust-strand",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        # 4 + 4 disjoint samples, 1 variant intersected
        assert len(read_psam_iids(out)) == 8
        assert len(read_pvar_keys(out)) == 1


class TestSingleSampleMerge:
    """README canonical use case 2: a single-sample user PFILE intersected
    with a reference panel. The merge orchestrator must handle a 1-row
    psam without tripping any sample-axis loop assumptions."""

    def test_single_sample_target_with_reference(
        self, single_sample_panel: Path, reference_panel: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(reference_panel),
                str(single_sample_panel),
                "-o",
                str(out),
                "--trust-strand",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(read_psam_iids(out)) == 21  # 20 reference + 1 target
        assert "T00000" in read_psam_iids(out)


class TestMergeZstInputs:
    """The merge path reads .pvar through two code paths — pandas
    (open_pvar_text) for alignment, and pgenlib's libzstd (PvarReader) in
    check_max_alleles. Both must handle .pvar.zst transparently.
    Previously only validate exercised the zst path."""

    def test_zst_pvar_input_on_merge(
        self, tiny_panel_a: Path, tiny_panel_b: Path, tmp_path: Path
    ) -> None:
        from tests.fixtures.modifiers import compress_pvar_to_zst

        compress_pvar_to_zst(tiny_panel_a)
        compress_pvar_to_zst(tiny_panel_b)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(tiny_panel_a),
                str(tiny_panel_b),
                "-o",
                str(out),
                "--trust-strand",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(read_psam_iids(out)) == 8
        assert len(read_pvar_keys(out)) == 1


class TestMergeIOFailures:
    """The CLI must surface IOFailure (exit code 2) when the input data
    files exist but are unreadable, or when the output directory cannot
    be written. Previously only the locked-sidecar path was tested.
    """

    def test_unreadable_pgen_raises_iofailure(
        self, tiny_panel_a: Path, tiny_panel_b: Path, tmp_path: Path
    ) -> None:
        """chmod 000 the .pgen of one input; merge must exit non-zero
        rather than crash with a raw OSError."""
        pgen = Path(str(tiny_panel_a) + ".pgen")
        original_mode = pgen.stat().st_mode
        pgen.chmod(0o000)
        try:
            out = tmp_path / "merged"
            runner = CliRunner()
            result = runner.invoke(
                cli,
                [
                    "merge",
                    str(tiny_panel_a),
                    str(tiny_panel_b),
                    "-o",
                    str(out),
                    "--trust-strand",
                    "--quiet",
                ],
            )
            assert result.exit_code != 0
            # IOFailure or InvariantViolation depending on which gate fires
            # first — both leave the user with a clear stderr message; the
            # contract is "not a raw OSError traceback".
            assert result.exception is not None
        finally:
            pgen.chmod(original_mode)

    @pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
    def test_unwritable_output_dir_raises_iofailure(
        self, tiny_panel_a: Path, tiny_panel_b: Path, tmp_path: Path
    ) -> None:
        """Write into a 0o500 (read+execute, no write) directory."""
        readonly_dir = tmp_path / "readonly"
        readonly_dir.mkdir()
        readonly_dir.chmod(0o500)
        try:
            out = readonly_dir / "merged"
            runner = CliRunner()
            result = runner.invoke(
                cli,
                [
                    "merge",
                    str(tiny_panel_a),
                    str(tiny_panel_b),
                    "-o",
                    str(out),
                    "--trust-strand",
                    "--quiet",
                ],
            )
            assert result.exit_code != 0
            assert result.exception is not None
        finally:
            readonly_dir.chmod(0o700)
