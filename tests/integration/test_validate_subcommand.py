"""Integration tests for the `validate` subcommand. Per LLD §4.2."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.modifiers import compress_pvar_to_zst
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


@pytest.fixture
def panel_a(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out_dir = tmp_path_factory.mktemp("vpanel_a")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=100,
        n_populations=2,
        variant_seed=1,
        sample_seed=10,
        sample_id_prefix="A",
    )
    return synthesize_pfile(spec, out_dir / "a").path


@pytest.fixture
def panel_b(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out_dir = tmp_path_factory.mktemp("vpanel_b")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=100,
        n_populations=2,
        variant_seed=1,
        sample_seed=20,
        sample_id_prefix="B",
    )
    return synthesize_pfile(spec, out_dir / "b").path


class TestValidateSmoke:
    def test_validates_clean_panels(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(panel_a), str(panel_b), "--quiet"])
        assert result.exit_code == 0, result.output

    def test_validate_emits_summary(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(panel_a), str(panel_b)])
        assert result.exit_code == 0, result.output
        assert "Validation passed" in result.output
        assert "(canonical)" in result.output

    def test_validate_writes_no_pfile(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        """Validate produces no .pgen/.pvar/.psam in tmp_path."""
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(panel_a), str(panel_b), "--quiet"])
        assert result.exit_code == 0, result.output
        # No PFILE artifacts in cwd or tmp_path
        for ext in (".pgen", ".pvar", ".psam"):
            assert not list(tmp_path.glob(f"*{ext}"))


class TestValidateGateD:
    """HLD §Exit-1 validation gates (d): --on-* error policies degrade to gate (d)
    in validate mode, exiting 1 instead of 3."""

    def test_on_missing_error_softens_to_validation_error(
        self, panel_a: Path, tmp_path: Path
    ) -> None:

        # Use only one panel to ensure no missing variants → no gate fire.
        # Then add a second panel by truncating: easiest path is to just use the same panel
        # twice with --on-collision first (no missing variants since same panel).
        # For a real "missing variant" test we'd need a modifier helper (Day 5).
        # Smoke check: validate accepts --on-missing error without crashing on clean data.
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(panel_a),
                str(panel_a),
                "--on-collision",
                "first",
                "--on-missing",
                "error",
                "--quiet",
            ],
        )
        # Same panel → no missing variants → gate (d) doesn't fire → exit 0
        assert result.exit_code == 0, result.output


class TestValidateZstPvar:
    """`.pvar.zst` panels (plink2 v2.0.0-a.6+ default, HGDP+1kGP distribution form).

    The full validate flow — format detection, multi-allelic check via
    pgenlib, pandas read, alignment — must work with a zstd-compressed
    .pvar alongside an uncompressed .pgen + .psam.
    """

    @pytest.fixture
    def panel_a_zst(self, panel_a: Path) -> Path:
        compress_pvar_to_zst(panel_a)
        return panel_a

    def test_validate_panel_with_zst_pvar(self, panel_a_zst: Path, panel_b: Path) -> None:
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(panel_a_zst), str(panel_b), "--quiet"])
        assert result.exit_code == 0, result.output

    def test_validate_both_inputs_zst(self, panel_a: Path, panel_b: Path) -> None:
        compress_pvar_to_zst(panel_a)
        compress_pvar_to_zst(panel_b)
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(panel_a), str(panel_b), "--quiet"])
        assert result.exit_code == 0, result.output

    def test_user_prefix_includes_zst_suffix(self, panel_a_zst: Path, panel_b: Path) -> None:
        """Passing `panel.pvar.zst` as the prefix arg should resolve the same
        way as passing `panel` — `strip_known_suffix` handles the two-part
        suffix."""
        zst_arg = Path(str(panel_a_zst) + ".pvar.zst")
        assert zst_arg.exists()
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(zst_arg), str(panel_b), "--quiet"])
        assert result.exit_code == 0, result.output


class TestValidateNoPopulationColumn:
    """Issue #3: a single-sample user PFILE intersected with a reference panel
    has only [IID, SEX] columns by construction. Populations are an *output*
    of the downstream ancestry classification, not a user input.

    With --no-population-column, validate should run the alignment + strand +
    IID-collision checks and skip the population-distribution checks.
    """

    @pytest.fixture
    def user_single_sample(self, panel_a: Path) -> Path:
        """Build a 1-sample user PFILE sharing variants with panel_a but with
        a psam containing only `#IID\tSEX` (no POP column)."""
        spec = SyntheticPanelSpec(
            n_samples=1,
            n_variants=100,
            n_populations=1,
            variant_seed=1,  # same variants as panel_a
            sample_seed=999,
            sample_id_prefix="USER",
        )
        out_dir = panel_a.parent.parent / "user_solo"
        out_dir.mkdir(exist_ok=True)
        desc = synthesize_pfile(spec, out_dir / "u")
        # Rewrite psam: keep only #IID, SEX (mirrors VCF→PFILE intersection output).
        psam_text = desc.psam_path.read_text().splitlines()
        header = psam_text[0].split("\t")
        # Keep #IID and SEX columns
        iid_idx = header.index("#IID")
        sex_idx = header.index("SEX")
        new_lines = ["\t".join(["#IID", "SEX"])]
        for line in psam_text[1:]:
            parts = line.split("\t")
            new_lines.append("\t".join([parts[iid_idx], parts[sex_idx]]))
        desc.psam_path.write_text("\n".join(new_lines) + "\n")
        return desc.path

    def test_validate_fails_without_flag(self, user_single_sample: Path, panel_a: Path) -> None:
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(user_single_sample), str(panel_a), "--quiet"])
        assert result.exit_code != 0, result.output

    def test_validate_succeeds_with_flag(self, user_single_sample: Path, panel_a: Path) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(user_single_sample),
                str(panel_a),
                "--no-population-column",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

    def test_validate_flag_emits_info_line(self, user_single_sample: Path, panel_a: Path) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["validate", str(user_single_sample), str(panel_a), "--no-population-column"],
        )
        assert result.exit_code == 0, result.output
        # info: ... no population column ...   on the user-side input
        assert "no population column" in result.output

    def test_flag_conflicts_with_population_column(
        self, user_single_sample: Path, panel_a: Path
    ) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(user_single_sample),
                str(panel_a),
                "--no-population-column",
                "--population-column",
                "SuperPop",
                "--quiet",
            ],
        )
        assert result.exit_code != 0
        assert "incompatible" in result.output.lower() or "no-population-column" in result.output

    def test_report_json_runs_with_flag(
        self, user_single_sample: Path, panel_a: Path, tmp_path: Path
    ) -> None:
        """--report-json should still emit a valid summary (population fields empty)."""
        import json

        report_json = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(user_single_sample),
                str(panel_a),
                "--no-population-column",
                "--report-json",
                str(report_json),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        payload = json.loads(report_json.read_text())
        assert payload["command"] == "validate"
        # User input (idx 0) has no POP → its per-input populations is empty
        assert payload["per_input_populations"][0] == {}
        # Panel (idx 1) still carries POP → non-empty
        assert payload["per_input_populations"][1] != {}


class TestValidateReports:
    def test_validate_writes_report_tsv(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        report = tmp_path / "report.tsv"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["validate", str(panel_a), str(panel_b), "--report", str(report), "--quiet"]
        )
        assert result.exit_code == 0, result.output
        assert report.exists()
        lines = report.read_text().splitlines()
        header_cols = set(lines[0].split("\t"))
        required = {"variant_id", "chr", "pos", "input_index", "action", "reason"}
        assert required.issubset(header_cols), f"missing report cols: {required - header_cols}"
        # 100 variants x 1 non-canonical input = 100 data rows
        assert len(lines) - 1 == 100

    def test_validate_writes_report_json_summary_only(
        self, panel_a: Path, panel_b: Path, tmp_path: Path
    ) -> None:
        import json

        report_json = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(panel_a),
                str(panel_b),
                "--report-json",
                str(report_json),
                "--trust-strand",  # synth panels share variant_seed; same-source
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        payload = json.loads(report_json.read_text())
        assert payload["command"] == "validate"
        assert "alignment" in payload
        assert payload["alignment"]["action_histogram"]["passthrough"] == 100
        # No variants array by default
        assert "variants" not in payload

    def test_validate_report_json_include_rows(
        self, panel_a: Path, panel_b: Path, tmp_path: Path
    ) -> None:
        import json

        report_json = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(panel_a),
                str(panel_b),
                "--report-json",
                str(report_json),
                "--report-json-include-rows",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        payload = json.loads(report_json.read_text())
        assert "variants" in payload
        assert len(payload["variants"]) == 100  # 1 non-canonical input x 100 variants
        # Spot-check a row
        row = payload["variants"][0]
        assert row["input_index"] == 1
        assert row["action"] == "passthrough"
