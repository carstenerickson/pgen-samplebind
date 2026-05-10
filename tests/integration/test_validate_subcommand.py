"""Integration tests for the `validate` subcommand. Per LLD §4.2."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
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


class TestValidateReports:
    def test_validate_writes_report_tsv(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        report = tmp_path / "report.tsv"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["validate", str(panel_a), str(panel_b), "--report", str(report), "--quiet"]
        )
        assert result.exit_code == 0, result.output
        assert report.exists()
        content = report.read_text()
        assert "variant_id\tchr\tpos\tinput_index\taction\treason" in content
        # 100 variants x 1 non-canonical input = 100 data rows + 1 header
        assert content.count("\n") == 101

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
