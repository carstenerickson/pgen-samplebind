"""Integration tests for --report and --report-json across merge and validate.

Covers HLD test 22 (Report-JSON default vs include-rows) at the integration
level. The summary-only-default and include-rows-opt-in contracts are pinned
in LLD §3.11.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


@pytest.fixture
def two_panels(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    out_dir = tmp_path_factory.mktemp("rpt")
    a = synthesize_pfile(
        SyntheticPanelSpec(
            n_samples=8,
            n_variants=50,
            n_populations=2,
            variant_seed=1,
            sample_seed=10,
            sample_id_prefix="A",
        ),
        out_dir / "a",
    ).path
    b = synthesize_pfile(
        SyntheticPanelSpec(
            n_samples=8,
            n_variants=50,
            n_populations=2,
            variant_seed=1,
            sample_seed=20,
            sample_id_prefix="B",
        ),
        out_dir / "b",
    ).path
    return a, b


class TestMergeReportTsv:
    def test_tsv_header_and_row_count(self, two_panels: tuple[Path, Path], tmp_path: Path) -> None:
        a, b = two_panels
        out = tmp_path / "merged"
        report = tmp_path / "report.tsv"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["merge", str(a), str(b), "-o", str(out), "--report", str(report), "--quiet"]
        )
        assert result.exit_code == 0, result.output
        assert report.exists()
        lines = report.read_text().splitlines()
        assert lines[0] == "variant_id\tchr\tpos\tinput_index\taction\treason"
        # 50 variants x 1 non-canonical input
        assert len(lines) == 51

    def test_tsv_actions_are_strings(self, two_panels: tuple[Path, Path], tmp_path: Path) -> None:
        a, b = two_panels
        out = tmp_path / "merged"
        report = tmp_path / "report.tsv"
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--report",
                str(report),
                "--trust-strand",  # synth panels share variant_seed; same-source
                "--quiet",
            ],
        )
        # All data rows should have action == "passthrough" (same variant set)
        rows = report.read_text().splitlines()[1:]
        for row in rows:
            fields = row.split("\t")
            assert fields[3] == "1"  # input_index
            assert fields[4] == "passthrough"


class TestMergeReportJsonSummaryOnly:
    """HLD test 22 case (i): --report-json default produces summary-only JSON."""

    def test_default_no_variants_array(self, two_panels: tuple[Path, Path], tmp_path: Path) -> None:
        a, b = two_panels
        out = tmp_path / "merged"
        report_json = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["merge", str(a), str(b), "-o", str(out), "--report-json", str(report_json), "--quiet"],
        )
        assert result.exit_code == 0, result.output

        payload = json.loads(report_json.read_text())
        assert "variants" not in payload
        # Sanity-check size: well under 100 KB for a small synth merge
        assert report_json.stat().st_size < 100 * 1024

    def test_summary_payload_shape(self, two_panels: tuple[Path, Path], tmp_path: Path) -> None:
        a, b = two_panels
        out = tmp_path / "merged"
        report_json = tmp_path / "report.json"
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(report_json),
                "--trust-strand",  # synth panels share variant_seed; same-source
                "--quiet",
            ],
        )
        payload = json.loads(report_json.read_text())

        # Top-level keys
        for key in (
            "tool",
            "tool_version",
            "command",
            "inputs",
            "policy",
            "alignment",
            "output",
            "pseudohaploid",
            "per_input_populations",
        ):
            assert key in payload, f"missing {key}"

        assert payload["command"] == "merge"
        assert len(payload["inputs"]) == 2
        assert payload["output"]["n_samples"] == 16  # 8 + 8
        assert payload["output"]["n_variants"] == 50
        assert payload["alignment"]["action_histogram"]["passthrough"] == 50

        # All 8 HLD-pinned action_histogram keys present
        expected = {
            "passthrough",
            "swap",
            "flip",
            "fill_missing",
            "dropped_ambiguous_strand",
            "dropped_allele_mismatch",
            "pre_alignment_filter_dropped",
            "drop",
        }
        assert set(payload["alignment"]["action_histogram"].keys()) == expected

        # v0.2: per-chrom 8-key breakdown also serialized. JSON keys are strings
        # (chrom ints stringified); each per-chrom value carries the same 8 keys
        # as the global histogram; per-chrom sums equal the global counts.
        per_chrom = payload["alignment"]["action_histogram_per_chrom"]
        assert isinstance(per_chrom, dict)
        assert all(k.isdigit() for k in per_chrom)  # chrom keys stringified
        assert all(set(v.keys()) == expected for v in per_chrom.values())
        for key in expected:
            per_chrom_sum = sum(v[key] for v in per_chrom.values())
            assert per_chrom_sum == payload["alignment"]["action_histogram"][key], (
                f"per_chrom sum != global for {key}"
            )


class TestMergeReportJsonIncludeRows:
    """HLD test 22 case (ii): --report-json-include-rows adds variants array."""

    def test_include_rows_adds_variants(
        self, two_panels: tuple[Path, Path], tmp_path: Path
    ) -> None:
        a, b = two_panels
        out = tmp_path / "merged"
        report_json = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(report_json),
                "--report-json-include-rows",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        payload = json.loads(report_json.read_text())
        assert "variants" in payload
        # 50 variants x 1 non-canonical input = 50 rows
        assert len(payload["variants"]) == 50
        row = payload["variants"][0]
        assert row["input_index"] == 1
        assert row["action"] == "passthrough"
        assert "variant_id" in row
        assert "chr" in row
        assert "pos" in row


class TestReportJsonSizeWarning:
    """HLD test 22 case (iii): large include-rows fires the >100 MB stderr warning.

    Constructed via stub: rather than building a 1M-variant fixture (slow),
    monkeypatch `_JSON_BYTES_PER_ROW` upward so the warning fires on a small
    synthetic merge. This verifies the warning path, not the magic number.
    """

    def test_warning_emitted_on_predicted_overflow(
        self,
        two_panels: tuple[Path, Path],
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        from pgen_samplebind import reporting

        # Force the per-row estimate way up so 50-row merge predicts > 100 MB.
        monkeypatch.setattr(reporting, "_JSON_BYTES_PER_ROW", 10_000_000)

        a, b = two_panels
        out = tmp_path / "merged"
        report_json = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(report_json),
                "--report-json-include-rows",
            ],
            # Don't pass --quiet so stderr warning is unsuppressed.
        )
        assert result.exit_code == 0, result.output
        # The warning goes to stderr; CliRunner captures both into result.output.
        assert "WARNING" in result.output
        assert "100 MB" in result.output
