"""Smoke test for `inspect` subcommand using the synthetic panel fixture.

Day 1 minimal coverage. Full HLD-test mapping comes in later passes.
"""

from __future__ import annotations

import json

from click.testing import CliRunner

from pgen_samplebind.cli import cli
from pgen_samplebind.types import InputDescriptor


def test_inspect_text_smoke(synth_panel_tiny: InputDescriptor) -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", str(synth_panel_tiny.path)])
    assert result.exit_code == 0, result.output
    assert "n_samples\t20" in result.output
    assert "n_populations\t3" in result.output
    assert "format\tpfile" in result.output


def test_inspect_json_smoke(synth_panel_tiny: InputDescriptor) -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", "--json", str(synth_panel_tiny.path)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_samples"] == 20
    assert payload["n_populations"] == 3
    assert payload["format"] == "pfile"
    assert "populations" in payload
    assert sum(payload["populations"].values()) == 20


def test_inspect_default_panel_runs(synth_panel_500x50k: InputDescriptor) -> None:
    """Confirm the default 500x50K spec also works end-to-end."""
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", "--json", str(synth_panel_500x50k.path)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_samples"] == 500
    assert payload["n_populations"] == 10


def test_inspect_missingness_histogram_populated(
    synth_panel_tiny: InputDescriptor,
) -> None:
    """Histogram has all 10 bins, counts sum to n_samples."""
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", "--json", str(synth_panel_tiny.path)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)

    hist = payload["missingness_histogram"]
    assert hist is not None
    assert "bins" in hist
    assert hist["n_variants_scanned"] == payload["n_variants_pre_filter"]

    bins = hist["bins"]
    expected_labels = {
        "0-10%",
        "10-20%",
        "20-30%",
        "30-40%",
        "40-50%",
        "50-60%",
        "60-70%",
        "70-80%",
        "80-90%",
        "90-100%",
    }
    assert set(bins.keys()) == expected_labels
    assert sum(bins.values()) == payload["n_samples"]


def test_inspect_missingness_at_default_5pct_lands_in_low_bin(
    synth_panel_tiny: InputDescriptor,
) -> None:
    """Synthesizer uses missing_rate=0.05 by default; per-sample missing rate
    should cluster well below 10%, so all 20 tiny-panel samples land in the
    0-10% bin."""
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", "--json", str(synth_panel_tiny.path)])
    payload = json.loads(result.output)
    bins = payload["missingness_histogram"]["bins"]
    assert bins["0-10%"] == payload["n_samples"]
    # Sanity: no samples in any higher bucket
    for label in (
        "10-20%",
        "20-30%",
        "30-40%",
        "40-50%",
        "50-60%",
        "60-70%",
        "70-80%",
        "80-90%",
        "90-100%",
    ):
        assert bins[label] == 0


def test_inspect_text_includes_histogram_section(
    synth_panel_tiny: InputDescriptor,
) -> None:
    """Text mode renders the histogram block."""
    runner = CliRunner()
    result = runner.invoke(cli, ["inspect", str(synth_panel_tiny.path)])
    assert result.exit_code == 0, result.output
    assert "missingness_histogram" in result.output
    assert "0-10%" in result.output
