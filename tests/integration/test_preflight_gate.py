"""Behavioral tests for the preflight gate (issue [#12](https://github.com/carstenerickson/pgen-samplebind/issues/12) step 4).

`evaluate_gate` is pure; the CLI wrapper does the output / exception side
effects. These tests cover both surfaces against the failure-mode corpus.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from pgen_samplebind.errors import ValidationError
from pgen_samplebind.formats import prepared_input
from pgen_samplebind.preflight import (
    GATE_FAILURE_LABELS,
    compute_preflight,
    evaluate_gate,
)
from pgen_samplebind.types import MergePolicy
from tests.fixtures.preflight_corpus import ALL_BUILDERS, CorpusPair


@pytest.fixture(scope="session")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, CorpusPair]:
    out: dict[str, CorpusPair] = {}
    for label, builder in ALL_BUILDERS.items():
        sub = tmp_path_factory.mktemp(f"gate_{label}")
        out[label] = builder(sub)
    return out


def _compute(pair: CorpusPair, *, preflight_policy: str, variant_key: str = "chr_pos"):
    with (
        prepared_input(pair.canonical, is_target=False, include_chrom=tuple(range(1, 23))) as a,
        prepared_input(pair.other, is_target=False, include_chrom=tuple(range(1, 23))) as b,
    ):
        policy = MergePolicy(
            variant_key=variant_key,  # type: ignore[arg-type]
            preflight_policy=preflight_policy,  # type: ignore[arg-type]
        )
        report = compute_preflight([a, b], policy, tool_version="test", command="merge")
        return evaluate_gate(report, policy), policy


# ---------------------------------------------------------------------------
# Pure-function tests
# ---------------------------------------------------------------------------


def test_gate_compatible_never_triggers(corpus: dict[str, CorpusPair]) -> None:
    report, _ = _compute(corpus["compatible"], preflight_policy="strict")
    assert report.gate["triggered"] is False
    assert report.gate["action"] == "none"


@pytest.mark.parametrize("label", ["build_mismatch", "disjoint"])
def test_gate_warn_emits_action_warn(corpus: dict[str, CorpusPair], label: str) -> None:
    """Failure-mode pairs under the default `warn` policy: triggered + warn."""
    report, _ = _compute(corpus[label], preflight_policy="warn")
    assert report.gate["triggered"] is True
    assert report.gate["action"] == "warn"
    assert report.gate["policy"] == "warn"
    assert len(report.gate["failing_inputs"]) == 1
    assert report.gate["failing_inputs"][0]["classification"] in GATE_FAILURE_LABELS


@pytest.mark.parametrize("label", ["build_mismatch", "disjoint"])
def test_gate_strict_emits_action_error(corpus: dict[str, CorpusPair], label: str) -> None:
    """Same pairs under `strict`: triggered + error."""
    report, _ = _compute(corpus[label], preflight_policy="strict")
    assert report.gate["triggered"] is True
    assert report.gate["action"] == "error"


@pytest.mark.parametrize("label", ["build_mismatch", "disjoint"])
def test_gate_off_never_acts(corpus: dict[str, CorpusPair], label: str) -> None:
    """`off` suppresses the gate even on failing classifications.

    `triggered` reflects whether the policy ACTED (false here); the
    classification-level signal lives in `would_trigger` and stays True
    so JSON consumers can still see what would have fired. See review
    finding #5.
    """
    report, _ = _compute(corpus[label], preflight_policy="off")
    assert report.gate["triggered"] is False  # policy did not act
    assert report.gate["would_trigger"] is True  # classification failed
    assert report.gate["action"] == "none"
    assert report.gate["policy"] == "off"
    # Diagnostics still populated under `off` so users can see why.
    assert len(report.gate["failing_inputs"]) >= 1


def test_gate_key_space_mismatch_triggers_under_wrong_key(
    corpus: dict[str, CorpusPair],
) -> None:
    """The key_space_mismatch fixture trips the gate only when variant_key=id
    (the user's mistake-shape); chr_pos sees it as compatible."""
    report_cp, _ = _compute(corpus["key_space_mismatch"], preflight_policy="warn")
    assert report_cp.gate["triggered"] is False

    report_id, _ = _compute(corpus["key_space_mismatch"], preflight_policy="warn", variant_key="id")
    assert report_id.gate["triggered"] is True
    assert report_id.gate["failing_inputs"][0]["classification"] == "key_space_mismatch"


# ---------------------------------------------------------------------------
# CLI integration: full merge invocation under each policy
# ---------------------------------------------------------------------------


def _run_merge_cli(canonical: Path, other: Path, out_prefix: Path, *, policy: str):
    """Invoke `pgen-samplebind merge` via CliRunner. Returns the full `result`.

    CliRunner merges stdout+stderr into `result.output`. Raised
    PgenSamplebindError instances surface as `result.exception` (since
    CliRunner doesn't go through `cli.main`'s exit-code routing).
    """
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "merge",
            str(canonical),
            str(other),
            "-o",
            str(out_prefix),
            "--preflight-policy",
            policy,
            "--on-collision",
            "first",  # corpus pairs share sample-id prefixes; sidestep collisions
        ],
        # catch_exceptions=True (default): PgenSamplebindError is captured as
        # result.exception. cli.main()'s exit-code routing isn't exercised here
        # (CliRunner bypasses it); tests/integration/test_exit_codes.py covers
        # the exit-code mapping via real subprocess invocation.
    )
    return result


def test_merge_cli_strict_unlinks_stale_triplet(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """Regression for review finding #4: a strict-policy failure must
    unlink any stale `.pgen`/`.pvar`/`.psam` left from a prior successful
    run at the same prefix, so workflow managers don't consume outputs
    that no longer match the current preflight report.
    """
    out_prefix = tmp_path / "out"
    # Plant a stale triplet from a previous successful run.
    Path(str(out_prefix) + ".pgen").write_bytes(b"stale-pgen")
    Path(str(out_prefix) + ".pvar").write_bytes(b"stale-pvar")
    Path(str(out_prefix) + ".psam").write_bytes(b"stale-psam")

    pair = corpus["build_mismatch"]
    result = _run_merge_cli(pair.canonical, pair.other, out_prefix, policy="strict")
    assert isinstance(result.exception, ValidationError)
    # Strict failure cleans up the stale triplet so on-disk state is
    # coherent: failing preflight.json + no genotype outputs.
    assert not Path(str(out_prefix) + ".pgen").exists()
    assert not Path(str(out_prefix) + ".pvar").exists()
    assert not Path(str(out_prefix) + ".psam").exists()
    assert (tmp_path / "out.preflight.json").exists()


def test_merge_cli_strict_raises_validation_error(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """End-to-end: `--preflight-policy strict` against build_mismatch
    panels raises ValidationError (→ exit 1 in `cli.main`) and the
    preflight JSON is still written so the user can inspect why.

    See `tests/integration/test_exit_codes.py` for the exit-code mapping
    test (subprocess-based); this test scopes to the in-process
    behavior CliRunner can observe.
    """
    pair = corpus["build_mismatch"]
    out_prefix = tmp_path / "out"
    result = _run_merge_cli(pair.canonical, pair.other, out_prefix, policy="strict")
    assert isinstance(result.exception, ValidationError)
    assert "Preflight gate triggered" in str(result.exception)
    assert (tmp_path / "out.preflight.json").exists()
    # Strict failure raises before pass-2 even starts → no .pgen written.
    assert not (tmp_path / "out.pgen").exists()


def test_merge_cli_warn_completes_with_output_warning(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """`warn` policy: output contains the warning, no exception raised."""
    pair = corpus["build_mismatch"]
    out_prefix = tmp_path / "out"
    result = _run_merge_cli(pair.canonical, pair.other, out_prefix, policy="warn")
    # Merge may still fail downstream (near-empty intersection) — what
    # matters is the gate emitted a warning and didn't raise here.
    assert "WARNING" in result.output
    assert "Preflight gate triggered" in result.output
    assert (tmp_path / "out.preflight.json").exists()


def test_merge_cli_off_silent_on_failure_mode(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """`off` policy: no preflight-gate warning in output."""
    pair = corpus["build_mismatch"]
    out_prefix = tmp_path / "out"
    result = _run_merge_cli(pair.canonical, pair.other, out_prefix, policy="off")
    assert "Preflight gate triggered" not in result.output
    assert (tmp_path / "out.preflight.json").exists()


def test_merge_cli_compatible_no_warning_under_warn(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """Negative control: compatible inputs under the default `warn` policy
    produce no preflight warning."""
    pair = corpus["compatible"]
    out_prefix = tmp_path / "out"
    result = _run_merge_cli(pair.canonical, pair.other, out_prefix, policy="warn")
    assert result.exit_code == 0
    assert "Preflight gate triggered" not in result.output
