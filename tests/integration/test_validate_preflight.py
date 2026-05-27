"""Behavioral tests for `validate` + preflight (issue #12 follow-up).

Validate is a dry-run for merge; the preflight gate fires under the same
classifications and the same `--preflight-policy` semantics, just without
writing a PFILE triplet. JSON emission is opt-in (`--preflight-json PATH`)
since validate has no output prefix to derive from.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.preflight_corpus import ALL_BUILDERS, CorpusPair


@pytest.fixture(scope="session")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, CorpusPair]:
    out: dict[str, CorpusPair] = {}
    for label, builder in ALL_BUILDERS.items():
        sub = tmp_path_factory.mktemp(f"validate_pf_{label}")
        out[label] = builder(sub)
    return out


def _run_validate_cli(
    canonical: Path,
    other: Path,
    *,
    policy: str = "warn",
    preflight_json: Path | None = None,
    extra_args: tuple[str, ...] = (),
):
    """Invoke `pgen-samplebind validate` via CliRunner.

    CliRunner captures PgenSamplebindError as `result.exception` rather
    than routing through cli.main()'s exit-code mapping (that path is
    covered by tests/integration/test_exit_codes.py via real subprocess
    invocation).
    """
    runner = CliRunner()
    args = [
        "validate",
        str(canonical),
        str(other),
        "--preflight-policy",
        policy,
        "--on-collision",
        "first",  # corpus pairs share sample-id prefixes
        *extra_args,
    ]
    if preflight_json is not None:
        args.extend(["--preflight-json", str(preflight_json)])
    return runner.invoke(cli, args)


def test_validate_strict_raises_on_build_mismatch(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """The whole point of validate-as-CI-dry-run: strict mode exits 1
    before the alignment table is even built."""
    pair = corpus["build_mismatch"]
    result = _run_validate_cli(pair.canonical, pair.other, policy="strict")
    assert result.exit_code != 0
    # ValidationError → exit 1 in cli.main, but CliRunner captures the
    # exception object rather than the routed code; check both surfaces.
    assert result.exception is not None
    assert "Preflight gate triggered" in str(result.exception)


def test_validate_warn_emits_stderr_continues(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """Warn-mode mirrors merge: stderr WARNING line, validate continues
    through the rest of pass 1 and exits 0 (alignment itself is fine —
    it's the input-compat picture that's bad)."""
    pair = corpus["build_mismatch"]
    result = _run_validate_cli(pair.canonical, pair.other, policy="warn")
    # CliRunner merges stdout + stderr into `result.output`.
    assert "WARNING: Preflight gate triggered" in result.output
    # Validate continues past the gate; failure (if any) is alignment-
    # level, not preflight-level. The corpus's build-mismatch pair has
    # zero coord overlap so pass 1 will then complain about extras —
    # which is the expected, NON-preflight failure mode. Exit-code
    # behavior under that downstream gate is covered elsewhere; here
    # we just check the preflight warning was emitted.


def test_validate_off_suppresses_gate(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """`off` policy: no stderr WARNING, no exception from the preflight
    gate. Downstream alignment behavior unchanged."""
    pair = corpus["build_mismatch"]
    out_json = tmp_path / "pf.json"
    result = _run_validate_cli(
        pair.canonical, pair.other, policy="off", preflight_json=out_json
    )
    assert "WARNING: Preflight gate triggered" not in result.output
    # The JSON file is still written (always-emit semantics) and its
    # `gate.would_trigger` still reflects the classification — only
    # `gate.triggered` / `action` are suppressed under `off`. See review
    # finding #5.
    assert out_json.exists()
    payload = json.loads(out_json.read_text())
    assert payload["gate"]["would_trigger"] is True
    assert payload["gate"]["triggered"] is False
    assert payload["gate"]["action"] == "none"


def test_validate_preflight_json_only_when_flagged(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """Without --preflight-json, no JSON file is written (validate has no
    `-o` to derive a default path from). The gate evaluator + stderr
    behavior runs regardless."""
    pair = corpus["compatible"]
    # No --preflight-json flag → no file should appear in tmp_path.
    files_before = set(tmp_path.iterdir())
    result = _run_validate_cli(pair.canonical, pair.other, policy="warn")
    assert result.exit_code == 0
    files_after = set(tmp_path.iterdir())
    assert files_before == files_after, (
        f"validate created unexpected files without --preflight-json: "
        f"{files_after - files_before}"
    )


def test_validate_preflight_json_schema_v1(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """The validate-emitted JSON conforms to the same schema v1 envelope
    that merge writes, with `command="validate"`."""
    pair = corpus["compatible"]
    out_json = tmp_path / "pf.json"
    result = _run_validate_cli(
        pair.canonical, pair.other, policy="warn", preflight_json=out_json
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(out_json.read_text())
    assert payload["schema_version"] == 1
    assert payload["command"] == "validate"
    assert payload["tool"] == "pgen-samplebind"
    assert "canonical" in payload
    assert payload["comparisons"][0]["classification"] == "compatible"


def test_validate_stdout_includes_preflight_block(
    corpus: dict[str, CorpusPair], tmp_path: Path
) -> None:
    """The validate-mode stdout summary surfaces the per-comparison
    classification + fraction so users see the picture without parsing
    JSON."""
    pair = corpus["compatible"]
    result = _run_validate_cli(pair.canonical, pair.other, policy="warn")
    assert result.exit_code == 0, result.output
    assert "Preflight (variant_key=chr_pos):" in result.output
    assert "classification=compatible" in result.output
    assert "gate: none" in result.output
