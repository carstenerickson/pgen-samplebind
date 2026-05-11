"""End-to-end exit-code harness via real subprocess invocation of cli.main().

Per LLD §3.16, exit codes are stable across versions:
    OK = 0
    VALIDATION_FAILURE = 1
    IO_FAILURE = 2
    INVARIANT_VIOLATION = 3
    USAGE_ERROR = 4

Day 3 noted that click's CliRunner doesn't go through cli.main(), so its
exit-code mapping isn't testable via runner.invoke. This file fills that
gap by spawning a real `python -m pgen_samplebind ...` subprocess for each
ExitCode case so the actual code-path that ships in the wheel is what's
under test.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from tests.fixtures.modifiers import drop_variants, subset_variants
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _run_cli(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    """Invoke `python -m pgen_samplebind ...` and capture exit code + streams."""
    return subprocess.run(
        [sys.executable, "-m", "pgen_samplebind", *args],
        capture_output=True,
        text=True,
        cwd=cwd,
        timeout=60,
    )


@pytest.fixture(scope="module")
def panel_a(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("exit_panel_a")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=80,
        n_populations=2,
        seed=0xE00C0DE,
        sample_id_prefix="A_",
    )
    desc = synthesize_pfile(spec, out / "A")
    return Path(str(desc.pgen_path)[:-5])  # strip .pgen


@pytest.fixture(scope="module")
def panel_b(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Same variants/seed as panel_a but disjoint sample IDs (B_*) so merge
    is a clean cross-panel bind."""
    out = tmp_path_factory.mktemp("exit_panel_b")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=80,
        n_populations=2,
        seed=0xE00C0DE,  # same variants
        sample_id_prefix="B_",
        sample_seed=0xB,
    )
    desc = synthesize_pfile(spec, out / "B")
    return Path(str(desc.pgen_path)[:-5])


def test_exit_0_happy_path(panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
    """Successful merge → exit 0."""
    out = tmp_path / "merged"
    result = _run_cli(["merge", str(panel_a), str(panel_b), "-o", str(out), "--quiet"])
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert out.with_suffix(".pgen").exists()


def test_exit_1_validation_failure_extras_threshold(panel_a: Path, tmp_path: Path) -> None:
    """Extras > 5% with --on-extra warn fires gate (a) → exit 1."""
    # Subset panel_a's variants down to 40 to use as canonical, panel_a stays
    # at 80 → 50% extras vs canonical → gate (a) fires.
    import numpy as np

    canonical_prefix = tmp_path / "canonical"
    subset_variants(panel_a, 10, 80, np.arange(40), canonical_prefix)
    out = tmp_path / "merged"
    # canonical_prefix shares panel_a's IIDs (subset copied .psam verbatim);
    # use --on-collision suffix so we hit gate (a) on extras instead of the
    # earlier collision exit-3 path.
    result = _run_cli(
        [
            "merge",
            str(canonical_prefix),
            str(panel_a),
            "-o",
            str(out),
            "--on-extra",
            "warn",
            "--on-collision",
            "suffix",
            "--quiet",
        ]
    )
    assert result.returncode == 1, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert "error:" in result.stderr.lower()


def test_nonexistent_input_is_usage_error(tmp_path: Path) -> None:
    """A missing prefix isn't an I/O failure (we never opened a file); it's a
    usage problem. Verifies the formats.detect_format → UsageError → exit 4
    path that ships in the wheel."""
    out = tmp_path / "merged"
    result = _run_cli(
        [
            "merge",
            str(tmp_path / "does_not_exist"),
            str(tmp_path / "also_missing"),
            "-o",
            str(out),
            "--quiet",
        ]
    )
    assert result.returncode == 4, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert "no PFILE/BFILE/EIGENSTRAT triplet" in result.stderr


def test_exit_2_io_failure_lock_held(panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
    """A pre-existing held lock at the output prefix → exit 2 with a useful
    message. Hold the lock from a child process so the merge subprocess
    actually contends with it."""
    import multiprocessing as mp
    import os

    out = tmp_path / "merged"

    def _child_holder(prefix_str: str, started_evt, release_evt) -> None:
        from pgen_samplebind.concurrency import output_lock

        with output_lock(Path(prefix_str)):
            started_evt.set()
            release_evt.wait(timeout=30.0)

    ctx = mp.get_context("fork" if os.name != "nt" else "spawn")
    started = ctx.Event()
    release = ctx.Event()
    proc = ctx.Process(target=_child_holder, args=(str(out), started, release))
    proc.start()
    try:
        assert started.wait(timeout=5.0), "child holder failed to start"
        result = _run_cli(["merge", str(panel_a), str(panel_b), "-o", str(out), "--quiet"])
    finally:
        release.set()
        proc.join(timeout=5.0)
        if proc.is_alive():
            proc.terminate()
            proc.join(timeout=2.0)

    assert result.returncode == 2, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert "lock" in result.stderr.lower()


def test_exit_3_invariant_violation_on_collision_error(panel_a: Path, tmp_path: Path) -> None:
    """Two copies of the same panel (same IIDs) with --on-collision error
    raises InvariantViolation → exit 3."""
    out = tmp_path / "merged"
    result = _run_cli(
        [
            "merge",
            str(panel_a),
            str(panel_a),
            "-o",
            str(out),
            "--on-collision",
            "error",
            "--quiet",
        ]
    )
    assert result.returncode == 3, f"stdout={result.stdout!r} stderr={result.stderr!r}"


def test_exit_4_usage_error_missing_required_arg() -> None:
    """Missing required --out triggers click usage error → exit 4."""
    result = _run_cli(["merge", "panel"])
    assert result.returncode == 4, f"stdout={result.stdout!r} stderr={result.stderr!r}"


def test_exit_4_usage_error_target_without_panel(tmp_path: Path, panel_a: Path) -> None:
    """--target without any positional INPUT → UsageError → exit 4."""
    out = tmp_path / "merged"
    result = _run_cli(["merge", "--target", str(panel_a), "-o", str(out), "--quiet"])
    # `merge` requires nargs=-1 with required=True positionals; click rejects
    # zero positionals as a usage error before our UsageError fires.
    assert result.returncode == 4, f"stdout={result.stdout!r} stderr={result.stderr!r}"


def test_validate_exit_0_alignment_ok(panel_a: Path, panel_b: Path) -> None:
    """validate on aligned inputs → exit 0."""
    result = _run_cli(["validate", str(panel_a), str(panel_b), "--quiet"])
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"


def test_validate_exit_1_softened_policy_error(
    panel_a: Path, panel_b: Path, tmp_path: Path
) -> None:
    """validate with --on-extra error and asymmetric extras → soft-recorded
    by build_alignment_table; gate (d) fires in evaluate_pass1_gates → exit 1.
    """
    # Drop 30 variants from panel_b → asymmetric (panel_a has 30 extras).
    import numpy as np

    canonical = tmp_path / "canonical"
    drop_variants(panel_b, 10, 80, np.arange(30), canonical)
    result = _run_cli(["validate", str(canonical), str(panel_a), "--on-extra", "error", "--quiet"])
    assert result.returncode == 1, f"stdout={result.stdout!r} stderr={result.stderr!r}"
