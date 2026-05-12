"""v0.2: tqdm-based progress bar over the pass-2 variant-block loop.

The orchestrator (`run_merge`) decides whether to enable it: on by default
when `sys.stderr.isatty()` AND `quiet=False`; off when piped or quiet so
workflow-manager stderr stays clean.

Direct-API tests here exercise the `merge_inputs` code path with
`show_progress` set both ways to confirm:
  (a) `show_progress=True` doesn't crash and completes the merge with the
      same output as the silent path.
  (b) `show_progress=False` (the existing default) is unchanged — no
      progress instrumentation in the output.

Visual rendering of the bar in a real terminal is verified manually; the
TTY-detection logic in the orchestrator is a single `sys.stderr.isatty()`
check that's exercised by CliRunner runs (no TTY) producing zero
progress output.
"""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _two_panels(tmp_path: Path) -> tuple[Path, Path]:
    a = synthesize_pfile(
        SyntheticPanelSpec(
            n_samples=6,
            n_variants=120,
            n_populations=1,
            variant_seed=1,
            sample_seed=2,
            sample_id_prefix="A",
        ),
        tmp_path / "a",
    ).path
    b = synthesize_pfile(
        SyntheticPanelSpec(
            n_samples=6,
            n_variants=120,
            n_populations=1,
            variant_seed=1,
            sample_seed=3,
            sample_id_prefix="B",
        ),
        tmp_path / "b",
    ).path
    return a, b


class TestProgressBarCodePath:
    """Verify the bar is gated correctly: --quiet always disables, piped
    stderr (no TTY) disables, identical .pgen output regardless. CliRunner
    runs under a non-TTY stderr by default, so it acts as the
    show_progress=False path."""

    def test_quiet_run_silent_and_produces_output(self, tmp_path: Path) -> None:
        a, b = _two_panels(tmp_path)
        out = tmp_path / "quiet"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["merge", str(a), str(b), "--trust-strand", "-o", str(out), "--quiet"],
        )
        assert result.exit_code == 0, result.output
        # --quiet → no progress bar AND no stdout summary.
        assert "Pass 2: streaming" not in result.output
        assert "Done:" not in result.output
        assert Path(str(out) + ".pgen").exists()

    def test_piped_run_no_progress_in_stderr(self, tmp_path: Path) -> None:
        """Non-TTY stderr (CliRunner's default) → progress bar suppressed
        even without --quiet. The stdout summary still appears."""
        a, b = _two_panels(tmp_path)
        out = tmp_path / "piped"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["merge", str(a), str(b), "--trust-strand", "-o", str(out)],
        )
        assert result.exit_code == 0, result.output
        # tqdm renders the "Pass 2: streaming genotypes" prefix; absence
        # confirms TTY detection correctly suppressed it.
        assert "Pass 2: streaming" not in result.output
        # Stdout summary IS still rendered (separate code path from progress).
        assert "Done:" in result.output

    def test_progress_path_byte_identical_output(self, tmp_path: Path) -> None:
        """Force-enable progress via the direct API and compare .pgen bytes
        against the silent path — bar is pure instrumentation, never affects
        the merged output."""
        from contextlib import ExitStack
        from dataclasses import replace

        from pgen_samplebind.formats import prepared_input
        from pgen_samplebind.merge import merge_inputs
        from pgen_samplebind.psam import (
            add_fid_from_pop,
            detect_population_column,
            read_psam,
            rename_to_pop,
            resolve_sample_identity,
        )
        from pgen_samplebind.pvar import count_raw_variants
        from pgen_samplebind.types import MergeContext, MergePolicy

        a, b = _two_panels(tmp_path)
        policy = MergePolicy(trust_strand=True, block_size=50)

        def _run(out_pgen: Path, out_pvar: Path, *, show_progress: bool) -> None:
            with ExitStack() as stack:
                descriptors = [stack.enter_context(prepared_input(p)) for p in (a, b)]
                psam_dfs = []
                for d in descriptors:
                    df = read_psam(d.psam_path)
                    df = rename_to_pop(df, detect_population_column(df, None))
                    psam_dfs.append(add_fid_from_pop(df))
                descriptors = [
                    replace(d, n_samples=len(df), n_variants=count_raw_variants(d.pvar_path))
                    for d, df in zip(descriptors, psam_dfs, strict=True)
                ]
                sample_plan = resolve_sample_identity(psam_dfs, policy)
                ctx = MergeContext(
                    policy=policy, sample_plan=sample_plan, show_progress=show_progress
                )
                merge_inputs(descriptors, out_pgen, out_pvar, ctx)

        silent_pgen = tmp_path / "silent.pgen"
        silent_pvar = tmp_path / "silent.pvar"
        progress_pgen = tmp_path / "progress.pgen"
        progress_pvar = tmp_path / "progress.pvar"
        _run(silent_pgen, silent_pvar, show_progress=False)
        _run(progress_pgen, progress_pvar, show_progress=True)

        assert silent_pgen.read_bytes() == progress_pgen.read_bytes()
        assert silent_pvar.read_bytes() == progress_pvar.read_bytes()
