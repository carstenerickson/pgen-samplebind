"""Row-order invariant survives `python -O`.

The merge orchestrator pairs three sample-axis structures (pass-2
counters, sample_plan.output_iids, merged_psam rows). If any one
diverges, the output .pgen carries misaligned genotypes against the
output .psam — silent data corruption.

The invariant assertion in merge_cmd.run_merge used to be gated on
`if __debug__:`, which `python -O` elides. v0.3.2 review follow-up
promoted it to an unconditional raise of InvariantViolation. This test
proves the check still fires under `-O` by spawning a subprocess with
`-O` and forcing the failure mode.

We construct the failure by monkey-patching `psam.resolve_sample_identity`
to return a SampleIdentityPlan whose output_iids disagree with what the
pass-2 counters produce. The merge must reject with InvariantViolation,
not produce a silent half-corrupt output.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from textwrap import dedent


def test_row_order_invariant_fires_under_python_o(tmp_path: Path) -> None:
    """Run a tiny merge with the row-order check forced to fail, inside a
    `python -O` subprocess. The check must still raise InvariantViolation
    rather than be silently elided by `-O`."""
    repo_root = Path(__file__).resolve().parents[2]

    # Stage a synth two-panel fixture in tmp_path/setup.py, then run a
    # second script under `python -O` that performs the merge with a
    # monkey-patched resolve_sample_identity that returns a mismatched
    # output_iids list.
    script = tmp_path / "force_invariant.py"
    script.write_text(
        dedent(
            f"""
            import sys, traceback
            sys.path.insert(0, {str(repo_root)!r})

            from pathlib import Path
            from click.testing import CliRunner

            from pgen_samplebind import psam
            from pgen_samplebind.cli import cli
            from pgen_samplebind.errors import InvariantViolation
            from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile

            assert not __debug__, "expected python -O; __debug__ is True"

            tmp = Path({str(tmp_path)!r})
            a = synthesize_pfile(
                SyntheticPanelSpec(
                    n_samples=4, n_variants=10, n_populations=1,
                    chromosomes=(1,), variant_seed=1, sample_seed=2,
                    sample_id_prefix="A",
                ),
                tmp / "a",
            ).path
            b = synthesize_pfile(
                SyntheticPanelSpec(
                    n_samples=4, n_variants=10, n_populations=1,
                    chromosomes=(1,), variant_seed=1, sample_seed=3,
                    sample_id_prefix="B",
                ),
                tmp / "b",
            ).path

            # Monkey-patch psam.merge_psams to return a psam whose IID
            # row order disagrees with sample_plan.output_iids. counters
            # and sample_plan are both derived from the unmodified plan,
            # so the mismatch lives only in the merged psam — exactly the
            # silent-misalignment failure mode the invariant guards.
            from pgen_samplebind.commands import merge_cmd
            original_merge_psams = psam.merge_psams
            def broken_merge_psams(psams, sample_plan):
                df = original_merge_psams(psams, sample_plan)
                df = df.copy()
                df.loc[0, "IID"] = "__INVARIANT_VIOLATION_PROBE__"
                return df
            # The orchestrator calls psam.merge_psams via the bound
            # module attribute, so patching the source module is enough.
            psam.merge_psams = broken_merge_psams
            merge_cmd.psam.merge_psams = broken_merge_psams

            runner = CliRunner()
            result = runner.invoke(
                cli,
                ["merge", str(a), str(b), "-o", str(tmp / "out"),
                 "--trust-strand", "--quiet"],
            )
            exc = result.exception
            if not isinstance(exc, InvariantViolation):
                print(f"UNEXPECTED: exit={{result.exit_code}} exc={{exc!r}}")
                if exc is not None:
                    traceback.print_exception(type(exc), exc, exc.__traceback__)
                sys.exit(2)
            if "row-order invariant violated" not in str(exc):
                print(f"WRONG MESSAGE: {{exc}}")
                sys.exit(3)
            print("OK")
            """
        )
    )

    result = subprocess.run(
        [sys.executable, "-O", str(script)],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"subprocess exit={result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "OK" in result.stdout


def test_assertion_helper_is_available() -> None:
    """Smoke check that the merge_cmd module imports cleanly. If the
    invariant check were accidentally moved back inside `if __debug__:`
    the source still imports; this is a sanity stub paired with the
    -O subprocess test above."""
    from pgen_samplebind.commands import merge_cmd  # noqa: F401
    # No-op assertion — the import is the test.
    assert True


if __name__ == "__main__":
    # Allow `python tests/integration/test_row_order_invariant.py` for ad-hoc
    # debugging (pytest runs both functions normally).
    shutil.rmtree("/tmp/row_order_smoke", ignore_errors=True)
    Path("/tmp/row_order_smoke").mkdir()
    test_row_order_invariant_fires_under_python_o(Path("/tmp/row_order_smoke"))
    print("standalone OK")
