"""HLD test 23: IID collision suffix scheme.

Per HLD v3.5 §IID collision handling. Four cases verifying:
- General-mode `_<input_idx>` suffix (case i, ii)
- Idempotent retry through layered slots when the suffixed name itself
  collides (case iii — iterative panel-build scenario)
- Target-mode `_target` semantic suffix coexisting with general scheme
  (case iv)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.modifiers import make_panel_with_iids


def _read_psam_iids(prefix: Path) -> list[str]:
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    iid_col = next(c for c in df.columns if c.lstrip("#") == "IID")
    return df[iid_col].tolist()


class TestHld23CollisionSuffixScheme:
    def test_case_i_two_inputs_share_iid(self, tmp_path: Path) -> None:
        """(i) Two inputs sharing IID `Sample1`, no `--target`:
        output `.psam` contains `Sample1` (from input[0]) and `Sample1_1`
        (from input[1])."""
        a = make_panel_with_iids(tmp_path / "a", ["Sample1"], seed=1)
        b = make_panel_with_iids(tmp_path / "b", ["Sample1"], seed=2)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        iids = _read_psam_iids(out)
        assert iids == ["Sample1", "Sample1_1"]

    def test_case_ii_three_inputs_all_share_iid(self, tmp_path: Path) -> None:
        """(ii) Three inputs all containing `Sample1`: output contains
        `Sample1`, `Sample1_1`, `Sample1_2`."""
        a = make_panel_with_iids(tmp_path / "a", ["Sample1"], seed=1)
        b = make_panel_with_iids(tmp_path / "b", ["Sample1"], seed=2)
        c = make_panel_with_iids(tmp_path / "c", ["Sample1"], seed=3)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                str(c),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        iids = _read_psam_iids(out)
        assert iids == ["Sample1", "Sample1_1", "Sample1_2"]

    def test_case_iii_iterative_build_idempotent_retry(self, tmp_path: Path) -> None:
        """(iii) Iterative build: input[0] is a panel that already contains
        `Sample1` AND `Sample1_1` (e.g., from a previous merge that renamed
        an input[1] sample), re-merged with a new input[1] that also has IID
        `Sample1`. The new sample wants `_1` but `Sample1_1` is taken
        (by input[0]), so it falls through to `Sample1_1_1`."""
        # input[0]: panel that already has both Sample1 and Sample1_1
        a = make_panel_with_iids(tmp_path / "a", ["Sample1", "Sample1_1"], seed=1)
        # input[1]: new panel with Sample1
        b = make_panel_with_iids(tmp_path / "b", ["Sample1"], seed=2)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        iids = _read_psam_iids(out)
        # input[0] contributes Sample1 + Sample1_1 unchanged (canonical never
        # suffixed); input[1]'s Sample1 collides → tries Sample1_1 (taken) →
        # falls through to Sample1_1_1.
        assert iids == ["Sample1", "Sample1_1", "Sample1_1_1"]

    def test_case_iv_target_mode_coexists_with_general_scheme(self, tmp_path: Path) -> None:
        """(iv) Target mode: panel containing `Sample1` plus target containing
        `Sample1`: output contains `Sample1` (panel) and `Sample1_target`
        (target). Confirms the `_<input_idx>` general scheme, the idempotent-
        retry numeric suffix, and the target-mode `_target` semantic suffix
        coexist correctly."""
        panel = make_panel_with_iids(tmp_path / "panel", ["Sample1"], seed=1)
        target = make_panel_with_iids(tmp_path / "target", ["Sample1"], seed=2)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel),
                "--target",
                str(target),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                # Disable gate (c) by setting min_call_rate to 0 — these
                # synthetic samples have no call-rate-meaningful denominator.
                "--target-min-call-rate",
                "0.0",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        iids = _read_psam_iids(out)
        # panel's Sample1 stays canonical; target's Sample1 → Sample1_target
        # (target-mode suffix takes precedence over the general `_<idx>`).
        assert iids == ["Sample1", "Sample1_target"]


class TestSuffixSchemeAdditionalCases:
    """Extra coverage beyond HLD test 23's 4 explicit cases."""

    def test_target_with_existing_target_suffix_falls_through(self, tmp_path: Path) -> None:
        """If panel already has `Sample1` AND `Sample1_target`, the target's
        Sample1 wants `_target` but it's taken → falls through to
        `Sample1_target_1`."""
        panel = make_panel_with_iids(tmp_path / "panel", ["Sample1", "Sample1_target"], seed=1)
        target = make_panel_with_iids(tmp_path / "target", ["Sample1"], seed=2)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel),
                "--target",
                str(target),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                "--target-min-call-rate",
                "0.0",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        assert _read_psam_iids(out) == ["Sample1", "Sample1_target", "Sample1_target_1"]

    def test_canonical_internal_duplicate_raises(self, tmp_path: Path) -> None:
        """input[0] with internal duplicate IID under --on-collision suffix
        raises InvariantViolation (canonical must have unique IIDs)."""
        from pgen_samplebind.errors import InvariantViolation

        a = make_panel_with_iids(tmp_path / "a", ["Sample1", "Sample1"], seed=1)
        b = make_panel_with_iids(tmp_path / "b", ["X1"], seed=2)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                "--quiet",
            ],
        )
        assert isinstance(result.exception, InvariantViolation)
        assert "input[0] contains duplicate IID" in str(result.exception)
