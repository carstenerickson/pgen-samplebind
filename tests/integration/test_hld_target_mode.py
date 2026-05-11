"""HLD test 13 (target+EIGENSTRAT compose) + gate (c) target call-rate tests.

Per HLD §Validation strategy / §Target mode:
- Test 13: --target user.eig + EIGENSTRAT panel positional; output AT2 f2
  matches mergeit reference within max_dev < 1e-9. (The AT2 f2 comparison
  requires AT2/mergeit binaries; only the end-to-end-runs-cleanly portion
  lands here. The f2 parity assertion is HLD test 17 territory and lives
  in a future test_mergeit_f2_parity.py marked external_tool.)
- Gate (c): target call rate below --target-min-call-rate exits 1
  (ValidationError). Denominator is canonical variant count per HLD
  §Target mode (prevents tiny-but-fully-called target from passing).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from pgen_samplebind.errors import ValidationError
from tests.fixtures import modifiers
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _plink2_available() -> bool:
    return shutil.which("plink2") is not None


def _read_psam_iids(prefix: Path) -> list[str]:
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    iid_col = next(c for c in df.columns if c.lstrip("#") == "IID")
    return df[iid_col].tolist()


# ---------- HLD test 13 (basic target + EIGENSTRAT compose) ------------------


@pytest.mark.eigenstrat
@pytest.mark.skipif(not _plink2_available(), reason="plink2 not on PATH")
class TestHld13TargetEigfileCompose:
    """Basic version: --target as EIGENSTRAT + EIGENSTRAT panel, end-to-end
    merge runs cleanly, target sample appears in output with `_target`
    suffix when colliding (or as-is otherwise). The AT2 f2 parity check
    (HLD test 13's full assertion) lives in the nightly external_tool
    test 17 path."""

    def test_target_eigfile_with_eigfile_panel_runs(self, tmp_path: Path) -> None:
        # Panel: synthesize PFILE with 6 samples x 30 variants, convert to EIGENSTRAT
        panel_pfile = tmp_path / "panel_pf"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=6,
                n_variants=30,
                n_populations=2,
                variant_seed=131,
                sample_seed=132,
                sample_id_prefix="P",
            ),
            panel_pfile,
        )
        panel_eig = tmp_path / "panel_eig"
        modifiers.pfile_to_eigenstrat(panel_pfile, panel_eig)

        # Target: synthesize PFILE with 1 sample sharing the panel's variant set,
        # convert to EIGENSTRAT. variant_seed matches so all variants align.
        target_pfile = tmp_path / "target_pf"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=1,
                n_variants=30,
                n_populations=1,
                variant_seed=131,
                sample_seed=232,
                sample_id_prefix="T",
            ),
            target_pfile,
        )
        target_eig = tmp_path / "target_eig"
        modifiers.pfile_to_eigenstrat(target_pfile, target_eig)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel_eig),
                "--target",
                str(target_eig),
                "-o",
                str(out),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        # Output has 6 panel samples + 1 target sample = 7
        iids = _read_psam_iids(out)
        assert len(iids) == 7
        # Target sample appears (no collision since prefixes differ → keeps T00000)
        assert "T00000" in iids

    def test_target_collision_gets_target_suffix(self, tmp_path: Path) -> None:
        """When target's IID collides with a panel IID, --on-collision suffix
        renames it with `_target` (per HLD test 23 case iv, also exercised
        end-to-end via the EIGENSTRAT path here)."""
        panel_pfile = tmp_path / "panel_pf"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=4,
                n_variants=20,
                n_populations=1,
                variant_seed=140,
                sample_seed=141,
                sample_id_prefix="X",
            ),
            panel_pfile,
        )
        panel_eig = tmp_path / "panel_eig"
        modifiers.pfile_to_eigenstrat(panel_pfile, panel_eig)

        # Target shares the same prefix → IID collision with panel's X00000
        target_pfile = tmp_path / "target_pf"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=1,
                n_variants=20,
                n_populations=1,
                variant_seed=140,
                sample_seed=241,
                sample_id_prefix="X",
            ),
            target_pfile,
        )
        target_eig = tmp_path / "target_eig"
        modifiers.pfile_to_eigenstrat(target_pfile, target_eig)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel_eig),
                "--target",
                str(target_eig),
                "-o",
                str(out),
                "--on-collision",
                "suffix",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        iids = _read_psam_iids(out)
        # 4 panel + 1 target = 5; target's colliding X00000 → X00000_target
        assert len(iids) == 5
        assert "X00000_target" in iids


# ---------- Gate (c) target call rate ----------------------------------------


class TestGateCTargetCallRate:
    """Gate (c) per HLD §Exit-1 validation gates (c) and §Target mode:
    target call rate below --target-min-call-rate exits 1.

    Synthesizer's `missing_rate` parameter controls per-genotype missingness.
    With missing_rate=0.95, the target's per-sample call rate is ~5%, well
    below the default 40% threshold → gate fires. With missing_rate=0.0
    (or default 0.05), call rate is ~95%, well above threshold → passes."""

    def test_low_call_rate_target_fires_gate_c(self, tmp_path: Path) -> None:
        # Panel: low missingness → high call rate, won't fail anything
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=4,
                n_variants=100,
                n_populations=1,
                variant_seed=900,
                sample_seed=901,
                sample_id_prefix="P",
                missing_rate=0.05,
            ),
            tmp_path / "panel",
        )
        # Target: 95% missingness → call rate ~5% → way below 40% default
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=1,
                n_variants=100,
                n_populations=1,
                variant_seed=900,
                sample_seed=910,
                sample_id_prefix="T",
                missing_rate=0.95,
            ),
            tmp_path / "target",
        )

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(tmp_path / "panel"),
                "--target",
                str(tmp_path / "target"),
                "-o",
                str(out),
                "--quiet",
            ],
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValidationError)
        assert "gate (c)" in str(result.exception).lower()
        assert "target call rate" in str(result.exception).lower()

    def test_low_call_rate_target_with_relaxed_threshold_passes(self, tmp_path: Path) -> None:
        """--target-min-call-rate 0.0 disables the gate; the same low-call-rate
        target merge succeeds."""
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=4,
                n_variants=100,
                n_populations=1,
                variant_seed=900,
                sample_seed=901,
                sample_id_prefix="P",
                missing_rate=0.05,
            ),
            tmp_path / "panel",
        )
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=1,
                n_variants=100,
                n_populations=1,
                variant_seed=900,
                sample_seed=910,
                sample_id_prefix="T",
                missing_rate=0.95,
            ),
            tmp_path / "target",
        )

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(tmp_path / "panel"),
                "--target",
                str(tmp_path / "target"),
                "-o",
                str(out),
                "--target-min-call-rate",
                "0.0",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

    def test_gate_c_failure_unlinks_output_triplet(self, tmp_path: Path) -> None:
        """LLD §4.1 fix #6 + §3.10 LLD pin: when gate (c) fires post-pass-2,
        the output cleanup wrapper unlinks the .pgen/.pvar/.psam triplet so
        the user doesn't see a half-built panel that's actually invalid."""
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=4,
                n_variants=100,
                n_populations=1,
                variant_seed=900,
                sample_seed=901,
                sample_id_prefix="P",
                missing_rate=0.05,
            ),
            tmp_path / "panel",
        )
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=1,
                n_variants=100,
                n_populations=1,
                variant_seed=900,
                sample_seed=910,
                sample_id_prefix="T",
                missing_rate=0.95,
            ),
            tmp_path / "target",
        )

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(tmp_path / "panel"),
                "--target",
                str(tmp_path / "target"),
                "-o",
                str(out),
                "--quiet",
            ],
        )
        assert result.exit_code != 0
        # No partial output files left on disk
        for ext in (".pgen", ".pvar", ".psam"):
            assert not (Path(str(out) + ext)).exists(), f"gate (c) failure left {out}{ext} on disk"


class TestTargetWithoutPanelRejected:
    """--target requires at least one positional INPUT (the panel).

    Click's `nargs=-1, required=True` on the positional INPUTS argument
    catches this BEFORE run_merge gets a chance to raise UsageError —
    the user gets click's "Missing argument 'INPUTS...'" message and
    exits with click's usage-error code. Either rejection path is fine
    from a UX perspective; the test just confirms the merge doesn't
    proceed with target-only input."""

    def test_target_alone_rejected(self, tmp_path: Path) -> None:
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=1,
                n_variants=10,
                n_populations=1,
                variant_seed=1,
                sample_seed=2,
                sample_id_prefix="T",
            ),
            tmp_path / "target",
        )
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["merge", "--target", str(tmp_path / "target"), "-o", str(out), "--quiet"],
        )
        assert result.exit_code != 0
        # Click rejects with SystemExit(2) (its default usage-error exit code);
        # the exact exit code mapping to ExitCode.USAGE_ERROR=4 happens in
        # cli.main() which CliRunner doesn't traverse. Either way: merge
        # doesn't run.
        assert "Missing argument" in result.output or "INPUTS" in result.output
