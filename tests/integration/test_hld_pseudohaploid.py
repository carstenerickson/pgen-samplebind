"""HLD test 9: pseudohaploid detection.

Per HLD §Validation strategy / LLD §5.3:

  9. Pseudohaploid detection: synthetic panel with 50% pseudohaploid +
     50% diploid samples; output `.psam` PSEUDOHAPLOID column matches
     synthesis ground truth.

The synthesizer writes the PSEUDOHAPLOID column to its output .psam as
the ground-truth label (what was requested). The merge orchestrator
ignores that column on input and re-derives it from genotypes via
pseudohaploid.update_block + classify_all (no information leaked between
input and output PSEUDOHAPLOID columns). So comparing input ground-truth
→ output detected is a true validation of the detection logic.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _read_psam_pseudohaploid_column(prefix: Path) -> list[str]:
    """Read PSEUDOHAPLOID column values from a .psam file, in psam order."""
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df["PSEUDOHAPLOID"].astype(str).tolist()


def _make_panel(
    tmp_path: Path,
    name: str,
    sample_id_prefix: str,
    sample_seed: int,
    pseudohaploid_fraction: float,
    n_samples: int = 20,
    n_variants: int = 200,
) -> Path:
    """Synthesize a panel with the given pseudohaploid_fraction."""
    spec = SyntheticPanelSpec(
        n_samples=n_samples,
        n_variants=n_variants,
        n_populations=2,
        variant_seed=901,
        sample_seed=sample_seed,
        sample_id_prefix=sample_id_prefix,
        pseudohaploid_fraction=pseudohaploid_fraction,
    )
    return synthesize_pfile(spec, tmp_path / name).path


class TestHld09PseudohaploidDetection:
    def test_50_50_mix_matches_ground_truth(self, tmp_path: Path) -> None:
        """Default 50/50 pseudohaploid_fraction → output PSEUDOHAPLOID column
        matches the synthesizer's ground truth per-sample."""
        a = _make_panel(tmp_path, "a", "A", sample_seed=10, pseudohaploid_fraction=0.5)
        b = _make_panel(tmp_path, "b", "B", sample_seed=20, pseudohaploid_fraction=0.5)

        gt_a = _read_psam_pseudohaploid_column(a)
        gt_b = _read_psam_pseudohaploid_column(b)
        ground_truth = gt_a + gt_b

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), str(b), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid_column(out)
        assert detected == ground_truth

        # Sanity: roughly 50/50 split. Bernoulli variance for 40 samples:
        # n_pseudo expected ~20 ± ~3. Allow a generous window.
        n_pseudo = ground_truth.count("1")
        n_diploid = ground_truth.count("0")
        assert n_pseudo + n_diploid == 40, "all samples should classify cleanly"
        assert 12 <= n_pseudo <= 28, f"50/50 split badly skewed: {n_pseudo}/40"

    def test_all_pseudohaploid_classifies_pseudohaploid(self, tmp_path: Path) -> None:
        """pseudohaploid_fraction=1.0 → all samples should classify as 1
        (no hets in synthesis → het_count == 0 → PSEUDOHAPLOID)."""
        a = _make_panel(tmp_path, "a", "A", sample_seed=11, pseudohaploid_fraction=1.0)
        b = _make_panel(tmp_path, "b", "B", sample_seed=21, pseudohaploid_fraction=1.0)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), str(b), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid_column(out)
        assert all(s == "1" for s in detected), f"expected all PSEUDOHAPLOID; got {set(detected)}"

    def test_all_diploid_classifies_diploid(self, tmp_path: Path) -> None:
        """pseudohaploid_fraction=0.0 → all samples should classify as 0
        (uniform 0/1/2 with ~33% het rate >> 5% threshold → DIPLOID)."""
        a = _make_panel(tmp_path, "a", "A", sample_seed=12, pseudohaploid_fraction=0.0)
        b = _make_panel(tmp_path, "b", "B", sample_seed=22, pseudohaploid_fraction=0.0)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), str(b), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid_column(out)
        # With 200 variants and ~33% expected het rate, the chance of any
        # diploid sample getting 0 hets is essentially zero.
        assert all(s == "0" for s in detected), f"expected all DIPLOID; got {set(detected)}"

    @pytest.mark.parametrize("fraction", [0.25, 0.75])
    def test_arbitrary_fraction_matches_ground_truth(self, tmp_path: Path, fraction: float) -> None:
        """Sanity check with non-50/50 fractions: detected still matches
        ground truth exactly."""
        a = _make_panel(
            tmp_path,
            "a",
            "A",
            sample_seed=int(13 + fraction * 100),
            pseudohaploid_fraction=fraction,
        )
        b = _make_panel(
            tmp_path,
            "b",
            "B",
            sample_seed=int(23 + fraction * 100),
            pseudohaploid_fraction=fraction,
        )

        ground_truth = _read_psam_pseudohaploid_column(a) + _read_psam_pseudohaploid_column(b)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), str(b), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid_column(out)
        assert detected == ground_truth
