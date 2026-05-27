"""Regression tests for `tests.fixtures.modifiers.shift_positions`.

Issue-12 step-2 review surfaced that the dict form silently defaulted
unmapped chroms to offset 0, letting test fixtures violate the
build_mismatch invariant without an obvious cause. The strict-rejection
behavior is pinned here so a future "just default it back to 0"
refactor surfaces in CI rather than as a mysterious classifier-test
flake.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tests.fixtures.modifiers import shift_positions
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


@pytest.fixture
def multi_chrom_panel(tmp_path: Path) -> Path:
    """Synthesize a panel spanning chr1 + chr2 + chr3."""
    spec = SyntheticPanelSpec(
        n_samples=4,
        n_variants=30,
        n_populations=1,
        chromosomes=(1, 2, 3),
        variant_seed=7,
        sample_seed=8,
    )
    return synthesize_pfile(spec, tmp_path / "in").path


def test_shift_positions_rejects_dict_missing_chroms(
    multi_chrom_panel: Path, tmp_path: Path
) -> None:
    """An offset dict that omits a chrom present in the pvar must raise."""
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="missing chromosomes"):
        shift_positions(multi_chrom_panel, out, per_chrom_offset={1: 1_000_000})


def test_shift_positions_dict_covering_all_chroms_works(
    multi_chrom_panel: Path, tmp_path: Path
) -> None:
    """A complete dict applies per-chrom offsets correctly."""
    out = tmp_path / "out"
    shift_positions(
        multi_chrom_panel,
        out,
        per_chrom_offset={1: 1_000_000, 2: 2_000_000, 3: 3_000_000},
    )
    df = pd.read_csv(Path(str(out) + ".pvar"), sep="\t")
    in_df = pd.read_csv(Path(str(multi_chrom_panel) + ".pvar"), sep="\t")
    # Every row's POS should be shifted by the matching chrom's offset.
    expected = {1: 1_000_000, 2: 2_000_000, 3: 3_000_000}
    for _, (c_in, p_in, c_out, p_out) in enumerate(
        zip(in_df["#CHROM"], in_df["POS"], df["#CHROM"], df["POS"], strict=True)
    ):
        assert int(c_in) == int(c_out)
        assert int(p_out) - int(p_in) == expected[int(c_in)]


def test_shift_positions_int_form_unchanged(multi_chrom_panel: Path, tmp_path: Path) -> None:
    """The int form (uniform shift) still applies to every row."""
    out = tmp_path / "out"
    shift_positions(multi_chrom_panel, out, per_chrom_offset=500_000)
    df = pd.read_csv(Path(str(out) + ".pvar"), sep="\t")
    in_df = pd.read_csv(Path(str(multi_chrom_panel) + ".pvar"), sep="\t")
    for p_in, p_out in zip(in_df["POS"], df["POS"], strict=True):
        assert int(p_out) - int(p_in) == 500_000
