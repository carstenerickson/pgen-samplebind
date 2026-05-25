"""Unit tests for alignment.evaluate_pass1_gates — -1 validation gates."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from pgen_samplebind.alignment import (
    build_alignment_table,
    evaluate_pass1_gates,
)
from pgen_samplebind.errors import ValidationError
from pgen_samplebind.types import AlignmentSummary, MergePolicy


def _pvar_df(
    chroms: list[int], positions: list[int], ids: list[str], refs: list[str], alts: list[str]
) -> pd.DataFrame:
    """Build a pvar DataFrame in the canonical schema produced by pvar.read_pvar
    (includes `_pgen_row` as uint32 — pgenlib's native variant_idx type —
    which build_alignment_table now reads from)."""
    df = pd.DataFrame(
        {"chrom": chroms, "pos": positions, "id": ids, "ref": refs, "alt": alts}
    ).astype({"chrom": "int8", "pos": "int64"})
    df["_pgen_row"] = np.arange(len(df), dtype=np.uint32)
    return df


@pytest.fixture
def policy() -> MergePolicy:
    return MergePolicy()


@pytest.fixture
def summary() -> AlignmentSummary:
    return AlignmentSummary()


class TestNoGateFires:
    def test_clean_alignment_passes(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        # No raise expected
        evaluate_pass1_gates(table, summary, policy, is_validate_mode=False)
        evaluate_pass1_gates(table, summary, policy, is_validate_mode=True)


class TestGateAExtrasThreshold:
    def test_extras_below_threshold_passes(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        # 100 canonical variants, 5 extras (5%) — under default 10%
        canonical = _pvar_df(
            [1] * 100, list(range(100)), [f"v{i}" for i in range(100)], ["A"] * 100, ["G"] * 100
        )
        other = _pvar_df(
            [1] * 105, list(range(105)), [f"v{i}" for i in range(105)], ["A"] * 105, ["G"] * 105
        )
        table = build_alignment_table(canonical, [other], policy, summary)
        evaluate_pass1_gates(table, summary, policy, is_validate_mode=False)

    def test_extras_above_threshold_raises(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        # 100 canonical, 20 extras (20%) — over default 10%
        canonical = _pvar_df(
            [1] * 100, list(range(100)), [f"v{i}" for i in range(100)], ["A"] * 100, ["G"] * 100
        )
        other = _pvar_df(
            [1] * 120, list(range(120)), [f"v{i}" for i in range(120)], ["A"] * 120, ["G"] * 120
        )
        table = build_alignment_table(canonical, [other], policy, summary)
        with pytest.raises(ValidationError, match="gate \\(a\\)"):
            evaluate_pass1_gates(table, summary, policy, is_validate_mode=False)


class TestGateBAmbiguousStrandThreshold:
    def test_below_intersection_threshold_passes(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """8% ambiguous-strand drops (below 10% default) → pass."""
        # 100 variants total; 8 are ambiguous AT/TA (drop); 92 pass through
        chroms = [1] * 100
        positions = list(range(1, 101))
        ids = [f"v{i}" for i in range(100)]
        canonical_refs = ["A"] * 8 + ["A"] * 92
        canonical_alts = ["T"] * 8 + ["G"] * 92
        # Other: ambiguous AT swapped (drop) for first 8, same as canonical for rest
        other_refs = ["T"] * 8 + ["A"] * 92
        other_alts = ["A"] * 8 + ["G"] * 92

        canonical = _pvar_df(chroms, positions, ids, canonical_refs, canonical_alts)
        other = _pvar_df(chroms, positions, ids, other_refs, other_alts)
        table = build_alignment_table(canonical, [other], policy, summary)
        evaluate_pass1_gates(table, summary, policy, is_validate_mode=False)

    def test_above_intersection_threshold_raises(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """25% ambiguous-strand drops (above 10% default) → ValidationError."""
        chroms = [1] * 100
        positions = list(range(1, 101))
        ids = [f"v{i}" for i in range(100)]
        canonical_refs = ["A"] * 25 + ["A"] * 75
        canonical_alts = ["T"] * 25 + ["G"] * 75
        other_refs = ["T"] * 25 + ["A"] * 75
        other_alts = ["A"] * 25 + ["G"] * 75

        canonical = _pvar_df(chroms, positions, ids, canonical_refs, canonical_alts)
        other = _pvar_df(chroms, positions, ids, other_refs, other_alts)
        table = build_alignment_table(canonical, [other], policy, summary)
        with pytest.raises(ValidationError, match="gate \\(b\\)"):
            evaluate_pass1_gates(table, summary, policy, is_validate_mode=False)

    def test_intersection_denominator_catches_wrong_panel(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """-1 validation gates (b) wrong-panel scenario.

        Small intersection x high ambiguous-drop rate = small canonical fraction.
        Intersection denominator should fire; canonical denominator wouldn't.
        """
        # 1000 canonical variants; only 10 overlap with other; 4 of the 10 are
        # ambiguous AT/TA swaps. Intersection = 10; ambig_drops = 4 → 40% > 10%.
        # Canonical fraction = 4 / 1000 = 0.4% (would silently pass canonical gate).
        chroms = [1] * 1000
        positions = list(range(1, 1001))
        ids = [f"v{i}" for i in range(1000)]
        canonical_refs = ["A"] * 1000
        canonical_alts = ["G"] * 996 + ["T"] * 4  # last 4 are A/T ambiguous
        canonical = _pvar_df(chroms, positions, ids, canonical_refs, canonical_alts)

        # Other has only 10 variants overlapping (positions 991-1000)
        other_chroms = [1] * 10
        other_positions = list(range(991, 1001))
        other_ids = [f"v{i}" for i in range(990, 1000)]
        # First 6 are passthrough (positions 991-996, A/G)
        # Last 4 are A/T canonical, swapped to T/A in other → ambiguous drop
        other_refs = ["A"] * 6 + ["T"] * 4
        other_alts = ["G"] * 6 + ["A"] * 4
        other = _pvar_df(other_chroms, other_positions, other_ids, other_refs, other_alts)

        table = build_alignment_table(canonical, [other], policy, summary)
        with pytest.raises(ValidationError, match="gate \\(b\\)"):
            evaluate_pass1_gates(table, summary, policy, is_validate_mode=False)


class TestGateDValidateModeOnly:
    def test_softened_on_missing_triggers_in_validate(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])  # missing variant b
        modified = replace(policy, on_missing="error")
        table = build_alignment_table(
            canonical, [other], modified, summary, soften_policy_errors=True
        )
        with pytest.raises(ValidationError, match="gate \\(d\\)"):
            evaluate_pass1_gates(table, summary, modified, is_validate_mode=True)

    def test_softened_on_missing_doesnt_trigger_in_merge_mode(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """Gate (d) only fires under is_validate_mode=True."""
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        modified = replace(policy, on_missing="error")
        table = build_alignment_table(
            canonical, [other], modified, summary, soften_policy_errors=True
        )
        # is_validate_mode=False → no raise
        evaluate_pass1_gates(table, summary, modified, is_validate_mode=False)

    def test_softened_on_extra_triggers(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        other = _pvar_df(
            [1, 1, 1], [100, 200, 300], ["a", "b", "c"], ["A", "A", "A"], ["G", "G", "G"]
        )
        modified = replace(policy, on_extra="error")
        # 2 extras (b, c) in 1-variant canonical = 200% → also fires gate (a),
        # but soften prevents the immediate raise during build_alignment_table.
        table = build_alignment_table(
            canonical, [other], modified, summary, soften_policy_errors=True
        )
        with pytest.raises(ValidationError):
            evaluate_pass1_gates(table, summary, modified, is_validate_mode=True)
