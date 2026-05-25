"""Unit tests for alignment.build_alignment_table — pass-1 pipeline.

Exercises:
- happy path (all passthrough)
- per-input action assignment (passthrough / swap / flip)
- missing-in-other handling (FILL_MISSING / drop_variant / error)
- extras detection
- count_kept_variants
- compute_intersection_size
- build_action_histogram (8-key contract)
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from pgen_samplebind.alignment import (
    build_action_histogram,
    build_action_histogram_per_chrom,
    build_alignment_table,
    compute_intersection_size,
    count_kept_variants,
)
from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.types import AlignmentSummary, DropReason, MergeAction, MergePolicy


def _pvar_df(
    chroms: list[int], positions: list[int], ids: list[str], refs: list[str], alts: list[str]
) -> pd.DataFrame:
    """Build a pvar DataFrame in the canonical schema produced by pvar.read_pvar.

    Includes `_pgen_row` (uint32, pgenlib's native variant_idx type — the
    original .pgen row position read_pvar stamps before its biallelic-SNP
    filter). These synthetic inputs have no pre-filter rows, so `_pgen_row`
    matches the row index.
    """
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


class TestHappyPath:
    def test_all_passthrough_two_inputs(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df(
            [1, 1, 2], [100, 200, 50], ["a", "b", "c"], ["A", "C", "G"], ["G", "T", "T"]
        )
        other = canonical.copy()
        table = build_alignment_table(canonical, [other], policy, summary)

        assert len(table) == 3
        assert (table["action_input_1"] == MergeAction.PASSTHROUGH.value).all()
        assert summary.n_passthrough == 3
        assert summary.n_extras_dropped == 0


class TestPerInputActions:
    def test_swap_detected(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1], [100], ["v1"], ["A"], ["G"])
        other = _pvar_df([1], [100], ["v1"], ["G"], ["A"])  # swapped
        table = build_alignment_table(canonical, [other], policy, summary)
        assert table["action_input_1"].iloc[0] == MergeAction.REF_ALT_SWAP.value
        assert summary.n_ref_alt_swap == 1

    def test_strand_flip_detected(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1], [100], ["v1"], ["A"], ["G"])
        other = _pvar_df([1], [100], ["v1"], ["T"], ["C"])  # complement
        table = build_alignment_table(canonical, [other], policy, summary)
        assert table["action_input_1"].iloc[0] == MergeAction.STRAND_FLIP.value
        assert summary.n_strand_flip == 1

    def test_ambiguous_drops_by_default(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1], [100], ["v1"], ["A"], ["T"])
        other = _pvar_df([1], [100], ["v1"], ["T"], ["A"])  # ambiguous swap
        table = build_alignment_table(canonical, [other], policy, summary)
        assert table["action_input_1"].iloc[0] == MergeAction.DROP.value
        assert table["drop_reason_input_1"].iloc[0] == DropReason.AMBIGUOUS_STRAND.value
        assert summary.n_dropped == 1
        assert summary.n_dropped_by_reason[DropReason.AMBIGUOUS_STRAND] == 1


class TestMissingInOther:
    def test_fill_missing_default(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])  # missing variant b
        table = build_alignment_table(canonical, [other], policy, summary)
        assert table["action_input_1"].iloc[0] == MergeAction.PASSTHROUGH.value
        assert table["action_input_1"].iloc[1] == MergeAction.FILL_MISSING.value
        assert summary.n_fill_missing == 1

    def test_drop_variant_policy(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        modified = replace(policy, on_missing="drop_variant")
        table = build_alignment_table(canonical, [other], modified, summary)
        assert table["action_input_1"].iloc[1] == MergeAction.DROP.value
        assert table["drop_reason_input_1"].iloc[1] == DropReason.ON_MISSING_DROP_VARIANT.value

    def test_error_policy_raises(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        modified = replace(policy, on_missing="error")
        with pytest.raises(InvariantViolation, match="--on-missing error"):
            build_alignment_table(canonical, [other], modified, summary)

    def test_error_policy_softens_in_validate_mode(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        modified = replace(policy, on_missing="error")
        table = build_alignment_table(
            canonical, [other], modified, summary, soften_policy_errors=True
        )
        # Doesn't raise; records the trigger
        assert summary.policy_error_triggers["on_missing_count"] == 1
        assert table["action_input_1"].iloc[1] == MergeAction.FILL_MISSING.value


class TestExtras:
    def test_extras_counted(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        other = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])  # extra: b
        build_alignment_table(canonical, [other], policy, summary)
        assert summary.n_extras_dropped == 1

    def test_extras_error_policy_raises(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        other = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        modified = replace(policy, on_extra="error")
        with pytest.raises(InvariantViolation, match="--on-extra error"):
            build_alignment_table(canonical, [other], modified, summary)


class TestCountKeptVariants:
    def test_no_drops(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        assert count_kept_variants(table) == 2

    def test_with_drops(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "T"], ["G", "C"])
        # First variant: A/T canonical with A/T → pass; T/A → drop ambiguous
        other = _pvar_df([1, 1], [100, 200], ["a", "b"], ["T", "A"], ["A", "G"])
        # Wait: variant a is A/T canonical, T/A other → DROP ambig
        # variant b is T/C canonical (non-ambig), A/G other = complement of T/C → STRAND_FLIP
        table = build_alignment_table(canonical, [other], policy, summary)
        assert count_kept_variants(table) == 1  # one drop, one keep


class TestComputeIntersectionSize:
    def test_full_intersection(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        assert compute_intersection_size(table) == 2

    def test_partial_intersection(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        table = build_alignment_table(canonical, [other], policy, summary)
        # variant a is in both (intersection); variant b is FILL_MISSING (not intersection)
        assert compute_intersection_size(table) == 1


class TestBuildActionHistogram:
    def test_histogram_keys_always_present(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        hist = build_action_histogram(table)
        expected_keys = {
            "passthrough",
            "swap",
            "flip",
            "fill_missing",
            "dropped_ambiguous_strand",
            "dropped_allele_mismatch",
            "dropped_on_strand",
            "pre_alignment_filter_dropped",
            "drop",
        }
        assert set(hist.keys()) == expected_keys

    def test_passthrough_counted(self, policy: MergePolicy, summary: AlignmentSummary) -> None:
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        hist = build_action_histogram(table)
        assert hist["passthrough"] == 2
        assert hist["swap"] == 0
        assert hist["flip"] == 0

    def test_drop_residual_for_on_missing_drop_variant(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """The bare 'drop' key is ON_MISSING_DROP_VARIANT residual"""
        canonical = _pvar_df([1, 1], [100, 200], ["a", "b"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        modified = replace(policy, on_missing="drop_variant")
        table = build_alignment_table(canonical, [other], modified, summary)
        hist = build_action_histogram(table)
        assert hist["drop"] == 1
        # Other dropped_* buckets should be 0 for this case
        assert hist["dropped_ambiguous_strand"] == 0
        assert hist["dropped_allele_mismatch"] == 0


class TestBuildActionHistogramPerChrom:
    """Per-chromosome 8-key breakdown (v0.2)."""

    def test_only_present_chroms_appear(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df(
            [1, 1, 6, 22], [100, 200, 300, 400], list("abcd"), list("ACAA"), list("GTGG")
        )
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        per_chrom = build_action_histogram_per_chrom(table)
        assert set(per_chrom.keys()) == {1, 6, 22}

    def test_per_chrom_sums_to_global_histogram(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """Sanity: summing each key across all chroms must reproduce the
        global `action_histogram`."""
        canonical = _pvar_df(
            [1, 1, 6, 22], [100, 200, 300, 400], list("abcd"), list("ACAA"), list("GTGG")
        )
        table = build_alignment_table(canonical, [canonical.copy()], policy, summary)
        global_hist = build_action_histogram(table)
        per_chrom = build_action_histogram_per_chrom(table)
        for key in global_hist:
            per_chrom_sum = sum(h[key] for h in per_chrom.values())
            assert per_chrom_sum == global_hist[key], f"key={key} mismatch"

    def test_concentrates_drops_on_correct_chrom(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        """A drop on chr 6 only must surface in per_chrom[6], not on other chroms.
        This is the motivating diagnostic: HLA-region drops localized to chr 6
        shouldn't average out into the global histogram."""
        # Canonical: chr 1 + chr 6 each have a variant. Other input has only chr 1.
        # With on_missing=drop_variant, chr 6 variant gets dropped on the OTHER side.
        canonical = _pvar_df([1, 6], [100, 200], ["a", "b"], ["A", "A"], ["G", "G"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        modified = replace(policy, on_missing="drop_variant")
        table = build_alignment_table(canonical, [other], modified, summary)
        per_chrom = build_action_histogram_per_chrom(table)
        assert per_chrom[1]["passthrough"] == 1
        assert per_chrom[1]["drop"] == 0
        assert per_chrom[6]["drop"] == 1
        assert per_chrom[6]["passthrough"] == 0

    def test_empty_table_returns_empty_dict(self) -> None:
        # Defensive: validate_cmd may produce an empty table on pre-filter failure.
        empty = pd.DataFrame({"chrom": pd.Series(dtype="int8"), "pos": pd.Series(dtype="int64")})
        assert build_action_histogram_per_chrom(empty) == {}


class TestUniqueKeyValidation:
    def test_duplicate_canonical_raises(
        self, policy: MergePolicy, summary: AlignmentSummary
    ) -> None:
        canonical = _pvar_df([1, 1], [100, 100], ["a", "a2"], ["A", "C"], ["G", "T"])
        other = _pvar_df([1], [100], ["a"], ["A"], ["G"])
        with pytest.raises(InvariantViolation, match="duplicate"):
            build_alignment_table(canonical, [other], policy, summary)
