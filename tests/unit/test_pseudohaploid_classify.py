"""Unit tests for pseudohaploid.classify cutoff boundaries.

Per  +  detection:
  het_count == 0                  → PSEUDOHAPLOID
  het_rate >= 5%                  → DIPLOID
  0 < het_rate < 5%               → UNKNOWN
  called_count == 0               → UNKNOWN  (boundary;  pin)

Denominator is `called_count` (non-missing autosomal calls), not total
sites — a sample with mostly-missing calls but observable hets in the
called subset gets classified on the called subset.
"""

from __future__ import annotations

import numpy as np
import pytest

from pgen_samplebind.pseudohaploid import classify, classify_all, init_counters, update_block
from pgen_samplebind.types import PseudohaploidStatus as P


class TestClassifyBoundaries:
    @pytest.mark.parametrize(
        "het, called, expected",
        [
            (0, 100, P.PSEUDOHAPLOID),  # zero hets → PSEUDOHAPLOID
            (0, 1, P.PSEUDOHAPLOID),  # zero hets, even with tiny called
            (5, 100, P.DIPLOID),  # 5% exactly → DIPLOID (>=)
            (10, 100, P.DIPLOID),  # well above threshold
            (50, 100, P.DIPLOID),  # very high het rate
            (1, 100, P.UNKNOWN),  # 1% → UNKNOWN
            (4, 100, P.UNKNOWN),  # 4% → UNKNOWN (just below 5%)
        ],
    )
    def test_cutoffs(self, het: int, called: int, expected: P) -> None:
        assert classify(het, called) == expected


class TestClassifyZeroCalledBoundary:
    """pin: called_count == 0 → UNKNOWN (no signal; honest answer)."""

    def test_zero_called_returns_unknown(self) -> None:
        assert classify(0, 0) == P.UNKNOWN

    def test_zero_called_does_not_divide_by_zero(self) -> None:
        # No exception even with het_count > 0 and called_count == 0
        # (which shouldn't happen in practice but defends against caller bugs).
        assert classify(5, 0) == P.UNKNOWN


class TestClassifyAllVectorized:
    def test_mixed_array_classification(self) -> None:
        het_counts = np.array([0, 0, 5, 50, 1, 4, 0], dtype=np.int64)
        called_counts = np.array([100, 1, 100, 100, 100, 100, 0], dtype=np.int64)
        statuses = classify_all(het_counts, called_counts)
        assert statuses.tolist() == [
            P.PSEUDOHAPLOID,  # 0 hets
            P.PSEUDOHAPLOID,  # 0 hets, called=1
            P.DIPLOID,  # 5% exactly
            P.DIPLOID,  # 50%
            P.UNKNOWN,  # 1%
            P.UNKNOWN,  # 4%
            P.UNKNOWN,  # called == 0 boundary
        ]


class TestUpdateBlockChromosomeFilter:
    """update_block is no-op for non-autosomal blocks (chrom > 22), per the
     detection autosome-only spec and  pin that
    blocks don't span chromosomes (caller-supplies-chromosome contract)."""

    def test_autosome_block_updates_counters(self) -> None:
        het_counts, called_counts = init_counters(n_samples=4)
        # All-het block on chrom 1 (autosome)
        block = np.array([[1, 1, 0, -9], [1, 0, 1, 1]], dtype=np.int8)
        update_block(block, chrom=1, het_counts=het_counts, called_counts=called_counts)
        # sample 0: 2 hets, 2 called
        # sample 1: 1 het, 2 called
        # sample 2: 1 het, 2 called
        # sample 3: 1 het, 1 called
        assert het_counts.tolist() == [2, 1, 1, 1]
        assert called_counts.tolist() == [2, 2, 2, 1]

    def test_sex_chromosome_block_is_noop(self) -> None:
        """chrom > 22 (X=23, Y=24, MT=26) → counters untouched."""
        het_counts, called_counts = init_counters(n_samples=4)
        block = np.array([[1, 1, 0, -9], [1, 0, 1, 1]], dtype=np.int8)
        update_block(block, chrom=23, het_counts=het_counts, called_counts=called_counts)
        assert het_counts.tolist() == [0, 0, 0, 0]
        assert called_counts.tolist() == [0, 0, 0, 0]

    def test_in_place_update_accumulates_across_blocks(self) -> None:
        het_counts, called_counts = init_counters(n_samples=2)
        block_1 = np.array([[1, 0]], dtype=np.int8)
        block_2 = np.array([[1, 1]], dtype=np.int8)
        update_block(block_1, chrom=1, het_counts=het_counts, called_counts=called_counts)
        update_block(block_2, chrom=1, het_counts=het_counts, called_counts=called_counts)
        assert het_counts.tolist() == [2, 1]
        assert called_counts.tolist() == [2, 2]
