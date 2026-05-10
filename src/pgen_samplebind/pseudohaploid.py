"""Per-sample heterozygosity tally during pass 2; classification at end.

Per LLD §3.8. Cutoffs match HLD §Pseudohaploid detection:
  het_count == 0 over non-missing autosomal calls → PSEUDOHAPLOID
  het_rate >= 5% → DIPLOID
  0 < het_rate < 5% → UNKNOWN
  called_count == 0 → UNKNOWN (boundary; review fix #3)

`update_block` trusts the caller-supplies-chromosome contract: blocks must
not span chromosome boundaries (merge.merge_inputs splits at boundaries).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .types import PseudohaploidStatus

_DIPLOID_HET_RATE_THRESHOLD = 0.05  # ≥ 5% het rate → diploid


def init_counters(n_samples: int) -> tuple[np.ndarray[Any, Any], np.ndarray[Any, Any]]:
    """(het_counts, called_counts) zero-init int64 arrays of shape (n_samples,)."""
    return (
        np.zeros(n_samples, dtype=np.int64),
        np.zeros(n_samples, dtype=np.int64),
    )


def update_block(
    block: np.ndarray[Any, Any],
    chrom: int,
    het_counts: np.ndarray[Any, Any],
    called_counts: np.ndarray[Any, Any],
) -> None:
    """Vectorized per-block tally update; in-place on the running counters.

    No-op for non-autosomal blocks (chrom > 22). The caller (merge.merge_inputs)
    guarantees blocks don't span chromosomes — see §3.10 chromosome-boundary pin.

    block shape: (block_variants, n_samples), int8.
    """
    if chrom > 22:
        return
    het_counts += (block == 1).sum(axis=0).astype(np.int64)
    called_counts += (block != -9).sum(axis=0).astype(np.int64)


def classify(het_count: int, called_count: int) -> PseudohaploidStatus:
    """Per-sample classification per HLD §Pseudohaploid detection.

    Boundary: called_count == 0 → UNKNOWN (no signal; honest answer).
    """
    if called_count == 0:
        return PseudohaploidStatus.UNKNOWN
    if het_count == 0:
        return PseudohaploidStatus.PSEUDOHAPLOID
    rate = het_count / called_count
    if rate >= _DIPLOID_HET_RATE_THRESHOLD:
        return PseudohaploidStatus.DIPLOID
    return PseudohaploidStatus.UNKNOWN


def classify_all(
    het_counts: np.ndarray[Any, Any], called_counts: np.ndarray[Any, Any]
) -> np.ndarray[Any, Any]:
    """Vectorized classification across all samples; returns object-dtype array."""
    return np.array(
        [
            classify(int(h), int(c))
            for h, c in zip(het_counts.tolist(), called_counts.tolist(), strict=True)
        ],
        dtype=object,
    )
