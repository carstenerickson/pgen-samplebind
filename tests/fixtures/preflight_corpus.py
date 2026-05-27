"""Failure-mode corpus for the preflight classifier (issue
[#12](https://github.com/carstenerickson/pgen-samplebind/issues/12)).

Each builder returns a `(canonical_prefix, other_prefix)` pair that
reproduces one specific failure shape: build mismatch, all-allele swap,
all-strand flip, disjoint panels, key-space (chr:pos vs rsID) mismatch,
and the happy-path compatible pair as a negative control.

The corpus is the long-term asset: step 3 of issue #12 will write the
classifier against these fixtures so the per-failure-mode evidence keys
are exercised in tests, not just code-reviewed. Each pair is intentionally
tiny (40 variants, 4 samples on two chromosomes) so the full corpus
builds in well under a second.

Pairs are session-scoped via pytest fixtures in `tests/conftest.py`'s
sibling — see `tests/integration/test_preflight_corpus.py` for usage.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .modifiers import (
    flip_strand,
    rename_ids_to_rsid,
    shift_positions,
    swap_ref_alt,
)
from .synthesize import SyntheticPanelSpec, synthesize_pfile

# Single shared shape across the corpus so per-fixture differences are
# isolated to the deliberate failure-mode mutation, not panel size.
_CORPUS_SHAPE = SyntheticPanelSpec(
    n_samples=4,
    n_variants=40,
    n_populations=1,
    chromosomes=(1, 2),
    variant_seed=20260527,
    sample_seed=20260528,
    sample_id_prefix="C",
)


@dataclass(frozen=True, slots=True)
class CorpusPair:
    """One (canonical, other) pair plus the expected-failure-mode label.

    The label is the *intent* of the fixture — what the classifier (step 3)
    should ultimately return. Step 2 tests only assert measurable
    properties (intersection size, per-chrom counts); step 3 tests use the
    label as ground truth.
    """

    canonical: Path
    other: Path
    label: str  # "compatible" | "build_mismatch" | ... | "key_space_mismatch"


def _build_canonical(out_dir: Path, *, sample_seed: int | None = None) -> Path:
    """The shared canonical panel: chr:pos-style IDs, no mutations applied."""
    spec = _CORPUS_SHAPE
    if sample_seed is not None:
        spec = SyntheticPanelSpec(
            n_samples=spec.n_samples,
            n_variants=spec.n_variants,
            n_populations=spec.n_populations,
            chromosomes=spec.chromosomes,
            seed=spec.seed,
            variant_seed=spec.variant_seed,
            sample_seed=sample_seed,
            sample_id_prefix=spec.sample_id_prefix,
        )
    return synthesize_pfile(spec, out_dir / "canonical").path


def build_compatible_pair(out_dir: Path) -> CorpusPair:
    """Negative control: identical variants, disjoint samples. Intersection ≈ 100%."""
    canonical = _build_canonical(out_dir)
    # Different sample_seed → different samples; variant_seed is shared so chr:pos sets match.
    other_spec = SyntheticPanelSpec(
        n_samples=_CORPUS_SHAPE.n_samples,
        n_variants=_CORPUS_SHAPE.n_variants,
        n_populations=_CORPUS_SHAPE.n_populations,
        chromosomes=_CORPUS_SHAPE.chromosomes,
        variant_seed=_CORPUS_SHAPE.variant_seed,
        sample_seed=_CORPUS_SHAPE.sample_seed + 1,
        sample_id_prefix="D",
    )
    other = synthesize_pfile(other_spec, out_dir / "other").path
    return CorpusPair(canonical=canonical, other=other, label="compatible")


def build_build_mismatch_pair(out_dir: Path) -> CorpusPair:
    """Build mismatch: same variants, every POS shifted by +1,000,000.

    Mimics an hg19 panel merged against an hg38 target. chr:pos intersection
    collapses to ~0; allele content unchanged at corresponding row indices.
    """
    canonical = _build_canonical(out_dir)
    pre_shift = out_dir / "_pre_shift"
    other = out_dir / "other"
    # Copy canonical → pre_shift, then shift onto `other`.
    pre_shift.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pgen", ".pvar", ".psam"):
        shutil.copy(Path(str(canonical) + ext), Path(str(pre_shift) + ext))
    shift_positions(pre_shift, other, per_chrom_offset=1_000_000)
    return CorpusPair(canonical=canonical, other=other, label="build_mismatch")


def build_allele_swap_pair(out_dir: Path) -> CorpusPair:
    """Every non-canonical variant has REF/ALT swapped (+ genotypes x→2-x).

    Intersection at chr:pos is full; allele orientation disagrees at every
    site. Classifier will flag a `dropped_allele_mismatch`-heavy histogram
    if the policy doesn't auto-recode.
    """
    canonical = _build_canonical(out_dir)
    other = out_dir / "other"
    n_variants = _CORPUS_SHAPE.n_variants
    n_samples = _CORPUS_SHAPE.n_samples
    swap_ref_alt(
        canonical,
        n_samples=n_samples,
        n_variants=n_variants,
        variant_indices=np.arange(n_variants),
        out_prefix=other,
    )
    return CorpusPair(canonical=canonical, other=other, label="allele_swap")


def build_strand_flip_pair(out_dir: Path) -> CorpusPair:
    """Every unambiguous variant has REF/ALT complemented.

    Intersection at chr:pos is full; the alignment table will accumulate
    `strand_flip` (or `dropped_allele_mismatch` under --on-strand drop).
    A/T+C/G variants are skipped because they're palindromic and the
    complement is the same pair.
    """
    canonical = _build_canonical(out_dir)
    other = out_dir / "other"
    # Flip every unambiguous variant; ambiguous ones (A/T, C/G) are no-ops
    # under complementation so leaving them flipped-by-helper would still
    # produce the same alleles. `flip_strand` happily complements both
    # but the resulting REF/ALT is unchanged for palindromes.
    flip_strand(canonical, np.arange(_CORPUS_SHAPE.n_variants), other)
    return CorpusPair(canonical=canonical, other=other, label="strand_flip")


def build_disjoint_pair(out_dir: Path) -> CorpusPair:
    """Two unrelated panels with *asymmetric chromosome coverage*.

    Canonical covers (1, 2); other covers (1, 3). Chr1 is shared but at
    independent positions; chr2 is canonical-only; chr3 is other-only.
    This asymmetric-chrom signature is what distinguishes "disjoint
    panels" from "build mismatch" — a build shift would preserve chrom
    presence on both sides, so the classifier can tell the two apart
    from the per-chrom coverage alone.
    """
    canonical = _build_canonical(out_dir)
    other_spec = SyntheticPanelSpec(
        n_samples=_CORPUS_SHAPE.n_samples,
        n_variants=_CORPUS_SHAPE.n_variants,
        n_populations=_CORPUS_SHAPE.n_populations,
        chromosomes=(1, 3),  # asymmetric vs. canonical's (1, 2)
        variant_seed=_CORPUS_SHAPE.variant_seed + 99_991,  # unrelated prime offset
        sample_seed=_CORPUS_SHAPE.sample_seed + 1,
        sample_id_prefix="D",
    )
    other = synthesize_pfile(other_spec, out_dir / "other").path
    return CorpusPair(canonical=canonical, other=other, label="disjoint")


def build_key_space_mismatch_pair(out_dir: Path) -> CorpusPair:
    """Canonical uses chr:pos-style IDs; other uses `rs<N>` IDs.

    chr:pos intersection is full (positions match by construction);
    ID intersection is zero. The classifier should distinguish a
    key-space mismatch from a genuine disjoint case by comparing chr:pos
    overlap against ID overlap.
    """
    canonical = _build_canonical(out_dir)
    pre_rename = out_dir / "_pre_rename"
    other = out_dir / "other"
    pre_rename.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pgen", ".pvar", ".psam"):
        shutil.copy(Path(str(canonical) + ext), Path(str(pre_rename) + ext))
    rename_ids_to_rsid(pre_rename, other, start=1)
    return CorpusPair(canonical=canonical, other=other, label="key_space_mismatch")


# Convenience: map label → builder, so test parametrization can iterate
# the full corpus without re-listing every fixture.
ALL_BUILDERS: dict[str, callable] = {  # type: ignore[type-arg]
    "compatible": build_compatible_pair,
    "build_mismatch": build_build_mismatch_pair,
    "allele_swap": build_allele_swap_pair,
    "strand_flip": build_strand_flip_pair,
    "disjoint": build_disjoint_pair,
    "key_space_mismatch": build_key_space_mismatch_pair,
}
