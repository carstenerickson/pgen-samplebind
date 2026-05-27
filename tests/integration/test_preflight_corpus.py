"""Behavioral tests against the preflight failure-mode corpus.

These pin observable preflight numbers (intersection size, per-chrom
breakdown) for each failure mode. Step 3 (classifier) will add
classification-label assertions using the same corpus; step 4 (gate)
will add policy-driven warn/error assertions. Asserting numerics today
guards against silent regressions while those follow-ups land.
"""

from __future__ import annotations

import pytest

from pgen_samplebind.formats import prepared_input
from pgen_samplebind.preflight import compute_preflight
from pgen_samplebind.types import MergePolicy
from tests.fixtures.preflight_corpus import ALL_BUILDERS, CorpusPair


@pytest.fixture(scope="session")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, CorpusPair]:
    """Build the full failure-mode corpus once per test session."""
    out: dict[str, CorpusPair] = {}
    for label, builder in ALL_BUILDERS.items():
        sub = tmp_path_factory.mktemp(f"preflight_{label}")
        out[label] = builder(sub)
    return out


def _run_preflight(pair: CorpusPair, *, variant_key: str = "chr_pos"):
    """Open both prefixes via the real format-detection path and run preflight."""
    with prepared_input(pair.canonical, is_target=False, include_chrom=tuple(range(1, 23))) as a, \
         prepared_input(pair.other, is_target=False, include_chrom=tuple(range(1, 23))) as b:
        policy = MergePolicy(variant_key=variant_key)  # type: ignore[arg-type]
        return compute_preflight([a, b], policy, tool_version="test", command="merge")


def test_corpus_compatible_full_intersection(corpus: dict[str, CorpusPair]) -> None:
    """Happy path: identical variants → ~100% chr:pos intersection."""
    report = _run_preflight(corpus["compatible"])
    pair = report.comparisons[0]
    assert pair.intersection_fraction_of_min == pytest.approx(1.0)
    assert pair.intersection == pair.n_variants


def test_corpus_build_mismatch_collapses_chr_pos(corpus: dict[str, CorpusPair]) -> None:
    """Position-shifted panel: chr:pos intersection collapses to ~0%."""
    report = _run_preflight(corpus["build_mismatch"])
    pair = report.comparisons[0]
    # Positions all shifted by 1M → zero overlap (no coincidental collisions
    # at 40 variants in a 100M-position space).
    assert pair.intersection == 0
    assert pair.intersection_fraction_of_min == 0.0
    # n_variants on the other side preserved (only POS column changed).
    assert pair.n_variants == report.canonical["n_variants"]


def test_corpus_allele_swap_preserves_chr_pos(corpus: dict[str, CorpusPair]) -> None:
    """REF/ALT swap leaves chr:pos sets untouched — the allele-level
    failure surfaces in alignment.action_histogram, not in preflight's
    key-space view. Step 3's classifier will need that downstream signal."""
    report = _run_preflight(corpus["allele_swap"])
    pair = report.comparisons[0]
    assert pair.intersection == pair.n_variants
    assert pair.intersection_fraction_of_min == pytest.approx(1.0)


def test_corpus_strand_flip_preserves_chr_pos(corpus: dict[str, CorpusPair]) -> None:
    """Complemented alleles also leave chr:pos sets untouched. Same caveat
    as allele_swap: signal lives in the alignment histogram."""
    report = _run_preflight(corpus["strand_flip"])
    pair = report.comparisons[0]
    assert pair.intersection == pair.n_variants


def test_corpus_disjoint_near_zero_intersection(corpus: dict[str, CorpusPair]) -> None:
    """Independent variant seeds: vanishingly small coincidental overlap."""
    report = _run_preflight(corpus["disjoint"])
    pair = report.comparisons[0]
    assert pair.intersection_fraction_of_min < 0.1


def test_corpus_key_space_mismatch_chr_pos_full_id_empty(
    corpus: dict[str, CorpusPair],
) -> None:
    """Same coordinates, rsID IDs on one side. chr:pos intersection is full;
    `--variant-key id` would see zero overlap. The classifier
    distinguishes key-space mismatch from disjoint by comparing both."""
    pair_fixture = corpus["key_space_mismatch"]
    chr_pos_report = _run_preflight(pair_fixture, variant_key="chr_pos")
    id_report = _run_preflight(pair_fixture, variant_key="id")

    assert chr_pos_report.comparisons[0].intersection_fraction_of_min == pytest.approx(1.0)
    assert id_report.comparisons[0].intersection == 0


@pytest.mark.parametrize(
    "label,variant_key,expected_classification",
    [
        ("compatible", "chr_pos", "compatible"),
        ("build_mismatch", "chr_pos", "build_mismatch"),
        # Allele-swap leaves chr:pos intersection at 100% → classifier (which
        # only sees key-space) labels these as `compatible`. The actual
        # allele-disagreement signal lives in alignment.action_histogram
        # (post-pass-1); a future commit may hoist a histogram-driven label
        # in by extending the classifier or adding a post-merge pass.
        ("allele_swap", "chr_pos", "compatible"),
        ("strand_flip", "chr_pos", "compatible"),
        ("disjoint", "chr_pos", "disjoint_panels"),
        # The key_space_mismatch fixture has matching chr:pos and divergent
        # IDs. Under variant_key=chr_pos the active key already works (full
        # intersection) → compatible. The mistake-shape — user specified
        # `--variant-key id` against panels with inconsistent IDs — surfaces
        # only when the active key is `id`, hence the second row below.
        ("key_space_mismatch", "chr_pos", "compatible"),
        ("key_space_mismatch", "id", "key_space_mismatch"),
    ],
)
def test_corpus_classifier_labels(
    corpus: dict[str, CorpusPair],
    label: str,
    variant_key: str,
    expected_classification: str,
) -> None:
    """Ground-truth classifier outputs for every failure mode in the corpus.

    Step-3 contract: each corpus pair has a known intent (the `label`
    field on `CorpusPair`); the classifier must agree on the four shapes
    it can actually distinguish from key-space data alone (compatible,
    build_mismatch, disjoint_panels, key_space_mismatch). allele_swap and
    strand_flip are deliberately not yet detected — see comment above.
    """
    report = _run_preflight(corpus[label], variant_key=variant_key)
    pair = report.comparisons[0]
    assert pair.classification == expected_classification
    # Evidence is always populated alongside classification.
    assert isinstance(pair.classification_evidence, dict)
    assert "active_key_fraction" in pair.classification_evidence
    assert "alternate_key_fraction" in pair.classification_evidence


def test_corpus_key_space_mismatch_evidence_keys_on_alternate(
    corpus: dict[str, CorpusPair],
) -> None:
    """The key_space_mismatch label only fires when alternate-key overlap
    is substantially better than active-key overlap. Pin the evidence so
    a future threshold tweak can't silently re-classify the fixture."""
    report = _run_preflight(corpus["key_space_mismatch"])
    pair = report.comparisons[0]
    ev = pair.classification_evidence
    assert ev is not None
    # Active key (chr_pos) is full → fraction near 1.0. Alternate key (id)
    # is empty → fraction 0.0. So this fixture should classify as
    # `compatible` (active is already great). The inverted shape — where
    # active is id and the user is hitting the rsID/chr:pos mismatch from
    # the opposite direction — fires the label.
    assert pair.classification == "compatible"

    # Now query with variant_key=id (the user's mistake): active drops to 0,
    # alternate is 1.0, lift >= 0.4 → key_space_mismatch.
    report_id = _run_preflight(corpus["key_space_mismatch"], variant_key="id")
    pair_id = report_id.comparisons[0]
    assert pair_id.classification == "key_space_mismatch"
    ev_id = pair_id.classification_evidence
    assert ev_id is not None
    assert ev_id["alternate_key_fraction"] - ev_id["active_key_fraction"] >= 0.4


def test_corpus_build_mismatch_per_chrom_signature(corpus: dict[str, CorpusPair]) -> None:
    """Build mismatch leaves each chromosome with equal-sized but disjoint
    key sets. That's the diagnostic signature: canonical_size > 0 AND
    other_size > 0 AND intersection == 0 on every chrom. Classifier will
    key on exactly this shape vs. the disjoint case (where chrom presence
    itself may differ)."""
    report = _run_preflight(corpus["build_mismatch"])
    pair = report.comparisons[0]
    assert len(pair.per_chrom) > 0
    for pc in pair.per_chrom:
        assert pc.canonical_size > 0
        assert pc.other_size > 0
        assert pc.intersection == 0
