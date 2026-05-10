"""Unit tests for alignment.resolve_alleles — HLD §Allele resolution truth table."""

from __future__ import annotations

import pytest

from pgen_samplebind.alignment import resolve_alleles
from pgen_samplebind.types import DropReason, MergeAction


class TestPassthrough:
    """Same alleles, same orientation."""

    @pytest.mark.parametrize(
        "ref,alt",
        [
            ("A", "C"),
            ("A", "G"),
            ("C", "T"),
            ("G", "T"),
            ("C", "A"),
            ("G", "A"),
            ("T", "C"),
            ("T", "G"),
        ],
    )
    def test_unambiguous_same(self, ref: str, alt: str) -> None:
        action, reason = resolve_alleles(ref, alt, ref, alt, trust_strand=False)
        assert action == MergeAction.PASSTHROUGH
        assert reason is None


class TestRefAltSwap:
    """Same alleles, swapped orientation, non-ambiguous."""

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [
            ("A", "G", "G", "A"),
            ("A", "C", "C", "A"),
            ("C", "T", "T", "C"),
            ("G", "T", "T", "G"),
        ],
    )
    def test_swap(self, c_ref: str, c_alt: str, o_ref: str, o_alt: str) -> None:
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=False)
        assert action == MergeAction.REF_ALT_SWAP
        assert reason is None


class TestStrandFlip:
    """Other alleles complement-of-canonical, same orientation."""

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [
            # A/G -> T/C (complement, no swap)
            ("A", "G", "T", "C"),
            # A/C -> T/G (complement, no swap)
            ("A", "C", "T", "G"),
            # C/T -> G/A
            ("C", "T", "G", "A"),
            # G/T -> C/A
            ("G", "T", "C", "A"),
        ],
    )
    def test_strand_flip(self, c_ref: str, c_alt: str, o_ref: str, o_alt: str) -> None:
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=False)
        assert action == MergeAction.STRAND_FLIP
        assert reason is None


class TestStrandFlipAndSwap:
    """Other alleles complement-of-canonical, swapped orientation."""

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [
            # A/G -> C/T (complement of G/A which is canonical swapped)
            ("A", "G", "C", "T"),
            # A/C -> G/T
            ("A", "C", "G", "T"),
            # C/T -> A/G
            ("C", "T", "A", "G"),
            # G/T -> A/C
            ("G", "T", "A", "C"),
        ],
    )
    def test_flip_and_swap(self, c_ref: str, c_alt: str, o_ref: str, o_alt: str) -> None:
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=False)
        assert action == MergeAction.STRAND_FLIP_AND_SWAP
        assert reason is None


class TestAmbiguousSnps:
    """A/T and C/G ambiguous SNPs — strand cannot be inferred from alleles alone."""

    @pytest.mark.parametrize("ref,alt", [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")])
    def test_ambig_same_passthrough(self, ref: str, alt: str) -> None:
        """Same canonical and other alleles → PASSTHROUGH even without trust_strand."""
        action, reason = resolve_alleles(ref, alt, ref, alt, trust_strand=False)
        assert action == MergeAction.PASSTHROUGH
        assert reason is None

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [("A", "T", "T", "A"), ("T", "A", "A", "T"), ("C", "G", "G", "C"), ("G", "C", "C", "G")],
    )
    def test_ambig_swapped_drops_by_default(
        self, c_ref: str, c_alt: str, o_ref: str, o_alt: str
    ) -> None:
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=False)
        assert action == MergeAction.DROP
        assert reason == DropReason.AMBIGUOUS_STRAND

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [("A", "T", "T", "A"), ("C", "G", "G", "C")],
    )
    def test_ambig_swapped_passes_with_trust_strand(
        self, c_ref: str, c_alt: str, o_ref: str, o_alt: str
    ) -> None:
        """Per HLD: --trust-strand returns PASSTHROUGH for ambiguous swapped pairs."""
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=True)
        assert action == MergeAction.PASSTHROUGH
        assert reason is None

    def test_ambig_canonical_with_nonambig_other_drops_mismatch(self) -> None:
        """A/T canonical with A/G other: different allele set → allele mismatch."""
        action, reason = resolve_alleles("A", "T", "A", "G", trust_strand=False)
        assert action == MergeAction.DROP
        assert reason == DropReason.ALLELE_MISMATCH


class TestAlleleMismatch:
    """Combinations that don't fit any resolution."""

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [
            ("A", "G", "A", "T"),  # different ALT, can't resolve
            ("A", "G", "G", "T"),  # different ALT
            ("A", "G", "C", "G"),  # different REF
            ("C", "T", "A", "T"),  # different REF
            ("A", "G", "A", "C"),  # different ALT
        ],
    )
    def test_mismatch(self, c_ref: str, c_alt: str, o_ref: str, o_alt: str) -> None:
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=False)
        assert action == MergeAction.DROP
        assert reason == DropReason.ALLELE_MISMATCH


class TestDefensiveAgainstNonAcgt:
    """Non-ACGT alleles (shouldn't reach resolve_alleles per pvar filter, but defend)."""

    @pytest.mark.parametrize(
        "c_ref,c_alt,o_ref,o_alt",
        [
            ("N", "G", "A", "G"),
            ("A", "G", "N", "G"),
            ("A", "0", "A", "G"),
            ("AT", "G", "A", "G"),  # multi-char (indel-ish)
        ],
    )
    def test_non_acgt_drops_as_mismatch(
        self, c_ref: str, c_alt: str, o_ref: str, o_alt: str
    ) -> None:
        action, reason = resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand=False)
        assert action == MergeAction.DROP
        assert reason == DropReason.ALLELE_MISMATCH
