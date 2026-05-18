"""Unit tests for hashing.canonicalize_pvar_bytes — canonical-hash spec.

Tests the canonicalization spec at the byte level (no SHA-256 involved).
Order invariance and ID invariance are exercised here at the unit level;
the integration version against on-disk files is in tests/integration/.
"""

from __future__ import annotations

import pandas as pd

from pgen_samplebind.hashing import canonicalize_pvar_bytes


def _make_pvar_df(
    chroms: list[int], positions: list[int], ids: list[str], refs: list[str], alts: list[str]
) -> pd.DataFrame:
    """Helper: build a DataFrame in the schema read_pvar produces."""
    return pd.DataFrame(
        {"chrom": chroms, "pos": positions, "id": ids, "ref": refs, "alt": alts}
    ).astype({"chrom": "int8", "pos": "int64"})


class TestCanonicalizationFormat:
    def test_single_variant(self) -> None:
        df = _make_pvar_df([1], [1000], ["rs1"], ["A"], ["G"])
        assert canonicalize_pvar_bytes(df) == b"1\t1000\tA\tG\n"

    def test_multiple_variants_emit_one_line_each(self) -> None:
        df = _make_pvar_df(
            [1, 1, 2], [100, 200, 50], ["rs1", "rs2", "rs3"], ["A", "C", "G"], ["G", "T", "T"]
        )
        out = canonicalize_pvar_bytes(df)
        # 3 lines, all newline-terminated
        assert out.count(b"\n") == 3
        # Tab-separated: 4 columns each line
        for line in out.strip().split(b"\n"):
            assert line.count(b"\t") == 3


class TestSortInvariance:
    """Same variant set in two different orderings → same hash."""

    def test_chromosome_order_irrelevant(self) -> None:
        df_a = _make_pvar_df([1, 2], [100, 200], ["rs1", "rs2"], ["A", "C"], ["G", "T"])
        df_b = _make_pvar_df([2, 1], [200, 100], ["rs2", "rs1"], ["C", "A"], ["T", "G"])
        assert canonicalize_pvar_bytes(df_a) == canonicalize_pvar_bytes(df_b)

    def test_position_order_irrelevant(self) -> None:
        df_a = _make_pvar_df(
            [1, 1, 1], [100, 200, 300], ["a", "b", "c"], ["A", "C", "G"], ["G", "T", "T"]
        )
        df_b = _make_pvar_df(
            [1, 1, 1], [300, 100, 200], ["c", "a", "b"], ["G", "A", "C"], ["T", "G", "T"]
        )
        assert canonicalize_pvar_bytes(df_a) == canonicalize_pvar_bytes(df_b)

    def test_numeric_not_lexicographic_sort(self) -> None:
        """Ensure (1, 100) sorts before (1, 9) numerically — lex would put '9' first."""
        df = _make_pvar_df([1, 1], [9, 100], ["a", "b"], ["A", "C"], ["G", "T"])
        out = canonicalize_pvar_bytes(df)
        first_line = out.split(b"\n")[0]
        assert first_line == b"1\t9\tA\tG"


class TestIDInvariance:
    """Same variant set with different ID conventions → same hash."""

    def test_rs_vs_chrposrefalt_id_irrelevant(self) -> None:
        df_a = _make_pvar_df([1, 2], [100, 200], ["rs1", "rs2"], ["A", "C"], ["G", "T"])
        df_b = _make_pvar_df([1, 2], [100, 200], ["1:100:A:G", "2:200:C:T"], ["A", "C"], ["G", "T"])
        assert canonicalize_pvar_bytes(df_a) == canonicalize_pvar_bytes(df_b)


class TestRefAltTiebreak:
    """When (chrom, pos) collide (e.g., a tri-allelic site stored as two
    biallelic rows), the canonical order must still be well-defined so the
    hash is invariant to the original physical row order across formats.
    """

    def test_same_locus_two_alts_order_invariant(self) -> None:
        df_a = _make_pvar_df([1, 1], [100, 100], ["v1", "v2"], ["A", "A"], ["C", "T"])
        df_b = _make_pvar_df([1, 1], [100, 100], ["v2", "v1"], ["A", "A"], ["T", "C"])
        assert canonicalize_pvar_bytes(df_a) == canonicalize_pvar_bytes(df_b)

    def test_same_locus_alt_then_ref_tiebreak(self) -> None:
        """ref tiebreaks before alt: rows with same (chrom, pos) emit in
        ref-ascending then alt-ascending order."""
        df = _make_pvar_df([1, 1], [100, 100], ["a", "b"], ["G", "A"], ["T", "C"])
        out = canonicalize_pvar_bytes(df)
        assert out == b"1\t100\tA\tC\n1\t100\tG\tT\n"
