"""Unit tests for pvar.normalize_chrom — chromosome string→int normalization."""

from __future__ import annotations

import pandas as pd
import pytest

from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.pvar import normalize_chrom, normalize_chrom_series


class TestNumericChromosomes:
    @pytest.mark.parametrize("input_value, expected", [(i, i) for i in range(1, 23)])
    def test_autosome_strings(self, input_value: int, expected: int) -> None:
        assert normalize_chrom(str(input_value)) == expected

    @pytest.mark.parametrize("input_value", [1, 22, 23, 26])
    def test_int_input_pass_through(self, input_value: int) -> None:
        assert normalize_chrom(input_value) == input_value


class TestSexChromosomesAndMito:
    @pytest.mark.parametrize("input_str", ["X", "x", "chrX", "chrx", "23"])
    def test_x_chromosome(self, input_str: str) -> None:
        assert normalize_chrom(input_str) == 23

    @pytest.mark.parametrize("input_str", ["Y", "y", "chrY", "24"])
    def test_y_chromosome(self, input_str: str) -> None:
        assert normalize_chrom(input_str) == 24

    @pytest.mark.parametrize("input_str", ["MT", "mt", "M", "m", "chrMT", "chrM", "26"])
    def test_mitochondrial(self, input_str: str) -> None:
        assert normalize_chrom(input_str) == 26


class TestChrPrefix:
    @pytest.mark.parametrize("input_str, expected", [("chr1", 1), ("CHR1", 1), ("chr22", 22)])
    def test_chr_prefix_stripped(self, input_str: str, expected: int) -> None:
        assert normalize_chrom(input_str) == expected


class TestInvalidInputs:
    @pytest.mark.parametrize("input_str", ["foo", "27", "-1", "0", "23.5", ""])
    def test_unparseable_raises(self, input_str: str) -> None:
        with pytest.raises(InvariantViolation):
            normalize_chrom(input_str)

    def test_int_out_of_range_raises(self) -> None:
        with pytest.raises(InvariantViolation):
            normalize_chrom(0)
        with pytest.raises(InvariantViolation):
            normalize_chrom(27)


class TestNormalizeChromSeries:
    """v0.3: vectorized variant for the read_pvar / ASCII-EIGENSTRAT hot
    paths. Same semantics as scalar `normalize_chrom`; single pandas pass."""

    def test_mixed_inputs_vectorized(self) -> None:
        s = pd.Series(["1", "22", "X", "chrX", "Y", "chrY", "MT", "chrM", "23", "26"])
        out = normalize_chrom_series(s)
        assert out.tolist() == [1, 22, 23, 23, 24, 24, 26, 26, 23, 26]
        assert out.dtype == "int8"

    def test_chr_prefix_case_insensitive(self) -> None:
        s = pd.Series(["chr1", "CHR1", "Chr1", "cHr22"])
        out = normalize_chrom_series(s)
        assert out.tolist() == [1, 1, 1, 22]

    def test_letter_case_insensitive(self) -> None:
        s = pd.Series(["x", "X", "y", "Y", "m", "M", "mt", "MT"])
        out = normalize_chrom_series(s)
        assert out.tolist() == [23, 23, 24, 24, 26, 26, 26, 26]

    def test_unparseable_raises_first_offender(self) -> None:
        s = pd.Series(["1", "22", "foo", "3"])
        with pytest.raises(InvariantViolation, match=r"unparseable.*'foo'"):
            normalize_chrom_series(s)

    @pytest.mark.parametrize("bad", ["0", "27", "-1", "23.5", "", "foo"])
    def test_individual_bad_values(self, bad: str) -> None:
        with pytest.raises(InvariantViolation):
            normalize_chrom_series(pd.Series([bad]))

    def test_whitespace_tolerated(self) -> None:
        s = pd.Series([" 1 ", "\tchrX\n", "  MT  "])
        out = normalize_chrom_series(s)
        assert out.tolist() == [1, 23, 26]

    def test_empty_series(self) -> None:
        s = pd.Series([], dtype=object)
        out = normalize_chrom_series(s)
        assert out.empty
        assert out.dtype == "int8"

    def test_matches_scalar_on_large_random_input(self) -> None:
        """Vectorized output must equal scalar `.map(normalize_chrom)` on the
        full case-space — the perf change must be a pure speedup, not a
        semantic change."""
        cases = ["1", "22", "X", "chrX", "chrx", "Y", "MT", "chrM", "23", "26", "chr5"] * 1000
        s = pd.Series(cases)
        vec_out = normalize_chrom_series(s)
        scalar_out = s.map(normalize_chrom).astype("int8")
        assert vec_out.tolist() == scalar_out.tolist()
