"""Unit tests for pvar.normalize_chrom — chromosome string→int normalization."""

from __future__ import annotations

import pytest

from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.pvar import normalize_chrom


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
