"""Native ASCII per-line EIGENSTRAT support (HLD §Format detection v2).

plink2 `--eigfile` only reads PACKEDANCESTRYMAP-format EIGENSTRAT (binary,
`GENO `/`TGENO ` header). The older ASCII per-line variant (one digit per
sample-variant cell, no header) is parsed natively by formats.py so users
don't need to pre-convert via convertf. Surfaced during the Track E
Phase 7 dogfood — Carsten's pileupCaller output is single-sample ASCII.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pgenlib
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from pgen_samplebind.formats import (
    InputFormat,
    _convert_ascii_eigenstrat,
    _is_packedancestrymap,
    detect_format,
    prepared_input,
)
from tests.fixtures.modifiers import write_ascii_eigenstrat


def _read_pgen_full(prefix: Path, n_samples: int, n_variants: int) -> np.ndarray:
    pgen_path = Path(str(prefix) + ".pgen")
    reader = pgenlib.PgenReader(str(pgen_path).encode(), raw_sample_ct=n_samples)
    try:
        buf = np.empty((n_variants, n_samples), dtype=np.int8)
        reader.read_range(0, n_variants, buf)
    finally:
        reader.close()
    return buf


class TestPackedancestrymapSniff:
    def test_packed_geno_starts_with_GENO(self, tmp_path: Path) -> None:
        p = tmp_path / "x.geno"
        p.write_bytes(b"GENO 5 10 abc...")
        assert _is_packedancestrymap(p)

    def test_tgeno_packed(self, tmp_path: Path) -> None:
        p = tmp_path / "x.geno"
        p.write_bytes(b"TGENO 5 10 abc")
        assert _is_packedancestrymap(p)

    def test_ascii_geno_does_not_match(self, tmp_path: Path) -> None:
        p = tmp_path / "x.geno"
        p.write_bytes(b"02192\n")
        assert not _is_packedancestrymap(p)

    def test_missing_file_returns_false(self, tmp_path: Path) -> None:
        assert not _is_packedancestrymap(tmp_path / "nope.geno")


class TestAsciiEigenstratNativeParser:
    """Synthesize a tiny ASCII EIGENSTRAT triplet, run it through
    `_convert_ascii_eigenstrat` directly, and verify the resulting PFILE
    has correct genotypes (with eig→plink REF/ALT count flip) and pvar."""

    def test_multisample_three_variants(self, tmp_path: Path) -> None:
        # 4 samples, 3 variants (all autosomal)
        # eig: 0=hom-REF, 1=het, 2=hom-ALT, 9=missing
        # Expected plink genotypes (count of ALT) = 2 - eig (with 9 → -9)
        eig_geno = np.array(
            [
                [0, 1, 2, 9],
                [2, 2, 0, 1],
                [9, 0, 1, 2],
            ],
            dtype=np.int8,
        )
        snp_rows = [
            ("rs1", 1, 0.001, 100, "A", "C"),
            ("rs2", 2, 0.500, 200, "G", "T"),
            ("rs3", 3, 1.000, 300, "C", "A"),
        ]
        ind_rows = [
            ("S1", "M", "POP1"),
            ("S2", "F", "POP1"),
            ("S3", "M", "POP2"),
            ("S4", "U", "POP2"),
        ]
        in_prefix = tmp_path / "ascii_in"
        write_ascii_eigenstrat(in_prefix, eig_geno, snp_rows, ind_rows)

        out_prefix = tmp_path / "out"
        _convert_ascii_eigenstrat(in_prefix, out_prefix, include_chrom=tuple(range(1, 23)))

        # Genotype check: eig 0/1/2/9 → plink 2/1/0/missing
        expected = np.array(
            [
                [2, 1, 0, -9],
                [0, 0, 2, 1],
                [-9, 2, 1, 0],
            ],
            dtype=np.int8,
        )
        actual = _read_pgen_full(out_prefix, n_samples=4, n_variants=3)
        np.testing.assert_array_equal(actual, expected)

        # pvar has CM column preserved
        pvar_lines = (out_prefix.with_suffix(".pvar")).read_text().splitlines()
        assert pvar_lines[0] == "#CHROM\tPOS\tID\tREF\tALT\tCM"
        assert pvar_lines[1] == "1\t100\trs1\tA\tC\t0.001"
        assert pvar_lines[2] == "2\t200\trs2\tG\tT\t0.5"
        assert pvar_lines[3] == "3\t300\trs3\tC\tA\t1.0"

        # psam has FID=0 (orchestrator adds FID=POP later) and PHENO1 = pop label
        psam_lines = (out_prefix.with_suffix(".psam")).read_text().splitlines()
        assert psam_lines[0] == "#FID\tIID\tSEX\tPHENO1"
        assert psam_lines[1] == "0\tS1\t1\tPOP1"  # M → 1
        assert psam_lines[2] == "0\tS2\t2\tPOP1"  # F → 2
        assert psam_lines[4] == "0\tS4\t0\tPOP2"  # U → 0

    def test_single_sample_pseudohaploid_style(self, tmp_path: Path) -> None:
        """Carsten-style: 1 sample, only 0/2/9 values (pseudohaploid)."""
        eig_geno = np.array([[0], [2], [2], [9], [0]], dtype=np.int8)
        snp_rows = [
            ("rs1", 1, 0.0, 100, "A", "C"),
            ("rs2", 1, 0.0, 200, "G", "T"),
            ("rs3", 5, 0.0, 300, "C", "A"),
            ("rs4", 7, 0.0, 400, "T", "G"),
            ("rs5", 22, 0.0, 500, "A", "G"),
        ]
        ind_rows = [("Carsten", "U", "Carsten_target")]
        in_prefix = tmp_path / "ascii_carsten"
        write_ascii_eigenstrat(in_prefix, eig_geno, snp_rows, ind_rows)

        out_prefix = tmp_path / "out"
        _convert_ascii_eigenstrat(in_prefix, out_prefix, include_chrom=tuple(range(1, 23)))

        expected = np.array([[2], [0], [0], [-9], [2]], dtype=np.int8)
        actual = _read_pgen_full(out_prefix, n_samples=1, n_variants=5)
        np.testing.assert_array_equal(actual, expected)

    def test_sex_chrom_filtered_out(self, tmp_path: Path) -> None:
        """include_chrom=(1..22) should drop chr 23 (X) variants at parse time."""
        eig_geno = np.array([[0, 2], [1, 1], [2, 0]], dtype=np.int8)
        snp_rows = [
            ("rs_auto", 5, 0.0, 100, "A", "G"),
            ("rs_x", 23, 0.0, 200, "C", "T"),  # to be dropped
            ("rs_auto2", 10, 0.0, 300, "T", "A"),
        ]
        ind_rows = [("S1", "M", "POP1"), ("S2", "F", "POP1")]
        in_prefix = tmp_path / "ascii_x"
        write_ascii_eigenstrat(in_prefix, eig_geno, snp_rows, ind_rows)

        out_prefix = tmp_path / "out"
        _convert_ascii_eigenstrat(in_prefix, out_prefix, include_chrom=tuple(range(1, 23)))

        # Only autosomes kept → 2 variants
        actual = _read_pgen_full(out_prefix, n_samples=2, n_variants=2)
        # eig rows 0,2 → plink: [2,0], [0,2]
        np.testing.assert_array_equal(actual, np.array([[2, 0], [0, 2]], dtype=np.int8))

    def test_geno_size_mismatch_raises_iofailure(self, tmp_path: Path) -> None:
        """Truncated .geno with wrong byte count → IOFailure."""
        from pgen_samplebind.errors import IOFailure

        in_prefix = tmp_path / "bad"
        # Write .ind and .snp matching 3 variants x 2 samples but truncated .geno
        write_ascii_eigenstrat(
            in_prefix,
            np.array([[0, 1], [2, 0], [1, 2]], dtype=np.int8),
            [
                ("rs1", 1, 0.0, 100, "A", "C"),
                ("rs2", 1, 0.0, 200, "G", "T"),
                ("rs3", 1, 0.0, 300, "C", "A"),
            ],
            [("S1", "M", "POP1"), ("S2", "F", "POP1")],
        )
        # Truncate .geno to 5 bytes (well short of 3 x 3 = 9)
        (in_prefix.with_suffix(".geno")).write_bytes(b"01\n2")

        out_prefix = tmp_path / "out"
        with pytest.raises(IOFailure, match="size mismatch"):
            _convert_ascii_eigenstrat(in_prefix, out_prefix, include_chrom=tuple(range(1, 23)))


class TestAsciiEigenstratEndToEnd:
    """Run the prepared_input → merge end-to-end on an ASCII EIGENSTRAT input."""

    def test_prepared_input_routes_ascii_through_native(self, tmp_path: Path) -> None:
        eig_geno = np.array([[0, 1], [2, 9], [1, 0]], dtype=np.int8)
        snp_rows = [
            ("rs1", 1, 0.0, 100, "A", "C"),
            ("rs2", 1, 0.0, 200, "G", "T"),
            ("rs3", 2, 0.0, 300, "C", "A"),
        ]
        ind_rows = [("S1", "M", "POP1"), ("S2", "F", "POP2")]
        in_prefix = tmp_path / "ascii_in"
        write_ascii_eigenstrat(in_prefix, eig_geno, snp_rows, ind_rows)

        # detect_format identifies it as EIGENSTRAT
        assert detect_format(in_prefix) is InputFormat.EIGENSTRAT

        with prepared_input(in_prefix) as desc:
            # Tempdir-routed PFILE should exist and be readable
            assert desc.fmt is InputFormat.EIGENSTRAT
            assert desc.pgen_path.exists()
            buf = _read_pgen_full(
                Path(str(desc.pgen_path).removesuffix(".pgen")),
                n_samples=2,
                n_variants=3,
            )
            np.testing.assert_array_equal(buf, np.array([[2, 1], [0, -9], [1, 2]], dtype=np.int8))

    def test_merge_two_ascii_inputs(self, tmp_path: Path) -> None:
        """Two ASCII inputs with disjoint samples + same variants → merge cleanly."""
        snp_rows = [
            ("rs1", 1, 0.001, 100, "A", "C"),
            ("rs2", 1, 0.005, 200, "G", "T"),
            ("rs3", 2, 0.010, 300, "C", "A"),
        ]
        # Panel A: 2 samples
        write_ascii_eigenstrat(
            tmp_path / "a",
            np.array([[0, 1], [2, 0], [1, 2]], dtype=np.int8),
            snp_rows,
            [("A1", "M", "POPA"), ("A2", "F", "POPA")],
        )
        # Panel B: 3 samples
        write_ascii_eigenstrat(
            tmp_path / "b",
            np.array([[2, 0, 1], [0, 1, 2], [9, 2, 0]], dtype=np.int8),
            snp_rows,
            [("B1", "M", "POPB"), ("B2", "F", "POPB"), ("B3", "U", "POPB")],
        )

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(tmp_path / "a"),
                str(tmp_path / "b"),
                "-o",
                str(out),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        # 5 samples (2 + 3, no collisions); 3 variants
        buf = _read_pgen_full(out, n_samples=5, n_variants=3)
        # Plink-encoded: 2-eig for non-missing; 9 → -9
        # A: [[2,1],[0,2],[1,0]];  B: [[0,2,1],[2,1,0],[-9,0,2]]
        # Concat by samples (A then B):
        expected = np.array(
            [
                [2, 1, 0, 2, 1],
                [0, 2, 2, 1, 0],
                [1, 0, -9, 0, 2],
            ],
            dtype=np.int8,
        )
        np.testing.assert_array_equal(buf, expected)

        # Output .pvar carries the CM column from the input .snp
        pvar_lines = (out.with_suffix(".pvar")).read_text().splitlines()
        assert pvar_lines[0] == "#CHROM\tPOS\tID\tREF\tALT\tCM"
        assert pvar_lines[1].endswith("\t0.001")
        assert pvar_lines[2].endswith("\t0.005")
        assert pvar_lines[3].endswith("\t0.01")
