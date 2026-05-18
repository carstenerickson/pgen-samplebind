"""Unit tests for parser robustness on malformed `.pvar` and `.psam` files.

These exercise the user-induced failure modes a hand-edited or
mis-converted input might produce — coverage gaps surfaced in the v0.3.2
review. The contract is: surface a clear `IOFailure` (or
`InvariantViolation` where appropriate) rather than a cryptic pandas /
numpy traceback.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pgen_samplebind.errors import IOFailure
from pgen_samplebind.psam import read_psam
from pgen_samplebind.pvar import read_pvar


class TestMalformedPvar:
    def test_missing_chrom_header_raises_iofailure(self, tmp_path: Path) -> None:
        """A .pvar with no `#CHROM` line at all (e.g., header line stripped)
        triggers the no-header diagnostic."""
        p = tmp_path / "bad.pvar"
        p.write_text("1\t1000\trs1\tA\tC\n2\t2000\trs2\tG\tT\n")
        with pytest.raises(IOFailure, match="no #CHROM header"):
            read_pvar(p)

    def test_missing_required_column_raises_iofailure(self, tmp_path: Path) -> None:
        """Header present but missing one of the required columns (REF here)."""
        p = tmp_path / "no_ref.pvar"
        p.write_text("#CHROM\tPOS\tID\tALT\n1\t1000\trs1\tC\n")
        with pytest.raises(IOFailure, match="missing required columns.*REF"):
            read_pvar(p)

    def test_extra_metadata_lines_tolerated(self, tmp_path: Path) -> None:
        """Multiple `##`-prefixed metadata lines before the header should be
        skipped, not rejected."""
        p = tmp_path / "with_metadata.pvar"
        p.write_text(
            "##fileformat=PVARv1.0\n"
            "##source=plink2\n"
            "##contig=<ID=1>\n"
            "#CHROM\tPOS\tID\tREF\tALT\n"
            "1\t1000\trs1\tA\tC\n"
        )
        df = read_pvar(p)
        assert len(df) == 1
        assert df["chrom"].tolist() == [1]

    def test_short_data_row_silently_filtered(self, tmp_path: Path) -> None:
        """Pinned current behavior, not necessarily desired: pandas with
        na_filter=False fills missing trailing fields with empty strings,
        and the biallelic-SNP filter (empty-REF/ALT) silently drops the
        row. Worth surfacing in a future fix; for now this test pins the
        behavior so anyone changing it has to update the assertion."""
        p = tmp_path / "short_row.pvar"
        p.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\n"
            "1\t1000\trs1\tA\tC\n"
            "1\t2000\trs2\n"  # short row — silently dropped
            "1\t3000\trs3\tG\tT\n"
        )
        df = read_pvar(p)
        assert len(df) == 2
        assert df["id"].tolist() == ["rs1", "rs3"]


class TestMalformedPsam:
    def test_missing_pound_header_raises_iofailure(self, tmp_path: Path) -> None:
        """psam without a `#`-prefixed header line is rejected with a
        diagnostic that names the expected form."""
        p = tmp_path / "no_header.psam"
        p.write_text("S1\t1\tPOP1\nS2\t2\tPOP2\n")
        with pytest.raises(IOFailure, match="no `#`-prefixed header"):
            read_psam(p)

    def test_metadata_only_file_raises_iofailure(self, tmp_path: Path) -> None:
        """psam that contains only `##` metadata (no real header) is also
        rejected. Reachable on an empty-after-comments edge case."""
        p = tmp_path / "metadata_only.psam"
        p.write_text("##fileformat=PSAMv1.0\n##created=2026-05-18\n")
        with pytest.raises(IOFailure, match="no `#`-prefixed header"):
            read_psam(p)

    def test_metadata_then_header_tolerated(self, tmp_path: Path) -> None:
        """Multiple `##` lines before the header should be skipped."""
        p = tmp_path / "with_metadata.psam"
        p.write_text(
            "##fileformat=PSAMv1.0\n"
            "##source=plink2\n"
            "#IID\tSEX\tPOP\n"
            "S1\t1\tPOP1\n"
        )
        df = read_psam(p)
        assert len(df) == 1
        assert df["IID"].tolist() == ["S1"]
