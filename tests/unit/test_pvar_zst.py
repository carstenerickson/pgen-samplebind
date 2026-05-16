"""Unit tests for `.pvar.zst` (zstd-compressed) input support.

Plink2 v2.0.0-a.6+ writes `.pvar.zst` by default; HGDP+1kGP and similar
reference panels distributed via HuggingFace / Dataverse arrive in this
form. The pvar reader path must accept both `.pvar` and `.pvar.zst`
transparently so users don't have to pre-decompress.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import zstandard

from pgen_samplebind.pvar import (
    _find_header_line,
    count_raw_variants,
    is_zst_path,
    open_pvar_text,
    read_pvar,
)

_SAMPLE_PVAR = (
    "##fileformat=PVARv1.0\n"
    "#CHROM\tPOS\tID\tREF\tALT\n"
    "1\t1000\trs1\tA\tC\n"
    "1\t2000\trs2\tA\tG\n"
    "1\t3000\trs3\tC\tT\n"
    "1\t4000\trs4\tG\tT\n"
    # one non-biallelic-SNP row (should be filtered out by read_pvar)
    "1\t5000\trs5\tAT\tA\n"
)


@pytest.fixture
def pvar_files(tmp_path: Path) -> tuple[Path, Path]:
    """Write the same .pvar to disk as both uncompressed and zst-compressed."""
    plain = tmp_path / "sample.pvar"
    plain.write_text(_SAMPLE_PVAR)
    zst = tmp_path / "sample2.pvar.zst"
    zst.write_bytes(zstandard.ZstdCompressor().compress(_SAMPLE_PVAR.encode()))
    return plain, zst


class TestIsZstPath:
    def test_zst_path(self) -> None:
        assert is_zst_path(Path("foo.pvar.zst")) is True

    def test_plain_pvar(self) -> None:
        assert is_zst_path(Path("foo.pvar")) is False


class TestOpenPvarText:
    def test_reads_plain_pvar(self, pvar_files: tuple[Path, Path]) -> None:
        plain, _ = pvar_files
        with open_pvar_text(plain) as f:
            assert f.readline() == "##fileformat=PVARv1.0\n"

    def test_reads_zst_transparently(self, pvar_files: tuple[Path, Path]) -> None:
        _, zst = pvar_files
        with open_pvar_text(zst) as f:
            lines = f.readlines()
        assert lines[0] == "##fileformat=PVARv1.0\n"
        assert lines[1].startswith("#CHROM\t")
        assert len(lines) == 7  # 2 header + 5 data


class TestCountRawVariants:
    def test_plain(self, pvar_files: tuple[Path, Path]) -> None:
        plain, _ = pvar_files
        assert count_raw_variants(plain) == 5

    def test_zst_matches_plain(self, pvar_files: tuple[Path, Path]) -> None:
        plain, zst = pvar_files
        assert count_raw_variants(zst) == count_raw_variants(plain) == 5


class TestFindHeaderLine:
    def test_skips_metadata_in_zst(self, pvar_files: tuple[Path, Path]) -> None:
        _, zst = pvar_files
        # Line 0 is `##fileformat=PVARv1.0`; the `#CHROM` header is line 1.
        assert _find_header_line(zst) == 1


class TestReadPvar:
    def test_zst_round_trips(self, pvar_files: tuple[Path, Path]) -> None:
        plain, zst = pvar_files
        df_plain = read_pvar(plain)
        df_zst = read_pvar(zst)
        # Biallelic-SNP filter drops rs5 (REF=AT); 4 rows survive on both paths
        assert len(df_plain) == 4
        assert len(df_zst) == 4
        pd.testing.assert_frame_equal(df_plain, df_zst)

    def test_dtypes_preserved_through_zst(self, pvar_files: tuple[Path, Path]) -> None:
        _, zst = pvar_files
        df = read_pvar(zst)
        assert df["chrom"].dtype == "int8"
        assert df["pos"].dtype == "int64"
        assert df["cm"].dtype == "float64"
