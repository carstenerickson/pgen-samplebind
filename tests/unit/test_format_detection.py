"""Unit tests for formats.detect_format — companion-file inference."""

from __future__ import annotations

from pathlib import Path

import pytest

from pgen_samplebind.errors import UsageError
from pgen_samplebind.formats import detect_format
from pgen_samplebind.types import InputFormat


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


class TestPfileDetection:
    def test_pfile_full_triplet_detected(self, tmp_path: Path) -> None:
        prefix = tmp_path / "data"
        for ext in (".pgen", ".pvar", ".psam"):
            _touch(prefix.with_suffix(ext))
        assert detect_format(prefix) is InputFormat.PFILE

    def test_pfile_with_extension_in_prefix(self, tmp_path: Path) -> None:
        """User passes `data.pgen` instead of `data` — still detects PFILE."""
        prefix = tmp_path / "data"
        for ext in (".pgen", ".pvar", ".psam"):
            _touch(prefix.with_suffix(ext))
        # Pass with .pgen extension
        assert detect_format(tmp_path / "data.pgen") is InputFormat.PFILE
        # Pass with .pvar extension
        assert detect_format(tmp_path / "data.pvar") is InputFormat.PFILE

    def test_pfile_missing_pvar_raises(self, tmp_path: Path) -> None:
        prefix = tmp_path / "data"
        _touch(prefix.with_suffix(".pgen"))
        _touch(prefix.with_suffix(".psam"))
        with pytest.raises(UsageError, match="missing companion"):
            detect_format(prefix)


class TestBfileDetection:
    def test_bfile_full_triplet_detected(self, tmp_path: Path) -> None:
        prefix = tmp_path / "data"
        for ext in (".bed", ".bim", ".fam"):
            _touch(prefix.with_suffix(ext))
        assert detect_format(prefix) is InputFormat.BFILE

    def test_bfile_missing_fam_raises(self, tmp_path: Path) -> None:
        prefix = tmp_path / "data"
        _touch(prefix.with_suffix(".bed"))
        _touch(prefix.with_suffix(".bim"))
        with pytest.raises(UsageError, match="missing companion"):
            detect_format(prefix)


class TestEigenstratDetection:
    def test_eigenstrat_full_triplet_detected(self, tmp_path: Path) -> None:
        prefix = tmp_path / "data"
        for ext in (".geno", ".snp", ".ind"):
            _touch(prefix.with_suffix(ext))
        assert detect_format(prefix) is InputFormat.EIGENSTRAT


class TestNoneFound:
    def test_no_triplet_raises(self, tmp_path: Path) -> None:
        with pytest.raises(UsageError, match="no PFILE/BFILE/EIGENSTRAT triplet"):
            detect_format(tmp_path / "nonexistent")
