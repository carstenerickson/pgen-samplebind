"""Unit tests for `pvar.check_pvar_pgen_row_count_consistent` — the
mis-paired-triplet guardrail added alongside the issue #10 _pgen_row fix.

The guardrail asserts that `count_raw_variants(.pvar)` matches
`pgenlib.PgenReader.get_variant_ct(.pgen)`. Without it, a row-count
mismatch silently over- or under-reads the .pgen and corrupts dosages.

These tests pin both directions of the mismatch, the happy path, the
`.pvar.zst` resolution branch, and the optional `n_pvar` parameter that
callers (merge_cmd, inspect_cmd) use to skip a redundant full-file
scan when they've already computed the count.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib
import pytest
import zstandard

from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.pvar import check_pvar_pgen_row_count_consistent


def _write_pfile_pair(prefix: Path, n_pvar_rows: int, n_pgen_rows: int, n_samples: int = 4) -> Path:
    """Write a (.pvar, .pgen) pair whose row counts can disagree. Returns
    the .pgen path (callers pass this to the guardrail; the .pvar sibling
    is auto-resolved). Genotypes are all-zero — only the row count matters
    for these tests, not the dosage values.
    """
    pvar_df = pd.DataFrame(
        [("1", i * 100, f"v{i}", "A", "C") for i in range(n_pvar_rows)],
        columns=["#CHROM", "POS", "ID", "REF", "ALT"],
    )
    pvar_path = Path(str(prefix) + ".pvar")
    pvar_df.to_csv(pvar_path, sep="\t", index=False, lineterminator="\n")

    pgen_path = Path(str(prefix) + ".pgen")
    geno = np.zeros((n_pgen_rows, n_samples), dtype=np.int8)
    writer = pgenlib.PgenWriter(str(pgen_path).encode(), n_samples, n_pgen_rows)
    try:
        writer.append_biallelic_batch(geno)
    finally:
        writer.close()
    return pgen_path


class TestHappyPath:
    """Matching row counts should return silently with the count."""

    def test_returns_pvar_line_count_when_consistent(self, tmp_path: Path) -> None:
        pgen_path = _write_pfile_pair(tmp_path / "ok", n_pvar_rows=7, n_pgen_rows=7)
        result = check_pvar_pgen_row_count_consistent(pgen_path)
        assert result == 7


class TestMismatchDetection:
    """Either direction of count mismatch must raise InvariantViolation
    with a message that names both counts so a user can grep for them."""

    def test_pvar_has_fewer_rows_than_pgen(self, tmp_path: Path) -> None:
        pgen_path = _write_pfile_pair(tmp_path / "fewer", n_pvar_rows=3, n_pgen_rows=5)
        with pytest.raises(InvariantViolation) as exc_info:
            check_pvar_pgen_row_count_consistent(pgen_path)
        msg = str(exc_info.value)
        assert "variant_ct (5)" in msg
        assert "data-line count (3)" in msg

    def test_pvar_has_more_rows_than_pgen(self, tmp_path: Path) -> None:
        pgen_path = _write_pfile_pair(tmp_path / "more", n_pvar_rows=10, n_pgen_rows=4)
        with pytest.raises(InvariantViolation) as exc_info:
            check_pvar_pgen_row_count_consistent(pgen_path)
        msg = str(exc_info.value)
        assert "variant_ct (4)" in msg
        assert "data-line count (10)" in msg


class TestNPvarParameter:
    """Callers that already computed `count_raw_variants` (merge_cmd's
    step-5 cache, inspect_cmd's n_pre_filter) can pass it through to
    skip a redundant full-file scan. The parameter must be honored
    verbatim — including a value that disagrees with the on-disk .pvar,
    so callers stay in control of the comparison.
    """

    def test_n_pvar_provided_matches_pgen_returns_value(self, tmp_path: Path) -> None:
        pgen_path = _write_pfile_pair(tmp_path / "match", n_pvar_rows=6, n_pgen_rows=6)
        assert check_pvar_pgen_row_count_consistent(pgen_path, n_pvar=6) == 6

    def test_n_pvar_provided_mismatches_pgen_raises(self, tmp_path: Path) -> None:
        # .pvar and .pgen both have 6 rows; caller passes a wrong n_pvar.
        # The guardrail must use the passed value (not re-read the .pvar)
        # so a stale cached count surfaces as a clear failure instead of
        # silently passing.
        pgen_path = _write_pfile_pair(tmp_path / "stale", n_pvar_rows=6, n_pgen_rows=6)
        with pytest.raises(InvariantViolation) as exc_info:
            check_pvar_pgen_row_count_consistent(pgen_path, n_pvar=99)
        msg = str(exc_info.value)
        assert "variant_ct (6)" in msg
        assert "data-line count (99)" in msg

    def test_n_pvar_overrides_disk_scan(self, tmp_path: Path) -> None:
        # On-disk .pvar has 3 rows, .pgen has 5. Without n_pvar, the
        # mismatch raises. With n_pvar=5 (caller asserts equality),
        # the guardrail trusts the caller and passes — the contract is
        # "verify against the supplied count," not "always re-scan."
        pgen_path = _write_pfile_pair(tmp_path / "override", n_pvar_rows=3, n_pgen_rows=5)
        assert check_pvar_pgen_row_count_consistent(pgen_path, n_pvar=5) == 5


class TestZstSiblingResolution:
    """When the .pvar sibling is gzstd-compressed (.pvar.zst), the
    guardrail must resolve to it correctly. Mirrors `check_max_alleles`'s
    sibling-resolution logic; pinned here so the two paths can't drift.
    """

    def test_zst_sibling_is_resolved_and_counted(self, tmp_path: Path) -> None:
        prefix = tmp_path / "compressed"
        n_samples = 4
        # Write a 4-row .pvar then compress and remove the plain copy.
        pvar_df = pd.DataFrame(
            [("1", i * 100, f"v{i}", "A", "C") for i in range(4)],
            columns=["#CHROM", "POS", "ID", "REF", "ALT"],
        )
        plain_pvar = Path(str(prefix) + ".pvar")
        pvar_df.to_csv(plain_pvar, sep="\t", index=False, lineterminator="\n")
        zst_pvar = Path(str(prefix) + ".pvar.zst")
        zst_pvar.write_bytes(zstandard.ZstdCompressor().compress(plain_pvar.read_bytes()))
        plain_pvar.unlink()

        pgen_path = Path(str(prefix) + ".pgen")
        geno = np.zeros((4, n_samples), dtype=np.int8)
        writer = pgenlib.PgenWriter(str(pgen_path).encode(), n_samples, 4)
        try:
            writer.append_biallelic_batch(geno)
        finally:
            writer.close()

        # No explicit n_pvar — exercises the count_raw_variants path
        # against the .zst sibling.
        assert check_pvar_pgen_row_count_consistent(pgen_path) == 4
