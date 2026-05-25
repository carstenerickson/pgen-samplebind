"""Regression tests for issue #10: silent panel GT corruption when the
input .pvar contains biallelic non-SNP rows (e.g., indels) that
`pvar.read_pvar`'s biallelic-SNP filter drops.

Pre-fix, dropped rows shifted all downstream `canonical_idx` /
`idx_input_<i>` values relative to the actual .pgen row layout, and
genotype reads landed on the wrong .pgen rows — producing silent
dosage corruption (no error, no warning, no `--report-json` flag).
Post-fix, `read_pvar` stamps `_pgen_row` before filtering, and
`alignment.build_alignment_table` uses it instead of `np.arange`, so
the read indices remain aligned with the .pgen.

Also covers the row-count guardrail
(`pvar.check_pvar_pgen_row_count_consistent`) for mis-paired triplets,
which is a related class of silent-corruption bug.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.helpers import read_pgen_full


def _write_mixed_pfile(
    prefix: Path,
    pvar_rows: list[tuple[str, int, str, str, str]],
    geno: np.ndarray,
    iids: list[str],
) -> None:
    """Write a PFILE triplet whose .pvar interleaves SNP and non-SNP rows
    (single-row population for simplicity). geno.shape == (n_variants,
    n_samples) and must include every .pvar row (incl. the non-SNP ones)
    so the .pgen row layout matches the .pvar 1:1."""
    n_samples = len(iids)
    assert geno.shape == (len(pvar_rows), n_samples)

    pvar_df = pd.DataFrame(pvar_rows, columns=["#CHROM", "POS", "ID", "REF", "ALT"])
    pvar_df.to_csv(str(prefix) + ".pvar", sep="\t", index=False, lineterminator="\n")

    psam_df = pd.DataFrame({"#IID": iids, "SEX": [1] * n_samples, "POP": ["p0"] * n_samples})
    psam_df.to_csv(str(prefix) + ".psam", sep="\t", index=False, lineterminator="\n")

    writer = pgenlib.PgenWriter(str(prefix).encode() + b".pgen", n_samples, len(pvar_rows))
    try:
        writer.append_biallelic_batch(geno)
    finally:
        writer.close()


@pytest.fixture
def panel_with_interleaved_indels(tmp_path: Path) -> tuple[Path, np.ndarray, np.ndarray]:
    """Panel whose .pvar has SNPs and biallelic indels interleaved.

    Returns (prefix, full_geno_matrix, kept_snp_row_mask). The .pgen has
    one row per .pvar line (including indel rows); read_pvar will drop
    the indel rows from the pandas DataFrame, so the surviving canonical
    indices must remap to the original .pgen rows via `_pgen_row` for
    byte-correct reads.
    """
    prefix = tmp_path / "mixed"
    # 9 .pvar rows, biallelic in the pgenlib sense (max_allele_ct == 2)
    # but only 6 pass read_pvar's single-char ACGT REF+ALT filter.
    pvar_rows = [
        ("1", 100, "v0_snp", "A", "C"),  # SNP keep — row 0
        ("1", 200, "v1_indel", "AT", "A"),  # INDEL drop — row 1
        ("1", 300, "v2_snp", "C", "G"),  # SNP keep — row 2
        ("1", 400, "v3_indel", "CG", "C"),  # INDEL drop — row 3
        ("1", 500, "v4_snp", "A", "G"),  # SNP keep — row 4
        ("1", 600, "v5_snp", "T", "C"),  # SNP keep — row 5
        ("1", 700, "v6_indel", "A", "AGTC"),  # INDEL drop — row 6
        ("1", 800, "v7_snp", "G", "T"),  # SNP keep — row 7
        ("1", 900, "v8_snp", "C", "A"),  # SNP keep — row 8
    ]
    iids = [f"S{i}" for i in range(6)]
    rng = np.random.default_rng(42)
    geno = rng.integers(0, 3, size=(len(pvar_rows), len(iids)), dtype=np.int8)

    _write_mixed_pfile(prefix, pvar_rows, geno, iids)
    # read_pvar keeps single-char ACGT REF + ALT only. Use explicit length
    # checks + a set lookup; `r in "ACGT"` would substring-match indel rows
    # like ("CG", "C") and silently miscount the expected kept count.
    bases = {"A", "C", "G", "T"}
    kept_mask = np.array(
        [len(r[3]) == 1 and len(r[4]) == 1 and r[3] in bases and r[4] in bases for r in pvar_rows],
        dtype=bool,
    )
    return prefix, geno, kept_mask


class TestCanonicalPassthroughByteIdentity:
    """Pre-fix, the canonical (input[0]) read path mis-indexed into the
    .pgen whenever read_pvar dropped any non-SNP rows above a kept SNP.
    This test exercises that path: self-merge a panel with interleaved
    indels and assert the canonical samples come back byte-identical at
    the surviving SNP rows.
    """

    def test_self_merge_preserves_canonical_bytes(
        self, panel_with_interleaved_indels: tuple[Path, np.ndarray, np.ndarray]
    ) -> None:
        prefix, full_geno, kept_mask = panel_with_interleaved_indels
        n_samples = full_geno.shape[1]
        out = prefix.parent / "merged"

        result = CliRunner().invoke(
            cli,
            [
                "merge",
                str(prefix),
                "-o",
                str(out),
                "--trust-strand",
                "--on-collision",
                "first",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)

        n_kept = int(kept_mask.sum())
        merged_geno = read_pgen_full(out, n_samples, n_kept)
        expected = full_geno[kept_mask]
        np.testing.assert_array_equal(merged_geno, expected)


class TestNonCanonicalReadByteIdentity:
    """Pre-fix, the non-canonical (input[1+]) read path had the same
    misindex bug via `idx_input_<i>`. Bind a canonical panel + a second
    panel with the same interleaved-indel layout but distinct IIDs (so
    both contribute samples to the output, no collision policy needed).
    Output sample axis is then S0..S5 (canonical) followed by B0..B5
    (non-canonical); both halves must be byte-identical to their source
    pgens at the surviving SNP rows — that's the visible signal that
    canonical_idx and idx_input_<i> both address the right .pgen rows.
    """

    def test_two_input_merge_canonical_block_is_byte_identical(
        self, panel_with_interleaved_indels: tuple[Path, np.ndarray, np.ndarray]
    ) -> None:
        prefix, full_geno, kept_mask = panel_with_interleaved_indels
        n_samples = full_geno.shape[1]
        # Build a second panel whose .pvar/.pgen layout matches the first
        # so we exercise the non-canonical read path. Use distinct IIDs so
        # both contribute to the output sample axis.
        second_prefix = prefix.parent / "mixed_b"
        pvar_b = [
            ("1", 100, "v0_snp", "A", "C"),
            ("1", 200, "v1_indel", "AT", "A"),
            ("1", 300, "v2_snp", "C", "G"),
            ("1", 400, "v3_indel", "CG", "C"),
            ("1", 500, "v4_snp", "A", "G"),
            ("1", 600, "v5_snp", "T", "C"),
            ("1", 700, "v6_indel", "A", "AGTC"),
            ("1", 800, "v7_snp", "G", "T"),
            ("1", 900, "v8_snp", "C", "A"),
        ]
        rng = np.random.default_rng(1337)
        geno_b = rng.integers(0, 3, size=(len(pvar_b), n_samples), dtype=np.int8)
        iids_b = [f"B{i}" for i in range(n_samples)]
        _write_mixed_pfile(second_prefix, pvar_b, geno_b, iids_b)

        out = prefix.parent / "merged_two"
        result = CliRunner().invoke(
            cli,
            [
                "merge",
                str(prefix),
                str(second_prefix),
                "-o",
                str(out),
                "--trust-strand",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)

        n_kept = int(kept_mask.sum())
        merged_geno = read_pgen_full(out, n_samples * 2, n_kept)
        # First n_samples columns are the canonical block (panel A); they
        # must match panel A's SNP rows byte-for-byte.
        np.testing.assert_array_equal(merged_geno[:, :n_samples], full_geno[kept_mask])
        # Next n_samples columns are panel B; same invariant for B.
        np.testing.assert_array_equal(merged_geno[:, n_samples:], geno_b[kept_mask])


class TestRowCountGuardrail:
    """The new `check_pvar_pgen_row_count_consistent` guardrail catches a
    related class of silent-corruption bug: mis-paired triplets where
    the .pgen and .pvar have different row counts (e.g., one was
    truncated, or the two came from different make-pgen runs). Without
    the guardrail, the merge would either over- or under-read the .pgen
    and produce dosage corruption that no other check catches.
    """

    def test_pvar_has_fewer_rows_than_pgen_is_rejected(self, tmp_path: Path) -> None:
        prefix = tmp_path / "mismatched"
        n_samples = 4
        # .pvar declares 3 variants; .pgen contains 5.
        pvar_df = pd.DataFrame(
            [("1", i * 100, f"v{i}", "A", "C") for i in range(3)],
            columns=["#CHROM", "POS", "ID", "REF", "ALT"],
        )
        pvar_df.to_csv(str(prefix) + ".pvar", sep="\t", index=False, lineterminator="\n")
        psam_df = pd.DataFrame(
            {
                "#IID": [f"S{i}" for i in range(n_samples)],
                "SEX": [1] * n_samples,
                "POP": ["p0"] * n_samples,
            }
        )
        psam_df.to_csv(str(prefix) + ".psam", sep="\t", index=False, lineterminator="\n")
        geno = np.zeros((5, n_samples), dtype=np.int8)
        writer = pgenlib.PgenWriter(str(prefix).encode() + b".pgen", n_samples, 5)
        try:
            writer.append_biallelic_batch(geno)
        finally:
            writer.close()

        out = tmp_path / "merged"
        result = CliRunner().invoke(
            cli, ["merge", str(prefix), "-o", str(out), "--trust-strand", "--quiet"]
        )
        assert result.exit_code != 0
        # The exception message names both counts so a user can grep for it.
        exc_msg = str(result.exception)
        assert "variant_ct (5)" in exc_msg
        assert "data-line count (3)" in exc_msg
