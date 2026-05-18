"""`afs` subcommand: PFILE → AT2-compatible AFS TSVs.

Bridges the PFILE-native pipeline gap until `pfile_to_afs()` lands in AT2
upstream. Tests verify per-population aggregation math + pseudohaploid
adjustment + the three-TSV-plus-manifest output shape.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from pgen_samplebind.afs import compute_afs
from pgen_samplebind.cli import cli
from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.formats import prepared_input
from tests.fixtures.helpers import read_pgen_full as _read_pgen_full
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


class TestAfsAggregationCorrectness:
    """Verify per-population allele-frequency math against a hand-computed
    reference. Synth panel is fully autosomal so include_chrom doesn't kick
    in; pseudohaploid_fraction=0 so all-diploid simplifies the math."""

    def test_freq_matches_hand_computed(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=12,
            n_variants=10,
            n_populations=3,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=1,
            sample_seed=2,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")

        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            result = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                adjust_pseudohaploid=False,
            )

        # Hand-compute the reference by reading PFILE directly.
        psam = pd.read_csv(desc.psam_path, sep="\t")
        psam.columns = [c.lstrip("#") for c in psam.columns]
        pops = psam["POP"].astype(str).to_numpy()
        pop_labels = sorted(set(pops.tolist()))

        # PFILE genotype convention: 0/1/2 = count of ALT; -9 = missing.
        # All-diploid + no missing → expected ALT freq per (variant, pop)
        # is just mean(g) / 2 over samples in pop.
        geno = _read_pgen_full(desc.path, 12, 10)  # (variants, samples)
        for v_idx in range(10):
            for p in pop_labels:
                sample_idx = np.where(pops == p)[0]
                values = geno[v_idx, sample_idx]
                expected_alt_count = float(values.sum())
                expected_called = int(2 * len(values))  # diploid, no missing
                expected_freq = expected_alt_count / expected_called
                # Find the corresponding row in the result
                row = result.freq.iloc[v_idx]
                row_counts = result.counts.iloc[v_idx]
                assert row[p] == pytest.approx(expected_freq), (
                    f"variant {v_idx}, pop {p}: expected freq {expected_freq}, got {row[p]}"
                )
                assert row_counts[p] == expected_called, (
                    f"variant {v_idx}, pop {p}: expected called {expected_called}, "
                    f"got {row_counts[p]}"
                )

    def test_missing_data_reduces_called_count(self, tmp_path: Path) -> None:
        """With missing genotypes, called_count should drop accordingly and
        freq should be computed only over called samples."""
        spec = SyntheticPanelSpec(
            n_samples=10,
            n_variants=5,
            n_populations=2,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.3,  # 30% missingness
            variant_seed=5,
            sample_seed=6,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")

        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            result = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                adjust_pseudohaploid=False,
            )

        # Cross-check against the raw genotype matrix.
        psam = pd.read_csv(desc.psam_path, sep="\t")
        psam.columns = [c.lstrip("#") for c in psam.columns]
        pops = psam["POP"].astype(str).to_numpy()
        pop_labels = sorted(set(pops.tolist()))
        geno = _read_pgen_full(desc.path, 10, 5)

        for v_idx in range(5):
            for p in pop_labels:
                sample_idx = np.where(pops == p)[0]
                values = geno[v_idx, sample_idx]
                non_missing = values != -9
                expected_called = int(2 * non_missing.sum())
                expected_alt = float(values[non_missing].sum())
                row_counts = result.counts.iloc[v_idx]
                row_freq = result.freq.iloc[v_idx]
                assert row_counts[p] == expected_called
                if expected_called > 0:
                    assert row_freq[p] == pytest.approx(expected_alt / expected_called)
                else:
                    assert pd.isna(row_freq[p])


class TestAfsPseudohaploidAdjustment:
    """With pseudohaploid samples in the input, called-allele counts should
    halve (1 per pseudohaploid sample instead of 2) and ALT counts should
    likewise scale by 1/2 — leaving the frequency UNCHANGED."""

    def test_pseudohap_halves_called_counts(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=10,
            n_variants=5,
            n_populations=2,
            pseudohaploid_fraction=1.0,  # ALL pseudohaploid
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=11,
            sample_seed=12,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")

        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            adjusted = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                adjust_pseudohaploid=True,
            )
            unadjusted = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                adjust_pseudohaploid=False,
            )

        # Adjusted: each pseudohap sample contributes 1 called allele.
        # Unadjusted: contributes 2. So adjusted counts == unadjusted / 2.
        adj_counts = adjusted.counts.drop(columns="variant_id").to_numpy()
        unadj_counts = unadjusted.counts.drop(columns="variant_id").to_numpy()
        np.testing.assert_array_equal(adj_counts * 2, unadj_counts)

        # Frequency unchanged (numerator and denominator both halve).
        adj_freq = adjusted.freq.drop(columns="variant_id").to_numpy()
        unadj_freq = unadjusted.freq.drop(columns="variant_id").to_numpy()
        np.testing.assert_allclose(adj_freq, unadj_freq, rtol=1e-10)

    def test_adjust_applied_flag_reflects_input(self, tmp_path: Path) -> None:
        """All-diploid input → adjust_pseudohaploid_applied should be False
        even if the flag is True, because no pseudohap samples were present."""
        spec = SyntheticPanelSpec(
            n_samples=8,
            n_variants=3,
            n_populations=2,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=21,
            sample_seed=22,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            result = compute_afs(pdesc, adjust_pseudohaploid=True)
        assert result.adjust_pseudohaploid_applied is False


class TestAfsPopulationSubset:
    def test_populations_filter_restricts_output_columns(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=15,
            n_variants=10,
            n_populations=5,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=31,
            sample_seed=32,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")

        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            result = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                populations=["pop_00", "pop_02"],
            )

        # Only the two requested pops appear as columns.
        freq_cols = [c for c in result.freq.columns if c != "variant_id"]
        assert freq_cols == ["pop_00", "pop_02"]
        assert result.populations == ["pop_00", "pop_02"]

    def test_unknown_population_raises(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=6,
            n_variants=3,
            n_populations=2,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=41,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        with (
            prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc,
            pytest.raises(InvariantViolation, match="not found"),
        ):
            compute_afs(
                descriptor=pdesc,
                populations=["nonexistent_population"],
            )

    def test_missing_population_column_raises(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=4,
            n_variants=3,
            n_populations=1,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        with (
            prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc,
            pytest.raises(InvariantViolation, match="not in"),
        ):
            compute_afs(
                descriptor=pdesc,
                population_column="DOES_NOT_EXIST",
            )


class TestAfsSubcommandEndToEnd:
    def test_three_tsvs_and_manifest_written(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=12,
            n_variants=20,
            n_populations=3,
            pseudohaploid_fraction=0.5,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.05,
            variant_seed=51,
            sample_seed=52,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        out_dir = tmp_path / "afs_out"

        runner = CliRunner()
        result = runner.invoke(cli, ["afs", str(desc.path), "-o", str(out_dir), "--quiet"])
        assert result.exit_code == 0, result.output

        # Four files: 3 TSVs + manifest
        files = {p.name for p in out_dir.iterdir()}
        assert files == {
            "afs_snp.tsv",
            "afs_freq.tsv",
            "afs_counts.tsv",
            "afs_manifest.json",
        }

        snp = pd.read_csv(out_dir / "afs_snp.tsv", sep="\t")
        assert list(snp.columns) == ["variant_id", "chrom", "pos", "ref", "alt", "cm"]
        assert len(snp) == 20

        freq = pd.read_csv(out_dir / "afs_freq.tsv", sep="\t")
        assert len(freq) == 20
        assert "variant_id" in freq.columns
        # 3 populations + variant_id
        assert len(freq.columns) == 4

        manifest = json.loads((out_dir / "afs_manifest.json").read_text())
        assert manifest["command"] == "afs"
        assert manifest["n_variants"] == 20
        assert manifest["n_populations"] == 3
        # Pseudohaploid mix present → adjustment applied.
        assert manifest["adjust_pseudohaploid_applied"] is True
        assert sum(manifest["n_samples_per_pop"].values()) == 12

    def test_no_pseudohaploid_adjust_flag(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=10,
            n_variants=10,
            n_populations=2,
            pseudohaploid_fraction=1.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=61,
            sample_seed=62,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        out_dir = tmp_path / "afs_out"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "afs",
                str(desc.path),
                "-o",
                str(out_dir),
                "--no-pseudohaploid-adjust",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        manifest = json.loads((out_dir / "afs_manifest.json").read_text())
        assert manifest["adjust_pseudohaploid_requested"] is False
        assert manifest["adjust_pseudohaploid_applied"] is False

    def test_population_filter_via_cli(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=15,
            n_variants=10,
            n_populations=5,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=71,
            sample_seed=72,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        out_dir = tmp_path / "afs_out"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "afs",
                str(desc.path),
                "-o",
                str(out_dir),
                "--populations",
                "pop_00",
                "--populations",
                "pop_02",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        manifest = json.loads((out_dir / "afs_manifest.json").read_text())
        assert manifest["populations"] == ["pop_00", "pop_02"]

        freq = pd.read_csv(out_dir / "afs_freq.tsv", sep="\t")
        assert [c for c in freq.columns if c != "variant_id"] == ["pop_00", "pop_02"]


class TestAfsIncludeSexChrom:
    """The `--include-sex-chrom` flag is the only sex-chromosome code path
    in the tool. Without it, autosomes (1-22) only; with it, X/Y/XY/MT
    (23-26) are included alongside. Previously untested at integration
    level despite being the documented option."""

    def test_autosomes_only_excludes_chr23(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=6,
            n_variants=12,
            n_populations=2,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            chromosomes=(1, 2, 23),  # mix autosomes + X
            variant_seed=991,
            sample_seed=992,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        out_dir = tmp_path / "afs_out"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["afs", str(desc.path), "-o", str(out_dir), "--quiet"]
        )
        assert result.exit_code == 0, result.output

        snp = pd.read_csv(out_dir / "afs_snp.tsv", sep="\t")
        # chrom 23 dropped → only autosomes remain
        assert set(snp["chrom"].astype(int)) <= {1, 2}
        assert 23 not in set(snp["chrom"].astype(int))

    def test_include_sex_chrom_admits_chr23(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=6,
            n_variants=12,
            n_populations=2,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            chromosomes=(1, 2, 23),
            variant_seed=991,
            sample_seed=992,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        out_dir = tmp_path / "afs_out"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["afs", str(desc.path), "-o", str(out_dir), "--include-sex-chrom", "--quiet"],
        )
        assert result.exit_code == 0, result.output

        snp = pd.read_csv(out_dir / "afs_snp.tsv", sep="\t")
        # chrom 23 retained when sex chroms included
        chrom_set = set(snp["chrom"].astype(int))
        assert 23 in chrom_set
