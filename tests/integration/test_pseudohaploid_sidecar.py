"""Integration tests for the pseudohaploid-sidecar override path (issue #2).

The sidecar is the authoritative per-sample pseudohaploid declaration from
upstream tools that know the answer by methodology (e.g., pileup-aadr's
single-BAM --randomDiploid → pseudohaploid by construction). When present,
it takes precedence over both the input `.psam` PSEUDOHAPLOID column AND
the heterozygosity-derived inference.

Coverage:
  - merge: override flips DIPLOID-by-genotype → PSEUDOHAPLOID in output .psam
  - merge: partial sidecar coverage (some IIDs covered, others fall back)
  - merge: orphan IID in sidecar → InvariantViolation (exit 3)
  - merge: regression — no sidecar → output unchanged (auto-detect path)
  - afs: sidecar drives adjust_pseudohaploid_applied=True even without .psam column
  - afs: orphan IID → InvariantViolation
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from pgen_samplebind.afs import compute_afs
from pgen_samplebind.cli import cli
from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.formats import prepared_input
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _read_psam_pseudohaploid(prefix: Path) -> dict[str, str]:
    """Return {IID: PSEUDOHAPLOID} from a .psam file."""
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    return dict(zip(df["IID"].astype(str), df["PSEUDOHAPLOID"].astype(str), strict=True))


def _write_sidecar(prefix: Path, samples: dict[str, int]) -> None:
    """Write a v1 sidecar mapping IID → 0/1 next to the given prefix."""
    payload = {
        "schema_version": 1,
        "samples": {iid: {"pseudohaploid": v} for iid, v in samples.items()},
    }
    Path(str(prefix) + ".pseudohaploid.json").write_text(json.dumps(payload))


def _diploid_panel(tmp_path: Path, name: str, *, sample_seed: int, prefix: str) -> Path:
    """Synthesize an all-DIPLOID panel (pseudohaploid_fraction=0.0)."""
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=200,
        n_populations=2,
        variant_seed=901,
        sample_seed=sample_seed,
        sample_id_prefix=prefix,
        pseudohaploid_fraction=0.0,  # all DIPLOID by genotype
    )
    return synthesize_pfile(spec, tmp_path / name).path


class TestMergeSidecarOverride:
    def test_sidecar_flips_diploid_to_pseudohaploid(self, tmp_path: Path) -> None:
        """Sidecar declares samples PSEUDOHAPLOID; output .psam matches the
        sidecar even though genotype heterozygosity would say DIPLOID."""
        a = _diploid_panel(tmp_path, "a", sample_seed=101, prefix="A")
        b = _diploid_panel(tmp_path, "b", sample_seed=102, prefix="B")

        # All input[0] (panel A) samples declared pseudohaploid via sidecar.
        a_iids = list(_read_psam_pseudohaploid(a).keys())
        _write_sidecar(a, {iid: 1 for iid in a_iids})
        # Panel B: no sidecar → falls back to heterozygosity inference (DIPLOID).

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), str(b), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid(out)
        # Panel A samples → PSEUDOHAPLOID (per sidecar override).
        for iid in a_iids:
            assert detected[iid] == "1", f"expected sidecar override for {iid}, got {detected[iid]}"
        # Panel B samples → DIPLOID (per heterozygosity auto-detect).
        b_iids = list(_read_psam_pseudohaploid(b).keys())
        for iid in b_iids:
            assert detected[iid] == "0", (
                f"expected auto-detect DIPLOID for {iid}, got {detected[iid]}"
            )

    def test_sidecar_partial_coverage(self, tmp_path: Path) -> None:
        """Sidecar with only a subset of IIDs: covered samples follow sidecar,
        uncovered fall back to heterozygosity inference."""
        a = _diploid_panel(tmp_path, "a", sample_seed=103, prefix="A")
        a_iids = list(_read_psam_pseudohaploid(a).keys())

        # First half declared PSEUDOHAPLOID via sidecar; rest absent.
        half = len(a_iids) // 2
        sidecar_overrides = {iid: 1 for iid in a_iids[:half]}
        _write_sidecar(a, sidecar_overrides)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid(out)
        for iid in a_iids[:half]:
            assert detected[iid] == "1", f"sidecar override missed {iid}"
        for iid in a_iids[half:]:
            assert detected[iid] == "0", f"expected DIPLOID auto-detect for {iid}"

    def test_sidecar_orphan_iid_raises(self, tmp_path: Path) -> None:
        """Sidecar with an IID that doesn't appear in the .psam → InvariantViolation
        (exits 3 via cli.main; CliRunner doesn't go through that path, so we
        check the raised exception type — same pattern as test_merge_subcommand
        for collision-error exits)."""
        a = _diploid_panel(tmp_path, "a", sample_seed=104, prefix="A")
        a_iids = list(_read_psam_pseudohaploid(a).keys())
        _write_sidecar(a, {a_iids[0]: 1, "GHOST_SAMPLE": 1})

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), "-o", str(out), "--quiet"])
        assert isinstance(result.exception, InvariantViolation), result.output
        assert "GHOST_SAMPLE" in str(result.exception)

    def test_no_sidecar_preserves_auto_detection(self, tmp_path: Path) -> None:
        """Regression: with no sidecar, the orchestrator's heterozygosity
        re-derivation continues to drive the output PSEUDOHAPLOID column."""
        a = _diploid_panel(tmp_path, "a", sample_seed=105, prefix="A")

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(a), "-o", str(out), "--quiet"])
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid(out)
        # All DIPLOID by synthesis → all auto-classified DIPLOID.
        assert set(detected.values()) == {"0"}

    def test_sidecar_with_target_mode_suffix_rename(self, tmp_path: Path) -> None:
        """When --on-collision suffix renames a sidecar'd sample (e.g., to
        `<iid>_target`), the override must follow the rename to the output IID."""
        a = _diploid_panel(tmp_path, "a", sample_seed=106, prefix="A")
        # Build a target panel that shares one IID with `a` so the suffix
        # scheme kicks in.
        shared_iid = next(iter(_read_psam_pseudohaploid(a)))
        target = _diploid_panel(tmp_path, "tgt", sample_seed=107, prefix="A")
        # Sidecar on the target: declare the shared (colliding) sample
        # pseudohaploid; pgen-samplebind will rename it to `<iid>_target`.
        _write_sidecar(target, {shared_iid: 1})

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                "--target",
                str(target),
                "--on-collision",
                "suffix",
                "-o",
                str(out),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        detected = _read_psam_pseudohaploid(out)
        renamed = f"{shared_iid}_target"
        assert renamed in detected, f"expected suffix-renamed {renamed} in output"
        assert detected[renamed] == "1", (
            f"sidecar override should follow suffix rename; got {detected[renamed]} for {renamed}"
        )


class TestAfsSidecarOverride:
    def test_sidecar_drives_pseudohaploid_adjustment(self, tmp_path: Path) -> None:
        """afs reads sidecar even when `.psam` has no PSEUDOHAPLOID column;
        adjust_pseudohaploid_applied flips True and called counts halve."""
        spec = SyntheticPanelSpec(
            n_samples=8,
            n_variants=50,
            n_populations=2,
            pseudohaploid_fraction=0.0,  # all DIPLOID by genotype + .psam column "0"
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=11,
            sample_seed=22,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        iids = list(_read_psam_pseudohaploid(desc.path).keys())

        # Sidecar marks all samples PSEUDOHAPLOID.
        _write_sidecar(desc.path, {iid: 1 for iid in iids})

        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            result = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                adjust_pseudohaploid=True,
            )

        assert result.adjust_pseudohaploid_applied, (
            "sidecar should drive adjust_pseudohaploid_applied=True"
        )

        # Called-allele math: 8 samples, all pseudohap → 8 alleles per pop-row
        # in total (split across 2 pops, 4 samples each = 4 alleles per pop).
        # Sum across pops per variant should equal the per-variant total.
        # Specifically: with adjust applied, each row's pop-summed count must
        # equal n_samples_per_pop[pop] (=4), not 2 * n_samples_per_pop[pop] (=8).
        pop_cols = [c for c in result.counts.columns if c != "variant_id"]
        first_row_counts = result.counts.iloc[0][pop_cols].astype(int).tolist()
        # Each pop has 4 samples; pseudohap → 4 called alleles. Diploid would give 8.
        assert first_row_counts == [4, 4], (
            f"pseudohaploid adjustment failed: got {first_row_counts} per pop, "
            f"expected [4, 4] (1 allele per sample x 4 samples)"
        )

    def test_afs_sidecar_orphan_iid_raises(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=6,
            n_variants=20,
            pseudohaploid_fraction=0.0,
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=33,
            sample_seed=44,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")
        _write_sidecar(desc.path, {"NOT_A_REAL_SAMPLE": 1})

        with (
            prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc,
            pytest.raises(InvariantViolation, match="not present in"),
        ):
            compute_afs(descriptor=pdesc, population_column="POP")

    def test_afs_sidecar_overrides_psam_column(self, tmp_path: Path) -> None:
        """Sidecar takes precedence over an existing `.psam` PSEUDOHAPLOID
        column. Setup: synth panel writes PSEUDOHAPLOID column "0" for every
        sample; sidecar flips a subset to "1"; only the sidecar-flipped
        samples should be adjusted."""
        spec = SyntheticPanelSpec(
            n_samples=8,
            n_variants=30,
            n_populations=1,
            pseudohaploid_fraction=0.0,  # .psam PSEUDOHAPLOID written "0"
            ambiguous_strand_fraction=0.0,
            missing_rate=0.0,
            variant_seed=55,
            sample_seed=66,
        )
        desc = synthesize_pfile(spec, tmp_path / "panel")

        iids = list(_read_psam_pseudohaploid(desc.path).keys())
        half = len(iids) // 2  # 4 samples
        _write_sidecar(desc.path, {iid: 1 for iid in iids[:half]})

        with prepared_input(desc.path, include_chrom=tuple(range(1, 23))) as pdesc:
            result = compute_afs(
                descriptor=pdesc,
                population_column="POP",
                adjust_pseudohaploid=True,
            )

        # Expected: 4 pseudohap (1 allele each) + 4 diploid (2 alleles each) = 12 total.
        pop_col = next(c for c in result.counts.columns if c != "variant_id")
        # Pick a variant with no missing data — the synth has missing_rate=0.0.
        first_row = result.counts.iloc[0][pop_col]
        assert int(first_row) == 4 * 1 + 4 * 2, (
            f"expected 12 alleles (4 pseudohap + 4 diploid), got {first_row}"
        )
