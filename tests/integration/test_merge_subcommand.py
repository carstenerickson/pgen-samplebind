"""Integration tests for the `merge` subcommand.

Day 3 minimal coverage: smoke test (two panels with shared variants but
disjoint samples concatenate correctly), --on-collision first behavior,
.psam output structure. HLD tests 3-6 (strand flip recovery, allele swap
recovery, ambiguous strand, missing variant default) wire later when the
genotype-modifier helpers land.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pgenlib
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _read_pgen_full(pgen_path: Path, n_samples: int, n_variants: int) -> np.ndarray:
    """Read a complete .pgen into an (n_variants, n_samples) int8 matrix."""
    reader = pgenlib.PgenReader(str(pgen_path).encode(), raw_sample_ct=n_samples)
    try:
        buf = np.empty((n_variants, n_samples), dtype=np.int8)
        reader.read_range(0, n_variants, buf, sample_maj=0)
        return buf
    finally:
        if hasattr(reader, "close"):
            reader.close()


def _read_psam_iids(psam_path: Path) -> list[str]:
    import pandas as pd

    df = pd.read_csv(psam_path, sep="\t")
    iid_col = next(c for c in df.columns if c.lstrip("#") == "IID")
    return df[iid_col].tolist()


def _read_pvar_keys(pvar_path: Path) -> list[tuple[int, int]]:
    import pandas as pd

    df = pd.read_csv(pvar_path, sep="\t")
    chrom_col = next(c for c in df.columns if c.lstrip("#") == "CHROM")
    return list(zip(df[chrom_col].astype(int), df["POS"].astype(int), strict=True))


@pytest.fixture
def panel_a(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Panel A: 10 samples (A00000..A00009), 100 variants from variant_seed=1."""
    out_dir = tmp_path_factory.mktemp("panel_a")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=100,
        n_populations=2,
        variant_seed=1,
        sample_seed=10,
        sample_id_prefix="A",
    )
    desc = synthesize_pfile(spec, out_dir / "a")
    return desc.path


@pytest.fixture
def panel_b(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Panel B: 10 samples (B00000..B00009), 100 variants from variant_seed=1
    (same variants as panel_a)."""
    out_dir = tmp_path_factory.mktemp("panel_b")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=100,
        n_populations=2,
        variant_seed=1,  # same variant set
        sample_seed=20,
        sample_id_prefix="B",
    )
    desc = synthesize_pfile(spec, out_dir / "b")
    return desc.path


@pytest.fixture
def panel_c(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Panel C: 10 samples (C00000..C00009), same 100 variants as A and B."""
    out_dir = tmp_path_factory.mktemp("panel_c")
    spec = SyntheticPanelSpec(
        n_samples=10,
        n_variants=100,
        n_populations=2,
        variant_seed=1,
        sample_seed=30,
        sample_id_prefix="C",
    )
    desc = synthesize_pfile(spec, out_dir / "c")
    return desc.path


class TestMultiInputMerge:
    """HLD project plan Day 4: multi-input support (3+ PFILEs in one bind)."""

    def test_three_panels_concatenate(
        self, panel_a: Path, panel_b: Path, panel_c: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel_a),
                str(panel_b),
                str(panel_c),
                "-o",
                str(out),
                "--trust-strand",  # synth panels share variant_seed; same-source
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        iids = _read_psam_iids(Path(str(out) + ".psam"))
        assert len(iids) == 30  # 10 + 10 + 10
        assert set(iids) == (
            {f"A{i:05d}" for i in range(10)}
            | {f"B{i:05d}" for i in range(10)}
            | {f"C{i:05d}" for i in range(10)}
        )

        # Variants intact (all-passthrough across both non-canonical inputs)
        keys = _read_pvar_keys(Path(str(out) + ".pvar"))
        assert len(keys) == 100

        # Pgen has the right shape
        buf = _read_pgen_full(Path(str(out) + ".pgen"), n_samples=30, n_variants=100)
        assert buf.shape == (100, 30)


class TestMergeSmoke:
    def test_two_disjoint_panels_concatenate(
        self, panel_a: Path, panel_b: Path, tmp_path: Path
    ) -> None:
        """Two panels with disjoint samples, same variants → 20 samples, 100 variants."""
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel_a),
                str(panel_b),
                "-o",
                str(out),
                "--trust-strand",  # synth panels share variant_seed; same-source
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        out_pgen = Path(str(out) + ".pgen")
        out_pvar = Path(str(out) + ".pvar")
        out_psam = Path(str(out) + ".psam")
        assert out_pgen.exists() and out_pvar.exists() and out_psam.exists()

        # 20 output samples (10 + 10, no collisions)
        iids = _read_psam_iids(out_psam)
        assert len(iids) == 20
        assert set(iids) == {f"A{i:05d}" for i in range(10)} | {f"B{i:05d}" for i in range(10)}

        # 100 output variants (full intersection, no drops since same variant_seed)
        keys = _read_pvar_keys(out_pvar)
        assert len(keys) == 100

        # Pgen has the right shape
        buf = _read_pgen_full(out_pgen, n_samples=20, n_variants=100)
        assert buf.shape == (100, 20)

    def test_round_trip_via_on_collision_first(self, panel_a: Path, tmp_path: Path) -> None:
        """Merge a panel with itself using --on-collision first → output equals input."""
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(panel_a),
                str(panel_a),
                "-o",
                str(out),
                "--on-collision",
                "first",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        # Output sample count == panel_a sample count (duplicates dropped)
        iids = _read_psam_iids(Path(str(out) + ".psam"))
        assert len(iids) == 10
        assert set(iids) == {f"A{i:05d}" for i in range(10)}


class TestMergePsamStructure:
    def test_psam_has_canonical_columns(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        out = tmp_path / "merged"
        runner = CliRunner()
        runner.invoke(cli, ["merge", str(panel_a), str(panel_b), "-o", str(out), "--quiet"])

        import pandas as pd

        df = pd.read_csv(Path(str(out) + ".psam"), sep="\t")
        cols = [c.lstrip("#") for c in df.columns]
        # Canonical column order per LLD §3.5: FID first (plink2 .psam spec
        # requires FID to be the first column when present, header `#FID...`).
        assert cols[:5] == ["FID", "IID", "SEX", "POP", "PSEUDOHAPLOID"]

    def test_fid_equals_pop(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        """HLD §Output PFILE: FID populated equal to POP for AT2 compatibility."""
        out = tmp_path / "merged"
        runner = CliRunner()
        runner.invoke(cli, ["merge", str(panel_a), str(panel_b), "-o", str(out), "--quiet"])

        import pandas as pd

        df = pd.read_csv(Path(str(out) + ".psam"), sep="\t")
        df.columns = [c.lstrip("#") for c in df.columns]
        assert (df["FID"] == df["POP"]).all()

    def test_pseudohaploid_column_populated(
        self, panel_a: Path, panel_b: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "merged"
        runner = CliRunner()
        runner.invoke(cli, ["merge", str(panel_a), str(panel_b), "-o", str(out), "--quiet"])

        import pandas as pd

        df = pd.read_csv(Path(str(out) + ".psam"), sep="\t")
        df.columns = [c.lstrip("#") for c in df.columns]
        assert "PSEUDOHAPLOID" in df.columns
        # Each value must be one of {0, 1, U}
        assert set(df["PSEUDOHAPLOID"].astype(str).unique()).issubset({"0", "1", "U"})


class TestMergeOnCollisionError:
    def test_collision_raises_under_default_policy(self, panel_a: Path, tmp_path: Path) -> None:
        """Default --on-collision error → InvariantViolation when both inputs share IIDs.

        CliRunner doesn't go through cli.main() (the exit-code mapping wrapper),
        so we assert exception type rather than exit_code 3. End-to-end exit
        code is exercised by tests/integration/test_cli_main_exit_codes.py
        (Day 4 or later).
        """
        from pgen_samplebind.errors import InvariantViolation

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["merge", str(panel_a), str(panel_a), "-o", str(out), "--quiet"]
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, InvariantViolation)
        assert "iid" in str(result.exception).lower()


class TestMergeRebindOwnOutput:
    """Issue #6: re-binding pgen-samplebind's own output with --population-column
    FID used to crash because the prior-output's POP column collided with the
    rename target.
    """

    def test_rebind_output_with_population_column_fid(
        self, panel_a: Path, panel_b: Path, tmp_path: Path
    ) -> None:
        out_1 = tmp_path / "merged_1"
        runner = CliRunner()
        result_1 = runner.invoke(
            cli,
            [
                "merge",
                str(panel_a),
                str(panel_b),
                "-o",
                str(out_1),
                "--trust-strand",
                "--quiet",
            ],
        )
        assert result_1.exit_code == 0, result_1.output

        out_2 = tmp_path / "merged_2"
        result_2 = runner.invoke(
            cli,
            [
                "merge",
                str(out_1),
                str(out_1),
                "-o",
                str(out_2),
                "--population-column",
                "FID",
                "--on-collision",
                "first",
                "--trust-strand",
                "--quiet",
            ],
        )
        assert result_2.exit_code == 0, result_2.output

        import pandas as pd

        df = pd.read_csv(Path(str(out_2) + ".psam"), sep="\t")
        df.columns = [c.lstrip("#") for c in df.columns]
        assert list(df.columns).count("POP") == 1
        assert (df["FID"] == df["POP"]).all()
        assert len(df) == 20


class TestMergeQuietAndOutput:
    def test_default_emits_summary(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        """Without --quiet, the simple summary block prints to stdout."""
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(cli, ["merge", str(panel_a), str(panel_b), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert "Read 2 inputs" in result.output
        assert "Wrote" in result.output

    def test_quiet_suppresses_summary(self, panel_a: Path, panel_b: Path, tmp_path: Path) -> None:
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["merge", str(panel_a), str(panel_b), "-o", str(out), "--quiet"]
        )
        assert result.exit_code == 0, result.output
        assert "Read 2 inputs" not in result.output
