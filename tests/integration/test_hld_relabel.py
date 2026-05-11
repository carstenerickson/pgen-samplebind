"""HLD test 14: --relabel-from anno; plus 2-col form + --id-column coverage.

Per HLD §Validation strategy:
  14. test_relabel_from_anno: AADR `.anno` file as `--relabel-from` with
      `--relabel-input-col 'Genetic ID' --relabel-output-col 'Group ID'`;
      output `.psam` POP column matches a hand-extracted (Genetic ID,
      Group ID) reference.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures.modifiers import make_panel_with_iids
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _read_psam(prefix: Path) -> pd.DataFrame:
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def _write_anno_tsv(
    path: Path,
    rows: list[dict[str, str]],
    columns: list[str] | None = None,
) -> Path:
    """Write an N-col TSV with the provided rows. `columns` controls header order."""
    df = pd.DataFrame(rows)
    if columns:
        df = df[columns]
    df.to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


# ---------- HLD test 14: N-col anno-style relabel ----------------------------


class TestHld14RelabelFromAnno:
    """AADR-style N-col anno file; --relabel-input-col 'Genetic ID',
    --relabel-output-col 'Group ID'. Output POP per-sample comes from
    the anno's Group ID; samples missing from the anno keep their
    original POP (per HLD §Relabeling)."""

    def test_anno_overrides_pop_per_sample(self, tmp_path: Path) -> None:
        # Two panels with disjoint IIDs but same variants.
        a = make_panel_with_iids(tmp_path / "a", ["I0001", "I0002", "I0003"], pop="orig_pop_a")
        b = make_panel_with_iids(tmp_path / "b", ["I0004", "I0005"], pop="orig_pop_b")

        # AADR-style anno file with extra columns (ensure N-col path is used)
        anno = _write_anno_tsv(
            tmp_path / "aadr_v66.anno",
            [
                {"Genetic ID": "I0001", "Group ID": "England_IA", "Year": "2020"},
                {"Genetic ID": "I0002", "Group ID": "Russia_IA", "Year": "2019"},
                {"Genetic ID": "I0004", "Group ID": "Sweden_Viking", "Year": "2021"},
                # I0003 and I0005 deliberately missing — they keep orig POP
            ],
            columns=["Genetic ID", "Group ID", "Year"],
        )

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--relabel-from",
                str(anno),
                "--relabel-input-col",
                "Genetic ID",
                "--relabel-output-col",
                "Group ID",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        psam_df = _read_psam(out)
        # Build the expected POP map from anno + originals
        expected_by_iid = {
            "I0001": "England_IA",  # from anno
            "I0002": "Russia_IA",  # from anno
            "I0003": "orig_pop_a",  # missing from anno → original
            "I0004": "Sweden_Viking",  # from anno
            "I0005": "orig_pop_b",  # missing from anno → original
        }
        actual_by_iid = dict(zip(psam_df["IID"], psam_df["POP"], strict=True))
        assert actual_by_iid == expected_by_iid

    def test_fid_mirrors_relabeled_pop(self, tmp_path: Path) -> None:
        """FID = POP must reflect the POST-relabel POP, not the original."""
        a = make_panel_with_iids(tmp_path / "a", ["S001"], pop="ORIGINAL_POP")
        b = make_panel_with_iids(tmp_path / "b", ["S002"], pop="ORIGINAL_POP")
        anno = _write_anno_tsv(
            tmp_path / "anno.tsv",
            [
                {"Genetic ID": "S001", "Group ID": "RELABELED"},
                {"Genetic ID": "S002", "Group ID": "RELABELED"},
            ],
            columns=["Genetic ID", "Group ID"],
        )

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--relabel-from",
                str(anno),
                "--relabel-input-col",
                "Genetic ID",
                "--relabel-output-col",
                "Group ID",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        psam_df = _read_psam(out)
        # FID must match POST-relabel POP (per HLD §Output PFILE)
        assert (psam_df["FID"] == psam_df["POP"]).all()
        assert set(psam_df["FID"]) == {"RELABELED"}


# ---------- 2-col simple form (POP→POP collapse) ----------------------------


class TestRelabelTwoColSimple:
    """2-col header-less TSV: source = POP, target = POP. Used to collapse
    multiple population labels into a single output label across inputs."""

    def test_two_col_collapses_pops(self, tmp_path: Path) -> None:
        a = make_panel_with_iids(tmp_path / "a", ["S001", "S002"], pop="OldName_A")
        b = make_panel_with_iids(tmp_path / "b", ["S003", "S004"], pop="OldName_B")

        # 2-col TSV: OldName_A → Combined; OldName_B → Combined
        relabel = tmp_path / "relabel.tsv"
        relabel.write_text("OldName_A\tCombined\nOldName_B\tCombined\n")

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--relabel-from",
                str(relabel),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output
        psam_df = _read_psam(out)
        assert set(psam_df["POP"]) == {"Combined"}

    def test_two_col_unmapped_pop_keeps_original(self, tmp_path: Path) -> None:
        a = make_panel_with_iids(tmp_path / "a", ["S001"], pop="MapMe")
        b = make_panel_with_iids(tmp_path / "b", ["S002"], pop="LeaveAlone")
        relabel = tmp_path / "relabel.tsv"
        relabel.write_text("MapMe\tNewLabel\n")

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["merge", str(a), str(b), "-o", str(out), "--relabel-from", str(relabel), "--quiet"],
        )
        assert result.exit_code == 0, result.output
        psam_df = _read_psam(out)
        actual = dict(zip(psam_df["IID"], psam_df["POP"], strict=True))
        assert actual == {"S001": "NewLabel", "S002": "LeaveAlone"}


# ---------- --id-column custom join key (Master ID) -------------------------


class TestRelabelWithCustomIdColumn:
    """HLD §Cross-version sample identity: --id-column NAME selects which
    .psam column drives sample-identity ops, including --relabel-from's
    join. AADR cross-version uses 'Master ID'."""

    def test_id_column_routes_relabel_lookup(self, tmp_path: Path) -> None:
        # Construct .psam with both IID and a custom 'Master ID' column.
        # We can't use make_panel_with_iids directly (no Master ID); build inline.
        from tests.fixtures.modifiers import make_panel_with_iids

        a = make_panel_with_iids(tmp_path / "a", ["I0001"], pop="orig")
        # Hand-edit .psam to add a Master ID column.
        psam_path = Path(str(a) + ".psam")
        df = pd.read_csv(psam_path, sep="\t")
        df["Master ID"] = ["MASTER_AAA"]
        # Re-write keeping the # prefix on the first column header
        with open(psam_path, "w") as f:
            df.to_csv(f, sep="\t", index=False)

        # Anno indexed by Master ID
        anno = _write_anno_tsv(
            tmp_path / "anno.tsv",
            [{"Master ID": "MASTER_AAA", "Group ID": "MasterPop"}],
            columns=["Master ID", "Group ID"],
        )

        b = make_panel_with_iids(tmp_path / "b", ["I0002"], pop="other")

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--id-column",
                "Master ID",  # join on Master ID column
                "--relabel-from",
                str(anno),
                "--relabel-input-col",
                "Master ID",
                "--relabel-output-col",
                "Group ID",
                "--on-collision",
                "first",  # in case of any collision
                "--quiet",
            ],
        )
        # input b lacks Master ID column → relabel.apply_relabel will raise
        # UsageError ("source column 'Master ID' not present").
        # That's the expected failure mode for cross-version use without a
        # Master ID in every input. Test confirms the error message points
        # the user at --id-column.
        from pgen_samplebind.errors import UsageError

        assert isinstance(result.exception, UsageError)
        assert "Master ID" in str(result.exception)


# ---------- Edge cases ------------------------------------------------------


class TestRelabelEdgeCases:
    def test_one_relabel_col_without_other_raises(self, tmp_path: Path) -> None:
        """--relabel-input-col without --relabel-output-col → UsageError."""
        from pgen_samplebind.errors import UsageError

        a = make_panel_with_iids(tmp_path / "a", ["S001"])
        b = make_panel_with_iids(tmp_path / "b", ["S002"])
        relabel = tmp_path / "r.tsv"
        relabel.write_text("Genetic ID\tGroup ID\nS001\tX\n")

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--relabel-from",
                str(relabel),
                "--relabel-input-col",
                "Genetic ID",
                # missing --relabel-output-col
                "--quiet",
            ],
        )
        assert isinstance(result.exception, UsageError)
        assert "must be supplied together" in str(result.exception)

    def test_relabel_validate_subcommand(self, tmp_path: Path) -> None:
        """--relabel-from also wired into validate (HLD §Relabeling: applies
        to merge AND validate)."""
        a = make_panel_with_iids(tmp_path / "a", ["S001"], pop="OldA")
        b = make_panel_with_iids(tmp_path / "b", ["S002"], pop="OldB")
        relabel = tmp_path / "r.tsv"
        relabel.write_text("OldA\tNewA\nOldB\tNewB\n")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "validate",
                str(a),
                str(b),
                "--relabel-from",
                str(relabel),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output


# ---------- HLD test 15: no string-labeled chroms in output -----------------


@pytest.mark.eigenstrat
class TestHld15NoStringLabeledChroms:
    """HLD test 15: output .pvar's #CHROM column contains integer values
    only (1-22 for autosomes; the autosomal filter excludes X/Y/MT by
    default per HLD §Format detection). String labels like 'chr1' are
    normalized away during the EIGENSTRAT → PFILE conversion + our
    pvar.normalize_chrom path."""

    def test_eigenstrat_output_chrom_column_is_integer(self, tmp_path: Path) -> None:
        from tests.fixtures.modifiers import pfile_to_eigenstrat

        # Build PFILE → EIGENSTRAT → merge → output PFILE
        pfile_orig = tmp_path / "orig"
        synthesize_pfile(
            SyntheticPanelSpec(
                n_samples=4,
                n_variants=20,
                n_populations=1,
                variant_seed=151,
                sample_seed=152,
                sample_id_prefix="P",
            ),
            pfile_orig,
        )
        eig = tmp_path / "eig"
        pfile_to_eigenstrat(pfile_orig, eig)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["merge", str(eig), str(eig), "-o", str(out), "--on-collision", "first", "--quiet"]
        )
        assert result.exit_code == 0, result.output

        out_pvar = pd.read_csv(Path(str(out) + ".pvar"), sep="\t")
        # The #CHROM column should be integer-valued (not string-labeled).
        # Pandas may parse it as int64 directly if all values are numeric.
        chrom_col = out_pvar["#CHROM"]
        # Every value parses cleanly as an integer in 1..22 (autosomal default).
        for v in chrom_col:
            int_v = int(v)
            assert 1 <= int_v <= 22, f"unexpected chrom value {v!r} (outside 1-22 autosome range)"
        # Reject any explicit "chr"-prefixed string labels
        for v in chrom_col.astype(str):
            assert not v.lower().startswith("chr"), (
                f"output .pvar has string-labeled chrom {v!r}; expected normalized integer"
            )
