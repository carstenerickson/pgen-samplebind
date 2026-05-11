"""HLD test 17: f2 parity vs `mergeit` reference.

Bind the same EIGENSTRAT inputs via (a) `mergeit` (AT2's reference
EIGENSOFT-style merger) and (b) `pgen-samplebind merge`. Run AT2's
`extract_f2` over both outputs and assert per-pair `max_dev < 1e-9`.

This test is marked `external_tool` + `slow`: it is skipped in the default
CI matrix and only runs in the nightly workflow against a custom container
that has both `mergeit` (from AT2 distribution) and `admixtools`/`admixr`
on PATH. Locally it auto-skips when the binaries are missing so day-to-day
dev doesn't bog down on infrastructure that isn't there.

Per LLD §6 nightly workflow + §5.5 markers.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from tests.fixtures.modifiers import pfile_to_eigenstrat
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile

pytestmark = [pytest.mark.external_tool, pytest.mark.slow]


def _which_or_skip(*names: str) -> dict[str, str]:
    """Return name→absolute-path for each binary, or pytest.skip if any
    is missing. Keeps the skip message specific so the nightly logs say
    exactly which tool is unavailable."""
    resolved = {n: shutil.which(n) for n in names}
    missing = [n for n, p in resolved.items() if p is None]
    if missing:
        pytest.skip(f"external tool(s) not on PATH: {', '.join(missing)}")
    return {n: p for n, p in resolved.items() if p is not None}


def _write_mergeit_par(par_path: Path, geno_a: Path, geno_b: Path, out_prefix: Path) -> None:
    """Minimal mergeit `parameter file` driving the merge. mergeit's
    parameter-file format is documented in EIGENSOFT/AT2 sources."""
    par_path.write_text(
        "\n".join(
            [
                f"geno1: {geno_a}.geno",
                f"snp1:  {geno_a}.snp",
                f"ind1:  {geno_a}.ind",
                f"geno2: {geno_b}.geno",
                f"snp2:  {geno_b}.snp",
                f"ind2:  {geno_b}.ind",
                f"genooutfilename: {out_prefix}.geno",
                f"snpoutfilename:  {out_prefix}.snp",
                f"indoutfilename:  {out_prefix}.ind",
                "outputformat: PACKEDANCESTRYMAP",
                "docheck: YES",
                "hashcheck: NO",
            ]
        )
        + "\n"
    )


def test_mergeit_f2_parity(tmp_path: Path) -> None:
    """End-to-end f2 parity. Skipped unless mergeit + AT2 extract_f2 are on
    PATH. The two converters MUST agree on every pairwise f2 to within
    1e-9 (AT2's documented numerical tolerance for identical genotype
    inputs)."""
    bins = _which_or_skip("plink2", "mergeit", "qpfstats")
    plink2 = bins["plink2"]
    mergeit = bins["mergeit"]
    qpfstats = bins["qpfstats"]

    # Two disjoint-sample panels sharing variants.
    spec_a = SyntheticPanelSpec(
        n_samples=40,
        n_variants=2_000,
        n_populations=4,
        seed=0xF21A,
        sample_id_prefix="A_",
    )
    spec_b = SyntheticPanelSpec(
        n_samples=40,
        n_variants=2_000,
        n_populations=4,
        seed=0xF21A,  # same variants
        sample_id_prefix="B_",
        sample_seed=0xF21B,
    )
    desc_a = synthesize_pfile(spec_a, tmp_path / "A")
    desc_b = synthesize_pfile(spec_b, tmp_path / "B")

    pfile_a = Path(str(desc_a.pgen_path)[:-5])
    pfile_b = Path(str(desc_b.pgen_path)[:-5])

    eig_a = pfile_to_eigenstrat(pfile_a, tmp_path / "eig_A")
    eig_b = pfile_to_eigenstrat(pfile_b, tmp_path / "eig_B")

    # Path 1: mergeit reference
    mergeit_out = tmp_path / "mergeit_out"
    par = tmp_path / "mergeit.par"
    _write_mergeit_par(par, eig_a, eig_b, mergeit_out)
    subprocess.run(
        [mergeit, "-p", str(par)],
        shell=False,
        check=True,
        capture_output=True,
        text=True,
    )

    # Path 2: pgen-samplebind merge over the same two EIGENSTRAT inputs
    samplebind_out = tmp_path / "samplebind_out"
    subprocess.run(
        [
            "pgen-samplebind",
            "merge",
            str(eig_a),
            str(eig_b),
            "-o",
            str(samplebind_out),
            "--quiet",
        ],
        shell=False,
        check=True,
        capture_output=True,
        text=True,
    )

    # Convert pgen-samplebind output back to EIGENSTRAT for an apples-to-apples
    # qpfstats input.
    samplebind_eig = pfile_to_eigenstrat(samplebind_out, tmp_path / "samplebind_eig")

    # f2 over both outputs; assert max_dev < 1e-9 across every pair.
    # qpfstats invocation: hand-written par per AT2 docs.
    def _run_qpfstats(eig_prefix: Path, label: str) -> Path:
        par_p = tmp_path / f"qpfstats_{label}.par"
        f2_out = tmp_path / f"f2_{label}.tsv"
        par_p.write_text(
            "\n".join(
                [
                    f"genotypename: {eig_prefix}.geno",
                    f"snpname:      {eig_prefix}.snp",
                    f"indivname:    {eig_prefix}.ind",
                    f"f2filename:   {f2_out}",
                ]
            )
            + "\n"
        )
        subprocess.run(
            [qpfstats, "-p", str(par_p)],
            shell=False,
            check=True,
            capture_output=True,
            text=True,
        )
        return f2_out

    f2_mergeit = _run_qpfstats(mergeit_out, "mergeit")
    f2_samplebind = _run_qpfstats(samplebind_eig, "samplebind")

    # Per-pair max_dev < 1e-9
    def _parse_f2(p: Path) -> dict[tuple[str, str], float]:
        out: dict[tuple[str, str], float] = {}
        for line in p.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 3:
                out[(parts[0], parts[1])] = float(parts[2])
        return out

    f2_a = _parse_f2(f2_mergeit)
    f2_b = _parse_f2(f2_samplebind)
    assert set(f2_a) == set(f2_b), (
        f"population pair sets differ: mergeit={set(f2_a)}, samplebind={set(f2_b)}"
    )
    for pair, val_a in f2_a.items():
        val_b = f2_b[pair]
        max_dev = abs(val_a - val_b)
        assert max_dev < 1e-9, (
            f"f2 mismatch for {pair}: mergeit={val_a}, samplebind={val_b}, max_dev={max_dev}"
        )

    # Touch plink2 just so the test doesn't claim a binary it never used.
    assert Path(plink2).exists()
