"""File-level mutate-and-write helpers for HLD integration tests.

Per LLD §5.2. Each helper reads a PFILE, mutates it according to the
operation, and writes a new PFILE. Materializing to disk lets failing
tests inspect the modified fixture post-mortem.

Operations:
- `flip_strand`: complement REF/ALT in .pvar at given variant indices;
  genotypes unchanged (strand-flip on biallelic SNPs is metadata-only).
- `swap_ref_alt`: swap REF/ALT in .pvar AND apply 2-x to genotypes at
  those variants (preserving -9). Used by HLD test 4 (allele-swap recovery).
- `drop_variants`: drop the specified variant indices; write the kept
  subset. Used by HLD test 6 (missing-variant default).
- `subset_variants`: keep only the specified indices. Inverse of drop.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _read_pvar_raw(prefix: Path) -> pd.DataFrame:
    """Read the .pvar without the biallelic-SNP filter (preserves all rows
    so the modifier helpers can operate by row index 1:1 with the .pgen)."""
    return pd.read_csv(Path(str(prefix) + ".pvar"), sep="\t")


def _write_pvar_raw(df: pd.DataFrame, prefix: Path) -> None:
    df.to_csv(Path(str(prefix) + ".pvar"), sep="\t", index=False, lineterminator="\n")


def _read_full_pgen(prefix: Path, n_samples: int, n_variants: int) -> np.ndarray:
    import pgenlib

    reader = pgenlib.PgenReader(str(Path(str(prefix) + ".pgen")).encode(), raw_sample_ct=n_samples)
    try:
        buf = np.empty((n_variants, n_samples), dtype=np.int8)
        reader.read_range(0, n_variants, buf, sample_maj=0)
        return buf
    finally:
        if hasattr(reader, "close"):
            reader.close()


def _write_full_pgen(geno: np.ndarray, prefix: Path) -> None:
    import pgenlib

    n_variants, n_samples = geno.shape
    writer = pgenlib.PgenWriter(str(Path(str(prefix) + ".pgen")).encode(), n_samples, n_variants)
    try:
        writer.append_biallelic_batch(geno)
    finally:
        writer.close()


def flip_strand(
    in_prefix: Path,
    variant_indices: np.ndarray,
    out_prefix: Path,
) -> None:
    """Complement REF and ALT at the given variant indices in the .pvar.

    Genotypes unchanged (per LLD §2.1: strand-flip on biallelic SNPs is
    metadata-only — the hardcall encoding stays the same once you accept
    that REF/ALT now refer to the complemented bases).
    """
    pvar_df = _read_pvar_raw(in_prefix)
    for idx in variant_indices:
        pvar_df.iat[int(idx), pvar_df.columns.get_loc("REF")] = _COMPLEMENT[
            pvar_df.iat[int(idx), pvar_df.columns.get_loc("REF")]
        ]
        pvar_df.iat[int(idx), pvar_df.columns.get_loc("ALT")] = _COMPLEMENT[
            pvar_df.iat[int(idx), pvar_df.columns.get_loc("ALT")]
        ]
    _write_pvar_raw(pvar_df, out_prefix)
    shutil.copy(Path(str(in_prefix) + ".pgen"), Path(str(out_prefix) + ".pgen"))
    shutil.copy(Path(str(in_prefix) + ".psam"), Path(str(out_prefix) + ".psam"))


def swap_ref_alt(
    in_prefix: Path,
    n_samples: int,
    n_variants: int,
    variant_indices: np.ndarray,
    out_prefix: Path,
) -> None:
    """Swap REF/ALT in .pvar AND apply x → 2-x to the genotype rows at the
    specified variants (preserving -9 for missing).

    Used by HLD test 4 (allele-swap recovery): merge should auto-recode the
    swapped genotypes back to the canonical encoding via REF_ALT_SWAP.
    """
    # .pvar: swap REF and ALT.
    pvar_df = _read_pvar_raw(in_prefix)
    ref_col = pvar_df.columns.get_loc("REF")
    alt_col = pvar_df.columns.get_loc("ALT")
    for idx in variant_indices:
        ridx = int(idx)
        old_ref = pvar_df.iat[ridx, ref_col]
        old_alt = pvar_df.iat[ridx, alt_col]
        pvar_df.iat[ridx, ref_col] = old_alt
        pvar_df.iat[ridx, alt_col] = old_ref
    _write_pvar_raw(pvar_df, out_prefix)

    # .pgen: 2-x at the specified rows, preserving -9.
    geno = _read_full_pgen(in_prefix, n_samples, n_variants)
    for idx in variant_indices:
        ridx = int(idx)
        not_missing = geno[ridx] != -9
        geno[ridx, not_missing] = 2 - geno[ridx, not_missing]
    _write_full_pgen(geno, out_prefix)

    shutil.copy(Path(str(in_prefix) + ".psam"), Path(str(out_prefix) + ".psam"))


def drop_variants(
    in_prefix: Path,
    n_samples: int,
    n_variants: int,
    drop_indices: np.ndarray,
    out_prefix: Path,
) -> None:
    """Drop the specified variant indices; write a new PFILE with the kept
    subset. Used by HLD test 6 (missing-variant default)."""
    keep_mask = np.ones(n_variants, dtype=bool)
    keep_mask[drop_indices] = False

    # Filter .pvar
    pvar_df = _read_pvar_raw(in_prefix)
    pvar_df = pvar_df.loc[keep_mask].reset_index(drop=True)
    _write_pvar_raw(pvar_df, out_prefix)

    # Filter .pgen
    geno = _read_full_pgen(in_prefix, n_samples, n_variants)
    _write_full_pgen(geno[keep_mask], out_prefix)

    shutil.copy(Path(str(in_prefix) + ".psam"), Path(str(out_prefix) + ".psam"))


def subset_variants(
    in_prefix: Path,
    n_samples: int,
    n_variants: int,
    keep_indices: np.ndarray,
    out_prefix: Path,
) -> None:
    """Keep only the specified variant indices; inverse of drop_variants."""
    keep_mask = np.zeros(n_variants, dtype=bool)
    keep_mask[keep_indices] = True

    pvar_df = _read_pvar_raw(in_prefix)
    pvar_df = pvar_df.loc[keep_mask].reset_index(drop=True)
    _write_pvar_raw(pvar_df, out_prefix)

    geno = _read_full_pgen(in_prefix, n_samples, n_variants)
    _write_full_pgen(geno[keep_mask], out_prefix)

    shutil.copy(Path(str(in_prefix) + ".psam"), Path(str(out_prefix) + ".psam"))


def find_ambiguous_variant_indices(in_prefix: Path) -> np.ndarray:
    """Return indices of A/T and C/G variants in the .pvar (1240k design has
    ~5-8% naturally; the synthesizer matches this fraction)."""
    pvar_df = _read_pvar_raw(in_prefix)
    is_at = ((pvar_df["REF"] == "A") & (pvar_df["ALT"] == "T")) | (
        (pvar_df["REF"] == "T") & (pvar_df["ALT"] == "A")
    )
    is_cg = ((pvar_df["REF"] == "C") & (pvar_df["ALT"] == "G")) | (
        (pvar_df["REF"] == "G") & (pvar_df["ALT"] == "C")
    )
    return np.where(is_at | is_cg)[0]


def find_unambiguous_variant_indices(in_prefix: Path) -> np.ndarray:
    """Return indices of variants whose REF/ALT pair is NOT A/T or C/G."""
    pvar_df = _read_pvar_raw(in_prefix)
    is_at = ((pvar_df["REF"] == "A") & (pvar_df["ALT"] == "T")) | (
        (pvar_df["REF"] == "T") & (pvar_df["ALT"] == "A")
    )
    is_cg = ((pvar_df["REF"] == "C") & (pvar_df["ALT"] == "G")) | (
        (pvar_df["REF"] == "G") & (pvar_df["ALT"] == "C")
    )
    return np.where(~(is_at | is_cg))[0]


def pfile_to_eigenstrat(pfile_prefix: Path, out_prefix: Path) -> Path:
    """Convert a PFILE to EIGENSTRAT (.geno/.snp/.ind) via plink2 shell-out.

    Used by HLD tests 7, 8, 10 to bootstrap valid PACKEDANCESTRYMAP
    EIGENSTRAT fixtures (plink2's --eigfile requires the binary header'd
    format; the simple ASCII per-line format isn't accepted).

    Returns out_prefix on success; raises subprocess.CalledProcessError
    on plink2 failure.
    """
    import shutil
    import subprocess

    plink2 = shutil.which("plink2")
    if plink2 is None:
        raise RuntimeError(
            "plink2 not on PATH; required for pfile_to_eigenstrat. "
            "Install plink2 v2.x or skip the EIGENSTRAT test (mark "
            "with @pytest.mark.eigenstrat)."
        )
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            plink2,
            "--pfile",
            str(pfile_prefix),
            "--export",
            "eig",
            "--out",
            str(out_prefix),
        ],
        shell=False,
        check=True,
        capture_output=True,
        text=True,
    )
    return out_prefix


def make_panel_with_iids(
    out_prefix: Path,
    iids: list[str],
    n_variants: int = 10,
    seed: int = 42,
    pop: str = "pop_a",
) -> Path:
    """Build a tiny PFILE with explicitly-specified IIDs.

    Used by HLD test 23 to construct collision scenarios: e.g., a panel
    that already contains both `Sample1` and `Sample1_1` (the iterative-
    build case), or two panels that share specific IIDs by construction.
    """
    import pgenlib

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    n_samples = len(iids)
    rng = np.random.default_rng(seed)
    geno = rng.integers(0, 3, size=(n_variants, n_samples), dtype=np.int8)

    out_pgen = Path(str(out_prefix) + ".pgen")
    out_pvar = Path(str(out_prefix) + ".pvar")
    out_psam = Path(str(out_prefix) + ".psam")

    writer = pgenlib.PgenWriter(str(out_pgen).encode(), n_samples, n_variants)
    try:
        writer.append_biallelic_batch(geno)
    finally:
        writer.close()

    pvar_df = pd.DataFrame(
        {
            "#CHROM": [1] * n_variants,
            "POS": list(range(1000, 1000 + n_variants)),
            "ID": [f"v{i}" for i in range(n_variants)],
            "REF": ["A"] * n_variants,
            "ALT": ["G"] * n_variants,
        }
    )
    pvar_df.to_csv(out_pvar, sep="\t", index=False, lineterminator="\n")

    psam_df = pd.DataFrame(
        {
            "#IID": iids,
            "SEX": [1] * n_samples,
            "POP": [pop] * n_samples,
        }
    )
    psam_df.to_csv(out_psam, sep="\t", index=False, lineterminator="\n")
    return out_prefix


def make_corrupt_eigenstrat(out_prefix: Path) -> Path:
    """Write a deliberately-malformed EIGENSTRAT triplet for HLD test 19.

    The .geno lacks the required `GENO ` header so plink2 --eigfile rejects
    it with a clear error — that's the failure path we want to exercise
    (subprocess fails → IOFailure with stderr surfaced → tempdir cleanup).
    """
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(out_prefix) + ".geno").write_bytes(b"this is not a valid GENO header\n")
    Path(str(out_prefix) + ".snp").write_text("rs1\t1\t0.0\t100\tA\tG\nrs2\t1\t0.0\t200\tC\tT\n")
    Path(str(out_prefix) + ".ind").write_text("S00000\tM\tpop_a\nS00001\tF\tpop_b\n")
    return out_prefix
