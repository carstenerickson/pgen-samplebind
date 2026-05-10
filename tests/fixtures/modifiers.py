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
