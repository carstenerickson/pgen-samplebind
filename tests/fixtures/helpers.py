"""Shared test helpers — readers for the PFILE triplet a merge writes.

Previously these were copy-pasted across 7 integration test files. All
prefix-based: `prefix.pgen` / `prefix.pvar` / `prefix.psam` resolved
internally so callers pass the same `Path` everywhere.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_pgen_full(prefix: Path, n_samples: int, n_variants: int) -> np.ndarray:
    """Read the entire `prefix.pgen` into an (n_variants, n_samples) int8 matrix."""
    import pgenlib

    pgen_path = Path(str(prefix) + ".pgen")
    reader = pgenlib.PgenReader(str(pgen_path).encode(), raw_sample_ct=n_samples)
    try:
        buf = np.empty((n_variants, n_samples), dtype=np.int8)
        reader.read_range(0, n_variants, buf, sample_maj=0)
        return buf
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()


def read_psam(prefix: Path) -> pd.DataFrame:
    """Read `prefix.psam` and strip the leading `#` from the first header cell."""
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_psam_iids(prefix: Path) -> list[str]:
    """Return the IID column of `prefix.psam` as a list."""
    df = read_psam(prefix)
    return df["IID"].tolist()


def read_pvar(prefix: Path) -> pd.DataFrame:
    """Read `prefix.pvar` and strip the leading `#` from the first header cell."""
    df = pd.read_csv(Path(str(prefix) + ".pvar"), sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_pvar_keys(prefix: Path) -> list[tuple[int, int]]:
    """Return the (chrom, pos) list from `prefix.pvar`."""
    df = read_pvar(prefix)
    return list(zip(df["CHROM"].astype(int), df["POS"].astype(int), strict=True))
