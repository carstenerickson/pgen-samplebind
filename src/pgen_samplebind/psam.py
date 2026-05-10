""".psam parsing, population-column auto-detect, basic column ops.

Per LLD §3.5. Day 1 implements only the read paths needed by `inspect`;
collision resolution, column union, and relabel are deferred to later days.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .errors import InvariantViolation, IOFailure, UsageError

_POP_COLUMN_FALLBACKS = ("POP", "PHENO", "PHENO1")


def read_psam(path: Path) -> pd.DataFrame:
    """Parse a .psam (or .fam-with-header) file.

    Plink2 .psam files use tab- or whitespace-delimited columns. The header
    line is prefixed with `#` (e.g., `#FID` or `#IID`). Column dtypes
    preserved as input.

    Raises:
        IOFailure: file unreadable.
    """
    try:
        # Read all lines; find the header line (starts with `#`).
        # Plink2 .psam may have ##-prefixed metadata lines (rare).
        with open(path, encoding="utf-8") as f:
            header_idx = -1
            for i, line in enumerate(f):
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header_idx = i
                    break
                # First non-comment line — no header present
                break
        if header_idx == -1:
            # No `#`-prefixed header. Plink2 spec allows this for legacy .fam-style files;
            # for now we require the explicit header form per HLD canonical.
            raise IOFailure(
                f"{path}: no `#`-prefixed header line found. pgen-samplebind requires "
                f"plink2 .psam with explicit header (e.g., '#IID\\tSEX\\tPOP')"
            )
        df = pd.read_csv(
            path,
            sep=r"\s+",  # whitespace tolerant; plink2 outputs tabs but legacy is space
            skiprows=header_idx,
            header=0,
            engine="python",  # \s+ requires python engine
            dtype=str,
            na_filter=False,
        )
    except (OSError, pd.errors.ParserError) as e:
        raise IOFailure(f"cannot parse {path}: {e}") from e

    # Strip the leading `#` from the first column (e.g., '#IID' → 'IID')
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def detect_population_column(df: pd.DataFrame, override: str | None) -> str:
    """Return the column holding population labels per HLD §Population labels.

    Try POP, then PHENO, then PHENO1; or use override if given.

    Raises:
        UsageError: override specified but column not present.
        InvariantViolation: no population column found and no override.
    """
    if override is not None:
        if override not in df.columns:
            raise UsageError(
                f"--population-column {override!r} not present in .psam columns: {list(df.columns)}"
            )
        return override
    for candidate in _POP_COLUMN_FALLBACKS:
        if candidate in df.columns:
            return candidate
    raise InvariantViolation(
        f".psam has no population column (tried {_POP_COLUMN_FALLBACKS}); "
        f"available columns: {list(df.columns)}. Use --population-column NAME to override."
    )


def rename_to_pop(df: pd.DataFrame, source_col: str) -> pd.DataFrame:
    """Rename `source_col` to 'POP'. No-op if already 'POP'."""
    if source_col == "POP":
        return df
    return df.rename(columns={source_col: "POP"})


def add_fid_from_pop(df: pd.DataFrame) -> pd.DataFrame:
    """Add or overwrite FID column = POP per HLD §Output PFILE.

    AT2 extract_f2 keys on FID, so emitting it ensures downstream interop.
    """
    out = df.copy()
    out["FID"] = out["POP"]
    return out
