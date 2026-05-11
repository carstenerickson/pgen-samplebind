""".psam parsing, population-column auto-detect, sample identity resolution,
column union, write.

Per LLD §3.5. Day 1: read paths. Day 3: resolve_sample_identity, merge_psams,
write_psam. --on-collision suffix and --relabel-from deferred to later days.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .errors import InvariantViolation, IOFailure, UsageError
from .types import MergePolicy, SampleIdentityPlan

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


def resolve_sample_identity(
    psams: list[pd.DataFrame],
    policy: MergePolicy,
    target_idx: int | None,
) -> SampleIdentityPlan:
    """Compute the SampleIdentityPlan from input psams + --on-collision policy.
    Per HLD §IID collision handling (v3.5) and LLD §3.5.

    Supports all three policies:
      - error: raises InvariantViolation on first collision (exit 3).
      - first: drops duplicates in input order; first occurrence wins.
      - suffix: renames duplicates per HLD v3.5:
          * General mode: input[N>0]'s colliding sample gets `_<input_idx>`.
          * Target mode (input_idx == target_idx): suffix is `_target`.
          * Idempotent retry: if the renamed slot is also taken, fall through
            to `<base>_<suffix>_1`, `<base>_<suffix>_2`, ... until free.
        Input[0]'s IIDs are never suffixed (canonical, preserved).

    Raises:
        InvariantViolation: collision under --on-collision error;
            input[0] contains its own internal duplicate (canonical must
            have unique IIDs even under --on-collision suffix).
    """
    if policy.on_collision == "suffix":
        return _resolve_with_suffix(psams, target_idx)

    output_iids: list[str] = []
    keep_masks: list[np.ndarray[Any, Any]] = []
    output_indices_per_input: list[np.ndarray[Any, Any]] = []
    seen: set[str] = set()

    for input_idx, psam in enumerate(psams):
        iids = psam["IID"].tolist()
        n = len(iids)
        keep = np.ones(n, dtype=bool)
        out_idx = np.full(n, -1, dtype=np.int64)

        for i, iid in enumerate(iids):
            if iid in seen:
                if policy.on_collision == "error":
                    raise InvariantViolation(
                        f"--on-collision error: IID {iid!r} appears in input[{input_idx}] "
                        f"and an earlier input. Use --on-collision first (drop duplicates) "
                        f"or suffix (rename duplicates) to allow merge to proceed."
                    )
                # first: drop this duplicate (mask out)
                keep[i] = False
            else:
                seen.add(iid)
                out_idx[i] = len(output_iids)
                output_iids.append(iid)

        keep_masks.append(keep)
        output_indices_per_input.append(out_idx)

    return SampleIdentityPlan(
        output_iids=tuple(output_iids),
        per_input_keep_mask=tuple(keep_masks),
        per_input_output_indices=tuple(output_indices_per_input),
    )


def _resolve_with_suffix(psams: list[pd.DataFrame], target_idx: int | None) -> SampleIdentityPlan:
    """--on-collision suffix scheme per HLD §IID collision handling (v3.5).

    Algorithm: sequential pass over inputs. Input[0] is canonical and never
    suffixed (raises if input[0] has internal duplicates). For input[N>0],
    each colliding IID gets `_<input_idx>` (general) or `_target` (when
    input_idx == target_idx) appended; if the suffixed name is also taken,
    fall through to `<base>_<suffix>_1`, `<base>_<suffix>_2`, ... until a
    free slot is found.
    """
    output_iids: list[str] = []
    keep_masks: list[np.ndarray[Any, Any]] = []
    output_indices_per_input: list[np.ndarray[Any, Any]] = []
    seen: set[str] = set()

    for input_idx, psam in enumerate(psams):
        iids = psam["IID"].tolist()
        n = len(iids)
        keep = np.ones(n, dtype=bool)
        out_idx = np.full(n, -1, dtype=np.int64)

        for i, iid in enumerate(iids):
            if input_idx == 0:
                # Canonical input — never suffixed. Internal duplicates here
                # are a bug in the input itself, not something --on-collision
                # suffix can resolve (which non-canonical sample would the
                # canonical's duplicate "rename to"?).
                if iid in seen:
                    raise InvariantViolation(
                        f"input[0] contains duplicate IID {iid!r}; canonical input "
                        f"must have unique IIDs. --on-collision suffix only renames "
                        f"non-canonical duplicates."
                    )
                seen.add(iid)
                out_idx[i] = len(output_iids)
                output_iids.append(iid)
                continue

            if iid not in seen:
                seen.add(iid)
                out_idx[i] = len(output_iids)
                output_iids.append(iid)
                continue

            # Collision — apply suffix scheme with idempotent retry.
            base_suffix = "_target" if input_idx == target_idx else f"_{input_idx}"
            candidate = iid + base_suffix
            retry = 1
            while candidate in seen:
                candidate = iid + base_suffix + f"_{retry}"
                retry += 1
            seen.add(candidate)
            out_idx[i] = len(output_iids)
            output_iids.append(candidate)

        keep_masks.append(keep)
        output_indices_per_input.append(out_idx)

    return SampleIdentityPlan(
        output_iids=tuple(output_iids),
        per_input_keep_mask=tuple(keep_masks),
        per_input_output_indices=tuple(output_indices_per_input),
    )


def merge_psams(
    psams: list[pd.DataFrame],
    sample_plan: SampleIdentityPlan,
) -> pd.DataFrame:
    """Concatenate psams using the resolved sample identity, with column union.

    Day 3: simple concat with NA fill across input column unions. Conflict
    detection (same column name, conflicting values for shared samples) is
    deferred to Day 4 reporting work.

    Returns a DataFrame with `len(sample_plan.output_iids)` rows in the
    output_iids order.
    """
    kept = []
    for psam, keep_mask in zip(psams, sample_plan.per_input_keep_mask, strict=True):
        kept.append(psam.iloc[keep_mask].reset_index(drop=True))
    merged = pd.concat(kept, ignore_index=True, sort=False)

    # Replace IIDs with the (possibly renamed-via-suffix-future) output_iids,
    # in case sample_plan introduces suffixes. Day 3 uses identity for first/
    # error policies, so this is a no-op here.
    merged["IID"] = list(sample_plan.output_iids)
    return merged


def write_psam(df: pd.DataFrame, path: Path) -> None:
    """Write canonical .psam: column order FID IID SEX POP PSEUDOHAPLOID then
    extras in alphabetical order. Tab-delimited, header line prefixed with '#'.

    Per plink2 .psam spec: when FID is present it MUST be the first column
    (the header line then starts with `#FID`). When FID is absent, IID is
    first (`#IID`). The merge orchestrator always emits FID (computed as
    FID = POP via add_fid_from_pop) so the produced .psam is the FID-first
    form.

    Raises:
        IOFailure: write failed.
    """
    canonical = ["FID", "IID", "SEX", "POP", "PSEUDOHAPLOID"]
    extras = sorted(c for c in df.columns if c not in canonical)
    present_canonical = [c for c in canonical if c in df.columns]
    column_order = present_canonical + extras

    out = df[column_order].copy()
    # Plink2 convention: header line starts with `#` on the first column.
    out.columns = ["#" + c if i == 0 else c for i, c in enumerate(out.columns)]
    try:
        out.to_csv(path, sep="\t", index=False, lineterminator="\n", na_rep="NA")
    except OSError as e:
        raise IOFailure(f"cannot write {path}: {e}") from e
