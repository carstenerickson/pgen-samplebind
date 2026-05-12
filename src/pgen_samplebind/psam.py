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


def read_relabel_tsv(
    path: Path,
    input_col: str | None,
    output_col: str | None,
) -> pd.DataFrame:
    """Read a relabel TSV in either 2-col or N-col form per HLD §Relabeling.

    Decision rule:
      - Both input_col and output_col absent → 2-col header-less TSV
        (rows are `<from_pop>\\t<to_pop>` pairs; source column is the
        sample's POP).
      - Both present → N-col header-required TSV (e.g., AADR .anno);
        the named columns are pulled out and renamed to "input" and
        "output". Source column is whatever id_column the orchestrator
        passes (default IID).
      - One present, the other absent → UsageError (ambiguous).

    Returns a DataFrame with two columns: "input" and "output".

    Raises:
        UsageError: only one of input_col/output_col supplied; or N-col
            named columns not present in the file.
        IOFailure: TSV unreadable or unparseable.
    """
    if (input_col is None) != (output_col is None):
        raise UsageError(
            "--relabel-input-col and --relabel-output-col must be supplied "
            "together (N-col form), or both omitted (2-col header-less form). "
            f"Got input_col={input_col!r}, output_col={output_col!r}."
        )

    try:
        if input_col is None:
            # 2-col header-less TSV
            df = pd.read_csv(
                path,
                sep="\t",
                header=None,
                names=["input", "output"],
                dtype=str,
                na_filter=False,
            )
            if df.shape[1] != 2:
                raise UsageError(
                    f"--relabel-from {path}: expected 2 tab-delimited columns "
                    f"(no header) but found {df.shape[1]}. For wider TSVs use "
                    f"--relabel-input-col and --relabel-output-col to pick "
                    f"two columns by name."
                )
            return df

        # N-col header-required TSV
        df = pd.read_csv(path, sep="\t", dtype=str, na_filter=False)
    except (OSError, pd.errors.ParserError) as e:
        raise IOFailure(f"cannot parse --relabel-from {path}: {e}") from e

    if input_col not in df.columns:
        raise UsageError(
            f"--relabel-input-col {input_col!r} not in {path}. "
            f"Available columns: {list(df.columns)}"
        )
    if output_col not in df.columns:
        raise UsageError(
            f"--relabel-output-col {output_col!r} not in {path}. "
            f"Available columns: {list(df.columns)}"
        )
    return df[[input_col, output_col]].rename(columns={input_col: "input", output_col: "output"})


def apply_relabel(
    psam: pd.DataFrame,
    relabel: pd.DataFrame,
    source_column: str,
    target_column: str = "POP",
) -> pd.DataFrame:
    """Map each sample's source_column value through the relabel; if found,
    set target_column to the relabel's output value. Otherwise leave
    target_column as-is. Per HLD §Relabeling.

    For 2-col relabel: source_column = "POP" (rows are POP→POP, collapse
    populations across inputs).
    For N-col relabel: source_column = id_column (default IID; rows are
    sample-id → POP, override per-sample).

    Returns a NEW DataFrame; input untouched.
    """
    if source_column not in psam.columns:
        raise UsageError(
            f"--relabel-from source column {source_column!r} not present in "
            f".psam columns: {list(psam.columns)}. Use --id-column NAME if your "
            f".psam uses a non-IID identifier."
        )
    if target_column not in psam.columns:
        raise UsageError(
            f"--relabel-from target column {target_column!r} not present in "
            f".psam columns: {list(psam.columns)} (expected POP after auto-detect)."
        )

    out = psam.copy()
    mapping: dict[str, str] = dict(
        zip(relabel["input"].astype(str), relabel["output"].astype(str), strict=True)
    )
    out[target_column] = [
        mapping.get(str(s), t) for s, t in zip(out[source_column], out[target_column], strict=True)
    ]
    return out


def resolve_sample_identity(
    psams: list[pd.DataFrame],
    policy: MergePolicy,
    target_idxs: tuple[int, ...] = (),
) -> SampleIdentityPlan:
    """Compute the SampleIdentityPlan from input psams + --on-collision policy.
    Per HLD §IID collision handling (v3.5) and LLD §3.5.

    `target_idxs` is the (possibly empty) tuple of input indexes marked as
    targets via --target. Zero or one target: target collision suffix is
    bare `_target`. Two or more targets: suffix is `_target_<input_idx>`
    so each target's renamed sample can be traced back unambiguously.

    Supports all three policies:
      - error: raises InvariantViolation on first collision (exit 3).
      - first: drops duplicates in input order; first occurrence wins.
      - suffix: renames duplicates per HLD v3.5:
          * General mode: input[N>0]'s colliding sample gets `_<input_idx>`.
          * Single-target mode (one entry in target_idxs): suffix is `_target`.
          * Multi-target mode (>=2 in target_idxs): suffix is `_target_<input_idx>`.
          * Idempotent retry: if the renamed slot is also taken, fall through
            to `<base>_<suffix>_1`, `<base>_<suffix>_2`, ... until free.
        Input[0]'s IIDs are never suffixed (canonical, preserved).

    Raises:
        InvariantViolation: collision under --on-collision error;
            input[0] contains its own internal duplicate (canonical must
            have unique IIDs even under --on-collision suffix).
    """
    if policy.on_collision == "suffix":
        return _resolve_with_suffix(psams, target_idxs)

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


def _resolve_with_suffix(
    psams: list[pd.DataFrame],
    target_idxs: tuple[int, ...],
) -> SampleIdentityPlan:
    """--on-collision suffix scheme per HLD §IID collision handling (v3.5).

    Algorithm: sequential pass over inputs. Input[0] is canonical and never
    suffixed (raises if input[0] has internal duplicates). For input[N>0],
    each colliding IID gets one of:
      - `_target`            — when input_idx is the single target
      - `_target_<input_idx>` — when input_idx is one of multiple targets
      - `_<input_idx>`        — non-target
    If the suffixed name is also taken, fall through to `<base>_<suffix>_1`,
    `<base>_<suffix>_2`, ... until a free slot is found.
    """
    target_idx_set = frozenset(target_idxs)
    single_target = len(target_idxs) == 1
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
            if input_idx in target_idx_set:
                base_suffix = "_target" if single_target else f"_target_{input_idx}"
            else:
                base_suffix = f"_{input_idx}"
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
