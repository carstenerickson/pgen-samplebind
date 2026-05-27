"""Pre-flight input compatibility report.

Always-emitted structured summary of how compatible the canonical input
is with each other (non-canonical) input. Written to
`<output_prefix>.preflight.json` before pass 2 of `merge`, so users (and
downstream pipelines) see the compatibility picture even on runs that
later fail mid-pass or that produce a misleadingly small merge.

This module owns the JSON schema v1 contract. Downstream pipelines may
assert against this file — e.g.,
`jq '.comparisons[0].intersection_fraction_of_min'`.

Schema-evolution policy: `classification` / `classification_evidence`
on each comparison and the top-level `gate` object are emitted with
placeholder values in v1 and will be populated by the failure-mode
classifier (issue [#12](https://github.com/carstenerickson/pgen-samplebind/issues/12),
step 3) and the policy gate (step 4) without bumping schema_version.
Reserve a schema_version bump for renames, removals, or
incompatible type changes.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import pandas as pd

from .errors import IOFailure, InvariantViolation
from .pvar import read_pvar
from .types import InputDescriptor, MergePolicy

PREFLIGHT_SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class PerChromCompat:
    """Per-chromosome variant-key counts and intersection for one (canonical, other) pair."""

    chrom: int
    canonical_size: int
    other_size: int
    intersection: int


@dataclass(frozen=True, slots=True)
class PairCompatibility:
    """How compatible one non-canonical input is with the canonical input.

    `alternate_key_*` fields carry intersection numbers under the
    non-active variant key (id when policy.variant_key=chr_pos and
    vice versa). The classifier compares the two to detect the case
    where the user picked the wrong --variant-key for their data:
    high alternate-key overlap with low active-key overlap is the
    `key_space_mismatch` signature.

    `classification` and `classification_evidence` are populated by
    `classify_pair`; downstream pipelines may key on them
    (`jq '.comparisons[0].classification'`).
    """

    input_index: int
    path: Path
    n_variants: int
    intersection: int
    intersection_fraction_of_min: float
    per_chrom: tuple[PerChromCompat, ...]
    alternate_key: str | None = None
    alternate_key_canonical_size: int | None = None
    alternate_key_other_size: int | None = None
    alternate_key_intersection: int | None = None
    alternate_key_fraction_of_min: float | None = None
    classification: str | None = None
    classification_evidence: dict[str, Any] | None = None


@dataclass(frozen=True, slots=True)
class PreflightReport:
    """Schema v1 envelope. The JSON file is a verbatim serialization of this shape."""

    schema_version: int
    tool: str
    tool_version: str
    command: str  # "merge" | "validate"
    variant_key: str  # mirrors policy.variant_key
    canonical: dict[str, Any]
    comparisons: tuple[PairCompatibility, ...]
    gate: dict[str, Any]


def compute_preflight(
    descriptors: list[InputDescriptor],
    policy: MergePolicy,
    *,
    tool_version: str,
    command: str,
    pvars: list[pd.DataFrame] | None = None,
) -> PreflightReport:
    """Build the preflight report by computing key-space intersection of
    each non-canonical input against the canonical (index 0).

    If `pvars` is supplied, it must be parallel to `descriptors` (same
    length, same order) and the reader skips the .pvar read entirely.
    `validate` passes the pvars it already read at its step 7 to avoid
    a double-read at 84M-variant scale; `merge` calls with `pvars=None`
    today (merge_inputs reads pvars internally for its own pass and we
    haven't yet plumbed that cache out).
    """
    if not descriptors:
        raise InvariantViolation("compute_preflight: no descriptors")
    if pvars is not None and len(pvars) != len(descriptors):
        raise InvariantViolation(
            f"compute_preflight: pvars length {len(pvars)} != descriptors length "
            f"{len(descriptors)} (callers must pass a list parallel to descriptors)"
        )

    canonical_desc = descriptors[0]
    canonical_df = pvars[0] if pvars is not None else read_pvar(canonical_desc.pvar_path)
    canonical_keys = _keys_for(canonical_df, policy.variant_key)
    canonical_by_chrom = _keys_by_chrom(canonical_df, policy.variant_key)
    alternate_key = _alternate_key(policy.variant_key)
    canonical_alt_keys = _keys_for(canonical_df, alternate_key)

    comparisons: list[PairCompatibility] = []
    for i, desc in enumerate(descriptors[1:], start=1):
        other_df = pvars[i] if pvars is not None else read_pvar(desc.pvar_path)
        other_keys = _keys_for(other_df, policy.variant_key)
        other_by_chrom = _keys_by_chrom(other_df, policy.variant_key)

        intersection = len(canonical_keys & other_keys)
        denom = min(len(canonical_keys), len(other_keys))
        fraction = (intersection / denom) if denom > 0 else 0.0

        per_chrom = tuple(
            PerChromCompat(
                chrom=chrom,
                canonical_size=len(canonical_by_chrom.get(chrom, set())),
                other_size=len(other_by_chrom.get(chrom, set())),
                intersection=len(
                    canonical_by_chrom.get(chrom, set()) & other_by_chrom.get(chrom, set())
                ),
            )
            for chrom in sorted(set(canonical_by_chrom) | set(other_by_chrom))
        )

        # Alternate-key view, for key_space_mismatch detection.
        other_alt_keys = _keys_for(other_df, alternate_key)
        alt_intersection = len(canonical_alt_keys & other_alt_keys)
        alt_denom = min(len(canonical_alt_keys), len(other_alt_keys))
        alt_fraction = (alt_intersection / alt_denom) if alt_denom > 0 else 0.0

        pair = PairCompatibility(
            input_index=i,
            path=desc.path,
            n_variants=len(other_keys),
            intersection=intersection,
            intersection_fraction_of_min=fraction,
            per_chrom=per_chrom,
            alternate_key=alternate_key,
            alternate_key_canonical_size=len(canonical_alt_keys),
            alternate_key_other_size=len(other_alt_keys),
            alternate_key_intersection=alt_intersection,
            alternate_key_fraction_of_min=alt_fraction,
        )
        label, evidence = classify_pair(pair)
        pair = replace(pair, classification=label, classification_evidence=evidence)
        comparisons.append(pair)

    return PreflightReport(
        schema_version=PREFLIGHT_SCHEMA_VERSION,
        tool="pgen-samplebind",
        tool_version=tool_version,
        command=command,
        variant_key=policy.variant_key,
        canonical={
            "path": str(canonical_desc.path),
            "index": 0,
            "n_variants": len(canonical_keys),
        },
        comparisons=tuple(comparisons),
        gate={"triggered": False, "action": "none", "threshold": None},
    )


# Labels that the gate treats as actionable failures. `compatible` is
# the happy path; `empty_input` is its own diagnostic the merge will
# surface anyway via downstream validation. Build/key-space/disjoint
# all silently produce near-empty merges today — they are the failure
# shapes the gate exists to catch.
GATE_FAILURE_LABELS: frozenset[str] = frozenset(
    {"build_mismatch", "key_space_mismatch", "disjoint_panels"}
)


def evaluate_gate(report: "PreflightReport", policy: MergePolicy) -> "PreflightReport":
    """Decide what action the policy mandates given the report's classifications.

    Returns a new `PreflightReport` with `gate` populated:
      - `triggered`: whether the policy DID take action (false under `off`
        even if classifications failed — `would_trigger` carries the
        classification-level signal independently).
      - `would_trigger`: any comparison classified into GATE_FAILURE_LABELS
        regardless of policy. Workflow consumers gating CI on the raw
        classification signal should key on this rather than `triggered`,
        so a user's explicit `--preflight-policy off` doesn't surface as
        a CI failure.
      - `action`: "none" | "warn" | "error" — what the caller should do.
      - `policy`: echo of policy.preflight_policy.
      - `threshold`: the classifier's compatibility threshold (echoed for
        consumers — the classifier itself owns the cutoff).
      - `failing_inputs`: list of `{input_index, classification}` for the
        offending comparisons. Populated whenever `would_trigger=true`,
        regardless of action, so diagnostics survive policy=off.

    Pure function — no I/O, no stderr — so callers can inspect the
    decision before deciding whether to raise / warn / continue. The CLI
    wrapper does the actual stderr emission and exception raising.

    Raises `InvariantViolation` if any comparison has `classification=None`
    — callers must run `classify_pair` (or use `compute_preflight`, which
    runs it automatically) before invoking the gate. The frozen-dataclass
    default allows None so unit tests can construct partial reports; the
    gate refuses to silently treat such cases as non-failing.
    """
    unclassified = [c.input_index for c in report.comparisons if c.classification is None]
    if unclassified:
        raise InvariantViolation(
            f"evaluate_gate called on un-classified comparisons "
            f"(input_index={unclassified}); run classify_pair first or call "
            f"compute_preflight, which classifies inline."
        )

    failing = [
        {"input_index": c.input_index, "classification": c.classification}
        for c in report.comparisons
        if c.classification in GATE_FAILURE_LABELS
    ]
    would_trigger = bool(failing)

    if not would_trigger:
        action = "none"
    elif policy.preflight_policy == "off":
        action = "none"
    elif policy.preflight_policy == "strict":
        action = "error"
    else:  # "warn"
        action = "warn"

    # `triggered` reflects whether the policy ACTED on the failure.
    # `would_trigger` is the classification-level signal, surfaced
    # separately so JSON consumers can tell apart "no failure" from
    # "failure suppressed by --preflight-policy off".
    new_gate = {
        "triggered": action != "none",
        "would_trigger": would_trigger,
        "action": action,
        "policy": policy.preflight_policy,
        "threshold": _COMPATIBLE_MIN_FRACTION,
        "failing_inputs": failing,
    }
    return replace(report, gate=new_gate)


def format_gate_message(report: "PreflightReport") -> str:
    """Render a human-readable one-paragraph summary of why the gate fired.

    Used by the CLI to format both the stderr warning (under `warn`) and
    the exception message (under `strict`). Includes the canonical path,
    every offending input's path + classification + intersection
    fraction so the user can diagnose without reading the JSON.
    """
    failing = report.gate.get("failing_inputs", [])
    if not failing:
        return "preflight gate not triggered"

    lines: list[str] = [
        f"Preflight gate triggered against canonical {report.canonical['path']!r}"
        f" (variant_key={report.variant_key}):"
    ]
    by_index = {c.input_index: c for c in report.comparisons}
    for entry in failing:
        idx = entry["input_index"]
        pair = by_index[idx]
        cls = entry["classification"]
        hint = _CLASSIFICATION_HINTS.get(cls, "")
        lines.append(
            f"  input[{idx}] {pair.path}: "
            f"classification={cls}, "
            f"intersection={pair.intersection:,}/{min(pair.n_variants, report.canonical['n_variants']):,} "
            f"({pair.intersection_fraction_of_min:.1%} of min)."
            f"{(' ' + hint) if hint else ''}"
        )
    lines.append("See <prefix>.preflight.json for full evidence.")
    return "\n".join(lines)


# One-line user-facing hints per classification label. Kept here (vs. on
# the classifier) because they're CLI-presentation concerns, not pure
# classification semantics.
_CLASSIFICATION_HINTS: dict[str, str] = {
    "build_mismatch": (
        "Likely cause: coordinate-build mismatch (hg19 vs hg38) OR two "
        "unrelated panels that happen to cover the same chromosomes. "
        "Verify panel/build identity; liftover one side if builds differ."
    ),
    "key_space_mismatch": (
        "The non-active variant key matches well. Try the other "
        "--variant-key value (chr_pos / id)."
    ),
    "disjoint_panels": (
        "Chromosome coverage differs between inputs. Verify you're "
        "merging the panels you intended."
    ),
}


# Classifier thresholds. The gate (step 4) consumes the classification
# label + its own policy thresholds — the classifier itself just
# distinguishes shapes, not severity.
_COMPATIBLE_MIN_FRACTION = 0.5
# `key_space_mismatch` requires the alternate-key fraction to be
# substantially better than the active key — a small gap could just be
# noise from different upstream filters.
_KEY_SPACE_ALTERNATE_LIFT = 0.4


def classify_pair(pair: PairCompatibility) -> tuple[str, dict[str, Any]]:
    """Classify the failure shape of one (canonical, other) comparison.

    Pure function over `PairCompatibility` — no I/O, no policy reads —
    so step 4's gate can re-invoke or re-evaluate without re-running
    pass 1. Tested directly against the failure-mode corpus.

    Returns `(label, evidence)`. Labels:
      - `"compatible"`: active-key intersection_fraction_of_min >= 0.5.
        The negative control; no gate action.
      - `"key_space_mismatch"`: active-key fraction low *and* alternate-key
        fraction substantially higher (lift >= 0.4). Signals the user
        picked the wrong `--variant-key` for their data (e.g., chr_pos
        against an rsID-only panel).
      - `"build_mismatch"`: active-key fraction low, alternate-key
        comparable, *and* every shared chromosome has canonical_size>0,
        other_size>0, intersection=0, with symmetric chrom presence on
        both sides. The hg19/hg38 fingerprint. **Caveat:** this label
        also fires for two genuinely unrelated panels that happen to
        cover the same chromosome set with zero coordinate overlap (the
        common case for two human panels both restricted to autosomes
        1-22). Distinguishing those from a true build mismatch requires
        signal we don't yet compute (e.g., per-chrom position-shift
        consistency) — the actionable advice ("verify panel/build
        identity") is the same in both cases, so the label is
        conservative on purpose.
      - `"disjoint_panels"`: active-key fraction low with asymmetric
        chrom presence (chroms exist in one side but not the other).
        Catch-all for unrelated panels whose chrom coverage itself
        differs — a strictly clearer signal than `build_mismatch`.
      - `"empty_input"`: at least one side has zero post-filter variants.

    Evidence carries the numbers the classifier keyed on so users can
    see WHY a label was chosen (and, when wrong, file a bug with a
    specific number to argue against).
    """
    active_frac = pair.intersection_fraction_of_min
    alt_frac = pair.alternate_key_fraction_of_min or 0.0
    n_canonical_chroms = sum(1 for pc in pair.per_chrom if pc.canonical_size > 0)
    n_other_chroms = sum(1 for pc in pair.per_chrom if pc.other_size > 0)
    n_shared_chroms = sum(
        1 for pc in pair.per_chrom if pc.canonical_size > 0 and pc.other_size > 0
    )
    n_shared_chroms_zero_overlap = sum(
        1
        for pc in pair.per_chrom
        if pc.canonical_size > 0 and pc.other_size > 0 and pc.intersection == 0
    )
    n_chroms_only_canonical = sum(
        1 for pc in pair.per_chrom if pc.canonical_size > 0 and pc.other_size == 0
    )
    n_chroms_only_other = sum(
        1 for pc in pair.per_chrom if pc.other_size > 0 and pc.canonical_size == 0
    )

    evidence: dict[str, Any] = {
        "active_key_fraction": active_frac,
        "alternate_key_fraction": alt_frac,
        "n_canonical_chroms": n_canonical_chroms,
        "n_other_chroms": n_other_chroms,
        "n_shared_chroms": n_shared_chroms,
        "n_shared_chroms_zero_overlap": n_shared_chroms_zero_overlap,
        "n_chroms_only_canonical": n_chroms_only_canonical,
        "n_chroms_only_other": n_chroms_only_other,
    }

    # Empty-input check first — protects downstream ratios from being
    # mistaken for "compatible" when one side actually has no variants.
    if pair.n_variants == 0 or (
        pair.alternate_key_canonical_size is not None
        and pair.alternate_key_canonical_size == 0
    ):
        return "empty_input", evidence

    if active_frac >= _COMPATIBLE_MIN_FRACTION:
        return "compatible", evidence

    # Active-key fraction is low. Three distinguishable causes follow.

    # 1. Key-space mismatch — alternate key would have made the merge work.
    if alt_frac - active_frac >= _KEY_SPACE_ALTERNATE_LIFT:
        return "key_space_mismatch", evidence

    # 2. Build mismatch — chrom coverage is *symmetric* (both sides cover
    #    the same chroms) and every shared chrom is fully disjoint at the
    #    coord level. The hg19/hg38 fingerprint: same panel, different
    #    coordinates. Symmetry matters because asymmetric chrom presence
    #    (e.g., canonical has chr2 but other doesn't) is the disjoint-
    #    panels signature, not a build issue.
    if (
        n_shared_chroms > 0
        and n_shared_chroms_zero_overlap == n_shared_chroms
        and n_chroms_only_canonical == 0
        and n_chroms_only_other == 0
    ):
        return "build_mismatch", evidence

    # 3. Otherwise: genuinely disjoint panels.
    return "disjoint_panels", evidence


def _alternate_key(active: str) -> str:
    """Return the non-active variant key for the alternate-key view."""
    if active == "chr_pos":
        return "id"
    if active == "id":
        return "chr_pos"
    raise InvariantViolation(f"unknown variant_key: {active!r}")


# Placeholder values commonly used in the .pvar ID column when no real
# variant ID is known. We filter these out before forming the id-keyed
# set: a real-world .pvar with all-'.' IDs would otherwise collapse to
# a single-element set on both sides, giving an alternate-key fraction
# of 1.0 and triggering a spurious `key_space_mismatch` classification.
# See review finding #1.
_PLACEHOLDER_VARIANT_IDS: frozenset[str] = frozenset({".", "", "0", "NA", "nan", "None"})


def _keys_for(df: pd.DataFrame, variant_key: str) -> set[Any]:
    if variant_key == "chr_pos":
        return set(zip(df["chrom"].tolist(), df["pos"].tolist(), strict=True))
    if variant_key == "id":
        return {
            vid
            for vid in df["id"].astype(str).tolist()
            if vid not in _PLACEHOLDER_VARIANT_IDS
        }
    raise InvariantViolation(f"unknown variant_key: {variant_key!r}")


def _keys_by_chrom(df: pd.DataFrame, variant_key: str) -> dict[int, set[Any]]:
    """Per-chrom key sets, used for the per_chrom breakdown.

    Placeholder IDs are filtered in the `id` branch (see
    `_PLACEHOLDER_VARIANT_IDS`) so per-chrom sizes match the
    deduplicated-key set produced by `_keys_for`.
    """
    out: dict[int, set[Any]] = {}
    if variant_key == "chr_pos":
        for chrom, pos in zip(df["chrom"].tolist(), df["pos"].tolist(), strict=True):
            out.setdefault(int(chrom), set()).add((chrom, pos))
        return out
    if variant_key == "id":
        for chrom, vid in zip(df["chrom"].tolist(), df["id"].astype(str).tolist(), strict=True):
            if vid in _PLACEHOLDER_VARIANT_IDS:
                continue
            out.setdefault(int(chrom), set()).add(vid)
        return out
    raise InvariantViolation(f"unknown variant_key: {variant_key!r}")


def write_preflight_json(report: PreflightReport, path: Path) -> None:
    """Serialize `report` to `path` per schema v1."""
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(_to_dict(report), f, indent=2)
            f.write("\n")
    except OSError as e:
        raise IOFailure(f"cannot write {path}: {e}") from e


def _to_dict(report: PreflightReport) -> dict[str, Any]:
    return {
        "schema_version": report.schema_version,
        "tool": report.tool,
        "tool_version": report.tool_version,
        "command": report.command,
        "variant_key": report.variant_key,
        "canonical": dict(report.canonical),
        "comparisons": [_pair_to_dict(c) for c in report.comparisons],
        "gate": dict(report.gate),
    }


def _pair_to_dict(c: PairCompatibility) -> dict[str, Any]:
    return {
        "input_index": int(c.input_index),
        "path": str(c.path),
        "n_variants": int(c.n_variants),
        "intersection": int(c.intersection),
        "intersection_fraction_of_min": float(c.intersection_fraction_of_min),
        "per_chrom": [
            {
                "chrom": int(pc.chrom),
                "canonical_size": int(pc.canonical_size),
                "other_size": int(pc.other_size),
                "intersection": int(pc.intersection),
            }
            for pc in c.per_chrom
        ],
        "alternate_key": c.alternate_key,
        "alternate_key_canonical_size": (
            None if c.alternate_key_canonical_size is None else int(c.alternate_key_canonical_size)
        ),
        "alternate_key_other_size": (
            None if c.alternate_key_other_size is None else int(c.alternate_key_other_size)
        ),
        "alternate_key_intersection": (
            None if c.alternate_key_intersection is None else int(c.alternate_key_intersection)
        ),
        "alternate_key_fraction_of_min": (
            None
            if c.alternate_key_fraction_of_min is None
            else float(c.alternate_key_fraction_of_min)
        ),
        "classification": c.classification,
        "classification_evidence": c.classification_evidence,
    }


__all__ = (
    "PREFLIGHT_SCHEMA_VERSION",
    "GATE_FAILURE_LABELS",
    "PerChromCompat",
    "PairCompatibility",
    "PreflightReport",
    "compute_preflight",
    "classify_pair",
    "evaluate_gate",
    "format_gate_message",
    "write_preflight_json",
)
