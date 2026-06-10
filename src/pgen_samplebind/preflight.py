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

import numpy as np
import pandas as pd

from .errors import InvariantViolation, IOFailure
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
    # Per-chrom position-shift consistency: distinguishes a true coordinate-
    # build mismatch (uniform shift within each chrom, MAD ~0) from two
    # unrelated panels that happen to cover the same chromosomes
    # (shifts look like noise). Computed only for shared chroms with
    # zero coord-level intersection AND >=5 variants on both sides;
    # `None` when no chroms qualify. See `_compute_build_shift_signature`.
    build_shift_signature: dict[str, Any] | None = None
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
    alternate_key = _alternate_key(policy.variant_key)
    # Build the canonical key universes once (active + alternate). For
    # chr_pos these are sorted-unique int64 code arrays; for id, str sets.
    canonical_univ = _key_universe(canonical_df, policy.variant_key)
    canonical_alt_univ = _key_universe(canonical_df, alternate_key)
    n_canonical = _key_count(canonical_univ)
    n_canonical_alt = _key_count(canonical_alt_univ)
    # Extract canonical chrom/pos numpy arrays once (used by every
    # _compute_build_shift_signature call in the loop). At 84M-variant
    # canonical scale these are ~336 MB each; materializing them once
    # avoids re-pulling per non-canonical input. For N=2 merges this is
    # a wash, but N>2 (multi-input cohort assembly) sees real wins.
    canonical_chrom_arr = canonical_df["chrom"].to_numpy()
    canonical_pos_arr = canonical_df["pos"].to_numpy()
    # Canonical per-chrom key structure, computed once (loop-invariant).
    # chr_pos → {chrom: unique-key-count}; id → {chrom: set[id]}. Was
    # rebuilt inside _per_chrom_compat on every non-canonical input — an
    # N-1x redundant np.unique / Python-bucketing pass for N-way merges.
    canonical_per_chrom = _prepare_canonical_per_chrom(
        policy.variant_key, canonical_df, canonical_univ
    )
    # Lazily-filled {chrom: sorted canonical positions} shared across all
    # build-shift calls so a given chrom's canonical column is masked and
    # sorted at most once for the whole merge, not once per non-canonical
    # input (the canonical positions don't change between comparisons).
    canonical_sorted_pos_cache: dict[int, np.ndarray[Any, Any]] = {}

    comparisons: list[PairCompatibility] = []
    for i, desc in enumerate(descriptors[1:], start=1):
        other_df = pvars[i] if pvars is not None else read_pvar(desc.pvar_path)
        other_univ = _key_universe(other_df, policy.variant_key)
        n_other = _key_count(other_univ)

        # Compute the active-key intersection ONCE and reuse it for both the
        # scalar count and the per-chrom breakdown. Previously the same
        # intersection was computed twice per comparison — once here for the
        # count, once inside _per_chrom_compat — doubling the dominant
        # np.intersect1d cost on every merge (not just N>2).
        inter = _intersect(canonical_univ, other_univ)
        intersection = _key_count(inter)
        denom = min(n_canonical, n_other)
        fraction = (intersection / denom) if denom > 0 else 0.0

        per_chrom = _per_chrom_compat(
            policy.variant_key, canonical_per_chrom, other_df, other_univ, inter
        )

        # Alternate-key view, for key_space_mismatch detection (count only).
        other_alt_univ = _key_universe(other_df, alternate_key)
        alt_intersection = _intersection_count(canonical_alt_univ, other_alt_univ)
        n_other_alt = _key_count(other_alt_univ)
        alt_denom = min(n_canonical_alt, n_other_alt)
        alt_fraction = (alt_intersection / alt_denom) if alt_denom > 0 else 0.0

        # Build-shift signature: only computed when the per-chrom shape
        # could plausibly indicate a build mismatch (shared chroms with
        # zero coord overlap). Cheap when nothing qualifies — early-out
        # inside the helper. Canonical arrays + the sorted-pos cache are
        # hoisted out of the loop (see above); other arrays are
        # per-iteration anyway.
        build_shift = _compute_build_shift_signature(
            canonical_chrom_arr,
            canonical_pos_arr,
            other_df["chrom"].to_numpy(),
            other_df["pos"].to_numpy(),
            per_chrom,
            canonical_sorted_pos_cache,
        )

        pair = PairCompatibility(
            input_index=i,
            path=desc.path,
            n_variants=n_other,
            intersection=intersection,
            intersection_fraction_of_min=fraction,
            per_chrom=per_chrom,
            alternate_key=alternate_key,
            alternate_key_canonical_size=n_canonical_alt,
            alternate_key_other_size=n_other_alt,
            alternate_key_intersection=alt_intersection,
            alternate_key_fraction_of_min=alt_fraction,
            build_shift_signature=build_shift,
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
            "n_variants": n_canonical,
        },
        comparisons=tuple(comparisons),
        gate={"triggered": False, "action": "none", "threshold": None},
    )


# Labels that the gate treats as actionable failures. `compatible` is
# the happy path. The four failure labels all produce near-empty merges
# without preflight intervention — `empty_input` (one side has zero
# post-filter variants) was previously excluded under the assumption
# that downstream alignment gates (a)-(d) would catch it, but they
# don't: gate (a) computes 0 extras against an empty other (silent),
# gate (b) is guarded by `if intersection_size > 0`, gate (c) only
# fires under --target, and gate (d) iterates the empty input.
# Including empty_input here means `--preflight-policy strict` actually
# protects against the most extreme near-empty-merge shape. Closes
# review finding #B6.
GATE_FAILURE_LABELS: frozenset[str] = frozenset(
    {"build_mismatch", "key_space_mismatch", "disjoint_panels", "empty_input"}
)


def evaluate_gate(report: PreflightReport, policy: MergePolicy) -> PreflightReport:
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

    if not would_trigger or policy.preflight_policy == "off":
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


def format_gate_message(report: PreflightReport) -> str:
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
        "Likely cause: coordinate-build mismatch (hg19 vs hg38) — per-chrom "
        "positions on both sides differ by a near-uniform shift, the "
        "signature of a coordinate remap. Liftover one side to the other's "
        "build (CrossMap or Picard LiftoverVcf) and re-run."
    ),
    "key_space_mismatch": (
        "The non-active variant key matches well. Try the other --variant-key value (chr_pos / id)."
    ),
    "disjoint_panels": (
        "Inputs appear to come from unrelated cohorts. Either chromosome "
        "coverage differs between inputs, or the chrom sets coincide but "
        "coordinates don't (random positions, not a uniform shift). Run "
        "`pgen-samplebind hash` on each input to verify panel identity; "
        "if both inputs SHOULD be related, suspect a coordinate-build "
        "mismatch on too few variants to detect (try liftover with "
        "CrossMap or Picard LiftoverVcf)."
    ),
    "empty_input": (
        "One side has zero post-filter variants. Common causes: input "
        "contains only multi-allelic / non-SNP rows (filtered by the "
        "biallelic-SNP gate), upstream filter emitted an empty .pvar, "
        "or `--include-chrom` excluded everything. Verify the input with "
        "`pgen-samplebind inspect`."
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
        comparable, every shared chromosome has canonical_size>0,
        other_size>0, intersection=0, with symmetric chrom presence on
        both sides, AND either (a) the per-chrom position-shift
        signature shows a consistent shift (has_consistent_shift=True
        — the hg19/hg38 fingerprint on panels large enough to compute
        the signature), OR (b) the signature couldn't be computed
        because no shared chrom has >=5 variants on both sides (small
        targeted panel — conservatively assume build mismatch since
        liftover is a cheap remediation to try and the alternative
        label would steer the user away from the real fix). Two
        unrelated panels with same chrom set, random independent
        positions, AND sufficient density compute a non-consistent
        signature and fall through to `disjoint_panels`.
      - `"disjoint_panels"`: active-key fraction low with either
        asymmetric chrom presence (chroms exist in one side but not
        the other) OR symmetric chroms with zero overlap but no
        consistent position shift (random positions, not a build
        remap). Covers genuinely unrelated panels regardless of
        whether their chrom sets coincide.
      - `"empty_input"`: at least one side has zero post-filter variants.

    Evidence carries the numbers the classifier keyed on so users can
    see WHY a label was chosen (and, when wrong, file a bug with a
    specific number to argue against).
    """
    active_frac = pair.intersection_fraction_of_min
    alt_frac = pair.alternate_key_fraction_of_min or 0.0
    n_canonical_chroms = sum(1 for pc in pair.per_chrom if pc.canonical_size > 0)
    n_other_chroms = sum(1 for pc in pair.per_chrom if pc.other_size > 0)
    n_shared_chroms = sum(1 for pc in pair.per_chrom if pc.canonical_size > 0 and pc.other_size > 0)
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
    # `canonical_active_size` is derived from per_chrom (which always
    # covers the active key) rather than `alternate_key_canonical_size`:
    # a canonical .pvar with all-placeholder IDs ('.', 'NA', etc.) under
    # the default --variant-key chr_pos has alternate_key_canonical_size=0
    # but real chr_pos data — keying the empty_input guard on the
    # alternate key would mis-label every such input as empty. Closes
    # review finding #C1.
    canonical_active_size = sum(int(pc.canonical_size) for pc in pair.per_chrom)
    evidence["canonical_active_size"] = canonical_active_size
    if pair.n_variants == 0 or canonical_active_size == 0:
        return "empty_input", evidence

    if active_frac >= _COMPATIBLE_MIN_FRACTION:
        return "compatible", evidence

    # Active-key fraction is low. Three distinguishable causes follow.

    # 1. Key-space mismatch — alternate key would have made the merge work.
    if alt_frac - active_frac >= _KEY_SPACE_ALTERNATE_LIFT:
        return "key_space_mismatch", evidence

    # 2. Build mismatch — chrom coverage is *symmetric* (both sides cover
    #    the same chroms), every shared chrom is fully disjoint at the
    #    coord level, AND the per-chrom position-shift signature shows
    #    a consistent translation. The shift-consistency check is what
    #    separates a true coordinate-build mismatch from two unrelated
    #    panels that happen to share a chrom set (random positions →
    #    no consistent shift → falls through to `disjoint_panels`).
    symmetric_zero_overlap = (
        n_shared_chroms > 0
        and n_shared_chroms_zero_overlap == n_shared_chroms
        and n_chroms_only_canonical == 0
        and n_chroms_only_other == 0
    )
    if symmetric_zero_overlap:
        sig = pair.build_shift_signature
        # Echo the shift evidence so users see WHY the classifier picked
        # `build_mismatch` vs `disjoint_panels` — the most-asked question
        # when this label fires.
        if sig is not None:
            evidence["build_shift_n_chroms_evaluated"] = len(sig["chroms_evaluated"])
            evidence["build_shift_n_consistent"] = sig["n_consistent_shift_chroms"]
            evidence["build_shift_has_consistent_shift"] = sig["has_consistent_shift"]
            if sig["has_consistent_shift"]:
                return "build_mismatch", evidence
            # Same-chrom no-overlap WITH signal but no consistent shift
            # → genuinely disjoint, fall through.
        else:
            # No shift signal could be computed (panel too small for
            # >=5 variants/chrom on every shared chrom). Conservatively
            # label as build_mismatch: the symmetric-zero-overlap shape
            # is the hg19/hg38 fingerprint, and recommending a liftover
            # check costs the user nothing if it turns out to be a
            # genuinely-disjoint small panel. The pre-sharpener default
            # for this case was build_mismatch — preserving it ensures
            # small targeted panels (e.g., 30-variant pharmacogenomic
            # panels) don't lose the liftover remediation hint. Closes
            # review finding #B1.
            evidence["build_shift_signature_available"] = False
            return "build_mismatch", evidence

    # 3. Otherwise (asymmetric chroms OR same-chrom-no-shift): disjoint.
    return "disjoint_panels", evidence


# --- Build-shift signature ----------------------------------------------
#
# Distinguishes a true hg19/hg38-style coordinate-build mismatch (uniform
# per-chrom position shift) from two unrelated panels that happen to
# cover the same chromosome set (positions look like noise). Only
# evaluated for shared chroms with zero coord-level intersection where
# both sides have enough variants for a stable median (>=5). For each
# qualifying chrom we sort positions on both sides, rank-align them, and
# look at:
#   - median(other_sorted[i] - canonical_sorted[i]) → the per-chrom shift
#   - MAD(shifts) / max(|median_shift|, 1) → "relative MAD"
# Build mismatch: relative MAD ≈ 0 (shift is uniform within a chrom).
# Disjoint panels: relative MAD is large (sorted-position differences
# are dominated by sampling noise, not a real translation).
_BUILD_SHIFT_MIN_VARIANTS_PER_CHROM = 5
_BUILD_SHIFT_MAX_RELATIVE_MAD = 0.1  # smaller = more uniform
_BUILD_SHIFT_MIN_MAGNITUDE = 1000  # bp; below this is sampling noise
_BUILD_SHIFT_MIN_FRACTION_CHROMS = 0.5  # at least half of evaluated chroms
# must show a consistent shift to call it build_mismatch


def _compute_build_shift_signature(
    canonical_chrom: np.ndarray[Any, Any],
    canonical_pos: np.ndarray[Any, Any],
    other_chrom: np.ndarray[Any, Any],
    other_pos: np.ndarray[Any, Any],
    per_chrom: tuple[PerChromCompat, ...],
    canonical_sorted_pos_cache: dict[int, np.ndarray[Any, Any]] | None = None,
) -> dict[str, Any] | None:
    """Per-chrom rank-aligned position-shift summary.

    Returns None when no chroms qualify (no shared zero-overlap chrom
    has >=_BUILD_SHIFT_MIN_VARIANTS_PER_CHROM on both sides). When some
    chroms qualify, returns a dict with:
      - chroms_evaluated: list[int]
      - median_shift_per_chrom: dict[str, float]  (int keys stringified
        for JSON consumers)
      - relative_mad_per_chrom: dict[str, float]
      - n_consistent_shift_chroms: int  (chroms passing both the
        magnitude and relative-MAD thresholds)
      - has_consistent_shift: bool  (>= _BUILD_SHIFT_MIN_FRACTION_CHROMS
        of evaluated chroms are consistent — the classifier's
        build-mismatch trigger)

    Pure helper. Takes pre-extracted numpy arrays for both sides (not
    pandas DataFrames) so the canonical arrays can be materialized once
    by `compute_preflight` and reused across every non-canonical input
    in an N-way merge — at 84M-variant canonical scale that's a
    ~336 MB array we avoid re-extracting from the pandas Series on each
    iteration.
    """
    qualifying: list[int] = [
        int(pc.chrom)
        for pc in per_chrom
        if pc.canonical_size >= _BUILD_SHIFT_MIN_VARIANTS_PER_CHROM
        and pc.other_size >= _BUILD_SHIFT_MIN_VARIANTS_PER_CHROM
        and pc.intersection == 0
    ]
    if not qualifying:
        return None

    median_shift_per_chrom: dict[str, float] = {}
    relative_mad_per_chrom: dict[str, float] = {}
    n_consistent = 0

    # Local aliases to keep the loop body readable.
    can_chrom = canonical_chrom
    can_pos = canonical_pos
    oth_chrom = other_chrom
    oth_pos = other_pos
    # The canonical sorted positions for a chrom are invariant across
    # non-canonical inputs; cache them (when the caller supplies a cache)
    # so each chrom's canonical mask+sort runs at most once per merge.
    cache = canonical_sorted_pos_cache

    for chrom in qualifying:
        if cache is not None and chrom in cache:
            can_sorted = cache[chrom]
        else:
            can_sorted = np.sort(can_pos[can_chrom == chrom])
            if cache is not None:
                cache[chrom] = can_sorted
        oth_sorted = np.sort(oth_pos[oth_chrom == chrom])
        n = min(len(can_sorted), len(oth_sorted))
        # rank-align: i-th smallest on each side
        shifts = oth_sorted[:n].astype(np.int64) - can_sorted[:n].astype(np.int64)
        median = float(np.median(shifts))
        mad = float(np.median(np.abs(shifts - median)))
        relative_mad = mad / max(abs(median), 1.0)
        median_shift_per_chrom[str(chrom)] = median
        relative_mad_per_chrom[str(chrom)] = relative_mad
        if (
            abs(median) >= _BUILD_SHIFT_MIN_MAGNITUDE
            and relative_mad <= _BUILD_SHIFT_MAX_RELATIVE_MAD
        ):
            n_consistent += 1

    # "Consistent shift across enough chroms" — protects against one
    # chrom coincidentally lining up while others stay noise. A pair
    # with a single qualifying chrom needs that chrom to be consistent.
    min_required = max(1, int(_BUILD_SHIFT_MIN_FRACTION_CHROMS * len(qualifying) + 0.5))
    has_consistent_shift = n_consistent >= min_required

    return {
        "chroms_evaluated": qualifying,
        "median_shift_per_chrom": median_shift_per_chrom,
        "relative_mad_per_chrom": relative_mad_per_chrom,
        "n_consistent_shift_chroms": n_consistent,
        "has_consistent_shift": has_consistent_shift,
    }


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


# (chrom, pos) packed into one int64: chrom in the high bits, pos in the
# low _CHROM_SHIFT bits. Lets the chr_pos key-space be represented as a
# sorted-unique int64 array so intersection / per-chrom counts run as
# vectorized numpy (np.unique / np.intersect1d) instead of building a
# Python set of 2-tuples — the latter was the dominant cost of the
# preflight pass at panel scale (issue #12 perf follow-up). 40-bit pos
# field holds any genomic coordinate (< 2^40 ≈ 1.1e12 bp) collision-free,
# and the max chrom (26) << 40 ≈ 2.9e13 stays well under the int64 ceiling.
_CHROM_SHIFT = 40


def _unique_chr_pos_codes(df: pd.DataFrame) -> np.ndarray[Any, Any]:
    """Encode each (chrom, pos) row as one int64 and return sorted-unique codes."""
    chrom = np.asarray(df["chrom"], dtype=np.int64)
    pos = np.asarray(df["pos"], dtype=np.int64)
    if pos.size and (pos.max() >= (1 << _CHROM_SHIFT) or pos.min() < 0):
        raise InvariantViolation(
            f"position out of range for chr_pos key encoding: "
            f"max={int(pos.max())}, min={int(pos.min())} (expected 0 <= pos < 2^{_CHROM_SHIFT})"
        )
    return np.unique((chrom << _CHROM_SHIFT) | pos)


def _counts_by_chrom_from_codes(codes: np.ndarray[Any, Any]) -> dict[int, int]:
    """Per-chrom unique-key counts from an int64 code array (chrom = code >> shift)."""
    if codes.size == 0:
        return {}
    uniq, counts = np.unique(codes >> _CHROM_SHIFT, return_counts=True)
    return dict(zip(uniq.tolist(), counts.tolist(), strict=True))


# A "key universe" is the deduplicated set of variant keys under one
# variant_key. chr_pos → sorted-unique int64 code array (vectorized);
# id → set[str] (arbitrary strings; placeholder IDs filtered). Both
# support _key_count / _intersection_count uniformly.
def _key_universe(df: pd.DataFrame, variant_key: str) -> Any:
    if variant_key == "chr_pos":
        return _unique_chr_pos_codes(df)
    if variant_key == "id":
        return {vid for vid in df["id"].astype(str).tolist() if vid not in _PLACEHOLDER_VARIANT_IDS}
    raise InvariantViolation(f"unknown variant_key: {variant_key!r}")


def _key_count(universe: Any) -> int:
    return int(universe.size) if isinstance(universe, np.ndarray) else len(universe)


def _intersect(a: Any, b: Any) -> Any:
    """Intersection in the universe's own representation: a sorted-unique
    int64 code array (chr_pos) or a set (id). Callers that only need the
    size go through `_intersection_count`; callers that also need the
    per-chrom breakdown reuse the returned object so the intersection is
    computed once."""
    if isinstance(a, np.ndarray):
        # Both arrays are sorted-unique → assume_unique is safe and faster.
        return np.intersect1d(a, b, assume_unique=True)
    return a & b


def _intersection_count(a: Any, b: Any) -> int:
    return _key_count(_intersect(a, b))


def _prepare_canonical_per_chrom(
    variant_key: str, canonical_df: pd.DataFrame, canonical_univ: Any
) -> Any:
    """Canonical-side per-chrom structure, computed once per merge (it is
    invariant across non-canonical inputs). chr_pos → {chrom: count} from
    the already-built code universe; id → {chrom: set[id]} buckets."""
    if variant_key == "chr_pos":
        return _counts_by_chrom_from_codes(canonical_univ)
    if variant_key == "id":
        return _ids_by_chrom(canonical_df)
    raise InvariantViolation(f"unknown variant_key: {variant_key!r}")


def _per_chrom_chr_pos(
    can_counts: dict[int, int],
    oth_codes: np.ndarray[Any, Any],
    inter_codes: np.ndarray[Any, Any],
) -> tuple[PerChromCompat, ...]:
    """Per-chrom breakdown for the chr_pos key. `can_counts` is the
    precomputed canonical {chrom: count}; `inter_codes` is the active-key
    intersection already computed by the caller (reused, not recomputed)."""
    oth_counts = _counts_by_chrom_from_codes(oth_codes)
    inter_counts = _counts_by_chrom_from_codes(inter_codes)
    return tuple(
        PerChromCompat(
            chrom=chrom,
            canonical_size=can_counts.get(chrom, 0),
            other_size=oth_counts.get(chrom, 0),
            intersection=inter_counts.get(chrom, 0),
        )
        for chrom in sorted(set(can_counts) | set(oth_counts))
    )


def _per_chrom_id(
    can_by: dict[int, set[str]], other_df: pd.DataFrame
) -> tuple[PerChromCompat, ...]:
    """Per-chrom breakdown for the id key. `can_by` is the precomputed
    canonical {chrom: set[id]}; the other side is bucketed and intersected
    per chrom — kept set-based because id keys are arbitrary strings, not
    the perf-critical path."""
    oth_by = _ids_by_chrom(other_df)
    return tuple(
        PerChromCompat(
            chrom=chrom,
            canonical_size=len(can_by.get(chrom, set())),
            other_size=len(oth_by.get(chrom, set())),
            intersection=len(can_by.get(chrom, set()) & oth_by.get(chrom, set())),
        )
        for chrom in sorted(set(can_by) | set(oth_by))
    )


def _ids_by_chrom(df: pd.DataFrame) -> dict[int, set[str]]:
    out: dict[int, set[str]] = {}
    for chrom, vid in zip(df["chrom"].tolist(), df["id"].astype(str).tolist(), strict=True):
        if vid in _PLACEHOLDER_VARIANT_IDS:
            continue
        out.setdefault(int(chrom), set()).add(vid)
    return out


def _per_chrom_compat(
    variant_key: str,
    canonical_per_chrom: Any,
    other_df: pd.DataFrame,
    oth_universe: Any,
    inter: Any,
) -> tuple[PerChromCompat, ...]:
    """Per-chrom compatibility for the active key, reusing the canonical
    structure (`canonical_per_chrom`) and the already-computed intersection
    (`inter`) from the caller. chr_pos uses the int64 code path; id falls
    back to the set-based path (where `inter` is unused — per-chrom id
    intersection is a per-bucket set-AND)."""
    if variant_key == "chr_pos":
        return _per_chrom_chr_pos(canonical_per_chrom, oth_universe, inter)
    if variant_key == "id":
        return _per_chrom_id(canonical_per_chrom, other_df)
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
        "build_shift_signature": c.build_shift_signature,
        "classification": c.classification,
        "classification_evidence": c.classification_evidence,
    }


__all__ = (
    "GATE_FAILURE_LABELS",
    "PREFLIGHT_SCHEMA_VERSION",
    "PairCompatibility",
    "PerChromCompat",
    "PreflightReport",
    "classify_pair",
    "compute_preflight",
    "evaluate_gate",
    "format_gate_message",
    "write_preflight_json",
)
