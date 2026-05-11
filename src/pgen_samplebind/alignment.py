"""Pass 1: variant matching, allele resolution, action assignment, gate evaluation.

Per LLD §3.9. The alignment table is the single source of truth for what pass 2
will do (`merge.merge_inputs`) and what reports will say (`reporting.write_report_*`).
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import pvar as pvar_module
from .errors import InvariantViolation, ValidationError
from .types import AlignmentSummary, DropReason, MergeAction, MergePolicy

# Watson-Crick complement.
_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _is_ambiguous_pair(ref: str, alt: str) -> bool:
    """True if the (ref, alt) pair is A/T or C/G — a strand-ambiguous SNP."""
    return _COMPLEMENT.get(ref) == alt


def resolve_alleles(
    canonical_ref: str,
    canonical_alt: str,
    other_ref: str,
    other_alt: str,
    trust_strand: bool,
) -> tuple[MergeAction, DropReason | None]:
    """Apply HLD §Allele resolution truth table for a single (canonical, other) pair.

    Returns (action, drop_reason). drop_reason is None except when action is DROP.

    Behavior summary (HLD truth table v2):
      same alleles + same orientation        → PASSTHROUGH
      same alleles + swapped orientation     → REF_ALT_SWAP
      complement-of-other matches canonical  → STRAND_FLIP
      complement-of-other matches swapped    → STRAND_FLIP_AND_SWAP
      ambiguous (A/T or C/G), any order      → DROP(ambiguous_strand)
                                               or PASSTHROUGH under trust_strand
      anything else                          → DROP(allele_mismatch)

    Note on `--trust-strand` for ambiguous SNPs: by default pgen-samplebind
    drops A/T and C/G ambiguous matches because strand cannot be verified
    even when alleles agree (an A/T pair on the forward strand looks
    identical to an A/T pair on the reverse strand — complementing A/T
    gives T/A which is the same pair). This matches mergeit's `strandcheck`
    behavior and is the safer default for cross-source merges. Pass
    `--trust-strand` for single-source merges where REF/ALT calls are
    guaranteed consistent across inputs even on ambiguous SNPs.
    """
    # Defensive: any non-ACGT input is an allele mismatch.
    if not (
        canonical_ref in _COMPLEMENT
        and canonical_alt in _COMPLEMENT
        and other_ref in _COMPLEMENT
        and other_alt in _COMPLEMENT
    ):
        return MergeAction.DROP, DropReason.ALLELE_MISMATCH

    if _is_ambiguous_pair(canonical_ref, canonical_alt):
        canonical_set = frozenset((canonical_ref, canonical_alt))
        other_set = frozenset((other_ref, other_alt))
        if canonical_set != other_set:
            return MergeAction.DROP, DropReason.ALLELE_MISMATCH
        # Both inputs ambiguous at same site — strand undetermined regardless
        # of REF/ALT order. Drop unless user explicitly trusts strand.
        if trust_strand:
            return MergeAction.PASSTHROUGH, None
        return MergeAction.DROP, DropReason.AMBIGUOUS_STRAND

    # Non-ambiguous canonical.
    if (canonical_ref, canonical_alt) == (other_ref, other_alt):
        return MergeAction.PASSTHROUGH, None
    if (canonical_ref, canonical_alt) == (other_alt, other_ref):
        return MergeAction.REF_ALT_SWAP, None

    other_ref_comp = _COMPLEMENT[other_ref]
    other_alt_comp = _COMPLEMENT[other_alt]
    if (canonical_ref, canonical_alt) == (other_ref_comp, other_alt_comp):
        return MergeAction.STRAND_FLIP, None
    if (canonical_ref, canonical_alt) == (other_alt_comp, other_ref_comp):
        return MergeAction.STRAND_FLIP_AND_SWAP, None

    return MergeAction.DROP, DropReason.ALLELE_MISMATCH


def _key_set(df: pd.DataFrame, variant_key: str) -> set[tuple[int, int] | str]:
    """Build a hashable key set from a pvar DataFrame for membership checks."""
    if variant_key == "chr_pos":
        return set(zip(df["chrom"].tolist(), df["pos"].tolist(), strict=True))
    return set(df["id"].tolist())


def _classify_per_input_action(
    c_ref: str,
    c_alt: str,
    o_ref: object,
    o_alt: object,
    trust_strand: bool,
) -> tuple[MergeAction, DropReason | None]:
    """Per-row classification handling the missing-in-other case.

    Returns the *raw* action — policy-error gating is applied by the caller.

    Detects "missing in other" by checking the dtype of the merge cell: after
    pandas left-join on an object-dtype string column, missing rows have NaN
    (a float). Concretely: anything not a `str` is treated as missing.
    """
    if not isinstance(o_ref, str) or not isinstance(o_alt, str):
        return MergeAction.FILL_MISSING, None
    return resolve_alleles(c_ref, c_alt, o_ref, o_alt, trust_strand)


def _apply_on_missing_policy(
    raw_action: MergeAction,
    raw_reason: DropReason | None,
    policy: MergePolicy,
    soften_policy_errors: bool,
    triggers: dict[str, int],
    input_idx: int,
) -> tuple[MergeAction, DropReason | None]:
    """For FILL_MISSING raw actions, apply --on-missing policy."""
    if raw_action != MergeAction.FILL_MISSING:
        return raw_action, raw_reason
    if policy.on_missing == "fill_missing":
        return MergeAction.FILL_MISSING, None
    if policy.on_missing == "drop_variant":
        return MergeAction.DROP, DropReason.ON_MISSING_DROP_VARIANT
    # error
    if soften_policy_errors:
        triggers["on_missing_count"] += 1
        return MergeAction.FILL_MISSING, None
    raise InvariantViolation(
        f"--on-missing error: variant absent in input[{input_idx}]. "
        f"Use --on-missing fill_missing (default; -9 for absent samples) or "
        f"drop_variant (drop the variant from output) to allow merge to proceed."
    )


def _apply_on_mismatch_policy(
    raw_action: MergeAction,
    raw_reason: DropReason | None,
    policy: MergePolicy,
    soften_policy_errors: bool,
    triggers: dict[str, int],
    input_idx: int,
) -> tuple[MergeAction, DropReason | None]:
    """For DROP/ALLELE_MISMATCH raw actions, apply --on-mismatch policy."""
    if raw_action != MergeAction.DROP or raw_reason != DropReason.ALLELE_MISMATCH:
        return raw_action, raw_reason
    if policy.on_mismatch == "drop":
        return raw_action, raw_reason
    # error
    if soften_policy_errors:
        triggers["on_mismatch_count"] += 1
        return raw_action, raw_reason
    raise InvariantViolation(
        f"--on-mismatch error: allele mismatch at variant in input[{input_idx}]. "
        f"Use --on-mismatch drop (default; drop mismatched variants) to allow "
        f"merge to proceed; the per-variant report (--report PATH) details which."
    )


def build_alignment_table(
    canonical_pvar: pd.DataFrame,
    other_pvars: list[pd.DataFrame],
    policy: MergePolicy,
    summary: AlignmentSummary,
    soften_policy_errors: bool = False,
) -> pd.DataFrame:
    """Pass 1 main entry. See LLD §3.9.

    Returns the alignment DataFrame with schema in LLD §2.5.
    """
    pvar_module.validate_unique_keys(canonical_pvar, policy.variant_key)

    # Initialize table from canonical (renamed `id` → `variant_id` per LLD §2.5).
    # `cm` carries genetic-position (centiMorgans) through the pipeline so
    # the output .pvar preserves it for Morgan-spaced jackknife consumers.
    cm_col = canonical_pvar["cm"] if "cm" in canonical_pvar.columns else 0.0
    table = canonical_pvar[["chrom", "pos", "id", "ref", "alt"]].copy()
    table["cm"] = cm_col
    table = table.rename(columns={"id": "variant_id"})
    table["canonical_idx"] = np.arange(len(table), dtype=np.int64)

    merge_keys = ["chrom", "pos"] if policy.variant_key == "chr_pos" else ["id"]
    canonical_key_set = _key_set(canonical_pvar, policy.variant_key)

    triggers = {"on_mismatch_count": 0, "on_missing_count": 0, "on_extra_count": 0}
    action_categories = [a.value for a in MergeAction]
    drop_categories = [r.value for r in DropReason]

    for i, other_pvar in enumerate(other_pvars):
        input_idx = i + 1  # input[0] is canonical; non-canonical are 1, 2, ...

        # Tag other rows with their original index for downstream pass 2.
        other_indexed = other_pvar.assign(_other_idx=np.arange(len(other_pvar), dtype=np.int64))
        merged = pd.merge(
            canonical_pvar[[*merge_keys, "ref", "alt"]],
            other_indexed[[*merge_keys, "_other_idx", "ref", "alt"]].rename(
                columns={"ref": "_other_ref", "alt": "_other_alt"}
            ),
            on=merge_keys,
            how="left",
        )

        actions: list[str] = []
        reasons: list[str | None] = []

        for c_ref, c_alt, o_ref, o_alt in zip(
            merged["ref"], merged["alt"], merged["_other_ref"], merged["_other_alt"], strict=True
        ):
            raw_action, raw_reason = _classify_per_input_action(
                c_ref, c_alt, o_ref, o_alt, policy.trust_strand
            )
            action, reason = _apply_on_missing_policy(
                raw_action, raw_reason, policy, soften_policy_errors, triggers, input_idx
            )
            action, reason = _apply_on_mismatch_policy(
                action, reason, policy, soften_policy_errors, triggers, input_idx
            )
            actions.append(action.value)
            reasons.append(reason.value if reason is not None else None)

        table[f"action_input_{input_idx}"] = pd.Categorical(actions, categories=action_categories)
        table[f"drop_reason_input_{input_idx}"] = pd.Categorical(
            reasons, categories=drop_categories
        )
        table[f"idx_input_{input_idx}"] = (
            merged["_other_idx"].astype("Int64").reset_index(drop=True)
        )

        # Extras: variants in input[N] absent from input[0].
        other_keys = _key_set(other_pvar, policy.variant_key)
        extras_count = len(other_keys - canonical_key_set)
        if extras_count > 0:
            if policy.on_extra == "error":
                if soften_policy_errors:
                    triggers["on_extra_count"] += extras_count
                else:
                    raise InvariantViolation(
                        f"--on-extra error: {extras_count} variants in input[{input_idx}] "
                        f"absent from input[0]"
                    )
            summary.n_extras_dropped += extras_count

        _tally_actions_to_summary(actions, reasons, summary)

    if soften_policy_errors:
        for k, v in triggers.items():
            summary.policy_error_triggers[k] = v

    # Reorder columns: canonical_idx first, then chrom/pos/variant_id/ref/alt/cm,
    # then per-input action/reason/idx triplets.
    leading = ["canonical_idx", "chrom", "pos", "variant_id", "ref", "alt", "cm"]
    others = [c for c in table.columns if c not in leading]
    return table[leading + others]


def _tally_actions_to_summary(
    actions: list[str], reasons: list[str | None], summary: AlignmentSummary
) -> None:
    """Update summary action-bucket counts in place."""
    for action_str in actions:
        if action_str == MergeAction.PASSTHROUGH.value:
            summary.n_passthrough += 1
        elif action_str == MergeAction.REF_ALT_SWAP.value:
            summary.n_ref_alt_swap += 1
        elif action_str in (
            MergeAction.STRAND_FLIP.value,
            MergeAction.STRAND_FLIP_AND_SWAP.value,
        ):
            summary.n_strand_flip += 1
        elif action_str == MergeAction.DROP.value:
            summary.n_dropped += 1
        elif action_str == MergeAction.FILL_MISSING.value:
            summary.n_fill_missing += 1
    for reason_str in reasons:
        if reason_str is not None:
            reason_enum = DropReason(reason_str)
            summary.n_dropped_by_reason[reason_enum] = (
                summary.n_dropped_by_reason.get(reason_enum, 0) + 1
            )


def count_kept_variants(table: pd.DataFrame) -> int:
    """Variants where no per-input action is DROP. Per LLD §3.9.

    Used by `merge.merge_inputs` to construct PgenWriter with the exact
    variant_ct (HLD §Two-pass merge).
    """
    action_cols = [c for c in table.columns if c.startswith("action_input_")]
    if not action_cols:
        return len(table)
    is_dropped = pd.concat([table[c] == MergeAction.DROP.value for c in action_cols], axis=1).any(
        axis=1
    )
    return int((~is_dropped).sum())


def compute_intersection_size(table: pd.DataFrame) -> int:
    """Variants where strand resolution was attempted at least once.

    Gate (b)'s denominator per HLD §Exit-1 validation gates: variants where
    at least one input had a matched (non-FILL_MISSING) variant.
    """
    action_cols = [c for c in table.columns if c.startswith("action_input_")]
    if not action_cols:
        return 0
    has_match = pd.concat(
        [table[c] != MergeAction.FILL_MISSING.value for c in action_cols], axis=1
    ).any(axis=1)
    return int(has_match.sum())


def warn_extras_threshold(
    n_extras: int,
    n_canonical: int,
    threshold_fraction: float,
    quiet: bool,
) -> None:
    """Emit stderr warning when extras > threshold per HLD §Extra-variant handling."""
    if quiet or n_canonical == 0:
        return
    if n_extras / n_canonical > threshold_fraction:
        print(
            f"WARNING: {n_extras} extra variants in input[N] absent from input[0] "
            f"(> {threshold_fraction:.0%} of canonical's {n_canonical} variants). "
            f"Input order may be reversed.",
            file=sys.stderr,
        )


def evaluate_pass1_gates(
    alignment_table: pd.DataFrame,
    summary: AlignmentSummary,
    policy: MergePolicy,
    is_validate_mode: bool,
) -> None:
    """Evaluate pass-1-checkable Exit-1 gates per HLD §Exit-1 validation gates.

    Called from inside `merge.merge_inputs` (between pass 1 and pass 2) and
    from `validate`'s orchestrator after `build_alignment_table`.

    Gates:
      (a) extras_count > policy.extras_warn_threshold * n_canonical
      (b) n_dropped_ambiguous_strand > policy.validate_strand_fail_pct/100 *
          intersection_size
      (d) is_validate_mode AND any --on-* error policy would have triggered

    Gate (c) target call rate is genotype-dependent and lives in
    `merge.merge_inputs` (post-pass-2).

    Raises:
        ValidationError: any gate fires.
    """
    n_canonical = len(alignment_table)
    if n_canonical == 0:
        return

    # Gate (a): extras above the --on-extra warn threshold.
    # The threshold is named "warn threshold" because it's the --on-extra warn
    # threshold; gate (a) only applies when the policy is `warn`. Under
    # --on-extra drop the user has explicitly opted into "extras are normal";
    # under --on-extra error the violation already raised in build_alignment_table.
    if policy.on_extra == "warn" and (
        summary.n_extras_dropped > policy.extras_warn_threshold * n_canonical
    ):
        raise ValidationError(
            f"gate (a): extras above warn threshold "
            f"({summary.n_extras_dropped} extras > "
            f"{policy.extras_warn_threshold:.0%} of {n_canonical} canonical variants). "
            f"Input order may be reversed; use --on-extra drop if extras are intentional."
        )

    # Gate (b): ambiguous-strand drops above intersection threshold.
    n_ambig_drops = summary.n_dropped_by_reason.get(DropReason.AMBIGUOUS_STRAND, 0)
    intersection_size = compute_intersection_size(alignment_table)
    if intersection_size > 0:
        ambig_pct = (n_ambig_drops / intersection_size) * 100.0
        if ambig_pct > policy.validate_strand_fail_pct:
            raise ValidationError(
                f"gate (b): ambiguous-strand drops above intersection threshold "
                f"({n_ambig_drops}/{intersection_size} = {ambig_pct:.1f}% > "
                f"{policy.validate_strand_fail_pct:.1f}%)"
            )

    # Gate (d): validate-mode soft-policy violations.
    if is_validate_mode:
        triggers = summary.policy_error_triggers
        msg_parts: list[str] = []
        if policy.on_mismatch == "error" and triggers.get("on_mismatch_count", 0) > 0:
            msg_parts.append(f"{triggers['on_mismatch_count']} on-mismatch trigger(s)")
        if policy.on_missing == "error" and triggers.get("on_missing_count", 0) > 0:
            msg_parts.append(f"{triggers['on_missing_count']} on-missing trigger(s)")
        if policy.on_extra == "error" and triggers.get("on_extra_count", 0) > 0:
            msg_parts.append(f"{triggers['on_extra_count']} on-extra trigger(s)")
        if msg_parts:
            raise ValidationError(
                f"gate (d): policy-error conditions would have fired in merge mode: "
                f"{'; '.join(msg_parts)}"
            )


def build_action_histogram(alignment_table: pd.DataFrame) -> dict[str, int]:
    """Build the 8-key action_histogram per LLD §2.10.

    Sums per-input action counts across all non-canonical inputs.
    All 8 keys always present (zero-valued if no variants matched) so the
    JSON schema stays stable for workflow consumers regardless of outcome.
    """
    histogram: dict[str, int] = {
        "passthrough": 0,
        "swap": 0,
        "flip": 0,
        "fill_missing": 0,
        "dropped_ambiguous_strand": 0,
        "dropped_allele_mismatch": 0,
        "pre_alignment_filter_dropped": 0,
        "drop": 0,  # residual: --on-missing drop_variant per LLD §2.10 mapping
    }

    action_cols = [c for c in alignment_table.columns if c.startswith("action_input_")]
    for action_col in action_cols:
        reason_col = action_col.replace("action_input_", "drop_reason_input_")
        actions = alignment_table[action_col]
        reasons = alignment_table[reason_col] if reason_col in alignment_table.columns else None

        histogram["passthrough"] += int((actions == MergeAction.PASSTHROUGH.value).sum())
        histogram["swap"] += int((actions == MergeAction.REF_ALT_SWAP.value).sum())
        histogram["flip"] += int(
            (
                (actions == MergeAction.STRAND_FLIP.value)
                | (actions == MergeAction.STRAND_FLIP_AND_SWAP.value)
            ).sum()
        )
        histogram["fill_missing"] += int((actions == MergeAction.FILL_MISSING.value).sum())

        if reasons is not None:
            histogram["dropped_ambiguous_strand"] += int(
                (reasons == DropReason.AMBIGUOUS_STRAND.value).sum()
            )
            histogram["dropped_allele_mismatch"] += int(
                (reasons == DropReason.ALLELE_MISMATCH.value).sum()
            )
            histogram["pre_alignment_filter_dropped"] += int(
                (
                    (reasons == DropReason.NON_BIALLELIC.value)
                    | (reasons == DropReason.NON_SNP.value)
                    | (reasons == DropReason.PRE_ALIGNMENT_OTHER.value)
                ).sum()
            )
            histogram["drop"] += int((reasons == DropReason.ON_MISSING_DROP_VARIANT.value).sum())

    return histogram
