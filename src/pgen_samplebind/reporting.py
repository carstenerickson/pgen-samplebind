"""Stdout summary, TSV per-variant report, JSON report.

Per  The TSV report is bulk-written from the alignment table
post-pass-2 (memory bounded by alignment_table size which we already
hold for the merge); the JSON report defaults to summary-only
(workflow-friendly, ~few KB) with `--report-json-include-rows` opting
into the per-variant array (~100 B/row x n_variants x n_inputs;
warns at >100 MB).
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from . import pseudohaploid as _pseudohaploid_mod
from .types import (
    InputDescriptor,
    MergeCounters,
    MergePolicy,
    PseudohaploidStatus,
    ReportRow,
)

# Conservative bytes-per-JSON-row estimate for the size warning.
# The buffered Python object footprint is higher than this; the warning
# fires before the user actually feels memory pressure (see  caveat).
_JSON_BYTES_PER_ROW = 100
_JSON_INCLUDE_ROWS_WARN_BYTES = 100 * 1024 * 1024  # 100 MB


def format_stdout_summary(
    counters: MergeCounters,
    descriptors: list[InputDescriptor],
    psam_dfs: list[pd.DataFrame],
    merged_psam: pd.DataFrame,
    output_paths: dict[str, Path],
    elapsed_s: float,
    quiet: bool,
) -> str:
    """Render the  summary block. Empty string if quiet=True."""
    if quiet:
        return ""

    lines: list[str] = []

    # Inputs block
    lines.append(f"Read {len(descriptors)} inputs:")
    for i, desc in enumerate(descriptors):
        canonical_marker = "  (canonical)" if i == 0 else ""
        lines.append(f"  [{i}] {desc.path}: {desc.n_samples:,} samples{canonical_marker}")
    lines.append("")

    # Variant alignment block (per HLD example)
    canonical_label = str(descriptors[0].path) if descriptors else "?"
    lines.append(f"Variant alignment (canonical = {canonical_label}):")
    hist = counters.action_histogram
    lines.append(f"  Passthrough:       {hist.get('passthrough', 0):>12,}")
    lines.append(f"  REF/ALT swapped:   {hist.get('swap', 0):>12,} (recoded)")
    lines.append(f"  Strand-flipped:    {hist.get('flip', 0):>12,} (complemented)")
    lines.append(
        f"  Filled missing:    {hist.get('fill_missing', 0):>12,} "
        f"(variants in [0] absent from [N>0])"
    )
    n_dropped_total = (
        hist.get("dropped_ambiguous_strand", 0)
        + hist.get("dropped_allele_mismatch", 0)
        + hist.get("pre_alignment_filter_dropped", 0)
        + hist.get("drop", 0)
    )
    lines.append(f"  Dropped:           {n_dropped_total:>12,}")
    if hist.get("dropped_ambiguous_strand", 0):
        lines.append(f"    Ambiguous strand:    {hist['dropped_ambiguous_strand']:,}")
    if hist.get("dropped_allele_mismatch", 0):
        lines.append(f"    Allele mismatch:     {hist['dropped_allele_mismatch']:,}")
    if hist.get("pre_alignment_filter_dropped", 0):
        lines.append(f"    Pre-filter dropped:  {hist['pre_alignment_filter_dropped']:,}")
    if hist.get("drop", 0):
        lines.append(f"    --on-missing drop:   {hist['drop']:,}")
    lines.append("")

    # Per-population breakdown
    lines.append("Per-population sample counts:")
    for i, (desc, psam_df) in enumerate(zip(descriptors, psam_dfs, strict=True)):
        if "POP" not in psam_df.columns:
            continue
        pop_counts = Counter(psam_df["POP"].tolist())
        largest = pop_counts.most_common(1)[0] if pop_counts else ("?", 0)
        lines.append(
            f"  [{i}] {desc.path}: {len(pop_counts)} populations, "
            f"{len(psam_df):,} samples (largest: {largest[0]} n={largest[1]:,})"
        )
    if "POP" in merged_psam.columns:
        out_pop_counts = Counter(merged_psam["POP"].tolist())
        lines.append(f"  Output: {len(out_pop_counts)} populations, {len(merged_psam):,} samples")
    lines.append("")

    # Pseudohaploid mix
    lines.append("Pseudohaploid mix (output):")
    pseudo: dict[str, int] = {"1": 0, "0": 0, "U": 0}
    for _, h, c in counters.per_sample_het:
        status = _pseudohaploid_mod.classify(h, c)
        pseudo[status.value] += 1
    lines.append(
        f"  {pseudo['1']:,} pseudohaploid, {pseudo['0']:,} diploid, {pseudo['U']:,} unknown"
    )
    lines.append("")

    # Output files
    lines.append("Wrote:")
    for kind, p in output_paths.items():
        size = p.stat().st_size if p.exists() else 0
        lines.append(f"  {p} ({_humanize_bytes(size)})  [{kind}]")
    lines.append("")

    lines.append(
        f"Done: {counters.n_output_samples:,} samples x "
        f"{counters.n_output_variants:,} variants in {elapsed_s:.2f}s."
    )
    return "\n".join(lines) + "\n"


def _humanize_bytes(n: int) -> str:
    """Render byte count in B/KB/MB/GB."""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}" if unit != "B" else f"{n:,} B"
        n //= 1024
    return f"{n} TB"


def write_report_tsv(
    alignment_table: pd.DataFrame,
    n_inputs: int,
    path: Path,
) -> None:
    """Per-variant action TSV.

    Header: variant_id\\tchr\\tpos\\tinput_index\\taction\\treason
    One row per (variant, non-canonical input_index). Canonical (input_index=0)
    rows are implicit PASSTHROUGH and skipped pin (saves
    O(V) rows on file size at no information loss).

    Bulk-write via pandas.to_csv per non-canonical input. Memory cost is
    bounded by alignment_table which we already hold; no separate buffer.
    """
    if n_inputs <= 1:
        # No non-canonical inputs to report.
        path.touch()
        with open(path, "w", encoding="utf-8") as f:
            f.write("variant_id\tchr\tpos\tinput_index\taction\treason\n")
        return

    # Build a long-form DataFrame: for each non-canonical input, expand into rows.
    # Done per-input to keep memory bounded by O(V) per chunk rather than O(V*N).
    with open(path, "w", encoding="utf-8") as f:
        f.write("variant_id\tchr\tpos\tinput_index\taction\treason\n")
        for i in range(1, n_inputs):
            action_col = f"action_input_{i}"
            reason_col = f"drop_reason_input_{i}"
            chunk = pd.DataFrame(
                {
                    "variant_id": alignment_table["variant_id"],
                    "chr": alignment_table["chrom"],
                    "pos": alignment_table["pos"],
                    "input_index": i,
                    "action": alignment_table[action_col].astype(str),
                    "reason": alignment_table[reason_col].astype(str).fillna("").replace("nan", ""),
                }
            )
            chunk.to_csv(f, sep="\t", index=False, header=False, lineterminator="\n")


def _serialize_variant_rows(rows: list[ReportRow]) -> list[dict[str, Any]]:
    """Render `ReportRow` objects as JSON-friendly dicts.

    Same field set as the TSV row format. Caller has already handled the
    canonical-row skip
    """
    return [
        {
            "variant_id": r.variant_id,
            "chr": int(r.chrom),
            "pos": int(r.pos),
            "input_index": int(r.input_index),
            "action": r.action.value,
            "reason": r.reason,
        }
        for r in rows
    ]


def build_variant_rows_from_alignment(
    alignment_table: pd.DataFrame, n_inputs: int
) -> list[ReportRow]:
    """Build `ReportRow` objects from the alignment table for the JSON
    include-rows path. Skips canonical (input_index=0)

    Used by `merge.merge_inputs` when ctx.collect_variant_rows is True,
    and by `validate`'s orchestrator (which builds rows from the alignment
    table directly since validate has no pass-2).
    """
    from .types import MergeAction

    rows: list[ReportRow] = []
    if n_inputs <= 1:
        return rows
    for i in range(1, n_inputs):
        action_col = f"action_input_{i}"
        reason_col = f"drop_reason_input_{i}"
        for vid, chrom, pos, action_str, reason_str in zip(
            alignment_table["variant_id"],
            alignment_table["chrom"],
            alignment_table["pos"],
            alignment_table[action_col].astype(str),
            alignment_table[reason_col].astype(str),
            strict=True,
        ):
            reason = "" if reason_str in ("nan", "None") else str(reason_str)
            rows.append(
                ReportRow(
                    variant_id=str(vid),
                    chrom=int(chrom),
                    pos=int(pos),
                    input_index=i,
                    action=MergeAction(action_str),
                    reason=reason,
                )
            )
    return rows


def _maybe_warn_json_size(n_output_variants: int, n_inputs: int, quiet: bool) -> None:
    """Stderr warning if predicted JSON-include-rows size > 100 MB.

    Threshold is the JSON byte size at write time; Python object footprint
    during merge is higher. Warning is conservative on purpose — fires
    before users see the memory wall.
    """
    if quiet:
        return
    n_rows = n_output_variants * max(n_inputs - 1, 0)
    predicted_bytes = n_rows * _JSON_BYTES_PER_ROW
    if predicted_bytes > _JSON_INCLUDE_ROWS_WARN_BYTES:
        print(
            f"WARNING: --report-json-include-rows predicted size "
            f"{_humanize_bytes(predicted_bytes)} (> 100 MB). For streaming "
            f"per-variant data at scale, prefer --report TSV (constant memory).",
            file=sys.stderr,
        )


def write_report_json(
    counters: MergeCounters,
    descriptors: list[InputDescriptor],
    psam_dfs: list[pd.DataFrame],
    output_paths: dict[str, Path],
    policy: MergePolicy,
    path: Path,
    include_rows: bool,
    quiet: bool,
    mode: str = "merge",
    tool_version: str = "",
) -> None:
    """Run-level JSON report.

    Default (include_rows=False): summary-only (~few KB). Workflow-friendly.
    include_rows=True: adds `variants` array from `counters.variant_rows`;
    warns at >100 MB predicted size.

    `mode` is "merge" or "validate"; validate mode lacks pass-2 fields
    (n_output_samples derived from psams; per_sample_het is empty).
    """
    n_inputs = len(descriptors)

    if include_rows:
        _maybe_warn_json_size(counters.n_output_variants, n_inputs, quiet)

    payload: dict[str, Any] = {
        "tool": "pgen-samplebind",
        "tool_version": tool_version,
        "command": mode,
        "inputs": [
            {
                "path": str(d.path),
                "format": d.fmt.value,
                "n_samples": int(d.n_samples),
                "n_variants": int(d.n_variants),
            }
            for d in descriptors
        ],
        "policy": {
            "on_mismatch": policy.on_mismatch,
            "on_missing": policy.on_missing,
            "on_extra": policy.on_extra,
            "on_strand": policy.on_strand,
            "on_collision": policy.on_collision,
            "trust_strand": policy.trust_strand,
            "variant_key": policy.variant_key,
            "validate_strand_fail_pct": policy.validate_strand_fail_pct,
            "extras_warn_threshold": policy.extras_warn_threshold,
        },
        "alignment": {
            "action_histogram": dict(counters.action_histogram),
            # Per-chromosome 9-key breakdown. JSON object keys must be strings,
            # so stringify the chrom int. Consumers read with `int(k)`. Added
            # in v0.2 as a diagnostic for chr-specific drop concentrations
            # (HLA strand artifacts, hg19/hg38 build mismatches, etc.).
            "action_histogram_per_chrom": {
                str(chrom): dict(hist)
                for chrom, hist in counters.action_histogram_per_chrom.items()
            },
            "intersection_size": int(counters.intersection_size),
            "extras_count": int(counters.extras_count),
        },
        "output": {
            "n_samples": int(counters.n_output_samples),
            "n_variants": int(counters.n_output_variants),
            "paths": {k: str(v) for k, v in output_paths.items()},
        },
        "pseudohaploid": _pseudohaploid_summary(counters),
        "per_input_populations": [
            dict(Counter(df["POP"].tolist())) if "POP" in df.columns else {} for df in psam_dfs
        ],
    }

    if include_rows and counters.variant_rows is not None:
        payload["variants"] = _serialize_variant_rows(counters.variant_rows)

    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
            f.write("\n")
    except OSError as e:
        from .errors import IOFailure

        raise IOFailure(f"cannot write {path}: {e}") from e


def _pseudohaploid_summary(counters: MergeCounters) -> dict[str, int]:
    """Compute {0,1,U} sample counts from counters.per_sample_het."""
    out = {
        PseudohaploidStatus.DIPLOID.value: 0,
        PseudohaploidStatus.PSEUDOHAPLOID.value: 0,
        PseudohaploidStatus.UNKNOWN.value: 0,
    }
    for _, h, c in counters.per_sample_het:
        status = _pseudohaploid_mod.classify(h, c)
        out[status.value] += 1
    return out
