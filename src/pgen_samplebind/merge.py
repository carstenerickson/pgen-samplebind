"""Pass 1 + Exit-1 gate evaluation + pass 2: PgenReader/Writer orchestration,
genotype recoding, output writing.

Per LLD §3.10 and HLD §Module orchestration. The HLD-pinned signature
`merge_inputs(inputs, out_pgen_path, out_pvar_path, options) -> MergeCounters`
encapsulates pass 1 + gates + pass 2 in a single call; the orchestrator
finalizes `.psam` from the returned MergeCounters.

Day 3 scope: end-to-end PFILE merge for the simple cases. Deferred:
- ctx.report_tsv_path streaming → Day 4
- ctx.collect_variant_rows → Day 4
- gate (c) target call rate is checked post-pass-2 (HLD §Exit-1 validation
  gates (c)). Detected by descriptor.is_target; uses canonical variant
  count as denominator and the per_sample_het called_count as numerator.
- Output cleanup wrapper → Day 6 (orchestrator-side)
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm import tqdm

from . import pseudohaploid, pvar, reporting
from .alignment import (
    build_action_histogram,
    build_action_histogram_per_chrom,
    build_alignment_table,
    compute_intersection_size,
    evaluate_pass1_gates,
    warn_extras_threshold,
)
from .errors import IOFailure, ValidationError
from .types import (
    AlignmentSummary,
    InputDescriptor,
    MergeAction,
    MergeContext,
    MergeCounters,
    SampleIdentityPlan,
)


def _check_target_call_rate(
    target_idx: int,
    sample_plan: SampleIdentityPlan,
    per_sample_het: list[tuple[str, int, int]],
    n_canonical_variants: int,
    min_call_rate: float,
) -> None:
    """Gate (c) per HLD §Exit-1 validation gates (c) and §Target mode.

    For each kept target sample, computes call_rate = called_count /
    n_canonical_variants and raises ValidationError if any sample falls
    below min_call_rate. The denominator is pinned to canonical (panel)
    variant count to prevent a tiny-but-fully-called target from spuriously
    passing the gate (HLD §Target mode).
    """
    if n_canonical_variants <= 0:
        return  # nothing to check
    target_keep_mask = sample_plan.per_input_keep_mask[target_idx]
    target_output_indices = sample_plan.per_input_output_indices[target_idx]

    failed: list[tuple[str, float]] = []
    for i, kept in enumerate(target_keep_mask.tolist()):
        if not kept:
            continue
        out_pos = int(target_output_indices[i])
        iid, _het, called = per_sample_het[out_pos]
        rate = called / n_canonical_variants
        if rate < min_call_rate:
            failed.append((iid, rate))

    if failed:
        head = ", ".join(f"{iid}={rate:.1%}" for iid, rate in failed[:5])
        more = "" if len(failed) <= 5 else f" (and {len(failed) - 5} more)"
        raise ValidationError(
            f"gate (c): target call rate below {min_call_rate:.0%} threshold "
            f"for {len(failed)} sample(s): {head}{more}. The target's missingness "
            f"is too high for reliable downstream analysis (e.g., qpAdm). "
            f"Override with --target-min-call-rate FLOAT to lower the threshold "
            f"if you accept the risk."
        )


def _kept_variants_table(alignment_table: pd.DataFrame) -> pd.DataFrame:
    """Return the subset of alignment_table where no per-input action is DROP."""
    action_cols = [c for c in alignment_table.columns if c.startswith("action_input_")]
    if not action_cols:
        return alignment_table.reset_index(drop=True)
    is_dropped = pd.concat(
        [alignment_table[c] == MergeAction.DROP.value for c in action_cols], axis=1
    ).any(axis=1)
    return alignment_table.loc[~is_dropped].reset_index(drop=True)


def _iter_blocks_chrom_aware(table: pd.DataFrame, block_size: int) -> Iterator[tuple[int, int]]:
    """Yield (start, end) row indices into `table` such that each block:
      - Has at most `block_size` rows.
      - Does not span a chromosome boundary.

    Per LLD §3.10 chromosome-boundary pin (lets pseudohaploid.update_block
    receive a single chrom value per block).
    """
    n = len(table)
    if n == 0:
        return
    chroms = table["chrom"].to_numpy()
    start = 0
    while start < n:
        end_max = min(start + block_size, n)
        first_chrom = chroms[start]
        end = start + 1
        while end < end_max and chroms[end] == first_chrom:
            end += 1
        yield start, end
        start = end


def _swap_genotypes_in_place(buf_row: np.ndarray[Any, Any]) -> None:
    """Apply x → 2-x to a single variant row, preserving -9 missing.

    REF_ALT_SWAP and STRAND_FLIP_AND_SWAP both reduce to this single op
    on hardcalls (LLD §2.1 action-collapse pin: strand-flip alone is
    metadata-only on biallelic SNPs).
    """
    not_missing = buf_row != -9
    buf_row[not_missing] = 2 - buf_row[not_missing]


def _write_pvar_tsv(kept_table: pd.DataFrame, out_pvar_path: Path) -> None:
    """Write the kept-variants .pvar (TSV, plink2 format with #CHROM header).

    Includes the `CM` column so plink2's downstream `--make-bed` / AT2's
    `extract_f2 blgsize` etc. see the original genetic-position values rather
    than zeros (which would collapse Morgan-spaced jackknife blocks).
    """
    cols = ["chrom", "pos", "variant_id", "ref", "alt"]
    rename = {
        "chrom": "#CHROM",
        "pos": "POS",
        "variant_id": "ID",
        "ref": "REF",
        "alt": "ALT",
    }
    if "cm" in kept_table.columns:
        cols.append("cm")
        rename["cm"] = "CM"
    out = kept_table[cols].rename(columns=rename)
    try:
        out.to_csv(out_pvar_path, sep="\t", index=False, lineterminator="\n")
    except OSError as e:
        raise IOFailure(f"cannot write {out_pvar_path}: {e}") from e


def _read_canonical_block(
    reader: object,  # pgenlib.PgenReader
    block_alignment: pd.DataFrame,
    n_samples_input: int,
    keep_mask: np.ndarray[Any, Any],
    out_indices: np.ndarray[Any, Any],
    output_buf: np.ndarray[Any, Any],
) -> None:
    """Read canonical (input[0]) for this block. Canonical is always passthrough
    by construction (it IS the reference); no action-column lookup needed.
    """
    canonical_idxs = block_alignment["canonical_idx"].to_numpy().astype(np.uint32)
    canonical_buf = np.empty((len(canonical_idxs), n_samples_input), dtype=np.int8)
    reader.read_list(canonical_idxs, canonical_buf, sample_maj=0)  # type: ignore[attr-defined]

    kept_samples = canonical_buf[:, keep_mask]
    kept_out_indices = out_indices[keep_mask]
    output_buf[:, kept_out_indices] = kept_samples


def _read_block_from_input(
    reader: object,  # pgenlib.PgenReader; typed as object to keep mypy happy w/o stubs
    block_alignment: pd.DataFrame,
    input_idx: int,
    n_samples_input: int,
) -> tuple[np.ndarray[Any, Any], np.ndarray[Any, Any]]:
    """Read this input's block of variants. Returns (read_buf, needs_read_mask).

    read_buf has shape (n_to_read, n_samples_input) — only rows that needed
    actual reading. needs_read_mask is a length-block_rows bool indicating
    which rows in the block were read (False = FILL_MISSING).
    """
    action_col = f"action_input_{input_idx}"
    idx_col = f"idx_input_{input_idx}"
    actions = block_alignment[action_col].to_numpy(dtype=object)
    needs_read = actions != MergeAction.FILL_MISSING.value
    if not needs_read.any():
        return np.empty((0, n_samples_input), dtype=np.int8), needs_read

    raw_idxs = block_alignment.loc[needs_read, idx_col]
    # Int64 nullable → uint32 numpy. needs_read excludes FILL_MISSING so all are non-null.
    read_idxs = raw_idxs.astype("int64").to_numpy().astype(np.uint32)

    read_buf = np.empty((len(read_idxs), n_samples_input), dtype=np.int8)
    reader.read_list(read_idxs, read_buf, sample_maj=0)  # type: ignore[attr-defined]
    return read_buf, needs_read


def _apply_actions_and_place(
    read_buf: np.ndarray[Any, Any],
    needs_read: np.ndarray[Any, Any],
    block_alignment: pd.DataFrame,
    input_idx: int,
    keep_mask: np.ndarray[Any, Any],
    out_indices: np.ndarray[Any, Any],
    output_buf: np.ndarray[Any, Any],
) -> None:
    """Apply per-row REF/ALT-swap recoding and place samples into output_buf.

    output_buf shape: (block_rows, n_output_samples). Pre-filled with -9.
    Updated in place for this input's contribution.
    """
    actions = block_alignment[f"action_input_{input_idx}"].to_numpy(dtype=object)
    actions_subset = actions[needs_read]

    # REF_ALT_SWAP and STRAND_FLIP_AND_SWAP both swap; STRAND_FLIP and PASSTHROUGH
    # are identity (LLD §2.1 action-collapse).
    swap_mask = (actions_subset == MergeAction.REF_ALT_SWAP.value) | (
        actions_subset == MergeAction.STRAND_FLIP_AND_SWAP.value
    )
    if swap_mask.any():
        for i in np.where(swap_mask)[0]:
            _swap_genotypes_in_place(read_buf[i])

    # Apply sample plan: keep mask + place into output positions
    kept_samples = read_buf[:, keep_mask]
    kept_out_indices = out_indices[keep_mask]
    block_rows_with_read = np.where(needs_read)[0]

    for write_idx, source_row in enumerate(block_rows_with_read):
        output_buf[source_row, kept_out_indices] = kept_samples[write_idx]


def merge_inputs(
    inputs: list[InputDescriptor],
    out_pgen_path: Path,
    out_pvar_path: Path,
    ctx: MergeContext,
) -> MergeCounters:
    """Pass 1 + Exit-1 gate evaluation + pass 2 per LLD §3.10.

    Internally:
      1. Read pvars; build alignment_table via alignment.build_alignment_table.
      2. evaluate_pass1_gates (a)/(b) — raises ValidationError on either.
      3. Filter to kept variants; open PgenWriter with exact variant_ct +
         sample_ct from sample_plan.
      4. Stream pass 2 in chrom-aware blocks: read each input via read_list,
         apply per-(variant, input) action vectorized, place into output via
         per_input_output_indices, write via append_biallelic_batch.
      5. Update per_sample_het counters per block via pseudohaploid.update_block.
      6. Write out_pvar_path (TSV) after the .pgen completes.

    Returns MergeCounters per HLD §Module orchestration.

    Raises:
        ValidationError: gate (a) or (b) fires.
        IOFailure: write failure or PgenWriter close mismatch.
    """
    import pgenlib

    # ----- Pass 1 -----
    pvars = [pvar.read_pvar(d.pvar_path) for d in inputs]
    canonical_pvar = pvars[0]
    other_pvars = pvars[1:]
    summary = AlignmentSummary()
    alignment_table = build_alignment_table(
        canonical_pvar, other_pvars, ctx.policy, summary, soften_policy_errors=False
    )

    # ----- Stderr warning + gates (a)/(b) -----
    # Warning fires before gate (a) so the user sees both signals (the
    # informational "extras count is high" line AND the structured
    # ValidationError that the gate then raises). Warning is suppressed
    # under --on-extra drop (user has opted into "extras are intentional").
    if ctx.policy.on_extra == "warn":
        warn_extras_threshold(
            summary.n_extras_dropped,
            len(canonical_pvar),
            ctx.policy.extras_warn_threshold,
            quiet=False,
        )
    evaluate_pass1_gates(alignment_table, summary, ctx.policy, is_validate_mode=False)

    # ----- Pass 2 setup -----
    kept_table = _kept_variants_table(alignment_table)
    n_kept = len(kept_table)
    n_output_samples = len(ctx.sample_plan.output_iids)
    n_samples_per_input = [d.n_samples for d in inputs]

    het_counts, called_counts = pseudohaploid.init_counters(n_output_samples)

    writer = pgenlib.PgenWriter(str(out_pgen_path).encode(), n_output_samples, n_kept)
    readers: list[object] = []
    try:
        for desc, n_samples in zip(inputs, n_samples_per_input, strict=True):
            readers.append(
                pgenlib.PgenReader(str(desc.pgen_path).encode(), raw_sample_ct=n_samples)
            )

        # ----- Pass 2: chrom-aware block loop -----
        # Progress bar wraps the block iterator when ctx.show_progress is set
        # (v0.2). Unit is variants (each block contributes block_size_actual);
        # tqdm renders rate as variants/s in the bar suffix. Off in piped or
        # --quiet contexts so workflow-manager stderr stays clean.
        block_iter: Iterator[tuple[int, int]] = _iter_blocks_chrom_aware(
            kept_table, ctx.policy.block_size
        )
        progress_bar: tqdm[Any] | None = None
        if ctx.show_progress:
            progress_bar = tqdm(
                total=n_kept,
                unit=" variants",
                desc="Pass 2: streaming genotypes",
                leave=False,
            )

        for block_start, block_end in block_iter:
            block_alignment = kept_table.iloc[block_start:block_end]
            block_size_actual = len(block_alignment)
            chrom_int = int(block_alignment["chrom"].iloc[0])

            output_buf = np.full((block_size_actual, n_output_samples), -9, dtype=np.int8)

            # Canonical (input[0]): always passthrough; no action column lookup.
            _read_canonical_block(
                readers[0],
                block_alignment,
                n_samples_per_input[0],
                ctx.sample_plan.per_input_keep_mask[0],
                ctx.sample_plan.per_input_output_indices[0],
                output_buf,
            )

            # Non-canonical inputs (1..N-1): action_input_<i> driven.
            for i, reader in enumerate(readers[1:], start=1):
                # i matches the per-input column suffix: action_input_<i>.
                read_buf, needs_read = _read_block_from_input(
                    reader, block_alignment, i, n_samples_per_input[i]
                )
                if read_buf.size == 0:
                    continue
                _apply_actions_and_place(
                    read_buf,
                    needs_read,
                    block_alignment,
                    i,
                    ctx.sample_plan.per_input_keep_mask[i],
                    ctx.sample_plan.per_input_output_indices[i],
                    output_buf,
                )

            pseudohaploid.update_block(output_buf, chrom_int, het_counts, called_counts)
            writer.append_biallelic_batch(output_buf)
            if progress_bar is not None:
                progress_bar.update(block_size_actual)

        if progress_bar is not None:
            progress_bar.close()

        try:
            writer.close()
        except RuntimeError as e:
            raise IOFailure(
                f"PgenWriter.close() raised on {out_pgen_path}: {e}. This indicates a "
                f"variant_ct mismatch — likely a bug in alignment.count_kept_variants."
            ) from e
    finally:
        for reader in readers:
            close = getattr(reader, "close", None)
            if callable(close):
                close()

    # ----- Pass 2 outputs -----
    _write_pvar_tsv(kept_table, out_pvar_path)

    # ----- Reports (per LLD §3.10 step 7-8) -----
    if ctx.report_tsv_path is not None:
        reporting.write_report_tsv(alignment_table, len(inputs), ctx.report_tsv_path)

    variant_rows = (
        reporting.build_variant_rows_from_alignment(alignment_table, len(inputs))
        if ctx.collect_variant_rows
        else None
    )

    per_sample_het = [
        (iid, int(het_counts[i]), int(called_counts[i]))
        for i, iid in enumerate(ctx.sample_plan.output_iids)
    ]

    # ----- Gate (c): target call rate -----
    # Per HLD §Exit-1 validation gates (c): each target's call rate (non-missing
    # genotypes / canonical variant count) must be >= policy.target_min_call_rate.
    # Genotype-dependent → checked here post-pass-2, not in evaluate_pass1_gates.
    # Per LLD §3.10: if it fires, the .pgen/.pvar exist on disk; the orchestrator's
    # try/except (LLD §4.1 fix #6) unlinks the partial triplet. Multi-target
    # mode evaluates the gate independently per target — strict semantics: any
    # failing target blocks the whole merge.
    for target_idx, desc in enumerate(inputs):
        if not desc.is_target:
            continue
        _check_target_call_rate(
            target_idx=target_idx,
            sample_plan=ctx.sample_plan,
            per_sample_het=per_sample_het,
            n_canonical_variants=len(canonical_pvar),
            min_call_rate=ctx.policy.target_min_call_rate,
        )

    return MergeCounters(
        action_histogram=build_action_histogram(alignment_table),
        action_histogram_per_chrom=build_action_histogram_per_chrom(alignment_table),
        intersection_size=compute_intersection_size(alignment_table),
        extras_count=summary.n_extras_dropped,
        per_sample_het=per_sample_het,
        n_output_samples=n_output_samples,
        n_output_variants=n_kept,
        variant_rows=variant_rows,
    )
