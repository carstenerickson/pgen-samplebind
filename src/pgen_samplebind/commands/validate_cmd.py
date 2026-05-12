"""`validate` subcommand orchestrator. Sequence in LLD §4.2.

Validate runs pass 1 only — no genotype reads, no .pgen output. Same
Exit-1 gates as merge (HLD §Exit-1 validation gates), but --on-* error
policies are softened to gate (d) per HLD §Exit-1 validation gates (d).
"""

from __future__ import annotations

import sys
import time
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path

from .. import __version__, psam, pvar, reporting
from ..alignment import (
    build_action_histogram,
    build_action_histogram_per_chrom,
    build_alignment_table,
    compute_intersection_size,
    count_kept_variants,
    evaluate_pass1_gates,
    warn_extras_threshold,
)
from ..pvar import check_max_alleles
from ..types import AlignmentSummary, MergeCounters, MergePolicy


def run_validate(
    input_paths: tuple[Path, ...],
    policy: MergePolicy,
    report_path: Path | None,
    report_json_path: Path | None,
    quiet: bool,
    relabel_from: Path | None = None,
    relabel_input_col: str | None = None,
    relabel_output_col: str | None = None,
) -> None:
    """Validate subcommand orchestrator. Per LLD §4.2.

    Day 4: works for PFILE inputs end-to-end. BFILE/EIGENSTRAT deferred to
    Day 6 (same shell-out path as merge).
    """
    started = time.perf_counter()

    with ExitStack() as stack:
        # Step 3: format detection per input
        from ..formats import prepared_input

        descriptors = [
            stack.enter_context(
                prepared_input(p, is_target=False, include_chrom=policy.include_chrom)
            )
            for p in input_paths
        ]

        # Step 4: per-input multi-allelic check
        for desc in descriptors:
            check_max_alleles(desc.pgen_path)

        # Step 5: read psams
        psam_dfs = []
        for desc in descriptors:
            df = psam.read_psam(desc.psam_path)
            pop_col = psam.detect_population_column(df, policy.population_column)
            df = psam.rename_to_pop(df, pop_col)
            psam_dfs.append(df)

        # --relabel-from per-input (HLD §Relabeling) — same logic as merge.
        if relabel_from is not None:
            relabel_df = psam.read_relabel_tsv(relabel_from, relabel_input_col, relabel_output_col)
            source_col = "POP" if relabel_input_col is None else policy.id_column
            psam_dfs = [
                psam.apply_relabel(df, relabel_df, source_column=source_col) for df in psam_dfs
            ]

        descriptors = [
            replace(d, n_samples=len(df)) for d, df in zip(descriptors, psam_dfs, strict=True)
        ]

        # Step 7: read pvars
        pvars = [pvar.read_pvar(d.pvar_path) for d in descriptors]
        descriptors = [
            replace(d, n_variants=len(pv)) for d, pv in zip(descriptors, pvars, strict=True)
        ]

        # Step 8: validate_unique_keys on canonical (delegated)
        # Step 9: build alignment table with soften_policy_errors=True (validate mode).
        summary = AlignmentSummary()
        alignment_table = build_alignment_table(
            pvars[0],
            pvars[1:],
            policy,
            summary,
            soften_policy_errors=True,
        )

        # Step 10: extras warning
        warn_extras_threshold(
            summary.n_extras_dropped,
            len(pvars[0]),
            policy.extras_warn_threshold,
            quiet,
        )

        # Step 11: evaluate_pass1_gates with is_validate_mode=True (gates a/b/d).
        evaluate_pass1_gates(alignment_table, summary, policy, is_validate_mode=True)

        # Step 12: collision detection only (no output written)
        psam.resolve_sample_identity(psam_dfs, policy)

        # Build a partial MergeCounters for reporting. Validate has no pass 2,
        # so per_sample_het is empty; n_output_samples derived from psams.
        n_output_samples = sum(len(df) for df in psam_dfs)
        n_kept = count_kept_variants(alignment_table)
        variant_rows = (
            reporting.build_variant_rows_from_alignment(alignment_table, len(descriptors))
            if (policy.report_json_include_rows and report_json_path is not None)
            else None
        )
        counters = MergeCounters(
            action_histogram=build_action_histogram(alignment_table),
            action_histogram_per_chrom=build_action_histogram_per_chrom(alignment_table),
            intersection_size=compute_intersection_size(alignment_table),
            extras_count=summary.n_extras_dropped,
            per_sample_het=[],  # no pass 2
            n_output_samples=n_output_samples,
            n_output_variants=n_kept,
            variant_rows=variant_rows,
        )

        # Step 13: report TSV
        if report_path is not None:
            reporting.write_report_tsv(alignment_table, len(descriptors), report_path)

        # Step 14: report JSON (mode="validate")
        if report_json_path is not None:
            reporting.write_report_json(
                counters=counters,
                descriptors=descriptors,
                psam_dfs=psam_dfs,
                output_paths={},  # validate writes no PFILE output
                policy=policy,
                path=report_json_path,
                include_rows=policy.report_json_include_rows,
                quiet=quiet,
                mode="validate",
                tool_version=__version__,
            )

    elapsed = time.perf_counter() - started

    # Step 15: stdout summary (validate-mode block — slimmer than merge's)
    if not quiet:
        sys.stdout.write(
            _format_validate_stdout(
                counters=counters,
                descriptors=descriptors,
                summary=summary,
                policy=policy,
                elapsed_s=elapsed,
            )
        )


def _format_validate_stdout(
    counters: MergeCounters,
    descriptors: list,  # type: ignore[type-arg]
    summary: AlignmentSummary,
    policy: MergePolicy,
    elapsed_s: float,
) -> str:
    """Validate-mode stdout summary. Subset of merge's HLD block (no output
    PFILE; no pseudohaploid mix since pass 2 didn't run)."""
    lines: list[str] = []
    lines.append(f"Read {len(descriptors)} inputs:")
    for i, d in enumerate(descriptors):
        canonical = "  (canonical)" if i == 0 else ""
        lines.append(
            f"  [{i}] {d.path}: {d.n_samples:,} samples, {d.n_variants:,} variants{canonical}"
        )
    lines.append("")

    lines.append("Variant alignment (pass 1 only):")
    hist = counters.action_histogram
    lines.append(f"  Passthrough:       {hist.get('passthrough', 0):>12,}")
    lines.append(f"  REF/ALT swapped:   {hist.get('swap', 0):>12,}")
    lines.append(f"  Strand-flipped:    {hist.get('flip', 0):>12,}")
    lines.append(f"  Filled missing:    {hist.get('fill_missing', 0):>12,}")
    n_dropped = (
        hist.get("dropped_ambiguous_strand", 0)
        + hist.get("dropped_allele_mismatch", 0)
        + hist.get("pre_alignment_filter_dropped", 0)
        + hist.get("drop", 0)
    )
    lines.append(f"  Dropped:           {n_dropped:>12,}")
    if hist.get("dropped_ambiguous_strand", 0):
        lines.append(f"    Ambiguous strand:    {hist['dropped_ambiguous_strand']:,}")
    if hist.get("dropped_allele_mismatch", 0):
        lines.append(f"    Allele mismatch:     {hist['dropped_allele_mismatch']:,}")
    lines.append("")

    triggers = summary.policy_error_triggers
    if any(triggers.values()):
        lines.append("Policy-error triggers (would have exited 3 in merge mode):")
        for k, v in triggers.items():
            if v:
                lines.append(f"  {k}: {v:,}")
        lines.append("")

    lines.append(
        f"Validation passed: alignment OK ({n_dropped:,} dropped, "
        f"{counters.intersection_size:,} intersection, {counters.extras_count:,} extras). "
        f"Done in {elapsed_s:.2f}s."
    )
    return "\n".join(lines) + "\n"
