"""`merge` subcommand orchestrator. Sequence in LLD §4.1.

Day 3-4 scope: end-to-end PFILE merge with reports + HLD-style stdout summary.
Deferred to later days:
- check_plink2_available (only needed for EIGENSTRAT/BFILE) — Day 6
- --relabel-from (psam relabel) — Day 9
- Output cleanup wrapper (try/except → unlink triplet) — Day 6
- Concurrency lock (output_lock) — Day 10
"""

from __future__ import annotations

import sys
import time
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path

import numpy as np

from .. import __version__, psam, pseudohaploid, reporting
from ..formats import prepared_input
from ..merge import merge_inputs
from ..pvar import check_max_alleles, count_raw_variants
from ..types import MergeContext, MergePolicy


def run_merge(
    input_paths: tuple[Path, ...],
    target_path: Path | None,
    output_prefix: Path,
    policy: MergePolicy,
    report_path: Path | None,
    report_json_path: Path | None,
    quiet: bool,
) -> None:
    """Merge subcommand orchestrator. Per LLD §4.1.

    Day 3-4: target_path is rejected (--target mode is Day 8).
    """
    if target_path is not None:
        raise NotImplementedError(
            "--target mode is deferred to project Day 8. For Day 3-4 use positional "
            "INPUT arguments only."
        )

    started = time.perf_counter()

    out_pgen_path = Path(str(output_prefix) + ".pgen")
    out_pvar_path = Path(str(output_prefix) + ".pvar")
    out_psam_path = Path(str(output_prefix) + ".psam")
    output_paths = {"pgen": out_pgen_path, "pvar": out_pvar_path, "psam": out_psam_path}

    with ExitStack() as stack:
        # Step 1: format detection per input via context manager
        descriptors = [
            stack.enter_context(
                prepared_input(p, is_target=False, include_chrom=policy.include_chrom)
            )
            for p in input_paths
        ]

        # Step 5: per-input multi-allelic startup check
        for desc in descriptors:
            check_max_alleles(desc.pgen_path)

        # Step 6: read psams; detect population column; rename → POP; add FID = POP
        psam_dfs = []
        for desc in descriptors:
            df = psam.read_psam(desc.psam_path)
            pop_col = psam.detect_population_column(df, policy.population_column)
            df = psam.rename_to_pop(df, pop_col)
            df = psam.add_fid_from_pop(df)
            psam_dfs.append(df)

        # Populate descriptor n_samples (from psam) and n_variants (cheap raw line
        # count from .pvar) — both used by reporting; merge_inputs reads pvars
        # internally and doesn't depend on n_variants here.
        descriptors = [
            replace(d, n_samples=len(df), n_variants=count_raw_variants(d.pvar_path))
            for d, df in zip(descriptors, psam_dfs, strict=True)
        ]

        # Step 10: resolve sample identity (collision policy applied)
        sample_plan = psam.resolve_sample_identity(psam_dfs, policy, target_idx=None)

        # Step 11: build MergeContext
        ctx = MergeContext(
            policy=policy,
            sample_plan=sample_plan,
            report_tsv_path=report_path,
            collect_variant_rows=(policy.report_json_include_rows and report_json_path is not None),
        )

        # Step 12: merge_inputs (pass 1 + gates + pass 2; writes .pgen + .pvar +
        # report TSV if ctx.report_tsv_path; populates counters.variant_rows
        # if ctx.collect_variant_rows).
        counters = merge_inputs(descriptors, out_pgen_path, out_pvar_path, ctx)

        # Step 13: psam finalization
        merged_psam = psam.merge_psams(psam_dfs, sample_plan)

        # Step 14: classify pseudohaploid from counters; assign PSEUDOHAPLOID column.
        # Row-order invariant assertion (LLD §4.1 fix #1):
        if __debug__:
            for i, (iid_in_counters, _, _) in enumerate(counters.per_sample_het):
                assert (
                    iid_in_counters == sample_plan.output_iids[i] == merged_psam.iloc[i]["IID"]
                ), (
                    f"row-order invariant violated at i={i}: "
                    f"counters={iid_in_counters!r}, "
                    f"plan={sample_plan.output_iids[i]!r}, "
                    f"psam={merged_psam.iloc[i]['IID']!r}"
                )

        het_array = np.array([h for _, h, _ in counters.per_sample_het], dtype=np.int64)
        called_array = np.array([c for _, _, c in counters.per_sample_het], dtype=np.int64)
        statuses = pseudohaploid.classify_all(het_array, called_array)
        merged_psam["PSEUDOHAPLOID"] = [s.value for s in statuses]

        # Step 15: write .psam
        psam.write_psam(merged_psam, out_psam_path)

    elapsed = time.perf_counter() - started

    # Step 16: report-JSON (TSV is already written inside merge_inputs)
    if report_json_path is not None:
        reporting.write_report_json(
            counters=counters,
            descriptors=descriptors,
            psam_dfs=psam_dfs,
            output_paths=output_paths,
            policy=policy,
            path=report_json_path,
            include_rows=policy.report_json_include_rows,
            quiet=quiet,
            mode="merge",
            tool_version=__version__,
        )

    # Step 17: stdout summary
    if not quiet:
        sys.stdout.write(
            reporting.format_stdout_summary(
                counters=counters,
                descriptors=descriptors,
                psam_dfs=psam_dfs,
                merged_psam=merged_psam,
                output_paths=output_paths,
                elapsed_s=elapsed,
                quiet=quiet,
            )
        )
