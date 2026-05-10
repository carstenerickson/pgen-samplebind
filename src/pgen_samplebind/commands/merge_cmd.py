"""`merge` subcommand orchestrator. Sequence in LLD §4.1.

Day 3 scope: end-to-end PFILE merge for the simple cases. Deferred to later days:
- check_plink2_available (only needed for EIGENSTRAT/BFILE) — Day 6
- --relabel-from (psam relabel) — Day 9
- Output cleanup wrapper (try/except → unlink triplet) — Day 6
- Concurrency lock (output_lock) — Day 10
- Detailed stdout summary formatting — Day 4
- Report TSV / report JSON — Day 4
"""

from __future__ import annotations

import sys
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path

from .. import psam, pseudohaploid
from ..formats import prepared_input
from ..merge import merge_inputs
from ..pvar import check_max_alleles
from ..types import InputDescriptor, MergeContext, MergeCounters, MergePolicy


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

    Day 3: target_path is rejected (--target mode is Day 8). report_path /
    report_json_path are accepted but not yet written (Day 4).
    """
    if target_path is not None:
        raise NotImplementedError(
            "--target mode is deferred to project Day 8. For Day 3 use positional "
            "INPUT arguments only."
        )

    out_pgen_path = Path(str(output_prefix) + ".pgen")
    out_pvar_path = Path(str(output_prefix) + ".pvar")
    out_psam_path = Path(str(output_prefix) + ".psam")

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

        # Populate descriptor n_samples for merge_inputs
        descriptors = [
            replace(d, n_samples=len(df)) for d, df in zip(descriptors, psam_dfs, strict=True)
        ]

        # Step 10: resolve sample identity (collision policy applied)
        sample_plan = psam.resolve_sample_identity(psam_dfs, policy, target_idx=None)

        # Step 11: build MergeContext
        ctx = MergeContext(
            policy=policy,
            sample_plan=sample_plan,
            report_tsv_path=report_path,  # Day 4 wires the actual streaming
            collect_variant_rows=policy.report_json_include_rows,
        )

        # Step 12: merge_inputs (pass 1 + gates + pass 2; writes .pgen + .pvar)
        counters = merge_inputs(descriptors, out_pgen_path, out_pvar_path, ctx)

        # Step 13: psam finalization
        merged_psam = psam.merge_psams(psam_dfs, sample_plan)

        # Step 14: classify pseudohaploid from counters; assign PSEUDOHAPLOID column
        # Row-order invariant assertion (LLD §4.1 fix #1):
        if __debug__:
            for i, (iid_in_counters, _, _) in enumerate(counters.per_sample_het):
                assert (
                    iid_in_counters == sample_plan.output_iids[i] == merged_psam.iloc[i]["IID"]
                ), (
                    f"row-order invariant violated at i={i}: "
                    f"counters={iid_in_counters!r}, plan={sample_plan.output_iids[i]!r}, "
                    f"psam={merged_psam.iloc[i]['IID']!r}"
                )

        het_array = [h for _, h, _ in counters.per_sample_het]
        called_array = [c for _, _, c in counters.per_sample_het]
        # Convert to numpy for classify_all
        import numpy as np

        statuses = pseudohaploid.classify_all(
            np.array(het_array, dtype=np.int64),
            np.array(called_array, dtype=np.int64),
        )
        merged_psam["PSEUDOHAPLOID"] = [s.value for s in statuses]

        # Step 15: write .psam
        psam.write_psam(merged_psam, out_psam_path)

    if not quiet:
        _print_simple_summary(
            counters,
            descriptors,
            out_pgen_path,
            out_pvar_path,
            out_psam_path,
        )


def _print_simple_summary(
    counters: MergeCounters,
    descriptors: list[InputDescriptor],
    out_pgen: Path,
    out_pvar: Path,
    out_psam: Path,
) -> None:
    """Stdout one-liner. Day 4 replaces this with the detailed HLD summary block."""
    print(f"Read {len(descriptors)} inputs:", file=sys.stdout)
    for i, d in enumerate(descriptors):
        print(f"  [{i}] {d.path}: {d.n_samples} samples", file=sys.stdout)
    print("", file=sys.stdout)
    print("Variant alignment counts (action_histogram):", file=sys.stdout)
    for k, v in counters.action_histogram.items():
        print(f"  {k}\t{v}", file=sys.stdout)
    print("", file=sys.stdout)
    pseudo_counts: dict[str, int] = {"0": 0, "1": 0, "U": 0}
    for _, h, c in counters.per_sample_het:
        if c == 0:
            pseudo_counts["U"] += 1
        elif h == 0:
            pseudo_counts["1"] += 1
        elif h / c >= 0.05:
            pseudo_counts["0"] += 1
        else:
            pseudo_counts["U"] += 1
    print(
        f"Pseudohaploid: {pseudo_counts['1']} pseudo, "
        f"{pseudo_counts['0']} diploid, {pseudo_counts['U']} unknown",
        file=sys.stdout,
    )
    print("", file=sys.stdout)
    print("Wrote:", file=sys.stdout)
    print(f"  {out_pgen} ({out_pgen.stat().st_size:,} B)", file=sys.stdout)
    print(f"  {out_pvar} ({out_pvar.stat().st_size:,} B)", file=sys.stdout)
    print(f"  {out_psam} ({out_psam.stat().st_size:,} B)", file=sys.stdout)
    print(
        f"Done: {counters.n_output_samples} samples x {counters.n_output_variants} variants",
        file=sys.stdout,
    )
