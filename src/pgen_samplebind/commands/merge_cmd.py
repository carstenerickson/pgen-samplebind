"""`merge` subcommand orchestrator. Sequence in """

from __future__ import annotations

import sys
import time
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__, psam, pseudohaploid, reporting
from ..concurrency import output_lock
from ..errors import InvariantViolation, PgenSamplebindError
from ..formats import prepared_input
from ..merge import merge_inputs
from ..pvar import check_max_alleles, count_raw_variants
from ..types import (
    InputDescriptor,
    MergeContext,
    MergePolicy,
    PseudohaploidStatus,
    SampleIdentityPlan,
)


def _collect_sidecar_overrides(
    descriptors: list[InputDescriptor],
    psam_dfs: list[pd.DataFrame],
    sample_plan: SampleIdentityPlan,
) -> dict[str, PseudohaploidStatus]:
    """Walk every input's pseudohaploid sidecar (if present) and produce a
    flat output-IID → status map.

    Threading semantics: the sidecar files keyed by INPUT IIDs; the merge
    orchestrator emits OUTPUT IIDs that may have been renamed by the
    `_target` / `_<input_idx>` suffix scheme under --on-collision suffix.
    We map input → output via `sample_plan.per_input_output_indices`,
    skipping samples that the collision plan dropped (keep_mask == False).

    Raises:
        InvariantViolation: a sidecar carries IIDs not present in its
            input's .psam — orphaned entries signal upstream/downstream
            disagreement and shouldn't be silently ignored.
    """
    overrides: dict[str, PseudohaploidStatus] = {}
    for input_idx, (desc, df) in enumerate(zip(descriptors, psam_dfs, strict=True)):
        sidecar = pseudohaploid.read_sidecar(desc.path)
        if sidecar is None:
            continue

        input_iids = df["IID"].astype(str).tolist()
        input_iid_set = set(input_iids)
        orphans = sorted(set(sidecar.keys()) - input_iid_set)
        if orphans:
            preview = orphans[:5]
            tail = f" ... +{len(orphans) - 5} more" if len(orphans) > 5 else ""
            raise InvariantViolation(
                f"input[{input_idx}] pseudohaploid sidecar at {desc.path} lists "
                f"{len(orphans)} sample(s) not present in the input .psam: "
                f"{preview}{tail}. Fix the upstream tool or the .psam to match."
            )

        keep_mask = sample_plan.per_input_keep_mask[input_idx]
        out_indices = sample_plan.per_input_output_indices[input_idx]
        for i, iid in enumerate(input_iids):
            if not keep_mask[i] or iid not in sidecar:
                continue
            out_iid = sample_plan.output_iids[int(out_indices[i])]
            overrides[out_iid] = sidecar[iid]
    return overrides


def _unlink_output_triplet(*paths: Path) -> None:
    """Best-effort unlink of partial output files. Per  fix #6:
    on PgenSamplebindError mid-pass-2 / psam-finalization, the orchestrator
    unlinks the .pgen / .pvar / .psam triplet so downstream pipelines never
    silently consume a half-built output. Atomic-rename across the triplet
    isn't actually atomic (three separate renames; NFS doesn't even
    guarantee single-file atomicity), so unlink-on-failure is the simpler
    and substantively equivalent choice.
    """
    for p in paths:
        try:
            if p.exists():
                p.unlink()
        except OSError:
            # Best-effort; if we can't unlink (permissions, fs error), the
            # original PgenSamplebindError is more important.
            pass


def run_merge(
    input_paths: tuple[Path, ...],
    target_paths: tuple[Path, ...],
    output_prefix: Path,
    policy: MergePolicy,
    report_path: Path | None,
    report_json_path: Path | None,
    quiet: bool,
    relabel_from: Path | None = None,
    relabel_input_col: str | None = None,
    relabel_output_col: str | None = None,
) -> None:
    """Merge subcommand orchestrator. Per 

    Target mode (--target): one or more targets are appended after the
    positional inputs. The canonical (input[0]) remains the first positional.
    Each target's descriptor is marked is_target=True so merge_inputs's
    gate (c) call-rate check fires per-target, and resolve_sample_identity
    uses target_idxs for the `_target` / `_target_<idx>` suffix scheme
    (bare `_target` for the single-target case; `_target_<input_idx>` when
    multiple targets are supplied).
    """
    if target_paths and not input_paths:
        from ..errors import UsageError

        raise UsageError("--target requires at least one positional INPUT (the panel).")

    all_input_paths: tuple[Path, ...] = (*input_paths, *target_paths)
    # Indexes (in all_input_paths) of every target. Empty when none.
    target_idxs: tuple[int, ...] = tuple(len(input_paths) + i for i in range(len(target_paths)))
    target_idx_set: frozenset[int] = frozenset(target_idxs)

    started = time.perf_counter()

    out_pgen_path = Path(str(output_prefix) + ".pgen")
    out_pvar_path = Path(str(output_prefix) + ".pvar")
    out_psam_path = Path(str(output_prefix) + ".psam")
    output_paths = {"pgen": out_pgen_path, "pvar": out_pvar_path, "psam": out_psam_path}

    with ExitStack() as stack:
        # Step 3: advisory output-prefix lock. Acquired BEFORE
        # any input read so a held-lock failure exits 2 without touching
        # inputs. Released on context exit. NFS/SMB/CIFS triggers a stderr
        # warning since flock semantics there are advisory-only-on-paper.
        stack.enter_context(output_lock(output_prefix))

        # Step 1: format detection per input via context manager. Target
        # descriptors are marked is_target=True so gate (c) finds them.
        descriptors = [
            stack.enter_context(
                prepared_input(
                    p,
                    is_target=(i in target_idx_set),
                    include_chrom=policy.include_chrom,
                )
            )
            for i, p in enumerate(all_input_paths)
        ]

        # Step 5: per-input multi-allelic startup check
        for desc in descriptors:
            check_max_alleles(desc.pgen_path)

        # Step 6: read psams; detect population column; rename → POP.
        # NOTE: add_fid_from_pop must run AFTER --relabel-from (Day 9), because
        # the relabel may change POP and FID should mirror the final POP.
        psam_dfs = []
        for desc in descriptors:
            df = psam.read_psam(desc.psam_path)
            pop_col = psam.detect_population_column(df, policy.population_column)
            df = psam.rename_to_pop(df, pop_col)
            psam_dfs.append(df)

        # Step 7: --relabel-from per-input. Applied per-input
        # before sample-bind so different inputs can have independent relabels.
        if relabel_from is not None:
            relabel_df = psam.read_relabel_tsv(relabel_from, relabel_input_col, relabel_output_col)
            # 2-col form (no input/output col flags) → source = POP (collapse).
            # N-col form (flags supplied) → source = id_column (per-sample override).
            source_col = "POP" if relabel_input_col is None else policy.id_column
            psam_dfs = [
                psam.apply_relabel(df, relabel_df, source_column=source_col) for df in psam_dfs
            ]

        # Now FID = POP (after any relabel applied)
        psam_dfs = [psam.add_fid_from_pop(df) for df in psam_dfs]

        # Populate descriptor n_samples (from psam) and n_variants (cheap raw line
        # count from .pvar) — both used by reporting; merge_inputs reads pvars
        # internally and doesn't depend on n_variants here.
        descriptors = [
            replace(d, n_samples=len(df), n_variants=count_raw_variants(d.pvar_path))
            for d, df in zip(descriptors, psam_dfs, strict=True)
        ]

        # Step 10: resolve sample identity (collision policy applied; target_idxs
        # drive the `_target` / `_target_<input_idx>` suffix scheme under
        # --on-collision suffix).
        sample_plan = psam.resolve_sample_identity(psam_dfs, policy, target_idxs=target_idxs)

        # Step 10b: per-input pseudohaploid sidecar (issue #2). Upstream tools
        # like pileup-aadr write `<prefix>.pseudohaploid.json` to assert
        # per-sample status by construction; we map sidecar IIDs through the
        # collision plan to output IIDs so Step 14 can override classification.
        sidecar_overrides = _collect_sidecar_overrides(descriptors, psam_dfs, sample_plan)

        # Step 11: build MergeContext. Progress bar is on by default when
        # stderr is a tty and the user didn't pass --quiet — workflow
        # managers (Snakemake, Nextflow) pipe stderr and so see a silent run.
        ctx = MergeContext(
            policy=policy,
            sample_plan=sample_plan,
            report_tsv_path=report_path,
            collect_variant_rows=(policy.report_json_include_rows and report_json_path is not None),
            show_progress=(not quiet and sys.stderr.isatty()),
            quiet=quiet,
        )

        # Steps 12-15: merge_inputs + psam finalization, wrapped in the
        # output-cleanup wrapper fix #6. On any PgenSamplebindError
        # (gate failure, IO failure, invariant violation), unlink the partial
        # triplet before re-raising so downstream pipelines never consume a
        # half-built output.
        try:
            # Step 12: merge_inputs (pass 1 + gates + pass 2; writes .pgen +
            # .pvar + report TSV if ctx.report_tsv_path; populates
            # counters.variant_rows if ctx.collect_variant_rows).
            counters = merge_inputs(descriptors, out_pgen_path, out_pvar_path, ctx)

            # Step 13: psam finalization
            merged_psam = psam.merge_psams(psam_dfs, sample_plan)

            # Step 14: classify pseudohaploid; assign PSEUDOHAPLOID column.
            # Row-order invariant: the IID at position i must agree across
            # the pass-2 counters (genotype-column order), the sample plan
            # (output-IID order), and the merged psam (row order). A
            # mismatch here means a refactor broke the sample-axis
            # alignment — the failure mode is silent genotype/sample
            # misassignment in the output, so we enforce unconditionally
            # rather than under `__debug__`.
            for i, (iid_in_counters, _, _) in enumerate(counters.per_sample_het):
                if not (
                    iid_in_counters == sample_plan.output_iids[i] == merged_psam.iloc[i]["IID"]
                ):
                    raise InvariantViolation(
                        f"row-order invariant violated at i={i}: "
                        f"counters={iid_in_counters!r}, "
                        f"plan={sample_plan.output_iids[i]!r}, "
                        f"psam={merged_psam.iloc[i]['IID']!r}"
                    )

            het_array = np.array([h for _, h, _ in counters.per_sample_het], dtype=np.int64)
            called_array = np.array([c for _, _, c in counters.per_sample_het], dtype=np.int64)
            statuses = pseudohaploid.classify_all(het_array, called_array)
            # Sidecar overrides (issue #2) take precedence over heterozygosity
            # inference when the upstream tool provided an authoritative status.
            merged_psam["PSEUDOHAPLOID"] = [
                sidecar_overrides.get(out_iid, status).value
                for out_iid, status in zip(sample_plan.output_iids, statuses, strict=True)
            ]

            # Step 15: write .psam
            psam.write_psam(merged_psam, out_psam_path)
        except PgenSamplebindError:
            _unlink_output_triplet(out_pgen_path, out_pvar_path, out_psam_path)
            raise

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
