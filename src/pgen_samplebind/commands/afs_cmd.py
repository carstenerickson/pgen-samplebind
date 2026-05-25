"""`afs` subcommand orchestrator. Bridges PFILE → AT2-compatible AFS TSVs."""

from __future__ import annotations

import json
import sys
import time
from contextlib import ExitStack
from pathlib import Path

from .. import __version__
from ..afs import compute_afs, write_afs_tsvs
from ..formats import prepared_input
from ..pvar import check_max_alleles, check_pvar_pgen_row_count_consistent


def run_afs(
    input_path: Path,
    output_dir: Path,
    population_column: str | None,
    populations_filter: tuple[str, ...] | None,
    adjust_pseudohaploid: bool,
    include_sex_chrom: bool,
    block_size: int,
    quiet: bool,
) -> None:
    """`afs` subcommand orchestrator.

    Reads one input (PFILE / BFILE / EIGENSTRAT — auto-detected), computes
    per-population allele frequencies + called-allele counts, writes three
    TSVs into output_dir matching AT2's `*_to_afs()` shape.

    A small R loader (`scripts/load_pgensb_afs.R` in the repo) turns the
    TSVs into the AT2-style three-data-frame list for downstream f-statistic
    work — bridge until `pfile_to_afs()` lands in admixtools upstream.
    """
    started = time.perf_counter()

    include_chrom: tuple[int, ...] = (
        tuple(range(1, 27)) if include_sex_chrom else tuple(range(1, 23))
    )

    with ExitStack() as stack:
        desc = stack.enter_context(
            prepared_input(input_path, is_target=False, include_chrom=include_chrom)
        )

        # Multi-allelic input would SIGSEGV inside pgenlib's read_range
        # (same C-layer crash check_max_alleles guards against in merge /
        # validate / inspect — afs hits the same C path). Row-count
        # consistency catches mis-paired triplets that would otherwise
        # over- or under-read the .pgen and silently corrupt frequencies.
        check_max_alleles(desc.pgen_path)
        check_pvar_pgen_row_count_consistent(desc.pgen_path)

        # compute_afs resolves --population-column via the same
        # POP/PHENO/PHENO1 auto-detect that merge and validate use when
        # None is passed, then reads the .psam once.
        result = compute_afs(
            descriptor=desc,
            population_column=population_column,
            populations=list(populations_filter) if populations_filter else None,
            adjust_pseudohaploid=adjust_pseudohaploid,
            include_chrom=include_chrom,
            block_size=block_size,
        )

    paths = write_afs_tsvs(result, output_dir)

    # Manifest with metadata for downstream consumers (R loader uses this).
    manifest = {
        "tool": "pgen-samplebind",
        "tool_version": __version__,
        "command": "afs",
        "input": str(input_path),
        "n_variants": len(result.snp),
        "n_populations": len(result.populations),
        "populations": result.populations,
        "n_samples_per_pop": result.n_samples_per_pop,
        "adjust_pseudohaploid_requested": adjust_pseudohaploid,
        "adjust_pseudohaploid_applied": result.adjust_pseudohaploid_applied,
        "include_chrom": list(include_chrom),
        "files": {k: v.name for k, v in paths.items()},
    }
    manifest_path = output_dir / "afs_manifest.json"
    try:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    except OSError as e:
        from ..errors import IOFailure

        raise IOFailure(f"cannot write {manifest_path}: {e}") from e

    elapsed = time.perf_counter() - started
    if not quiet:
        sys.stdout.write(
            f"AFS computed: {len(result.snp):,} variants x {len(result.populations)} populations "
            f"in {elapsed:.2f}s.\n"
            f"  pseudohaploid adjustment: "
            f"{'applied' if result.adjust_pseudohaploid_applied else 'not applied'}\n"
            f"  output: {output_dir}/\n"
            f"    afs_snp.tsv     ({len(result.snp):,} rows)\n"
            f"    afs_freq.tsv    ({len(result.freq):,} rows x {len(result.populations) + 1} cols)\n"
            f"    afs_counts.tsv  ({len(result.counts):,} rows x {len(result.populations) + 1} cols)\n"
            f"    afs_manifest.json\n"
        )
