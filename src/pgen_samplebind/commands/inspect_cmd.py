"""`inspect` subcommand orchestrator.

Emits format / sample count / variant count / populations / sex distribution
/ per-sample missingness histogram via PgenReader iteration.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from ..errors import InvariantViolation
from ..formats import prepared_input
from ..psam import detect_population_column, read_psam, rename_to_pop
from ..pvar import (
    check_max_alleles,
    check_pvar_pgen_row_count_consistent,
    count_raw_variants,
    read_pvar,
)

# Per-sample-missing-rate bin edges: 10 deciles spanning [0, 1.0001) so that
# rate == 1.0 (100% missing) lands in the final bin rather than escaping it.
_HISTOGRAM_BIN_EDGES = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001])
_HISTOGRAM_BIN_LABELS = (
    "0-10%",
    "10-20%",
    "20-30%",
    "30-40%",
    "40-50%",
    "50-60%",
    "60-70%",
    "70-80%",
    "80-90%",
    "90-100%",
)

# Block size for the genotype iteration. Per  step 5: not user-tunable
# for inspect; matches the merge default.
_INSPECT_BLOCK_SIZE = 2048


def _compute_missingness_histogram(
    pgen_path: Path,
    n_samples: int,
    n_variants: int,
    block_size: int = _INSPECT_BLOCK_SIZE,
) -> dict[str, Any] | None:
    """Per-sample missing-genotype histogram via PgenReader.read_range iteration.

    Returns a dict with `bins` (label → count) and `n_variants_scanned`, or
    None if scanning would be trivially empty.
    """
    if n_samples == 0 or n_variants == 0:
        return None

    import pgenlib

    missing_counts = np.zeros(n_samples, dtype=np.int64)

    # PgenReader requires raw_sample_ct; variant_ct is optional but passing it
    # avoids a silent header-read on every method call.
    reader = pgenlib.PgenReader(str(pgen_path).encode(), raw_sample_ct=n_samples)
    try:
        for start in range(0, n_variants, block_size):
            end = min(start + block_size, n_variants)
            buf = np.empty((end - start, n_samples), dtype=np.int8)
            reader.read_range(start, end, buf, sample_maj=0)
            # Per HLD: missing genotypes are encoded as -9 in pgenlib's int8 buffers.
            missing_counts += (buf == -9).sum(axis=0).astype(np.int64)
    finally:
        # Pgenlib readers expose close(); fall back to GC if not present.
        if hasattr(reader, "close"):
            reader.close()

    # Per-sample missing rate = missing_count / total_variants.
    missing_rates = missing_counts.astype(np.float64) / float(n_variants)
    counts, _ = np.histogram(missing_rates, bins=_HISTOGRAM_BIN_EDGES)

    bins = {label: int(c) for label, c in zip(_HISTOGRAM_BIN_LABELS, counts.tolist(), strict=True)}
    return {"bins": bins, "n_variants_scanned": int(n_variants)}


def _build_summary(input_path: Path) -> dict[str, Any]:
    """Build the structured inspect summary."""
    with prepared_input(input_path) as desc:
        # Sample-side
        psam_df = read_psam(desc.psam_path)
        pop_col = detect_population_column(psam_df, override=None)
        psam_df = rename_to_pop(psam_df, pop_col)
        n_samples = len(psam_df)

        pop_counts = Counter(psam_df["POP"].tolist())
        sex_counts: dict[str, int] = {}
        if "SEX" in psam_df.columns:
            sex_counts = dict(Counter(psam_df["SEX"].tolist()))

        # Variant-side
        n_pre_filter = count_raw_variants(desc.pvar_path)
        pvar_df = read_pvar(desc.pvar_path)
        n_post_filter = len(pvar_df)

        chrom_counts = dict(Counter(pvar_df["chrom"].tolist()))
        chrom_counts_sorted = {int(k): int(v) for k, v in sorted(chrom_counts.items())}

        # Missingness histogram via PgenReader iteration.
        # Multi-allelic input would SIGSEGV inside read_range, and a
        # pvar/pgen row-count mismatch would either over- or under-read
        # the .pgen and produce a garbage histogram; both soft-skip with
        # a structured reason rather than failing the whole inspect.
        missingness: dict[str, Any] | None
        try:
            check_max_alleles(desc.pgen_path)
            check_pvar_pgen_row_count_consistent(desc.pgen_path)
        except InvariantViolation as e:
            missingness = {"status": "skipped", "reason": str(e)}
        else:
            missingness = _compute_missingness_histogram(desc.pgen_path, n_samples, n_pre_filter)

        return {
            "input_path": str(input_path),
            "format": desc.fmt.value,
            "n_samples": n_samples,
            "n_variants_pre_filter": n_pre_filter,
            "n_variants_post_filter_biallelic_snp": n_post_filter,
            "n_populations": len(pop_counts),
            "populations": dict(pop_counts.most_common()),
            "sex_distribution": sex_counts,
            "variants_per_chrom": chrom_counts_sorted,
            "missingness_histogram": missingness,
        }


def _format_text(summary: dict[str, Any]) -> str:
    """TSV-style human-readable output."""
    lines = []
    lines.append(f"input_path\t{summary['input_path']}")
    lines.append(f"format\t{summary['format']}")
    lines.append(f"n_samples\t{summary['n_samples']}")
    lines.append(f"n_variants_pre_filter\t{summary['n_variants_pre_filter']}")
    lines.append(
        f"n_variants_post_filter_biallelic_snp\t{summary['n_variants_post_filter_biallelic_snp']}"
    )
    lines.append(f"n_populations\t{summary['n_populations']}")
    lines.append("")
    lines.append("populations:")
    for pop, n in summary["populations"].items():
        lines.append(f"  {pop}\t{n}")
    if summary["sex_distribution"]:
        lines.append("")
        lines.append("sex_distribution:")
        for sex, n in sorted(summary["sex_distribution"].items()):
            lines.append(f"  {sex}\t{n}")
    lines.append("")
    lines.append("variants_per_chrom:")
    for chrom, n in summary["variants_per_chrom"].items():
        lines.append(f"  {chrom}\t{n}")

    histogram = summary["missingness_histogram"]
    lines.append("")
    if histogram is None:
        lines.append("missingness_histogram\t(no variants to scan)")
    elif histogram.get("status") == "skipped":
        lines.append(f"missingness_histogram\tskipped: {histogram['reason']}")
    else:
        n_scanned = histogram["n_variants_scanned"]
        lines.append(f"missingness_histogram (per-sample missing rate over {n_scanned} variants):")
        for label, count in histogram["bins"].items():
            lines.append(f"  {label}\t{count}")
    return "\n".join(lines) + "\n"


def run_inspect(input_path: Path, json_output: bool) -> None:
    """Build a structured summary of one input and emit it."""
    summary = _build_summary(input_path)
    if json_output:
        json.dump(summary, sys.stdout, indent=2, default=str)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(_format_text(summary))
