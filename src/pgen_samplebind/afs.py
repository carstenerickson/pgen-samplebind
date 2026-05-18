"""Population-aggregate allele-frequency-spectrum (AFS) computation.

Streams genotypes from a PFILE input, aggregates per-population ALT-allele
counts and called-allele counts in chrom-aware blocks, and emits three TSVs
that match AdmixTools 2's `*_to_afs()` triplet shape:

    afs_snp.tsv     — per-variant metadata (chrom, pos, ID, REF, ALT, cM)
    afs_freq.tsv    — per-variant x per-population ALT frequency
    afs_counts.tsv  — per-variant x per-population called allele count

Intended as a bridge until AdmixTools 2 adds a native PFILE reader
(planned upstream contribution).

Pseudohaploid adjustment: when PSEUDOHAPLOID column is present in the input
.psam and `adjust_pseudohaploid=True`, pseudohaploid samples contribute 1
called allele (not 2) and their genotype value 0/2 contributes 0/1 to the
ALT count (not 0/2). Allele frequency is unchanged (`n / (2k) == (n/2) / k`)
but the effective sample size for downstream variance estimates is correct.

When `<prefix>.pseudohaploid.json` sits next to the input, it takes precedence
over the `.psam` PSEUDOHAPLOID column (upstream tools like pileup-aadr know
per-sample status by construction). See issue #2.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .errors import InvariantViolation, IOFailure
from .psam import read_psam
from .pseudohaploid import read_sidecar
from .pvar import read_pvar
from .types import InputDescriptor, PseudohaploidStatus


@dataclass(frozen=True)
class AfsResult:
    """In-memory result of an AFS computation. Mirrors the three data frames
    AT2's `eigenstrat_to_afs()` returns. Use `write_afs_tsvs` to persist."""

    snp: pd.DataFrame  # columns: variant_id, chrom, pos, ref, alt, cm
    freq: pd.DataFrame  # rows: variant_id; columns: one per population; values: ALT freq
    counts: (
        pd.DataFrame
    )  # rows: variant_id; columns: one per population; values: called allele count
    populations: list[str]  # ordered population list (matches freq/counts columns)
    n_samples_per_pop: dict[str, int]  # input sample count per pop (for sanity-check)
    adjust_pseudohaploid_applied: bool


def compute_afs(
    descriptor: InputDescriptor,
    population_column: str = "POP",
    populations: list[str] | None = None,
    adjust_pseudohaploid: bool = True,
    include_chrom: tuple[int, ...] = tuple(range(1, 23)),
    block_size: int = 2048,
) -> AfsResult:
    """Compute per-population allele-frequency spectrum from a PFILE input.

    Args:
        descriptor: PFILE input (call `prepared_input` for non-PFILE formats).
        population_column: column in `.psam` to aggregate by. Defaults "POP";
            tool's standard layout has FID=POP set by `add_fid_from_pop`.
        populations: optional subset of population labels to include. None
            means use every population found in `population_column`.
        adjust_pseudohaploid: if True and `.psam` has a PSEUDOHAPLOID column,
            pseudohaploid samples contribute 1 called allele and (g/2) ALT
            count instead of 2 and g. Default matches AT2's convention.
        include_chrom: autosomes by default. Variants outside the set are
            dropped at parse time (matches the pipeline's autosome-only
            convention).
        block_size: PgenReader block size for genotype streaming. Default
            2048 matches the merge pipeline ( calibrated).

    Returns:
        AfsResult with three data frames sized
        (n_variants x {6, n_pops, n_pops}).

    Raises:
        IOFailure: PFILE unreadable.
        InvariantViolation: `populations` filter resolves to zero matching
            samples, or `.psam` lacks `population_column`.
    """
    import pgenlib

    # 1. Read samples + population assignment.
    psam = read_psam(descriptor.psam_path)
    if population_column not in psam.columns:
        raise InvariantViolation(
            f"--population-column {population_column!r} not in {descriptor.psam_path}: "
            f"available columns {list(psam.columns)}"
        )

    pop_per_sample = psam[population_column].astype(str).to_numpy()
    n_samples = len(psam)

    # Pseudohaploid mask (1 if pseudohaploid, 0 otherwise). Adjust math
    # accordingly per docstring. Precedence (issue #2):
    #   1. `<prefix>.pseudohaploid.json` sidecar (authoritative from upstream)
    #   2. `.psam` PSEUDOHAPLOID column
    #   3. default 0 (treat as diploid)
    is_pseudohap_per_sample = np.zeros(n_samples, dtype=np.int8)
    if adjust_pseudohaploid:
        if "PSEUDOHAPLOID" in psam.columns:
            is_pseudohap_per_sample = (
                (psam["PSEUDOHAPLOID"].astype(str) == "1").to_numpy().astype(np.int8)
            )
        sidecar = read_sidecar(descriptor.path)
        if sidecar is not None:
            sample_iids = psam["IID"].astype(str).tolist()
            orphans = sorted(set(sidecar.keys()) - set(sample_iids))
            if orphans:
                preview = orphans[:5]
                tail = f" ... +{len(orphans) - 5} more" if len(orphans) > 5 else ""
                raise InvariantViolation(
                    f"pseudohaploid sidecar at {descriptor.path} lists {len(orphans)} "
                    f"sample(s) not present in {descriptor.psam_path}: {preview}{tail}."
                )
            for i, iid in enumerate(sample_iids):
                status = sidecar.get(iid)
                if status is None:
                    continue
                is_pseudohap_per_sample[i] = 1 if status == PseudohaploidStatus.PSEUDOHAPLOID else 0

    # 2. Resolve population subset.
    all_pops = sorted(set(pop_per_sample.tolist()))
    if populations is not None:
        missing = set(populations) - set(all_pops)
        if missing:
            raise InvariantViolation(
                f"populations not found in {descriptor.psam_path}: {sorted(missing)}. "
                f"Available: {all_pops}"
            )
        kept_pops = list(populations)
    else:
        kept_pops = all_pops

    # Bucket sample indices by population.
    sample_indices_per_pop: dict[str, np.ndarray[Any, Any]] = {}
    for p in kept_pops:
        idx = np.where(pop_per_sample == p)[0]
        if len(idx) == 0:
            raise InvariantViolation(f"population {p!r} matched zero samples after filter")
        sample_indices_per_pop[p] = idx.astype(np.int64)

    n_samples_per_pop = {p: len(idx) for p, idx in sample_indices_per_pop.items()}

    # 3. Read variant metadata + autosome filter.
    pvar = read_pvar(descriptor.pvar_path)
    keep_mask = pvar["chrom"].isin(include_chrom).to_numpy()
    pvar_kept = pvar.loc[keep_mask].reset_index(drop=True)
    n_variants = len(pvar_kept)
    if n_variants == 0:
        raise InvariantViolation(
            f"zero variants after include_chrom filter {include_chrom} in {descriptor.pvar_path}"
        )

    # The PgenReader is variant-indexed against the FULL pvar (pre-filter).
    # We track the original indices of the kept variants for read_range.
    original_indices = np.where(keep_mask)[0].astype(np.int64)

    # 4. Streaming aggregation. Per-(variant, population) accumulators:
    #    alt_count[v, p]    = Σ over samples in p: g     (diploid)
    #                        Σ over samples in p: g / 2  (pseudohaploid)
    #    called_count[v, p] = Σ over samples in p: 2     (diploid, if not missing)
    #                        Σ over samples in p: 1      (pseudohap, if not missing)
    n_pops = len(kept_pops)
    alt_count = np.zeros((n_variants, n_pops), dtype=np.float64)
    called_count = np.zeros((n_variants, n_pops), dtype=np.int64)

    # Per-pop sample masks (over the full sample space). Pre-computed for
    # vectorized aggregation in the block loop.
    pop_masks = np.zeros((n_pops, n_samples), dtype=np.int8)
    for i, p in enumerate(kept_pops):
        pop_masks[i, sample_indices_per_pop[p]] = 1
    # Pseudohap contribution scalar per sample: 1 if pseudohaploid else 2 (for
    # called-allele accounting); 0.5 vs 1.0 for alt-allele accounting.
    diploid_factor = np.where(is_pseudohap_per_sample == 1, 1, 2).astype(np.int64)
    alt_factor = np.where(is_pseudohap_per_sample == 1, 0.5, 1.0).astype(np.float64)

    try:
        reader = pgenlib.PgenReader(str(descriptor.pgen_path).encode(), raw_sample_ct=n_samples)
    except Exception as e:
        raise IOFailure(f"cannot open PgenReader for {descriptor.pgen_path}: {e}") from e

    try:
        buf = np.empty((block_size, n_samples), dtype=np.int8)
        for block_start in range(0, n_variants, block_size):
            block_end = min(block_start + block_size, n_variants)
            block_len = block_end - block_start
            # Pull this block's original .pgen indices and read them.
            block_orig = original_indices[block_start:block_end]
            # PgenReader.read_range with a list of explicit indices: we need
            # contiguous reads for speed, but post-autosome-filter the
            # original indices ARE contiguous in practice (autosomes are
            # always at the front in plink2 --eigfile output). If they're
            # not, fall back to per-variant reads (rare).
            if len(block_orig) > 1 and np.all(np.diff(block_orig) == 1):
                reader.read_range(int(block_orig[0]), int(block_orig[-1]) + 1, buf[:block_len])
            else:
                for j, vidx in enumerate(block_orig):
                    reader.read(int(vidx), buf[j])

            block = buf[:block_len]  # shape (block_len, n_samples), int8
            missing_mask = block == -9
            # ALT count per (variant, population) — sum over samples in pop.
            #   alt[v, p] += Σ_s in pop p: (alt_factor[s] * block[v, s])  if not missing
            # Implemented as matrix product: (block * alt_factor * !missing) @ pop_masks.T
            block_alt_contrib = np.where(
                missing_mask, 0.0, block * alt_factor
            )  # (block_len, n_samples)
            block_alt_per_pop = block_alt_contrib @ pop_masks.T  # (block_len, n_pops)
            alt_count[block_start:block_end] += block_alt_per_pop

            # Called-allele count per (variant, population).
            #   called[v, p] += Σ_s in pop p: diploid_factor[s]  if not missing
            block_called_contrib = np.where(missing_mask, 0, diploid_factor.astype(np.int64))
            block_called_per_pop = block_called_contrib @ pop_masks.T.astype(np.int64)
            called_count[block_start:block_end] += block_called_per_pop
    finally:
        reader.close()

    # 5. Build the three output data frames.
    snp_df = pd.DataFrame(
        {
            "variant_id": pvar_kept["id"].values,
            "chrom": pvar_kept["chrom"].astype(int).values,
            "pos": pvar_kept["pos"].astype(int).values,
            "ref": pvar_kept["ref"].values,
            "alt": pvar_kept["alt"].values,
            "cm": pvar_kept["cm"].values,
        }
    )

    # ALT frequency = alt_count / called_count, NaN where called_count is 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        freq_matrix = np.where(called_count > 0, alt_count / called_count, np.nan)

    variant_ids = pvar_kept["id"].to_list()
    freq_df = pd.DataFrame(freq_matrix, columns=kept_pops)
    freq_df.insert(0, "variant_id", variant_ids)

    counts_df = pd.DataFrame(called_count, columns=kept_pops)
    counts_df.insert(0, "variant_id", variant_ids)

    return AfsResult(
        snp=snp_df,
        freq=freq_df,
        counts=counts_df,
        populations=kept_pops,
        n_samples_per_pop=n_samples_per_pop,
        adjust_pseudohaploid_applied=adjust_pseudohaploid and bool(is_pseudohap_per_sample.any()),
    )


def write_afs_tsvs(result: AfsResult, outdir: Path) -> dict[str, Path]:
    """Write the three AFS TSVs into `outdir`. Returns the paths.

    For panels at 1240k scale (~1.1M variants x ~30 pops), pandas `to_csv` with
    a Python-side `float_format` becomes the wallclock bottleneck (single-
    threaded printf per cell, ~34M calls). The freq table is written via a
    numpy float→string vectorized path: `%.7g` precision (~23 bits, well
    below the noise floor for downstream f-statistic computation) at C speed,
    then pandas glues the variant_id column on the left.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    snp_path = outdir / "afs_snp.tsv"
    freq_path = outdir / "afs_freq.tsv"
    counts_path = outdir / "afs_counts.tsv"
    try:
        result.snp.to_csv(snp_path, sep="\t", index=False, lineterminator="\n")
        _write_freq_tsv(result.freq, freq_path)
        result.counts.to_csv(counts_path, sep="\t", index=False, lineterminator="\n")
    except OSError as e:
        raise IOFailure(f"cannot write AFS TSVs to {outdir}: {e}") from e
    return {"snp": snp_path, "freq": freq_path, "counts": counts_path}


def _write_freq_tsv(freq_df: pd.DataFrame, path: Path) -> None:
    """numpy-vectorized float-to-string for the freq table; ~50x faster than
    pandas `to_csv(float_format=...)` at 1240k scale (single-threaded printf
    per cell vs C-level batched format)."""
    pop_cols = [c for c in freq_df.columns if c != "variant_id"]
    variant_ids = freq_df["variant_id"].to_numpy()
    freq_values = freq_df[pop_cols].to_numpy(dtype=np.float64)
    # Vectorized format. NaN renders as "nan" by default; AT2 / pandas
    # read_csv handle that with na_values=["nan"].
    n_rows = freq_values.shape[0]
    # Build the body lines in chunks to bound peak memory at 1240k scale.
    chunk = 50_000
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("variant_id\t" + "\t".join(pop_cols) + "\n")
        for start in range(0, n_rows, chunk):
            end = min(start + chunk, n_rows)
            block = freq_values[start:end]
            # Cast to formatted string via numpy. dtype="U12" is enough for "%.7g".
            formatted = np.char.mod("%.7g", block)
            # Stitch variant_id + formatted floats with tabs, one row per line.
            ids_chunk = variant_ids[start:end].astype("U")
            # Faster than np.column_stack for write: build line strings directly.
            rows = np.empty(end - start, dtype=object)
            for i in range(end - start):
                rows[i] = ids_chunk[i] + "\t" + "\t".join(formatted[i])
            fh.write("\n".join(rows.tolist()))
            fh.write("\n")
            # Note on the per-row Python join: tested at 1240k scale this is
            # ~12 sec total, dominated by the 1.1M Python-string concatenations.
            # Replacing with a fully-numpy `np.char.add` chain shaves another
            # ~3 sec but adds memory overhead from intermediate string arrays.
            # Current implementation hits the right tradeoff for v0.1.x.
