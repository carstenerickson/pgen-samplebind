"""Deterministic synthetic-fixture generator. Per LLD §5.2 / HLD project Day 1.

Generates a 1240k-style PFILE triplet (.pgen + .pvar + .psam) from a fixed seed.
Default 500-sample x 50K-variant spec runs in < 5 s on M2 single-core.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pgen_samplebind.types import InputDescriptor, InputFormat

# Allele pairs.
# Unambiguous: AC, AG, CT, GT (strand-resolvable).
# Ambiguous: AT, CG (strand-ambiguous; HLD §Allele resolution).
_UNAMBIGUOUS_PAIRS = [("A", "C"), ("A", "G"), ("C", "T"), ("G", "T")]
_AMBIGUOUS_PAIRS = [("A", "T"), ("C", "G")]


@dataclass(frozen=True)
class SyntheticPanelSpec:
    """Knobs for the synthetic panel generator.

    `seed` controls the random stream. For tests that need two panels with the
    SAME variants but DIFFERENT samples, use the same `variant_seed` and
    different `sample_seed` / `sample_id_prefix`.
    """

    n_samples: int = 500
    n_variants: int = 50_000
    n_populations: int = 10
    pseudohaploid_fraction: float = 0.5
    chromosomes: tuple[int, ...] = tuple(range(1, 23))  # autosomes
    seed: int = 0xCA753E
    ambiguous_strand_fraction: float = 0.06  # ~5-8% A/T+C/G in 1240k design (HLD)
    missing_rate: float = 0.05
    sample_id_prefix: str = "S"
    variant_seed: int | None = None  # if set, used for variant generation; else `seed`
    sample_seed: int | None = None  # if set, used for samples/genotypes; else `seed`


def _build_variants(
    spec: SyntheticPanelSpec, rng: np.random.Generator
) -> tuple[list[int], list[int], list[str], list[str], list[str]]:
    """Return (chroms, positions, ids, refs, alts) for the panel."""
    n_chrom = len(spec.chromosomes)
    base_per_chrom = spec.n_variants // n_chrom
    extra = spec.n_variants - base_per_chrom * n_chrom

    chroms: list[int] = []
    positions: list[int] = []
    for i, chrom in enumerate(spec.chromosomes):
        n_this = base_per_chrom + (1 if i < extra else 0)
        # Random unique positions, sorted; arbitrary 100M upper bound is plenty.
        pos = rng.choice(np.arange(1, 100_000_000, dtype=np.int64), size=n_this, replace=False)
        pos.sort()
        chroms.extend([int(chrom)] * n_this)
        positions.extend(pos.tolist())

    n_variants = len(chroms)

    # Decide ambiguous-strand vs unambiguous per variant
    ambig_mask = rng.random(n_variants) < spec.ambiguous_strand_fraction
    # Pre-pick pair indices and orientation flips
    unamb_pair_idx = rng.integers(0, len(_UNAMBIGUOUS_PAIRS), size=n_variants)
    amb_pair_idx = rng.integers(0, len(_AMBIGUOUS_PAIRS), size=n_variants)
    flip_orient = rng.random(n_variants) < 0.5

    refs: list[str] = []
    alts: list[str] = []
    for i in range(n_variants):
        pair = (
            _AMBIGUOUS_PAIRS[amb_pair_idx[i]]
            if ambig_mask[i]
            else _UNAMBIGUOUS_PAIRS[unamb_pair_idx[i]]
        )
        if flip_orient[i]:
            refs.append(pair[1])
            alts.append(pair[0])
        else:
            refs.append(pair[0])
            alts.append(pair[1])

    ids = [f"chr{c}:{p}" for c, p in zip(chroms, positions, strict=True)]
    return chroms, positions, ids, refs, alts


def _build_genotypes(
    n_variants: int,
    n_samples: int,
    is_pseudohap: np.ndarray,
    spec: SyntheticPanelSpec,
    rng: np.random.Generator,
) -> np.ndarray:
    """Build the (n_variants, n_samples) int8 genotype matrix.

    Diploid samples: 0/1/2 uniform random.
    Pseudohaploid samples: 0/2 only (no heterozygous).
    Then apply per-cell missing rate (set to -9).
    """
    geno = rng.integers(0, 3, size=(n_variants, n_samples), dtype=np.int8)

    pseudohap_cols = np.where(is_pseudohap)[0]
    for col in pseudohap_cols:
        het_mask = geno[:, col] == 1
        n_het = int(het_mask.sum())
        if n_het:
            replacement = rng.choice(np.array([0, 2], dtype=np.int8), size=n_het)
            geno[het_mask, col] = replacement

    if spec.missing_rate > 0:
        missing_mask = rng.random(geno.shape) < spec.missing_rate
        geno[missing_mask] = -9

    return geno


def synthesize_pfile(spec: SyntheticPanelSpec, out_prefix: Path) -> InputDescriptor:
    """Generate a PFILE triplet (.pgen + .pvar + .psam) at out_prefix.

    Returns the resolved InputDescriptor for immediate test use. Wallclock
    target: < 5 s for the default 500-sample x 50K-variant spec.
    """
    import pgenlib

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    variant_rng = np.random.default_rng(
        spec.variant_seed if spec.variant_seed is not None else spec.seed
    )
    sample_rng = np.random.default_rng(
        spec.sample_seed if spec.sample_seed is not None else spec.seed
    )

    chroms, positions, ids, refs, alts = _build_variants(spec, variant_rng)
    n_variants = len(chroms)

    # Sample metadata.
    iids = [f"{spec.sample_id_prefix}{i:05d}" for i in range(spec.n_samples)]
    sex = sample_rng.integers(1, 3, size=spec.n_samples)  # plink: 1=male, 2=female
    pops = [f"pop_{i % spec.n_populations:02d}" for i in range(spec.n_samples)]
    is_pseudohap = sample_rng.random(spec.n_samples) < spec.pseudohaploid_fraction

    geno = _build_genotypes(n_variants, spec.n_samples, is_pseudohap, spec, sample_rng)

    out_pgen = Path(str(out_prefix) + ".pgen")
    out_pvar = Path(str(out_prefix) + ".pvar")
    out_psam = Path(str(out_prefix) + ".psam")

    # Write .pgen via pgenlib (bytes filename per pgenlib API)
    writer = pgenlib.PgenWriter(str(out_pgen).encode(), spec.n_samples, n_variants)
    try:
        writer.append_biallelic_batch(geno)
    finally:
        writer.close()

    # Write .pvar (TSV, plink2 format with #CHROM header)
    pvar_df = pd.DataFrame(
        {"#CHROM": chroms, "POS": positions, "ID": ids, "REF": refs, "ALT": alts}
    )
    pvar_df.to_csv(out_pvar, sep="\t", index=False, lineterminator="\n")

    # Write .psam (TSV, plink2 format with #IID header)
    psam_df = pd.DataFrame(
        {
            "#IID": iids,
            "SEX": sex,
            "POP": pops,
            "PSEUDOHAPLOID": np.where(is_pseudohap, "1", "0"),
        }
    )
    psam_df.to_csv(out_psam, sep="\t", index=False, lineterminator="\n")

    return InputDescriptor(
        path=out_prefix,
        pgen_path=out_pgen,
        pvar_path=out_pvar,
        psam_path=out_psam,
        fmt=InputFormat.PFILE,
        n_samples=spec.n_samples,
        n_variants=n_variants,
        is_target=False,
        eigfile_tempdir=None,
    )
