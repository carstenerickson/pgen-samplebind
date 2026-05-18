"""Canonical variant-set hash for cross-format panel-identity verification.

The canonicalization spec is the contract: two formats representing the same
panel must produce identical hashes.
"""

from __future__ import annotations

import hashlib
from io import StringIO

import pandas as pd

from .types import InputDescriptor, VariantHash


def canonicalize_pvar_bytes(df: pd.DataFrame) -> bytes:
    """Apply the canonical-variant-hash spec to a pre-filtered .pvar
    DataFrame (output of pvar.read_pvar — biallelic-SNP-filtered, chrom
    normalized to int, REF/ALT uppercased).

    Spec steps applied here (1-4 happen upstream in read_pvar):
      6. Sort by (chrom_int, pos_int, ref, alt) ascending — numeric on the
         first two, bytewise on the last two. The ref/alt tiebreak makes
         the canonical order well-defined even when the upstream caller
         did not enforce key uniqueness (e.g., a tri-allelic site stored
         as two biallelic rows). Without it, stable sort preserves
         original row order, and two formats that load the same variants
         in different physical orders produce different hashes.
      7. One line per variant: `{chrom}\\t{pos}\\t{ref}\\t{alt}\\n`.
      Returns the UTF-8 bytestream that step 8 will SHA-256.

    Variant ID is NOT part of the hash (step 5) — different sources name the
    same variant differently; the hash should be invariant.
    """
    df_sorted = df.sort_values(["chrom", "pos", "ref", "alt"], kind="mergesort").reset_index(
        drop=True
    )
    # Use StringIO + to_csv for speed (~5x faster than Python-level join at 1240k scale)
    buf = StringIO()
    df_sorted[["chrom", "pos", "ref", "alt"]].to_csv(
        buf, sep="\t", header=False, index=False, lineterminator="\n"
    )
    return buf.getvalue().encode("utf-8")


def hash_input(desc: InputDescriptor) -> VariantHash:
    """Read .pvar, canonicalize, SHA-256, return wrapped VariantHash."""
    h, _ = hash_input_with_canonical(desc)
    return h


def hash_input_with_canonical(desc: InputDescriptor) -> tuple[VariantHash, bytes]:
    """As hash_input plus the canonical bytestream for --emit-canonical."""
    from . import pvar as pvar_module  # avoid circular at module load

    n_pre_filter = pvar_module.count_raw_variants(desc.pvar_path)
    pvar_df = pvar_module.read_pvar(desc.pvar_path)
    canonical = canonicalize_pvar_bytes(pvar_df)
    sha = hashlib.sha256(canonical).hexdigest()
    return (
        VariantHash(
            sha256_hex=sha,
            n_variants_hashed=len(pvar_df),
            n_variants_pre_filter=n_pre_filter,
            canonical_form_bytes=len(canonical),
        ),
        canonical,
    )
