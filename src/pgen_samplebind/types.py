"""Shared dataclasses and enums.

The contract surface between modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, IntEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    import numpy as np


# -----------------------------------------------------------------------------
# Enums
# -----------------------------------------------------------------------------


class InputFormat(Enum):
    PFILE = "pfile"
    BFILE = "bfile"
    EIGENSTRAT = "eigenstrat"


class MergeAction(Enum):
    """Per-(variant, input) action assigned in pass 1, applied in pass 2.

    Six values, but pass 2 collapses to two genotype operations:
      - identity (PASSTHROUGH, STRAND_FLIP)
      - swap, i.e. x -> 2-x with -9 preserved (REF_ALT_SWAP, STRAND_FLIP_AND_SWAP)
    DROP and FILL_MISSING are control-flow, not arithmetic.
    """

    PASSTHROUGH = "passthrough"
    REF_ALT_SWAP = "ref_alt_swap"
    STRAND_FLIP = "strand_flip"
    STRAND_FLIP_AND_SWAP = "strand_flip_and_swap"
    DROP = "drop"
    FILL_MISSING = "fill_missing"


class DropReason(Enum):
    """Why a variant was dropped in pass 1. Surfaced in the per-variant report."""

    AMBIGUOUS_STRAND = "ambiguous_strand"
    ALLELE_MISMATCH = "allele_mismatch"
    ON_MISSING_DROP_VARIANT = "on_missing_drop"
    ON_STRAND_DROP = "on_strand_drop"
    NON_BIALLELIC = "non_biallelic"
    NON_SNP = "non_snp"
    PRE_ALIGNMENT_OTHER = "pre_alignment_other"


class PseudohaploidStatus(Enum):
    DIPLOID = "0"
    PSEUDOHAPLOID = "1"
    UNKNOWN = "U"


class ExitCode(IntEnum):
    """Stable across versions."""

    OK = 0
    VALIDATION_FAILURE = 1
    IO_FAILURE = 2
    INVARIANT_VIOLATION = 3
    USAGE_ERROR = 4


# -----------------------------------------------------------------------------
# Policy literal types
# -----------------------------------------------------------------------------

OnMismatch = Literal["drop", "error"]
OnMissing = Literal["fill_missing", "drop_variant", "error"]
OnExtra = Literal["warn", "drop", "error"]
OnStrand = Literal["drop", "flip", "error"]
OnCollision = Literal["error", "first", "suffix"]
VariantKey = Literal["chr_pos", "id"]
PreflightPolicy = Literal["warn", "strict", "off"]


# -----------------------------------------------------------------------------
# Input descriptor
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class InputDescriptor:
    """Resolved metadata for one input."""

    path: Path
    pgen_path: Path
    pvar_path: Path
    psam_path: Path
    fmt: InputFormat
    n_samples: int = 0
    n_variants: int = 0
    is_target: bool = False
    eigfile_tempdir: Path | None = None


# -----------------------------------------------------------------------------
# Policy
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MergePolicy:
    """All user-tunable behavior flags for `merge` and `validate`."""

    on_mismatch: OnMismatch = "drop"
    on_missing: OnMissing = "fill_missing"
    on_extra: OnExtra = "warn"
    on_strand: OnStrand = "flip"
    on_collision: OnCollision = "error"
    trust_strand: bool = False
    variant_key: VariantKey = "chr_pos"
    target_min_call_rate: float = 0.40
    include_chrom: tuple[int, ...] = tuple(range(1, 23))
    population_column: str | None = None
    no_population_column: bool = False
    id_column: str = "IID"
    block_size: int = 2048
    extras_warn_threshold: float = 0.10
    validate_strand_fail_pct: float = 10.0
    report_json_include_rows: bool = False
    # Preflight gate: warn (default — stderr + continue), strict (raise
    # ValidationError on any non-compatible classification), or off
    # (compute + emit JSON but never warn or fail). See issue
    # [#12](https://github.com/carstenerickson/pgen-samplebind/issues/12) step 4.
    preflight_policy: PreflightPolicy = "warn"


# -----------------------------------------------------------------------------
# Run summary
# -----------------------------------------------------------------------------


@dataclass(slots=True)
class AlignmentSummary:
    """Pass-1 outcome counts."""

    n_passthrough: int = 0
    n_ref_alt_swap: int = 0
    n_strand_flip: int = 0
    n_fill_missing: int = 0
    n_dropped: int = 0
    n_dropped_by_reason: dict[DropReason, int] = field(default_factory=dict)
    n_extras_dropped: int = 0
    n_pre_alignment_filter_dropped: int = 0
    policy_error_triggers: dict[str, int] = field(default_factory=dict)


# -----------------------------------------------------------------------------
# Per-variant report row
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ReportRow:
    """One row of the --report TSV."""

    variant_id: str
    chrom: int
    pos: int
    input_index: int
    action: MergeAction
    reason: str = ""


# -----------------------------------------------------------------------------
# Variant hash output
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VariantHash:
    """Output of the hash subcommand. Wraps a SHA-256 with provenance."""

    sha256_hex: str
    n_variants_hashed: int
    n_variants_pre_filter: int
    canonical_form_bytes: int

    def render(self) -> str:
        """Default emission: `sha256:<hex>`."""
        return f"sha256:{self.sha256_hex}"


# -----------------------------------------------------------------------------
# Orchestration types
# -----------------------------------------------------------------------------


@dataclass(slots=True)
class MergeCounters:
    """Returned by merge.merge_inputs after pass 2; consumed by run_merge.

    action_histogram has a fixed 9-key schema — always all of:
    passthrough, swap, flip, fill_missing, dropped_ambiguous_strand,
    dropped_allele_mismatch, dropped_on_strand,
    pre_alignment_filter_dropped, drop. All keys are always present in
    the emitted JSON so workflow consumers can rely on the schema
    regardless of merge outcome — don't drop or rename a key without
    auditing the report-JSON consumers.
    """

    action_histogram: dict[str, int] = field(default_factory=dict)
    action_histogram_per_chrom: dict[int, dict[str, int]] = field(default_factory=dict)
    intersection_size: int = 0
    extras_count: int = 0
    per_sample_het: list[tuple[str, int, int]] = field(default_factory=list)
    n_output_samples: int = 0
    n_output_variants: int = 0
    variant_rows: list[ReportRow] | None = None


@dataclass(frozen=True, slots=True)
class SampleIdentityPlan:
    """Resolved sample identity precomputed by psam.resolve_sample_identity."""

    output_iids: tuple[str, ...]
    per_input_keep_mask: tuple[np.ndarray[Any, Any], ...]
    per_input_output_indices: tuple[np.ndarray[Any, Any], ...]


@dataclass(frozen=True)
class MergeContext:
    """Bundle passed to merge_inputs as the fourth `options` arg per HLD
    §Module orchestration."""

    policy: MergePolicy
    sample_plan: SampleIdentityPlan
    report_tsv_path: Path | None = None
    collect_variant_rows: bool = False
    # v0.2: tqdm-driven progress bar on the pass-2 variant-block loop.
    # Set by `run_merge` to `not quiet and sys.stderr.isatty()` so workflow
    # managers (piped stderr) and `--quiet` runs stay silent.
    show_progress: bool = False
    # Passed from the CLI --quiet flag so merge_inputs can suppress advisory
    # stderr warnings (e.g., extras-count threshold) consistently.
    quiet: bool = False
