"""Shared dataclasses and enums.

The contract surface between modules. See LLD §2 for the design rationale on each type.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, IntEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    import numpy as np


# -----------------------------------------------------------------------------
# Enums (LLD §2.1)
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
    NON_BIALLELIC = "non_biallelic"
    NON_SNP = "non_snp"
    PRE_ALIGNMENT_OTHER = "pre_alignment_other"


class PseudohaploidStatus(Enum):
    DIPLOID = "0"
    PSEUDOHAPLOID = "1"
    UNKNOWN = "U"


class ExitCode(IntEnum):
    """Stable across versions per HLD §Exit codes."""

    OK = 0
    VALIDATION_FAILURE = 1
    IO_FAILURE = 2
    INVARIANT_VIOLATION = 3
    USAGE_ERROR = 4


# -----------------------------------------------------------------------------
# Policy literal types (LLD §2.4)
# -----------------------------------------------------------------------------

OnMismatch = Literal["drop", "error"]
OnMissing = Literal["fill_missing", "drop_variant", "error"]
OnExtra = Literal["warn", "drop", "error"]
OnStrand = Literal["drop", "flip", "error"]
OnCollision = Literal["error", "first", "suffix"]
VariantKey = Literal["chr_pos", "id"]


# -----------------------------------------------------------------------------
# Input descriptor (LLD §2.3)
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
# Policy (LLD §2.4)
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
    id_column: str = "IID"
    threads: int = 1
    block_size: int = 2048
    extras_warn_threshold: float = 0.10
    validate_strand_fail_pct: float = 10.0
    report_json_include_rows: bool = False


# -----------------------------------------------------------------------------
# Variant + sample views (LLD §2.5, §2.6)
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VariantRow:
    """Typed handle on a single canonical variant. Bulk pass-1 representation
    is a pandas DataFrame; this dataclass wraps a single row when a typed
    handle is needed."""

    canonical_idx: int
    chrom: int
    pos: int
    variant_id: str
    ref: str
    alt: str


@dataclass(slots=True)
class SampleRecord:
    """Logical view of one .psam row."""

    iid: str
    fid: str
    sex: int  # plink convention: 0=unknown, 1=male, 2=female
    pop: str
    pseudohaploid: PseudohaploidStatus = PseudohaploidStatus.UNKNOWN
    extra_cols: dict[str, str] = field(default_factory=dict)
    source_input_idx: int = -1


# -----------------------------------------------------------------------------
# Run summary (LLD §2.7)
# -----------------------------------------------------------------------------


@dataclass(slots=True)
class PerInputSummary:
    descriptor: InputDescriptor
    n_samples: int
    n_variants_pre_filter: int
    n_variants_post_filter: int
    n_populations: int
    largest_population: tuple[str, int]
    pseudohaploid_counts: dict[PseudohaploidStatus, int] = field(default_factory=dict)
    sex_counts: dict[int, int] = field(default_factory=dict)


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


@dataclass(slots=True)
class RunSummary:
    """Top-level run summary; serialized to JSON by --report-json."""

    inputs: list[PerInputSummary] = field(default_factory=list)
    alignment: AlignmentSummary = field(default_factory=AlignmentSummary)
    output_n_samples: int = 0
    output_n_variants: int = 0
    output_populations: dict[str, int] = field(default_factory=dict)
    output_pseudohaploid: dict[PseudohaploidStatus, int] = field(default_factory=dict)
    iid_collisions: list[str] = field(default_factory=list)
    psam_column_conflicts: list[tuple[str, int]] = field(default_factory=list)
    output_paths: dict[str, Path] = field(default_factory=dict)
    elapsed_s: float = 0.0
    pgen_samplebind_version: str = ""
    plink2_version: str | None = None


# -----------------------------------------------------------------------------
# Per-variant report row (LLD §2.8)
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ReportRow:
    """One row of the --report TSV (HLD §Reports)."""

    variant_id: str
    chrom: int
    pos: int
    input_index: int
    action: MergeAction
    reason: str = ""


# -----------------------------------------------------------------------------
# Variant hash output (LLD §2.9)
# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VariantHash:
    """Output of the hash subcommand. Wraps a SHA-256 with provenance."""

    sha256_hex: str
    n_variants_hashed: int
    n_variants_pre_filter: int
    canonical_form_bytes: int

    def render(self) -> str:
        """Default emission: `sha256:<hex>` per HLD §Variant hash step 9."""
        return f"sha256:{self.sha256_hex}"


# -----------------------------------------------------------------------------
# Orchestration types (LLD §2.10)
# -----------------------------------------------------------------------------


@dataclass(slots=True)
class MergeCounters:
    """Returned by merge.merge_inputs after pass 2; consumed by run_merge.

    Field set per HLD §Module orchestration. action_histogram keys per LLD
    §2.10 8-key mapping pin (passthrough, swap, flip, fill_missing,
    dropped_ambiguous_strand, dropped_allele_mismatch,
    pre_alignment_filter_dropped, drop).
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
