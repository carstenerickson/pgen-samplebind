"""Per-sample heterozygosity tally during pass 2; classification at end.

Per LLD §3.8. Cutoffs match HLD §Pseudohaploid detection:
  het_count == 0 over non-missing autosomal calls → PSEUDOHAPLOID
  het_rate >= 5% → DIPLOID
  0 < het_rate < 5% → UNKNOWN
  called_count == 0 → UNKNOWN (boundary; review fix #3)

`update_block` trusts the caller-supplies-chromosome contract: blocks must
not span chromosome boundaries (merge.merge_inputs splits at boundaries).

Sidecar reader (issue #2): `read_sidecar` parses an upstream tool's
`<prefix>.pseudohaploid.json` file when present. Sibling tools (e.g.,
pileup-aadr) know per-sample pseudohaploid status by construction;
honoring the sidecar lets that authoritative signal flow through to
the output `.psam` PSEUDOHAPLOID column instead of being re-inferred
from heterozygosity counts.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .errors import IOFailure, UsageError
from .types import PseudohaploidStatus

_DIPLOID_HET_RATE_THRESHOLD = 0.05  # ≥ 5% het rate → diploid

# Sidecar schema versions this reader accepts. Bump when the on-disk
# schema gains incompatible required fields.
_SIDECAR_SUPPORTED_SCHEMA_VERSIONS = frozenset({1})
_SIDECAR_SUFFIX = ".pseudohaploid.json"


def init_counters(n_samples: int) -> tuple[np.ndarray[Any, Any], np.ndarray[Any, Any]]:
    """(het_counts, called_counts) zero-init int64 arrays of shape (n_samples,)."""
    return (
        np.zeros(n_samples, dtype=np.int64),
        np.zeros(n_samples, dtype=np.int64),
    )


def update_block(
    block: np.ndarray[Any, Any],
    chrom: int,
    het_counts: np.ndarray[Any, Any],
    called_counts: np.ndarray[Any, Any],
) -> None:
    """Vectorized per-block tally update; in-place on the running counters.

    No-op for non-autosomal blocks (chrom > 22). The caller (merge.merge_inputs)
    guarantees blocks don't span chromosomes — see §3.10 chromosome-boundary pin.

    block shape: (block_variants, n_samples), int8.
    """
    if chrom > 22:
        return
    het_counts += (block == 1).sum(axis=0).astype(np.int64)
    called_counts += (block != -9).sum(axis=0).astype(np.int64)


def classify(het_count: int, called_count: int) -> PseudohaploidStatus:
    """Per-sample classification per HLD §Pseudohaploid detection.

    Boundary: called_count == 0 → UNKNOWN (no signal; honest answer).
    """
    if called_count == 0:
        return PseudohaploidStatus.UNKNOWN
    if het_count == 0:
        return PseudohaploidStatus.PSEUDOHAPLOID
    rate = het_count / called_count
    if rate >= _DIPLOID_HET_RATE_THRESHOLD:
        return PseudohaploidStatus.DIPLOID
    return PseudohaploidStatus.UNKNOWN


def classify_all(
    het_counts: np.ndarray[Any, Any], called_counts: np.ndarray[Any, Any]
) -> np.ndarray[Any, Any]:
    """Vectorized classification across all samples; returns object-dtype array."""
    return np.array(
        [
            classify(int(h), int(c))
            for h, c in zip(het_counts.tolist(), called_counts.tolist(), strict=True)
        ],
        dtype=object,
    )


def read_sidecar(prefix: Path) -> dict[str, PseudohaploidStatus] | None:
    """Read `<prefix>.pseudohaploid.json` if present; otherwise return None.

    The sidecar is written by upstream tools (e.g., pileup-aadr) to assert
    per-sample pseudohaploid status by construction, taking precedence over
    heterozygosity-derived inference. Missing sidecar is silent (not an
    error) — callers fall back to existing behavior.

    Accepts either bare prefix (`/data/sample`) or one carrying a recognized
    format suffix (`/data/sample.geno`); the suffix is stripped to locate
    the sidecar at the canonical base.

    Schema v1 (per pileup-aadr LLD §output.py):
        {
          "schema_version": 1,
          "samples": {
            "<iid>": {"pseudohaploid": 0 | 1, ...optional fields...}
          }
        }

    Returns:
        dict mapping IID → PseudohaploidStatus (DIPLOID or PSEUDOHAPLOID;
        UNKNOWN never appears since the sidecar schema is binary).
        None when no sidecar file exists.

    Raises:
        IOFailure: sidecar file exists but is unreadable.
        UsageError: JSON parse failure; missing/unsupported schema_version;
            missing `samples`; per-sample entry missing the `pseudohaploid`
            field or carrying a value outside {0, 1}.
    """
    # Late import to avoid a circular dependency (formats imports nothing
    # from pseudohaploid, but keeping the import local is cheaper than
    # restructuring).
    from .formats import strip_known_suffix

    base = strip_known_suffix(prefix)
    sidecar_path = Path(str(base) + _SIDECAR_SUFFIX)
    if not sidecar_path.exists():
        return None

    try:
        raw_text = sidecar_path.read_text(encoding="utf-8")
    except OSError as e:
        raise IOFailure(f"cannot read pseudohaploid sidecar {sidecar_path}: {e}") from e

    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError as e:
        raise UsageError(
            f"pseudohaploid sidecar {sidecar_path} is not valid JSON: {e.msg} "
            f"(line {e.lineno}, col {e.colno})"
        ) from e

    if not isinstance(data, dict):
        raise UsageError(
            f"pseudohaploid sidecar {sidecar_path}: top level must be a JSON object, "
            f"got {type(data).__name__}"
        )

    version = data.get("schema_version")
    if version is None:
        raise UsageError(
            f"pseudohaploid sidecar {sidecar_path}: missing required `schema_version` field. "
            f"Supported versions: {sorted(_SIDECAR_SUPPORTED_SCHEMA_VERSIONS)}."
        )
    if version not in _SIDECAR_SUPPORTED_SCHEMA_VERSIONS:
        raise UsageError(
            f"pseudohaploid sidecar {sidecar_path}: unsupported schema_version={version!r}. "
            f"This pgen-samplebind build accepts: {sorted(_SIDECAR_SUPPORTED_SCHEMA_VERSIONS)}."
        )

    samples = data.get("samples")
    if not isinstance(samples, dict):
        raise UsageError(
            f"pseudohaploid sidecar {sidecar_path}: missing or non-object `samples` field "
            f"(got {type(samples).__name__})."
        )

    out: dict[str, PseudohaploidStatus] = {}
    for iid, entry in samples.items():
        if not isinstance(entry, dict):
            raise UsageError(
                f"pseudohaploid sidecar {sidecar_path}: samples[{iid!r}] must be an object, "
                f"got {type(entry).__name__}"
            )
        if "pseudohaploid" not in entry:
            raise UsageError(
                f"pseudohaploid sidecar {sidecar_path}: samples[{iid!r}] missing required "
                f"`pseudohaploid` field."
            )
        value = entry["pseudohaploid"]
        # Reject booleans and floats explicitly: `True == 1` and `1.0 == 1` are
        # both True in Python, but the schema is strict {0, 1} integers.
        if not isinstance(value, int) or isinstance(value, bool) or value not in (0, 1):
            raise UsageError(
                f"pseudohaploid sidecar {sidecar_path}: samples[{iid!r}].pseudohaploid must "
                f"be 0 or 1, got {value!r}."
            )
        out[str(iid)] = (
            PseudohaploidStatus.PSEUDOHAPLOID if value == 1 else PseudohaploidStatus.DIPLOID
        )

    return out
