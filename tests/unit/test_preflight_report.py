"""Unit tests for preflight report (schema v1 contract).

Locks the JSON envelope and field types so downstream pipelines that
assert against `<prefix>.preflight.json` see a stable contract. The
classifier (step 3) and policy gate (step 4) will populate the currently-
placeholder fields without bumping schema_version — those follow-ups
must preserve every field this test pins.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pgen_samplebind.errors import InvariantViolation
from pgen_samplebind.preflight import (
    PREFLIGHT_SCHEMA_VERSION,
    PairCompatibility,
    PreflightReport,
    classify_pair,
    compute_preflight,
    evaluate_gate,
    write_preflight_json,
)
from pgen_samplebind.types import MergePolicy
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


def _panel(tmp_path: Path, name: str, *, variant_seed: int, n_variants: int = 50) -> Path:
    """One tiny synthetic panel (4 samples, configurable variants on chr1+chr2)."""
    spec = SyntheticPanelSpec(
        n_samples=4,
        n_variants=n_variants,
        n_populations=1,
        chromosomes=(1, 2),
        variant_seed=variant_seed,
        sample_seed=variant_seed + 1,
        sample_id_prefix=name.upper(),
    )
    return synthesize_pfile(spec, tmp_path / name).path


def _descriptor(prefix: Path):
    """Build an InputDescriptor by running the same format-detection path
    `run_merge` uses, so the test exercises the real read path."""
    from pgen_samplebind.formats import prepared_input

    with prepared_input(prefix, is_target=False, include_chrom=tuple(range(1, 23))) as desc:
        return desc


def test_compute_preflight_schema_v1_shape(tmp_path: Path) -> None:
    """Canonical + one other input → comparisons[0] populated; envelope fields fixed."""
    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=11)  # same variant_seed → identical variants

    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(
        descriptors, MergePolicy(), tool_version="test", command="merge"
    )

    assert report.schema_version == PREFLIGHT_SCHEMA_VERSION == 1
    assert report.tool == "pgen-samplebind"
    assert report.command == "merge"
    assert report.variant_key == "chr_pos"
    assert report.canonical["index"] == 0
    assert report.canonical["n_variants"] > 0
    assert len(report.comparisons) == 1
    pair = report.comparisons[0]
    assert pair.input_index == 1
    # Identical variant_seed → identical (chrom, pos) sets → full intersection.
    assert pair.intersection == pair.n_variants
    assert pair.intersection_fraction_of_min == pytest.approx(1.0)
    # Classifier (step 3) populates these. `gate` is the pre-evaluate
    # default; populated by `evaluate_gate` in step-4 tests below.
    assert pair.classification == "compatible"
    assert isinstance(pair.classification_evidence, dict)
    assert report.gate == {"triggered": False, "action": "none", "threshold": None}
    # Per-chrom breakdown covers both chromosomes used by the fixture.
    chroms_seen = {pc.chrom for pc in pair.per_chrom}
    assert chroms_seen == {1, 2}


def test_placeholder_ids_dont_inflate_alt_key_fraction(tmp_path: Path) -> None:
    """Regression for review finding #1: real-world .pvar files with '.' /
    empty / 'NA' IDs on every variant must not collapse the id-keyed set
    to a single element on both sides, which would trigger a spurious
    key_space_mismatch classification on a real chr_pos failure.
    """
    import pandas as pd

    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=999)

    # Rewrite both panels' ID columns to all-'.' placeholders to mimic a
    # production .pvar with no variant IDs assigned.
    for prefix in (a, b):
        pvar_path = Path(str(prefix) + ".pvar")
        df = pd.read_csv(pvar_path, sep="\t")
        df["ID"] = "."
        df.to_csv(pvar_path, sep="\t", index=False, lineterminator="\n")

    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(
        descriptors, MergePolicy(), tool_version="test", command="merge"
    )
    pair = report.comparisons[0]
    # Placeholder IDs filtered → alternate-key intersection is 0 (no
    # real IDs on either side), so the lift can't manufacture a
    # key_space_mismatch label. Without the fix, alt fraction would be
    # 1.0 (both sides collapse to {'.'}) and lift would exceed 0.4.
    assert pair.alternate_key_intersection == 0
    assert pair.classification != "key_space_mismatch"


def test_evaluate_gate_rejects_unclassified_comparisons(tmp_path: Path) -> None:
    """Regression for review finding #7: evaluate_gate must refuse to
    silently treat `classification=None` as not-failed."""
    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=11)
    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(
        descriptors, MergePolicy(), tool_version="test", command="merge"
    )
    # Strip the classification from the first comparison to simulate a
    # caller that forgot to run classify_pair.
    from dataclasses import replace

    unclassified = replace(
        report.comparisons[0], classification=None, classification_evidence=None
    )
    broken_report = replace(report, comparisons=(unclassified,))
    with pytest.raises(InvariantViolation, match="un-classified"):
        evaluate_gate(broken_report, MergePolicy())


def test_compute_preflight_low_intersection_same_chroms(tmp_path: Path) -> None:
    """Different variant seeds on the same chromosomes → near-zero intersection.

    The classifier labels this `build_mismatch` (not `disjoint_panels`) because
    chromosome coverage is symmetric on both sides — the classifier can't
    distinguish 'same panel, different build' from 'two unrelated panels on
    the same chromosomes' without signal it doesn't currently compute. The
    actionable advice is the same in both cases; see classify_pair's
    docstring for the caveat. The `disjoint_panels` label requires
    *asymmetric* chrom presence, which the corpus's disjoint pair exercises.
    """
    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=999)

    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(
        descriptors, MergePolicy(), tool_version="test", command="merge"
    )

    pair = report.comparisons[0]
    # Random positions in a 1..100M space across 50 variants per panel: collisions vanishingly rare.
    assert pair.intersection_fraction_of_min < 0.1
    # Pin the classifier label so any future refinement that distinguishes
    # the build-vs-disjoint cases here surfaces deliberately.
    assert pair.classification == "build_mismatch"


def test_compute_preflight_single_input(tmp_path: Path) -> None:
    """Canonical only → comparisons is empty but envelope is well-formed."""
    a = _panel(tmp_path, "a", variant_seed=11)
    descriptors = [_descriptor(a)]
    report = compute_preflight(
        descriptors, MergePolicy(), tool_version="test", command="merge"
    )
    assert report.comparisons == ()
    assert report.canonical["n_variants"] > 0


def test_write_preflight_json_roundtrip(tmp_path: Path) -> None:
    """Serialized JSON contains exactly the v1 envelope keys at top level.

    Mirrors the production write order: `compute_preflight` then
    `evaluate_gate` then `write_preflight_json`. Earlier versions of
    this test skipped `evaluate_gate` and pinned only the pre-gate
    3-key `gate` shape, missing the production 6-key shape — see
    review finding #3.
    """
    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=11)
    descriptors = [_descriptor(a), _descriptor(b)]
    policy = MergePolicy()
    report = compute_preflight(descriptors, policy, tool_version="test", command="merge")
    report = evaluate_gate(report, policy)  # production write order

    out = tmp_path / "preflight.json"
    write_preflight_json(report, out)
    payload = json.loads(out.read_text())

    assert set(payload.keys()) == {
        "schema_version",
        "tool",
        "tool_version",
        "command",
        "variant_key",
        "canonical",
        "comparisons",
        "gate",
    }
    assert payload["schema_version"] == 1
    # Comparison row shape — pin every key so step 3/4 additions are deliberate.
    # alternate_key_* added in step 3 (classifier); classification now populated.
    comp = payload["comparisons"][0]
    assert set(comp.keys()) == {
        "input_index",
        "path",
        "n_variants",
        "intersection",
        "intersection_fraction_of_min",
        "per_chrom",
        "alternate_key",
        "alternate_key_canonical_size",
        "alternate_key_other_size",
        "alternate_key_intersection",
        "alternate_key_fraction_of_min",
        "classification",
        "classification_evidence",
    }
    per_chrom_row = comp["per_chrom"][0]
    assert set(per_chrom_row.keys()) == {
        "chrom",
        "canonical_size",
        "other_size",
        "intersection",
    }
    # Pin the post-`evaluate_gate` gate-dict shape. Production always
    # runs evaluate_gate before write, so the 6 keys below are the
    # actual schema-v1 contract — not the pre-gate 3-key default that
    # earlier versions of this test mistakenly pinned.
    assert set(payload["gate"].keys()) == {
        "triggered",
        "would_trigger",
        "action",
        "policy",
        "threshold",
        "failing_inputs",
    }
