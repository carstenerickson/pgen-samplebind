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
    _intersection_count,
    _key_universe,
    _unique_chr_pos_codes,
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
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")

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
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
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
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
    # Strip the classification from the first comparison to simulate a
    # caller that forgot to run classify_pair.
    from dataclasses import replace

    unclassified = replace(report.comparisons[0], classification=None, classification_evidence=None)
    broken_report = replace(report, comparisons=(unclassified,))
    with pytest.raises(InvariantViolation, match="un-classified"):
        evaluate_gate(broken_report, MergePolicy())


def test_compute_preflight_low_intersection_same_chroms(tmp_path: Path) -> None:
    """Different variant seeds on the same chromosomes → near-zero intersection
    AND no consistent per-chrom position shift → `disjoint_panels`.

    Before the build-shift sharpener (commit-TBD), this case was forced into
    the `build_mismatch` bucket purely on the symmetric-chrom signal. The
    sharpener splits build_mismatch (uniform per-chrom position shift, the
    real hg19/hg38 fingerprint) from disjoint-same-chroms (random positions
    on overlapping chroms). Random variant_seeds with no shift→ disjoint.
    """
    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=999)

    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")

    pair = report.comparisons[0]
    # Random positions in a 1..100M space across 50 variants per panel: collisions vanishingly rare.
    assert pair.intersection_fraction_of_min < 0.1
    # The sharpened classifier: random positions on same chroms → not a
    # consistent shift → disjoint, not build_mismatch.
    assert pair.classification == "disjoint_panels"
    # And the signature itself: computed, but flags no consistent shift.
    sig = pair.build_shift_signature
    assert sig is not None
    assert sig["has_consistent_shift"] is False
    # Evidence echoes the shift decision so users can audit it.
    assert pair.classification_evidence is not None
    assert pair.classification_evidence["build_shift_has_consistent_shift"] is False


def test_compute_preflight_single_input(tmp_path: Path) -> None:
    """Canonical only → comparisons is empty but envelope is well-formed."""
    a = _panel(tmp_path, "a", variant_seed=11)
    descriptors = [_descriptor(a)]
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
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
        "build_shift_signature",
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


def test_build_shift_sharpener_uniform_shift_keeps_build_mismatch_label(
    tmp_path: Path,
) -> None:
    """A uniform +1M per-chrom POS shift (the hg19/hg38 fingerprint) must
    still classify as `build_mismatch` after the sharpener. Catches a
    regression where tightening the threshold accidentally drops the
    label for the very case it's designed to flag.
    """
    import pandas as pd

    a = _panel(tmp_path, "a", variant_seed=11, n_variants=200)
    b_prefix = tmp_path / "b" / "shifted"
    b_prefix.parent.mkdir(parents=True, exist_ok=True)
    # Build b as a uniform-shifted copy of a's .pvar (keep .pgen and .psam
    # untouched — the classifier only reads pvar). +1.5M per chrom is
    # large enough to clear the magnitude floor and the +0/+1 noise.
    import shutil

    shutil.copy(Path(str(a) + ".pgen"), Path(str(b_prefix) + ".pgen"))
    shutil.copy(Path(str(a) + ".psam"), Path(str(b_prefix) + ".psam"))
    pvar = pd.read_csv(Path(str(a) + ".pvar"), sep="\t")
    pvar["POS"] = pvar["POS"] + 1_500_000
    pvar["ID"] = pvar.apply(lambda r: f"chr{r['#CHROM']}:{r['POS']}", axis=1)
    pvar.to_csv(Path(str(b_prefix) + ".pvar"), sep="\t", index=False, lineterminator="\n")

    descriptors = [_descriptor(a), _descriptor(b_prefix)]
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
    pair = report.comparisons[0]
    assert pair.classification == "build_mismatch"
    sig = pair.build_shift_signature
    assert sig is not None
    assert sig["has_consistent_shift"] is True
    # Every evaluated chrom should show the +1.5M shift cleanly.
    for chrom_str, median in sig["median_shift_per_chrom"].items():
        assert abs(median - 1_500_000) < 1, (
            f"chr{chrom_str}: expected uniform shift ~1.5M, got median={median}"
        )
        assert sig["relative_mad_per_chrom"][chrom_str] < 0.01, (
            f"chr{chrom_str}: shifts are uniform, relative MAD should be ~0"
        )


def test_placeholder_canonical_ids_chr_pos_classifies_compatible(tmp_path: Path) -> None:
    """Regression for review finding #C1: canonical .pvar with all-placeholder
    IDs ('.', 'NA', etc.) under the default --variant-key chr_pos must NOT
    be mis-labeled as `empty_input`. The chr_pos data is intact — placeholder
    IDs only affect the alternate-key view, which the empty_input guard
    must not key on.
    """
    import pandas as pd

    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=11)  # same seed → identical chr_pos

    # Strip canonical's IDs to all-'.' (the most-common real-world placeholder).
    pvar_path = Path(str(a) + ".pvar")
    df = pd.read_csv(pvar_path, sep="\t")
    df["ID"] = "."
    df.to_csv(pvar_path, sep="\t", index=False, lineterminator="\n")

    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
    pair = report.comparisons[0]
    # Active-key (chr_pos) data is intact → compatible. Pre-fix this was
    # `empty_input` because alternate_key_canonical_size==0 tripped the guard.
    assert pair.classification == "compatible"
    assert pair.alternate_key_canonical_size == 0  # placeholders filtered
    # Evidence carries the active-key canonical size used by the new guard.
    assert pair.classification_evidence["canonical_active_size"] > 0


def test_empty_other_input_classifies_empty_input_and_is_gated(tmp_path: Path) -> None:
    """Regression for review finding #B6: an `other` input with zero post-
    filter variants must classify as `empty_input` AND must trip the gate
    under --preflight-policy strict. Pre-fix, empty_input was excluded
    from GATE_FAILURE_LABELS, so strict-mode users got no protection.
    """
    import pandas as pd

    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=22)

    # Wipe other's .pvar to header-only (zero variant rows). The biallelic-
    # SNP filter would catch a multi-allelic-only panel too; we use a
    # header-only file for test determinism.
    pvar_path = Path(str(b) + ".pvar")
    df = pd.read_csv(pvar_path, sep="\t")
    df.iloc[0:0].to_csv(pvar_path, sep="\t", index=False, lineterminator="\n")

    descriptors = [_descriptor(a), _descriptor(b)]
    strict_policy = MergePolicy(preflight_policy="strict")
    report = compute_preflight(descriptors, strict_policy, tool_version="test", command="merge")
    assert report.comparisons[0].classification == "empty_input"

    gated = evaluate_gate(report, strict_policy)
    assert gated.gate["would_trigger"] is True
    assert gated.gate["triggered"] is True
    assert gated.gate["action"] == "error"
    assert gated.gate["failing_inputs"][0]["classification"] == "empty_input"


def test_small_panel_build_mismatch_keeps_build_mismatch_label(tmp_path: Path) -> None:
    """Regression for review finding #B1: a small targeted panel (~3
    variants/chrom — below `_BUILD_SHIFT_MIN_VARIANTS_PER_CHROM = 5`)
    with symmetric chroms, zero coord overlap, and a hg19/hg38-style
    shift must still classify as `build_mismatch`. Pre-fix, the
    signature couldn't be computed (no qualifying chroms) and the
    classifier demoted to `disjoint_panels`, hiding the liftover hint.
    """
    import shutil

    import pandas as pd

    # Tiny panel: 3 chroms x 3 variants = 9 total. Each chrom has only
    # 3 variants, well under the 5-variant signature threshold.
    spec = SyntheticPanelSpec(
        n_samples=4,
        n_variants=9,
        n_populations=1,
        chromosomes=(1, 2, 3),
        variant_seed=701,
        sample_seed=702,
    )
    a = synthesize_pfile(spec, tmp_path / "small_a").path

    # Build b as a uniformly-shifted copy of a's .pvar (+1.5M, all chroms).
    b_prefix = tmp_path / "small_b"
    shutil.copy(Path(str(a) + ".pgen"), Path(str(b_prefix) + ".pgen"))
    shutil.copy(Path(str(a) + ".psam"), Path(str(b_prefix) + ".psam"))
    pvar = pd.read_csv(Path(str(a) + ".pvar"), sep="\t")
    pvar["POS"] = pvar["POS"] + 1_500_000
    pvar["ID"] = pvar.apply(lambda r: f"chr{r['#CHROM']}:{r['POS']}", axis=1)
    pvar.to_csv(Path(str(b_prefix) + ".pvar"), sep="\t", index=False, lineterminator="\n")

    descriptors = [_descriptor(a), _descriptor(b_prefix)]
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
    pair = report.comparisons[0]
    # Below-threshold panel: signature is None, but conservative fallback
    # preserves the build_mismatch label so the user sees the liftover hint.
    assert pair.classification == "build_mismatch"
    assert pair.build_shift_signature is None
    assert pair.classification_evidence["build_shift_signature_available"] is False


def test_disjoint_panels_hint_acknowledges_both_subshapes(tmp_path: Path) -> None:
    """Regression for review finding #B2: the disjoint_panels hint must
    acknowledge BOTH sub-shapes the label now covers (asymmetric chroms
    OR symmetric chroms without consistent shift). The old hint said
    'Chromosome coverage differs' which is wrong for the symmetric case.
    """
    from pgen_samplebind.preflight import _CLASSIFICATION_HINTS

    hint = _CLASSIFICATION_HINTS["disjoint_panels"]
    # Both sub-shape descriptions must appear so users in either case
    # see something true about their situation.
    assert "coverage differs" in hint.lower() or "different" in hint.lower()
    assert "coordinates don't" in hint or "coordinate" in hint
    # And the cross-reference to hash for panel identity should be there.
    assert "hash" in hint


def test_build_shift_signature_absent_when_no_qualifying_chroms(tmp_path: Path) -> None:
    """When the pair is compatible (high intersection), no chrom has
    zero-overlap, so the signature isn't computed — `None` rather than
    an empty dict. Pins the cost guarantee: identical-data pairs pay no
    shift-computation overhead."""
    a = _panel(tmp_path, "a", variant_seed=11)
    b = _panel(tmp_path, "b", variant_seed=11)
    descriptors = [_descriptor(a), _descriptor(b)]
    report = compute_preflight(descriptors, MergePolicy(), tool_version="test", command="merge")
    pair = report.comparisons[0]
    assert pair.classification == "compatible"
    assert pair.build_shift_signature is None


def test_chr_pos_code_encoding_matches_tuple_set_semantics() -> None:
    """The int64 (chrom, pos) encoding must be collision-free and produce
    exactly the same dedup + intersection counts as the prior
    `set(zip(chrom, pos))` representation it replaced (issue #12 perf
    follow-up — vectorizing the dominant cost must not change results).
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(7)
    chrom_a = rng.integers(1, 27, size=5000)
    pos_a = rng.integers(1, 250_000_000, size=5000)
    a = pd.DataFrame({"chrom": chrom_a, "pos": pos_a, "id": ["."] * 5000})
    # b shares the first 2000 rows of a, plus 3000 fresh ones.
    chrom_b = np.concatenate([chrom_a[:2000], rng.integers(1, 27, size=3000)])
    pos_b = np.concatenate([pos_a[:2000], rng.integers(1, 250_000_000, size=3000)])
    b = pd.DataFrame({"chrom": chrom_b, "pos": pos_b, "id": ["."] * 5000})

    # Reference: the old set-of-tuples semantics.
    set_a = set(zip(chrom_a.tolist(), pos_a.tolist(), strict=True))
    set_b = set(zip(chrom_b.tolist(), pos_b.tolist(), strict=True))

    codes_a = _unique_chr_pos_codes(a)
    codes_b = _unique_chr_pos_codes(b)
    assert codes_a.size == len(set_a)
    assert codes_b.size == len(set_b)
    assert _intersection_count(codes_a, codes_b) == len(set_a & set_b)


def test_chr_pos_code_encoding_rejects_out_of_range_pos() -> None:
    """A non-physical position (>= 2^40) would collide in the packed
    encoding, so the encoder raises rather than silently miscounting."""
    import pandas as pd

    bad = pd.DataFrame({"chrom": [1], "pos": [1 << 40], "id": ["."]})
    with pytest.raises(InvariantViolation, match="out of range"):
        _unique_chr_pos_codes(bad)


def test_key_universe_id_filters_placeholders() -> None:
    """The id universe stays set-based (faster than numpy for strings) and
    filters placeholder IDs — the vectorization only touched chr_pos."""
    import pandas as pd

    df = pd.DataFrame(
        {"chrom": [1, 1, 2, 2], "pos": [10, 20, 30, 40], "id": ["rs1", ".", "rs2", "NA"]}
    )
    universe = _key_universe(df, "id")
    assert universe == {"rs1", "rs2"}
