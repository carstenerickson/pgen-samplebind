"""HLD-named integration tests for correctness invariants.

Maps to LLD §5.3 / HLD §Validation strategy. Day 5 lands tests 1-6 + 20:
- Test 1: round-trip identity (via --on-collision first; suffix is Day 8)
- Test 2: three-way associativity
- Test 3: strand-flip recovery
- Test 4: allele-swap recovery
- Test 5: ambiguous-strand handling (default drop vs --trust-strand pass)
- Test 6: missing-variant default (fill_missing vs drop_variant)
- Test 20: asymmetric extras handling (warn / drop / error)

Tests 7, 8, 13, 15, 16 (EIGENSTRAT) deferred to Day 6.
Tests 9 (pseudohaploid detection) deferred to Day 7.
Tests 11, 12 (hash invariance) already covered as unit tests in
test_hash_canonicalization.py.
Test 17 (mergeit f2 parity) requires AT2/mergeit binaries; nightly only.
Test 18 (perf benchmark) gets its own file; not on Day 5.
Tests 21, 22 already covered: 21 mostly via test_alignment_evaluate_gates,
22 via test_reports.
Test 23 (collision suffix) → Day 8 alongside --target.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib
import pytest
from click.testing import CliRunner

from pgen_samplebind.cli import cli
from tests.fixtures import modifiers
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile

# ---------- helpers ----------------------------------------------------------


def _read_pgen(prefix: Path, n_samples: int, n_variants: int) -> np.ndarray:
    pgen_path = Path(str(prefix) + ".pgen")
    reader = pgenlib.PgenReader(str(pgen_path).encode(), raw_sample_ct=n_samples)
    try:
        buf = np.empty((n_variants, n_samples), dtype=np.int8)
        reader.read_range(0, n_variants, buf, sample_maj=0)
        return buf
    finally:
        if hasattr(reader, "close"):
            reader.close()


def _read_psam_iids(prefix: Path) -> list[str]:
    df = pd.read_csv(Path(str(prefix) + ".psam"), sep="\t")
    iid_col = next(c for c in df.columns if c.lstrip("#") == "IID")
    return df[iid_col].tolist()


def _read_pvar_keys(prefix: Path) -> list[tuple[int, int]]:
    df = pd.read_csv(Path(str(prefix) + ".pvar"), sep="\t")
    chrom_col = next(c for c in df.columns if c.lstrip("#") == "CHROM")
    return list(zip(df[chrom_col].astype(int), df["POS"].astype(int), strict=True))


# ---------- HLD test 1: round-trip identity ----------------------------------


class TestHld01RoundTripIdentity:
    """Bind a panel with itself and verify nothing is lost.

    HLD test 1 specifies --on-collision suffix (Day 8). For Day 5 we use
    --on-collision first, which keeps the canonical samples and drops the
    duplicates — output equals input panel. Same code path exercised
    (full alignment + pass 2 + psam finalization)."""

    def test_self_merge_first_equals_input(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=10,
            n_variants=80,
            n_populations=2,
            variant_seed=11,
            sample_seed=12,
            sample_id_prefix="A",
        )
        a = synthesize_pfile(spec, tmp_path / "a")

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a.path),
                str(a.path),
                "-o",
                str(out),
                "--on-collision",
                "first",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        # IIDs unchanged (duplicates dropped → 10 samples kept)
        assert _read_psam_iids(out) == _read_psam_iids(a.path)
        # All variants kept (passthrough across the second input → no drops)
        assert _read_pvar_keys(out) == _read_pvar_keys(a.path)
        # Genotypes byte-identical
        original = _read_pgen(a.path, 10, 80)
        merged = _read_pgen(out, 10, 80)
        np.testing.assert_array_equal(merged, original)


# ---------- HLD test 2: three-way associativity ------------------------------


class TestHld02ThreeWayAssociativity:
    """merge(merge(A, B), C) byte-equal to merge(A, B, C) modulo sample-axis
    ordering. We verify both the variant ordering and the genotype matrix
    after sample-IID alignment."""

    def test_associative_for_three_panels(self, tmp_path: Path) -> None:
        common = dict(n_samples=8, n_variants=60, n_populations=2, variant_seed=20)
        a = synthesize_pfile(
            SyntheticPanelSpec(sample_seed=1, sample_id_prefix="A", **common),
            tmp_path / "a",
        )
        b = synthesize_pfile(
            SyntheticPanelSpec(sample_seed=2, sample_id_prefix="B", **common),
            tmp_path / "b",
        )
        c = synthesize_pfile(
            SyntheticPanelSpec(sample_seed=3, sample_id_prefix="C", **common),
            tmp_path / "c",
        )

        runner = CliRunner()

        # Path 1: merge(A, B), then merge with C
        ab = tmp_path / "ab"
        r1 = runner.invoke(cli, ["merge", str(a.path), str(b.path), "-o", str(ab), "--quiet"])
        assert r1.exit_code == 0, r1.output
        abc1 = tmp_path / "abc1"
        r2 = runner.invoke(cli, ["merge", str(ab), str(c.path), "-o", str(abc1), "--quiet"])
        assert r2.exit_code == 0, r2.output

        # Path 2: merge(A, B, C) directly
        abc2 = tmp_path / "abc2"
        r3 = runner.invoke(
            cli,
            ["merge", str(a.path), str(b.path), str(c.path), "-o", str(abc2), "--quiet"],
        )
        assert r3.exit_code == 0, r3.output

        # Variants identical and in same order
        assert _read_pvar_keys(abc1) == _read_pvar_keys(abc2)

        # Sample IIDs identical and in same order (canonical-first ordering
        # preserved through both paths since A is always input[0])
        assert _read_psam_iids(abc1) == _read_psam_iids(abc2)

        # Genotype matrices identical
        n_samples = len(_read_psam_iids(abc1))
        n_variants = len(_read_pvar_keys(abc1))
        m1 = _read_pgen(abc1, n_samples, n_variants)
        m2 = _read_pgen(abc2, n_samples, n_variants)
        np.testing.assert_array_equal(m1, m2)


# ---------- HLD test 3: strand-flip recovery ---------------------------------


class TestHld03StrandFlipRecovery:
    """A + flipped(A) merge: --on-strand flip (default) → passthrough-equivalent
    output (since strand-flip on biallelic SNPs is metadata-only on hardcalls
    per LLD §2.1 action-collapse)."""

    def test_strand_flip_default_recovers(self, tmp_path: Path) -> None:
        # Use unambiguous variants only by setting ambiguous_strand_fraction=0.0
        # so we can flip without falling into ambiguous-strand drops.
        spec = SyntheticPanelSpec(
            n_samples=8,
            n_variants=50,
            n_populations=2,
            variant_seed=30,
            sample_seed=31,
            sample_id_prefix="A",
            ambiguous_strand_fraction=0.0,
        )
        a = synthesize_pfile(spec, tmp_path / "a")

        # B: same data as A but with sample IDs prefixed differently AND strand
        # flipped at a subset of variants.
        spec_b = SyntheticPanelSpec(
            n_samples=8,
            n_variants=50,
            n_populations=2,
            variant_seed=30,
            sample_seed=32,
            sample_id_prefix="B",
            ambiguous_strand_fraction=0.0,
        )
        b_unflipped = synthesize_pfile(spec_b, tmp_path / "b_unflipped")

        # Flip strand at first 20 variants (all unambiguous since fraction=0.0)
        flip_indices = np.arange(20)
        b = tmp_path / "b"
        modifiers.flip_strand(b_unflipped.path, flip_indices, b)

        # Merge with default policy
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a.path),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(tmp_path / "r.json"),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        import json

        payload = json.loads((tmp_path / "r.json").read_text())
        hist = payload["alignment"]["action_histogram"]
        # 20 variants strand-flipped (input[1] complement matches input[0])
        assert hist["flip"] == 20
        # 30 passthrough (the remaining unflipped variants)
        assert hist["passthrough"] == 30
        # No drops
        assert hist["dropped_ambiguous_strand"] == 0
        assert hist["dropped_allele_mismatch"] == 0


# ---------- HLD test 4: allele-swap recovery ---------------------------------


class TestHld04AlleleSwapRecovery:
    """A + swap_ref_alt(A) merge: default policy auto-recodes the swapped
    genotypes via REF_ALT_SWAP. Output should match A's encoding for the
    swapped variants."""

    def test_swap_default_recodes(self, tmp_path: Path) -> None:
        spec_a = SyntheticPanelSpec(
            n_samples=6,
            n_variants=40,
            n_populations=2,
            variant_seed=40,
            sample_seed=41,
            sample_id_prefix="A",
            ambiguous_strand_fraction=0.0,
        )
        a = synthesize_pfile(spec_a, tmp_path / "a")

        spec_b = SyntheticPanelSpec(
            n_samples=6,
            n_variants=40,
            n_populations=2,
            variant_seed=40,
            sample_seed=42,
            sample_id_prefix="B",
            ambiguous_strand_fraction=0.0,
        )
        b_unswapped = synthesize_pfile(spec_b, tmp_path / "b_unswapped")

        swap_indices = np.arange(10)
        b = tmp_path / "b"
        modifiers.swap_ref_alt(b_unswapped.path, 6, 40, swap_indices, b)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a.path),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(tmp_path / "r.json"),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        import json

        payload = json.loads((tmp_path / "r.json").read_text())
        hist = payload["alignment"]["action_histogram"]
        assert hist["swap"] == 10  # the swapped variants
        assert hist["passthrough"] == 30
        assert hist["dropped_allele_mismatch"] == 0

        # Verify B's samples in the output have their swapped variants recoded
        # back to A's encoding (i.e., merged[swapped_indices, B_samples] should
        # equal b_unswapped's original genotypes, not b's).
        original_b = _read_pgen(b_unswapped.path, 6, 40)
        merged = _read_pgen(out, 12, 40)
        # B samples are columns 6..11 in the output
        b_cols_in_output = merged[:, 6:12]
        np.testing.assert_array_equal(b_cols_in_output[:40], original_b)


# ---------- HLD test 5: ambiguous-strand handling ----------------------------


class TestHld05AmbiguousStrand:
    """Default policy drops A/T+C/G ambiguous swap pairs; --trust-strand
    passes them through (per HLD spec; biological footgun documented in
    resolve_alleles)."""

    @pytest.fixture
    def panels_with_ambig_swap(self, tmp_path: Path) -> tuple[Path, Path, int]:
        """Build A + B where B has REF/ALT swapped at A's ambiguous positions only."""
        spec_a = SyntheticPanelSpec(
            n_samples=6,
            n_variants=100,
            n_populations=2,
            variant_seed=50,
            sample_seed=51,
            sample_id_prefix="A",
            ambiguous_strand_fraction=0.20,  # ~20% ambig for the test
        )
        a = synthesize_pfile(spec_a, tmp_path / "a")

        spec_b = SyntheticPanelSpec(
            n_samples=6,
            n_variants=100,
            n_populations=2,
            variant_seed=50,
            sample_seed=52,
            sample_id_prefix="B",
            ambiguous_strand_fraction=0.20,
        )
        b_unmodified = synthesize_pfile(spec_b, tmp_path / "b_unmodified")

        ambig_indices = modifiers.find_ambiguous_variant_indices(b_unmodified.path)
        b = tmp_path / "b"
        modifiers.swap_ref_alt(b_unmodified.path, 6, 100, ambig_indices, b)
        return a.path, b, len(ambig_indices)

    def test_default_drops_ambiguous(
        self, panels_with_ambig_swap: tuple[Path, Path, int], tmp_path: Path
    ) -> None:
        a, b, n_ambig = panels_with_ambig_swap
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(tmp_path / "r.json"),
                "--validate-strand-fail-pct",
                "100",  # disable gate (b)
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        import json

        payload = json.loads((tmp_path / "r.json").read_text())
        hist = payload["alignment"]["action_histogram"]
        assert hist["dropped_ambiguous_strand"] == n_ambig
        # Output variant count = canonical - ambiguous drops
        assert payload["output"]["n_variants"] == 100 - n_ambig

    def test_trust_strand_passes_through(
        self, panels_with_ambig_swap: tuple[Path, Path, int], tmp_path: Path
    ) -> None:
        a, b, _ = panels_with_ambig_swap
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a),
                str(b),
                "-o",
                str(out),
                "--trust-strand",
                "--report-json",
                str(tmp_path / "r.json"),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        import json

        payload = json.loads((tmp_path / "r.json").read_text())
        hist = payload["alignment"]["action_histogram"]
        # No ambiguous drops (trust-strand passes them as PASSTHROUGH per HLD)
        assert hist["dropped_ambiguous_strand"] == 0
        assert payload["output"]["n_variants"] == 100


# ---------- HLD test 6: missing-variant default ------------------------------


class TestHld06MissingVariantDefault:
    """A + subset(A) merge: default fill_missing keeps all of A's variants;
    --on-missing drop_variant drops the missing-in-input[1] variants."""

    def test_default_fill_missing_keeps_all(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=6,
            n_variants=100,
            n_populations=2,
            variant_seed=60,
            sample_seed=61,
            sample_id_prefix="A",
        )
        a = synthesize_pfile(spec, tmp_path / "a")

        spec_b = SyntheticPanelSpec(
            n_samples=6,
            n_variants=100,
            n_populations=2,
            variant_seed=60,
            sample_seed=62,
            sample_id_prefix="B",
        )
        b_full = synthesize_pfile(spec_b, tmp_path / "b_full")

        # B drops the last 25 variants
        drop = np.arange(75, 100)
        b = tmp_path / "b"
        modifiers.drop_variants(b_full.path, 6, 100, drop, b)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a.path),
                str(b),
                "-o",
                str(out),
                "--report-json",
                str(tmp_path / "r.json"),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        import json

        payload = json.loads((tmp_path / "r.json").read_text())
        hist = payload["alignment"]["action_histogram"]
        assert hist["fill_missing"] == 25
        assert hist["passthrough"] == 75
        # All 100 canonical variants kept (default fill_missing)
        assert payload["output"]["n_variants"] == 100

        # B's samples should be -9 at the missing variants in the output
        merged = _read_pgen(out, 12, 100)
        # Output variant order matches canonical (a). Missing-in-B variants are
        # those with positions matching the dropped indices in b_full.
        # Here we used the same variant_seed for a and b_full so positions match;
        # the dropped indices map directly to canonical row indices for the same
        # positions. Sample columns 6..11 are B's samples.
        np.testing.assert_array_equal(merged[75:100, 6:12], np.full((25, 6), -9, dtype=np.int8))

    def test_drop_variant_policy_drops_them(self, tmp_path: Path) -> None:
        spec = SyntheticPanelSpec(
            n_samples=6,
            n_variants=100,
            n_populations=2,
            variant_seed=60,
            sample_seed=61,
            sample_id_prefix="A",
        )
        a = synthesize_pfile(spec, tmp_path / "a")

        spec_b = SyntheticPanelSpec(
            n_samples=6,
            n_variants=100,
            n_populations=2,
            variant_seed=60,
            sample_seed=62,
            sample_id_prefix="B",
        )
        b_full = synthesize_pfile(spec_b, tmp_path / "b_full")

        drop = np.arange(75, 100)
        b = tmp_path / "b"
        modifiers.drop_variants(b_full.path, 6, 100, drop, b)

        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge",
                str(a.path),
                str(b),
                "-o",
                str(out),
                "--on-missing",
                "drop_variant",
                "--report-json",
                str(tmp_path / "r.json"),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.output

        import json

        payload = json.loads((tmp_path / "r.json").read_text())
        hist = payload["alignment"]["action_histogram"]
        # The 25 missing variants now drop instead of fill_missing
        assert hist["drop"] == 25
        assert hist["fill_missing"] == 0
        assert payload["output"]["n_variants"] == 75


# ---------- HLD test 20: asymmetric extras handling --------------------------


class TestHld20AsymmetricExtras:
    """100-variant panel as input[0] + 1.15M-variant panel as input[1] (the
    "input order reversed" failure mode). Default --on-extra warn fires a
    stderr warning naming the count; --on-extra error exits 3; --on-extra
    drop produces no warning."""

    @pytest.fixture
    def small_canonical_large_other(self, tmp_path: Path) -> tuple[Path, Path, int]:
        """Build a small canonical (50 variants) and a large other (500
        variants superset). 450 extras → 900% of canonical, well above the
        10% default warn threshold."""
        spec_large = SyntheticPanelSpec(
            n_samples=4,
            n_variants=500,
            n_populations=2,
            variant_seed=70,
            sample_seed=72,
            sample_id_prefix="B",
        )
        b = synthesize_pfile(spec_large, tmp_path / "b")

        # Build small canonical = first 50 variants of the large panel
        spec_small_full = SyntheticPanelSpec(
            n_samples=4,
            n_variants=500,
            n_populations=2,
            variant_seed=70,
            sample_seed=71,
            sample_id_prefix="A",
        )
        a_full = synthesize_pfile(spec_small_full, tmp_path / "a_full")
        keep = np.arange(50)
        a = tmp_path / "a"
        modifiers.subset_variants(a_full.path, 4, 500, keep, a)
        n_extras = 500 - 50  # 450
        return a, b.path, n_extras

    def test_default_warn_emits_stderr(
        self, small_canonical_large_other: tuple[Path, Path, int], tmp_path: Path
    ) -> None:
        a, b, _ = small_canonical_large_other
        out = tmp_path / "merged"
        runner = CliRunner()
        # The merge should fail at gate (a) since extras (450) >> 10% of canonical (5).
        # That's exit 1 (ValidationError), not exit 0 with a warning. Test that.
        result = runner.invoke(cli, ["merge", str(a), str(b), "-o", str(out), "--on-extra", "warn"])
        from pgen_samplebind.errors import ValidationError

        # CliRunner returns exception type; our gate (a) test elsewhere
        # confirms exit code 1 mapping. Here just confirm gate (a) fires.
        assert isinstance(result.exception, ValidationError)
        assert "gate (a)" in str(result.exception).lower()
        # Warning should ALSO have been written to stderr before the gate fired
        assert "WARNING" in result.output

    def test_on_extra_error_raises(
        self, small_canonical_large_other: tuple[Path, Path, int], tmp_path: Path
    ) -> None:
        a, b, _ = small_canonical_large_other
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["merge", str(a), str(b), "-o", str(out), "--on-extra", "error", "--quiet"]
        )
        from pgen_samplebind.errors import InvariantViolation

        assert isinstance(result.exception, InvariantViolation)
        assert "on-extra error" in str(result.exception).lower()

    def test_on_extra_drop_succeeds_silently(
        self, small_canonical_large_other: tuple[Path, Path, int], tmp_path: Path
    ) -> None:
        """--on-extra drop is the explicit "smaller canonical, extras intentional"
        mode: gate (a) is bypassed (its threshold is the *warn* threshold);
        merge succeeds with no warning."""
        a, b, _ = small_canonical_large_other
        out = tmp_path / "merged"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["merge", str(a), str(b), "-o", str(out), "--on-extra", "drop", "--quiet"]
        )
        assert result.exit_code == 0, result.output
        # No warning emitted (drop suppresses warn_extras_threshold via the
        # `on_extra == "warn"` check in alignment.evaluate_pass1_gates and the
        # warning is gated by the same policy semantics in merge_inputs).
        assert "WARNING" not in result.output
