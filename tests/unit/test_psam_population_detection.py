"""Unit tests for psam.detect_population_column / rename_to_pop / add_fid_from_pop.

Per HLD §Population labels: auto-detect tries POP, then PHENO, then PHENO1
in that order. plink2's --eigfile writes PHENO1 (after our PHENO1 → POP
rename in the orchestrator); modern hand-built psam files use POP;
legacy plink1.9-converted files use PHENO. Override with --population-column NAME.
"""

from __future__ import annotations

import pandas as pd
import pytest

from pgen_samplebind.errors import InvariantViolation, UsageError
from pgen_samplebind.psam import (
    add_fid_from_pop,
    detect_population_column,
    rename_to_pop,
    try_detect_population_column,
)


class TestDetectPopulationColumnFallbackOrder:
    """POP > PHENO > PHENO1 priority per HLD §Population labels."""

    def test_pop_takes_precedence_over_pheno(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"], "PHENO": ["y"]})
        assert detect_population_column(df, override=None) == "POP"

    def test_pop_takes_precedence_over_pheno1(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"], "PHENO1": ["y"]})
        assert detect_population_column(df, override=None) == "POP"

    def test_pheno_takes_precedence_over_pheno1(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "PHENO": ["x"], "PHENO1": ["y"]})
        assert detect_population_column(df, override=None) == "PHENO"

    def test_only_pop_present(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"]})
        assert detect_population_column(df, override=None) == "POP"

    def test_only_pheno_present(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "PHENO": ["x"]})
        assert detect_population_column(df, override=None) == "PHENO"

    def test_only_pheno1_present(self) -> None:
        """plink2 --eigfile --make-pgen writes PHENO1 — this is the post-
        EIGENSTRAT-conversion case."""
        df = pd.DataFrame({"IID": ["a"], "PHENO1": ["x"]})
        assert detect_population_column(df, override=None) == "PHENO1"


class TestDetectPopulationColumnOverride:
    def test_explicit_override(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"], "MyPop": ["custom"]})
        assert detect_population_column(df, override="MyPop") == "MyPop"

    def test_override_takes_precedence_over_fallback(self) -> None:
        """Even when POP is present, override wins."""
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"], "MyPop": ["custom"]})
        assert detect_population_column(df, override="MyPop") == "MyPop"

    def test_override_missing_raises_usage_error(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"]})
        with pytest.raises(UsageError, match="not present"):
            detect_population_column(df, override="DoesNotExist")


class TestDetectPopulationColumnNoMatch:
    def test_no_population_column_raises_invariant(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "SEX": ["1"]})
        with pytest.raises(InvariantViolation, match="no population column"):
            detect_population_column(df, override=None)


class TestTryDetectPopulationColumn:
    """Non-raising lookup used by `validate --no-population-column` (issue #3)."""

    def test_returns_pop_when_present(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"]})
        assert try_detect_population_column(df, override=None) == "POP"

    def test_returns_none_when_no_fallback(self) -> None:
        df = pd.DataFrame({"IID": ["GFX0442453"], "SEX": ["1"]})
        assert try_detect_population_column(df, override=None) is None

    def test_bad_override_still_raises(self) -> None:
        """A typoed --population-column is a user error, not a missing column."""
        df = pd.DataFrame({"IID": ["a"], "SEX": ["1"]})
        with pytest.raises(UsageError, match="not present"):
            try_detect_population_column(df, override="NoSuchColumn")


class TestRenameToPop:
    def test_rename_pheno1_to_pop(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "PHENO1": ["x"]})
        out = rename_to_pop(df, source_col="PHENO1")
        assert "POP" in out.columns
        assert "PHENO1" not in out.columns
        assert out["POP"].tolist() == ["x"]

    def test_rename_pop_is_noop(self) -> None:
        df = pd.DataFrame({"IID": ["a"], "POP": ["x"]})
        out = rename_to_pop(df, source_col="POP")
        assert "POP" in out.columns
        assert out["POP"].tolist() == ["x"]

    def test_rebind_with_existing_pop_equal_to_source_drops_duplicate(self) -> None:
        """Round-tripping pgen-samplebind's own output: prior merge wrote both
        FID and POP equal to the population label. Re-binding with
        --population-column FID would otherwise produce two POP columns
        after rename (issue #6).
        """
        df = pd.DataFrame(
            {
                "FID": ["pop_a", "pop_b"],
                "IID": ["s1", "s2"],
                "SEX": ["1", "2"],
                "POP": ["pop_a", "pop_b"],
                "PSEUDOHAPLOID": ["0", "1"],
            }
        )
        out = rename_to_pop(df, source_col="FID")
        assert list(out.columns).count("POP") == 1
        assert "FID" not in out.columns
        assert out["POP"].tolist() == ["pop_a", "pop_b"]

    def test_rebind_with_distinct_pop_raises(self) -> None:
        """A genuine ambiguity — input has two candidate population columns
        with different values — should raise rather than silently drop one.
        """
        df = pd.DataFrame(
            {
                "FID": ["fam_a", "fam_b"],
                "IID": ["s1", "s2"],
                "POP": ["pop_x", "pop_y"],
            }
        )
        with pytest.raises(UsageError, match="distinct 'POP' column"):
            rename_to_pop(df, source_col="FID")

    def test_duplicate_pop_columns_raises_invariant(self) -> None:
        """Guard against duplicate POP columns arriving at rename_to_pop.
        Not reachable via the normal read_psam path (pandas mangles
        duplicate headers), but the guard turns a confusing pandas error
        into a clear InvariantViolation if some other caller ever
        constructs such a DataFrame.
        """
        df = pd.DataFrame([["a", "x", "y"]], columns=["IID", "POP", "POP"])
        with pytest.raises(InvariantViolation, match="2 'POP' columns"):
            rename_to_pop(df, source_col="IID")

    def test_rebind_with_nan_in_pop_raises_with_clear_diagnostic(self) -> None:
        """NaN in either candidate column breaks the equality short-circuit
        (str(nan) compares equal across columns coincidentally). Surface
        the real cause instead of a misleading 'distinct POP' error.
        """
        df = pd.DataFrame(
            {
                "FID": ["pop_a", "pop_b"],
                "IID": ["s1", "s2"],
                "POP": ["pop_a", None],
            }
        )
        with pytest.raises(UsageError, match="missing values"):
            rename_to_pop(df, source_col="FID")


class TestAddFidFromPop:
    def test_fid_equals_pop(self) -> None:
        df = pd.DataFrame({"IID": ["a", "b"], "POP": ["pop_x", "pop_y"]})
        out = add_fid_from_pop(df)
        assert "FID" in out.columns
        assert out["FID"].tolist() == ["pop_x", "pop_y"]

    def test_fid_overwrites_existing(self) -> None:
        """When .fam-derived input has FID already, add_fid_from_pop
        overwrites with POP per HLD §Output PFILE (AT2 keys on FID)."""
        df = pd.DataFrame({"FID": ["old"], "IID": ["a"], "POP": ["new"]})
        out = add_fid_from_pop(df)
        assert out["FID"].tolist() == ["new"]

    def test_duplicate_pop_columns_raises_invariant(self) -> None:
        """Defense in depth: rename_to_pop now guarantees a single POP, but
        the same guard in add_fid_from_pop ensures a direct caller passing
        duplicates gets a clear error instead of pandas' "Cannot set a
        DataFrame with multiple columns to FID" (issue #6's symptom).
        """
        df = pd.DataFrame([["a", "x", "y"]], columns=["IID", "POP", "POP"])
        with pytest.raises(InvariantViolation, match="2 'POP' columns"):
            add_fid_from_pop(df)
