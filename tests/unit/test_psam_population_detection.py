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
