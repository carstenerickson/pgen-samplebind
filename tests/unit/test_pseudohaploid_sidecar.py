"""Unit tests for `pseudohaploid.read_sidecar` (issue #2).

The sidecar is an authoritative per-sample pseudohaploid declaration written
by sibling tools (e.g., pileup-aadr). See HLD §Pseudohaploid sidecar and the
schema v1 spec in pileup-aadr's LLD §output.py.

Coverage:
  - Missing sidecar → silent None (caller falls back to existing behavior).
  - Valid v1 → IID → PseudohaploidStatus mapping with 0/1 → DIPLOID/PSEUDOHAPLOID.
  - Suffix tolerance: `<prefix>` and `<prefix>.geno` both locate the same file.
  - Schema enforcement: missing/unknown version, missing `samples`, malformed
    per-sample entries, out-of-range `pseudohaploid` values, malformed JSON.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pgen_samplebind.errors import IOFailure, UsageError
from pgen_samplebind.pseudohaploid import read_sidecar
from pgen_samplebind.types import PseudohaploidStatus as P


def _write_sidecar(prefix: Path, payload: object) -> Path:
    path = Path(str(prefix) + ".pseudohaploid.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class TestReadSidecarMissing:
    def test_missing_sidecar_returns_none(self, tmp_path: Path) -> None:
        prefix = tmp_path / "nothing_here"
        assert read_sidecar(prefix) is None

    def test_unrelated_files_do_not_trigger(self, tmp_path: Path) -> None:
        # Other companion files alongside but no .pseudohaploid.json.
        prefix = tmp_path / "panel"
        Path(str(prefix) + ".geno").write_text("")
        Path(str(prefix) + ".ind").write_text("")
        Path(str(prefix) + ".snp").write_text("")
        assert read_sidecar(prefix) is None


class TestReadSidecarValid:
    def test_v1_binary_mapping(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(
            prefix,
            {
                "schema_version": 1,
                "samples": {
                    "sampleA": {"pseudohaploid": 1, "note": "single-BAM randomDiploid"},
                    "sampleB": {"pseudohaploid": 0},
                },
            },
        )
        result = read_sidecar(prefix)
        assert result == {"sampleA": P.PSEUDOHAPLOID, "sampleB": P.DIPLOID}

    def test_extra_metadata_fields_ignored(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(
            prefix,
            {
                "schema_version": 1,
                "tool": "pileup-aadr",
                "tool_version": "0.1.0",
                "samples": {
                    "carsten": {
                        "pseudohaploid": 1,
                        "het_count": 0,
                        "het_rate": 0.0,
                        "note": "...",
                    },
                },
            },
        )
        assert read_sidecar(prefix) == {"carsten": P.PSEUDOHAPLOID}

    def test_empty_samples_dict_returns_empty_map(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1, "samples": {}})
        assert read_sidecar(prefix) == {}


class TestReadSidecarSuffixTolerance:
    """Caller may pass either bare prefix or one with a format suffix."""

    def test_bare_prefix_locates_sidecar(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1, "samples": {"A": {"pseudohaploid": 1}}})
        assert read_sidecar(prefix) == {"A": P.PSEUDOHAPLOID}

    @pytest.mark.parametrize("suffix", [".geno", ".pgen", ".pvar", ".psam", ".bed"])
    def test_format_suffix_stripped(self, tmp_path: Path, suffix: str) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1, "samples": {"A": {"pseudohaploid": 0}}})
        # User passes panel.geno (or similar) — same sidecar is located.
        assert read_sidecar(Path(str(prefix) + suffix)) == {"A": P.DIPLOID}


class TestReadSidecarSchemaErrors:
    def test_missing_schema_version(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"samples": {"A": {"pseudohaploid": 1}}})
        with pytest.raises(UsageError, match="missing required `schema_version`"):
            read_sidecar(prefix)

    def test_unsupported_schema_version(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 99, "samples": {}})
        with pytest.raises(UsageError, match="unsupported schema_version=99"):
            read_sidecar(prefix)

    def test_missing_samples_field(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1})
        with pytest.raises(UsageError, match="missing or non-object `samples`"):
            read_sidecar(prefix)

    def test_samples_must_be_object(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1, "samples": ["A", "B"]})
        with pytest.raises(UsageError, match="non-object `samples`"):
            read_sidecar(prefix)

    def test_per_sample_entry_must_be_object(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1, "samples": {"A": 1}})
        with pytest.raises(UsageError, match=r"samples\['A'\] must be an object"):
            read_sidecar(prefix)

    def test_missing_pseudohaploid_field(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        _write_sidecar(prefix, {"schema_version": 1, "samples": {"A": {"note": "..."}}})
        with pytest.raises(UsageError, match="missing required `pseudohaploid`"):
            read_sidecar(prefix)

    @pytest.mark.parametrize("bad_value", [2, -1, "1", 1.5, None, True, False])
    def test_pseudohaploid_must_be_0_or_1(self, tmp_path: Path, bad_value: object) -> None:
        # bool slips through the `value == 1` test in Python (True == 1),
        # so confirm we reject True/False explicitly via the parametrization.
        prefix = tmp_path / "panel"
        _write_sidecar(
            prefix,
            {"schema_version": 1, "samples": {"A": {"pseudohaploid": bad_value}}},
        )
        with pytest.raises(UsageError, match="must be 0 or 1"):
            read_sidecar(prefix)

    def test_top_level_not_object(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        path = Path(str(prefix) + ".pseudohaploid.json")
        path.write_text(json.dumps(["not", "an", "object"]))
        with pytest.raises(UsageError, match="top level must be a JSON object"):
            read_sidecar(prefix)


class TestReadSidecarMalformed:
    def test_malformed_json(self, tmp_path: Path) -> None:
        prefix = tmp_path / "panel"
        path = Path(str(prefix) + ".pseudohaploid.json")
        path.write_text("{this is not: valid JSON")
        with pytest.raises(UsageError, match="not valid JSON"):
            read_sidecar(prefix)

    def test_unreadable_file_raises_iofailure(self, tmp_path: Path) -> None:
        # Simulate read error by writing an empty file then chmod 000.
        prefix = tmp_path / "panel"
        path = Path(str(prefix) + ".pseudohaploid.json")
        path.write_text("{}")
        path.chmod(0o000)
        try:
            with pytest.raises((IOFailure, UsageError)):
                # On macOS / Linux, chmod 000 yields PermissionError → IOFailure.
                # If the test runs as root (CI containers occasionally do), the
                # read succeeds and we get UsageError (missing schema_version).
                read_sidecar(prefix)
        finally:
            path.chmod(0o644)  # restore so pytest cleanup works
