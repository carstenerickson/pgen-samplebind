"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

from pgen_samplebind.types import InputDescriptor

from .fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile


@pytest.fixture(scope="session")
def synth_panel_500x50k(tmp_path_factory: pytest.TempPathFactory) -> InputDescriptor:
    """Default synthetic panel; built once per session and reused."""
    out_dir = tmp_path_factory.mktemp("synth_500x50k")
    return synthesize_pfile(SyntheticPanelSpec(), out_dir / "panel")


@pytest.fixture(scope="session")
def synth_panel_tiny(tmp_path_factory: pytest.TempPathFactory) -> InputDescriptor:
    """Tiny synthetic panel for fast unit-style integration tests."""
    out_dir = tmp_path_factory.mktemp("synth_tiny")
    spec = SyntheticPanelSpec(n_samples=20, n_variants=200, n_populations=3)
    return synthesize_pfile(spec, out_dir / "tiny")


@pytest.fixture(scope="function")
def tmp_run_dir(tmp_path: Path) -> Path:
    """Per-function tempdir for output triplets and reports."""
    return tmp_path
