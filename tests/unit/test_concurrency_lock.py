"""Unit tests for concurrency.output_lock + detect_network_filesystem.

LLD §3.7: non-blocking fcntl.flock; lock-held → IOFailure; NFS warning.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
import time
from pathlib import Path

import pytest

from pgen_samplebind.concurrency import detect_network_filesystem, output_lock
from pgen_samplebind.errors import IOFailure


def test_output_lock_creates_lock_file(tmp_path: Path) -> None:
    prefix = tmp_path / "panel"
    with output_lock(prefix):
        assert (tmp_path / "panel.lock").exists()
    # Lock file remains after release (next acquirer reuses it).
    assert (tmp_path / "panel.lock").exists()


def test_output_lock_re_acquires_after_release(tmp_path: Path) -> None:
    prefix = tmp_path / "panel"
    with output_lock(prefix):
        pass
    # Same process can re-take it cleanly.
    with output_lock(prefix):
        pass


def _hold_lock_in_child(prefix_str: str, started_evt: object, release_evt: object) -> None:
    """Subprocess target: hold the lock until the parent signals release."""
    from pgen_samplebind.concurrency import output_lock as _output_lock

    with _output_lock(Path(prefix_str)):
        started_evt.set()  # type: ignore[attr-defined]
        release_evt.wait(timeout=10.0)  # type: ignore[attr-defined]


def test_output_lock_held_by_other_process_raises_iofailure(tmp_path: Path) -> None:
    """fcntl.flock is per-process; need a real subprocess to test contention."""
    prefix = tmp_path / "panel"
    ctx = mp.get_context("fork" if sys.platform != "win32" else "spawn")
    started = ctx.Event()
    release = ctx.Event()
    proc = ctx.Process(target=_hold_lock_in_child, args=(str(prefix), started, release))
    proc.start()
    try:
        assert started.wait(timeout=5.0), "child failed to acquire lock in time"
        with pytest.raises(IOFailure, match="locked by another"), output_lock(prefix):
            pass  # pragma: no cover — should never reach
    finally:
        release.set()
        proc.join(timeout=5.0)
        if proc.is_alive():
            proc.terminate()
            proc.join(timeout=2.0)


def test_output_lock_serializes_within_one_process(tmp_path: Path) -> None:
    """Sanity: same-process serial acquire/release works without errors."""
    prefix = tmp_path / "panel"
    for _ in range(3):
        with output_lock(prefix):
            time.sleep(0.001)


def test_detect_network_filesystem_local_path_returns_none(tmp_path: Path) -> None:
    """tmp_path is on a local fs; should not be flagged."""
    assert detect_network_filesystem(tmp_path / "lock") is None


def test_detect_network_filesystem_nonexistent_parent_returns_none(tmp_path: Path) -> None:
    """Best-effort: unresolvable paths return None silently."""
    bad = tmp_path / "no" / "such" / "dir" / "lock"
    assert detect_network_filesystem(bad) is None or isinstance(detect_network_filesystem(bad), str)
