"""HLD test 18: throughput benchmark gate.

Measures wallclock + throughput (genotypes/sec) + peak RSS for a fixed-size
synthetic merge. Fails if throughput regresses below `threshold_pct_of_baseline`
of the recorded baseline (default 80% of HLD's 100 M genotypes/sec target).

Per LLD §5.7: Linux x86_64 only — macOS GitHub runners have inconsistent CPU
profiles that produce flaky regressions. Marked `slow`; runs in the dedicated
bench cell of the CI matrix (linux + py3.12 + --runslow).

Baseline lives in tests/integration/perf_baseline.json and updates manually
(PR with bench rationale) when an intentional perf change lands.
"""

from __future__ import annotations

import json
import platform
import resource
import sys
import time
from pathlib import Path

import pytest

from pgen_samplebind.commands.merge_cmd import run_merge
from pgen_samplebind.types import MergePolicy
from tests.fixtures.synthesize import SyntheticPanelSpec, synthesize_pfile

pytestmark = [pytest.mark.slow]

BASELINE_PATH = Path(__file__).parent / "perf_baseline.json"


def _is_linux_x86_64() -> bool:
    return sys.platform == "linux" and platform.machine() in {"x86_64", "AMD64"}


def _peak_rss_mb() -> float:
    """Peak resident set size in MB. ru_maxrss is bytes on macOS, KB on Linux."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def test_perf_benchmark(tmp_path: Path) -> None:
    """Measure throughput on a fixed two-panel merge; gate against baseline."""
    if not _is_linux_x86_64():
        pytest.skip(
            f"perf bench gate is linux-x86_64 only "
            f"(detected: {sys.platform}/{platform.machine()}); "
            "macOS CPU profiles produce flaky regressions per LLD §5.7"
        )

    baseline = json.loads(BASELINE_PATH.read_text())
    fixture_spec = baseline["fixture"]
    n_samples = fixture_spec["n_samples_per_panel"]
    n_variants = fixture_spec["n_variants"]
    assert fixture_spec["n_panels"] == 2, (
        "this benchmark fixture is hard-wired to two panels; if you change "
        "the baseline n_panels, update the test fixture too"
    )

    spec_a = SyntheticPanelSpec(
        n_samples=n_samples,
        n_variants=n_variants,
        n_populations=10,
        seed=0xBEC0DE,
        sample_id_prefix="A_",
    )
    spec_b = SyntheticPanelSpec(
        n_samples=n_samples,
        n_variants=n_variants,
        n_populations=10,
        seed=0xBEC0DE,  # same variants
        sample_id_prefix="B_",
        sample_seed=0xB,
    )
    desc_a = synthesize_pfile(spec_a, tmp_path / "A")
    desc_b = synthesize_pfile(spec_b, tmp_path / "B")

    pfile_a = Path(str(desc_a.pgen_path)[:-5])
    pfile_b = Path(str(desc_b.pgen_path)[:-5])
    out = tmp_path / "merged"

    policy = MergePolicy()

    started = time.perf_counter()
    run_merge(
        input_paths=(pfile_a, pfile_b),
        target_path=None,
        output_prefix=out,
        policy=policy,
        report_path=None,
        report_json_path=None,
        quiet=True,
    )
    elapsed_s = time.perf_counter() - started

    total_genotypes = 2 * n_samples * n_variants  # both panels combined
    throughput = total_genotypes / elapsed_s
    peak_rss_mb = _peak_rss_mb()

    baseline_throughput = baseline["throughput_genotypes_per_sec_baseline"]
    threshold = baseline["threshold_pct_of_baseline"] * baseline_throughput

    print(
        f"\n[perf] elapsed={elapsed_s:.2f}s, throughput={throughput:,.0f} g/s, "
        f"peak_rss={peak_rss_mb:.0f}MB, "
        f"baseline={baseline_throughput:,.0f} g/s, threshold={threshold:,.0f} g/s",
        flush=True,
    )

    assert elapsed_s < baseline["wallclock_s_ceiling"], (
        f"wallclock {elapsed_s:.1f}s exceeds ceiling "
        f"{baseline['wallclock_s_ceiling']}s — likely a fixture-size mismatch"
    )
    assert peak_rss_mb < baseline["max_rss_mb_ceiling"], (
        f"peak RSS {peak_rss_mb:.0f}MB exceeds ceiling "
        f"{baseline['max_rss_mb_ceiling']}MB — likely a memory regression"
    )
    assert throughput >= threshold, (
        f"throughput regression: {throughput:,.0f} g/s < {threshold:,.0f} g/s "
        f"({100 * throughput / baseline_throughput:.0f}% of baseline; gate at "
        f"{100 * baseline['threshold_pct_of_baseline']:.0f}%). "
        "Update perf_baseline.json with a PR rationale if this is intentional."
    )
