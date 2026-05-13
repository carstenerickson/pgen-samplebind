# Contributing to pgen-samplebind

The project is small and opinionated. Substantive scope changes (multi-allelic merges, dosage data, BFILE-only output, etc.) should start with a design discussion on a [GitHub issue](https://github.com/carstenerickson/pgen-samplebind/issues) before implementation — the scope is deliberately bounded to the sample-bind operation and we want to keep it that way.

For bug reports, small fixes, and additive improvements, a PR is the right next step.

## Dev setup

```bash
git clone https://github.com/carstenerickson/pgen-samplebind.git
cd pgen-samplebind
python3.11 -m venv .venv          # or 3.12 / 3.13 / 3.14
source .venv/bin/activate
pip install -e '.[dev]'
```

The `[dev]` extra pulls in `pytest`, `pytest-cov`, `ruff`, `mypy`, `pandas-stubs`, `build`, and `twine`.

EIGENSTRAT and BFILE input paths require **plink2 v2.0.0-a.7.1 or newer** on `PATH`. Pure-PFILE development needs no plink2. See the main [README §Install](README.md#install) for the install snippet.

## Running tests

```bash
# Default suite — fast, no external tools required (~60s)
pytest tests/ -m "not slow and not external_tool and not dogfood_full"

# Full suite — runs everything that's available locally. Tests gated on
# external tools (mergeit, AdmixTools, R) auto-skip if the tool isn't
# installed; you don't need to set them all up just to run the suite.
pytest tests/
```

### pytest markers

| Marker | Means |
|---|---|
| `slow` | > 30s wallclock; e.g., `tests/integration/test_mergeit_f2_parity.py`, `tests/integration/test_perf_benchmark.py` |
| `external_tool` | Requires AdmixTools 2 / mergeit binaries on PATH |
| `eigenstrat` | Requires plink2 with EIGENSTRAT support |
| `network` | Requires internet (AADR Dataverse fetch) |
| `dogfood` | AADR-derivative regression tests; no external deps for the default tier |
| `dogfood_plink2` | Dogfood subset needing plink2 on PATH |
| `dogfood_full` | Dogfood subset needing R + admixtools (full qpAdm shootout) |

### Dogfood — published-research-workflow regression test

The dogfood suite is the project's primary **trust artifact**: anyone can clone, run, and verify pgen-samplebind reproduces the established `mergeit + plink2 + AdmixTools 2` pipeline on a published-research-shape workload, without trusting maintainer claims.

A 44-sample × 50K-variant fixture (Patterson 7-source + 4 English target pops + 1 individual target, drawn from AADR v66 under fair-use for non-commercial scholarly verification — see [`tests/dogfood/README.md`](tests/dogfood/README.md)) flows through `pgen-samplebind merge` and is verified against a vendored mergeit-pipeline reference qpAdm TSV.

```bash
pytest tests/dogfood/         # runs all 6 tests; auto-skips tier-2/3 if tools missing
```

The dogfood tests gate via `pytest.skipif` against `HAS_PLINK2` / `HAS_ADMIXTOOLS` probes in `tests/dogfood/conftest.py`, not via the marker selector. So bare `pytest tests/dogfood/` runs everything that's available on your machine and silently skips the rest:

| Tier | Marker | Requires | Verifies |
|---|---|---|---|
| 1 | `dogfood` (all six) | nothing (pgen-samplebind only) | panel shape, cM preservation, FID=POP, PSEUDOHAPLOID column populated |
| 2 | `dogfood_plink2` | `plink2` on PATH | PFILE → BED conversion preserves cM end-to-end |
| 3 | `dogfood_full` | `plink2` + R + `admixtools` | extract_f2 + qpAdm shootout numerically matches the vendored mergeit reference (tolerance 1e-6 on weights, 1e-4 on p_tail) |

To force the full tier-3 shootout locally, install plink2 + R + admixtools first; CI's dedicated `dogfood` job sets up all three on every commit and runs the whole set. The `dogfood_plink2` / `dogfood_full` markers exist for CI status-check naming and intent documentation; they don't subset what `pytest tests/dogfood/` runs.

## Code quality

All three must pass before a PR can be merged. CI enforces them on every push:

```bash
ruff check src/ tests/
ruff format --check src/ tests/
mypy src/
```

- **`ruff`**: linting (E/F/W/I/B/UP/SIM/RUF rule sets per `pyproject.toml`) and formatting. Run `ruff format src/ tests/` to auto-fix formatting.
- **`mypy`**: strict mode over `src/pgen_samplebind` only. The scope is set via `files = ["src/pgen_samplebind"]` in `pyproject.toml`'s `[tool.mypy]` section, so `mypy` with no args (the form CI runs as `mypy src/`) checks just the library code. `pgenlib` has no stubs upstream, so it's marked `ignore_missing_imports`.

Running `mypy tests/` explicitly would type-check the test suite under strict mode and surface a different (noisier) set of findings — the strict gate is intentionally narrow to keep test-writing friction low. New tests should still type-annotate fixtures and helpers, but they aren't gated by CI.

## Commit + PR conventions

- **Conventional-commit style** for the subject line: `feat(scope): …`, `fix(scope): …`, `perf(scope): …`, `chore(deps): …`, `ci(scope): …`, `docs: …`, `test(scope): …`. Scopes mirror module names where applicable (`merge`, `afs`, `pvar`, `alignment`, `pseudohaploid`, …) or use `cli`/`reporting`/`packaging` for cross-cutting.
- **Body explains the why**, not just the what — what bug it fixes, why this approach, what alternatives were ruled out. Future maintainers (often us in six months) read commit bodies, not blame lines.
- **One logical change per commit** when feasible. The v0.3 perf pass landed five poles as five separate commits so each could be reverted independently if a regression slipped through.
- **PR title = top commit subject**; PR body summarizes the change set, links any related issue, and lists the verification path (tests added, dogfood tier verified, benchmark numbers).

## Release process

`pgen-samplebind` ships to PyPI via OIDC trusted publishing on tag push:

1. Land all changes on `main` via PR. Ensure CI is green.
2. Bump the version in `pyproject.toml` and `src/pgen_samplebind/__init__.py` (`__version__ = "X.Y.Z"`).
3. Update `CHANGELOG.md`: close `[Unreleased]` → open `[X.Y.Z] - YYYY-MM-DD` with the change set. Update the bottom-of-file diff links.
4. Update the `## Status` section in `README.md` with a one-paragraph summary.
5. Commit + push: `chore(release): vX.Y.Z — version bump + CHANGELOG + README finalize`.
6. Tag: `git tag -a vX.Y.Z -m "vX.Y.Z" && git push origin vX.Y.Z`.

The tag push fires `.github/workflows/release.yml`, which:

- Builds sdist + wheel.
- Validates metadata with `twine check --strict`.
- Smoke-tests the built wheel across the full 8-cell matrix (Python 3.11/3.12/3.13/3.14 × ubuntu/macos).
- Publishes to PyPI via OIDC trusted publishing — no API token in repo secrets.

Manual workflow dispatch (without a tag) runs the build + smoke-test path without publishing — useful for dry-running the release before tagging.

## Performance baseline maintenance

The CI bench gate (`tests/integration/test_perf_benchmark.py`) reads `tests/integration/perf_baseline.json` and fails if measured throughput drops below `threshold_pct_of_baseline × baseline` on the recorded fixture.

When an intentional perf change lands:

1. Push the change to `main` and observe the CI bench-job log for the `[perf] elapsed=… throughput=… peak_rss=…` line. (Surfaced on PASS thanks to the bench step's `-s` flag in `.github/workflows/ci.yml`.)
2. If the new throughput is significantly higher than the current baseline (≥ ~20% above), update `perf_baseline.json`:
   - Bump `throughput_genotypes_per_sec_baseline` to ~0.92x of the measured value (leaves ~8% headroom for runner variance — single-digit-percent in practice).
   - Add an entry to `_calibration_notes` recording the date, the measured number, and the rationale.
3. Commit as a separate `ci(perf): calibrate baseline to NM g/s` change. Do not silently roll the baseline into a feature commit — the calibration history is the only durable record of when each number was set.

Concrete example — the v0.3.1 calibration: the new fixture measured **70.84 M g/s** on `ubuntu-latest`; the baseline was bumped to **65 M g/s** (91.8% of measured, ≈0.92x), gating at the existing 0.80 threshold-pct = **52 M g/s** = 73% of current. That tolerates the observed ~8.5% run-to-run runner variance while catching a real ~27% regression. See the matching `_calibration_notes` entry in `perf_baseline.json` for the audit trail.

Do **not** auto-update the baseline. Silent drift defeats the gate's purpose.

## Design philosophy

Three principles, in tension order:

1. **Correctness first.** Byte-equality with the mergeit reference on the dogfood fixture is the regression bar. Every feature commit must keep the dogfood green; perf and ergonomics matter only after that.
2. **Bounded scope.** The tool does sample-bind. Not variant-union, not dosage, not multi-allelic, not BFILE-only output. Scope creep here adds maintenance debt that doesn't differentiate the project from plink2 once `--pmerge` lands the non-concatenating case upstream.
3. **Vectorized hot paths.** Python loops over n_variants or n_samples are a code smell. The 5-pole v0.3 perf pass (see `CHANGELOG.md` v0.3.0 entry) replaced every such loop in the merge hot path with numpy / pandas C-level passes; new code should follow the same pattern. Per-row `apply` / `iterrows` / Python `for` over series will be flagged in review.

When the right call is ambiguous, file an issue first — `pgen-samplebind` errs toward shipping a smaller surface than the user's first ask, because the project's value is in being one well-tested step rather than five plausible ones.
