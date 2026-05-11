# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **`--trust-strand` semantic extended** to cover the matching-allele ambiguous case. By default A/T and C/G SNPs are now dropped whenever strand cannot be verified — including when canonical and other inputs both have the same ambiguous pair (previously: passthrough). Matches mergeit's `strandcheck: YES` convention and is the safe default for cross-source merges. Pass `--trust-strand` for single-source pipelines where REF/ALT calls are guaranteed consistent.
- **cM column preserved end-to-end** through the merge pipeline. Output `.pvar` now carries the genetic-position values from input `.pvar`/`.snp` files (default 0.0 if absent). Previously cM was dropped, which propagated to `.bim` after `plink2 --make-bed` and broke Morgan-spaced jackknife block partitioning in downstream tools (e.g., AT2 `extract_f2 blgsize=0.05`).

### Added

- **Native ASCII per-line EIGENSTRAT input.** plink2 `--eigfile` only reads PACKEDANCESTRYMAP-format EIGENSTRAT (binary `GENO ` header). Older ASCII per-line files (one digit per sample-variant cell, no header) are now parsed natively — no convertf pre-conversion required. Format is auto-detected from the `.geno` header bytes.
- **PyPI release pipeline.** `.github/workflows/release.yml` triggers on `v*` tag pushes, builds sdist + wheel, validates with `twine check --strict`, smoke-tests the built wheel against the unit suite on `ubuntu-latest` + `macos-latest` × Python 3.11/3.12, and publishes to PyPI via OIDC trusted publishing (no long-lived API token in repo secrets). Manual `workflow_dispatch` runs build + test only as a dry-run.
- **`py.typed` marker** ships with the wheel so type-checkers honor the package's inline type hints. `Typing :: Typed` classifier added.

### Packaging

- pyproject.toml polished: added `Repository` and `Changelog` URLs, expanded classifiers (Linux/macOS, Python 3-only, `Console`, `Typing :: Typed`), added `twine` to dev extras.

- `merge` subcommand: bind two or more PFILE/BFILE/EIGENSTRAT inputs sharing a variant set into one output PFILE.
  - Two-pass architecture: pass 1 alignment via pandas, pass 2 genotype streaming via pgenlib `read_range` / `append_biallelic_batch` with chromosome-aware block iteration.
  - Allele resolution: passthrough, REF/ALT swap, strand flip, strand-flip-and-swap, drop, fill-missing.
  - Policy flags: `--on-mismatch`, `--on-missing`, `--on-extra`, `--on-strand`, `--on-collision`, `--trust-strand`.
  - `--target` mode for single-sample / small-cohort append: asymmetric strand-check, per-sample call-rate gate (`--target-min-call-rate`).
  - `--on-collision suffix`: `_<input_idx>` general suffix scheme (`_target` in target mode), with idempotent numeric retry layered as `<base>_<suffix>_<n>`.
  - `--relabel-from`: 2-col header-less POP→POP collapse, or N-col header-required per-sample override (e.g., AADR `Master ID`→`Group ID`).
  - Pseudohaploid detection per-sample via heterozygosity tally; written to output `.psam` `PSEUDOHAPLOID` column.
  - Population label auto-detection (POP / PHENO / PHENO1 fallback); FID mirrors POST-relabel POP.
- `validate` subcommand: pass-1 only, no genotype reads, no output written. Same Exit-1 gates as merge; `--on-* error` policies softened to gate (d) per HLD §Exit-1 validation gates.
- `hash` subcommand: canonical variant-set hash for cross-format panel-identity verification (PFILE/BFILE/EIGENSTRAT yield identical hashes for the same variant set).
- `inspect` subcommand: structured summary of one input — format, samples, variants, populations, sex distribution, missingness histogram (TSV / `--json`).
- EIGENSTRAT and BFILE input via `plink2 --eigfile`/`--bfile --make-pgen` shell-out into a per-invocation `tempfile.TemporaryDirectory`. Requires plink2 v2.0.0-a.7.1+ on PATH.
- Reports: `--report` per-variant TSV (streamed; constant memory) and `--report-json` summary JSON (default ~few KB; `--report-json-include-rows` opt-in with >100 MB stderr warning).
- Concurrency: `fcntl.flock`-based advisory output-prefix lock (`{prefix}.lock`); raises exit 2 on contention; NFS/SMB/CIFS detection emits stderr warning.
- Stable exit codes: 0 success / 1 validation failure / 2 I/O failure / 3 invariant violation / 4 usage error.
- CI matrix: ubuntu-latest + macos-latest × Python 3.11 + 3.12, plus a dedicated linux-x86_64 perf-benchmark gate (HLD test 18) that fails on regression below 80% of the recorded baseline (currently 25 M genotypes/sec end-to-end).
- Test coverage: 248 tests including 23 HLD correctness/regression tests (1-15, 19-23) and full subprocess exit-code harness; HLD #17 (mergeit f2 parity) and #16 (Phase 6/7 panel hash invariance) gated to nightly external_tool / dogfood runs.

### Documentation

- README expanded with the four canonical use cases (panel extension, target append, cross-source merge, AADR cross-version cohort assembly), full CLI reference, exit-code table, validation-gate explanation, plink2 a7.x troubleshooting.

[Unreleased]: https://github.com/carstenerickson/pgen-samplebind/compare/v0.1.0...HEAD
