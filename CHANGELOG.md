# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-05-12

Quality-of-life and upstream-integration release. Four feature areas plus broadened dependency support: pseudohaploid sidecar reader for sibling-tool integration, repeatable `--target`, per-chromosome diagnostic histogram, and a live progress bar during the pass-2 genotype stream.

### Added

- **Pseudohaploid sidecar reader** (closes [#2](https://github.com/carstenerickson/pgen-samplebind/issues/2)). `merge` and `afs` now honor an optional `<prefix>.pseudohaploid.json` sidecar next to any input. Schema v1 maps `IID → {pseudohaploid: 0|1}`; missing sidecar is silent (current behavior preserved). Precedence: **sidecar > `.psam` PSEUDOHAPLOID column > heterozygosity inference**. The sidecar is written by sibling tools that know per-sample status by methodology (e.g., pileup-aadr's single-BAM `--randomDiploid` output is pseudohaploid by construction); honoring it lets that authoritative signal flow through to the output `.psam` PSEUDOHAPLOID column instead of being re-derived. Under `--on-collision suffix`, overrides correctly follow the input → output IID rename. Orphaned sidecar entries (IIDs not present in the input `.psam`) raise `InvariantViolation` (exit 3). 35 new tests (27 unit covering schema enforcement + suffix tolerance, 8 integration covering merge + afs + target-mode rename + regression).
- **`--target` is repeatable** for multi-sample append. Pass `--target` more than once to bind several targets in a single merge invocation; each is appended after the positional inputs and marked `is_target=True`. Gate (c) (`--target-min-call-rate`) fires per-target — strict semantics: any failing target blocks the merge. Collision-suffix scheme: bare `_target` is preserved when exactly one target is supplied (back-compat with v0.1); two or more targets switch to `_target_<input_idx>` so each renamed sample traces unambiguously to its source target. Motivating workflow: pileup-aadr produces one EIGENSTRAT triplet per ancient BAM, so a researcher running it on N samples collapses N loop iterations into a single `merge --target s1 --target s2 ... --target sN` invocation. 4 new integration tests in `TestMultiTarget`.
- **Per-chromosome action histogram.** `MergeCounters` gains `action_histogram_per_chrom: dict[int, dict[str, int]]` — the same 8 keys as the global `action_histogram`, bucketed by canonical chromosome. Surfaced in `--report-json` under `alignment.action_histogram_per_chrom` (JSON object keys stringified; consumers read with `int(k)`). Diagnostic for chr-specific drop concentrations that the autosome-wide average smooths out: chr 6 at 20% ambiguous-strand drops while autosomes average 5% signals an HLA-region strand artifact in one source; chr 22 at 80% allele_mismatch drops usually means hg19/hg38 build disagreement. Strictly additive — global `action_histogram` stays as the sum across chroms, no breaking change for existing JSON consumers. 4 new unit tests + JSON-shape integration assertion.
- **Pass-2 progress bar** (`tqdm`). The variant-block loop in `merge.merge_inputs` is now wrapped with a live progress bar when `sys.stderr.isatty()` AND `--quiet` is not passed. Interactive terminal runs see `Pass 2: streaming genotypes [████░░░░] 60% · 720k/1.2M variants · 23M geno/s · ETA 2m14s`; piped contexts (workflow managers, `--quiet`, redirected stderr) stay silent so CI logs and Snakemake/Nextflow runs aren't polluted. Bar is pure instrumentation: byte-identical output regardless of whether it's displayed. Adds `tqdm>=4.66,<5` to dependencies (~80 KB wheel). 3 new integration tests.

### Changed

- **`requires-python = ">=3.11,<3.15"`** — bumped the Python cap from `<3.13` to admit Python 3.13 and 3.14. The previous cap was tied to numpy 1.26's lack of 3.13 support; with the numpy upper bound now at `<3` (admits 2.x), 3.13 and 3.14 work cleanly. Verified with the full 278-test suite on both Python 3.13 + numpy 2.4.4 and Python 3.14 + numpy 2.4.4. CI matrix now runs `[3.11, 3.12, 3.13, 3.14] × [ubuntu, macos]` = 8 cells (up from 4); the release-pipeline smoke-test matrix matches.
- **Bumped `numpy` and `pandas` upper bounds** to allow current major versions:
  - `numpy>=1.26,<2` → `numpy>=1.26,<3` (admits numpy 2.x)
  - `pandas>=2.2,<3` → `pandas>=2.2,<4` (admits pandas 3.x)
  Pgenlib 0.94.0 was compiled against the stable NumPy C API, so it works against both 1.x and 2.x ABIs without rebuild. Validated locally with the full test suite (278 passed) on numpy 2.4.4 + pandas 3.0.3 + pgenlib 0.94.0. The AADR dogfood `dogfood_full` tier in CI runs the same numerical-parity check against the mergeit reference TSV on every commit, so future numpy / pandas / pgenlib bumps are auto-validated end-to-end through to qpAdm.


### Added

- **AADR-derivative dogfood regression test** (`tests/dogfood/`). A 44-sample × 50K-variant fixture drawn from AADR v66 under fair-use for non-commercial scholarly verification, with a vendored mergeit-pipeline reference qpAdm TSV. Six tests across three tiers (default / `dogfood_plink2` / `dogfood_full`) verify pgen-samplebind reproduces the established `mergeit + plink2 + AdmixTools 2` pipeline end-to-end. Includes provenance script `build_fixture.py` so the fixture can be regenerated from AADR + auxiliary panels. Default-tier tests run in ~10s with no external dependencies; full-tier (with R + admixtools installed) verifies qpAdm output numerically matches the reference within 1e-6 on weights and 1e-4 on p_tail (cross-architecture float-arithmetic noise floor is ~1e-12). See `tests/dogfood/README.md` for fair-use rationale and run instructions.


### Fixed

- **`pgen-samplebind afs` write-path 15-min → ~2-min on 1240k-scale panels.** Previous implementation used `pandas.to_csv(float_format="%.10g")` for the freq matrix, which is single-threaded printf-per-cell and dominated wallclock at scale (1.1M variants × 30 populations = 34M format calls). Replaced with a numpy-vectorized `np.char.mod("%.7g", block)` + chunked join in `afs._write_freq_tsv_fast`. ~50× speedup on large panels; below-noise-floor precision impact on downstream f-statistics (`%.7g` ≈ 23 bits, well below the float64 arithmetic-order noise in AT2's jackknife chain).

### Added

- **`scripts/pgensb_afs_to_at2_f2_cache.R`** — end-to-end R bridge that reads a pgen-samplebind AFS bundle and writes an AdmixTools 2-compatible f2 cache directory. Applies the `discard_from_aftable(maxmiss=0, ...)` filter that `extract_f2` silently applies (without which downstream `afs_to_f2` produces divergent f2; this was discovered during the Phase 7 dogfood-2). Drives `afs_to_f2()` for both `type='f2'` (with `poly_only=TRUE`) and `type='ap'` (with `poly_only=FALSE`) to match `extract_f2(qpfstats=FALSE)`'s behavior. Output is consumable directly by `f2_from_precomp(afprod=TRUE)` → `qpadm()` etc. The previous `scripts/load_pgensb_afs.R` is now documented as a raw loader for inspection — not the AT2 entry point.
- README AFS section updated to reference the new end-to-end bridge and document the qpfstats-bypass limitation.


### Added

- **`afs` subcommand** — bridge for the PFILE-native pipeline. Streams genotypes from a PFILE/BFILE/EIGENSTRAT input and emits three TSVs + a manifest matching AdmixTools 2's `*_to_afs()` shape (per-variant ALT frequencies, called-allele counts, SNP metadata). With pseudohaploid adjustment by default (`--no-pseudohaploid-adjust` to skip). Includes `scripts/load_pgensb_afs.R` for direct loading into the AT2 three-data-frame format. Removes the need for a `plink2 --make-bed` last-mile conversion before f-statistic computation. Bridge until `pfile_to_afs()` lands in admixtools upstream.

## [0.1.0] - 2026-05-10

Initial public release. The missing `plink2 --pmerge` non-concatenating case for ancient-DNA / population-genetics workflows.

### Added

- **`merge` subcommand** — bind two or more PFILE/BFILE/EIGENSTRAT inputs sharing a variant set into one output PFILE.
  - Two-pass architecture: pass 1 alignment via pandas, pass 2 genotype streaming via pgenlib `read_range` / `append_biallelic_batch` with chromosome-aware block iteration.
  - Allele resolution: passthrough, REF/ALT swap, strand flip, strand-flip-and-swap, drop, fill-missing.
  - Policy flags: `--on-mismatch`, `--on-missing`, `--on-extra`, `--on-strand`, `--on-collision`, `--trust-strand`.
  - `--target` mode for single-sample / small-cohort append: asymmetric strand-check, per-sample call-rate gate (`--target-min-call-rate`).
  - `--on-collision suffix`: `_<input_idx>` general suffix scheme (`_target` in target mode), with idempotent numeric retry layered as `<base>_<suffix>_<n>`.
  - `--relabel-from`: 2-col header-less POP→POP collapse, or N-col header-required per-sample override (e.g., AADR `Master ID`→`Group ID`).
  - Pseudohaploid detection per-sample via heterozygosity tally; written to output `.psam` `PSEUDOHAPLOID` column.
  - Population label auto-detection (POP / PHENO / PHENO1 fallback); FID auto-mirrors POST-relabel POP.
  - cM (genetic-position) column preserved end-to-end so downstream Morgan-spaced jackknife consumers (e.g., AT2 `extract_f2 blgsize=0.05`) get correct block partitioning.
- **`validate` subcommand** — pass-1 only, no genotype reads, no output written. Same Exit-1 gates as merge; `--on-* error` policies softened to gate (d) per HLD §Exit-1 validation gates.
- **`hash` subcommand** — canonical variant-set hash for cross-format panel-identity verification. PFILE/BFILE/EIGENSTRAT yield identical hashes for the same variant set.
- **`inspect` subcommand** — structured summary of one input: format, samples, variants, populations, sex distribution, missingness histogram. TSV by default, `--json` for machine-readable.
- **EIGENSTRAT input in both flavors**: PACKEDANCESTRYMAP (binary; converted via `plink2 --eigfile --make-pgen`) and ASCII per-line (parsed natively, no convertf pre-conversion). Format auto-detected from the `.geno` header bytes.
- **BFILE input** via `plink2 --bfile --make-pgen` shell-out into a per-invocation `tempfile.TemporaryDirectory`.
- **Reports**: `--report` per-variant TSV (streamed; constant memory) and `--report-json` summary JSON (default ~few KB; `--report-json-include-rows` opt-in with >100 MB stderr warning).
- **Concurrency**: `fcntl.flock`-based advisory output-prefix lock (`{prefix}.lock`); raises exit 2 on contention; NFS/SMB/CIFS detection emits a stderr warning.
- **Stable exit codes**: 0 success / 1 validation failure / 2 I/O failure / 3 invariant violation / 4 usage error.
- **CI matrix**: ubuntu-latest + macos-latest × Python 3.11/3.12, plus a dedicated linux-x86_64 perf-benchmark gate (HLD test 18) that fails on regression below 80% of the recorded baseline.
- **PyPI release pipeline** (`.github/workflows/release.yml`): tag-triggered (`v*`); builds sdist + wheel; validates with `twine check --strict`; smoke-tests the built wheel against the unit suite on the full CI matrix; publishes to PyPI via OIDC trusted publishing (no long-lived API token in repo secrets).
- **`py.typed` marker** ships with the wheel so type-checkers honor the package's inline type hints.

### Defaults worth noting

- `--trust-strand` is **off** by default. A/T and C/G ambiguous SNPs are dropped wherever strand cannot be verified — including the case where canonical and other inputs both have the same ambiguous pair (e.g., both have A/T at the same position). Strand cannot be proven the same because complementing A/T gives T/A which is the same pair. This matches `mergeit`'s `strandcheck: YES` convention and is the safe default for cross-source merges. Pass `--trust-strand` for single-source pipelines (same AADR release, same processing) where REF/ALT calls are guaranteed consistent.
- `--on-extra warn` (default) fires gate (a) → exit 1 when extras exceed 10% of canonical's variant count. Catches the input-order-reversed failure mode (smaller panel placed first; larger panel's distinct variants silently dropped).
- `--validate-strand-fail-pct 10.0` (default) fires gate (b) → exit 1 when ambiguous-strand drops exceed 10% of the alignment intersection. Intersection denominator catches the wrong-panel failure mode at small intersection sizes.

### Test coverage

262 tests at release: 165 unit + ~95 integration + dogfood end-to-end byte-equality proof against the `mergeit + plink2 + awk` reference pipeline on a real-world Reich-Lab-style ancient-DNA panel build (md5-identical proximal qpAdm shootout output). HLD test 17 (`mergeit f2` parity) and HLD test 18 (perf benchmark) are gated to dedicated CI cells; HLD test 16 (panel hash invariance) is the manual nightly path.

### Known limitations

- Biallelic SNPs only. Multi-allelic input is rejected at startup with a clear preprocessing recipe.
- Phase / dosage data not supported (out of scope; AdmixTools is hardcall-based).
- Variant union with sample-bind not supported (use `plink1.9 --bmerge` for that).
- BFILE-only output not supported (use `plink2 --pfile out --make-bed` if needed).
- EIGENSTRAT/BFILE input requires `plink2 v2.0.0-a.7.1+` on PATH (the `--eigfile`/`--make-pgen` path); pure-PFILE workflows have no plink2 dependency.

[Unreleased]: https://github.com/carstenerickson/pgen-samplebind/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/carstenerickson/pgen-samplebind/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/carstenerickson/pgen-samplebind/releases/tag/v0.1.0
