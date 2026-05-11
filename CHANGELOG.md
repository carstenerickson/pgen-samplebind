# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/carstenerickson/pgen-samplebind/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/carstenerickson/pgen-samplebind/releases/tag/v0.1.0
