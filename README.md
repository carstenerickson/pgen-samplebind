# pgen-samplebind

A focused CLI for the one thing `plink2 --pmerge` doesn't yet do: bind two or more PFILE/BFILE/EIGENSTRAT genotype datasets that share a variant set but contain **different samples**. Targeted at the ancient-DNA / population-genetics community working with the 1240k SNP capture panel and downstream AdmixTools 2 / qpAdm pipelines.

```bash
pgen-samplebind merge panel_a panel_b -o merged
```

That's it for the happy path. Variant alignment, strand canonicalization, ambiguous-strand policy, missing-call fill, IID collision handling, pseudohaploid detection, population-label preservation, and report emission all run with one command.

## Why this exists

`plink2 --pmerge` errors with "under development" for the non-concatenating case (same variants, different samples). Today's workarounds — `plink1.9 --bmerge`, EIGENSOFT's `mergeit`, or VCF round-trips through `bcftools merge` — each carry friction (deprecated tool, format-conversion overhead, single-thread C engine, multi-format hops). The AdmixTools / ancient-DNA community spends nontrivial time on this every panel build.

`pgen-samplebind` is a stopgap until plink2 lands the feature natively; it's deliberately scoped to the bind operation only — not a plink2 replacement.

## Install

```bash
pip install pgen-samplebind
```

Requires Python 3.11 through 3.14. CI tests all four versions on both Ubuntu and macOS. For EIGENSTRAT or BFILE inputs, also install **plink2 v2.0.0-a.7.1 or newer** and put it on PATH. Pick the right asset for your platform:

```bash
# Linux x86_64
ASSET=plink2_linux_x86_64.zip
# macOS Apple Silicon
# ASSET=plink2_mac_arm64.zip
# macOS Intel
# ASSET=plink2_mac.zip

curl -fsSL "https://github.com/chrchang/plink-ng/releases/download/v2.0.0-a.7.1/${ASSET}" -o /tmp/plink2.zip
unzip /tmp/plink2.zip -d ~/bin
plink2 --version  # confirm v2.0.0-a.7.1 (or newer) on PATH
```

The `v2.0.0-a.7.x` line introduced the `--eigfile`/`--export eig` support pgen-samplebind shells out to; older versions including `v2.00a5.x` (still shipped by Homebrew at the time of this writing) will silently lack EIGENSTRAT support. PFILE-only workflows have no plink2 dependency.

## Canonical use cases

Four workflows the tool is shaped around. Each is one command.

### 1. Panel extension

Existing AADR-derived panel (~700 samples × ~1.15M 1240k variants), add new populations from a newer AADR release. Variants identical, alleles canonical (same source).

```bash
pgen-samplebind merge \
    /data/aadr_v66_phase6 \
    /data/aadr_v66_new_pops \
    -o /data/aadr_v66_phase7 \
    --report phase7_bind_report.tsv
```

### 2. Single-sample target append

A user's WGS pseudohaploid genotypes (from `pileupCaller --randomDiploid`) merged into an AADR-derived panel for downstream qpAdm. The user is missing some panel sites (low-coverage regions) — the tool emits missing-call codes for them, doesn't drop the variants. `--target` activates the asymmetric strand-check + per-sample call-rate gate (default 0.40).

```bash
pgen-samplebind merge \
    --target /data/carsten_pileupcaller \
    /data/aadr_v66_subset \
    -o /data/aadr_with_carsten
```

### 3. Cross-source merge (different formats, different strand conventions)

Two cohort releases from different labs — one in EIGENSTRAT, one in PFILE — both built on 1240k. Strand canonicalization may differ; ambiguous A/T and C/G SNPs may need explicit policy.

```bash
# cohort_a is EIGENSTRAT (cohort_a.geno + .snp + .ind); cohort_b is PFILE
# (cohort_b.pgen + .pvar + .psam). Format is auto-detected from the prefix.
pgen-samplebind merge \
    /data/cohort_a \
    /data/cohort_b \
    -o /data/cohort_ab \
    --on-strand flip \
    --validate-strand-fail-pct 10
```

### 4. AADR cross-version cohort assembly

Recreating a published cohort using AADR v66 genotypes (more SNPs) requires mapping v44.3 sample IDs to v66 sample IDs through the Master-ID join — AADR de-anonymized many samples between releases (`I0001` became `Loschbour.AG`, etc.). Sample-bind plus identity remapping, one step.

```bash
pgen-samplebind merge \
    /data/aadr_v44_3_panel \
    /data/aadr_v66_extras \
    --id-column 'Master ID' \
    --relabel-from /data/aadr_v66.anno \
    --relabel-input-col 'Master ID' \
    --relabel-output-col 'Group ID' \
    -o /data/cross_version_panel
```

## Subcommands

```
pgen-samplebind merge    INPUT [INPUT ...] -o OUTPUT  [options]
pgen-samplebind validate INPUT [INPUT ...]            [options]
pgen-samplebind afs      INPUT             -o OUTDIR  [options]
pgen-samplebind hash     INPUT                        [options]
pgen-samplebind inspect  INPUT                        [options]
```

Inputs are PFILE/BFILE/EIGENSTRAT prefixes; format auto-detected from companion files (`.pgen`+`.pvar`+`.psam` / `.bed`+`.bim`+`.fam` / `.geno`+`.snp`+`.ind`).

### `merge` — bind inputs into one output PFILE

| Option | Default | Purpose |
|---|---|---|
| `-o, --out PREFIX` | required | Output PFILE prefix |
| `--target PATH` | none | Single-sample / small-cohort mode; activates asymmetric strand-check + per-target call-rate gate. Targets are appended after the positional inputs; canonical remains input[0]. Repeatable — pass `--target` multiple times to append several targets in one merge |
| `--variant-key {chr_pos\|id}` | `chr_pos` | Match key |
| `--on-mismatch {drop\|error}` | `drop` | Allele mismatch beyond resolution |
| `--on-missing {fill_missing\|drop_variant\|error}` | `fill_missing` | Variant in input[0] absent in input N |
| `--on-extra {warn\|drop\|error}` | `warn` | Variant in input N absent from input[0] |
| `--on-strand {drop\|flip\|error}` | `flip` | Strand mismatch on unambiguous SNPs |
| `--trust-strand` | off | Pass A/T and C/G ambiguous SNPs without flagging (footgun for cross-source panels) |
| `--on-collision {error\|first\|suffix}` | `error` | IID collision policy. `suffix` uses `_<input_idx>` (general); `_target` (single-target mode); `_target_<input_idx>` (multi-target mode, ≥ 2 targets, so renames disambiguate). Idempotent numeric retry on further collisions |
| `--id-column NAME` | `IID` | `.psam` column for identity ops (e.g., `'Master ID'` for AADR anno) |
| `--population-column NAME` | auto | Holds population labels (default: POP / PHENO / PHENO1 fallback) |
| `--target-min-call-rate FLOAT` | `0.40` | Target-mode per-sample minimum call rate before exit-1 |
| `--validate-strand-fail-pct N` | `10.0` | Exit 1 if ambiguous-strand drops exceed N% of intersection |
| `--relabel-from FILE` | none | TSV-driven POP relabel. 2-col header-less form maps POP→POP; N-col form (with `--relabel-input-col` / `--relabel-output-col`) joins per-sample on `--id-column` |
| `--relabel-input-col NAME` | none | For N-col TSVs: which column matches `--id-column` |
| `--relabel-output-col NAME` | none | For N-col TSVs: which column becomes the new POP value |
| `--report PATH` | none | Per-variant action TSV (streamed; constant memory) |
| `--report-json PATH` | none | Run-level summary JSON (~few KB; rows excluded by default) |
| `--report-json-include-rows` | off | Include per-variant rows in JSON (buffered; warns at >100 MB predicted size) |
| `--quiet` | off | Suppress progress to stdout |
| `--block-size N` | `2048` | Variants per pgenlib read block |

### `validate` — check alignment without writing

Same alignment / strand options as `merge`. No output written; reports go to stdout plus optional `--report` / `--report-json`. Exits 0 if alignment OK, 1 if any of the gates below fires.

### `hash` — emit canonical variant-set hash

```bash
# Hashing the same variant set via two different formats yields the same
# digest: PFILE on /data/aadr_v66_subset.{pgen,pvar,psam} and EIGENSTRAT
# on /data/aadr_v66_subset_eig.{geno,snp,ind}. Format is auto-detected
# from the companion files present at the prefix.
pgen-samplebind hash /data/aadr_v66_subset
# 7c4f8e...  PFILE
pgen-samplebind hash /data/aadr_v66_subset_eig
# 7c4f8e...  EIGENSTRAT  ← same hash → format-equivalent panels
```

`--emit-canonical` prints the canonicalized bytestream that's hashed (for diagnosis when two inputs *should* match but don't).

### `afs` — per-population allele-frequency-spectrum TSVs

```bash
pgen-samplebind afs panel -o panel_afs/
```

Streams genotypes from a PFILE/BFILE/EIGENSTRAT input and emits three TSVs + a manifest JSON matching the shape AdmixTools 2's `*_to_afs()` family returns:

```
panel_afs/
├── afs_snp.tsv       (variant_id, chrom, pos, ref, alt, cm)
├── afs_freq.tsv      (variant_id × population, ALT-allele frequency)
├── afs_counts.tsv    (variant_id × population, called-allele counts)
└── afs_manifest.json (tool version, sample counts per pop, parameters)
```

Useful for the **PFILE-native pipeline**: skip the `plink2 --pfile … --make-bed` last-mile step before AT2's non-qpfstats f2 / qpAdm path. Bridge until `pfile_to_afs()` lands in admixtools upstream.

**Feeding AFS into AT2** — use the end-to-end bridge script:

```bash
Rscript scripts/pgensb_afs_to_at2_f2_cache.R panel_afs/ panel_at2_f2_cache/
```

This loads the AFS bundle, applies the `discard_from_aftable(maxmiss=0, …)` filter that AT2's `extract_f2` silently applies before writing its cache, then calls `afs_to_f2()` twice (once with `poly_only=TRUE` for `type='f2'`, once with `poly_only=FALSE` for `type='ap'`) to produce an AT2-ready f2 cache. From R:

```r
library(admixtools)
f2 <- f2_from_precomp("panel_at2_f2_cache/", pops = my_pops, afprod = TRUE)
qpadm(f2, left = ..., right = ..., target = ...)
```

The sibling `scripts/load_pgensb_afs.R` is a **raw loader** that returns the AFS as in-memory data frames without filtering — use it for inspection / debugging, not as the AT2 entry point (raw AFS fed into `afs_to_f2()` produces divergent f2 because `extract_f2`'s `maxmiss=0` filter isn't applied).

**Limitation**: AT2's `extract_f2(qpfstats=TRUE)` path reads genotypes directly and bypasses the AFS layer entirely. For workflows that need qpfstats (e.g., ancient-DNA with high missingness), the PFILE→BED last-mile remains. The byte-equal-to-mergeit-`qpfstats` proof comes from the PFILE→BED→qpfstats path, not this AFS bridge.

Key flags:
- `--populations POP` (repeatable) — restrict to a subset of populations
- `--no-pseudohaploid-adjust` — skip pseudohaploid called-allele adjustment (treats all samples as diploid). Default applies the adjustment when the input has a `PSEUDOHAPLOID` column — pseudohaploid samples then contribute 1 called allele (not 2) for correct effective sample sizes in downstream variance estimates. Allele frequencies are unchanged.
- `--include-sex-chrom` — extend to chr 23/24/25/26 (default: autosomes only)
- `--population-column NAME` — aggregate by a `.psam` column other than `POP`

### `inspect` — structured summary of one input

Format, sample count, variant count, populations, pseudohaploid mix, sex distribution, missingness histogram. TSV by default; `--json` for machine-readable.

## Validation gates (exit 1)

`merge` applies soft-validation gates (a)-(c); `validate` applies all four (a)-(d). Any firing gate exits 1.

- **(a) Extras above warn threshold.** Variants in input[N] absent from input[0] exceed the `--on-extra` warn threshold (10% of input[0] by default). Catches the input-order-reversed failure mode (smaller panel placed first; larger panel's distinct variants silently dropped).
- **(b) Ambiguous-strand drops above intersection threshold.** Drops due to A/T or C/G allele ambiguity exceed `--validate-strand-fail-pct` of the alignment intersection (10% by default). Intersection denominator is deliberate: catches the wrong-panel failure mode (tiny intersection × 30% drop rate) that a canonical denominator would silently hide.
- **(c) Target call rate below threshold.** In `--target` mode, target's non-missing call count over canonical variants is below `--target-min-call-rate`.
- **(d) `--on-* error` policy would have triggered.** Validate-mode only — `validate` softens the `error` policies into a count + gate-(d) trigger so it can report the full picture rather than aborting at the first hit. In `merge` mode those policies aren't softened: they raise an `InvariantViolation` and exit 3 (explicit invariant enforcement is honored).

For `merge`, gates (a)-(c) run between pass 1 (alignment) and pass 2 (genotype streaming) — failing fast saves the pass-2 wallclock for an alignment that wouldn't have validated.

## Exit codes

Stable across versions; safe to script against.

| Code | Meaning |
|---|---|
| `0` | Success |
| `1` | Validation failure — gate (a)/(b)/(c)/(d) fired |
| `2` | I/O failure — read/write failure, plink2 subprocess failed, advisory output-prefix lock held |
| `3` | Invariant violation — multi-allelic input, duplicate canonical variant keys, `--on-* error` policy triggered, `--on-collision error` on duplicate IIDs |
| `4` | Usage error — bad CLI argument combination, unknown input prefix |

## Concurrency

- **Same input prefixes**: safe, read-only.
- **Same output prefix**: tool advisory-locks `{out}.lock` via `fcntl.flock`; fails clearly with exit 2 if held; cleans up on success or signal.
- **Cache directory**: not managed by this tool — caller's problem.

The lock prevents two concurrent invocations from corrupting each other's output. It does not synchronize across machines (file lock is local-fs only). On Linux, the tool detects NFS/SMB/CIFS at lock-acquire time (via `/proc/self/mountinfo`) and emits a stderr warning since `fcntl.flock` semantics over network filesystems are implementation-defined and effectively no-op; the lock is still attempted. On macOS the detection is suppressed (no `/proc`-equivalent fs-type API) to avoid false positives — the network-fs caveat applies the same way, you just won't get the diagnostic warning.

## Performance

Expected throughput **50-70 M genotypes/sec end-to-end** on a typical Linux x86_64 machine (one core, default `--block-size 2048`), post-v0.3 vectorization. The v0.3.1 CI baseline measured 70.84 M g/s on a 100M-genotype synthetic fixture (`ubuntu-latest`); the perf benchmark gates against regression below 80% of that recorded baseline. For a 700-sample × 1.15M-variant 1240k panel, plan on roughly 15-25 seconds wallclock including pass-1 alignment, pass-2 genotype rewrite, and `.psam` finalization.

## Troubleshooting

### "plink2 not found, required for EIGENSTRAT input"

Install plink2 v2.0.0-a.7.1+ (see [Install](#install)). PFILE-only workflows have no plink2 dependency.

### "FID column header on line 1 is not at the beginning"

Symptom of an older plink2 reading our `.psam` output. We emit FID-first column order per the plink2 spec; if you see this from a downstream tool, verify the consumer is `plink2 v2.0.0-a.7.x` or newer. Earlier plink2 lines require `--no-fid` or an `.fam`-style fallback.

### `--target` mode call-rate gate fires unexpectedly

Default threshold is 0.40 (40% of canonical variants called for the target sample). For very-low-coverage targets this may be too strict. Pass `--target-min-call-rate 0.20` (or lower) if appropriate; with `0.0` the gate is fully disabled.

### "ASCII per-line EIGENSTRAT" inputs

Both flavors of EIGENSTRAT input are supported: PACKEDANCESTRYMAP (binary, `GENO `/`TGENO ` header — converted via plink2 `--eigfile`) and the older ASCII per-line variant (one digit per sample-variant cell — parsed natively, no plink2 dependency for the ASCII path). Format is auto-detected from the `.geno` header bytes.

### Cross-source merges drop more variants than expected

By default A/T and C/G ambiguous SNPs are dropped wherever strand cannot be verified — including the case where both inputs have the same allele pair (e.g., both have A/T at the same position). The strand cannot be proven the same because complementing A/T gives T/A which is the same pair. This matches mergeit's `strandcheck: YES` convention and is the safe default for cross-source merges (different cohorts, different processing pipelines).

For **single-source merges** where REF/ALT calls are guaranteed consistent across inputs (same AADR release, same conversion pipeline), pass `--trust-strand` to passthrough matching ambiguous variants. The 1240k panel has ~1.3% A/T+C/G matching ambiguous SNPs so this typically recovers 10-20K additional variants on a full 1.15M-variant panel.

### "ambiguous-strand drops exceed 10% of intersection"

Gate (b) fired. Three common causes:

1. **Wrong panel pairing**: tiny intersection × normal drop rate computes as a high *fraction* of intersection. Run `pgen-samplebind hash` on each input — if the hashes differ unexpectedly, the panels aren't what you thought they were.
2. **Different strand conventions**: cross-source merges where one source already strand-flipped at A/T+C/G sites.
3. **Legitimate but high A/T+C/G fraction**: 1240k naturally has ~5-8% ambiguous; a panel restricted to ambiguous-only sites would push above 10%. Raise the threshold with `--validate-strand-fail-pct 20` if that's your case, or use `--trust-strand` for explicit pass-through (footgun warning: silently passes potentially-flipped genotypes).

### Half-built output files after a failure

Don't happen. The merge orchestrator wraps pass 2 + `.psam` finalization in a try/except that unlinks the `.pgen`/`.pvar`/`.psam` triplet on any tool-internal exception, so downstream pipelines never silently consume a partial output. (Atomic-rename across the triplet isn't actually atomic on NFS — three separate file ops — so unlink-on-failure is the safer default.)

### Output prefix is locked by another pgen-samplebind process

Another invocation holds the lock at `{prefix}.lock`. If that process has died (kill -9, OOM, lost terminal), the lock file persists and you must remove it manually:

```bash
ls -la /data/merged.lock     # check lock-file age vs. expected runtime
rm /data/merged.lock         # only after confirming no live process holds it
```

## Verification

End-to-end byte-equal qpAdm parity against the established `mergeit + plink2 + AdmixTools 2` pipeline is verified on every commit via a vendored AADR-derivative regression test (Patterson 7-source + 4 English target pops + 1 individual target, drawn from AADR v66 under fair-use). The result is the project's trust artifact — anyone can clone, run, and confirm pgen-samplebind reproduces the established pipeline on a published-research-shape workload without trusting the maintainer's claims. See [CONTRIBUTING.md §Dogfood](CONTRIBUTING.md#dogfood--published-research-workflow-regression-test) for run instructions and the three-tier breakdown.

## Status

v0.3.0 — performance + cleanup release. Five vectorizations across the merge hot path (pass-1 alignment truth-table lookup, pass-2 block placement, chromosome-cast, duplicate-key check, action histogram) cut Python-level loops from the long poles; estimated 20-40% wallclock reduction at 1240k scale. Drops the never-implemented `--threads` flag. End-to-end byte-equal qpAdm parity continues to be proven on every commit against the `mergeit + plink2 + AdmixTools 2` reference pipeline via the AADR-derivative dogfood CI.

See [CHANGELOG.md](CHANGELOG.md) for the full feature list and known limitations.

## Contributing

Issues and pull requests welcome at <https://github.com/carstenerickson/pgen-samplebind/issues>. Dev setup, test-runner conventions, commit + release process, and design philosophy are in [CONTRIBUTING.md](CONTRIBUTING.md). The project is small and opinionated — substantive scope changes (e.g., multi-allelic merges, dosage data, BFILE-only output) should start with a design discussion before implementation.

## License

MIT — see [LICENSE](LICENSE).
