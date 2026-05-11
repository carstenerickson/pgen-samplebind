# Dogfood test — AADR v66 panel build replication

This directory contains a regression test that exercises `pgen-samplebind`
end-to-end on a real Allen Ancient DNA Resource (AADR) v66 derivative,
reproducing the established `mergeit + plink2 + AdmixTools 2` pipeline. The
output qpAdm verdicts must be byte-equal to the vendored reference.

## What this test exercises

Replicates a Patterson-2022-style ancestry decomposition workflow:

1. **Multi-input merge with cross-source ambiguous-strand drops**:
   28 Patterson-source samples (v66) + 12 English target samples (c29i) +
   1 single-sample target → 41-sample, 12-population panel.

2. **`--target` mode with pseudohaploid sample**: one
   Patterson_England_IA individual (`I17014`) is appended via `--target`
   with the pop-label `dogfood_target` so it gets its own qpAdm row.

3. **cM preservation through the merge pipeline**: panel `.bim` must
   carry the genetic-position values from the source `.snp` files (a
   real bug found and fixed during the Phase 7 dogfood; the fixture is
   one of the regression guards against that bug recurring).

4. **AT2 qpAdm parity (when available)**: `pgen-samplebind`'s output
   panel produces md5-equal qpAdm results to the
   `mergeit + plink2 + AT2` reference pipeline on the same fixture
   data. The vendored `qpadm_reference.tsv` is the reference; the
   test runs the pipeline and asserts byte-equality.

## Source attribution

- **Dataset**: AADR v66.0 (April 2026 release).
- **Maintainer**: David Reich Lab, Harvard Medical School.
- **Original data**: <https://doi.org/10.7910/DVN/FFIDCW>.
- **License**: CC-BY-NC-ND 4.0.

If you're doing actual research with AADR data, fetch the full release
from the David Reich Lab Dataverse. This fixture is intentionally too
small for meaningful population-genetics inference (50,000 variants vs
the full ~1.15M autosomal 1240k panel; 44 samples vs ~16,000 in v66).

## Fair-use rationale

This derivative is used under fair-use doctrine for non-commercial
scholarly software-verification purposes:

1. **Non-commercial scholarly use** — open-source software regression
   test. The repository, the test, and the derived data are all freely
   distributed without charge.
2. **Factual source material** — genotype calls are uncopyrightable
   biological facts; the ND clause applies to the compilation, not the
   underlying facts. A small derivative for testing purposes does not
   substitute for the compilation.
3. **Small fraction of original** — 44 samples (~0.3 % of v66's ~16K) ×
   50,000 variants (~4 % of 1.15M autosomal 1240k) = ~0.012 % of the
   original dataset by cell count.
4. **No commercial harm** — this fixture is free and does not compete
   with the original distribution. Anyone wanting AADR data for real
   research goes to the primary source.
5. **Transformative use** — the data has been subsetted, repackaged,
   relabeled, and bound to test code. The fixture's purpose is
   software verification, not the original purpose of population-
   genetics inference.

## Fixture contents

| File | Size | Description |
|---|---|---|
| `panel_v66_subset.{geno,snp,ind}` | ~4 MB | 28 samples × 50K variants — Patterson 7 source pops (WHGA, WHGB, Balkan_N, OldSteppe, OldAfrica, Turkey_N, Russia_Afanasievo), 4 per pop, drawn from AADR v66 via the Master-ID join used in our internal Phase 7 panel build |
| `brit_subset_subset.{geno,snp,ind}` | ~4 MB | 12 samples × 50K variants — 4 English target pops (England_C_EBA, _MBA, _LBA, _IA), 3 per pop, drawn from a v44-era cohort that's also AADR-derived |
| `target_individual.{geno,snp,ind}` | ~4 MB | 1 sample (`I17014`, originally Patterson_England_IA in v44) × 50K variants — relabeled `dogfood_target` so the test can name it as a `--target`-mode target |
| `qpadm_reference.tsv` | <1 KB | mergeit-pipeline reference: 5 target rows × 8 cols (3 source weights + 3 SEs + sum_w + p_tail). md5 = `0e654786c60ab3cefdb61664727bab73`. Pre-computed once on this fixture data using `mergeit + plink2 --make-bed + AT2 extract_f2(qpfstats=FALSE) + qpadm` |
| `build_fixture.py` | 9 KB | The cloud-side fixture builder. Not run by the test; committed for provenance and so this fixture is reproducible from AADR + brit_subset if someone needs to regenerate or extend it |

Total fixture size: ~12 MB raw, ~2.4 MB if gzipped.

## Running the dogfood test

### Default CI (no plink2 or AT2 needed)

```bash
pytest tests/dogfood/ -v
```

Always-runnable subset: validates `pgen-samplebind merge` produces
the expected panel shape (41 samples × 49,382 variants after
1,236 strand-ambiguous drops), correct FID=POP, populated
PSEUDOHAPLOID column, preserved cM column in the output `.pvar`.

### With plink2 on PATH

```bash
pytest tests/dogfood/ -v -m "dogfood or dogfood_plink2"
```

Adds: PFILE → BED conversion + `.bim` schema check (cM column non-zero
distinct values, FIDs match populations).

### Full pipeline including AT2 R package (external_tool)

```bash
pytest tests/dogfood/ -v -m "dogfood_full"
```

Adds: full qpAdm shootout + md5-compare against
`qpadm_reference.tsv`. Requires `R` with `admixtools` installed
(typically via the `uqrmaie1/admixtools` GitHub repo).

## Why this fixture and not synthesized data

Synthetic data would let us avoid the AADR license question entirely.
We chose the AADR derivative because:

1. **Realistic ancient-DNA characteristics**: pseudohaploid samples
   (adaptive pulldown methodology), cross-source strand ambiguity
   (~6% A/T+C/G dropped by default), cM-preserving merge requirements.
   Synthetic data would have to mimic these characteristics; using real
   data is the cleanest demonstration.

2. **Compatibility with published research**: the qpAdm shootout uses
   the Patterson 3-source model (WHG/EEF/Steppe), the canonical
   ancient-DNA decomposition. A maintainer or contributor can sanity-
   check the reference numbers (e.g., Patterson_England_IA at WHG≈0.10,
   EEF≈0.33, Steppe≈0.56) against published values.

3. **Trust signal for the community**: the test reproduces the
   established pipeline on real ancient-DNA data. That's a stronger
   claim than "we passed unit tests on synthetic fixtures."

## Regenerating the fixture

If AADR changes or new pops should be added: `build_fixture.py` walks
the v62 → v66 Master-ID join, picks samples deterministically (seed
`0xD06F00D`), and subsets variants. Run on a host that has AADR v66 +
the v44-era brit_subset locally; not part of any user-facing flow.

The vendored `qpadm_reference.tsv` was generated once with:
- `mergeit` (EIGENSOFT, from AT2's Docker image)
- `plink2 v2.0.0-a.7.1`
- `admixtools` R package commit hash from May 2026

Any future-AT2-version regeneration that produces different reference
numbers should be a deliberate, documented update to the fixture, not
a silent change.
