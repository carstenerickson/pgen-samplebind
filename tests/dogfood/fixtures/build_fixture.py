#!/usr/bin/env python3
"""Build the public dogfood fixture from AADR v66 + c29i brit_subset.

Outputs four EIGENSTRAT triplets (ASCII, since plink2 --eigfile needs
PACKEDANCESTRYMAP; we'll convert via plink2 in the next step):
  panel_v66_subset.{geno,snp,ind}    -- 28 Patterson source samples, 50K variants
  brit_subset_subset.{geno,snp,ind}  -- 16 English target samples, same 50K
  target_individual.{geno,snp,ind}   -- 1 Patterson_England_IA pulled aside
  keep_variants.tsv                  -- the chosen 50K variant IDs

Deterministic seed: 0xD06F00D.

Subsetting strategy:
- 4 samples per Patterson pop (7 source pops x 4 = 28), seeded random selection
- 3 samples per English target pop (3 of 4 pops x 3 = 9 brit_subset samples),
  plus 1 Patterson_England_IA pulled aside as the target sample (4 of 4 pops covered)
- 50K random autosomal variants (seeded), evenly distributed across chr 1-22

Configuration via environment variables (no hardcoded local paths):

    PGENSB_DOGFOOD_V66_PREFIX    AADR v66.0 PACKEDANCESTRYMAP prefix
                                 (.geno + .snp + .ind, sample-major)
    PGENSB_DOGFOOD_V62_ANNO      AADR v62 .anno TSV (for the v62→Master ID
                                 join; used to resolve a small number of
                                 sample-ID changes between v62 and v66)
    PGENSB_DOGFOOD_V62_PATCHED   v62-patched .ind file with corrected
                                 Patterson_England_IA pop labels
    PGENSB_DOGFOOD_BRIT_PREFIX   c29i English-targets EIGENSTRAT prefix
                                 (3 of 4 Patterson_England_* target pops)
    PGENSB_DOGFOOD_OUTDIR        Where to write the four fixture triplets +
                                 keep_variants.tsv. Defaults to ./build/.

The fixture-build is a provenance / regeneration step, not part of CI —
the vendored `panel_v66_subset.*`, `brit_subset_subset.*`,
`target_individual.*` triplets in this directory are the artifacts CI
consumes. This script exists so anyone with the source AADR + c29i data
can independently verify that those vendored triplets follow the
documented subsetting strategy. Set the five environment variables to
point at your local AADR v66 + v62 anno + c29i build, then run.
"""

import os
import random
import sys
from collections import defaultdict
from pathlib import Path

SEED = 0xD06F00D
TARGET_VARIANTS = 50_000
SAMPLES_PER_SOURCE_POP = 4
SAMPLES_PER_TARGET_POP = 3

_REQUIRED_ENV = {
    "PGENSB_DOGFOOD_V66_PREFIX": "AADR v66.0 PACKEDANCESTRYMAP prefix",
    "PGENSB_DOGFOOD_V62_ANNO": "AADR v62 .anno TSV for v62→Master ID join",
    "PGENSB_DOGFOOD_V62_PATCHED": "v62-patched .ind with corrected pop labels",
    "PGENSB_DOGFOOD_BRIT_PREFIX": "c29i English-targets EIGENSTRAT prefix",
}
_missing = [(k, desc) for k, desc in _REQUIRED_ENV.items() if not os.environ.get(k)]
if _missing:
    msg = "missing required environment variable(s):\n"
    for k, desc in _missing:
        msg += f"  {k}  ({desc})\n"
    msg += "See module docstring for the full configuration interface."
    sys.exit(msg)

V66_PREFIX = os.environ["PGENSB_DOGFOOD_V66_PREFIX"]
V62_ANNO = os.environ["PGENSB_DOGFOOD_V62_ANNO"]
V62_PATCHED = os.environ["PGENSB_DOGFOOD_V62_PATCHED"]
BRIT_PREFIX = os.environ["PGENSB_DOGFOOD_BRIT_PREFIX"]
OUTDIR = Path(os.environ.get("PGENSB_DOGFOOD_OUTDIR", "./build"))

PATTERSON_SOURCES = [
    "Patterson_WHGA",
    "Patterson_WHGB",
    "Patterson_Balkan_N",
    "Patterson_OldSteppe",
    "Patterson_OldAfrica",
    "Patterson_Turkey_N",
    "Patterson_Russia_Afanasievo",
]
PATTERSON_TARGETS = [
    "Patterson_England_C_EBA",
    "Patterson_England_MBA",
    "Patterson_England_LBA",
    "Patterson_England_IA",
]

# ---------------------------- Master-ID join ---------------------------


def build_v62_to_master():
    """v62 GeneticID -> v62 MasterID, from v62 anno col 1 -> col 2."""
    m = {}
    with open(V62_ANNO) as f:
        next(f)  # header
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                m[cols[0]] = cols[1]
    return m


def build_master_to_v66():
    """v66 IndividualID -> v66 GeneticID. v66 anno col 3 -> col 1.
    Priority: AG > DG > SG (matches rung-8's preference)."""
    PRIORITY = {".AG": 3, ".DG": 2, ".SG": 1}
    m = {}
    seen_prio = {}
    with open(f"{V66_PREFIX}.anno") as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            gid, iid = cols[0], cols[2]
            prio = 0
            for suf, p in PRIORITY.items():
                if suf in gid:
                    prio = p
                    break
            if iid not in seen_prio or prio > seen_prio[iid]:
                m[iid] = gid
                seen_prio[iid] = prio
    return m


def build_v66_iid_set():
    """Set of v66 GeneticIDs that actually appear in v66.ind."""
    s = set()
    with open(f"{V66_PREFIX}.ind") as f:
        for line in f:
            parts = line.split()
            if parts:
                s.add(parts[0])
    return s


# ---------------------------- Sample selection -------------------------


def find_patterson_samples_in_v66():
    """Walk v62_patched.ind; for each Patterson_* row, find the v66 GeneticID
    via Master-ID join. Returns dict[pop] -> list of v66 GeneticIDs."""
    v62_to_master = build_v62_to_master()
    master_to_v66 = build_master_to_v66()
    v66_ind_set = build_v66_iid_set()

    by_pop = defaultdict(list)
    with open(V62_PATCHED) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            v62_gid, _sex, pop = parts[0], parts[1], parts[2]
            if pop not in PATTERSON_SOURCES:
                continue
            master = v62_to_master.get(v62_gid)
            if not master:
                continue
            v66_gid = master_to_v66.get(master)
            if v66_gid and v66_gid in v66_ind_set:
                by_pop[pop].append(v66_gid)
    return by_pop


def select_samples():
    """Deterministic sample selection. Returns three lists:
    (panel_samples, brit_samples, target_sample) where each panel/brit row is
    (sample_id, sex, pop) — sex is filled with 'U' since we don't read it
    here (plink2 picks it up from the source .ind / .anno)."""
    rng = random.Random(SEED)

    # Sources: pick from v66 via Master-ID join
    by_pop = find_patterson_samples_in_v66()
    panel_samples = []
    for pop in PATTERSON_SOURCES:
        pool = sorted(by_pop.get(pop, []))
        rng.shuffle(pool)
        chosen = pool[:SAMPLES_PER_SOURCE_POP]
        for iid in chosen:
            panel_samples.append((iid, "U", pop))
        print(f"  {pop}: {len(pool)} available, picked {len(chosen)}")

    # English targets: pull from brit_subset; pick 3 per pop for brit_subset
    # plus 1 Patterson_England_IA as the standalone --target
    brit_samples_by_pop = defaultdict(list)
    with open(f"{BRIT_PREFIX}.ind") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            iid, sex, pop = parts[0], parts[1], parts[2]
            if pop in PATTERSON_TARGETS:
                brit_samples_by_pop[pop].append((iid, sex, pop))

    brit_samples = []
    target_sample = None
    for pop in PATTERSON_TARGETS:
        pool = sorted(brit_samples_by_pop[pop], key=lambda r: r[0])
        rng.shuffle(pool)
        if pop == "Patterson_England_IA":
            # Pull one out as the target, take the next 3 for brit_subset.
            target_sample = pool[0]
            chosen = pool[1 : 1 + SAMPLES_PER_TARGET_POP]
        else:
            chosen = pool[:SAMPLES_PER_TARGET_POP]
        brit_samples.extend(chosen)
        print(
            f"  {pop}: {len(pool)} available, picked {len(chosen)}"
            + (f" + 1 target ({target_sample[0]})" if pop == "Patterson_England_IA" else "")
        )

    return panel_samples, brit_samples, target_sample


# ---------------------------- Variant selection ------------------------


def select_variants():
    """50K random autosomal variants, deterministic. Returns sorted list of
    SNP IDs."""
    rng = random.Random(SEED + 1)
    auto_variants = []
    with open(f"{V66_PREFIX}.snp") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                chrom = int(parts[1])
            except ValueError:
                continue
            if 1 <= chrom <= 22:
                auto_variants.append(parts[0])
    print(f"  v66 autosomal variants: {len(auto_variants):,}")

    # Also intersect with brit_subset variants
    brit_variants = set()
    with open(f"{BRIT_PREFIX}.snp") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                chrom = int(parts[1])
            except ValueError:
                continue
            if 1 <= chrom <= 22:
                brit_variants.add(parts[0])
    print(f"  brit_subset autosomal variants: {len(brit_variants):,}")

    intersection = [v for v in auto_variants if v in brit_variants]
    print(f"  intersection: {len(intersection):,}")

    rng.shuffle(intersection)
    chosen = sorted(intersection[:TARGET_VARIANTS])
    print(f"  chosen: {len(chosen):,}")
    return chosen


# ---------------------------- Write outputs ----------------------------


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("=== Sample selection ===")
    panel_samples, brit_samples, target_sample = select_samples()
    print(f"\n  panel: {len(panel_samples)} samples")
    print(f"  brit:  {len(brit_samples)} samples")
    print(f"  target: 1 sample ({target_sample[0]})")

    print("\n=== Variant selection ===")
    chosen_variants = select_variants()

    # Write keeplist files
    with open(OUTDIR / "keep_variants.tsv", "w") as f:
        for v in chosen_variants:
            f.write(v + "\n")

    with open(OUTDIR / "panel_keep_samples.tsv", "w") as f:
        for iid, sex, pop in panel_samples:
            f.write(f"{iid}\t{sex}\t{pop}\n")

    with open(OUTDIR / "brit_keep_samples.tsv", "w") as f:
        for iid, sex, pop in brit_samples:
            f.write(f"{iid}\t{sex}\t{pop}\n")

    with open(OUTDIR / "target_keep_sample.tsv", "w") as f:
        iid, sex, pop = target_sample
        f.write(f"{iid}\t{sex}\t{pop}\n")

    print(f"\n=== Wrote keeplist files to {OUTDIR} ===")
    for fn in (
        "keep_variants.tsv",
        "panel_keep_samples.tsv",
        "brit_keep_samples.tsv",
        "target_keep_sample.tsv",
    ):
        with open(OUTDIR / fn) as fh:
            n = sum(1 for _ in fh)
        print(f"  {fn}: {n} rows")


if __name__ == "__main__":
    main()
