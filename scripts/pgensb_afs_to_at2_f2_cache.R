#!/usr/bin/env Rscript
# pgensb_afs_to_at2_f2_cache.R — end-to-end conversion from a pgen-samplebind
# `afs` output bundle to an AdmixTools 2 f2-cache directory, ready for
# `f2_from_precomp()` / `qpadm()` / etc.
#
# This is the script to use when you want PFILE-native AT2 work without the
# `plink2 --make-bed` last-mile step. The lighter-weight `load_pgensb_afs.R`
# helper in this directory only loads the AFS into memory; it does NOT
# apply the filter `extract_f2` silently applies (`discard_from_aftable
# (maxmiss=0, minmaf=0, maxmaf=0.5, ...)`) before writing its cache. Feeding
# an unfiltered AFS into `afs_to_f2()` produces divergent downstream f2 and
# qpAdm — see the Phase 7 dogfood-2 writeup.
#
# Usage:
#   Rscript pgensb_afs_to_at2_f2_cache.R <afs_bundle_dir> <out_f2_cache_dir>
#
# Then from R:
#   library(admixtools)
#   f2 <- f2_from_precomp("<out_f2_cache_dir>", pops = my_pops, afprod = TRUE)
#   qpadm(f2, left = ..., right = ..., target = ...)
#
# Limitation: AT2's `extract_f2(qpfstats=TRUE)` reads genotypes directly and
# bypasses the AFS layer entirely. This bridge covers the non-qpfstats path
# only. If you need qpfstats (e.g., ancient-DNA with high missingness), the
# PFILE → BED conversion (`plink2 --pfile <pref> --make-bed`) remains
# necessary for that single AT2 call.

suppressPackageStartupMessages({
  library(admixtools)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: pgensb_afs_to_at2_f2_cache.R <afs_bundle_dir> <out_f2_cache_dir>")
}
AFS_BUNDLE   <- normalizePath(args[1])
F2_CACHE_DIR <- args[2]

# Intermediate afdir (the layout AT2's extract_afs produces and afs_to_f2 reads).
# We keep it as a sister dir to the f2 cache so the user can inspect / re-run.
AT2_AFDIR    <- paste0(F2_CACHE_DIR, ".afdir")
dir.create(F2_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AT2_AFDIR,    recursive = TRUE, showWarnings = FALSE)
F2_CACHE_DIR <- normalizePath(F2_CACHE_DIR)
AT2_AFDIR    <- normalizePath(AT2_AFDIR)

cat(sprintf("[afs->f2] %s\n  -> afdir:   %s\n  -> f2cache: %s\n",
            AFS_BUNDLE, AT2_AFDIR, F2_CACHE_DIR))

# ---- 1. Load AFS bundle ----------------------------------------------------

snp_df    <- read_tsv(file.path(AFS_BUNDLE, "afs_snp.tsv"),    col_types = cols())
freq_df   <- read_tsv(file.path(AFS_BUNDLE, "afs_freq.tsv"),   col_types = cols())
counts_df <- read_tsv(file.path(AFS_BUNDLE, "afs_counts.tsv"), col_types = cols())

snpdat <- snp_df %>%
  transmute(
    SNP = as.character(variant_id),
    CHR = as.character(chrom),
    cm  = as.numeric(cm),
    POS = as.numeric(pos),
    A1  = as.character(ref),
    A2  = as.character(alt)
  )
pop_cols   <- setdiff(colnames(freq_df), "variant_id")
afs_mat    <- as.matrix(freq_df[, pop_cols])
counts_mat <- as.matrix(counts_df[, pop_cols])
rownames(afs_mat)    <- snpdat$SNP
rownames(counts_mat) <- snpdat$SNP
colnames(afs_mat)    <- pop_cols
colnames(counts_mat) <- pop_cols
cat(sprintf("  loaded %d variants, %d populations\n", nrow(snpdat), length(pop_cols)))

# ---- 2. Apply AT2's hidden `extract_f2` filter -----------------------------
# `extract_f2` calls `discard_from_aftable(maxmiss=0, minmaf=0, maxmaf=0.5,
# minac2=FALSE, auto_only=TRUE)` before writing its cache. The strict
# `maxmiss=0` filter drops any variant where ANY population has zero
# called alleles. We inline it here rather than calling
# `admixtools:::discard_from_aftable` since that function isn't exported.

keep_no_miss <- rowSums(counts_mat == 0) == 0
afs_mat    <- afs_mat[keep_no_miss, , drop = FALSE]
counts_mat <- counts_mat[keep_no_miss, , drop = FALSE]
snpdat     <- snpdat[keep_no_miss, , drop = FALSE]
cat(sprintf("  after maxmiss=0 filter: %d variants (dropped %d)\n",
            nrow(snpdat), sum(!keep_no_miss)))

# Polymorphic flag — `extract_f2` records this in snpdat and afs_to_f2 uses
# it when `poly_only=TRUE`.
snpdat$poly <- apply(afs_mat, 1, function(x) any(x > 0 & x < 1, na.rm = TRUE))
cat(sprintf("  polymorphic: %d / %d\n", sum(snpdat$poly), length(snpdat$poly)))

# ---- 3. Write the afdir layout afs_to_f2 expects ---------------------------

write_tsv(snpdat, file.path(AT2_AFDIR, "snpdat.tsv.gz"))
saveRDS(afs_mat,    file = file.path(AT2_AFDIR, "afs1.rds"))
saveRDS(counts_mat, file = file.path(AT2_AFDIR, "counts1.rds"))

# ---- 4. Build the f2 cache -------------------------------------------------
# Mirror `extract_f2`'s default `poly_only = c('f2')`:
#   - type='f2'  uses polymorphic-only variants
#   - type='ap'  uses all variants
# Both are needed by `f2_from_precomp(afprod=TRUE)`.

cat("  computing type='f2' (polymorphic only)\n")
afs_to_f2(AT2_AFDIR, F2_CACHE_DIR, chunk1 = 1, chunk2 = 1, blgsize = 0.05,
          type = "f2", poly_only = TRUE,  apply_corr = TRUE,
          overwrite = FALSE, verbose = TRUE)

cat("  computing type='ap' (all variants; required for afprod=TRUE)\n")
afs_to_f2(AT2_AFDIR, F2_CACHE_DIR, chunk1 = 1, chunk2 = 1, blgsize = 0.05,
          type = "ap", poly_only = FALSE, apply_corr = TRUE,
          overwrite = FALSE, verbose = TRUE)

# Stamp the cache so downstream AT2 helpers (which expect extract_f2's
# `.cached_pops.txt` manifest) can detect it as complete.
writeLines(sort(pop_cols), file.path(F2_CACHE_DIR, ".cached_pops.txt"))
file.create(file.path(F2_CACHE_DIR, ".extract_f2_done"))

cat(sprintf("[afs->f2] done. f2 cache ready at %s\n", F2_CACHE_DIR))
cat("Next from R:\n")
cat("  library(admixtools)\n")
cat(sprintf("  f2 <- f2_from_precomp(\"%s\", pops = ..., afprod = TRUE)\n",
            F2_CACHE_DIR))
cat("  qpadm(f2, left = ..., right = ..., target = ...)\n")
