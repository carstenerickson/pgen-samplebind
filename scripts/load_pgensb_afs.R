#!/usr/bin/env Rscript
# load_pgensb_afs.R — read pgen-samplebind `afs` subcommand output into a
# three-data-frame list (freq, counts, snp).
#
# This is a RAW loader: it returns the AFS exactly as pgen-samplebind
# wrote it. It does NOT apply the filter AT2's `extract_f2()` silently
# applies (`discard_from_aftable(maxmiss=0, minmaf=0, maxmaf=0.5, ...)`)
# before writing its cache. Feeding this raw output into AT2's
# `afs_to_f2()` will produce f2 / qpAdm results that diverge from
# `extract_f2()` on the same panel (concrete failure mode documented in
# the Phase 7 dogfood-2 writeup).
#
# For feeding pgen-samplebind AFS into AT2's f2 / qpAdm chain, use the
# sibling script `pgensb_afs_to_at2_f2_cache.R` which applies the right
# filter and drives `afs_to_f2` end-to-end into an AT2-ready cache.
#
# Use this raw loader for inspecting / diffing / debugging AFS contents
# directly, not as the AT2 entry point.
#
# Usage from R:
#   source("load_pgensb_afs.R")
#   afs <- load_pgensb_afs("path/to/afs_output_dir")
#   # afs$freq, afs$counts, afs$snp are unfiltered data frames.

suppressPackageStartupMessages({
  library(jsonlite)
})

load_pgensb_afs <- function(dir) {
  manifest_path <- file.path(dir, "afs_manifest.json")
  if (!file.exists(manifest_path)) {
    stop(sprintf("manifest not found at %s — is this a pgen-samplebind afs output directory?",
                 manifest_path))
  }
  manifest <- jsonlite::fromJSON(manifest_path)

  snp_path    <- file.path(dir, manifest$files$snp)
  freq_path   <- file.path(dir, manifest$files$freq)
  counts_path <- file.path(dir, manifest$files$counts)

  snp    <- read.delim(snp_path,    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  freq   <- read.delim(freq_path,   sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)
  counts <- read.delim(counts_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

  # Move variant_id to row names for compatibility with AT2's matrix-style
  # AFS shape (rows = variants, cols = populations).
  rownames(freq)   <- freq$variant_id
  rownames(counts) <- counts$variant_id
  freq$variant_id   <- NULL
  counts$variant_id <- NULL

  list(
    freq      = freq,
    counts    = counts,
    snp       = snp,
    manifest  = manifest
  )
}

# When invoked as a script with a dir arg, print a summary.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  afs <- load_pgensb_afs(args[1])
  cat(sprintf("Loaded AFS from %s:\n", args[1]))
  cat(sprintf("  Variants:    %d\n", nrow(afs$snp)))
  cat(sprintf("  Populations: %d (%s)\n",
              ncol(afs$freq),
              paste(colnames(afs$freq), collapse = ", ")))
  cat(sprintf("  Tool version: %s\n", afs$manifest$tool_version))
  cat(sprintf("  Pseudohaploid adjustment: %s\n",
              ifelse(afs$manifest$adjust_pseudohaploid_applied, "applied", "not applied")))
}
