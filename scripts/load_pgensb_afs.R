#!/usr/bin/env Rscript
# load_pgensb_afs.R — read pgen-samplebind `afs` subcommand output into the
# three-data-frame list shape that AdmixTools 2's `*_to_afs()` family returns.
#
# Usage from R:
#   source("load_pgensb_afs.R")
#   afs <- load_pgensb_afs("path/to/afs_output_dir")
#   # afs$freq, afs$counts, afs$snp now match AT2's eigenstrat_to_afs() shape
#
# Bridge until `pfile_to_afs()` lands in admixtools upstream.

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
