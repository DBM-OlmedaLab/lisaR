#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Maintainer check. This file is excluded from the source tarball.
# It verifies that installed runtime helpers have one source copy only.

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) args[[1L]] else normalizePath(".", mustWork = TRUE)
canonical_dir <- file.path(repo, "inst", "scripts")
legacy_dir <- file.path(repo, "scripts")

if (!dir.exists(canonical_dir)) {
  stop("Canonical runtime-script directory is missing: ", canonical_dir,
       call. = FALSE)
}

canonical <- list.files(canonical_dir, pattern = "\\.[Rr]$", full.names = FALSE)
if (!length(canonical)) {
  stop("Canonical runtime-script directory contains no R scripts.", call. = FALSE)
}

if (dir.exists(legacy_dir)) {
  duplicates <- intersect(
    canonical,
    list.files(legacy_dir, pattern = "\\.[Rr]$", full.names = FALSE)
  )
  if (length(duplicates)) {
    stop(
      "Runtime scripts have duplicate source copies under top-level scripts/: ",
      paste(sort(duplicates), collapse = ", "),
      call. = FALSE
    )
  }
}

absolute_hits <- character()
for (path in file.path(canonical_dir, canonical)) {
  text <- readLines(path, warn = FALSE)
  if (any(grepl("(?:/home|/Users)/[[:alnum:]_.-]+/", text, fixed = FALSE))) {
    absolute_hits <- c(absolute_hits, path)
  }
}
if (length(absolute_hits)) {
  stop(
    "Canonical runtime scripts contain development-only absolute paths: ",
    paste(absolute_hits, collapse = ", "),
    call. = FALSE
  )
}

vignette_dir <- file.path(repo, "vignettes")
vignettes <- list.files(vignette_dir, pattern = "[.]Rmd$")
if (!length(vignettes)) stop("No R Markdown vignette sources found.")
generated_dirs <- file.path(repo, c("docs/vignettes", "inst/vignettes-md"))
generated_files <- unlist(lapply(generated_dirs, list.files, recursive = TRUE))
if (length(generated_files) ||
    length(list.files(vignette_dir, pattern = "[.]md$"))) {
  stop("Keep vignette sources in vignettes/; do not commit generated copies.")
}

cat("SOURCE_LAYOUT_GATE=PASS\n")
cat("canonical_runtime_scripts=", length(canonical), "\n", sep = "")
cat("vignette_sources=", length(vignettes), "\n", sep = "")
