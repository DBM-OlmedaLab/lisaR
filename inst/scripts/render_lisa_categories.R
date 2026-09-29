#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


args <- commandArgs(trailingOnly = TRUE)
values <- list(source_run = "", output_dir = "", selection = "full")
while (length(args)) {
  if (length(args) < 2L || !startsWith(args[[1]], "--")) stop("Arguments must be --flag value pairs.", call. = FALSE)
  key <- gsub("-", "_", sub("^--", "", args[[1]]))
  if (!key %in% names(values)) stop("Unknown argument: ", args[[1]], call. = FALSE)
  values[[key]] <- args[[2]]; args <- args[-c(1L, 2L)]
}
if (!nzchar(values$source_run)) stop("Required argument: --source-run", call. = FALSE)
if (!nzchar(values$output_dir)) stop("Required argument: --output-dir", call. = FALSE)
lisaR::render_lisa_categories(values$source_run, values$selection, values$output_dir)
