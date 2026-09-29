#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Reproduce one NES annotation variant from its exact shared source and settings.
# Usage: Rscript reproduce_lisa_category_nes.R source.tsv settings.json output.png
# Settings declare the presentation and device context. New figures use one
# fresh pinned R process per output; archived settings retain legacy rendering.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: reproduce_lisa_category_nes.R source.tsv settings.json output.png", call. = FALSE)
if (!requireNamespace("lisaR", quietly = TRUE)) stop("The matching lisaR package must be installed.", call. = FALSE)
getFromNamespace("lisa_reproduce_category_nes_variant", "lisaR")(args[[1L]], args[[2L]], args[[3L]])
