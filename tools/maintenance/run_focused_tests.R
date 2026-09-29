# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Run focused lisaR testthat files from source with a pinned Rscript runtime.
# Usage: Rscript --vanilla tools/maintenance/run_focused_tests.R 'workers|run-management'

filter <- commandArgs(trailingOnly = TRUE)
filter <- if (length(filter)) filter[[1L]] else NULL
testthat::test_local(filter = filter, reporter = "summary")
