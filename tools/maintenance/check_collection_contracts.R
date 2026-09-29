# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Execute the focused collection and resource validation tests.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
root <- normalizePath(root, mustWork = TRUE)

old <- Sys.getenv("LISAR_RUN_PENDING_CONTRACTS", unset = NA_character_)
old_max_fails <- Sys.getenv("TESTTHAT_MAX_FAILS", unset = NA_character_)
old_options <- options(
  testthat.progress.max_fails = 100L,
  testthat.summary.max_reports = 100L
)
on.exit({
  options(old_options)
  if (is.na(old)) {
    Sys.unsetenv("LISAR_RUN_PENDING_CONTRACTS")
  } else {
    Sys.setenv(LISAR_RUN_PENDING_CONTRACTS = old)
  }
  if (is.na(old_max_fails)) {
    Sys.unsetenv("TESTTHAT_MAX_FAILS")
  } else {
    Sys.setenv(TESTTHAT_MAX_FAILS = old_max_fails)
  }
}, add = TRUE)
Sys.setenv(LISAR_RUN_PENDING_CONTRACTS = "true")
Sys.setenv(TESTTHAT_MAX_FAILS = "100")

results <- testthat::test_local(
  root,
  filter = "collection-contracts",
  reporter = "summary"
)

if (!testthat:::all_passed(results)) {
  quit(save = "no", status = 1L, runLast = FALSE)
}
