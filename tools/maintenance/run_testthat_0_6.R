# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
if (length(args) >= 2L) .libPaths(c(normalizePath(args[[2]], mustWork = TRUE), .libPaths()))
library(testthat)
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, quiet = TRUE, export_all = TRUE)
} else {
  library(lisaR)
}
target <- if (length(args) >= 3L) file.path(root, "tests", "testthat", args[[3]]) else
  file.path(root, "tests", "testthat")
result <- if (file.info(target)$isdir) {
  testthat::test_dir(target, reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE)
} else {
  testthat::test_file(target, reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE)
}
failed <- vapply(result, function(item) {
  any(vapply(item$results, inherits, logical(1), what = c("expectation_failure", "expectation_error")))
}, logical(1))
if (any(failed)) quit(status = 1L)
