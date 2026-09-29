# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
files <- c(
  list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
  list.files(file.path(root, "inst", "scripts"), pattern = "[.]R$", full.names = TRUE)
)
failures <- vapply(files, function(path) {
  tryCatch({ parse(path); "" }, error = conditionMessage)
}, character(1))
failures <- failures[nzchar(failures)]
if (length(failures)) {
  for (path in names(failures)) message(path, ": ", failures[[path]])
  quit(status = 1L)
}
message("Parsed ", length(files), " R files successfully.")
