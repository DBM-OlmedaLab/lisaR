# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

write_lisa_manifest <- function(rows, output_file) {
  if (!is.data.frame(rows)) {
    stop("rows must be a data.frame", call. = FALSE)
  }
  if (!all(c("stage", "path", "status") %in% names(rows))) {
    stop("manifest rows must contain stage, path and status columns", call. = FALSE)
  }
  write_lisa_tsv(rows, output_file)
}
