# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Shared logic for the source-tarball retired-dependency scan, factored out
# of tests/verify_source_tarball.R so it can be exercised by a fast testthat
# unit test without paying for a full `R CMD build`.
#
# NEWS.md and inst/THIRD_PARTY_NOTICES.md are the canonical places to
# disclose that the package named by retired_layout_package (GPL-2 only)
# still reaches an installation as a transitive dependency of `fgsea`;
# GPL-3 licence compliance requires naming it there (see the "Licence of the
# lisaR package" section). That is not a regression of the 2026-09-03
# refactor that removed it as a direct GeneCard layout dependency, so those
# two files are exempt from the scan below while every other shipped text
# file remains covered.
lisa_retired_layout_exempt_basenames <- function() {
  c("NEWS.md", "THIRD_PARTY_NOTICES.md")
}

lisa_retired_layout_hits <- function(files, retired_layout_package = paste0("cow", "plot"),
                                      exempt_basenames = lisa_retired_layout_exempt_basenames()) {
  candidates <- files[!basename(files) %in% exempt_basenames]
  candidates[vapply(candidates, function(path) {
    text <- paste(suppressWarnings(readLines(path, warn = FALSE)), collapse = "\n")
    text <- iconv(text, from = "", to = "UTF-8", sub = "byte")
    grepl(tolower(retired_layout_package), tolower(text), fixed = TRUE)
  }, logical(1))]
}
