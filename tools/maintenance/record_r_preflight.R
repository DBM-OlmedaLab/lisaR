# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Record the pinned R runtime used for reproducible lisaR maintenance checks.
# Run with the absolute Rscript executable from the selected Conda environment.

namespaces <- c(
  "lisaR", "testthat", "DESeq2", "limma", "fgsea", "AnnotationDbi",
  "KEGGREST", "org.Hs.eg.db", "org.Mm.eg.db", "rmarkdown", "yaml",
  "jsonlite", "knitr"
)
commands <- c("salmon", "samtools", "bedtools", "fastqc", "multiqc", "pandoc")

rscript_suffix <- if (.Platform$OS.type == "windows") ".exe" else ""
rscript_default <- file.path(R.home("bin"), paste0("Rscript", rscript_suffix))
cat("Rscript=", Sys.getenv("LISAR_RSCRIPT", unset = rscript_default), "\n", sep = "")
cat("R.version=", R.version.string, "\n", sep = "")
cat("R.home=", R.home(), "\n", sep = "")
cat("CONDA_PREFIX=", Sys.getenv("CONDA_PREFIX", unset = "unset"), "\n", sep = "")
cat("libPaths=", paste(.libPaths(), collapse = " | "), "\n", sep = "")
for (package in namespaces) {
  available <- requireNamespace(package, quietly = TRUE)
  version <- if (available) as.character(utils::packageVersion(package)) else "unavailable"
  cat("namespace.", package, "=", version, "\n", sep = "")
}
for (command in commands) {
  path <- Sys.which(command)
  cat("command.", command, "=", if (nzchar(path)) path else "unavailable", "\n", sep = "")
}
