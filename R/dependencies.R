# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Dependency diagnostics for optional lisaR capabilities.

#' Report availability of optional lisaR dependencies
#'
#' @param packages Character vector of package names.
#'
#' @return A data frame with package availability and installed versions.
#' @keywords internal
lisa_optional_dependency_status <- function(packages) {
  packages <- as.character(packages)
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  versions <- vapply(packages, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
  }, character(1))
  data.frame(package = packages, available = available, version = versions, stringsAsFactors = FALSE)
}

#' Require an optional dependency with an actionable repair message
#'
#' @param package Package name.
#' @param feature Human-readable feature name.
#'
#' @return Invisibly `TRUE` when the dependency is available.
#' @keywords internal
lisa_require_optional <- function(package, feature) {
  if (requireNamespace(package, quietly = TRUE)) return(invisible(TRUE))
  manager <- if (package %in% lisa_bioconductor_packages()) {
    sprintf("BiocManager::install(%s)", shQuote(package))
  } else {
    sprintf("install.packages(%s)", shQuote(package))
  }
  stop(
    sprintf(
      "Optional package '%s' is required for %s. Install it with `%s`, or rerun the standalone lisaR installer with `--profile=full`.",
      package,
      feature,
      manager
    ),
    call. = FALSE
  )
}

lisa_bioconductor_packages <- function() {
  c(
    "AnnotationDbi", "DESeq2", "fgsea", "KEGGREST",
    "org.Hs.eg.db", "org.Mm.eg.db", "SummarizedExperiment"
  )
}

lisa_config_dependency_requirements <- function(
  pipeline,
  species,
  report = list(),
  available = function(package) requireNamespace(package, quietly = TRUE)
) {
  packages <- c("fgsea", "ggplot2", "hommel", "jsonlite", "openxlsx", "yaml")
  reasons <- rep("core lisaR workflow", length(packages))

  add <- function(package, reason) {
    packages <<- c(packages, package)
    reasons <<- c(reasons, reason)
  }

  run_gene_level <- identical(lisa_config_get(report, "mode", "standard"), "full")
  run_reports <- TRUE
  run_kegg <- lisa_config_bool(lisa_config_get(pipeline, "run_kegg_maps", FALSE), "pipeline.run_kegg_maps")

  if (run_gene_level) {
    for (package in c("ggrepel", "gridExtra", "patchwork")) {
      add(package, "gene-level figures and cards")
    }
    # The gene-level cards and overlays are saved through ggplot2::ggsave(),
    # whose SVG device is svglite. The category figures of a standard run use
    # the cairo device in grDevices instead and do not need it. Declaring the
    # requirement here turns a missing svglite into a readiness failure before
    # anything is written, rather than an error raised after the whole
    # canonical run has already completed.
    svg_requested <- isTRUE(lisa_config_bool(
      lisa_config_get(lisa_config_get(report, "formats", list()), "svg", FALSE),
      "report.formats.svg"
    ))
    if (svg_requested) {
      add("svglite", "SVG gene-level cards and overlays in the full report")
    }
  }
  if (run_reports) {
    add("openxlsx", "report and XLSX output")
  }
  if (run_kegg) {
    for (package in c("AnnotationDbi", "KEGGREST", "png")) {
      add(package, "KEGG pathway figures")
    }
    annotation_package <- switch(
      tolower(as.character(species)),
      "homo sapiens" = "org.Hs.eg.db",
      "mus musculus" = "org.Mm.eg.db",
      ""
    )
    if (nzchar(annotation_package)) {
      add(annotation_package, paste("gene annotation for", species))
    }
  }

  out <- data.frame(
    package = packages,
    reason = reasons,
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$package), , drop = FALSE]
  out$manager <- ifelse(
    out$package %in% lisa_bioconductor_packages(),
    "Bioconductor",
    "CRAN"
  )
  out$available <- vapply(out$package, available, logical(1))
  row.names(out) <- NULL
  out
}

lisa_assert_config_dependencies <- function(
  pipeline,
  species,
  report = list(),
  available = function(package) requireNamespace(package, quietly = TRUE)
) {
  status <- lisa_config_dependency_requirements(
    pipeline,
    species,
    report = report,
    available = available
  )
  missing <- status[!status$available, , drop = FALSE]
  if (!nrow(missing)) return(status)

  cran <- missing$package[missing$manager == "CRAN"]
  bioc <- missing$package[missing$manager == "Bioconductor"]
  repair <- character()
  if (length(cran)) {
    repair <- c(
      repair,
      sprintf(
        "install.packages(c(%s))",
        paste(shQuote(cran), collapse = ", ")
      )
    )
  }
  if (length(bioc)) {
    repair <- c(
      repair,
      "install.packages('BiocManager')",
      sprintf(
        "BiocManager::install(c(%s))",
        paste(shQuote(bioc), collapse = ", ")
      )
    )
  }
  stop(
    "LISA-DEPENDENCY-001 missing packages for the requested configuration: ",
    paste(missing$package, collapse = ", "),
    ". Repair before running: ",
    paste(repair, collapse = "; "),
    ". The standalone kit can install the complete set with `--profile=full`.",
    call. = FALSE
  )
}
