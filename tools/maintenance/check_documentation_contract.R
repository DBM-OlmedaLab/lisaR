# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

fail <- function(...) stop(paste0(...), call. = FALSE)
read_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
contains_ci <- function(needle, haystack) {
  grepl(tolower(needle), tolower(haystack), fixed = TRUE)
}

schema_path <- file.path(root, "inst", "schema", "lisa-config.schema.json")
schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
if (!identical(schema$properties$pipeline$properties$schema_version$const, "1.0.0")) {
  fail("Documentation gate requires configuration schema 1.0.0.")
}
if (!identical(schema$properties$source_data$items$additionalProperties, FALSE)) {
  fail("source_data rows are not strict in the JSON Schema.")
}
if (!identical(
  schema$properties$pipeline$properties$evidence$additionalProperties,
  FALSE
)) {
  fail("pipeline.evidence is not strict in the JSON Schema.")
}

files <- c(
  file.path(root, "README.md"),
  list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$", full.names = TRUE),
  file.path(
    root, "inst", "extdata", "dictionaries",
    "TECHNICAL_REBUILD_PROVENANCE.md"
  )
)
missing <- files[!file.exists(files)]
if (length(missing)) fail("Missing documentation source: ", paste(missing, collapse = ", "))
text <- stats::setNames(vapply(files, read_text, character(1)), basename(files))
example_readme_path <- file.path(
  root, "inst", "examples", "riaz-gse91061", "README.md"
)
if (!file.exists(example_readme_path)) fail("Missing Riaz example README.")
example_readme <- read_text(example_readme_path)
all_text <- paste(c(text, example_readme), collapse = "\n")

forbidden_everywhere <- c(
  "delta_A_minus_B",
  "mean_NES_A - mean_NES_B",
  "mean_NES_A-mean_NES_B",
  "schema_version: \"0.5.0\"",
  "build_lisa_report_package()",
  "approved_lisa_dictionary@1.0.0",
  "approved_term2gene@2026.1",
  "approved_category_map@1.0.0"
)
for (needle in forbidden_everywhere) {
  if (grepl(needle, all_text, fixed = TRUE)) {
    fail("Obsolete documentation contract found: ", needle)
  }
}

riaz <- text[["riaz-gse91061-worked-example.Rmd"]]
for (needle in c(
  "Rscript standalone/install_lisaR.R",
  "tables\", \"response_dynamics.tsv", "figures\", \"response_dynamics.png",
  "extdata/riaz-gse91061", "gene_category_contributions.tsv"
)) {
  if (grepl(needle, riaz, fixed = TRUE)) {
    fail("The Riaz guide still depends on an obsolete or non-R/RStudio route: ", needle)
  }
}
for (needle in c("build_lisa_report_package", "compact, precomputed")) {
  if (grepl(needle, example_readme, fixed = TRUE)) {
    fail("The installed Riaz README still exposes an obsolete route: ", needle)
  }
}

# Parse every R chunk, including eval=FALSE. Purl can comment those chunks,
# which would otherwise make an invalid executable example pass this gate.
for (source in files[grepl("[.]Rmd$", files)]) {
  lines <- readLines(source, warn = FALSE)
  inside <- FALSE
  code <- character()
  for (line in lines) {
    if (grepl("^```\\{r([ ,}]|$)", line)) { inside <- TRUE; next }
    if (inside && grepl("^```[[:space:]]*$", line)) {
      inside <- FALSE; code <- c(code, ""); next
    }
    if (inside) code <- c(code, line)
  }
  if (inside) fail("Unclosed R chunk in ", basename(source))
  tryCatch(
    parse(text = code, keep.source = FALSE),
    error = function(error) fail(
      "R code in ", basename(source), " does not parse: ",
      conditionMessage(error)
    )
  )
}

report <- text[["scientific-report.Rmd"]]
required_ab <- c(
  "mean_NES_contextual_A", "display_mean_NES_A", "endpoint_source_A",
  "has_significant_support_A", "plot_has_any_significant_support",
  "gsea_padj_cutoff", "delta_mean_NES"
)
absent_ab <- required_ab[!vapply(
  required_ab, grepl, logical(1), x = report, fixed = TRUE
)]
if (length(absent_ab)) {
  fail("Scientific report omits A/B contract field(s): ", paste(absent_ab, collapse = ", "))
}

readme <- text[["README.md"]]
for (needle in c(
  "install_github", "prepare_lisa_msigdb_resource", "install_lisa_example_bundle",
  "run_lisa", "responders_vs_pd_on", "documentation.html", "lisa-logo.svg"
)) {
  if (!grepl(needle, readme, fixed = TRUE)) {
    fail("README omits required user route: ", needle)
  }
}

# Detailed validation/configuration routes belong in the advanced guides,
# rather than being required in the short repository entry point.
advanced_text <- paste(text[["getting-started.Rmd"]],
                       text[["function-reference.Rmd"]],
                       text[["advanced.Rmd"]], collapse = "\n")
for (needle in c("validate_lisa_config", "plan_lisa_outputs",
                 "verify_lisa_run", "resources-and-provenance", "troubleshooting")) {
  if (!grepl(needle, advanced_text, fixed = TRUE))
    fail("Advanced route omits required content: ", needle)
}
for (stem in c("documentation", "riaz-on-quick-start", "use-your-own-de",
               "read-your-report", "advanced", "statistical-analysis")) {
  if (!paste0(stem, ".Rmd") %in% names(text))
    fail("Missing guided documentation page: ", stem)
}

semantic_text <- paste(
  readme,
  text[["how-lisar-works.Rmd"]],
  text[["resources-and-provenance.Rmd"]],
  text[["TECHNICAL_REBUILD_PROVENANCE.md"]],
  sep = "\n"
)
for (needle in c(
  "biological categor",
  "supercategor",
  "large language models",
  "provenance archive",
  "remain external"
)) {
  if (!contains_ci(needle, semantic_text)) {
    fail("Semantic provenance documentation omits: ", needle)
  }
}
for (needle in c("no longer available", "some of those older records are no")) {
  if (grepl(needle, semantic_text, fixed = TRUE)) {
    fail("Obsolete semantic provenance claim found: ", needle)
  }
}

collection_ids <- c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS", "HALLMARKS")
collection_documents <- c(
  "how-lisar-works.Rmd", "configuration-reference.Rmd",
  "function-reference.Rmd", "custom-dictionaries-and-supercategories.Rmd"
)
missing_collection_documents <- setdiff(collection_documents, names(text))
if (length(missing_collection_documents)) {
  fail(
    "Missing collection-facing documentation: ",
    paste(missing_collection_documents, collapse = ", ")
  )
}
for (document in collection_documents) {
  absent <- collection_ids[!vapply(
    collection_ids, grepl, logical(1), x = text[[document]], fixed = TRUE
  )]
  if (length(absent)) {
    fail(
      document, " does not expose all five collection options: ",
      paste(absent, collapse = ", ")
    )
  }
}

custom <- text[["custom-dictionaries-and-supercategories.Rmd"]]
for (needle in c(
  "lisa_dictionary@2", "category_map@1", "term2gene@1",
  "macrogroup_id", "macrogroup_name", "tier",
  "install_lisa_resource(", "validate_lisa_custom_resources("
)) {
  if (!grepl(needle, custom, fixed = TRUE)) {
    fail("Custom-resource guide omits required contract text: ", needle)
  }
}
for (needle in c("large language model", "multi-LLM", "prompt")) {
  if (contains_ci(needle, custom)) {
    fail("Custom-resource guide must not teach LLM-based creation: ", needle)
  }
}

resource_docs <- paste(
  readme,
  text[["getting-started.Rmd"]],
  text[["configuration-reference.Rmd"]],
  text[["function-reference.Rmd"]],
  text[["resources-and-provenance.Rmd"]],
  sep = "\n"
)
for (needle in c(
  "lisa_dictionary_core@1.0.0", "msigdb_term2gene@2026.1",
  "lisa_category_map@1.0.0"
)) {
  if (!grepl(needle, resource_docs, fixed = TRUE)) {
    fail("Resource documentation omits the current default: ", needle)
  }
}
bundle <- utils::read.delim(file.path(root, "inst", "extdata", "dictionaries",
  "RESOURCE_BUNDLE_MANIFEST.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("lisa_dictionary_core", "lisa_dictionary_expanded", "lisa_category_map") %in% bundle$logical_id) ||
    !all(file.exists(file.path(root, "inst", bundle$artifact)))) {
  fail("Bundled dictionary/hierarchy documentation disagrees with the manifest.")
}
if (!contains_ci("bundled", resource_docs) ||
    !contains_ci("external", resource_docs) ||
    !grepl("prepare_lisa_msigdb_resource", resource_docs, fixed = TRUE)) {
  fail("Document bundled dictionaries/hierarchy and the external membership acquisition route.")
}
for (needle in c("Scientific dictionaries are not bundled",
                 "Scientific resources are not bundled")) {
  if (contains_ci(needle, resource_docs)) fail("Obsolete blanket resource claim: ", needle)
}
for (stem in c("installation", "cptac-ccrcc-proteomics")) {
  if (!file.exists(file.path(root, "vignettes", paste0(stem, ".Rmd"))))
    fail("Missing approved user route: ", stem)
}
active_user_text <- paste(text[grepl("[.]Rmd$", names(text))], readme, collapse = "\n")
if (grepl("interpreting-evidence", active_user_text, fixed = TRUE) ||
    grepl("install-for-evaluators", active_user_text, fixed = TRUE)) {
  fail("A retired vignette is still in the active user route.")
}


message(
  "Documentation contract checks passed for ", length(files) + 1L,
  " source files."
)
