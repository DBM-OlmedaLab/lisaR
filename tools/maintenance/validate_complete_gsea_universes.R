#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Validate the permanent complete-universe GSEA contract against exact
# TERM2GENE, dictionary and active category-map resources.

parse_args <- function(args) {
  allowed <- c(
    "package-root", "term2gene", "core", "expanded", "category-map",
    "output-dir"
  )
  index <- seq.int(1L, length(args), 2L)
  if (length(args) != 2L * length(allowed) ||
      any(!startsWith(args[index], "--"))) {
    stop(
      paste(
        "Usage: validate_complete_gsea_universes.R",
        "--package-root DIR --term2gene FILE --core FILE --expanded FILE",
        "--category-map FILE --output-dir DIR"
      ),
      call. = FALSE
    )
  }
  keys <- sub("^--", "", args[index])
  if (anyDuplicated(keys) || !setequal(keys, allowed)) {
    stop("All universe-validator arguments must be supplied exactly once.",
         call. = FALSE)
  }
  stats::setNames(args[index + 1L], keys)
}

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character()
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "", fileEncoding = "UTF-8"
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
for (key in c("package-root", "term2gene", "core", "expanded", "category-map")) {
  cli[[key]] <- normalizePath(cli[[key]], mustWork = TRUE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The maintenance validator requires pkgload.", call. = FALSE)
}
pkgload::load_all(cli[["package-root"]], quiet = TRUE, export_all = TRUE)

output_dir <- normalizePath(cli[["output-dir"]], mustWork = FALSE)
if (dir.exists(output_dir) &&
    length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be absent or empty: ", output_dir,
       call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

term2gene <- read_tsv(cli[["term2gene"]])
dictionaries <- list(
  core = read_tsv(cli[["core"]]),
  expanded = read_tsv(cli[["expanded"]])
)
category_map <- read_tsv(cli[["category-map"]])
collections <- c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
c2_sources <- c(
  "BIOCARTA", "KEGG_LEGACY", "KEGG_MEDICUS", "PID", "REACTOME",
  "WIKIPATHWAYS", "CP_OTHER"
)

dictionary_categories <- sort(unique(unlist(lapply(
  dictionaries, function(x) as.character(x$category_id)
), use.names = FALSE)))
assert(!anyDuplicated(category_map$category_id),
       "Active category map contains duplicate category IDs.")
assert(setequal(category_map$category_id, dictionary_categories),
       "Active category map differs from assigned dictionary categories.")
assert(all(nzchar(category_map$macrogroup_id)),
       "Active category map contains an empty supercategory ID.")
assert(all(table(category_map$macrogroup_id) >= 1L),
       "An active supercategory has no category.")

coverage_rows <- list()
residual_rows <- list()
for (tier in names(dictionaries)) {
  dictionary <- dictionaries[[tier]]
  for (universe in collections) {
    selected <- lisaR:::select_lisa_term2gene_universe(
      term2gene_raw = term2gene,
      universe = universe,
      c2_sources = c2_sources,
      dictionary_gene_sets = unique(dictionary$gene_set_id[
        dictionary$universe == universe
      ])
    )
    evaluated <- sort(unique(as.character(selected$gs_name)))
    assigned <- sort(unique(as.character(dictionary$gene_set_id[
      dictionary$universe == universe
    ])))
    absent_from_source <- setdiff(assigned, evaluated)
    assert(
      !length(absent_from_source),
      paste0(tier, "/", universe,
             " dictionary sets are absent from the selected source space.")
    )
    residual <- setdiff(evaluated, assigned)
    assert(
      !length(intersect(assigned, residual)) &&
        setequal(evaluated, c(assigned, residual)),
      paste0(tier, "/", universe,
             " classified/residual partition is not exact.")
    )
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      tier = tier,
      universe = universe,
      evaluated_gene_sets = length(evaluated),
      classified_gene_sets = length(assigned),
      unclassified_gene_sets = length(residual),
      biological_coverage_pct = round(
        100 * length(assigned) / length(evaluated), 2
      ),
      gsea_accounting_coverage_pct = 100,
      residual_table_policy = if (universe %in%
                                  c("GOBP-C2", "GOMF", "GOCC")) {
        "dedicated_table"
      } else {
        "ledger_only"
      },
      plot_policy = "classified_only",
      stringsAsFactors = FALSE
    )
    residual_rows[[length(residual_rows) + 1L]] <- data.frame(
      tier = tier,
      universe = universe,
      gene_set_id = residual,
      classification_status = "unclassified",
      classification_bucket = "OTHER_UNCLASSIFIED",
      stringsAsFactors = FALSE
    )
  }
}

coverage <- do.call(rbind, coverage_rows)
residuals <- do.call(rbind, residual_rows)
row.names(coverage) <- NULL
row.names(residuals) <- NULL
write_tsv(coverage, file.path(output_dir, "UNIVERSE_COVERAGE.tsv"))
write_tsv(residuals, file.path(output_dir, "RESIDUAL_GENE_SETS.tsv"))
write_tsv(data.frame(
  metric = c(
    "active_categories", "active_supercategories",
    "zero_assignment_categories", "empty_supercategories"
  ),
  value = c(
    nrow(category_map), length(unique(category_map$macrogroup_id)), 0L, 0L
  ),
  status = "PASS",
  stringsAsFactors = FALSE
), file.path(output_dir, "CATEGORY_MAP_VALIDATION.tsv"))

print(coverage, row.names = FALSE)
cat(sprintf(
  "COMPLETE_GSEA_UNIVERSE_VALIDATION=PASS categories=%d supercategories=%d\n",
  nrow(category_map), length(unique(category_map$macrogroup_id))
))
