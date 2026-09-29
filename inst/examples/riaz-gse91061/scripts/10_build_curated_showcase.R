#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Build a small, reviewable showcase from a verified standard Riaz run.
# This script copies selected existing outputs. It does not rerun differential
# expression, enrichment or LISA annotation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste(
      "Usage: 10_build_curated_showcase.R",
      "<verified-lisa-run> <new-output-directory>"
    ),
    call. = FALSE
  )
}
if (!requireNamespace("lisaR", quietly = TRUE)) {
  stop("lisaR 1.x (>= 1.0.0 and < 2.0.0) is required.", call. = FALSE)
}
active_version <- utils::packageVersion("lisaR")
if (active_version < "1.0.0" || active_version >= "2.0.0") {
  stop(
    "lisaR 1.x (>= 1.0.0 and < 2.0.0) is required; found ", active_version, ".",
    call. = FALSE
  )
}

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

source_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_root <- riaz_new_output_path(args[[2]], "showcase destination")
if (riaz_path_entry_exists(output_root)) {
  stop(
    "The output path already exists; choose a new directory: ",
    output_root,
    call. = FALSE
  )
}
if (riaz_path_is_within(output_root, source_root)) {
  stop(
    paste(
      "The showcase destination must be outside the verified source run;",
      "choose a new sibling directory."
    ),
    call. = FALSE
  )
}

verification <- lisaR::verify_lisa_run(source_root)
if (!identical(verification$gate, "PASS") ||
    !identical(verification$artifact_type, "scientific_run")) {
  stop(
    "The source must verify as a PASS scientific_run before selection.",
    call. = FALSE
  )
}

analyses <- c("differential_longitudinal_response", "responders_on_vs_pre")
collections <- c("PATHWAYS", "GOBP-C2")
max_files <- 2000L
max_bytes <- 250 * 1024^2

staging <- paste0(output_root, ".staging")
if (riaz_path_entry_exists(staging)) {
  stop("Staging path already exists: ", staging, call. = FALSE)
}
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
completed <- FALSE
on.exit(
  if (!completed && dir.exists(staging)) {
    unlink(staging, recursive = TRUE, force = TRUE)
  },
  add = TRUE
)

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
}

sha256_file <- function(path) {
  get(
    "lisa_sha256_file", envir = asNamespace("lisaR"), inherits = FALSE
  )(path)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

find_one <- function(directory, pattern, label) {
  hits <- sort(list.files(directory, pattern = pattern, full.names = TRUE))
  if (length(hits) != 1L) {
    stop(
      label, " must resolve to exactly one file; found ", length(hits),
      " under ", directory, ".",
      call. = FALSE
    )
  }
  hits[[1]]
}

de_index_path <- file.path(source_root, "config", "de_index.tsv")
contrast_index_path <- file.path(source_root, "config", "contrast_index.tsv")
contrast_status_path <- file.path(source_root, "contrast_status.tsv")
required_indexes <- c(de_index_path, contrast_index_path, contrast_status_path)
if (!all(file.exists(required_indexes))) {
  stop("The verified run is missing a required analysis/contrast index.", call. = FALSE)
}
de_index <- read_tsv(de_index_path)
contrast_index <- read_tsv(contrast_index_path)
contrast_status <- read_tsv(contrast_status_path)
if (!"analysis_id" %in% names(de_index) ||
    !all(analyses %in% de_index$analysis_id)) {
  stop("The selected showcase analyses are absent from config/de_index.tsv.", call. = FALSE)
}
required_contrast_columns <- c("contrast_id", "output_id")
if (!all(required_contrast_columns %in% names(contrast_index)) ||
    !all(c(required_contrast_columns, "collection", "status") %in%
      names(contrast_status))) {
  stop("The contrast index/status contract is incomplete.", call. = FALSE)
}

summary_paths <- list()
selection <- list()
for (analysis in analyses) {
  for (collection in collections) {
    summary_dir <- file.path(
      source_root, "outputs", "single_de", analysis,
      paste0("collection_", collection), "lisa_tables"
    )
    summary_path <- find_one(
      summary_dir, "_GSEA_category_summary[.]tsv$", "Category summary"
    )
    key <- paste(analysis, collection, sep = "\r")
    summary_paths[[key]] <- summary_path
    tab <- read_tsv(summary_path)
    required <- c(
      "category_id", "category_display_name", "n_genesets_mapped",
      "n_genesets_evaluable", "n_genesets_significant", "mean_NES",
      "min_padj"
    )
    if (!all(required %in% names(tab))) {
      stop("Category summary lacks the lisaR 0.6 support fields: ", summary_path,
           call. = FALSE)
    }
    tab <- tab[
      is.finite(tab$min_padj) & tab$min_padj <= 0.05 &
        is.finite(tab$mean_NES) & tab$n_genesets_evaluable > 0L &
        tab$n_genesets_significant > 0L,
      , drop = FALSE
    ]
    choose_direction <- function(direction) {
      part <- tab[direction * tab$mean_NES > 0, , drop = FALSE]
      part <- part[
        order(part$min_padj, -abs(part$mean_NES), part$category_id),
        , drop = FALSE
      ]
      utils::head(part, 2L)
    }
    chosen <- rbind(choose_direction(1), choose_direction(-1))
    if (nrow(chosen) == 0L) {
      stop(
        "No supported positive or negative category found for ", analysis,
        " / ", collection, ".",
        call. = FALSE
      )
    }
    chosen$analysis_id <- analysis
    chosen$collection <- collection
    chosen$selection_rule <- paste(
      "up to two positive and two negative categories; min_padj<=0.05;",
      "n_genesets_evaluable>0; n_genesets_significant>0; ordered by",
      "min_padj, abs(mean_NES), category_id"
    )
    selection[[length(selection) + 1L]] <- chosen
  }
}
selection <- do.call(rbind, selection)
selection <- selection[, c(
  "analysis_id", "collection", "category_id", "category_display_name",
  "n_genesets_mapped", "n_genesets_evaluable", "n_genesets_significant",
  "mean_NES", "min_padj", "selection_rule"
)]
write_tsv(selection, file.path(staging, "selection.tsv"))
write_tsv(
  data.frame(
    source_run = source_root,
    gate = verification$gate,
    artifact_type = verification$artifact_type,
    stringsAsFactors = FALSE
  ),
  file.path(staging, "source_verification.tsv")
)

copied <- list()
relative_to_source <- function(path) {
  prefix <- paste0(source_root, "/")
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!startsWith(normalized, prefix)) {
    stop(
      "A selected showcase source escaped the verified run: ", normalized,
      call. = FALSE
    )
  }
  substring(normalized, nchar(prefix) + 1L)
}

copy_one <- function(source, destination, kind, analysis = "",
                     collection = "", category = "", contrast = "") {
  if (!file.exists(source)) return(invisible(FALSE))
  destination <- riaz_safe_relative_path(
    destination, "showcase inventory destination"
  )
  target <- file.path(staging, destination)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, target, overwrite = FALSE, copy.mode = TRUE,
                 copy.date = TRUE)) {
    stop("Could not copy: ", source, call. = FALSE)
  }
  copied[[length(copied) + 1L]] <<- data.frame(
    path = destination,
    kind = kind,
    analysis_id = analysis,
    contrast_id = contrast,
    collection = collection,
    category_id = category,
    source_path = relative_to_source(source),
    stringsAsFactors = FALSE
  )
  invisible(TRUE)
}

copy_matching <- function(directory, pattern, destination_dir, kind, analysis,
                          collection, category = "", contrast = "") {
  if (!dir.exists(directory)) return(0L)
  matches <- sort(list.files(directory, pattern = pattern, full.names = TRUE))
  for (source in matches) {
    copy_one(
      source, file.path(destination_dir, basename(source)), kind,
      analysis, collection, category, contrast
    )
  }
  length(matches)
}

for (key in names(summary_paths)) {
  parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
  copy_one(
    summary_paths[[key]],
    file.path("analyses", parts[[1]], parts[[2]], "category_summary.tsv"),
    "category_summary", analysis = parts[[1]], collection = parts[[2]]
  )
}

for (i in seq_len(nrow(selection))) {
  analysis <- selection$analysis_id[[i]]
  collection <- selection$collection[[i]]
  category <- selection$category_id[[i]]
  plot_root <- file.path(
    source_root, "outputs", "single_de", analysis,
    paste0("collection_", collection), "plots"
  )
  member_dirs <- list.dirs(plot_root, recursive = FALSE, full.names = TRUE)
  member_dirs <- member_dirs[grepl("_GSEA_category_pathways$", member_dirs)]
  destination <- file.path(
    "analyses", analysis, collection, "categories", category,
    "member_gene_sets"
  )
  pattern <- paste0(
    "^[0-9]+_", regex_escape(category), "(_source[.]tsv|[.]png)$"
  )
  copied_members <- sum(vapply(member_dirs, function(directory) {
    copy_matching(
      directory, pattern, destination, "member_gene_set", analysis,
      collection, category
    )
  }, integer(1)))
  if (copied_members == 0L) {
    stop(
      "No standard member-gene-set output found for ", analysis, " / ",
      collection, " / ", category, ".",
      call. = FALSE
    )
  }
}

completed_contrasts <- contrast_status[
  contrast_status$status == "completed" &
    contrast_status$collection %in% collections,
  c("contrast_id", "output_id", "collection"), drop = FALSE
]
completed_contrasts <- unique(completed_contrasts)
if (nrow(completed_contrasts) == 0L) {
  stop("No completed configured contrast is available for the showcase.", call. = FALSE)
}
known_pairs <- paste(
  contrast_index$contrast_id, contrast_index$output_id, sep = "\r"
)
status_pairs <- paste(
  completed_contrasts$contrast_id, completed_contrasts$output_id, sep = "\r"
)
if (!all(status_pairs %in% known_pairs)) {
  stop("contrast_status.tsv contains an entry absent from contrast_index.tsv.",
       call. = FALSE)
}

for (i in seq_len(nrow(completed_contrasts))) {
  row <- completed_contrasts[i, , drop = FALSE]
  output_name <- paste(row$contrast_id[[1]], row$output_id[[1]], sep = "_")
  contrast_base <- file.path(
    source_root, "outputs", "category_contrasts", output_name,
    paste0("collection_", row$collection[[1]])
  )
  table_path <- find_one(
    file.path(contrast_base, "lisa_tables"),
    "_GSEA_category_contrast[.]tsv$", "Contrast table"
  )
  plot_path <- find_one(
    file.path(contrast_base, "plots"),
    "_GSEA_contrast_dumbbell_all_annotated[.]png$", "Annotated contrast plot"
  )
  plot_source <- find_one(
    file.path(contrast_base, "plots"),
    "_GSEA_contrast_dumbbell_all_annotated_source[.]tsv$",
    "Annotated contrast source table"
  )
  destination <- file.path(
    "contrasts", row$contrast_id[[1]], row$output_id[[1]],
    row$collection[[1]]
  )
  for (path in c(table_path, plot_path, plot_source)) {
    copy_one(
      path, file.path(destination, basename(path)),
      if (grepl("[.]png$", path)) "contrast_plot" else "contrast_source",
      collection = row$collection[[1]], contrast = row$contrast_id[[1]]
    )
  }
}

if (length(copied) == 0L) {
  stop("No showcase files were copied.", call. = FALSE)
}
inventory <- do.call(rbind, copied)
absolute <- file.path(staging, inventory$path)
inventory$bytes <- as.numeric(file.info(absolute)$size)
inventory$sha256 <- vapply(absolute, sha256_file, character(1))
inventory <- inventory[order(inventory$path), , drop = FALSE]
write_tsv(inventory, file.path(staging, "inventory.tsv"))

png_rows <- inventory[grepl("[.]png$", inventory$path), , drop = FALSE]
cards <- vapply(seq_len(nrow(png_rows)), function(i) {
  rel <- png_rows$path[[i]]
  label <- paste(
    Filter(nzchar, c(
      png_rows$analysis_id[[i]], png_rows$contrast_id[[i]],
      png_rows$collection[[i]], png_rows$category_id[[i]]
    )),
    collapse = " · "
  )
  paste0(
    "<article><a href=\"", html_escape(rel), "\"><img src=\"",
    html_escape(rel), "\" alt=\"", html_escape(label),
    "\" loading=\"lazy\"></a><p><strong>",
    html_escape(png_rows$kind[[i]]), "</strong><br>",
    html_escape(label), "</p></article>"
  )
}, character(1))
html <- c(
  "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
  "<title>lisaR Riaz curated showcase</title><style>body{font-family:system-ui,sans-serif;margin:0;color:#20242b;background:#f5f3ef}header,main{max-width:1200px;margin:auto;padding:24px}header{background:#fff;border-bottom:1px solid #ddd}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:18px}article{background:#fff;border:1px solid #ddd;border-radius:10px;padding:12px}img{width:100%;height:260px;object-fit:contain;background:#fff}p{line-height:1.45}</style></head><body>",
  "<header><h1>lisaR Riaz curated showcase</h1>",
  "<p>This is a deterministic selection copied from a verified lisaR 0.6 standard scientific run. No analysis was recalculated.</p>",
  "<p>It is a local review candidate, not an approved public distribution. See <a href=\"source_verification.tsv\">source verification</a>, <a href=\"selection.tsv\">selection</a>, and <a href=\"inventory.tsv\">inventory and hashes</a>.</p></header>",
  "<main><div class=\"grid\">", cards, "</div></main></body></html>"
)
writeLines(html, file.path(staging, "index.html"), useBytes = TRUE)

all_files <- list.files(
  staging, recursive = TRUE, full.names = TRUE, all.files = TRUE,
  no.. = TRUE
)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
total_bytes <- sum(as.numeric(file.info(all_files)$size))
if (length(all_files) > max_files || total_bytes > max_bytes) {
  stop(
    "Showcase budget exceeded: ", length(all_files), " files and ",
    total_bytes, " bytes.",
    call. = FALSE
  )
}

if (riaz_path_entry_exists(output_root)) {
  stop(
    "The showcase destination appeared before promotion; staging is preserved.",
    call. = FALSE
  )
}
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
if (!file.rename(staging, output_root)) {
  stop("Could not promote the completed showcase.", call. = FALSE)
}
completed <- TRUE
cat(
  "Created ", output_root, " with ", length(all_files), " files and ",
  total_bytes, " bytes.\nRIAZ_SHOWCASE_GATE=PASS\n",
  sep = ""
)
