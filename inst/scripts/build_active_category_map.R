#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Derive the active lisaR category map from an immutable source map and the
# exact core/expanded dictionaries. Categories with no assignment in either
# dictionary are removed. Because supercategories exist only through their
# category rows, this also removes any supercategory with no retained category.

parse_args <- function(args) {
  allowed <- c(
    "core", "expanded", "source-map", "source-map-sha256",
    "output", "audit"
  )
  keys_at <- seq.int(1L, length(args), 2L)
  if (length(args) != 2L * length(allowed) ||
      any(!startsWith(args[keys_at], "--"))) {
    stop(
      paste(
        "Usage: build_active_category_map.R",
        "--core FILE --expanded FILE --source-map FILE",
        "--source-map-sha256 HASH --output FILE --audit FILE"
      ),
      call. = FALSE
    )
  }
  keys <- sub("^--", "", args[keys_at])
  if (anyDuplicated(keys) || !setequal(keys, allowed)) {
    stop("All active-map builder arguments must be supplied exactly once.",
         call. = FALSE)
  }
  stats::setNames(args[keys_at + 1L], keys)
}

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character()
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "", fileEncoding = "UTF-8"
  )
}

sha256_file <- function(path) {
  unname(digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ))
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
for (key in c("core", "expanded", "source-map")) {
  cli[[key]] <- normalizePath(cli[[key]], mustWork = TRUE)
}
expected_map_sha256 <- tolower(cli[["source-map-sha256"]])
assert(
  grepl("^[0-9a-f]{64}$", expected_map_sha256),
  "--source-map-sha256 must be one reviewed lower-case SHA-256."
)
assert(
  identical(sha256_file(cli[["source-map"]]), expected_map_sha256),
  "Source category-map SHA-256 differs from the reviewed input."
)

core <- read_tsv(cli[["core"]])
expanded <- read_tsv(cli[["expanded"]])
source_map <- read_tsv(cli[["source-map"]])
required_dictionary <- c("universe", "gene_set_id", "category_id")
required_map <- c(
  "category_id", "display_name", "macrogroup_id", "macrogroup_name",
  "macrogroup_order", "category_order_within_macrogroup"
)
assert(
  !length(setdiff(required_dictionary, names(core))) &&
    !length(setdiff(required_dictionary, names(expanded))),
  "Core or expanded dictionary lacks required identity columns."
)
assert(
  !length(setdiff(required_map, names(source_map))),
  "Source category map lacks required hierarchy columns."
)
assert(!anyDuplicated(source_map$category_id),
       "Source category map contains duplicate category IDs.")

dictionary_categories <- sort(unique(c(
  as.character(core$category_id), as.character(expanded$category_id)
)))
dictionary_categories <- dictionary_categories[
  !is.na(dictionary_categories) & nzchar(dictionary_categories)
]
missing_from_map <- setdiff(dictionary_categories, source_map$category_id)
assert(
  !length(missing_from_map),
  paste0(
    "A dictionary category is absent from the source map: ",
    paste(missing_from_map, collapse = ", ")
  )
)

retained <- source_map$category_id %in% dictionary_categories
active_map <- source_map[retained, , drop = FALSE]
row.names(active_map) <- NULL
assert(nrow(active_map) == length(dictionary_categories),
       "Active map does not contain exactly one row per assigned category.")
assert(all(nzchar(active_map$macrogroup_id)),
       "A retained category lacks a supercategory.")

# A supercategory is represented only by category rows. Verify that every
# retained supercategory therefore has at least one retained category.
supercategory_sizes <- table(active_map$macrogroup_id)
assert(length(supercategory_sizes) > 0L && all(supercategory_sizes >= 1L),
       "A retained supercategory has zero categories.")

audit <- source_map[c("category_id", "display_name", "macrogroup_id")]
audit$source_map_resource <- "lisa_category_map@0.1.0"
audit$active_map_resource <- "lisa_category_map@1.0.0"
audit$core_assignment_rows <- vapply(
  audit$category_id,
  function(id) sum(as.character(core$category_id) == id), integer(1L)
)
audit$expanded_assignment_rows <- vapply(
  audit$category_id,
  function(id) sum(as.character(expanded$category_id) == id), integer(1L)
)
audit$retained <- retained
audit$disposition <- ifelse(
  retained, "retained_assigned_category", "removed_zero_assignment_category"
)

write_tsv(active_map, cli[["output"]])
write_tsv(audit, cli[["audit"]])

cat(sprintf(
  paste0(
    "active category map: %d/%d categories; %d supercategories; ",
    "removed=%d; sha256=%s\n"
  ),
  nrow(active_map), nrow(source_map), length(supercategory_sizes),
  sum(!retained), sha256_file(cli[["output"]])
))
