#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Rebuild the complete private lisaR 0.1 dictionary resources from frozen,
# already-scored ledgers. This script does not reproduce the upstream
# scoring/model-response phase, grant a licence, or publish its outputs.
# BioCarta, KEGG LEGACY and KEGG MEDICUS assignments are retained exactly as
# reviewed in the historical 0.1 resources. The category-map 0.1 bytes remain
# immutable.

parse_args <- function(args) {
  allowed <- c(
    "dictionary-build-dir", "category-palette", "pathways-hierarchy",
    "category-map-sha256", "output-dir", "manifest"
  )
  if (length(args) != 2L * length(allowed) || any(!startsWith(args[seq.int(1L, length(args), 2L)], "--"))) {
    stop(
      paste(
        "Usage: build_compact_dictionary_resources.R",
        "--dictionary-build-dir DIR --category-palette FILE",
        "--pathways-hierarchy FILE --category-map-sha256 HASH",
        "--output-dir DIR --manifest FILE"
      ),
      call. = FALSE
    )
  }
  keys <- sub("^--", "", args[seq.int(1L, length(args), 2L)])
  if (anyDuplicated(keys) || !setequal(keys, allowed)) {
    stop("All compact-builder arguments must be supplied exactly once.", call. = FALSE)
  }
  stats::setNames(args[seq.int(2L, length(args), 2L)], keys)
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
build_dir <- normalizePath(cli[["dictionary-build-dir"]], mustWork = TRUE)
palette_path <- normalizePath(cli[["category-palette"]], mustWork = TRUE)
pathways_meta_path <- normalizePath(cli[["pathways-hierarchy"]], mustWork = TRUE)
expected_category_map_sha256 <- cli[["category-map-sha256"]]
if (!grepl("^[0-9a-f]{64}$", expected_category_map_sha256)) {
  stop(
    "--category-map-sha256 must be the reviewed lower-case digest from the private input manifest.",
    call. = FALSE
  )
}
output_dir <- normalizePath(cli[["output-dir"]], mustWork = FALSE)
manifest_path <- normalizePath(cli[["manifest"]], mustWork = FALSE)

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
}

sha256_file <- function(path) {
  if (!requireNamespace("lisaR", quietly = TRUE)) {
    stop(
      "LISA-SHA256-003 the installed lisaR namespace is required for SHA-256.",
      call. = FALSE
    )
  }
  get("lisa_sha256_file", envir = asNamespace("lisaR"), inherits = FALSE)(path)
}

# Two distinct contracts, deliberately kept apart :
#
#   `required_ledger`  - what the frozen historical ledgers must contain. They
#                        still carry `LISA_score` and are never rewritten, so
#                        this list is unchanged and the historical score check
#                        below still runs against them.
#   `required_runtime` - what the runtime resource emits. `LISA_score` is gone;
#                        it was an exact redundant encoding of `tier`.
required_ledger <- c(
  "universe", "gene_set_id", "gene_set_name", "source_id",
  "category_id", "category_display_name", "LISA_score", "tier"
)
required_runtime <- setdiff(required_ledger, "LISA_score")

source_family <- function(gene_set_id) {
  x <- as.character(gene_set_id)
  out <- rep("CP_OTHER", length(x))
  out[grepl("^GOBP_", x)] <- "GO:BP"
  out[grepl("^GOMF_", x)] <- "GO:MF"
  out[grepl("^GOCC_", x)] <- "GO:CC"
  out[grepl("^BIOCARTA_", x)] <- "BIOCARTA"
  out[grepl("^PID_", x)] <- "PID"
  out[grepl("^REACTOME_", x)] <- "REACTOME"
  out[grepl("^WP_", x)] <- "WIKIPATHWAYS"
  out[grepl("^KEGG_", x)] <- "KEGG_LEGACY"
  out[grepl("^KEGG_MEDICUS_|^N[0-9]", x)] <- "KEGG_MEDICUS"
  out
}

family_order <- c(
  "GO:BP", "GO:MF", "GO:CC", "REACTOME", "WIKIPATHWAYS",
  "KEGG_MEDICUS", "PID", "BIOCARTA", "KEGG_LEGACY", "CP_OTHER"
)
# C6 emits the score-free runtime contract under its own artifact names. Those
# contents are now published as the first LISA dictionaries,
# `lisa_dictionary_core@1.0.0` and `lisa_dictionary_expanded@1.0.0`, and the
# artifact names follow that identity. The v0_1 artifacts are retained
# unmodified in the package as the historical record; this builder no longer
# writes them, so a published identity can never be overwritten.
#
# These names are the bytes the active bundle manifest pins. Changing one here
# without changing RESOURCE_BUNDLE_MANIFEST.tsv would make the rebuild
# unverifiable, so the two are reviewed together.
output_artifact <- c(
  core = "lisa_dictionary_core_runtime_v1_0.tsv",
  expanded = "lisa_dictionary_expanded_runtime_v1_0.tsv"
)

source_files <- character()
runtime_tables <- list()
inventory_tables <- list()
build_tier <- function(tier, stage_dir) {
  suffix <- if (identical(tier, "expanded")) "expanded" else "core"
  base <- file.path(build_dir, paste0("lisa_dictionary_", suffix, "_v0_1.tsv"))
  pathways <- file.path(build_dir, paste0("lisa_pathways_dictionary_", suffix, "_v0_1.tsv"))
  missing_files <- c(base, pathways)[!file.exists(c(base, pathways))]
  if (length(missing_files)) {
    stop("Missing frozen ledger(s): ", paste(basename(missing_files), collapse = ", "), call. = FALSE)
  }
  source_files <<- c(source_files, base, pathways)
  x <- rbind(read_tsv(base), read_tsv(pathways))
  missing_columns <- setdiff(required_ledger, names(x))
  if (length(missing_columns)) {
    stop("Missing ledger column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  x <- x[, required_ledger, drop = FALSE]
  for (column in c("universe", "gene_set_id", "category_id", "tier")) {
    if (any(is.na(x[[column]]) | !nzchar(trimws(as.character(x[[column]]))))) {
      stop("Empty required values in ", column, " for tier ", tier, call. = FALSE)
    }
  }
  # Historical check, retained: the frozen ledger must still carry the scores it
  # was reviewed with. C6 removes the score from the runtime resource, not from
  # the construction record, and never re-scores anything upstream.
  expected_score <- if (identical(tier, "core")) 4 else c(3, 4)
  if (any(!as.numeric(x$LISA_score) %in% expected_score)) {
    stop("Unexpected LISA_score in ", tier, " ledger.", call. = FALSE)
  }
  if (identical(tier, "core") && any(as.character(x$tier) != "core")) {
    stop("Core ledger contains a non-core tier.", call. = FALSE)
  }
  if (identical(tier, "expanded") && any(!as.character(x$tier) %in% c("core", "expanded"))) {
    stop("Expanded ledger contains an unknown tier.", call. = FALSE)
  }
  # The score must be a function of tier, or dropping it would lose information
  # and the runtime projection below would not be equivalence-preserving.
  coherent <- as.numeric(x$LISA_score) ==
    ifelse(as.character(x$tier) == "core", 4, 3)
  if (any(!coherent)) {
    stop("LISA_score is not a function of tier in the ", tier,
         " ledger; the runtime projection would lose information.", call. = FALSE)
  }

  # De-duplicate and order on the LEDGER columns, exactly as before this
  # migration, so row identity and order are decided by the same data they
  # always were. Only then project to the runtime contract. Doing it in this
  # order makes the projection provably incapable of merging rows or permuting
  # them: it removes a column from an already-fixed table.
  x <- unique(x)
  x <- x[order(x$universe, x$gene_set_id, x$category_id), , drop = FALSE]
  rownames(x) <- NULL
  rows_before_projection <- nrow(x)
  x <- x[, required_runtime, drop = FALSE]
  if (!identical(nrow(x), rows_before_projection) || anyDuplicated(x)) {
    stop("Runtime projection changed row identity for tier ", tier,
         "; refusing to emit a resource that is not an exact projection.",
         call. = FALSE)
  }

  families <- source_family(x$gene_set_id)
  present_families <- family_order[family_order %in% unique(families)]
  unexpected_families <- sort(setdiff(unique(families), family_order))
  present_families <- c(present_families, unexpected_families)
  inventory_tables[[tier]] <<- do.call(rbind, lapply(
    present_families,
    function(family) {
      before <- families == family
      data.frame(
        tier = tier,
        source_family = family,
        rows_before = sum(before),
        unique_gene_set_ids_before = length(unique(x$gene_set_id[before])),
        rows_after = sum(before),
        unique_gene_set_ids_after = length(unique(x$gene_set_id[before])),
        excluded_rows = 0L,
        excluded_unique_gene_set_ids = 0L,
        disposition = "retained_in_complete_v0_1_resource",
        stringsAsFactors = FALSE
      )
    }
  ))
  out <- file.path(stage_dir, output_artifact[[tier]])
  write_tsv(x, out)
  runtime_tables[[tier]] <<- x
  out
}

parent <- dirname(output_dir)
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
stage_dir <- tempfile("lisa-compact-build-", tmpdir = parent)
dir.create(stage_dir)
on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)

core_out <- build_tier("core", stage_dir)
expanded_out <- build_tier("expanded", stage_dir)
core_key <- do.call(paste, c(runtime_tables$core[required_runtime], sep = "\r"))
expanded_key <- do.call(paste, c(runtime_tables$expanded[required_runtime], sep = "\r"))
if (!all(core_key %in% expanded_key)) {
  stop("Expanded ledger does not contain every accepted core assignment.", call. = FALSE)
}

category_map <- read_tsv(palette_path)
required_palette <- c(
  "category_id", "display_name", "macrogroup_id", "macrogroup_name",
  "macrogroup_order", "category_order_within_macrogroup", "notes"
)
if (length(setdiff(required_palette, names(category_map)))) {
  stop("Category palette lacks required hierarchy columns.", call. = FALSE)
}
category_map$color <- NULL

pathways <- read_tsv(pathways_meta_path)
required_pathways <- c("family_id", "pathway_id", "pathway_name")
if (length(setdiff(required_pathways, names(pathways)))) {
  stop("PATHWAYS hierarchy lacks required columns.", call. = FALSE)
}
if ("include_in_strict_PATHWAYS" %in% names(pathways)) {
  keep <- toupper(trimws(pathways$include_in_strict_PATHWAYS)) %in% c("SI", "YES", "TRUE", "1")
  pathways <- pathways[keep, , drop = FALSE]
}
pathways <- pathways[!duplicated(pathways$pathway_id), , drop = FALSE]
family_levels <- unique(pathways$family_id)
pathways_map <- data.frame(
  category_id = pathways$pathway_id,
  display_name = pathways$pathway_name,
  macrogroup_id = pathways$family_id,
  macrogroup_name = pathways$family_id,
  macrogroup_order = match(pathways$family_id, family_levels),
  category_order_within_macrogroup = ave(
    seq_len(nrow(pathways)), pathways$family_id, FUN = seq_along
  ),
  notes = paste0(
    "External PATHWAYS hierarchy from approved local source ",
    basename(pathways_meta_path)
  ),
  family_id = pathways$family_id,
  pathway_id = pathways$pathway_id,
  stringsAsFactors = FALSE
)
all_columns <- union(names(category_map), names(pathways_map))
for (column in setdiff(all_columns, names(category_map))) category_map[[column]] <- NA
for (column in setdiff(all_columns, names(pathways_map))) pathways_map[[column]] <- NA
category_map <- rbind(
  category_map[, all_columns, drop = FALSE],
  pathways_map[, all_columns, drop = FALSE]
)
category_map <- category_map[!duplicated(category_map$category_id), , drop = FALSE]
category_map <- category_map[order(
  as.numeric(category_map$macrogroup_order),
  as.numeric(category_map$category_order_within_macrogroup),
  category_map$category_id
), , drop = FALSE]
category_out <- file.path(stage_dir, "lisa_category_map_runtime_v0_1.tsv")
write_tsv(category_map, category_out)
if (!identical(sha256_file(category_out), expected_category_map_sha256)) {
  stop("Category-map 0.1 bytes changed; refusing to build a new version implicitly.", call. = FALSE)
}

inventory <- do.call(rbind, inventory_tables)
row.names(inventory) <- NULL
inventory_out <- file.path(stage_dir, "SOURCE_FAMILY_INVENTORY.tsv")
write_tsv(inventory, inventory_out)

outputs <- c(core_out, expanded_out, category_out)
evidence_outputs <- inventory_out
sources <- unique(c(source_files, palette_path, pathways_meta_path))
manifest <- rbind(
  data.frame(
    role = "source_ledger", artifact = basename(sources),
    rows = vapply(sources, function(path) nrow(read_tsv(path)), integer(1L)),
    bytes = unname(file.info(sources)$size),
    sha256 = vapply(sources, sha256_file, character(1L)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    role = "private_resource_candidate", artifact = basename(outputs),
    rows = vapply(outputs, function(path) nrow(read_tsv(path)), integer(1L)),
    bytes = unname(file.info(outputs)$size),
    sha256 = vapply(outputs, sha256_file, character(1L)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    role = "build_evidence", artifact = basename(evidence_outputs),
    rows = vapply(evidence_outputs, function(path) nrow(read_tsv(path)), integer(1L)),
    bytes = unname(file.info(evidence_outputs)$size),
    sha256 = vapply(evidence_outputs, sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
)

if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be absent or empty: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (path in c(outputs, evidence_outputs)) {
  if (!file.copy(path, file.path(output_dir, basename(path)), copy.mode = FALSE)) {
    stop("Could not promote rebuilt artifact: ", basename(path), call. = FALSE)
  }
}
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(manifest, manifest_path)

cat(sprintf(
  paste0(
    "LISA_COMPACT_DICTIONARY_PASS core=%d expanded=%d categories=%d ",
    "biocarta_core_rows=%d kegg_legacy_core_rows=%d bytes=%d\n"
  ),
  nrow(runtime_tables$core), nrow(runtime_tables$expanded), nrow(category_map),
  sum(inventory$rows_after[
    inventory$tier == "core" & inventory$source_family == "BIOCARTA"
  ]),
  sum(inventory$rows_after[
    inventory$tier == "core" & inventory$source_family == "KEGG_LEGACY"
  ]),
  sum(file.info(file.path(output_dir, basename(outputs)))$size)
))
