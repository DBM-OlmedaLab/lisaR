#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Regenerate the C6 rows of RESOURCE_MIGRATION_MANIFEST.tsv.
#
# Why this is a separate tool from build_resource_migration.R
# -----------------------------------------------------------
# The C5 builder asserts that no artifact's bytes differ from the recorded
# manifest, because C5 was a pure identity relabel. C6 deliberately transforms
# content (it removes the runtime LISA_score column), so that assertion is now
# expected to fail and the C5 builder must not be re-run against a post-C6 tree.
# Its rows are historical and are preserved verbatim by this script.
#
# What a C6 row means
# -------------------
#   change_kind      = "supersede_transform"
#   content_transformed = TRUE
#   old_state        = "retired_transformed"
#
# These rows are deliberately invisible to lisa_superseded_resource_identities(),
# which filters on change_kind == "relabel". They therefore create NO redirect:
# supported exact old pins use the separate code-authenticated historical
# registry and retain their original bytes; these transform rows authorize no
# redirect. Different-byte alias targets are refused by LISA-RESOURCE-030.
#
# Usage:
#   Rscript --vanilla build_resource_migration_c6.R --package-root DIR

parse_args <- function(args) {
  if (length(args) != 2L || !identical(args[[1L]], "--package-root")) {
    stop("usage: build_resource_migration_c6.R --package-root DIR", call. = FALSE)
  }
  normalizePath(args[[2L]], mustWork = TRUE)
}

package_root <- parse_args(commandArgs(trailingOnly = TRUE))
dictionaries <- file.path(package_root, "inst", "extdata", "dictionaries")

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = character()
  )
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
}
sha256_file <- function(path) {
  if (!file.exists(path)) stop("Cannot digest a missing artifact: ", path, call. = FALSE)
  digest::digest(file = path, algo = "sha256")
}

migration_path <- file.path(dictionaries, "RESOURCE_MIGRATION_MANIFEST.tsv")
migration <- read_tsv(migration_path)

required_columns <- c(
  "change_kind", "logical_id", "old_resource_id", "new_resource_id",
  "old_version", "new_version", "species", "modality", "schema",
  "old_sha256", "new_sha256", "old_artifact", "new_artifact",
  "artifact_renamed", "content_transformed", "old_state", "new_state",
  "approved_origin", "rationale"
)
if (!identical(names(migration), required_columns)) {
  stop("RESOURCE_MIGRATION_MANIFEST.tsv does not use the documented schema.",
       call. = FALSE)
}

# Preserve every non-C6 row exactly as found. C5's decisions are history and
# are never recomputed here.
historical <- migration[migration$change_kind != "supersede_transform", , drop = FALSE]
row.names(historical) <- NULL

rationale <- paste(
  "C6 removed the runtime LISA_score column, which was an exact redundant",
  "encoding of tier (core=4, expanded=3, verified over every row). The new",
  "resource is the old one with that single column deleted: same rows, same",
  "assignments, same keys, same multiplicity, same order, same tier. Because the",
  "bytes differ this is NOT a relabel and is deliberately not redirected; the old",
  "identity is retired and its artifact preserved unmodified, so no published",
  "identity can ever resolve to different content."
)

# FROZEN HISTORICAL GENERATOR. This script reproduces the C6 ledger rows as
# they were written at the C6 transition, and is retained as that record.
#
# It is deliberately NOT updated to the first-publication release. Two facts
# have moved since, and the live RESOURCE_MIGRATION_MANIFEST.tsv, not this
# script, is authoritative for both:
#
#   * `new_artifact` below names the files as they were called at C6. Those
#     bytes now live at lisa_dictionary_core_runtime_v1_0.tsv and
#     lisa_dictionary_expanded_runtime_v1_0.tsv. Nothing about the bytes
#     changed; only the filenames did.
#   * `new_version` 2.0.0 records the identity C6 assigned. Those identities
#     are now superseded compatibility aliases for
#     lisa_dictionary_core@1.0.0 and lisa_dictionary_expanded@1.0.0.
#
# Re-running this script would therefore regenerate the C6-era rows, not the
# current ledger. Do that only to audit C6 itself.
c6_map <- data.frame(
  logical_id = c("lisa_core", "lisa_expanded"),
  old_version = c("1.0.0", "1.0.0"),
  new_version = c("2.0.0", "2.0.0"),
  old_artifact = c(
    "extdata/dictionaries/lisa_core_runtime_v0_1.tsv",
    "extdata/dictionaries/lisa_expanded_runtime_v0_1.tsv"
  ),
  new_artifact = c(
    "extdata/dictionaries/lisa_core_runtime_v2_0.tsv",
    "extdata/dictionaries/lisa_expanded_runtime_v2_0.tsv"
  ),
  stringsAsFactors = FALSE
)

artifact_path <- function(artifact) file.path(package_root, "inst", artifact)

# The retired artifact must still be present and unmodified: the whole point of
# retiring rather than deleting is that the historical bytes remain auditable.
rows <- lapply(seq_len(nrow(c6_map)), function(i) {
  row <- c6_map[i, , drop = FALSE]
  old_path <- artifact_path(row$old_artifact)
  new_path <- artifact_path(row$new_artifact)
  old_sha <- sha256_file(old_path)
  new_sha <- sha256_file(new_path)
  if (identical(old_sha, new_sha)) {
    stop(sprintf(
      "'%s' old and new artifacts are byte-identical; record a relabel, not a transform.",
      row$logical_id
    ), call. = FALSE)
  }

  # Assert the transformation really is the declared one: the new table must be
  # the old table with exactly the LISA_score column removed, same row order.
  old_table <- read_tsv(old_path)
  new_table <- read_tsv(new_path)
  expected <- setdiff(names(old_table), "LISA_score")
  if (!"LISA_score" %in% names(old_table) ||
      !identical(names(new_table), expected) ||
      !identical(nrow(new_table), nrow(old_table))) {
    stop(sprintf("'%s' is not an exact LISA_score projection.", row$logical_id),
         call. = FALSE)
  }
  for (column in expected) {
    if (!identical(new_table[[column]], old_table[[column]])) {
      stop(sprintf("'%s' column '%s' changed during the projection.",
                   row$logical_id, column), call. = FALSE)
    }
  }

  data.frame(
    change_kind = "supersede_transform",
    logical_id = row$logical_id,
    old_resource_id = paste0(row$logical_id, "@", row$old_version),
    new_resource_id = paste0(row$logical_id, "@", row$new_version),
    old_version = row$old_version, new_version = row$new_version,
    species = "Homo sapiens", modality = "all",
    schema = "lisa_dictionary@2",
    old_sha256 = old_sha, new_sha256 = new_sha,
    old_artifact = row$old_artifact, new_artifact = row$new_artifact,
    artifact_renamed = "TRUE", content_transformed = "TRUE",
    old_state = "retired_transformed", new_state = "active",
    approved_origin = paste(
      "LISA project candidate; CSIC-owned classification;",
      "profile-independent gene-set/category relationships; RESOURCE_PROVENANCE.tsv"
    ),
    rationale = rationale,
    stringsAsFactors = FALSE
  )
})

out <- rbind(historical, do.call(rbind, rows))
row.names(out) <- NULL
write_tsv(out, migration_path)

cat("C6 migration rows regenerated:", nrow(out) - nrow(historical),
    "| preserved historical rows:", nrow(historical), "\n")
cat("manifest:", migration_path, "\n")
