#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# C5 resource-version migration builder.
#
# Purpose: derive the migrated `RESOURCE_BUNDLE_MANIFEST.tsv` and the
# machine-readable `RESOURCE_MIGRATION_MANIFEST.tsv` from the declared version
# map plus the bytes that are actually on disk. SHA-256 values are always
# recomputed from the artifacts, never copied from a previous manifest and
# never typed by hand, so a mistyped digest cannot enter the registry.
#
# The script is deliberately idempotent: running it again on an already
# migrated tree reproduces byte-identical outputs.
#
# Usage:
#   Rscript --vanilla tools/maintenance/build_resource_migration.R \
#     --package-root . [--check]
#
#   --check  Do not write; fail if the on-disk manifests differ from the
#            derived ones. Used by the installed contract test.

args <- commandArgs(trailingOnly = TRUE)

cli_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[[1L]] == length(args)) {
    stop(sprintf("Missing value for %s", flag), call. = FALSE)
  }
  args[[hit[[1L]] + 1L]]
}

package_root <- cli_value("--package-root", ".")
check_only <- "--check" %in% args

dictionaries <- file.path(package_root, "inst", "extdata", "dictionaries")
if (!dir.exists(dictionaries)) {
  stop("Not a lisaR source tree: ", dictionaries, call. = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = character()
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
    na = ""
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot digest a missing artifact: ", path, call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256")
}

# ---------------------------------------------------------------------------
# Declared C5 version map.
#
# `change_kind` is the auditable distinction the migration contract requires:
#   relabel  - the same bytes keep their content and gain a new resource
#              version. No file is renamed and no content is transformed.
#   transform- the bytes themselves change. No C5 resource uses this; it is
#              declared so the schema can express it without a later edit.
#   hold     - the resource deliberately keeps its current version. The
#              `rationale` records why, including refused collisions.
#
# Every `relabel` row registers its previous identity as a superseded
# compatibility identity: old pins keep resolving, to the very same bytes,
# and nobody else may re-mint that identity with different content.
# ---------------------------------------------------------------------------
version_map <- rbind(
  data.frame(
    change_kind = "relabel", logical_id = "lisa_core",
    old_version = "0.1.0", new_version = "1.0.0",
    rationale = paste(
      "First-publication version for the LISA project's own core dictionary.",
      "lisa_core has carried exactly one version (0.1.0) in the whole recorded",
      "history and lisa_core@1.0.0 has never existed, so 1.0.0 is free for this",
      "logical identity. Bytes, assignments, tier and LISA_score are unchanged."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    change_kind = "relabel", logical_id = "lisa_expanded",
    old_version = "0.1.0", new_version = "1.0.0",
    rationale = paste(
      "First-publication version for the expanded dictionary. Same identity",
      "history argument as lisa_core; lisa_expanded@1.0.0 has never existed.",
      "Bytes and assignments unchanged."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    change_kind = "relabel", logical_id = "lisa_category_map",
    old_version = "0.1.1", new_version = "1.0.0",
    rationale = paste(
      "First-publication version for the active category map. The historical",
      "external map lisa_category_map@0.1.0 stays private, byte-preserved and",
      "unreferenced as a bundled resource; lisa_category_map@1.0.0 has never",
      "existed. Category labels, order and supercategories are unchanged."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    change_kind = "hold", logical_id = "lisa_quickstart_dictionary",
    old_version = "3.0.0", new_version = "3.0.0",
    rationale = paste(
      "REFUSED 1.0.0. lisa_quickstart_dictionary@1.0.0 already exists for this",
      "same logical identity with different bytes",
      "(extdata/quick-start/example_dictionary.tsv). Re-minting it would point a",
      "published identity at different content. The fixture also teaches",
      "version migration, so its ascending versions are meaningful, not cosmetic."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    change_kind = "hold", logical_id = "lisa_quickstart_term2gene",
    old_version = "2.0.0", new_version = "2.0.0",
    rationale = paste(
      "REFUSED 1.0.0. lisa_quickstart_term2gene@1.0.0 already exists for this",
      "same logical identity with different bytes",
      "(extdata/quick-start/example_term2gene.tsv)."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    change_kind = "hold", logical_id = "lisa_quickstart_category_map",
    old_version = "1.0.0", new_version = "1.0.0",
    rationale = paste(
      "Already at its first-publication version; nothing to migrate."
    ),
    stringsAsFactors = FALSE
  )
)

bundle_path <- file.path(dictionaries, "RESOURCE_BUNDLE_MANIFEST.tsv")
migration_path <- file.path(dictionaries, "RESOURCE_MIGRATION_MANIFEST.tsv")

# C6 guard. This script is the historical C5 builder and asserts, below, that no
# artifact's bytes ever change. C6 deliberately transformed the dictionary
# artifacts (it removed the runtime LISA_score column), so re-running this
# script against a post-C6 tree would both fail that assertion and silently drop
# the C6 rows from the manifest. Refuse explicitly instead of failing obscurely.
if (file.exists(migration_path)) {
  existing <- read_tsv(migration_path)
  if ("change_kind" %in% names(existing) &&
      any(existing$change_kind == "supersede_transform")) {
    stop(
      paste(
        "This tree already contains C6 transform rows. build_resource_migration.R",
        "is the C5 identity-relabel builder and must not be re-run here; it would",
        "drop them. Use tools/maintenance/build_resource_migration_c6.R."
      ),
      call. = FALSE
    )
  }
}

bundle <- read_tsv(bundle_path)
required_columns <- c(
  "logical_id", "version", "species", "modality", "schema", "sha256",
  "compatibility", "approved_origin", "artifact"
)
if (!identical(names(bundle), required_columns)) {
  stop(
    "RESOURCE_BUNDLE_MANIFEST.tsv does not use the closed registry schema.",
    call. = FALSE
  )
}

package_dir <- file.path(package_root, "inst")
artifact_path <- function(artifact) {
  file.path(package_dir, sub("^extdata/", "extdata/", artifact))
}

# The migration is expressed over identities, not over row order. Resolve each
# manifest row against the declared map; an unmapped row is an error rather
# than a silent pass-through, so a newly added resource cannot slip through
# the migration without an explicit decision.
resolve_new_version <- function(logical_id, version) {
  hit <- which(version_map$logical_id == logical_id)
  if (!length(hit)) {
    stop(sprintf(
      "Resource '%s' has no C5 decision. Add an explicit row to version_map.",
      logical_id
    ), call. = FALSE)
  }
  row <- version_map[hit[[1L]], , drop = FALSE]
  if (identical(version, row$new_version)) return(version)      # idempotent
  if (identical(version, row$old_version)) return(row$new_version)
  # A historical variant that the map does not move (for example the retained
  # quick-start 2.0.0 dictionary) keeps its version untouched.
  version
}

migrated <- bundle
migrated$version <- vapply(seq_len(nrow(bundle)), function(i) {
  resolve_new_version(bundle$logical_id[[i]], bundle$version[[i]])
}, character(1))

# Recompute every digest from the artifact that the row points at. A relabel
# must not change any digest; that invariant is asserted, not assumed.
observed <- vapply(migrated$artifact, function(a) sha256_file(artifact_path(a)),
                   character(1), USE.NAMES = FALSE)
changed <- which(observed != bundle$sha256)
if (length(changed)) {
  stop(sprintf(
    paste0(
      "Artifact bytes differ from the recorded manifest for '%s@%s'. C5 is an ",
      "identity migration and must not transform content; investigate before ",
      "rebuilding."
    ),
    bundle$logical_id[[changed[[1L]]]], bundle$version[[changed[[1L]]]]
  ), call. = FALSE)
}
migrated$sha256 <- observed
migrated <- migrated[order(
  migrated$logical_id, migrated$version, migrated$species, migrated$modality,
  method = "radix"
), , drop = FALSE]
row.names(migrated) <- NULL

# ---------------------------------------------------------------------------
# Machine-readable migration manifest.
# ---------------------------------------------------------------------------
manifest_rows <- lapply(seq_len(nrow(version_map)), function(i) {
  row <- version_map[i, , drop = FALSE]
  # Before the migration the row still carries the old version; afterwards it
  # carries the new one. Accept either so the script stays idempotent and can
  # re-verify an already migrated tree.
  source_rows <- bundle[
    bundle$logical_id == row$logical_id &
      bundle$version %in% c(row$old_version, row$new_version), ,
    drop = FALSE
  ]
  if (!nrow(source_rows)) {
    stop(sprintf(
      "Declared migration for '%s@%s' has no manifest row.",
      row$logical_id, row$old_version
    ), call. = FALSE)
  }
  if (nrow(source_rows) > 1L) {
    # Both identities registered at once would mean two competing active
    # variants, which the migration explicitly rules out.
    stop(sprintf(
      "Resource '%s' has both %s and %s registered; the migration expects one active variant.",
      row$logical_id, row$old_version, row$new_version
    ), call. = FALSE)
  }
  source_row <- source_rows[1L, , drop = FALSE]
  sha <- sha256_file(artifact_path(source_row$artifact))
  relabel <- identical(row$change_kind, "relabel")
  data.frame(
    change_kind = row$change_kind,
    logical_id = row$logical_id,
    old_resource_id = paste0(row$logical_id, "@", row$old_version),
    new_resource_id = paste0(row$logical_id, "@", row$new_version),
    old_version = row$old_version,
    new_version = row$new_version,
    species = source_row$species,
    modality = source_row$modality,
    schema = source_row$schema,
    old_sha256 = sha,
    new_sha256 = sha,
    old_artifact = source_row$artifact,
    new_artifact = source_row$artifact,
    artifact_renamed = "FALSE",
    content_transformed = "FALSE",
    old_state = if (relabel) "superseded_compatibility" else "active",
    new_state = "active",
    approved_origin = source_row$approved_origin,
    rationale = row$rationale,
    stringsAsFactors = FALSE
  )
})
migration <- do.call(rbind, manifest_rows)
migration <- migration[order(migration$logical_id, method = "radix"), ,
                       drop = FALSE]
row.names(migration) <- NULL

if (any(migration$old_sha256 != migration$new_sha256)) {
  stop("A C5 relabel changed bytes; refusing to write.", call. = FALSE)
}

if (check_only) {
  compare <- function(actual, expected, label) {
    tmp_a <- tempfile(); tmp_b <- tempfile()
    on.exit(unlink(c(tmp_a, tmp_b)), add = TRUE)
    write_tsv(actual, tmp_a)
    write_tsv(expected, tmp_b)
    if (!identical(sha256_file(tmp_a), sha256_file(tmp_b))) {
      stop(sprintf("%s is not the derived migration output.", label),
           call. = FALSE)
    }
  }
  compare(bundle[order(bundle$logical_id, bundle$version, bundle$species,
                       bundle$modality, method = "radix"), , drop = FALSE],
          migrated, "RESOURCE_BUNDLE_MANIFEST.tsv")
  compare(read_tsv(migration_path), migration, "RESOURCE_MIGRATION_MANIFEST.tsv")
  cat("LISA_RESOURCE_MIGRATION_CHECK_PASS\n")
  quit(save = "no", status = 0L)
}

write_tsv(migrated, bundle_path)
write_tsv(migration, migration_path)

cat(sprintf(
  "LISA_RESOURCE_MIGRATION_PASS registry_rows=%d migrated_identities=%d relabelled=%d held=%d\n",
  nrow(migrated), nrow(migration),
  sum(migration$change_kind == "relabel"),
  sum(migration$change_kind == "hold")
))
