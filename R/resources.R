# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Versioned resource contracts for dictionaries and immutable KEGG snapshots.

lisa_dictionary_registry_columns <- function() {
  c("logical_id", "version", "species", "modality", "schema", "sha256",
    "compatibility", "approved_origin", "artifact")
}

lisa_resource_schemas <- function() {
  c(lisa_dictionary_schemas(), "term2gene@1", "category_map@1")
}

# Dictionary schema identities, oldest first.
#
# `lisa_dictionary@1` is the historical eight-column contract whose seventh
# column, `LISA_score`, was a redundant re-encoding of `tier` inherited from the
# upstream construction project (core <-> 4, expanded <-> 3, verified over both
# reviewed resources). `lisa_dictionary@2` is the canonical runtime contract and
# omits it.
#
# `@1` stays listed and fully supported on the read path: the meaning of a
# published schema identity is never rewritten in place. Resources declaring it
# are validated under their own historical rules and then normalised to the
# score-free runtime projection, so both schemas yield identical analysis.
lisa_dictionary_schemas <- function() {
  c("lisa_dictionary@1", "lisa_dictionary@2")
}

lisa_dictionary_runtime_columns <- function() {
  c("universe", "gene_set_id", "gene_set_name", "source_id",
    "category_id", "category_display_name", "tier")
}

# The single point where a dictionary becomes runtime data. Everything after
# this call sees the canonical score-free contract regardless of which schema
# the resource on disk declared.
#
# The historical column is read, validated by the caller and then discarded; it
# never re-acquires analytical meaning. Dropping it cannot merge rows: the
# reviewed resources have no two rows differing only in `LISA_score` (it is a
# function of `tier`), and `(universe, gene_set_id, category_id)` uniqueness is
# enforced independently by lisa_validate_dictionary_shape().
lisa_dictionary_runtime_projection <- function(x) {
  if ("LISA_score" %in% names(x)) x <- x[, setdiff(names(x), "LISA_score"), drop = FALSE]
  x
}

lisa_dictionary_has_legacy_score <- function(x) "LISA_score" %in% names(x)

lisa_resource_modalities <- function() {
  c("all", "transcriptomic/genomic", "global proteomic", "targeted", "custom")
}

lisa_empty_resource_registry <- function(internal = FALSE) {
  columns <- lisa_dictionary_registry_columns()
  out <- as.data.frame(
    stats::setNames(rep(list(character()), length(columns)), columns),
    stringsAsFactors = FALSE
  )
  if (isTRUE(internal)) {
    out$.storage <- character()
    out$.path <- character()
  }
  out
}

# The active scientific dictionaries. Their contents are the first published
# LISA dictionaries, so they carry version 1.0.0 under their own
# first-publication logical namespace (`lisa_dictionary_*`).
#
# The older `lisa_core`/`lisa_expanded` namespace is NOT reused for them.
# `lisa_core@1.0.0` and `lisa_expanded@1.0.0` are immutable published pins to
# different (score-carrying) bytes, and a published identity is never re-minted
# for other content. The old namespace stays readable through the closed
# current-release aliases below.
lisa_dictionary_tier_resources <- function() {
  c(core = "lisa_dictionary_core@1.0.0",
    expanded = "lisa_dictionary_expanded@1.0.0")
}

# Logical IDs that denote a *standard* LISA dictionary tier, as opposed to a
# custom or synthetic one. Both the current first-publication namespace and the
# retired scientific namespace are listed: a resource that was standard when it
# was published does not become custom because the active identity moved.
#
# Kept in one place because two independent decisions read it: which validator
# install_lisa_resource() applies, and which tier a configured
# `dictionary_resource` implies.
lisa_standard_dictionary_tiers <- function() {
  c(
    lisa_dictionary_core = "core",
    lisa_dictionary_expanded = "expanded",
    lisa_core = "core",
    lisa_expanded = "expanded"
  )
}

lisa_dictionary_logical_id_for_tier <- function(tier = c("core", "expanded")) {
  tier <- match.arg(tier)
  resource_id <- lisa_dictionary_tier_resources()[[tier]]
  sub("@.*$", "", resource_id)
}

# Closed identity authority for the score-carrying public dictionaries.
#
# These rows are intentionally code-authenticated rather than reconstructed
# from RESOURCE_MIGRATION_MANIFEST.tsv. That manifest is an audit ledger for the
# C5/C6 transition, not a trust root for runtime resolution: changing a
# transform row must never be able to redirect an old public pin to other bytes,
# another schema or another path. Historical rows also stay outside the active
# registry so unversioned/default selection remains on the score-free 2.0.0
# resources.
lisa_historical_dictionary_definitions <- function() {
  origin <- paste0(
    "LISA project candidate; CSIC-owned classification; profile-independent ",
    "gene-set/category relationships; RESOURCE_PROVENANCE.tsv; immutable ",
    "historical C5 public dictionary"
  )
  data.frame(
    logical_id = c("lisa_core", "lisa_expanded"),
    version = c("1.0.0", "1.0.0"),
    species = c("Homo sapiens", "Homo sapiens"),
    modality = c("all", "all"),
    schema = c("lisa_dictionary@1", "lisa_dictionary@1"),
    sha256 = c(
      "32d5ccafaec18fd6bdc51d3a23d80c21b7e996f5fe061a510c1b24ef3fafc876",
      "f8dd8bdbf61278d080c6349cebe9fd74993dc40bb07051164412e783358edaf4"
    ),
    compatibility = c("lisaR>=0.6.0", "lisaR>=0.6.0"),
    approved_origin = c(origin, origin),
    artifact = c(
      "extdata/dictionaries/lisa_core_runtime_v0_1.tsv",
      "extdata/dictionaries/lisa_expanded_runtime_v0_1.tsv"
    ),
    stringsAsFactors = FALSE
  )
}

# C5's 0.1.0 -> 1.0.0 compatibility aliases are also closed here. In
# particular, neither the relabel nor the later transform row in the migration
# ledger supplies runtime path/hash/schema identity.
lisa_historical_dictionary_aliases <- function() {
  definitions <- lisa_historical_dictionary_definitions()
  data.frame(
    logical_id = definitions$logical_id,
    old_resource_id = paste0(definitions$logical_id, "@0.1.0"),
    new_resource_id = paste0(definitions$logical_id, "@", definitions$version),
    old_version = rep("0.1.0", nrow(definitions)),
    new_version = definitions$version,
    species = definitions$species,
    modality = definitions$modality,
    schema = definitions$schema,
    old_sha256 = definitions$sha256,
    new_sha256 = definitions$sha256,
    stringsAsFactors = FALSE
  )
}

lisa_historical_dictionary_alias_match <- function(
  logical_id, version, species = NULL, modality = NULL
) {
  if (is.null(version)) return(NULL)
  aliases <- lisa_historical_dictionary_aliases()
  hit <- aliases$logical_id == logical_id & aliases$old_version == version
  if (!is.null(species)) hit <- hit & aliases$species == species
  if (!is.null(modality)) {
    hit <- hit & (aliases$modality == as.character(modality)[[1L]] |
                    aliases$modality == "all")
  }
  if (sum(hit) != 1L) return(NULL)
  aliases[which(hit)[[1L]], , drop = FALSE]
}

# Closed cross-namespace compatibility table for the first-publication rename.
#
# Each row says: a reference written against the identity a previous lisaR
# release resolved reads, today, exactly the same bytes under the current
# first-publication identity. `old_version = NA` covers the version-less form,
# which previously resolved through the registry and must keep resolving
# deterministically.
#
# Like lisa_historical_dictionary_definitions(), this is code-authenticated on
# purpose. RESOURCE_MIGRATION_MANIFEST.tsv documents the same mapping as an
# audit ledger, but it is not a trust root: editing a ledger row must never be
# able to redirect a published reference to different bytes. The declared
# digest and schema are re-checked against the resolved artifact before the
# resolver returns, so a future content edit cannot silently inherit one of
# these names.
#
# Deliberately absent: `lisa_core@1.0.0`, `lisa_core@0.1.0`,
# `lisa_expanded@1.0.0` and `lisa_expanded@0.1.0`. Those are historical pins to
# the score-carrying bytes and keep their existing authenticated path, hash and
# schema. Only identities that already pointed at the *current* contents are
# listed here, so no alias can ever redirect a historical pin.
lisa_current_release_dictionary_aliases <- function() {
  data.frame(
    old_logical_id = c(
      "lisa_core", "lisa_core",
      "lisa_expanded", "lisa_expanded",
      "lisa_quickstart_dictionary",
      "lisa_example_custom_dictionary"
    ),
    old_version = c(
      "2.0.0", NA_character_,
      "2.0.0", NA_character_,
      "4.0.0",
      "2.0.0"
    ),
    new_logical_id = c(
      "lisa_dictionary_core", "lisa_dictionary_core",
      "lisa_dictionary_expanded", "lisa_dictionary_expanded",
      "lisa_dictionary_quickstart",
      "lisa_dictionary_custom_example"
    ),
    new_version = rep("1.0.0", 6L),
    species = rep("Homo sapiens", 6L),
    modality = c(
      "all", "all",
      "all", "all",
      "transcriptomic/genomic",
      "transcriptomic/genomic"
    ),
    schema = rep("lisa_dictionary@2", 6L),
    sha256 = c(
      "f24b5bd8d9ac66b6d64c1a91c56ab567aa64a37f05dee5ce0ed3011b2bed8994",
      "f24b5bd8d9ac66b6d64c1a91c56ab567aa64a37f05dee5ce0ed3011b2bed8994",
      "0914c2d2e40369041a2313e913898f9c8e4065779e3973b076b70ea473045736",
      "0914c2d2e40369041a2313e913898f9c8e4065779e3973b076b70ea473045736",
      "72f36fe1a82cece837435e09939450e7bf8b90e57d179c218fd907c909235d0d",
      "2f1b4161824e700f8252374dd3f878088ae89f54d7f384d16a712000e1de7026"
    ),
    stringsAsFactors = FALSE
  )
}

# Match one requested identity against the closed table above. Returns NULL
# whenever the request is not a superseded first-publication reference, so
# callers fall through to ordinary resolution unchanged.
#
# `version = NULL` matches only the rows that declare the version-less form;
# an explicit version must match exactly. Modality is matched the way the
# registry matches it: a profile-independent "all" row answers any profile.
lisa_current_release_alias_match <- function(logical_id, version = NULL,
                                             species = NULL, modality = NULL) {
  aliases <- lisa_current_release_dictionary_aliases()
  hit <- aliases$old_logical_id == logical_id
  hit <- hit & if (is.null(version)) {
    is.na(aliases$old_version)
  } else {
    !is.na(aliases$old_version) & aliases$old_version == version
  }
  if (!is.null(species)) hit <- hit & aliases$species == species
  if (!is.null(modality)) {
    hit <- hit & (aliases$modality == as.character(modality)[[1L]] |
                    aliases$modality == "all")
  }
  if (sum(hit) != 1L) return(NULL)
  aliases[which(hit)[[1L]], , drop = FALSE]
}

lisa_pipeline_resource_defaults <- function(tier = c("core", "expanded")) {
  tier <- match.arg(tier)
  c(
    # The LISA project's own resources carry their first-publication version.
    # `msigdb_term2gene@2026.1` is deliberately untouched: it names the upstream
    # MSigDB release, which is provenance, not a version of our resource.
    dictionary_resource = unname(lisa_dictionary_tier_resources()[[tier]]),
    term2gene_resource = "msigdb_term2gene@2026.1",
    category_map_resource = "lisa_category_map@1.0.0"
  )
}

lisa_resource_migration_manifest_path <- function() {
  path <- system.file(
    "extdata", "dictionaries", "RESOURCE_MIGRATION_MANIFEST.tsv",
    package = "lisaR"
  )
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) ||
      !isTRUE(utils::file_test("-f", path)) || lisa_path_is_link(path)) {
    stop(
      "LISA-RESOURCE-029 installed resource migration manifest is absent or unsafe. Repair: reinstall the exact lisaR package candidate.",
      call. = FALSE
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

lisa_resource_migration_columns <- function() {
  c("change_kind", "logical_id", "old_resource_id", "new_resource_id",
    "old_version", "new_version", "species", "modality", "schema",
    "old_sha256", "new_sha256", "old_artifact", "new_artifact",
    "artifact_renamed", "content_transformed", "old_state", "new_state",
    "approved_origin", "rationale")
}

# Identities that a previous lisaR release resolved and that this release has
# renamed to a first-publication version. They are intentionally NOT registry
# rows: keeping them out of the registry preserves one ordinary active variant
# per resource, so a version-less reference such as "lisa_core" still resolves
# deterministically. They are instead a closed compatibility ledger, read from
# the same machine-readable migration manifest that documents the migration.
lisa_superseded_resource_identities <- function() {
  migration <- lisa_read_resource_migration_manifest()
  superseded <- migration[
    migration$change_kind == "relabel" &
      migration$old_state == "superseded_compatibility", ,
    drop = FALSE
  ]
  row.names(superseded) <- NULL
  superseded
}

# Exact all-modality keys that remain reserved even though they are not active
# defaults. The core/expanded rows come from code-authenticated historical
# definitions; category-map compatibility remains the unchanged C5 relabel.
lisa_reserved_historical_identities <- function() {
  superseded <- lisa_superseded_resource_identities()
  from_relabel <- data.frame(
    logical_id = superseded$logical_id,
    version = superseded$old_version,
    species = superseded$species,
    modality = superseded$modality,
    resource_id = superseded$old_resource_id,
    replacement_id = superseded$new_resource_id,
    stringsAsFactors = FALSE
  )
  definitions <- lisa_historical_dictionary_definitions()
  aliases <- lisa_historical_dictionary_aliases()
  authenticated <- rbind(
    data.frame(
      logical_id = definitions$logical_id,
      version = definitions$version,
      species = definitions$species,
      modality = definitions$modality,
      resource_id = paste0(definitions$logical_id, "@", definitions$version),
      replacement_id = paste0(definitions$logical_id, "@", definitions$version),
      stringsAsFactors = FALSE
    ),
    data.frame(
      logical_id = aliases$logical_id,
      version = aliases$old_version,
      species = aliases$species,
      modality = aliases$modality,
      resource_id = aliases$old_resource_id,
      replacement_id = aliases$new_resource_id,
      stringsAsFactors = FALSE
    )
  )
  # Identities that were ordinary built-in registry rows until the
  # first-publication rename. They lost the LISA-RESOURCE-020 built-in
  # collision guard when they left the bundle manifest, so reserve them
  # explicitly: an external all-modality row must not be able to re-mint
  # `lisa_core@2.0.0` with different bytes now that it is served by an alias.
  #
  # Only the version-carrying rows are reserved. The version-less alias rows
  # are not identities at all, and the synthetic custom-example fixture was
  # never a built-in registry row, so neither is reserved here.
  release_aliases <- lisa_current_release_dictionary_aliases()
  release_aliases <- release_aliases[
    !is.na(release_aliases$old_version) &
      release_aliases$old_logical_id != "lisa_example_custom_dictionary", ,
    drop = FALSE
  ]
  superseded_current <- data.frame(
    logical_id = release_aliases$old_logical_id,
    version = release_aliases$old_version,
    species = release_aliases$species,
    modality = release_aliases$modality,
    resource_id = lisa_resource_id(
      release_aliases$old_logical_id, release_aliases$old_version
    ),
    replacement_id = lisa_resource_id(
      release_aliases$new_logical_id, release_aliases$new_version
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(from_relabel, authenticated, superseded_current)
  key <- paste(out$logical_id, out$version, out$species, out$modality, sep = "\r")
  out[!duplicated(key), , drop = FALSE]
}

lisa_read_resource_migration_manifest <- function(
  path = lisa_resource_migration_manifest_path()
) {
  migration <- tryCatch(
    utils::read.delim(
      path, sep = "\t", header = TRUE, quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
      na.strings = character()
    ),
    error = function(error) stop(
      "LISA-RESOURCE-029 resource migration manifest could not be read: ",
      conditionMessage(error),
      ". Repair: reinstall the exact lisaR package candidate.",
      call. = FALSE
    )
  )
  required <- lisa_resource_migration_columns()
  if (!identical(names(migration), required)) {
    stop(
      "LISA-RESOURCE-029 resource migration manifest does not use the documented schema. Repair: rebuild it with tools/maintenance/build_resource_migration.R.",
      call. = FALSE
    )
  }
  if (!nrow(migration)) return(migration)
  # A relabel that claims transformed content or a renamed artifact would make
  # the compatibility guarantee below unverifiable, so refuse it outright.
  #
  # Subset first, then test. The previous form combined a full-length mask with
  # already-subset vectors, so R recycled the shorter operand: the guard did not
  # actually test the relabel rows it named, and emitted a recycling warning as
  # soon as the row count stopped being a multiple of the relabel count. C6
  # exposed that by adding non-relabel rows.
  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]
  if (nrow(relabel) && any(relabel$content_transformed != "FALSE" |
                           relabel$artifact_renamed != "FALSE" |
                           relabel$old_sha256 != relabel$new_sha256)) {
    stop(
      "LISA-RESOURCE-029 resource migration manifest declares a relabel that changes bytes. Repair: record a transformation with a new logical identity instead.",
      call. = FALSE
    )
  }
  if (!lisa_sha256_all_valid(migration$old_sha256, nrow(migration)) ||
      !lisa_sha256_all_valid(migration$new_sha256, nrow(migration))) {
    stop(
      "LISA-RESOURCE-029 resource migration manifest has a malformed SHA-256. Repair: rebuild it with tools/maintenance/build_resource_migration.R.",
      call. = FALSE
    )
  }
  migration
}

# Resolve a requested identity through the compatibility ledger. Returns NULL
# when the request is not a superseded identity, so callers fall through to
# ordinary registry resolution unchanged.
lisa_superseded_identity_match <- function(logical_id, version, species = NULL,
                                           modality = NULL) {
  if (is.null(version)) return(NULL)
  superseded <- lisa_superseded_resource_identities()
  if (!nrow(superseded)) return(NULL)
  hit <- superseded$logical_id == logical_id &
    superseded$old_version == version
  if (!is.null(species)) hit <- hit & superseded$species == species
  # Modality is matched permissively: a superseded row registered for the
  # profile-independent modality "all" answers any requested profile, exactly
  # as the registry itself does.
  if (!is.null(modality)) {
    hit <- hit & (superseded$modality == as.character(modality)[[1L]] |
                    superseded$modality == "all")
  }
  if (sum(hit) != 1L) return(NULL)
  superseded[which(hit)[[1L]], , drop = FALSE]
}

lisa_builtin_resource_manifest_path <- function() {
  path <- system.file(
    "extdata", "dictionaries", "RESOURCE_BUNDLE_MANIFEST.tsv",
    package = "lisaR"
  )
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) ||
      !isTRUE(utils::file_test("-f", path)) || lisa_path_is_link(path)) {
    stop(
      "LISA-RESOURCE-009 installed package resource manifest is absent or unsafe. Repair: reinstall the exact lisaR package candidate.",
      call. = FALSE
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

lisa_builtin_resource_definitions <- function() {
  lisa_read_dictionary_registry(lisa_builtin_resource_manifest_path())
}

lisa_builtin_resource_registry <- function() {
  entries <- lisa_builtin_resource_definitions()
  paths <- vapply(entries$artifact, function(artifact) {
    system.file(artifact, package = "lisaR")
  }, character(1))
  regular <- nzchar(paths) & file.exists(paths) & !dir.exists(paths) &
    vapply(paths, function(path) isTRUE(utils::file_test("-f", path)), logical(1)) &
    !vapply(paths, lisa_path_is_link, logical(1))
  if (!all(regular)) {
    stop(sprintf(
      "LISA-RESOURCE-009 installed package resource artifact is absent or unsafe: %s. Repair: reinstall the exact lisaR package candidate.",
      entries$artifact[[which(!regular)[[1L]]]]
    ), call. = FALSE)
  }
  paths <- vapply(paths, normalizePath, character(1), winslash = "/", mustWork = TRUE)
  observed <- vapply(paths, lisa_sha256_file, character(1))
  mismatch <- which(observed != entries$sha256)
  if (length(mismatch)) {
    stop(sprintf(
      "LISA-RESOURCE-007 installed package resource '%s' SHA-256 mismatch. Repair: reinstall the exact lisaR package candidate.",
      lisa_resource_id(entries$logical_id[[mismatch[[1L]]]], entries$version[[mismatch[[1L]]]])
    ), call. = FALSE)
  }
  entries <- entries[, lisa_dictionary_registry_columns(), drop = FALSE]
  entries$.storage <- "builtin"
  entries$.path <- paths
  row.names(entries) <- NULL
  entries
}

# Resolve and authenticate the closed historical rows without making them
# active/default registry entries. This deliberately mirrors the regular
# built-in integrity checks: exact code-owned metadata, regular non-symlink
# package paths, and SHA-256 verification before the resolver may return a row.
lisa_historical_dictionary_registry <- function() {
  entries <- lisa_read_dictionary_registry(
    lisa_historical_dictionary_definitions()
  )
  paths <- vapply(entries$artifact, function(artifact) {
    system.file(artifact, package = "lisaR")
  }, character(1))
  regular <- nzchar(paths) & file.exists(paths) & !dir.exists(paths) &
    vapply(paths, function(path) isTRUE(utils::file_test("-f", path)), logical(1)) &
    !vapply(paths, lisa_path_is_link, logical(1))
  if (!all(regular)) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-009 installed historical resource artifact is absent or unsafe: %s. ",
        "Repair: reinstall the exact lisaR package candidate."
      ),
      entries$artifact[[which(!regular)[[1L]]]]
    ), call. = FALSE)
  }
  paths <- vapply(
    paths, normalizePath, character(1), winslash = "/", mustWork = TRUE
  )
  observed <- vapply(paths, lisa_sha256_file, character(1))
  mismatch <- which(observed != entries$sha256)
  if (length(mismatch)) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-007 installed historical resource '%s' SHA-256 mismatch. ",
        "Repair: reinstall the exact lisaR package candidate."
      ),
      lisa_resource_id(
        entries$logical_id[[mismatch[[1L]]]],
        entries$version[[mismatch[[1L]]]]
      )
    ), call. = FALSE)
  }
  entries <- entries[, lisa_dictionary_registry_columns(), drop = FALSE]
  entries$.storage <- "historical_builtin"
  entries$.path <- paths
  row.names(entries) <- NULL
  entries
}

# Backward-compatible internal name. This now describes only synthetic files
# that are genuinely installed with lisaR; scientific defaults are registered
# external resources and never placeholder built-ins.
lisa_builtin_dictionary_registry <- function() {
  registry <- lisa_builtin_resource_registry()
  registry[, lisa_dictionary_registry_columns(), drop = FALSE]
}

lisa_dictionary_cache_root <- function() {
  root <- getOption("lisaR.dictionary_cache_root", Sys.getenv("LISAR_DICTIONARY_CACHE", ""))
  if (!nzchar(root)) root <- file.path(path.expand("~"), ".local", "share", "lisaR", "dictionaries")
  # Every other managed path in this package is canonical and forward-slashed
  # (see lisa_resource_cache_root below). Plain normalizePath() returned
  # backslashes on Windows and left a not-yet-created cache root in the caller
  # spelling, so paths derived from here could not be compared with the
  # canonical registry path that install_lisa_resource reports.
  lisa_path_canonical(root)
}

lisa_builtin_registered_resource <- function(logical_id, version = NULL,
                                             species = NULL) {
  entries <- lisa_builtin_resource_registry()
  selected <- entries[entries$logical_id == logical_id, , drop = FALSE]
  if (!is.null(version)) selected <- selected[selected$version == version, , drop = FALSE]
  if (!is.null(species)) selected <- selected[selected$species == species, , drop = FALSE]
  if (nrow(selected) != 1L) return(NULL)
  path <- selected$.path[[1L]]
  c(
    as.list(selected[1, lisa_dictionary_registry_columns(), drop = FALSE]),
    list(
      path = path,
      resource_id = lisa_resource_id(
        selected$logical_id[[1]], selected$version[[1]]
      )
    )
  )
}

lisa_builtin_dictionary_path <- function(tier = c("core", "expanded")) {
  tier <- match.arg(tier)
  # Read the active version from the tier map rather than hard-coding it, so a
  # future resource migration cannot leave this helper pointing at a retired
  # identity, as C6 would otherwise have done.
  resource_id <- lisa_dictionary_tier_resources()[[tier]]
  # Derive both halves from the tier map. Deriving the logical ID as
  # paste0("lisa_", tier) silently assumed the retired namespace and would have
  # pointed this helper at an identity the registry no longer serves.
  resource <- lisa_builtin_registered_resource(
    sub("@.*$", "", resource_id),
    sub("^.*@", "", resource_id), "Homo sapiens"
  )
  if (is.null(resource)) {
    stop(sprintf(
      "LISA-RESOURCE-013 installed scientific %s dictionary is unavailable. Repair: reinstall the exact lisaR package candidate.",
      tier
    ), call. = FALSE)
  }
  resource$path
}

lisa_builtin_category_map_path <- function() {
  resource <- lisa_builtin_registered_resource(
    "lisa_category_map", "1.0.0", "Homo sapiens"
  )
  if (is.null(resource)) {
    stop(
      "LISA-RESOURCE-014 installed scientific category map is unavailable. Repair: reinstall the exact lisaR package candidate.",
      call. = FALSE
    )
  }
  resource$path
}

lisa_registry_artifact_is_safe <- function(artifact) {
  if (length(artifact) != 1L || is.na(artifact) ||
      !nzchar(trimws(as.character(artifact)))) return(FALSE)
  artifact <- as.character(artifact)
  if (!identical(artifact, trimws(artifact)) ||
      grepl("[[:cntrl:]\\\\]", artifact) || startsWith(artifact, "/") ||
      grepl("^[A-Za-z]:", artifact)) return(FALSE)
  parts <- strsplit(artifact, "/", fixed = TRUE)[[1L]]
  length(parts) > 0L && all(nzchar(parts)) && !any(parts %in% c(".", ".."))
}

lisa_read_dictionary_registry <- function(registry = NULL) {
  if (is.null(registry)) return(lisa_empty_resource_registry())
  if (is.character(registry) && length(registry) == 1L) {
    if (is.na(registry) || !nzchar(trimws(registry)) ||
        !file.exists(registry) || dir.exists(registry) ||
        !isTRUE(utils::file_test("-f", registry)) ||
        lisa_path_is_link(registry)) {
      stop(
        "LISA-RESOURCE-001 configured registry must be one existing regular non-symlink resource_registry.tsv. Repair: configure the reviewed registry file directly.",
        call. = FALSE
      )
    }
    registry <- tryCatch(
      utils::read.delim(
        registry, sep = "\t", header = TRUE, quote = "",
        comment.char = "", check.names = FALSE,
        stringsAsFactors = FALSE, colClasses = "character",
        na.strings = character()
      ),
      error = function(error) stop(
        "LISA-RESOURCE-001 configured registry TSV could not be read: ",
        conditionMessage(error),
        ". Repair: provide a valid tab-separated resource_registry.tsv.",
        call. = FALSE
      )
    )
  }
  registry <- as.data.frame(registry, stringsAsFactors = FALSE)
  required <- lisa_dictionary_registry_columns()
  missing <- setdiff(required, names(registry))
  if (length(missing)) stop(sprintf("LISA-RESOURCE-001 registry is missing required column(s): %s. Repair: provide %s.", paste(missing, collapse = ", "), paste(required, collapse = ", ")), call. = FALSE)
  extra <- setdiff(names(registry), required)
  if (length(extra)) stop(sprintf(
    paste0(
      "LISA-RESOURCE-001 registry contains unsupported column(s): %s. ",
      "Repair: use exactly the documented registry schema; keep administrative metadata in a separate ledger."
    ),
    paste(extra, collapse = ", ")
  ), call. = FALSE)
  registry <- registry[, required, drop = FALSE]
  if (!nrow(registry)) return(registry)
  registry[] <- lapply(registry, as.character)
  nonempty <- required
  has_empty <- vapply(registry[nonempty], function(value) {
    any(is.na(value) | !nzchar(trimws(value)))
  }, logical(1))
  if (any(has_empty)) {
    stop("LISA-RESOURCE-002 registry contains an empty required value. Repair: complete every identity, integrity and compatibility field.", call. = FALSE)
  }
  noncanonical <- vapply(registry, function(value) {
    any(value != trimws(value) | grepl("[[:cntrl:]]", value))
  }, logical(1))
  if (any(noncanonical)) {
    stop(
      "LISA-RESOURCE-002 registry contains whitespace-delimited or control-character values. Repair: use canonical one-line TSV fields without leading or trailing whitespace.",
      call. = FALSE
    )
  }
  safe_ids <- tryCatch({
    vapply(registry$logical_id, lisa_safe_id, character(1), field = "registry logical_id")
    vapply(registry$version, lisa_safe_id, character(1), field = "registry version")
    TRUE
  }, error = function(error) FALSE)
  if (!isTRUE(safe_ids) ||
      !all(vapply(registry$artifact, lisa_registry_artifact_is_safe, logical(1)))) {
    stop("LISA-RESOURCE-002 registry contains an unsafe logical_id, version or artifact path. Repair: use portable IDs and a relative artifact path without traversal.", call. = FALSE)
  }
  valid_species <- vapply(registry$species, function(value) {
    tryCatch(
      identical(lisa_species_contract(value)$scientific_name, value),
      error = function(error) FALSE
    )
  }, logical(1))
  if (!all(valid_species) ||
      !all(registry$modality %in% lisa_resource_modalities()) ||
      !all(registry$schema %in% lisa_resource_schemas()) ||
      !lisa_sha256_all_valid(registry$sha256, nrow(registry))) {
    stop(
      "LISA-RESOURCE-002 registry contains an unsupported species, modality, schema or SHA-256 value. Repair: use the documented canonical registry contract.",
      call. = FALSE
    )
  }
  valid_compatibility <- grepl(
    "^lisaR\\s*(>=|<=|==|>|<)\\s*[0-9]+(?:\\.[0-9]+){1,3}$",
    registry$compatibility, perl = TRUE
  )
  if (!all(valid_compatibility)) {
    stop(
      "LISA-RESOURCE-022 registry compatibility is malformed. Repair: use one closed constraint such as lisaR>=0.6.0.",
      call. = FALSE
    )
  }
  registry_key <- paste(
    registry$logical_id, registry$version, registry$species,
    registry$modality, sep = "@"
  )
  if (anyDuplicated(registry_key)) {
    stop("LISA-RESOURCE-003 registry has duplicate logical ID/version/species/modality rows. Repair: retain one approved row per compatible resource.", call. = FALSE)
  }
  registry
}

lisa_resource_id <- function(logical_id, version) paste0(logical_id, "@", version)

lisa_active_resource_registry <- function(registry = NULL) {
  builtin <- lisa_builtin_resource_registry()
  external <- lisa_read_dictionary_registry(registry)
  if (nrow(external) && nrow(builtin)) {
    builtin_key <- paste(
      builtin$logical_id, builtin$version, builtin$species, builtin$modality,
      sep = "\r"
    )
    external_key <- paste(
      external$logical_id, external$version, external$species,
      external$modality, sep = "\r"
    )
    collision <- match(external_key, builtin_key, nomatch = 0L)
    if (any(collision > 0L)) {
      row <- which(collision > 0L)[[1L]]
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-020 external registry redefines built-in resource '%s' for species '%s' and modality '%s'. ",
          "Repair: remove that row; built-in fixture identities cannot be hidden or replaced."
        ),
        lisa_resource_id(external$logical_id[[row]], external$version[[row]]),
        external$species[[row]], external$modality[[row]]
      ), call. = FALSE)
    }
  }
  if (nrow(external)) {
    # Published built-in identities stay reserved even when they are served by
    # the closed historical path rather than the active/default registry.
    # Otherwise an external all-modality row could re-mint an old public pin
    # with different bytes. Profile-scoped rows remain supported exactly as in
    # C5 and are resolved/checked below.
    reserved <- lisa_reserved_historical_identities()
    if (nrow(reserved)) {
      reserved_key <- paste(
        reserved$logical_id, reserved$version, reserved$species,
        reserved$modality, sep = "\r"
      )
      external_key <- paste(
        external$logical_id, external$version, external$species,
        external$modality, sep = "\r"
      )
      conflict <- which(external_key %in% reserved_key)
      if (length(conflict)) {
        row <- conflict[[1L]]
        historical <- reserved[
          match(external_key[[row]], reserved_key), , drop = FALSE
        ]
        stop(sprintf(
          paste0(
            "LISA-RESOURCE-030 external registry redefines historical built-in identity '%s'. ",
            "Repair: remove that all-modality row and reference '%s'; a published resource identity cannot be re-minted with different content."
          ),
          historical$resource_id[[1L]], historical$replacement_id[[1L]]
        ), call. = FALSE)
      }
    }
  }
  if (nrow(external)) {
    external$.storage <- "cache"
    external$.path <- ""
  } else {
    external <- lisa_empty_resource_registry(internal = TRUE)
  }
  out <- rbind(builtin, external)
  row.names(out) <- NULL
  out
}

lisa_validate_resource_compatibility <- function(
  compatibility,
  package_version = as.character(utils::packageVersion("lisaR"))
) {
  if (!is.character(compatibility) || length(compatibility) != 1L ||
      is.na(compatibility)) {
    stop(
      "LISA-RESOURCE-022 resource compatibility is malformed. Repair: use a constraint such as lisaR>=0.6.0.",
      call. = FALSE
    )
  }
  match <- regexec(
    "^lisaR\\s*(>=|<=|==|>|<)\\s*([0-9]+(?:\\.[0-9]+){1,3})$",
    trimws(compatibility), perl = TRUE
  )
  parts <- regmatches(trimws(compatibility), match)[[1L]]
  if (length(parts) != 3L) {
    stop(sprintf(
      "LISA-RESOURCE-022 compatibility '%s' is malformed. Repair: use one closed constraint such as lisaR>=0.6.0.",
      compatibility
    ), call. = FALSE)
  }
  comparison <- tryCatch(
    utils::compareVersion(package_version, parts[[3L]]),
    error = function(error) NA_integer_
  )
  satisfied <- switch(parts[[2L]],
    ">=" = comparison >= 0L,
    "<=" = comparison <= 0L,
    "==" = comparison == 0L,
    ">" = comparison > 0L,
    "<" = comparison < 0L,
    FALSE
  )
  if (is.na(comparison) || !isTRUE(satisfied)) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-022 compatibility constraint '%s' is not satisfied by active lisaR '%s'. ",
        "Repair: install a compatible resource version or use the required lisaR version."
      ), compatibility, package_version
    ), call. = FALSE)
  }
  invisible(TRUE)
}

lisa_validate_hpc_dictionary_root <- function(root) {
  root <- lisa_path_walk_lexically(root)$path
  if (lisa_path_is_link(root)) {
    stop("LISA-RESOURCE-005 configured shared dictionary root cannot be a symbolic link. Repair: configure the reviewed directory directly.", call. = FALSE)
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!dir.exists(root) || file.access(root, 4L) != 0L) stop("LISA-RESOURCE-004 configured shared dictionary root is not readable. Repair: configure a readable administrator-managed root.", call. = FALSE)
  if (file.access(root, 2L) == 0L) stop("LISA-RESOURCE-005 configured shared dictionary root must be read-only. Repair: mount or configure the shared registry root read-only.", call. = FALSE)
  root
}

lisa_resource_cache_root <- function(root) {
  root <- lisa_path_walk_lexically(root)$path
  if (lisa_path_is_link(root)) {
    stop(
      "LISA-RESOURCE-008 configured resource cache root cannot be a symbolic link. Repair: configure the managed cache directory directly.",
      call. = FALSE
    )
  }
  if (dir.exists(root)) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  }
  root
}

lisa_resource_artifact_path <- function(root, relative) {
  root <- lisa_resource_cache_root(root)
  candidate <- lisa_path_walk_lexically(file.path(root, relative))$path
  if (!lisa_path_within(candidate, root)) {
    stop(
      "LISA-RESOURCE-002 registered artifact escapes its managed resource root. Repair: use a relative artifact path without traversal.",
      call. = FALSE
    )
  }
  walked <- lisa_path_walk_lexically(candidate)
  managed_components <- walked$visited[vapply(
    walked$visited, lisa_path_within, logical(1), root = root
  )]
  linked <- managed_components[vapply(
    managed_components, lisa_path_is_link, logical(1)
  )]
  if (length(linked)) {
    stop(sprintf(
      "LISA-RESOURCE-008 registered artifact path contains a symbolic link inside the managed root: %s. Repair: install regular directories and files directly in the cache.",
      linked[[1L]]
    ), call. = FALSE)
  }
  if (file.exists(candidate)) {
    canonical <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
    if (!lisa_path_within(canonical, root)) {
      stop(
        "LISA-RESOURCE-008 registered artifact resolves outside its managed resource root. Repair: remove linked cache components and reinstall the resource.",
        call. = FALSE
      )
    }
    candidate <- canonical
  }
  candidate
}

lisa_dictionary_resource_path <- function(entry, cache_root = lisa_dictionary_cache_root(), shared_root = NULL) {
  root <- if (!is.null(shared_root) && nzchar(shared_root)) {
    lisa_validate_hpc_dictionary_root(shared_root)
  } else {
    lisa_resource_cache_root(cache_root)
  }
  if (!lisa_registry_artifact_is_safe(entry$artifact[[1L]])) {
    stop("LISA-RESOURCE-002 registered artifact path is unsafe. Repair: use one relative cache path without traversal.", call. = FALSE)
  }
  lisa_resource_artifact_path(
    root,
    file.path(
      entry$logical_id[[1]], entry$version[[1]], entry$artifact[[1]]
    )
  )
}

# Select the registry rows eligible for one requested identity: logical ID,
# then version, species and finally profile, preferring an exact modality row
# before the profile-independent "all" row. Kept separate from the resolver so
# the identical selection can be attempted for a superseded pin and, when that
# leaves nothing eligible, for the identity the compatibility ledger renamed it
# to. `before_modality` is returned so the caller can still tell "no such
# version or species" apart from "registered, but not for this profile".
lisa_select_resource_rows <- function(registry, logical_id, version = NULL,
                                      species = NULL, modality = NULL) {
  selected <- registry[registry$logical_id == logical_id, , drop = FALSE]
  if (!is.null(version)) selected <- selected[selected$version == version, , drop = FALSE]
  if (!is.null(species)) selected <- selected[selected$species == species, , drop = FALSE]
  before_modality <- selected
  if (!is.null(modality)) {
    exact <- selected[selected$modality == modality, , drop = FALSE]
    selected <- if (nrow(exact)) {
      exact
    } else {
      selected[selected$modality == "all", , drop = FALSE]
    }
  } else if (nrow(selected) > 1L) {
    # A profile-independent row is the only deterministic choice when this
    # low-level resolver is called without a requested modality. Pipeline
    # resolution always supplies its profile and therefore prefers an exact
    # row before falling back to `all`.
    independent <- selected[selected$modality == "all", , drop = FALSE]
    if (nrow(independent) == 1L) selected <- independent
  }
  list(selected = selected, before_modality = before_modality)
}

lisa_resolve_dictionary_resource <- function(logical_id, version = NULL, species = NULL,
                                             modality = NULL, registry = NULL,
                                             cache_root = lisa_dictionary_cache_root(), shared_root = NULL) {
  logical_id <- lisa_safe_id(logical_id, "resource logical_id")
  if (!is.null(version)) version <- lisa_safe_id(version, "resource version")
  if (!is.null(species)) {
    species <- lisa_species_contract(species)$scientific_name
  }
  registry <- lisa_active_resource_registry(registry)
  if (!is.null(modality)) modality <- as.character(modality)[[1L]]
  requested_logical_id <- logical_id
  requested_resource_id <- if (is.null(version)) {
    NULL
  } else {
    lisa_resource_id(logical_id, version)
  }
  # Historical rows participate only for an exact supported version. They are
  # never appended for an unversioned/default request, so the active 2.0.0
  # registry remains the sole default-selection authority.
  resolution_registry <- registry
  if (!is.null(version)) {
    historical_definitions <- lisa_historical_dictionary_definitions()
    exact_historical <- historical_definitions$logical_id == logical_id &
      historical_definitions$version == version
    if (any(exact_historical)) {
      resolution_registry <- rbind(
        registry,
        lisa_historical_dictionary_registry()[exact_historical, , drop = FALSE]
      )
    }
  }
  selection <- lisa_select_resource_rows(
    resolution_registry, logical_id, version,
    species = species, modality = modality
  )
  # An old pin keeps working. The compatibility ledger is consulted when the
  # requested identity has no *eligible* row for the requested species and
  # profile (exact modality first, then "all") — not merely when the version is
  # absent from the registry entirely. An explicitly installed resource still
  # wins for the profile it was registered under, but a profile-scoped override
  # of a superseded identity no longer disables the historical all-modality
  # fallback under every other profile. A row that does match is never replaced
  # silently: it is selected here and, if it is ambiguous, corrupt or otherwise
  # invalid, it fails loudly below.
  superseded <- NULL
  authenticated_historical_alias <- NULL
  current_release_alias <- NULL
  if (!is.null(version) && !nrow(selection$selected)) {
    authenticated_historical_alias <- lisa_historical_dictionary_alias_match(
      logical_id, version, species = species, modality = modality
    )
    candidate <- authenticated_historical_alias
    if (is.null(candidate)) {
      candidate <- lisa_superseded_identity_match(
        logical_id, version, species = species, modality = modality
      )
    }
    if (!is.null(candidate)) {
      redirect_registry <- registry
      if (!is.null(authenticated_historical_alias)) {
        redirect_registry <- rbind(
          registry, lisa_historical_dictionary_registry()
        )
      }
      redirected <- lisa_select_resource_rows(
        redirect_registry, logical_id, candidate$new_version[[1L]],
        species = species, modality = modality
      )
      # Adopt the redirect only when it yields a row. Otherwise keep the
      # original outcome, whose diagnostics name the identity that was pinned.
      if (nrow(redirected$selected)) {
        superseded <- candidate
        version <- candidate$new_version[[1L]]
        selection <- redirected
      }
    }
  }
  # Last, and only if nothing above resolved: the cross-namespace
  # first-publication rename. Ordering matters. Every same-namespace route —
  # the active registry, the exact historical rows and the C5 relabel ledger —
  # has already been tried, so this alias can only ever answer a reference that
  # would otherwise have failed. It cannot shadow a historical pin.
  #
  # Unlike the branches above it also runs for a version-less reference, which
  # is exactly the case the rename would otherwise have broken.
  current_release_alias_requested_id <- requested_resource_id %||%
    requested_logical_id
  if (is.null(superseded) && !nrow(selection$selected)) {
    candidate <- lisa_current_release_alias_match(
      logical_id, version, species = species, modality = modality
    )
    if (!is.null(candidate)) {
      redirected <- lisa_select_resource_rows(
        registry, candidate$new_logical_id[[1L]], candidate$new_version[[1L]],
        species = species, modality = modality
      )
      if (nrow(redirected$selected)) {
        current_release_alias <- candidate
        logical_id <- candidate$new_logical_id[[1L]]
        version <- candidate$new_version[[1L]]
        selection <- redirected
      }
    }
  }
  selected <- selection$selected
  if (!is.null(modality) && !nrow(selected) && nrow(selection$before_modality)) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-021 resource '%s' is not registered for modality '%s'. ",
        "Repair: register a species- and modality-matched resource or select the matching profile."
      ),
      if (is.null(version)) logical_id else lisa_resource_id(logical_id, version),
      modality
    ), call. = FALSE)
  }
  if (nrow(selected) != 1L) stop(sprintf("LISA-RESOURCE-006 resource '%s' did not resolve to exactly one registered version for species '%s'. Repair: declare an approved logical ID, version, species and modality.", if (is.null(version)) logical_id else lisa_resource_id(logical_id, version), species %||% "unspecified"), call. = FALSE)
  lisa_validate_resource_compatibility(selected$compatibility[[1L]])
  if (!lisa_sha256_all_valid(selected$sha256, nrow(selected))) stop(sprintf("LISA-RESOURCE-007 resource '%s' has no approved lower-case SHA-256. Repair: install a reviewed registry entry before analysis.", lisa_resource_id(selected$logical_id[[1]], selected$version[[1]])), call. = FALSE)
  path <- if (selected$.storage[[1L]] %in% c("builtin", "historical_builtin")) {
    selected$.path[[1L]]
  } else {
    lisa_dictionary_resource_path(selected, cache_root, shared_root)
  }
  if (!file.exists(path) || dir.exists(path) ||
      !isTRUE(utils::file_test("-f", path)) || lisa_path_is_link(path)) stop(sprintf("LISA-RESOURCE-008 registered resource '%s' is absent or is not one regular non-symlink file at %s. Repair: install the approved resource into the canonical cache; do not supply a per-run path.", lisa_resource_id(selected$logical_id[[1]], selected$version[[1]]), path), call. = FALSE)
  observed <- lisa_sha256_file(path)
  if (!identical(observed, selected$sha256[[1]])) stop(sprintf("LISA-RESOURCE-009 SHA-256 mismatch for '%s'. Repair: replace the corrupt cache artifact with the approved version.", lisa_resource_id(selected$logical_id[[1]], selected$version[[1]])), call. = FALSE)
  if (!selected$.storage[[1L]] %in% c("builtin", "historical_builtin")) {
    checked_path <- lisa_dictionary_resource_path(
      selected, cache_root = cache_root, shared_root = shared_root
    )
    if (!identical(checked_path, path)) {
      stop(
        "LISA-RESOURCE-008 registered resource path changed during validation. Repair: secure the managed cache and retry.",
        call. = FALSE
      )
    }
  }
  if (!is.null(superseded)) {
    # The whole point of the compatibility ledger is that an old pin reads the
    # same science it always read. Prove it against the recorded digest rather
    # than trusting the rename, so a future content edit cannot silently
    # inherit a published identity.
    schema_mismatch <- !is.null(authenticated_historical_alias) &&
      !identical(selected$schema[[1L]], authenticated_historical_alias$schema[[1L]])
    if (!identical(observed, superseded$old_sha256[[1L]]) || schema_mismatch) {
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-030 superseded resource '%s' no longer resolves to its recorded bytes. ",
          "Repair: restore the migrated artifact or register a new logical identity; a published identity must not be reused for different content."
        ), superseded$old_resource_id[[1L]]
      ), call. = FALSE)
    }
  }
  if (!is.null(current_release_alias)) {
    # Same obligation as the compatibility ledger above: prove that the old
    # name still reads the science it always read. The alias declares the
    # digest and schema it promises, so an unrelated content edit under the new
    # identity cannot quietly acquire a published name.
    if (!identical(observed, current_release_alias$sha256[[1L]]) ||
        !identical(selected$schema[[1L]], current_release_alias$schema[[1L]])) {
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-030 superseded resource '%s' no longer resolves to its recorded bytes. ",
          "Repair: restore the migrated artifact or register a new logical identity; a published identity must not be reused for different content."
        ),
        current_release_alias_requested_id
      ), call. = FALSE)
    }
  }
  out <- c(as.list(selected[1, lisa_dictionary_registry_columns(), drop = FALSE]), list(path = normalizePath(path, winslash = "/", mustWork = TRUE), resource_id = lisa_resource_id(selected$logical_id[[1]], selected$version[[1]])))
  # Provenance keeps both identities: what the run asked for and what it read.
  out$requested_resource_id <- requested_resource_id %||% out$resource_id
  # The logical ID the caller actually wrote. A version-less reference through
  # the first-publication rename has no requested *resource* ID, so without
  # this the receipt would lose the fact that `lisa_core` was asked for.
  out$requested_logical_id <- requested_logical_id
  out$superseded_resource_id <- if (!is.null(superseded)) {
    superseded$old_resource_id[[1L]]
  } else if (!is.null(current_release_alias)) {
    current_release_alias_requested_id
  } else {
    NA_character_
  }
  out
}

lisa_pipeline_resource_spec <- function() {
  data.frame(
    config_key = c(
      "dictionary_resource", "term2gene_resource", "category_map_resource"
    ),
    # The canonical schema for each role, used when reporting an unresolved
    # resource. Acceptance is governed by lisa_role_accepted_schemas(), which is
    # wider for dictionaries because the legacy schema is still readable.
    expected_schema = c(
      "lisa_dictionary@2", "term2gene@1", "category_map@1"
    ),
    stringsAsFactors = FALSE
  )
}

# Schemas a resolved resource may legitimately declare for a given role.
# Dictionaries accept both contracts; the reader normalises `@1` to the
# score-free runtime projection. Other roles are unchanged by C6.
lisa_role_accepted_schemas <- function(config_key) {
  switch(
    config_key,
    dictionary_resource = lisa_dictionary_schemas(),
    term2gene_resource = "term2gene@1",
    category_map_resource = "category_map@1",
    stop(sprintf("unknown resource role '%s'", config_key), call. = FALSE)
  )
}

lisa_active_resource_context <- function() {
  cache_root <- lisa_dictionary_cache_root()
  registry <- getOption("lisaR.dictionary_registry", NULL)
  if (is.null(registry)) {
    default_registry <- file.path(cache_root, "resource_registry.tsv")
    if (file.exists(default_registry)) registry <- default_registry
  }
  shared_root <- getOption(
    "lisaR.shared_dictionary_root",
    Sys.getenv("LISAR_SHARED_DICTIONARY_ROOT", "")
  )
  if (!is.character(shared_root) || length(shared_root) != 1L ||
      is.na(shared_root) || !nzchar(shared_root)) {
    shared_root <- NULL
  }
  list(
    registry = registry,
    cache_root = cache_root,
    shared_root = shared_root
  )
}

lisa_resolve_low_level_resource_paths <- function(
  species, lisa_dictionary, universes,
  lisa_dictionary_path = NULL, term2gene_path = NULL,
  category_map_path = NULL,
  profile = "transcriptomic/genomic"
) {
  universes <- normalize_lisa_collections(universes)
  has_semantic_collection <- any(universes != "HALLMARKS")
  paths <- list(
    dictionary_resource = lisa_dictionary_path,
    term2gene_resource = term2gene_path,
    category_map_resource = category_map_path
  )

  if (has_semantic_collection) {
    supplied <- !vapply(paths, is.null, logical(1))
    if (all(!supplied)) {
      tier <- as.character(lisa_dictionary)[[1L]]
      if (!tier %in% c("core", "expanded")) {
        stop(
          paste0(
            "LISA-RESOURCE-010 the low-level API can resolve registered ",
            "defaults only for the core or expanded tier. Repair: provide ",
            "all three explicit resource paths for an external custom tier."
          ),
          call. = FALSE
        )
      }
      resolved <- lisa_resolve_pipeline_resources(
        pipeline = list(profile = profile, lisa_dictionary = tier),
        species = species,
        stop_on_error = TRUE
      )
      return(list(
        lisa_dictionary_path = resolved$resolved$dictionary_resource$path,
        term2gene_path = resolved$resolved$term2gene_resource$path,
        category_map_path = resolved$resolved$category_map_resource$path,
        registered_category_map = TRUE,
        source = "registered_defaults",
        resources = resolved$status
      ))
    }
    if (!all(supplied)) {
      missing <- names(paths)[!supplied]
      stop(
        paste0(
          "LISA-RESOURCE-027 incomplete low-level resource bundle; missing ",
          paste(missing, collapse = ", "),
          ". Repair: omit all three paths to use registered defaults, or ",
          "provide lisa_dictionary_path, term2gene_path and category_map_path together."
        ),
        call. = FALSE
      )
    }
    return(list(
      lisa_dictionary_path = lisa_dictionary_path,
      term2gene_path = term2gene_path,
      category_map_path = category_map_path,
      registered_category_map = FALSE,
      source = "explicit_paths",
      resources = NULL
    ))
  }

  if (is.null(term2gene_path)) {
    selection <- lisa_select_pipeline_resource_references(list(profile = profile))
    reference <- selection$references$term2gene_resource
    context <- lisa_active_resource_context()
    resource <- lisa_resolve_dictionary_resource(
      logical_id = reference$logical_id,
      version = reference$version,
      species = species,
      modality = profile,
      registry = context$registry,
      cache_root = context$cache_root,
      shared_root = context$shared_root
    )
    if (!identical(resource$schema, "term2gene@1")) {
      stop(
        "LISA-RESOURCE-017 the default HALLMARKS membership resource is not term2gene@1. Repair: install the registered TERM2GENE resource for this species.",
        call. = FALSE
      )
    }
    term2gene <- lisa_read_custom_resource_tsv(
      resource, "term2gene_resource"
    )
    lisa_validate_custom_term2gene(term2gene, resource$resource_id)
    term2gene_path <- resource$path
    resources <- data.frame(
      config_key = "term2gene_resource",
      resource_id = resource$resource_id,
      species = resource$species,
      modality = resource$modality,
      schema = resource$schema,
      sha256 = resource$sha256,
      path = resource$path,
      stringsAsFactors = FALSE
    )
  } else {
    resources <- NULL
  }
  list(
    lisa_dictionary_path = lisa_dictionary_path,
    term2gene_path = term2gene_path,
    category_map_path = category_map_path,
    registered_category_map = FALSE,
    source = if (is.null(resources)) "explicit_paths" else "registered_term2gene",
    resources = resources
  )
}

lisa_parse_resource_reference <- function(value, config_key,
                                          require_version = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop(sprintf(
      "LISA-RESOURCE-014 pipeline.%s is required. Repair: provide a registered logical ID such as name@version.",
      config_key
    ), call. = FALSE)
  }
  parts <- strsplit(trimws(value), "@", fixed = TRUE)[[1L]]
  if (!length(parts) || length(parts) > 2L || any(!nzchar(parts))) {
    stop(sprintf(
      "LISA-RESOURCE-016 pipeline.%s has malformed resource ID '%s'. Repair: use one registered name or name@version reference.",
      config_key, value
    ), call. = FALSE)
  }
  reference <- list(
    logical_id = parts[[1L]],
    version = if (length(parts) == 2L) parts[[2L]] else NULL,
    configured_id = trimws(value)
  )
  tryCatch({
    reference$logical_id <- lisa_safe_id(
      reference$logical_id, paste0("pipeline.", config_key, " logical_id")
    )
    if (!is.null(reference$version)) {
      reference$version <- lisa_safe_id(
        reference$version, paste0("pipeline.", config_key, " version")
      )
    }
  }, error = function(error) {
    stop(sprintf(
      "LISA-RESOURCE-016 pipeline.%s has unsafe resource ID '%s'. Repair: use portable logical_id@version syntax.",
      config_key, value
    ), call. = FALSE)
  })
  if (isTRUE(require_version) && is.null(reference$version)) {
    stop(sprintf(
      "LISA-RESOURCE-019 explicit resource '%s' has no version. Repair: configure pipeline.%s as logical_id@version.",
      reference$logical_id, config_key
    ), call. = FALSE)
  }
  reference
}

lisa_dictionary_tier_for_resource <- function(resource_id) {
  logical_id <- strsplit(as.character(resource_id)[[1L]], "@", fixed = TRUE)[[1L]][[1L]]
  tiers <- lisa_standard_dictionary_tiers()
  if (!logical_id %in% names(tiers)) return(NULL)
  unname(tiers[[logical_id]])
}

lisa_select_pipeline_resource_references <- function(pipeline) {
  if (is.null(pipeline)) pipeline <- list()
  if (!is.list(pipeline)) {
    stop(
      "LISA-RESOURCE-016 pipeline resource configuration must be a named list. Repair: provide a valid pipeline object.",
      call. = FALSE
    )
  }
  tier_is_explicit <- !is.null(pipeline$lisa_dictionary)
  tier <- if (tier_is_explicit) {
    value <- tolower(trimws(as.character(pipeline$lisa_dictionary)[[1L]]))
    if (!value %in% c("core", "expanded")) {
      stop(sprintf(
        "LISA-RESOURCE-018 pipeline.lisa_dictionary '%s' is unsupported. Repair: use core or expanded, or omit it for a custom dictionary_resource.",
        value
      ), call. = FALSE)
    }
    value
  } else {
    "core"
  }
  defaults <- lisa_pipeline_resource_defaults(tier)
  references <- stats::setNames(vector("list", length(defaults)), names(defaults))

  dictionary_value <- pipeline$dictionary_resource
  if (is.null(dictionary_value)) {
    reference <- lisa_parse_resource_reference(
      defaults[["dictionary_resource"]], "dictionary_resource",
      require_version = TRUE
    )
    reference$configured_id <- if (tier_is_explicit) {
      paste0("lisa_dictionary:", tier)
    } else {
      ""
    }
    reference$selection_source <- if (tier_is_explicit) "tier" else "default"
    effective_tier <- tier
  } else {
    reference <- lisa_parse_resource_reference(
      dictionary_value, "dictionary_resource", require_version = TRUE
    )
    resource_id <- lisa_resource_id(reference$logical_id, reference$version)
    derived_tier <- lisa_dictionary_tier_for_resource(resource_id)
    if (tier_is_explicit &&
        (is.null(derived_tier) || !identical(derived_tier, tier))) {
      qualifier <- if (is.null(derived_tier)) "custom " else ""
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-018 pipeline.lisa_dictionary '%s' conflicts with %sdictionary_resource '%s'. ",
          "Repair: remove lisa_dictionary for a custom resource, or select the resource matching that tier."
        ), tier, qualifier, resource_id
      ), call. = FALSE)
    }
    effective_tier <- if (is.null(derived_tier)) "custom" else derived_tier
    reference$selection_source <- "explicit"
  }
  reference$resource_id <- lisa_resource_id(
    reference$logical_id, reference$version
  )
  references[["dictionary_resource"]] <- reference

  for (key in c("term2gene_resource", "category_map_resource")) {
    configured <- pipeline[[key]]
    if (is.null(configured)) {
      item <- lisa_parse_resource_reference(
        defaults[[key]], key, require_version = TRUE
      )
      item$configured_id <- ""
      item$selection_source <- "default"
    } else {
      item <- lisa_parse_resource_reference(
        configured, key, require_version = TRUE
      )
      item$selection_source <- "explicit"
    }
    item$resource_id <- lisa_resource_id(item$logical_id, item$version)
    references[[key]] <- item
  }
  list(references = references, dictionary_tier = effective_tier)
}

lisa_resolve_pipeline_resources <- function(
  pipeline,
  species,
  registry = lisa_active_resource_context()$registry,
  cache_root = lisa_active_resource_context()$cache_root,
  shared_root = lisa_active_resource_context()$shared_root,
  stop_on_error = TRUE
) {
  spec <- lisa_pipeline_resource_spec()
  selection <- lisa_select_pipeline_resource_references(pipeline)
  profile <- lisa_config_get(pipeline, "profile", NULL)
  if (!is.character(profile) || length(profile) != 1L || is.na(profile) ||
      !nzchar(trimws(profile))) {
    stop(
      "LISA-RESOURCE-021 pipeline.profile is required to select modality-matched resources. Repair: declare the normalized scientific profile.",
      call. = FALSE
    )
  }
  profile <- tolower(trimws(profile))
  effective_default_tier <- if (selection$dictionary_tier %in%
                                c("core", "expanded")) {
    selection$dictionary_tier
  } else {
    "core"
  }
  effective_defaults <- lisa_pipeline_resource_defaults(effective_default_tier)
  resolved <- stats::setNames(vector("list", nrow(spec)), spec$config_key)
  status <- vector("list", nrow(spec))

  resolve_one <- function(i) {
    key <- spec$config_key[[i]]
    reference <- selection$references[[key]]
    resource <- tryCatch(
      lisa_resolve_dictionary_resource(
        reference$logical_id,
        version = reference$version,
        species = species,
        modality = profile,
        registry = registry,
        cache_root = cache_root,
        shared_root = shared_root
      ),
      error = function(error) {
        message <- conditionMessage(error)
        is_effective_default <- identical(
          reference$resource_id,
          unname(effective_defaults[[key]])
        )
        if (is_effective_default &&
            grepl("^LISA-RESOURCE-(006|008)", message)) {
          label <- if (identical(key, "term2gene_resource")) {
            "default TERM2GENE"
          } else {
            paste("default", key)
          }
          stop(sprintf(
            paste0(
              "LISA-RESOURCE-023 %s '%s' for species '%s' and modality '%s' is not installed. ",
              "Repair: run lisaR::install_lisa_resource(...) once with the authorised local TSV, or configure the administrator registry/cache, then rerun validate_lisa_config()."
            ), label, reference$resource_id,
            lisa_species_contract(species)$scientific_name, profile
          ), call. = FALSE)
        }
        stop(error)
      }
    )
    accepted <- lisa_role_accepted_schemas(key)
    if (!resource$schema %in% accepted) {
      stop(sprintf(
        "LISA-RESOURCE-017 pipeline.%s resolved schema '%s', expected %s. Repair: select the registered resource for this role.",
        key, resource$schema,
        paste(sprintf("'%s'", accepted), collapse = " or ")
      ), call. = FALSE)
    }
    resource$configured_id <- reference$configured_id
    resource$selection_source <- reference$selection_source
    resolved[[key]] <<- resource
    data.frame(
      config_key = key,
      configured_id = reference$configured_id,
      selection_source = reference$selection_source,
      resource_id = resource$resource_id,
      species = as.character(resource$species),
      modality = as.character(resource$modality),
      schema = as.character(resource$schema),
      compatibility = as.character(resource$compatibility),
      sha256 = as.character(resource$sha256),
      path = as.character(resource$path),
      checked = TRUE,
      ready = TRUE,
      status = "ready",
      message = "",
      stringsAsFactors = FALSE
    )
  }

  for (i in seq_len(nrow(spec))) {
    key <- spec$config_key[[i]]
    if (isTRUE(stop_on_error)) {
      status[[i]] <- resolve_one(i)
      next
    }
    status[[i]] <- tryCatch(
      resolve_one(i),
      error = function(error) {
        reference <- selection$references[[key]]
        data.frame(
          config_key = key,
          configured_id = reference$configured_id,
          selection_source = reference$selection_source,
          resource_id = reference$resource_id,
          species = lisa_species_contract(species)$scientific_name,
          modality = profile,
          schema = spec$expected_schema[[i]],
          compatibility = NA_character_,
          sha256 = NA_character_,
          path = NA_character_,
          checked = TRUE,
          ready = FALSE,
          status = "not_ready",
          message = conditionMessage(error),
          stringsAsFactors = FALSE
        )
      }
    )
  }
  status <- do.call(rbind, status)
  row.names(status) <- NULL
  result <- list(
    ready = nrow(status) == nrow(spec) && all(status$ready),
    status = status,
    resolved = resolved,
    dictionary_tier = selection$dictionary_tier,
    selection = selection$references
  )
  if (isTRUE(result$ready)) {
    content_validation <- tryCatch(
      lisa_validate_resolved_bundle(
        result$resolved, result$status,
        lisa_species_contract(species)$scientific_name,
        dictionary_tier = result$dictionary_tier
      ),
      error = function(error) {
        if (isTRUE(stop_on_error)) stop(error)
        error
      }
    )
    if (inherits(content_validation, "error")) {
      result$ready <- FALSE
      result$status$ready <- FALSE
      result$status$status <- "not_ready"
      result$status$message <- conditionMessage(content_validation)
    } else {
      result$content_validation <- content_validation
      if (identical(result$dictionary_tier, "custom")) {
        result$custom_validation <- content_validation
      }
    }
  }
  result
}

lisa_custom_abort <- function(code, role, resource_id, reason, repair) {
  stop(sprintf(
    "%s role=%s resource=%s; %s Repair: %s",
    code, role, resource_id, reason, repair
  ), call. = FALSE)
}

lisa_read_custom_resource_tsv <- function(resource, role) {
  tryCatch(
    utils::read.delim(
      resource$path,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character",
      na.strings = character()
    ),
    error = function(error) {
      lisa_custom_abort(
        "LISA-CUSTOM-001", role, resource$resource_id,
        paste0("TSV could not be read: ", conditionMessage(error), "."),
        "write one non-empty UTF-8 tab-separated table with a header."
      )
    }
  )
}

lisa_custom_assert_table <- function(x, required, role, resource_id,
                                     allow_extra = FALSE) {
  missing <- setdiff(required, names(x))
  extra <- setdiff(names(x), required)
  duplicated_names <- unique(names(x)[duplicated(names(x))])
  if (!is.data.frame(x) || !nrow(x) || length(missing) ||
      length(duplicated_names) || (!isTRUE(allow_extra) && length(extra))) {
    details <- c(
      if (!is.data.frame(x) || !nrow(x)) "table is empty" else character(),
      if (length(missing)) paste0("missing columns: ", paste(missing, collapse = ", ")) else character(),
      if (length(extra)) paste0("unsupported columns: ", paste(extra, collapse = ", ")) else character(),
      if (length(duplicated_names)) paste0("duplicate column names: ", paste(duplicated_names, collapse = ", ")) else character()
    )
    lisa_custom_abort(
      "LISA-CUSTOM-001", role, resource_id,
      paste0(paste(details, collapse = "; "), "."),
      paste0("provide the schema columns: ", paste(required, collapse = ", "), ".")
    )
  }
  invisible(TRUE)
}

lisa_custom_assert_nonempty <- function(x, columns, role, resource_id) {
  for (column in columns) {
    value <- as.character(x[[column]])
    bad <- is.na(value) | !nzchar(trimws(value)) | value != trimws(value)
    if (any(bad)) {
      lisa_custom_abort(
        "LISA-CUSTOM-001", role, resource_id,
        sprintf("required field '%s' contains an empty or non-canonical value.", column),
        sprintf("fill '%s' without leading or trailing whitespace.", column)
      )
    }
  }
  invisible(TRUE)
}

lisa_custom_assert_safe_ids <- function(x, columns, role, resource_id) {
  for (column in columns) {
    values <- unique(as.character(x[[column]]))
    valid <- vapply(values, function(value) {
      tryCatch({
        lisa_safe_id(value, paste0(role, " ", column))
        TRUE
      }, error = function(error) FALSE)
    }, logical(1))
    if (!all(valid)) {
      lisa_custom_abort(
        "LISA-CUSTOM-001", role, resource_id,
        sprintf("technical field '%s' contains a non-portable ID.", column),
        paste0(
          "use ASCII letters, numbers, dot, hyphen or underscore in '",
          column, "', and keep human-readable text in its label column."
        )
      )
    }
  }
  invisible(TRUE)
}

lisa_validate_dictionary_shape <- function(x, resource_id, schema = NULL) {
  role <- "dictionary_resource"
  # Transitional reader: accept both dictionary schemas. `lisa_dictionary@2` is
  # the canonical seven-column contract; `lisa_dictionary@1` is the same table
  # plus the historical `LISA_score`. Presence of that one column is the only
  # permitted difference, and it is still exact -- no other extra column is
  # tolerated, because `allow_extra` stays FALSE either way.
  #
  # When the caller knows the declared schema (the installer and the registry
  # do), the table must match it. Otherwise a resource could be registered as
  # `@2` while carrying `@1` bytes, leaving the registry telling a reviewer
  # something the file contradicts. When `schema` is NULL the shape is inferred,
  # which is what direct validation helpers and ad-hoc tables need.
  runtime <- lisa_dictionary_runtime_columns()
  legacy <- lisa_dictionary_has_legacy_score(x)
  if (!is.null(schema)) {
    declared_legacy <- identical(schema, "lisa_dictionary@1")
    if (!identical(legacy, declared_legacy)) {
      lisa_custom_abort(
        "LISA-CUSTOM-001", role, resource_id,
        sprintf(
          "table declares schema '%s' but %s a LISA_score column.",
          schema, if (legacy) "carries" else "omits"
        ),
        if (declared_legacy) {
          "declare schema 'lisa_dictionary@2' for a score-free dictionary, or restore the legacy column."
        } else {
          "declare schema 'lisa_dictionary@1' for a legacy dictionary, or delete the LISA_score column."
        }
      )
    }
    legacy <- declared_legacy
  }
  required <- if (legacy) {
    # Keep the historical column order so a legacy resource is reported against
    # the schema it actually declares.
    append(runtime, "LISA_score", after = match("category_display_name", runtime))
  } else {
    runtime
  }
  lisa_custom_assert_table(x, required, role, resource_id)
  lisa_custom_assert_nonempty(
    x,
    c("universe", "gene_set_id", "category_id", "category_display_name", "tier"),
    role, resource_id
  )
  lisa_custom_assert_safe_ids(x, "category_id", role, resource_id)
  semantic_universes <- c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
  bad_universes <- setdiff(unique(x$universe), semantic_universes)
  if (length(bad_universes)) {
    lisa_custom_abort(
      "LISA-CUSTOM-003", role, resource_id,
      paste0("unsupported dictionary universe(s): ", paste(bad_universes, collapse = ", "), "."),
      paste0("use only ", paste(semantic_universes, collapse = ", "), "; HALLMARKS is direct TERM2GENE content.")
    )
  }
  assignment <- paste(x$universe, x$gene_set_id, x$category_id, sep = "\r")
  if (anyDuplicated(assignment)) {
    lisa_custom_abort(
      "LISA-CUSTOM-002", role, resource_id,
      "duplicate (universe, gene_set_id, category_id) assignment detected.",
      "retain each exact assignment once; many-to-many relations remain supported."
    )
  }
  invisible(TRUE)
}

lisa_validate_custom_dictionary <- function(x, resource_id, schema = NULL) {
  role <- "dictionary_resource"
  lisa_validate_dictionary_shape(x, resource_id, schema = schema)
  if (any(x$tier != "custom")) {
    lisa_custom_abort(
      "LISA-CUSTOM-003", role, resource_id,
      "tier must be 'custom'; core and expanded are reserved reviewed tiers.",
      "set every custom dictionary tier value to custom."
    )
  }
  # The canonical contract has no score column at all, so a new custom
  # dictionary must not carry one -- not even an empty one. A legacy
  # `lisa_dictionary@1` custom resource keeps its original rule unchanged: the
  # column must be present and empty, and a populated value is still refused.
  if (lisa_dictionary_has_legacy_score(x)) {
    score <- as.character(x$LISA_score)
    if (any(!is.na(score) & nzchar(trimws(score)))) {
      lisa_custom_abort(
        "LISA-CUSTOM-003", role, resource_id,
        "LISA_score must be empty/NA for a legacy lisa_dictionary@1 custom dictionary and is never an analytical weight.",
        "leave LISA_score empty for every custom row, or migrate the resource to lisa_dictionary@2 by deleting the column."
      )
    }
  }
  invisible(TRUE)
}

lisa_validate_standard_dictionary <- function(x, resource_id,
                                              tier = c("core", "expanded"),
                                              schema = NULL) {
  tier <- match.arg(tier)
  lisa_validate_dictionary_shape(x, resource_id, schema = schema)
  # Tier is the contract. The retired `LISA_score` was a redundant encoding of
  # it (core <-> 4, expanded <-> 3) across both reviewed resources, so checking
  # tier accepts exactly the resources the score check accepted -- an
  # equivalence, not a relaxation.
  valid <- if (identical(tier, "core")) {
    x$tier == "core"
  } else {
    x$tier %in% c("core", "expanded")
  }
  if (any(!valid)) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-024 '%s' violates its reviewed %s tier contract. ",
        "Repair: install the exact versioned candidate."
      ), resource_id, tier
    ), call. = FALSE)
  }
  # A legacy lisa_dictionary@1 resource is additionally held to the historical
  # tier/score coherence, so an old resource is never checked less strictly than
  # before this migration.
  if (lisa_dictionary_has_legacy_score(x)) {
    score <- suppressWarnings(as.numeric(x$LISA_score))
    expected <- ifelse(x$tier == "core", 4, ifelse(x$tier == "expanded", 3, NA_real_))
    if (any(!is.finite(score) | is.na(expected) | score != expected)) {
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-024 legacy lisa_dictionary@1 resource '%s' violates the historical %s tier/score coherence. ",
          "Repair: install the exact versioned candidate, or migrate it to lisa_dictionary@2 by deleting the LISA_score column."
        ), resource_id, tier
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

lisa_validate_category_map <- function(x, resource_id) {
  standard <- identical(sub("@.*$", "", resource_id), "lisa_category_map")
  # The reviewed standard map combines the semantic plotting hierarchy and
  # the strict PATHWAYS hierarchy. Their order values are local to those two
  # scopes; neither the IDs nor the stored order values may be renumbered.
  scope <- rep.int("semantic", nrow(x))
  if (standard && any(c("family_id", "pathway_id") %in% names(x))) {
    if (!all(c("family_id", "pathway_id") %in% names(x))) {
      lisa_custom_abort("LISA-CUSTOM-004", "category_map_resource", resource_id,
        "standard hierarchy scope requires both family_id and pathway_id.",
        "retain the paired hierarchy identity columns from the reviewed map.")
    }
    family <- !is.na(x$family_id) & nzchar(trimws(x$family_id))
    pathway <- !is.na(x$pathway_id) & nzchar(trimws(x$pathway_id))
    if (any(family != pathway) ||
        any(family & (is.na(x$macrogroup_id) | is.na(x$category_id) |
          x$family_id != x$macrogroup_id | x$pathway_id != x$category_id))) {
      lisa_custom_abort("LISA-CUSTOM-004", "category_map_resource", resource_id,
        "standard family_id/pathway_id do not match their hierarchy identities.",
        "retain family_id=macrogroup_id and pathway_id=category_id for PATHWAYS rows; leave both empty otherwise.")
    }
    scope[family] <- "PATHWAYS"
  }
  lisa_validate_custom_category_map(x, resource_id, order_scope = scope)
}

lisa_validate_custom_category_map <- function(x, resource_id,
                                              order_scope = rep.int("global", nrow(x))) {
  role <- "category_map_resource"
  required <- c(
    "category_id", "display_name", "macrogroup_id", "macrogroup_name",
    "macrogroup_order", "category_order_within_macrogroup"
  )
  lisa_custom_assert_table(
    x, required, role, resource_id, allow_extra = TRUE
  )
  lisa_custom_assert_nonempty(
    x, required, role, resource_id
  )
  lisa_custom_assert_safe_ids(
    x, c("category_id", "macrogroup_id"), role, resource_id
  )
  duplicated_category <- duplicated(x$category_id) |
    duplicated(x$category_id, fromLast = TRUE)
  if (any(duplicated_category)) {
    groups <- split(x[duplicated_category, , drop = FALSE],
                    x$category_id[duplicated_category])
    conflicting <- vapply(groups, function(rows) {
      length(unique(paste(
        rows$macrogroup_id, rows$macrogroup_name, sep = "\r"
      ))) > 1L
    }, logical(1))
    if (any(conflicting)) {
      lisa_custom_abort(
        "LISA-CUSTOM-004", role, resource_id,
        "one category_id maps to multiple supercategories.",
        "retain exactly one macrogroup_id/macrogroup_name row per category_id."
      )
    }
    lisa_custom_abort(
      "LISA-CUSTOM-002", role, resource_id,
      "duplicate category_id row detected.",
      "retain exactly one registered row per category_id."
    )
  }
  order_columns <- c("macrogroup_order", "category_order_within_macrogroup")
  for (column in order_columns) {
    value <- as.character(x[[column]])
    if (any(!grepl("^[1-9][0-9]*$", value))) {
      lisa_custom_abort(
        "LISA-CUSTOM-004", role, resource_id,
        sprintf("'%s' contains a non-positive or non-whole number.", column),
        sprintf("encode '%s' as positive whole numbers.", column)
      )
    }
  }
  group_rows <- split(seq_len(nrow(x)), x$macrogroup_id)
  inconsistent_group <- names(group_rows)[vapply(group_rows, function(index) {
    length(unique(x$macrogroup_name[index])) != 1L ||
      length(unique(x$macrogroup_order[index])) != 1L ||
      length(unique(order_scope[index])) != 1L
  }, logical(1))]
  if (length(inconsistent_group)) {
    lisa_custom_abort(
      "LISA-CUSTOM-004", role, resource_id,
      paste0(
        "macrogroup_id has inconsistent name or order: ",
        paste(inconsistent_group, collapse = ", "), "."
      ),
      "use one macrogroup_name and one macrogroup_order for each macrogroup_id."
    )
  }
  group_identity <- vapply(group_rows, function(index) {
    paste(order_scope[index[[1L]]], x$macrogroup_order[index[[1L]]], sep = "\r")
  }, character(1))
  if (anyDuplicated(group_identity)) {
    lisa_custom_abort(
      "LISA-CUSTOM-004", role, resource_id,
      "multiple macrogroup_id values share one macrogroup_order within an ordering scope.",
      "assign one unique positive macrogroup_order to each supercategory within its ordering scope."
    )
  }
  within_group_order <- paste(
    x$macrogroup_id, x$category_order_within_macrogroup, sep = "\r"
  )
  if (anyDuplicated(within_group_order)) {
    lisa_custom_abort(
      "LISA-CUSTOM-004", role, resource_id,
      "category_order_within_macrogroup is duplicate within a macrogroup_id.",
      "assign one unique positive category order within each supercategory."
    )
  }
  invisible(TRUE)
}

lisa_validate_custom_term2gene <- function(x, resource_id) {
  role <- "term2gene_resource"
  required <- c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source",
    "gene_symbol"
  )
  lisa_custom_assert_table(x, required, role, resource_id)
  lisa_custom_assert_nonempty(
    x, c("gs_collection", "gs_name", "gene_symbol"), role, resource_id
  )
  key_values <- lapply(x[required], function(value) {
    value <- as.character(value)
    value[is.na(value)] <- ""
    value
  })
  membership <- do.call(paste, c(key_values, list(sep = "\r")))
  if (anyDuplicated(membership)) {
    lisa_custom_abort(
      "LISA-CUSTOM-002", role, resource_id,
      "TERM2GENE contains a duplicate exact membership row.",
      "retain each exact membership row once."
    )
  }
  metadata_columns <- c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source"
  )
  metadata <- unique(x[metadata_columns])
  conflicting_names <- names(which(table(metadata$gs_name) > 1L))
  if (length(conflicting_names)) {
    lisa_custom_abort(
      "LISA-CUSTOM-006", role, resource_id,
      paste0(
        "gs_name has inconsistent collection/provenance metadata: ",
        paste(conflicting_names, collapse = ", "), "."
      ),
      paste0(
        "use one gs_collection, gs_subcollection and gs_exact_source triple ",
        "per gs_name while retaining all of its gene_symbol memberships."
      )
    )
  }
  invisible(TRUE)
}

lisa_validate_custom_resource_relations <- function(dictionary, category_map,
                                                    term2gene, resource_ids) {
  category_labels <- split(
    dictionary$category_display_name, dictionary$category_id
  )
  inconsistent <- names(category_labels)[vapply(
    category_labels, function(value) length(unique(value)) != 1L, logical(1)
  )]
  if (length(inconsistent)) {
    lisa_custom_abort(
      "LISA-CUSTOM-005", "dictionary_resource",
      resource_ids[["dictionary_resource"]],
      paste0("category_id has inconsistent display labels: ", paste(inconsistent, collapse = ", "), "."),
      "use one identical category display label throughout the dictionary."
    )
  }
  missing_categories <- setdiff(
    unique(dictionary$category_id), category_map$category_id
  )
  if (length(missing_categories)) {
    lisa_custom_abort(
      "LISA-CUSTOM-005", "category_map_resource",
      resource_ids[["category_map_resource"]],
      paste0("orphan dictionary category_id(s): ", paste(missing_categories, collapse = ", "), "."),
      "add one registered category-map row for every dictionary category_id."
    )
  }
  dictionary_label <- vapply(category_labels, `[[`, character(1), 1L)
  map_label <- stats::setNames(
    category_map$display_name, category_map$category_id
  )
  mismatched_labels <- names(dictionary_label)[
    dictionary_label != unname(map_label[names(dictionary_label)])
  ]
  if (length(mismatched_labels)) {
    lisa_custom_abort(
      "LISA-CUSTOM-005", "category_map_resource",
      resource_ids[["category_map_resource"]],
      paste0("category display name differs from the dictionary for: ", paste(mismatched_labels, collapse = ", "), "."),
      "make dictionary category_display_name and category-map display_name identical."
    )
  }
  missing_gene_sets <- setdiff(
    unique(dictionary$gene_set_id), unique(term2gene$gs_name)
  )
  if (length(missing_gene_sets)) {
    lisa_custom_abort(
      "LISA-CUSTOM-006", "term2gene_resource",
      resource_ids[["term2gene_resource"]],
      paste0("dictionary gene_set_id is absent from TERM2GENE: ", paste(missing_gene_sets, collapse = ", "), "."),
      "add species-matched memberships for every dictionary gene_set_id."
    )
  }
  invisible(TRUE)
}

lisa_validate_resolved_bundle <- function(resolved, resource_status,
                                          canonical_species,
                                          dictionary_tier = "custom") {
  modalities <- unique(vapply(resolved, function(resource) {
    as.character(resource$modality)[[1L]]
  }, character(1)))
  exact_modalities <- setdiff(modalities, "all")
  if (!length(modalities) || any(!nzchar(modalities)) ||
      length(exact_modalities) > 1L) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-021 resource bundle has incompatible modalities: %s. ",
        "Repair: use one scientific profile, with optional profile-independent 'all' resources."
      ),
      paste(modalities, collapse = ", ")
    ), call. = FALSE)
  }

  dictionary <- lisa_read_custom_resource_tsv(
    resolved$dictionary_resource, "dictionary_resource"
  )
  category_map <- lisa_read_custom_resource_tsv(
    resolved$category_map_resource, "category_map_resource"
  )
  term2gene <- lisa_read_custom_resource_tsv(
    resolved$term2gene_resource, "term2gene_resource"
  )
  # The registry's declared schema is authoritative here: a resolved resource
  # that does not match what it claims must fail preflight, not be silently
  # reinterpreted.
  dictionary_schema <- as.character(resolved$dictionary_resource$schema)
  if (identical(dictionary_tier, "custom")) {
    lisa_validate_custom_dictionary(
      dictionary, resolved$dictionary_resource$resource_id,
      schema = dictionary_schema
    )
  } else {
    lisa_validate_standard_dictionary(
      dictionary, resolved$dictionary_resource$resource_id,
      tier = dictionary_tier, schema = dictionary_schema
    )
  }
  lisa_validate_category_map(
    category_map, resolved$category_map_resource$resource_id
  )
  lisa_validate_custom_term2gene(
    term2gene, resolved$term2gene_resource$resource_id
  )
  resource_ids <- vapply(resolved, `[[`, character(1), "resource_id")
  lisa_validate_custom_resource_relations(
    dictionary, category_map, term2gene, resource_ids
  )
  row.names(resource_status) <- NULL
  list(
    valid = TRUE,
    species = canonical_species,
    resources = resource_status,
    summary = list(
      n_dictionary_assignments = as.integer(nrow(dictionary)),
      n_categories = as.integer(nrow(category_map)),
      n_gene_sets = as.integer(length(unique(dictionary$gene_set_id))),
      n_term2gene_memberships = as.integer(nrow(term2gene)),
      n_genes = as.integer(length(unique(term2gene$gene_symbol)))
    )
  )
}

#' Validate a registered custom LISA resource bundle
#'
#' Resolves a versioned custom dictionary, category map and TERM2GENE resource
#' from the active lisaR registry and managed cache. It verifies registry
#' identity, species, modality, schema, compatibility and SHA-256 before
#' validating the three TSV tables and their cross-resource relations. It
#' never downloads data and does not accept direct per-run file paths.
#'
#' This standalone validator has no pipeline profile argument. Each reference
#' must therefore resolve unambiguously for the requested species: a unique
#' exact-modality row is accepted, and a profile-independent `modality =
#' "all"` row is preferred when several modality rows share the same
#' reference. The resolved bundle may contain `all` resources plus at most one
#' exact modality.
#'
#' @param dictionary_resource Versioned registered `lisa_dictionary@2`
#'   reference in `logical_id@version` form.
#' @param category_map_resource Versioned registered `category_map@1`
#'   reference in `logical_id@version` form.
#' @param term2gene_resource Versioned registered `term2gene@1` reference in
#'   `logical_id@version` form.
#' @param species Species represented by the registered resources: `"Homo
#'   sapiens"` or `"Mus musculus"`.
#'
#' @return A list containing `valid = TRUE`, canonical `species`, a resource
#'   status table and a compact row/set summary. Invalid bundles stop with an
#'   actionable `LISA-CUSTOM-*` or `LISA-RESOURCE-*` error.
#' @export
validate_lisa_custom_resources <- function(
  dictionary_resource,
  category_map_resource,
  term2gene_resource,
  species
) {
  canonical_species <- lisa_species_contract(species)$scientific_name
  context <- lisa_active_resource_context()
  references <- list(
    dictionary_resource = dictionary_resource,
    category_map_resource = category_map_resource,
    term2gene_resource = term2gene_resource
  )
  resolved <- stats::setNames(vector("list", length(references)), names(references))
  rows <- vector("list", length(references))
  for (i in seq_along(references)) {
    role <- names(references)[[i]]
    reference <- lisa_parse_resource_reference(
      references[[role]], role, require_version = TRUE
    )
    resource <- lisa_resolve_dictionary_resource(
      reference$logical_id,
      version = reference$version,
      species = canonical_species,
      modality = NULL,
      registry = context$registry,
      cache_root = context$cache_root,
      shared_root = context$shared_root
    )
    accepted <- lisa_role_accepted_schemas(role)
    if (!resource$schema %in% accepted) {
      stop(sprintf(
        "LISA-RESOURCE-017 %s resolved schema '%s', expected %s. Repair: register the resource for its declared role.",
        role, resource$schema,
        paste(sprintf("'%s'", accepted), collapse = " or ")
      ), call. = FALSE)
    }
    resource$configured_id <- reference$configured_id
    resource$selection_source <- "explicit"
    resolved[[role]] <- resource
    rows[[i]] <- data.frame(
      config_key = role,
      configured_id = reference$configured_id,
      selection_source = "explicit",
      resource_id = resource$resource_id,
      species = as.character(resource$species),
      modality = as.character(resource$modality),
      schema = as.character(resource$schema),
      compatibility = as.character(resource$compatibility),
      sha256 = as.character(resource$sha256),
      path = as.character(resource$path),
      checked = TRUE,
      ready = TRUE,
      status = "ready",
      message = "",
      stringsAsFactors = FALSE
    )
  }
  resource_status <- do.call(rbind, rows)
  lisa_validate_resolved_bundle(
    resolved, resource_status, canonical_species,
    dictionary_tier = "custom"
  )
}

lisa_species_contract <- function(species, msigdb_mode = NULL, term2gene_target_species = NULL,
                                  kegg_code = NULL, annotation_db = NULL) {
  aliases <- c("human" = "Homo sapiens", "homo sapiens" = "Homo sapiens", "mouse" = "Mus musculus", "mus musculus" = "Mus musculus")
  key <- tolower(trimws(as.character(species[[1]])))
  scientific_name <- unname(aliases[key])
  if (!key %in% names(aliases) || is.na(scientific_name)) stop("LISA-SPECIES-001 species is missing or unsupported. Repair: use 'Homo sapiens' or 'Mus musculus'.", call. = FALSE)
  expected <- if (identical(scientific_name, "Homo sapiens")) list(msigdb_mode = "human", namespace = "HGNC symbol", annotation_db = "org.Hs.eg.db", kegg_code = "hsa") else list(msigdb_mode = "mouse_native", namespace = "MGI symbol", annotation_db = "org.Mm.eg.db", kegg_code = "mmu")
  supplied <- list(msigdb_mode = msigdb_mode, kegg_code = kegg_code, annotation_db = annotation_db)
  for (field in names(supplied)) if (!is.null(supplied[[field]]) && nzchar(supplied[[field]]) && !identical(as.character(supplied[[field]]), expected[[field]])) stop(sprintf("LISA-SPECIES-002 %s '%s' conflicts with species '%s'. Repair: use %s.", field, supplied[[field]], scientific_name, expected[[field]]), call. = FALSE)
  if (!is.null(term2gene_target_species) && nzchar(term2gene_target_species) && !identical(as.character(term2gene_target_species), scientific_name)) stop(sprintf("LISA-SPECIES-003 TERM2GENE target species '%s' conflicts with '%s'. Repair: install a matching resource.", term2gene_target_species, scientific_name), call. = FALSE)
  c(list(scientific_name = scientific_name), expected)
}

lisa_kegg_snapshot_path <- function(cache_root, snapshot_id, organism, resource_type, key) {
  snapshot_id <- lisa_safe_id(snapshot_id, "snapshot_id")
  organism <- lisa_safe_id(organism, "organism")
  resource_type <- lisa_safe_id(resource_type, "resource_type")
  key_component <- lisa_external_key_path_component(key, "KEGG resource key")
  file.path(
    normalizePath(cache_root, mustWork = FALSE), "kegg", snapshot_id,
    organism, resource_type, key_component
  )
}

lisa_read_kegg_snapshot <- function(cache_root, snapshot_id, organism, resource_type, key) {
  path <- lisa_kegg_snapshot_path(cache_root, snapshot_id, organism, resource_type, key)
  metadata <- paste0(path, ".metadata.tsv")
  if (!file.exists(path) || !file.exists(metadata)) stop(sprintf("LISA-KEGG-001 immutable KEGG resource is absent: %s. Repair: explicitly materialize snapshot '%s' before cache_only analysis.", path, snapshot_id), call. = FALSE)
  meta <- read_lisa_tsv(metadata)
  required <- c("snapshot_id", "organism", "resource_type", "key", "sha256", "retrieval_status")
  if (!all(required %in% names(meta)) || nrow(meta) != 1L ||
      !identical(meta$retrieval_status[[1]], "validated") ||
      !lisa_sha256_all_valid(meta$sha256, 1L)) stop("LISA-KEGG-002 KEGG snapshot metadata is invalid. Repair: refresh the resource into a new snapshot.", call. = FALSE)
  if (!identical(lisa_sha256_file(path), meta$sha256[[1]])) stop("LISA-KEGG-003 KEGG snapshot SHA-256 mismatch. Repair: discard the corrupt snapshot and explicitly refresh it.", call. = FALSE)
  list(path = path, metadata = meta)
}

lisa_store_kegg_snapshot <- function(source, cache_root, snapshot_id, organism, resource_type, key,
                                     query = "local_fixture", retrieved_at = as.character(Sys.time())) {
  if (!file.exists(source)) stop("LISA-KEGG-004 staged KEGG resource does not exist. Repair: provide a validated staged response.", call. = FALSE)
  cache_root <- lisa_run_root(cache_root)
  snapshot_id <- lisa_safe_id(snapshot_id, "snapshot_id")
  organism <- lisa_safe_id(organism, "organism")
  resource_type <- lisa_safe_id(resource_type, "resource_type")
  target <- lisa_kegg_snapshot_path(cache_root, snapshot_id, organism, resource_type, key)
  if (file.exists(target)) stop("LISA-KEGG-005 immutable KEGG snapshot already exists. Repair: use a new snapshot ID for refresh.", call. = FALSE)
  lisa_guarded_dir_create(dirname(target), cache_root)
  staged <- paste0(target, ".staging-", Sys.getpid())
  lisa_guarded_copy(source, staged, overwrite = FALSE, run_root = cache_root)
  sha <- lisa_sha256_file(staged)
  metadata <- data.frame(snapshot_id = snapshot_id, organism = organism, resource_type = resource_type,
    key = key, sha256 = sha, query = query, retrieved_at = retrieved_at, retrieval_status = "validated", stringsAsFactors = FALSE)
  meta_stage <- paste0(staged, ".metadata.tsv")
  lisa_guarded_write(
    meta_stage,
    function(path) utils::write.table(
      metadata, file = path, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = TRUE, na = ""
    ),
    run_root = cache_root
  )
  lisa_guarded_rename(staged, target, run_root = cache_root)
  lisa_guarded_rename(
    meta_stage, paste0(target, ".metadata.tsv"), run_root = cache_root
  )
  lisa_read_kegg_snapshot(cache_root, snapshot_id, organism, resource_type, key)
}

lisa_kegg_binary_response <- function(fetch) {
  encoding_note <- NA_character_
  value <- withCallingHandlers(
    fetch(),
    warning = function(w) {
      message <- conditionMessage(w)
      if (grepl(
        "^No encoding supplied: defaulting to UTF-8\\.?$",
        message
      )) {
        encoding_note <<- paste(
          "KEGG binary PNG response declared no text encoding;",
          "encoding is not applicable to the decoded image array."
        )
        invokeRestart("muffleWarning")
      }
    }
  )
  list(value = value, encoding_note = encoding_note)
}

lisa_fetch_kegg_resource <- function(organism, resource_type, key) {
  lisa_require_optional("KEGGREST", "retrieving KEGG pathway resources")
  if (!organism %in% c("hsa", "mmu")) {
    stop("LISA-KEGG-020 unsupported KEGG organism code.", call. = FALSE)
  }
  if (!resource_type %in% c("pathway_links", "pathway_list", "kgml", "image")) {
    stop("LISA-KEGG-021 unsupported KEGG resource type.", call. = FALSE)
  }

  staged <- tempfile(
    paste0("lisaR-kegg-", resource_type, "-"),
    fileext = if (identical(resource_type, "image")) ".png" else
      if (identical(resource_type, "kgml")) ".xml" else ".rds"
  )
  query <- ""
  if (identical(resource_type, "pathway_links")) {
    value <- KEGGREST::keggLink("pathway", organism)
    saveRDS(value, staged)
    query <- sprintf("KEGGREST::keggLink('pathway', '%s')", organism)
  } else if (identical(resource_type, "pathway_list")) {
    value <- KEGGREST::keggList("pathway", organism)
    saveRDS(value, staged)
    query <- sprintf("KEGGREST::keggList('pathway', '%s')", organism)
  } else {
    if (!grepl(paste0("^", organism, "[0-9]{5}$"), key)) {
      stop("LISA-KEGG-022 pathway key does not match the declared organism.",
           call. = FALSE)
    }
    if (identical(resource_type, "kgml")) {
      value <- KEGGREST::keggGet(key, option = "kgml")
      if (!is.character(value) || length(value) != 1L || !nzchar(value)) {
        stop("LISA-KEGG-023 KEGGREST returned invalid KGML.", call. = FALSE)
      }
      writeLines(value, staged, useBytes = TRUE)
    } else {
      lisa_require_optional("png", "writing a retrieved KEGG pathway image")
      binary <- lisa_kegg_binary_response(
        function() KEGGREST::keggGet(key, option = "image")
      )
      value <- binary$value
      if (!is.array(value) || length(dim(value)) < 2L) {
        stop("LISA-KEGG-024 KEGGREST returned an invalid pathway image.",
             call. = FALSE)
      }
      png::writePNG(value, staged)
      if (!is.na(binary$encoding_note)) {
        attr(staged, "transport_note") <- binary$encoding_note
      }
    }
    query <- sprintf("KEGGREST::keggGet('%s', option='%s')",
                     key, resource_type)
    if (identical(resource_type, "image") &&
        !is.null(attr(staged, "transport_note"))) {
      query <- paste(query, attr(staged, "transport_note"), sep = "; ")
    }
  }
  attr(staged, "query") <- query
  staged
}

lisa_resolve_kegg_resource <- function(mode = c("cache_only", "prefer_cache", "refresh"), cache_root,
                                       snapshot_id, organism, resource_type, key, fetch = NULL) {
  mode <- match.arg(mode)
  existing <- tryCatch(lisa_read_kegg_snapshot(cache_root, snapshot_id, organism, resource_type, key), error = function(e) e)
  if (!inherits(existing, "error") && mode != "refresh") return(existing)
  if (identical(mode, "cache_only")) stop(conditionMessage(existing), call. = FALSE)
  if (is.null(fetch) || !is.function(fetch)) stop("LISA-KEGG-008 live retrieval requires an explicit fetch function. Repair: use cache_only with a frozen snapshot or explicitly provide approved retrieval code.", call. = FALSE)
  staged <- fetch()
  if (!is.character(staged) || length(staged) != 1L) stop("LISA-KEGG-009 fetch did not return one staged file.", call. = FALSE)
  lisa_store_kegg_snapshot(
    staged, cache_root, snapshot_id, organism, resource_type, key,
    query = attr(staged, "query") %||% "explicit external fetch"
  )
}
