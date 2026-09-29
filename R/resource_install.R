# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

#' Install and register one local LISA scientific resource
#'
#' Copies one already-obtained resource into lisaR's managed local cache and
#' records its immutable identity in `resource_registry.tsv`. The helper never
#' downloads data, accepts licence terms, or grants redistribution rights.
#' Obtain the source through an authorised route and record that route in
#' `approved_origin` before calling this function.
#'
#' @param source Existing local TSV file. Symbolic links are rejected.
#' @param logical_id Stable registry identifier without a version suffix.
#' @param version Immutable resource version.
#' @param species Resource species, `"Homo sapiens"` or `"Mus musculus"`.
#' @param modality Scientific profile for which the resource was prepared.
#' @param schema One of `"lisa_dictionary@2"`, legacy
#'   `"lisa_dictionary@1"`, `"term2gene@1"` or `"category_map@1"`.
#' @param approved_origin Non-empty auditable description of the authorised
#'   source and local review. This is not a licence grant.
#' @param expected_sha256 Reviewed lower-case SHA-256 expected for `source`.
#'   The helper refuses self-certification by computing and accepting an
#'   otherwise unknown digest implicitly.
#' @param compatibility Closed lisaR version constraint.
#' @param cache_root Managed cache root. The default is shared with resource
#'   resolution.
#' @param registry_path Registry TSV to update. It must remain inside
#'   `cache_root`; the default is loaded automatically by subsequent sessions.
#'
#' @return The installed registry row plus canonical `path` and `resource_id`.
#' @export
install_lisa_resource <- function(
  source,
  logical_id,
  version,
  species,
  modality,
  schema,
  approved_origin,
  expected_sha256,
  compatibility = "lisaR>=0.6.0",
  cache_root = lisa_dictionary_cache_root(),
  registry_path = file.path(cache_root, "resource_registry.tsv")
) {
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !file.exists(source) || dir.exists(source) ||
      !isTRUE(utils::file_test("-f", source)) || lisa_path_is_link(source)) {
    stop(
      "LISA-RESOURCE-024 source must be one existing regular non-symlink TSV. Repair: provide the authorised local resource file.",
      call. = FALSE
    )
  }
  source <- normalizePath(source, winslash = "/", mustWork = TRUE)
  logical_id <- lisa_safe_id(logical_id, "resource logical_id")
  version <- lisa_safe_id(version, "resource version")
  species <- lisa_species_contract(species)$scientific_name
  if (!is.character(modality) || length(modality) != 1L || is.na(modality) ||
      !nzchar(trimws(modality)) || !identical(modality, trimws(modality)) ||
      !modality %in% lisa_resource_modalities()) {
    stop(
      "LISA-RESOURCE-024 modality must be 'all' or one canonical scientific profile. Repair: record the exact compatible pipeline.profile, or 'all' only for a reviewed profile-independent resource.",
      call. = FALSE
    )
  }
  # Both dictionary schemas install to the same cached artifact name: the file
  # layout is keyed by role, not by schema version, so a legacy resource and a
  # migrated one occupy the same slot for their own logical_id/version.
  schemas <- c(
    "lisa_dictionary@1" = "dictionary.tsv",
    "lisa_dictionary@2" = "dictionary.tsv",
    "term2gene@1" = "term2gene.tsv",
    "category_map@1" = "category_map.tsv"
  )
  if (!is.character(schema) || length(schema) != 1L || is.na(schema) ||
      !schema %in% names(schemas)) {
    stop(
      "LISA-RESOURCE-024 schema is unsupported. Repair: use lisa_dictionary@2 (or legacy lisa_dictionary@1), term2gene@1 or category_map@1.",
      call. = FALSE
    )
  }
  if (!is.character(approved_origin) || length(approved_origin) != 1L ||
      is.na(approved_origin) || !nzchar(trimws(approved_origin)) ||
      !identical(approved_origin, trimws(approved_origin)) ||
      grepl("[[:cntrl:]]", approved_origin)) {
    stop(
      "LISA-RESOURCE-024 approved_origin is required. Repair: record the authorised source and local review; do not use it as a licence field.",
      call. = FALSE
    )
  }
  if (!is.character(expected_sha256) || length(expected_sha256) != 1L ||
      is.na(expected_sha256) ||
      !lisa_sha256_all_valid(expected_sha256, 1L)) {
    stop(
      "LISA-RESOURCE-025 expected_sha256 must be one reviewed lower-case SHA-256. Repair: obtain the digest from the approved resource manifest; do not self-certify unknown bytes.",
      call. = FALSE
    )
  }
  source_sha256 <- lisa_sha256_file(source)
  if (!identical(source_sha256, expected_sha256)) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-025 source SHA-256 does not match expected_sha256 for '%s'. ",
        "Repair: obtain the reviewed immutable file or correct the expected digest from its approved manifest."
      ), lisa_resource_id(logical_id, version)
    ), call. = FALSE)
  }
  lisa_validate_resource_compatibility(compatibility)

  table <- tryCatch(
    utils::read.delim(
      source, sep = "\t", header = TRUE, quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
      na.strings = character()
    ),
    error = function(error) stop(
      "LISA-RESOURCE-024 source TSV could not be read: ",
      conditionMessage(error), ". Repair: provide a valid tab-separated file.",
      call. = FALSE
    )
  )
  resource_id <- lisa_resource_id(logical_id, version)
  if (schema %in% lisa_dictionary_schemas()) {
    standard_tiers <- lisa_standard_dictionary_tiers()
    if (logical_id %in% names(standard_tiers)) {
      lisa_validate_standard_dictionary(
        table, resource_id,
        tier = unname(standard_tiers[[logical_id]]),
        schema = schema
      )
    } else {
      lisa_validate_custom_dictionary(table, resource_id, schema = schema)
    }
  } else if (identical(schema, "category_map@1")) {
    lisa_validate_category_map(table, resource_id)
  } else {
    lisa_validate_custom_term2gene(table, resource_id)
  }
  if (!identical(lisa_sha256_file(source), source_sha256)) {
    stop(
      "LISA-RESOURCE-025 source changed while it was being validated. Repair: use an immutable local source and retry.",
      call. = FALSE
    )
  }

  cache_root <- lisa_path_walk_lexically(cache_root)$path
  lisa_guarded_dir_create(cache_root, run_root = NULL)
  cache_root <- lisa_resource_cache_root(cache_root)
  registry_path <- lisa_guarded_path(registry_path, cache_root)
  lisa_guarded_dir_create(dirname(registry_path), cache_root)
  species_artifact <- c(
    "Homo sapiens" = "homo_sapiens",
    "Mus musculus" = "mus_musculus"
  )[[species]]
  modality_artifact <- c(
    all = "all",
    `transcriptomic/genomic` = "transcriptomic-genomic",
    `global proteomic` = "global-proteomic",
    targeted = "targeted",
    custom = "custom"
  )[[modality]]
  row <- data.frame(
    logical_id = logical_id,
    version = version,
    species = species,
    modality = modality,
    schema = schema,
    sha256 = source_sha256,
    compatibility = compatibility,
    approved_origin = trimws(approved_origin),
    artifact = paste(
      unname(species_artifact), unname(modality_artifact),
      unname(schemas[[schema]]), sep = "/"
    ),
    stringsAsFactors = FALSE
  )
  # Validate the row before any write and reject attempts to replace a
  # synthetic resource installed with the package.
  lisa_read_dictionary_registry(row)
  lisa_active_resource_registry(row)

  lock_path <- lisa_guarded_path(
    paste0(registry_path, ".lock"), cache_root
  )
  if (lisa_path_entry_exists(lock_path) ||
      !dir.create(lock_path, recursive = FALSE, showWarnings = FALSE,
                  mode = "0700")) {
    stop(sprintf(
      paste0(
        "LISA-RESOURCE-026 resource registry is locked: %s. ",
        "Repair: wait for the active installer; remove the lock only after verifying that no installer is running."
      ), lock_path
    ), call. = FALSE)
  }
  lisa_set_private_mode(lock_path, directory = TRUE)
  lock_owned <- TRUE
  on.exit(if (lock_owned && lisa_path_entry_exists(lock_path)) {
    try(
      lisa_guarded_delete(lock_path, recursive = TRUE, run_root = cache_root),
      silent = TRUE
    )
  }, add = TRUE)

  registry <- if (file.exists(registry_path)) {
    lisa_read_dictionary_registry(registry_path)
  } else {
    lisa_empty_resource_registry()
  }
  key <- paste(
    row$logical_id, row$version, row$species, row$modality, sep = "\r"
  )
  existing_key <- paste(
    registry$logical_id, registry$version, registry$species,
    registry$modality, sep = "\r"
  )
  hit <- which(existing_key == key)
  if (length(hit)) {
    columns <- lisa_dictionary_registry_columns()
    identical_row <- all(vapply(columns, function(column) {
      identical(
        as.character(registry[[column]][hit[[1L]]]),
        as.character(row[[column]][[1L]])
      )
    }, logical(1)))
    if (!isTRUE(identical_row)) {
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-024 registry identity '%s' already exists with different metadata or SHA-256. ",
          "Repair: use a new version for changed bytes or metadata."
        ), resource_id
      ), call. = FALSE)
    }
  }

  target <- lisa_dictionary_resource_path(row, cache_root, shared_root = NULL)
  lisa_guarded_dir_create(dirname(target), cache_root)
  created_target <- FALSE
  committed <- FALSE
  on.exit(if (created_target && !committed && lisa_path_entry_exists(target)) {
    try(lisa_guarded_delete(target, recursive = FALSE, run_root = cache_root),
        silent = TRUE)
  }, add = TRUE)
  if (!identical(lisa_sha256_file(source), source_sha256)) {
    stop(
      "LISA-RESOURCE-025 source changed before installation. Repair: use an immutable local source and retry.",
      call. = FALSE
    )
  }
  if (lisa_path_entry_exists(target)) {
    if (dir.exists(target) || lisa_path_is_link(target) ||
        !identical(lisa_sha256_file(target), row$sha256[[1L]])) {
      stop(sprintf(
        paste0(
          "LISA-RESOURCE-024 immutable target '%s' already exists with a different identity. ",
          "Repair: use a new resource version; never overwrite registered bytes."
        ), resource_id
      ), call. = FALSE)
    }
  } else {
    lisa_guarded_copy(source, target, overwrite = FALSE, run_root = cache_root)
    created_target <- TRUE
  }
  if (!identical(lisa_sha256_file(target), source_sha256)) {
    stop(
      "LISA-RESOURCE-025 installed bytes do not match the reviewed SHA-256. Repair: secure the source and cache, then retry.",
      call. = FALSE
    )
  }

  if (!length(hit)) {
    registry <- rbind(registry, row)
    registry <- registry[order(
      registry$logical_id, registry$version, registry$species,
      registry$modality,
      method = "radix"
    ), , drop = FALSE]
    row.names(registry) <- NULL
    lisa_guarded_write(
      registry_path,
      function(path) utils::write.table(
        registry, path, sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = TRUE, na = ""
      ),
      run_root = cache_root,
      overwrite = TRUE
    )
  }
  committed <- TRUE
  c(
    as.list(row[1L, , drop = FALSE]),
    list(
      path = normalizePath(target, winslash = "/", mustWork = TRUE),
      resource_id = resource_id,
      registry_path = normalizePath(
        registry_path, winslash = "/", mustWork = TRUE
      )
    )
  )
}

lisa_msigdb_2026_1_term2gene_contract <- function() {
  list(
    source_url = "https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1",
    source_zip_sha256 = "84c8ff3270db946ced2a7a4beacf866a3d4ce138af2cbdb1dae76eacb4d4c576",
    output_sha256 = "716576959ec2941c1a7ecf63ff38e765667e3cd2eea3312499c1c2a05211d4a0",
    output_bytes = 125742425,
    max_archive_bytes = 2 * 1024^3,
    download_timeout_seconds = 900L,
    collections = c("H", "C2", "C5"),
    rds_files = c(
      "msigdb.2026.1.Hs.H.rds",
      "msigdb.2026.1.Hs.C2.rds",
      "msigdb.2026.1.Hs.C5.rds"
    )
  )
}

lisa_msigdb_2026_1_acquisition_terms <- function() {
  paste0(
    "MSigDB terms: https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp; ",
    "KEGG terms: https://www.kegg.jp/kegg/legal.html; BioCarta disclaimer/",
    "conditions are linked from the MSigDB terms."
  )
}

lisa_msigdb_2026_1_require_terms <- function(accept_terms) {
  if (isTRUE(accept_terms)) return(invisible(TRUE))
  stop(
    paste0(
      "LISA-RESOURCE-028 acquisition is not started because the recipient has not explicitly accepted the applicable MSigDB, KEGG and BioCarta terms. ",
      lisa_msigdb_2026_1_acquisition_terms(),
      " Review them for your intended use, then run prepare_lisa_msigdb_resource(accept_terms = TRUE). ",
      "This assertion is not a licence, permission grant, or acceptance on anyone else's behalf; alternatively supply an already-authorised local ZIP via source_zip."
    ),
    call. = FALSE
  )
}

lisa_msigdb_2026_1_term2gene_origin <- function(contract) {
  paste0(
    "MSigDB 2026.1 user-supplied local ZIP; ZIP SHA-256=",
    contract$source_zip_sha256,
    "; exact H,C2,C5 Homo sapiens RDS extraction; deterministic ",
    "distinct/order(gs_name,gene_symbol)/readr TSV transform; ",
    "canonical SHA-256=", contract$output_sha256,
    "; msigdbr 26.1.0 / Zenodo record 18968178 provenance; ",
    "see lisaR inst/THIRD_PARTY_NOTICES.md; no redistribution right granted"
  )
}

lisa_msigdb_2026_1_existing_resource <- function(registry_path, cache_root,
                                                   contract) {
  if (!file.exists(registry_path)) return(NULL)
  registry <- lisa_read_dictionary_registry(registry_path)
  hit <- registry[
    registry$logical_id == "msigdb_term2gene" &
      registry$version == "2026.1" &
      registry$species == "Homo sapiens" & registry$modality == "all",
    , drop = FALSE
  ]
  if (!nrow(hit)) return(NULL)
  expected <- data.frame(
    logical_id = "msigdb_term2gene", version = "2026.1",
    species = "Homo sapiens", modality = "all", schema = "term2gene@1",
    sha256 = contract$output_sha256, compatibility = "lisaR>=0.6.0",
    approved_origin = lisa_msigdb_2026_1_term2gene_origin(contract),
    artifact = "homo_sapiens/all/term2gene.tsv", stringsAsFactors = FALSE
  )
  columns <- lisa_dictionary_registry_columns()
  same <- all(vapply(columns, function(column) {
    identical(as.character(hit[[column]][[1L]]),
              as.character(expected[[column]][[1L]]))
  }, logical(1)))
  if (!isTRUE(same)) {
    stop(
      "LISA-RESOURCE-024 registry identity 'msigdb_term2gene@2026.1' already exists with different metadata or SHA-256. Repair: use a new version for changed bytes or metadata; never overwrite the pinned identity.",
      call. = FALSE
    )
  }
  lisa_resolve_dictionary_resource(
    "msigdb_term2gene", "2026.1", "Homo sapiens", modality = "all",
    registry = registry_path, cache_root = cache_root
  )
}

lisa_assert_msigdb_2026_1_zip <- function(source_zip, contract) {
  if (!is.character(source_zip) || length(source_zip) != 1L ||
      is.na(source_zip) || !file.exists(source_zip) || dir.exists(source_zip) ||
      !isTRUE(utils::file_test("-f", source_zip)) || lisa_path_is_link(source_zip)) {
    stop(
      "LISA-RESOURCE-028 source_zip must be one existing regular non-symlink MSigDB 2026.1 ZIP. Repair: supply the reviewed local ZIP; this helper never downloads it.",
      call. = FALSE
    )
  }
  source_zip <- normalizePath(source_zip, winslash = "/", mustWork = TRUE)
  observed <- lisa_sha256_file(source_zip)
  if (!identical(observed, contract$source_zip_sha256)) {
    stop(
      "LISA-RESOURCE-028 source_zip SHA-256 does not match the pinned MSigDB 2026.1 ZIP. Repair: use the exact reviewed local ZIP; do not substitute a later or differently packaged release.",
      call. = FALSE
    )
  }

  listing <- tryCatch(
    utils::unzip(source_zip, list = TRUE),
    error = function(error) stop(
      "LISA-RESOURCE-028 source_zip could not be listed safely: ",
      conditionMessage(error), call. = FALSE
    )
  )
  if (!is.data.frame(listing) || !all(c("Name", "Length") %in% names(listing)) ||
      anyDuplicated(listing$Name)) {
    stop(
      "LISA-RESOURCE-028 source_zip listing is malformed or contains duplicate members. Repair: obtain the exact reviewed ZIP.",
      call. = FALSE
    )
  }
  wanted <- match(contract$rds_files, listing$Name)
  if (anyNA(wanted)) {
    stop(
      "LISA-RESOURCE-028 source_zip does not contain the required MSigDB 2026.1 human H/C2/C5 RDS members. Repair: obtain the exact reviewed ZIP.",
      call. = FALSE
    )
  }
  lengths <- suppressWarnings(as.numeric(listing$Length[wanted]))
  # The fixed archive hash anchors identity; these independent bounds ensure
  # that only the three required members can be expanded into the private
  # staging directory and prevent an accidental unbounded extraction path.
  if (any(!is.finite(lengths)) || any(lengths <= 0) || any(lengths > 512 * 1024^2) ||
      sum(lengths) > 1024 * 1024^2) {
    stop(
      "LISA-RESOURCE-028 required ZIP members exceed the bounded extraction policy. Repair: obtain the exact reviewed ZIP.",
      call. = FALSE
    )
  }
  source_zip
}

lisa_msigdb_2026_1_archive_path <- function(cache_root) {
  file.path(cache_root, "archives", "msigdb.2026.1.zip")
}

lisa_msigdb_2026_1_download_file <- function(url, destination, contract) {
  old_timeout <- getOption("timeout")
  options(timeout = contract$download_timeout_seconds)
  on.exit(options(timeout = old_timeout), add = TRUE)
  status <- tryCatch(
    utils::download.file(url, destination, method = "libcurl", mode = "wb", quiet = TRUE),
    error = function(error) {
      stop(
        "LISA-RESOURCE-028 MSigDB 2026.1 download failed: ",
        conditionMessage(error),
        ". Repair: check connectivity and retry; existing registered resources were not changed.",
        call. = FALSE
      )
    }
  )
  if (!identical(as.integer(status), 0L)) {
    stop(
      "LISA-RESOURCE-028 MSigDB 2026.1 download returned a non-zero status. Repair: check connectivity and retry; existing registered resources were not changed.",
      call. = FALSE
    )
  }
  invisible(destination)
}

lisa_msigdb_2026_1_download_zip <- function(cache_root, contract) {
  archive <- lisa_msigdb_2026_1_archive_path(cache_root)
  if (lisa_path_entry_exists(archive)) {
    verified <- tryCatch(
      lisa_assert_msigdb_2026_1_zip(archive, contract),
      error = identity
    )
    if (!inherits(verified, "error")) return(verified)
    if (lisa_path_is_link(archive) || dir.exists(archive)) {
      stop(
        "LISA-RESOURCE-028 cached MSigDB archive is not one regular file. Repair: remove the invalid cache entry after review, then retry.",
        call. = FALSE
      )
    }
  }
  lisa_guarded_write(archive, function(path) {
    lisa_msigdb_2026_1_download_file(contract$source_url, path, contract)
    bytes <- unname(file.info(path)$size)
    if (!is.finite(bytes) || bytes <= 0 || bytes > contract$max_archive_bytes) {
      stop(
        "LISA-RESOURCE-028 downloaded MSigDB archive exceeds the bounded acquisition policy. Repair: retry the pinned source; no bytes were registered.",
        call. = FALSE
      )
    }
    lisa_assert_msigdb_2026_1_zip(path, contract)
  }, run_root = cache_root, overwrite = TRUE)
  lisa_assert_msigdb_2026_1_zip(archive, contract)
}

lisa_msigdb_2026_1_term2gene_table <- function(extract_dir, contract) {
  required <- c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source",
    "db_gene_symbol"
  )
  output_columns <- c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source",
    "gene_symbol"
  )
  tables <- lapply(contract$rds_files, function(member) {
    path <- file.path(extract_dir, member)
    if (!file.exists(path) || dir.exists(path) ||
        !isTRUE(utils::file_test("-f", path)) || lisa_path_is_link(path)) {
      stop(
        "LISA-RESOURCE-028 bounded ZIP extraction did not produce one regular required RDS member. Repair: obtain the exact reviewed ZIP.",
        call. = FALSE
      )
    }
    raw <- tryCatch(readRDS(path), error = function(error) stop(
      "LISA-RESOURCE-028 required MSigDB RDS member could not be read: ",
      conditionMessage(error), call. = FALSE
    ))
    if (!is.data.frame(raw) || !all(required %in% names(raw))) {
      stop(
        "LISA-RESOURCE-028 required MSigDB RDS member has an unexpected schema. Repair: obtain the exact reviewed ZIP.",
        call. = FALSE
      )
    }
    out <- as.data.frame(raw[, required, drop = FALSE], stringsAsFactors = FALSE)
    names(out) <- output_columns
    out
  })
  term2gene <- unique(do.call(rbind, tables))
  row.names(term2gene) <- NULL
  term2gene <- term2gene[order(
    term2gene$gs_name, term2gene$gene_symbol, method = "radix"
  ), output_columns, drop = FALSE]
  row.names(term2gene) <- NULL
  if (!identical(names(term2gene), output_columns) || !nrow(term2gene) ||
      anyNA(term2gene)) {
    stop(
      "LISA-RESOURCE-028 canonical TERM2GENE schema/content invariant failed. Repair: obtain the exact reviewed ZIP.",
      call. = FALSE
    )
  }
  term2gene
}

#' Prepare and register the pinned MSigDB 2026.1 TERM2GENE resource
#'
#' Converts the exact user-supplied local `msigdb.2026.1.zip` into the fixed
#' Homo sapiens H/C2/C5 `term2gene@1` resource, then delegates immutable cache
#' and registry installation to [install_lisa_resource()]. With no local ZIP,
#' it can acquire the one pinned public `msigdbr` 26.1.0 archive only after the
#' recipient explicitly asserts `accept_terms = TRUE`. That assertion is not a
#' licence, permission grant or acceptance for anyone else. Existing exact
#' installed resources resolve first, without a ZIP, terms assertion or network
#' access. The ZIP digest, canonical output digest, output size, schema and
#' transformation contract are internal; callers do not supply a hash. A new
#' download is written privately, bounded, verified before promotion, and then
#' transformed; failure leaves the existing registered resource untouched.
#'
#' @param source_zip Optional existing local MSigDB 2026.1 ZIP with the pinned
#'   digest. Symbolic links are rejected. This local route never asks for a new
#'   assertion of terms.
#' @param accept_terms Set to `TRUE` only after the recipient has reviewed and
#'   accepted the applicable MSigDB, KEGG and BioCarta terms for their intended
#'   use. It permits this call's public acquisition only when `source_zip` is
#'   `NULL`; it is not a licence or permission grant.
#' @param cache_root Managed cache root shared with resource resolution.
#' @param registry_path Registry TSV to update. It must remain inside
#'   `cache_root`; the default is used by subsequent normal resolution.
#'
#' @return The immutable registered row plus canonical `path`, `resource_id`,
#'   `registry_path`, source ZIP SHA-256 and transformation provenance.
#' @export
prepare_lisa_msigdb_resource <- function(
  source_zip = NULL,
  accept_terms = FALSE,
  cache_root = lisa_dictionary_cache_root(),
  registry_path = file.path(cache_root, "resource_registry.tsv")
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "LISA-RESOURCE-028 readr is required to reproduce the pinned canonical TSV bytes. Repair: reinstall lisaR with its required runtime dependencies; this helper will not install packages.",
      call. = FALSE
    )
  }
  contract <- lisa_msigdb_2026_1_term2gene_contract()
  cache_root <- lisa_resource_cache_root(
    lisa_path_walk_lexically(cache_root)$path
  )
  registry_path <- lisa_path_walk_lexically(registry_path)$path
  # `lisa_resource_cache_root()` returns the canonical spelling of the managed
  # cache, while `registry_path` is still the caller spelling. Compare
  # filesystem identity so a platform alias of the same directory (macOS
  # `/private/var`, a Windows 8.3 short name) is not read as an escape.
  if (is.null(lisa_path_relative_within(registry_path, cache_root))) {
    stop(
      "LISA-RESOURCE-028 registry_path must remain inside cache_root. Repair: use the canonical managed resource registry.",
      call. = FALSE
    )
  }
  existing <- lisa_msigdb_2026_1_existing_resource(
    registry_path, cache_root, contract
  )
  if (!is.null(existing)) {
    return(c(existing, list(
      source_zip_sha256 = contract$source_zip_sha256,
      source_url = contract$source_url,
      acquired = FALSE,
      transformation = "H,C2,C5 human RDS; distinct; order(gs_name,gene_symbol); readr::write_tsv"
    )))
  }
  acquired <- is.null(source_zip)
  if (acquired) {
    lisa_msigdb_2026_1_require_terms(accept_terms)
    source_zip <- lisa_msigdb_2026_1_download_zip(cache_root, contract)
  } else {
    source_zip <- lisa_assert_msigdb_2026_1_zip(source_zip, contract)
  }
  staging <- tempfile("lisar-msigdb-2026.1-")
  if (!dir.create(staging, mode = "0700", showWarnings = FALSE) ||
      lisa_path_is_link(staging)) {
    stop("LISA-RESOURCE-028 could not create a private TERM2GENE staging directory.", call. = FALSE)
  }
  lisa_set_private_mode(staging, directory = TRUE)
  on.exit(if (lisa_path_entry_exists(staging)) {
    try(lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL), silent = TRUE)
  }, add = TRUE)
  extracted <- tryCatch(
    utils::unzip(source_zip, files = contract$rds_files, exdir = staging),
    error = function(error) stop(
      "LISA-RESOURCE-028 bounded ZIP extraction failed: ", conditionMessage(error),
      call. = FALSE
    )
  )
  if (length(extracted) != length(contract$rds_files)) {
    stop("LISA-RESOURCE-028 bounded ZIP extraction did not return all required RDS members.", call. = FALSE)
  }
  term2gene <- lisa_msigdb_2026_1_term2gene_table(staging, contract)
  source_after_read <- lisa_sha256_file(source_zip)
  if (!identical(source_after_read, contract$source_zip_sha256)) {
    stop("LISA-RESOURCE-028 source_zip changed while it was being transformed. Repair: use an immutable local ZIP and retry.", call. = FALSE)
  }
  canonical_tsv <- file.path(staging, "msigdb_term2gene_2026.1.tsv")
  lisa_guarded_write(canonical_tsv, function(path) {
    readr::write_tsv(term2gene, path)
  }, run_root = staging, overwrite = FALSE)
  observed_sha256 <- lisa_sha256_file(canonical_tsv)
  observed_bytes <- unname(file.info(canonical_tsv)$size)
  if (!identical(observed_sha256, contract$output_sha256) ||
      !identical(as.numeric(observed_bytes), as.numeric(contract$output_bytes))) {
    stop(
      "LISA-RESOURCE-028 canonical MSigDB TERM2GENE digest or byte-size mismatch. Repair: use the exact reviewed ZIP and compatible canonical writer; do not register these bytes.",
      call. = FALSE
    )
  }
  installed <- install_lisa_resource(
    source = canonical_tsv,
    logical_id = "msigdb_term2gene",
    version = "2026.1",
    species = "Homo sapiens",
    modality = "all",
    schema = "term2gene@1",
    approved_origin = lisa_msigdb_2026_1_term2gene_origin(contract),
    expected_sha256 = contract$output_sha256,
    cache_root = cache_root,
    registry_path = registry_path
  )
  c(installed, list(
    source_zip_sha256 = contract$source_zip_sha256,
    source_url = contract$source_url,
    acquired = acquired,
    transformation = "H,C2,C5 human RDS; distinct; order(gs_name,gene_symbol); readr::write_tsv"
  ))
}
