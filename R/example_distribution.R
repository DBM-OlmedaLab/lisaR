# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# General acquisition route for the reviewed example distribution archives.
#
# The archives are prepared and reviewed outside this package; their identity
# (name, SHA-256, byte size, member contract) is pinned here.  The recipient
# never supplies a digest and is never asked to point at a previously reviewed
# directory: this entry obtains one fixed archive, verifies it, expands it into
# a private staging tree, re-verifies every scientific file against the bundle
# manifest that travels inside the archive, and only then promotes the bundle
# atomically.  The existing per-example preparers consume that bundle
# directory unchanged.

lisa_example_distribution_contracts <- function() {
  list(
    `riaz-gse91061` = list(
      example_id = "riaz-gse91061",
      version = "1.0.0",
      archive_name = "riaz-gse91061-lisa-examples-v1.0.0.tar.gz",
      archive_sha256 = "983efbe6fce2bfbbd0a3df691d60603c53b452d1e2cdc3d5ab08642df808b144",
      archive_bytes = 19320799,
      # Independent bounds; the pinned digest anchors identity, these keep an
      # unexpected archive from being expanded at all.
      max_archive_bytes = 64 * 1024^2,
      max_members = 64L,
      max_total_extracted_bytes = 256 * 1024^2,
      download_timeout_seconds = 900L,
      bundle_manifest = "prepared_bundle_manifest.tsv",
      required_files = c(
        "prepared_bundle_manifest.tsv", "prepared_bundle_provenance.tsv",
        "ATTRIBUTION_AND_LICENSE.md", "CHECKSUMS.sha256", "VERSION.txt"
      ),
      manifest_rows = 17L,
      manifest_prefixes = c(
        "outputs/design/", "outputs/models/", "outputs/lisa_inputs/"
      ),
      # Published versioned data asset; its exact bytes remain pinned above.
      source_url = "https://github.com/DBM-OlmedaLab/lisaR-example-data/releases/download/v1.0.0/riaz-gse91061-lisa-examples-v1.0.0.tar.gz",
      preparer_env = "LISAR_RIAZ_PREPARED_BUNDLE",
      preparer_script = "examples/riaz-gse91061/scripts/11_prepare_prepared_project.R",
      project_env = "LISAR_RIAZ_PROJECT_DIR",
      # This preparer stops at the promoted project and prints its own next
      # stage; it never calls run_lisa(), so it needs no prepare-only switch.
      prepare_only_env = NA_character_,
      citation = paste0(
        "Riaz N, Havel JJ, Makarov V, et al. Tumor and Microenvironment ",
        "Evolution during Immunotherapy with Nivolumab. Cell 2017;171:934-949. ",
        "GEO accession GSE91061."
      ),
      provenance = paste0(
        "riaz-prepared-bundle@2; DESeq2 results and paired design exported by ",
        "the reviewed canonical Riaz project; redistributed derivatives only, ",
        "no raw sequencing data; see ATTRIBUTION_AND_LICENSE.md inside the ",
        "archive. Distribution of these bytes is not itself a licence grant."
      )
    ),
    `cptac-ccrcc` = list(
      example_id = "cptac-ccrcc",
      version = "1.0.0",
      archive_name = "cptac-ccrcc-lisa-prepared-v1.0.0.tar.gz",
      archive_sha256 = "70b5715d9cfeb5e29cda65828b1c9ca1bba7b2d53405f2452d79b9c44fc98857",
      archive_bytes = 371383,
      max_archive_bytes = 16 * 1024^2,
      max_members = 32L,
      max_total_extracted_bytes = 64 * 1024^2,
      download_timeout_seconds = 900L,
      bundle_manifest = "cptac_ccrcc_bundle_manifest.tsv",
      required_files = c(
        "cptac_ccrcc_bundle_manifest.tsv", "CHECKSUMS.sha256", "README.md"
      ),
      manifest_rows = 6L,
      manifest_prefixes = c("de/", "audit/"),
      source_url = "https://github.com/DBM-OlmedaLab/lisaR-example-data/releases/download/v1.0.0/cptac-ccrcc-lisa-prepared-v1.0.0.tar.gz",
      preparer_env = "LISAR_CPTAC_CCRCC_BUNDLE",
      preparer_script = "examples/cptac-ccrcc/prepare_and_run_cptac_ccrcc.R",
      project_env = "LISAR_CPTAC_CCRCC_PROJECT_DIR",
      # This preparer continues into run_lisa() unless it is told otherwise, so
      # the command this entry prints must pin its prepare-only route: obtaining
      # a bundle never starts a scientific run on the recipient's behalf.
      prepare_only_env = "LISAR_CPTAC_CCRCC_PREPARE_ONLY",
      citation = paste0(
        "Clark DJ, Dhanasekaran SM, Petralia F, et al. Integrated Proteogenomic ",
        "Characterization of Clear Cell Renal Cell Carcinoma. Cell 2019;179:964-983. ",
        "CPTAC ccRCC (PDC000127)."
      ),
      provenance = paste0(
        "cptac-ccrcc-prepared-bundle@1; accepted paired differential ",
        "expression over 6482 genes plus five design/QC/provenance tables; ",
        "derived tables only, not the original abundance matrix. The source ",
        "protein-level table is not asserted to be byte-identical to any PDC ",
        "supplementary file, and no particular data licence is claimed for it."
      )
    )
  )
}

lisa_example_distribution_contract <- function(example) {
  contracts <- lisa_example_distribution_contracts()
  if (!is.character(example) || length(example) != 1L || is.na(example) ||
      !example %in% names(contracts)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-030 unknown example '%s'. ",
        "Repair: use one of: %s."
      ),
      paste(as.character(example), collapse = ", "),
      paste(names(contracts), collapse = ", ")
    ), call. = FALSE)
  }
  contracts[[example]]
}

# Resolve the one source this call may read.  A published pinned URL is used
# automatically; otherwise the caller must name an explicit reviewed source.
# The absence of a published endpoint is reported as a gate, never worked
# around with a guessed release address.
lisa_example_distribution_source <- function(source, contract) {
  if (is.null(source)) {
    if (!is.na(contract$source_url) && nzchar(contract$source_url)) {
      return(list(url = contract$source_url, pinned = TRUE))
    }
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-031 no public distribution endpoint is published yet for '%s'. ",
        "The reviewed archive '%s' (SHA-256 %s) exists but has not been released, ",
        "so there is no default download address and lisaR will not invent one. ",
        "Repair: wait for the reviewed publication, or pass source= naming one ",
        "explicit local file or file://, http:// or https:// address that serves ",
        "exactly that reviewed archive."
      ),
      contract$example_id, contract$archive_name, contract$archive_sha256
    ), call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(trimws(source)) || !identical(source, trimws(source)) ||
      grepl("[[:cntrl:]]", source)) {
    stop(
      "LISA-EXAMPLE-031 source must be one non-empty address or local archive path. Repair: pass a single file://, http:// or https:// URL, or one existing local archive file.",
      call. = FALSE
    )
  }
  if (grepl("^(https|http|file)://", source)) {
    return(list(url = source, pinned = FALSE))
  }
  if (!file.exists(source) || dir.exists(source) ||
      !isTRUE(utils::file_test("-f", source)) || lisa_path_is_link(source)) {
    stop(
      "LISA-EXAMPLE-031 local source must be one existing regular non-symlink archive file. Repair: supply the reviewed archive or an explicit file://, http:// or https:// address.",
      call. = FALSE
    )
  }
  list(
    url = paste0(
      "file://", normalizePath(source, winslash = "/", mustWork = TRUE)
    ),
    pinned = FALSE
  )
}

lisa_example_distribution_download <- function(url, destination, contract) {
  old_timeout <- getOption("timeout")
  options(timeout = contract$download_timeout_seconds)
  on.exit(options(timeout = old_timeout), add = TRUE)
  status <- tryCatch(
    utils::download.file(
      url, destination, method = "libcurl", mode = "wb", quiet = TRUE
    ),
    error = function(error) stop(sprintf(
      paste0(
        "LISA-EXAMPLE-032 example distribution download failed for '%s': %s. ",
        "Repair: check connectivity and the address, then retry; the ",
        "destination project and any existing bundle were not changed."
      ), contract$example_id, conditionMessage(error)
    ), call. = FALSE)
  )
  if (!identical(as.integer(status), 0L)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-032 example distribution download returned a non-zero status for '%s'. ",
        "Repair: check connectivity and the address, then retry; nothing was changed."
      ), contract$example_id
    ), call. = FALSE)
  }
  invisible(destination)
}

lisa_example_distribution_assert_archive <- function(archive, contract) {
  if (!file.exists(archive) || dir.exists(archive) ||
      !isTRUE(utils::file_test("-f", archive)) || lisa_path_is_link(archive)) {
    stop(
      "LISA-EXAMPLE-032 obtained archive is not one regular non-symlink file. Repair: retry the reviewed source.",
      call. = FALSE
    )
  }
  bytes <- unname(file.info(archive)$size)
  if (!is.finite(bytes) || bytes <= 0 || bytes > contract$max_archive_bytes) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-032 obtained archive for '%s' is empty or exceeds the bounded ",
        "acquisition policy. Repair: retry the reviewed source; no bytes were installed."
      ), contract$example_id
    ), call. = FALSE)
  }
  if (!identical(as.numeric(bytes), as.numeric(contract$archive_bytes))) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-033 obtained archive for '%s' has %.0f bytes but the reviewed ",
        "archive has %.0f. Repair: obtain the exact reviewed archive; nothing was installed."
      ), contract$example_id, as.numeric(bytes),
      as.numeric(contract$archive_bytes)
    ), call. = FALSE)
  }
  observed <- lisa_sha256_file(archive)
  if (!identical(observed, contract$archive_sha256)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-033 obtained archive SHA-256 does not match the reviewed archive for '%s'. ",
        "Repair: obtain the exact reviewed archive; a partial, mirrored or re-packaged ",
        "copy is refused and nothing was installed."
      ), contract$example_id
    ), call. = FALSE)
  }
  normalizePath(archive, winslash = "/", mustWork = TRUE)
}

lisa_example_distribution_archive_path <- function(cache_root, contract) {
  short <- switch(contract$example_id, `riaz-gse91061` = "riaz", `cptac-ccrcc` = "cptac")
  file.path(cache_root, "example-archives", paste0(short, "-", contract$version, ".tar.gz"))
}

# Obtain the archive into the managed cache.  A cached archive that already
# matches the pinned identity is reused without touching the network; a new
# download is written through the guarded temporary-plus-rename path and is
# verified before it can be reused.
lisa_example_distribution_obtain <- function(cache_root, contract, source) {
  archive <- lisa_example_distribution_archive_path(cache_root, contract)
  if (lisa_path_entry_exists(archive)) {
    verified <- tryCatch(
      lisa_example_distribution_assert_archive(archive, contract),
      error = identity
    )
    if (!inherits(verified, "error")) {
      return(list(path = verified, acquired = FALSE, url = NA_character_))
    }
    if (lisa_path_is_link(archive) || dir.exists(archive)) {
      stop(
        "LISA-EXAMPLE-032 cached example archive is not one regular file. Repair: remove the invalid cache entry after review, then retry.",
        call. = FALSE
      )
    }
  }
  resolved <- lisa_example_distribution_source(source, contract)
  lisa_guarded_dir_create(dirname(archive), cache_root)
  lisa_guarded_write(archive, function(path) {
    lisa_example_distribution_download(resolved$url, path, contract)
    lisa_example_distribution_assert_archive(path, contract)
  }, run_root = cache_root, overwrite = TRUE)
  list(
    path = lisa_example_distribution_assert_archive(archive, contract),
    acquired = TRUE,
    url = resolved$url
  )
}

# Reject an unsafe or defectively packaged archive from its listing, before a
# single member is written to disk.  Duplicate names are refused explicitly:
# a hard-linked or auto-recursed tar expands the same path more than once and
# defeats a post-extraction tree check.
lisa_example_distribution_members <- function(archive, contract) {
  members <- tryCatch(
    utils::untar(archive, list = TRUE, tar = "internal"),
    error = function(error) stop(sprintf(
      paste0(
        "LISA-EXAMPLE-034 example archive for '%s' could not be listed safely: %s. ",
        "Repair: obtain the exact reviewed archive."
      ), contract$example_id, conditionMessage(error)
    ), call. = FALSE)
  )
  if (!is.character(members) || !length(members) || anyNA(members)) {
    stop(sprintf(
      "LISA-EXAMPLE-034 example archive for '%s' has an unreadable member listing. Repair: obtain the exact reviewed archive.",
      contract$example_id
    ), call. = FALSE)
  }
  if (length(members) > contract$max_members) {
    stop(sprintf(
      "LISA-EXAMPLE-034 example archive for '%s' exceeds the bounded member policy. Repair: obtain the exact reviewed archive.",
      contract$example_id
    ), call. = FALSE)
  }
  normalised <- sub("/+$", "", members)
  if (anyDuplicated(normalised)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-034 example archive for '%s' lists the same path more than once. ",
        "Repair: obtain a deterministically packaged archive without duplicate or ",
        "hard-linked entries; nothing was extracted."
      ), contract$example_id
    ), call. = FALSE)
  }
  unsafe <- !nzchar(normalised) |
    grepl("^[/~]", normalised) |
    grepl("^[A-Za-z]:", normalised) |
    grepl("(^|/)\\.\\.(/|$)", normalised) |
    grepl("(^|/)\\.(/|$)", normalised) |
    grepl("[[:cntrl:]]", normalised) |
    grepl("\\\\", normalised)
  if (any(unsafe)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-034 example archive for '%s' contains an absolute, parent-relative ",
        "or otherwise unsafe member path. Repair: obtain the exact reviewed archive; ",
        "nothing was extracted."
      ), contract$example_id
    ), call. = FALSE)
  }
  files <- members[!grepl("/$", members)]
  if (!length(files)) {
    stop(sprintf(
      "LISA-EXAMPLE-034 example archive for '%s' contains no regular member. Repair: obtain the exact reviewed archive.",
      contract$example_id
    ), call. = FALSE)
  }
  missing <- setdiff(contract$required_files, files)
  if (length(missing)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-034 example archive for '%s' is missing required bundle-root ",
        "file(s): %s. Repair: obtain the exact reviewed archive; a bundle without its ",
        "manifest and provenance cannot be prepared."
      ), contract$example_id, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  list(all = members, files = files)
}

lisa_example_distribution_read_manifest <- function(path, contract) {
  manifest <- tryCatch(
    utils::read.delim(
      path, sep = "\t", header = TRUE, quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
      na.strings = character()
    ),
    error = function(error) stop(sprintf(
      "LISA-EXAMPLE-035 bundle manifest of '%s' could not be read: %s. Repair: obtain the exact reviewed archive.",
      contract$example_id, conditionMessage(error)
    ), call. = FALSE)
  )
  required <- c("path", "bytes", "sha256", "role")
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      !identical(nrow(manifest), as.integer(contract$manifest_rows))) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-035 bundle manifest of '%s' is malformed or does not list the ",
        "reviewed %d scientific file(s). Repair: obtain the exact reviewed archive."
      ), contract$example_id, as.integer(contract$manifest_rows)
    ), call. = FALSE)
  }
  manifest <- manifest[, required, drop = FALSE]
  bytes <- suppressWarnings(as.numeric(manifest$bytes))
  safe_path <- !grepl("^[/~]", manifest$path) &
    !grepl("(^|/)\\.\\.?(/|$)", manifest$path) &
    !grepl("[[:cntrl:]\\\\]", manifest$path) &
    nzchar(manifest$path)
  if (anyDuplicated(manifest$path) || !all(safe_path) ||
      !lisa_sha256_all_valid(manifest$sha256, nrow(manifest)) ||
      any(!is.finite(bytes)) || any(bytes <= 0)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-035 bundle manifest of '%s' has duplicate, unsafe or invalid rows. ",
        "Repair: obtain the exact reviewed archive."
      ), contract$example_id
    ), call. = FALSE)
  }
  if (!all(vapply(manifest$path, function(entry) {
    any(startsWith(entry, contract$manifest_prefixes))
  }, logical(1)))) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-035 bundle manifest of '%s' lists a file outside its reviewed ",
        "scientific directories. Repair: obtain the exact reviewed archive."
      ), contract$example_id
    ), call. = FALSE)
  }
  manifest$bytes <- bytes
  manifest
}

# Verify the expanded tree against the manifest that travelled inside the
# archive.  Every listed scientific file must exist as a regular non-symlink
# file with the manifest's exact size and digest, and the tree may not contain
# a file the archive listing did not declare.
lisa_example_distribution_verify_tree <- function(staging, members, contract) {
  lisa_assert_run_tree_safe(staging)
  tree <- lisa_scan_run_tree(staging)
  present <- substring(tree$path[!tree$isdir], nchar(staging) + 2L)
  unexpected <- setdiff(present, members$files)
  if (length(unexpected)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-036 expanded bundle of '%s' contains file(s) the archive listing ",
        "did not declare: %s. Repair: obtain the exact reviewed archive; nothing was installed."
      ), contract$example_id, paste(utils::head(unexpected, 5L), collapse = ", ")
    ), call. = FALSE)
  }
  missing <- setdiff(members$files, present)
  if (length(missing)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-036 bounded extraction of '%s' did not produce declared member(s): %s. ",
        "Repair: obtain the exact reviewed archive; nothing was installed."
      ), contract$example_id, paste(utils::head(missing, 5L), collapse = ", ")
    ), call. = FALSE)
  }
  manifest_path <- file.path(staging, contract$bundle_manifest)
  manifest <- lisa_example_distribution_read_manifest(manifest_path, contract)
  scientific <- file.path(staging, manifest$path)
  regular <- vapply(scientific, function(path) {
    file.exists(path) && !dir.exists(path) &&
      isTRUE(utils::file_test("-f", path)) && !lisa_path_is_link(path)
  }, logical(1))
  if (!all(regular)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-036 expanded bundle of '%s' has a missing or non-regular scientific ",
        "file. Repair: obtain the exact reviewed archive; nothing was installed."
      ), contract$example_id
    ), call. = FALSE)
  }
  observed_bytes <- as.numeric(file.info(scientific)$size)
  observed_sha256 <- unname(vapply(scientific, lisa_sha256_file, character(1)))
  if (!identical(observed_bytes, as.numeric(manifest$bytes)) ||
      !identical(observed_sha256, as.character(manifest$sha256))) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-036 expanded bundle of '%s' does not match its own reviewed manifest. ",
        "Repair: obtain the exact reviewed archive; nothing was installed."
      ), contract$example_id
    ), call. = FALSE)
  }
  total <- sum(as.numeric(file.info(file.path(staging, present))$size),
               na.rm = TRUE)
  if (!is.finite(total) || total > contract$max_total_extracted_bytes) {
    stop(sprintf(
      "LISA-EXAMPLE-036 expanded bundle of '%s' exceeds the bounded extraction policy. Repair: obtain the exact reviewed archive.",
      contract$example_id
    ), call. = FALSE)
  }
  list(
    manifest = manifest,
    manifest_sha256 = lisa_sha256_file(manifest_path),
    files = length(present),
    bytes = total
  )
}

# Build the command that consumes the promoted bundle.  It must be runnable
# exactly as printed, so every value is a concrete quoted path: an unquoted
# `<placeholder>` would be read by the shell as a redirection and would leave
# the project variable empty.  The project directory is a new sibling of the
# bundle; the preparers refuse to overwrite an existing one, so naming it here
# cannot destroy anything, and the recipient can edit it before running.
lisa_example_distribution_next_command <- function(contract, bundle) {
  template <- system.file(contract$preparer_script, package = "lisaR")
  if (!nzchar(template)) template <- contract$preparer_script
  project <- paste0(bundle, "-project")
  assignments <- c(
    sprintf("%s=%s", contract$preparer_env, shQuote(bundle)),
    sprintf("%s=%s", contract$project_env, shQuote(project))
  )
  # A preparer that would otherwise continue into run_lisa() is pinned to its
  # prepare-only route: this entry promises no scientific computation.
  if (!is.na(contract$prepare_only_env)) {
    assignments <- c(assignments, sprintf("%s=1", contract$prepare_only_env))
  }
  paste(c(assignments, "Rscript", "--vanilla", shQuote(template)),
        collapse = " ")
}

#' Install one reviewed lisaR example input bundle from its distribution archive
#'
#' Obtains the single reviewed distribution archive for one installed example,
#' verifies it against lisaR's pinned identity, expands it inside a private
#' staging tree, re-verifies every scientific file against the bundle manifest
#' carried in the archive, and promotes the verified bundle to `destination`
#' with one atomic rename. The recipient never supplies a digest and is never
#' asked to point at an already-reviewed directory: the archive identity, the
#' member contract and the manifest contract are internal.
#'
#' A network or integrity failure at any point leaves `destination` and any
#' previously installed bundle untouched; nothing is promoted until the whole
#' expanded tree has been verified. An archive already present in the managed
#' cache with the pinned identity is reused without network access.
#'
#' The resulting directory is the input of the example's existing preparer, so
#' this function performs no scientific computation: it neither refits models
#' nor runs LISA. The returned `next_command` is the preparer call, runnable
#' exactly as printed: it carries no placeholder, it names a concrete new
#' project directory beside the bundle, and where a preparer would otherwise
#' continue into `run_lisa()` it pins that preparer's prepare-only route. Edit
#' the project directory before running it if another location is wanted;
#' starting the scientific run afterwards stays an explicit separate decision.
#'
#' By default the reviewed archive is downloaded from its versioned public
#' release in `DBM-OlmedaLab/lisaR-example-data`, without credentials. Its checksum
#' is verified before use. Pass `source` to use a mirror or local archive
#' containing exactly the same bytes; verified cached archives are reused offline.
#'
#' @param example Installed example identifier, `"riaz-gse91061"` or
#'   `"cptac-ccrcc"`.
#' @param destination New directory to create for the verified bundle. It must
#'   not already exist; lisaR never overwrites it.
#' @param source Optional single explicit source for the reviewed archive: one
#'   `file://`, `http://` or `https://` address, or one existing local archive
#'   file. It selects where the pinned bytes are read from; it can never change
#'   which bytes are accepted.
#' @param cache_root Managed cache root used for the verified archive. The
#'   default is shared with resource resolution, so a second example install
#'   reuses an already-verified archive.
#'
#' @return Invisibly, a list with the example identifier and version, the
#'   promoted `bundle` directory, the cached `archive` and its SHA-256, whether
#'   the archive was `acquired` in this call, the verified `files` count and
#'   `bundle_manifest_sha256`, the `public_endpoint` status, the citation and
#'   provenance, and the `next_command` for the example's preparer.
#' @export
install_lisa_example_bundle <- function(
  example,
  destination,
  source = NULL,
  cache_root = lisa_dictionary_cache_root()
) {
  contract <- lisa_example_distribution_contract(example)
  target <- lisa_managed_destination(destination, create_parent = TRUE)
  if (lisa_path_entry_exists(target$path)) {
    stop(sprintf(
      paste0(
        "LISA-EXAMPLE-030 destination already exists and lisaR will not overwrite it: %s. ",
        "Repair: name a new directory, or use the existing verified bundle."
      ), target$path
    ), call. = FALSE)
  }
  cache_root <- lisa_resource_cache_root(
    lisa_path_walk_lexically(cache_root)$path
  )
  lisa_guarded_dir_create(cache_root, run_root = NULL)

  archive <- lisa_example_distribution_obtain(cache_root, contract, source)
  members <- lisa_example_distribution_members(archive$path, contract)

  staging <- file.path(
    target$parent,
    paste0(".", basename(target$path), ".bundle-staging-", Sys.getpid())
  )
  if (lisa_path_entry_exists(staging)) {
    stop(sprintf(
      "LISA-EXAMPLE-036 bundle staging path already exists: %s. Repair: inspect it and retry with a new destination.",
      staging
    ), call. = FALSE)
  }
  if (!dir.create(staging, mode = "0700", showWarnings = FALSE) ||
      lisa_path_is_link(staging)) {
    stop(
      "LISA-EXAMPLE-036 could not create a private bundle staging directory.",
      call. = FALSE
    )
  }
  lisa_set_private_mode(staging, directory = TRUE)
  staging <- normalizePath(staging, winslash = "/", mustWork = TRUE)
  promoted <- FALSE
  on.exit(if (!promoted && lisa_path_entry_exists(staging)) {
    try(lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL),
        silent = TRUE)
  }, add = TRUE)

  extraction <- tryCatch(
    utils::untar(
      archive$path, files = members$files, exdir = staging, tar = "internal"
    ),
    error = function(error) stop(sprintf(
      "LISA-EXAMPLE-036 bounded extraction of '%s' failed: %s. Repair: obtain the exact reviewed archive.",
      contract$example_id, conditionMessage(error)
    ), call. = FALSE)
  )
  if (!identical(as.integer(extraction), 0L)) {
    stop(sprintf(
      "LISA-EXAMPLE-036 bounded extraction of '%s' returned a non-zero status. Repair: obtain the exact reviewed archive; nothing was installed.",
      contract$example_id
    ), call. = FALSE)
  }
  verified <- lisa_example_distribution_verify_tree(staging, members, contract)
  # The verified archive must still be the pinned one after the whole read.
  if (!identical(lisa_sha256_file(archive$path), contract$archive_sha256)) {
    stop(sprintf(
      "LISA-EXAMPLE-033 example archive of '%s' changed while it was being expanded. Repair: retry with an immutable source; nothing was installed.",
      contract$example_id
    ), call. = FALSE)
  }

  lisa_promote_managed_directory(staging, target$path)
  promoted <- TRUE
  bundle <- normalizePath(target$path, winslash = "/", mustWork = TRUE)

  result <- list(
    example = contract$example_id,
    version = contract$version,
    bundle = bundle,
    archive = archive$path,
    archive_sha256 = contract$archive_sha256,
    archive_bytes = as.numeric(contract$archive_bytes),
    acquired = archive$acquired,
    source_url = archive$url,
    files = verified$files,
    scientific_files = nrow(verified$manifest),
    bundle_manifest = contract$bundle_manifest,
    bundle_manifest_sha256 = verified$manifest_sha256,
    public_endpoint = if (is.na(contract$source_url)) {
      "pending_reviewed_stable_distribution"
    } else {
      contract$source_url
    },
    citation = contract$citation,
    provenance = contract$provenance,
    next_command = lisa_example_distribution_next_command(contract, bundle)
  )
  invisible(result)
}
