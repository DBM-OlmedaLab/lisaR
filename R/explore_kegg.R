# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# Native KEGG pathway maps: resources and preparation.
#
# Two facts drive everything in this file.
#
# 1. The map's inputs are not in the run. A completed STANDARD run keeps its
#    DE tables and its category evidence, but a painted KEGG map also needs a
#    KGML diagram, a base image and a gene->pathway link table. Those live in a
#    cached snapshot, and an accepted run may well declare none: the run this
#    milestone is built against has `run_kegg_maps: false`. So the snapshot is
#    an EXPLICIT declaration made once against the exploration workspace, never
#    an inference, and it is always resolved `cache_only` -- a missing snapshot
#    is an error that names what is absent, never a download or an install.
#
# 2. Enumerating which pathways exist is expensive. Deciding that a category
#    can offer `hsa04110` means mapping every DE symbol to Entrez and joining
#    the snapshot's pathway links -- far too much to do on the navigation path
#    that H2_SCOPE requires to stay free of work. So enumeration is an explicit,
#    reusable PREPARATION step that writes a ranked index outside the source
#    run, and navigation only reads that table.
#
# The preparation is not a reanalysis. It runs the accepted painter's own
# summary code, in a mode that refuses to paint, and stores the ranking that
# painter would have used. No DE, GSEA or ORA is invoked, and nothing is written
# inside the source run.
# ---------------------------------------------------------------------------

lisa_explore_kegg_declaration_path <- function(ws) {
  file.path(ws$root, "prepared", "kegg-resources.json")
}

lisa_explore_kegg_index_path <- function(ws) {
  file.path(ws$root, "prepared", "kegg-pathway-index.tsv")
}

lisa_explore_kegg_preparation_path <- function(ws) {
  file.path(ws$root, "prepared", "kegg-pathway-preparation.json")
}

lisa_explore_kegg_empty_units <- function() {
  data.frame(analysis_id = character(), collection = character(),
             stringsAsFactors = FALSE)
}

# Normalise the JSON representation without relying on jsonlite simplifying a
# one-row array in the same way as a many-row array. The declaration is internal,
# but this defensive reader also makes a partially written or hand-edited file
# fail closed instead of turning it into invented preparation provenance.
lisa_explore_kegg_unit_table <- function(value) {
  empty <- lisa_explore_kegg_empty_units()
  if (is.null(value)) return(empty)
  if (is.data.frame(value)) {
    if (!all(c("analysis_id", "collection") %in% names(value))) return(empty)
    out <- value[, c("analysis_id", "collection"), drop = FALSE]
  } else if (is.list(value) &&
             all(c("analysis_id", "collection") %in% names(value))) {
    out <- data.frame(analysis_id = as.character(value$analysis_id),
                      collection = as.character(value$collection),
                      stringsAsFactors = FALSE)
  } else if (is.list(value)) {
    rows <- lapply(value, function(unit) {
      if (!is.list(unit) ||
          !all(c("analysis_id", "collection") %in% names(unit))) return(NULL)
      data.frame(analysis_id = as.character(unit$analysis_id)[[1L]],
                 collection = as.character(unit$collection)[[1L]],
                 stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(empty)
    out <- do.call(rbind, rows)
  } else {
    return(empty)
  }
  out$analysis_id <- lisa_explore_blank(out$analysis_id)
  out$collection <- lisa_explore_blank(out$collection)
  out <- unique(out[nzchar(out$analysis_id) & nzchar(out$collection), , drop = FALSE])
  out[order(out$analysis_id, out$collection), , drop = FALSE]
}

lisa_explore_kegg_preparation <- function(ws) {
  record <- lisa_explore_read_json(lisa_explore_kegg_preparation_path(ws))
  if (is.null(record) ||
      !identical(as.character(record$schema),
                 "lisa-explore-kegg-preparation/1") ||
      !is.list(record$context) || is.null(record$index_sha256)) return(NULL)
  record$prepared_units <- lisa_explore_kegg_unit_table(record$prepared_units)
  record
}

lisa_explore_kegg_declaration <- function(ws) {
  declaration <- lisa_explore_read_json(lisa_explore_kegg_declaration_path(ws))
  if (is.null(declaration)) return(NULL)
  needed <- c("cache_root", "snapshot_id", "species", "organism")
  if (!all(needed %in% names(declaration))) return(NULL)
  declaration
}

lisa_explore_kegg_index <- function(ws) {
  path <- lisa_explore_kegg_index_path(ws)
  if (!file.exists(path)) return(NULL)
  index <- read_lisa_tsv(path)
  if (!nrow(index)) return(NULL)
  declaration <- lisa_explore_kegg_declaration(ws)
  if (is.null(declaration) ||
      !lisa_explore_kegg_index_is_current(ws, index, declaration)) return(NULL)
  # New indexes have a separately persisted prepared-unit declaration. Bind it
  # to the exact TSV bytes so the two-file update fails closed during the brief
  # interval between their guarded atomic writes. A legacy index has no such
  # record and remains readable for its positive rows, but missing units are
  # never inferred from it (see lisa_explore_unavailable_reason()).
  if (file.exists(lisa_explore_kegg_preparation_path(ws))) {
    preparation <- lisa_explore_kegg_preparation_status(ws, declaration)
    if (!identical(preparation$state, "current")) return(NULL)
  }
  index
}

# A prepared index is reusable only for the source and resource declaration that
# produced it. Workspaces persist across sessions, so merely finding a TSV is not
# enough: a reconfigured snapshot or a changed source manifest must fail closed
# and require explicit preparation again rather than silently offering stale
# pathway/category associations under a new paint context.
lisa_explore_kegg_index_context <- function(ws, declaration) {
  list(
    source_manifest_hash = as.character(ws$source_manifest_hash),
    kegg_cache_root = as.character(declaration$cache_root),
    kegg_snapshot_id = as.character(declaration$snapshot_id),
    kegg_species = as.character(declaration$species),
    kegg_organism = as.character(declaration$organism)
  )
}

lisa_explore_kegg_index_is_current <- function(ws, index, declaration) {
  expected <- lisa_explore_kegg_index_context(ws, declaration)
  all(vapply(names(expected), function(field) {
    field %in% names(index) && nrow(index) > 0L &&
      all(lisa_explore_blank(index[[field]]) == expected[[field]])
  }, logical(1)))
}

lisa_explore_kegg_same_context <- function(left, right) {
  fields <- names(right)
  is.list(left) && all(fields %in% names(left)) &&
    all(vapply(fields, function(field) {
      identical(as.character(left[[field]]), as.character(right[[field]]))
    }, logical(1)))
}

# State, rather than a bare TRUE/FALSE, is important to user-facing reasons:
# absent is a legacy/unprepared declaration, stale names a different source or
# snapshot, and invalid includes a TSV changed independently of its provenance.
lisa_explore_kegg_preparation_status <- function(ws, declaration = NULL) {
  path <- lisa_explore_kegg_preparation_path(ws)
  if (!file.exists(path)) {
    return(list(state = "absent", units = lisa_explore_kegg_empty_units()))
  }
  record <- lisa_explore_kegg_preparation(ws)
  if (is.null(record)) {
    return(list(state = "invalid", units = lisa_explore_kegg_empty_units()))
  }
  if (is.null(declaration)) declaration <- lisa_explore_kegg_declaration(ws)
  if (is.null(declaration) ||
      !lisa_explore_kegg_same_context(
        record$context, lisa_explore_kegg_index_context(ws, declaration))) {
    return(list(state = "stale", units = lisa_explore_kegg_empty_units(),
                record = record))
  }
  index_path <- lisa_explore_kegg_index_path(ws)
  if (!file.exists(index_path) ||
      !identical(as.character(record$index_sha256),
                 as.character(lisa_sha256_file(index_path)))) {
    return(list(state = "invalid", units = lisa_explore_kegg_empty_units(),
                record = record))
  }
  list(state = "current", units = record$prepared_units, record = record)
}

lisa_explore_kegg_unit_selected <- function(units, selected) {
  if (!nrow(units) || !nrow(selected)) return(rep(FALSE, nrow(units)))
  paste(units$analysis_id, units$collection, sep = "\r") %in%
    paste(selected$analysis_id, selected$collection, sep = "\r")
}

lisa_explore_kegg_empty_index <- function(context = NULL) {
  result <- data.frame(
    analysis_id = character(), collection = character(),
    category_id = character(), kegg_id = character(),
    kegg_title = character(), rank = integer(),
    n_evidence_genes = integer(), resources_available = logical(),
    kegg_resource_digest = character(), stringsAsFactors = FALSE)
  if (!is.null(context)) {
    for (field in names(context)) result[[field]] <- character()
  }
  result
}

lisa_explore_kegg_preparation_record <- function(context, units, index_path) {
  unit_rows <- lapply(seq_len(nrow(units)), function(index) {
    list(analysis_id = as.character(units$analysis_id[[index]]),
         collection = as.character(units$collection[[index]]))
  })
  list(schema = "lisa-explore-kegg-preparation/1",
       context = context, prepared_units = unname(unit_rows),
       index_sha256 = lisa_sha256_file(index_path),
       updated_at = lisa_explore_now())
}

# The declaration the catalog and the renderer both bind. Kept in one function so
# the reuse key, the catalog row and the painter arguments cannot drift apart.
lisa_explore_kegg_settings <- function(declaration) {
  list(
    kegg_cache_root = as.character(declaration$cache_root),
    kegg_snapshot_id = as.character(declaration$snapshot_id),
    kegg_species = as.character(declaration$species),
    max_abs_log2fc = as.character(declaration$max_abs_log2fc %||% "0.5"),
    color_power = as.character(declaration$color_power %||% "1")
  )
}

lisa_explore_kegg_setting <- function(value, field) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0) {
    stop("LISA-EXPLORE-067 ", field,
         " must be one finite number greater than zero.", call. = FALSE)
  }
  format(value, scientific = FALSE, trim = TRUE, digits = 15)
}

# Return the validated content identity recorded by the immutable-snapshot
# contract. `lisa_read_kegg_snapshot()` verifies both the metadata and the bytes
# before the digest is accepted; callers never infer identity from a path or a
# snapshot label alone.
lisa_explore_kegg_resource_identity <- function(cache_root, snapshot_id,
                                                 organism, resource_type, key) {
  resource <- lisa_read_kegg_snapshot(cache_root, snapshot_id, organism,
                                      resource_type, key)
  metadata <- resource$metadata
  paste(resource_type, key, as.character(metadata$sha256[[1L]]), sep = "=")
}

lisa_explore_kegg_pathway_resource_digest <- function(cache_root, snapshot_id,
                                                       organism, kegg_id) {
  resources <- c(
    lisa_explore_kegg_resource_identity(
      cache_root, snapshot_id, organism, "pathway_list", "all"),
    lisa_explore_kegg_resource_identity(
      cache_root, snapshot_id, organism, "pathway_links", "all"),
    lisa_explore_kegg_resource_identity(
      cache_root, snapshot_id, organism, "kgml", kegg_id),
    lisa_explore_kegg_resource_identity(
      cache_root, snapshot_id, organism, "image", kegg_id)
  )
  lisa_sha256_text(paste(resources, collapse = "\n"))
}

#' Declare the cached KEGG snapshot native pathway maps may use
#'
#' Native KEGG maps are painted over a KGML diagram and a base image that live in
#' an immutable local snapshot. This records which snapshot an exploration
#' workspace may read. The snapshot must already exist: it is verified here and
#' always resolved in cache-only mode, so lisaR never downloads or installs one.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param cache_root Root of the existing local KEGG snapshot cache.
#' @param snapshot_id Identifier of the immutable snapshot inside `cache_root`.
#' @param species `"Homo sapiens"` or `"Mus musculus"`.
#' @param max_abs_log2fc Node-fill cap, the painter's own default when omitted.
#' @param color_power Node-fill colour power, the painter's own default when
#'   omitted.
#' @return Invisibly, the recorded declaration.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_configure_kegg(ws, cache_root, "snapshot-2026", "Homo sapiens")
#' }
lisa_explore_configure_kegg <- function(ws, cache_root, snapshot_id, species,
                                        max_abs_log2fc = 0.5, color_power = 1) {
  lisa_explore_assert_workspace(ws)
  cache_root <- lisa_explore_scalar(cache_root, "cache_root")
  snapshot_id <- lisa_safe_id(lisa_explore_scalar(snapshot_id, "snapshot_id"),
                              "snapshot_id")
  species <- lisa_explore_scalar(species, "species")
  if (!dir.exists(cache_root)) {
    stop("LISA-EXPLORE-060 the KEGG snapshot cache root does not exist: ",
         cache_root, ". Native pathway maps read an existing local snapshot; ",
         "lisaR does not create or download one.", call. = FALSE)
  }
  cache_root <- normalizePath(cache_root, winslash = "/", mustWork = TRUE)
  organism <- switch(tolower(trimws(species)),
    "homo sapiens" = , "human" = "hsa",
    "mus musculus" = , "mouse" = "mmu",
    stop("LISA-EXPLORE-061 species must be Homo sapiens or Mus musculus; ",
         "received: ", species, ".", call. = FALSE))
  species <- if (identical(organism, "hsa")) "Homo sapiens" else "Mus musculus"
  max_abs_log2fc <- lisa_explore_kegg_setting(max_abs_log2fc, "max_abs_log2fc")
  color_power <- lisa_explore_kegg_setting(color_power, "color_power")
  # Verify the snapshot is really there before recording it, so a typo is caught
  # now rather than by a worker halfway through a render.
  organism_root <- file.path(cache_root, "kegg", snapshot_id, organism)
  missing <- c("pathway_list", "pathway_links", "kgml", "image")
  missing <- missing[!dir.exists(file.path(organism_root, missing))]
  if (length(missing)) {
    stop("LISA-EXPLORE-062 the declared KEGG snapshot is incomplete at ",
         organism_root, "; absent resource directories: ",
         paste(missing, collapse = ", "),
         ". Supply a complete local snapshot; no resource is fetched.",
         call. = FALSE)
  }
  # The two collection-wide inputs are required by every preparation and every
  # painted map. Validate their metadata and bytes now, not just their parent
  # directories, so a corrupt declaration is never persisted as usable.
  invisible(lisa_explore_kegg_resource_identity(
    cache_root, snapshot_id, organism, "pathway_list", "all"))
  invisible(lisa_explore_kegg_resource_identity(
    cache_root, snapshot_id, organism, "pathway_links", "all"))
  declaration <- list(
    schema = "lisa-explore-kegg/2", cache_root = cache_root,
    snapshot_id = snapshot_id, species = species, organism = organism,
    max_abs_log2fc = max_abs_log2fc,
    color_power = color_power,
    access_mode = "cache_only", declared_at = lisa_explore_now())
  lisa_guarded_dir_create(file.path(ws$root, "prepared"), run_root = ws$root)
  lisa_explore_write_json(declaration, lisa_explore_kegg_declaration_path(ws),
                          run_root = ws$root)
  invisible(declaration)
}

# The exact-product declaration handed to the catalog and, through the worker, to
# `plan_lisa_extension()`. Reading it touches two prepared files and nothing else.
lisa_explore_exact_products <- function(ws) {
  products <- list(heatmap_variants = lisa_explore_heatmap_variants(),
                   native_single_de = TRUE,
                   contrast_products = lisa_explore_contrast_products(ws))
  declaration <- lisa_explore_kegg_declaration(ws)
  index <- if (is.null(declaration)) NULL else lisa_explore_kegg_index(ws)
  if (!is.null(index)) {
    settings <- lisa_explore_kegg_settings(declaration)
    # Only the pathways the declared snapshot can actually paint become catalog
    # rows. The rest stay explicitly flagged in the prepared index. An exact API
    # request can obtain the missing-resource reason; the current category UI
    # does not claim a separate notice for every unpaintable pathway when that
    # category also has paintable pathways.
    if ("resources_available" %in% names(index)) {
      index <- index[as.logical(index$resources_available) %in% TRUE, , drop = FALSE]
    }
    if (nrow(index)) {
      # One catalog row per (analysis, collection, associated category, pathway).
      # The pathway is the entity; the category is the navigation association,
      # and it is deliberately NOT part of the artifact's identity -- see
      # `lisa_explore_identity_fields()`.
      rows <- index[, c("analysis_id", "collection", "category_id", "kegg_id",
                        "kegg_title", "rank", "kegg_resource_digest"),
                    drop = FALSE]
      for (name in names(settings)) rows[[name]] <- settings[[name]]
      products$kegg_maps <- rows
    }
  }
  products
}

# --- preparation -----------------------------------------------------------

lisa_explore_kegg_prepare_root <- function(ws) {
  file.path(ws$root, "prepared", "kegg-work")
}

#' Prepare the native KEGG pathway index for an exploration workspace
#'
#' Computes, once, which KEGG pathways each collection's saved evidence genes map
#' to and how the accepted painter ranks them, so an interface can offer an exact
#' pathway selector without any generator running on the navigation path. It
#' reads saved evidence and the declared local snapshot, writes only inside the
#' exploration workspace, and paints nothing.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param analyses Optional analysis identifiers to restrict preparation to.
#' @param collections Optional collection identifiers to restrict preparation to.
#' @return Invisibly, the prepared index as a data frame.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_prepare_kegg_index(ws)
#' }
lisa_explore_prepare_kegg_index <- function(ws, analyses = NULL,
                                            collections = NULL) {
  lisa_explore_assert_workspace(ws)
  declaration <- lisa_explore_kegg_declaration(ws)
  if (is.null(declaration)) {
    stop("LISA-EXPLORE-063 no KEGG snapshot is declared for this workspace. ",
         "Call lisa_explore_configure_kegg() first.", call. = FALSE)
  }
  settings <- lisa_explore_kegg_settings(declaration)
  # Deliberately the UNEXPANDED catalog: preparation must not depend on the
  # index it is about to write.
  catalog <- lisa_extension_discover_catalog(ws$source_run)
  units <- unique(catalog[catalog$unit_type == "single_de" &
                          catalog$product %in% c("gene_cards", "volcano"),
                          c("analysis_id", "collection", "gsea_path", "de_path",
                            "summary_path"), drop = FALSE])
  if (!is.null(analyses)) {
    units <- units[units$analysis_id %in% as.character(analyses), , drop = FALSE]
  }
  if (!is.null(collections)) {
    units <- units[units$collection %in% as.character(collections), , drop = FALSE]
  }
  if (!nrow(units)) {
    stop("LISA-EXPLORE-064 no gene-level collection matched; there is nothing ",
         "to prepare.", call. = FALSE)
  }
  selected_units <- unique(units[, c("analysis_id", "collection"), drop = FALSE])
  selected_units <- selected_units[order(selected_units$analysis_id,
                                           selected_units$collection), , drop = FALSE]
  context <- lisa_explore_kegg_index_context(ws, declaration)
  package_dir <- lisa_resolve_package_dir()
  prepared <- list()
  for (index in seq_len(nrow(units))) {
    unit <- units[index, , drop = FALSE]
    prepared[[length(prepared) + 1L]] <- lisa_explore_prepare_one_kegg_index(
      ws, unit, settings, package_dir,
      as.character(declaration$organism))
  }
  prepared <- Filter(function(rows) !is.null(rows) && nrow(rows), prepared)
  result <- if (length(prepared)) do.call(rbind, prepared) else {
    lisa_explore_kegg_empty_index(context)
  }
  if (nrow(result)) {
    for (field in names(context)) result[[field]] <- context[[field]]
  }
  if (nrow(result)) {
    # Reported, not hidden: a reader who expected a map that the snapshot does
    # not carry should learn it here, once, rather than from a failed render.
    pathways <- unique(result[, c("kegg_id", "resources_available")])
    message(sprintf(
      "Prepared %d pathway(s) for %d category association(s); %d paintable from snapshot '%s', %d absent from it.",
      nrow(pathways), nrow(result), sum(pathways$resources_available),
      settings$kegg_snapshot_id, sum(!pathways$resources_available)))
  }
  # Enumeration happens outside the shared index lock. Promotion is short and
  # atomic per file: under the lock we re-check context, merge only provenance-
  # backed units from that exact context, replace every selected unit (including
  # a now-empty result), then bind the declaration to the final TSV digest.
  final <- lisa_explore_with_lock(ws, {
    current_declaration <- lisa_explore_kegg_declaration(ws)
    if (is.null(current_declaration) ||
        !lisa_explore_kegg_same_context(
          context, lisa_explore_kegg_index_context(ws, current_declaration))) {
      stop("LISA-EXPLORE-068 the source or KEGG snapshot declaration changed ",
           "during preparation; the prepared result was not promoted.",
           call. = FALSE)
    }

    previous_status <- lisa_explore_kegg_preparation_status(
      ws, current_declaration)
    previous_rows <- lisa_explore_kegg_empty_index(context)
    previous_units <- lisa_explore_kegg_empty_units()
    if (identical(previous_status$state, "current")) {
      previous_units <- previous_status$units
      readable <- lisa_explore_kegg_index(ws)
      if (!is.null(readable)) previous_rows <- readable
    } else if (file.exists(lisa_explore_kegg_index_path(ws))) {
      message("LISA-EXPLORE-069 the existing KEGG pathway index has no current ",
              "prepared-unit provenance; unverifiable legacy/stale rows are ",
              "not merged. Only the selected units will be recorded.")
    }

    keep_rows <- !lisa_explore_kegg_unit_selected(previous_rows, selected_units)
    kept_rows <- previous_rows[keep_rows, , drop = FALSE]
    combined <- rbind(kept_rows, result)
    if (nrow(combined)) {
      combined <- unique(combined)
      combined <- combined[order(combined$analysis_id, combined$collection,
                                 combined$rank, combined$kegg_id,
                                 combined$category_id), , drop = FALSE]
      rownames(combined) <- NULL
    }
    kept_units <- previous_units[
      !lisa_explore_kegg_unit_selected(previous_units, selected_units), ,
      drop = FALSE]
    combined_units <- unique(rbind(kept_units, selected_units))
    combined_units <- combined_units[order(combined_units$analysis_id,
                                           combined_units$collection), , drop = FALSE]
    rownames(combined_units) <- NULL

    lisa_guarded_dir_create(file.path(ws$root, "prepared"), run_root = ws$root)
    target <- lisa_explore_kegg_index_path(ws)
    old_run_root <- getOption("lisaR.run_root", NULL)
    options(lisaR.run_root = ws$root)
    on.exit(options(lisaR.run_root = old_run_root), add = TRUE)
    write_lisa_tsv(combined, target)
    lisa_explore_write_json(
      lisa_explore_kegg_preparation_record(context, combined_units, target),
      lisa_explore_kegg_preparation_path(ws), run_root = ws$root)
    combined
  })
  invisible(final)
}

# One (analysis, collection). The painter's summary is produced by the painter
# itself, in `--index-only true` mode, so the ranking offered to the reader is by
# construction the ranking a paint would use -- not a reimplementation of it that
# could drift.
lisa_explore_prepare_one_kegg_index <- function(ws, unit, settings, package_dir,
                                                declaration_organism) {
  analysis <- as.character(unit$analysis_id[[1L]])
  collection <- as.character(unit$collection[[1L]])
  work_root <- file.path(lisa_explore_kegg_prepare_root(ws),
                         paste0(analysis, "__", collection))
  if (lisa_path_entry_exists(work_root)) {
    lisa_guarded_delete(work_root, recursive = TRUE, run_root = ws$root)
  }
  lisa_guarded_dir_create(work_root, run_root = ws$root)
  work_root <- lisa_run_root(work_root)
  old_run_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = work_root)
  on.exit(options(lisaR.run_root = old_run_root), add = TRUE)

  stage <- c("config/de_index.tsv", "config/contrast_index.tsv")
  stage <- stage[file.exists(file.path(ws$source_run, stage))]
  stage <- c(stage, as.character(unit$gsea_path[[1L]]),
             as.character(unit$de_path[[1L]]),
             as.character(unit$summary_path[[1L]]))
  lisa_extension_stage_files(ws$source_run, work_root, unique(stage[nzchar(stage)]))
  prefix <- lisa_extension_file_label_prefix(as.character(unit$gsea_path[[1L]]),
                                             analysis)
  built <- lisa_build_single_gene_level_tables(analysis, collection,
    file.path(work_root, "outputs"), prefix)
  if (!identical(built$status, "completed")) {
    stop("LISA-EXPLORE-065 the gene-level evidence needed to enumerate pathways ",
         "could not be built for ", analysis, "/", collection, ": ",
         built$message, call. = FALSE)
  }
  emitted <- file.path(work_root, "pathway-index.tsv")
  code_ledger <- lisa_write_code_identity_ledger(
    package_dir, "build_single_de_kegg_pathway_painter.R",
    file.path(work_root, "code_identity.tsv"),
    executable_recipes = "reproduce_lisa_figure.R")
  result <- lisa_run_post_script(package_dir,
    "build_single_de_kegg_pathway_painter.R",
    c("--project-dir", work_root, "--analysis-id", analysis,
      "--universe", collection, "--species", settings$kegg_species,
      "--kegg-cache-root", settings$kegg_cache_root,
      "--kegg-snapshot-id", settings$kegg_snapshot_id,
      "--kegg-access-mode", "cache_only",
      "--max-abs-log2fc", settings$max_abs_log2fc,
      "--color-power", settings$color_power,
      "--emit-pathway-index", emitted,
      # The whole point: enumerate without painting.
      "--index-only", "true"),
    trusted_run_root = work_root, code_ledger = code_ledger)
  if (!identical(result$status, "completed")) {
    stop("LISA-EXPLORE-066 pathway enumeration failed for ", analysis, "/",
         collection, ": ", result$message, call. = FALSE)
  }
  if (!file.exists(emitted)) return(NULL)
  summary <- read_lisa_tsv(emitted)
  if (!nrow(summary)) return(NULL)
  # The painter records the categories a pathway's evidence genes came from as a
  # semicolon-joined field. Expanded here into one row per association, which is
  # what the catalog needs and what makes the map reachable from each category
  # sheet without duplicating the map itself.
  rows <- lapply(seq_len(nrow(summary)), function(index) {
    row <- summary[index, , drop = FALSE]
    categories <- unlist(strsplit(as.character(row$categories %||% ""), ";",
                                  fixed = TRUE), use.names = FALSE)
    categories <- trimws(categories)
    categories <- categories[nzchar(categories)]
    if (!length(categories)) return(NULL)
    data.frame(analysis_id = analysis, collection = collection,
               category_id = categories, kegg_id = as.character(row$kegg_id),
               kegg_title = as.character(row$kegg_title),
               rank = as.integer(row$rank),
               n_evidence_genes = as.integer(row$n_evidence_genes),
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  rows <- do.call(rbind, rows)

  organism <- as.character(declaration_organism)
  common <- vapply(c("pathway_list", "pathway_links"), function(resource_type) {
    lisa_explore_kegg_resource_identity(
      settings$kegg_cache_root, settings$kegg_snapshot_id, organism,
      resource_type, "all")
  }, character(1))
  pathway_ids <- unique(as.character(rows$kegg_id))
  digests <- vapply(pathway_ids, function(kegg_id) {
    tryCatch({
      specific <- vapply(c("kgml", "image"), function(resource_type) {
        lisa_explore_kegg_resource_identity(
          settings$kegg_cache_root, settings$kegg_snapshot_id, organism,
          resource_type, kegg_id)
      }, character(1))
      lisa_sha256_text(paste(c(common, specific), collapse = "\n"))
    }, error = function(error) "")
  }, character(1))
  rows$kegg_resource_digest <- unname(digests[match(rows$kegg_id, pathway_ids)])
  rows$resources_available <- nzchar(rows$kegg_resource_digest)
  rows
}
