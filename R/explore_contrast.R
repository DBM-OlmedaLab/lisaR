# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


lisa_explore_contrast_prepared_root <- function(ws) {
  file.path(ws$root, "prepared", "contrast")
}

# The `<context>` path segment. Binding the directory name to the source
# manifest makes "never silently mix prepared rows from different snapshots" a
# property of the filesystem layout rather than only of the receipt check.
lisa_explore_contrast_context_id <- function(ws) {
  substr(as.character(ws$source_manifest_hash), 1L, 16L)
}

lisa_explore_contrast_bundle_dir <- function(ws, contrast_id, collection) {
  file.path(lisa_explore_contrast_prepared_root(ws),
            lisa_explore_contrast_context_id(ws),
            lisa_safe_id(contrast_id, "contrast_id"),
            lisa_safe_id(collection, "collection"))
}

lisa_explore_contrast_receipt_path <- function(bundle_dir) {
  file.path(bundle_dir, "receipt.json")
}

lisa_explore_contrast_bundle_files <- function() {
  c(paired = "paired_gene_evidence.tsv",
    summary = "category_gene_summary.tsv")
}

lisa_explore_contrast_pathway_index_file <- function() "contrast_pathway_index.tsv"

# The effective fixed settings of the preparation and of every product drawn
# from it. They are the accepted builders' own defaults, named here so they
# enter the receipt -- and therefore the preparation digest and the artifact
# identity -- instead of being implicit. This is not a new control surface:
# there is no way to ask for a different value.
lisa_explore_contrast_settings <- function() {
  list(top_categories = "0", top_genes_per_category = "28", lfc_cap = "1.5",
       network_variant = "all_categories_top7",
       network_top_genes_per_category = "7",
       heatmap_top_categories = "6")
}

# The preparation's own code identity. A change to any script that produced the
# bundle invalidates it, which is what makes "fail closed and require
# re-preparation if any bound code changes" true rather than aspirational.
lisa_explore_contrast_preparation_scripts <- function() {
  c("build_contrast_gene_level_product.R",
    "build_contrast_kegg_pathway_painter.R",
    "build_contrast_gene_category_network.R")
}

lisa_explore_contrast_code_digest <- function(package_dir = NULL) {
  if (is.null(package_dir)) package_dir <- lisa_resolve_package_dir()
  rows <- lisa_code_identity_rows(package_dir,
                                  lisa_explore_contrast_preparation_scripts())
  text <- if (!nrow(rows)) "" else paste(
    apply(rows[order(rows$relative_path),
               c("relative_path", "sha256"), drop = FALSE], 1L, paste,
          collapse = ":"), collapse = "\n")
  lisa_sha256_text(text)
}

lisa_explore_contrast_context <- function(ws, declaration = NULL) {
  context <- list(source_manifest_hash = as.character(ws$source_manifest_hash),
                  kegg_cache_root = "", kegg_snapshot_id = "",
                  kegg_species = "", kegg_organism = "")
  if (!is.null(declaration)) {
    context$kegg_cache_root <- as.character(declaration$cache_root)
    context$kegg_snapshot_id <- as.character(declaration$snapshot_id)
    context$kegg_species <- as.character(declaration$species)
    context$kegg_organism <- as.character(declaration$organism)
  }
  context
}

# One digest over everything the bundle is bound to. It is what
# `lisa_explore_resource_digest()` folds into the artifact key, so a figure drawn
# from one preparation can never be served for another.
lisa_explore_contrast_preparation_digest <- function(receipt) {
  flatten <- function(rows) {
    if (is.null(rows) || !length(rows)) return("")
    paste(vapply(rows, function(row)
      paste(as.character(row$path), as.character(row$sha256), sep = "="),
      character(1)), collapse = "\n")
  }
  named <- function(values) {
    if (is.null(values) || !length(values)) return("")
    keys <- sort(names(values))
    paste(paste(keys, vapply(keys, function(key) as.character(values[[key]]),
                             character(1)), sep = "="), collapse = "\n")
  }
  lisa_sha256_text(paste(c(
    paste0("schema=", as.character(receipt$schema)),
    paste0("contrast_id=", as.character(receipt$contrast_id)),
    paste0("collection=", as.character(receipt$collection)),
    named(receipt$context), named(receipt$settings),
    paste0("code=", as.character(receipt$code_sha256)),
    paste0("kegg_resource_digest=", as.character(receipt$kegg_resource_digest)),
    flatten(receipt$inputs), flatten(receipt$files)), collapse = "\n"))
}

lisa_explore_contrast_read_receipt <- function(bundle_dir) {
  receipt <- lisa_explore_read_json(lisa_explore_contrast_receipt_path(bundle_dir))
  if (is.null(receipt) ||
      !identical(as.character(receipt$schema),
                 "lisa-explore-contrast-preparation/1") ||
      !is.list(receipt$context) || !is.list(receipt$settings)) return(NULL)
  receipt
}

# --- fail-closed bundle verification ---------------------------------------
#
# A bundle is usable only when every one of these still holds. Each branch names
# the concrete thing that changed, because "not applicable" is not an answer a
# reader can act on.
lisa_explore_contrast_bundle_status <- function(ws, contrast_id, collection,
                                                 declaration = NULL,
                                                 code_digest = NULL) {
  bundle_dir <- lisa_explore_contrast_bundle_dir(ws, contrast_id, collection)
  absent <- list(state = "absent", dir = bundle_dir, receipt = NULL)
  if (!file.exists(lisa_explore_contrast_receipt_path(bundle_dir))) return(absent)
  receipt <- lisa_explore_contrast_read_receipt(bundle_dir)
  if (is.null(receipt)) {
    return(list(state = "invalid", dir = bundle_dir, receipt = NULL,
                reason = "the preparation receipt is unreadable"))
  }
  if (is.null(declaration)) declaration <- lisa_explore_kegg_declaration(ws)
  expected <- lisa_explore_contrast_context(ws, declaration)
  if (!lisa_explore_kegg_same_context(receipt$context, expected)) {
    return(list(state = "stale", dir = bundle_dir, receipt = receipt,
                reason = paste0("the recorded preparation belongs to a ",
                                "different source run or KEGG snapshot context")))
  }
  if (is.null(code_digest)) code_digest <- lisa_explore_contrast_code_digest()
  if (!identical(as.character(receipt$code_sha256), as.character(code_digest))) {
    return(list(state = "stale", dir = bundle_dir, receipt = receipt,
                reason = paste0("the preparation code changed after this bundle ",
                                "was written")))
  }
  # Prepared files are verified against their OWN receipt, inside their own
  # contained workspace path. They are deliberately not looked for in the run
  # manifest: they are not source data.
  for (file in receipt$files) {
    relative <- lisa_explore_presentation_path(as.character(file$path),
                                               "prepared file")
    path <- file.path(bundle_dir, relative)
    if (!file.exists(path) ||
        !identical(lisa_sha256_file(path), as.character(file$sha256))) {
      return(list(state = "invalid", dir = bundle_dir, receipt = receipt,
                  reason = paste0("a prepared file is missing or changed on ",
                                  "disk: ", relative)))
    }
  }
  # Source inputs are verified the other way round: against the run manifest,
  # exactly as every other scientific input of this package is.
  manifest <- read_lisa_tsv(file.path(ws$source_run, "run_manifest.tsv"))
  for (input in receipt$inputs) {
    relative <- as.character(input$path)
    matched <- match(relative, manifest$path)
    if (is.na(matched) ||
        !identical(as.character(manifest$sha256[[matched]]),
                   as.character(input$sha256))) {
      return(list(state = "stale", dir = bundle_dir, receipt = receipt,
                  reason = paste0("a bound scientific input changed in the ",
                                  "source run: ", relative)))
    }
  }
  list(state = "current", dir = bundle_dir, receipt = receipt,
       digest = lisa_explore_contrast_preparation_digest(receipt))
}

# --- discovery of preparable contrast units ---------------------------------

lisa_explore_contrast_endpoints <- function(ws, contrast_id) {
  index <- read_lisa_tsv(file.path(ws$source_run, "config", "contrast_index.tsv"))
  if (!nrow(index)) return(NULL)
  names <- paste(as.character(index$contrast_id), as.character(index$output_id),
                 sep = "_")
  keep <- as.character(index$contrast_id) == contrast_id | names == contrast_id
  row <- index[keep, , drop = FALSE]
  if (nrow(row) != 1L) return(NULL)
  list(contrast_id = as.character(row$contrast_id[[1L]]),
       output_id = as.character(row$output_id[[1L]]),
       analysis_a = as.character(row$contrast_a[[1L]]),
       analysis_b = as.character(row$contrast_b[[1L]]))
}

# Every contrast/collection this run can prepare, with the exact source inputs
# each one binds. The catalog is read UNEXPANDED: preparation must not depend on
# the catalog rows it is about to make available.
lisa_explore_contrast_units <- function(ws, contrasts = NULL, collections = NULL,
                                        catalog = NULL) {
  if (is.null(catalog)) catalog <- lisa_extension_discover_catalog(ws$source_run)
  contrast_rows <- unique(catalog[catalog$unit_type == "contrast",
                                  c("contrast_id", "collection", "contrast_path"),
                                  drop = FALSE])
  if (!is.null(contrasts)) {
    contrast_rows <- contrast_rows[contrast_rows$contrast_id %in%
                                     as.character(contrasts), , drop = FALSE]
  }
  if (!is.null(collections)) {
    contrast_rows <- contrast_rows[contrast_rows$collection %in%
                                     as.character(collections), , drop = FALSE]
  }
  if (!nrow(contrast_rows)) return(list())
  single <- unique(catalog[catalog$unit_type == "single_de",
                           c("analysis_id", "collection", "gsea_path", "de_path",
                             "summary_path"), drop = FALSE])
  units <- list()
  for (index in seq_len(nrow(contrast_rows))) {
    row <- contrast_rows[index, , drop = FALSE]
    contrast_id <- as.character(row$contrast_id[[1L]])
    collection <- as.character(row$collection[[1L]])
    endpoints <- lisa_explore_contrast_endpoints(ws, contrast_id)
    if (is.null(endpoints)) next
    sides <- lapply(c(endpoints$analysis_a, endpoints$analysis_b),
                    function(analysis) {
      single[single$analysis_id == analysis & single$collection == collection, ,
             drop = FALSE]
    })
    # A contrast whose endpoints have no gene-level evidence in this collection
    # cannot be prepared. It is skipped rather than half-prepared.
    if (any(vapply(sides, nrow, integer(1)) != 1L)) next
    inputs <- c("config/contrast_index.tsv", "config/de_index.tsv",
                as.character(row$contrast_path[[1L]]),
                unlist(lapply(sides, function(side)
                  c(as.character(side$gsea_path[[1L]]),
                    as.character(side$de_path[[1L]]),
                    as.character(side$summary_path[[1L]]))), use.names = FALSE))
    units[[length(units) + 1L]] <- list(
      contrast_id = contrast_id, collection = collection,
      endpoints = endpoints, sides = sides,
      inputs = sort(unique(inputs[nzchar(inputs)])))
  }
  units
}

# Manifest-bound, by construction. An input that is absent from the canonical
# manifest, or whose bytes differ from it, is refused here rather than silently
# prepared from -- which is what makes the receipt's `inputs` a real provenance
# record instead of a list of paths that happened to exist.
lisa_explore_contrast_verify_inputs <- function(ws, inputs) {
  manifest <- read_lisa_tsv(file.path(ws$source_run, "run_manifest.tsv"))
  matched <- match(inputs, manifest$path)
  if (anyNA(matched)) {
    stop("LISA-EXPLORE-070 a contrast preparation input is absent from the ",
         "canonical manifest: ", paste(inputs[is.na(matched)], collapse = ", "),
         ". Preparation reads only manifest-bound saved evidence.", call. = FALSE)
  }
  observed <- vapply(file.path(ws$source_run, inputs), lisa_sha256_file,
                     character(1))
  recorded <- as.character(manifest$sha256[matched])
  differing <- inputs[unname(observed) != recorded]
  if (length(differing)) {
    stop("LISA-EXPLORE-071 a contrast preparation input differs from the ",
         "canonical manifest: ", paste(differing, collapse = ", "), ".",
         call. = FALSE)
  }
  lapply(seq_along(inputs), function(index) {
    list(path = inputs[[index]], sha256 = unname(observed[[index]]))
  })
}

# --- the preparation itself -------------------------------------------------

#' Prepare the shared native contrast products for an exploration workspace
#'
#' Builds, once per contrast and collection, the paired gene evidence and the
#' complete ranked category summary that the native contrast card, paired
#' heatmap, gene-category network and pathway map all read, plus -- when a local
#' KEGG snapshot is declared -- the complete ranked contrast pathway index.
#'
#' This is data preparation, not analysis. It reads manifest-bound saved
#' evidence of the source run and the declared local snapshot only; it runs no
#' DE, GSEA, ORA or LISA step, writes nothing inside the source run, fetches no
#' resource, and draws no figure: every generator it calls is invoked in a mode
#' that refuses to plot.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param contrasts Optional contrast identifiers to restrict preparation to.
#' @param collections Optional collection identifiers to restrict preparation to.
#' @return Invisibly, a data frame with one row per prepared contrast bundle.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_prepare_contrast_products(ws)
#' }
lisa_explore_prepare_contrast_products <- function(ws, contrasts = NULL,
                                                   collections = NULL) {
  lisa_explore_assert_workspace(ws)
  units <- lisa_explore_contrast_units(ws, contrasts, collections)
  if (!length(units)) {
    stop("LISA-EXPLORE-072 no contrast/collection matched; there is nothing to ",
         "prepare.", call. = FALSE)
  }
  declaration <- lisa_explore_kegg_declaration(ws)
  context <- lisa_explore_contrast_context(ws, declaration)
  settings <- lisa_explore_contrast_settings()
  package_dir <- lisa_resolve_package_dir()
  code_digest <- lisa_explore_contrast_code_digest(package_dir)
  rows <- lapply(units, function(unit) {
    lisa_explore_prepare_one_contrast(ws, unit, declaration, context, settings,
                                      package_dir, code_digest)
  })
  rows <- Filter(Negate(is.null), rows)
  result <- if (!length(rows)) {
    data.frame(contrast_id = character(), collection = character(),
               prepared_dir = character(), preparation_digest = character(),
               n_categories = integer(), n_pathways = integer(),
               stringsAsFactors = FALSE)
  } else do.call(rbind, rows)
  if (nrow(result)) {
    message(sprintf(
      "Prepared %d contrast/collection bundle(s): %d ranked categor(y/ies) and %d ranked pathway(s) in total. No figure was drawn.",
      nrow(result), sum(result$n_categories), sum(result$n_pathways)))
  }
  invisible(result)
}

lisa_explore_prepare_one_contrast <- function(ws, unit, declaration, context,
                                              settings, package_dir,
                                              code_digest) {
  contrast_id <- unit$contrast_id
  collection <- unit$collection
  inputs <- lisa_explore_contrast_verify_inputs(ws, unit$inputs)

  # An isolated private work root. Staging is what keeps the source run
  # read-only: the builders write beside their inputs, so their inputs are
  # copies, never the run's own files.
  work_root <- file.path(ws$root, "prepared", "contrast-work",
                         paste0(lisa_safe_id(contrast_id, "contrast_id"), "__",
                                lisa_safe_id(collection, "collection")))
  if (lisa_path_entry_exists(work_root)) {
    lisa_guarded_delete(work_root, recursive = TRUE, run_root = ws$root)
  }
  lisa_guarded_dir_create(work_root, run_root = ws$root)
  work_root <- lisa_run_root(work_root)
  old_run_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = work_root)
  on.exit(options(lisaR.run_root = old_run_root), add = TRUE)

  lisa_extension_stage_files(ws$source_run, work_root, unit$inputs)

  # The endpoints' gene-level derivative tables. This is the accepted
  # single-DE derivative builder, run over staged copies; it derives tables
  # from saved evidence and performs no enrichment.
  for (side in unit$sides) {
    analysis <- as.character(side$analysis_id[[1L]])
    prefix <- lisa_extension_file_label_prefix(as.character(side$gsea_path[[1L]]),
                                               analysis)
    built <- lisa_build_single_gene_level_tables(analysis, collection,
      file.path(work_root, "outputs"), prefix)
    if (!identical(built$status, "completed")) {
      stop("LISA-EXPLORE-073 the endpoint gene-level evidence needed to pair ",
           contrast_id, " / ", collection, " could not be built for ", analysis,
           ": ", built$message, call. = FALSE)
    }
  }

  code_ledger <- lisa_write_code_identity_ledger(
    package_dir, lisa_explore_contrast_preparation_scripts(),
    file.path(work_root, "code_identity.tsv"),
    executable_recipes = "reproduce_lisa_figure.R")
  run <- function(script, args) {
    result <- lisa_run_post_script(package_dir, script,
      c("--project-dir", work_root, "--contrast-id", contrast_id,
        "--universe", collection, args),
      trusted_run_root = work_root, code_ledger = code_ledger)
    if (!identical(result$status, "completed")) {
      stop("LISA-EXPLORE-074 contrast preparation failed for ", contrast_id,
           " / ", collection, " in ", script, ": ", result$message, call. = FALSE)
    }
  }

  # The whole point: prepare the two tables and draw nothing.
  run("build_contrast_gene_level_product.R",
      c("--tables-only", "true",
        "--top-categories", settings$top_categories,
        "--top-genes-per-category", settings$top_genes_per_category,
        "--lfc-cap", settings$lfc_cap))

  produced <- file.path(work_root, "outputs", "gene_level", "category_contrasts",
                        paste(unit$endpoints$contrast_id,
                              unit$endpoints$output_id, sep = "_"),
                        paste0("collection_", collection))
  prefix <- paste(unit$endpoints$output_id, collection, "contrast_gene_level",
                  sep = "_")
  paired_path <- file.path(produced, paste0(prefix, "_paired_gene_evidence.tsv"))
  summary_path <- file.path(produced,
                            paste0(prefix, "_contrast_category_gene_summary.tsv"))
  if (!file.exists(paired_path) || !file.exists(summary_path)) {
    stop("LISA-EXPLORE-075 the plot-free preparation produced no paired ",
         "evidence for ", contrast_id, " / ", collection, ".", call. = FALSE)
  }
  summary <- read_lisa_tsv(summary_path)

  # The contrast pathway index, from the SAME paired evidence. Index and data
  # preparation only: no map is parsed or drawn, and no resource is fetched.
  pathway_index <- NULL
  kegg_resource_digest <- ""
  if (!is.null(declaration)) {
    kegg_settings <- lisa_explore_kegg_settings(declaration)
    emitted <- file.path(work_root, "contrast-pathway-index.tsv")
    run("build_contrast_kegg_pathway_painter.R",
        c("--species", kegg_settings$kegg_species,
          "--kegg-cache-root", kegg_settings$kegg_cache_root,
          "--kegg-snapshot-id", kegg_settings$kegg_snapshot_id,
          "--kegg-access-mode", "cache_only",
          "--max-abs-log2fc", kegg_settings$max_abs_log2fc,
          "--color-power", kegg_settings$color_power,
          "--paired-input", paired_path, "--summary-input", summary_path,
          "--emit-pathway-index", emitted, "--index-only", "true"))
    if (file.exists(emitted)) {
      pathway_index <- read_lisa_tsv(emitted)
      if (nrow(pathway_index)) {
        organism <- as.character(declaration$organism)
        digests <- vapply(as.character(pathway_index$kegg_id), function(kegg_id) {
          tryCatch(lisa_explore_kegg_pathway_resource_digest(
            kegg_settings$kegg_cache_root, kegg_settings$kegg_snapshot_id,
            organism, kegg_id), error = function(error) "")
        }, character(1))
        pathway_index$kegg_resource_digest <- unname(digests)
        pathway_index$resources_available <- nzchar(pathway_index$kegg_resource_digest)
        kegg_resource_digest <- lisa_sha256_text(paste(
          paste(as.character(pathway_index$kegg_id),
                pathway_index$kegg_resource_digest, sep = "="),
          collapse = "\n"))
      } else pathway_index <- NULL
    }
  }

  # --- promote the bundle --------------------------------------------------
  bundle_dir <- lisa_explore_contrast_bundle_dir(ws, contrast_id, collection)
  if (lisa_path_entry_exists(bundle_dir)) {
    lisa_guarded_delete(bundle_dir, recursive = TRUE, run_root = ws$root)
  }
  lisa_guarded_dir_create(bundle_dir, run_root = ws$root)
  bundle_names <- lisa_explore_contrast_bundle_files()
  lisa_guarded_copy(paired_path, file.path(bundle_dir, bundle_names[["paired"]]),
                    overwrite = FALSE, run_root = ws$root)
  lisa_guarded_copy(summary_path, file.path(bundle_dir, bundle_names[["summary"]]),
                    overwrite = FALSE, run_root = ws$root)
  bundle_files <- unname(bundle_names)
  if (!is.null(pathway_index)) {
    target <- file.path(bundle_dir, lisa_explore_contrast_pathway_index_file())
    previous <- getOption("lisaR.run_root", NULL)
    options(lisaR.run_root = ws$root)
    write_lisa_tsv(pathway_index, target)
    options(lisaR.run_root = previous)
    bundle_files <- c(bundle_files, lisa_explore_contrast_pathway_index_file())
  }
  receipt <- list(
    schema = "lisa-explore-contrast-preparation/1",
    contrast_id = contrast_id, collection = collection,
    output_id = unit$endpoints$output_id,
    analysis_a = unit$endpoints$analysis_a,
    analysis_b = unit$endpoints$analysis_b,
    context = context, settings = settings, code_sha256 = code_digest,
    kegg_resource_digest = kegg_resource_digest,
    inputs = unname(inputs),
    files = unname(lapply(bundle_files, function(relative) {
      list(path = relative,
           sha256 = lisa_sha256_file(file.path(bundle_dir, relative)))
    })),
    prepared_at = lisa_explore_now())
  lisa_explore_write_json(receipt,
                          lisa_explore_contrast_receipt_path(bundle_dir),
                          run_root = ws$root)

  # The work root is scratch, not evidence. Removing it keeps a prepared
  # workspace from accumulating a second copy of every staged input.
  lisa_guarded_delete(work_root, recursive = TRUE, run_root = ws$root)

  data.frame(contrast_id = contrast_id, collection = collection,
             prepared_dir = bundle_dir,
             preparation_digest = lisa_explore_contrast_preparation_digest(receipt),
             n_categories = nrow(summary),
             n_pathways = if (is.null(pathway_index)) 0L else nrow(pathway_index),
             stringsAsFactors = FALSE)
}

# --- the exact-product declaration the catalog reads -------------------------
#
# Reading this touches the prepared receipts and two prepared tables per bundle.
# No generator is reachable from here, which is what keeps navigation free of
# renders.
lisa_explore_contrast_products <- function(ws) {
  root <- lisa_explore_contrast_prepared_root(ws)
  context_dir <- file.path(root, lisa_explore_contrast_context_id(ws))
  if (!dir.exists(context_dir)) return(NULL)
  declaration <- lisa_explore_kegg_declaration(ws)
  code_digest <- lisa_explore_contrast_code_digest()
  units <- list(); categories <- list(); maps <- list()
  for (contrast_id in sort(list.files(context_dir))) {
    contrast_dir <- file.path(context_dir, contrast_id)
    if (!dir.exists(contrast_dir)) next
    for (collection in sort(list.files(contrast_dir))) {
      status <- tryCatch(
        lisa_explore_contrast_bundle_status(ws, contrast_id, collection,
                                            declaration, code_digest),
        error = function(error) list(state = "invalid"))
      if (!identical(status$state, "current")) next
      bundle_names <- lisa_explore_contrast_bundle_files()
      summary <- read_lisa_tsv(file.path(status$dir, bundle_names[["summary"]]))
      if (!nrow(summary) || !"category_id" %in% names(summary)) next
      units[[length(units) + 1L]] <- data.frame(
        contrast_id = contrast_id, collection = collection,
        prepared_dir = status$dir, preparation_digest = status$digest,
        # The endpoints the preparation ACTUALLY paired, read from its receipt.
        # A contrast map binds its A/B DE tables from these, so a row can never
        # name an endpoint the bundle was not built from.
        analysis_a = as.character(status$receipt$analysis_a),
        analysis_b = as.character(status$receipt$analysis_b),
        stringsAsFactors = FALSE)
      categories[[length(categories) + 1L]] <- data.frame(
        contrast_id = contrast_id, collection = collection,
        category_id = as.character(summary$category_id),
        preparation_digest = status$digest, stringsAsFactors = FALSE)
      index_path <- file.path(status$dir,
                              lisa_explore_contrast_pathway_index_file())
      if (!file.exists(index_path) || is.null(declaration)) next
      index <- read_lisa_tsv(index_path)
      if (!nrow(index)) next
      if ("resources_available" %in% names(index)) {
        index <- index[as.logical(index$resources_available) %in% TRUE, ,
                       drop = FALSE]
      }
      if (!nrow(index)) next
      rows <- data.frame(
        contrast_id = contrast_id, collection = collection,
        kegg_id = as.character(index$kegg_id),
        kegg_title = as.character(index$kegg_title),
        rank = as.integer(index$rank),
        kegg_resource_digest = as.character(index$kegg_resource_digest),
        preparation_digest = status$digest, stringsAsFactors = FALSE)
      kegg_settings <- lisa_explore_kegg_settings(declaration)
      for (setting in names(kegg_settings)) rows[[setting]] <- kegg_settings[[setting]]
      maps[[length(maps) + 1L]] <- rows
    }
  }
  if (!length(units)) return(NULL)
  list(units = do.call(rbind, units),
       categories = do.call(rbind, categories),
       maps = if (length(maps)) do.call(rbind, maps) else NULL)
}

lisa_explore_contrast_unavailable_reason <- function(ws, request, generic) {
  product <- as.character(request$product)
  if (identical(product, "de_recurrent_genes")) {
    return(paste0(generic, " A recurrent gene screen needs this analysis and ",
                  "collection to carry gene-level evidence: a GSEA table with ",
                  "a leadingEdge column and a standardized DE table with ",
                  "symbol, log2FoldChange and an adjusted p-value."))
  }
  status <- tryCatch(
    lisa_explore_contrast_bundle_status(ws, as.character(request$contrast_id),
                                        as.character(request$collection)),
    error = function(error) list(state = "invalid",
                                 reason = conditionMessage(error)))
  prepare <- paste0("Run lisa_explore_prepare_contrast_products() for this ",
                    "contrast and collection; it reads manifest-bound saved ",
                    "evidence and the declared local snapshot only, and draws ",
                    "nothing.")
  if (identical(status$state, "absent")) {
    return(paste0("the shared contrast preparation has not been built for ",
                  request$contrast_id, " / ", request$collection, " yet. ",
                  prepare))
  }
  if (!identical(status$state, "current")) {
    return(paste0("the recorded contrast preparation for ", request$contrast_id,
                  " / ", request$collection, " is no longer usable: ",
                  as.character(status$reason %||% status$state), ". ", prepare))
  }
  if (identical(product, "contrast_gene_card")) {
    return(paste0("the prepared ranked summary for ", request$contrast_id,
                  " / ", request$collection, " contains no category ",
                  request$category_id, "."))
  }
  if (!identical(product, "contrast_kegg_map")) return(generic)

  declaration <- lisa_explore_kegg_declaration(ws)
  if (is.null(declaration)) {
    return(paste0("native contrast pathway maps need a local cached KEGG ",
                  "snapshot, and this workspace declares none. Declare one ",
                  "with lisa_explore_configure_kegg() and prepare again; lisaR ",
                  "never downloads or installs a snapshot on your behalf."))
  }
  index_path <- file.path(status$dir, lisa_explore_contrast_pathway_index_file())
  if (!file.exists(index_path)) {
    return(paste0("the contrast pathway index was not prepared for ",
                  request$contrast_id, " / ", request$collection,
                  ", so no map can be offered for it. ", prepare))
  }
  index <- read_lisa_tsv(index_path)
  if (!nrow(index)) {
    return(paste0("the contrast pathway index was prepared for ",
                  request$contrast_id, " / ", request$collection,
                  " and found no KEGG pathway for its paired evidence genes ",
                  "in the declared snapshot."))
  }
  entity <- lisa_explore_blank(request$entity)[[1L]]
  if (!nzchar(entity)) return(generic)
  hit <- index[as.character(index$kegg_id) == entity, , drop = FALSE]
  if (!nrow(hit)) {
    return(paste0("pathway ", entity, " is not among the ", nrow(index),
                  " pathways this contrast's paired evidence genes map to."))
  }
  if ("resources_available" %in% names(hit) &&
      !any(as.logical(hit$resources_available) %in% TRUE)) {
    return(paste0("pathway ", entity, " ranks in this contrast, but the ",
                  "declared KEGG snapshot carries no KGML diagram and base ",
                  "image validated for it. Materialise them into a new ",
                  "snapshot to paint it; lisaR does not download one."))
  }
  generic
}
