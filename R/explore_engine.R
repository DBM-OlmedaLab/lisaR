# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# The exploration engine: three strictly separated operations.
#
#   plan      read saved evidence, decide what a figure needs. Renders nothing.
#   submit    run exactly one figure in a separate process, validate, promote.
#   assemble  read validated artifacts for display or export. Never renders.
#
# Merely navigating produces zero renders. This is enforced by construction:
# nothing on the planning path can reach a generator, because the planning path
# only ever touches the catalog.
# ---------------------------------------------------------------------------

# --- planning (read-only) --------------------------------------------------

lisa_explore_catalog_cache_path <- function(ws, exact_products = NULL) {
  # The cache is keyed on the source manifest hash AND on the exact-product
  # declaration. Without the second half, preparing a KEGG pathway index would
  # leave the catalog serving a cached copy that still knows nothing about the
  # native maps it just made available.
  exact_digest <- substr(lisa_explore_exact_digest(exact_products), 1L, 12L)
  file.path(ws$root, "prepared",
            paste0("catalog-", substr(ws$source_manifest_hash, 1L, 12L),
                   "-", exact_digest, ".tsv"))
}

lisa_explore_exact_digest <- function(exact_products) {
  if (is.null(exact_products)) return(lisa_sha256_text(""))
  maps <- exact_products$kegg_maps
  map_text <- if (is.data.frame(maps) && nrow(maps)) {
    paste(apply(maps[order(maps$analysis_id, maps$collection, maps$category_id,
                           maps$kegg_id), , drop = FALSE], 1L, paste,
                collapse = ":"), collapse = "\n")
  } else ""
  contrast <- exact_products$contrast_products
  contrast_text <- if (!is.null(contrast) && is.data.frame(contrast$units) &&
                       nrow(contrast$units)) {
    units <- contrast$units[order(contrast$units$contrast_id,
                                  contrast$units$collection), , drop = FALSE]
    paste(paste(units$contrast_id, units$collection, units$preparation_digest,
                sep = ":"), collapse = "\n")
  } else ""
  lisa_sha256_text(paste(c(
    paste(sort(as.character(exact_products$heatmap_variants)), collapse = ","),
    paste0("native_single_de=", as.character(isTRUE(exact_products$native_single_de))),
    map_text, contrast_text), collapse = "\n"))
}

lisa_explore_catalog_character_fields <- function(catalog) {
  fields <- intersect(c(
    "unit_type", "analysis_id", "contrast_id", "collection", "category_id",
    "supercategory_id", "product", "entity", "variant", "summary_path",
    "gsea_path", "de_path", "contrast_path", "prepared_dir", "de_a_path",
    "de_b_path", "kegg_id", "kegg_title", "kegg_cache_root",
    "kegg_snapshot_id", "kegg_species", "kegg_resource_digest",
    "preparation_digest"), names(catalog))
  for (field in fields) catalog[[field]] <- as.character(catalog[[field]])
  catalog
}

lisa_explore_read_catalog <- function(path) {
  if (!file.exists(path) || !is.finite(file.info(path)$size) ||
      file.info(path)$size == 0) return(data.frame())
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
                    comment.char = "", check.names = FALSE,
                    colClasses = "character", stringsAsFactors = FALSE)
}

#' Enumerate every figure that can be requested
#'
#' Reads the saved evidence of the source run and returns one row per
#' requestable figure. No generator is invoked and nothing is written inside the
#' source run. A catalog cache may be written in the separate exploration workspace.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param refresh Recompute the catalog instead of reusing the cached copy.
#' @return A data frame with one row per requestable figure.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   ws <- lisa_explore_open(run_dir)
#'   head(lisa_explore_catalog(ws))
#' }
lisa_explore_catalog <- function(ws, refresh = FALSE) {
  lisa_explore_assert_workspace(ws)
  exact_products <- lisa_explore_exact_products(ws)
  cache <- lisa_explore_catalog_cache_path(ws, exact_products)
  # The cache is keyed on the source manifest hash, so a changed source run can
  # never be served from a stale catalog.
  if (!isTRUE(refresh) && file.exists(cache)) {
    cached <- lisa_explore_catalog_character_fields(lisa_explore_read_catalog(cache))
    if (nrow(cached)) return(cached)
  }
  catalog <- lisa_explore_catalog_character_fields(
    lisa_extension_discover_catalog(ws$source_run, exact_products))
  write_lisa_tsv(catalog, cache)
  catalog
}

lisa_explore_blank <- function(value) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  value
}

lisa_explore_catalog_row <- function(ws, request, catalog = NULL) {
  if (is.null(catalog)) catalog <- lisa_explore_catalog(ws)
  keep <- catalog$unit_type == request$unit_type &
    catalog$collection == request$collection &
    lisa_explore_blank(catalog$category_id) ==
      lisa_explore_blank(request$category_id) &
    catalog$product == request$product
  for (field in c("entity", "variant")) {
    catalog_value <- if (field %in% names(catalog)) {
      lisa_explore_blank(catalog[[field]])
    } else rep("", nrow(catalog))
    keep <- keep & catalog_value == lisa_explore_blank(request[[field]])
  }
  owner <- if (identical(request$unit_type, "single_de")) {
    catalog$analysis_id == request$analysis_id
  } else {
    catalog$contrast_id == request$contrast_id
  }
  selected <- keep & owner
  # Data-frame subscripting with NA retains an all-NA row. A malformed catalogue
  # value must therefore fail closed, never masquerade as one exact match.
  selected[is.na(selected)] <- FALSE
  catalog[selected, , drop = FALSE]
}

# Identity of the code and resources that would draw this figure. Combined with
# the source manifest hash and the presentation policy, this defines a reuse
# key that does not depend solely on a PNG existing at a path.
lisa_explore_engine_digest <- function(units) {
  package_dir <- lisa_resolve_package_dir()
  scripts <- lisa_extension_planned_scripts(units)
  rows <- lisa_code_identity_rows(package_dir, scripts)
  text <- if (!nrow(rows)) "" else paste(
    apply(rows[order(rows$relative_path),
               c("package", "package_version", "relative_path", "sha256")],
          1L, paste, collapse = ":"),
    collapse = "\n")
  lisa_sha256_text(paste(c(text, as.character(utils::packageVersion("lisaR"))),
                         collapse = "\n"))
}

# Finding 6 of the parent review. Everything in the reuse key except the request
# itself and the planned-script set is invariant across the whole catalog, yet
# `lisa_explore_status()` recomputed all of it once per row: a file read for the
# GSEA cutoff, a policy digest, and -- by far the worst -- a package directory
# scan plus SHA-256 of every code-identity row for `lisa_explore_engine_digest()`.
# On a 160-row catalog that is 160 rescans of the package, on the exact path the
# UI calls on every navigation.
#
# The context below computes each invariant once and memoises the engine digest
# per planned-script set, which is the only part that genuinely varies. The key
# VALUE is unchanged -- the same inputs are hashed in the same order -- so reuse
# keys minted before this change still match. Nothing is weakened: no component
# is dropped from the fingerprint.
lisa_explore_key_context <- function(ws, policy = lisa_explore_report_policy()) {
  structure(list(
    gsea_padj_cutoff = format(lisa_extension_gsea_padj_cutoff(ws$source_run),
                              digits = 17),
    source_manifest_hash = ws$source_manifest_hash,
    policy_digest = lisa_explore_policy_digest(policy),
    engine_digests = new.env(parent = emptyenv())
  ), class = "lisa_explore_key_context")
}

lisa_explore_context_engine_digest <- function(context, row) {
  scripts <- tryCatch(lisa_extension_planned_scripts(row),
                      error = function(error) character())
  # Hash the planned-script set so the cache key is always a single, valid,
  # non-empty name -- including when a row plans no scripts at all.
  cache_key <- lisa_sha256_text(paste(sort(unique(as.character(scripts))),
                                      collapse = "\n"))
  if (exists(cache_key, envir = context$engine_digests, inherits = FALSE)) {
    return(get(cache_key, envir = context$engine_digests, inherits = FALSE))
  }
  digest <- lisa_explore_engine_digest(row)
  assign(cache_key, digest, envir = context$engine_digests)
  digest
}

lisa_explore_resource_digest <- function(row) {
  columns <- c("kegg_species", "kegg_snapshot_id", "kegg_cache_root",
               "max_abs_log2fc", "color_power", "kegg_resource_digest",
               "preparation_digest")
  values <- vapply(columns, function(column) {
    if (!column %in% names(row)) return("")
    lisa_explore_blank(row[[column]])[[1L]]
  }, character(1))
  if (!any(nzchar(values))) return("")
  paste(paste(columns, values, sep = "="), collapse = "\n")
}

lisa_explore_key <- function(ws, request, row,
                             policy = lisa_explore_report_policy(),
                             context = NULL) {
  if (is.null(context)) context <- lisa_explore_key_context(ws, policy)
  # The IDENTITY vector, not the whole request vector. They are the same thing
  # for every product except the native pathway map, where the category is a
  # navigation association and must not split one map into two artifacts. See
  # `lisa_explore_identity_fields()`.
  values <- lisa_explore_identity_vector(request)
  resources <- lisa_explore_resource_digest(row)
  fingerprint <- paste(c(
    paste(paste(names(values), values, sep = "="), collapse = "\n"),
    context$source_manifest_hash,
    context$gsea_padj_cutoff,
    context$policy_digest,
    lisa_explore_context_engine_digest(context, row),
    if (nzchar(resources)) resources else NULL
  ), collapse = "\n")
  paste0("fig-", substr(lisa_sha256_text(fingerprint), 1L, 24L))
}

lisa_explore_unavailable_reason <- function(ws, request) {
  # Scope-aware wording. Saying "for category <NA>" about a figure that has no
  # category is worse than unhelpful: it suggests the reader chose the wrong
  # category when the product does not have one at all.
  scope <- if (identical(lisa_explore_product_scope(request$product), "category")) {
    paste0(" for category ", request$category_id)
  } else ""
  generic <- paste0("the saved evidence of this run offers no ", request$product,
                    scope, " in collection ", request$collection, ".")
  if (request$product %in% lisa_explore_h3_products()) {
    return(lisa_explore_contrast_unavailable_reason(ws, request, generic))
  }
  if (!identical(request$product, "kegg_pathway_map")) {
    if (identical(request$product, "heatmap") && !is.na(request$variant)) {
      return(paste0(generic, " The ", request$variant, " scale needs a saved ",
                    "expression matrix and sample annotation for this analysis."))
    }
    return(generic)
  }
  declaration <- lisa_explore_kegg_declaration(ws)
  if (is.null(declaration)) {
    return(paste0("native KEGG pathway maps need a local cached KEGG snapshot, ",
                  "and this workspace declares none. Declare one with ",
                  "lisa_explore_configure_kegg(); lisaR never downloads or ",
                  "installs a snapshot on your behalf."))
  }
  preparation <- lisa_explore_kegg_preparation_status(ws, declaration)
  index <- lisa_explore_kegg_index(ws)
  if (identical(preparation$state, "stale")) {
    return(paste0("the recorded native KEGG preparation belongs to a different ",
                  "source run or snapshot context. Run ",
                  "lisa_explore_prepare_kegg_index() for this context; it reads ",
                  "saved evidence and local snapshot files only, and paints ",
                  "nothing."))
  }
  if (identical(preparation$state, "invalid")) {
    return(paste0("the native KEGG preparation provenance does not match its ",
                  "local pathway index. Re-run ",
                  "lisa_explore_prepare_kegg_index(); lisaR does not download ",
                  "or paint anything during preparation."))
  }
  if (identical(preparation$state, "absent") &&
      file.exists(lisa_explore_kegg_index_path(ws))) {
    return(paste0("this workspace has a legacy native KEGG pathway index that ",
                  "does not declare which analysis/collection units were ",
                  "prepared. Whether ", request$analysis_id, " / ",
                  request$collection, " was attempted cannot be determined; ",
                  "re-run lisa_explore_prepare_kegg_index() for the needed ",
                  "units. It uses saved evidence and local snapshot files only ",
                  "and paints nothing."))
  }
  if (identical(preparation$state, "absent")) {
    return(paste0("the native KEGG pathway index has not been prepared for this ",
                  "workspace yet. Run lisa_explore_prepare_kegg_index() once; it ",
                  "reads saved evidence and the declared snapshot only, and ",
                  "paints nothing."))
  }
  prepared_units <- preparation$units
  prepared_here <- nrow(prepared_units) && any(
    prepared_units$analysis_id == request$analysis_id &
      prepared_units$collection == request$collection)
  if (!prepared_here) {
    return(paste0("the native KEGG pathway index has not been prepared for ",
                  request$analysis_id, " / ", request$collection,
                  ". Run lisa_explore_prepare_kegg_index() for this unit; it ",
                  "reads saved evidence and the declared local snapshot only, ",
                  "and paints nothing."))
  }
  if (is.null(index)) index <- lisa_explore_kegg_empty_index()
  owned <- index[as.character(index$analysis_id) == request$analysis_id &
                 as.character(index$collection) == request$collection, ,
                 drop = FALSE]
  if (!nrow(owned)) {
    return(paste0("the native KEGG pathway index was prepared for ",
                  request$analysis_id, " / ", request$collection,
                  " and found no pathway associations for its saved evidence ",
                  "in the declared snapshot."))
  }
  entity <- lisa_explore_blank(request$entity)[[1L]]
  if (!nzchar(entity)) {
    # Asked about the product as a whole rather than about one pathway: the
    # reader wants to know why this category offers no map at all.
    associated <- owned[as.character(owned$category_id) == request$category_id,
                        , drop = FALSE]
    if (nrow(associated) && "resources_available" %in% names(associated) &&
        !any(as.logical(associated$resources_available) %in% TRUE)) {
      return(paste0("the prepared index contains ",
                    length(unique(as.character(associated$kegg_id))),
                    " known pathway association(s) for category ",
                    request$category_id, ", but the declared local snapshot ",
                    "has no validated KGML diagram and base image for them. ",
                    "Materialise those resources into a new snapshot to paint ",
                    "a map; lisaR does not download them."))
    }
    return(paste0("the prepared index contains no pathway association for ",
                  "category ", request$category_id,
                  " in this analysis/collection."))
  }
  if (!entity %in% as.character(owned$kegg_id)) {
    return(paste0("pathway ", entity, " is not among the ",
                  length(unique(as.character(owned$kegg_id))),
                  " pathways this collection's evidence genes map to."))
  }
  paintable <- owned[as.character(owned$kegg_id) == entity, , drop = FALSE]
  if ("resources_available" %in% names(paintable) &&
      !any(as.logical(paintable$resources_available) %in% TRUE)) {
    # The pathway is real and ranked; the declared snapshot simply does not carry
    # its diagram. Naming that precisely is the difference between a usable
    # message and "not applicable".
    return(paste0("pathway ", entity, " ranks in this collection, but the ",
                  "declared KEGG snapshot carries no KGML diagram and base ",
                  "image validated for it. Materialise them into a new snapshot ",
                  "to paint it; ",
                  "lisaR does not download one."))
  }
  paste0("pathway ", request$entity, " is not associated with category ",
         request$category_id, " in this collection's prepared index.")
}

#' Plan one exact figure without rendering it
#'
#' Determines whether a figure is applicable, which saved inputs it needs and
#' what its current state is. Nothing is rendered and nothing is scheduled.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param request A [lisa_figure_request()].
#' @return A list describing the request, its state and its inputs.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_plan(ws, request)
#' }
lisa_explore_plan <- function(ws, request) {
  lisa_explore_assert_workspace(ws)
  lisa_explore_assert_request(request)
  row <- lisa_explore_catalog_row(ws, request)
  request_id <- lisa_request_id(request)
  if (!nrow(row)) {
    return(list(request = request, request_id = request_id, key = NA_character_,
                applicable = FALSE, state = "not_applicable",
                reason = lisa_explore_unavailable_reason(ws, request),
                inputs = character(), expected_files = 0L))
  }
  if (nrow(row) != 1L) {
    # Defensive: the catalog is deduplicated upstream, so this cannot normally
    # happen. If it ever does we refuse rather than render an ambiguous figure.
    stop("LISA-EXPLORE-020 the request resolved to ", nrow(row),
         " catalog rows; exactly one is required.", call. = FALSE)
  }
  path_columns <- grep("_path$", names(row), value = TRUE)
  inputs <- sort(unique(unlist(row[path_columns], use.names = FALSE)))
  inputs <- inputs[nzchar(inputs)]
  policy <- lisa_explore_report_policy()
  key <- lisa_explore_key(ws, request, row, policy)
  entry <- lisa_explore_index_find(lisa_explore_read_index(ws), key)
  state <- lisa_explore_entry_state(ws, entry)
  list(request = request, request_id = request_id, key = key,
       applicable = TRUE, state = state$state, reason = state$reason,
       inputs = inputs,
       expected_files = sum(unlist(policy$formats)) +
         as.integer(isTRUE(policy$source_data)) +
         as.integer(isTRUE(policy$recipes)),
       artifact_dir = state$artifact_dir, job_id = state$job_id)
}

# --- artifact validity -----------------------------------------------------

# Complete-download validation. Every recorded file must be present, non-empty
# and hash-identical. A partially promoted or corrupted artifact is never served
# as available.
lisa_explore_validate_artifact <- function(ws, entry) {
  if (is.null(entry) || is.null(entry$artifact_dir)) {
    return(list(ok = FALSE, reason = "no artifact recorded"))
  }
  dir <- file.path(ws$root, "extensions", as.character(entry$artifact_dir))
  if (!dir.exists(dir)) {
    return(list(ok = FALSE, reason = "the artifact directory is absent"))
  }
  files <- entry$files
  if (is.null(files) || !length(files)) {
    return(list(ok = FALSE, reason = "no files recorded for the artifact"))
  }
  for (file in files) {
    path <- file.path(dir, as.character(file$path))
    if (!file.exists(path)) {
      return(list(ok = FALSE,
                  reason = paste0("a recorded file is missing: ", file$path)))
    }
    size <- as.numeric(file.info(path)$size)
    if (!is.finite(size) || size <= 0) {
      return(list(ok = FALSE,
                  reason = paste0("a recorded file is empty: ", file$path)))
    }
    if (!identical(lisa_sha256_file(path), as.character(file$sha256))) {
      return(list(ok = FALSE,
                  reason = paste0("a recorded file changed on disk: ", file$path)))
    }
  }
  list(ok = TRUE, reason = "")
}

lisa_explore_entry_state <- function(ws, entry) {
  empty <- list(state = "ungenerated", reason = "", artifact_dir = NULL,
                job_id = NULL)
  if (is.null(entry)) return(empty)
  state <- as.character(entry$state)
  job_id <- if (is.null(entry$job_id)) NULL else as.character(entry$job_id)
  if (identical(state, "available")) {
    check <- lisa_explore_validate_artifact(ws, entry)
    if (!isTRUE(check$ok)) {
      # A recorded artifact that no longer validates is not available. It is
      # reported as ungenerated with the concrete reason, never served as if it
      # were still good.
      return(list(state = "ungenerated", reason = check$reason,
                  artifact_dir = NULL, job_id = job_id))
    }
    return(list(state = "available", reason = "",
                artifact_dir = as.character(entry$artifact_dir), job_id = job_id))
  }
  if (state %in% c("queued", "running")) {
    job <- if (is.null(job_id)) NULL else lisa_explore_read_job(ws, job_id)
    reconciled <- lisa_explore_reconcile_job(ws, job)
    return(list(state = reconciled$state, reason = reconciled$message,
                artifact_dir = NULL, job_id = job_id))
  }
  list(state = state,
       reason = if (is.null(entry$message)) "" else as.character(entry$message),
       artifact_dir = NULL, job_id = job_id)
}

# A job whose owning process has gone is an orphan. It is reported failed. It is
# never reported available and never silently relaunched.
lisa_explore_reconcile_job <- function(ws, job) {
  if (is.null(job)) {
    return(list(state = "failed", message = "the job record is missing"))
  }
  state <- as.character(job$state)
  subject <- if (identical(lisa_explore_job_kind(job), "export"))
    "export" else "figure"
  if (!state %in% c("queued", "running")) {
    return(list(state = state,
                message = if (is.null(job$message)) "" else as.character(job$message)))
  }
  if (isTRUE(job$cancel_requested) && identical(state, "queued")) {
    return(list(state = "failed", message = "cancelled before starting"))
  }
  if (identical(state, "running") &&
      !lisa_explore_identity_alive(list(host = job$host, pid = job$pid,
                                        start_time = job$start_time))) {
    return(list(state = "failed",
                message = paste("the worker process ended without completing the",
                                subject)))
  }
  if (identical(state, "queued")) {
    # A job that has not been dispatched yet is genuinely waiting for the single
    # worker slot. It has no owning process, so there is nothing to orphan.
    if (!isTRUE(job$launched)) {
      return(list(state = "queued",
                  message = paste("waiting for the single worker slot to run the",
                                  subject)))
    }
    # A background job that never reached `running` because its worker died at
    # startup is an orphan too, and must not sit in the queue for ever.
    launcher <- lisa_explore_read_json(
      lisa_explore_launcher_path(ws, as.character(job$job_id)))
    if (!is.null(launcher) && !lisa_explore_identity_alive(launcher)) {
      return(list(state = "failed",
                  message = paste("the worker process ended before starting the",
                                  subject)))
    }
    if (is.null(launcher) && !is.null(job$dispatcher) &&
        !lisa_explore_identity_alive(job$dispatcher)) {
      # Marked for dispatch, but the session that claimed the slot died before
      # the child existed. Without this the slot would stay occupied for ever.
      return(list(state = "failed",
                  message = "the dispatching session ended before the worker started"))
    }
  }
  list(state = state,
       message = if (is.null(job$stage)) "" else as.character(job$stage))
}

# --- the single worker slot ------------------------------------------------
#
# Finding 2 of the parent review. Deduplication by key stops the *same* figure
# being rendered twice, but the previous submit path released the index lock and
# then launched immediately, so two *different* requests -- or two sessions over
# the same workspace -- each started their own `callr` worker. The approved
# policy is one costly worker at a time, and a disabled button in one UI cannot
# enforce it for a second session.
#
# The slot is therefore derived from the job records and claimed *inside the
# index lock*, which is the only thing both sessions serialise on. A submission
# that cannot claim the slot is recorded as queued-but-undispatched and is
# started later by `lisa_explore_pump()`; it is never dropped and never run
# concurrently.

# The job currently occupying the worker slot, or NULL. Only a job that has
# actually been dispatched can hold the slot.
lisa_explore_slot_holder <- function(ws, jobs = NULL) {
  if (is.null(jobs)) jobs <- lisa_explore_all_jobs(ws)
  for (job in jobs) {
    state <- as.character(job$state)
    if (!state %in% c("queued", "running")) next
    if (!isTRUE(job$launched)) next
    if (identical(lisa_explore_reconcile_job(ws, job)$state, "failed")) next
    return(job)
  }
  NULL
}

# The oldest dispatchable job: explicitly submitted, still queued, not
# cancelled, and not yet launched. Ordered by creation so the queue is fair.
lisa_explore_next_dispatchable <- function(ws, jobs = NULL) {
  if (is.null(jobs)) jobs <- lisa_explore_all_jobs(ws)
  waiting <- Filter(function(job) {
    identical(as.character(job$state), "queued") && !isTRUE(job$launched) &&
      !isTRUE(job$cancel_requested) && isTRUE(job$background)
  }, jobs)
  if (!length(waiting)) return(NULL)
  created <- vapply(waiting, function(job) as.character(job$created_at), character(1))
  waiting[[order(created, seq_along(created))[[1L]]]]
}

# Claim the slot for `job_id` if it is free. Must be called with the index lock
# held: the read of the slot and the write that claims it have to be atomic with
# respect to the other session.
lisa_explore_claim_slot <- function(ws, job_id) {
  jobs <- lisa_explore_all_jobs(ws)
  holder <- lisa_explore_slot_holder(ws, jobs)
  if (!is.null(holder) && !identical(as.character(holder$job_id),
                                     as.character(job_id))) {
    return(list(claimed = FALSE, blocked_by = as.character(holder$job_id)))
  }
  job <- lisa_explore_read_job(ws, job_id)
  if (is.null(job)) return(list(claimed = FALSE, blocked_by = NA_character_))
  if (isTRUE(job$cancel_requested)) {
    return(list(claimed = FALSE, blocked_by = NA_character_))
  }
  job$launched <- TRUE
  job$dispatcher <- lisa_explore_process_identity()
  job$launched_at <- lisa_explore_now()
  lisa_explore_write_job(ws, job)
  list(claimed = TRUE, blocked_by = NA_character_)
}

#' Start the next queued figure if the single worker slot is free
#'
#' Exploration runs one costly worker at a time. A submission made while another
#' figure is rendering waits in the queue; this function dispatches the oldest
#' waiting job once the slot frees. It only ever starts work that was already
#' explicitly submitted, so it never turns browsing into rendering.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param lock_timeout_seconds Maximum time to wait for the short index-lock
#'   claim. Interactive polling can set zero so a busy or unverifiable lock is
#'   reported without blocking the Shiny event loop.
#' @return Invisibly, the dispatched job identifier, or `NULL`.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_pump(ws)
#' }
lisa_explore_pump <- function(ws, lock_timeout_seconds = 30) {
  lisa_explore_assert_workspace(ws)
  if (!requireNamespace("callr", quietly = TRUE)) return(invisible(NULL))
  candidate <- lisa_explore_with_lock(ws, {
    jobs <- lisa_explore_all_jobs(ws)
    if (!is.null(lisa_explore_slot_holder(ws, jobs))) {
      NULL
    } else {
      next_job <- lisa_explore_next_dispatchable(ws, jobs)
      if (is.null(next_job)) {
        NULL
      } else {
        claim <- lisa_explore_claim_slot(ws, next_job$job_id)
        if (isTRUE(claim$claimed)) as.character(next_job$job_id) else NULL
      }
    }
  }, timeout_seconds = lock_timeout_seconds)
  if (is.null(candidate)) return(invisible(NULL))
  # dispatch path. Same reasoning as in `lisa_explore_submit()`: the claim
  # is inside the lock, the launch is not, so a launch failure has to release the
  # slot by marking the job failed. Without this the pump would hand the slot to
  # a job that cannot start and then never hand it to anything else.
  launched <- tryCatch({
    lisa_explore_launch_background(ws, candidate)
    TRUE
  }, error = function(error) {
    lisa_explore_mark_job_failed(ws, candidate,
      paste0("the worker process could not be started: ",
             conditionMessage(error)))
    FALSE
  })
  if (!isTRUE(launched)) return(invisible(NULL))
  invisible(candidate)
}

# --- submission ------------------------------------------------------------

lisa_explore_new_job_id <- function() {
  paste0("job-", format(Sys.time(), "%Y%m%dT%H%M%S"), "-",
         substr(lisa_sha256_text(paste(Sys.getpid(), Sys.time(),
                                       stats::runif(1))), 1L, 8L))
}

#' Submit one exact figure for generation
#'
#' Schedules exactly one figure. Repeating a submission that is already
#' validated, queued or running returns the existing result or job instead of
#' starting a second one.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param request A [lisa_figure_request()].
#' @param background Run in a separate R process. Requires the optional `callr`
#'   package; set `FALSE` to run in the calling process.
#' @return A list with the job identifier, the reuse key and the resulting state.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_submit(ws, request)
#' }
lisa_explore_submit <- function(ws, request, background = TRUE) {
  lisa_explore_assert_workspace(ws)
  lisa_explore_assert_request(request)
  if (!request$product %in% lisa_explore_generatable_products()) {
    stop("LISA-EXPLORE-021 this engine can currently generate only: ",
         paste(lisa_explore_generatable_products(), collapse = ", "),
         ". The request named ", request$product,
         ", which is planned for a later milestone.", call. = FALSE)
  }
  plan <- lisa_explore_plan(ws, request)
  if (!isTRUE(plan$applicable)) {
    stop("LISA-EXPLORE-022 the figure is not applicable to this run: ",
         plan$reason, call. = FALSE)
  }
  if (isTRUE(background) && !requireNamespace("callr", quietly = TRUE)) {
    stop("LISA-EXPLORE-023 background generation needs the optional package ",
         "'callr'. Install it, or call with background = FALSE.", call. = FALSE)
  }

  decision <- lisa_explore_with_lock(ws, {
    entries <- lisa_explore_read_index(ws)
    entry <- lisa_explore_index_find(entries, plan$key)
    state <- lisa_explore_entry_state(ws, entry)

    # Record the category this request navigates from. For a shared native
    # pathway map the second category adds itself here and reuses the artifact:
    # one map, one paint, two places it can be read from. `associations` never
    # takes part in the reuse key, so adding one cannot invalidate the artifact.
    if (!is.null(entry)) {
      merged <- lisa_explore_union_associations(entry$associations,
                                                request$category_id)
      if (!identical(merged, lisa_explore_union_associations(entry$associations))) {
        entry$associations <- merged
        entries <- lisa_explore_index_upsert(entries, entry)
        lisa_explore_write_index(ws, entries)
      }
    }

    if (identical(state$state, "available")) {
      list(action = "reuse", job_id = state$job_id, state = "available")
    } else if (state$state %in% c("queued", "running")) {
      # Double click, reload or a second session join the job in flight; they
      # never start a duplicate.
      list(action = "join", job_id = state$job_id, state = state$state)
    } else if (!isTRUE(background) &&
               !is.null(foreground_blocker <- lisa_explore_slot_holder(ws))) {
      # A *foreground* caller is itself the worker, so a submission it
      # cannot start cannot be deferred either. Previously the job record and
      # the `queued` index entry were written first and the refusal happened
      # afterwards, which left behind a job with `background = FALSE` and
      # `launched = FALSE`. `lisa_explore_next_dispatchable()` only ever
      # dispatches background jobs, so the pump could never start it, and every
      # later submit for the same request -- background or not -- found the
      # `queued` entry and *joined* the ghost instead of starting a real render.
      # The figure became ungeneratable until someone cancelled a job they had
      # no reason to know existed.
      #
      # Deciding the claim before anything is persisted removes the failure
      # mode rather than cleaning up after it: a refused foreground submit
      # writes no job file and no index entry at all, so the very next submit
      # for the same figure starts a real render as soon as the slot is free.
      list(action = "refuse", job_id = NA_character_, state = "ungenerated",
           blocked_by = as.character(foreground_blocker$job_id))
    } else {
      job_id <- lisa_explore_new_job_id()
      job <- list(
        job_id = job_id, key = plan$key, request_id = plan$request_id,
        request = as.list(lisa_explore_request_vector(request)),
        source_run = ws$source_run, root = ws$root,
        state = "queued", stage = "queued", cancel_requested = FALSE,
        background = isTRUE(background), launched = FALSE,
        host = lisa_explore_host(), pid = NA_integer_,
        start_time = NA_character_,
        created_at = lisa_explore_now(), message = ""
      )
      lisa_explore_write_job(ws, job)
      entries <- lisa_explore_index_upsert(entries, list(
        key = plan$key, request_id = plan$request_id,
        request = as.list(lisa_explore_request_vector(request)),
        associations = lisa_explore_union_associations(request$category_id),
        state = "queued", job_id = job_id,
        source_manifest_hash = ws$source_manifest_hash,
        updated_at = lisa_explore_now(), message = "",
        artifact_dir = NULL, files = list()
      ))
      lisa_explore_write_index(ws, entries)

      # Claim the single worker slot here, while the lock is still held. This is
      # what stops two *different* requests -- or two sessions -- from running
      # two costly workers at once.
      claim <- lisa_explore_claim_slot(ws, job_id)
      if (!isTRUE(claim$claimed) && !isTRUE(background)) {
        lisa_explore_mark_job_failed(ws, job_id,
          "the worker slot was taken before this foreground submission started")
        list(action = "refuse", job_id = NA_character_, state = "failed",
             blocked_by = claim$blocked_by)
      } else {
        list(action = if (isTRUE(claim$claimed)) "start" else "wait",
             job_id = job_id, state = "queued", blocked_by = claim$blocked_by)
      }
    }
  })

  if (identical(decision$action, "refuse")) {
    stop("LISA-EXPLORE-036 another figure is already being generated (job ",
         decision$blocked_by, "). Exploration runs one worker at a time; wait ",
         "for it to finish, or submit with background = TRUE to queue this ",
         "figure.", call. = FALSE)
  }

  if (identical(decision$action, "start")) {
    if (isTRUE(background)) {
      # The slot is claimed inside the lock but the child is started
      # outside it, so a launch failure used to leave the job `queued` with
      # `launched = TRUE` and a live dispatcher. `lisa_explore_reconcile_job()`
      # only fails such a job once the launcher or the dispatcher is dead, so
      # while the Shiny session stayed alive every later request queued behind a
      # job that could never run. Marking it failed here -- atomically, inside
      # the index lock, by the same path a worker failure takes -- frees the
      # slot immediately, so the next valid job runs without restarting R.
      launched <- tryCatch({
        lisa_explore_launch_background(ws, decision$job_id)
        TRUE
      }, error = function(error) {
        lisa_explore_mark_job_failed(ws, decision$job_id,
          paste0("the worker process could not be started: ",
                 conditionMessage(error)))
        FALSE
      })
      if (!isTRUE(launched)) {
        # The slot is free again: hand it to whatever was waiting behind this
        # job rather than leaving the queue stalled by a failure.
        lisa_explore_pump(ws)
      }
    } else {
      lisa_explore_execute_job(ws$root, decision$job_id)
    }
    job <- lisa_explore_read_job(ws, decision$job_id)
    decision$state <- lisa_explore_reconcile_job(ws, job)$state
    # The slot is free again the moment a foreground job finishes, so honour any
    # request that queued behind it rather than leaving it stranded.
    if (!isTRUE(background)) lisa_explore_pump(ws)
  }

  list(job_id = decision$job_id, request_id = plan$request_id, key = plan$key,
       state = decision$state, reused = !decision$action %in% c("start", "wait"),
       queued_behind = if (identical(decision$action, "wait"))
         decision$blocked_by else NA_character_)
}

# Mark one job failed and record the failure in the index, atomically under the
# index lock. A failed job is skipped by `lisa_explore_slot_holder()`, so this is
# also how a claimed-but-unstartable job releases the single worker slot without
# any separate "release" step that could itself fail half-way.
lisa_explore_mark_job_failed <- function(ws, job_id, message) {
  lisa_explore_with_lock(ws, {
    job <- lisa_explore_read_job(ws, job_id)
    if (is.null(job)) return(invisible(NULL))
    finished <- lisa_explore_job_set(ws, job, state = "failed", stage = "failed",
                                     message = message,
                                     launched = FALSE,
                                     finished_at = lisa_explore_now())
    if (identical(lisa_explore_job_kind(job), "figure")) {
      entries <- lisa_explore_index_upsert(lisa_explore_read_index(ws),
                                           lisa_explore_failed_entry(ws, job, message))
      lisa_explore_write_index(ws, entries)
    }
    invisible(finished)
  })
}

#' Report the state of exploration jobs
#'
#' Reconciles each job with its owning process before reporting, so an orphaned
#' job is reported as failed rather than as running or available.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param job_id Optional job identifier; all jobs are reported when omitted.
#' @param dispatch Hand the free worker slot to the oldest already-submitted
#'   figure that is still waiting. Set `FALSE` to report state without
#'   dispatching anything.
#' @return A data frame with one row per job.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_poll(ws)
#' }
lisa_explore_poll <- function(ws, job_id = NULL, dispatch = TRUE) {
  lisa_explore_assert_workspace(ws)
  # Release handles of workers that have already finished, then hand the free
  # slot to the oldest figure that was explicitly submitted while it was busy.
  # Polling never *creates* work: `lisa_explore_pump()` can only start a job that
  # some earlier `lisa_explore_submit()` already put in the queue.
  lisa_explore_reap_processes()
  if (isTRUE(dispatch)) lisa_explore_pump(ws)
  jobs <- if (is.null(job_id)) {
    lisa_explore_all_jobs(ws)
  } else {
    job <- lisa_explore_read_job(ws, lisa_explore_scalar(job_id, "job_id"))
    if (is.null(job)) list() else list(job)
  }
  if (!length(jobs)) {
    return(data.frame(job_id = character(), key = character(),
                      request_id = character(), state = character(),
                      stage = character(), message = character(),
                      stringsAsFactors = FALSE))
  }
  rows <- lapply(jobs, function(job) {
    reconciled <- lisa_explore_reconcile_job(ws, job)
    data.frame(job_id = as.character(job$job_id), key = as.character(job$key),
               request_id = as.character(job$request_id),
               state = reconciled$state,
               stage = if (is.null(job$stage)) "" else as.character(job$stage),
               message = reconciled$message, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Request cancellation of an exploration job
#'
#' A queued job is cancelled at once. A running job is asked to stop and honours
#' the request at its next safe point; no process is ever forcibly terminated.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param job_id The job identifier to cancel.
#' @return A list describing the cancellation outcome.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_cancel(ws, job_id)
#' }
lisa_explore_cancel <- function(ws, job_id) {
  lisa_explore_assert_workspace(ws)
  job_id <- lisa_explore_scalar(job_id, "job_id")
  lisa_explore_with_lock(ws, {
    job <- lisa_explore_read_job(ws, job_id)
    if (is.null(job)) {
      stop("LISA-EXPLORE-024 unknown job: ", job_id, call. = FALSE)
    }
    if (!as.character(job$state) %in% c("queued", "running")) {
      return(list(job_id = job_id, state = as.character(job$state),
                  cancelled = FALSE,
                  message = "the job had already finished"))
    }
    job$cancel_requested <- TRUE
    if (identical(as.character(job$state), "queued")) {
      job$state <- "failed"
      job$message <- "cancelled before starting"
      lisa_explore_write_job(ws, job)
      if (identical(lisa_explore_job_kind(job), "figure")) {
        entries <- lisa_explore_index_upsert(
          lisa_explore_read_index(ws),
          lisa_explore_failed_entry(ws, job, "cancelled before starting"))
        lisa_explore_write_index(ws, entries)
      }
      return(list(job_id = job_id, state = "failed", cancelled = TRUE,
                  message = "cancelled before starting"))
    }
    lisa_explore_write_job(ws, job)
    list(job_id = job_id, state = "running", cancelled = FALSE,
         message = paste("cancellation requested; the worker will stop at its",
                         "next safe point"))
  })
}

lisa_explore_failed_entry <- function(ws, job, message) {
  list(key = as.character(job$key), request_id = as.character(job$request_id),
       request = job$request, state = "failed",
       job_id = as.character(job$job_id),
       source_manifest_hash = ws$source_manifest_hash,
       updated_at = lisa_explore_now(), message = message,
       artifact_dir = NULL, files = list())
}

# --- entry currency (finding 5) --------------------------------------------
#
# An index entry is *current* only when the key it would be minted with today is
# the key it carries. That single comparison covers every drift the review names,
# because all of them are inputs to the key:
#
#   * the source run changed          -> source_manifest_hash differs
#   * the presentation policy changed -> policy digest differs
#   * the GSEA cutoff changed         -> cutoff differs
#   * lisaR or a generator changed    -> engine/code identity digest differs
#
# and separately: the source run may no longer offer the figure at all, in which
# case there is no catalog row to key against.

lisa_explore_entry_currency <- function(ws, entry, catalog = NULL, context = NULL) {
  if (is.null(entry)) return(list(current = FALSE, reason = "no index entry"))
  recorded_hash <- as.character(entry$source_manifest_hash)
  if (length(recorded_hash) == 1L && !is.na(recorded_hash) &&
      !identical(recorded_hash, ws$source_manifest_hash)) {
    return(list(current = FALSE,
                reason = "the source run changed after this figure was produced"))
  }
  request <- tryCatch(
    lisa_explore_request_from_row(as.data.frame(entry$request,
                                                stringsAsFactors = FALSE)),
    error = function(error) NULL)
  if (is.null(request)) {
    return(list(current = FALSE, reason = "the recorded request is unreadable"))
  }
  if (is.null(catalog)) catalog <- lisa_explore_catalog(ws)
  row <- lisa_explore_catalog_row(ws, request, catalog)
  if (nrow(row) != 1L) {
    return(list(current = FALSE,
                reason = "the source run no longer offers this figure"))
  }
  if (is.null(context)) context <- lisa_explore_key_context(ws)
  current_key <- lisa_explore_key(ws, request, row, context = context)
  if (!identical(current_key, as.character(entry$key))) {
    return(list(current = FALSE,
                reason = paste0("the engine or resources changed after this ",
                                "figure was produced (current key ", current_key,
                                ", recorded ", as.character(entry$key), ")")))
  }
  list(current = TRUE, reason = "")
}

lisa_explore_entry_is_current <- function(ws, entry, catalog = NULL,
                                          context = NULL) {
  isTRUE(lisa_explore_entry_currency(ws, entry, catalog, context)$current)
}

#' Report index entries that are no longer current
#'
#' Lists artifacts whose own files still validate but whose reuse key no longer
#' matches the key the same request would be minted with today, because the
#' source run, the presentation policy or the engine has changed. These entries
#' are deliberately hidden from [lisa_explore_artifacts()] and from every export;
#' this function exists so their exclusion is visible rather than silent.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @return A data frame of stale entries with the reason each is stale.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_stale_artifacts(ws)
#' }
lisa_explore_stale_artifacts <- function(ws) {
  lisa_explore_assert_workspace(ws)
  entries <- lisa_explore_read_index(ws)
  catalog <- lisa_explore_catalog(ws)
  context <- lisa_explore_key_context(ws)
  rows <- lapply(entries, function(entry) {
    if (!identical(lisa_explore_entry_state(ws, entry)$state, "available")) {
      return(NULL)
    }
    currency <- lisa_explore_entry_currency(ws, entry, catalog, context)
    if (isTRUE(currency$current)) return(NULL)
    data.frame(key = as.character(entry$key),
               request_id = as.character(entry$request_id),
               reason = currency$reason, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(key = character(), request_id = character(),
                      reason = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

# --- assembly (read-only) --------------------------------------------------

# The live shell and static exporter must attach an available artifact to the
# same category set. The catalogue-derived status already carries the complete
# reuse key, so exact-key equality gives every currently valid association while
# preventing a pathway from leaking into another source, analysis, collection,
# settings or resource context. Recorded submit associations are deliberately
# not the authority: they may be only the subset of categories that clicked.
lisa_explore_catalogue_associations <- function(status, key,
                                                 fallback = character()) {
  values <- character()
  if (is.data.frame(status) && nrow(status) &&
      all(c("key", "category_id") %in% names(status))) {
    values <- lisa_explore_blank(
      status$category_id[lisa_explore_blank(status$key) == as.character(key)])
  }
  values <- values[nzchar(values)]
  if (!length(values)) values <- lisa_explore_blank(fallback)
  sort(unique(values[nzchar(values)]))
}

#' List validated exploration artifacts
#'
#' Reads the index and reports only artifacts that still validate. Nothing is
#' rendered.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param request Optional [lisa_figure_request()] to report a single figure.
#' @return A data frame of validated artifacts.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_artifacts(ws)
#' }
lisa_explore_artifacts <- function(ws, request = NULL) {
  lisa_explore_assert_workspace(ws)
  entries <- lisa_explore_read_index(ws)
  if (!is.null(request)) {
    plan <- lisa_explore_plan(ws, request)
    entries <- Filter(function(entry) identical(as.character(entry$key), plan$key),
                      entries)
  }
  catalog <- lisa_explore_catalog(ws)
  context <- lisa_explore_key_context(ws)
  rows <- lapply(entries, function(entry) {
    state <- lisa_explore_entry_state(ws, entry)
    if (!identical(state$state, "available")) return(NULL)
    # Finding 5 of the parent review. File hashes matching the *recorded* hashes
    # only proves the artifact has not been corrupted since it was written. It
    # says nothing about whether the source run, the presentation policy or the
    # engine that produced it are still the current ones. Recompute the key that
    # this entry's request would get *now* and require it to match; a stale
    # entry is therefore never displayed or exported as current, even though its
    # own files still validate perfectly.
    if (!isTRUE(lisa_explore_entry_is_current(ws, entry, catalog, context))) {
      return(NULL)
    }
    fields <- entry$request
    data.frame(
      key = as.character(entry$key),
      request_id = as.character(entry$request_id),
      unit_type = as.character(fields$unit_type),
      analysis_id = as.character(fields$analysis_id),
      contrast_id = as.character(fields$contrast_id),
      collection = as.character(fields$collection),
      category_id = as.character(fields$category_id),
      product = as.character(fields$product),
      entity = lisa_explore_blank(fields$entity),
      variant = lisa_explore_blank(fields$variant),
      # Semicolon-joined because this is a data frame column, not a key. The
      # exporter splits it again to attach one artifact at every category view
      # it belongs to.
      associations = paste(unlist(lisa_explore_union_associations(
        entry$associations, fields$category_id)), collapse = ";"),
      artifact_dir = as.character(entry$artifact_dir),
      n_files = length(entry$files),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(key = character(), request_id = character(),
                      unit_type = character(), analysis_id = character(),
                      contrast_id = character(), collection = character(),
                      category_id = character(), product = character(),
                      entity = character(), variant = character(),
                      associations = character(),
                      artifact_dir = character(), n_files = integer(),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Report the state of every requestable figure
#'
#' Joins the read-only catalog with the index so an interface can render the six
#' product states without scheduling anything.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @return The catalog with a `state` column added.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   table(lisa_explore_status(ws)$state)
#' }
lisa_explore_status <- function(ws) {
  lisa_explore_assert_workspace(ws)
  catalog <- lisa_explore_catalog(ws)
  entries <- lisa_explore_read_index(ws)
  # One context for the whole catalog instead of one per row (finding 6).
  context <- lisa_explore_key_context(ws)
  resolved <- lapply(seq_len(nrow(catalog)), function(index) {
    row <- catalog[index, , drop = FALSE]
    request <- lisa_explore_request_from_row(row)
    key <- lisa_explore_key(ws, request, row, context = context)
    state <- lisa_explore_entry_state(ws, lisa_explore_index_find(entries, key))
    c(key = key, state = state$state, reason = state$reason)
  })
  catalog$state <- vapply(resolved, `[[`, character(1), "state")
  # The key is reported alongside the state because it, not the request id, is
  # the artifact's identity. Two categories sharing one native pathway map have
  # two request ids and one key, and an interface that matched on request id
  # would show the second category "ungenerated" for a map already on disk.
  catalog$key <- vapply(resolved, `[[`, character(1), "key")
  catalog$reason <- vapply(resolved, `[[`, character(1), "reason")
  catalog
}
