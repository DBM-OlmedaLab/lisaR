# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Background export jobs share the existing exploration worker slot with figure
# renders. They are ordinary persisted jobs, but they never write a figure index
# entry: their immutable payload is the exact set of validated artifact keys that
# existed when the user requested the snapshot.

lisa_explore_job_kind <- function(job) {
  kind <- if (is.null(job$kind)) "figure" else as.character(job$kind)
  if (!kind %in% c("figure", "export")) "unknown" else kind
}

lisa_explore_export_signature <- function(ws, output_dir, format,
                                          allow_absolute_references,
                                          snapshot_keys) {
  text <- paste(c(
    lisa_path_key(output_dir), as.character(format),
    as.character(isTRUE(allow_absolute_references)),
    ws$source_manifest_hash, sort(unique(as.character(snapshot_keys)))),
    collapse = "\n")
  paste0("export-", substr(lisa_sha256_text(text), 1L, 24L))
}

lisa_explore_snapshot_rows <- function(value) {
  if (is.data.frame(value)) return(value)
  if (is.null(value) || !length(value)) return(data.frame())
  rows <- lapply(value, function(row) {
    if (is.data.frame(row)) row else
      as.data.frame(lapply(row, function(cell) {
        cell <- unlist(cell, use.names = FALSE)
        if (!length(cell)) "" else as.character(cell[[1L]])
      }), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

lisa_explore_latest_export_job <- function(ws) {
  jobs <- Filter(function(job) identical(lisa_explore_job_kind(job), "export"),
                 lisa_explore_all_jobs(ws))
  if (!length(jobs)) return(NULL)
  created <- vapply(jobs, function(job) as.character(job$created_at), character(1L))
  jobs[[order(created, seq_along(created), decreasing = TRUE)[[1L]]]]
}

lisa_explore_export_job_status <- function(ws, job) {
  if (is.null(job)) return("")
  reconciled <- lisa_explore_reconcile_job(ws, job)
  state <- as.character(reconciled$state)
  stage <- if (is.null(job$stage)) "" else as.character(job$stage)
  detail <- if (is.null(job$message)) "" else as.character(job$message)
  destination <- if (is.null(job$destination)) "" else as.character(job$destination)
  if (identical(state, "completed")) {
    return(paste("Saved expanded report to", destination))
  }
  if (identical(state, "failed")) {
    message <- if (nzchar(reconciled$message)) reconciled$message else detail
    return(paste("Export failed:", message))
  }
  prefix <- if (identical(state, "queued")) "Export queued" else "Export running"
  cancel_pending <- isTRUE(job$cancel_requested) &&
    state %in% c("queued", "running")
  pieces <- c(prefix, if (nzchar(stage)) paste0("stage: ", stage),
              if (nzchar(detail) && !identical(detail, stage)) detail,
              if (cancel_pending)
                "cancellation requested; stopping at the next safe point")
  paste(pieces, collapse = " - ")
}

lisa_explore_export_owner_path <- function(output_dir, staging_id = NULL) {
  suffix <- if (is.null(staging_id)) "" else {
    staging_id <- lisa_explore_scalar(staging_id, "staging_id")
    if (!grepl("^[0-9a-f]{24}$", staging_id)) {
      stop("Unsafe export staging identity.", call. = FALSE)
    }
    paste0("-", staging_id)
  }
  file.path(dirname(output_dir), paste0(".", basename(output_dir),
                                       ".export-owner", suffix, ".json"))
}

lisa_explore_claim_export_staging <- function(ws, job, output_dir) {
  legacy_staging <- lisa_explore_export_staging_path(output_dir)
  legacy_owner <- lisa_explore_export_owner_path(output_dir)
  prior <- lisa_explore_read_json(legacy_owner)
  if (lisa_path_entry_exists(legacy_staging)) {
    if (is.null(prior)) {
      stop("LISA-EXPLORE-080 export staging exists without a verifiable owner. ",
           "It was retained for manual inspection and no destination was ",
           "published: ", legacy_staging, call. = FALSE)
    }
    if (lisa_explore_identity_alive(prior)) {
      stop("LISA-EXPLORE-081 another live export owner is using staging for ",
           output_dir, ". It was left untouched.", call. = FALSE)
    }
    lisa_guarded_delete(legacy_staging, recursive = TRUE, run_root = NULL)
    lisa_guarded_delete(legacy_owner, recursive = FALSE, run_root = NULL)
  } else if (!is.null(prior)) {
    if (lisa_explore_identity_alive(prior)) {
      stop("LISA-EXPLORE-081 another live export owner is preparing ", output_dir,
           ". It was left untouched.", call. = FALSE)
    }
    lisa_guarded_delete(legacy_owner, recursive = FALSE, run_root = NULL)
  }
  token <- lisa_explore_new_lock_token()
  staging <- lisa_explore_export_staging_path(output_dir, token)
  owner_path <- lisa_explore_export_owner_path(output_dir, token)
  if (!dir.create(staging, showWarnings = FALSE, mode = "0700")) {
    stop("LISA-EXPLORE-081 an isolated export staging claim could not be ",
         "created. Existing entries were left untouched: ", staging,
         call. = FALSE)
  }
  owner <- c(lisa_explore_process_identity(), list(
    token = token, job_id = as.character(job$job_id),
    destination = output_dir, claimed_at = lisa_explore_now()))
  written <- tryCatch({
    lisa_explore_write_json(owner, owner_path, run_root = NULL)
    TRUE
  }, error = function(error) error)
  if (inherits(written, "error")) {
    if (!length(list.files(staging, all.files = TRUE, no.. = TRUE))) {
      lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
    }
    stop(written)
  }
  token
}

lisa_explore_release_export_staging <- function(output_dir, token) {
  owner_path <- lisa_explore_export_owner_path(output_dir, token)
  owner <- lisa_explore_read_json(owner_path)
  if (is.null(owner) || is.null(owner$token) ||
      !identical(as.character(owner$token), as.character(token))) {
    return(invisible(FALSE))
  }
  lisa_guarded_delete(owner_path, recursive = FALSE, run_root = NULL)
  invisible(TRUE)
}

lisa_explore_finish_export_failed <- function(ws, job, message,
                                              elapsed = NA_real_) {
  lisa_explore_with_lock(ws, {
    current <- lisa_explore_read_job(ws, as.character(job$job_id))
    if (is.null(current)) current <- job
    lisa_explore_job_set(ws, current, state = "failed", stage = "failed",
                         message = message, launched = FALSE,
                         finished_at = lisa_explore_now(),
                         elapsed_seconds = elapsed)
  })
}

lisa_explore_execute_export_job <- function(root, job_id) {
  root <- lisa_existing_run_root(root)
  job_id <- lisa_explore_scalar(job_id, "job_id")
  job <- lisa_explore_read_json(file.path(root, "jobs", paste0(job_id, ".json")))
  if (is.null(job)) stop("LISA-EXPLORE-024 unknown job: ", job_id, call. = FALSE)
  if (!identical(lisa_explore_job_kind(job), "export")) {
    stop("LISA-EXPLORE-082 job ", job_id, " is not an export job.", call. = FALSE)
  }
  ws <- lisa_explore_open(as.character(job$source_run), root)
  job <- lisa_explore_read_job(ws, job_id)
  if (isTRUE(job$cancel_requested)) {
    return(invisible(lisa_explore_finish_export_failed(
      ws, job, "cancelled before starting")))
  }

  started <- Sys.time()
  job <- lisa_explore_with_lock(ws, {
    identity <- lisa_explore_process_identity()
    lisa_explore_job_set(ws, job, state = "running",
      stage = "validating snapshot", message = "validating the requested snapshot",
      pid = identity$pid, host = identity$host, start_time = identity$start_time,
      started_at = lisa_explore_now())
  })
  output_dir <- as.character(job$destination)
  token <- NULL
  result <- tryCatch({
    token <- lisa_explore_claim_export_staging(ws, job, output_dir)
    on.exit(if (!is.null(token))
      lisa_explore_release_export_staging(output_dir, token), add = TRUE)
    progress <- function(stage, message = stage) {
      lisa_explore_job_progress(ws, job_id, stage = stage, message = message)
      invisible(NULL)
    }
    should_cancel <- function() {
      current <- lisa_explore_read_job(ws, job_id)
      isTRUE(current$cancel_requested)
    }
    lisa_explore_export(ws, output_dir, format = as.character(job$format),
      allow_absolute_references = isTRUE(job$allow_absolute_references),
      snapshot_artifacts = lisa_explore_snapshot_rows(job$snapshot_artifacts),
      progress = progress, should_cancel = should_cancel,
      staging_id = token)
  }, error = function(error) error)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(result, "error")) {
    return(invisible(lisa_explore_finish_export_failed(
      ws, lisa_explore_read_job(ws, job_id), conditionMessage(result), elapsed)))
  }
  invisible(lisa_explore_with_lock(ws, {
    current <- lisa_explore_read_job(ws, job_id)
    lisa_explore_job_set(ws, current, state = "completed", stage = "completed",
      message = paste("saved", result$standard_files, "STANDARD files and",
                      result$extension_files, "extension files"),
      launched = FALSE, finished_at = lisa_explore_now(),
      elapsed_seconds = elapsed, standard_files = result$standard_files,
      extension_files = result$extension_files, extensions = result$extensions,
      archive = result$archive)
  }))
}

lisa_explore_execute_queued_job <- function(root, job_id) {
  root <- lisa_existing_run_root(root)
  job <- lisa_explore_read_json(file.path(root, "jobs", paste0(job_id, ".json")))
  if (identical(lisa_explore_job_kind(job), "export")) {
    lisa_explore_execute_export_job(root, job_id)
  } else {
    lisa_explore_execute_job(root, job_id)
  }
}

lisa_explore_submit_export <- function(ws, destination, format = c("dir", "zip"),
                                       allow_absolute_references = FALSE,
                                       background = TRUE) {
  lisa_explore_assert_workspace(ws)
  format <- match.arg(format)
  destination <- lisa_explore_scalar(destination, "destination")
  resolved <- lisa_explore_export_destination(ws, destination)
  output_dir <- resolved$path
  if (isTRUE(background) && !requireNamespace("callr", quietly = TRUE)) {
    stop("LISA-EXPLORE-023 background export needs the optional package 'callr'.",
         call. = FALSE)
  }

  decision <- lisa_explore_with_lock(ws, {
    jobs <- lisa_explore_all_jobs(ws)
    active_exports <- Filter(function(job) {
      identical(lisa_explore_job_kind(job), "export") &&
        lisa_explore_reconcile_job(ws, job)$state %in% c("queued", "running")
    }, jobs)
    if (length(active_exports)) {
      existing <- active_exports[[1L]]
      if (!identical(lisa_path_key(as.character(existing$destination)),
                     lisa_path_key(output_dir))) {
        stop("LISA-EXPLORE-083 export job ", as.character(existing$job_id),
             " is already queued or running. Wait for it to finish before ",
             "requesting a later snapshot at a new destination.", call. = FALSE)
      }
      return(list(action = "join", job_id = as.character(existing$job_id),
                  state = lisa_explore_reconcile_job(ws, existing)$state,
                  key = as.character(existing$key),
                  snapshot_count = as.integer(existing$snapshot_count),
                  destination = output_dir))
    }
    if (!isTRUE(background) &&
        !is.null(foreground_blocker <- lisa_explore_slot_holder(ws))) {
      stop("LISA-EXPLORE-084 another figure or export already owns the single ",
           "worker slot (job ", as.character(foreground_blocker$job_id),
           "). No export was queued; wait for it to finish, or submit with ",
           "background = TRUE.", call. = FALSE)
    }
    # This is the snapshot instant. It is serialized with index publication, and
    # artifact validation happens before the keys are persisted.
    available <- lisa_explore_artifacts(ws)
    keys <- sort(unique(as.character(available$key)))
    if (nrow(available)) available <- available[match(keys, available$key), , drop = FALSE]
    snapshot_rows <- lapply(seq_len(nrow(available)), function(index) {
      as.list(available[index, , drop = FALSE])
    })
    signature <- lisa_explore_export_signature(
      ws, output_dir, format, allow_absolute_references, keys)
    active <- Filter(function(job) {
      identical(lisa_explore_job_kind(job), "export") &&
        identical(as.character(job$key), signature) &&
        lisa_explore_reconcile_job(ws, job)$state %in% c("queued", "running")
    }, jobs)
    if (length(active)) {
      existing <- active[[1L]]
      list(action = "join", job_id = as.character(existing$job_id),
           state = lisa_explore_reconcile_job(ws, existing)$state,
           key = signature, snapshot_count = length(keys))
    } else {
      job_id <- lisa_explore_new_job_id()
      job <- list(
        kind = "export", job_id = job_id, key = signature,
        request_id = signature, request = list(), source_run = ws$source_run,
        source_manifest_hash = ws$source_manifest_hash, root = ws$root,
        destination = output_dir, format = format,
        allow_absolute_references = isTRUE(allow_absolute_references),
        snapshot_keys = as.list(keys), snapshot_artifacts = snapshot_rows,
        snapshot_count = length(keys),
        state = "queued", stage = "queued",
        message = paste("snapshot contains", length(keys), "available figure(s)"),
        cancel_requested = FALSE, background = isTRUE(background),
        launched = FALSE, host = lisa_explore_host(), pid = NA_integer_,
        start_time = NA_character_, created_at = lisa_explore_now())
      lisa_explore_write_job(ws, job)
      claim <- lisa_explore_claim_slot(ws, job_id)
      list(action = if (isTRUE(claim$claimed)) "start" else "wait",
           job_id = job_id, state = "queued", key = signature,
           snapshot_count = length(keys), blocked_by = claim$blocked_by)
    }
  })

  if (identical(decision$action, "start")) {
    if (isTRUE(background)) {
      launched <- tryCatch({
        lisa_explore_launch_background(ws, decision$job_id)
        TRUE
      }, error = function(error) {
        lisa_explore_mark_job_failed(ws, decision$job_id,
          paste0("the export worker process could not be started: ",
                 conditionMessage(error)))
        FALSE
      })
      if (!isTRUE(launched)) lisa_explore_pump(ws)
    } else {
      lisa_explore_execute_export_job(ws$root, decision$job_id)
      decision$state <- lisa_explore_reconcile_job(
        ws, lisa_explore_read_job(ws, decision$job_id))$state
      lisa_explore_pump(ws)
    }
  }
  decision$destination <- output_dir
  decision
}
