# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# The generation worker.
#
# One worker renders one figure. It receives the pinned R environment and its
# library paths explicitly, so it does not depend on the working directory or on
# an interactive session.
#
# Promotion happens only after the produced set of files has been validated.
# A failure leaves a failed attempt, never an available figure.
# ---------------------------------------------------------------------------

# --- live process registry -------------------------------------------------
#
# Process handles must stay reachable until their workers finish:
#
#   * `callr::r_bg()` collects `...` into `options$extra` and builds a
#     `processx::process` whose `cleanup` argument keeps its default TRUE.
#     `cleanup` is NOT reachable through r_bg -- passing `cleanup = FALSE` is
#     silently swallowed. So every callr child is registered for cleanup.
#   * cleanup = TRUE installs a finaliser on the handle. If the handle becomes
#     unreachable and the garbage collector runs, THE CHILD IS KILLED. Probe 1
#     reproduced exactly this: handle dropped, `gc()` forced, the worker died
#     before writing its marker.
#   * Retaining a strong reference makes the child survive `gc()` (probe 2, Q5:
#     marker written, exit status 0).
#
# The previous implementation returned the handle `invisible()` and
# `lisa_explore_submit()` discarded it, so every background render was one
# garbage collection away from being killed mid-figure. The registry below is
# the fix: a strong reference lives in the package namespace for as long as the
# job runs.

.lisa_explore_processes <- new.env(parent = emptyenv())

lisa_explore_register_process <- function(job_id, process) {
  assign(as.character(job_id), process, envir = .lisa_explore_processes)
  invisible(process)
}

lisa_explore_process_handle <- function(job_id) {
  job_id <- as.character(job_id)
  if (!exists(job_id, envir = .lisa_explore_processes, inherits = FALSE)) return(NULL)
  get(job_id, envir = .lisa_explore_processes, inherits = FALSE)
}

lisa_explore_forget_process <- function(job_id) {
  job_id <- as.character(job_id)
  if (exists(job_id, envir = .lisa_explore_processes, inherits = FALSE)) {
    rm(list = job_id, envir = .lisa_explore_processes)
  }
  invisible(NULL)
}

# Drop handles whose process has finished. Called from the polling path so the
# registry cannot grow without bound in a long-lived session. A handle is only
# released once the process is no longer alive, so releasing can never kill a
# running worker.
lisa_explore_reap_processes <- function() {
  for (job_id in ls(envir = .lisa_explore_processes, all.names = TRUE)) {
    handle <- lisa_explore_process_handle(job_id)
    finished <- tryCatch(!isTRUE(handle$is_alive()), error = function(error) TRUE)
    if (finished) lisa_explore_forget_process(job_id)
  }
  invisible(NULL)
}

# A thin seam over `callr::r_bg()`. Its only purpose is that a test can make the
# launch fail without rendering anything: injecting a failure here exercises the
# whole claim/launch/release path at zero render cost.
lisa_explore_start_process <- function(...) {
  callr::r_bg(...)
}

lisa_explore_launch_background <- function(ws, job_id) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("LISA-EXPLORE-023 background generation needs the optional package 'callr'.",
         call. = FALSE)
  }
  log_dir <- file.path(ws$root, "jobs")
  process <- lisa_explore_start_process(
    func = function(root, job_id) {
      lisaR:::lisa_explore_execute_queued_job(root, job_id)
    },
    args = list(root = ws$root, job_id = job_id),
    libpath = .libPaths(),
    stdout = file.path(log_dir, paste0(job_id, ".out")),
    stderr = file.path(log_dir, paste0(job_id, ".err")),
    # supervise = FALSE: the supervisor exists to kill children when the parent
    # dies, which processx's own cleanup already does here. It bought nothing and
    # added a third process to the budget. Probe 2 Q6 showed parent exit ends the
    # child with supervise either way, so this is not a behaviour change.
    supervise = FALSE
  )
  # Keep a strong reference for the lifetime of the job. Without this the child
  # is killed by the next garbage collection (probe 1, Q2).
  lisa_explore_register_process(job_id, process)

  # Record the child PID in its own file. Writing it back into the job record
  # would race with the child, which may already have advanced that record to
  # `running`; the parent would then overwrite it back to `queued`.
  lisa_explore_write_json(
    c(lisa_explore_process_identity(process$get_pid()),
      list(launched_at = lisa_explore_now())),
    lisa_explore_launcher_path(ws, job_id), run_root = ws$root)
  invisible(process)
}

lisa_explore_launcher_path <- function(ws, job_id) {
  file.path(ws$root, "jobs", paste0(job_id, ".launcher.json"))
}


lisa_explore_staging_prefix <- function(key) paste0(".", key, ".staging-")

lisa_explore_staging_owner_path <- function(ws, key) {
  file.path(ws$root, "extensions", paste0(".", key, ".staging-owner.json"))
}

lisa_explore_staging_directories <- function(ws, key) {
  extensions <- file.path(ws$root, "extensions")
  if (!dir.exists(extensions)) return(character())
  entries <- list.files(extensions, all.files = TRUE, no.. = TRUE)
  entries <- entries[startsWith(entries, lisa_explore_staging_prefix(key))]
  paths <- file.path(extensions, entries)
  paths[dir.exists(paths)]
}

lisa_explore_clear_staging_owner <- function(ws, key) {
  marker <- lisa_explore_staging_owner_path(ws, key)
  if (file.exists(marker)) {
    lisa_guarded_delete(marker, recursive = FALSE, run_root = ws$root)
  }
  invisible(NULL)
}

lisa_explore_recover_staging <- function(ws, job, key) {
  job_id <- as.character(job$job_id)
  stale <- lisa_explore_staging_directories(ws, key)
  if (length(stale)) {
    owner <- lisa_explore_read_json(lisa_explore_staging_owner_path(ws, key))
    if (!is.null(owner) && !identical(as.character(owner$job_id), job_id) &&
        lisa_explore_identity_alive(owner)) {
      stop("LISA-EXPLORE-039 another live worker (job ",
           as.character(owner$job_id), ") is rendering this figure; its ",
           "staging directory was left untouched. Wait for it to finish.",
           call. = FALSE)
    }
    for (path in stale) {
      lisa_guarded_delete(path, recursive = TRUE, run_root = ws$root)
    }
  }
  lisa_explore_clear_staging_owner(ws, key)
  lisa_explore_write_json(
    c(lisa_explore_process_identity(),
      list(job_id = job_id, key = key, claimed_at = lisa_explore_now())),
    lisa_explore_staging_owner_path(ws, key), run_root = ws$root)
  invisible(length(stale))
}

lisa_explore_job_set <- function(ws, job, ...) {
  updates <- list(...)
  for (name in names(updates)) job[[name]] <- updates[[name]]
  job$updated_at <- lisa_explore_now()
  lisa_explore_write_job(ws, job)
  job
}

# Progress and cancellation update the same persisted job record. Serialize the
# read-modify-write under the workspace lock so a worker can never write an old
# `cancel_requested = FALSE` value over a concurrent cancellation request.
lisa_explore_job_progress <- function(ws, job_id, stage, message = stage) {
  lisa_explore_with_lock(ws, {
    current <- lisa_explore_read_job(ws, job_id)
    if (is.null(current)) {
      stop("LISA-EXPLORE-024 unknown job: ", job_id, call. = FALSE)
    }
    if (isTRUE(current$cancel_requested)) {
      stop("cancelled at the safe point before ", stage, call. = FALSE)
    }
    lisa_explore_job_set(ws, current, stage = stage, message = message)
  })
}

lisa_explore_peak_memory_mb <- function() {
  gc_mb <- tryCatch(sum(as.numeric(gc(full = TRUE)[, 6L])),
                    error = function(error) NA_real_)
  rss_mb <- tryCatch({
    if (requireNamespace("ps", quietly = TRUE)) {
      as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^2
    } else NA_real_
  }, error = function(error) NA_real_)
  list(gc_max_mb = gc_mb, rss_mb = rss_mb)
}

#' Execute one queued exploration job
#'
#' Renders exactly one figure, validates the produced files and promotes the
#' result. This is the worker entry point; it is called in a separate process by
#' [lisa_explore_submit()] and is not normally invoked directly.
#'
#' @param root The exploration workspace root.
#' @param job_id The queued job identifier.
#' @return Invisibly, the final job record.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_execute_job(ws$root, job_id)
#' }
lisa_explore_execute_job <- function(root, job_id) {
  root <- lisa_existing_run_root(root)
  job_id <- lisa_explore_scalar(job_id, "job_id")
  job_file <- file.path(root, "jobs", paste0(job_id, ".json"))
  job <- lisa_explore_read_json(job_file)
  if (is.null(job)) stop("LISA-EXPLORE-024 unknown job: ", job_id, call. = FALSE)
  ws <- lisa_explore_open(as.character(job$source_run), root)

  # Safe point one: honour a cancellation that arrived while queued.
  job <- lisa_explore_read_job(ws, job_id)
  if (isTRUE(job$cancel_requested)) {
    return(invisible(lisa_explore_finish_failed(ws, job, "cancelled before starting")))
  }

  started <- Sys.time()
  gc(reset = TRUE, full = TRUE)
  job <- lisa_explore_with_lock(ws, {
    # Record the full process identity, not only the number, so a recycled PID
    # cannot make this job look alive after the worker is gone.
    identity <- lisa_explore_process_identity()
    updated <- lisa_explore_job_set(ws, job, state = "running", stage = "planning",
                                    pid = identity$pid, host = identity$host,
                                    start_time = identity$start_time,
                                    started_at = lisa_explore_now())
    entries <- lisa_explore_index_upsert(lisa_explore_read_index(ws), list(
      key = as.character(job$key), request_id = as.character(job$request_id),
      request = job$request, state = "running", job_id = job_id,
      source_manifest_hash = ws$source_manifest_hash,
      updated_at = lisa_explore_now(), message = "",
      artifact_dir = NULL, files = list()))
    lisa_explore_write_index(ws, entries)
    updated
  })

  result <- tryCatch({
    request <- lisa_explore_request_from_row(as.data.frame(job$request,
                                                           stringsAsFactors = FALSE))
    policy <- lisa_explore_report_policy()
    selection <- lisa_explore_selection(request, policy)

    # Verify the plan resolves to exactly one unit BEFORE rendering. This is the
    # structural guarantee that no sibling figure can be drawn: the renderer
    # receives a single row, so its `--categories` argument carries one category
    # and its product list carries one product.
    # The same exact-product declaration the catalog was built from, so the plan
    # sees the very row the reader chose -- this scale, this pathway -- and not a
    # broader family that would then have to be filtered after rendering.
    exact_products <- lisa_explore_exact_products(ws)
    extension_plan <- plan_lisa_extension(ws$source_run, selection, exact_products)
    if (nrow(extension_plan$units) != 1L) {
      stop("LISA-EXPLORE-025 the request expanded to ",
           nrow(extension_plan$units),
           " render units; exactly one is required. Refusing to render.",
           call. = FALSE)
    }

    # Safe point two: the render itself is one indivisible unit.
    job <- lisa_explore_read_job(ws, job_id)
    if (isTRUE(job$cancel_requested)) {
      stop("LISA-EXPLORE-026 cancelled at the safe point before rendering.",
           call. = FALSE)
    }
    job <- lisa_explore_job_progress(ws, job_id, stage = "rendering")

    artifact_dir <- file.path(ws$root, "extensions", as.character(job$key))
    if (lisa_path_entry_exists(artifact_dir)) {
      # A previous attempt left an artifact that no longer validates. Remove it
      # inside our own workspace before rendering again.
      lisa_guarded_delete(artifact_dir, recursive = TRUE, run_root = ws$root)
    }
    # `render_lisa_categories()` stages under a *deterministic* hidden
    # sibling of the artifact directory, and refuses with LISA-EXTENSION-022 if
    # that directory already exists. A worker killed mid-render -- session exit,
    # OOM, power loss, a cooperative memory pause that ends the process -- never
    # runs its `on.exit` cleanup, so the staging directory survives and every
    # later retry of that same figure fails with an error that names no recovery.
    # The worker's own pre-render cleanup above removed only `extensions/<key>`,
    # never the staging sibling.
    lisa_explore_recover_staging(ws, job, as.character(job$key))
    on.exit(lisa_explore_clear_staging_owner(ws, as.character(job$key)),
            add = TRUE)
    render_lisa_categories(ws$source_run, selection, artifact_dir, exact_products)

    job <- lisa_explore_job_progress(ws, job_id, stage = "validating")
    files <- lisa_explore_collect_artifact_files(artifact_dir)
    lisa_explore_assert_complete_product(files, request, policy)
    list(artifact_dir = artifact_dir, files = files)
  }, error = function(error) {
    structure(list(message = conditionMessage(error)), class = "lisa_explore_failure")
  })

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  memory <- lisa_explore_peak_memory_mb()

  if (inherits(result, "lisa_explore_failure")) {
    return(invisible(lisa_explore_finish_failed(ws, lisa_explore_read_job(ws, job_id),
                                                result$message, elapsed, memory)))
  }

  invisible(lisa_explore_with_lock(ws, {
    finished <- lisa_explore_job_set(ws, lisa_explore_read_job(ws, job_id),
                                     state = "completed", stage = "completed",
                                     message = "", finished_at = lisa_explore_now(),
                                     elapsed_seconds = elapsed,
                                     gc_max_mb = memory$gc_max_mb,
                                     rss_mb = memory$rss_mb)
    entries <- lisa_explore_index_upsert(lisa_explore_read_index(ws), list(
      key = as.character(job$key), request_id = as.character(job$request_id),
      request = job$request, state = "available", job_id = job_id,
      source_manifest_hash = ws$source_manifest_hash,
      updated_at = lisa_explore_now(), message = "",
      artifact_dir = basename(result$artifact_dir),
      files = result$files,
      elapsed_seconds = elapsed))
    lisa_explore_write_index(ws, entries)
    finished
  }))
}

lisa_explore_finish_failed <- function(ws, job, message, elapsed = NA_real_,
                                       memory = list(gc_max_mb = NA_real_,
                                                     rss_mb = NA_real_)) {
  lisa_explore_with_lock(ws, {
    finished <- lisa_explore_job_set(ws, job, state = "failed", stage = "failed",
                                     message = message,
                                     finished_at = lisa_explore_now(),
                                     elapsed_seconds = elapsed,
                                     gc_max_mb = memory$gc_max_mb,
                                     rss_mb = memory$rss_mb)
    entries <- lisa_explore_index_upsert(lisa_explore_read_index(ws),
                                         lisa_explore_failed_entry(ws, job, message))
    lisa_explore_write_index(ws, entries)
    finished
  })
}

lisa_explore_collect_artifact_files <- function(artifact_dir) {
  tree <- lisa_scan_run_tree(artifact_dir)
  paths <- sort(tree$path[!tree$isdir])
  relative <- substring(paths, nchar(normalizePath(artifact_dir, winslash = "/",
                                                   mustWork = TRUE)) + 2L)
  lapply(seq_along(paths), function(index) {
    list(path = relative[[index]],
         bytes = as.numeric(file.info(paths[[index]])$size),
         sha256 = lisa_sha256_file(paths[[index]]))
  })
}

# The file name a product's source data actually carries.
#
# The accepted base assumed `_source.tsv` for everything. That held while volcano
# was the only generatable product, but the leading-edge heatmap builder has
# always written its source table as `<stem>_matrix.tsv` -- so the first heatmap
# ever requested would have rendered correctly and then been refused as
# incomplete. Naming the real convention per product fixes that without relaxing
# the check: something must still match, and it must be the right something.
lisa_explore_source_data_pattern <- function(product) {
  switch(as.character(product),
    heatmap = "_matrix[.]tsv$",
    contrast_kegg_map = "_contrast_painted_source[.]tsv$",
    "_source[.]tsv$")
}

lisa_explore_recipe_pattern <- function(product) {
  switch(as.character(product),
    kegg_pathway_map = "_painted_recipe[.]R$",
    contrast_kegg_map = "_contrast_painted_recipe[.]R$",
    "_recipe[.]R$")
}

lisa_explore_product_directory <- function(request) {
  if (identical(request$unit_type, "single_de")) {
    paste0("artifacts/", request$analysis_id, "/", request$collection, "/",
           request$product, "/")
  } else {
    paste0("artifacts/contrasts/", request$contrast_id, "/",
           request$collection, "/", request$product, "/")
  }
}

# Complete-product validation. A figure is one figure delivered as PNG, PDF, its
# source table and its recipe. An attempt that produced only some of them is
# incomplete and must not be promoted.
lisa_explore_assert_complete_product <- function(files, request, policy) {
  paths <- vapply(files, function(file) as.character(file$path), character(1))
  artifacts <- paths[startsWith(paths, "artifacts/")]
  if (!length(artifacts)) {
    stop("LISA-EXPLORE-027 the render produced no artifact files.", call. = FALSE)
  }
  expected <- character()
  if (isTRUE(policy$formats$png)) expected <- c(expected, "[.]png$")
  if (isTRUE(policy$formats$pdf)) expected <- c(expected, "[.]pdf$")
  if (isTRUE(policy$source_data)) {
    expected <- c(expected, lisa_explore_source_data_pattern(request$product))
  }
  if (isTRUE(policy$recipes)) {
    expected <- c(expected, lisa_explore_recipe_pattern(request$product))
  }
  missing <- expected[!vapply(expected, function(pattern) {
    any(grepl(pattern, artifacts))
  }, logical(1))]
  if (length(missing)) {
    stop("LISA-EXPLORE-028 the figure is incomplete; no file matched: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  product_dir <- lisa_explore_product_directory(request)
  drawn <- artifacts[grepl("[.](png|pdf|svg)$", artifacts)]
  foreign <- drawn[!startsWith(drawn, product_dir)]
  if (length(foreign)) {
    stop("LISA-EXPLORE-029 the render produced figures outside the requested ",
         "product: ", paste(utils::head(foreign, 5L), collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}
