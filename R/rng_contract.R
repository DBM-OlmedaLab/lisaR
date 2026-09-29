# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Deterministic random-number contract for lisaR execution.

lisa_rng_snapshot <- function() {
  list(
    kind = RNGkind(),
    has_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  )
}

lisa_restore_rng <- function(snapshot) {
  do.call(RNGkind, as.list(snapshot$kind))
  if (isTRUE(snapshot$has_seed)) {
    assign(".Random.seed", snapshot$seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

lisa_preserve_rng <- function(code) {
  snapshot <- lisa_rng_snapshot()
  on.exit(lisa_restore_rng(snapshot), add = TRUE)
  force(code)
}

#' Derive a deterministic effective seed for one task
#'
#' @param root_seed Integer root seed. Defaults to 1729.
#' @param task_id Stable task identifier.
#'
#' @return An integer seed suitable for `set.seed()`.
#' @keywords internal
lisa_task_seed <- function(root_seed = 1729L, task_id) {
  root_seed <- as.integer(root_seed[[1]])
  if (is.na(root_seed)) stop("root_seed must be a finite integer.", call. = FALSE)
  text <- enc2utf8(as.character(task_id[[1]]))
  bytes <- utf8ToInt(text)
  hash <- 0
  for (byte in bytes) hash <- (hash * 131 + byte) %% 2147483646
  as.integer((abs(as.numeric(root_seed)) + hash) %% 2147483646 + 1)
}

#' Evaluate code using a deterministic L'Ecuyer-CMRG task seed
#'
#' @param root_seed Integer root seed.
#' @param task_id Stable task identifier.
#' @param expr Expression to evaluate.
#'
#' @return The value of `expr`.
#' @keywords internal
lisa_with_task_seed <- function(root_seed = 1729L, task_id, expr) {
  expr <- substitute(expr)
  environment <- parent.frame()
  lisa_preserve_rng({
    RNGkind("L'Ecuyer-CMRG")
    set.seed(lisa_task_seed(root_seed, task_id))
    eval(expr, envir = environment)
  })
}

lisa_validate_workers <- function(workers = 4L) {
  if (length(workers) != 1L || is.logical(workers) || is.na(workers) ||
      !is.numeric(workers) || !is.finite(workers) || workers < 1 ||
      workers > .Machine$integer.max || workers != floor(workers)) {
    stop("pipeline.workers must be one positive integer.", call. = FALSE)
  }
  as.integer(workers)
}

lisa_runtime_platform <- function(os_type = .Platform$OS.type,
                                  sysname = unname(Sys.info()[["sysname"]])) {
  os_type <- if (length(os_type)) tolower(as.character(os_type[[1L]])) else ""
  sysname <- if (length(sysname)) tolower(as.character(sysname[[1L]])) else ""
  if (is.na(os_type)) os_type <- ""
  if (is.na(sysname)) sysname <- ""
  if (identical(os_type, "windows") || grepl("windows", sysname, fixed = TRUE)) return("windows")
  if (identical(sysname, "linux")) return("linux")
  if (identical(sysname, "darwin")) return("macos")
  "other"
}

lisa_positive_worker_limit <- function(value) {
  if (!length(value)) return(NA_integer_)
  value <- suppressWarnings(as.numeric(as.character(value[[1L]])))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value < 1 ||
      value > .Machine$integer.max || value != floor(value)) return(NA_integer_)
  as.integer(value)
}

lisa_visible_worker_limit <- function(
  detected_cores = suppressWarnings(parallel::detectCores(logical = TRUE)),
  scheduler = Sys.getenv(
    c("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS", "NCPUS", "LSB_DJOB_NUMPROC"),
    unset = ""
  )
) {
  detected <- lisa_positive_worker_limit(detected_cores)
  scheduler_limits <- vapply(as.list(scheduler), lisa_positive_worker_limit, integer(1))
  scheduler_limits <- scheduler_limits[!is.na(scheduler_limits)]
  candidates <- c(visible_cores = detected, scheduler_allocation = scheduler_limits)
  candidates <- candidates[!is.na(candidates)]
  if (!length(candidates)) return(list(workers = NA_integer_, reason = ""))
  smallest <- min(candidates)
  reason <- if (any(scheduler_limits == smallest)) "scheduler_allocation" else "visible_cores"
  list(workers = as.integer(smallest), reason = reason)
}

lisa_parallel_plan <- function(
  workers,
  tasks,
  platform = lisa_runtime_platform(),
  check_limit = Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""),
  allocation = lisa_visible_worker_limit(),
  backend_requested = "auto"
) {
  workers <- lisa_validate_workers(workers)
  if (length(tasks) != 1L || is.logical(tasks) || is.na(tasks) ||
      !is.numeric(tasks) || !is.finite(tasks) || tasks < 0 ||
      tasks > .Machine$integer.max || tasks != floor(tasks)) {
    stop("tasks must be one non-negative integer.", call. = FALSE)
  }
  tasks <- as.integer(tasks)
  platform <- match.arg(as.character(platform[[1L]]), c("linux", "windows", "macos", "other"))
  backend_requested <- match.arg(as.character(backend_requested[[1L]]), c("auto", "sequential", "multicore"))
  effective <- min(workers, max(1L, tasks))
  reasons <- character()
  if (tasks < workers) reasons <- c(reasons, "task_count")

  check_limit <- tolower(trimws(as.character(check_limit[[1L]])))
  if (check_limit %in% c("true", "yes", "1") && effective > 2L) {
    effective <- 2L
    reasons <- c(reasons, "r_check_limit_cores")
  }

  if (identical(platform, "linux") && !is.null(allocation)) {
    allocation_workers <- lisa_positive_worker_limit(allocation$workers %||% NA_integer_)
    if (!is.na(allocation_workers) && effective > allocation_workers) {
      effective <- allocation_workers
      reasons <- c(reasons, as.character(allocation$reason %||% "visible_allocation"))
    }
  }

  if (identical(backend_requested, "sequential") && effective > 1L) {
    effective <- 1L
    reasons <- c(reasons, "requested_serial")
  }
  unsupported_request <- !identical(platform, "linux") &&
    (workers > 1L || identical(backend_requested, "multicore")) &&
    !identical(backend_requested, "sequential")
  if (unsupported_request) {
    if (effective > 1L) effective <- 1L
    reasons <- c(reasons, "unsupported_platform_serial")
  }
  backend_effective <- if (identical(platform, "linux") && effective > 1L &&
                              !identical(backend_requested, "sequential")) "multicore" else "sequential"

  structure(list(
    platform = platform,
    backend_requested = backend_requested,
    backend_effective = backend_effective,
    workers_requested = workers,
    workers_effective = as.integer(effective),
    task_count = tasks,
    cap_reason = if (length(reasons)) paste(unique(reasons), collapse = ";") else "none",
    inner_threads = 1L
  ), class = "lisa_parallel_plan")
}

lisa_effective_workers <- function(workers, tasks = .Machine$integer.max,
                                   platform = lisa_runtime_platform(),
                                   check_limit = Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""),
                                   allocation = lisa_visible_worker_limit()) {
  lisa_parallel_plan(
    workers = workers, tasks = tasks, platform = platform,
    check_limit = check_limit, allocation = allocation
  )$workers_effective
}

lisa_inner_thread_variables <- function() {
  c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS"
  )
}

lisa_thread_environment <- function() {
  variables <- lisa_inner_thread_variables()
  stats::setNames(Sys.getenv(variables, unset = ""), variables)
}

lisa_with_inner_threads <- function(threads = 1L, code) {
  threads <- lisa_validate_workers(threads)
  variables <- lisa_inner_thread_variables()
  previous <- Sys.getenv(variables, unset = NA_character_)
  on.exit({
    present <- !is.na(previous)
    if (any(present)) do.call(Sys.setenv, as.list(previous[present]))
    if (any(!present)) Sys.unsetenv(variables[!present])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep(as.character(threads), length(variables)), variables)))
  force(code)
}

lisa_task_runtime_context <- function(plan, root_seed, task_id) {
  list(
    task_id = as.character(task_id),
    root_seed = if (is.null(root_seed)) NA_integer_ else as.integer(root_seed),
    task_seed = if (is.null(root_seed)) NA_integer_ else lisa_task_seed(root_seed, task_id),
    backend_requested = plan$backend_requested,
    backend_effective = plan$backend_effective,
    workers_requested = plan$workers_requested,
    workers_effective = plan$workers_effective,
    cap_reason = plan$cap_reason,
    inner_threads = plan$inner_threads,
    platform = plan$platform
  )
}

lisa_with_task_runtime_context <- function(context, code) {
  previous <- getOption("lisaR.outer_parallel_context", NULL)
  options(lisaR.outer_parallel_context = context)
  on.exit(options(lisaR.outer_parallel_context = previous), add = TRUE)
  lisa_with_inner_threads(context$inner_threads, force(code))
}

lisa_task_envelope <- function(task_id, evaluator) {
  tryCatch({
    value <- evaluator(task_id)
    if (is.null(value)) {
      return(list(
        task_id = task_id, ok = FALSE, value = NULL,
        error = list(class = "lisa_null_task_result", message = "task returned NULL")
      ))
    }
    list(task_id = task_id, ok = TRUE, value = value, error = NULL)
  }, error = function(error) {
    list(
      task_id = task_id, ok = FALSE, value = NULL,
      error = list(class = class(error)[[1L]], message = conditionMessage(error))
    )
  })
}

lisa_multicore_apply <- function(task_ids, worker, workers, preschedule = FALSE) {
  parallel::mclapply(
    task_ids, worker, mc.cores = workers, mc.preschedule = preschedule,
    mc.set.seed = FALSE, mc.allow.recursive = FALSE
  )
}

lisa_reconcile_task_envelopes <- function(envelopes, task_ids) {
  if (!is.list(envelopes) || length(envelopes) != length(task_ids)) {
    stop(
      sprintf(
        "LISA-PARALLEL-001 worker result cardinality mismatch: expected %d, received %d.",
        length(task_ids), if (is.list(envelopes)) length(envelopes) else 0L
      ),
      call. = FALSE
    )
  }
  if (is.null(names(envelopes)) || !identical(names(envelopes), task_ids)) {
    stop("LISA-PARALLEL-002 worker result names do not match the planned task IDs.", call. = FALSE)
  }

  required <- c("task_id", "ok", "value", "error")
  malformed <- vapply(seq_along(task_ids), function(index) {
    envelope <- envelopes[[index]]
    !is.list(envelope) || anyDuplicated(names(envelope)) ||
      !setequal(names(envelope), required) ||
      length(envelope$task_id) != 1L || !is.character(envelope$task_id) ||
      !identical(envelope$task_id, task_ids[[index]]) ||
      length(envelope$ok) != 1L || !is.logical(envelope$ok) || is.na(envelope$ok) ||
      (isTRUE(envelope$ok) && (is.null(envelope$value) || !is.null(envelope$error))) ||
      (!isTRUE(envelope$ok) && (
        !is.null(envelope$value) || !is.list(envelope$error) ||
          !identical(sort(names(envelope$error)), c("class", "message")) ||
          length(envelope$error$class) != 1L || !nzchar(as.character(envelope$error$class)) ||
          length(envelope$error$message) != 1L || !nzchar(as.character(envelope$error$message))
      ))
  }, logical(1))
  if (any(malformed)) {
    stop(
      "LISA-PARALLEL-003 malformed worker envelope for task(s): ",
      paste(task_ids[malformed], collapse = ", "), ".",
      call. = FALSE
    )
  }

  failed <- !vapply(envelopes, function(envelope) isTRUE(envelope$ok), logical(1))
  if (any(failed)) {
    details <- vapply(envelopes[failed], function(envelope) {
      sprintf("%s [%s]: %s", envelope$task_id, envelope$error$class, envelope$error$message)
    }, character(1))
    stop("LISA-PARALLEL-004 worker task failed: ", paste(details, collapse = " | "), call. = FALSE)
  }

  values <- lapply(envelopes, `[[`, "value")
  names(values) <- task_ids
  values
}

lisa_map_task_values <- function(task_ids, fun, plan, root_seed = NULL,
                                 preschedule = FALSE) {
  task_ids <- as.character(task_ids)
  if (!is.function(fun)) stop("fun must be a function.", call. = FALSE)
  if (!length(task_ids)) return(stats::setNames(list(), character()))
  if (anyDuplicated(task_ids)) stop("task_ids must be unique.", call. = FALSE)

  named_task_ids <- stats::setNames(task_ids, task_ids)
  worker <- function(task_id) {
    context <- lisa_task_runtime_context(plan, root_seed, task_id)
    lisa_task_envelope(task_id, function(id) {
      lisa_with_task_runtime_context(context, {
        if (is.null(root_seed)) fun(id) else lisa_with_task_seed(root_seed, id, fun(id))
      })
    })
  }
  envelopes <- lisa_preserve_rng({
    if (identical(plan$backend_effective, "multicore")) {
      lisa_multicore_apply(
        named_task_ids, worker, workers = plan$workers_effective,
        preschedule = preschedule
      )
    } else {
      lapply(named_task_ids, worker)
    }
  })
  values <- lisa_reconcile_task_envelopes(envelopes, task_ids)
  values
}

#' Apply a function with deterministic per-task random streams
#'
#' Linux multicore is the only parallel backend. Other platforms execute the
#' same task contract serially.
#'
#' @param task_ids Stable task identifiers.
#' @param fun Function accepting one task identifier.
#' @param root_seed Integer root seed.
#' @param backend Either `"sequential"` or Linux-only `"multicore"`.
#' @param workers Number of requested workers.
#' @param platform Injectable platform family used by contract tests.
#'
#' @return A list named by `task_ids`.
#' @keywords internal
lisa_map_with_task_seeds <- function(task_ids, fun, root_seed = 1729L,
                                     backend = c("sequential", "multicore"), workers = 1L,
                                     platform = lisa_runtime_platform()) {
  backend <- match.arg(backend)
  task_ids <- as.character(task_ids)
  if (!length(task_ids)) return(stats::setNames(list(), character()))
  plan <- lisa_parallel_plan(
    workers = workers, tasks = length(task_ids), platform = platform,
    backend_requested = backend
  )
  lisa_map_task_values(task_ids, fun, plan, root_seed = root_seed)
}

lisa_parallel_phase <- function(task_ids) {
  phase <- sub(":.*$", "", task_ids)
  if (length(unique(phase)) == 1L) unique(phase) else "mixed"
}

lisa_record_parallel_execution <- function(plan, task_ids, root_seed,
                                           run_root = getOption("lisaR.run_root", NULL)) {
  if (is.null(run_root) || !length(task_ids) || !dir.exists(run_root)) return(invisible(FALSE))
  phase <- lisa_parallel_phase(task_ids)
  thread_values <- paste0(lisa_inner_thread_variables(), "=", plan$inner_threads)
  execution <- data.frame(
    phase = phase,
    platform = plan$platform,
    backend_requested = plan$backend_requested,
    backend_effective = plan$backend_effective,
    workers_requested = plan$workers_requested,
    workers_effective = plan$workers_effective,
    task_count = plan$task_count,
    cap_reason = plan$cap_reason,
    root_seed = as.integer(root_seed),
    inner_threads = plan$inner_threads,
    thread_environment = paste(thread_values, collapse = ";"),
    locale = Sys.getlocale(),
    timezone = Sys.timezone(),
    stringsAsFactors = FALSE
  )
  seeds <- data.frame(
    phase = phase,
    task_id = task_ids,
    root_seed = as.integer(root_seed),
    effective_seed = vapply(task_ids, function(task_id) lisa_task_seed(root_seed, task_id), integer(1)),
    stringsAsFactors = FALSE
  )
  execution_path <- file.path(run_root, "parallel_execution.tsv")
  seeds_path <- file.path(run_root, "parallel_task_seeds.tsv")
  if (file.exists(execution_path)) {
    execution <- rbind(read_lisa_tsv(execution_path), execution)
    execution <- execution[!duplicated(execution$phase, fromLast = TRUE), , drop = FALSE]
  }
  if (file.exists(seeds_path)) {
    seeds <- rbind(read_lisa_tsv(seeds_path), seeds)
    seeds <- seeds[!duplicated(seeds$task_id, fromLast = TRUE), , drop = FALSE]
  }
  write_lisa_tsv(execution, execution_path)
  write_lisa_tsv(seeds, seeds_path)
  invisible(TRUE)
}

lisa_pipeline_map <- function(task_ids, fun, workers = 1L, root_seed = 1729L,
                              platform = lisa_runtime_platform()) {
  workers <- lisa_validate_workers(workers)
  task_ids <- as.character(task_ids)
  if (!length(task_ids)) return(stats::setNames(list(), character()))
  plan <- lisa_parallel_plan(
    workers = workers, tasks = length(task_ids), platform = platform,
    backend_requested = "auto"
  )
  values <- lisa_map_task_values(task_ids, fun, plan, root_seed = root_seed)
  lisa_record_parallel_execution(plan, task_ids, root_seed)
  values
}

lisa_find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(20L)) {
    if (file.exists(file.path(current, "DESCRIPTION"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  normalizePath(start, winslash = "/", mustWork = FALSE)
}

#' Capture deterministic execution and environment identity
#'
#' @param root_seed Integer root seed.
#' @param task_ids Stable task identifiers.
#' @param backend Inner execution backend.
#' @param workers Number of inner workers or threads.
#' @param outer_context Outer pipeline task context. Defaults to the context
#'   installed by `lisa_pipeline_map()`.
#' @param project_root Package source root containing lock material.
#'
#' @return A `lisa_runtime_manifest` list.
#' @keywords internal
lisa_runtime_manifest <- function(root_seed = 1729L, task_ids = character(), backend = "sequential",
                                  workers = 1L,
                                  outer_context = getOption("lisaR.outer_parallel_context", NULL),
                                  project_root = lisa_find_project_root()) {
  task_ids <- as.character(task_ids)
  lock_files <- c("renv.lock", "envs/lisa-core.yml", "envs/lisa-full.yml", "envs/LOCK_GENERATION_BLOCKER.md")
  lock_paths <- file.path(project_root, lock_files)
  hashes <- vapply(lock_paths, function(path) if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_, character(1))
  declared <- c("lisaR", "jsonlite", "yaml", "fgsea", "hommel", "openxlsx", "ggplot2", "msigdbr", "DESeq2")
  installed <- vapply(declared, requireNamespace, logical(1), quietly = TRUE)
  versions <- vapply(seq_along(declared), function(index) {
    if (!installed[[index]]) return(NA_character_)
    as.character(utils::packageVersion(declared[[index]]))
  }, character(1))
  if (is.null(outer_context)) {
    outer_context <- list(
      task_id = "",
      root_seed = NA_integer_,
      task_seed = NA_integer_,
      backend_requested = backend,
      backend_effective = backend,
      workers_requested = as.integer(workers),
      workers_effective = as.integer(workers),
      cap_reason = "none",
      inner_threads = as.integer(workers),
      platform = lisa_runtime_platform()
    )
  }
  structure(list(
    root_seed = as.integer(root_seed),
    task_seeds = data.frame(task_id = task_ids, effective_seed = vapply(task_ids, function(task_id) lisa_task_seed(root_seed, task_id), integer(1)), stringsAsFactors = FALSE),
    rng_kind = c(kind = "L'Ecuyer-CMRG", normal.kind = RNGkind()[[2]], sample.kind = RNGkind()[[3]]),
    backend = as.character(outer_context$backend_effective),
    workers = as.integer(outer_context$workers_effective),
    backend_requested = as.character(outer_context$backend_requested),
    backend_effective = as.character(outer_context$backend_effective),
    workers_requested = as.integer(outer_context$workers_requested),
    workers_effective = as.integer(outer_context$workers_effective),
    cap_reason = as.character(outer_context$cap_reason),
    outer_task_id = as.character(outer_context$task_id),
    outer_root_seed = as.integer(outer_context$root_seed),
    outer_task_seed = as.integer(outer_context$task_seed),
    inner_backend = as.character(backend),
    inner_workers = as.integer(workers),
    inner_threads = as.integer(outer_context$inner_threads),
    platform = as.character(outer_context$platform),
    locale = Sys.getlocale(),
    timezone = Sys.timezone(),
    thread_environment = lisa_thread_environment(),
    versions = data.frame(package = c("R", declared), version = c(R.version.string, versions), stringsAsFactors = FALSE),
    lock_hashes = data.frame(file = lock_files, md5 = unname(hashes), stringsAsFactors = FALSE)
  ), class = "lisa_runtime_manifest")
}
