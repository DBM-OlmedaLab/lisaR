# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# The exploration workspace.
#
# The STANDARD run is opened read-only. Everything this engine writes lives in a
# separate root, which is asserted to be outside the source run before any
# directory is created. Inventories and manifests are never written under
# the source tree.
# ---------------------------------------------------------------------------

lisa_explore_schema <- "lisa-explore-workspace/1"

lisa_explore_host <- function() {
  host <- tryCatch(as.character(Sys.info()[["nodename"]]), error = function(e) "")
  if (!length(host) || is.na(host) || !nzchar(host)) "unknown-host" else host
}

lisa_explore_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
}

# --- process identity ------------------------------------------------------
#
# A host+PID pair is not an identity: PIDs are recycled, so a dead owner whose
# number has been reissued to an unrelated process looks alive for ever, and a
# reclaimed lock can be deleted out from under its new owner. We therefore carry
# the process *start time* as well. host+pid+start_time is unique on a host for
# as long as the process exists, which is exactly the lifetime we need.
#
# `start_time` is best-effort: `ps` when available, otherwise field 22 of
# /proc/<pid>/stat (start time in jiffies since boot). When neither is available
# the field is NA and we degrade to the old PID-only test, but we say so instead
# of pretending the check was strong.

lisa_explore_proc_start_time <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (is.na(pid) || pid <= 0L) return(NA_character_)
  if (requireNamespace("ps", quietly = TRUE)) {
    value <- tryCatch(
      format(as.numeric(ps::ps_create_time(ps::ps_handle(pid))), digits = 17),
      error = function(error) NA_character_)
    if (!is.na(value)) return(value)
  }
  stat_file <- file.path("/proc", pid, "stat")
  if (file.exists(stat_file)) {
    return(tryCatch({
      line <- readLines(stat_file, warn = FALSE)[[1L]]
      # The comm field may contain spaces and parentheses; everything after the
      # final ") " is positional, and start time is the 20th of those fields.
      fields <- strsplit(sub("^.*\\) ", "", line), " ", fixed = TRUE)[[1L]]
      as.character(fields[[20L]])
    }, error = function(error) NA_character_))
  }
  NA_character_
}

lisa_explore_process_identity <- function(pid = Sys.getpid()) {
  list(host = lisa_explore_host(), pid = as.integer(pid),
       start_time = lisa_explore_proc_start_time(pid))
}

# Fail-closed liveness. When we cannot judge a process -- a different host, or no
# usable process interface -- we report it as alive. That way we never steal a
# lock from a healthy worker and never declare a running job dead. The cost is
# that a genuinely orphaned job on a foreign host stays `running`; it is never
# promoted to `available`, which is the guarantee that actually matters.
lisa_explore_pid_alive <- function(pid, host) {
  pid <- suppressWarnings(as.integer(pid))
  if (is.na(pid) || pid <= 0L) return(FALSE)
  if (!identical(as.character(host), lisa_explore_host())) return(TRUE)
  if (requireNamespace("ps", quietly = TRUE)) {
    return(tryCatch({
      handle <- ps::ps_handle(pid)
      isTRUE(ps::ps_is_running(handle))
    }, error = function(error) FALSE))
  }
  if (dir.exists("/proc")) return(dir.exists(file.path("/proc", pid)))
  TRUE
}

# Liveness of a recorded identity rather than of a bare number. A PID that is
# alive but whose start time differs from the recorded one has been *reused*: the
# original process is gone, so the identity is dead even though the number is
# busy. This is the check that makes stale-lock reclaim and orphan detection
# correct rather than merely usually-correct.
lisa_explore_identity_alive <- function(identity) {
  if (is.null(identity)) return(FALSE)
  host <- as.character(identity$host)
  pid <- suppressWarnings(as.integer(identity$pid))
  if (is.na(pid) || pid <= 0L) return(FALSE)
  if (!identical(host, lisa_explore_host())) return(TRUE)   # fail-closed
  if (!lisa_explore_pid_alive(pid, host)) return(FALSE)
  recorded <- as.character(identity$start_time)
  if (length(recorded) != 1L || is.na(recorded) || !nzchar(recorded)) {
    # No start time was recorded (older record, or no process interface). We can
    # only answer the weaker PID question, and we already know it is alive.
    return(TRUE)
  }
  current <- lisa_explore_proc_start_time(pid)
  if (is.na(current)) return(TRUE)                          # fail-closed
  identical(current, recorded)
}

lisa_explore_write_json <- function(object, path, run_root) {
  text <- jsonlite::toJSON(object, auto_unbox = TRUE, pretty = TRUE,
                           null = "null", na = "null")
  lisa_guarded_write(path, function(target) {
    writeLines(as.character(text), target, useBytes = TRUE)
  }, run_root = run_root)
  invisible(path)
}

lisa_explore_read_json <- function(path) {
  if (!file.exists(path)) return(NULL)
  jsonlite::fromJSON(paste(readLines(path, warn = FALSE), collapse = "\n"),
                     simplifyVector = TRUE, simplifyDataFrame = FALSE)
}

# --- default workspace location -------------------------------------
#
# The previous default was `<dirname(run)>/<run>_explore`, i.e. a sibling of the
# run *inside the results tree that holds it*. For an accepted run such as the
# Riaz STANDARD run, this places generated and exported material inside the
# accepted reports tree. The source run itself was never written
# to, but the folder boundary was still crossed by default.
#
# The default is now a persistent, user-owned location outside every run:
# `tools::R_user_dir("lisaR", "data")/explore/<run>-<digest>`. It is the standard
# CRAN-policy location for package data, it is not project-private, and it is
# stable across sessions so "close and reopen reuses the artifacts" still holds.
# The digest of the *absolute* run path keeps two identically named runs apart.
#
# When no safe default can be established the workspace is not invented
# somewhere arbitrary: `lisa_explore_open()` refuses and asks for an explicit
# `extension_root`. An explicit root keeps going through `lisa_managed_destination()`
# and the guarded IO layer exactly as before; nothing about it is relaxed here.

lisa_explore_user_data_dir <- function() {
  tryCatch(tools::R_user_dir("lisaR", which = "data"),
           error = function(error) NA_character_)
}

lisa_explore_default_root <- function(run_dir) {
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  base <- lisa_explore_user_data_dir()
  if (length(base) != 1L || is.na(base) || !nzchar(base)) {
    stop("LISA-EXPLORE-016 no safe default workspace location could be ",
         "established on this machine. Pass extension_root = \"<a writable ",
         "directory outside the run>\" to choose where exploration output is ",
         "kept.", call. = FALSE)
  }
  file.path(base, "explore",
            paste0(basename(run_dir), "-",
                   substr(lisa_sha256_text(lisa_path_key(run_dir)), 1L, 10L)))
}

# A default root must be creatable. If it is not -- read-only home, no HOME at
# all, a file in the way -- say so once, with the action that fixes it, rather
# than failing later inside the guarded IO layer with a path the user never chose.
lisa_explore_assert_default_root_usable <- function(root) {
  parent <- dirname(root)
  created <- tryCatch({
    if (!dir.exists(parent)) {
      dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    }
    dir.exists(parent) && file.access(parent, mode = 2L) == 0L
  }, error = function(error) FALSE)
  if (!isTRUE(created)) {
    stop("LISA-EXPLORE-016 the default workspace location is not writable: ",
         parent, ". Pass extension_root = \"<a writable directory outside the ",
         "run>\" to choose where exploration output is kept.", call. = FALSE)
  }
  invisible(root)
}

#' Open an exploration workspace over a completed run
#'
#' Opens a saved STANDARD run read-only and prepares a separate writable
#' workspace outside it. The source run is never written to.
#'
#' @param run_dir A completed lisaR run directory, used read-only.
#' @param extension_root Optional workspace location. It must lie outside
#'   `run_dir`. When omitted, a persistent user-owned directory under
#'   [tools::R_user_dir()] is used, so no generated or exported material is ever
#'   placed inside the results tree that holds an accepted run. Pass an explicit
#'   path to keep the workspace somewhere else.
#' @return An object of class `lisa_explore_workspace`.
#' @export
#'
#' @examples
#' # Opening a workspace requires a completed run; the sample project is
#' # offline, so the prerequisite run is kept interactive.
#' if (interactive()) {
#'   ws <- lisa_explore_open(run_dir)
#'   ws$root
#' }
lisa_explore_open <- function(run_dir, extension_root = NULL) {
  source_run <- lisa_assert_run_tree_safe(run_dir)
  manifest_path <- file.path(source_run, "run_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("LISA-EXPLORE-010 run_dir is not a completed run: run_manifest.tsv is absent.",
         call. = FALSE)
  }
  source_manifest_hash <- lisa_sha256_file(manifest_path)

  requested <- if (is.null(extension_root)) {
    lisa_explore_assert_default_root_usable(lisa_explore_default_root(source_run))
  } else {
    lisa_explore_scalar(extension_root, "extension_root")
  }
  # Refuse before creating anything: the workspace may never live inside the
  # immutable source run.
  probe <- lisa_managed_destination(requested, create_parent = FALSE)$path
  if (lisa_path_within(probe, source_run) ||
      identical(lisa_path_key(probe), lisa_path_key(source_run))) {
    stop("LISA-EXPLORE-011 extension_root must be outside the immutable source run.",
         call. = FALSE)
  }

  root <- lisa_run_root(requested)
  for (child in c("jobs", "prepared", "staging", "extensions", "exports")) {
    lisa_guarded_dir_create(file.path(root, child), run_root = root)
  }

  session_path <- file.path(root, "session.json")
  existing <- lisa_explore_read_json(session_path)
  if (!is.null(existing)) {
    if (!identical(as.character(existing$schema), lisa_explore_schema)) {
      stop("LISA-EXPLORE-012 the workspace was written by an incompatible engine version: ",
           as.character(existing$schema), call. = FALSE)
    }
    if (!identical(as.character(existing$source_run_basename), basename(source_run))) {
      stop("LISA-EXPLORE-013 this workspace belongs to a different run: ",
           as.character(existing$source_run_basename), call. = FALSE)
    }
  } else {
    lisa_explore_write_json(list(
      schema = lisa_explore_schema,
      created_at = lisa_explore_now(),
      source_run_basename = basename(source_run),
      source_manifest_hash = source_manifest_hash,
      lisaR_version = as.character(utils::packageVersion("lisaR")),
      r_version = R.version.string
    ), session_path, run_root = root)
  }

  index_path <- file.path(root, "index.json")
  if (!file.exists(index_path)) {
    lisa_explore_write_json(list(schema = "lisa-explore-index/1", entries = list()),
                            index_path, run_root = root)
  }

  structure(list(
    root = root,
    source_run = source_run,
    source_manifest_hash = source_manifest_hash,
    index_path = index_path,
    session_path = session_path
  ), class = "lisa_explore_workspace")
}

lisa_explore_assert_workspace <- function(ws) {
  if (!inherits(ws, "lisa_explore_workspace")) {
    stop("LISA-EXPLORE-014 a lisa_explore_workspace is required. Open it with ",
         "lisa_explore_open().", call. = FALSE)
  }
  ws
}

#' @export
print.lisa_explore_workspace <- function(x, ...) {
  cat("<lisa_explore_workspace>\n")
  cat("  source run (read-only):", x$source_run, "\n")
  cat("  workspace             :", x$root, "\n")
  invisible(x)
}

# --- single index owner ----------------------------------------------------
#
# A lock directory is created atomically by the filesystem. The owner file names
# the host and PID so a lock left by a crashed process can be reclaimed, while a
# lock held by a live process is respected. A second session over the same
# workspace therefore never becomes a second index writer.

lisa_explore_lock_path <- function(ws) file.path(ws$root, ".index-lock")

# A stale/released sibling is retained only when a transition could not restore
# ownership safely. Its presence is therefore a deliberate fail-closed barrier,
# not litter to ignore: another caller that obtained `.index-lock` in the narrow
# rename window must relinquish it before doing index work. This turns the
# documented three-contender race into an actionable manual-recovery state
# instead of two usable writers.
lisa_explore_lock_quarantines <- function(ws) {
  root <- lisa_existing_run_root(ws$root)
  entries <- list.files(root, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  entries[grepl("^[.]index-lock[.](stale|released)-[0-9a-f]{24}$",
                basename(entries))]
}

# Normal release and reclaim verification happens after the lock directory is
# atomically renamed. These short-lived names are not quarantine: an acquirer
# waits for them to finish and retries. If verification cannot safely restore
# ownership, the transition is renamed to one of the quarantine names above.
lisa_explore_lock_transitions <- function(ws) {
  root <- lisa_existing_run_root(ws$root)
  entries <- list.files(root, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  entries[grepl("^[.]index-lock[.](releasing|reclaiming)-[0-9a-f]{24}$",
                basename(entries))]
}

lisa_explore_lock_transition_owner_verified <- function(evidence) {
  if (!isTRUE(evidence$verified) || is.null(evidence$owner)) return(FALSE)
  owner <- evidence$owner
  scalar_text <- function(value) {
    value <- as.character(value)
    length(value) == 1L && !is.na(value) && nzchar(value)
  }
  pid <- suppressWarnings(as.integer(owner$pid))
  all(vapply(owner[c("host", "start_time", "token")], scalar_text,
             logical(1L))) && length(pid) == 1L && !is.na(pid) && pid > 0L
}

# Complete an interrupted release/reclaim only when the transition still
# carries one verified local owner identity and that exact identity is dead.
# The transition is first moved back to the canonical lock path atomically;
# the normal token-checked stale-lock reclaimer then re-verifies it before any
# evidence is removed. Live, foreign and ambiguous transitions are never moved.
lisa_explore_recover_orphaned_lock_transitions <- function(ws) {
  lock <- lisa_explore_lock_path(ws)
  transitions <- lisa_explore_lock_transitions(ws)
  # More than one transition is conflicting ownership evidence. Even if every
  # recorded PID appears dead, choosing an order would be an unsafe guess.
  if (length(transitions) != 1L) return(invisible(0L))
  recovered <- 0L
  for (transition in transitions) {
    evidence <- lisa_explore_lock_evidence(transition,
                                            allow_marker_only = TRUE)
    if (!lisa_explore_lock_transition_owner_verified(evidence) ||
        lisa_explore_identity_alive(evidence$owner) ||
        lisa_path_entry_exists(lock)) next
    if (!file.rename(transition, lock)) next
    if (isTRUE(lisa_explore_reclaim_stale_lock(
      lock, as.character(evidence$owner$token)))) {
      recovered <- recovered + 1L
    }
  }
  invisible(recovered)
}

lisa_explore_assert_no_lock_quarantine <- function(ws) {
  quarantines <- lisa_explore_lock_quarantines(ws)
  if (!length(quarantines)) return(invisible(TRUE))
  stop("LISA-EXPLORE-075 index-lock ownership evidence is quarantined (",
       paste(basename(utils::head(quarantines, 3L)), collapse = ", "),
       "). No new writer was admitted. Inspect the retained ownership records ",
       "and recover them manually only after every owner is confirmed closed.",
       call. = FALSE)
}

# A token minted per acquisition. It is what makes release *ownership-checked*:
# releasing compares the token on disk with the one we hold, so a caller can
# never delete a lock that some other process acquired in the meantime.
lisa_explore_new_lock_token <- function() {
  substr(lisa_sha256_text(paste(Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%OS6"),
                                stats::runif(1))), 1L, 24L)
}

lisa_explore_lock_marker <- function(lock, token) {
  file.path(lock, paste0(".lisa-owner-", token))
}

# Read only ownership evidence that is strong enough to act on. A guarded write
# interrupted after its payload was flushed but before promotion leaves exactly
# one `.lisa-write-*.json`; that complete record is equivalent ownership
# evidence. Empty, malformed, multiple or otherwise surprising entries are
# deliberately unverifiable and therefore never reclaimed.
lisa_explore_lock_evidence <- function(lock, allow_marker_only = FALSE) {
  invalid <- function(reason) list(verified = FALSE, owner = NULL,
                                   source = "unverified", reason = reason)
  if (!dir.exists(lock)) return(invalid("the lock directory disappeared"))
  entries <- list.files(lock, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  owner_file <- file.path(lock, "owner.json")
  partials <- entries[grepl("^[.]lisa-write-[[:alnum:]]+[.]json$",
                            basename(entries))]
  markers <- entries[grepl("^[.]lisa-owner-[0-9a-f]{24}$", basename(entries))]
  candidate <- if (lisa_path_entry_exists(owner_file)) owner_file else {
    if (length(partials) != 1L) {
      if (isTRUE(allow_marker_only) && length(entries) == 1L &&
          length(markers) == 1L && dir.exists(markers[[1L]])) {
        return(list(verified = TRUE,
                    owner = list(token = sub("^[.]lisa-owner-", "",
                                             basename(markers[[1L]]))),
                    source = "initialization marker", reason = ""))
      }
      reason <- if (!length(partials)) {
        "owner.json is absent and there is no complete partial owner record"
      } else {
        "owner.json is absent and multiple partial owner records are present"
      }
      return(invalid(reason))
    }
    partials[[1L]]
  }
  if (lisa_path_is_link(candidate) || dir.exists(candidate) ||
      !isTRUE(utils::file_test("-f", candidate))) {
    return(invalid("the ownership record is not one regular, non-link file"))
  }
  owner <- tryCatch(lisa_explore_read_json(candidate), error = identity)
  if (inherits(owner, "error") || !is.list(owner)) {
    return(invalid("the ownership record is malformed JSON"))
  }
  required <- c("host", "pid", "start_time", "token", "acquired_at")
  if (!all(required %in% names(owner))) {
    return(invalid(paste0("the ownership record lacks: ",
                          paste(setdiff(required, names(owner)), collapse = ", "))))
  }
  scalar_text <- function(value) {
    value <- as.character(value)
    length(value) == 1L && !is.na(value) && nzchar(value)
  }
  pid <- suppressWarnings(as.integer(owner$pid))
  identity_fields <- c("host", "token", "acquired_at")
  partial <- !identical(candidate, owner_file)
  if (!all(vapply(owner[identity_fields], scalar_text, logical(1L))) ||
      (partial && !scalar_text(owner$start_time)) || length(pid) != 1L ||
      is.na(pid) || pid <= 0L) {
    return(invalid("the ownership record has incomplete identity or token fields"))
  }
  expected_marker <- basename(lisa_explore_lock_marker(lock,
                                                       as.character(owner$token)))
  allowed <- c(basename(candidate), expected_marker)
  unexpected <- setdiff(basename(entries), allowed)
  if (length(unexpected) || length(markers) > 1L ||
      (length(markers) == 1L && !identical(basename(markers), expected_marker))) {
    return(invalid("the lock directory contains ambiguous ownership entries"))
  }
  list(verified = TRUE, owner = owner,
       source = if (!partial) "owner.json" else
         "partial owner record",
       reason = "")
}

# Reclaim is done by *renaming* the lock directory first. file.rename is atomic,
# so if several sessions decide simultaneously that a lock is stale exactly one
# of them succeeds in moving it; the losers see the rename fail and simply retry
# the normal acquisition path. Without this, two reclaimers could both unlink and
# both then create the lock, producing two index owners -- the precise race the
# parent review asked about.
lisa_explore_reclaim_stale_lock <- function(lock, expected_token) {
  # Re-check in place immediately before the rename. Most stale inspectors lose
  # their race here, without ever moving the replacement owner's live lock.
  current <- lisa_explore_lock_evidence(lock)
  if (!isTRUE(current$verified) || is.null(expected_token) ||
      !identical(as.character(current$owner$token),
                 as.character(expected_token))) return(FALSE)
  suffix <- lisa_explore_new_lock_token()
  transition <- paste0(lock, ".reclaiming-", suffix)
  quarantine <- paste0(lock, ".stale-", suffix)
  if (!file.rename(lock, transition)) return(FALSE)
  evidence <- lisa_explore_lock_evidence(transition)
  matches <- isTRUE(evidence$verified) && !is.null(expected_token) &&
    identical(as.character(evidence$owner$token), as.character(expected_token))
  if (!matches) {
    # It changed between inspection and rename. Restore only if the original
    # path is still free; otherwise leave the quarantined directory intact for
    # inspection rather than deleting ownership evidence we cannot verify.
    if (!lisa_path_entry_exists(lock) && file.rename(transition, lock)) {
      return(FALSE)
    }
    if (!file.rename(transition, quarantine) &&
        lisa_path_entry_exists(transition)) {
      # Even if the evidence rename itself is refused, leave an atomic sibling
      # barrier so no new writer can mistake this abnormal state for completion.
      dir.create(quarantine, showWarnings = FALSE, mode = "0700")
    }
    return(FALSE)
  }
  unlink(transition, recursive = TRUE, force = TRUE)
  !lisa_path_entry_exists(transition)
}

lisa_explore_acquire_lock <- function(ws, timeout_seconds = 30) {
  lock <- lisa_explore_lock_path(ws)
  owner_file <- file.path(lock, "owner.json")
  token <- lisa_explore_new_lock_token()
  deadline <- Sys.time() + timeout_seconds
  last_evidence <- NULL
  repeat {
    lisa_explore_assert_no_lock_quarantine(ws)
    lisa_explore_recover_orphaned_lock_transitions(ws)
    transitions <- lisa_explore_lock_transitions(ws)
    if (length(transitions)) {
      transition_evidence <- lapply(
        transitions, lisa_explore_lock_evidence, allow_marker_only = TRUE)
      still_present <- vapply(transitions, lisa_path_entry_exists, logical(1L))
      transitions <- transitions[still_present]
      transition_evidence <- transition_evidence[still_present]
      if (!length(transitions)) next
      owner_verified <- vapply(
        transition_evidence, lisa_explore_lock_transition_owner_verified,
        logical(1L))
      if (any(!owner_verified)) {
        failed <- which(!owner_verified)
        details <- vapply(failed, function(index) paste0(
          basename(transitions[[index]]), ": ",
          if (isTRUE(transition_evidence[[index]]$verified))
            "complete process identity is absent" else
              transition_evidence[[index]]$reason), character(1L))
        stop("LISA-EXPLORE-075 index-lock transition ownership could not be ",
             "verified (", paste(details, collapse = "; "), "). The retained ",
             "evidence was not moved or deleted. Confirm every owner is closed ",
             "before manual recovery.", call. = FALSE)
      }
      if (Sys.time() >= deadline) {
        dead <- vapply(transition_evidence, function(item)
          !lisa_explore_identity_alive(item$owner), logical(1L))
        if (any(dead)) {
          stop("LISA-EXPLORE-075 a verified dead index-lock transition could ",
               "not be reconciled while another lock path was present. All ",
               "ownership evidence was retained; inspect the lock records and ",
               "retry only after the active owner is closed.", call. = FALSE)
        }
        stop("LISA-EXPLORE-015 another session is completing an index-lock ",
             "ownership transition. Its verified owner is live or is on another ",
             "host and cannot be proven closed. Retry after the current writer ",
             "releases it; if this persists, inspect the retained owner record.",
             call. = FALSE)
      }
      Sys.sleep(0.05)
      next
    }
    if (dir.create(lock, showWarnings = FALSE, mode = "0700")) {
      process_identity <- lisa_explore_process_identity()
      marker <- lisa_explore_lock_marker(lock, token)
      if (!dir.create(marker, showWarnings = FALSE, mode = "0700")) {
        unlink(lock, recursive = TRUE, force = TRUE)
        stop("LISA-EXPLORE-074 could not reserve the index-lock ownership marker.",
             call. = FALSE)
      }
      initialization_error <- tryCatch({
        lisa_explore_write_json(c(process_identity, list(token = token,
                                                         acquired_at = lisa_explore_now())),
                                owner_file, run_root = ws$root)
        # The marker is an empty, token-named directory inside the lock. Remove
        # that exact entry directly after re-checking it: guarded recursive
        # deletion first renamed it to an unlabelled `.lisa-delete-*` directory,
        # so a crash between rename and unlink made otherwise complete ownership
        # permanently ambiguous. A crash here now leaves either the expected
        # token marker or no marker, both of which evidence() understands.
        if (lisa_path_is_link(marker) || !dir.exists(marker) ||
            length(list.files(marker, all.files = TRUE, no.. = TRUE))) {
          stop("LISA-EXPLORE-076 the index-lock ownership marker is not the ",
               "expected empty token directory.", call. = FALSE)
        }
        if (!identical(unlink(marker, recursive = TRUE, force = FALSE), 0L) ||
            lisa_path_entry_exists(marker)) {
          stop("LISA-EXPLORE-076 the index-lock ownership marker could not be ",
               "removed safely.", call. = FALSE)
        }
        NULL
      }, error = base::identity)
      if (inherits(initialization_error, "error")) {
        lisa_explore_release_lock(ws, token)
        stop(initialization_error)
      }
      # A competing rename/reclaim may have left retained evidence between our
      # pre-create check and this acquisition. Relinquish our own token before it
      # can be used for index work, then refuse all writers until reconciliation.
      quarantine_error <- tryCatch({
        lisa_explore_assert_no_lock_quarantine(ws)
        NULL
      }, error = base::identity)
      if (inherits(quarantine_error, "error")) {
        lisa_explore_release_lock(ws, token)
        stop(quarantine_error)
      }
      # A release/reclaim can begin after the pre-create check but before this
      # owner is fully initialised. Relinquish this token and retry after the
      # ordinary transition completes; unlike retained quarantine, this is not
      # a manual-recovery error.
      if (length(lisa_explore_lock_transitions(ws))) {
        lisa_explore_release_lock(ws, token)
        if (Sys.time() >= deadline) {
          stop("LISA-EXPLORE-015 another session is completing an index-lock ",
               "ownership transition. Retry after the current writer releases it.",
               call. = FALSE)
        }
        Sys.sleep(0.05)
        next
      }
      return(invisible(token))
    }
    evidence <- lisa_explore_lock_evidence(lock)
    last_evidence <- evidence
    if (isTRUE(evidence$verified) &&
        !lisa_explore_identity_alive(evidence$owner)) {
      # The previous owner is gone -- and gone by identity, not merely by PID
      # number. A complete guarded-write temporary record is accepted here too,
      # closing the interruption window between payload flush and promotion.
      lisa_explore_reclaim_stale_lock(lock, evidence$owner$token)
      next
    }
    if (Sys.time() >= deadline) {
      detail <- if (!isTRUE(last_evidence$verified)) {
        paste0(" Lock recovery was refused because ", last_evidence$reason,
               ". Inspect .index-lock and close only a confirmed owner before retrying.")
      } else if (identical(last_evidence$source, "partial owner record")) {
        " Its partial owner record belongs to a live or unverifiable process; recovery was refused."
      } else ""
      stop("LISA-EXPLORE-015 another session holds the workspace index lock. ",
           "Only one index owner is permitted.", detail, call. = FALSE)
    }
    Sys.sleep(0.05)
  }
}

# Ownership-checked release. We remove the lock only when the token on disk is
# still ours. If it is not, another session legitimately owns the lock now and
# deleting it would hand the index to two writers at once.
lisa_explore_release_lock <- function(ws, token = NULL) {
  lock <- lisa_explore_lock_path(ws)
  if (!dir.exists(lock)) return(invisible(FALSE))
  if (is.null(token)) return(invisible(FALSE))
  # Never move a directory that is not already proven to carry our token. The
  # post-rename verification below is retained to catch replacement after this
  # check; together they prevent release by a caller with the wrong token.
  current <- lisa_explore_lock_evidence(lock, allow_marker_only = TRUE)
  if (!isTRUE(current$verified) ||
      !identical(as.character(current$owner$token),
                 as.character(token))) return(invisible(FALSE))
  suffix <- lisa_explore_new_lock_token()
  transition <- paste0(lock, ".releasing-", suffix)
  quarantine <- paste0(lock, ".released-", suffix)
  if (!file.rename(lock, transition)) return(invisible(FALSE))
  evidence <- lisa_explore_lock_evidence(transition, allow_marker_only = TRUE)
  matches <- isTRUE(evidence$verified) &&
    identical(as.character(evidence$owner$token), as.character(token))
  if (!matches) {
    if (!lisa_path_entry_exists(lock) && file.rename(transition, lock)) {
      return(invisible(FALSE))
    }
    if (!file.rename(transition, quarantine) &&
        lisa_path_entry_exists(transition)) {
      dir.create(quarantine, showWarnings = FALSE, mode = "0700")
    }
    return(invisible(FALSE))
  }
  unlink(transition, recursive = TRUE, force = TRUE)
  invisible(!lisa_path_entry_exists(transition))
}

lisa_explore_with_lock <- function(ws, expr, timeout_seconds = 30) {
  token <- lisa_explore_acquire_lock(ws, timeout_seconds = timeout_seconds)
  on.exit(lisa_explore_release_lock(ws, token), add = TRUE)
  force(expr)
}

# --- index -----------------------------------------------------------------

lisa_explore_read_index <- function(ws) {
  index <- lisa_explore_read_json(ws$index_path)
  if (is.null(index)) return(list())
  entries <- index$entries
  if (is.null(entries)) list() else entries
}

lisa_explore_write_index <- function(ws, entries) {
  lisa_explore_write_json(list(schema = "lisa-explore-index/1",
                               updated_at = lisa_explore_now(),
                               entries = unname(entries)),
                          ws$index_path, run_root = ws$root)
}

lisa_explore_index_find <- function(entries, key) {
  for (entry in entries) {
    if (identical(as.character(entry$key), as.character(key))) return(entry)
  }
  NULL
}

# The category views an artifact is reachable from. For a native KEGG pathway map
# one artifact is legitimately associated with several categories -- the map's
# identity is the pathway and its paint context, and the category is navigation
# (H2_SCOPE). Kept sorted and deduplicated so the list is a set, not a log.
lisa_explore_union_associations <- function(...) {
  values <- as.character(unlist(list(...), use.names = FALSE))
  values <- values[!is.na(values) & nzchar(values)]
  as.list(sort(unique(values)))
}

lisa_explore_index_upsert <- function(entries, entry) {
  replaced <- FALSE
  out <- lapply(entries, function(existing) {
    if (identical(as.character(existing$key), as.character(entry$key))) {
      replaced <<- TRUE
      # Associations are carried forward here rather than at each call site: the
      # worker rewrites the whole entry when it promotes a figure, so an
      # association recorded at submission time would otherwise be lost exactly
      # when the artifact becomes available.
      entry$associations <- lisa_explore_union_associations(
        existing$associations, entry$associations)
      return(entry)
    }
    existing
  })
  if (!replaced) out <- c(out, list(entry))
  out
}

# --- jobs ------------------------------------------------------------------

lisa_explore_job_path <- function(ws, job_id) {
  file.path(ws$root, "jobs", paste0(job_id, ".json"))
}

lisa_explore_read_job <- function(ws, job_id) {
  lisa_explore_read_json(lisa_explore_job_path(ws, job_id))
}

lisa_explore_write_job <- function(ws, job) {
  lisa_explore_write_json(job, lisa_explore_job_path(ws, job$job_id),
                          run_root = ws$root)
}

lisa_explore_all_jobs <- function(ws) {
  files <- list.files(file.path(ws$root, "jobs"), pattern = "[.]json$",
                      full.names = TRUE)
  # `<job>.launcher.json` records the child PID; it is not a job record.
  files <- files[!grepl("[.]launcher[.]json$", files)]
  jobs <- lapply(sort(files), lisa_explore_read_json)
  Filter(Negate(is.null), jobs)
}
