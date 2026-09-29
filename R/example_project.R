# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

#' Create the installed lisaR sample project
#'
#' Copies a bounded, synthetic and fully offline sample project from the
#' installed package into a new directory. The destination must not already
#' exist.
#'
#' @param path New directory to create.
#'
#' @return Invisibly, the normalized project directory.
#' @export
#'
#' @examples
#' project <- lisa_init_project(tempfile("lisa-quick-start-"))
#' file.exists(file.path(project, "study.yml"))
lisa_init_project <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(trimws(path))) {
    stop("`path` must name one new project directory.", call. = FALSE)
  }
  destination <- lisa_managed_destination(path, create_parent = TRUE)
  target <- destination$path
  parent <- destination$parent
  if (lisa_path_entry_exists(target)) {
    stop("LISA-INIT-001 destination already exists; lisaR will not overwrite it: ",
         target, call. = FALSE)
  }
  template <- system.file("examples", "quick-start", package = "lisaR")
  if (!nzchar(template) || !dir.exists(template)) {
    stop("LISA-INIT-004 installed sample template is missing; reinstall lisaR.",
         call. = FALSE)
  }
  staging <- file.path(parent, paste0(".", basename(target), ".staging-", Sys.getpid()))
  if (lisa_path_entry_exists(staging)) {
    stop("LISA-INIT-005 staging destination already exists: ", staging,
         call. = FALSE)
  }
  staging <- lisa_run_root(staging)
  promoted <- FALSE
  recovery_required <- FALSE
  on.exit(if (!promoted && !recovery_required &&
              lisa_path_entry_exists(staging)) {
    lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  }, add = TRUE)
  template <- lisa_assert_run_tree_safe(template)
  template_tree <- lisa_scan_run_tree(template)
  if (!nrow(template_tree)) {
    stop("LISA-INIT-006 installed sample project is empty.", call. = FALSE)
  }
  relative <- substring(template_tree$path, nchar(template) + 2L)
  directories <- which(template_tree$isdir)
  for (i in directories) {
    lisa_guarded_dir_create(file.path(staging, relative[[i]]), staging)
  }
  files <- which(!template_tree$isdir)
  tryCatch({
    for (i in files) {
      lisa_guarded_copy(
        template_tree$path[[i]], file.path(staging, relative[[i]]),
        overwrite = FALSE, run_root = staging
      )
    }
  },
    lisa_filesystem_recovery_error = function(error) {
      recovery_required <<- TRUE
      stop(error)
    }
  )
  if (!length(files) ||
      !file.exists(file.path(staging, "study.yml"))) {
    stop("LISA-INIT-006 failed to copy the installed sample project.",
         call. = FALSE)
  }
  lisa_assert_run_tree_safe(staging)
  promotion_error <- NULL
  promoted_ok <- tryCatch({
    lisa_promote_managed_directory(staging, target)
    TRUE
  }, error = function(error) {
    promotion_error <<- error
    if (inherits(error, "lisa_filesystem_recovery_error")) {
      recovery_required <<- TRUE
    }
    FALSE
  })
  if (!promoted_ok) {
    if (inherits(promotion_error, "lisa_filesystem_recovery_error")) {
      stop(promotion_error)
    }
    stop(
      "LISA-INIT-007 failed to promote the sample project atomically: ",
      conditionMessage(promotion_error), call. = FALSE
    )
  }
  promoted <- TRUE
  invisible(normalizePath(target, winslash = "/", mustWork = TRUE))
}
