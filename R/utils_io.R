# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

read_lisa_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) {
      stop("Required TSV file does not exist: ", path, call. = FALSE)
    }
    return(data.frame())
  }
  if (file.info(path)$size == 0) {
    return(data.frame())
  }
  first_lines <- readLines(path, n = 2L, warn = FALSE)
  if (length(first_lines) == 0 || !any(nzchar(trimws(first_lines)))) {
    return(data.frame())
  }
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_lisa_tsv <- function(x, path) {
  lisa_guarded_write(path, function(target) utils::write.table(x, file = target, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""))
  invisible(path)
}

lisa_require_columns <- function(x, columns, label = deparse(substitute(x))) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lisa_is_absolute_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(as.character(path))) {
    return(FALSE)
  }
  path <- gsub("\\\\", "/", path.expand(as.character(path)))
  if (startsWith(path, "//?/") || startsWith(path, "//./")) {
    stop("Windows device-path namespaces are not supported.", call. = FALSE)
  }
  if (grepl("^[A-Za-z]:(?:$|[^/])", path, perl = TRUE)) {
    stop("Windows drive-relative paths are not supported.", call. = FALSE)
  }
  # Three or more leading slashes still denote the POSIX filesystem root.
  # Test this before the UNC branch so such paths are never rebased as if they
  # were relative input.
  if (startsWith(path, "///")) {
    return(TRUE)
  }
  if (startsWith(path, "//")) {
    if (!grepl("^//[^/]+/[^/]+(?:/|$)", path, perl = TRUE)) {
      stop("Incomplete UNC paths are not supported.", call. = FALSE)
    }
    return(TRUE)
  }
  startsWith(path, "/") || grepl("^[A-Za-z]:/", path)
}

lisa_rscript_executable <- function() {
  suffix <- if (.Platform$OS.type == "windows") ".exe" else ""
  path <- file.path(R.home("bin"), paste0("Rscript", suffix))
  if (!file.exists(path) || dir.exists(path)) {
    stop("The Rscript executable is unavailable.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

lisa_norm_path <- function(path, base_dir = getwd()) {
  if (is.na(path) || !nzchar(path)) {
    return(path)
  }
  if (lisa_is_absolute_path(path)) {
    return(path)
  }
  file.path(base_dir, path)
}

lisa_empty_df <- function(columns) {
  stats::setNames(
    as.data.frame(rep(list(character()), length(columns)), stringsAsFactors = FALSE),
    columns
  )
}
