#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

report_route_label <- "standard"
report_renderer_sha256 <- NA_character_
report_project_dir <- NULL
report_gallery_assets <- list()

parse_args <- function(args) {
  out <- list(
    project_dir = NA_character_,
    title = "LISA scientific report", presentation_config = "",
    refresh_category_links = "false", complete_missing = "false", kegg_maps = "false", kegg_cache_root = "",
    kegg_snapshot_id = "", kegg_access_mode = "cache_only", report_mode = "",
    lisa_internal_renderer_sha256 = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Usage: build_LISA_report.R --project-dir PATH [--title TITLE] [--report-mode standard|full]", call. = FALSE)
    }
    j <- i + 1
    while (j < length(args) && !startsWith(args[[j + 1]], "--")) {
      j <- j + 1
    }
    val <- paste(args[(i + 1):j], collapse = " ")
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    out[[key]] <- val
    i <- j + 1
  }
  if (is.na(out$project_dir)) stop("--project-dir is required", call. = FALSE)
  out
}

report_executable_package_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Cannot identify the verified build_LISA_report.R executable.", call. = FALSE)
  }
  script <- normalizePath(
    sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE
  )
  script_dir <- dirname(script)
  package_dir <- if (identical(basename(dirname(script_dir)), "inst")) {
    dirname(dirname(script_dir))
  } else {
    dirname(script_dir)
  }
  package_dir <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(package_dir, "DESCRIPTION"))) {
    stop("Verified report executable is outside a lisaR package tree.", call. = FALSE)
  }
  package_dir
}

esc <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  project_dir <- get0("report_project_dir", inherits = TRUE, ifnotfound = NULL)
  if (!is.null(project_dir) && length(project_dir) == 1L &&
      !is.na(project_dir) && nzchar(project_dir)) {
    x <- report_portable_text(x, project_dir)
  }
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

slug <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  x <- gsub("(^-|-$)", "", x)
  ifelse(nzchar(x), x, "section")
}

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.frame())
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) == 0 || !nzchar(trimws(first_line[[1]]))) return(data.frame())
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

# Native KEGG maps are rendered inside the run whenever run_kegg_maps was
# requested (STANDARD or FULL); their layer is shown in both report modes.
report_kegg_maps_requested <- function(project_dir) {
  contract <- read_tsv(file.path(project_dir, "contract_manifest.tsv"))
  if (!nrow(contract) || !all(c("key", "value") %in% names(contract))) return(FALSE)
  any(contract$key == "run_kegg_maps" & tolower(contract$value) %in% c("true", "yes", "1"))
}

# D5: an unavailable KEGG snapshot is announced on the report cover, never
# silently omitted.
report_kegg_maps_notice <- function(project_dir) {
  status <- read_tsv(file.path(project_dir, "kegg_maps_status.tsv"))
  if (!nrow(status) || !"status" %in% names(status)) return("")
  problems <- status[status$status != "completed", , drop = FALSE]
  if (!nrow(problems)) return("")
  message <- unique(as.character(problems$message))
  where <- if (all(c("kind", "owner", "collection") %in% names(problems)) &&
      any(nzchar(problems$owner))) paste0(" Affected: ", paste(unique(paste(problems$owner[nzchar(problems$owner)],
      problems$collection[nzchar(problems$owner)])), collapse = "; "), ".") else ""
  notice_html("warn", "KEGG maps unavailable",
    paste0(paste(utils::head(message, 3L), collapse = " "), where))
}

report_requested <- function(product, default = TRUE) {
  policy <- get0("report_output_policy", inherits = TRUE, ifnotfound = NULL)
  if (is.null(policy) || !is.data.frame(policy) || !all(c("product", "requested") %in% names(policy))) {
    return(isTRUE(default))
  }
  row <- policy[as.character(policy$product) == product, , drop = FALSE]
  if (nrow(row) != 1L) return(isTRUE(default))
  value <- row$requested[[1]]
  if (is.logical(value)) return(isTRUE(value))
  tolower(as.character(value)) %in% c("true", "t", "1", "yes")
}

rbind_fill_local <- function(xs) {
  xs <- xs[!vapply(xs, function(x) is.null(x) || nrow(x) == 0, logical(1))]
  if (length(xs) == 0) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(cols, names(x))
    for (col in missing) x[[col]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

report_path_components <- function(path) {
  raw_path <- report_slash_path(path.expand(as.character(path)))
  if (!report_is_native_absolute_path(
    raw_path, reject_drive_relative = TRUE
  )) {
    stop("Cannot resolve an absolute report path.", call. = FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (startsWith(path, "//")) {
    components <- strsplit(sub("^//+", "", path), "/", fixed = TRUE)[[1L]]
    components <- components[nzchar(components)]
    if (length(components) < 2L) {
      stop("Cannot resolve an incomplete UNC path.", call. = FALSE)
    }
    return(list(
      root = paste0("//", components[[1L]], "/", components[[2L]]),
      parts = components[-c(1L, 2L)],
      case_insensitive = TRUE
    ))
  }
  if (grepl("^[A-Za-z]:($|/)", path)) {
    return(list(
      root = substr(path, 1L, 2L),
      parts = Filter(nzchar, strsplit(sub("^[A-Za-z]:/?", "", path), "/",
                                     fixed = TRUE)[[1L]]),
      case_insensitive = TRUE
    ))
  }
  if (startsWith(path, "/")) {
    return(list(
      root = "/",
      parts = Filter(nzchar, strsplit(sub("^/+", "", path), "/",
                                     fixed = TRUE)[[1L]]),
      case_insensitive = identical(.Platform$OS.type, "windows")
    ))
  }
  stop("Cannot resolve an absolute report path.", call. = FALSE)
}

report_relative_path <- function(path, from_dir) {
  target <- report_path_components(path)
  base <- report_path_components(from_dir)
  case_insensitive <- isTRUE(target$case_insensitive) ||
    isTRUE(base$case_insensitive)
  compare <- function(x) if (case_insensitive) tolower(x) else x
  if (!identical(compare(target$root), compare(base$root))) {
    stop("Cannot form a relative report path across filesystem roots.",
         call. = FALSE)
  }
  common <- 0L
  limit <- min(length(target$parts), length(base$parts))
  if (limit > 0L) {
    for (index in seq_len(limit)) {
      if (!identical(compare(target$parts[[index]]),
                     compare(base$parts[[index]]))) break
      common <- index
    }
  }
  upward <- rep("..", length(base$parts) - common)
  downward <- if (common < length(target$parts)) {
    target$parts[seq.int(common + 1L, length(target$parts))]
  } else {
    character()
  }
  relative <- c(upward, downward)
  if (!length(relative)) "." else paste(relative, collapse = "/")
}

rel_path <- function(path, from_file) {
  esc(report_relative_path(path, dirname(from_file)))
}

report_root_for_page <- function(page_file) {
  page_dir <- dirname(normalizePath(page_file, mustWork = FALSE))
  if (basename(page_dir) == "report_pages") return(dirname(page_dir))
  page_dir
}

short_media_copy <- function(path, page_file) {
  if (!file.exists(path)) return("")
  root <- report_root_for_page(page_file)
  media_dir <- file.path(root, "report_media")
  dir.create(media_dir, recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(path))
  hash <- tryCatch(unname(tools::md5sum(path)), error = function(e) "")
  if (!nzchar(hash)) hash <- slug(basename(path))
  dest <- file.path(media_dir, paste0(substr(hash, 1, 16), ".", ext))
  if (!file.exists(dest)) file.copy(path, dest, overwrite = TRUE)
  # One canonical media copy is enough. Pages below report_pages/ link back to
  # ../report_media instead of creating a second hardlink/copy subtree.
  rel_path(dest, page_file)
}

gzip_copy <- function(source, destination, compression = 9L) {
  input <- file(source, open = "rb")
  output <- gzfile(destination, open = "wb", compression = compression)
  ok <- FALSE
  on.exit({
    close(input)
    close(output)
    if (!ok && file.exists(destination)) unlink(destination)
  }, add = TRUE)
  repeat {
    block <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
  ok <- TRUE
  invisible(destination)
}

# The complete run is the private provenance record.  These helpers construct
# only the portable *view* used by the HTML report.  They never rewrite the
# canonical source table: internal paths become run-relative and external
# paths lose their private parent while keeping the basename.  Other values
# remain byte-for-byte equivalent when a delimited table is parsed.
report_slash_path <- function(x) gsub("\\\\", "/", enc2utf8(as.character(x)))

report_is_file_uri <- function(x) {
  grepl(
    "^file:(?://|/|[A-Za-z]:/)", trimws(as.character(x)),
    ignore.case = TRUE, perl = TRUE
  )
}

report_file_uri_path <- function(x) {
  path <- report_slash_path(trimws(as.character(x)))
  path <- sub("^file:", "", path, ignore.case = TRUE)
  path <- utils::URLdecode(path)
  if (grepl("^///[A-Za-z]:/", path)) {
    path <- substring(path, 4L)
  } else if (startsWith(path, "///")) {
    # A POSIX file URI has three slashes: two delimit the empty authority and
    # the third is the filesystem root. Preserve that root slash.
    path <- substring(path, 3L)
  }
  path
}

report_is_native_absolute_path <- function(x, reject_drive_relative = FALSE) {
  x <- report_slash_path(trimws(as.character(x)))
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(FALSE)
  if (startsWith(x, "//?/") || startsWith(x, "//./")) {
    stop(
      "LISA-REPORT-PATH-001 Windows device-path namespaces are not supported.",
      call. = FALSE
    )
  }
  if (isTRUE(reject_drive_relative) &&
      grepl("^[A-Za-z]:(?:$|[^/])", x, perl = TRUE)) {
    stop(
      "LISA-REPORT-PATH-002 Windows drive-relative paths are not supported.",
      call. = FALSE
    )
  }
  if (startsWith(x, "///")) return(TRUE)
  if (startsWith(x, "//")) {
    if (!grepl("^//[^/]+/[^/]+(?:/|$)", x, perl = TRUE)) {
      stop("LISA-REPORT-PATH-003 incomplete UNC path.", call. = FALSE)
    }
    return(TRUE)
  }
  startsWith(x, "/") || grepl("^[A-Za-z]:/", x)
}

report_is_absolute_path <- function(x) {
  x <- report_slash_path(trimws(as.character(x)))
  if (report_is_file_uri(x)) x <- report_file_uri_path(x)
  report_is_native_absolute_path(x)
}

report_is_path_column <- function(name) {
  grepl(
    "(^|_)(path|paths|file|files|dir|directory|root|folder|subdir|png|svg|pdf|html|tsv|csv|rds|xlsx|recipe)(_|$)",
    tolower(as.character(name)), perl = TRUE
  )
}

report_path_basename <- function(path) {
  parts <- strsplit(sub("/+$", "", report_slash_path(path)), "/",
                    fixed = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  if (length(parts)) parts[[length(parts)]] else "external-path"
}

report_normalize_existing_ancestor <- function(path) {
  path <- report_slash_path(path.expand(as.character(path)))
  invisible(report_is_native_absolute_path(
    path, reject_drive_relative = TRUE
  ))
  native_absolute <- startsWith(path, "/") ||
    (identical(.Platform$OS.type, "windows") &&
       grepl("^[A-Za-z]:/|^//", path))
  if (!native_absolute) return(path)

  suffix <- character()
  cursor <- path
  repeat {
    if (file.exists(cursor) || dir.exists(cursor)) {
      cursor <- normalizePath(cursor, winslash = "/", mustWork = TRUE)
      break
    }
    parent <- report_slash_path(dirname(cursor))
    if (identical(parent, cursor)) {
      cursor <- normalizePath(cursor, winslash = "/", mustWork = FALSE)
      break
    }
    suffix <- c(basename(cursor), suffix)
    cursor <- parent
  }
  if (length(suffix)) {
    cursor <- do.call(file.path, as.list(c(cursor, suffix)))
  }
  report_slash_path(cursor)
}

report_path_within_root <- function(path, project_dir) {
  path <- report_slash_path(path)
  project_dir <- report_normalize_existing_ancestor(project_dir)
  candidate <- if (grepl("^[A-Za-z]:/|^/", path)) {
    report_normalize_existing_ancestor(path)
  } else {
    path
  }
  insensitive <- identical(.Platform$OS.type, "windows") ||
    grepl("^[A-Za-z]:/", candidate) || grepl("^[A-Za-z]:/", project_dir)
  compare <- function(x) if (insensitive) tolower(x) else x
  candidate_cmp <- compare(candidate)
  root_cmp <- compare(sub("/+$", "", project_dir))
  if (identical(candidate_cmp, root_cmp)) return(".")
  prefix <- paste0(root_cmp, "/")
  if (startsWith(candidate_cmp, prefix)) {
    return(substr(candidate, nchar(prefix) + 1L, nchar(candidate)))
  }
  NA_character_
}

report_portable_path <- function(path, project_dir) {
  path <- report_slash_path(trimws(as.character(path)))
  if (!nzchar(path)) return(path)
  if (report_is_absolute_path(path)) {
    relative <- report_path_within_root(path, project_dir)
    if (!is.na(relative)) return(relative)
    return(report_path_basename(path))
  }
  path <- sub("^[.]/", "", path)
  if (identical(path, "..") || startsWith(path, "../")) {
    return(report_path_basename(path))
  }
  path
}

report_absolute_token_pattern <- paste0(
  "(?<![A-Za-z0-9:./<])",
  "(?:[A-Za-z]:/|//|/)",
  "[^\\t\\r\\n <>\\\"']+"
)

report_file_uri_pattern <- paste0(
  "(?i)(?<![A-Za-z0-9+.-])file:(?://|/|[A-Za-z]:/)",
  "[^\\t\\r\\n <>\\\"']+"
)

report_portable_text_one <- function(value, project_dir, force_path = FALSE,
                                     normalized_root = NULL) {
  if (length(value) == 0L || is.na(value)) return(value)
  value <- enc2utf8(as.character(value))
  if (!nzchar(value)) return(value)
  slash <- report_slash_path(value)
  if (isTRUE(force_path)) {
    invisible(report_is_native_absolute_path(
      trimws(slash), reject_drive_relative = TRUE
    ))
  }
  # This file also runs as a standalone Rscript without attached package helpers.
  root <- if (is.null(normalized_root)) {
    normalizePath(project_dir, winslash = "/", mustWork = FALSE)
  } else normalized_root
  has_root <- grepl(root, slash, fixed = TRUE)
  has_file_uri <- grepl(report_file_uri_pattern, slash, perl = TRUE)
  has_absolute <- report_is_absolute_path(trimws(slash)) ||
    grepl(report_absolute_token_pattern, slash, perl = TRUE)
  if (!force_path && !has_root && !has_absolute && !has_file_uri) return(value)
  if (report_is_file_uri(trimws(slash)) && identical(trimws(slash), slash)) {
    return(report_portable_path(report_file_uri_path(slash), project_dir))
  }
  if (report_is_absolute_path(trimws(slash)) &&
      identical(trimws(slash), slash)) {
    return(report_portable_path(slash, project_dir))
  }
  if (has_root) slash <- gsub(root, ".", slash, fixed = TRUE)
  file_uri_matches <- gregexpr(report_file_uri_pattern, slash, perl = TRUE)
  file_uri_tokens <- regmatches(slash, file_uri_matches)
  file_uri_replacements <- lapply(file_uri_tokens, function(group) {
    vapply(group, function(token) {
      report_portable_path(report_file_uri_path(token), project_dir)
    }, character(1))
  })
  if (length(file_uri_tokens) && any(lengths(file_uri_tokens))) {
    regmatches(slash, file_uri_matches) <- file_uri_replacements
  }
  matches <- gregexpr(report_absolute_token_pattern, slash, perl = TRUE)
  tokens <- regmatches(slash, matches)
  replacements <- lapply(tokens, function(group) {
    vapply(group, report_portable_path, character(1), project_dir = project_dir)
  })
  if (length(tokens) && any(lengths(tokens))) {
    regmatches(slash, matches) <- replacements
  }
  if (force_path) report_portable_path(slash, project_dir) else slash
}

report_portable_text <- function(value, project_dir, force_path = FALSE) {
  # Vector screening avoids filesystem/path work for every numeric cell in an
  # expression matrix. Read-as-text downloads retain their original precision.
  value <- as.character(value)
  if (!length(value)) return(value)
  root <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)
  slash <- report_slash_path(value)
  candidate <- !is.na(value) & nzchar(value) & (
    force_path | grepl(root, slash, fixed = TRUE) |
    grepl(report_file_uri_pattern, slash, perl = TRUE) |
    grepl(report_absolute_token_pattern, slash, perl = TRUE) |
    grepl("^[[:space:]]*(/|[A-Za-z]:/)", slash)
  )
  if (any(candidate)) {
    distinct <- unique(value[candidate])
    converted <- vapply(distinct, report_portable_text_one, character(1),
      project_dir = project_dir, force_path = force_path,
      normalized_root = root, USE.NAMES = FALSE)
    value[candidate] <- converted[match(value[candidate], distinct)]
  }
  value
}

report_portable_data_frame <- function(df, project_dir) {
  if (is.null(df) || !is.data.frame(df) || !length(names(df))) return(df)
  for (index in seq_along(df)) {
    # Numeric/logical columns cannot contain local paths. Preserve their types,
    # including NA/NaN/Inf, instead of converting them cell by cell.
    if (is.numeric(df[[index]]) || is.logical(df[[index]])) next
    original <- as.character(df[[index]])
    portable <- report_portable_text(
      original, project_dir,
      force_path = report_is_path_column(names(df)[[index]])
    )
    changed <- (is.na(original) != is.na(portable)) |
      (!is.na(original) & !is.na(portable) & original != portable)
    if (any(changed)) df[[index]] <- portable
  }
  df
}

report_read_delimited_text <- function(path) {
  lower <- tolower(path)
  ext <- if (grepl("[.]csv([.]gz)?$", lower)) "csv" else "tsv"
  separator <- if (identical(ext, "csv")) "," else "\t"
  quote <- if (identical(ext, "csv")) "\"" else ""
  utils::read.table(
    path, sep = separator, header = TRUE, quote = quote,
    comment.char = "", check.names = FALSE, colClasses = "character",
    na.strings = NULL, fill = TRUE, stringsAsFactors = FALSE
  )
}

report_write_delimited_gzip <- function(df, destination, extension) {
  output <- gzfile(destination, open = "wt", compression = 9L)
  ok <- FALSE
  on.exit({
    close(output)
    if (!ok && file.exists(destination)) unlink(destination)
  }, add = TRUE)
  utils::write.table(
    df, output,
    sep = if (identical(extension, "csv")) "," else "\t",
    quote = identical(extension, "csv"), row.names = FALSE, na = ""
  )
  ok <- TRUE
  invisible(destination)
}

report_portable_table_copy <- function(source, destination, project_dir) {
  extension <- tolower(tools::file_ext(source))
  if (!extension %in% c("tsv", "csv")) {
    stop("LISA-REPORT-PRIVACY-001 portable downloads accept TSV/CSV only.",
         call. = FALSE)
  }
  first_line <- readLines(source, n = 1L, warn = FALSE)
  if (!length(first_line) || !nzchar(trimws(first_line[[1L]]))) {
    gzip_copy(source, destination)
    return(invisible(destination))
  }
  table <- tryCatch(
    report_read_delimited_text(source),
    error = function(error) stop(
      "LISA-REPORT-PRIVACY-002 cannot inspect portable table ",
      basename(source), ": ", conditionMessage(error), call. = FALSE
    )
  )
  portable <- report_portable_data_frame(table, project_dir)
  if (identical(table, portable)) {
    gzip_copy(source, destination)
  } else {
    report_write_delimited_gzip(portable, destination, extension)
  }
  invisible(destination)
}

# Only transformations performed by this invocation may be reused. A previous
# report's unsanitized copy is never trusted merely because its name matches.
report_download_cache <- new.env(parent = emptyenv())

short_file_copy <- function(path, page_file) {
  if (!file.exists(path)) return("")
  root <- report_root_for_page(page_file)
  files_dir <- file.path(root, "report_files")
  dir.create(files_dir, recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("tsv", "csv")) {
    stop("LISA-REPORT-PRIVACY-003 only TSV/CSV tables enter report_files.",
         call. = FALSE)
  }
  hash <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  dest_ext <- paste0(ext, ".gz")
  dest <- file.path(files_dir, paste0(hash, ".", dest_ext))
  cache_key <- paste(normalizePath(root, winslash = "/", mustWork = TRUE),
                     hash, ext, "portable-v2", sep = "|")
  prior <- report_download_cache[[cache_key]]
  valid <- !is.null(prior) && identical(prior$path, dest) && file.exists(dest) &&
    identical(as.numeric(file.info(dest)$size), prior$bytes) &&
    identical(digest::digest(file = dest, algo = "sha256", serialize = FALSE), prior$sha256)
  if (!valid) {
    report_portable_table_copy(path, dest, root)
    report_download_cache[[cache_key]] <- list(path = dest,
      bytes = as.numeric(file.info(dest)$size),
      sha256 = digest::digest(file = dest, algo = "sha256", serialize = FALSE))
  }
  rel_path(dest, page_file)
}

# Global native indexes are downloads, not category-specific products. Use the
# same portable table surface as every other report table, rather than linking
# an unlisted original under outputs/.
report_unclassified_product_href <- function(path, page_file) {
  if (tolower(tools::file_ext(path)) %in% c("tsv", "csv"))
    short_file_copy(path, page_file)
  else rel_path(path, page_file)
}

report_unclassified_product_link <- function(path, page_file) {
  if (tolower(tools::file_ext(path)) %in% c("txt", "md")) {
    # Native README metadata is readable in-place, without introducing a
    # dangling download or leaking original filesystem paths.
    text <- report_portable_text(readLines(path, warn = FALSE), report_root_for_page(page_file))
    return(paste0('<details><summary>', esc(basename(path)),
      '</summary><pre>', esc(paste(text, collapse = "\n")), '</pre></details>'))
  }
  paste0('<a class="file-link" download href="',
    esc(report_unclassified_product_href(path, page_file)), '">', esc(basename(path)), '</a>')
}

report_shareable_manifest_name <- "shareable_report_manifest.tsv"

report_shareable_state <- function(path) {
  state <- tryCatch(
    fs::file_info(path, follow = FALSE),
    error = function(error) {
      stop(
        "LISA-REPORT-SHARE-001 could not inspect a shareable component: ",
        conditionMessage(error), call. = FALSE
      )
    }
  )
  if (nrow(state) != 1L) {
    stop(
      "LISA-REPORT-SHARE-001 shareable-component inspection was inconclusive.",
      call. = FALSE
    )
  }
  type <- as.character(state$type[[1L]])
  if (identical(type, "symlink")) {
    stop("LISA-REPORT-SHARE-001 shareable components cannot be symlinks.",
         call. = FALSE)
  }
  if (!is.na(type) && !type %in% c("file", "directory")) {
    stop(
      "LISA-REPORT-SHARE-011 shareable components must be regular files or directories.",
      call. = FALSE
    )
  }
  list(type = if (is.na(type)) "missing" else type, info = state)
}

report_assert_shareable_not_link <- function(path) {
  report_shareable_state(path)
  invisible(path)
}

report_shareable_lexical_parts <- function(path, project_dir) {
  candidate <- report_slash_path(as.character(path))
  root <- sub("/+$", "", report_slash_path(as.character(project_dir)))
  insensitive <- identical(.Platform$OS.type, "windows") ||
    grepl("^[A-Za-z]:/", candidate) || grepl("^[A-Za-z]:/", root)
  compare <- function(value) if (insensitive) tolower(value) else value
  candidate_cmp <- compare(candidate)
  root_cmp <- compare(root)
  prefix <- paste0(root_cmp, "/")
  if (!startsWith(candidate_cmp, prefix)) {
    stop("LISA-REPORT-SHARE-002 shareable component escaped the run root.",
         call. = FALSE)
  }
  relative <- substring(candidate, nchar(root) + 2L)
  raw_parts <- strsplit(relative, "/", fixed = TRUE)[[1L]]
  parts <- character()
  for (part in raw_parts) {
    if (!nzchar(part) || identical(part, ".")) next
    if (identical(part, "..")) {
      if (!length(parts)) {
        stop("LISA-REPORT-SHARE-002 shareable component escaped the run root.",
             call. = FALSE)
      }
      parts <- parts[-length(parts)]
    } else {
      parts <- c(parts, part)
    }
  }
  if (!length(parts)) {
    stop("LISA-REPORT-SHARE-002 invalid lexical shareable component.",
         call. = FALSE)
  }
  parts
}

report_assert_shareable_file <- function(path, project_dir) {
  root_state <- report_shareable_state(project_dir)
  if (!identical(root_state$type, "directory")) {
    stop("LISA-REPORT-SHARE-011 shareable root is not a directory.",
         call. = FALSE)
  }
  parts <- report_shareable_lexical_parts(path, project_dir)
  current <- project_dir
  final_state <- NULL
  for (index in seq_along(parts)) {
    current <- file.path(current, parts[[index]])
    state <- report_shareable_state(current)
    expected <- if (identical(index, length(parts))) "file" else "directory"
    if (!identical(state$type, expected)) {
      stop(
        "LISA-REPORT-SHARE-011 shareable path changed or has the wrong type: ",
        current, call. = FALSE
      )
    }
    final_state <- state$info
  }
  final_state
}

report_shareable_file_receipt <- function(path, project_dir) {
  before <- report_assert_shareable_file(path, project_dir)
  digest <- lisaR:::lisa_sha256_file(path)
  after <- report_assert_shareable_file(path, project_dir)
  stable_fields <- intersect(
    c("size", "modification_time", "change_time", "device_id", "inode"),
    intersect(names(before), names(after))
  )
  stable <- vapply(stable_fields, function(field) {
    identical(as.character(before[[field]][[1L]]),
              as.character(after[[field]][[1L]]))
  }, logical(1))
  if (length(stable) && !all(stable)) {
    stop("LISA-REPORT-SHARE-012 shareable file changed while it was hashed.",
         call. = FALSE)
  }
  list(bytes = as.character(as.numeric(after$size[[1L]])), sha256 = digest)
}

report_shareable_tree_files <- function(path) {
  state <- report_shareable_state(path)
  if (identical(state$type, "missing")) return(character())
  if (!identical(state$type, "directory")) {
    stop("LISA-REPORT-SHARE-011 shareable subtree is not a directory.",
         call. = FALSE)
  }
  entries <- list.files(
    path, recursive = FALSE, full.names = TRUE, all.files = TRUE, no.. = TRUE
  )
  if (!identical(report_shareable_state(path)$type, "directory")) {
    stop("LISA-REPORT-SHARE-012 shareable subtree changed during discovery.",
         call. = FALSE)
  }
  files <- character()
  for (entry in entries) {
    entry_state <- report_shareable_state(entry)
    if (identical(entry_state$type, "directory")) {
      files <- c(files, report_shareable_tree_files(entry))
    } else if (identical(entry_state$type, "file")) {
      if (!identical(report_shareable_state(entry)$type, "file")) {
        stop("LISA-REPORT-SHARE-012 shareable file changed during discovery.",
             call. = FALSE)
      }
      files <- c(files, entry)
    } else if (identical(entry_state$type, "missing")) {
      stop("LISA-REPORT-SHARE-012 shareable entry disappeared during discovery.",
           call. = FALSE)
    } else {
      stop("LISA-REPORT-SHARE-011 invalid shareable component type.",
           call. = FALSE)
    }
  }
  files
}

report_embedded_json_payload <- function(path, data_id) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  marker <- paste0('<script type="application/json" id="', data_id, '">')
  start <- regexpr(marker, html, fixed = TRUE)[[1L]]
  if (start < 1L) return(NULL)
  from <- start + nchar(marker)
  tail <- substr(html, from, nchar(html))
  finish <- regexpr("</script>", tail, fixed = TRUE)[[1L]]
  if (finish < 1L) {
    stop("LISA-REPORT-PRODUCT-001 unterminated category-product payload: ",
         path, call. = FALSE)
  }
  tryCatch(
    jsonlite::fromJSON(
      substr(tail, 1L, finish - 1L), simplifyVector = FALSE
    ),
    error = function(error) stop(
      "LISA-REPORT-PRODUCT-001 unreadable category-product payload: ",
      path, call. = FALSE
    )
  )
}

report_category_product_reference_rows <- function(project_dir) {
  page_roots <- file.path(
    project_dir, "report_pages", c("evidence", "contrast_evidence")
  )
  pages <- unlist(lapply(page_roots[dir.exists(page_roots)], function(root) {
    list.files(
      root, pattern = "^index[.]html$", recursive = TRUE,
      full.names = TRUE
    )
  }), use.names = FALSE)
  gallery_pages <- file.path(project_dir, "report_pages", c("single_de.html", "contrasts.html"))
  pages <- c(pages, gallery_pages[file.exists(gallery_pages)])
  if (!length(pages)) return(data.frame())
  rows <- list()
  for (page in sort(pages)) {
    contrast <- grepl("/report_pages/contrast_evidence/", report_slash_path(page), fixed = TRUE)
    payload <- report_embedded_json_payload(page,
      if (page %in% gallery_pages) "report-full-gallery-data" else
        if (contrast) "contrast-evidence-data" else "evidence-data")
    if (is.null(payload) || is.null(payload$category_products)) next
    if (!is.list(payload$category_products)) {
      stop("LISA-REPORT-PRODUCT-002 category_products must be a list: ",
           page, call. = FALSE)
    }
    assets <- unlist(lapply(payload$category_products, function(product) {
      if (!is.list(product) || is.null(product$assets)) return(list())
      if (!is.list(product$assets)) {
        stop("LISA-REPORT-PRODUCT-002 product assets must be a list: ",
             page, call. = FALSE)
      }
      product$assets
    }), recursive = FALSE)
    if (!length(assets)) next
    required <- c("format", "href", "source_path", "source_name", "sha256")
    for (asset in assets) {
      if (!is.list(asset) || !all(required %in% names(asset))) {
        stop("LISA-REPORT-PRODUCT-002 incomplete category-product reference: ",
             page, call. = FALSE)
      }
      values <- vapply(required, function(field) {
        value <- asset[[field]]
        if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
          stop("LISA-REPORT-PRODUCT-002 invalid category-product field ",
               field, ": ", page, call. = FALSE)
        }
        as.character(value)
      }, character(1L))
      source_path <- report_slash_path(values[["source_path"]])
      parts <- strsplit(source_path, "/", fixed = TRUE)[[1L]]
      if (!length(parts) || !parts[[1L]] %in% c("outputs", "artifacts") ||
          any(!nzchar(parts) | parts %in% c(".", "..")) ||
          report_is_absolute_path(source_path) ||
          grepl("\\\\|[[:cntrl:]]", source_path)) {
        stop("LISA-REPORT-PRODUCT-003 category product is outside canonical roots: ",
             source_path, call. = FALSE)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        page = report_relative_path(page, project_dir),
        href = values[["href"]],
        source_path = source_path,
        source_name = values[["source_name"]],
        format = tolower(values[["format"]]),
        sha256 = values[["sha256"]],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  if (anyDuplicated(paste(out$page, out$source_path, sep = "\r"))) {
    stop("LISA-REPORT-PRODUCT-004 duplicate category-product reference.",
         call. = FALSE)
  }
  out
}

report_category_product_files <- function(project_dir) {
  references <- report_category_product_reference_rows(project_dir)
  if (!nrow(references)) return(character())
  file.path(project_dir, sort(unique(references$source_path)))
}

report_shareable_files <- function(project_dir, include_manifest = TRUE) {
  if (!identical(report_shareable_state(project_dir)$type, "directory")) {
    stop("LISA-REPORT-SHARE-011 shareable root is not a directory.",
         call. = FALSE)
  }
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  files <- character()
  root_index <- file.path(project_dir, "report_index.html")
  root_index_state <- report_shareable_state(root_index)
  if (identical(root_index_state$type, "file")) {
    files <- c(files, root_index)
  } else if (!identical(root_index_state$type, "missing")) {
    stop("LISA-REPORT-SHARE-011 report_index.html is not a regular file.",
         call. = FALSE)
  }
  for (directory in c(
    "report_pages", "report_assets", "report_media", "report_files",
    "report_figure_data"
  )) {
    path <- file.path(project_dir, directory)
    path_state <- report_shareable_state(path)
    if (identical(path_state$type, "directory")) {
      files <- c(files, report_shareable_tree_files(path))
    } else if (!identical(path_state$type, "missing")) {
      stop("LISA-REPORT-SHARE-011 shareable subtree has the wrong type: ",
           path, call. = FALSE)
    }
  }
  files <- c(files, report_category_product_files(project_dir))
  if (include_manifest) {
    manifest <- file.path(project_dir, report_shareable_manifest_name)
    manifest_state <- report_shareable_state(manifest)
    if (identical(manifest_state$type, "file")) {
      files <- c(files, manifest)
    } else if (!identical(manifest_state$type, "missing")) {
      stop("LISA-REPORT-SHARE-011 manifest is not a regular file.",
           call. = FALSE)
    }
  }
  files <- unique(files)
  if (!length(files)) return(character())
  lexical_files <- files
  invisible(lapply(
    lexical_files, report_assert_shareable_file, project_dir = project_dir
  ))
  files <- normalizePath(lexical_files, winslash = "/", mustWork = TRUE)
  invisible(lapply(
    lexical_files, report_assert_shareable_file, project_dir = project_dir
  ))
  relative <- vapply(
    files, report_relative_path, character(1), from_dir = project_dir,
    USE.NAMES = FALSE
  )
  if (any(relative == ".." | startsWith(relative, "../")) ||
      any(grepl("\\\\|(^|/)[.][.]?($|/)", relative))) {
    stop("LISA-REPORT-SHARE-002 shareable component escaped the run root.",
         call. = FALSE)
  }
  files[order(relative)]
}

report_shareable_role <- function(relative_path) {
  if (identical(relative_path, report_shareable_manifest_name)) return("manifest")
  if (identical(relative_path, "report_index.html")) return("entrypoint")
  prefix <- sub("/.*$", "", relative_path)
  switch(prefix,
    report_pages = "page",
    report_assets = "asset",
    report_media = "media",
    report_files = "portable_table",
    report_figure_data = "figure_contract",
    outputs = "canonical_product",
    artifacts = "canonical_product",
    "component"
  )
}

write_shareable_report_manifest <- function(project_dir) {
  report_assert_shareable_not_link(project_dir)
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  manifest_path <- file.path(project_dir, report_shareable_manifest_name)
  manifest_state <- report_shareable_state(manifest_path)
  if (!manifest_state$type %in% c("missing", "file")) {
    stop("LISA-REPORT-SHARE-011 manifest destination has the wrong type.",
         call. = FALSE)
  }
  files <- report_shareable_files(project_dir, include_manifest = FALSE)
  relative <- vapply(
    files, report_relative_path, character(1), from_dir = project_dir,
    USE.NAMES = FALSE
  )
  receipts <- lapply(
    files, report_shareable_file_receipt, project_dir = project_dir
  )
  rows <- data.frame(
    relative_path = relative,
    role = vapply(relative, report_shareable_role, character(1)),
    bytes = vapply(receipts, `[[`, character(1), "bytes"),
    sha256 = vapply(receipts, `[[`, character(1), "sha256"),
    integrity = "sha256",
    stringsAsFactors = FALSE
  )
  # A file cannot contain its own stable digest.  It is still part of the exact
  # copy set, but the explicit SELF marker is the sole self-reference exception.
  rows <- rbind(rows, data.frame(
    relative_path = report_shareable_manifest_name,
    role = "manifest", bytes = "SELF", sha256 = "SELF",
    integrity = "self-reference-not-hashed", stringsAsFactors = FALSE
  ))
  rows <- rows[order(rows$relative_path), , drop = FALSE]
  lisaR:::lisa_guarded_write(
    manifest_path,
    function(target) utils::write.table(
      rows, target, sep = "\t", quote = FALSE,
      row.names = FALSE, na = ""
    ),
    run_root = project_dir,
    overwrite = TRUE
  )
  report_assert_shareable_file(manifest_path, project_dir)
  invisible(manifest_path)
}

read_shareable_report_manifest <- function(project_dir) {
  path <- file.path(project_dir, report_shareable_manifest_name)
  state <- report_shareable_state(path)
  if (identical(state$type, "missing") ||
      (identical(state$type, "file") &&
       identical(as.numeric(state$info$size[[1L]]), 0))) {
    stop("LISA-REPORT-SHARE-003 shareable report manifest is missing.",
         call. = FALSE)
  }
  if (!identical(state$type, "file")) {
    stop("LISA-REPORT-SHARE-011 manifest is not a regular file.",
         call. = FALSE)
  }
  report_assert_shareable_file(path, project_dir)
  manifest <- utils::read.delim(
    path, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    colClasses = "character", na.strings = NULL, stringsAsFactors = FALSE
  )
  report_assert_shareable_file(path, project_dir)
  manifest
}

report_read_portable_text_file <- function(path) {
  connection <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  paste(readLines(connection, warn = FALSE), collapse = "\n")
}

report_assert_no_private_paths <- function(project_dir, files) {
  text_files <- files[grepl(
    "[.]html$|[.](tsv|csv)[.]gz$", files, ignore.case = TRUE
  )]
  if (!length(text_files)) return(invisible(TRUE))
  root <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  root_backslash <- gsub("/", "\\\\", root, fixed = TRUE)
  for (path in text_files) {
    report_assert_shareable_file(path, project_dir)
    text <- report_read_portable_text_file(path)
    report_assert_shareable_file(path, project_dir)
    normalized <- report_slash_path(text)
    if (grepl(root, normalized, fixed = TRUE) ||
        grepl(root_backslash, text, fixed = TRUE) ||
        grepl(report_file_uri_pattern, normalized, perl = TRUE) ||
        grepl(report_absolute_token_pattern, normalized, perl = TRUE)) {
      stop(
        "LISA-REPORT-PRIVACY-004 absolute/private path remains in portable component ",
        report_relative_path(path, project_dir), call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

report_html_links <- function(path) {
  text <- report_read_portable_text_file(path)
  matches <- gregexpr(
    "(?:href|src)=[\\\"'][^\\\"']*[\\\"']", text,
    perl = TRUE, ignore.case = TRUE
  )
  attributes <- regmatches(text, matches)[[1L]]
  if (!length(attributes)) return(character())
  links <- sub("^[^=]+=[\\\"']", "", attributes, perl = TRUE)
  sub("[\\\"']$", "", links, perl = TRUE)
}

report_validate_shareable_links <- function(project_dir, manifest) {
  declared <- as.character(manifest$relative_path)
  html <- declared[grepl("[.]html$", declared, ignore.case = TRUE)]
  for (relative_html in html) {
    html_path <- file.path(project_dir, relative_html)
    report_assert_shareable_file(html_path, project_dir)
    links <- report_html_links(html_path)
    report_assert_shareable_file(html_path, project_dir)
    links <- links[!grepl(
      "^(#|[A-Za-z][A-Za-z0-9+.-]*:|//)", links, perl = TRUE
    )]
    links <- sub("[?#].*$", "", links)
    links <- links[nzchar(links)]
    for (link in links) {
      target <- file.path(dirname(html_path), utils::URLdecode(link))
      target_state <- report_shareable_state(target)
      if (!identical(target_state$type, "file")) {
        stop(
          "LISA-REPORT-SHARE-004 portable HTML link is missing: ", link,
          " from ", relative_html, call. = FALSE
        )
      }
      report_assert_shareable_file(target, project_dir)
      target <- normalizePath(target, winslash = "/", mustWork = TRUE)
      relative_target <- report_relative_path(target, project_dir)
      if (!relative_target %in% declared) {
        stop(
          "LISA-REPORT-SHARE-005 portable HTML link is outside the manifest: ",
          relative_target, call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

report_validate_category_product_references <- function(project_dir, manifest) {
  references <- report_category_product_reference_rows(project_dir)
  if (!nrow(references)) return(invisible(TRUE))
  declared <- as.character(manifest$relative_path)
  for (index in seq_len(nrow(references))) {
    row <- references[index, , drop = FALSE]
    page <- file.path(project_dir, row$page[[1L]])
    target <- file.path(project_dir, row$source_path[[1L]])
    report_assert_shareable_file(page, project_dir)
    report_assert_shareable_file(target, project_dir)
    expected <- lisaR:::lisa_category_product_reference_asset(
      project_dir, dirname(page), target,
      allowed_roots = c("outputs", "artifacts"),
      label = "FULL category product"
    )
    fields <- c(
      "format", "href", "source_path", "source_name", "sha256"
    )
    observed <- unname(vapply(
      fields, function(field) as.character(row[[field]][[1L]]), character(1L)
    ))
    wanted <- unname(vapply(
      fields, function(field) as.character(expected[[field]]), character(1L)
    ))
    if (!identical(observed, wanted)) {
      stop("LISA-REPORT-PRODUCT-005 category-product reference failed reconciliation: ",
           row$source_path[[1L]], call. = FALSE)
    }
    if (!row$source_path[[1L]] %in% declared) {
      stop("LISA-REPORT-SHARE-005 dynamic category-product link is outside the manifest: ",
           row$source_path[[1L]], call. = FALSE)
    }
    manifest_row <- manifest[manifest$relative_path == row$source_path[[1L]],
                             , drop = FALSE]
    if (nrow(manifest_row) != 1L ||
        !identical(as.character(manifest_row$sha256[[1L]]),
                   row$sha256[[1L]])) {
      stop("LISA-REPORT-PRODUCT-006 dynamic category-product hash is not manifest-bound: ",
           row$source_path[[1L]], call. = FALSE)
    }
    resolved <- normalizePath(
      file.path(dirname(page), utils::URLdecode(row$href[[1L]])),
      winslash = "/", mustWork = TRUE
    )
    if (!identical(resolved,
                   normalizePath(target, winslash = "/", mustWork = TRUE))) {
      stop("LISA-REPORT-PRODUCT-005 category-product href does not resolve to source_path: ",
           row$source_path[[1L]], call. = FALSE)
    }
  }
  invisible(TRUE)
}

validate_shareable_report_manifest <- function(project_dir) {
  report_assert_shareable_not_link(project_dir)
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  manifest <- read_shareable_report_manifest(project_dir)
  expected_names <- c("relative_path", "role", "bytes", "sha256", "integrity")
  if (!identical(names(manifest), expected_names) || !nrow(manifest)) {
    stop("LISA-REPORT-SHARE-006 invalid shareable manifest schema.",
         call. = FALSE)
  }
  paths <- as.character(manifest$relative_path)
  if (anyNA(paths) || any(!nzchar(paths)) || anyDuplicated(paths) ||
      !identical(paths, sort(paths)) ||
      any(grepl("^/|^[A-Za-z]:/|\\\\|(^|/)[.][.]?($|/)", paths))) {
    stop("LISA-REPORT-SHARE-007 manifest paths are not unique and portable.",
         call. = FALSE)
  }
  self <- paths == report_shareable_manifest_name
  if (sum(self) != 1L || !identical(manifest$role[self], "manifest") ||
      !identical(manifest$bytes[self], "SELF") ||
      !identical(manifest$sha256[self], "SELF") ||
      !identical(manifest$integrity[self], "self-reference-not-hashed")) {
    stop("LISA-REPORT-SHARE-008 invalid manifest self-reference exception.",
         call. = FALSE)
  }
  files <- report_shareable_files(project_dir, include_manifest = TRUE)
  observed <- vapply(
    files, report_relative_path, character(1), from_dir = project_dir,
    USE.NAMES = FALSE
  )
  if (!identical(paths, sort(observed))) {
    stop("LISA-REPORT-SHARE-009 manifest and portable surface differ.",
         call. = FALSE)
  }
  component_paths <- file.path(project_dir, paths[!self])
  receipts <- lapply(
    component_paths, report_shareable_file_receipt,
    project_dir = project_dir
  )
  component_sha256 <- vapply(receipts, `[[`, character(1), "sha256")
  component_bytes <- vapply(receipts, `[[`, character(1), "bytes")
  if (any(!grepl("^[0-9a-f]{64}$", manifest$sha256[!self])) ||
      !identical(
        manifest$sha256[!self],
        component_sha256
      ) ||
      !identical(
        manifest$bytes[!self],
        component_bytes
      )) {
    stop("LISA-REPORT-SHARE-010 portable component digest/size mismatch.",
         call. = FALSE)
  }
  report_assert_no_private_paths(project_dir, files[!self])
  report_validate_shareable_links(project_dir, manifest)
  report_validate_category_product_references(project_dir, manifest)
  invisible(manifest)
}

fmt <- function(x, digits = 3) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return("")
  y <- suppressWarnings(as.numeric(x))
  if (is.na(y)) return(as.character(x))
  if (abs(y) >= 100) return(sprintf("%.0f", y))
  if (abs(y) > 0 && abs(y) < 0.001) return(sprintf("%.2e", y))
  format(round(y, digits), trim = TRUE, scientific = FALSE)
}

metric_value <- function(df, metric, default = NA_integer_) {
  if (is.null(df) || nrow(df) == 0 || !"metric" %in% names(df) || !"value" %in% names(df)) return(default)
  hit <- df$value[df$metric == metric]
  if (length(hit) == 0) return(default)
  suppressWarnings(as.integer(hit[[1]]))
}

pretty_label <- function(path) {
  label <- tools::file_path_sans_ext(basename(path))
  display_labels <- get0("report_display_labels", inherits = TRUE, ifnotfound = character())
  for (id in names(display_labels)[order(nchar(names(display_labels)), decreasing = TRUE)]) {
    if (!is.na(display_labels[[id]]) && nzchar(display_labels[[id]]) && startsWith(label, paste0(id, "_"))) {
      label <- paste0(display_labels[[id]], " · ", substring(label, nchar(id) + 2L))
      break
    }
  }
  label <- gsub("_gene_level_top_recurrent_genes$", " — recurrent genes across categories", label)
  label <- gsub("^[0-9]+_", "", label)
  label <- gsub("_[^_]+_GSEA_", " GSEA ", label)
  label <- gsub("_[^_]+_ORA_", " ORA ", label)
  label <- gsub("_category_gene_card$", "", label)
  label <- gsub("_category_volcano_overlay$", "", label)
  label <- gsub("_macrogroup_heatmap$", "", label)
  label <- gsub("_contrast_category_card$", "", label)
  label <- gsub("_contrast_macrogroup_heatmap$", "", label)
  label <- gsub("_contrast_painted$", "", label)
  label <- gsub("_paired_gene_heatmap$", " paired gene heatmap", label)
  label <- gsub("_", " ", label)
  label <- gsub("\\s+", " ", label)
  trimws(label)
}

notice_html <- function(kind, title, text) {
  sprintf(
    '<div class="notice %s"><strong>%s</strong><span>%s</span></div>',
    esc(kind), esc(title), esc(text)
  )
}

runtime_for_path <- function(runtime, path) {
  if (is.null(runtime) || nrow(runtime) == 0 || !"path" %in% names(runtime)) return(data.frame())
  normalized <- normalizePath(path, mustWork = FALSE)
  runtime[startsWith(runtime$path, normalized), , drop = FALSE]
}

runtime_notice <- function(runtime, path) {
  rows <- runtime_for_path(runtime, path)
  if (nrow(rows) == 0 || !"status" %in% names(rows)) return("")
  failed <- rows[rows$status == "failed", , drop = FALSE]
  if (nrow(failed) == 0) return('<span class="pill status-ok">complete</span>')
  present <- vapply(failed$path, function(p) {
    if (!dir.exists(p)) return(FALSE)
    length(list.files(p, pattern = "\\.(png|svg)$", recursive = TRUE, full.names = TRUE)) > 0
  }, logical(1))
  if (all(present)) return('<span class="pill status-warn">files present; post-process flagged failed</span>')
  sprintf('<span class="pill status-warn">%d failed layer(s); %d missing visual layer(s)</span>', nrow(failed), sum(!present))
}

toc_html <- function(items) {
  if (length(items) == 0) return("")
  links <- paste(vapply(items, function(item) {
    sprintf('<a href="#%s">%s</a>', esc(item[["id"]]), esc(item[["label"]]))
  }, character(1)), collapse = "")
  paste0('<details class="toc"><summary>On this page (', length(items), ' sections)</summary><div>', links, '</div></details>')
}

first_collection <- function(path) {
  dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[dir.exists(dirs)]
  if (length(dirs) == 0) return("")
  sub("^collection_", "", basename(dirs[[1]]))
}

analysis_label_lookup <- function(index, id_col, id, fallback = id) {
  if (is.null(index) || nrow(index) == 0 || !id_col %in% names(index)) return(fallback)
  rows <- index[index[[id_col]] == id, , drop = FALSE]
  if (nrow(rows) == 0) return(fallback)
  for (col in c("label", "contrast_label", "output_id", "biological_question", "question")) {
    if (col %in% names(rows) && nzchar(as.character(rows[[col]][[1]]))) return(as.character(rows[[col]][[1]]))
  }
  fallback
}

contrast_output_id <- function(index, contrast_dir_id) {
  if (!is.null(index) && nrow(index) > 0L &&
      all(c("contrast_id", "output_id") %in% names(index))) {
    composite_ids <- paste(index$contrast_id, index$output_id, sep = "_")
    matched <- which(composite_ids == contrast_dir_id)
    if (length(matched) == 1L) return(as.character(index$output_id[[matched]]))
  }
  sub("^[^_]+_", "", contrast_dir_id)
}

contrast_short_title <- function(index, contrast_dir_id) {
  if (!is.null(index) && nrow(index) > 0L &&
      all(c("contrast_id", "output_id") %in% names(index))) {
    composite_ids <- paste(index$contrast_id, index$output_id, sep = "_")
    matched <- which(composite_ids == contrast_dir_id)
    if (length(matched) == 1L && "label" %in% names(index) && nzchar(index$label[[matched]])) return(index$label[[matched]])
    if (length(matched) == 1L &&
        identical(
          as.character(index$contrast_id[[matched]]),
          as.character(index$output_id[[matched]])
        )) {
      return(short_title(as.character(index$contrast_id[[matched]])))
    }
  }
  short_title(contrast_dir_id)
}

nav_anchor_for_de <- function(project_dir, analysis_id, layer = "single_de") {
  lisa_dir <- file.path(project_dir, "outputs", "lisa", "single_de", analysis_id)
  if (!dir.exists(lisa_dir)) lisa_dir <- file.path(project_dir, "outputs", "single_de", analysis_id)
  gene_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis_id)
  collection <- if (layer == "single_de") first_collection(lisa_dir) else first_collection(gene_dir)
  if (!nzchar(collection)) collection <- "GOBP-C2"
  suffix <- switch(layer,
    single_de = "LISA categories",
    gene_cards = "GeneCards",
    volcano = "Volcano overlays",
    recurrent = "Recurrent genes",
    "LISA categories"
  )
  slug(paste(short_title(analysis_id), collection, suffix, sep = " - "))
}

nav_anchor_for_contrast <- function(project_dir, contrast_id, layer = "category",
                                    contrast_index = NULL) {
  lisa_dir <- file.path(project_dir, "outputs", "lisa", "contrast", contrast_id)
  if (!dir.exists(lisa_dir)) lisa_dir <- file.path(project_dir, "outputs", "category_contrasts", contrast_id)
  gene_dir <- file.path(project_dir, "outputs", "gene_level", "category_contrasts", contrast_id)
  collection <- if (layer == "category") first_collection(lisa_dir) else first_collection(gene_dir)
  if (!nzchar(collection)) collection <- "GOBP-C2"
  suffix <- switch(layer,
    category = "LISA category shifts",
    cards = "Contrast GeneCards",
    heatmap = "Heatmaps",
    network = "KEGG maps / painted pathways",
    "LISA category shifts"
  )
  label <- paste(contrast_short_title(contrast_index, contrast_id), collection, sep = " - ")
  if (nzchar(suffix)) label <- paste(label, suffix, sep = " - ")
  slug(label)
}

collection_order <- c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS", "HALLMARKS")
single_layer_names <- c("LISA categories", "Category evidence", "GeneCards", "Volcano overlays", "Recurrent genes", "Heatmaps", "KEGG maps / painted pathways")
contrast_layer_names <- c("LISA category shifts", "Contrast evidence", "Contrast GeneCards", "Gene-category network", "Paired gene heatmaps", "KEGG maps / painted pathways")

order_collection_dirs <- function(dirs) {
  if (length(dirs) == 0) return(dirs)
  collections <- sub("^collection_", "", basename(dirs))
  idx <- match(collections, collection_order)
  dirs[order(is.na(idx), idx, collections)]
}

single_collection_dirs <- function(project_dir, analysis_id) {
  root <- file.path(project_dir, "outputs", "lisa", "single_de", analysis_id)
  if (!dir.exists(root)) root <- file.path(project_dir, "outputs", "single_de", analysis_id)
  order_collection_dirs(list.dirs(root, recursive = FALSE, full.names = TRUE))
}

contrast_collection_dirs <- function(project_dir, contrast_id) {
  root <- file.path(project_dir, "outputs", "lisa", "contrast", contrast_id)
  if (!dir.exists(root)) root <- file.path(project_dir, "outputs", "category_contrasts", contrast_id)
  order_collection_dirs(list.dirs(root, recursive = FALSE, full.names = TRUE))
}

layer_anchor <- function(base_label, collection, layer_name) {
  slug(paste(base_label, collection, layer_name, sep = " - "))
}

nav_tree_section_html <- function(label, groups) {
  if (length(groups) == 0) return("")
  body <- paste(vapply(groups, function(group) {
    collection_html <- paste(vapply(group$collections, function(collection) {
      links <- paste(vapply(collection$layers, function(layer) {
        sprintf('<a href="%s">%s</a>', esc(layer$href), esc(layer$label))
      }, character(1)), collapse = "")
      sprintf(
        '<details class="nav-collection" open><summary>%s</summary><div class="nav-sub">%s</div></details>',
        esc(collection$label), links
      )
    }, character(1)), collapse = "")
    sprintf(
      '<details class="nav-tree" open><summary><a class="%s" href="%s">%s</a></summary><div class="nav-collections">%s</div></details>',
      if (isTRUE(group$active)) "active" else "", esc(group$href), esc(group$label), collection_html
    )
  }, character(1)), collapse = "\n")
  paste0('<div class="nav-section"><span class="nav-label">', esc(label), '</span>', body, '</div>')
}

nav_section_html <- function(label, items) {
  if (length(items) == 0) return("")
  body <- paste(vapply(items, function(item) {
    cls <- if (isTRUE(item$active)) "active" else ""
    subs <- item$subs
    sub_html <- if (length(subs) == 0) "" else paste(vapply(subs, function(sub) {
      sprintf('<a href="%s">%s</a>', esc(sub[[2]]), esc(sub[[1]]))
    }, character(1)), collapse = "")
    sprintf(
      '<div class="nav-group"><a class="%s" href="%s">%s</a>%s</div>',
      cls, esc(item$href), esc(item$label),
      if (nzchar(sub_html)) paste0('<div class="nav-sub">', sub_html, '</div>') else ""
    )
  }, character(1)), collapse = "\n")
  paste0('<div class="nav-section"><span class="nav-label">', esc(label), '</span>', body, '</div>')
}

copy_assets <- function(package_dir, project_dir) {
  asset_dir <- file.path(project_dir, "report_assets")
  dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
  logo_candidates <- c(
    file.path(package_dir, "report_assets"),
    file.path(package_dir, "inst", "report_assets"),
    file.path(package_dir, "LISA_Logo")
  )
  logo_dir <- logo_candidates[dir.exists(logo_candidates)][1]
  if (is.na(logo_dir)) stop("LISA report logo assets are not installed.", call. = FALSE)
  logos <- c(
    "LISA_logo_A1_muted_red_S_automated_annotation_final.svg",
    "LISA_logo_C_compact_icon_muted_red_S.svg",
    "LISA_logo_A1_muted_red_S_automated_annotation_final.png",
    "LISA_logo_C_compact_icon_muted_red_S.png"
  )
  for (logo in logos) {
    src <- file.path(logo_dir, logo)
    if (file.exists(src)) file.copy(src, file.path(asset_dir, logo), overwrite = TRUE)
  }
  theme <- file.path(logo_dir, "lisa_pastel.css")
  if (!file.exists(theme) || !file.copy(theme,
      file.path(asset_dir, "lisa_pastel.css"), overwrite = TRUE))
    stop("LISA report theme is not installed.", call. = FALSE)
  invisible(asset_dir)
}

# Apply only to assembled report pages; plots and their scientific colours are
# not restyled. Keeping this last in the head also covers the evidence viewers.
report_apply_theme <- function(html, page_file, project_dir) {
  relative <- report_relative_path(page_file, project_dir)
  scope <- if (identical(relative, "report_index.html")) "overview" else
    if (relative == "report_pages/contrasts.html" ||
        startsWith(relative, "report_pages/contrast_evidence/")) "contrast" else
    if (relative == "report_pages/single_de.html" ||
        startsWith(relative, "report_pages/evidence/")) "analysis" else "neutral"
  if (!grepl('data-style-scope="', html, fixed = TRUE))
    html <- sub("<body\\b", paste0('<body data-style-scope="', scope, '"'), html, perl = TRUE)
  if (!grepl('data-lisa-theme="pastel"', html, fixed = TRUE)) {
    href <- report_relative_path(file.path(project_dir, "report_assets", "lisa_pastel.css"), dirname(page_file))
    html <- sub("</head>", paste0('<link rel="stylesheet" data-lisa-theme="pastel" href="',
      esc(href), '"></head>'), html, fixed = TRUE)
  }
  html
}

report_result_card <- function(label, href, kind) {
  # A human-supplied label may have a short heading before the first colon.
  # Otherwise show it unchanged; never infer a contrast direction from IDs.
  parts <- strsplit(label, ":", fixed = TRUE)[[1L]]
  title <- if (length(parts)) trimws(parts[[1L]]) else label
  description <- if (length(parts) > 1L) trimws(paste(parts[-1L], collapse = ":")) else ""
  paste0('<a class="style-result-card style-result-card--', kind, '" href="', esc(href), '">',
    '<span class="style-card-kind">', if (kind == "analysis") "DE analysis" else "Contrast", '</span>',
    '<strong>', esc(title), '</strong>',
    if (nzchar(description)) paste0('<span>', esc(description), '</span>') else "",
    '<span class="style-card-action">Open results <span aria-hidden="true">&#8594;</span></span></a>')
}

report_overview_entry <- function(title, mode, de_index, contrast_index, project_dir,
                                  analysis_map, details) {
  cards <- character()
  if (nrow(de_index) && "analysis_id" %in% names(de_index)) {
    cards <- vapply(as.character(de_index$analysis_id), function(id) {
      report_result_card(analysis_label_lookup(de_index, "analysis_id", id, short_title(id)),
        paste0("report_pages/single_de.html#", nav_anchor_for_de(project_dir, id, "single_de")), "analysis")
    }, character(1))
  }
  if (nrow(contrast_index) && "output_id" %in% names(contrast_index)) {
    cards <- c(cards, vapply(seq_len(nrow(contrast_index)), function(i) {
      output_id <- as.character(contrast_index$output_id[[i]])
      id <- if ("contrast_id" %in% names(contrast_index))
        paste(contrast_index$contrast_id[[i]], output_id, sep = "_") else output_id
      report_result_card(analysis_label_lookup(contrast_index, "output_id", output_id, short_title(output_id)),
        paste0("report_pages/contrasts.html#", nav_anchor_for_contrast(project_dir, id, "category", contrast_index)), "contrast")
    }, character(1)))
  }
  paste0('<section class="style-overview-title"><div><span class="style-mode-label">',
    esc(toupper(mode)), '</span><h1>', esc(title), '</h1></div><p class="style-counts">',
    '<span><strong>', nrow(de_index), '</strong> ', if (nrow(de_index) == 1L) "DE analysis" else "DE analyses", '</span>',
    '<span><strong>', nrow(contrast_index), '</strong> ', if (nrow(contrast_index) == 1L) "contrast" else "contrasts", '</span></p></section>',
    if (length(cards)) paste0('<nav aria-label="Study results" class="style-summary-grid">', paste(cards, collapse = ""), '</nav>') else "",
    analysis_map,
    '<details class="style-report-details" id="style-report-details"><summary>Report details</summary>',
    '<div class="style-report-details-body">', details, '</div></details>')
}

write_assets <- function(project_dir) {
  asset_dir <- file.path(project_dir, "report_assets")
  css <- c(
    ":root{--ink:#172026;--muted:#66717b;--line:#d9dee2;--soft:#f5f7f8;--panel:#fff;--red:#9f2f2f;--red2:#c85b4c;--blue:#2c6c8f;--green:#26705b;--nav:#11181d}",
    "*{box-sizing:border-box}body{margin:0;font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:var(--ink);background:#f1f3f4;line-height:1.5}",
    "a{color:#244860;text-decoration:none}a:hover{text-decoration:underline}.content{padding:24px 28px 48px;max-width:1580px;width:100%;min-width:0;margin:0 auto}.hero{background:#fff;border:1px solid var(--line);padding:18px 20px;border-radius:8px;margin-bottom:20px}.hero h1{margin:0 0 6px;font-size:27px}.hero p{margin:0;color:var(--muted)}.report-method-links{padding:14px 0}.report-method-links .links{gap:16px;margin:10px 0}.category-navigator{overflow-x:auto;margin:14px 0;max-width:100%}.category-navigator svg{display:block;width:100%;min-width:840px;height:auto}.nes-variants{padding:0;margin:0}.nes-controls{display:flex;gap:20px;flex-wrap:wrap;align-items:end;padding:8px 0}.nes-controls label{display:grid;gap:4px;font-size:13px;font-weight:700}.nes-downloads{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:12px 0;padding:10px 12px;background:var(--soft);border-radius:6px}.nes-downloads strong{font-size:12px}.category-complement{border-top:1px solid var(--line);margin-top:18px;padding-top:12px}.category-complement>summary{font-weight:700;cursor:pointer}.nes-view-description{font-size:13px;max-width:90ch}.nes-variants select{font:inherit;padding:6px 10px;max-width:100%;border:1px solid var(--line);border-radius:5px}.nes-variant-image{display:block;max-width:100%;height:auto}.nes-variant-panel[hidden]{display:none}.report-scope-context{font-size:13px;color:var(--muted);margin:8px 0;overflow-wrap:anywhere}",
    ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}.metric{background:#fff;border:1px solid var(--line);border-radius:8px;padding:14px}.metric span{display:block;color:var(--muted);font-size:12px;text-transform:uppercase;letter-spacing:.06em}.metric strong{font-size:25px}.panel{background:#fff;border:1px solid var(--line);border-radius:8px;padding:18px;margin:16px 0}.section-head{display:flex;justify-content:space-between;gap:12px;align-items:flex-start;border-bottom:1px solid var(--line);padding-bottom:10px;margin-bottom:14px}.section-head h2,.panel h2{margin:0;font-size:20px}.eyebrow{margin:0 0 4px;color:var(--red);font-size:12px;font-weight:700;text-transform:uppercase;letter-spacing:.08em}.muted{color:var(--muted)}.pill{display:inline-flex;align-items:center;border:1px solid var(--line);border-radius:999px;padding:4px 9px;font-size:12px;color:#39444c;background:#fafafa;margin:2px}.status-ok{color:var(--green);font-weight:700}.status-warn{color:#9a6a00;font-weight:700}.links{display:flex;gap:8px;flex-wrap:wrap}.btn,.file-link{display:inline-flex;align-items:center;gap:6px;border:1px solid var(--line);background:#fff;border-radius:6px;padding:6px 9px;font-size:12px;color:#27323a;cursor:pointer}.btn:hover,.file-link:hover{border-color:#b8c1c8;text-decoration:none;background:#f9fafb}.notice{display:flex;gap:10px;align-items:flex-start;border-radius:8px;padding:12px 14px;margin:14px 0;border:1px solid #ead7a8;background:#fff8e8}.notice strong{min-width:150px}.notice.warn{border-color:#ead7a8;background:#fff8e8}.notice.error{border-color:#e8b4b4;background:#fff1f1}.notice.info{border-color:#bfd7e6;background:#eef7fb}.toc{position:relative;background:#fff;border:1px solid var(--line);border-radius:8px;padding:10px 12px;margin:0 0 14px}.toc summary{cursor:pointer;font-size:12px;text-transform:uppercase;color:var(--muted);letter-spacing:.06em;font-weight:800}.toc[open] div{display:flex;flex-wrap:wrap;gap:8px;margin-top:8px;padding-right:4px}.toc a{border:1px solid var(--line);border-radius:999px;background:#fff;padding:4px 8px;font-size:12px}.detail-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px}.detail-card{display:block;background:#fff;border:1px solid var(--line);border-radius:8px;padding:13px;color:var(--ink)}.detail-card:hover{text-decoration:none;border-color:#b8c1c8;background:#fbfcfc}.detail-card strong{display:block;color:var(--ink);font-size:15px}.detail-card span{display:block;color:var(--muted);font-size:12px;margin-top:3px}.caption-meta{color:var(--muted);font-size:11px;margin-top:2px;overflow-wrap:anywhere}",
    ".table-tools{display:flex;justify-content:space-between;gap:8px;align-items:center;margin:8px 0}.table-search{border:1px solid var(--line);border-radius:6px;padding:7px 9px;min-width:220px}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:8px}table{border-collapse:collapse;width:100%;font-size:13px;background:#fff}th,td{border-bottom:1px solid #e8ecef;padding:7px 9px;text-align:left;vertical-align:top}th{background:#f7f8f9;position:sticky;top:0;z-index:1}tr:hover td{background:#fbfcfc}",
    ".analysis-block,.collection-block{background:#fff;border:1px solid var(--line);border-radius:8px;margin:16px 0;overflow:hidden}.analysis-block>summary,.collection-block>summary{cursor:pointer;list-style:none;padding:14px 16px;border-bottom:1px solid var(--line);display:flex;align-items:center;justify-content:space-between;gap:12px}.analysis-block>summary::-webkit-details-marker,.collection-block>summary::-webkit-details-marker{display:none}.analysis-block>summary{background:#f8fafb}.analysis-block>summary strong{font-size:19px}.collection-block>summary strong{font-size:17px}.analysis-block>summary span,.collection-block>summary span{color:var(--muted);font-size:12px}.analysis-inner{padding:4px 14px 16px;background:#fbfcfc}.collection-inner{padding:4px 16px 14px}.top-tabs{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 18px}.top-tabs a{border:1px solid var(--line);background:#fff;border-radius:6px;padding:8px 10px;color:#27323a}.top-tabs a.active{background:#263139;color:#fff;border-color:#263139;text-decoration:none}",
    ".image-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:14px}.image-grid.compact{grid-template-columns:repeat(auto-fit,minmax(210px,1fr))}.fig-card{border:1px solid var(--line);border-radius:8px;background:#fff;overflow:hidden}.fig-card button{display:block;border:0;background:#fff;padding:0;width:100%;cursor:zoom-in}.fig-card img{width:100%;height:220px;object-fit:contain;background:#fafafa;border-bottom:1px solid var(--line)}.fig-caption{padding:9px 10px}.fig-caption strong{display:block;font-size:13px;white-space:normal;line-height:1.25}.fig-actions{display:flex;gap:6px;flex-wrap:wrap;margin-top:6px}.modal{position:fixed;inset:0;background:rgba(10,16,20,.82);display:none;align-items:center;justify-content:center;z-index:1000;padding:24px}.modal.open{display:flex}.modal-inner{background:#fff;border-radius:8px;max-width:96vw;max-height:94vh;padding:12px}.modal img{max-width:92vw;max-height:78vh;display:block}.modal-head{display:flex;justify-content:space-between;gap:12px;align-items:center;margin-bottom:8px}.modal-actions{display:flex;gap:8px;align-items:center}.close{border:1px solid var(--line);background:#fff;border-radius:6px;padding:5px 9px;cursor:pointer}",
    ".hommel-support-view .fig-card img{width:100%;height:auto;max-height:none}",
    ".full-products-cover>p{margin:0 0 10px;max-width:96ch}.full-product-previews{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:14px 0}.full-product-preview{margin:0;border:1px solid var(--line);border-radius:8px;overflow:hidden;background:#fff}.full-product-thumb{display:block;aspect-ratio:4/3;background:#fafafa;border-bottom:1px solid var(--line)}.full-product-thumb img{display:block;width:100%;height:100%;object-fit:contain}.full-product-preview figcaption{padding:9px 10px;display:grid;gap:4px}.full-product-preview figcaption strong{font-size:13px}.full-product-preview figcaption>span{color:var(--muted);font-size:11px;overflow-wrap:anywhere}.full-product-preview .links{gap:6px;margin-top:2px}",
    "@media(max-width:900px){.content{padding:16px}.section-head{display:block}.analysis-inner{padding:0 8px 10px}.collection-inner{padding:0 8px 10px}.panel{padding:12px}.analysis-block>summary,.collection-block>summary{flex-wrap:wrap}.table-search{min-width:0;width:100%}}@media print{.report-method-links,.toc,.table-tools,.fig-actions{display:none}.content{padding:0}.nes-variant-image{max-height:90vh;object-fit:contain}}"
  )
  js <- c(
    "function openFigure(src,title,svg){const m=document.getElementById('modal');m.querySelector('img').src=src;m.querySelector('strong').textContent=title||src;const png=m.querySelector('[data-download-png]');const svgl=m.querySelector('[data-download-svg]');png.href=src;svgl.href=svg||src;svgl.style.display=svg?'inline-flex':'none';m.classList.add('open')}",
    "function closeFigure(){document.getElementById('modal').classList.remove('open')}",
    "function copyTable(id){const t=document.getElementById(id);if(!t)return;let rows=[...t.querySelectorAll('tr')].map(r=>[...r.children].map(c=>c.innerText).join('\\t')).join('\\n');navigator.clipboard.writeText(rows)}",
    "function filterTable(input,id){const q=input.value.toLowerCase();document.querySelectorAll('#'+id+' tbody tr').forEach(r=>{r.style.display=r.innerText.toLowerCase().includes(q)?'':'none'})}",
    "document.addEventListener('keydown',e=>{if(e.key==='Escape')closeFigure()});",
    "document.querySelectorAll('[data-nes-variants]').forEach(viewer=>{const view=viewer.querySelector('[data-nes-select]'),subset=viewer.querySelector('[data-nes-subset]'),panels=[...viewer.querySelectorAll('[data-nes-panel]')];function update(save){const set=subset?subset.value:panels[0].dataset.nesPlotSet;[...view.options].forEach(o=>{o.disabled=!panels.some(p=>p.dataset.nesPanel===o.value&&p.dataset.nesPlotSet===set)});if(view.selectedOptions[0]?.disabled)view.value=[...view.options].find(o=>!o.disabled)?.value||'';panels.forEach(p=>{p.hidden=p.dataset.nesPanel!==view.value||p.dataset.nesPlotSet!==set});const description=viewer.querySelector('.nes-view-description');if(description)description.textContent=view.value==='hommel_support'?'Stars show adjusted category P from robust Hommel multiple-testing correction; minimum support d/N is in Category evidence. Mean NES is descriptive.':'Mean NES is retained in every view. Percentages use significant member sets; P25-P75 describes dispersion, not a confidence interval.';if(save){const u=new URL(location.href);u.searchParams.set('nes_view',view.value);u.searchParams.set('nes_set',set);history.replaceState(history.state,'',u.href);document.dispatchEvent(new CustomEvent('lisa:nes-view',{detail:{view:view.value,subset:set}}))}}function select(next){if(subset&&[...subset.options].some(o=>o.value===next.subset))subset.value=next.subset;if([...view.options].some(o=>o.value===next.view))view.value=next.view;update(false)}function restore(){const q=new URLSearchParams(location.search);select({view:q.get('nes_view'),subset:q.get('nes_set')})}view.addEventListener('change',()=>update(true));if(subset)subset.addEventListener('change',()=>update(true));document.addEventListener('lisa:nes-view',e=>select(e.detail));window.addEventListener('popstate',restore);restore()});"
  )
  writeLines(css, file.path(asset_dir, "lisa_report.css"), useBytes = TRUE)
  writeLines(js, file.path(asset_dir, "lisa_report.js"), useBytes = TRUE)
}

table_html <- function(df, cols = NULL, id, source = NULL, limit = 40) {
  if (is.null(df) || nrow(df) == 0) return('<p class="muted">No table rows available.</p>')
  if (!is.null(cols)) cols <- intersect(cols, names(df)) else cols <- names(df)
  df <- df[seq_len(min(nrow(df), limit)), cols, drop = FALSE]
  head <- paste(sprintf("<th>%s</th>", esc(cols)), collapse = "")
  body <- paste(apply(df, 1, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", esc(vapply(row, fmt, character(1)))), collapse = ""), "</tr>")
  }), collapse = "\n")
  paste0(
    '<div class="table-tools"><input class="table-search" placeholder="Filter table" oninput="filterTable(this,\'', id, '\')">',
    '<div><button class="btn" onclick="copyTable(\'', id, '\')">Copy table</button></div></div>',
    '<div class="table-wrap"><table id="', id, '"><thead><tr>', head, '</tr></thead><tbody>', body, '</tbody></table></div>'
  )
}

download_table_html <- function(files, page_file, id = "download_table") {
  if (length(files) == 0) return('<p class="muted">No downloadable tables found.</p>')
  rows <- paste(vapply(files, function(path) {
    sprintf(
      '<tr><td>%s</td><td>%s</td><td><a class="file-link" download href="%s">Download</a></td></tr>',
      esc(basename(path)),
      esc(dirname(path)),
      short_file_copy(path, page_file)
    )
  }, character(1)), collapse = "\n")
  paste0(
    '<div class="table-tools"><input class="table-search" placeholder="Filter downloads" oninput="filterTable(this,\'', id, '\')">',
    '<div><button class="btn" onclick="copyTable(\'', id, '\')">Copy table</button></div></div>',
    '<div class="table-wrap"><table id="', id, '"><thead><tr><th>file</th><th>folder</th><th>download</th></tr></thead><tbody>', rows, '</tbody></table></div>'
  )
}

report_image_formats <- function() {
  formats <- c("png", "svg", "pdf")
  formats[vapply(formats, function(format) report_requested(format, format == "png"), logical(1))]
}

report_asset_requested <- function(path) {
  ext <- tolower(tools::file_ext(path))
  key <- if (ext %in% c("png", "svg", "pdf")) ext else if (ext == "r") "recipes" else
    if (ext %in% c("tsv", "csv", "gz")) "source_data" else if (ext == "json") "recipes" else ""
  !nzchar(key) || report_requested(key, TRUE)
}

report_figure_preview <- function(href, title, format) {
  if (tolower(format) == "pdf") {
    sprintf('<object class="figure-preview" data="%s" type="application/pdf" aria-label="%s"><a href="%s">Open PDF</a></object>',
      esc(href), esc(title), esc(href))
  } else sprintf('<img src="%s" loading="lazy" alt="%s">', esc(href), esc(title))
}

figure_card <- function(path, page_file, title = basename(path)) {
  stem <- tools::file_path_sans_ext(path)
  requested <- report_image_formats()
  if (!length(requested)) return("")
  formats <- requested[file.exists(paste0(stem, ".", requested))]
  missing <- setdiff(requested, formats)
  if (length(missing)) stop("Requested figure format missing: ", stem, " (",
    paste(missing, collapse = ", "), ")", call. = FALSE)
  hrefs <- vapply(paste0(stem, ".", formats), short_media_copy, character(1), page_file = page_file)
  links <- paste(vapply(seq_along(formats), function(i) sprintf(
    '<a class="file-link" download href="%s">%s</a>', esc(hrefs[[i]]), toupper(formats[[i]])), character(1)), collapse = "")
  figure_contract <- figure_source_contract(path, page_file)
  source_link <- if (!report_requested("source_data", TRUE)) "" else if (nzchar(figure_contract$source_tsv)) {
    sprintf('<a class="file-link" download href="%s">Source data</a>', figure_contract$source_tsv)
  } else '<span class="pill status-warn">source data missing</span>'
  recipe_link <- if (!report_requested("recipes", TRUE)) "" else if (nzchar(figure_contract$recipe_r)) {
    sprintf('<a class="file-link" download href="%s">R script</a>', figure_contract$recipe_r)
  } else '<span class="pill status-warn">R script missing</span>'
  display_title <- pretty_label(path)
  paste0('<figure class="fig-card"><a class="figure-png" target="_blank" rel="noopener" href="', esc(hrefs[[1L]]), '">',
    report_figure_preview(hrefs[[1L]], display_title, formats[[1L]]),
    '</a><figcaption class="fig-caption"><strong>', esc(display_title),
    '</strong><div class="fig-actions">', links, source_link, recipe_link, '</div></figcaption></figure>')
}

source_data_for_figure <- function(path) {
  if (!file.exists(path)) return("")
  dir <- dirname(path)
  stem <- tools::file_path_sans_ext(basename(path))
  exact <- c(
    file.path(dir, paste0(stem, "_source.tsv")),
    file.path(dir, paste0(stem, ".tsv")),
    file.path(dir, sub("_zscore$", "_matrix.tsv", stem)),
    file.path(dir, sub("_contrast_painted$", "_nodes.tsv", stem)),
    file.path(dir, sub("_painted$", "_nodes.tsv", stem))
  )
  exact <- exact[file.exists(exact)]
  if (length(exact) > 0) {
    return(normalizePath(
      exact[[1L]], winslash = "/", mustWork = TRUE
    ))
  }
  source_candidates <- character()
  if (grepl("GSEA_.*lollipop|GSEA_lollipop", basename(path), ignore.case = TRUE)) {
    source_candidates <- c(
      source_candidates,
      list.files(file.path(dirname(dir), "lisa_tables"), pattern = "GSEA.*category_summary[.]tsv$", full.names = TRUE)
    )
  }
  if (grepl("ORA_.*barplot|ORA_barplot", basename(path), ignore.case = TRUE)) {
    source_candidates <- c(
      source_candidates,
      list.files(file.path(dirname(dir), "lisa_tables"), pattern = "ORA.*category_summary[.]tsv$", full.names = TRUE)
    )
  }
  if (grepl("GSEA_contrast_dumbbell|^cx_[a-f0-9]{12}_", basename(path), ignore.case = TRUE)) {
    source_candidates <- c(
      source_candidates,
      list.files(file.path(dirname(dir), "lisa_tables"), pattern = "GSEA.*category_contrast[.]tsv$", full.names = TRUE)
    )
  }
  if (grepl("recurrent_gene_screen[.]png$", basename(path), ignore.case = TRUE)) {
    source_candidates <- c(source_candidates,
      list.files(dir, pattern = "_top_recurrent_genes[.]tsv$", full.names = TRUE))
  }
  if (grepl("_paired_gene_heatmap[.]png$", basename(path), ignore.case = TRUE)) {
    source_candidates <- c(source_candidates,
      list.files(dir, pattern = "_paired_gene_evidence[.]tsv$", full.names = TRUE))
  }
  if (grepl("_network[.]png$", basename(path), ignore.case = TRUE)) {
    source_candidates <- c(source_candidates,
      list.files(dir, pattern = "_network_edges[.]tsv$", full.names = TRUE))
  }
  source_candidates <- source_candidates[file.exists(source_candidates)]
  if (length(source_candidates) > 0) {
    return(normalizePath(
      source_candidates[[1L]], winslash = "/", mustWork = TRUE
    ))
  }
  ""
}

# Build one immutable, figure-specific source-data contract.  The report never
# points a card at a shared table, manifest, or index: every target carries a
# unique figure_id, source row order and explicit defaults for visual flags.
stable_text_hash <- function(text) {
  ints <- utf8ToInt(enc2utf8(text))
  if (!length(ints)) return("0000000000000000")
  h1 <- sum((seq_along(ints) %% 65521) * ints) %% 2147483647
  h2 <- sum(((seq_along(ints) * 131) %% 65521) * rev(ints)) %% 2147483647
  paste0(sprintf("%08x", as.integer(h1)), sprintf("%08x", as.integer(h2)))
}

figure_id_for_path <- function(path, project_dir) {
  rel <- report_relative_path(path, project_dir)
  paste0("F_", stable_text_hash(rel))
}

figure_type_for_path <- function(path) {
  b <- basename(path)
  if (grepl("GSEA_category_pathways", path)) return("lisa_category_gene_sets")
  if (grepl("category_volcano_overlay", b)) return("volcano_overlay")
  if (grepl("category_gene_card", b) && !grepl("contrast", b)) return("gene_card")
  if (grepl("contrast_category_card", b)) return("contrast_gene_card")
  if (grepl("contrast_painted", b)) return("kegg_contrast")
  if (grepl("_painted[.](png|svg|pdf)$", b)) return("kegg_single")
  if (grepl("heatmap", b)) return("heatmap")
  if (grepl("volcano", b)) return("volcano")
  if (grepl("network", b)) return("network")
  if (grepl("lollipop", b)) return("lisa_lollipop")
  if (grepl("barplot", b)) return("lisa_barplot")
  if (grepl("dumbbell|^cx_[a-f0-9]{12}_", b)) return("lisa_dumbbell")
  "generic"
}

figure_source_contract <- function(path, page_file) {
  source <- source_data_for_figure(path)
  source_requested <- report_requested("source_data", TRUE)
  recipe_requested <- report_requested("recipes", TRUE)
  if (!source_requested && !recipe_requested) return(list(source_tsv = "", recipe_r = "", figure_id = ""))
  if (!nzchar(source) || !file.exists(source)) return(list(source_tsv = "", recipe_r = "", figure_id = ""))
  root <- report_root_for_page(page_file)
  data_dir <- file.path(root, "report_figure_data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  figure_id <- figure_id_for_path(path, root)
  # Figure-specific volcano and GeneCard tables can contain tens of thousands
  # of rows. Keep the one-file-per-figure contract, but store it as a
  # transparently readable gzip stream so a complete auditable report does not
  # multiply into tens of gigabytes of repeated plain text.
  source_dest <- file.path(data_dir, paste0(figure_id, ".tsv.gz"))
  recipe_dest <- file.path(data_dir, paste0(figure_id, "_recipe.R"))
  if (source_requested && !file.exists(source_dest)) {
    x <- lisaR:::lisa_read_figure_source_tsv(source)
    if (!nrow(x) && file.info(source)$size > 0) stop("LISA-REPORT-SOURCE-001 unable to read figure TSV: ", source, call. = FALSE)
    x$figure_id <- figure_id
    # A source emitted by the original renderer is authoritative (including
    # the single-DE 'dumbbell' alias, which is a lollipop, not an A/B plot).
    if (!"figure_type" %in% names(x)) x$figure_type <- figure_type_for_path(path)
    x$source_row_order <- seq_len(nrow(x))
    if (!"selected_for_plot" %in% names(x)) x$selected_for_plot <- TRUE
    if (!"highlighted" %in% names(x)) x$highlighted <- x$selected_for_plot
    if (!"labelled" %in% names(x)) x$labelled <- FALSE
    x$source_file <- basename(source)
    if (grepl("^kegg_", unique(x$figure_type)[[1]])) {
      base_name <- if ("base_image_file" %in% names(x)) unique(as.character(x$base_image_file))[[1]] else ""
      base_source <- if (nzchar(base_name)) file.path(dirname(source), base_name) else ""
      if (!nzchar(base_source) || !file.exists(base_source)) stop("LISA-REPORT-SOURCE-005 KEGG base image missing: ", path, call. = FALSE)
      base_dest <- file.path(data_dir, paste0(figure_id, "_kegg_base.png"))
      if (!file.copy(base_source, base_dest, overwrite = TRUE)) stop("LISA-REPORT-SOURCE-006 cannot copy KEGG base image.", call. = FALSE)
      x$base_image_file <- basename(base_dest)
    }
    # Figure contracts belong to the portable report surface. Scientific
    # values are retained; only fields containing filesystem paths are rebased
    # or reduced to their external basename.
    x <- report_portable_data_frame(x, root)
    if (identical(unique(x$figure_type), "lisa_category_inference") || "support_schema_version" %in% names(x))
      x <- lisaR:::lisa_inference_exact_columns(x)
    con <- gzfile(source_dest, open = "wt", compression = 9)
    tryCatch(lisaR:::lisa_write_figure_source_tsv(x, con), finally = close(con))
  }
  if (recipe_requested && !file.exists(recipe_dest)) {
    # Executable recipes always come from the verified package closure.  Do
    # not prefer a sibling recipe in the scientific output tree: that would
    # allow mutable run data to select code copied into the final report.
    renderer <- file.path(report_package_dir, "scripts", "reproduce_lisa_figure.R")
    if (!file.exists(renderer)) {
      renderer <- file.path(report_package_dir, "inst", "scripts", "reproduce_lisa_figure.R")
    }
    if (!file.exists(renderer)) {
      stop("LISA-REPORT-SOURCE-007 executable reproduction recipe missing for: ", path, call. = FALSE)
    }
    lisaR:::lisa_copy_verified_figure_recipe(
      renderer, recipe_dest, report_renderer_sha256, run_root = root
    )
  }
  list(
    source_tsv = if (source_requested) rel_path(source_dest, page_file) else "",
    recipe_r = if (recipe_requested) rel_path(recipe_dest, page_file) else "",
    figure_id = figure_id
  )
}

validate_figure_source_contracts <- function(project_dir) {
  data_dir <- file.path(project_dir, "report_figure_data")
  tsvs <- list.files(data_dir, pattern = "^F_[0-9a-f]{16}[.]tsv[.]gz$", full.names = TRUE)
  recipes <- list.files(data_dir, pattern = "^F_[0-9a-f]{16}_recipe[.]R$", full.names = TRUE)
  if (!report_requested("source_data", TRUE) || !length(tsvs)) return(invisible(data.frame()))
  rows <- lapply(tsvs, function(path) {
    x <- read_tsv(path)
    required <- c("figure_id", "source_row_order", "selected_for_plot", "highlighted", "labelled")
    if (!all(required %in% names(x)) || length(unique(x$figure_id)) != 1L) {
      stop("LISA-REPORT-SOURCE-002 invalid figure-specific source schema: ", path, call. = FALSE)
    }
    id <- unique(x$figure_id)[[1]]
    recipe <- file.path(data_dir, paste0(id, "_recipe.R"))
    # The external filename is keyed by the deterministic figure ID, never by
    # a generic table role.  It must resolve to a companion recipe.
    if (report_requested("recipes", TRUE)) {
      if (!file.exists(recipe)) stop("LISA-REPORT-SOURCE-003 figure recipe missing: ", id, call. = FALSE)
      tryCatch(parse(recipe), error = function(e) stop("LISA-REPORT-SOURCE-008 invalid figure recipe ", id, ": ", conditionMessage(e), call. = FALSE))
    }
    figure_type <- unique(as.character(x$figure_type))[[1]]
    required_by_type <- list(
      lisa_category_inference = c("category_id", "category_order", "category_p_adjusted", "significant", "minimum_enriched_sets", "minimum_enriched_pct", "NES"),
      lisa_category_gene_sets = c("plot_x_nes", "plot_point_size_neg_log10_fdr", "plotted_order", "pathway_label"),
      volcano_overlay = c("symbol", "plot_x_log2FC", "plot_y_neg_log10_fdr", "threshold_abs_log2FC", "threshold_de_fdr"),
      gene_card = c("symbol", "log2FC", "padj", "gene_contribution_score", "plotted_order"),
      contrast_gene_card = c("symbol", "log2FC_A", "log2FC_B", "plotted_order"),
      contrast_heatmap = c("symbol", "category_display_name", "log2FC_A", "log2FC_B"),
      kegg_single = c("base_image_file", "x", "y", "width", "height", "log2FC"),
      kegg_contrast = c("base_image_file", "x", "y", "width", "height", "log2FC_A", "log2FC_B")
    )
    type_required <- required_by_type[[figure_type]]
    if (!is.null(type_required) && !all(type_required %in% names(x))) {
      stop("LISA-REPORT-SOURCE-009 incomplete figure data for ", id, " (", figure_type, ").", call. = FALSE)
    }
    if (grepl("^kegg_", figure_type)) {
      base <- file.path(data_dir, unique(as.character(x$base_image_file))[[1]])
      if (!file.exists(base)) stop("LISA-REPORT-SOURCE-010 packaged KEGG base image missing for ", id, call. = FALSE)
    }
    data.frame(figure_id = id, figure_type = figure_type, source_tsv = basename(path),
      recipe_r = if (report_requested("recipes", TRUE)) basename(recipe) else "not_requested",
      rows = nrow(x), recipe_parse = if (report_requested("recipes", TRUE)) "PASS" else "not_requested",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (anyDuplicated(out$figure_id) || anyDuplicated(out$source_tsv) ||
      (report_requested("recipes", TRUE) && length(recipes) != nrow(out))) {
    stop("LISA-REPORT-SOURCE-004 figure source targets must be unique one-to-one contracts.", call. = FALSE)
  }
  out
}

validate_category_member_plot_coverage <- function(project_dir) {
  manifests <- list.files(project_dir, pattern = "plot_manifest[.]tsv$", recursive = TRUE, full.names = TRUE)
  manifests <- manifests[grepl("GSEA_category_pathways", manifests, fixed = TRUE)]
  if (!length(manifests)) return(data.frame())
  rows <- lapply(manifests, function(manifest_path) {
    d <- dirname(manifest_path)
    manifest <- read_tsv(manifest_path)
    plots <- lisaR:::lisa_plot_files(d, formats = report_image_formats())
    plots <- plots[!grepl("_regenerated[.]png$|_kegg_base[.]png$", plots)]
    source_ok <- vapply(plots, function(p) file.exists(paste0(tools::file_path_sans_ext(p), "_source.tsv")), logical(1))
    recipe_ok <- vapply(plots, function(p) file.exists(paste0(tools::file_path_sans_ext(p), "_recipe.R")), logical(1))
    source_requested <- report_requested("source_data", TRUE)
    recipe_requested <- report_requested("recipes", TRUE)
    expected <- if (length(report_image_formats())) nrow(manifest) else 0L
    observed <- length(plots)
    status <- if (expected == observed &&
        (!source_requested || all(source_ok)) &&
        (!recipe_requested || all(recipe_ok))) "PASS" else "FAIL"
    data.frame(directory = d, expected_category_plots = expected, observed_category_plots = observed,
      source_contracts = if (source_requested) sum(source_ok) else "not_requested",
      recipes = if (recipe_requested) sum(recipe_ok) else "not_requested",
      status = status, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (any(out$status != "PASS")) {
    stop("LISA-REPORT-CATEGORY-001 category member-plot coverage or reproduction contracts are incomplete.", call. = FALSE)
  }
  out
}

# Inventory comes from the sections actually assembled below, including full
# products. It is data, never executable code or a hand-maintained menu.
report_navigation_json <- function(page, contexts = list()) {
  json <- jsonlite::toJSON(list(version = 1L, page = page, contexts = contexts),
    auto_unbox = TRUE, null = "null", digits = NA)
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  paste0('<script type="application/json" id="lisa-report-navigation">', json, '</script>')
}

report_navigation_section <- function(title, label) {
  list(id = slug(title), label = label, kind = switch(label,
    "LISA category shifts" = "lisa-categories", "LISA categories" = "lisa-categories",
    "Category evidence" = "category-evidence", "Contrast evidence" = "category-evidence",
    # The FULL extras section needs its own filterable kind: without it the
    # label falls through to slug(label) and the section dropdown cannot be
    # told apart from an ordinary layer, either by a reader or by a test.
    "FULL figures" = "full-products",
    slug(label)), href = paste0("#", slug(title)))
}

report_write_navigation_inventory <- function(project_dir, analyses, contrasts) {
  rebase <- function(contexts, page) lapply(contexts, function(context) {
    context$collections <- lapply(context$collections, function(collection) {
      collection$sections <- lapply(collection$sections, function(section) {
        section$href <- paste0("report_pages/", page, "#", section$id)
        section
      })
      collection
    })
    context
  })
  inventory <- list(version = 1L, contexts = c(rebase(analyses, "single_de.html"),
    rebase(contrasts, "contrasts.html")))
  target <- file.path(project_dir, "report_pages", "navigation_inventory.json")
  # A raw in-place write here previously corrupted any tree sharing this path
  # via a hard link (e.g. a reused/reassembled project directory), because it
  # mutated the shared inode instead of replacing the directory entry. Use
  # the same atomic guarded writer already used for every other generated
  # JSON/HTML artifact in this script.
  lisaR:::lisa_guarded_write(target,
    function(staged) jsonlite::write_json(inventory, staged, auto_unbox = TRUE, pretty = TRUE, null = "null"))
  invisible(target)
}

# Evidence pages are written once by the scientific evidence-generation
# pipeline and are treated as pre-existing, read-only inputs by this
# generator -- but their static executable/style assets (not scientific
# data) must always match the currently installed package, including when
# the evidence directory was inherited via a hard link from another
# project tree. A stale copy is invisible to any content check (the
# embedded JSON payload is untouched) but can silently break client-side
# rendering. Refresh unconditionally via the atomic guarded writer so a
# shared inode is replaced, never mutated in place.
report_refresh_evidence_static_assets <- function(evidence_dir, kind) {
  asset_root <- system.file(kind, package = "lisaR")
  if (!nzchar(asset_root)) stop("Installed ", kind, " assets are missing.", call. = FALSE)
  asset_dir <- file.path(evidence_dir, "assets")
  if (!dir.exists(asset_dir)) return(invisible(NULL))
  for (extension in c("css", "js")) {
    from <- file.path(asset_root, paste0("viewer.", extension))
    if (!file.exists(from)) next
    destination <- file.path(asset_dir, paste0(kind, ".", extension))
    if (!file.exists(destination)) next
    lisaR:::lisa_guarded_write(destination, function(target) {
      if (!file.copy(from, target, overwrite = TRUE)) stop("Cannot refresh ", kind, " static asset.", call. = FALSE)
    })
  }
  invisible(asset_dir)
}

# Main report pages and standalone evidence share one fixed navigation shell.
page_shell <- function(title, active, body, project_dir, page_file, root = FALSE,
                       navigation = list()) {
  gallery_assets <- report_gallery_assets[[report_relative_path(page_file, project_dir)]]
  if (length(gallery_assets)) {
    gallery_json <- jsonlite::toJSON(list(category_products = list(list(assets = gallery_assets))),
      auto_unbox = TRUE, null = "null", digits = 17)
    gallery_json <- gsub("<", "\\u003c", gallery_json, fixed = TRUE)
    body <- paste0(body, '<script type="application/json" id="report-full-gallery-data">',
      gallery_json, '</script>')
  }
  ns <- asNamespace("lisaR")
  routes <- get("lisa_presentation_routes", ns)(project_dir, dirname(page_file), report_build = TRUE)
  route <- switch(active, single_de = "analyses", downloads = "methods", qc = "methods", active)
  if (!route %in% names(routes)) route <- NULL
  assets <- get(".lisa_copy_report_shell_assets", ns)(project_dir,
    asset_subdir = "report_assets/lisa-shell")
  for (key in intersect(c("css_href", "js_href", "logo_href", "logo_compact_href"), names(assets))) {
    assets[[key]] <- report_relative_path(file.path(project_dir, assets[[key]]), dirname(page_file))
  }
  shell <- get(".lisa_report_shell", ns)(routes = routes, active = route,
    context = list(study = title), assets = assets, main_id = "lisa-main")
  asset_prefix <- report_relative_path(file.path(project_dir, "report_assets"), dirname(page_file))
  qc <- rel_path(file.path(project_dir, "report_pages", "qc.html"), page_file)
  methods <- rel_path(file.path(project_dir, "report_pages", "downloads.html"), page_file)
  method_links <- paste0('<details class="report-method-links"><summary>Methods and checks</summary><div class="links">',
    '<a href="', methods, '">Data and methods</a>',
    '<a href="', qc, '#report-validation">Report validation</a>',
    if (file.exists(file.path(project_dir, "post_lisa_status.tsv"))) paste0(
      '<a href="', qc, '#qc-post-lisa-status-tsv">Additional outputs</a>') else "", '</div></details>')
  paste0(
    '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>', esc(title), '</title><link rel="stylesheet" href="', esc(asset_prefix), '/lisa_report.css">',
    shell$head, '</head><body class="', shell$body_class, '">', shell$header,
    report_navigation_json(active, navigation),
    '<main id="lisa-main" class="content" tabindex="-1">', body, method_links, '</main>',
    lisaR:::lisa_hommel_help_html(),
    '<div id="modal" class="modal" onclick="if(event.target.id===\'modal\')closeFigure()"><div class="modal-inner"><div class="modal-head"><strong></strong><div class="modal-actions"><a class="file-link" data-download-png download href="#">PNG</a><a class="file-link" data-download-svg download href="#">SVG</a><button class="close" onclick="closeFigure()">Close</button></div></div><img alt=""></div></div>',
    '<script src="', esc(asset_prefix), '/lisa_report.js"></script></body></html>'
  )
}

report_scope_context <- function(id, label, collection, metadata_path = NULL, contrast = FALSE) {
  metadata <- list()
  if (!is.null(metadata_path) && file.exists(metadata_path)) {
    tab <- utils::read.delim(metadata_path, sep = "\t", quote = "", comment.char = "",
      colClasses = "character", na.strings = NULL, check.names = FALSE)
    if (!all(c("key", "value") %in% names(tab)) || anyDuplicated(tab$key))
      stop("Invalid evidence metadata table.", call. = FALSE)
    metadata <- as.list(stats::setNames(tab$value, tab$key))
  }
  exact_id <- if (contrast && !is.null(metadata$contrast_id)) metadata$contrast_id else id
  ids <- stats::setNames(list(exact_id), if (contrast) "contrast_id" else "analysis_id")
  if (contrast) for (key in c("analysis_a", "analysis_b")) if (!is.null(metadata[[key]])) ids[[key]] <- metadata[[key]]
  ids$collection <- collection
  if (!is.null(metadata$tier)) ids$tier <- metadata$tier
  context <- list(analysis = label, collection = collection, selection = "",
    direction = if (contrast) "Profiles: A and B; difference products: A minus B" else if (!is.null(metadata$positive_contrast)) metadata$positive_contrast else "",
    cutoff = if (!is.null(metadata$gsea_padj_cutoff)) paste("GSEA FDR <=", metadata$gsea_padj_cutoff) else "",
    exact_ids = ids)
  esc(jsonlite::toJSON(context, auto_unbox = TRUE, digits = NA, null = "null"))
}

# Inline the generated scientific SVG, not its standalone HTML document. Rebase
# links and namespace accessibility IDs so all scopes remain distinct offline.
report_inline_category_navigation <- function(navigation_dir, page_file) {
  path <- file.path(navigation_dir, "category_navigation.svg")
  if (!file.exists(path)) return("")
  svg <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!grepl("^\\s*<svg\\b", svg, perl = TRUE))
    stop("Category navigator is not an SVG fragment.", call. = FALSE)
  prefix <- paste0("nav-", stable_text_hash(report_relative_path(path, dirname(page_file))), "-")
  for (id in c("navigation-title", "navigation-description")) svg <- gsub(id, paste0(prefix, id), svg, fixed = TRUE)
  matches <- gregexpr('(?:xlink:)?href="[^"]+"', svg, perl = TRUE)
  links <- regmatches(svg, matches)[[1L]]
  if (length(links)) {
    replacements <- vapply(links, function(attribute) {
      href <- sub('^[^=]+="', "", sub('"$', "", attribute))
      if (startsWith(href, "#")) return(attribute)
      if (grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^[/\\\\]", href))
        stop("Navigator links must remain report-relative.", call. = FALSE)
      target <- sub("[?#].*$", "", href)
      suffix <- substring(href, nchar(target) + 1L)
      target_path <- file.path(navigation_dir, utils::URLdecode(target))
      if (!file.exists(target_path)) stop("Navigator target is missing: ", target, call. = FALSE)
      prefix_attr <- sub('=".*$', "", attribute)
      paste0(prefix_attr, '="', rel_path(target_path, page_file), suffix, '"')
    }, character(1L))
    regmatches(svg, matches) <- list(replacements)
  }
  paste0('<div class="category-navigator" aria-label="Category navigator map">', svg, '</div>')
}

report_category_evidence_section <- function(title, navigation_dir, evidence_index, page_file,
                                             unclassified = character(), collection_dir = NULL) {
  if (!file.exists(evidence_index)) return("")
  nav_svg <- file.path(navigation_dir, "category_navigation.svg")
  source <- file.path(navigation_dir, "tables", "category_navigation_source.tsv")
  paste0('<section class="panel category-evidence-cover" id="', esc(slug(title)), '">',
    '<h2>Category evidence</h2>',
    '<p>Select a category to inspect its gene sets and genes.</p>',
    if (file.exists(source) && "support_schema_version" %in% names(read_tsv(source))) lisaR:::lisa_hommel_help_button() else "",
    report_inline_category_navigation(navigation_dir, page_file),
    '<div class="links"><a class="file-link" href="', rel_path(evidence_index, page_file), '">Explore a category</a>',
    if (report_requested("svg", TRUE) && file.exists(nav_svg)) paste0('<a class="file-link" download href="', rel_path(nav_svg, page_file), '">Download navigator SVG</a>') else "",
    if (report_requested("source_data", TRUE) && file.exists(source)) paste0('<a class="file-link" download href="', rel_path(source, page_file), '">Download category data</a>') else "",
    '</div>',
    if (!is.null(collection_dir)) report_category_inference_downloads(collection_dir, page_file) else "",
    if (length(unclassified)) paste0('<details class="report-unclassified-products"><summary>Unclassified saved product files (',
      length(unclassified), ')</summary><p class="muted">These files could not be attached to one selected category from explicit source metadata. They remain available here and are not assigned by filename.</p><div class="links">',
      paste(vapply(sort(unique(unclassified)), function(path) report_unclassified_product_link(path, page_file), character(1L)), collapse = ""), '</div></details>') else "",
    '</section>')
}

# ---------------------------------------------------------------------------
# FULL extras: make the attached saved products reachable from the report's
# own navigation.
#
# The saved FULL products are one-category products: they belong inside each
# category's evidence sheet and must never be duplicated as a second top-level
# gallery, so the generator removes the corresponding standalone layers for a
# selected owner. Until now it inserted nothing in their place, and the result
# was that a selected owner's page, its table of contents and the section
# dropdown were label-identical to an unselected owner's: the extras existed,
# were correctly attached and correctly referenced, but were unreachable -- and
# the only FULL-specific content left in the main navigation was the honest
# "not generated by this engine" notices, which communicate the opposite of the
# truth. These helpers add the missing positive signal: per-family counts using
# each product's own recorded label verbatim, canonical download links and
# per-category deep links, built entirely from the already-computed attachment
# payload. No figure is rendered, no directory is scanned a second time, no
# byte is copied and nothing is written under `artifacts/` or `outputs/`.
# ---------------------------------------------------------------------------

report_evidence_scope_values <- function(evidence_index) {
  metadata_path <- file.path(dirname(evidence_index), "tables", "metadata.tsv")
  tab <- tryCatch(read_tsv(metadata_path), error = function(error) data.frame())
  if (!all(c("key", "value") %in% names(tab)) || anyDuplicated(tab$key)) return(character())
  stats::setNames(as.character(tab$value), as.character(tab$key))
}

# The two evidence viewers deliberately parse deep links differently: the
# single-analysis viewer reads a query parameter and its category key is
# `category` (category-evidence/viewer.js), while the contrast viewer reads a
# URL fragment and its category key is `category_id`
# (contrast-evidence/viewer.js). One naive builder would emit contrast links
# that load the page and silently select the FIRST category instead of the
# requested one, with no error anywhere. Emit exactly the form each viewer
# parses, and take the scope keys from the page's own recorded metadata so a
# link can never disagree with the page it points at.
report_evidence_deeplink <- function(evidence_index, page_file, kind = c("category", "contrast"),
                                     owner = NULL, collection = NULL, tier = NULL,
                                     category_id = NULL) {
  kind <- match.arg(kind)
  base <- report_relative_path(evidence_index, dirname(page_file))
  scope <- report_evidence_scope_values(evidence_index)
  pick <- function(key, supplied = NULL) {
    value <- if (length(supplied) == 1L && !is.na(supplied) && nzchar(as.character(supplied))) {
      as.character(supplied)
    } else if (key %in% names(scope)) scope[[key]] else ""
    if (length(value) != 1L || is.na(value)) "" else value
  }
  pairs <- if (identical(kind, "category")) {
    c(category = pick("category_id", category_id), analysis_id = pick("analysis_id", owner),
      collection = pick("collection", collection), tier = pick("tier", tier))
  } else {
    c(contrast_id = pick("contrast_id", owner), analysis_a = pick("analysis_a"),
      analysis_b = pick("analysis_b"), collection = pick("collection", collection),
      tier = pick("tier", tier), category_id = pick("category_id", category_id))
  }
  pairs <- pairs[nzchar(pairs)]
  if (!length(pairs)) return(base)
  query <- paste(paste0(names(pairs), "=",
    vapply(pairs, function(value) utils::URLencode(value, reserved = TRUE), character(1L))),
    collapse = "&")
  paste0(base, if (identical(kind, "category")) "?" else "#", query)
}

report_full_gallery_section <- function(title, layers, page_file, project_dir) {
  blocks <- character()
  assets <- list()
  for (name in names(layers)) {
    spec <- layers[[name]]
    if (!dir.exists(spec$path)) next
    paths <- list.files(spec$path, recursive = TRUE, full.names = TRUE)
    paths <- paths[tolower(tools::file_ext(paths)) %in% c("png", "svg", "pdf", "tsv", "csv", "r", "json")]
    paths <- Filter(report_asset_requested, paths)
    if (!length(paths)) next
    links <- vapply(paths, function(path) {
      asset <- lisaR:::lisa_category_product_reference_asset(project_dir, dirname(page_file),
        path, allowed_roots = c("artifacts"), label = "FULL gallery product")
      assets[[length(assets) + 1L]] <<- asset
      paste0('<a class="file-link" download href="', esc(asset$href), '">',
        esc(basename(path)), '</a>')
    }, character(1L))
    blocks <- c(blocks, paste0('<details><summary>', esc(name), '</summary><div class="links">',
      paste(links, collapse = ""), '</div></details>'))
  }
  if (!length(assets)) return("")
  key <- report_relative_path(page_file, project_dir)
  combined <- c(report_gallery_assets[[key]], assets)
  ids <- vapply(combined, function(x) x$source_path, character(1L))
  report_gallery_assets[[key]] <<- combined[!duplicated(ids)]
  paste0('<section class="panel full-products-gallery" id="', esc(slug(title)),
    '"><h2>FULL figures</h2><p>Additional figures for this collection are shown in the galleries below. ',
    'Download their files by figure family.</p>', paste(blocks, collapse = ""), '</section>')
}

report_full_products_section <- function(title, attachment, evidence_index, page_file,
                                         kind = c("category", "contrast"), owner = NULL,
                                         collection = NULL, tier = NULL, project_dir = NULL,
                                         preview_limit = Inf, deeplink_limit = 24L, related_sections = character(), native_layers = list()) {
  kind <- match.arg(kind)
  products <- if (is.list(attachment)) attachment$category_products else NULL
  if (!is.list(products) || !length(products)) return("")
  products <- Filter(function(p) !identical(p$product, "kegg"), products)
  if (!length(products)) return("")
  if (!file.exists(evidence_index)) return("")
  if (is.null(project_dir)) {
    project_dir <- get0("report_project_dir", inherits = TRUE, ifnotfound = NULL)
  }
  if (is.null(project_dir) || length(project_dir) != 1L || is.na(project_dir) ||
      !nzchar(project_dir)) {
    stop("FULL product navigation requires the canonical report root.", call. = FALSE)
  }
  page_dir <- dirname(page_file)
  # A displayed link points at the one canonical artifact the payload already
  # recorded, resolved from this page. Never a copy, never a second inventory
  # scan and never a path outside the allowed canonical roots.
  canonical <- function(asset) {
    if (!is.list(asset)) return(NULL)
    source_path <- as.character(asset$source_path)
    if (length(source_path) != 1L || is.na(source_path) || !nzchar(source_path)) return(NULL)
    parts <- strsplit(source_path, "/", fixed = TRUE)[[1L]]
    if (!length(parts) || !parts[[1L]] %in% c("artifacts", "outputs")) return(NULL)
    target <- file.path(project_dir, source_path)
    if (!file.exists(target)) return(NULL)
    format <- as.character(asset$format)
    list(format = if (length(format) == 1L && !is.na(format) && nzchar(format)) toupper(format) else "FILE",
      href = report_relative_path(target, page_dir))
  }
  # A family name is never synthesised: it is the label the inventory already
  # recorded for that product, verbatim. "KEGG gene sets" (an artifact gene-set
  # plot) and "KEGG maps" (an in-run painted pathway diagram) are different
  # products and must never be merged or relabelled into each other.
  family_of <- function(product) {
    label <- as.character(product$label)
    if (length(label) == 1L && !is.na(label) && nzchar(label)) return(label)
    as.character(product$product)[[1L]]
  }
  labels <- vapply(products, family_of, character(1L))
  categories <- vapply(products, function(product)
    as.character(product$category_id)[[1L]], character(1L))
  unique_categories <- sort(unique(categories))
  family_pills <- vapply(sort(unique(labels)), function(label) {
    n <- length(unique(categories[labels == label]))
    sprintf('<span class="pill">%s: %d categor%s</span>', esc(label), n,
      if (n == 1L) "y" else "ies")
  }, character(1L))
  deeplink <- function(category_id) report_evidence_deeplink(evidence_index, page_file, kind,
    owner, collection, tier, category_id)
  # A bounded preview, in a stable order. It exists so the section is visibly a
  # figure section rather than a list of words; the complete set stays in the
  # category sheets, which is the one place it is not duplicated. Families are
  # interleaved rather than taken in block order, so a reader sees every family
  # that the counts above claim, not `preview_limit` copies of the first one.
  by_family <- split(order(labels, categories, method = "radix"),
    labels[order(labels, categories, method = "radix")])
  by_family <- by_family[sort(names(by_family))]
  # One representative per family; full category inventories remain available.
  interleaved <- vapply(by_family, function(group) group[[1L]], integer(1L))
  preview <- character()
  for (index in interleaved) {
    if (length(preview) >= preview_limit) break
    assets <- lapply(products[[index]]$assets, canonical)
    assets <- assets[!vapply(assets, is.null, logical(1L))]
    image <- Filter(function(asset) asset$format %in% c("PNG", "SVG", "PDF"), assets)
    if (!length(image)) next
    downloads <- paste(vapply(assets, function(asset) sprintf(
      '<a class="file-link" download href="%s">%s</a>', esc(asset$href), esc(if (asset$format == "R") "R script" else asset$format)),
      character(1L)), collapse = "")
    asset_formats <- vapply(assets, function(asset) asset$format, character(1L))
    if (report_requested("svg", FALSE) && "PNG" %in% asset_formats && !"SVG" %in% asset_formats) {
      downloads <- paste0(downloads, if (grepl("kegg_pathway_map|contrast_kegg_map|painted",
        paste(products[[index]]$product, paste(vapply(assets, `[[`, "", "href"), collapse = " ")))) {
        '<span class="pill status-warn">SVG not available for this product type</span>'
      } else '<span class="pill status-warn">SVG missing</span>')
    }
    preview <- c(preview, paste0(
      '<figure class="full-product-preview"><a class="full-product-thumb" target="_blank" rel="noopener" href="', esc(image[[1L]]$href), '">',
      report_figure_preview(image[[1L]]$href, paste(labels[[index]], "for category", categories[[index]]), image[[1L]]$format),
      '</a><figcaption><strong>', esc(labels[[index]]), '</strong><span>', esc(categories[[index]]),
      '</span><span class="links">', downloads, '</span></figcaption></figure>'))
  }
  shown <- utils::head(unique_categories, deeplink_limit)
  category_links <- paste(vapply(shown, function(category_id) sprintf(
    '<a class="file-link" href="%s">%s</a>', esc(deeplink(category_id)), esc(category_id)),
    character(1L)), collapse = "")
  paste0('<section class="panel full-products-cover" id="', esc(slug(title)), '">',
    '<h2>FULL figures</h2>',
    '<p>', length(products), ' saved figure product',
    if (length(products) == 1L) " is" else "s are", ' attached to ',
    length(unique_categories), ' categor', if (length(unique_categories) == 1L) "y" else "ies",
    ' in this collection. One example of each figure type is shown below. ',
    'Select an image to open its full-resolution figure, or select a category to explore all its figures.</p>',
    '<div class="links">', paste(family_pills, collapse = ""), '</div>',
    if (length(related_sections)) paste0('<p>More extended figures in this collection:</p><div class="links">',
      paste(vapply(related_sections,function(label) sprintf('<a class="file-link" href="#%s">%s</a>',
        esc(slug(label)),esc(sub("^.* - ","",label))),character(1)),collapse=""),'</div>') else "",
    if (length(native_layers)) paste0('<div class="full-native-previews image-grid compact">',
      paste(vapply(names(native_layers),function(label) {
        spec <- native_layers[[label]]
        if(!length(spec$files)) return("")
        paste0('<div data-native-family="',esc(label),'">',
          '<h3>',esc(label),'</h3>',figure_card(spec$files[[1]],page_file),'</div>')
      },character(1)),collapse=""),'</div>') else "",
    if (length(preview)) paste0('<div class="full-product-previews">',
      paste(preview, collapse = ""), '</div>') else "",
    '<p class="muted">Open the figures for one category directly:</p>',
    '<div class="links">', category_links,
    if (length(unique_categories) > length(shown)) sprintf(
      '<span class="pill">%d further categories</span>',
      length(unique_categories) - length(shown)) else "",
    '<a class="file-link" href="', esc(deeplink(NULL)), '">Open the category explorer</a>',
    '</div></section>')
}

# Full-report attachments are presentation-only.  A saved product enters a
# category card only when its own source/matrix sidecar states one exact
# category_id; directory scope supplies the already selected context and
# collection.  Filename text never decides category membership.
report_category_product_family <- function(path) {
  stem <- sub("_(source|matrix)[.]tsv$|_settings[.]json$|_recipe[.]R$|[.](png|svg|pdf)$", "",
    path, ignore.case = TRUE)
  # Heatmap images encode their display scale while their matrix/recipe does
  # not.  This is a format-family normalization, never category inference.
  sub("(_leading_edge_gene_heatmap)_(zscore|log2|none)$", "\\1", stem)
}

# Build the exact-ID attachment inventory without using filenames as evidence.
# The caller has already selected one context and collection through canonical
# output directories; the sidecar must still state exactly one category ID.
report_category_product_inventory <- function(specs) {
  valid_path <- function(path) file.exists(path) && !grepl("[[:cntrl:]]", path)
  products <- list(); unclassified <- character()
  for (spec in specs) {
    files <- if (dir.exists(spec$path)) list.files(spec$path, recursive = TRUE, full.names = TRUE) else character()
    files <- files[vapply(files, valid_path, logical(1L))]
    sidecars <- files[grepl("_(source|matrix)[.]tsv$", files, ignore.case = TRUE)]
    grouped <- split(files, report_category_product_family(files))
    classified <- character()
    for (sidecar in sidecars) {
      source <- tryCatch(read_tsv(sidecar), error = function(error) data.frame())
      if (!"category_id" %in% names(source)) next
      ids <- unique(as.character(source$category_id)); ids <- ids[!is.na(ids) & nzchar(ids)]
      if (length(ids) != 1L) next
      assets <- grouped[[report_category_product_family(sidecar)]]
      if (!length(assets)) next
      key <- paste(ids[[1L]], spec$product, sep = "\r")
      if (is.null(products[[key]])) {
        products[[key]] <- list(category_id = ids[[1L]], product = spec$product,
          label = spec$label, assets = sort(assets))
      } else products[[key]]$assets <- sort(unique(c(products[[key]]$assets, assets)))
      classified <- c(classified, assets)
    }
    unclassified <- c(unclassified, setdiff(files, unique(classified)))
  }
  list(category_products = unname(products), unclassified = sort(unique(unclassified)))
}

report_evidence_attachment_scope <- function(evidence_index, owner, collection, contrast = FALSE) {
  metadata_path <- file.path(dirname(evidence_index), "tables", "metadata.tsv")
  categories_path <- file.path(dirname(evidence_index), "tables", "categories.tsv")
  metadata <- tryCatch(read_tsv(metadata_path), error = function(error) data.frame())
  categories <- tryCatch(read_tsv(categories_path), error = function(error) data.frame())
  owner_key <- if (contrast) "contrast_id" else "analysis_id"
  if (!all(c("key", "value") %in% names(metadata)) || anyDuplicated(metadata$key) ||
      !"category_id" %in% names(categories)) {
    return(list(valid = FALSE, category_ids = character(), reason = "evidence metadata or category table is incomplete"))
  }
  values <- stats::setNames(as.character(metadata$value), as.character(metadata$key))
  declared_ids <- as.character(categories$category_id)
  declared_ids <- declared_ids[!is.na(declared_ids) & nzchar(declared_ids)]
  ids <- unique(declared_ids)
  if (!all(c(owner_key, "collection") %in% names(values)) ||
      !identical(unname(values[[owner_key]]), as.character(owner)) ||
      !identical(unname(values[["collection"]]), as.character(collection)) || !length(ids) || anyDuplicated(declared_ids)) {
    return(list(valid = FALSE, category_ids = character(), reason = "evidence context does not match the selected owner and collection"))
  }
  list(valid = TRUE, category_ids = ids, reason = "")
}

report_move_category_products_to_end <- function(html) {
  # Saved evidence from an earlier package can carry the former early panel
  # placement. Normalize that installed-route input before writing the updated
  # payload, so retained FULL extras follow—not precede—the standard evidence.
  pattern <- '<section id="category-products"[^>]*>.*?</section>'
  match <- regexpr(pattern, html, perl = TRUE)
  if (match[[1L]] < 1L) return(html)
  product_section <- regmatches(html, match)
  without_products <- paste0(
    substr(html, 1L, match[[1L]] - 1L),
    substr(html, match[[1L]] + attr(match, "match.length")[[1L]], nchar(html))
  )
  main_close <- regexpr("</main>", without_products, fixed = TRUE)[[1L]]
  if (main_close < 1L) stop("Evidence page has no closing main element for category products.", call. = FALSE)
  paste0(substr(without_products, 1L, main_close - 1L), product_section,
    substr(without_products, main_close, nchar(without_products)))
}

# The saved-products panel stays at the end of the evidence page: the standard
# evidence must be read first, and the presentation-layout tests pin that
# order. The consequence is that a reader arriving at the page has no way of
# knowing it carries FULL figures at all -- they sit ~96% of the way down.
# Inject a small static announcement immediately after the page intro instead
# of moving the panel: plain HTML with a real fragment link, so it survives
# `--dump-dom` and still works with JavaScript disabled. Idempotent by
# construction -- any previously injected block is removed before a new one is
# written -- and an empty payload leaves the page exactly as it was.
report_evidence_products_anchor_pattern <-
  '(?s)<div id="evidence-products-anchor"[^>]*>.*?</div>'

report_inject_evidence_products_anchor <- function(html, products) {
  html <- gsub(report_evidence_products_anchor_pattern, "", html, perl = TRUE)
  if (!is.list(products) || !length(products)) return(html)
  categories <- unique(vapply(products, function(product)
    as.character(product$category_id)[[1L]], character(1L)))
  block <- paste0(
    '<div id="evidence-products-anchor" class="evidence-products-anchor"',
    ' data-lisa-products-total="', length(products), '"',
    ' data-lisa-products-categories="', length(categories), '">',
    '<p><strong>FULL figures on this page.</strong> ',
    '<span data-lisa-products-summary>', length(products),
    ' saved figure product', if (length(products) == 1L) "" else "s",
    ' across ', length(categories), ' categor',
    if (length(categories) == 1L) "y" else "ies", '.</span> ',
    '<span data-lisa-products-current></span> ',
    '<a class="evidence-products-jump" href="#category-products">',
    'Jump to the saved figures for this category</a></p></div>')
  marker <- '<div class="page-intro">'
  intro <- regexpr(marker, html, fixed = TRUE)[[1L]]
  if (intro < 1L) {
    # A page with no intro is not a rendered evidence sheet (both installed
    # templates have one) -- it is a reduced input such as a test stub or a
    # pre-template archive. There is nothing to anchor to, so leave it exactly
    # as it is. But a page that DOES carry the saved-products panel and yet has
    # no intro is a real evidence sheet built from an unexpected template:
    # fail closed rather than ship it silently missing its announcement.
    if (regexpr('id="category-products"', html, fixed = TRUE)[[1L]] > 0L) {
      stop("Evidence page carries saved products but has no page intro to anchor them.",
           call. = FALSE)
    }
    return(html)
  }
  tail <- substr(html, intro + nchar(marker), nchar(html))
  # Neither installed template nests a <div> inside the page intro; verify that
  # rather than assuming it, so a future template change fails loudly instead
  # of silently anchoring the block in the wrong place.
  close <- regexpr("</div>", tail, fixed = TRUE)[[1L]]
  open <- regexpr("<div", tail, fixed = TRUE)[[1L]]
  if (close < 1L || (open > 0L && open < close)) {
    stop("Evidence page intro is not a flat block; cannot anchor saved products.", call. = FALSE)
  }
  at <- intro + nchar(marker) + close + 5L
  paste0(substr(html, 1L, at - 1L), block, substr(html, at, nchar(html)))
}

# Refresh presentation labels in HTML markup only, never scientific JSON.
report_script_link_labels <- function(markup) {
  if (!report_requested("recipes", TRUE)) {
    markup <- gsub('<a\\b[^>]*href=["\'][^"\']+[.]R([?#][^"\']*)?["\'][^>]*>[^<]*</a>',
      "", markup, perl = TRUE, ignore.case = TRUE)
  }
  gsub('(<a\\b[^>]*href=["\'][^"\']+[.]R([?#][^"\']*)?["\'][^>]*>)[^<]*(</a>)',
    "\\1R script\\3", markup, perl = TRUE, ignore.case = TRUE)
}

report_replace_evidence_payload <- function(evidence_index, products, data_id = "evidence-data") {
  html <- paste(readLines(evidence_index, warn = FALSE), collapse = "\n")
  marker <- paste0('<script type="application/json" id="', data_id, '">')
  start <- regexpr(marker, html, fixed = TRUE)[[1L]]
  if (start < 1L) stop("Evidence page has no data payload to attach products.", call. = FALSE)
  json_start <- start + nchar(marker)
  tail <- substr(html, json_start, nchar(html))
  end <- regexpr("</script>", tail, fixed = TRUE)[[1L]]
  if (end < 1L) stop("Evidence page has an unterminated data payload.", call. = FALSE)
  payload <- jsonlite::fromJSON(substr(tail, 1L, end - 1L), simplifyVector = FALSE)
  if (is.null(products)) products <- payload$category_products
  payload$category_products <- products
  # An old standalone sheet may lack flags. The containing run is the final
  # authority on which presentation exports the user requested.
  payload$source_data <- report_requested("source_data", TRUE)
  payload$recipes <- report_requested("recipes", TRUE) &&
    (is.null(payload$recipes) || isTRUE(payload$recipes))
  if (is.list(payload$member_figure_index)) {
    payload$member_figure_index <- lapply(payload$member_figure_index, function(entry) {
      if (!payload$source_data) entry$source_tsv <- ""
      if (!payload$recipes) entry$recipe_r <- ""
      entry
    })
  }
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", na = "null", null = "null", digits = 17)
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  replacement <- paste0(report_script_link_labels(substr(html, 1L, json_start - 1L)), json, "</script>",
    report_script_link_labels(substr(tail, end + nchar("</script>"), nchar(tail))))
  replacement <- report_move_category_products_to_end(replacement)
  replacement <- report_inject_evidence_products_anchor(replacement, products)
  lisaR:::lisa_guarded_write(evidence_index, function(target) writeLines(replacement, target, useBytes = TRUE))
  invisible(products)
}

# Attach saved full products to the existing evidence sheet, never to a second
# gallery or copy layer. Every href resolves from the page to the one canonical
# artifact already present under the run's `outputs/` or `artifacts/` root.
report_attach_category_products <- function(evidence_index, specs, owner, collection, contrast = FALSE,
                                            exclude_products = character(), project_dir = NULL) {
  inventory <- report_category_product_inventory(specs)
  # Member-gene-set charts are the canonical standard evidence already present
  # in every category sheet. A FULL inventory can carry another rendering of
  # that same family; exclude it before staging so it cannot become a duplicate
  # card beside the standard evidence.
  exclude_products <- unique(as.character(exclude_products))
  exclude_products <- exclude_products[!is.na(exclude_products) & nzchar(exclude_products)]
  if (length(exclude_products) && length(inventory$category_products)) {
    inventory$category_products <- inventory$category_products[!vapply(
      inventory$category_products,
      function(product) product$product %in% exclude_products,
      logical(1L)
    )]
  }
  inventory$unclassified <- Filter(report_asset_requested, inventory$unclassified)
  if (!file.exists(evidence_index)) return(inventory)
  scope <- report_evidence_attachment_scope(evidence_index, owner, collection, contrast)
  if (!scope$valid) {
    inventory$unclassified <- sort(unique(c(inventory$unclassified,
      unlist(lapply(inventory$category_products, `[[`, "assets"), use.names = FALSE))))
    inventory$category_products <- list()
    return(inventory)
  }
  accepted <- inventory$category_products[vapply(inventory$category_products, function(product)
    product$category_id %in% scope$category_ids, logical(1L))]
  rejected <- inventory$category_products[!vapply(inventory$category_products, function(product)
    product$category_id %in% scope$category_ids, logical(1L))]
  inventory$unclassified <- Filter(report_asset_requested, sort(unique(c(inventory$unclassified,
    unlist(lapply(rejected, `[[`, "assets"), use.names = FALSE)))))
  evidence_dir <- dirname(evidence_index)
  if (is.null(project_dir)) {
    project_dir <- get0("report_project_dir", inherits = TRUE, ifnotfound = NULL)
  }
  if (is.null(project_dir) || length(project_dir) != 1L || is.na(project_dir) ||
      !nzchar(project_dir)) {
    stop("Category-product assembly requires the canonical report root.",
         call. = FALSE)
  }
  payload <- lapply(accepted, function(product) list(category_id = product$category_id,
    product = product$product, label = product$label, assets = lapply(Filter(report_asset_requested, product$assets), function(path) {
      # Native builders may predate the installed reproduction closure. Keep
      # their immutable outputs and source recipes, but serve the verified
      # package recipe for native products. Extension artifacts keep their
      # original, independently verified recipes unchanged.
      relative <- report_relative_path(path, project_dir)
      if (product$product %in% c("gene_cards", "contrast_gene_cards") &&
          grepl("^outputs/gene_level/", relative) &&
          grepl("_recipe[.]R$", basename(path))) {
        recipe_dir <- file.path(project_dir, "outputs", "report_recipes")
        dir.create(recipe_dir, recursive = TRUE, showWarnings = FALSE)
        renderer <- file.path(report_package_dir, "scripts", "reproduce_lisa_figure.R")
        if (!file.exists(renderer)) renderer <- file.path(report_package_dir,
          "inst", "scripts", "reproduce_lisa_figure.R")
        recipe <- file.path(recipe_dir, paste0(figure_id_for_path(path, project_dir), "_recipe.R"))
        lisaR:::lisa_copy_verified_figure_recipe(renderer, recipe,
          report_renderer_sha256, run_root = project_dir)
        path <- recipe
      }
      lisaR:::lisa_category_product_reference_asset(
        project_dir, evidence_dir, path,
        allowed_roots = c("outputs", "artifacts"),
        label = "FULL category product"
      )
    })))
  report_replace_evidence_payload(evidence_index, payload,
    if (contrast) "contrast-evidence-data" else "evidence-data")
  list(category_products = payload, unclassified = inventory$unclassified)
}

# A category-scoped FULL product may replace an otherwise duplicate top-level
# gallery only when at least one exact-ID card has the complete displayed
# contract. If attachment/scope/source/recipe is incomplete, the caller keeps
# the ordinary required layer so the final missing-layer gate still fails.
report_has_complete_category_product <- function(attachment, product,
                                                 required_formats = c(report_image_formats(),
                                                   if (report_requested("source_data", TRUE)) "tsv",
                                                   if (report_requested("recipes", TRUE)) "r")) {
  products <- attachment$category_products
  if (!is.list(products) || !length(products)) return(FALSE)
  matching <- products[vapply(products, function(entry)
    is.list(entry) && identical(as.character(entry$product), product), logical(1L))]
  if (!length(matching)) return(FALSE)
  all(vapply(matching, function(entry) {
    assets <- entry$assets
    if (!is.list(assets)) return(FALSE)
    if (!length(required_formats)) return(TRUE)
    if (!length(assets)) return(FALSE)
    formats <- vapply(assets, function(asset) {
      if (!is.list(asset) || is.null(asset$format) || length(asset$format) != 1L) return(NA_character_)
      tolower(as.character(asset$format))
    }, character(1L))
    !anyNA(formats) && all(required_formats %in% formats)
  }, logical(1L)))
}

# A FULL run's `report_mode` is one report-wide switch, but the selected
# extras (from a fresh full-pipeline run or a reused saved-figure
# `artifacts/` subtree) may legitimately cover only some analyses/contrasts.
# An owner never selected for extras must be treated exactly like a STANDARD
# owner (its optional galleries are "not_requested", never a hard "missing"
# failure). An owner that IS selected must fail closed as "missing" if its
# products never materialized -- even when that means no directory exists
# for it at all (e.g. a bug, a partial pipeline failure, or a reassembly
# that only hard-links some owners' artifacts while declaring intent for
# others). `contract_manifest.tsv`'s `full_scope_analyses`/
# `full_scope_contrasts` keys (written once per run by config_pipeline.R)
# are the authoritative selection record when present, because directory
# existence alone cannot distinguish "never selected" from "selected but
# its output never showed up" -- both look identical on disk (no
# directory). When the contract does not declare these keys (an older
# project tree from before this record existed, or genuinely no
# information), fall back to the previous directory-existence heuristic
# unchanged, for full backward compatibility.
report_full_scope_owner <- function(project_dir, owner, contrast = FALSE) {
  # A selected report is a first-class render view, not a modified scientific
  # run contract. The entry point records the exact scope before assembly.
  selected_file <- file.path(project_dir, "report_selection.tsv")
  if (file.exists(selected_file)) {
    selected <- read_tsv(selected_file)
    key <- if (contrast) "contrast_owners" else "analyses"
    value <- selected$value[selected$key == key]
    if (length(value) != 1L || is.na(value)) stop("Invalid report selection.", call. = FALSE)
    return(owner %in% strsplit(value, ";", fixed = TRUE)[[1L]])
  }
  contract_file <- file.path(project_dir, "contract_manifest.tsv")
  declared_key <- if (contrast) "full_scope_contrasts" else "full_scope_analyses"
  if (file.exists(contract_file)) {
    ec <- read_tsv(contract_file)
    if (all(c("key", "value") %in% names(ec)) && any(ec$key == declared_key)) {
      declared_value <- as.character(ec$value[ec$key == declared_key][[1L]])
      declared_owners <- if (is.na(declared_value) || !nzchar(declared_value)) {
        character(0)
      } else {
        strsplit(declared_value, ";", fixed = TRUE)[[1L]]
      }
      return(owner %in% declared_owners)
    }
  }
  if (contrast) {
    dir.exists(file.path(project_dir, "artifacts", "contrasts", owner)) ||
      dir.exists(file.path(project_dir, "outputs", "gene_level", "category_contrasts", owner))
  } else {
    dir.exists(file.path(project_dir, "artifacts", owner)) ||
      dir.exists(file.path(project_dir, "outputs", "gene_level", "single_de", owner))
  }
}

# Same original bitmap; only transparent, keyboard-accessible links are added.
report_link_existing_category_figure <- function(image, href, source, settings,
    evidence_index, page_file, title) {
  if(!file.exists(evidence_index)) stop('Category evidence target is missing.')
  cfg <- jsonlite::read_json(settings,simplifyVector=TRUE)
  cats <- read_tsv(file.path(dirname(evidence_index),'tables/categories.tsv'))
  meta <- read_tsv(file.path(dirname(evidence_index),'tables/metadata.tsv'))
  if(!all(c('key','value')%in%names(meta))) stop('Missing evidence identity.')
  root <- report_root_for_page(page_file)
  key <- stable_text_hash(paste(digest::digest(file=image,algo='sha256'),
    digest::digest(file=source,algo='sha256'),digest::digest(file=settings,algo='sha256'),
    'original-head-links-v1',sep=':'))
  cache <- file.path(root,'report_figure_data/category_links',paste0(key,'.tsv'))
  if(!file.exists(cache)) {
    regions <- lisaR:::lisa_category_figure_regions(source,settings)
    dir.create(dirname(cache),recursive=TRUE,showWarnings=FALSE)
    utils::write.table(regions,cache,sep='\t',quote=TRUE,row.names=FALSE)
  } else regions <- utils::read.delim(cache,sep="\t",quote='"',check.names=FALSE,comment.char="")
  if(nrow(regions)&&(!'category_id'%in%names(cats)||any(!regions$category_id%in%cats$category_id)))
    stop('Original figure categories do not match evidence.')
  identity <- meta[meta$key%in%c('contrast_id','analysis_id','analysis_a','analysis_b','collection','tier'),]
  query <- paste(paste0(identity$key,'=',vapply(as.character(identity$value),utils::URLencode,
    character(1L),reserved=TRUE)),collapse='&')
  links <- vapply(seq_len(nrow(regions)),function(i) {
    r <- regions[i,];label <- cats$category_display_name[match(r$category_id,cats$category_id)]
    if(!length(label)||is.na(label))label<-r$category_id
    url <- paste0(rel_path(evidence_index,page_file),'?category_id=',
      utils::URLencode(r$category_id,reserved=TRUE),'&',query)
    sprintf('<a class="category-figure-hit category-figure-%s" data-category-id="%s" data-side="%s" href="%s" aria-label="%s" title="%s" style="position:absolute;left:%.9f%%;top:%.9f%%;width:%.9f%%;height:%.9f%%;border-radius:3px;box-sizing:border-box;"></a>',
      esc(r$kind),esc(r$category_id),esc(r$side),esc(url),esc(paste(label,r$side,'— evidence')),
      esc(paste(label,r$side,'— evidence')),100*r$left,100*r$top,100*r$width,100*r$height)
  },character(1L))
  paste0('<style>.category-figure-hit:hover,.category-figure-hit:focus-visible{outline:2px solid #16598a;background:#16598a18}.original-category-figure{line-height:0}.original-category-figure img{display:block!important;width:100%!important;height:auto!important;margin:0!important}</style>',
    '<div class="original-category-figure" data-original-sha256="',digest::digest(file=image,algo='sha256'),
    '" style="position:relative;width:100%;aspect-ratio:',cfg$width,'/',cfg$height,';">',
    '<img class="nes-variant-image" loading="lazy" decoding="async" src="',href,'" alt="',esc(title),'">',
    paste(links,collapse=''),'</div>')
}

report_existing_category_evidence <- function(collection_dir,evidence_index,page_file,title) {
  if (!length(report_image_formats())) return(paste0(
    '<p>No image exports were requested. The category and gene evidence remains available below.</p>',
    '<a class="file-link" href="', rel_path(evidence_index,page_file), '">Explore contrast categories and genes</a>'))
  diagram <- report_nes_variants_section(collection_dir,page_file,title,evidence_index,id_suffix='-evidence')
  if(!nzchar(diagram)) stop('Original LISA categories figure unavailable; cannot substitute a different plot.')
  paste0('<p>The original LISA categories figure is shown unchanged. Select a symbol (head) or category label to open its evidence.</p>',
    diagram,'<div class="links"><a class="file-link" href="',rel_path(evidence_index,page_file),
    '">Explore contrast categories and genes</a></div>')
}

# Replace a generated HTML element without reserializing other sections or
# embedded JSON. This bounded mode shares the normal assembly components.
report_replace_element <- function(html,start,tag,replacement) {
  tail <- substr(html,start,nchar(html))
  pattern <- paste0('</?',tag,'\\b[^>]*>')
  hits <- gregexpr(pattern,tail,perl=TRUE)[[1L]];lens<-attr(hits,'match.length')
  if(hits[1L]!=1L)stop('Element boundary is not exact.')
  depth<-0L;last<-NA_integer_
  for(i in seq_along(hits)) {
    token<-substr(tail,hits[i],hits[i]+lens[i]-1L)
    depth<-depth+if(startsWith(token,'</'))-1L else 1L
    if(depth==0L){last<-hits[i]+lens[i]-1L;break}
  }
  if(is.na(last))stop('Unclosed generated element.')
  paste0(substr(html,1L,start-1L),replacement,substr(tail,last+1L,nchar(tail)))
}

report_refresh_category_links <- function(project_dir) {
  page <- file.path(project_dir,'report_pages/contrasts.html')
  index <- read_tsv(file.path(project_dir,'config','contrast_index.tsv'))
  if(!nrow(index)) {message('No contrasts: no category-link change required.');return(invisible(character()))}
  html <- paste(readLines(page,warn=FALSE),collapse='\n')
  original <- html
  dirs <- list.dirs(file.path(project_dir,'outputs/category_contrasts'),recursive=FALSE,full.names=TRUE)
  for(owner in basename(dirs)) for(collection_dir in contrast_collection_dirs(project_dir,owner)) {
    collection <- sub('^collection_','',basename(collection_dir))
    evidence <- file.path(project_dir,'report_pages/contrast_evidence',owner,collection,'index.html')
    manifest <- list.files(file.path(collection_dir,'qc'),pattern='_category_nes_variants[.]tsv$',full.names=TRUE)
    if(length(manifest)!=1L)stop('Original NES manifest is required for refresh.')
    family <- paste0('N_',stable_text_hash(report_relative_path(manifest,project_dir)))
    # Original categories section: update only its variant container.
    needle <- paste0('<div class="nes-variants" data-nes-variants="',family,'">')
    pos<-regexpr(needle,html,fixed=TRUE)[1L]
    if(pos<1L)stop('Original LISA categories section not found.')
    title <- paste(owner,collection) # alt fallback; bitmap text stays untouched.
    replacement <- report_nes_variants_section(collection_dir,page,title,evidence)
    html <- report_replace_element(html,pos,'div',replacement)
    # Replace only the rejected evidence cover, preserving its navigation ID.
    starts<-gregexpr('<section class="panel contrast-evidence-cover" id="[^"]+">',html,perl=TRUE)[[1L]]
    found<-FALSE
    for(pos in starts) {
      tail<-substr(html,pos,nchar(html));end<-regexpr('</section>',tail,fixed=TRUE)[1L]+9L
      block<-substr(tail,1,end)
      if(!grepl(rel_path(evidence,page),block,fixed=TRUE))next
      opening<-regmatches(tail,regexpr('^<section[^>]+>',tail))
      replacement<-paste0(opening,'<h2>Contrast evidence</h2>',
        report_existing_category_evidence(collection_dir,evidence,page,title),'</section>')
      html<-report_replace_element(html,pos,'section',replacement);found<-TRUE;break
    }
    if(!found)stop('Contrast evidence cover not found.')
  }
  # Only page + region tables changed. Update their entries in the existing
  # manifest; retain every other row/hash, no recompression or full traversal.
  write_page(page,html)
  manifest <- read_shareable_report_manifest(project_dir)
  changed <- c('report_pages/contrasts.html',substring(list.files(file.path(project_dir,'report_figure_data/category_links'),
    full.names=TRUE),nchar(project_dir)+2L))
  for(rel in changed) {
    row<-data.frame(relative_path=rel,role=if(grepl('[.]html$',rel))'page' else 'figure_contract',
      bytes=file.info(file.path(project_dir,rel))$size,
      sha256=digest::digest(file=file.path(project_dir,rel),algo='sha256'),integrity='sha256')
    where<-match(rel,manifest$relative_path)
    if(is.na(where))manifest<-rbind(manifest,row) else manifest[where,]<-row
  }
  manifest<-manifest[order(manifest$relative_path),]
  utils::write.table(manifest,file.path(project_dir,report_shareable_manifest_name),sep='\t',quote=FALSE,row.names=FALSE,na='')
  message('CATEGORY_LINK_REFRESH_PASS: contrasts page and ',length(changed)-1L,' region tables; no scientific output changed.')
  invisible(changed)
}

# Runtime manifests name exact source/settings files. Keep their bytes together
# because NES reproduction verifies the stored source hash; do not run these
# through the generic, augmented figure-source recipe.
report_nes_variants_section <- function(collection_dir, page_file, title, evidence_index = NULL, id_suffix = "") {
  manifests <- list.files(file.path(collection_dir, "qc"),
    pattern = "_category_nes_variants[.]tsv$", full.names = TRUE)
  if (!length(manifests)) return(report_hommel_only_section(collection_dir, page_file, id_suffix))
  if (length(manifests) != 1L) stop("Ambiguous NES-variant manifest.", call. = FALSE)
  x <- read_tsv(manifests[[1L]])
  required <- c("variant", "status", "format", "path", "source", "settings", "recipe")
  if (!all(required %in% names(x))) stop("Invalid NES-variant manifest schema.", call. = FALSE)
  x <- x[x$status == "rendered", , drop = FALSE]
  if (!nrow(x)) return(report_hommel_only_section(collection_dir, page_file, id_suffix))
  # Archived manifests predate the independent contrast subset selector.
  if (!"plot_set" %in% names(x)) x$plot_set <- "all"
  if (any(!x$variant %in% c("clean", "percentages", "direction", "dispersion")) ||
      any(!x$plot_set %in% c("all", "same_direction", "opposite_direction")) ||
      any(!x$format %in% c("png", "svg", "pdf")) ||
      anyDuplicated(paste(x$plot_set, x$variant, x$format))) stop("Invalid NES-variant identities.", call. = FALSE)
  resolve <- function(relative) {
    if (is.na(relative) || !nzchar(relative) ||
        grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^[/\\\\]|(^|[/\\\\])[.][.]([/\\\\]|$)", relative))
      stop("NES-variant artifacts must be collection-relative.", call. = FALSE)
    target <- file.path(collection_dir, relative)
    if (!file.exists(target)) stop("Missing NES-variant artifact: ", relative, call. = FALSE)
    target
  }
  for (key in c("path", "source", "settings", "recipe")) x[[key]] <- vapply(x[[key]], resolve, character(1L))
  root <- report_root_for_page(page_file)
  family_id <- paste0("N_", stable_text_hash(report_relative_path(manifests[[1L]], root)))
  dom_id <- paste0(family_id, id_suffix)
  hommel_panel <- report_hommel_view_panel(collection_dir, page_file, dom_id)
  bundle <- file.path(root, "report_figure_data", "nes_variants", family_id)
  dir.create(bundle, recursive = TRUE, showWarnings = FALSE)
  copy_exact <- function(path, name = basename(path)) {
    target <- file.path(bundle, name)
    if (file.exists(target) && !identical(digest::digest(file = path, algo = "sha256"),
        digest::digest(file = target, algo = "sha256")))
      stop("Conflicting NES-variant artifact names.", call. = FALSE)
    if (!file.copy(path, target, overwrite = TRUE)) stop("Cannot copy NES-variant artifact.", call. = FALSE)
    if (!identical(digest::digest(file = path, algo = "sha256"), digest::digest(file = target, algo = "sha256")))
      stop("NES-variant copy failed integrity verification.", call. = FALSE)
    target
  }
  package <- get0("report_package_dir", inherits = TRUE, ifnotfound = "")
  recipe_candidates <- c(file.path(package, "scripts", "reproduce_lisa_category_nes.R"),
    file.path(package, "inst", "scripts", "reproduce_lisa_category_nes.R"))
  recipe <- recipe_candidates[file.exists(recipe_candidates)][1L]
  if (is.na(recipe)) stop("Installed NES reproduction recipe is missing.", call. = FALSE)
  recipe_copy <- if (report_requested("recipes", TRUE)) copy_exact(recipe, "reproduce_category_nes.R") else ""
  available <- c("clean", "percentages", "direction", "dispersion")
  available <- available[available %in% x$variant]
  plot_sets <- c("all", "same_direction", "opposite_direction")
  plot_sets <- plot_sets[plot_sets %in% x$plot_set]
  labels <- c(clean = "Clean", percentages = "Percentages", direction = "Mean + median", dispersion = "Percentiles")
  subset_labels <- c(all = "All", same_direction = "Same direction", opposite_direction = "Opposite direction")
  identities <- unique(x[, c("plot_set", "variant"), drop = FALSE])
  identities <- identities[order(match(identities$plot_set, plot_sets), match(identities$variant, available)), , drop = FALSE]
  panels <- vapply(seq_len(nrow(identities)), function(j) {
    variant <- identities$variant[[j]]; plot_set <- identities$plot_set[[j]]
    rows <- x[x$variant == variant & x$plot_set == plot_set, , drop = FALSE]
    source <- unique(rows$source); settings <- unique(rows$settings)
    if (length(source) != 1L || length(settings) != 1L) stop("NES variant has ambiguous data contracts.", call. = FALSE)
    cfg <- jsonlite::read_json(settings, simplifyVector = TRUE)
    cfg_plot_set <- if (is.null(cfg$plot_set)) "all" else cfg$plot_set
    if (!identical(cfg$variant, variant) || !identical(cfg_plot_set, plot_set) ||
        !identical(cfg$source_sha256, digest::digest(file = source, algo = "sha256")))
      stop("NES variant settings disagree with its source or identity.", call. = FALSE)
    source_copy <- if (report_requested("source_data", TRUE)) copy_exact(source) else ""
    settings_copy <- if (report_requested("recipes", TRUE)) copy_exact(settings) else ""
    formats <- c("png", "svg", "pdf")
    preview_format <- formats[formats %in% rows$format][1L]
    preview <- rows$path[rows$format == preview_format][[1L]]
    preview_href <- short_media_copy(preview, page_file)
    downloads <- paste(vapply(seq_len(nrow(rows)), function(i) paste0('<a class="file-link" download href="',
      short_media_copy(rows$path[[i]], page_file), '">', toupper(esc(rows$format[[i]])), '</a>'), character(1L)), collapse = "")
    paste0('<div id="', dom_id, '-', plot_set, '-', variant, '" class="nes-variant-panel" data-nes-plot-set="', plot_set,
      '" data-nes-panel="', variant, '"', if (j != 1L || nzchar(hommel_panel)) ' hidden' else "", '>',
      '<div class="nes-downloads"><strong>Download this view</strong><div class="links">', downloads,
      if (report_requested("source_data", TRUE)) paste0('<a class="file-link" download href="', rel_path(source_copy, page_file), '">Source data</a>') else "",
      if (report_requested("recipes", TRUE)) paste0('<a class="file-link" download href="', rel_path(settings_copy, page_file), '">View settings</a>',
        '<a class="file-link" download href="', rel_path(recipe_copy, page_file), '">R script</a>') else "",
      '</div></div>',
      if (preview_format %in% c("png", "svg") && !is.null(evidence_index)) report_link_existing_category_figure(
        preview, preview_href, source, settings, evidence_index, page_file,
        paste(title, labels[[variant]], subset_labels[[plot_set]], sep=" - ")) else
      if (preview_format %in% c("png", "svg")) paste0('<img class="nes-variant-image" loading="lazy" decoding="async" src="', preview_href, '" alt="',
        esc(paste(title, labels[[variant]], subset_labels[[plot_set]], sep = " - ")), '">') else '<p>Download this view as PDF.</p>',
      '</div>')
  }, character(1L))
  paste0('<div class="nes-variants" data-nes-variants="', dom_id, '"><div class="nes-controls">',
    '<label for="', dom_id, '-select">View <select id="', dom_id, '-select" data-nes-select>',
    if (nzchar(hommel_panel)) '<option value="hommel_support" selected>Category significance</option>' else "",
    paste(vapply(available, function(v) paste0('<option value="', v, '">', labels[[v]], '</option>'), character(1L)), collapse = ""),
    '</select></label>',
    if (length(plot_sets) > 1L) paste0('<label for="', dom_id, '-subset">Categories <select id="', dom_id, '-subset" data-nes-subset>',
      paste(vapply(plot_sets, function(v) paste0('<option value="', v, '">', subset_labels[[v]], '</option>'), character(1L)), collapse = ""), '</select></label>') else "",
    '</div><p class="muted nes-view-description" aria-live="polite">', if (nzchar(hommel_panel)) 'Stars show adjusted category P from robust Hommel multiple-testing correction; minimum support d/N is in Category evidence. Mean NES is descriptive.' else 'Mean NES is retained in every view. Percentages use significant member sets; P25-P75 describes dispersion, not a confidence interval.', '</p>',
    hommel_panel, paste(panels, collapse = ""), '</div>')
}

# New presentation product, not a mutation of historical NES-variant manifests.
report_hommel_view_panel <- function(collection_dir, page_file, dom_id) {
  if (identical(basename(collection_dir), "collection_HALLMARKS") ||
      !grepl("/single_de/", gsub("\\\\", "/", collection_dir), fixed = TRUE)) return("")
  formats <- report_image_formats()
  if (!length(formats)) return("")
  paths <- list.files(file.path(collection_dir, "plots"),
    pattern = "_GSEA_lollipop_hommel_support[.](png|svg|pdf)$", full.names = TRUE)
  paths <- paths[tools::file_ext(paths) %in% formats]
  if (!length(paths)) return("")
  stems <- unique(sub("[.](png|svg|pdf)$", "", paths))
  if (length(stems) != 1L) stop("Ambiguous Hommel-support figure.", call. = FALSE)
  path <- paths[order(match(tools::file_ext(paths), c("png", "svg", "pdf")))][[1L]]
  paste0('<div id="', dom_id, '-all-hommel_support" class="nes-variant-panel hommel-support-view" ',
    'data-nes-plot-set="all" data-nes-panel="hommel_support">',
    '<p class="hommel-support-subtitle"><strong>Category significance &middot; adjusted P (robust Hommel)</strong></p>',
    lisaR:::lisa_hommel_help_button(),
    '<p class="muted">Stars show adjusted category P: * &le; 0.05, ** &le; 0.01, *** &le; 0.001, **** &le; 0.0001; no star means not significant or not evaluable. ',
    'See Category evidence for the exact P and minimum support d/N at 95% simultaneous confidence. ',
    'The separate NES panel describes the individually significant sets.</p>',
    figure_card(path, page_file, "Category significance"), lisaR:::lisa_hommel_reference_html(), '</div>')
}

report_hommel_only_section <- function(collection_dir, page_file, id_suffix = "") {
  id <- paste0("hommel-", stable_text_hash(report_relative_path(collection_dir, report_root_for_page(page_file))), id_suffix)
  panel <- report_hommel_view_panel(collection_dir, page_file, id)
  if (!nzchar(panel)) return("")
  paste0('<div class="nes-variants" data-nes-variants="', id, '"><div class="nes-controls"><label for="', id,
    '-select">View <select id="', id, '-select" data-nes-select><option value="hommel_support" selected>Category significance</option>',
    '</select></label></div>', panel, '</div>')
}

# Support decoration is a separate, neutral column. Keep both original
# historical lollipop variants downloadable without conflating their
# significant-member denominator with the category-inference denominator.
report_support_lollipops <- function(collection_dir, page_file) {
  formats <- report_image_formats()
  if (!length(formats)) return("")
  sources <- list.files(file.path(collection_dir, "plots"), pattern = "_GSEA_lollipop(_direction_stats)?_support_source[.]tsv$")
  if (!length(sources)) return('<p class="support-lollipop-unavailable muted">Support lollipop unavailable: this run has no complete historical lollipop source. The complete support table and member-set NES figure remain available.</p>')
  card <- function(suffix, css, title, expanded = FALSE) {
    paths <- list.files(file.path(collection_dir, "plots"),
      pattern = paste0(suffix, "[.](png|svg|pdf)$"), full.names = TRUE)
    paths <- paths[tools::file_ext(paths) %in% formats]
    if (!length(paths)) {
      source <- list.files(file.path(collection_dir, "plots"), pattern = paste0(suffix, "_source[.]tsv$"))
      if (length(source)) stop("Required support lollipop is missing: ", suffix, call. = FALSE)
      return("")
    }
    path <- paths[order(match(tools::file_ext(paths), c("png", "svg", "pdf")))][[1L]]
    paste0('<details class="', css, '"', if (expanded) ' open' else '', '><summary>',
      esc(title), '</summary>', figure_card(path, page_file, title), '</details>')
  }
  paste0('<details class="support-historical-lollipops"><summary>Historical lollipops (unchanged)</summary>',
    '<p>These earlier views remain available. Their NES points and sizes describe individually significant gene sets. ',
    'The S1 to S5 support column uses all evaluated sets, including those that do not pass the individual cutoff. ',
    'In the direction view, percentages describe the individually significant sets, not all evaluated sets.</p>',
    card("_GSEA_lollipop_support", "support-lollipop", "Historical lollipop with minimum support"),
    card("_GSEA_lollipop_direction_stats_support", "support-direction-lollipop", "Lollipop with minimum support and observed direction"),
    card("_GSEA_lollipop", "historical-lollipop", "Historical mean NES lollipop"),
    card("_GSEA_lollipop_direction_stats", "historical-direction-lollipop", "Historical lollipop with direction statistics"),
    '</details>')
}

# A single reading surface replaces the historical lollipop/dumbbell gallery
# when versioned variants exist. The pathway dotplot is a gene-set product,
# not a duplicate category summary; optional ORA retains its own identity.
report_category_inference_downloads <- function(collection_dir, page_file) {
  if (identical(basename(collection_dir), "collection_HALLMARKS")) return("")
  paths <- list.files(file.path(collection_dir, "lisa_tables"), pattern = "_LISA_category_inference[.]tsv$", full.names = TRUE)
  if (!length(paths)) return("")
  if (length(paths) != 1L) stop("Ambiguous category-inference table.", call. = FALSE)
  x <- read_tsv(paths[[1L]])
  images <- list.files(file.path(collection_dir, "plots"), pattern = "_LISA_category_inference[.](png|svg|pdf)$", full.names = TRUE)
  formats <- report_image_formats()
  images <- images[tools::file_ext(images) %in% formats]
  if (any(x$n_sets_evaluable > 0L) && length(formats) && !length(images))
    stop("Required category-inference figure is missing.", call. = FALSE)
  image <- if (length(images)) images[order(match(tools::file_ext(images), c("png", "svg", "pdf")))][[1L]] else ""
  figure_links <- if (nzchar(image)) {
    stem <- tools::file_path_sans_ext(image)
    if (any(!file.exists(paste0(stem, ".", formats))))
      stop("Requested category-inference figure format is missing.", call. = FALSE)
    contract <- figure_source_contract(image, page_file)
    paste0('<span data-figure-kind="member-set-nes-distribution">Member-set NES distribution (not the LISA category plot):</span>',
      paste(vapply(formats, function(format) paste0('<a class="file-link" download href="',
        esc(short_media_copy(paste0(stem, ".", format), page_file)), '">', toupper(format), '</a>'), character(1L)), collapse = ""),
      if (nzchar(contract$source_tsv)) paste0('<a class="file-link" download href="', esc(contract$source_tsv), '">Figure data</a>') else "",
      if (nzchar(contract$recipe_r)) paste0('<a class="file-link" download href="', esc(contract$recipe_r), '">Figure R script</a>') else "")
  } else ""
  table_link <- if (report_requested("source_data", TRUE)) paste0('<a class="file-link" data-category-results-tsv download href="',
    esc(short_file_copy(paths[[1L]], page_file)), '">Complete category results (TSV)</a>') else ""
  if (!nzchar(table_link) && !nzchar(figure_links)) return("")
  paste0('<div class="links category-results-downloads" data-download-context="category-evidence">', table_link, figure_links, '</div>')
}

report_lisa_categories_section <- function(title, collection_dir, files, page_file,
    status_html = "", note = "", evidence_index = NULL) {
  has_inference <- !identical(basename(collection_dir), "collection_HALLMARKS") &&
    length(list.files(file.path(collection_dir, "lisa_tables"), pattern = "_LISA_category_inference[.]tsv$")) > 0L
  historical <- if (has_inference) report_support_lollipops(collection_dir, page_file) else ""
  downloads <- if (is.null(evidence_index) || !file.exists(evidence_index))
    report_category_inference_downloads(collection_dir, page_file) else ""
  variants <- report_nes_variants_section(collection_dir, page_file, title, evidence_index)
  if (!nzchar(variants)) return(paste0(layer_section(title, files, page_file, note, status_html,
    empty_state = if (!any(vapply(c("png", "svg", "pdf"), report_requested, logical(1L), default = TRUE)))
      "not_requested" else "missing"), historical, downloads))
  complementary <- files[grepl("pathway_dotplot", basename(files), ignore.case = TRUE)]
  ora <- files[grepl("ORA", basename(files), fixed = TRUE)]
  gallery <- function(paths, label, description) {
    if (!length(paths)) return("")
    paste0('<details class="category-complement"><summary>', esc(label), '</summary><p class="muted">',
      esc(description), '</p><div class="image-grid">',
      paste(vapply(paths, figure_card, character(1L), page_file = page_file), collapse = ""), '</div></details>')
  }
  paste0('<section id="', esc(slug(title)), '" class="panel lisa-categories-section"><div class="section-head">',
    '<div><h2>LISA categories</h2><p class="muted">Switch views without changing the category order or scale. ',
    'Supercategory bands group the categories in the figures and downloads.</p></div>', status_html, '</div>',
    variants,
    gallery(complementary, "Gene-set overview", "A complementary view of individual gene sets, rather than category averages."),
    gallery(ora, "Over-representation analysis (ORA)", "Optional ORA figures; separate from the GSEA category summaries."),
    historical, downloads, '</section>')
}

completed_like_statuses <- c("completed", "skipped_existing", "done", "success", "ok")

status_count <- function(path, status_col = "status") {
  df <- read_tsv(path)
  if (nrow(df) == 0 || !status_col %in% names(df)) return(c(total = 0, completed = 0, failed = 0))
  status <- tolower(as.character(df[[status_col]]))
  c(
    total = nrow(df),
    completed = sum(status %in% completed_like_statuses),
    failed = sum(status == "failed")
  )
}

status_count_any <- function(paths, status_col = "status") {
  for (path in paths) {
    df <- read_tsv(path)
    if (nrow(df) > 0 && status_col %in% names(df)) {
      status <- tolower(as.character(df[[status_col]]))
      return(c(
        total = nrow(df),
        completed = sum(status %in% completed_like_statuses),
        failed = sum(status == "failed")
      ))
    }
  }
  c(total = 0, completed = 0, failed = 0)
}

biological_de_status <- function(de_index) {
  if (is.null(de_index) || nrow(de_index) == 0) return(c(total = 0, completed = 0, failed = 0))
  total <- nrow(de_index)
  if ("de_path" %in% names(de_index)) {
    present <- vapply(de_index$de_path, file.exists, logical(1))
    return(c(total = total, completed = sum(present), failed = sum(!present)))
  }
  c(total = total, completed = total, failed = 0)
}

biological_contrast_status <- function(contrast_index, contrast_status_df = data.frame()) {
  if (is.null(contrast_index) || nrow(contrast_index) == 0) return(c(total = 0, completed = 0, failed = 0))
  total <- nrow(contrast_index)
  if (nrow(contrast_status_df) == 0 || !"status" %in% names(contrast_status_df)) {
    return(c(total = total, completed = total, failed = 0))
  }
  id_col <- if ("contrast_id" %in% names(contrast_status_df) && "contrast_id" %in% names(contrast_index)) {
    "contrast_id"
  } else if ("output_id" %in% names(contrast_status_df) && "output_id" %in% names(contrast_index)) {
    "output_id"
  } else {
    ""
  }
  if (!nzchar(id_col)) return(c(total = total, completed = total, failed = 0))
  completed <- 0L
  failed <- 0L
  for (id in as.character(contrast_index[[id_col]])) {
    rows <- contrast_status_df[as.character(contrast_status_df[[id_col]]) == id, , drop = FALSE]
    status <- tolower(as.character(rows$status))
    if (length(status) > 0 && all(status %in% completed_like_statuses)) completed <- completed + 1L
    if (length(status) > 0 && any(status == "failed")) failed <- failed + 1L
  }
  c(total = total, completed = completed, failed = failed)
}

read_first_status_table <- function(paths, status_col = "status") {
  for (path in paths) {
    df <- read_tsv(path)
    if (nrow(df) > 0 && status_col %in% names(df)) return(df)
  }
  data.frame()
}

pngs <- function(path, pattern = "[.]png$") {
  files <- lisaR:::lisa_plot_files(path, pattern, formats = report_image_formats())
  files[!is_report_excluded_path(files)]
}

is_report_excluded_path <- function(path) {
  grepl(
    "GSEA_pathway_dotplot|pathway_dotplot|kegg_map_layer|contrast_kegg_map_layer|kegg_gene_heatmap",
    path,
    ignore.case = TRUE
  )
}

has_kegg_report_outputs <- function(project_dir) {
  length(pngs(
    file.path(project_dir, "outputs", "gene_level", "category_contrasts"),
    "contrast_painted\\.png$"
  )) > 0
}

kegg_report_state <- function(project_dir, kind, owner_id, collection, painter_dir, files) {
  stage <- if (identical(kind, "single")) "single_de_kegg_pathway_painter" else "contrast_kegg_pathway_painter"
  status <- read_tsv(file.path(project_dir, "post_lisa_status.tsv"))
  rows <- status[status$stage == stage & status$collection == collection, , drop = FALSE]
  if (identical(kind, "single") && "analysis_id" %in% names(rows)) rows <- rows[rows$analysis_id == owner_id, , drop = FALSE]
  if (identical(kind, "contrast") && "contrast_id" %in% names(rows)) rows <- rows[rows$contrast_id == owner_id, , drop = FALSE]
  index_files <- list.files(painter_dir, pattern = "_kegg_pathway_painter_index[.]tsv$", full.names = TRUE)
  index <- if (length(index_files)) read_tsv(index_files[[1]]) else data.frame()
  node_files <- list.files(painter_dir, pattern = "_kegg_pathway_painter_nodes[.]tsv$", full.names = TRUE)
  status_value <- if (nrow(rows) && "status" %in% names(rows)) tolower(as.character(rows$status[[nrow(rows)]])) else ""
  reason <- if (nrow(rows) && "message" %in% names(rows)) as.character(rows$message[[nrow(rows)]]) else ""
  if (length(files)) {
    if (!length(index_files) || !length(node_files)) stop("LISA-REPORT-KEGG-001 rendered KEGG map lacks required index or node layer: ", painter_dir, call. = FALSE)
    output_col <- if ("output_png" %in% names(index)) "output_png" else ""
    indexed <- if (nzchar(output_col)) as.character(index[[output_col]]) else character()
    indexed <- indexed[nzchar(indexed)]
    local_indexed <- file.path(painter_dir, basename(indexed))
    if (!nzchar(output_col) || !length(indexed) ||
        !all(file.exists(files)) ||
        !setequal(tools::file_path_sans_ext(basename(files)), tools::file_path_sans_ext(basename(local_indexed)))) {
      stop("LISA-REPORT-KEGG-002 KEGG index has broken rendered-map links: ", painter_dir, call. = FALSE)
    }
    return(list(state = "rendered", reason = "Rendered KEGG maps have index, node layer and local output links."))
  }
  if (!length(report_image_formats())) return(list(state = "not_requested", reason = "No image exports were requested."))
  if (length(index_files) && nrow(index) == 0) return(list(state = "inapplicable", reason = "No significant KEGG pathway candidates were selected by the existing enrichment evidence."))
  if (status_value %in% c("skipped_no_kegg_rows", "skipped_species_mismatch", "skipped_missing_input", "not_applicable")) {
    return(list(state = "skipped", reason = if (nzchar(reason)) reason else status_value))
  }
  if (status_value %in% c("completed", "failed") || length(index_files)) {
    stop("LISA-REPORT-KEGG-003 KEGG layer was expected but is neither rendered nor explicitly inapplicable: ", painter_dir, call. = FALSE)
  }
  list(state = "inapplicable", reason = "KEGG maps were not requested for this analysis/collection.")
}

kegg_layer_section <- function(title, files, page_file, note, state, status_html = "", subtype_html = "") {
  if (identical(state$state, "rendered")) return(section_figures(title, files, page_file, max(length(files), 1), note, status_html, subtype_html))
  section_id <- slug(title)
  sprintf('<section id="%s" class="panel"><div class="section-head"><div><p class="eyebrow">%s</p><h2>%s</h2><p class="muted">%s</p></div><span class="pill status-warn">%s</span></div></section>',
    esc(section_id), esc(title), esc(title), esc(state$reason), esc(state$state))
}

section_figures <- function(title, files, page_file, limit = 24, note = "", status_html = "", subtype_html = "") {
  files <- files[file.exists(files)]
  if (length(files) == 0) return("")
  section_id <- slug(title)
  cards <- paste(vapply(head(files, limit), function(f) figure_card(f, page_file, basename(f)), character(1)), collapse = "\n")
  shown <- min(length(files), limit)
  truncated <- if (length(files) > limit) sprintf('<span class="pill status-warn">%d hidden by page limit</span>', length(files) - limit) else ""
  sprintf(
    '<section id="%s" class="panel"><div class="section-head"><div><h2>%s</h2>%s</div><div>%s<span class="pill">%d/%d figures shown</span>%s</div></div>%s<div class="image-grid compact">%s</div></section>',
    esc(section_id), esc(sub("^.* - ", "", title)),
    if (nzchar(note)) paste0('<p class="muted">', esc(note), '</p>') else "",
    status_html, shown, length(files), truncated, subtype_html, cards
  )
}

missing_layer_section <- function(title, page_file, note = "No reportable files were found for this expected layer.") {
  section_id <- slug(title)
  parts <- strsplit(title, " - ", fixed = TRUE)[[1L]]
  family <- tail(parts, 1L)
  context <- paste(head(parts, -1L), collapse = " - ")
  sprintf(
    '<section id="%s" class="panel" data-lisa-layer-family="%s" data-lisa-layer-context="%s"><div class="section-head"><div><p class="eyebrow">%s</p><h2>%s</h2><p class="muted">%s</p></div><span class="pill status-warn">missing expected layer</span></div></section>',
    esc(section_id), esc(family), esc(context), esc(title), esc(title), esc(note)
  )
}

# Preserve the fail-closed missing-layer gate while reporting the actual report
# family and scope. A contrast Heatmaps placeholder is not a KEGG diagnosis.
report_expected_layer_diagnostics <- function(report_html) {
  marker <- 'data-lisa-layer-family="([^"]*)" data-lisa-layer-context="([^"]*)"'
  matches <- gregexpr(marker, report_html, perl = TRUE)
  values <- regmatches(report_html, matches)[[1L]]
  if (!length(values) || identical(values, "")) return(character())
  unique(vapply(values, function(value) {
    family <- sub(marker, "\\1", value, perl = TRUE)
    context <- sub(marker, "\\2", value, perl = TRUE)
    paste0("family=", family, "; context=", context)
  }, character(1L)))
}

empty_completed_section <- function(title, note, status_html = "") {
  section_id <- slug(title)
  sprintf(
    '<section id="%s" class="panel"><div class="section-head"><div><p class="eyebrow">%s</p><h2>%s</h2><p class="muted">%s</p></div><div>%s<span class="pill">0 figures generated</span></div></div></section>',
    esc(section_id), esc(title), esc(title), esc(note), status_html
  )
}

inapplicable_layer_section <- function(title, note, status_html = "") {
  section_id <- slug(title)
  sprintf(
    '<section id="%s" class="panel"><div class="section-head"><div><p class="eyebrow">%s</p><h2>%s</h2><p class="muted">%s</p></div><div>%s<span class="pill status-warn">not generated by this engine</span></div></div></section>',
    esc(section_id), esc(title), esc(title), esc(note), status_html
  )
}

not_requested_layer_section <- function(title, note, status_html = "") {
  section_id <- slug(title)
  sprintf(
    '<section id="%s" class="panel"><div class="section-head"><div><p class="eyebrow">%s</p><h2>%s</h2><p class="muted">%s</p></div><div>%s<span class="pill">not_requested</span></div></div></section>',
    esc(section_id), esc(title), esc(title), esc(note), status_html
  )
}

layer_section <- function(title, files, page_file, note, status_html = "", subtype_html = "",
                          empty_state = c("missing", "inapplicable", "completed_empty", "not_requested")) {
  empty_state <- match.arg(empty_state)
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    if (!length(report_image_formats())) return(not_requested_layer_section(title, "No image exports were requested.", status_html))
    if (identical(empty_state, "inapplicable")) return(inapplicable_layer_section(title, note, status_html))
    if (identical(empty_state, "completed_empty")) return(empty_completed_section(title, note, status_html))
    if (identical(empty_state, "not_requested")) return(not_requested_layer_section(title, note, status_html))
    return(missing_layer_section(title, page_file, note))
  }
  section_figures(title, files, page_file, max(length(files), 1), note, status_html, subtype_html)
}

subtype_badges <- function(files, specs) {
  if (length(files) == 0) return("")
  badges <- vapply(names(specs), function(label) {
    n <- sum(grepl(specs[[label]], files, ignore.case = TRUE))
    cls <- if (n > 0) "pill" else "pill status-warn"
    sprintf('<span class="%s">%s: %d</span>', cls, esc(label), n)
  }, character(1))
  paste0('<div class="links">', paste(badges, collapse = ""), '</div>')
}

collection_chips <- function(collection_dirs, page, base_label) {
  if (length(collection_dirs) == 0) return("")
  links <- paste(vapply(collection_dirs, function(collection_dir) {
    collection <- sub("^collection_", "", basename(collection_dir))
    href <- paste0(page, "#", layer_anchor(base_label, collection, "LISA categories"))
    sprintf('<a class="file-link" href="%s">%s</a>', esc(href), esc(collection))
  }, character(1)), collapse = "")
  paste0('<div class="links">', links, '</div>')
}

contrast_collection_chips <- function(collection_dirs, page, base_label) {
  if (length(collection_dirs) == 0) return("")
  links <- paste(vapply(collection_dirs, function(collection_dir) {
    collection <- sub("^collection_", "", basename(collection_dir))
    href <- paste0(page, "#", layer_anchor(base_label, collection, "LISA category shifts"))
    sprintf('<a class="file-link" href="%s">%s</a>', esc(href), esc(collection))
  }, character(1)), collapse = "")
  paste0('<div class="links">', links, '</div>')
}

top_gsea_rows <- function(project_dir, limit = 10) {
  files <- list.files(
    file.path(project_dir, "outputs", "lisa", "single_de"),
    pattern = "_GSEA_category_summary.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    files <- list.files(
      file.path(project_dir, "outputs", "single_de"),
      pattern = "_GSEA_category_summary.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
  }
  if (length(files) == 0) return(data.frame())
  rows <- rbind_fill_local(lapply(files, function(path) {
    df <- read_tsv(path)
    if (nrow(df) == 0) return(data.frame())
    parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
    collection <- sub("^collection_", "", parts[grep("^collection_", parts)[1]])
    single_idx <- which(parts == "single_de")
    analysis <- if (length(single_idx) > 0) parts[single_idx[1] + 1] else ""
    df$analysis_id <- analysis
    df$collection <- collection
    df
  }))
  if (is.null(rows) || nrow(rows) == 0) return(data.frame())
  rows$min_padj_num <- suppressWarnings(as.numeric(rows$min_padj))
  rows$mean_NES_abs <- abs(suppressWarnings(as.numeric(rows$mean_NES)))
  rows <- rows[order(rows$min_padj_num, -rows$mean_NES_abs, na.last = TRUE), , drop = FALSE]
  rows[seq_len(min(nrow(rows), limit)), , drop = FALSE]
}

short_title <- function(id) {
  labels <- get0("report_analysis_labels", inherits=TRUE, ifnotfound=character())
  if (id %in% names(labels) && nzchar(labels[[id]])) return(labels[[id]])
  id <- gsub("^P[0-9]+_", "", id)
  id <- gsub("GSE91061_", "", id)
  id <- gsub("on_treatment_vs_pre_all", "On-treatment vs pre", id)
  id <- gsub("pre_responder_vs_nonresponder", "Pre-treatment responder vs non-responder", id)
  id <- gsub("_vs_", " vs ", id)
  id <- gsub("_", " ", id)
  id
}

write_page <- function(file, html) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, file, useBytes = TRUE)
}

main <- function() {
  report_gallery_assets <<- list()
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  invisible(report_is_native_absolute_path(
    args$project_dir, reject_drive_relative = TRUE
  ))
  project_dir <- normalizePath(args$project_dir, mustWork = TRUE)
  package_dir <- report_executable_package_dir()
  report_package_dir <<- package_dir
  report_renderer_sha256 <<- args$lisa_internal_renderer_sha256
  report_project_dir <<- project_dir
  policy <- read_tsv(file.path(project_dir, "report_output_policy.tsv"))
  if (!nrow(policy)) policy <- data.frame(product = c("source_data", "recipes"), requested = TRUE)
  report_output_policy <<- policy
  if (identical(tolower(args$refresh_category_links), "true")) {
    report_refresh_category_links(project_dir)
    return(invisible(project_dir))
  }
  metadata <- read_tsv(file.path(project_dir, "lisa_pipeline_metadata.tsv"))
  report_mode <- if (nrow(metadata) && all(c("key", "value") %in% names(metadata)) &&
      any(metadata$key == "report_mode")) as.character(metadata$value[metadata$key == "report_mode"][[1]]) else "full"
  selected_path <- file.path(project_dir, "report_selection.tsv")
  if (file.exists(selected_path)) {
    selected <- read_tsv(selected_path)
    if (!all(c("key", "value") %in% names(selected)) || anyDuplicated(selected$key))
      stop("Invalid report selection manifest.", call. = FALSE)
    values <- stats::setNames(as.character(selected$value), selected$key)
    if (!all(c("schema_version", "report_mode", "analyses", "contrasts", "contrast_owners") %in% names(values)) ||
        values[["schema_version"]] != "1" || values[["report_mode"]] != "full")
      stop("Unsupported report selection manifest.", call. = FALSE)
    exact <- function(key) strsplit(values[[key]], ";", fixed = TRUE)[[1L]]
    saved_de <- read_tsv(file.path(project_dir, "config", "de_index.tsv"))
    saved_co <- read_tsv(file.path(project_dir, "config", "contrast_index.tsv"))
    de_dirs <- basename(list.dirs(file.path(project_dir, "outputs", "single_de"), recursive = FALSE))
    co_dirs <- basename(list.dirs(file.path(project_dir, "outputs", "category_contrasts"), recursive = FALSE))
    if (!identical(as.character(saved_de$analysis_id), exact("analyses")) ||
        !identical(as.character(saved_co$contrast_id), exact("contrasts")) ||
        !setequal(de_dirs, exact("analyses")) || !setequal(co_dirs, exact("contrast_owners")))
      stop("Report inputs do not match the declared selection.", call. = FALSE)
    report_mode <- values[["report_mode"]]
  }
  # --report-mode is the explicit request of the single run generator
  # (run_lisa(mode = ...)). Without it the historical resolution above applies,
  # so older run trees rebuild exactly as before.
  if (nzchar(args$report_mode)) {
    if (!args$report_mode %in% c("standard", "full"))
      stop("--report-mode must be standard or full.", call. = FALSE)
    if (file.exists(selected_path) && !identical(args$report_mode, report_mode))
      stop("--report-mode conflicts with the declared report selection manifest.", call. = FALSE)
    report_mode <- args$report_mode
  }
  report_route_label <<- report_mode
  standard_report <- identical(report_mode, "standard")
  kegg_maps_requested <- report_kegg_maps_requested(project_dir)
  if (args$complete_missing == "true") {
    if (standard_report) stop("Missing FULL products cannot be requested for a standard report.")
    lisaR:::lisa_complete_full_products(project_dir, args$kegg_maps == "true",
      args$kegg_cache_root, args$kegg_snapshot_id, args$kegg_access_mode)
  }
  copy_assets(package_dir, project_dir)
  write_assets(project_dir)
  pages_dir <- file.path(project_dir, "report_pages")
  dir.create(pages_dir, recursive = TRUE, showWarnings = FALSE)

  de_index <- read_tsv(file.path(project_dir, "config", "de_index.tsv"))
  contrast_index <- read_tsv(file.path(project_dir, "config", "contrast_index.tsv"))
  # Presentation is an example configuration, not edits to scientific results.
  presentation_path <- if (nzchar(args$presentation_config)) args$presentation_config else file.path(project_dir,"config","report_presentation.tsv")
  presentation <- if (file.exists(presentation_path)) read_tsv(presentation_path) else data.frame()
  if (nrow(presentation)) {
    if (!all(c("role","id","label","description") %in% names(presentation)) || anyDuplicated(paste(presentation$role,presentation$id))) stop("Invalid report presentation configuration.")
    for (role in c("analysis","contrast")) {
      ix <- if (role == "analysis") de_index else contrast_index
      key <- if (role == "analysis") "analysis_id" else "contrast_id"
      rows <- presentation[presentation$role == role,,drop=FALSE]
      if (nrow(rows) && !all(rows$id %in% ix[[key]])) stop("Presentation IDs outside selected scope.")
      for (i in seq_len(nrow(rows))) {
        match_row <- match(rows$id[[i]], ix[[key]])
        ix$label[match_row] <- rows$label[[i]]
        ix$comparison[match_row] <- rows$description[[i]]
      }
      if (role == "analysis") de_index <- ix else contrast_index <- ix
    }
    study <- presentation[presentation$role == "study",,drop=FALSE]
    if (nrow(study) == 1L) args$title <- study$label[[1L]]
  }
  report_analysis_labels <<- stats::setNames(if ("label" %in% names(de_index)) de_index$label else de_index$analysis_id, de_index$analysis_id)
  report_display_labels <<- report_analysis_labels
  if (nrow(contrast_index) && all(c("contrast_id", "output_id", "label") %in% names(contrast_index))) {
    report_display_labels <<- c(report_display_labels, stats::setNames(contrast_index$label,
      paste(contrast_index$contrast_id, contrast_index$output_id, sep = "_")))
  }
  fig_manifest <- read_tsv(file.path(project_dir, "manifests", "figure_format_manifest.tsv"))
  de_status <- biological_de_status(de_index)
  single_de_collection_status <- status_count_any(c(
    file.path(project_dir, "single_de_status.tsv"),
    file.path(project_dir, "reports", "lisa_single_de_status.tsv"),
    file.path(project_dir, "state", "de_status.tsv"),
    file.path(project_dir, "metrics", "de_runtime_size.tsv")
  ))
  contrast_status_paths <- c(
    file.path(project_dir, "contrast_status.tsv"),
    file.path(project_dir, "reports", "lisa_contrast_status.tsv"),
    file.path(project_dir, "state", "contrast_status.tsv"),
    file.path(project_dir, "metrics", "contrast_runtime_size.tsv")
  )
  contrast_status_df <- read_first_status_table(contrast_status_paths)
  contrast_status <- biological_contrast_status(contrast_index, contrast_status_df)
  contrast_collection_status <- status_count_any(contrast_status_paths)
  visual_status <- status_count_any(c(
    file.path(project_dir, "post_lisa_status.tsv"),
    file.path(project_dir, "metrics", "visual_runtime_size.tsv")
  ))
  runtime_df <- read_tsv(file.path(project_dir, "metrics", "visual_runtime_size.tsv"))
  validation_summary <- read_tsv(file.path(project_dir, "manifests", "report_product_validation_summary.tsv"))
  all_png <- list.files(file.path(project_dir, "outputs"), pattern = "[.]png$", recursive = TRUE, full.names = TRUE)
  all_svg <- list.files(file.path(project_dir, "outputs"), pattern = "\\.svg$", recursive = TRUE, full.names = TRUE)
  all_pdf <- list.files(file.path(project_dir, "outputs"), pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
  has_kegg_outputs <- has_kegg_report_outputs(project_dir)
  png_count <- metric_value(validation_summary, "png_files", length(all_png))
  svg_count <- metric_value(validation_summary, "svg_files", length(all_svg))
  pdf_count <- metric_value(validation_summary, "pdf_files_in_outputs", length(all_pdf))
  checked_links <- metric_value(validation_summary, "checked_local_links", NA_integer_)
  missing_links <- metric_value(validation_summary, "missing_local_links", NA_integer_)
  failed_layers <- if (nrow(runtime_df) > 0 && "status" %in% names(runtime_df)) runtime_df[runtime_df$status == "failed", , drop = FALSE] else data.frame()
  top_rows <- top_gsea_rows(project_dir, 8)
  de_cards <- if (nrow(de_index) > 0 && "analysis_id" %in% names(de_index)) {
    paste(vapply(seq_len(nrow(de_index)), function(i) {
      analysis_id <- as.character(de_index$analysis_id[[i]])
      label <- analysis_label_lookup(de_index, "analysis_id", analysis_id, short_title(analysis_id))
      chips <- collection_chips(single_collection_dirs(project_dir, analysis_id), "report_pages/single_de.html", short_title(analysis_id))
      model_note <- if ("model_note" %in% names(de_index)) {
        as.character(de_index$model_note[[i]])
      } else {
        ""
      }
      if (is.na(model_note) || !nzchar(trimws(model_note))) {
        model_note <- "Model design/adjustment not supplied."
      }
      sprintf(
        '<div class="detail-card"><a href="report_pages/single_de.html#%s"><strong>%s</strong></a><span>%s</span><span>%s</span>%s</div>',
        esc(nav_anchor_for_de(project_dir, analysis_id, "single_de")),
        esc(label), esc(model_note),
        if (standard_report) "Inspect the LISA categories and their member gene sets." else "Inspect the LISA categories and all applicable expanded evidence figures.",
        chips
      )
    }, character(1)), collapse = "")
  } else {
    '<div class="detail-card"><strong>No DE analyses indexed</strong><span>No rows were found in config/de_index.tsv.</span></div>'
  }
  contrast_cards <- if (nrow(contrast_index) > 0 && "output_id" %in% names(contrast_index)) {
    paste(vapply(seq_len(nrow(contrast_index)), function(i) {
      output_id <- as.character(contrast_index$output_id[[i]])
      contrast_id <- if ("contrast_id" %in% names(contrast_index)) {
        paste(as.character(contrast_index$contrast_id[[i]]), output_id, sep = "_")
      } else {
        output_id
      }
      label <- analysis_label_lookup(contrast_index, "output_id", output_id, short_title(output_id))
      contrast_title <- contrast_short_title(contrast_index, contrast_id)
      chips <- contrast_collection_chips(contrast_collection_dirs(project_dir, contrast_id), "report_pages/contrasts.html", contrast_title)
      sprintf(
        '<div class="detail-card"><a href="report_pages/contrasts.html#%s"><strong>%s</strong></a><span>%s</span>%s</div>',
        esc(nav_anchor_for_contrast(project_dir, contrast_id, "category", contrast_index)), esc(label),
        if (standard_report) "Inspect the LISA category shifts for this DE-vs-DE comparison." else "Inspect category shifts and all applicable expanded contrast figures.",
        chips
      )
    }, character(1)), collapse = "")
  } else {
    '<div class="detail-card"><strong>No DE-vs-DE contrasts</strong><span>This project has a single DE analysis or no planned DE-vs-DE contrast.</span></div>'
  }

  root_index <- file.path(project_dir, "report_index.html")
  warning_banner <- if (nrow(failed_layers) > 0) {
    notice_html(
      "warn",
      sprintf("%d post-LISA layers flagged failed", nrow(failed_layers)),
      "The runtime table flags post-LISA layers as failed, but recovered visual files are shown when present. QC should be read before biological interpretation."
    )
  } else {
    notice_html("info", "No failed visual layers", "All recorded post-LISA visual layers completed.")
  }
  validation_banner <- if (is.na(checked_links) || is.na(missing_links)) {
    notice_html("info", "Final link check follows report construction",
      "The verification gate records the final checked and missing-link counts after these pages have been written.")
  } else if (missing_links == 0) {
    notice_html("info", "Local links validated", sprintf("%s local links checked; no missing local target detected.", checked_links))
  } else {
    notice_html("warn", "Link validation requires review", sprintf("%s missing links reported in manifest.", missing_links))
  }
  analysis_map <- paste0(
    '<section id="analysis-map" class="panel"><div class="section-head"><div><p class="eyebrow">Navigation</p><h2>Analysis map</h2></div><a class="file-link" href="report_pages/qc.html#report-validation">Open validation/QC</a></div>',
    '<h3>DE analyses</h3><div class="detail-grid">', de_cards, '</div>',
    '<h3>DE-vs-DE contrasts</h3><div class="detail-grid">', contrast_cards, '</div></section>'
  )
  overview_body <- paste0(
    '<section class="hero"><h1>', esc(args$title), '</h1><p>',
    if (standard_report) "Standard LISA report with complete scientific tables, category summaries and member-gene-set figures. Optional expanded products have not been created." else "Full LISA report with the complete scientific tables and every applicable expanded product.",
    ' Start with the top category table, then open the collection pages for figures and source tables.</p></section>',
    notice_html(
      "info", "Portable report surface",
      paste0(
        "The complete run remains the canonical private provenance record. ",
        "Share only the files named in shareable_report_manifest.tsv; other ",
        "run folders and receipts are not part of this portable report."
      )
    ),
    report_kegg_maps_notice(project_dir),
    '<div class="grid">',
    sprintf('<div class="metric"><span>DE analyses</span><strong>%d/%d</strong></div>', de_status[["completed"]], de_status[["total"]]),
    sprintf('<div class="metric"><span>Single-DE collection runs</span><strong>%d/%d</strong></div>', single_de_collection_status[["completed"]], single_de_collection_status[["total"]]),
    sprintf('<div class="metric"><span>Biological contrasts</span><strong>%d/%d</strong></div>', contrast_status[["completed"]], contrast_status[["total"]]),
    sprintf('<div class="metric"><span>Contrast collection runs</span><strong>%d/%d</strong></div>', contrast_collection_status[["completed"]], contrast_collection_status[["total"]]),
    sprintf('<div class="metric"><span>Post-LISA layers</span><strong>%d/%d</strong></div>', visual_status[["completed"]], visual_status[["total"]]),
    sprintf('<div class="metric"><span>PNG figures</span><strong>%d</strong></div>', png_count),
    sprintf('<div class="metric"><span>SVG figures</span><strong>%d</strong></div>', svg_count),
    sprintf('<div class="metric"><span>PDF in outputs</span><strong>%d</strong></div>', pdf_count),
    '</div>',
    '<section id="top-category-signals" class="panel"><div class="section-head"><div><p class="eyebrow">Interpretation start point</p><h2>Top LISA category signals</h2></div><a class="file-link" href="report_pages/single_de.html">Open category figures</a></div>',
    '<p class="muted">These are the strongest category-level GSEA summaries across analyses and ontology collections, ranked by the minimum member-gene-set FDR and descriptive mean NES (not a category-level FDR). Use them as triage, not as a replacement for the collection-specific GSEA/ORA tables.</p>',
    table_html(top_rows, c("analysis_id", "collection", "category_display_name", "macrogroup_name", "mean_NES", "min_padj", "same_direction_pct", "mean_NES_direction"), "top_gsea_table", NULL, 8),
    '</section>',
    '<section class="panel"><div class="section-head"><div><p class="eyebrow">Run structure</p><h2>Analyses and contrasts</h2></div></div>',
    table_html(de_index, c("analysis_id", "label", "priority", "block", "dataset_id", "paper", "contrast_label", "species", "notes"), "de_index_table", file.path(project_dir, "config", "de_index.tsv"), 20),
    table_html(contrast_index, c("contrast_id", "output_id", "priority", "contrast_a", "contrast_b", "label_a", "label_b", "biological_question", "question"), "contrast_index_table", file.path(project_dir, "config", "contrast_index.tsv"), 20),
    '</section>'
  )
  if (!standard_report) {
    # Make the extra figures directly visible from the entry page, rather than
    # requiring a reader to discover a section dropdown on a standard-looking
    # page. The targets are the SAME integrated sections, not a second gallery.
    links <- character()
    for (analysis in as.character(de_index$analysis_id)) {
      if (!report_full_scope_owner(project_dir, analysis)) next
      for (dir in single_collection_dirs(project_dir, analysis)) {
        collection <- sub("^collection_", "", basename(dir))
        label <- paste(short_title(analysis), collection, "FULL figures", sep = " - ")
        links <- c(links, sprintf('<a class="file-link" href="report_pages/single_de.html#%s">%s</a>',
          esc(slug(label)), esc(label)))
      }
    }
    for (i in seq_len(nrow(contrast_index))) {
      owner <- paste(contrast_index$contrast_id[[i]], contrast_index$output_id[[i]], sep = "_")
      if (!report_full_scope_owner(project_dir, owner, TRUE)) next
      for (dir in contrast_collection_dirs(project_dir, owner)) {
        collection <- sub("^collection_", "", basename(dir))
        label <- paste(contrast_short_title(contrast_index, owner), collection, "FULL figures", sep = " - ")
        links <- c(links, sprintf('<a class="file-link" href="report_pages/contrasts.html#%s">%s</a>',
          esc(slug(label)), esc(label)))
      }
    }
    if (length(links)) overview_body <- paste0(
      '<section id="full-figures-entry" class="panel"><h2>FULL figures — open the extra figures</h2>',
      '<p>Heatmaps, volcano plots and the other saved category products are incorporated into each analysis and contrast below.</p>',
      '<div class="links">', paste(links, collapse = ""), '</div></section>', overview_body)
  }
  evidence_indexes <- list.files(file.path(project_dir, "report_pages", "evidence"),
    pattern = "^index[.]html$", recursive = TRUE, full.names = TRUE)
  if (length(evidence_indexes)) {
    evidence_links <- vapply(sort(evidence_indexes), function(path) {
      label <- paste(basename(dirname(dirname(path))), basename(dirname(path)), sep = " · ")
      paste0('<li><a href="', esc(rel_path(path, root_index)), '">', esc(label), '</a></li>')
    }, character(1))
    overview_body <- paste0(overview_body,
      '<section class="panel"><h2>Category evidence: gene sets and genes</h2>',
      '<p>Inspect exact support and leading-edge overlap without a new prioritization score.</p><ul>',
      paste(evidence_links, collapse = ""), '</ul></section>')
  }
  gene_evidence_index <- file.path(project_dir, "report_pages", "gene_evidence", "index.html")
  if (file.exists(gene_evidence_index)) overview_body <- paste0(overview_body,
    '<section class="panel"><h2>Evidence for your gene</h2>',
    '<p>Search a symbol or identifier to inspect DE measurements and every significant member gene set. Full set membership and leading-edge membership are reported separately, for individual analyses.</p>',
    '<a class="file-link" href="', esc(rel_path(gene_evidence_index, root_index)), '">Open gene evidence search</a></section>')
  overview_body <- report_overview_entry(args$title, report_mode, de_index, contrast_index,
    project_dir, analysis_map, overview_body)
  write_page(root_index, page_shell(args$title, "overview", overview_body, project_dir, root_index, root = TRUE))

  single_file <- file.path(pages_dir, "single_de.html")
  single_blocks <- character()
  single_toc <- list()
  single_navigation <- list()
  single_body <- paste0(
    '<section class="hero"><h1>Analyses</h1><p>',
    if (standard_report) "Each analysis contains its LISA category summary and the gene sets that support each applicable category." else "Each analysis contains its LISA category results and all applicable expanded evidence figures.",
    '</p></section>'
  )
  single_ids <- if ("analysis_id" %in% names(de_index)) as.character(de_index$analysis_id) else basename(list.dirs(file.path(project_dir, "outputs", "single_de"), recursive = FALSE, full.names = TRUE))
  for (analysis in single_ids) {
    analysis_label <- analysis_label_lookup(de_index, "analysis_id", analysis, short_title(analysis))
    # A FULL run may be scoped to selected owners only (STANDARD shell plus
    # extra figures for a subset of analyses/contrasts): an owner with no
    # full-only directory anywhere (neither the in-run
    # `outputs/gene_level/...` product tree nor a saved `artifacts/...`
    # subtree) was never selected for extras and must degrade exactly like a
    # STANDARD owner, not fail the required-layer gate below.
    analysis_full_scope <- !standard_report &&
      report_full_scope_owner(project_dir, analysis, contrast = FALSE)
    collection_blocks <- character()
    navigation_collections <- list()
    for (collection_dir in single_collection_dirs(project_dir, analysis)) {
      collection <- sub("^collection_", "", basename(collection_dir))
      gene_collection_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis, paste0("collection_", collection))
      base_title <- paste(short_title(analysis), collection, sep = " - ")
      lisa_imgs <- pngs(collection_dir)
      lisa_imgs <- lisa_imgs[grepl("GSEA|ORA|lollipop|barplot", lisa_imgs) &
        !grepl("GSEA_category_pathways|/category_nes/", lisa_imgs)]
      layers <- list(
        "LISA categories" = list(path = collection_dir, files = lisa_imgs, note = "LISA category summaries: mean-NES and direction/support lollipops, plus the gene-set overview. ORA figures appear only when explicitly requested."),
        "GeneCards" = list(path = file.path(gene_collection_dir, "category_gene_cards"), files = pngs(file.path(gene_collection_dir, "category_gene_cards"), "\\.png$"), note = "Category GeneCards: top supporting genes, log2FC, DE FDR, GSEA/ORA evidence and recurrence per LISA category.", empty_state = if (!analysis_full_scope) "not_requested" else "missing"),
        "Volcano overlays" = list(path = file.path(gene_collection_dir, "category_volcano_overlays"), files = pngs(file.path(gene_collection_dir, "category_volcano_overlays"), "\\.png$"), note = "Category-specific volcano overlays: DE landscape with LISA-supporting genes highlighted and labelled.", empty_state = if (!analysis_full_scope) "not_requested" else "missing"),
        "Recurrent genes" = list(path = file.path(gene_collection_dir, "recurrent_gene_screen"), files = pngs(file.path(gene_collection_dir, "recurrent_gene_screen"), "\\.png$"), note = if (identical(collection, "HALLMARKS")) "Not applicable to HALLMARKS: recurrence across semantic LISA categories is undefined for a direct Hallmark gene-set collection." else "Recurrent gene analysis: genes repeatedly supporting semantic categories within this analysis and collection.", empty_state = if (!analysis_full_scope) "not_requested" else if (identical(collection, "HALLMARKS")) "inapplicable" else "missing"),
        "Heatmaps" = list(path = file.path(gene_collection_dir, "leading_edge_gene_heatmaps"), files = pngs(file.path(gene_collection_dir, "leading_edge_gene_heatmaps"), "\\.png$"), note = "Leading-edge/supporting-gene heatmaps from the registered auxiliary expression matrix.", empty_state = if (!analysis_full_scope) "not_requested" else "missing"),
        "KEGG maps / painted pathways" = list(path = file.path(gene_collection_dir, "kegg_painter"), files = pngs(file.path(gene_collection_dir, "kegg_painter"), "_painted\\.png$"), note = "Canonical KEGG pathway diagrams painted from this collection's post-LISA evidence genes.", empty_state = if (!analysis_full_scope) "not_requested" else "missing")
      )
      category_product_specs <- list(
        list(product = "member_gene_sets", label = "Member gene sets", path = file.path(collection_dir, "plots", "GSEA_category_pathways")),
        list(product = "gene_cards", label = "Prioritized genes", path = file.path(gene_collection_dir, "category_gene_cards")),
        list(product = "volcano", label = "Volcano overlays", path = file.path(gene_collection_dir, "category_volcano_overlays")),
        list(product = "heatmap", label = "Gene heatmaps", path = file.path(gene_collection_dir, "leading_edge_gene_heatmaps")),
        # Immutable full-extension artifacts use the same exact sidecar
        # contract.  Absent directories are simply empty, so this does not
        # manufacture products for ordinary pipeline reports.
        list(product = "member_gene_sets", label = "Member gene sets", path = file.path(project_dir, "artifacts", analysis, collection, "member_gene_sets")),
        list(product = "gene_cards", label = "Prioritized genes", path = file.path(project_dir, "artifacts", analysis, collection, "gene_cards")),
        list(product = "volcano", label = "Volcano overlays", path = file.path(project_dir, "artifacts", analysis, collection, "volcano")),
        list(product = "heatmap", label = "Gene heatmaps", path = file.path(project_dir, "artifacts", analysis, collection, "heatmap")),
        list(product = "kegg", label = "KEGG gene sets", path = file.path(project_dir, "artifacts", analysis, collection, "kegg"))
      )
      # These are one-category products.  They are rendered once in the
      # category cards below rather than duplicated as independent galleries.
      # An owner outside the selected FULL scope degrades exactly like a
      # STANDARD owner (only "LISA categories"), never an unattached
      # "missing" gallery.
      single_evidence_available <- file.exists(file.path(project_dir, "report_pages",
        "evidence", analysis, collection, "index.html"))
      if (analysis_full_scope && !single_evidence_available) {
        # A collection without a category-evidence page (HALLMARKS) cannot attach
        # its FULL products to category cards; they are shown as complete
        # galleries so that no FULL figure is left without a place in the report.
        full_art <- file.path(project_dir, "artifacts", analysis, collection)
        full_galleries <- c("Member gene sets" = "member_gene_sets", "GeneCards" = "gene_cards",
          "Volcano overlays" = "volcano", "Heatmaps" = "heatmap", "KEGG gene sets" = "kegg")
        for (gallery_name in names(full_galleries)) {
          gallery_dir <- file.path(full_art, full_galleries[[gallery_name]])
          gallery_files <- pngs(gallery_dir, "\\.png$")
          if (length(gallery_files)) layers[[gallery_name]] <- list(path = gallery_dir,
            files = gallery_files, note = paste0(gallery_name,
              ": every category of this collection. This collection has no category-evidence page, so its FULL figures are listed here."),
            empty_state = "missing")
        }
      } else if (analysis_full_scope) layers <- layers[!names(layers) %in% c("GeneCards", "Volcano overlays", "Heatmaps")]
      if (!analysis_full_scope) layers <- layers[intersect(c("LISA categories",
        if (kegg_maps_requested) "KEGG maps / painted pathways"), names(layers))]
      # Category evidence follows category figures as one continuous journey.
      evidence_index <- file.path(project_dir, "report_pages", "evidence", analysis, collection, "index.html")
      navigation_dir <- file.path(project_dir, "report_pages", "category_navigation", analysis, collection)
      nav_meta_path <- file.path(navigation_dir,"metadata.json")
      if (file.exists(nav_meta_path)) {
        nm <- jsonlite::read_json(nav_meta_path,simplifyVector=TRUE)
        tab <- read_tsv(file.path(navigation_dir,"tables","category_navigation_source.tsv"))
        nm$display_label <- short_title(analysis)
        one_nav <- structure(list(metadata=nm,categories=tab),class="lisa_category_navigation")
        support <- lisaR:::lisa_hommel_read_results(collection_dir, analysis, collection)
        if (!is.null(support)) one_nav <- lisaR:::lisa_hommel_attach(one_nav, support)
        lisaR:::render_lisa_category_navigation(one_nav,navigation_dir,
          download_policy = list(svg = report_requested("svg", TRUE),
            source_data = report_requested("source_data", TRUE),
            recipes = report_requested("recipes", TRUE)))
      }
      product_attachment <- if (analysis_full_scope) {
        # An explicitly selected saved-result report must account for its
        # entire supplied artifact inventory. Keep the saved FULL member-set
        # rendering (with its own PDF/source/recipe), not a second standard
        # output rendering merged into the same asset family.
        # run_lisa(mode = "full") renders its own saved FULL products into this
        # run's artifacts/, exactly like a selected saved-result report; they are
        # attached (member gene sets included), never left unreferenced.
        selected_view <- file.exists(selected_path) || identical(args$report_mode, "full")
        if (selected_view && dir.exists(file.path(project_dir,"artifacts",analysis,collection,"member_gene_sets")))
          category_product_specs <- Filter(function(spec)
            spec$product != "member_gene_sets" || startsWith(spec$path, file.path(project_dir, "artifacts")),
            category_product_specs)
        report_attach_category_products(evidence_index, category_product_specs, analysis, collection,
          exclude_products = if (selected_view) character() else "member_gene_sets", project_dir = project_dir)
      } else list(category_products = list(), unclassified = character())
      block_sections <- character()
      navigation_sections <- list()
      if (file.exists(evidence_index)) {
        support <- lisaR:::lisa_hommel_read_results(collection_dir, analysis, collection)
        if (!is.null(support)) lisaR:::lisa_refresh_hommel_evidence_page(evidence_index, support)
        report_refresh_evidence_static_assets(dirname(evidence_index), "category-evidence")
        report_replace_evidence_payload(evidence_index, NULL, "evidence-data")
        contract_file <- file.path(project_dir, "contract_manifest.tsv")
        legacy_enabled <- FALSE
        if (file.exists(contract_file)) {
          ec <- read_tsv(contract_file)
          if (all(c("key", "value") %in% names(ec))) {
            legacy_enabled <- any(ec$key == "report_legacy_gene_products" & tolower(ec$value) == "true")
          }
        }
        # FULL retains prioritized and recurrent genes; category evidence is not a substitute.
      }
      if (analysis_full_scope && !file.exists(evidence_index)) {
        full_title <- paste(base_title, "FULL figures", sep = " - ")
        full_section <- report_full_gallery_section(full_title,
          Filter(function(z) startsWith(z$path, file.path(project_dir, "artifacts")), layers),
          single_file, project_dir)
        if (nzchar(full_section)) {
          single_toc[[length(single_toc) + 1L]] <- c(id = slug(full_title), label = full_title)
          navigation_sections[[length(navigation_sections) + 1L]] <- report_navigation_section(full_title, "FULL figures")
          block_sections <- c(block_sections, full_section)
        }
      }
      for (layer_name in names(layers)) {
        spec <- layers[[layer_name]]
        title <- paste(base_title, layer_name, sep = " - ")
        if (identical(layer_name, "LISA categories") && file.exists(evidence_index)) {
          evidence_title <- paste(base_title, "Category evidence", sep = " - ")
          single_toc[[length(single_toc) + 1]] <- c(id = slug(evidence_title), label = evidence_title)
          navigation_sections[[length(navigation_sections) + 1L]] <- report_navigation_section(evidence_title, "Category evidence")
          block_sections <- c(block_sections, report_category_evidence_section(evidence_title,
            navigation_dir, evidence_index, single_file, product_attachment$unclassified,
            collection_dir = collection_dir))
          # The attached extras were previously subtracted from this page and
          # never surfaced. Announce them here -- after the evidence cover, so
          # the standard evidence still comes first, and before the LISA
          # categories layer, so a reader meets them without scrolling past the
          # standard tables. Unreachable in STANDARD mode: analysis_full_scope
          # is hard-wired FALSE whenever standard_report is TRUE.
          full_title <- paste(base_title, "FULL figures", sep = " - ")
          full_section <- if (analysis_full_scope) {
            report_full_products_section(full_title, product_attachment, evidence_index,
              single_file, "category", analysis, collection, project_dir = project_dir,
              related_sections = paste(base_title,names(layers)[vapply(layers,function(z) length(z$files)>0L,logical(1)) &
                names(layers) %in% c("Recurrent genes","KEGG maps / painted pathways")],sep=" - "),
              native_layers=layers[intersect(c("KEGG maps / painted pathways","Recurrent genes"),names(layers))])
          } else ""
          if (nzchar(full_section)) {
            single_toc[[length(single_toc) + 1]] <- c(id = slug(full_title), label = full_title)
            navigation_sections[[length(navigation_sections) + 1L]] <-
              report_navigation_section(full_title, "FULL figures")
            block_sections <- c(block_sections, full_section)
          }
        }
        single_toc[[length(single_toc) + 1]] <- c(id = slug(title), label = title)
        navigation_sections[[length(navigation_sections) + 1L]] <- report_navigation_section(title, layer_name)
        if (identical(layer_name, "LISA categories")) {
          block_sections <- c(block_sections, report_lisa_categories_section(title, collection_dir,
            spec$files, single_file, runtime_notice(runtime_df, spec$path), spec$note, evidence_index))
        } else if (identical(layer_name, "KEGG maps / painted pathways") && standard_report && !kegg_maps_requested) {
          block_sections <- c(block_sections, layer_section(title, character(), single_file,
            spec$note, runtime_notice(runtime_df, spec$path), "", empty_state = "not_requested"))
        } else if (identical(layer_name, "KEGG maps / painted pathways")) {
          state <- kegg_report_state(project_dir, "single", analysis, collection, spec$path, spec$files)
          block_sections <- c(block_sections, kegg_layer_section(title, spec$files, single_file, spec$note, state,
            runtime_notice(runtime_df, spec$path), subtype_badges(spec$files, c("PNG figures" = "\\.png$"))))
        } else block_sections <- c(block_sections, layer_section(
          title, spec$files, single_file, spec$note, runtime_notice(runtime_df, spec$path),
          subtype_badges(spec$files, c("PNG figures" = "\\.png$")),
          empty_state = if (!is.null(spec$empty_state)) spec$empty_state else "missing"))
      }
      collection_blocks <- c(collection_blocks, sprintf(
        '<details class="collection-block"%s data-lisa-nav-collection="%s" data-lisa-context-json="%s"><summary><strong>%s</strong><span>%s sections</span></summary><div class="collection-inner">%s</div></details>',
        if (!length(collection_blocks)) " open" else "", esc(collection),
        report_scope_context(analysis, analysis_label, collection,
          file.path(dirname(evidence_index), "tables", "metadata.tsv")),
        esc(collection), length(block_sections), paste(block_sections, collapse = "\n")
      ))
      navigation_collections[[length(navigation_collections) + 1L]] <- list(
        id = collection, label = collection, sections = navigation_sections)
    }
    single_blocks <- c(single_blocks, sprintf(
      '<details id="%s" class="analysis-block"%s data-lisa-nav-context="%s"><summary><strong>%s</strong><span>%d collections</span></summary><div class="analysis-inner">%s</div></details>',
      esc(slug(short_title(analysis))), if (!length(single_blocks)) " open" else "", esc(analysis),
      esc(analysis_label), length(collection_blocks), paste(collection_blocks, collapse = "\n")
    ))
    single_navigation[[length(single_navigation) + 1L]] <- list(
      id = analysis, scientific_id = analysis, label = analysis_label, kind = "analysis", collections = navigation_collections)
  }
  single_body <- paste0(single_body, toc_html(single_toc), paste(single_blocks, collapse = "\n"))
  write_page(single_file, page_shell(args$title, "single_de", single_body, project_dir, single_file,
    navigation = single_navigation))

  contrasts_file <- file.path(pages_dir, "contrasts.html")
  contrast_body <- paste0('<section class="hero"><h1>Contrasts</h1><p>',
    if (standard_report) "Standard reports show the LISA category shifts between each compatible pair of analyses." else "Full reports show category shifts and every applicable expanded contrast product.",
    '</p></section>')
  contrast_blocks <- character()
  contrast_toc <- list()
  contrast_navigation <- list()
  contrast_dirs <- list.dirs(file.path(project_dir, "outputs", "lisa", "contrast"), recursive = FALSE, full.names = TRUE)
  if (length(contrast_dirs) == 0) contrast_dirs <- list.dirs(file.path(project_dir, "outputs", "category_contrasts"), recursive = FALSE, full.names = TRUE)
  cg_root <- file.path(project_dir, "outputs", "gene_level", "category_contrasts")
  contrast_ids <- basename(contrast_dirs)
  for (contrast_name in contrast_ids) {
    output_id <- contrast_output_id(contrast_index, contrast_name)
    contrast_title <- contrast_short_title(contrast_index, contrast_name)
    contrast_label <- analysis_label_lookup(contrast_index, "output_id", output_id, contrast_title)
    contrast_evidence_owner <- contrast_name
    if (all(c("contrast_id", "output_id") %in% names(contrast_index))) {
      matched <- which(paste(contrast_index$contrast_id, contrast_index$output_id, sep = "_") == contrast_name)
      if (length(matched) == 1L) contrast_evidence_owner <- as.character(contrast_index$contrast_id[[matched]])
    }
    # See the matching analysis-side comment: a contrast outside the selected
    # FULL scope has neither an `outputs/gene_level/...` product tree nor a
    # saved `artifacts/contrasts/...` subtree, and must degrade like STANDARD.
    contrast_full_scope <- !standard_report &&
      report_full_scope_owner(project_dir, contrast_name, contrast = TRUE)
    collection_blocks <- character()
    navigation_collections <- list()
    for (collection_dir in contrast_collection_dirs(project_dir, contrast_name)) {
      collection <- sub("^collection_", "", basename(collection_dir))
      gene_collection_dir <- file.path(cg_root, contrast_name, paste0("collection_", collection))
      base_title <- paste(contrast_title, collection, sep = " - ")
      layers <- list(
        "LISA category shifts" = list(path = collection_dir, files = pngs(collection_dir)[!grepl("/category_nes/", pngs(collection_dir))], note = "Dumbbell plots summarize category-level movement between the two DE analyses. Prioritize annotated all-category and opposite-direction views."),
        "Contrast GeneCards" = list(path = file.path(gene_collection_dir, "contrast_category_cards"), files = pngs(file.path(gene_collection_dir, "contrast_category_cards"), "\\.png$"), note = "Contrast category cards: interpretable A/B category-level gene drivers.", empty_state = if (!contrast_full_scope) "not_requested" else "missing"),
        "Gene-category network" = list(path = file.path(gene_collection_dir, "contrast_gene_category_network"), files = pngs(file.path(gene_collection_dir, "contrast_gene_category_network"), "[.]png$"), note = "Supporting genes across categories, with paired A/B effects.", empty_state = if (!contrast_full_scope) "not_requested" else "missing"),
        "Paired gene heatmaps" = list(path = gene_collection_dir, files = pngs(gene_collection_dir, "paired_gene_heatmap\\.png$"), note = "Paired gene heatmap: direct A/B log2FC comparison for supporting genes.", empty_state = if (!contrast_full_scope) "not_requested" else "missing"),
        "KEGG maps / painted pathways" = list(path = file.path(gene_collection_dir, "contrast_kegg_pathway_painter"), files = pngs(file.path(gene_collection_dir, "contrast_kegg_pathway_painter"), "contrast_painted\\.png$"), note = "KEGG painted pathway diagrams with source tables when generated.")
      )
      category_product_specs <- list(
        list(product = "contrast_gene_cards", label = "Prioritized genes (paired)", path = file.path(gene_collection_dir, "contrast_category_cards")),
        list(product = "contrast_profile", label = "Category profiles", path = file.path(project_dir, "artifacts", "contrasts", contrast_name, collection, "contrast_profile")),
        list(product = "contrast_heatmap", label = "Category heatmaps", path = file.path(project_dir, "artifacts", "contrasts", contrast_name, collection, "contrast_heatmap"))
      )
      evidence_index <- file.path(project_dir, "report_pages", "contrast_evidence", contrast_name, collection, "index.html")
      navigation_dir <- file.path(project_dir, "report_pages", "category_navigation", contrast_name, collection)
      nav_meta_path <- file.path(navigation_dir, "metadata.json")
      if (file.exists(nav_meta_path)) {
        nm <- jsonlite::read_json(nav_meta_path, simplifyVector = TRUE)
        tab <- read_tsv(file.path(navigation_dir, "tables", "contrast_navigation_source.tsv"))
        one_nav <- structure(list(metadata = nm, categories = tab), class = "lisa_contrast_navigation")
        lisaR:::render_lisa_contrast_navigation(one_nav, navigation_dir,
          download_policy = list(svg = report_requested("svg", TRUE),
            source_data = report_requested("source_data", TRUE),
            recipes = report_requested("recipes", TRUE)))
      }
      if (file.exists(evidence_index)) {
        report_refresh_evidence_static_assets(dirname(evidence_index), "contrast-evidence")
        report_replace_evidence_payload(evidence_index, NULL, "contrast-evidence-data")
      }
      product_attachment <- if (contrast_full_scope) {
        report_attach_category_products(evidence_index, category_product_specs,
          contrast_evidence_owner, collection, contrast = TRUE,
          project_dir = project_dir)
      } else list(category_products = list(), unclassified = character())
      if (contrast_full_scope) {
        # Contrast heatmaps in a FULL saved-input extension are exact-category
        # cards, just like single-analysis heatmaps. Do not also require or
        # relabel them as a nonexistent multicategory paired-gene gallery.
        represented_in_cards <- file.exists(evidence_index) &&
          report_has_complete_category_product(product_attachment, "contrast_heatmap")
        attached_layers <- "Contrast GeneCards"
        if (file.exists(evidence_index)) {
          layers <- layers[!names(layers) %in% attached_layers]
        } else {
          # No contrast-evidence page for this collection (HALLMARKS): show its
          # FULL contrast products as galleries instead of attaching them.
          full_art <- file.path(project_dir, "artifacts", "contrasts", contrast_name, collection)
          for (gallery in list(c("Category profiles", "contrast_profile"), c("Category heatmaps", "contrast_heatmap"))) {
            gallery_dir <- file.path(full_art, gallery[[2L]])
            gallery_files <- pngs(gallery_dir, "\\.png$")
            if (length(gallery_files)) layers[[gallery[[1L]]]] <- list(path = gallery_dir,
              files = gallery_files, note = paste0(gallery[[1L]],
                ": every category of this collection. This collection has no contrast-evidence page, so its FULL figures are listed here."),
              empty_state = "missing")
          }
        }
      }
      if (!contrast_full_scope) layers <- layers[intersect(c("LISA category shifts",
        if (kegg_maps_requested) "KEGG maps / painted pathways"), names(layers))]
      network_files <- pngs(file.path(gene_collection_dir, "contrast_gene_category_network"), "\\.png$")
      if (length(network_files) > 0) {
        layers[["Gene-category network"]] <- list(path = file.path(gene_collection_dir, "contrast_gene_category_network"), files = network_files, note = "Gene-category networks: category-gene connectivity for contrast interpretation.")
      }
      block_sections <- character()
      navigation_sections <- list()
      if (contrast_full_scope && !file.exists(evidence_index)) {
        full_title <- paste(base_title, "FULL figures", sep = " - ")
        full_section <- report_full_gallery_section(full_title,
          Filter(function(z) startsWith(z$path, file.path(project_dir, "artifacts")), layers),
          contrasts_file, project_dir)
        if (nzchar(full_section)) {
          contrast_toc[[length(contrast_toc) + 1L]] <- c(id = slug(full_title), label = full_title)
          navigation_sections[[length(navigation_sections) + 1L]] <- report_navigation_section(full_title, "FULL figures")
          block_sections <- c(block_sections, full_section)
        }
      }
      for (layer_name in names(layers)) {
        spec <- layers[[layer_name]]
        title <- paste(base_title, layer_name, sep = " - ")
        if (identical(layer_name, "LISA category shifts") && file.exists(evidence_index)) {
          evidence_title <- paste(base_title, "Contrast evidence", sep = " - ")
          contrast_toc[[length(contrast_toc) + 1]] <- c(id = slug(evidence_title), label = evidence_title)
          navigation_sections[[length(navigation_sections) + 1L]] <- report_navigation_section(evidence_title, "Contrast evidence")
          block_sections <- c(block_sections, paste0('<section class="panel contrast-evidence-cover" id="',
            esc(slug(evidence_title)), '"><h2>Contrast evidence</h2>',
            report_existing_category_evidence(collection_dir, evidence_index, contrasts_file, base_title),
            if (length(product_attachment$unclassified)) paste0('<details class="report-unclassified-products"><summary>Unclassified saved product files (',
              length(product_attachment$unclassified), ')</summary><p class="muted">These files could not be attached to one selected category from explicit source metadata. They remain available here and are not assigned by filename.</p><div class="links">',
              paste(vapply(sort(unique(product_attachment$unclassified)), function(path) report_unclassified_product_link(path, contrasts_file), character(1L)), collapse = ""),
              '</div></details>') else "", '</section>'))
          # Native contrast products and saved category products share this report.
          full_title <- paste(base_title, "FULL figures", sep = " - ")
          full_section <- if (contrast_full_scope) {
            report_full_products_section(full_title, product_attachment, evidence_index,
              contrasts_file, "contrast", contrast_evidence_owner, collection,
              project_dir = project_dir, related_sections = paste(base_title,
                names(layers)[vapply(layers,function(z) length(z$files)>0L,logical(1)) &
                  names(layers) %in% c("Recurrent genes","Gene-category network","Paired gene heatmaps","KEGG maps / painted pathways")], sep=" - "),
              native_layers=layers[intersect(c("KEGG maps / painted pathways","Paired gene heatmaps","Gene-category network"),names(layers))])
          } else ""
          if (nzchar(full_section)) {
            contrast_toc[[length(contrast_toc) + 1]] <- c(id = slug(full_title), label = full_title)
            navigation_sections[[length(navigation_sections) + 1L]] <-
              report_navigation_section(full_title, "FULL figures")
            block_sections <- c(block_sections, full_section)
          }
        }
        contrast_toc[[length(contrast_toc) + 1]] <- c(id = slug(title), label = title)
        navigation_sections[[length(navigation_sections) + 1L]] <- report_navigation_section(title,
          if (identical(layer_name, "LISA category shifts")) "LISA categories" else layer_name)
        if (identical(layer_name, "LISA category shifts")) {
          block_sections <- c(block_sections, report_lisa_categories_section(title, collection_dir,
            spec$files, contrasts_file, runtime_notice(runtime_df, spec$path), spec$note, evidence_index))
        } else if (identical(layer_name, "KEGG maps / painted pathways") && standard_report && !kegg_maps_requested) {
          block_sections <- c(block_sections, layer_section(title, character(), contrasts_file,
            spec$note, runtime_notice(runtime_df, spec$path), "", empty_state = "not_requested"))
        } else if (identical(layer_name, "KEGG maps / painted pathways")) {
          state <- kegg_report_state(project_dir, "contrast", contrast_name, collection, spec$path, spec$files)
          block_sections <- c(block_sections, kegg_layer_section(title, spec$files, contrasts_file, spec$note, state,
            runtime_notice(runtime_df, spec$path), subtype_badges(spec$files, c("PNG figures" = "\\.png$"))))
        } else block_sections <- c(block_sections, layer_section(
          title, spec$files, contrasts_file, spec$note, runtime_notice(runtime_df, spec$path),
          subtype_badges(spec$files, c("PNG figures" = "\\.png$")),
          empty_state = if (!is.null(spec$empty_state)) spec$empty_state else "missing"))
      }
      collection_blocks <- c(collection_blocks, sprintf(
        '<details class="collection-block"%s data-lisa-nav-collection="%s" data-lisa-context-json="%s"><summary><strong>%s</strong><span>%s sections</span></summary><div class="collection-inner">%s</div></details>',
        if (!length(collection_blocks)) " open" else "", esc(collection),
        report_scope_context(contrast_name, contrast_label, collection,
          file.path(project_dir, "report_pages", "contrast_evidence", contrast_name, collection, "tables", "metadata.tsv"), contrast = TRUE),
        esc(collection), length(block_sections), paste(block_sections, collapse = "\n")
      ))
      navigation_collections[[length(navigation_collections) + 1L]] <- list(
        id = collection, label = collection, sections = navigation_sections)
    }
    contrast_blocks <- c(contrast_blocks, sprintf(
      '<details id="%s" class="analysis-block"%s data-lisa-nav-context="%s"><summary><strong>%s</strong><span>%d collections</span></summary><div class="analysis-inner">%s</div></details>',
      esc(slug(contrast_title)), if (!length(contrast_blocks)) " open" else "", esc(contrast_name),
      esc(contrast_label), length(collection_blocks), paste(collection_blocks, collapse = "\n")
    ))
    scientific_id <- contrast_name
    if (all(c("contrast_id", "output_id") %in% names(contrast_index))) {
      matching <- which(paste(contrast_index$contrast_id, contrast_index$output_id, sep = "_") == contrast_name)
      if (length(matching) == 1L) scientific_id <- as.character(contrast_index$contrast_id[[matching]])
    }
    contrast_navigation[[length(contrast_navigation) + 1L]] <- list(
      id = contrast_name, scientific_id = scientific_id, label = contrast_label, kind = "contrast", collections = navigation_collections)
  }
  contrast_body <- paste0(contrast_body, toc_html(contrast_toc), paste(contrast_blocks, collapse = "\n"))
  write_page(contrasts_file, page_shell(args$title, "contrasts", contrast_body, project_dir, contrasts_file,
    navigation = contrast_navigation))
  report_write_navigation_inventory(project_dir, single_navigation, contrast_navigation)

  for (stale_page in file.path(pages_dir, c("single_gene.html", "contrast_gene.html", "kegg.html"))) {
    if (file.exists(stale_page)) unlink(stale_page)
  }

  downloads_file <- file.path(pages_dir, "downloads.html")
  tsvs <- list.files(
    project_dir, pattern = "\\.(tsv|csv)$", recursive = TRUE,
    full.names = TRUE
  )
  tsvs <- tsvs[
    !is_report_excluded_path(tsvs) &
      basename(tsvs) != report_shareable_manifest_name
  ]
  if (!report_requested("source_data", TRUE)) {
    # Keep analytical tables, inputs and QC downloadable. Figure sidecars are
    # internal evidence, not exports when the caller disabled source data.
    relative <- vapply(tsvs, report_relative_path, character(1), from_dir = project_dir)
    tsvs <- tsvs[!grepl("^artifacts/|^report_figure_data/|_(source|matrix)[.]tsv$", relative)]
  }
  downloads_body <- paste0(
    '<section class="hero"><h1>Downloads</h1><p>Portable TSV/CSV copies generated from the canonical private run.</p></section>',
    notice_html(
      "info", "Sharing boundary",
      paste0(
        "Copy only report_index.html and the components listed in ",
        "shareable_report_manifest.tsv. The rest of the run can contain ",
        "private paths, source receipts or unpublished inputs and is not ",
        "directly shareable. The manifest lists itself with SELF because a ",
        "file cannot contain its own stable SHA-256 digest."
      )
    ),
    '<p><a class="file-link" download href="../',
    report_shareable_manifest_name,
    '">Download shareable report manifest</a></p>',
    download_table_html(tsvs, downloads_file, "download_table")
  )
  write_page(downloads_file, page_shell(args$title, "downloads", downloads_body, project_dir, downloads_file))

  qc_file <- file.path(pages_dir, "qc.html")
  status_files <- c(
    file.path(project_dir, "single_de_status.tsv"),
    file.path(project_dir, "contrast_status.tsv"),
    file.path(project_dir, "post_lisa_status.tsv"),
    file.path(project_dir, "lisa_pipeline_plan.tsv"),
    file.path(project_dir, "lisa_pipeline_metadata.tsv"),
    file.path(project_dir, "reports", "de_status.tsv"),
    file.path(project_dir, "reports", "lisa_single_de_status.tsv"),
    file.path(project_dir, "reports", "lisa_contrast_status.tsv"),
    file.path(project_dir, "state", "de_status.tsv"),
    file.path(project_dir, "state", "contrast_status.tsv"),
    file.path(project_dir, "metrics", "visual_runtime_size.tsv"),
    file.path(project_dir, "manifests", "report_product_validation_summary.tsv"),
    file.path(project_dir, "manifests", "report_contract_validation.tsv"),
    file.path(project_dir, "manifests", "report_link_validation.tsv"),
    file.path(project_dir, "manifests", "figure_format_summary.tsv"),
    file.path(project_dir, "manifests", "figure_format_manifest.tsv")
  )
  qc_body <- paste0(
    '<section class="hero"><h1>QC and manifests</h1><p>Execution status, figure format policy and validation tables. This is the only place where runtime/link-validation banners are shown so they do not distract from biological interpretation.</p></section>',
    '<section id="report-validation" class="panel"><div class="section-head"><div><p class="eyebrow">Report validation</p><h2>Validation summary</h2></div></div>',
    warning_banner,
    validation_banner,
    '<div class="grid">',
    sprintf('<div class="metric"><span>Checked local links</span><strong>%s</strong></div>', if (is.na(checked_links)) "pending final gate" else checked_links),
    sprintf('<div class="metric"><span>Missing local links</span><strong>%s</strong></div>', if (is.na(missing_links)) "pending final gate" else missing_links),
    sprintf('<div class="metric"><span>PNG figures</span><strong>%d</strong></div>', png_count),
    sprintf('<div class="metric"><span>SVG figures</span><strong>%d</strong></div>', svg_count),
    sprintf('<div class="metric"><span>PDF in outputs</span><strong>%d</strong></div>', pdf_count),
    '</div></section>'
  )
  for (sf in status_files[file.exists(status_files)]) {
    df <- read_tsv(sf)
    section_id <- paste0("qc-", slug(basename(sf)))
    qc_body <- paste0(qc_body, '<section id="', esc(section_id), '" class="panel"><div class="section-head"><div><p class="eyebrow">QC</p><h2>', esc(basename(sf)), '</h2></div><a class="file-link" download href="', short_file_copy(sf, qc_file), '">Download</a></div>', table_html(df, NULL, paste0("qc_", slug(basename(sf))), sf, 80), '</section>')
  }
  write_page(qc_file, page_shell(args$title, "qc", qc_body, project_dir, qc_file))
  if (exists("lisa_refresh_evidence_navigation", envir = asNamespace("lisaR"), inherits = FALSE))
    get("lisa_refresh_evidence_navigation", envir = asNamespace("lisaR"))(project_dir)

  # Theme every generated HTML page before computing the portable manifest.
  themed_pages <- c(root_index, list.files(pages_dir, pattern = "[.]html$",
    recursive = TRUE, full.names = TRUE))
  for (page in themed_pages) {
    html <- paste(readLines(page, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    themed <- report_apply_theme(html, page, project_dir)
    if (!identical(html, themed)) write_page(page, themed)
  }

  # Fail closed: a report must never ship a broken expected-layer placeholder
  # nor a generic source-data target. Figure data contracts are checked after
  # all cards have been rendered, so coverage is derived from the report itself.
  report_html <- list.files(project_dir, pattern = "[.]html$", recursive = TRUE, full.names = TRUE)
  report_lines <- unlist(lapply(report_html, readLines, warn = FALSE))
  expected_layer_diagnostics <- report_expected_layer_diagnostics(paste(report_lines, collapse = "\n"))
  if (length(expected_layer_diagnostics)) {
    stop("LISA-REPORT-LAYER-004 required report layer missing: ",
      paste(expected_layer_diagnostics, collapse = " | "), call. = FALSE)
  }
  missing_contract_pattern <- paste(c(
    if (report_requested("source_data", TRUE)) "source data missing" else character(),
    if (report_requested("recipes", TRUE)) "recipe missing" else character()
  ), collapse = "|")
  if (nzchar(missing_contract_pattern) && any(grepl(missing_contract_pattern, report_lines, ignore.case = TRUE))) {
    stop("LISA-REPORT-SOURCE-011 every displayed figure must have figure-specific source data and an executable recipe.", call. = FALSE)
  }
  figure_contracts <- validate_figure_source_contracts(project_dir)
  dir.create(file.path(project_dir, "audit"), recursive = TRUE, showWarnings = FALSE)
  if (nrow(figure_contracts)) utils::write.table(figure_contracts,
    file.path(project_dir, "audit", "figure_source_contracts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  category_coverage <- validate_category_member_plot_coverage(project_dir)
  if (nrow(category_coverage)) utils::write.table(category_coverage,
    file.path(project_dir, "audit", "category_member_plot_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  completion_path <- file.path(project_dir,"audit/full_products/coverage.tsv")
  if (!standard_report && file.exists(completion_path)) {
    completion <- read_tsv(completion_path)
    required <- lisaR:::lisa_full_product_inventory(project_dir,
      kegg_maps = args$kegg_maps == "true" || any(completion$family == "kegg_maps"))
    attached <- report_category_product_reference_rows(project_dir)
    attached_paths <- as.character(attached$source_path)
    rendered_html <- paste(report_lines,collapse="\n")
    for (p in required$path[nzchar(required$path)]) {
      if (p %in% attached_paths) next
      media <- paste0(substr(unname(tools::md5sum(file.path(project_dir,p))),1,16),
        ".", tolower(tools::file_ext(p)))
      if (!grepl(paste0("report_media/",media),rendered_html,fixed=TRUE,useBytes=TRUE))
        stop("Required FULL figure has no report attachment or gallery: ",p,call.=FALSE)
    }
  }

  write_shareable_report_manifest(project_dir)
  validate_shareable_report_manifest(project_dir)

  audit_dir <- file.path(project_dir, "audit")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "# Professional report build",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Root report:", root_index),
    paste("PNG figures:", png_count),
    paste("SVG figures:", svg_count),
    paste("PDF files in outputs:", pdf_count)
  ), file.path(audit_dir, "professional_report_build.md"), useBytes = TRUE)
  if (!lisaR:::lisa_finalize_support_grades_derivation(project_dir))
    lisaR:::lisa_finalize_category_inference_derivation(project_dir)
  message("Wrote root report: ", root_index)
}

if (sys.nframe() == 0L) main()
