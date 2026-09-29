# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Shared presentation context. Scientific keys remain exact; labels are display only.
# Both operands must be canonicalized the same way. `normalizePath()` alone
# leaves a planned (not yet created) route unresolved while the existing page
# directory is resolved, so on macOS `/var` versus `/private/var` and on
# Windows `RUNNER~1` versus the long user name share no common prefix and the
# link degenerates into a `../../..`-chain plus an absolute path.
lisa_presentation_relative <- function(path, page_dir) {
  parts <- function(x) {
    pieces <- strsplit(gsub("\\\\", "/", lisa_path_canonical(x)), "/",
                       fixed = TRUE)[[1L]]
    pieces[nzchar(pieces)]
  }
  to <- parts(path); from <- parts(page_dir)
  n <- 0L
  while (n < min(length(to), length(from)) && identical(to[[n + 1L]], from[[n + 1L]])) n <- n + 1L
  tail <- if (n < length(to)) to[seq.int(n + 1L, length(to))] else character()
  paste(c(rep("..", length(from) - n), tail), collapse = "/")
}

lisa_category_product_href_safe <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^[\\\\/]", value) ||
      grepl("\\\\|[[:cntrl:]]", value)) return(FALSE)
  parts <- strsplit(value, "/", fixed = TRUE)[[1L]]
  safe_part <- function(part) {
    if (!nzchar(part) || part %in% c(".", "..")) return(FALSE)
    decoded <- tryCatch(utils::URLdecode(part), error = function(error) NA_character_)
    length(decoded) == 1L && !is.na(decoded) && nzchar(decoded) &&
      !decoded %in% c(".", "..") &&
      !grepl("[/\\\\]|[[:cntrl:]]", decoded)
  }
  if (all(vapply(parts, safe_part, logical(1L)))) return(TRUE)
  extension_artifact <- length(parts) > 4L &&
    all(parts[seq_len(3L)] == "..") && identical(parts[[4L]], "artifacts") &&
    all(vapply(parts[-seq_len(4L)], safe_part, logical(1L)))
  run_artifact <- length(parts) > 5L &&
    all(parts[seq_len(4L)] == "..") &&
    parts[[5L]] %in% c("outputs", "artifacts") &&
    all(vapply(parts[-seq_len(5L)], safe_part, logical(1L)))
  extension_artifact || run_artifact
}

lisa_presentation_root <- function(output_dir) {
  candidate <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  for (i in seq_len(7L)) {
    if (file.exists(file.path(candidate, "config", "de_index.tsv"))) return(candidate)
    parent <- dirname(candidate)
    if (identical(candidate, parent)) break
    candidate <- parent
  }
  NULL
}

lisa_presentation_routes <- function(project_dir, page_dir, report_build = FALSE) {
  if (is.null(project_dir)) return(list())
  files <- c(overview = "report_index.html", analyses = "report_pages/single_de.html",
    contrasts = "report_pages/contrasts.html", genes = "report_pages/gene_evidence/index.html",
    methods = "report_pages/downloads.html")
  # Evidence pages are written before the gene search and assembled report.
  # Use the persisted product contract, not directory creation order, for
  # future links. With no contract, standalone pages link only to real HTML.
  expected_path <- file.path(project_dir, "expected_artifacts.tsv")
  plan_path <- file.path(project_dir, "lisa_pipeline_plan.tsv")
  expected <- if (file.exists(expected_path)) lisa_evidence_read(expected_path) else data.frame()
  plan <- if (file.exists(plan_path)) lisa_evidence_read(plan_path) else data.frame()
  known <- all(c("stage", "expectation") %in% names(expected)) || "stage" %in% names(plan)
  planned <- function(stage) {
    if (all(c("stage", "expectation") %in% names(expected)))
      return(any(expected$stage == stage & expected$expectation == "required", na.rm = TRUE))
    if ("stage" %in% names(plan)) return(any(plan$stage == stage, na.rm = TRUE))
    FALSE
  }
  available <- stats::setNames(file.exists(file.path(project_dir, files)), names(files))
  report_routes <- c("overview", "analyses", "contrasts", "methods")
  if (known) {
    available[report_routes] <- FALSE
    if (planned("root_html_report")) {
      available[c("overview", "analyses", "methods")] <- TRUE
      available[["contrasts"]] <- planned("category_contrasts")
    }
    available[["genes"]] <- available[["genes"]] || planned("single_de_gene_evidence")
  }
  # A manual assembler has explicit intent to write its pages even when no
  # pipeline contract exists. Merely rendering an evidence page has no such
  # implication. An empty contrast-index file does not promise contrasts.
  if (isTRUE(report_build)) {
    available[c("overview", "analyses", "methods")] <- TRUE
    index_path <- file.path(project_dir, "config", "contrast_index.tsv")
    index <- if (file.exists(index_path) && file.info(index_path)$size > 2)
      lisa_evidence_read(index_path) else data.frame()
    available[["contrasts"]] <- available[["contrasts"]] ||
      ("contrast_id" %in% names(index) && nrow(index) > 0L)
  }
  files <- files[available]
  if (!length(files)) return(list())
  as.list(vapply(files, function(x) lisa_presentation_relative(file.path(project_dir, x), page_dir), character(1L)))
}

lisa_presentation_label <- function(project_dir, id, contrast = FALSE) {
  if (is.null(project_dir) || is.null(id) || length(id) != 1L || !nzchar(id)) return(id)
  path <- file.path(project_dir, "config", if (contrast) "contrast_index.tsv" else "de_index.tsv")
  if (!file.exists(path)) return(id)
  tab <- lisa_evidence_read(path)
  key <- if (contrast) {
    if ("contrast_id" %in% names(tab)) "contrast_id" else "output_id"
  } else "analysis_id"
  if (!key %in% names(tab)) return(id)
  row <- tab[as.character(tab[[key]]) == id, , drop = FALSE]
  if (nrow(row) != 1L) return(id)
  label_columns <- c("display_title", if (contrast) "contrast_label", "label", "title", "analysis_label", "comparison_name")
  for (name in label_columns) {
    if (name %in% names(row) && !is.na(row[[name]][[1L]]) && nzchar(row[[name]][[1L]])) return(as.character(row[[name]][[1L]]))
  }
  id
}

lisa_present_evidence_html <- function(html, output_dir, metadata = list(), active = "analyses", main_id = NULL) {
  html <- paste(html, collapse = "\n")
  project <- lisa_presentation_root(output_dir)
  id <- metadata$analysis_id
  if (identical(active, "contrasts")) id <- metadata$contrast_id
  exact_ids <- stats::setNames(list(id), if (identical(active, "contrasts")) "contrast_id" else "analysis_id")
  for (key in c("analysis_a", "analysis_b", "collection", "tier"))
    if (!is.null(metadata[[key]])) exact_ids[[key]] <- metadata[[key]]
  context <- list(analysis = lisa_presentation_label(project, id, identical(active, "contrasts")),
    collection = metadata$collection, direction = metadata$positive_contrast,
    cutoff = if (!is.null(metadata$gsea_padj_cutoff)) paste("GSEA FDR <=", metadata$gsea_padj_cutoff) else NULL,
    exact_ids = exact_ids)
  routes <- lisa_presentation_routes(project, output_dir)
  if (is.null(main_id)) {
    tag <- regmatches(html, regexpr("<main\\b[^>]*>", html, perl = TRUE))
    id_match <- regexec('\\bid=["\x27]([^"\x27]+)["\x27]', tag, perl = TRUE)
    found <- regmatches(tag, id_match)[[1L]]
    main_id <- if (length(found) > 1L) found[[2L]] else "lisa-main"
  }
  shell <- .lisa_report_shell(routes, active = if (active %in% names(routes)) active else NULL,
    context = context, assets = .lisa_copy_report_shell_assets(output_dir), main_id = main_id)
  html <- .lisa_inject_report_shell(html, shell, add_main = FALSE)
  if (!is.null(metadata$hommel_support_schema) || grepl("data-hommel-help", html, fixed = TRUE))
    html <- lisa_hommel_install_help(html)
  lisa_attach_evidence_navigation(html, project, output_dir, metadata, active)
}

# The assembler publishes root-relative routes after it knows actual sections.
# Standalone viewers consume the same inventory, never a guessed second menu.
lisa_attach_evidence_navigation <- function(html, project, page_dir,
    metadata = list(), active = "analyses") {
  if (is.null(project)) return(html)
  path <- file.path(project, "report_pages", "navigation_inventory.json")
  if (!file.exists(path)) return(html)
  nav <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!identical(nav$version, 1L) && !identical(nav$version, 1))
    stop("Unsupported report navigation inventory version.", call. = FALSE)
  for (i in seq_along(nav$contexts)) {
    for (j in seq_along(nav$contexts[[i]]$collections)) {
      sections <- nav$contexts[[i]]$collections[[j]]$sections
      for (k in seq_along(sections)) {
        href <- sections[[k]]$href
        if (!is.character(href) || length(href) != 1L ||
            grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^/|(^|/)\\.\\.(/|$)", href))
          stop("Report navigation must use contained relative routes.", call. = FALSE)
        pieces <- strsplit(href, "#", fixed = TRUE)[[1L]]
        sections[[k]]$href <- paste0(lisa_presentation_relative(
          file.path(project, pieces[[1L]]), page_dir),
          if (length(pieces) > 1L) paste0("#", paste(pieces[-1L], collapse = "#")) else "")
      }
      nav$contexts[[i]]$collections[[j]]$sections <- sections
    }
  }
  kind <- if (identical(active, "contrasts")) "contrast" else "analysis"
  id <- if (kind == "contrast") metadata$contrast_id else metadata$analysis_id
  matched <- which(vapply(nav$contexts, function(x) identical(x$kind, kind) &&
    (identical(x$id, id) || identical(x$scientific_id, id)), logical(1L)))
  nav$mode <- "external"
  nav$current <- list(context = if (length(matched) == 1L) nav$contexts[[matched]]$id else "",
    collection = if (is.null(metadata$collection)) "" else metadata$collection,
    section = "category-evidence")
  nav$page <- active
  encoded <- jsonlite::toJSON(nav, auto_unbox = TRUE, null = "null", digits = NA)
  encoded <- gsub("<", "\\u003c", encoded, fixed = TRUE)
  html <- gsub('(?s)<script[^>]*id="lisa-report-navigation"[^>]*>.*?</script>', "", html, perl = TRUE)
  marker <- regexpr("(?i)</body\\s*>", html, perl = TRUE)
  if (marker[[1L]] < 0L) stop("Evidence HTML body is missing.", call. = FALSE)
  payload <- paste0('<script id="lisa-report-navigation" type="application/json">', encoded, '</script>')
  paste0(substr(html, 1L, marker[[1L]] - 1L), payload, substr(html, marker[[1L]], nchar(html)))
}

lisa_refresh_evidence_navigation <- function(project_dir) {
  folders <- file.path(project_dir, "report_pages",
    c("evidence", "category_navigation", "contrast_evidence", "gene_evidence"))
  pages <- unlist(lapply(folders[dir.exists(folders)], list.files,
    pattern = "^index\\.html$", recursive = TRUE, full.names = TRUE), use.names = FALSE)
  for (page in pages) {
    relative <- lisa_presentation_relative(page, file.path(project_dir, "report_pages"))
    parts <- strsplit(relative, "/", fixed = TRUE)[[1L]]
    active <- if (parts[[1L]] == "contrast_evidence") "contrasts" else if (parts[[1L]] == "gene_evidence") "genes" else "analyses"
    metadata <- list()
    if (length(parts) == 4L) {
      metadata[[if (active == "contrasts") "contrast_id" else "analysis_id"]] <- parts[[2L]]
      metadata$collection <- parts[[3L]]
    }
    html <- paste(readLines(page, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    updated <- lisa_attach_evidence_navigation(html, project_dir, dirname(page), metadata, active)
    if (!identical(html, updated)) lisa_guarded_write(page, function(target)
      writeLines(updated, target, useBytes = TRUE))
  }
  general <- c(file.path(project_dir, "report_index.html"),
    list.files(file.path(project_dir, "report_pages"), pattern = "[.]html$", full.names = TRUE))
  general <- general[file.exists(general) & !basename(general) %in% c("single_de.html", "contrasts.html")]
  for (page in general) {
    html <- paste(readLines(page, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (!grepl("data-lisa-shell", html, fixed = TRUE)) next
    updated <- lisa_attach_evidence_navigation(html, project_dir, dirname(page), active = "overview")
    if (!identical(html, updated)) lisa_guarded_write(page, function(target)
      writeLines(updated, target, useBytes = TRUE))
  }
  invisible(pages)
}
