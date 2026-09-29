# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# Static export.
#
# The export is a pure read-and-assemble operation. It copies the complete
# STANDARD report plus every extension that was validated at the moment the
# export started, and it never invokes a generator or `--complete-missing`.
# The bundle is placed so the STANDARD report keeps its own
# relative links: the run tree lands at the bundle root, exactly as it was, and
# the extensions live beside it under `explore_extensions/`.
#
# This deliberately does NOT build an alternative gallery. The source report
# remains the presentation; `report_presentation.tsv` records the exact route
# and anchor at which each extension attaches, for the interface to use.
# ---------------------------------------------------------------------------

# Where a product attaches in the existing report.
#
# The previous answer was the per-collection navigator,
# `report_pages/category_navigation/<owner>/<collection>/index.html`, with the
# category id recorded as an anchor that nothing ever used for placement. That
# page is a *diagram of all the categories in a collection*: each
# `<a class="category-row" data-category-id="...">` links onward to the sheet the
# reader actually opens. Attaching there stacked every category's figure and
# control at the foot of one shared page, so the sheet for category A was just as
# likely to be sitting under the diagram next to B's, and the category the reader
# had opened received nothing at all.
#
# The reader-facing page is the evidence sheet the navigator links to:
#
#   report_pages/evidence/<analysis>/<collection>/index.html?category=<ID>
#   report_pages/contrast_evidence/<contrast>/<collection>/index.html?category=<ID>
#
# which is exactly the href the navigator emits, and which the sheet's own script
# mirrors back into `?category=` plus `#category-<ID>` as the reader switches
# category. So route + query + anchor below are the report's real routing, not a
# convention invented here.
#
# The navigator remains the fallback for a run whose evidence sheet is absent, so
# a figure is never silently dropped; `surface` records which of the two was
# used, and the presentation table carries the query and anchor so the shell and
# the static bundle scope a block to one category rather than to a page.
lisa_explore_attachment_surface <- function(unit_type) {
  if (identical(unit_type, "single_de")) "evidence" else "contrast_evidence"
}

lisa_explore_attachment_page <- function(unit_type) {
  if (identical(unit_type, "single_de")) "report_pages/single_de.html"
  else "report_pages/contrasts.html"
}

lisa_explore_attachment <- function(request, source_run = NULL) {
  owner <- if (identical(request$unit_type, "single_de")) {
    request$analysis_id
  } else {
    request$contrast_id
  }
  collection <- as.character(request$collection)
  if (identical(lisa_explore_product_scope(request$product), "collection")) {
    route <- lisa_explore_attachment_page(request$unit_type)
    if (!is.null(source_run) && !file.exists(file.path(source_run, route))) {
      # This run has no native collection page. Reported rather than invented:
      # the row keeps its route in `report_presentation.tsv` and the exporter
      # skips a page it does not have.
      return(list(route = route, anchor = "", query = "", category_id = "",
                  surface = "native_collection", owner = owner,
                  scope = "collection", context = owner,
                  collection = collection))
    }
    return(list(route = route, anchor = "", query = "", category_id = "",
                surface = "native_collection", owner = owner,
                scope = "collection", context = owner, collection = collection))
  }
  surface <- lisa_explore_attachment_surface(request$unit_type)
  route <- file.path("report_pages", surface, owner, request$collection,
                     "index.html")
  if (!is.null(source_run) && !file.exists(file.path(source_run, route))) {
    # This run has no evidence sheet for the collection. Fall back to the
    # navigator rather than naming a page that is not in the bundle.
    surface <- "category_navigation"
    route <- file.path("report_pages", surface, owner, request$collection,
                       "index.html")
  }
  list(route = route, anchor = paste0("category-", request$category_id),
       query = paste0("category=", utils::URLencode(as.character(request$category_id),
                                                    reserved = TRUE)),
       category_id = as.character(request$category_id),
       surface = surface, owner = owner,
       scope = "category", context = owner, collection = collection)
}

lisa_explore_relative_paths <- function(root) {
  tree <- lisa_scan_run_tree(root)
  paths <- sort(tree$path[!tree$isdir])
  normalized <- normalizePath(root, winslash = "/", mustWork = TRUE)
  substring(paths, nchar(normalized) + 2L)
}

# --- visible static attachment (finding 1) ---------------------------------
#
# The previous export wrote `report_presentation.tsv` and stopped there. The
# bundle's HTML was the verbatim STANDARD copy, so a user opening the exported
# report saw no new figure at all -- it was reachable only by reading a TSV and
# following paths by hand. A route table is not a report.
#
# This hook edits the copied category pages *inside staging*, before the manifest
# and receipt are hashed, so the manifest describes what the user will actually
# open. It is a pure text transformation over already-validated artifacts:
# `lisa_explore_presentation_attach_html()` renders only rows handed to it, and
# nothing here plans, renders, or completes a missing product.

lisa_explore_export_stylesheet <- "report_assets/lisa-explore.css"
lisa_explore_export_script <- "report_assets/lisa-explore.js"

lisa_explore_attach_head <- function(html, route, scoped) {
  if (grepl('data-lisa-explore-style="true"', html, fixed = TRUE)) return(html)
  head <- paste0(
    '<link rel="stylesheet" href="',
    lisa_explore_presentation_relative_href(lisa_explore_export_stylesheet, route),
    '" data-lisa-explore-style="true">',
    # The reveal script is added only where blocks are category-scoped, so a page
    # that shows every category at once carries no script it does not need.
    if (isTRUE(scoped)) paste0(
      '<script defer src="',
      lisa_explore_presentation_relative_href(lisa_explore_export_script, route),
      '" data-lisa-explore-script="true"></script>') else "")
  marker <- regexpr("(?i)</head\\s*>", html, perl = TRUE)
  if (marker[[1L]] < 0L) return(paste0(head, html))
  paste0(substr(html, 1L, marker[[1L]] - 1L), head,
         substr(html, marker[[1L]], nchar(html)))
}

lisa_explore_bundle_asset <- function(staging, asset_name, relative) {
  asset <- system.file("shiny_assets", asset_name, package = "lisaR")
  if (!nzchar(asset) || !file.exists(asset)) {
    stop("LISA-EXPLORE-038 the exploration asset '", asset_name, "' is missing ",
         "from the installed package; the attached figures would not render ",
         "correctly.", call. = FALSE)
  }
  target <- file.path(staging, relative)
  if (lisa_path_entry_exists(target)) return(invisible(relative))
  lisa_guarded_dir_create(dirname(target), run_root = staging)
  lisa_guarded_copy(asset, target, overwrite = FALSE, run_root = staging)
  invisible(relative)
}

lisa_explore_attach_presentation <- function(staging, presentation) {
  attached <- character()
  if (!nrow(presentation)) {
    return(list(pages = attached, stylesheet = NA_character_,
                script = NA_character_))
  }
  lisa_explore_bundle_asset(staging, "explore.css", lisa_explore_export_stylesheet)

  script_used <- FALSE
  routes <- unique(presentation$attach_route[nzchar(presentation$attach_route)])
  for (route in routes) {
    page <- file.path(staging, route)
    if (!file.exists(page)) {
      # The STANDARD report does not contain this route. We do not invent a page
      # for it; the row stays in report_presentation.tsv and is reported.
      next
    }
    original <- paste(readLines(page, warn = FALSE), collapse = "\n")
    scoped <- lisa_explore_presentation_scoped(original)
    updated <- lisa_explore_presentation_attach_html(original, presentation, route)
    updated <- lisa_explore_attach_head(updated, route, scoped)
    if (identical(updated, original)) next
    if (isTRUE(scoped) && !script_used) {
      lisa_explore_bundle_asset(staging, "explore-static.js",
                                lisa_explore_export_script)
      script_used <- TRUE
    }
    lisa_guarded_write(page, function(target) {
      writeLines(updated, target, useBytes = TRUE)
    }, run_root = staging, overwrite = TRUE)
    attached <- c(attached, route)
  }
  list(pages = attached, stylesheet = lisa_explore_export_stylesheet,
       script = if (script_used) lisa_explore_export_script else NA_character_)
}


lisa_explore_absolute_reference_pattern <- function() paste0(
  "(?i)(?:href|src|srcset|action|data|poster)\\s*=\\s*[\"']\\s*",
  "(?:[A-Za-z][A-Za-z0-9+.-]*:|//|/|\\\\\\\\)")

lisa_explore_absolute_css_pattern <- function() paste0(
  "(?i)url\\(\\s*[\"']?\\s*(?:[A-Za-z][A-Za-z0-9+.-]*:|//|/|\\\\\\\\)")

# Absolute values that are not local paths, so relocation does not break them.
lisa_explore_allowed_reference_pattern <- function() paste0(
  "(?i)(?:href|src|srcset|action|data|poster)\\s*=\\s*[\"']\\s*",
  "(?:https?:|mailto:|tel:|data:|#)")

lisa_explore_allowed_css_pattern <- function() "(?i)url\\(\\s*[\"']?\\s*(?:https?:|data:)"

lisa_explore_absolute_references <- function(text, css = FALSE) {
  pattern <- if (isTRUE(css)) lisa_explore_absolute_css_pattern() else
    lisa_explore_absolute_reference_pattern()
  allowed <- if (isTRUE(css)) lisa_explore_allowed_css_pattern() else
    lisa_explore_allowed_reference_pattern()
  # Remove the accepted external forms first so they cannot satisfy the generic
  # scheme branch of the pattern, then ask whether anything absolute remains.
  remaining <- gsub(allowed, "", text, perl = TRUE)
  any(grepl(pattern, remaining, perl = TRUE))
}

# Enumerated with list.files rather than lisa_scan_run_tree, and evaluated BEFORE
# the guarded scan is used for anything else. That ordering is what makes the
# symlink answer truthful: the guarded scan throws on the first symlink it meets,
# so while the report ran after it the symlink branch was unreachable and
# "symlinks are enforced" held only because the guarded copy layer cannot create
# one. Running the report first makes it a real check -- a symlink that reached
# staging by any route is named here and refused by the caller -- and the caller
# below no longer needs the guarded scan to have succeeded in order to refuse.
lisa_explore_portability_report <- function(staging) {
  entries <- list.files(staging, recursive = TRUE, all.files = TRUE,
                        no.. = TRUE, include.dirs = TRUE)
  targets <- Sys.readlink(file.path(staging, entries))
  links <- entries[!is.na(targets) & nzchar(targets)]

  absolute <- character()
  external <- character()
  for (relative in setdiff(entries, links)) {
    is_html <- grepl("[.]x?html?$", relative, perl = TRUE)
    is_css <- grepl("[.]css$", relative, perl = TRUE)
    if (!is_html && !is_css) next
    if (dir.exists(file.path(staging, relative))) next
    text <- readLines(file.path(staging, relative), warn = FALSE)
    if (lisa_explore_absolute_references(text, css = is_css)) {
      absolute <- c(absolute, relative)
    }
    if (any(grepl("(?i)(?:href|src)\\s*=\\s*[\"']\\s*https?:|(?i)url\\(\\s*[\"']?https?:",
                  text, perl = TRUE))) {
      external <- c(external, relative)
    }
  }
  list(absolute_reference_files = absolute, symlinks = links,
       external_reference_files = external)
}

lisa_explore_export_destination <- function(ws, destination) {
  destination <- lisa_explore_scalar(destination, "destination")
  resolved <- lisa_managed_destination(destination, create_parent = TRUE)
  output_dir <- resolved$path
  if (lisa_path_within(output_dir, ws$source_run) ||
      identical(lisa_path_key(output_dir), lisa_path_key(ws$source_run))) {
    stop("LISA-EXPLORE-030 the export must be written outside the immutable source run.",
         call. = FALSE)
  }
  if (lisa_path_entry_exists(output_dir)) {
    stop("LISA-EXPLORE-031 the export destination already exists: ", output_dir,
         call. = FALSE)
  }
  resolved
}

lisa_explore_export_staging_path <- function(output_dir, staging_id = NULL) {
  suffix <- if (is.null(staging_id)) "" else {
    staging_id <- lisa_explore_scalar(staging_id, "staging_id")
    if (!grepl("^[0-9a-f]{24}$", staging_id)) {
      stop("Unsafe export staging identity.", call. = FALSE)
    }
    paste0("-", staging_id)
  }
  lisa_short_staging_path(output_dir, paste0("export", suffix))
}

lisa_explore_export_callback <- function(callback, ...) {
  if (is.function(callback)) callback(...)
  invisible(NULL)
}

lisa_explore_export_cancelled <- function(callback, stage) {
  if (is.function(callback) && isTRUE(callback())) {
    stop("LISA-EXPLORE-079 export cancelled at the safe point before ", stage,
         "; no destination was published.", call. = FALSE)
  }
  invisible(FALSE)
}

# Copy one already-validated regular tree in bounded vectorized chunks. The old
# exporter called lisa_guarded_copy() once per file (and redundantly created each
# parent before that call); on reports with thousands of tiny assets, repeated
# canonicalization and tree-boundary walks dominated minutes of wall time. This
# helper keeps the same safety properties at tree granularity:
#
# * source and destination roots are guarded and symlink-free;
# * target parents are created through guarded IO before any batch;
# * each source is hashed before copying, and source + destination hashes must
#   both still match after copying;
# * a final guarded scan rejects any unexpected entry or symlink.
lisa_explore_guarded_copy_tree <- function(source_root, staging,
                                           destination_prefix = "",
                                           progress = NULL,
                                           should_cancel = NULL,
                                           stage = "copying files",
                                           chunk_size = 128L) {
  source_root <- lisa_existing_run_root(source_root)
  staging <- lisa_existing_run_root(staging)
  tree <- lisa_scan_run_tree(source_root)
  files <- tree$path[!tree$isdir]
  root_key <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
  relative <- if (length(files))
    substring(files, nchar(root_key) + 2L) else character()
  destination_relative <- file.path(destination_prefix, relative)
  targets <- file.path(staging, destination_relative)
  if (anyDuplicated(destination_relative)) {
    stop("LISA-EXPLORE-077 a guarded export tree contains duplicate paths.",
         call. = FALSE)
  }
  blocked <- targets[vapply(targets, lisa_path_entry_exists, logical(1L))]
  if (length(blocked)) {
    stop("LISA-EXPLORE-077 export copy would overwrite an existing path: ",
         blocked[[1L]], call. = FALSE)
  }
  directories <- unique(dirname(targets))
  directories <- directories[order(nchar(directories), directories)]
  for (directory in directories) {
    lisa_guarded_dir_create(directory, run_root = staging)
  }
  if (!length(files)) return(data.frame(
    path = character(), bytes = numeric(), sha256 = character(),
    stringsAsFactors = FALSE))

  chunks <- split(seq_along(files), ceiling(seq_along(files) / as.integer(chunk_size)))
  records <- vector("list", length(chunks))
  copied_count <- 0L
  for (chunk_index in seq_along(chunks)) {
    lisa_explore_export_cancelled(should_cancel, stage)
    selected <- chunks[[chunk_index]]
    from <- files[selected]
    to <- targets[selected]
    if (any(vapply(from, lisa_path_is_link, logical(1L))) ||
        any(!file.exists(from)) || any(dir.exists(from)) ||
        !all(vapply(from, function(path) utils::file_test("-f", path),
                    logical(1L)))) {
      stop("LISA-EXPLORE-077 guarded export source changed before copying.",
           call. = FALSE)
    }
    hashes <- vapply(from, lisa_sha256_file, character(1L))
    bytes <- as.numeric(file.info(from)$size)
    copied <- file.copy(from, to, overwrite = FALSE,
                        copy.mode = FALSE, copy.date = FALSE)
    if (length(copied) != length(from) || !all(copied)) {
      failed <- relative[selected][!copied]
      stop("LISA-EXPLORE-077 guarded export copy failed for: ",
           paste(utils::head(failed, 3L), collapse = ", "), call. = FALSE)
    }
    if (any(vapply(from, lisa_path_is_link, logical(1L))) ||
        any(vapply(to, lisa_path_is_link, logical(1L))) ||
        !all(vapply(to, function(path) utils::file_test("-f", path),
                    logical(1L)))) {
      stop("LISA-EXPLORE-077 guarded export copy produced or observed a ",
           "non-regular entry.", call. = FALSE)
    }
    source_after <- vapply(from, lisa_sha256_file, character(1L))
    target_after <- vapply(to, lisa_sha256_file, character(1L))
    if (!identical(unname(hashes), unname(source_after)) ||
        !identical(unname(hashes), unname(target_after))) {
      stop("LISA-EXPLORE-078 source or copied bytes changed during export; ",
           "the snapshot was not published.", call. = FALSE)
    }
    records[[chunk_index]] <- data.frame(
      path = destination_relative[selected], bytes = bytes,
      sha256 = unname(hashes), stringsAsFactors = FALSE)
    copied_count <- copied_count + length(selected)
    lisa_explore_export_callback(progress, stage,
      paste("copied and verified", copied_count, "of", length(files), "files"))
  }
  lisa_scan_run_tree(staging)
  do.call(rbind, records)
}

#' Export a static report bundle
#'
#' Copies the complete STANDARD report plus every exploration artifact that is
#' validated when the export begins. No figure is generated and no missing
#' product is completed. The result is self-contained: it opens from any
#' location with R and Shiny closed.
#'
#' @param ws A workspace from [lisa_explore_open()].
#' @param destination Directory to create for the bundle. It must not exist.
#' @param format `"dir"` for a folder, `"zip"` to also produce a ZIP archive.
#' @param allow_absolute_references Promote the bundle even when it contains
#'   absolute HTML references or symlinks. The default `FALSE` refuses instead,
#'   because such a bundle does not open from another location.
#' @param snapshot_artifacts Optional validated artifact snapshot used by the
#'   background exporter. Leave as `NULL` to snapshot currently valid artifacts.
#' @param progress Optional callback accepting a stage and message; used to
#'   report export progress without changing the selected artifact snapshot.
#' @param should_cancel Optional zero-argument callback returning `TRUE` to
#'   request cooperative cancellation at an export checkpoint.
#' @param staging_id Internal per-job staging identity used by the background
#'   exporter. Interactive callers should leave this as `NULL`.
#' @return Invisibly, a list describing the exported bundle.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   lisa_explore_export(ws, file.path(ws$root, "exports", "bundle-001"))
#' }
lisa_explore_export <- function(ws, destination, format = c("dir", "zip"),
                                allow_absolute_references = FALSE,
                                snapshot_artifacts = NULL,
                                progress = NULL,
                                should_cancel = NULL,
                                staging_id = NULL) {
  lisa_explore_assert_workspace(ws)
  format <- match.arg(format)
  resolved <- lisa_explore_export_destination(ws, destination)
  output_dir <- resolved$path

  # Snapshot taken once, now. Figures that finish later belong to the next
  # export, never to a half-written bundle.
  available_now <- lisa_explore_artifacts(ws)
  catalogue_status <- NULL
  if (is.null(snapshot_artifacts)) {
    available <- available_now
    catalogue_status <- lisa_explore_status(ws)
  } else {
    snapshot <- lisa_explore_snapshot_rows(snapshot_artifacts)
    if (!nrow(snapshot)) {
      available <- available_now[0, , drop = FALSE]
    } else {
      required <- c("key", "request_id", "unit_type", "analysis_id",
                    "contrast_id", "collection", "category_id", "product",
                    "entity", "variant", "associations", "artifact_dir")
      missing <- setdiff(required, names(snapshot))
      if (length(missing) || anyDuplicated(as.character(snapshot$key))) {
        stop("LISA-EXPLORE-078 the persisted export snapshot is malformed; ",
             "no destination was published.", call. = FALSE)
      }
      matched <- match(as.character(snapshot$key), as.character(available_now$key))
      if (anyNA(matched)) {
        stop("LISA-EXPLORE-078 one or more snapshot artifacts are no longer ",
             "valid; no partial export was published: ",
             paste(as.character(snapshot$key[is.na(matched)]), collapse = ", "),
             call. = FALSE)
      }
      current <- available_now[matched, , drop = FALSE]
      identity <- c("request_id", "unit_type", "analysis_id", "contrast_id",
                    "collection", "category_id", "product", "entity", "variant",
                    "artifact_dir")
      consistent <- vapply(identity, function(column) {
        identical(lisa_explore_blank(snapshot[[column]]),
                  lisa_explore_blank(current[[column]]))
      }, logical(1L))
      if (!all(consistent)) {
        stop("LISA-EXPLORE-078 snapshot artifact identity changed; no destination ",
             "was published.", call. = FALSE)
      }
      available <- snapshot
    }
  }
  started_at <- lisa_explore_now()

  manifest_path <- file.path(ws$source_run, "run_manifest.tsv")
  if (!file.exists(manifest_path) ||
      !identical(lisa_sha256_file(manifest_path), ws$source_manifest_hash)) {
    stop("LISA-EXPLORE-078 the STANDARD run manifest changed after this ",
         "workspace was opened; no snapshot was published.", call. = FALSE)
  }
  lisa_explore_export_callback(progress, "validating snapshot",
                               paste("validated", nrow(available),
                                     "available figure(s) in the snapshot"))
  lisa_explore_export_cancelled(should_cancel, "STANDARD copy")

  staging <- lisa_explore_export_staging_path(output_dir, staging_id)
  if (is.null(staging_id) && lisa_path_entry_exists(staging)) {
    stop("LISA-EXPLORE-032 the deterministic export staging directory already exists.",
         call. = FALSE)
  }
  if (is.null(staging_id)) {
    staging <- lisa_run_root(staging)
  } else {
    owner <- lisa_explore_read_json(
      lisa_explore_export_owner_path(output_dir, staging_id))
    if (!lisa_path_entry_exists(staging) || is.null(owner) ||
        !identical(as.character(owner$token), as.character(staging_id)) ||
        !identical(lisa_path_key(as.character(owner$destination)),
                   lisa_path_key(output_dir)) ||
        length(list.files(staging, all.files = TRUE, no.. = TRUE))) {
      stop("LISA-EXPLORE-080 the isolated export staging claim could not be ",
           "verified. It was retained and no destination was published: ",
           staging, call. = FALSE)
    }
    staging <- lisa_existing_run_root(staging)
  }
  promoted <- FALSE
  on.exit(if (!promoted && lisa_path_entry_exists(staging)) {
    lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  }, add = TRUE)

  # --- the complete STANDARD, verbatim -------------------------------------
  standard_records <- lisa_explore_guarded_copy_tree(
    ws$source_run, staging, progress = progress, should_cancel = should_cancel,
    stage = "copying STANDARD")
  standard_files <- standard_records$path

  # --- the extensions validated at snapshot time ---------------------------
  presentation <- list()
  extension_files <- 0L
  if (nrow(available)) {
    for (index in seq_len(nrow(available))) {
      row <- available[index, , drop = FALSE]
      source_dir <- file.path(ws$root, "extensions", row$artifact_dir)
      lisa_explore_export_cancelled(should_cancel, "extension copy")
      extension_records <- lisa_explore_guarded_copy_tree(
        source_dir, staging, file.path("explore_extensions", row$key),
        progress = progress, should_cancel = should_cancel,
        stage = "copying available extensions")
      relatives <- substring(extension_records$path,
        nchar(file.path("explore_extensions", row$key)) + 2L)
      extension_files <- extension_files + length(relatives)
      bundle_relative <- file.path("explore_extensions", row$key, relatives)
      # Product-aware selection, shared with the live shell. For a native pathway
      # map this is what keeps the preserved unpainted base diagram from being
      # exported as though it were the figure.
      pick <- function(role) {
        lisa_explore_pick_product_file(bundle_relative, row$product, role)
      }
      associations <- if (identical(lisa_explore_product_scope(row$product),
                                    "collection")) "" else if (!is.null(snapshot_artifacts)) {
        values <- strsplit(lisa_explore_blank(row$associations)[[1L]], ";",
                           fixed = TRUE)[[1L]]
        sort(unique(values[nzchar(values)]))
      } else {
        lisa_explore_catalogue_associations(catalogue_status, row$key,
                                            fallback = row$category_id)
      }
      for (category in associations) {
        request <- lisa_explore_request_from_row(row)
        request$category_id <- if (nzchar(category)) category else NA_character_
        # Resolved against the *staged* report, so the route named in the table is
        # a page that is actually in the bundle and the navigator fallback is only
        # used when the evidence sheet genuinely is not there.
        attachment <- lisa_explore_attachment(request, staging)
        presentation[[length(presentation) + 1L]] <- data.frame(
          key = row$key, unit_type = row$unit_type,
          analysis_id = row$analysis_id, contrast_id = row$contrast_id,
          collection = row$collection, category_id = category,
          product = row$product,
          entity = lisa_explore_blank(row$entity),
          variant = lisa_explore_blank(row$variant),
          attach_route = attachment$route, attach_anchor = attachment$anchor,
          attach_query = attachment$query, attach_surface = attachment$surface,
          attach_scope = attachment$scope, attach_context = attachment$context,
          attach_collection = attachment$collection,
          png = pick("png"), pdf = pick("pdf"),
          source_data = pick("source_data"), recipe = pick("recipe"),
          stringsAsFactors = FALSE)
      }
    }
  }

  presentation <- if (length(presentation)) {
    do.call(rbind, presentation)
  } else {
    data.frame(key = character(), unit_type = character(),
               analysis_id = character(), contrast_id = character(),
               collection = character(), category_id = character(),
               product = character(), entity = character(),
               variant = character(), attach_route = character(),
               attach_anchor = character(), attach_query = character(),
               attach_surface = character(), attach_scope = character(),
               attach_context = character(), attach_collection = character(),
               png = character(), pdf = character(),
               source_data = character(), recipe = character(),
               stringsAsFactors = FALSE)
  }
  old_run_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = staging)
  on.exit(options(lisaR.run_root = old_run_root), add = TRUE)
  write_lisa_tsv(presentation, file.path(staging, "report_presentation.tsv"))
  lisa_explore_export_callback(progress, "attaching figures",
                               paste("attaching", nrow(presentation),
                                     "figure placement(s)"))
  lisa_explore_export_cancelled(should_cancel, "figure attachment")

  # --- make the figures actually visible in the report ---------------------
  # Runs before the manifest and receipt are hashed, so they describe the pages
  # the user will open. Touches only `staging`; `ws$source_run` is never opened
  # for writing here or anywhere else in this function.
  attachment <- lisa_explore_attach_presentation(staging, presentation)

  # --- portability checks, enforced rather than merely counted -------------
  # Evaluated before `lisa_explore_relative_paths()`, which uses the guarded scan
  # and throws on a symlink: running the report first is what lets the refusal
  # below actually name a symlink instead of dying on it.
  lisa_explore_export_callback(progress, "checking portability",
                               "checking local references and links")
  lisa_explore_export_cancelled(should_cancel, "portability checks")
  portability <- lisa_explore_portability_report(staging)
  absolute_reference_files <- portability$absolute_reference_files
  if (!isTRUE(allow_absolute_references) &&
      (length(absolute_reference_files) || length(portability$symlinks))) {
    # Refuse to promote. `on.exit` removes the staging tree, so a bundle that
    # would not survive relocation never appears at the destination at all.
    stop("LISA-EXPLORE-037 the bundle is not relocatable and was not written. ",
         if (length(absolute_reference_files))
           paste0(length(absolute_reference_files),
                  " HTML file(s) reference absolute paths (",
                  paste(utils::head(absolute_reference_files, 3L), collapse = ", "),
                  "). ") else "",
         if (length(portability$symlinks))
           paste0(length(portability$symlinks), " symlink(s) remain (",
                  paste(utils::head(portability$symlinks, 3L), collapse = ", "),
                  "). ") else "",
         "Fix the report, or pass allow_absolute_references = TRUE to accept a ",
         "bundle that only opens from its original machine.", call. = FALSE)
  }

  lisa_explore_export_callback(progress, "building manifest",
                               "hashing the complete staged bundle")
  lisa_explore_export_cancelled(should_cancel, "manifest hashing")
  bundle_files <- lisa_explore_relative_paths(staging)
  manifest <- data.frame(
    path = bundle_files,
    bytes = as.numeric(file.info(file.path(staging, bundle_files))$size),
    sha256 = vapply(file.path(staging, bundle_files), lisa_sha256_file,
                    character(1)),
    stringsAsFactors = FALSE)
  write_lisa_tsv(manifest, file.path(staging, "explore_export_manifest.tsv"))

  receipt <- data.frame(
    schema = "lisa-explore-export/1",
    source_run = basename(ws$source_run),
    source_manifest_hash = ws$source_manifest_hash,
    started_at = started_at, completed_at = lisa_explore_now(),
    standard_files = length(standard_files),
    extension_files = extension_files,
    extensions = nrow(available),
    generators_invoked = 0L,
    complete_missing = "false",
    attached_pages = length(attachment$pages),
    html_files_with_absolute_references = length(absolute_reference_files),
    symlinks = length(portability$symlinks),
    # Recorded, not refused: an external URL is a network dependency, not a
    # relocation failure. Naming it keeps `relocatable` honest about what it does
    # and does not promise.
    files_with_external_references = length(portability$external_reference_files),
    relocatable = as.character(!length(absolute_reference_files) &&
                                 !length(portability$symlinks)),
    absolute_references_allowed = as.character(isTRUE(allow_absolute_references)),
    lisaR_version = as.character(utils::packageVersion("lisaR")),
    stringsAsFactors = FALSE)
  write_lisa_tsv(receipt, file.path(staging, "explore_export_receipt.tsv"))
  options(lisaR.run_root = old_run_root)

  lisa_explore_export_callback(progress, "final verification",
                               "verifying the complete bundle before atomic publication")
  lisa_explore_export_cancelled(should_cancel, "publication")
  lisa_promote_managed_directory(staging, output_dir)
  promoted <- TRUE

  archive <- NA_character_
  if (identical(format, "zip")) {
    archive <- lisa_explore_zip_bundle(output_dir)
  }

  invisible(list(output_dir = output_dir, archive = archive,
                 standard_files = length(standard_files),
                 extension_files = extension_files,
                 extensions = nrow(available),
                 presentation = presentation,
                 attached_pages = attachment$pages,
                 stylesheet = attachment$stylesheet,
                 script = attachment$script,
                 external_reference_files = portability$external_reference_files,
                 html_files_with_absolute_references = absolute_reference_files,
                 symlinks = portability$symlinks,
                 manifest = file.path(output_dir, "explore_export_manifest.tsv")))
}

lisa_explore_zip_bundle <- function(output_dir) {
  zip_program <- Sys.getenv("R_ZIPCMD", "zip")
  if (!nzchar(zip_program)) zip_program <- "zip"
  if (!nzchar(Sys.which(zip_program))) {
    stop("LISA-EXPLORE-033 no zip program is available; export the bundle as a ",
         "directory instead.", call. = FALSE)
  }
  archive <- paste0(output_dir, ".zip")
  if (lisa_path_entry_exists(archive)) {
    stop("LISA-EXPLORE-034 the archive already exists: ", archive, call. = FALSE)
  }
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(dirname(output_dir))
  # `zip =` must be passed explicitly: `utils::zip()`'s own default is the same
  # `Sys.getenv("R_ZIPCMD", "zip")` expression, so it inherits the empty-string
  # problem and fails with "'zip' must be a non-empty character string".
  status <- utils::zip(archive, basename(output_dir), flags = "-r9Xq",
                       zip = zip_program)
  if (!identical(as.integer(status), 0L)) {
    stop("LISA-EXPLORE-035 the zip program failed with status ", status,
         call. = FALSE)
  }
  archive
}
