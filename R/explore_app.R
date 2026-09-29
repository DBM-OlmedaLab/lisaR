# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# The native report remains the presentation. This shell injects an exact
# per-category control into that report; it is deliberately not a dashboard.

lisa_explore_shiny_scalar <- function(value, field) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("LISA-EXPLORE-050 ", field, " must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

lisa_explore_shiny_resource_path <- function(prefix, path) {
  path <- lisa_explore_presentation_path(path, "workspace asset")
  paste0("/", prefix, "/", path)
}

lisa_explore_shiny_presentation <- function(ws) {
  artifacts <- lisa_explore_artifacts(ws)
  empty <- data.frame(key = character(), category_id = character(), product = character(),
    attach_route = character(), attach_anchor = character(), attach_query = character(),
    png = character(), pdf = character(), source_data = character(),
    recipe = character(), entity = character(), variant = character(),
    attach_scope = character(), attach_context = character(),
    attach_collection = character(),
    stringsAsFactors = FALSE)
  if (!nrow(artifacts)) return(empty)
  entries <- lisa_explore_read_index(ws)
  rows <- lapply(seq_len(nrow(artifacts)), function(index) {
    artifact <- artifacts[index, , drop = FALSE]
    entry <- lisa_explore_index_find(entries, artifact$key)
    if (is.null(entry) || is.null(entry$files)) return(NULL)
    files <- vapply(entry$files, function(file) as.character(file$path), character(1L))
    # Product-aware, not "first file ending in .png": a native pathway map also
    # ships the preserved unpainted base diagram, and publishing that as the
    # figure would be the wrong image with the right extension.
    pick <- function(role) {
      hit <- lisa_explore_pick_product_file(files, artifact$product, role)
      if (!nzchar(hit)) "" else file.path("extensions", artifact$artifact_dir, hit)
    }
    attachment <- lisa_explore_attachment(lisa_explore_request_from_row(artifact),
                                          ws$source_run)
    data.frame(key = artifact$key, category_id = artifact$category_id,
      product = artifact$product, attach_route = attachment$route,
      attach_anchor = attachment$anchor, attach_query = attachment$query,
      png = pick("png"), pdf = pick("pdf"),
      source_data = pick("source_data"), recipe = pick("recipe"),
      entity = lisa_explore_blank(artifact$entity),
      variant = lisa_explore_blank(artifact$variant),
      attach_scope = attachment$scope, attach_context = attachment$context,
      attach_collection = attachment$collection,
      stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) empty else do.call(rbind, rows)
}

lisa_explore_shiny_control_payload <- function(ws, workspace_prefix) {
  status <- lisa_explore_status(ws)
  unavailable <- lisa_explore_shiny_unavailable(ws)
  status <- status[status$product %in% lisa_explore_generatable_products(), ,
                   drop = FALSE]
  if (!nrow(status)) return(unavailable)
  presentation <- lisa_explore_shiny_presentation(ws)
  artifacts <- lisa_explore_artifacts(ws)
  offered <- lapply(seq_len(nrow(status)), function(index) {
    row <- status[index, , drop = FALSE]
    request <- lisa_explore_request_from_row(row)
    # Resolved against the source run so the route is a page that exists, and so
    # it is the category sheet the reader opens rather than the collection
    # navigator.
    attachment <- lisa_explore_attachment(request, ws$source_run)
    request_id <- lisa_request_id(request)
    # Matched on the reuse KEY, not the request id. A native pathway map shared
    # by two categories has one key and two request ids; matching on request id
    # would leave the second category offering to generate a map that is already
    # on disk, and then quietly reusing it -- a button that lies about its work.
    key <- as.character(row$key)
    collection_scope <- identical(lisa_explore_product_scope(row$product),
                                  "collection")
    associations <- lisa_explore_catalogue_associations(status, key)
    available <- presentation[presentation$key %in% key, , drop = FALSE]
    block <- ""
    if (nrow(available) &&
        (collection_scope || as.character(row$category_id) %in% associations)) {
      # Re-attach at THIS category's route, not at the route recorded for
      # whichever category happened to request the map first.
      one <- available[1L, , drop = FALSE]
      one$category_id <- as.character(row$category_id)
      one$attach_anchor <- attachment$anchor
      one$attach_scope <- attachment$scope
      one$attach_context <- attachment$context
      one$attach_collection <- attachment$collection
      block <- lisa_explore_presentation_block(one, attachment$route,
        href = function(path)
          lisa_explore_shiny_resource_path(workspace_prefix, path))
    }
    list(route = attachment$route, query = attachment$query,
      anchor = attachment$anchor, category_id = as.character(row$category_id),
      # Carried explicitly so the live shell places a block the same way the
      # static exporter does, rather than inferring it from the product name.
      scope = attachment$scope, context = attachment$context,
      nav_collection = attachment$collection,
      product = as.character(row$product), state = as.character(row$state),
      reason = lisa_explore_blank(row$reason),
      label = lisa_explore_action_label(
        row$product, row$entity, row$variant,
        if ("kegg_title" %in% names(row)) row$kegg_title else ""),
      entity = lisa_explore_blank(row$entity),
      entity_title = if ("kegg_title" %in% names(row))
        lisa_explore_blank(row$kegg_title) else "",
      entity_rank = if ("rank" %in% names(row))
        lisa_explore_blank(row$rank) else "",
      variant = lisa_explore_blank(row$variant),
      key = key, request_id = request_id, unit_type = as.character(row$unit_type),
      analysis_id = as.character(row$analysis_id), contrast_id = as.character(row$contrast_id),
      collection = as.character(row$collection), block = block)
  })
  c(offered, unavailable)
}

# Export stage writes happen frequently while large trees are copied. They do
# not change any figure control, so they must not invalidate the expensive
# catalog/presentation payload on Shiny's single event loop. Export status is
# read separately by the polling observer below.
lisa_explore_control_fingerprint <- function(ws) {
  job_files <- list.files(file.path(ws$root, "jobs"), pattern = "[.]json$",
                          full.names = TRUE)
  figure_files <- Filter(function(path) {
    job <- lisa_explore_read_json(path)
    !is.null(job) && identical(lisa_explore_job_kind(job), "figure")
  }, job_files)
  files <- c(ws$index_path, unlist(figure_files, use.names = FALSE))
  files <- files[file.exists(files)]
  if (!length(files)) return("")
  info <- file.info(files)
  paste(basename(files), info$size, format(info$mtime, "%Y%m%d%H%M%OS3"),
        collapse = "|")
}

lisa_explore_shiny_unavailable <- function(ws, catalog = NULL) {
  if (is.null(catalog)) catalog <- lisa_explore_catalog(ws)
  if (!nrow(catalog)) return(list())
  products <- setdiff(lisa_explore_generatable_products(), "kegg")
  notices <- list()
  notice <- function(pseudo, request_id) {
    attachment <- lisa_explore_attachment(pseudo, ws$source_run)
    list(route = attachment$route, query = attachment$query,
         anchor = attachment$anchor, scope = attachment$scope,
         context = attachment$context, nav_collection = attachment$collection,
         category_id = lisa_explore_blank(pseudo$category_id)[[1L]],
         product = pseudo$product, state = "not_applicable",
         reason = lisa_explore_unavailable_reason(ws, pseudo),
         label = "", entity = "", variant = "", key = "",
         request_id = request_id,
         unit_type = pseudo$unit_type,
         analysis_id = lisa_explore_blank(pseudo$analysis_id)[[1L]],
         contrast_id = lisa_explore_blank(pseudo$contrast_id)[[1L]],
         collection = as.character(pseudo$collection), block = "")
  }
  for (unit_type in c("single_de", "contrast")) {
    rows <- catalog[catalog$unit_type == unit_type, , drop = FALSE]
    if (!nrow(rows)) next
    owner_column <- if (identical(unit_type, "single_de")) "analysis_id" else
      "contrast_id"
    kinds <- Filter(function(product)
      identical(lisa_explore_product_unit_type(product), unit_type) ||
        (identical(unit_type, "single_de") &&
           is.na(lisa_explore_product_unit_type(product))), products)
    if (!length(kinds)) next
    category_kinds <- Filter(lisa_explore_product_needs_category, kinds)
    collection_kinds <- setdiff(kinds, category_kinds)

    if (length(category_kinds)) {
      scopes <- unique(rows[, c(owner_column, "collection", "category_id"),
                            drop = FALSE])
      scopes <- scopes[nzchar(scopes$category_id), , drop = FALSE]
      for (index in seq_len(nrow(scopes))) {
        scope <- scopes[index, , drop = FALSE]
        owner <- as.character(scope[[owner_column]][[1L]])
        in_scope <- rows[[owner_column]] == owner &
          rows$collection == scope$collection &
          rows$category_id == scope$category_id
        for (product in setdiff(category_kinds,
                                unique(as.character(rows$product[in_scope])))) {
          pseudo <- list(unit_type = unit_type,
            analysis_id = if (identical(unit_type, "single_de")) owner else NA_character_,
            contrast_id = if (identical(unit_type, "contrast")) owner else NA_character_,
            collection = as.character(scope$collection),
            category_id = as.character(scope$category_id),
            product = product, entity = NA_character_, variant = NA_character_)
          notices[[length(notices) + 1L]] <- notice(pseudo,
            paste0("na-", product, "-", owner, "-", scope$collection, "-",
                   scope$category_id))
        }
      }
    }
    if (length(collection_kinds)) {
      scopes <- unique(rows[, c(owner_column, "collection"), drop = FALSE])
      for (index in seq_len(nrow(scopes))) {
        scope <- scopes[index, , drop = FALSE]
        owner <- as.character(scope[[owner_column]][[1L]])
        in_scope <- rows[[owner_column]] == owner &
          rows$collection == scope$collection
        for (product in setdiff(collection_kinds,
                                unique(as.character(rows$product[in_scope])))) {
          pseudo <- list(unit_type = unit_type,
            analysis_id = if (identical(unit_type, "single_de")) owner else NA_character_,
            contrast_id = if (identical(unit_type, "contrast")) owner else NA_character_,
            collection = as.character(scope$collection),
            category_id = NA_character_,
            product = product, entity = NA_character_, variant = NA_character_)
          notices[[length(notices) + 1L]] <- notice(pseudo,
            paste0("na-", product, "-", owner, "-", scope$collection))
        }
      }
    }
  }
  notices
}

lisa_explore_action_label <- function(product, entity = "", variant = "",
                                      entity_title = "") {
  entity <- lisa_explore_blank(entity)[[1L]]
  variant <- lisa_explore_blank(variant)[[1L]]
  entity_title <- lisa_explore_blank(entity_title)[[1L]]
  base <- switch(as.character(product),
    volcano = "Generate volcano",
    gene_cards = "Generate gene card",
    heatmap = "Generate gene heatmap",
    kegg_pathway_map = "Generate KEGG pathway map",
    de_recurrent_genes = "Generate recurrent genes figure",
    contrast_gene_card = "Generate contrast gene card",
    contrast_paired_heatmap = "Generate paired gene heatmap",
    contrast_gene_category_network = "Generate gene-category network",
    contrast_kegg_map = "Generate contrast KEGG map",
    "Generate figure")
  entity_label <- if (nzchar(entity_title) && nzchar(entity)) {
    paste0(entity_title, " (", entity, ")")
  } else if (nzchar(entity)) entity else character()
  detail <- c(entity_label,
              if (nzchar(variant)) paste0(variant, " scale"))
  if (!length(detail)) base else paste0(base, " (", paste(detail, collapse = ", "), ")")
}

lisa_explore_shiny_request <- function(payload) {
  if (!is.list(payload)) stop("LISA-EXPLORE-051 malformed figure action.", call. = FALSE)
  optional <- function(value) if (is.null(value) || is.na(value) || !nzchar(value)) NULL else value
  lisa_figure_request(unit_type = payload$unit_type, analysis_id = optional(payload$analysis_id),
    contrast_id = optional(payload$contrast_id), collection = payload$collection,
    category_id = optional(payload$category_id), product = payload$product,
    entity = optional(payload$entity), variant = optional(payload$variant))
}

#' Explore a completed lisaR run interactively
#'
#' Opens the existing report in a local Shiny shell. Navigation reads the report
#' and engine catalog only; a figure is submitted only after its explicit native
#' per-category product button is pressed. Completed artifacts appear
#' at that category route with their PNG, PDF, source-data and R-script downloads.
#' Saving the expanded report is explicit and never generates missing products.
#'
#' @param run_dir Completed STANDARD run, kept read-only.
#' @param extension_root Optional external exploration workspace.
#' @param launch.browser Whether Shiny should launch a browser.
#' @return Invisibly, the Shiny application result.
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   explore_lisa_run(run_dir)
#' }
explore_lisa_run <- function(run_dir, extension_root = NULL, launch.browser = TRUE) {
  launch.browser <- lisa_explore_shiny_scalar(launch.browser, "launch.browser")
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("LISA-EXPLORE-052 explore_lisa_run() needs optional package 'shiny'. ",
      "Install it in the selected project environment; lisaR never installs optional packages automatically.",
      call. = FALSE)
  }
  ws <- lisa_explore_open(run_dir, extension_root)
  nonce <- substr(lisa_sha256_text(paste(ws$root, Sys.getpid(), Sys.time())), 1L, 12L)
  report_prefix <- paste0("lisa-explore-report-", nonce)
  workspace_prefix <- paste0("lisa-explore-workspace-", nonce)
  assets_prefix <- paste0("lisa-explore-assets-", nonce)
  shiny::addResourcePath(report_prefix, ws$source_run)
  shiny::addResourcePath(workspace_prefix, ws$root)
  asset_dir <- system.file("shiny_assets", package = "lisaR")
  if (!nzchar(asset_dir) || !dir.exists(asset_dir)) {
    stop("LISA-EXPLORE-053 installed shiny assets are missing.", call. = FALSE)
  }
  shiny::addResourcePath(assets_prefix, asset_dir)

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$link(id = "lisa-explore-stylesheet", rel = "stylesheet",
                       href = paste0("/", assets_prefix, "/explore.css")),
      shiny::tags$script(src = paste0("/", assets_prefix, "/explore.js"))),
    shiny::tags$div(class = "lisa-explore-toolbar",
      shiny::tags$div(shiny::tags$strong("Explore this report"),
        shiny::tags$span("Controls appear in the category evidence sheet you open."),
        # Engine-level failures (an index lock held elsewhere, a worker that
        # could not start) are not about one category, so they are reported
        # here rather than attached to an arbitrary one. Per-request refusals
        # are shown in the category the reader clicked.
        shiny::tags$span(class = "lisa-explore-status", role = "alert",
                         shiny::textOutput("explore_engine_status", inline = TRUE))),
      shiny::textInput("explore_export_path", "Save expanded report to", value = ""),
      shiny::actionButton("explore_export", "Save expanded report"),
      shiny::actionButton("explore_export_cancel", "Cancel export"),
      # Must be a real Shiny output binding. A bare <output id> element carries
      # no `shiny-text-output` class, so `renderText()` would never reach it and
      # the export result -- including its error message -- would be invisible.
      shiny::tags$span(class = "lisa-explore-status",
                       shiny::textOutput("explore_export_status", inline = TRUE))),
    shiny::tags$iframe(id = "lisa-explore-report", class = "lisa-explore-report",
      src = paste0("/", report_prefix, "/report_index.html"), title = "lisaR report"))

  shiny::runApp(shiny::shinyApp(
    ui = ui, server = lisa_explore_shiny_server(ws, workspace_prefix)),
    launch.browser = launch.browser)
}

# The server, built separately from `explore_lisa_run()` so its behaviour can be
# driven headlessly by `shiny::testServer()`. The paths under test are the two
# observers; asserting
# that a refused submission surfaces its message and leaves the session alive
# requires a running reactive context, not a source-text check.
lisa_explore_export_feedback_record <- function(ws, kind, message,
                                                followed_job_id = NULL) {
  record <- list(kind = kind, message = message)
  if (is.null(followed_job_id)) return(record)
  followed_job <- lisa_explore_read_job(ws, followed_job_id)
  if (is.null(followed_job)) return(record)
  record$followed_job_id <- as.character(followed_job_id)
  record$followed_state <- as.character(
    lisa_explore_reconcile_job(ws, followed_job)$state)
  record
}

lisa_explore_shiny_server <- function(ws, workspace_prefix) {
  function(input, output, session) {
    # Rebuilding the whole payload every second was too expensive to do blindly:
    # it walks the catalog, mints a reuse key per row and reads the index, which
    # is exactly the path finding 6 identified as the hot one. The payload is a
    # pure function of the index file and the job records, so we recompute only
    # when their fingerprint changes and otherwise resend the cached value.
    cached <- shiny::reactiveValues(fingerprint = NULL, payload = list())
    # Both observers below call into the engine, and the engine can
    # legitimately raise: LISA-EXPLORE-015 when another session holds the index
    # lock for 30 s, a launch failure, LISA-EXPLORE-021/022 for a request this
    # milestone will not accept, LISA-EXPLORE-036 when a worker is already busy,
    # LISA-EXPLORE-051 on a malformed payload. Shiny closes the session on an
    # unhandled observer error, so an unwrapped failure here did not merely lose
    # a message -- it ended the app, and any queued or completed work became
    # unreachable until the user restarted R.
    #
    # Both are wrapped. The polling observer stays alive after a transient
    # failure and reports it; a refused submission is reported against the
    # request the reader clicked, in the category they are looking at, and the
    # button becomes usable again. Neither path reports success for work that
    # was rejected.
    engine_status <- shiny::reactiveVal("")
    consecutive_lock_failures <- shiny::reactiveVal(0L)
    feedback <- shiny::reactiveVal(list())
    previous_export <- lisa_explore_latest_export_job(ws)
    export_job_id <- shiny::reactiveVal(
      if (is.null(previous_export)) NULL else as.character(previous_export$job_id))
    export_status <- shiny::reactiveVal(
      lisa_explore_export_job_status(ws, previous_export))
    # User-action feedback is intentionally distinct from polled job status.
    # A refused request stays visible until a later successful action; a running
    # cancellation stays pinned only until that job reaches a terminal state.
    export_feedback <- shiny::reactiveVal(NULL)
    output$explore_export_status <- shiny::renderText(export_status())
    shiny::observe({
      # Scheduled first, so the timer survives whatever the body does.
      shiny::invalidateLater(1000, session)
      result <- tryCatch({
        # Hand the free worker slot to anything queued behind a finished figure.
        # This only ever starts a figure some earlier click already submitted;
        # opening or reopening the app never creates a request of its own.
        # Polling is read-mostly and runs on Shiny's single event loop. If a
        # different session owns the short write lock, report that state now
        # and try again next tick instead of blocking all normal interaction for
        # the engine's 30-second batch timeout.
        lisa_explore_pump(ws, lock_timeout_seconds = 0)
        fingerprint <- lisa_explore_control_fingerprint(ws)
        if (!identical(fingerprint, shiny::isolate(cached$fingerprint))) {
          cached$fingerprint <- fingerprint
          cached$payload <- lisa_explore_shiny_control_payload(ws, workspace_prefix)
        }
        NULL
      }, error = function(error) error)
      if (inherits(result, "error")) {
        # Keep the last good payload on screen: a transient lock or read failure
        # must not blank out figures that are perfectly valid.
        message <- conditionMessage(result)
        if (startsWith(message, "LISA-EXPLORE-015")) {
          failures <- shiny::isolate(consecutive_lock_failures()) + 1L
          consecutive_lock_failures(failures)
          # A healthy index owner holds this lock for milliseconds. Do not flash
          # that ordinary transition as a toolbar fault; the same message becomes
          # actionable after it persists across two one-second observations.
          if (failures >= 2L) engine_status(message)
        } else {
          consecutive_lock_failures(0L)
          engine_status(message)
        }
        # Force a recompute on the next tick rather than trusting a fingerprint
        # that may have been recorded before the failure.
        cached$fingerprint <- NULL
      } else {
        consecutive_lock_failures(0L)
        if (nzchar(shiny::isolate(engine_status()))) engine_status("")
      }
      session$sendCustomMessage("lisa-explore-controls", list(
        rows = shiny::isolate(cached$payload),
        feedback = shiny::isolate(feedback())))
      active_export_id <- shiny::isolate(export_job_id())
      if (!is.null(active_export_id)) {
        export_job <- lisa_explore_read_job(ws, active_export_id)
        if (!is.null(export_job)) {
          action_feedback <- shiny::isolate(export_feedback())
          reconciled_state <- as.character(
            lisa_explore_reconcile_job(ws, export_job)$state)
          cancel_finished <- !is.null(action_feedback) &&
            identical(action_feedback$kind, "cancel_pending") &&
            !(reconciled_state %in% c("queued", "running"))
          followed_state_changed <- !is.null(action_feedback) &&
            !identical(action_feedback$kind, "cancel_pending") &&
            !is.null(action_feedback$followed_job_id) &&
            identical(as.character(action_feedback$followed_job_id),
                      as.character(active_export_id)) &&
            !identical(as.character(action_feedback$followed_state),
                       reconciled_state)
          if (cancel_finished || followed_state_changed) {
            export_feedback(NULL)
            action_feedback <- NULL
          }
          if (is.null(action_feedback)) {
            export_status(lisa_explore_export_job_status(ws, export_job))
          }
        }
      }
    })
    output$explore_engine_status <- shiny::renderText(engine_status())
    shiny::observeEvent(input$explore_generate, {
      payload <- input$explore_generate
      request <- tryCatch(lisa_explore_shiny_request(payload),
                          error = function(error) error)
      request_id <- if (inherits(request, "error")) NA_character_ else
        lisa_request_id(request)
      outcome <- if (inherits(request, "error")) request else
        tryCatch({
          lisa_explore_submit(ws, request, background = TRUE)
          NULL
        }, error = function(error) error)
      if (inherits(outcome, "error")) {
        notes <- shiny::isolate(feedback())
        id <- if (is.na(request_id)) "request" else request_id
        notes[[id]] <- list(ok = FALSE, message = conditionMessage(outcome))
        feedback(notes)
      } else {
        notes <- shiny::isolate(feedback())
        notes[[request_id]] <- NULL
        feedback(notes)
      }
      # The payload is derived from the index and job files, which this
      # submission has just changed; drop the fingerprint so the next tick
      # reports the new state instead of the pre-click one.
      cached$fingerprint <- NULL
    })
    shiny::observeEvent(input$explore_export, {
      destination <- trimws(as.character(input$explore_export_path))
      if (!nzchar(destination)) {
        message <- "Choose a new destination before saving."
        export_feedback(lisa_explore_export_feedback_record(
          ws, "error", message, shiny::isolate(export_job_id())))
        export_status(message)
      } else {
        result <- tryCatch(
          lisa_explore_submit_export(ws, destination, background = TRUE),
          error = function(error) error)
        if (inherits(result, "error")) {
          message <- conditionMessage(result)
          export_feedback(lisa_explore_export_feedback_record(
            ws, "error", message, shiny::isolate(export_job_id())))
          export_status(message)
        } else {
          export_feedback(NULL)
          export_job_id(result$job_id)
          export_job <- lisa_explore_read_job(ws, result$job_id)
          export_status(lisa_explore_export_job_status(ws, export_job))
        }
      }
    })
    shiny::observeEvent(input$explore_export_cancel, {
      job_id <- shiny::isolate(export_job_id())
      if (is.null(job_id)) {
        message <- "No export is queued or running."
        export_feedback(lisa_explore_export_feedback_record(
          ws, "error", message, job_id))
        export_status(message)
      } else {
        outcome <- tryCatch(lisa_explore_cancel(ws, job_id),
                            error = function(error) error)
        if (inherits(outcome, "error")) {
          message <- conditionMessage(outcome)
          export_feedback(lisa_explore_export_feedback_record(
            ws, "error", message, job_id))
          export_status(message)
        } else {
          kind <- if (identical(as.character(outcome$state), "running"))
            "cancel_pending" else "action"
          export_feedback(lisa_explore_export_feedback_record(
            ws, kind, outcome$message, job_id))
          export_status(outcome$message)
        }
      }
    })
  }
}
