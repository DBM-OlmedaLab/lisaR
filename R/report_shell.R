# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Shared, portable page chrome.  This module changes navigation and display
# context only: exact scientific identifiers remain owned by the caller.
# It can be sourced by installed report scripts without extra dependencies.

.lisa_shell_escape <- function(x) {
  x <- as.character(x)
  for (pair in list(c("&", "&amp;"), c("<", "&lt;"), c(">", "&gt;"),
      c('"', "&quot;"), c("'", "&#39;"))) {
    x <- gsub(pair[[1L]], pair[[2L]], x, fixed = TRUE)
  }
  x
}

.lisa_shell_text <- function(x, name, empty = TRUE) {
  if (is.null(x)) x <- ""
  if (length(x) != 1L || is.na(x) || !is.atomic(x))
    stop(name, " must be one non-missing text value.", call. = FALSE)
  x <- as.character(x)
  if (!empty && !nzchar(x)) stop(name, " must not be empty.", call. = FALSE)
  x
}

# Routes are deliberately explicit.  Do not guess a home page that may not
# exist when a viewer is opened independently of an assembled report.
.lisa_shell_href <- function(x, name) {
  x <- .lisa_shell_text(x, name, empty = FALSE)
  if (grepl("^[[:space:]]|[[:space:]]$|[[:cntrl:]]", x) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^/|\\\\", x))
    stop(name, " must be a relative offline URL, optionally with query/fragment.", call. = FALSE)
  x
}

# Return links relative to the page output directory. asset_dir refers to the
# installed report_assets directory; callers sourcing a checkout pass its
# inst/report_assets path explicitly.  The existing logo is copied unchanged.
.lisa_copy_report_shell_assets <- function(output_dir, asset_dir = NULL,
    asset_subdir = "lisa-shell") {
  output_dir <- .lisa_shell_text(output_dir, "output_dir", empty = FALSE)
  if (is.null(asset_dir)) asset_dir <- system.file("report_assets", package = "lisaR")
  asset_dir <- .lisa_shell_text(asset_dir, "asset_dir", empty = FALSE)
  asset_subdir <- .lisa_shell_href(asset_subdir, "asset_subdir")
  components <- strsplit(asset_subdir, "/", fixed = TRUE)[[1L]]
  if (any(components %in% c("", ".", "..")) || grepl("[?#%]", asset_subdir))
    stop("asset_subdir must name a local descendant directory.", call. = FALSE)
  filenames <- c(css_href = "lisa_shell.css", js_href = "lisa_shell.js",
    logo_href = "LISA_logo_A1_muted_red_S_automated_annotation_final.svg",
    logo_compact_href = "LISA_logo_C_compact_icon_muted_red_S.svg")
  source <- file.path(asset_dir, unname(filenames))
  if (!all(file.exists(source))) stop("Missing shared shell assets: ",
    paste(unname(filenames)[!file.exists(source)], collapse = ", "), call. = FALSE)
  destination <- file.path(output_dir, asset_subdir)
  if (!dir.exists(destination) && !dir.create(destination, recursive = TRUE))
    stop("Cannot create shared shell asset directory.", call. = FALSE)
  paths <- file.path(destination, unname(filenames))
  for (i in seq_along(source)) {
    # Re-rendering into the same directory should not attempt to copy a file
    # onto itself.  Ordinary generated assets may otherwise be refreshed.
    same <- file.exists(paths[[i]]) && identical(
      normalizePath(source[[i]], winslash = "/"), normalizePath(paths[[i]], winslash = "/"))
    if (!same && !file.copy(source[[i]], paths[[i]], overwrite = TRUE))
      stop("Failed to copy shared shell asset: ", unname(filenames[[i]]), call. = FALSE)
  }
  links <- lapply(unname(filenames), function(x) paste0(asset_subdir, "/", x))
  names(links) <- names(filenames)
  links$files <- paths
  links
}

# routes: named list of overview, analyses, contrasts, genes, methods.  Each
# provided entry is an href string or list(href = ..., label = ...). NULL/""
# entries are omitted. active is the ONE current route key, or NULL.
# context: optional study/analysis/collection/selection/direction/cutoff labels
# and exact_ids = named list/vector. Human labels never replace exact IDs.
.lisa_report_shell <- function(routes = list(), active = NULL, context = list(),
    assets, main_id = "lisa-main", brand_label = "lisaR") {
  labels <- c(overview = "Overview", analyses = "Analyses", contrasts = "Contrasts",
    genes = "Gene search", methods = "Data & methods")
  if (!is.list(routes) || (length(routes) && (is.null(names(routes)) ||
      any(!nzchar(names(routes))) || anyDuplicated(names(routes)) ||
      any(!names(routes) %in% names(labels)))))
    stop("routes must be a uniquely named list of supported route keys.", call. = FALSE)
  if (!is.list(context)) stop("context must be a list.", call. = FALSE)
  main_id <- .lisa_shell_text(main_id, "main_id", empty = FALSE)
  if (!grepl("^[A-Za-z][A-Za-z0-9_.:-]*$", main_id))
    stop("main_id must be a simple HTML identifier.", call. = FALSE)
  brand_label <- .lisa_shell_text(brand_label, "brand_label", empty = FALSE)
  if (!is.list(assets) || !all(c("css_href", "js_href", "logo_href") %in% names(assets)))
    stop("assets must supply css_href, js_href and logo_href.", call. = FALSE)
  hrefs <- lapply(c("css_href", "js_href", "logo_href"), function(key)
    .lisa_shell_href(assets[[key]], paste0("assets$", key)))
  names(hrefs) <- c("css_href", "js_href", "logo_href")
  hrefs$logo_compact_href <- if (is.null(assets$logo_compact_href)) hrefs$logo_href else
    .lisa_shell_href(assets$logo_compact_href, "assets$logo_compact_href")
  resolved <- list()
  for (key in names(labels)) {
    entry <- routes[[key]]
    if (is.null(entry) || identical(entry, "")) next
    if (is.character(entry)) entry <- list(href = entry)
    if (!is.list(entry)) stop("Route ", key, " must be a URL or route object.", call. = FALSE)
    if (is.null(entry$href) || identical(entry$href, "")) next
    resolved[[key]] <- list(href = .lisa_shell_href(entry$href, paste0("routes$", key)),
      label = .lisa_shell_text(if (is.null(entry$label)) labels[[key]] else entry$label,
        paste0("route label ", key), empty = FALSE))
  }
  if (!is.null(active)) {
    active <- .lisa_shell_text(active, "active", empty = FALSE)
    if (!active %in% names(resolved))
      stop("active must identify one supplied, non-empty route.", call. = FALSE)
  }
  esc <- .lisa_shell_escape
  nav <- vapply(names(resolved), function(key) paste0('<a data-lisa-route="', key,
    '" href="', esc(resolved[[key]]$href), '"',
    if (identical(active, key)) ' aria-current="page"' else "", '>',
    esc(resolved[[key]]$label), '</a>'), character(1L))
  brand_inner <- paste0('<picture><source media="(max-width: 600px)" srcset="',
    esc(hrefs$logo_compact_href), '"><img class="lisa-shell-logo" src="', esc(hrefs$logo_href),
    '" alt="LISA - Geneset analysis, automated annotation and biological interpretation" ',
    'width="410" height="128"></picture><span class="lisa-shell-brand-short">',
    esc(brand_label), '</span>')
  brand <- if ("overview" %in% names(resolved)) paste0('<a class="lisa-shell-brand" href="',
    esc(resolved$overview$href), '" aria-label="', esc(paste0(brand_label, " overview")),
    '">', brand_inner, '</a>') else paste0('<span class="lisa-shell-brand">', brand_inner, '</span>')
  context_labels <- c(study = "Study", analysis = "Analysis", collection = "Collection",
    selection = "Viewing", direction = "Positive direction", cutoff = "Cutoff")
  context_values <- vapply(names(context_labels), function(key)
    .lisa_shell_text(context[[key]], paste0("context$", key)), character(1L))
  fields <- vapply(names(context_labels), function(key) paste0(
    '<span class="lisa-shell-context-field" data-lisa-context-field="', key, '"',
    if (!nzchar(context_values[[key]])) ' hidden' else "", '>',
    '<span class="lisa-shell-context-label">', context_labels[[key]], ': </span>',
    '<span data-lisa-context="', key, '">', esc(context_values[[key]]), '</span></span>'), character(1L))
  ids <- context$exact_ids
  if (is.null(ids)) ids <- list()
  if (length(ids) && (is.null(names(ids)) || any(!nzchar(names(ids))) || anyDuplicated(names(ids))))
    stop("context$exact_ids must have unique nonempty names.", call. = FALSE)
  id_values <- vapply(seq_along(ids), function(i)
    .lisa_shell_text(ids[[i]], paste0("context$exact_ids$", names(ids)[[i]])), character(1L))
  exact <- paste(vapply(seq_along(ids), function(i) paste0('<div><dt>', esc(names(ids)[[i]]),
    '</dt><dd><code>', esc(id_values[[i]]), '</code></dd></div>'), character(1L)), collapse = "")
  has_context <- any(nzchar(context_values)) || length(ids) > 0L
  header <- paste0('<a class="lisa-shell-skip" href="#', esc(main_id), '">Skip to content</a>',
    '<header class="lisa-shared-shell" data-lisa-shell>',
    '<div class="lisa-shell-bar">', brand,
    if (length(resolved)) paste0('<button type="button" class="lisa-shell-menu" aria-expanded="false" ',
      'aria-controls="lisa-shell-navigation">Menu</button>',
      '<nav id="lisa-shell-navigation" class="lisa-shell-nav" aria-label="Main navigation">',
      paste(nav, collapse = ""), '</nav>') else "", '</div>',
    '<div class="lisa-shell-explorer" data-lisa-explorer hidden ',
      'role="navigation" aria-label="Explore report results">',
    '<label class="lisa-shell-picker" for="lisa-nav-context">',
      '<span data-lisa-nav-context-label>Analysis / contrast</span>',
      '<select id="lisa-nav-context" data-lisa-nav-select="context"></select></label>',
    '<label class="lisa-shell-picker" for="lisa-nav-collection"><span>Collection</span>',
      '<select id="lisa-nav-collection" data-lisa-nav-select="collection"></select></label>',
    '<label class="lisa-shell-picker" for="lisa-nav-section"><span>Section</span>',
      '<select id="lisa-nav-section" data-lisa-nav-select="section"></select></label>',
    '<span class="lisa-shell-location" data-lisa-nav-location aria-live="polite"></span>',
    '</div>',
    '<div class="lisa-shell-context"', if (!has_context) ' hidden' else "", '>',
    '<div class="lisa-shell-context-fields">', paste(fields, collapse = ""), '</div>',
    '<details class="lisa-shell-identifiers"', if (!length(ids)) ' hidden' else "", '>',
    '<summary>Exact identifiers</summary><dl data-lisa-exact-ids>', exact, '</dl></details>',
    '</div></header>')
  list(head = paste0('<link rel="stylesheet" href="', esc(hrefs$css_href),
      '"><script defer src="', esc(hrefs$js_href), '"></script>'),
    header = header, body_class = "lisa-shell-page", main_id = main_id,
    main_open = paste0('<main id="', esc(main_id), '" class="lisa-shell-main" tabindex="-1">'),
    main_close = '</main>')
}

# Injection is optional; an assembler can instead use shell$head/$header and
# shell$main_open/$main_close directly. Existing page headers are NOT removed:
# the caller owns its content and should drop redundant legacy branding first.
.lisa_inject_report_shell <- function(html, shell, add_main = TRUE) {
  if (!is.character(html) || !length(html) || anyNA(html))
    stop("html must contain non-missing HTML text.", call. = FALSE)
  html <- paste(html, collapse = "\n")
  if (grepl('data-lisa-shell(?:[ >]|=)', html, perl = TRUE))
    stop("This page already contains a shared shell.", call. = FALSE)
  if (!is.list(shell) || !all(c("head", "header", "body_class", "main_id", "main_open", "main_close") %in% names(shell)))
    stop("shell must be produced by .lisa_report_shell().", call. = FALSE)
  if (!is.logical(add_main) || length(add_main) != 1L || is.na(add_main))
    stop("add_main must be TRUE or FALSE.", call. = FALSE)
  insert <- function(text, pattern, value, after = FALSE) {
    matches <- gregexpr(pattern, text, perl = TRUE)[[1L]]
    if (length(matches) != 1L || matches[[1L]] < 0L)
      stop("Expected exactly one HTML element matching ", pattern, ".", call. = FALSE)
    position <- matches[[1L]] + if (after) attr(matches, "match.length")[[1L]] else 0L
    paste0(substr(text, 1L, position - 1L), value, substr(text, position, nchar(text)))
  }
  if (add_main && grepl("(?i)<main\\b", html, perl = TRUE))
    stop("Page already has a main element; supply its ID and use add_main = FALSE.", call. = FALSE)
  # Existing main content must not be nested. Give an unlabelled <main> our
  # target ID, but preserve an existing ID rather than breaking its deep links.
  if (!add_main) {
    targets <- c(paste0('id="', shell$main_id, '"'), paste0("id='", shell$main_id, "'"))
    target_present <- any(vapply(targets, function(target) grepl(target, html, fixed = TRUE), logical(1L)))
    if (!target_present) {
      main <- gregexpr("(?i)<main\\b[^>]*>", html, perl = TRUE)[[1L]]
      if (length(main) != 1L || main[[1L]] < 0L)
        stop("Provide one existing main element or the shell main_id target.", call. = FALSE)
      opening <- substr(html, main[[1L]], main[[1L]] + attr(main, "match.length")[[1L]] - 1L)
      if (grepl("(?i)\\bid\\s*=", opening, perl = TRUE))
        stop("Existing main has another ID; pass that main_id to .lisa_report_shell().", call. = FALSE)
      replacement <- paste0(substr(opening, 1L, nchar(opening) - 1L), ' id="', shell$main_id, '"',
        if (grepl("(?i)\\btabindex\\s*=", opening, perl = TRUE)) "" else ' tabindex="-1"', '>')
      html <- paste0(substr(html, 1L, main[[1L]] - 1L), replacement,
        substr(html, main[[1L]] + attr(main, "match.length")[[1L]], nchar(html)))
    }
  }
  html <- insert(html, "(?i)</head\\s*>", shell$head)
  # Preserve the body attributes and append only our namespaced page class.
  body <- regexpr("(?i)<body\\b[^>]*>", html, perl = TRUE)
  if (body[[1L]] < 0L) stop("Expected an HTML body element.", call. = FALSE)
  opening <- regmatches(html, body)
  cls <- regexpr('(?i)\\bclass\\s*=\\s*(["\x27])[^"\x27]*["\x27]', opening, perl = TRUE)
  if (cls[[1L]] > 0L) {
    old <- regmatches(opening, cls)
    regmatches(opening, cls) <- paste0(substr(old, 1L, nchar(old) - 1L), " ", shell$body_class,
      substr(old, nchar(old), nchar(old)))
  } else opening <- paste0(substr(opening, 1L, nchar(opening) - 1L),
    ' class="', shell$body_class, '">')
  regmatches(html, body) <- opening
  html <- insert(html, "(?i)<body\\b[^>]*>", paste0(shell$header,
    if (add_main) shell$main_open else ""), after = TRUE)
  if (add_main) html <- insert(html, "(?i)</body\\s*>", shell$main_close)
  html
}
