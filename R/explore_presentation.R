# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# Presentation of already validated exploration artifacts.
#
# These helpers intentionally do not open a workspace, plan, submit, poll, or
# export. Given a presentation table and existing HTML they return transformed
# HTML only. This lets the Shiny shell and the static exporter share one safe
# attachment representation without turning browsing/export into rendering.
# ---------------------------------------------------------------------------

lisa_explore_presentation_columns <- function() {
  c("key", "category_id", "product", "attach_route", "attach_anchor",
    "attach_query", "png", "pdf", "source_data", "recipe")
}

# --- where a block goes on the page ---------------------------------
#
# One placement rule, stated once, honoured by both surfaces that render a block:
# the static exporter below (over HTML text) and `inst/shiny_assets/explore.js`
# (over the live iframe DOM). Writing it twice in two languages is the reason the
# exported bundle and the running app have to agree, so the rule is deliberately
# simple enough to express identically in both:
#
#   insert before the first `section.panel` inside `main`; if there is none,
#   append at the end of `main`.
#
# On an evidence sheet that lands the block directly under the category heading
# and its counts, above the analytical panels -- inside the category the reader
# opened, not at the foot of a collection. On the navigator fallback, which has
# no `.panel`, it degrades to the previous end-of-main behaviour.
#
# `lisa_explore_presentation_scoped()` decides whether the page routes by
# category at all. An evidence sheet does: it renders one category at a time and
# mirrors the choice into `?category=` and `#category-<id>`. Blocks there are
# emitted `hidden` and revealed by the active category, which is what stops
# category A from displaying B's artifact as its own. The navigator shows every
# category at once, so blocks there are not hidden.


lisa_explore_presentation_panel_pattern <- function() {
  "(?i)<section\\b[^>]*\\bclass\\s*=\\s*([\"'])[^\"']*\\bpanel\\b[^\"']*\\1[^>]*>"
}

lisa_explore_presentation_scoped <- function(html) {
  grepl('id="category-info"', html, fixed = TRUE)
}

# The character span of one owner's analysis block: from its
# `data-lisa-nav-context` attribute to the start of the next one, or to the end
# of the document. Bounding it this way is what stops a collection attribute
# belonging to a DIFFERENT owner -- every owner has a `GOCC` block -- from being
# matched, which would attach one contrast's figure under another's.
lisa_explore_presentation_context_span <- function(html, context) {
  marker <- paste0('data-lisa-nav-context="',
                   lisa_explore_presentation_escape(context, TRUE), '"')
  at <- regexpr(marker, html, fixed = TRUE)
  if (at[[1L]] < 0L) return(NULL)
  from <- at[[1L]]
  rest <- substr(html, from + attr(at, "match.length")[[1L]], nchar(html))
  nxt <- regexpr("data-lisa-nav-context=", rest, fixed = TRUE)
  to <- if (nxt[[1L]] < 0L) nchar(html) else
    from + attr(at, "match.length")[[1L]] + nxt[[1L]] - 2L
  list(from = from, to = to)
}

lisa_explore_presentation_collection_insert_at <- function(html, context,
                                                            collection) {
  span <- lisa_explore_presentation_context_span(html, context)
  if (is.null(span)) return(NA_integer_)
  body <- substr(html, span$from, span$to)
  marker <- paste0('data-lisa-nav-collection="',
                   lisa_explore_presentation_escape(collection, TRUE), '"')
  at <- regexpr(marker, body, fixed = TRUE)
  if (at[[1L]] < 0L) return(NA_integer_)
  inner <- regexpr('(?i)<div\\b[^>]*\\bclass\\s*=\\s*(["\'])[^"\']*\\bcollection-inner\\b[^"\']*\\1[^>]*>',
                   substr(body, at[[1L]], nchar(body)), perl = TRUE)
  if (inner[[1L]] < 0L) return(NA_integer_)
  after_inner <- at[[1L]] + inner[[1L]] - 1L + attr(inner, "match.length")[[1L]]
  tail <- substr(body, after_inner, nchar(body))
  panel <- regexpr(lisa_explore_presentation_panel_pattern(), tail, perl = TRUE)
  offset <- if (panel[[1L]] < 0L) after_inner else after_inner + panel[[1L]] - 1L
  as.integer(span$from + offset - 1L)
}

# Absolute character position at which a block must be inserted, or NA when the
# page has no usable `main` element.
lisa_explore_presentation_insert_at <- function(html) {
  closes <- gregexpr("(?i)</main\\s*>", html, perl = TRUE)[[1L]]
  if (length(closes) != 1L || closes[[1L]] < 0L) return(NA_integer_)
  close_at <- closes[[1L]]
  opens <- gregexpr("(?i)<main\\b[^>]*>", html, perl = TRUE)[[1L]]
  if (length(opens) != 1L || opens[[1L]] < 0L) return(as.integer(close_at))
  body_from <- opens[[1L]] + attr(opens, "match.length")[[1L]]
  if (body_from >= close_at) return(as.integer(close_at))
  body <- substr(html, body_from, close_at - 1L)
  panel <- regexpr(lisa_explore_presentation_panel_pattern(), body, perl = TRUE)
  if (panel[[1L]] < 0L) return(as.integer(close_at))
  as.integer(body_from + panel[[1L]] - 1L)
}

# Add the native section to the page's EXISTING "On this page" index rather than
# creating a parallel navigation, and keep its stated section count honest.
lisa_explore_presentation_add_toc <- function(html, anchor, label) {
  marker <- paste0('href="#', lisa_explore_presentation_escape(anchor, TRUE), '"')
  if (grepl(marker, html, fixed = TRUE)) return(html)
  toc <- regexpr('(?i)<details\\b[^>]*\\bclass\\s*=\\s*(["\'])[^"\']*\\btoc\\b[^"\']*\\1[^>]*>',
                 html, perl = TRUE)
  if (toc[[1L]] < 0L) return(html)
  from <- toc[[1L]] + attr(toc, "match.length")[[1L]]
  close <- regexpr("(?i)</details\\s*>", substr(html, from, nchar(html)), perl = TRUE)
  if (close[[1L]] < 0L) return(html)
  block_end <- from + close[[1L]] - 1L
  div_close <- gregexpr("(?i)</div\\s*>", substr(html, from, block_end - 1L),
                        perl = TRUE)[[1L]]
  if (div_close[[1L]] < 0L) return(html)
  at <- from + utils::tail(div_close, 1L)[[1L]] - 1L
  link <- paste0('<a href="#', lisa_explore_presentation_escape(anchor, TRUE),
                 '" data-lisa-explore-toc="true">',
                 lisa_explore_presentation_escape(label), "</a>")
  html <- paste0(substr(html, 1L, at - 1L), link, substr(html, at, nchar(html)))
  # The summary advertises a count; a section added without updating it makes
  # the page contradict itself. The count is read from the matched text rather
  # than by anchoring a pattern at the start of the document, because the page
  # is one long string in which `.` does not cross the line breaks.
  pattern <- "(?i)<summary>On this page \\(([0-9]+)(\\s+sections?\\))"
  found <- regexpr(pattern, html, perl = TRUE)
  if (found[[1L]] < 0L) return(html)
  matched <- regmatches(html, found)
  digits <- regmatches(matched, regexpr("[0-9]+", matched, perl = TRUE))
  if (!length(digits)) return(html)
  count <- suppressWarnings(as.integer(digits[[1L]]))
  if (is.na(count)) return(html)
  replacement <- sub("[0-9]+", as.character(count + 1L), matched, perl = TRUE)
  paste0(substr(html, 1L, found[[1L]] - 1L), replacement,
         substr(html, found[[1L]] + attr(found, "match.length")[[1L]],
                nchar(html)))
}

lisa_explore_presentation_escape <- function(value, attribute = FALSE) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  if (isTRUE(attribute)) {
    value <- gsub('"', "&quot;", value, fixed = TRUE)
    value <- gsub("'", "&#39;", value, fixed = TRUE)
  }
  value
}

lisa_explore_presentation_path <- function(path, field, allow_empty = FALSE) {
  path <- as.character(path)
  if (length(path) != 1L || is.na(path) || (!allow_empty && !nzchar(path))) {
    stop("LISA-EXPLORE-040 ", field, " must be one contained relative path.",
         call. = FALSE)
  }
  if (!nzchar(path) && isTRUE(allow_empty)) return("")
  path <- gsub("\\\\", "/", path)
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^/|^~|(^|/)\\.\\.(/|$)|[\r\n]", path,
            perl = TRUE)) {
    stop("LISA-EXPLORE-041 ", field,
         " must be a contained relative path, not a URL or traversal.", call. = FALSE)
  }
  pieces <- strsplit(path, "/", fixed = TRUE)[[1L]]
  if (!length(pieces) || any(!nzchar(pieces)) || any(pieces == ".")) {
    stop("LISA-EXPLORE-041 ", field, " must be a normal contained relative path.",
         call. = FALSE)
  }
  path
}

# A query route, not a path: `category=<encoded id>`. It is recorded so a reader
# can be sent straight to the category a figure belongs to, and so it is
# reachable from the exported bundle by URL. Anything that is not a plain
# `name=value` pair is refused rather than pasted into a link.
lisa_explore_presentation_query <- function(query) {
  query <- as.character(query)
  if (length(query) != 1L || is.na(query)) return("")
  if (!nzchar(query)) return("")
  if (!grepl("^[A-Za-z0-9_.-]+=[A-Za-z0-9%_.~-]*$", query)) {
    stop("LISA-EXPLORE-047 attach_query must be one simple name=value route.",
         call. = FALSE)
  }
  query
}

lisa_explore_presentation_validate <- function(presentation) {
  if (!is.data.frame(presentation)) {
    stop("LISA-EXPLORE-042 presentation must be a data frame.", call. = FALSE)
  }
  needed <- lisa_explore_presentation_columns()
  missing <- setdiff(needed, names(presentation))
  if (length(missing)) {
    stop("LISA-EXPLORE-043 presentation is missing columns: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  optional <- intersect(c("entity", "variant", "attach_scope", "attach_context",
                          "attach_collection"), names(presentation))
  presentation <- presentation[, c(needed, optional), drop = FALSE]
  for (index in seq_len(nrow(presentation))) {
    presentation$attach_route[[index]] <- lisa_explore_presentation_path(
      presentation$attach_route[[index]], "attach_route")
    presentation$attach_query[[index]] <- lisa_explore_presentation_query(
      presentation$attach_query[[index]])
    for (field in c("png", "pdf", "source_data", "recipe")) {
      presentation[[field]][[index]] <- lisa_explore_presentation_path(
        presentation[[field]][[index]], field, allow_empty = TRUE)
    }
  }
  presentation
}

lisa_explore_presentation_relative_href <- function(target, page_route) {
  target <- lisa_explore_presentation_path(target, "target")
  page_route <- lisa_explore_presentation_path(page_route, "page_route")
  page_dir <- dirname(page_route)
  from <- if (identical(page_dir, ".")) character() else
    strsplit(page_dir, "/", fixed = TRUE)[[1L]]
  to <- strsplit(target, "/", fixed = TRUE)[[1L]]
  common <- 0L
  while (common < min(length(from), length(to)) &&
         identical(from[[common + 1L]], to[[common + 1L]])) common <- common + 1L
  remainder <- if (common == length(to)) character() else
    to[seq.int(common + 1L, length(to))]
  pieces <- c(rep("..", length(from) - common), remainder)
  paste(pieces, collapse = "/")
}

lisa_explore_presentation_title <- function(product) {
  labels <- c(volcano = "Volcano", gene_cards = "Prioritized gene card",
              heatmap = "Gene heatmap",
              # `kegg` is the legacy GSEA gene-set chart family and keeps its own
              # wording; `kegg_pathway_map` is the native painted diagram. The two
              # must never read as the same product.
              kegg = "KEGG gene sets",
              kegg_pathway_map = "KEGG pathway map",
              contrast_profile = "Contrast profile",
              contrast_heatmap = "Paired gene heatmap",
              de_recurrent_genes = "Recurrent genes",
              contrast_gene_card = "Contrast gene card",
              contrast_paired_heatmap = "Paired gene heatmaps",
              contrast_gene_category_network = "Gene-category network",
              contrast_kegg_map = "KEGG maps / painted pathways")
  label <- labels[[as.character(product)]]
  if (is.null(label)) "Extended figure" else label
}

# The presentation scope of one row, asked of the central product-scope helper
# so the page, the live shell and the exporter cannot disagree about whether a
# block belongs to a category or to a collection.
lisa_explore_presentation_row_scope <- function(row) {
  declared <- as.character(row$attach_scope %||% "")
  if (length(declared) == 1L && !is.na(declared) && nzchar(declared)) return(declared)
  lisa_explore_product_scope(row$product)
}

# The element id a native collection block carries, so the page's existing
# section index can link to it. It is derived from the artifact key, which makes
# it stable across renders and unique per figure.
lisa_explore_presentation_anchor_id <- function(key) {
  paste0("lisa-explore-", as.character(key))
}

# A human label for the exact selector a request carried, appended to the title
# so two variants of one category's figure are distinguishable on the page.
lisa_explore_presentation_selector <- function(row) {
  entity <- as.character(row$entity %||% "")
  variant <- as.character(row$variant %||% "")
  parts <- c(if (length(entity) && !is.na(entity) && nzchar(entity)) entity,
             if (length(variant) && !is.na(variant) && nzchar(variant))
               paste0(variant, " scale"))
  if (!length(parts)) "" else paste(parts, collapse = ", ")
}

lisa_explore_product_file_patterns <- function(product) {
  switch(as.character(product),
    kegg_pathway_map = list(
      png = "_painted[.]png$", pdf = "_painted[.]pdf$",
      source_data = "_painted_source[.]tsv$", recipe = "_painted_recipe[.]R$"),
    contrast_kegg_map = list(
      png = "_contrast_painted[.]png$", pdf = "_contrast_painted[.]pdf$",
      source_data = "_contrast_painted_source[.]tsv$",
      recipe = "_contrast_painted_recipe[.]R$"),
    heatmap = list(
      png = "[.]png$", pdf = "[.]pdf$",
      source_data = "_matrix[.]tsv$", recipe = "_recipe[.]R$"),
    list(png = "[.]png$", pdf = "[.]pdf$",
         source_data = "_source[.]tsv$", recipe = "_recipe[.]R$"))
}

lisa_explore_pick_product_file <- function(files, product, role) {
  pattern <- lisa_explore_product_file_patterns(product)[[role]]
  if (is.null(pattern)) return("")
  hit <- sort(files[grepl(pattern, files)])
  if (!length(hit)) "" else hit[[1L]]
}

lisa_explore_presentation_block <- function(row, page_route, href = NULL,
                                            scoped = FALSE) {
  key <- as.character(row$key)
  if (length(key) != 1L || is.na(key) || !nzchar(key)) {
    stop("LISA-EXPLORE-044 presentation key must be non-empty.", call. = FALSE)
  }
  anchor <- as.character(row$attach_anchor)
  category <- as.character(row$category_id)
  collection_scope <- identical(lisa_explore_presentation_row_scope(row),
                                "collection")
  title <- lisa_explore_presentation_title(row$product)
  # The exact selector is part of the figure's name on the page. Without it a
  # z-score heatmap and a raw-scale heatmap of one category, or two pathway maps
  # reachable from one category, would present as the same unnamed "Gene heatmap"
  # or "KEGG pathway map".
  selector <- lisa_explore_presentation_selector(row)
  if (nzchar(selector)) title <- paste0(title, " (", selector, ")")
  esc <- lisa_explore_presentation_escape
  href_for <- if (is.null(href)) function(path)
    lisa_explore_presentation_relative_href(path, page_route) else href
  link <- function(field, label) {
    value <- as.character(row[[field]])
    if (!nzchar(value)) return("")
    paste0('<a download href="', esc(href_for(value), TRUE),
           '">', esc(label), "</a>")
  }
  downloads <- c(link("png", "Download PNG"), link("pdf", "Download PDF"),
                 link("source_data", "Source data"), link("recipe", "R script"))
  downloads <- downloads[nzchar(downloads)]
  # A collection-wide figure is described by its owner and collection, never by
  # a category: calling it "Recurrent genes for <category>" would assert exactly
  # the category crop H3_SCOPE forbids.
  subject <- if (collection_scope)
    paste(as.character(row$attach_context %||% ""),
          as.character(row$attach_collection %||% row$collection %||% "")) else
    category
  image <- if (nzchar(as.character(row$png))) paste0(
    '<img class="lisa-explore-attachment-image" src="',
    esc(href_for(row$png), TRUE),
    '" alt="', esc(paste(title, "for", trimws(subject)), TRUE), '">') else ""
  # Stated where the download actually is, because it is easy to assume the
  # opposite: the figure, its PDF and its source table are plain files and open
  # anywhere, while the recipe is an R script and needs an R with lisaR and its
  # dependencies installed. That is the accepted base's design for recipes,
  # inherited unchanged here; it is not a defect of this figure.
  note <- if (nzchar(as.character(row$recipe)))
    paste0('<p class="lisa-explore-note">The PNG, PDF and source table open ',
           'without R. The R script re-draws this figure and needs an R ',
           'installation with lisaR and its dependencies.</p>') else ""
  # A collection-wide block is never hidden by the active category, and carries
  # the id the page's section index links to. A category block keeps the
  # accepted behaviour exactly.
  hidden <- isTRUE(scoped) && !collection_scope
  paste0('<section class="lisa-explore-attachment" data-lisa-explore-attachment="true" ',
         if (collection_scope) paste0('id="',
           esc(lisa_explore_presentation_anchor_id(key), TRUE), '" ') else "",
         'data-lisa-explore-key="', esc(key, TRUE), '" data-lisa-explore-category="',
         esc(category, TRUE), '" data-lisa-explore-scope="',
         esc(if (collection_scope) "collection" else "category", TRUE),
         '" data-lisa-attach-anchor="',
         esc(anchor, TRUE), '"', if (hidden) " hidden" else "",
         '><h2>', esc(title), if (nzchar(trimws(subject)))
           paste0(' \u00b7 ', esc(trimws(subject))) else "",
         '</h2><p>', esc(if (collection_scope)
           "Available extended figure for this whole collection."
           else "Available extended figure for this category."), '</p>', image,
         '<nav class="lisa-explore-downloads" aria-label="', esc(paste(title, "downloads"), TRUE),
         '">', paste(downloads, collapse = " "), "</nav>", note, "</section>")
}

#' Attach available figures to their existing report route
#'
#' A pure HTML transformation shared by the interactive shell and static export.
#' It adds only already validated artifacts from `presentation`; it never plans or
#' renders a figure. `route` is the report-relative page being assembled.
#'
#' @param html Existing page HTML.
#' @param presentation A validated `report_presentation.tsv` data frame.
#' @param route Report-relative route of `html`.
#' @return HTML with the route's available figure blocks appended inside `main`.
lisa_explore_presentation_attach_html <- function(html, presentation, route) {
  if (!is.character(html) || !length(html) || anyNA(html)) {
    stop("LISA-EXPLORE-045 html must be non-missing text.", call. = FALSE)
  }
  route <- lisa_explore_presentation_path(route, "route")
  presentation <- lisa_explore_presentation_validate(presentation)
  rows <- presentation[presentation$attach_route == route, , drop = FALSE]
  if (!nrow(rows)) return(paste(html, collapse = "\n"))
  existing <- paste(html, collapse = "\n")
  duplicate <- vapply(seq_len(nrow(rows)), function(index) {
    marker <- paste0('data-lisa-explore-key="',
                     lisa_explore_presentation_escape(rows$key[[index]], TRUE),
                     '" data-lisa-explore-category="',
                     lisa_explore_presentation_escape(rows$category_id[[index]], TRUE),
                     '"')
    grepl(marker, existing, fixed = TRUE)
  }, logical(1L))
  rows <- rows[!duplicate, , drop = FALSE]
  if (!nrow(rows)) return(existing)
  scoped <- lisa_explore_presentation_scoped(existing)

  placements <- lapply(seq_len(nrow(rows)), function(index) {
    row <- rows[index, , drop = FALSE]
    collection_scope <- identical(lisa_explore_presentation_row_scope(row),
                                  "collection")
    position <- if (collection_scope) {
      lisa_explore_presentation_collection_insert_at(existing,
        as.character(row$attach_context %||% ""),
        as.character(row$attach_collection %||% row$collection %||% ""))
    } else lisa_explore_presentation_insert_at(existing)
    if (!collection_scope && is.na(position)) {
      stop("LISA-EXPLORE-046 report page must contain exactly one main element.",
           call. = FALSE)
    }
    # A native section whose owner/collection context is genuinely absent from
    # this page is reported by staying in `report_presentation.tsv`; it is never
    # dropped into the first panel on the page instead, which would attach a
    # global figure to somebody else's collection.
    if (is.na(position)) return(NULL)
    list(position = as.integer(position),
         html = lisa_explore_presentation_block(row, route, scoped = scoped),
         collection_scope = collection_scope,
         anchor = lisa_explore_presentation_anchor_id(row$key),
         label = paste(lisa_explore_presentation_title(row$product),
                       trimws(paste(as.character(row$attach_context %||% ""),
                                    as.character(row$attach_collection %||%
                                                   row$collection %||% "")))))
  })
  placements <- Filter(Negate(is.null), placements)
  if (!length(placements)) return(existing)
  order_by <- order(vapply(placements, `[[`, integer(1), "position"),
                    decreasing = TRUE)
  for (placement in placements[order_by]) {
    existing <- paste0(substr(existing, 1L, placement$position - 1L),
                       placement$html,
                       substr(existing, placement$position, nchar(existing)))
  }
  # The section index is updated after every block is in place, so each anchor
  # it names is really on the page.
  for (placement in placements) {
    if (!isTRUE(placement$collection_scope)) next
    existing <- lisa_explore_presentation_add_toc(existing, placement$anchor,
                                                  trimws(placement$label))
  }
  existing
}
