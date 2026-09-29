# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# ---------------------------------------------------------------------------
# Exact figure requests.
#
# The whole point of this file is that a request cannot express more than one
# figure. `lisa_extension_filter()` in report_architecture.R applies `categories`
# and `products` as independent vectors, so a two-category two-product selection
# expands to a Cartesian product. That behaviour is correct for the existing
# "render this whole selection" API and is deliberately left untouched.
#
# The exploration engine needs the opposite guarantee: one request is one figure.
# We obtain it structurally -- every field is validated as a single scalar and a
# vector is a hard error. Unrequested figures are never rendered.
# ---------------------------------------------------------------------------

# Products the shared catalog can enumerate. Requests may name any of them so the
# interface can display their state; only `lisa_explore_h1_products()` can be
# submitted for rendering in this milestone.
lisa_explore_catalog_products <- function() {
  c("member_gene_sets", "gene_cards", "volcano", "heatmap", "kegg",
    "kegg_pathway_map", "contrast_profile", "contrast_heatmap",
    lisa_explore_h3_products())
}

lisa_explore_h1_products <- function() "volcano"

lisa_explore_h2_products <- function() {
  c("volcano", "gene_cards", "heatmap", "kegg_pathway_map")
}

lisa_explore_h3_products <- function() {
  c("de_recurrent_genes", "contrast_gene_card", "contrast_paired_heatmap",
    "contrast_gene_category_network", "contrast_kegg_map")
}

lisa_explore_generatable_products <- function() {
  c(lisa_explore_h2_products(), lisa_explore_h3_products())
}

lisa_explore_collection_products <- function() {
  c("de_recurrent_genes", "contrast_paired_heatmap",
    "contrast_gene_category_network", "contrast_kegg_map")
}

lisa_explore_product_scope <- function(product) {
  if (as.character(product) %in% lisa_explore_collection_products()) "collection"
  else "category"
}

lisa_explore_product_needs_category <- function(product) {
  identical(lisa_explore_product_scope(product), "category")
}

lisa_explore_product_unit_type <- function(product) {
  switch(as.character(product),
    de_recurrent_genes = "single_de",
    contrast_gene_card = ,
    contrast_paired_heatmap = ,
    contrast_gene_category_network = ,
    contrast_kegg_map = "contrast",
    NA_character_)
}

# The explicit scales the native gene heatmap offers. These are the builder's
# own `--scale` values, not new science: `zscore` is per-gene z-score across
# samples and `raw` plots the saved matrix values unscaled.
lisa_explore_heatmap_variants <- function() c("zscore", "raw")

# The single accepted exploration variant of the gene-category network. It maps
# onto the network builder's own accepted settings -- `top_categories = 0`
# (every category) and `top_genes_per_category = 7` -- so naming it selects the
# accepted figure rather than introducing a new scientific control. There is
# deliberately no second variant: H3_SCOPE forbids a category-local crop and any
# new scientific setting control.
lisa_explore_network_variants <- function() "all_categories_top7"

# Products whose request carries a concrete entity, and products whose request
# carries an explicit variant. Anything not listed here must carry neither, so a
# stray entity or variant cannot quietly widen a request or, worse, make two
# different requests collide on one key.
lisa_explore_product_entity <- function(product) {
  as.character(product) %in% c("kegg_pathway_map", "contrast_kegg_map")
}

lisa_explore_product_variants <- function(product) {
  switch(as.character(product),
    heatmap = lisa_explore_heatmap_variants(),
    contrast_gene_category_network = lisa_explore_network_variants(),
    character())
}

# A KEGG pathway identifier is three lowercase organism letters and five digits
# (`hsa04110`). Validated here so a malformed id is refused when the request is
# built, not after a worker has been started for it.
lisa_explore_valid_pathway_id <- function(entity) {
  grepl("^[a-z]{3}[0-9]{5}$", as.character(entity))
}

lisa_explore_states <- function() {
  c("ungenerated", "queued", "running", "available", "not_applicable", "failed")
}

lisa_explore_scalar <- function(value, field, required = TRUE) {
  if (is.null(value)) {
    if (isTRUE(required)) {
      stop("LISA-EXPLORE-001 ", field, " is required and must name exactly one value.",
           call. = FALSE)
    }
    return(NA_character_)
  }
  if (length(value) != 1L) {
    stop("LISA-EXPLORE-002 ", field, " must name exactly one value; ",
         length(value), " were given. A figure request is one figure: build one ",
         "request per figure instead of passing a vector, which would otherwise ",
         "expand into a cross-product.", call. = FALSE)
  }
  value <- as.character(value)
  if (is.na(value) || !nzchar(trimws(value))) {
    stop("LISA-EXPLORE-003 ", field, " must be a non-empty value.", call. = FALSE)
  }
  trimws(value)
}

#' Name one exact figure to explore
#'
#' Builds the request object used by every exploration entry point. A request
#' always denotes a single figure: passing a vector to any field is an error, so
#' a selection can never silently widen into a Cartesian product.
#'
#' @param unit_type `"single_de"` or `"contrast"`.
#' @param collection The collection identifier the figure belongs to.
#' @param product The figure product, for example `"volcano"`.
#' @param category_id The category the figure is drawn for. Required for a
#'   category product and refused for a native collection-wide product such as
#'   `"contrast_paired_heatmap"`, which is one figure for the whole collection.
#' @param analysis_id Analysis identifier; required when `unit_type` is
#'   `"single_de"`.
#' @param contrast_id Contrast identifier; required when `unit_type` is
#'   `"contrast"`.
#' @param entity Optional concrete entity, such as a KEGG map identifier.
#' @param variant Optional explicit variant, such as a heatmap scale.
#' @return An object of class `lisa_figure_request`.
#' @export
#'
#' @examples
#' request <- lisa_figure_request(
#'   unit_type = "single_de",
#'   analysis_id = "demo_analysis",
#'   collection = "H",
#'   category_id = "SYN_SIGNAL",
#'   product = "volcano"
#' )
#' lisa_request_id(request)
lisa_figure_request <- function(unit_type, collection, product,
                                category_id = NULL,
                                analysis_id = NULL, contrast_id = NULL,
                                entity = NULL, variant = NULL) {
  unit_type <- lisa_explore_scalar(unit_type, "unit_type")
  if (!unit_type %in% c("single_de", "contrast")) {
    stop("LISA-EXPLORE-004 unit_type must be single_de or contrast.", call. = FALSE)
  }
  collection <- lisa_explore_scalar(collection, "collection")
  product <- lisa_explore_scalar(product, "product")
  if (!product %in% lisa_explore_catalog_products()) {
    stop("LISA-EXPLORE-005 unknown product: ", product, ". Known products: ",
         paste(lisa_explore_catalog_products(), collapse = ", "), ".", call. = FALSE)
  }
  expected_unit <- lisa_explore_product_unit_type(product)
  if (!is.na(expected_unit) && !identical(unit_type, expected_unit)) {
    stop("LISA-EXPLORE-015 ", product, " is a ", expected_unit,
         " product; it cannot be requested with unit_type = ", unit_type, ".",
         call. = FALSE)
  }

  needs_category <- lisa_explore_product_needs_category(product)
  category_id <- lisa_explore_scalar(category_id, "category_id",
                                     required = needs_category)
  if (!needs_category && !is.na(category_id)) {
    stop("LISA-EXPLORE-014 ", product, " is a collection-wide figure: it is ",
         "drawn once for the whole collection, so category_id does not apply. ",
         "Omit it; do not pass the browsed category or a placeholder such as ",
         "ALL or GLOBAL.", call. = FALSE)
  }

  # The owner field is mandatory for its unit type and forbidden for the other,
  # so a request can never be ambiguous about which unit owns the figure.
  if (identical(unit_type, "single_de")) {
    analysis_id <- lisa_explore_scalar(analysis_id, "analysis_id")
    if (!is.null(contrast_id)) {
      stop("LISA-EXPLORE-006 contrast_id does not apply to a single_de request.",
           call. = FALSE)
    }
    contrast_id <- NA_character_
  } else {
    contrast_id <- lisa_explore_scalar(contrast_id, "contrast_id")
    if (!is.null(analysis_id)) {
      stop("LISA-EXPLORE-006 analysis_id does not apply to a contrast request.",
           call. = FALSE)
    }
    analysis_id <- NA_character_
  }

  entity <- lisa_explore_scalar(entity, "entity", required = FALSE)
  variant <- lisa_explore_scalar(variant, "variant", required = FALSE)

  if (lisa_explore_product_entity(product)) {
    if (is.na(entity)) {
      stop("LISA-EXPLORE-008 ", product, " needs an explicit entity naming one ",
           "pathway, for example entity = \"hsa04110\". A request for ",
           "\"whichever map\" is not one figure.", call. = FALSE)
    }
    if (!lisa_explore_valid_pathway_id(entity)) {
      stop("LISA-EXPLORE-009 entity must be one KEGG pathway identifier such as ",
           "hsa04110; received: ", entity, ".", call. = FALSE)
    }
  } else if (!is.na(entity)) {
    stop("LISA-EXPLORE-010 entity does not apply to product ", product, ".",
         call. = FALSE)
  }

  allowed_variants <- lisa_explore_product_variants(product)
  if (length(allowed_variants)) {
    if (is.na(variant)) {
      stop("LISA-EXPLORE-011 ", product, " needs one explicit variant; ",
           "available: ", paste(allowed_variants, collapse = ", "),
           ". Rendering every variant would be a cross-product.", call. = FALSE)
    }
    if (!variant %in% allowed_variants) {
      stop("LISA-EXPLORE-012 unknown ", product, " variant: ", variant,
           ". Available: ", paste(allowed_variants, collapse = ", "), ".",
           call. = FALSE)
    }
  } else if (!is.na(variant)) {
    stop("LISA-EXPLORE-013 variant does not apply to product ", product, ".",
         call. = FALSE)
  }

  request <- list(
    unit_type = unit_type,
    analysis_id = analysis_id,
    contrast_id = contrast_id,
    collection = collection,
    category_id = category_id,
    product = product,
    entity = entity,
    variant = variant
  )
  structure(request, class = "lisa_figure_request")
}

lisa_explore_is_request <- function(x) inherits(x, "lisa_figure_request")

lisa_explore_assert_request <- function(request) {
  if (!lisa_explore_is_request(request)) {
    stop("LISA-EXPLORE-007 a lisa_figure_request is required. Build it with ",
         "lisa_figure_request().", call. = FALSE)
  }
  request
}

lisa_explore_request_fields <- function() {
  c("unit_type", "analysis_id", "contrast_id", "collection", "category_id",
    "product", "entity", "variant")
}

lisa_explore_request_vector <- function(request) {
  lisa_explore_assert_request(request)
  values <- vapply(lisa_explore_request_fields(), function(field) {
    value <- request[[field]]
    if (is.null(value) || is.na(value)) "" else as.character(value)
  }, character(1))
  values
}

lisa_explore_identity_fields <- function(product) {
  fields <- lisa_explore_request_fields()
  product <- as.character(product)
  if (identical(product, "kegg_pathway_map") ||
      identical(lisa_explore_product_scope(product), "collection")) {
    return(setdiff(fields, "category_id"))
  }
  fields
}

lisa_explore_identity_vector <- function(request) {
  lisa_explore_assert_request(request)
  values <- lisa_explore_request_vector(request)
  values[lisa_explore_identity_fields(request$product)]
}

#' Stable identifier of an exact figure request
#'
#' Derive a deterministic identifier from the canonical fields of a figure
#' request. The identifier labels the request; it does not render a figure
#' or certify that an artifact has been generated.
#'
#' @param request A [lisa_figure_request()].
#' @return A short stable identifier for the requested figure.
#' @export
#'
#' @examples
#' lisa_request_id(lisa_figure_request(
#'   unit_type = "single_de", analysis_id = "demo", collection = "H",
#'   category_id = "SYN_SIGNAL", product = "volcano"
#' ))
lisa_request_id <- function(request) {
  values <- lisa_explore_request_vector(request)
  canonical <- paste(paste(names(values), values, sep = "="), collapse = "\n")
  paste0("req-", substr(lisa_sha256_text(canonical), 1L, 16L))
}

#' Describe an exact figure request as one row
#'
#' @param request A [lisa_figure_request()].
#' @return A one-row data frame describing the request.
#' @export
#'
#' @examples
#' as.data.frame(lisa_figure_request(
#'   unit_type = "single_de", analysis_id = "demo", collection = "H",
#'   category_id = "SYN_SIGNAL", product = "volcano"
#' ))
as.data.frame.lisa_figure_request <- function(x, ...) {
  values <- lisa_explore_request_vector(x)
  row <- as.data.frame(as.list(values), stringsAsFactors = FALSE)
  row$request_id <- lisa_request_id(x)
  row
}

#' @export
print.lisa_figure_request <- function(x, ...) {
  cat("<lisa_figure_request>", lisa_request_id(x), "\n")
  values <- lisa_explore_request_vector(x)
  for (field in names(values)) {
    if (nzchar(values[[field]])) {
      cat(sprintf("  %-12s %s\n", field, values[[field]]))
    }
  }
  invisible(x)
}

lisa_explore_request_from_row <- function(row) {
  owner <- if (identical(as.character(row$unit_type), "single_de")) {
    list(analysis_id = as.character(row$analysis_id), contrast_id = NULL)
  } else {
    list(analysis_id = NULL, contrast_id = as.character(row$contrast_id))
  }
  optional <- function(value) {
    if (is.null(value) || is.na(value) || !nzchar(as.character(value))) NULL
    else as.character(value)
  }
  lisa_figure_request(
    unit_type = as.character(row$unit_type),
    collection = as.character(row$collection),
    product = as.character(row$product),
    category_id = optional(row$category_id),
    analysis_id = owner$analysis_id,
    contrast_id = owner$contrast_id,
    entity = optional(row$entity),
    variant = optional(row$variant)
  )
}

# The presentation policy for an exploration artifact requires
# PNG and PDF to be two formats of the same figure, plus its source data and its
# recipe; they are never separate scientific requests.
lisa_explore_report_policy <- function() {
  list(mode = "selected",
       formats = list(png = TRUE, svg = FALSE, pdf = TRUE),
       source_data = TRUE,
       recipes = TRUE)
}

lisa_explore_policy_digest <- function(policy = lisa_explore_report_policy()) {
  flat <- c(
    mode = policy$mode,
    png = as.character(isTRUE(policy$formats$png)),
    svg = as.character(isTRUE(policy$formats$svg)),
    pdf = as.character(isTRUE(policy$formats$pdf)),
    source_data = as.character(isTRUE(policy$source_data)),
    recipes = as.character(isTRUE(policy$recipes))
  )
  lisa_sha256_text(paste(paste(names(flat), flat, sep = "="), collapse = "\n"))
}

# Translate one exact request into the selection understood by the existing
# extension API. Only scalars are ever placed in the selection, which is what
# keeps `lisa_extension_filter()` from producing a cross-product.
lisa_explore_selection <- function(request, policy = lisa_explore_report_policy()) {
  lisa_explore_assert_request(request)
  selection <- list(
    mode = "selected",
    collections = request$collection,
    products = request$product
  )
  if (!is.na(request$category_id)) selection$categories <- request$category_id
  if (identical(request$unit_type, "single_de")) {
    selection$analyses <- request$analysis_id
  } else {
    selection$contrasts <- request$contrast_id
  }
  if (!is.na(request$entity)) selection$entities <- request$entity
  if (!is.na(request$variant)) selection$variants <- request$variant
  list(selection = selection, report = policy)
}
