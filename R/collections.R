# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

lisa_builtin_collection_registry <- function() {
  data.frame(
    analysis_collection = c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS", "HALLMARKS"),
    dictionary_id = c("LISA_GOBP_C2", "LISA_GOMF", "LISA_GOCC", "LISA_PATHWAYS", "MSIGDB_HALLMARKS"),
    dictionary_kind = c("semantic_lisa", "semantic_lisa", "semantic_lisa", "semantic_lisa", "msigdb_direct"),
    post_lisa_profile = c("full", "full", "full", "full", "hallmarks_minimal"),
    run_gsea = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    # ORA is available for an explicit opt-in but is not part of the default
    # ranked-GSEA workflow.
    run_ora = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    run_contrasts = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    run_category_pathway_plots = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    run_gene_cards = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    run_volcano_overlays = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    run_recurrent_gene_screen = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    run_enrichmentmap = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    run_kegg_layers = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    allow_missing_dictionary = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

lisa_default_collections <- function() {
  lisa_builtin_collection_registry()$analysis_collection
}

normalize_lisa_collections <- function(collections) {
  if (missing(collections) || is.null(collections)) {
    collections <- lisa_default_collections()
  }
  collections <- as.character(collections)
  if (length(collections) == 1 && tolower(collections) == "all") {
    return(lisa_default_collections())
  }

  aliases <- c(
    gobp = "GOBP-C2",
    gobpc2 = "GOBP-C2",
    `gobp-c2` = "GOBP-C2",
    c2 = "GOBP-C2",
    bp = "GOBP-C2",
    gomf = "GOMF",
    mf = "GOMF",
    gocc = "GOCC",
    cc = "GOCC",
    pathways = "PATHWAYS",
    pathway = "PATHWAYS",
    lisa_pathways = "PATHWAYS",
    lisa_patways = "PATHWAYS",
    hallmarks = "HALLMARKS",
    hallmark = "HALLMARKS",
    h = "HALLMARKS",
    msigdb_h = "HALLMARKS"
  )

  key <- tolower(collections)
  out <- collections
  out[key %in% names(aliases)] <- unname(aliases[key[key %in% names(aliases)]])

  valid <- lisa_builtin_collection_registry()$analysis_collection
  bad <- setdiff(out, valid)
  if (length(bad) > 0) {
    stop(
      "Unknown LISA analysis collection(s): ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  unique(out)
}

lisa_collection_registry <- function(collections = lisa_default_collections()) {
  collections <- normalize_lisa_collections(collections)
  registry <- lisa_builtin_collection_registry()
  registry <- registry[match(collections, registry$analysis_collection), , drop = FALSE]
  rownames(registry) <- NULL
  registry
}

lisa_collection_dictionary_ids <- function(collections = lisa_default_collections()) {
  registry <- lisa_collection_registry(collections)
  stats::setNames(registry$dictionary_id, registry$analysis_collection)
}
