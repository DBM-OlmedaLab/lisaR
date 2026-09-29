# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Command-line interfaces for the canonical post-processing scripts.
#
# Keeping this contract in the package lets the orchestrator reject an
# unsupported argument before it starts a child R process. Tests compare this
# table with every script's own parser, so the two sides cannot drift silently.
lisa_post_script_interfaces <- function() {
  list(
    "build_category_evidence.R" = c(
      "project-dir", "analysis-id", "universe", "tier", "gsea-padj-cutoff",
      "de-padj-cutoff", "positive-contrast", "max-sets", "max-genes", "formats",
      "output-dir", "category-dictionary", "categories", "source-data", "recipes", "overlap-export"
    ),
    "build_category_navigation.R" = c(
      "project-dir", "analysis-id", "universe", "tier", "gsea-padj-cutoff",
      "positive-contrast", "output-dir"
    ),
    "build_gene_evidence.R" = c(
      "project-dir", "analysis-ids", "term2gene", "tier", "gsea-padj-cutoff",
      "de-padj-cutoff", "output-dir"
    ),
    "build_contrast_evidence.R" = c(
      "project-dir", "contrast-id", "analysis-a", "analysis-b", "label", "universe", "tier",
      "evidence-a", "evidence-b", "contrast-summary", "output-dir", "max-sets", "max-genes",
      "formats", "variants", "report-href", "gene-explorer-href"
    ),
    "build_contrast_navigation.R" = c(
      "project-dir", "contrast-id", "scope", "analysis-a", "analysis-b", "universe", "tier",
      "gsea-padj-cutoff", "label", "contrast-summary", "output-dir"
    ),
    # --report-mode is how run_lisa() asks its one report generator for a
    # standard or full report; omitted, the saved run metadata decides.
    "build_LISA_report.R" = c(
      "project-dir", "title", "presentation-config", "report-mode"
    ),
    # Public reusable entry point that assembles a scoped report from saved
    # results. --analyses/--contrasts carry the exact selection; "none" is the
    # accepted zero-contrast form. --plan-only stops at the selection summary,
    # before any copy. --assemble false prepares the selected input view and its
    # category links, then stops without invoking the report generator.
    # --complete-missing adds only absent presentation products, never science.
    "build_selected_LISA_report.R" = c(
      "source-dir", "output-dir", "analyses", "contrasts",
      "artifacts-dir", "artifact-manifest", "title",
      "plan-only", "complete-missing", "assemble",
      "kegg-maps", "kegg-cache-root", "kegg-snapshot-id", "kegg-access-mode"
    ),
    "build_active_category_map.R" = c(
      "core", "expanded", "source-map", "source-map-sha256",
      "output", "audit"
    ),
    "build_compact_dictionary_resources.R" = c(
      "dictionary-build-dir", "category-palette", "pathways-hierarchy",
      "category-map-sha256", "output-dir", "manifest"
    ),
    "build_contrast_gene_category_network.R" = c(
      "project-dir", "contrast-id", "universe", "top-categories",
      "top-genes-per-category",
      "variant", "paired-input", "summary-input"
    ),
    # `tables-only` is the plot-free shared preparation; `exact-product` plus
    # `category-id` select exactly one figure; `paired-input`/`summary-input`
    # reuse a prepared bundle. Omitting all five reproduces the established
    # default: prepare the tables, draw every card, draw the paired heatmap.
    "build_contrast_gene_level_product.R" = c(
      "project-dir", "contrast-id", "universe", "top-categories",
      "top-genes-per-category", "lfc-cap",
      "tables-only", "exact-product", "category-id",
      "paired-input", "summary-input"
    ),
    "build_contrast_kegg_map_layer.R" = c(
      "project-dir", "contrast-id", "universe", "top-pathways",
      "top-genes", "lfc-cap"
    ),
    "build_contrast_kegg_pathway_painter.R" = c(
      "project-dir", "contrast-id", "universe", "species",
      "kegg-cache-root", "kegg-snapshot-id", "kegg-access-mode",
      "top-pathways", "top-genes", "max-abs-log2fc", "color-power",
      "kegg-id", "emit-pathway-index", "index-only",
      "paired-input", "summary-input"
    ),
    "build_contrast_macrogroup_heatmaps.R" = c(
      "project-dir", "contrast-id", "universe", "max-macrogroups",
      "max-categories-per-macrogroup", "max-genes", "lfc-cap",
      "macrogroups"
    ),
    "build_single_de_category_count_heatmaps.R" = c(
      "project-dir", "universe", "max-categories"
    ),
    "build_single_de_category_gene_cards.R" = c(
      "project-dir", "analysis-id", "universe", "max-categories",
      "top-genes", "categories", "formats"
    ),
    "build_single_de_category_volcano_overlays.R" = c(
      "project-dir", "analysis-id", "universe", "max-categories",
      "label-genes", "categories", "formats", "de-padj-cutoff", "lfc-cutoff"
    ),
    "build_single_de_enrichmentmap.R" = c(
      "project-dir", "analysis-id", "universe", "min-jaccard",
      "max-nodes", "max-edges", "top-genes-per-category",
      "render-graphs"
    ),
    # `kegg-id`, `emit-pathway-index` and `index-only` are the exact-selector
    # flags. Omitting all three reproduces the established default: rank the
    # collection, paint the top `top-pathways`.
    "build_single_de_kegg_pathway_painter.R" = c(
      "project-dir", "analysis-id", "universe", "species",
      "kegg-cache-root", "kegg-snapshot-id", "kegg-access-mode",
      "top-pathways", "max-abs-log2fc", "color-power",
      "kegg-id", "emit-pathway-index", "index-only"
    ),
    "build_single_de_leading_edge_gene_heatmaps.R" = c(
      "project-dir", "analysis-id", "universe", "expression-matrix",
      "expression-matrix-column", "scale", "top-genes", "categories", "formats"
    ),
    "build_single_de_recurrent_gene_screen.R" = c(
      "project-dir", "analysis-id", "universe", "top-genes",
      "top-categories", "min-recurrent-categories",
      "min-high-categories", "min-high-macrogroups",
      "de-padj-cutoff", "lfc-cutoff"
    )
  )
}

lisa_validate_post_script_args <- function(script_name, args) {
  interfaces <- lisa_post_script_interfaces()
  allowed <- interfaces[[script_name]]
  if (is.null(allowed)) {
    stop(
      "No registered command-line interface for post-processing script: ",
      script_name,
      call. = FALSE
    )
  }
  if (length(args) %% 2L != 0L) {
    stop(
      "Post-processing arguments must be flag/value pairs for ",
      script_name,
      call. = FALSE
    )
  }
  flags <- args[seq.int(1L, length(args), by = 2L)]
  if (any(!startsWith(flags, "--"))) {
    stop("Invalid post-processing flag for ", script_name, call. = FALSE)
  }
  names <- sub("^--", "", flags)
  duplicated_names <- unique(names[duplicated(names)])
  if (length(duplicated_names)) {
    stop(
      "Duplicate post-processing argument(s) for ",
      script_name,
      ": ",
      paste(paste0("--", duplicated_names), collapse = ", "),
      call. = FALSE
    )
  }
  unknown <- setdiff(names, allowed)
  if (length(unknown)) {
    stop(
      "Unsupported argument(s) for ",
      script_name,
      ": ",
      paste(paste0("--", unknown), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
