# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

run_lisa_single_de_heatmaps <- function(project_dir,
                                        analysis_id,
                                        universe = "GOBP-C2",
                                        scale = "zscore",
                                        expression_matrix = NULL) {
  package_dir <- lisa_resolve_package_dir()
  args <- c(
    "--project-dir", project_dir,
    "--analysis-id", analysis_id,
    "--universe", universe,
    "--scale", scale
  )
  if (!is.null(expression_matrix) && nzchar(as.character(expression_matrix))) {
    args <- c(args, "--expression-matrix", as.character(expression_matrix))
  }
  lisa_run_post_script(package_dir, "build_single_de_leading_edge_gene_heatmaps.R", args)
}

run_lisa_single_de_kegg_maps <- function(project_dir,
                                         analysis_id,
                                         universe = "GOBP-C2",
                                         species,
                                         kegg_cache_root = getOption("lisaR.kegg_cache_root", Sys.getenv("LISAR_KEGG_CACHE", "")),
                                         kegg_snapshot_id = getOption("lisaR.kegg_snapshot_id", Sys.getenv("LISAR_KEGG_SNAPSHOT", ""))) {
  package_dir <- lisa_resolve_package_dir()
  lisa_run_post_script(
    package_dir,
    "build_single_de_kegg_pathway_painter.R",
    c("--project-dir", project_dir, "--analysis-id", analysis_id, "--universe", universe,
      "--species", species, "--kegg-cache-root", kegg_cache_root,
      "--kegg-snapshot-id", kegg_snapshot_id)
  )
}

run_lisa_single_de_recurrent_gene_screen <- function(project_dir,
                                                     analysis_id,
                                                     universe = "GOBP-C2",
                                                     top_genes = 40,
                                                     min_recurrent_categories = 3) {
  package_dir <- lisa_resolve_package_dir()
  lisa_run_post_script(
    package_dir,
    "build_single_de_recurrent_gene_screen.R",
    c(
      "--project-dir", project_dir,
      "--analysis-id", analysis_id,
      "--universe", universe,
      "--top-genes", as.character(top_genes),
      "--min-recurrent-categories", as.character(min_recurrent_categories)
    )
  )
}

run_lisa_contrast_heatmaps <- function(project_dir,
                                       contrast_id,
                                       universe = "GOBP-C2") {
  package_dir <- lisa_resolve_package_dir()
  lisa_run_post_script(
    package_dir,
    "build_contrast_macrogroup_heatmaps.R",
    c("--project-dir", project_dir, "--contrast-id", contrast_id, "--universe", universe)
  )
}

run_lisa_contrast_kegg_maps <- function(project_dir,
                                        contrast_id,
                                        universe = "GOBP-C2",
                                        species,
                                        kegg_cache_root = getOption("lisaR.kegg_cache_root", Sys.getenv("LISAR_KEGG_CACHE", "")),
                                        kegg_snapshot_id = getOption("lisaR.kegg_snapshot_id", Sys.getenv("LISAR_KEGG_SNAPSHOT", ""))) {
  package_dir <- lisa_resolve_package_dir()
  lisa_run_post_script(
    package_dir,
    "build_contrast_kegg_pathway_painter.R",
    c("--project-dir", project_dir, "--contrast-id", contrast_id, "--universe", universe,
      "--species", species, "--kegg-cache-root", kegg_cache_root,
      "--kegg-snapshot-id", kegg_snapshot_id)
  )
}

run_lisa_contrast_kegg_pathway_painter <- function(project_dir,
                                                   contrast_id,
                                                   universe = "GOBP-C2",
                                                   species,
                                                   kegg_cache_root = getOption("lisaR.kegg_cache_root", Sys.getenv("LISAR_KEGG_CACHE", "")),
                                                   kegg_snapshot_id = getOption("lisaR.kegg_snapshot_id", Sys.getenv("LISAR_KEGG_SNAPSHOT", "")),
                                                   max_abs_log2fc = 0.5,
                                                   color_power = 1.0) {
  package_dir <- lisa_resolve_package_dir()
  lisa_run_post_script(
    package_dir,
    "build_contrast_kegg_pathway_painter.R",
    c(
      "--project-dir", project_dir,
      "--contrast-id", contrast_id,
      "--universe", universe,
      "--species", species,
      "--kegg-cache-root", kegg_cache_root,
      "--kegg-snapshot-id", kegg_snapshot_id,
      "--max-abs-log2fc", as.character(max_abs_log2fc),
      "--color-power", as.character(color_power)
    )
  )
}
