# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Helpers for immutable, offline KEGG painter resources.

kegg_species_contract <- function(species) {
  key <- tolower(trimws(as.character(species[[1]])))
  if (key %in% c("human", "homo sapiens")) {
    return(list(scientific_name = "Homo sapiens", kegg_code = "hsa", annotation_db = "org.Hs.eg.db"))
  }
  if (key %in% c("mouse", "mus musculus")) {
    return(list(scientific_name = "Mus musculus", kegg_code = "mmu", annotation_db = "org.Mm.eg.db"))
  }
  stop("LISA-SPECIES-001 species is missing or unsupported. Repair: use Homo sapiens or Mus musculus.", call. = FALSE)
}

kegg_snapshot_resource <- function(cfg, resource_type, key) {
  access_mode <- if (!is.null(cfg$kegg_access_mode)) cfg$kegg_access_mode else "external"
  mode <- if (identical(access_mode, "external")) "prefer_cache" else "cache_only"
  lisaR:::lisa_resolve_kegg_resource(
    mode = mode, cache_root = cfg$kegg_cache_root,
    snapshot_id = cfg$kegg_snapshot_id, organism = cfg$kegg_code,
    resource_type = resource_type, key = key,
    fetch = if (identical(mode, "prefer_cache")) {
      function() lisaR:::lisa_fetch_kegg_resource(
        cfg$kegg_code, resource_type, key
      )
    } else NULL
  )$path
}

kegg_snapshot_rds <- function(cfg, resource_type, key) {
  readRDS(kegg_snapshot_resource(cfg, resource_type, key))
}

kegg_annotation_database <- function(cfg) {
  if (!requireNamespace(cfg$annotation_db, quietly = TRUE)) {
    stop(sprintf("Missing species annotation package: %s", cfg$annotation_db), call. = FALSE)
  }
  getExportedValue(cfg$annotation_db, cfg$annotation_db)
}
