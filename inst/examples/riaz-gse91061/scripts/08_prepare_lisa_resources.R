#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Prepare the ordinary C1 resources used by the Riaz example.
#
# lisa_core and lisa_category_map resolve from installed lisaR. TERM2GENE is
# deliberately not copied into this project: the normal C1 helper verifies the
# pinned MSigDB 2026.1 ZIP and registers it in LISAR_RIAZ_RESOURCE_CACHE (or
# lisaR's normal cache). A local ZIP grants no redistribution right. Without a
# ZIP, the helper requires the caller's own explicit TRUE assertion before its
# fixed public acquisition route; this script never accepts terms silently.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

riaz_require_packages(c("data.table", "lisaR"))
config_path <- file.path(riaz_project_dir(), "config", "riaz-gse91061.yml")
riaz_assert(file.exists(config_path), "Run 07_prepare_lisa_config.R first.")

resource_cache <- riaz_configure_normal_resources()
source_zip <- Sys.getenv("LISAR_RIAZ_MSIGDB_ZIP", unset = "")
accept_terms <- identical(
  Sys.getenv("LISAR_RIAZ_ACCEPT_MSIGDB_TERMS", unset = ""), "TRUE"
)
prepared <- lisaR::prepare_lisa_msigdb_resource(
  source_zip = if (nzchar(source_zip)) source_zip else NULL,
  accept_terms = accept_terms,
  cache_root = resource_cache
)
validation <- lisaR::validate_lisa_config(
  config_path, check_files = TRUE, strict = TRUE
)
riaz_assert(isTRUE(validation$execution_ready),
            "Normal C1 resources did not make the Riaz configuration ready.")
resolved <- validation$resources
riaz_assert(
  nrow(resolved) == 3L && all(resolved$checked) && all(resolved$ready),
  "Expected three checked, ready ordinary C1 resources."
)
resolved$resource_cache <- resource_cache
resolved$term2gene_source_zip_sha256 <- prepared$source_zip_sha256
resolved$term2gene_acquired_this_call <- prepared$acquired
resolved$status <- "ordinary_c1_resources_verified"
riaz_write_tsv(resolved, file.path(
  riaz_project_dir(), "logs", "lisa_resource_resolution.tsv"
))
print(resolved, row.names = FALSE)
cat("LISA_RESOURCE_GATE=PASS\n")
