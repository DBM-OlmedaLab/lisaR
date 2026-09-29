# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Versioned non-scientific contracts for lisaR development.
#
# These constants identify the package and report contracts only. They do not
# select resources, alter analytical parameters, or change scientific output.

lisa_contract_versions <- function() {
  list(
    package_development_version = "1.0.0",
    package_r_version = "1.0.0",
    pipeline_schema_version = lisa_pipeline_schema_version(),
    report_schema_version = "7.0.0",
    contract_manifest_version = "1.0.0"
  )
}

#' Create a minimal lisaR contract-manifest record
#'
#' This Initial scaffold creates an in-memory contract identity. Run-time
#' capture, resource identities, and file checksums are intentionally deferred
#' to later approved milestones.
#'
#' @param source_tree A source commit or source-tree identifier, when known.
#' @param configuration_schema_version A configuration-schema version, when
#'   known. It is `NA` until the versioned configuration schema is implemented.
#'
#' @return A `lisa_contract_manifest` list.
#' @keywords internal
new_lisa_contract_manifest <- function(
    source_tree = NA_character_,
    configuration_schema_version = NA_character_) {
  versions <- lisa_contract_versions()

  structure(
    list(
      contract_manifest_version = versions$contract_manifest_version,
      package_development_version = versions$package_development_version,
      package_r_version = versions$package_r_version,
      pipeline_schema_version = versions$pipeline_schema_version,
      report_schema_version = versions$report_schema_version,
      configuration_schema_version = as.character(configuration_schema_version),
      source_tree = as.character(source_tree)
    ),
    class = "lisa_contract_manifest"
  )
}
