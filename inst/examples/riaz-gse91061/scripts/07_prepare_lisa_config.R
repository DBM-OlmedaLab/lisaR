#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Materialise the canonical YAML configuration and an equivalent JSON file.
#
# YAML is the tutorial format because comments can explain every decision.
# JSON is generated from YAML rather than maintained independently, preventing
# the two examples from drifting apart.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

riaz_require_packages(c("jsonlite", "yaml"))

template <- file.path(
  dirname(riaz_example_dir()),
  "config",
  "riaz-gse91061.yml"
)
riaz_assert(file.exists(template), paste("Missing YAML template:", template))

config_dir <- file.path(riaz_project_dir(), "config")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
yaml_path <- file.path(config_dir, "riaz-gse91061.yml")
json_path <- file.path(config_dir, "riaz-gse91061.json")

copied <- file.copy(template, yaml_path, overwrite = TRUE, copy.mode = TRUE)
riaz_assert(copied, paste("Could not write runtime YAML:", yaml_path))

configuration <- yaml::read_yaml(yaml_path)
jsonlite::write_json(
  configuration,
  json_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null",
  digits = NA
)

# Parse both representations and require semantic identity before lisaR sees
# either file. This checks values and structure; YAML comments intentionally do
# not exist in JSON.
yaml_object <- yaml::read_yaml(yaml_path)
json_object <- jsonlite::read_json(json_path, simplifyVector = FALSE)
riaz_assert(
  identical(
    jsonlite::toJSON(
      yaml_object, auto_unbox = TRUE, null = "null", digits = NA
    ),
    jsonlite::toJSON(
      json_object, auto_unbox = TRUE, null = "null", digits = NA
    )
  ),
  "Generated JSON is not semantically identical to the canonical YAML."
)

manifest <- data.frame(
  format = c("YAML", "JSON"),
  # Keep this receipt valid after the enclosing prepared-project staging
  # directory is atomically renamed by stage 11.
  path = file.path("config", basename(c(yaml_path, json_path))),
  sha256 = c(riaz_sha256(yaml_path), riaz_sha256(json_path)),
  stringsAsFactors = FALSE
)
riaz_write_tsv(manifest, file.path(config_dir, "config_manifest.tsv"))

print(manifest, row.names = FALSE)
cat("CONFIG_EQUIVALENCE_GATE=PASS\n")
