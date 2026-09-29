# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

lisa_pipeline_plan <- function(
  de_index,
  contrast_index = NULL,
  collections = lisa_default_collections(),
  run_gene_level = FALSE,
  run_reports = TRUE,
  run_kegg_maps = FALSE,
  run_hallmarks = TRUE,
  run_ora = FALSE,
  run_kegg = NULL,
  category_evidence = TRUE
) {
  if (!is.null(run_kegg)) {
    warning("run_kegg is deprecated; migrate to run_kegg_maps. It is not interpreted automatically.", call. = FALSE)
  }
  de_index <- validate_lisa_de_index(de_index, require_files = FALSE)
  contrast_index <- validate_lisa_contrast_index(contrast_index, de_index = de_index)
  registry <- lisa_collection_registry(collections)
  if (!isTRUE(run_hallmarks)) registry <- registry[registry$analysis_collection != "HALLMARKS", , drop = FALSE]
  registry$run_ora <- registry$run_ora & isTRUE(run_ora)

  has_full_collections <- any(registry$post_lisa_profile == "full")
  has_hallmarks <- any(registry$analysis_collection == "HALLMARKS")
  has_contrast_collections <- any(registry$run_contrasts)
  has_kegg_collections <- any(registry$run_kegg_layers)
  full_collections <- registry$analysis_collection[registry$post_lisa_profile == "full"]
  hallmark_collections <- registry$analysis_collection[registry$analysis_collection == "HALLMARKS"]
  contrast_collections <- registry$analysis_collection[registry$run_contrasts]
  kegg_collections <- registry$analysis_collection[registry$run_kegg_layers]
  enrichmentmap_collections <- registry$analysis_collection[registry$run_enrichmentmap]

  stage_rows <- list()
  add_stage <- function(stage, stage_collections, n_analyses = nrow(de_index), n_contrasts = nrow(contrast_index)) {
    if (length(stage_collections) == 0) {
      return(invisible(NULL))
    }
    stage_registry <- lisa_collection_registry(stage_collections)
    stage_rows[[length(stage_rows) + 1L]] <<- data.frame(
      stage_order = length(stage_rows) + 1L,
      stage = stage,
      n_analyses = n_analyses,
      n_contrasts = n_contrasts,
      collections = paste(stage_registry$analysis_collection, collapse = ","),
      dictionaries = paste(stage_registry$dictionary_id, collapse = ","),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  add_stage("single_de_lisa", registry$analysis_collection)
  add_stage("single_de_category_inference", setdiff(registry$analysis_collection, "HALLMARKS"), n_contrasts = 0L)
  add_stage("single_de_enrichmentmap", enrichmentmap_collections)
  if (isTRUE(category_evidence)) {
    add_stage("single_de_category_evidence", full_collections, n_contrasts = 0L)
    add_stage("single_de_category_navigation", full_collections, n_contrasts = 0L)
    add_stage("single_de_gene_evidence", full_collections, n_contrasts = 0L)
  }
  if (nrow(contrast_index) > 0 && isTRUE(has_contrast_collections)) {
    add_stage("category_contrasts", contrast_collections)
    if (isTRUE(category_evidence)) add_stage("contrast_category_evidence", setdiff(contrast_collections, "HALLMARKS"), n_analyses = 0L)
  }
  if (isTRUE(run_gene_level)) {
    if (isTRUE(has_full_collections)) {
      add_stage("single_de_gene_level_tables", full_collections)
      add_stage("single_de_gene_cards", full_collections)
      add_stage("single_de_volcano_overlays", full_collections)
      add_stage("single_de_recurrent_gene_screen", full_collections)
      add_stage("single_de_leading_edge_gene_heatmaps", full_collections)
    }
    if (isTRUE(has_hallmarks)) {
      add_stage("single_de_hallmarks_lollipop", hallmark_collections)
      add_stage("single_de_hallmarks_gene_cards", hallmark_collections)
      add_stage("single_de_hallmarks_volcano_overlays", hallmark_collections)
      add_stage("single_de_hallmarks_leading_edge_gene_heatmaps", hallmark_collections)
    }
    if (nrow(contrast_index) > 0 && isTRUE(has_contrast_collections)) {
      add_stage("contrast_gene_level_tables", contrast_collections)
      add_stage("contrast_category_cards", contrast_collections)
      add_stage("contrast_paired_gene_heatmaps", contrast_collections)
      add_stage("contrast_gene_category_networks", contrast_collections)
    }
  }
  if (isTRUE(run_kegg_maps) && isTRUE(has_kegg_collections)) {
    add_stage("kegg_pathway_layers", kegg_collections)
    add_stage("kegg_pathway_painter", kegg_collections)
  }
  if (isTRUE(run_reports)) {
    add_stage("root_html_report", registry$analysis_collection)
    add_stage("manifest_validation", registry$analysis_collection)
  }

  do.call(rbind, stage_rows)
}

lisa_expected_artifact_row <- function(
  artifact_id, stage, expectation, status_source,
  analysis_id = "", contrast_id = "", collection = "",
  witnesses = character(), policy = "", reason = "",
  validation_scope = "canonical_material", external_contract = ""
) {
  expectation <- match.arg(
    expectation, c("required", "optional", "not_applicable")
  )
  data.frame(
    artifact_id = artifact_id,
    stage = stage,
    analysis_id = analysis_id,
    contrast_id = contrast_id,
    collection = collection,
    expectation = expectation,
    status_source = status_source,
    witness_rule = if (length(witnesses)) "all_material" else "none",
    witnesses = paste(witnesses, collapse = "|"),
    policy = policy,
    reason = reason,
    validation_scope = validation_scope,
    external_contract = external_contract,
    stringsAsFactors = FALSE
  )
}

lisa_pipeline_requested_report_mode <- function(
  output_dir, fallback = "standard"
) {
  fallback <- match.arg(fallback, c("standard", "full"))
  path <- file.path(output_dir, "contract_manifest.tsv")
  if (!file.exists(path)) return(fallback)
  contract <- read_lisa_tsv(path)
  if (!all(c("key", "value") %in% names(contract))) {
    stop(
      "LISA-ARTIFACT-002 contract_manifest.tsv lacks key/value columns.",
      call. = FALSE
    )
  }
  value <- as.character(contract$value[contract$key == "report_mode"])
  if (length(value) != 1L || !value %in% c("standard", "full")) {
    stop(
      "LISA-ARTIFACT-002 contract_manifest.tsv must declare exactly one ",
      "report_mode equal to standard or full.",
      call. = FALSE
    )
  }
  value
}

lisa_expected_artifacts <- function(
  de_index, contrast_index, registry, output_dir,
  file_label_prefix = "semantic", run_reports = TRUE, run_ora = FALSE,
  plot_formats = "png", export_formats = "tsv", source_data = TRUE,
  recipes = FALSE, requested_report_mode = "standard", category_evidence = FALSE,
  category_nes_variants = c("clean", "percentages", "direction", "dispersion"),
  integrated_full = FALSE
) {
  requested_report_mode <- match.arg(
    requested_report_mode, c("standard", "full")
  )
  prefix <- safe_file_label(file_label_prefix %||% "semantic")
  category_nes_variants <- lisa_nes_variants_validate(category_nes_variants)
  rows <- list()
  add <- function(...) {
    rows[[length(rows) + 1L]] <<- lisa_expected_artifact_row(...)
  }

  policy_text <- paste0(
    "plot_formats=", paste(sort(unique(plot_formats)), collapse = ","),
    ";export_formats=", paste(sort(unique(export_formats)), collapse = ","),
    ";source_data=", tolower(as.character(isTRUE(source_data))),
    ";recipes=", tolower(as.character(isTRUE(recipes)))
  )
  add(
    "contract:report_output_policy", "report_output_policy", "required",
    "contract", witnesses = "file:report_output_policy.tsv",
    policy = policy_text,
    reason = "The canonical output policy must be materialized before computation."
  )

  registry_capabilities <- lisa_collection_registry(
    registry$analysis_collection
  )
  for (i in seq_len(nrow(de_index))) {
    analysis_id <- as.character(de_index$analysis_id[[i]])
    analysis_file_id <- safe_file_label(analysis_id)
    for (j in seq_len(nrow(registry))) {
      collection <- as.character(registry$analysis_collection[[j]])
      collection_dir <- file.path(
        "outputs", "single_de", analysis_id,
        paste0("collection_", collection)
      )
      gsea_witnesses <- c(
        paste0("dir:", collection_dir),
        paste0(
          "file:", file.path(
            collection_dir, "enrichment",
            paste0(analysis_file_id, "_GSEA_", prefix, "_annotated.tsv")
          )
        ),
        paste0(
          "file:", file.path(
            collection_dir, "lisa_tables",
            paste0(
              analysis_file_id, "_", prefix,
              "_GSEA_category_summary.tsv"
            )
          )
        ),
        paste0(
          "file:", file.path(
            collection_dir, "qc",
            paste0(analysis_file_id, "_GSEA_universe_ledger.tsv")
          )
        )
      )
      if (collection %in% c("GOBP-C2", "GOMF", "GOCC")) {
        gsea_witnesses <- c(
          gsea_witnesses,
          paste0(
            "file:", file.path(
              collection_dir, "lisa_tables",
              paste0(
                analysis_file_id, "_", prefix,
                "_GSEA_OTHER_UNCLASSIFIED.tsv"
              )
            )
          )
        )
      }
      if (!identical(collection, "HALLMARKS")) gsea_witnesses <- c(gsea_witnesses,
        lisa_nes_expected_witnesses(collection_dir, analysis_file_id, prefix,
          category_nes_variants, plot_formats))
      add(
        paste("single_de_lisa", analysis_id, collection, sep = ":"),
        "single_de_lisa", "required", "single_de_status",
        analysis_id = analysis_id, collection = collection,
        witnesses = gsea_witnesses,
        reason = paste0(
          "Every configured analysis/collection requires a completed GSEA ",
          "run and its canonical LISA summary."
        )
      )

      capability <- registry_capabilities[
        registry_capabilities$analysis_collection == collection,
        , drop = FALSE
      ]
      if (!identical(collection, "HALLMARKS")) add(
        paste("single_de_category_inference", analysis_id, collection, sep = ":"),
        "single_de_category_inference", "required", "post_lisa_status",
        analysis_id = analysis_id, collection = collection,
        witnesses = paste0("file:", file.path(collection_dir, "lisa_tables",
          paste0(analysis_file_id, "_LISA_category_inference.tsv"))),
        reason = "Simultaneous category inference across applicable collections; HALLMARKS excluded.")
      if (isTRUE(category_evidence) && !identical(collection, "HALLMARKS")) {
        evidence_files <- "index.html"
        # Interactive evidence remains required when every static format is
        # disabled. Its renderer deliberately omits all member-figure files.
        if (length(intersect(plot_formats, c("png", "pdf", "svg")))) {
          evidence_files <- c(evidence_files, "figures/category_members_index.tsv",
            "figures/reproduce_category_member_evidence.R")
        }
        add(paste("single_de_category_evidence", analysis_id, collection, sep = ":"),
          "single_de_category_evidence", "required", "post_lisa_status",
          analysis_id = analysis_id, collection = collection,
          witnesses = paste0("file:", file.path("report_pages", "evidence",
            analysis_id, collection, evidence_files)),
          reason = "Category evidence is the primary gene-set/gene interpretation product.")
        add(paste("single_de_category_navigation", analysis_id, collection, sep = ":"),
          "single_de_category_navigation", "required", "post_lisa_status",
          analysis_id = analysis_id, collection = collection,
          witnesses = paste0("file:", file.path("report_pages", "category_navigation",
            analysis_id, collection, "index.html")),
          reason = "Individual category overview links to exact scoped evidence.")
      }
      ora_supported <- nrow(capability) == 1L &&
        isTRUE(capability$run_ora[[1L]])
      ora_expectation <- if (!ora_supported) {
        "not_applicable"
      } else if (isTRUE(run_ora)) {
        "required"
      } else {
        "optional"
      }
      ora_witnesses <- c(
        paste0(
          "file:", file.path(
            collection_dir, "enrichment",
            paste0(analysis_file_id, "_ORA_all_tested_annotated.tsv")
          )
        ),
        paste0(
          "file:", file.path(
            collection_dir, "enrichment",
            paste0(analysis_file_id, "_ORA_", prefix, "_annotated.tsv")
          )
        ),
        paste0(
          "file:", file.path(
            collection_dir, "lisa_tables",
            paste0(
              analysis_file_id, "_", prefix,
              "_ORA_category_summary.tsv"
            )
          )
        )
      )
      add(
        paste("single_de_ora", analysis_id, collection, sep = ":"),
        "single_de_ora", ora_expectation,
        if (identical(ora_expectation, "required")) {
          "single_de_status"
        } else {
          "none"
        },
        analysis_id = analysis_id, collection = collection,
        witnesses = if (identical(ora_expectation, "required")) {
          ora_witnesses
        } else {
          character()
        },
        policy = paste0("run_ora=", tolower(as.character(isTRUE(run_ora)))),
        reason = if (!ora_supported) {
          "ORA is not applicable to this collection registry entry."
        } else if (isTRUE(run_ora)) {
          "ORA was explicitly requested and all tested, selected and LISA tables are required."
        } else {
          "ORA is supported but was not requested."
        }
      )
    }
  }

  if (isTRUE(category_evidence) && any(registry$analysis_collection != "HALLMARKS")) {
    add("single_de_gene_evidence:global", "single_de_gene_evidence", "required", "post_lisa_status",
      witnesses = "file:report_pages/gene_evidence/index.html",
      reason = "Global gene evidence includes full memberships and individual-analysis statistics, not contrasts.")
  }

  if (is.null(contrast_index) || !nrow(contrast_index)) {
    add(
      "category_contrasts:none", "category_contrasts",
      "not_applicable", "none",
      reason = "No category-profile contrasts were configured."
    )
  } else {
    for (i in seq_len(nrow(contrast_index))) {
      contrast_row <- contrast_index[i, , drop = FALSE]
      contrast_id <- as.character(contrast_row$contrast_id[[1L]])
      output_id <- lisa_contrast_output_id(contrast_row)
      contrast_name <- lisa_contrast_name(contrast_row)
      for (j in seq_len(nrow(registry))) {
        collection <- as.character(registry$analysis_collection[[j]])
        supported <- isTRUE(registry$run_contrasts[[j]])
        collection_dir <- file.path(
          "outputs", "category_contrasts", contrast_name,
          paste0("collection_", collection)
        )
        witnesses <- c(
          paste0("dir:", collection_dir),
          paste0(
            "file:", file.path(
              collection_dir, "lisa_tables",
              paste0(
                safe_file_label(output_id), "_", prefix,
                "_GSEA_category_contrast.tsv"
              )
            )
          )
        )
        if (supported && !identical(collection, "HALLMARKS")) witnesses <- c(witnesses,
          lisa_nes_expected_witnesses(collection_dir, safe_file_label(output_id), prefix,
            category_nes_variants, plot_formats, contrast = TRUE))
        if (supported && isTRUE(category_evidence) && !identical(collection, "HALLMARKS")) {
          evidence_dir <- file.path("report_pages", "contrast_evidence", contrast_name, collection)
          evidence_files <- c("index.html", file.path("tables", paste0(c(
            "categories", "sets", "genes", "leading_edges", "conservation",
            "metadata", "provenance", "figure_index"), ".tsv")),
            "assets/contrast-evidence.css", "assets/contrast-evidence.js")
          add(paste("contrast_category_evidence", contrast_id, collection, sep = ":"),
            "contrast_category_evidence", "required", "post_lisa_status",
            contrast_id = contrast_id, collection = collection,
            witnesses = paste0("file:", file.path(evidence_dir, evidence_files)),
            reason = "A/B evidence requires exact scoped individual evidence and original contrast endpoints; missing support is not zero.")
        }
        add(
          paste("category_contrasts", contrast_id, collection, sep = ":"),
          "category_contrasts",
          if (supported) "required" else "not_applicable",
          if (supported) "contrast_status" else "none",
          contrast_id = contrast_id, collection = collection,
          witnesses = if (supported) witnesses else character(),
          reason = if (supported) {
            "The configured contrast is supported by this collection and requires its canonical contrast table."
          } else {
            "Category contrasts are not applicable to this collection registry entry."
          }
        )
      }
    }
  }

  add(
    "report:root_html_report", "root_html_report",
    if (isTRUE(run_reports)) "required" else "optional",
    if (isTRUE(run_reports)) "post_lisa_status" else "none",
    witnesses = if (isTRUE(run_reports)) {
      "file:report_index.html"
    } else {
      character()
    },
    policy = paste0("run_reports=", tolower(as.character(isTRUE(run_reports)))),
    reason = if (isTRUE(run_reports)) {
      "The canonical run requires a non-empty root HTML report."
    } else {
      "The root HTML report is supported but was not requested."
    }
  )

  full_requested <- identical(requested_report_mode, "full")
  if (full_requested && isTRUE(integrated_full)) {
    # run_lisa(mode = "full"): the FULL products are rendered inside this run
    # (artifacts/ and outputs/gene_level/) before its single report and its
    # single promotion, so they are canonical material of this manifest.
    add(
      "report:full_products", "full_report_products", "required",
      "contract",
      witnesses = c("file:contract_manifest.tsv", "file:product_plan.tsv",
        "file:artifacts/full_products_inventory.tsv"),
      policy = "requested_report_mode=full",
      reason = paste0(
        "Full products are rendered into artifacts/ of this same run and ",
        "integrated by the one report generator before promotion."
      )
    )
    return(do.call(rbind, rows))
  }
  add(
    "extension:full_report_contract", "full_report_extension_contract",
    if (full_requested) "required" else "optional",
    if (full_requested) "extension_contract" else "none",
    witnesses = if (full_requested) {
      c("file:contract_manifest.tsv", "file:product_plan.tsv")
    } else {
      character()
    },
    policy = paste0("requested_report_mode=", requested_report_mode),
    reason = if (full_requested) {
      paste0(
        "Full expansion is required by configuration and is delegated to the ",
        "separately validated post-promotion extension contract."
      )
    } else {
      "Full expansion is supported but was not requested."
    },
    validation_scope = "external_extension",
    external_contract = paste(
      c(
        "extension_receipt.tsv", "extension_inventory.tsv",
        "extension_manifest.tsv"
      ),
      collapse = "|"
    )
  )

  do.call(rbind, rows)
}

# Supplemental views have their own material witnesses, never canonical summary columns.
lisa_nes_expected_witnesses <- function(collection_dir, comparison, prefix,
    variants, formats, contrast = FALSE) {
  status <- file.path(collection_dir, "qc", paste0(comparison, "_category_nes_variants.tsv"))
  formats <- intersect(formats, c("png", "pdf", "svg"))
  paths <- status
  if (length(formats)) {
    stem <- lisa_nes_runtime_stem(comparison, prefix, contrast)
    stems <- if (contrast) c(stem, paste0(stem, c("_same_direction", "_opposite_direction"))) else stem
    files <- c(paste0(stem, "_recipe.R"), unlist(lapply(stems, function(selected_stem) c(
      paste0(selected_stem, "_source.tsv"),
      paste0(selected_stem, "_", variants, "_settings.json"),
      unlist(lapply(variants, function(v) paste0(selected_stem, "_", v, ".", formats)), use.names = FALSE))), use.names = FALSE))
    paths <- c(paths, file.path(collection_dir, "plots", "category_nes", files))
  }
  paste0("file:", paths)
}

lisa_artifact_status_rows <- function(
  expected_row, single_status, contrast_status, post_lisa_status,
  output_dir
) {
  source <- as.character(expected_row$status_source[[1L]])
  if (identical(source, "none")) return(data.frame())
  if (identical(source, "contract")) {
    return(data.frame(status = "completed", stringsAsFactors = FALSE))
  }
  if (identical(source, "extension_contract")) {
    observed_mode <- lisa_pipeline_requested_report_mode(
      output_dir, fallback = "standard"
    )
    return(data.frame(
      status = if (identical(observed_mode, "full")) {
        "completed"
      } else {
        "missing_contract"
      },
      stringsAsFactors = FALSE
    ))
  }
  status <- switch(
    source,
    single_de_status = single_status,
    contrast_status = contrast_status,
    post_lisa_status = post_lisa_status,
    stop("LISA-ARTIFACT-003 unknown status source: ", source, call. = FALSE)
  )
  if (is.null(status) || !is.data.frame(status) || !nrow(status) ||
      !"status" %in% names(status)) {
    return(data.frame())
  }
  match_fields <- switch(
    source,
    single_de_status = c("analysis_id", "collection"),
    contrast_status = c("contrast_id", "collection"),
    post_lisa_status = c(
      "stage", "analysis_id", "contrast_id", "collection"
    )
  )
  keep <- rep(TRUE, nrow(status))
  for (field in match_fields) {
    expected_value <- as.character(expected_row[[field]][[1L]])
    if (!nzchar(expected_value)) next
    if (!field %in% names(status)) return(data.frame())
    keep <- keep & !is.na(status[[field]]) &
      as.character(status[[field]]) == expected_value
  }
  status[keep, , drop = FALSE]
}

lisa_artifact_witness_validation <- function(expected_row, output_dir) {
  encoded <- as.character(expected_row$witnesses[[1L]])
  if (!nzchar(encoded)) {
    return(list(expected = 0L, material = 0L, missing = character()))
  }
  witnesses <- strsplit(encoded, "|", fixed = TRUE)[[1L]]
  material <- logical(length(witnesses))
  for (i in seq_along(witnesses)) {
    witness <- witnesses[[i]]
    if (startsWith(witness, "file:")) {
      relative <- substring(witness, 6L)
      path <- file.path(output_dir, relative)
      info <- file.info(path)
      material[[i]] <- file.exists(path) && !dir.exists(path) &&
        nrow(info) == 1L && is.finite(info$size[[1L]]) && info$size[[1L]] > 0
    } else if (startsWith(witness, "dir:")) {
      relative <- substring(witness, 5L)
      material[[i]] <- dir.exists(file.path(output_dir, relative))
    } else {
      stop(
        "LISA-ARTIFACT-004 unknown witness type in expected ledger: ",
        witness,
        call. = FALSE
      )
    }
  }
  list(
    expected = length(witnesses),
    material = sum(material),
    missing = witnesses[!material]
  )
}

lisa_reconcile_expected_artifacts <- function(
  expected, single_status = data.frame(), contrast_status = data.frame(),
  post_lisa_status = data.frame(), output_dir
) {
  required_columns <- c(
    "artifact_id", "stage", "analysis_id", "contrast_id", "collection",
    "expectation", "status_source", "witnesses", "validation_scope"
  )
  if (!is.data.frame(expected) || !nrow(expected) ||
      !all(required_columns %in% names(expected))) {
    stop(
      "LISA-ARTIFACT-005 expected artifact ledger is empty or malformed.",
      call. = FALSE
    )
  }
  rows <- lapply(seq_len(nrow(expected)), function(i) {
    row <- expected[i, , drop = FALSE]
    status_rows <- lisa_artifact_status_rows(
      row, single_status, contrast_status, post_lisa_status, output_dir
    )
    status_count <- nrow(status_rows)
    observed_status <- if (status_count == 1L) {
      tolower(as.character(status_rows$status[[1L]]))
    } else if (status_count > 1L) {
      "ambiguous"
    } else {
      "missing"
    }
    witness <- lisa_artifact_witness_validation(row, output_dir)
    expectation <- as.character(row$expectation[[1L]])
    validation_scope <- as.character(row$validation_scope[[1L]])
    validation_status <- if (identical(expectation, "required")) {
      if (status_count == 0L) {
        "missing_status"
      } else if (status_count > 1L) {
        "ambiguous_status"
      } else if (!identical(observed_status, "completed")) {
        "status_not_completed"
      } else if (witness$expected == 0L) {
        "missing_witness_contract"
      } else if (witness$material != witness$expected) {
        "missing_witness"
      } else if (identical(validation_scope, "external_extension")) {
        "external_contract_validated"
      } else {
        "validated"
      }
    } else if (identical(expectation, "optional")) {
      if (status_count == 0L) {
        "optional_not_requested"
      } else if (status_count == 1L && identical(observed_status, "completed") &&
                 (witness$expected == 0L || witness$material == witness$expected)) {
        "optional_completed"
      } else {
        "optional_incomplete"
      }
    } else {
      if (status_count == 0L && witness$material == 0L) {
        "not_applicable"
      } else {
        "not_applicable_but_present"
      }
    }
    data.frame(
      artifact_id = as.character(row$artifact_id[[1L]]),
      stage = as.character(row$stage[[1L]]),
      analysis_id = as.character(row$analysis_id[[1L]]),
      contrast_id = as.character(row$contrast_id[[1L]]),
      collection = as.character(row$collection[[1L]]),
      expectation = expectation,
      validation_scope = validation_scope,
      status_source = as.character(row$status_source[[1L]]),
      observed_status = observed_status,
      observed_status_rows = status_count,
      observed_message = if (status_count > 0L && "message" %in% names(status_rows)) {
        paste(unique(as.character(status_rows$message[!is.na(status_rows$message)])), collapse = " | ")
      } else "",
      witnesses_expected = witness$expected,
      witnesses_material = witness$material,
      missing_witnesses = paste(witness$missing, collapse = "|"),
      validation_status = validation_status,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lisa_assert_required_artifacts <- function(validation) {
  passed <- validation$validation_status == "validated" |
    (validation$validation_scope == "external_extension" &
      validation$validation_status == "external_contract_validated")
  failed <- validation[
    validation$expectation == "required" & !passed,
    , drop = FALSE
  ]
  if (!nrow(failed)) return(invisible(TRUE))
  detail <- vapply(seq_len(nrow(failed)), function(i) {
    row <- failed[i, , drop = FALSE]
    paste0(
      row$artifact_id[[1L]], "=", row$validation_status[[1L]],
      if ("observed_message" %in% names(row) && nzchar(row$observed_message[[1L]])) {
        paste0(": ", row$observed_message[[1L]])
      } else "",
      if (nzchar(row$missing_witnesses[[1L]])) {
        paste0(" (", row$missing_witnesses[[1L]], ")")
      } else {
        ""
      }
    )
  }, character(1))
  condition <- structure(
    list(
      message = paste0(
        "LISA-ARTIFACT-001 required pipeline artifacts did not validate; ",
        "the run cannot be promoted: ", paste(detail, collapse = " | ")
      ),
      call = NULL,
      validation = failed
    ),
    class = c(
      "lisa_artifact_validation_error", "lisa_error", "error", "condition"
    )
  )
  stop(condition)
}

run_lisa_pipeline <- function(
  de_index,
  contrast_index = NULL,
  dictionary_dir,
  term2gene,
  dictionary_path = NULL,
  category_map_path = NULL,
  output_dir,
  collections = lisa_default_collections(),
  base_dir = getwd(),
  lisa_project_root = NULL,
  lisa_dictionary = "core",
  file_label_prefix = "semantic",
  plot_formats = c("png"),
  export_formats = c("tsv", "xlsx"),
  run_gene_level = FALSE,
  run_reports = TRUE,
  run_kegg_maps = FALSE,
  kegg_access_mode = getOption("lisaR.kegg_access_mode", Sys.getenv("LISAR_KEGG_ACCESS_MODE", "external")),
  kegg_cache_root = getOption("lisaR.kegg_cache_root", Sys.getenv("LISAR_KEGG_CACHE", "")),
  kegg_snapshot_id = getOption("lisaR.kegg_snapshot_id", Sys.getenv("LISAR_KEGG_SNAPSHOT", "")),
  run_hallmarks = TRUE,
  run_ora = FALSE,
  gsea_padj_cutoff = 0.25,
  run_kegg = NULL,
  report_title = "LISA report",
  report_mode = c("standard", "full"),
  source_data = TRUE,
  recipes = FALSE,
  registered_category_map = FALSE,
  workers = 4L,
  dry_run = TRUE,
  .test_package_dir = NULL,
  category_evidence = TRUE,
  evidence_max_sets = 25L, evidence_max_genes = 40L,
  cache_dir = NULL, cache_mode = "off", cache_max_bytes = 512 * 1024^2,
  category_nes_variants = c("clean", "percentages", "direction", "dispersion"),
  .presentation = NULL
) {
  run_gene_level <- lisa_config_bool(run_gene_level, "run_gene_level")
  if (run_gene_level) {
    stop(
      "LISA-GENE-LEVEL-001 `run_gene_level = TRUE` is not executed by the canonical pipeline. Use `plan_lisa_extension()` and `render_lisa_categories()` after the canonical run completes.",
      call. = FALSE
    )
  }
  workers <- lisa_validate_workers(workers)
  category_nes_variants <- lisa_nes_variants_validate(category_nes_variants)
  gsea_padj_cutoff <- lisa_validate_gsea_padj_cutoff(gsea_padj_cutoff)
  report_mode <- match.arg(report_mode)
  source_data <- lisa_report_bool(source_data, "report.source_data", TRUE)
  recipes <- lisa_report_bool(recipes, "report.recipes", FALSE)
  if (!identical(report_mode, "standard")) {
    stop("LISA-REPORT-MODE-001 selected/full products require a completed canonical run. Use run_lisa_pipeline_from_config() for full mode, or plan_lisa_extension() and render_lisa_categories() explicitly.", call. = FALSE)
  }
  # The scientific computation is always the standard one. `.presentation` is
  # the internal route by which run_lisa(mode = "full") and KEGG maps render
  # their presentation products inside this same run, before its one report.
  presentation_mode <- if (is.list(.presentation)) {
    match.arg(.presentation$report_mode, c("standard", "full"))
  } else "standard"
  run_root <- lisa_run_root(output_dir)
  old_run_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = run_root)
  on.exit(options(lisaR.run_root = old_run_root), add = TRUE)
  if (!is.null(run_kegg)) warning("run_kegg is deprecated; migrate to run_kegg_maps. It is not interpreted automatically.", call. = FALSE)
  kegg_access_mode <- match.arg(kegg_access_mode, c("external", "cache_only"))
  if (isTRUE(run_kegg_maps) && identical(kegg_access_mode, "external")) {
    lisa_require_optional("KEGGREST", "retrieving optional KEGG pathway maps")
    lisa_require_optional("png", "rendering optional KEGG pathway maps")
    if (!nzchar(kegg_cache_root)) kegg_cache_root <- tempfile("lisaR-kegg-runtime-")
    if (!nzchar(kegg_snapshot_id)) {
      kegg_snapshot_id <- paste0("runtime-", format(Sys.time(), "%Y%m%dT%H%M%S"))
    }
  }
  if (isTRUE(run_kegg_maps) && identical(kegg_access_mode, "cache_only") &&
      (!nzchar(kegg_cache_root) || !nzchar(kegg_snapshot_id))) {
    stop("LISA-KEGG-010 cache_only KEGG maps require an explicit cache root and immutable snapshot ID.", call. = FALSE)
  }
  if (is.character(de_index) && length(de_index) == 1) {
    base_dir <- dirname(de_index)
  }
  de_index_valid <- validate_lisa_de_index(de_index, base_dir = base_dir, require_files = !isTRUE(dry_run))
  if (!"species" %in% names(de_index_valid) || any(!nzchar(as.character(de_index_valid$species)))) {
    stop("LISA-SPECIES-004 every analysis must declare species. Repair: add a species column with Homo sapiens or Mus musculus.", call. = FALSE)
  }
  species_contracts <- do.call(rbind, lapply(seq_len(nrow(de_index_valid)), function(i) {
    contract <- lisa_species_contract(de_index_valid$species[[i]],
      msigdb_mode = if ("msigdb_mode" %in% names(de_index_valid)) de_index_valid$msigdb_mode[[i]] else NULL,
      term2gene_target_species = if ("term2gene_target_species" %in% names(de_index_valid)) de_index_valid$term2gene_target_species[[i]] else NULL)
    data.frame(analysis_id = de_index_valid$analysis_id[[i]], as.data.frame(contract, stringsAsFactors = FALSE), stringsAsFactors = FALSE)
  }))
  contrast_index_valid <- validate_lisa_contrast_index(contrast_index, de_index = de_index_valid)
  registry <- lisa_collection_registry(collections)
  if (!isTRUE(run_hallmarks)) registry <- registry[registry$analysis_collection != "HALLMARKS", , drop = FALSE]
  registry$run_ora <- registry$run_ora & isTRUE(run_ora)
  plan <- lisa_pipeline_plan(
    de_index_valid,
    contrast_index_valid,
    collections = registry$analysis_collection,
    run_gene_level = run_gene_level,
    category_evidence = category_evidence,
    run_reports = run_reports,
    run_kegg_maps = run_kegg_maps,
    run_hallmarks = run_hallmarks,
    run_ora = run_ora
  )

  lisa_guarded_dir_create(output_dir, run_root)
  outputs_root <- file.path(output_dir, "outputs")
  lisa_guarded_dir_create(outputs_root, run_root)
  lisa_write_pipeline_config(output_dir, de_index_valid, contrast_index_valid)
  write_lisa_tsv(species_contracts, file.path(output_dir, "species_contracts.tsv"))
  write_lisa_tsv(plan, file.path(output_dir, "lisa_pipeline_plan.tsv"))
  write_lisa_tsv(registry, file.path(output_dir, "lisa_collection_registry.tsv"))
  output_policy <- lisa_output_policy_table(plot_formats, source_data, recipes)
  write_lisa_tsv(output_policy, file.path(output_dir, "report_output_policy.tsv"))
  requested_report_mode <- lisa_pipeline_requested_report_mode(
    output_dir, fallback = report_mode
  )
  expected_artifacts <- lisa_expected_artifacts(
    de_index = de_index_valid,
    contrast_index = contrast_index_valid,
    registry = registry,
    output_dir = output_dir,
    file_label_prefix = file_label_prefix,
    run_reports = run_reports,
    run_ora = run_ora,
    plot_formats = plot_formats,
    export_formats = export_formats,
    source_data = source_data,
    recipes = recipes,
    requested_report_mode = requested_report_mode,
    category_evidence = category_evidence,
    category_nes_variants = category_nes_variants,
    integrated_full = identical(presentation_mode, "full")
  )
  write_lisa_tsv(
    expected_artifacts, file.path(output_dir, "expected_artifacts.tsv")
  )
  artifact_validation <- data.frame()

  metadata <- data.frame(
    key = c("dictionary_dir", "term2gene", "dry_run", "workers", "gsea_padj_cutoff", "report_mode", "category_nes_variants", "package_stage", "collections", "dictionaries", "outputs_root"),
    value = c(
      dictionary_dir,
      term2gene,
      as.character(dry_run),
      as.character(workers),
      lisa_config_text(gsea_padj_cutoff),
      presentation_mode,
      paste(category_nes_variants, collapse = ","),
      "single_de_and_category_contrast_engine_wired",
      paste(registry$analysis_collection, collapse = ","),
      paste(registry$dictionary_id, collapse = ","),
      outputs_root
    ),
    stringsAsFactors = FALSE
  )
  write_lisa_tsv(metadata, file.path(output_dir, "lisa_pipeline_metadata.tsv"))

  if (!isTRUE(dry_run)) {
    single_status <- run_lisa_pipeline_single_de(
      de_index = de_index_valid,
      registry = registry,
      base_dir = base_dir,
      output_dir = outputs_root,
      term2gene = term2gene,
      dictionary_path = dictionary_path,
      category_map_path = category_map_path,
      lisa_project_root = lisa_project_root,
      lisa_dictionary = lisa_dictionary,
      file_label_prefix = file_label_prefix,
      plot_formats = plot_formats,
      export_formats = export_formats,
      report_mode = report_mode,
      source_data = source_data,
      recipes = recipes,
      registered_category_map = registered_category_map,
      gsea_padj_cutoff = gsea_padj_cutoff,
      cache_dir = cache_dir, cache_mode = cache_mode, cache_max_bytes = cache_max_bytes,
      category_evidence = category_evidence,
      category_nes_variants = category_nes_variants,
      workers = workers
    )
    write_lisa_tsv(single_status, file.path(output_dir, "single_de_status.tsv"))

    contrast_status <- run_lisa_pipeline_contrasts(
      de_index = de_index_valid,
      contrast_index = contrast_index_valid,
      registry = registry,
      output_dir = outputs_root,
      single_status = single_status,
      file_label_prefix = file_label_prefix,
      plot_formats = plot_formats,
      export_formats = export_formats,
      source_data = source_data,
      recipes = recipes,
      gsea_padj_cutoff = gsea_padj_cutoff,
      category_nes_variants = category_nes_variants,
      workers = workers
    )
    write_lisa_tsv(contrast_status, file.path(output_dir, "contrast_status.tsv"))

    post_lisa_status <- run_lisa_pipeline_post_lisa(
      de_index = de_index_valid,
      contrast_index = contrast_index_valid,
      registry = registry,
      output_dir = output_dir,
      outputs_root = outputs_root,
      single_status = single_status,
      contrast_status = contrast_status,
      file_label_prefix = file_label_prefix,
      run_gene_level = run_gene_level,
      run_reports = run_reports,
      run_kegg = run_kegg_maps,
      kegg_access_mode = kegg_access_mode,
      kegg_cache_root = kegg_cache_root,
      kegg_snapshot_id = kegg_snapshot_id,
      species_contracts = species_contracts,
      report_title = report_title,
      workers = workers,
      .test_package_dir = .test_package_dir,
      category_evidence = category_evidence, evidence_max_sets = evidence_max_sets,
      evidence_max_genes = evidence_max_genes, lisa_dictionary = lisa_dictionary,
      gsea_padj_cutoff = gsea_padj_cutoff, plot_formats = plot_formats,
      source_data = source_data, recipes = recipes, term2gene = term2gene,
      category_nes_variants = category_nes_variants,
      .presentation = .presentation
    )
    write_lisa_tsv(post_lisa_status, file.path(output_dir, "post_lisa_status.tsv"))

    artifact_validation <- lisa_reconcile_expected_artifacts(
      expected_artifacts,
      single_status = single_status,
      contrast_status = contrast_status,
      post_lisa_status = post_lisa_status,
      output_dir = output_dir
    )
    write_lisa_tsv(
      artifact_validation, file.path(output_dir, "artifact_validation.tsv")
    )
    lisa_assert_required_artifacts(artifact_validation)
    lisa_assert_no_failed_stages(single_status, contrast_status, post_lisa_status)

    pending_stages <- setdiff(
      plan$stage,
      c(
        "single_de_lisa", "category_contrasts", "single_de_category_evidence", "single_de_category_inference",
        "single_de_category_navigation", "single_de_gene_evidence", "contrast_category_evidence",
        "single_de_category_count_heatmaps", "single_de_enrichmentmap",
        "single_de_gene_level_tables", "single_de_gene_cards",
        "single_de_volcano_overlays", "single_de_recurrent_gene_screen",
        "single_de_hallmarks_lollipop", "single_de_hallmarks_gene_cards",
        "single_de_hallmarks_volcano_overlays",
        "single_de_leading_edge_gene_heatmaps",
        "single_de_hallmarks_leading_edge_gene_heatmaps",
        "contrast_gene_level_tables", "contrast_category_cards",
        "contrast_paired_gene_heatmaps", "contrast_gene_category_networks",
        "kegg_pathway_layers", "kegg_pathway_painter",
        "root_html_report", "manifest_validation"
      )
    )
    if (length(pending_stages) > 0) {
      pending <- data.frame(
        stage = pending_stages,
        status = "pending_migration",
        note = "Post-LISA product generation is planned but not wired in lisaR yet.",
        stringsAsFactors = FALSE
      )
      write_lisa_tsv(pending, file.path(output_dir, "pending_migration_stages.tsv"))
    }
  }

  invisible(list(
    plan = plan,
    metadata = metadata,
    expected_artifacts = expected_artifacts,
    artifact_validation = artifact_validation
  ))
}

lisa_assert_no_failed_stages <- function(...) {
  failures <- unlist(lapply(list(...), function(status) {
    if (is.null(status) || nrow(status) == 0L || !"status" %in% names(status)) {
      return(character())
    }
    failed <- status[tolower(as.character(status$status)) == "failed", , drop = FALSE]
    if (nrow(failed) == 0L) return(character())
    vapply(seq_len(nrow(failed)), function(i) {
      identifiers <- intersect(
        c("stage", "analysis_id", "contrast_id", "collection"),
        names(failed)
      )
      detail <- paste(
        sprintf(
          "%s=%s",
          identifiers,
          as.character(failed[i, identifiers, drop = TRUE])
        ),
        collapse = ", "
      )
      message <- if ("message" %in% names(failed)) {
        as.character(failed$message[[i]])
      } else {
        ""
      }
      if (is.na(message)) message <- ""
      paste0(detail, if (nzchar(message)) paste0(": ", message) else "")
    }, character(1))
  }), use.names = FALSE)
  if (length(failures)) {
    stop(
      "LISA-RUN-FAILED: one or more pipeline stages failed; ",
      "the run cannot be promoted as PASS. ",
      paste(failures, collapse = " | "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

run_lisa_pipeline_single_de <- function(de_index, registry, base_dir, output_dir, term2gene, dictionary_path, category_map_path,
                                        lisa_project_root, lisa_dictionary, file_label_prefix,
                                        plot_formats, export_formats, report_mode = "standard",
                                        source_data = TRUE, recipes = FALSE,
                                        registered_category_map = FALSE,
                                        gsea_padj_cutoff = 0.25,
                                        workers = 4L, category_evidence = FALSE, cache_dir = NULL,
                                        cache_mode = "off", cache_max_bytes = 512 * 1024^2,
                                        category_nes_variants = c("clean", "percentages", "direction", "dispersion")) {
  # Collection tasks for one analysis share the same parent. Create those
  # parents serially before forking so workers never race while establishing a
  # managed-directory boundary.
  single_de_root <- file.path(output_dir, "single_de")
  lisa_guarded_dir_create(single_de_root, output_dir)
  for (analysis_id in unique(as.character(de_index$analysis_id))) {
    lisa_guarded_dir_create(
      file.path(single_de_root, analysis_id), output_dir
    )
  }
  # Establish the shared cache directory before collection workers start.
  cache_options <- lisa_content_cache_options(cache_dir, cache_mode, cache_max_bytes)
  if (cache_options$mode %in% c("readwrite", "refresh")) {
    lisa_cache_slot(cache_options, paste(rep("0", 64L), collapse = ""), create = TRUE)
  }
  tasks <- expand.grid(
    registry_row = seq_len(nrow(registry)),
    de_index_row = seq_len(nrow(de_index)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  task_ids <- sprintf(
    "single_de:%s:%s",
    de_index$analysis_id[tasks$de_index_row],
    registry$analysis_collection[tasks$registry_row]
  )
  values <- lisa_pipeline_map(task_ids, function(task_id) {
    task_index <- match(task_id, task_ids)
    i <- tasks$de_index_row[[task_index]]
    j <- tasks$registry_row[[task_index]]
    de_row <- de_index[i, , drop = FALSE]
    analysis_id <- de_row$analysis_id[[1]]
    de_path <- lisa_norm_path(de_row$de_path[[1]], base_dir)
    collection <- registry$analysis_collection[[j]]
    out_dir <- file.path(output_dir, "single_de", analysis_id, paste0("collection_", collection))
    status <- "completed"
    message <- ""
    started_at <- as.character(Sys.time())
    result <- tryCatch({
      run_LISA_DE(
        input = de_path,
        output_dir = out_dir,
        comparison_name = analysis_id,
        display_title = lisa_de_label(de_index, analysis_id),
        comparison_subtitle = if ("comparison" %in% names(de_row)) de_row$comparison[[1]] else analysis_id,
        positive_direction = if ("positive_direction" %in% names(de_row)) de_row$positive_direction[[1]] else "",
        species = de_row$species[[1]],
        msigdb_mode = lisa_species_contract(de_row$species[[1]], if ("msigdb_mode" %in% names(de_row)) de_row$msigdb_mode[[1]] else NULL)$msigdb_mode,
        symbol_col = if ("symbol_col" %in% names(de_row)) de_row$symbol_col[[1]] else NULL,
        rank_col = if ("rank_col" %in% names(de_row)) de_row$rank_col[[1]] else NULL,
        logfc_col = if ("logfc_col" %in% names(de_row)) de_row$logfc_col[[1]] else NULL,
        padj_col = if ("padj_col" %in% names(de_row)) de_row$padj_col[[1]] else NULL,
        pvalue_col = if ("pvalue_col" %in% names(de_row)) de_row$pvalue_col[[1]] else NULL,
        gene_id_col = if ("gene_id_col" %in% names(de_row)) de_row$gene_id_col[[1]] else NULL,
        lisa_dictionary = lisa_dictionary,
        lisa_dictionary_path = dictionary_path,
        term2gene_path = term2gene,
        category_map_path = category_map_path,
        universes = collection,
        gsea_padj_cutoff = gsea_padj_cutoff,
        run_ora = isTRUE(registry$run_ora[[j]]),
        plots = lisa_plots_for_collection(registry[j, , drop = FALSE], report_mode = report_mode,
          category_evidence = category_evidence),
        export_formats = export_formats,
        plot_formats = plot_formats,
        figure_source_data = source_data,
        figure_recipes = recipes,
        category_nes_variants = category_nes_variants,
        registered_category_map = registered_category_map,
        file_label_prefix = file_label_prefix,
        lisa_project_root = lisa_project_root,
        allow_missing_dictionary = isTRUE(registry$allow_missing_dictionary[[j]]),
        n_threads = 1L,
        cache_dir = cache_dir, cache_mode = cache_mode, cache_max_bytes = cache_max_bytes,
        verbose = TRUE
      )
    }, error = function(e) {
      status <<- "failed"
      message <<- conditionMessage(e)
      NULL
    })
    if (!is.null(result) && is.list(result) && identical(result$status, "pending_dictionary")) {
      status <- "pending_dictionary"
    }
    data.frame(
      analysis_id = analysis_id,
      collection = collection,
      dictionary_id = registry$dictionary_id[[j]],
      status = status,
      output_dir = out_dir,
      message = message,
      started_at = started_at,
      ended_at = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
  }, workers = workers)
  do.call(rbind, unname(values))
}

run_lisa_pipeline_contrasts <- function(de_index, contrast_index, registry, output_dir, single_status,
                                        file_label_prefix, plot_formats, export_formats,
                                        source_data = TRUE, recipes = FALSE,
                                        gsea_padj_cutoff = 0.25, workers = 4L,
                                        category_nes_variants = c("clean", "percentages", "direction", "dispersion")) {
  if (is.null(contrast_index) || nrow(contrast_index) == 0) {
    return(data.frame())
  }
  contrast_registry <- registry[registry$run_contrasts, , drop = FALSE]
  if (nrow(contrast_registry) == 0) {
    return(data.frame())
  }

  # As above, every collection for a contrast shares a parent directory.
  # Establish it in the parent process before multicore execution.
  contrast_root <- file.path(output_dir, "category_contrasts")
  lisa_guarded_dir_create(contrast_root, output_dir)
  contrast_names <- vapply(
    seq_len(nrow(contrast_index)),
    function(i) lisa_contrast_name(contrast_index[i, , drop = FALSE]),
    character(1)
  )
  for (contrast_name in unique(contrast_names)) {
    lisa_guarded_dir_create(
      file.path(contrast_root, contrast_name), output_dir
    )
  }

  tasks <- expand.grid(
    registry_row = seq_len(nrow(contrast_registry)),
    contrast_index_row = seq_len(nrow(contrast_index)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  task_ids <- sprintf(
    "contrast:%s:%s",
    contrast_index$contrast_id[tasks$contrast_index_row],
    contrast_registry$analysis_collection[tasks$registry_row]
  )
  values <- lisa_pipeline_map(task_ids, function(task_id) {
    task_index <- match(task_id, task_ids)
    i <- tasks$contrast_index_row[[task_index]]
    j <- tasks$registry_row[[task_index]]
    contrast_row <- contrast_index[i, , drop = FALSE]
    contrast_id <- contrast_row$contrast_id[[1]]
    output_id <- lisa_contrast_output_id(contrast_row)
    contrast_name <- lisa_contrast_name(contrast_row)
    a_id <- contrast_row$analysis_a[[1]]
    b_id <- contrast_row$analysis_b[[1]]
    label <- if ("contrast_label" %in% names(contrast_row)) contrast_row$contrast_label[[1]] else contrast_id
    label_a <- lisa_de_label(de_index, a_id)
    label_b <- lisa_de_label(de_index, b_id)
    collection <- contrast_registry$analysis_collection[[j]]
    out_dir <- file.path(output_dir, "category_contrasts", contrast_name, paste0("collection_", collection))
    status <- "completed"
    message <- ""
    started_at <- as.character(Sys.time())
    ready <- lisa_single_de_collection_ready(single_status, a_id, b_id, collection)
    if (!isTRUE(ready$ready)) {
      status <- "skipped_single_de_not_completed"
      message <- ready$message
    } else {
      tryCatch({
        run_LISA_contrast(
          contrast_a = file.path(output_dir, "single_de", a_id),
          contrast_b = file.path(output_dir, "single_de", b_id),
          output_dir = out_dir,
          comparison_name = output_id,
          contrast_a_label = label_a,
          contrast_b_label = label_b,
          contrast_a_title = label_a,
          contrast_b_title = label_b,
          comparison_subtitle = if ("comparison" %in% names(contrast_row)) contrast_row$comparison[[1]] else NULL,
          positive_direction = if ("positive_direction" %in% names(contrast_row)) contrast_row$positive_direction[[1]] else NULL,
          universes = collection,
          gsea_padj_cutoff = gsea_padj_cutoff,
          export_formats = export_formats,
          plot_formats = plot_formats,
          figure_source_data = source_data,
          figure_recipes = recipes,
          category_nes_variants = category_nes_variants,
          file_label_prefix = file_label_prefix,
          plot_title_prefix = label,
          verbose = TRUE
        )
      }, error = function(e) {
        status <<- "failed"
        message <<- conditionMessage(e)
      })
    }
    data.frame(
      contrast_id = contrast_id,
      output_id = output_id,
      contrast_name = contrast_name,
      collection = collection,
      dictionary_id = contrast_registry$dictionary_id[[j]],
      status = status,
      output_dir = out_dir,
      message = message,
      started_at = started_at,
      ended_at = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
  }, workers = workers)
  do.call(rbind, unname(values))
}

lisa_single_de_collection_ready <- function(single_status, analysis_a, analysis_b, collection) {
  if (is.null(single_status) || nrow(single_status) == 0) {
    return(list(ready = FALSE, message = "single_de_status is empty"))
  }
  needed <- single_status[
    single_status$analysis_id %in% c(analysis_a, analysis_b) &
      single_status$collection == collection,
    ,
    drop = FALSE
  ]
  if (nrow(needed) != 2) {
    return(list(
      ready = FALSE,
      message = sprintf("missing completed single-DE outputs for %s/%s collection %s", analysis_a, analysis_b, collection)
    ))
  }
  bad <- needed[needed$status != "completed", , drop = FALSE]
  if (nrow(bad) > 0) {
    return(list(
      ready = FALSE,
      message = sprintf(
        "single-DE outputs not completed for collection %s: %s",
        collection,
        paste(sprintf("%s=%s", bad$analysis_id, bad$status), collapse = ", ")
      )
    ))
  }
  list(ready = TRUE, message = "")
}

lisa_write_pipeline_config <- function(output_dir, de_index, contrast_index) {
  config_dir <- file.path(output_dir, "config")
  lisa_guarded_dir_create(config_dir)
  write_lisa_tsv(de_index, file.path(config_dir, "de_index.tsv"))
  if (is.null(contrast_index) || nrow(contrast_index) == 0) {
    write_lisa_tsv(data.frame(), file.path(config_dir, "contrast_index.tsv"))
    return(invisible(config_dir))
  }
  rows <- lapply(seq_len(nrow(contrast_index)), function(i) {
    row <- contrast_index[i, , drop = FALSE]
    data.frame(
      contrast_id = row$contrast_id[[1]],
      output_id = lisa_contrast_output_id(row),
      analysis_a = row$analysis_a[[1]],
      analysis_b = row$analysis_b[[1]],
      contrast_a = row$analysis_a[[1]],
      contrast_b = row$analysis_b[[1]],
      label_a = lisa_de_label(de_index, row$analysis_a[[1]]),
      label_b = lisa_de_label(de_index, row$analysis_b[[1]]),
      contrast_label = if ("contrast_label" %in% names(row)) row$contrast_label[[1]] else row$contrast_id[[1]],
      stringsAsFactors = FALSE
    )
  })
  write_lisa_tsv(do.call(rbind, rows), file.path(config_dir, "contrast_index.tsv"))
  invisible(config_dir)
}

lisa_de_label <- function(de_index, analysis_id) {
  row <- de_index[de_index$analysis_id == analysis_id, , drop = FALSE]
  for (col in c("label", "analysis_label", "title")) {
    if (nrow(row) == 1 && col %in% names(row) && nzchar(as.character(row[[col]][[1]]))) {
      return(as.character(row[[col]][[1]]))
    }
  }
  analysis_id
}

lisa_contrast_output_id <- function(row) {
  if ("output_id" %in% names(row) && nzchar(as.character(row$output_id[[1]]))) {
    return(lisa_safe_id(as.character(row$output_id[[1]]), "output_id"))
  }
  lisa_safe_id(as.character(row$contrast_id[[1]]), "contrast_id")
}

lisa_contrast_name <- function(row) {
  paste(lisa_safe_id(as.character(row$contrast_id[[1]]), "contrast_id"), lisa_contrast_output_id(row), sep = "_")
}

lisa_plots_for_collection <- function(registry_row, report_mode = "standard", category_evidence = FALSE) {
  report_mode <- match.arg(report_mode, c("standard", "full"))
  # Preserve category overview families independently of the new evidence view.
  # Detailed member-set panels remain optional; ORA is generated only on request.
  standard_plots <- if (identical(registry_row$post_lisa_profile[[1]], "hallmarks_minimal")) {
    "lollipop"
  } else c("lollipop", "direction_lollipop", "dotplot")
  if (isTRUE(registry_row$run_ora[[1]])) standard_plots <- c(standard_plots, "barplot")
  if (isTRUE(registry_row$run_category_pathway_plots[[1]]) && !isTRUE(category_evidence)) {
    standard_plots <- c(standard_plots, "category_pathways")
  }
  if (identical(report_mode, "standard")) return(standard_plots)
  if (identical(registry_row$post_lisa_profile[[1]], "hallmarks_minimal")) {
    return(c("lollipop"))
  }
  plots <- c("lollipop", "direction_lollipop", "barplot", "dotplot")
  if (isTRUE(registry_row$run_category_pathway_plots[[1]])) {
    plots <- c(plots, "category_pathways")
  }
  plots
}

lisa_post_lisa_planned_scripts <- function(
  contrast_index, registry, single_status, contrast_status,
  run_gene_level, run_reports, run_kegg, species_contracts, category_evidence = FALSE
) {
  scripts <- character()
  completed_single <- single_status[single_status$status == "completed", , drop = FALSE]
  for (i in seq_len(nrow(completed_single))) {
    row <- completed_single[i, , drop = FALSE]
    registry_row <- registry[
      registry$analysis_collection == row$collection[[1L]], , drop = FALSE
    ]
    if (!nrow(registry_row)) next
    if (isTRUE(registry_row$run_enrichmentmap[[1L]])) {
      scripts <- c(scripts, "build_single_de_enrichmentmap.R")
    }
    if (isTRUE(run_gene_level)) {
      if (isTRUE(registry_row$run_gene_cards[[1L]])) {
        scripts <- c(scripts, "build_single_de_category_gene_cards.R")
      }
      if (isTRUE(registry_row$run_volcano_overlays[[1L]])) {
        scripts <- c(scripts, "build_single_de_category_volcano_overlays.R")
      }
      if (isTRUE(registry_row$run_recurrent_gene_screen[[1L]])) {
        scripts <- c(scripts, "build_single_de_recurrent_gene_screen.R")
      }
      scripts <- c(scripts, "build_single_de_leading_edge_gene_heatmaps.R")
      if (isTRUE(run_kegg) && isTRUE(registry_row$run_kegg_layers[[1L]])) {
        scripts <- c(scripts, "build_single_de_kegg_pathway_painter.R")
      }
    }
  }

  completed_contrasts <- contrast_status[
    contrast_status$status == "completed", , drop = FALSE
  ]
  if (isTRUE(run_gene_level) && nrow(completed_contrasts)) {
    scripts <- c(
      scripts, "build_contrast_gene_level_product.R",
      "build_contrast_gene_category_network.R"
    )
    if (isTRUE(run_kegg)) {
      for (i in seq_len(nrow(completed_contrasts))) {
        row <- completed_contrasts[i, , drop = FALSE]
        registry_row <- registry[
          registry$analysis_collection == row$collection[[1L]], , drop = FALSE
        ]
        if (!nrow(registry_row) || !isTRUE(registry_row$run_kegg_layers[[1L]])) next
        contrast_def <- contrast_index[
          contrast_index$contrast_id == row$contrast_id[[1L]], , drop = FALSE
        ]
        if (length(lisa_contrast_species(contrast_def, species_contracts)) == 1L) {
          scripts <- c(scripts, "build_contrast_kegg_pathway_painter.R")
        }
      }
    }
  }
  if (isTRUE(category_evidence) && any(completed_single$collection != "HALLMARKS")) {
    scripts <- c(scripts, "build_category_evidence.R", "build_category_navigation.R", "build_gene_evidence.R")
    if (any(completed_contrasts$collection != "HALLMARKS"))
      scripts <- c(scripts, "build_contrast_evidence.R", "build_contrast_navigation.R")
  }
  if (isTRUE(run_reports)) scripts <- c(scripts, "build_LISA_report.R")
  sort(unique(scripts))
}

run_lisa_pipeline_post_lisa <- function(de_index, contrast_index, registry, output_dir, outputs_root,
                                        single_status, contrast_status, file_label_prefix,
                                        run_gene_level, run_reports, run_kegg,
                                        report_title = "LISA report", species_contracts = data.frame(),
                                        kegg_access_mode = "external",
                                        kegg_cache_root = "", kegg_snapshot_id = "",
                                        workers = 4L, .test_package_dir = NULL,
                                        category_evidence = FALSE, evidence_max_sets = 25L,
                                        evidence_max_genes = 40L, lisa_dictionary = "core",
                                        gsea_padj_cutoff = 0.25, plot_formats = "png",
                                        source_data = TRUE, recipes = FALSE, term2gene = NULL,
                                        category_nes_variants = c("clean", "percentages", "direction", "dispersion"),
                                        .presentation = NULL) {
  rows <- list()
  add_row <- function(stage, analysis_id = "", contrast_id = "", collection = "", status = "completed",
                      output_path = "", message = "") {
    rows[[length(rows) + 1L]] <<- lisa_post_status_row(
      stage, analysis_id, contrast_id, collection, status, output_path, message
    )
    try(write_lisa_tsv(do.call(rbind, rows), file.path(output_dir, "post_lisa_status.tsv")), silent = TRUE)
  }

  inference_status <- lisa_run_category_inference(output_dir, term2gene,
    single_status = single_status, plot_formats = plot_formats)
  if (nrow(inference_status)) rows[[length(rows) + 1L]] <- inference_status

  package_dir <- lisa_resolve_package_dir(
    .test_package_dir = .test_package_dir
  )
  trusted_run_root <- lisa_assert_run_tree_safe(output_dir)
  planned_scripts <- lisa_post_lisa_planned_scripts(
    contrast_index, registry, single_status, contrast_status,
    run_gene_level, run_reports, run_kegg, species_contracts, category_evidence
  )
  code_ledger <- lisa_write_code_identity_ledger(
    package_dir, planned_scripts, file.path(output_dir, "code_identity.tsv")
  )
  completed_single <- single_status[single_status$status == "completed", , drop = FALSE]
  single_ids <- sprintf("post_single:%s:%s", completed_single$analysis_id, completed_single$collection)
  single_values <- lisa_pipeline_map(single_ids, function(task_id) {
    i <- match(task_id, single_ids)
    lisa_run_single_post_unit(
      completed_single[i, , drop = FALSE], registry, output_dir, outputs_root,
      file_label_prefix, package_dir, run_gene_level, trusted_run_root,
      code_ledger, category_evidence = category_evidence,
      evidence_max_sets = evidence_max_sets, evidence_max_genes = evidence_max_genes,
      lisa_dictionary = lisa_dictionary, gsea_padj_cutoff = gsea_padj_cutoff,
      plot_formats = plot_formats, source_data = source_data, recipes = recipes,
      de_index = de_index
    )
  }, workers = workers)
  for (value in unname(single_values)) {
    if (nrow(value)) rows[[length(rows) + 1L]] <- value
  }

  if (isTRUE(category_evidence) && any(completed_single$collection != "HALLMARKS")) {
    if (!is.character(term2gene) || length(term2gene) != 1L || !file.exists(term2gene)) {
      stop("LISA-GENE-EVIDENCE-INPUT: full TERM2GENE is required for gene membership evidence.", call. = FALSE)
    }
    evidence_ids <- unique(completed_single$analysis_id[completed_single$collection != "HALLMARKS"])
    res <- lisa_run_post_script(package_dir, "build_gene_evidence.R",
      c("--project-dir", output_dir, "--analysis-ids", paste(evidence_ids, collapse = ","),
        "--term2gene", term2gene, "--tier", lisa_dictionary,
        "--gsea-padj-cutoff", format(gsea_padj_cutoff, digits = 17),
        "--output-dir", file.path(output_dir, "report_pages", "gene_evidence")),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger)
    add_row("single_de_gene_evidence", status = res$status, output_path = res$output_path, message = res$message)
  }

  completed_contrasts <- contrast_status[contrast_status$status == "completed", , drop = FALSE]
  # All individual evidence builders above have joined before any A/B builder
  # starts. A failed side is explicit and never replaced by the other side.
  evidence_contrasts <- completed_contrasts[completed_contrasts$collection != "HALLMARKS", , drop = FALSE]
  if (isTRUE(category_evidence) && nrow(evidence_contrasts)) {
    evidence_status <- if (length(rows)) do.call(rbind, rows) else data.frame()
    contrast_ids <- sprintf("post_contrast_evidence:%s:%s", evidence_contrasts$contrast_id, evidence_contrasts$collection)
    values <- lisa_pipeline_map(contrast_ids, function(task_id) {
      i <- match(task_id, contrast_ids)
      lisa_run_contrast_evidence_unit(evidence_contrasts[i, , drop = FALSE],
        contrast_index, evidence_status, output_dir, package_dir, trusted_run_root,
        code_ledger, lisa_dictionary, evidence_max_sets, evidence_max_genes,
        plot_formats, category_nes_variants, gsea_padj_cutoff)
    }, workers = workers)
    for (value in unname(values)) if (nrow(value)) rows[[length(rows) + 1L]] <- value
  }
  if (isTRUE(run_gene_level)) {
    contrast_ids <- sprintf("post_contrast:%s:%s", completed_contrasts$contrast_id, completed_contrasts$collection)
    contrast_values <- lisa_pipeline_map(contrast_ids, function(task_id) {
      i <- match(task_id, contrast_ids)
      lisa_run_contrast_post_unit(
        completed_contrasts[i, , drop = FALSE], output_dir, package_dir,
        trusted_run_root, code_ledger
      )
    }, workers = workers)
    for (value in unname(contrast_values)) {
      if (nrow(value)) rows[[length(rows) + 1L]] <- value
    }
  }

  # KEGG cache population is deliberately serial: immutable external snapshots
  # use check-then-store semantics and must never be written by competing tasks.
  if (isTRUE(run_gene_level) && isTRUE(run_kegg)) {
    for (i in seq_len(nrow(completed_single))) {
      row <- completed_single[i, , drop = FALSE]
      registry_row <- registry[registry$analysis_collection == row$collection[[1]], , drop = FALSE]
      built <- do.call(rbind, rows)
      ready <- built$stage == "single_de_gene_level_tables" & built$analysis_id == row$analysis_id[[1]] &
        built$collection == row$collection[[1]] & built$status == "completed"
      if (any(ready) && isTRUE(registry_row$run_kegg_layers[[1]])) {
        species_row <- species_contracts[species_contracts$analysis_id == row$analysis_id[[1]], , drop = FALSE]
        if (nrow(species_row) != 1L) stop("LISA-SPECIES-005 missing inherited species contract for KEGG painter.", call. = FALSE)
        res <- lisa_run_post_script(package_dir, "build_single_de_kegg_pathway_painter.R",
          c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]], "--universe", row$collection[[1]],
            "--species", species_row$scientific_name[[1]], "--kegg-cache-root", kegg_cache_root,
            "--kegg-snapshot-id", kegg_snapshot_id, "--kegg-access-mode", kegg_access_mode),
          trusted_run_root = trusted_run_root, code_ledger = code_ledger)
        add_row("single_de_kegg_pathway_painter", row$analysis_id[[1]], "", row$collection[[1]], res$status, res$output_path, res$message)
      }
    }
    for (i in seq_len(nrow(completed_contrasts))) {
      row <- completed_contrasts[i, , drop = FALSE]
      registry_row <- registry[registry$analysis_collection == row$collection[[1]], , drop = FALSE]
      built <- do.call(rbind, rows)
      ready <- built$stage == "contrast_gene_level_tables" & built$contrast_id == row$contrast_id[[1]] &
        built$collection == row$collection[[1]] & built$status == "completed"
      if (any(ready) && isTRUE(registry_row$run_kegg_layers[[1]])) {
        contrast_def <- contrast_index[contrast_index$contrast_id == row$contrast_id[[1]], , drop = FALSE]
        contrast_species <- lisa_contrast_species(contrast_def, species_contracts)
        if (length(contrast_species) != 1L) {
          add_row("contrast_kegg_pathway_painter", "", row$contrast_id[[1]], row$collection[[1]],
            "skipped_species_mismatch", "", paste0(
              "KEGG contrast map omitted because the two analyses use different organisms: ",
              paste(contrast_species, collapse = " vs "),
              ". Category-level LISA contrast remains valid; a shared organism map would be misleading."))
        } else {
          painter <- lisa_run_post_script(package_dir, "build_contrast_kegg_pathway_painter.R",
            c("--project-dir", output_dir, "--contrast-id", row$contrast_id[[1]], "--universe", row$collection[[1]],
              "--species", contrast_species[[1]], "--kegg-cache-root", kegg_cache_root,
              "--kegg-snapshot-id", kegg_snapshot_id, "--kegg-access-mode", kegg_access_mode),
            trusted_run_root = trusted_run_root, code_ledger = code_ledger)
          add_row("contrast_kegg_pathway_painter", "", row$contrast_id[[1]], row$collection[[1]], painter$status, painter$output_path, painter$message)
        }
      }
    }
  }

  presentation_mode <- "standard"
  if (is.list(.presentation)) {
    presentation_mode <- .presentation$report_mode
    lisa_render_run_presentation(.presentation, output_dir, outputs_root,
      de_index, contrast_index, registry, species_contracts,
      completed_single, completed_contrasts, package_dir, workers, add_row)
  }
  trusted_run_root <- lisa_assert_run_tree_safe(output_dir)
  if (isTRUE(run_reports)) {
    res <- lisa_run_post_script(package_dir, "build_LISA_report.R",
      c("--project-dir", output_dir, "--title", report_title,
        "--report-mode", presentation_mode),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger
    )
    add_row("root_html_report", "", "", paste(registry$analysis_collection, collapse = ","), res$status, res$output_path, res$message)
  }

  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

lisa_post_status_row <- function(stage, analysis_id = "", contrast_id = "", collection = "",
                                 status = "completed", output_path = "", message = "") {
  message <- paste(as.character(message), collapse = "\n")
  message <- gsub("[\r\n\t]+", " | ", message)
  data.frame(
    stage = stage, analysis_id = analysis_id, contrast_id = contrast_id,
    collection = collection, status = status, output_path = output_path,
    message = message, timestamp = as.character(Sys.time()), stringsAsFactors = FALSE
  )
}

lisa_run_single_post_unit <- function(row, registry, output_dir, outputs_root,
                                      file_label_prefix, package_dir, run_gene_level,
                                      trusted_run_root, code_ledger,
                                      category_evidence = FALSE, evidence_max_sets = 25L,
                                      evidence_max_genes = 40L, lisa_dictionary = "core",
                                      gsea_padj_cutoff = 0.25, plot_formats = "png",
                                      source_data = TRUE, recipes = FALSE,
                                      de_index = data.frame()) {
  rows <- list()
  add <- function(stage, result) {
    rows[[length(rows) + 1L]] <<- lisa_post_status_row(
      stage, row$analysis_id[[1]], "", row$collection[[1]],
      result$status, result$output_path, result$message
    )
  }
  registry_row <- registry[registry$analysis_collection == row$collection[[1]], , drop = FALSE]
  if (nrow(registry_row) > 0 && isTRUE(registry_row$run_enrichmentmap[[1]])) {
    add("single_de_enrichmentmap", lisa_run_post_script(package_dir, "build_single_de_enrichmentmap.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]], "--universe", row$collection[[1]],
        "--render-graphs", "true"), trusted_run_root = trusted_run_root,
      code_ledger = code_ledger))
  }
  if (isTRUE(category_evidence) && !identical(row$collection[[1]], "HALLMARKS")) {
    direction <- "not recorded"
    if (all(c("analysis_id", "positive_direction") %in% names(de_index))) {
      d <- as.character(de_index$positive_direction[de_index$analysis_id == row$analysis_id[[1]]])
      if (length(d) == 1L && !is.na(d) && nzchar(d)) direction <- d
    }
    target <- file.path(output_dir, "report_pages", "evidence", row$analysis_id[[1]], row$collection[[1]])
    add("single_de_category_evidence", lisa_run_post_script(package_dir, "build_category_evidence.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]],
        "--universe", row$collection[[1]], "--tier", lisa_dictionary,
        "--gsea-padj-cutoff", format(gsea_padj_cutoff, digits = 17),
        "--positive-contrast", direction, "--max-sets", as.character(evidence_max_sets),
        "--max-genes", as.character(evidence_max_genes), "--formats", paste(plot_formats, collapse = ","),
        "--source-data", tolower(as.character(source_data)),
        "--recipes", tolower(as.character(recipes)), "--output-dir", target),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger))
    add("single_de_category_navigation", lisa_run_post_script(package_dir, "build_category_navigation.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]],
        "--universe", row$collection[[1]], "--tier", lisa_dictionary,
        "--gsea-padj-cutoff", format(gsea_padj_cutoff, digits = 17),
        "--positive-contrast", direction, "--output-dir",
        file.path(output_dir, "report_pages", "category_navigation", row$analysis_id[[1]], row$collection[[1]])),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger))
  }
  if (!isTRUE(run_gene_level)) return(if (length(rows)) do.call(rbind, rows) else data.frame())
  built <- lisa_build_single_gene_level_tables(
    row$analysis_id[[1]], row$collection[[1]], outputs_root, file_label_prefix
  )
  add("single_de_gene_level_tables", built)
  if (identical(built$status, "completed") && isTRUE(registry_row$run_gene_cards[[1]])) {
    add("single_de_gene_cards", lisa_run_post_script(package_dir, "build_single_de_category_gene_cards.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]], "--universe", row$collection[[1]]),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger))
  }
  if (identical(built$status, "completed") && isTRUE(registry_row$run_volcano_overlays[[1]])) {
    add("single_de_volcano_overlays", lisa_run_post_script(package_dir, "build_single_de_category_volcano_overlays.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]], "--universe", row$collection[[1]]),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger))
  }
  if (identical(built$status, "completed") && isTRUE(registry_row$run_recurrent_gene_screen[[1]])) {
    add("single_de_recurrent_gene_screen", lisa_run_post_script(package_dir, "build_single_de_recurrent_gene_screen.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]], "--universe", row$collection[[1]]),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger))
  }
  if (identical(built$status, "completed")) {
    add("single_de_leading_edge_gene_heatmaps", lisa_run_post_script(package_dir, "build_single_de_leading_edge_gene_heatmaps.R",
      c("--project-dir", output_dir, "--analysis-id", row$analysis_id[[1]], "--universe", row$collection[[1]], "--scale", "zscore"),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger))
  }
  do.call(rbind, rows)
}

# Builds presentation only; canonical single/contrast tables remain immutable.
lisa_run_contrast_evidence_unit <- function(row, contrast_index, evidence_status,
    output_dir, package_dir, trusted_run_root, code_ledger,
    lisa_dictionary = "core", max_sets = 25L, max_genes = 40L,
    plot_formats = "png", variants = c("clean", "percentages", "direction", "dispersion"),
    gsea_padj_cutoff = 0.25) {
  id <- as.character(row$contrast_id[[1L]])
  collection <- as.character(row$collection[[1L]])
  definition <- contrast_index[contrast_index$contrast_id == id, , drop = FALSE]
  if (nrow(definition) != 1L) stop("LISA-CONTRAST-EVIDENCE-001 contrast identity is ambiguous.", call. = FALSE)
  a <- as.character(definition$analysis_a[[1L]]); b <- as.character(definition$analysis_b[[1L]])
  scope <- lisa_contrast_name(definition)
  target <- file.path(output_dir, "report_pages", "contrast_evidence", scope, collection)
  ready <- if (all(c("stage", "analysis_id", "collection", "status") %in% names(evidence_status))) {
    vapply(c(a, b), function(id) {
      selected <- evidence_status[evidence_status$stage == "single_de_category_evidence" &
        evidence_status$analysis_id == id & evidence_status$collection == collection, , drop = FALSE]
      nrow(selected) == 1L && identical(selected$status[[1L]], "completed")
    }, logical(1L))
  } else c(FALSE, FALSE)
  if (!all(ready)) return(lisa_post_status_row("contrast_category_evidence", "", id,
    collection, "skipped_individual_evidence_not_completed", target,
    paste("Exact individual evidence is not completed for", paste(c(a, b)[!ready], collapse = ", "))))
  res <- lisa_run_post_script(package_dir, "build_contrast_evidence.R",
    c("--project-dir", output_dir, "--contrast-id", id,
      "--analysis-a", a, "--analysis-b", b,
      "--universe", collection, "--tier", lisa_dictionary,
      "--evidence-a", file.path(output_dir, "report_pages", "evidence", a, collection),
      "--evidence-b", file.path(output_dir, "report_pages", "evidence", b, collection),
      "--output-dir", target, "--max-sets", as.character(max_sets),
      "--max-genes", as.character(max_genes), "--formats", paste(plot_formats, collapse = ","),
      "--variants", paste(variants, collapse = ","),
      "--report-href", "../../../contrasts.html",
      "--gene-explorer-href", "../../../gene_evidence/index.html"),
    trusted_run_root = trusted_run_root, code_ledger = code_ledger)
  rows <- list(lisa_post_status_row("contrast_category_evidence", "", id, collection,
    res$status, res$output_path, res$message))
  if (identical(res$status, "completed")) {
    label <- if ("contrast_label" %in% names(definition) && !is.na(definition$contrast_label[[1L]]) &&
      nzchar(as.character(definition$contrast_label[[1L]]))) as.character(definition$contrast_label[[1L]]) else id
    nav_target <- file.path(output_dir, "report_pages", "category_navigation", scope, collection)
    nav <- lisa_run_post_script(package_dir, "build_contrast_navigation.R",
      c("--project-dir", output_dir, "--contrast-id", id, "--scope", scope,
        "--analysis-a", a, "--analysis-b", b, "--universe", collection, "--tier", lisa_dictionary,
        "--gsea-padj-cutoff", format(gsea_padj_cutoff, digits = 17),
        "--label", label, "--output-dir", nav_target),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger)
    rows[[length(rows) + 1L]] <- lisa_post_status_row("contrast_category_navigation", "", id, collection,
      nav$status, nav$output_path, nav$message)
  }
  do.call(rbind, rows)
}

lisa_run_contrast_post_unit <- function(row, output_dir, package_dir, trusted_run_root,
                                        code_ledger) {
  rows <- list()
  res <- lisa_run_post_script(package_dir, "build_contrast_gene_level_product.R",
    c("--project-dir", output_dir, "--contrast-id", row$contrast_id[[1]], "--universe", row$collection[[1]]),
    trusted_run_root = trusted_run_root, code_ledger = code_ledger)
  rows[[1L]] <- lisa_post_status_row("contrast_gene_level_tables", "", row$contrast_id[[1]],
    row$collection[[1]], res$status, res$output_path, res$message)
  if (identical(res$status, "completed")) {
    net <- lisa_run_post_script(package_dir, "build_contrast_gene_category_network.R",
      c("--project-dir", output_dir, "--contrast-id", row$contrast_id[[1]], "--universe", row$collection[[1]]),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger)
    rows[[2L]] <- lisa_post_status_row("contrast_gene_category_networks", "", row$contrast_id[[1]],
      row$collection[[1]], net$status, net$output_path, net$message)
  }
  do.call(rbind, rows)
}

lisa_contrast_species <- function(contrast_def, species_contracts) {
  if (nrow(contrast_def) != 1L ||
      !all(c("analysis_a", "analysis_b") %in% names(contrast_def))) {
    stop("LISA-SPECIES-006 contrast definition must contain one analysis_a and one analysis_b.", call. = FALSE)
  }
  analysis_ids <- as.character(c(
    contrast_def$analysis_a[[1]],
    contrast_def$analysis_b[[1]]
  ))
  matched <- match(analysis_ids, as.character(species_contracts$analysis_id))
  if (anyNA(matched)) {
    missing <- analysis_ids[is.na(matched)]
    stop(
      "LISA-SPECIES-007 missing inherited species contract for contrast analysis: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  species <- as.character(species_contracts$scientific_name[matched])
  if (anyNA(species) || any(!nzchar(species))) {
    stop("LISA-SPECIES-008 inherited contrast species must not be empty.", call. = FALSE)
  }
  unique(species)
}

lisa_build_single_gene_level_tables <- function(analysis_id, collection, outputs_root, file_label_prefix) {
  collection_dir <- file.path(outputs_root, "single_de", analysis_id, paste0("collection_", collection))
  gene_dir <- file.path(outputs_root, "gene_level", "single_de", analysis_id, paste0("collection_", collection))
  de_path <- file.path(collection_dir, "inputs", paste0(analysis_id, "_standardized_DE.tsv"))
  gsea_path <- file.path(collection_dir, "enrichment", paste0(analysis_id, "_GSEA_", file_label_prefix, "_annotated.tsv"))
  ora_path <- file.path(collection_dir, "enrichment", paste0(analysis_id, "_ORA_", file_label_prefix, "_annotated.tsv"))
  if (!file.exists(de_path)) {
    return(list(status = "skipped_missing_input", output_path = gene_dir, message = paste("missing standardized DE:", de_path)))
  }
  if (!file.exists(gsea_path)) {
    return(list(status = "skipped_missing_input", output_path = gene_dir, message = paste("missing GSEA table:", gsea_path)))
  }
  msg <- ""
  status <- "completed"
  tryCatch({
    build_lisa_gene_level_tables(
      de_table = de_path,
      gsea_table = gsea_path,
      ora_table = if (file.exists(ora_path)) ora_path else NULL,
      output_dir = gene_dir,
      symbol_col = "symbol",
      logfc_col = "log2FoldChange",
      padj_col = "padj",
      analysis_id = analysis_id,
      collection = collection
    )
  }, error = function(e) {
    status <<- "failed"
    msg <<- conditionMessage(e)
  })
  list(status = status, output_path = gene_dir, message = msg)
}

lisa_code_identity_abort <- function(code, message, files = character()) {
  condition <- structure(
    list(
      message = paste(code, message), call = NULL, code = code,
      files = as.character(files)
    ),
    class = c("lisa_code_identity_error", "lisa_error", "error", "condition")
  )
  stop(condition)
}

lisa_resolve_package_dir <- function(.test_package_dir = NULL) {
  # Production has exactly one authority for executable code: the installed
  # lisaR namespace.  Source checkouts are available only through an explicit
  # internal test/development injection and can never arrive from YAML/JSON.
  test_injection <- !is.null(.test_package_dir)
  injection_kind <- if (isTRUE(test_injection)) "explicit" else "none"

  if (isTRUE(test_injection)) {
    if (!is.character(.test_package_dir) || length(.test_package_dir) != 1L ||
        is.na(.test_package_dir) || !nzchar(.test_package_dir)) {
      lisa_code_identity_abort(
        "LISA-CODE-ROOT-002",
        "the internal .test_package_dir injection must name one package tree."
      )
    }
    package_dir <- normalizePath(.test_package_dir, winslash = "/", mustWork = TRUE)
    candidates <- c(
      file.path(package_dir, "inst", "scripts"),
      file.path(package_dir, "scripts")
    )
    candidates <- candidates[dir.exists(candidates)]
    description <- tryCatch(
      read.dcf(file.path(package_dir, "DESCRIPTION"), fields = "Package"),
      error = function(error) NULL
    )
    if (is.null(description) || nrow(description) != 1L ||
        !identical(as.character(description[1L, "Package"]), "lisaR") ||
        !length(candidates)) {
      lisa_code_identity_abort(
        "LISA-CODE-ROOT-002",
        "the explicit development injection is not a lisaR source or installed package tree.",
        .test_package_dir
      )
    }
    script_dir <- normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
  } else {
    raw_package_dir <- system.file(package = "lisaR")
    raw_script_dir <- system.file("scripts", package = "lisaR")
    if (!nzchar(raw_package_dir) || !nzchar(raw_script_dir) || !dir.exists(raw_script_dir)) {
      lisa_code_identity_abort(
        "LISA-CODE-ROOT-001",
        "cannot locate the scripts in the installed lisaR package; install lisaR before running."
      )
    }
    package_dir <- normalizePath(raw_package_dir, winslash = "/", mustWork = TRUE)
    script_dir <- normalizePath(raw_script_dir, winslash = "/", mustWork = TRUE)
    # pkgload maps system.file(package = "lisaR") to the source `inst/`
    # directory. Treat that layout as an explicit development injection only
    # when its parent is demonstrably the lisaR source tree. Installed packages
    # keep DESCRIPTION directly under package_dir and never take this branch.
    if (!file.exists(file.path(package_dir, "DESCRIPTION"))) {
      source_root <- normalizePath(dirname(package_dir), winslash = "/",
                                   mustWork = TRUE)
      description <- tryCatch(
        read.dcf(file.path(source_root, "DESCRIPTION"), fields = "Package"),
        error = function(error) NULL
      )
      expected_scripts <- normalizePath(
        file.path(source_root, "inst", "scripts"), winslash = "/",
        mustWork = FALSE
      )
      if (is.null(description) || nrow(description) != 1L ||
          !identical(as.character(description[1L, "Package"]), "lisaR") ||
          !identical(script_dir, expected_scripts)) {
        lisa_code_identity_abort(
          "LISA-CODE-ROOT-001",
          "cannot bind the loaded lisaR namespace to one package tree."
        )
      }
      package_dir <- source_root
      test_injection <- TRUE
      injection_kind <- "pkgload"
    }
  }

  if (!identical(script_dir, package_dir) &&
      !startsWith(script_dir, paste0(package_dir, "/"))) {
    lisa_code_identity_abort(
      "LISA-CODE-ROOT-003",
      "the resolved script directory is outside the selected lisaR package tree.",
      script_dir
    )
  }
  attr(package_dir, "lisaR.script_dir") <- script_dir
  attr(package_dir, "lisaR.test_injection") <- isTRUE(test_injection)
  attr(package_dir, "lisaR.test_injection_kind") <- injection_kind
  package_dir
}

lisa_package_script_dir <- function(package_dir) {
  script_dir <- attr(package_dir, "lisaR.script_dir", exact = TRUE)
  if (is.null(script_dir)) {
    installed <- lisa_resolve_package_dir(.test_package_dir = NULL)
    supplied <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
    if (!identical(supplied, as.character(installed)[[1L]])) {
      lisa_code_identity_abort(
        "LISA-CODE-ROOT-004",
        "an unverified package directory reached the executable-code boundary.",
        supplied
      )
    }
    package_dir <- installed
    script_dir <- attr(installed, "lisaR.script_dir", exact = TRUE)
  }

  package_root <- gsub("\\\\", "/", as.character(package_dir)[[1L]])
  script_dir <- gsub("\\\\", "/", as.character(script_dir)[[1L]])
  if (!dir.exists(package_root) || lisa_path_is_link(package_root)) {
    lisa_code_identity_abort(
      "LISA-CODE-ROOT-005",
      "the selected lisaR package root disappeared or became a symbolic link.",
      package_root
    )
  }
  resolved_root <- normalizePath(
    package_root, winslash = "/", mustWork = TRUE
  )
  if (!identical(lisa_path_key(resolved_root), lisa_path_key(package_root)) ||
      !lisa_path_within(script_dir, package_root)) {
    lisa_code_identity_abort(
      "LISA-CODE-ROOT-005",
      "the selected lisaR code root changed after it was resolved.",
      c(package_root, script_dir)
    )
  }

  assert_descendant_components <- function() {
    relative <- substring(script_dir, nchar(package_root) + 2L)
    parts <- strsplit(relative, "/", fixed = TRUE)[[1L]]
    current <- package_root
    for (part in parts[nzchar(parts)]) {
      current <- paste0(sub("/+$", "", current), "/", part)
      if (lisa_path_is_link(current)) {
        lisa_code_identity_abort(
          "LISA-CODE-ROOT-005",
          "the selected lisaR scripts subtree contains a symbolic link.",
          current
        )
      }
    }
    invisible(TRUE)
  }
  assert_descendant_components()
  if (!dir.exists(script_dir)) {
    lisa_code_identity_abort(
      "LISA-CODE-ROOT-005",
      "the selected lisaR scripts directory is no longer available.",
      script_dir
    )
  }
  resolved_script_dir <- normalizePath(
    script_dir, winslash = "/", mustWork = TRUE
  )
  assert_descendant_components()
  if (!identical(lisa_path_key(resolved_script_dir), lisa_path_key(script_dir)) ||
      !lisa_path_within(resolved_script_dir, resolved_root)) {
    lisa_code_identity_abort(
      "LISA-CODE-ROOT-005",
      "the selected lisaR scripts directory escaped its package root.",
      resolved_script_dir
    )
  }
  resolved_script_dir
}

lisa_post_script_dependencies <- function() {
  # Keys are exactly the Rscript entry points launched by lisaR. Values are
  # sibling R files sourced or copied as executable recipes by that entry
  # point. A renderer copied into an output bundle is executable code too, so
  # it must be covered by the same pre-launch ledger and re-hash as its builder.
  list(
    "build_LISA_report.R" = c("reproduce_lisa_figure.R", "reproduce_lisa_category_nes.R"),
    "build_category_evidence.R" = character(),
    "build_category_navigation.R" = character(),
    "build_gene_evidence.R" = character(),
    "build_contrast_evidence.R" = character(),
    "build_contrast_navigation.R" = character(),
    "build_contrast_gene_category_network.R" = character(),
    "build_contrast_gene_level_product.R" = "reproduce_lisa_figure.R",
    "build_contrast_kegg_pathway_painter.R" = c(
      "kegg_snapshot_helpers.R", "build_contrast_kegg_map_layer.R",
      "reproduce_lisa_figure.R"
    ),
    "build_contrast_macrogroup_heatmaps.R" = character(),
    "build_single_de_category_gene_cards.R" = "reproduce_lisa_figure.R",
    "build_single_de_category_volcano_overlays.R" = c(
      "lisa_plot_metadata.R", "reproduce_lisa_figure.R"
    ),
    "build_single_de_enrichmentmap.R" = "lisa_plot_metadata.R",
    "build_single_de_kegg_pathway_painter.R" = c(
      "kegg_snapshot_helpers.R", "lisa_plot_metadata.R",
      "reproduce_lisa_figure.R"
    ),
    "build_single_de_leading_edge_gene_heatmaps.R" = c(
      "lisa_plot_metadata.R", "lisa_leading_edge_heatmap_native.R"
    ),
    "build_single_de_recurrent_gene_screen.R" = "lisa_plot_metadata.R"
  )
}

lisa_post_script_path <- function(package_dir, script_name) {
  if (!is.character(script_name) || length(script_name) != 1L ||
      is.na(script_name) || !nzchar(script_name) ||
      !identical(basename(script_name), script_name) ||
      grepl("[[:cntrl:]/\\\\]", script_name)) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-003",
      "a post-processing script name must be one plain file name."
    )
  }
  script_dir <- lisa_package_script_dir(package_dir)
  script_path <- file.path(script_dir, script_name)
  if (lisa_path_is_link(script_path)) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-001",
      paste0("required installed script/helper must not be a symlink: ", script_name, "."),
      script_path
    )
  }
  regular <- file.exists(script_path) && !dir.exists(script_path) &&
    isTRUE(utils::file_test("-f", script_path))
  if (!isTRUE(regular)) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-001",
      paste0("required installed script/helper is missing: ", script_name, "."),
      script_path
    )
  }
  resolved <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  if (!startsWith(resolved, paste0(script_dir, "/"))) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-001",
      paste0("required installed script/helper resolves outside its package tree: ", script_name, "."),
      resolved
    )
  }
  resolved
}

lisa_post_script_identity <- function(package_dir, script_name) {
  dependencies <- lisa_post_script_dependencies()
  if (!script_name %in% names(dependencies)) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-003",
      paste0("no executable-code closure is registered for ", script_name, ".")
    )
  }
  relative <- unique(c(script_name, dependencies[[script_name]]))
  paths <- vapply(
    relative,
    function(name) lisa_post_script_path(package_dir, name),
    character(1), USE.NAMES = FALSE
  )
  relative <- file.path("scripts", relative)
  if (script_name %in% c("build_category_evidence.R", "build_gene_evidence.R", "build_contrast_evidence.R")) {
    asset_name <- switch(script_name, "build_category_evidence.R" = "category-evidence",
      "build_gene_evidence.R" = "gene-evidence", "build_contrast_evidence.R" = "contrast-evidence")
    asset_root <- file.path(dirname(lisa_package_script_dir(package_dir)), asset_name)
    assets <- file.path(asset_root, c("viewer.html", "viewer.css", "viewer.js"))
    for (asset in assets) {
      lisa_assert_regular_managed_file(asset, dirname(asset_root), "Category evidence asset")
    }
    paths <- c(paths, normalizePath(assets, winslash = "/", mustWork = TRUE))
    relative <- c(relative, file.path(asset_name, basename(assets)))
  }
  if (script_name %in% c("build_LISA_report.R", "build_category_evidence.R",
      "build_category_navigation.R", "build_gene_evidence.R", "build_contrast_evidence.R",
      "build_contrast_navigation.R")) {
    root <- file.path(dirname(lisa_package_script_dir(package_dir)), "report_assets")
    assets <- file.path(root, c("lisa_shell.css", "lisa_shell.js",
      "LISA_logo_C_compact_icon_muted_red_S.svg"))
    for (asset in assets) lisa_assert_regular_managed_file(asset, dirname(root), "Shared report asset")
    paths <- c(paths, normalizePath(assets, winslash = "/", mustWork = TRUE))
    relative <- c(relative, file.path("report_assets", basename(assets)))
  }
  hashes <- vapply(paths, lisa_sha256_file, character(1), USE.NAMES = FALSE)
  structure(
    list(
      package_dir = package_dir,
      script_name = script_name,
      files = data.frame(
        relative_path = relative,
        absolute_path = paths,
        sha256 = hashes,
        stringsAsFactors = FALSE
      )
    ),
    class = "lisa_post_script_identity"
  )
}

lisa_verify_post_script_identity <- function(identity) {
  if (!inherits(identity, "lisa_post_script_identity")) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-003",
      "the executable-code identity receipt is malformed."
    )
  }
  current <- lisa_post_script_identity(identity$package_dir, identity$script_name)
  if (!identical(current$files, identity$files)) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-002",
      paste0(
        "the installed script/helper closure changed after inventory for ",
        identity$script_name, "; reinstall lisaR and start a new run."
      ),
      identity$files$absolute_path
    )
  }
  invisible(current)
}

lisa_code_identity_rows <- function(package_dir, script_names) {
  script_names <- sort(unique(as.character(script_names)))
  script_names <- script_names[nzchar(script_names)]
  description <- read.dcf(
    file.path(as.character(package_dir), "DESCRIPTION"),
    fields = c("Package", "Version")
  )
  package_name <- as.character(description[1L, "Package"])
  package_version <- as.character(description[1L, "Version"])
  package_root <- if (isTRUE(attr(package_dir, "lisaR.test_injection", exact = TRUE))) {
    "internal-test:lisaR-source"
  } else {
    "installed:lisaR"
  }
  empty <- data.frame(
    package = character(), package_version = character(),
    package_root = character(), entrypoint = character(),
    relative_path = character(), sha256 = character(),
    stringsAsFactors = FALSE
  )
  if (!length(script_names)) return(empty)
  rows <- lapply(script_names, function(script_name) {
    identity <- lisa_post_script_identity(package_dir, script_name)
    data.frame(
      package = package_name,
      package_version = package_version,
      package_root = package_root,
      entrypoint = script_name,
      relative_path = identity$files$relative_path,
      sha256 = identity$files$sha256,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lisa_write_code_identity_ledger <- function(package_dir, script_names, path,
                                            executable_recipes = character()) {
  entries <- lisa_code_identity_rows(package_dir, script_names)
  executable_recipes <- sort(unique(as.character(executable_recipes)))
  executable_recipes <- executable_recipes[nzchar(executable_recipes)]
  if (length(executable_recipes)) {
    description <- read.dcf(
      file.path(as.character(package_dir), "DESCRIPTION"),
      fields = c("Package", "Version")
    )
    package_root <- if (isTRUE(attr(package_dir, "lisaR.test_injection", exact = TRUE))) {
      "internal-test:lisaR-source"
    } else {
      "installed:lisaR"
    }
    recipe_rows <- lapply(executable_recipes, function(recipe) {
      recipe_path <- lisa_post_script_path(package_dir, recipe)
      data.frame(
        package = as.character(description[1L, "Package"]),
        package_version = as.character(description[1L, "Version"]),
        package_root = package_root,
        entrypoint = paste0("copy:", recipe),
        relative_path = file.path("scripts", recipe),
        sha256 = lisa_sha256_file(recipe_path),
        stringsAsFactors = FALSE
      )
    })
    entries <- rbind(entries, do.call(rbind, recipe_rows))
  }
  write_lisa_tsv(entries, path)
  structure(
    list(
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      sha256 = lisa_sha256_file(path),
      entries = entries
    ),
    class = "lisa_code_identity_ledger"
  )
}

lisa_verify_code_recipe_ledger <- function(ledger, package_dir, recipe_name) {
  if (!inherits(ledger, "lisa_code_identity_ledger") || !file.exists(ledger$path)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-001",
      "the executable-code ledger is absent or malformed."
    )
  }
  if (!identical(lisa_sha256_file(ledger$path), ledger$sha256)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-002",
      "code_identity.tsv changed after its serial pre-execution snapshot.",
      ledger$path
    )
  }
  recipe_path <- lisa_post_script_path(package_dir, recipe_name)
  expected <- ledger$entries[
    ledger$entries$entrypoint == paste0("copy:", recipe_name),
    c("relative_path", "sha256"), drop = FALSE
  ]
  observed <- data.frame(
    relative_path = file.path("scripts", recipe_name),
    sha256 = lisa_sha256_file(recipe_path),
    stringsAsFactors = FALSE
  )
  rownames(expected) <- NULL
  if (nrow(expected) != 1L || !identical(expected, observed)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-003",
      paste0("the executable recipe does not match code_identity.tsv: ", recipe_name, "."),
      recipe_path
    )
  }
  invisible(list(path = recipe_path, sha256 = observed$sha256[[1L]]))
}

lisa_copy_verified_code_file <- function(source, target, expected_sha256,
                                         run_root = getOption("lisaR.run_root", NULL),
                                         .copy = NULL) {
  if ((!is.null(.copy) && !is.function(.copy)) ||
      !lisa_sha256_all_valid(expected_sha256, 1L)) {
    lisa_code_identity_abort(
      "LISA-CODE-COPY-001",
      "an executable-code copy requires one valid expected SHA-256 digest.",
      source
    )
  }
  source_before <- lisa_sha256_file(source)
  if (!identical(source_before, as.character(expected_sha256)[[1L]])) {
    lisa_code_identity_abort(
      "LISA-CODE-COPY-001",
      "the executable source changed immediately before copy.",
      source
    )
  }
  copy_root <- lisa_run_root(if (is.null(run_root)) dirname(target) else run_root)
  target <- lisa_guarded_path(target, copy_root)
  copied <- if (is.null(.copy)) {
    lisa_guarded_copy(
      source, target, overwrite = TRUE, run_root = copy_root
    )
    TRUE
  } else {
    lisa_guarded_dir_create(dirname(target), copy_root)
    isTRUE(.copy(source, target, overwrite = TRUE))
  }
  valid_copy <- FALSE
  if (copied && file.exists(target)) {
    valid_copy <- tryCatch(
      identical(lisa_sha256_file(source), source_before) &&
        identical(lisa_sha256_file(target), source_before),
      error = function(error) FALSE
    )
  }
  if (!isTRUE(valid_copy)) {
    if (file.exists(target)) {
      lisa_guarded_delete(target, run_root = copy_root)
    }
    lisa_code_identity_abort(
      "LISA-CODE-COPY-001",
      "the executable source or copied target failed post-copy identity verification; the target was removed.",
      c(source, target)
    )
  }
  invisible(target)
}

# Install the portable renderer from the code closure already verified by the
# parent process. Direct, standalone script use has no parent receipt, so it
# takes one immediate source snapshot and then delegates to the same guarded
# copy primitive. A supplied but malformed digest is never treated as absent.
lisa_copy_verified_figure_recipe <- function(
    source, target, expected_sha256 = NA_character_,
    run_root = getOption("lisaR.run_root", NULL)) {
  digest_absent <- is.null(expected_sha256) || length(expected_sha256) == 0L ||
    (length(expected_sha256) == 1L && is.na(expected_sha256))
  if (isTRUE(digest_absent)) {
    expected_sha256 <- lisa_sha256_file(source)
  } else if (!lisa_sha256_all_valid(expected_sha256, 1L)) {
    lisa_code_identity_abort(
      "LISA-CODE-COPY-001",
      "the private renderer handoff requires one valid SHA-256 digest.",
      source
    )
  }
  lisa_copy_verified_code_file(
    source, target, as.character(expected_sha256)[[1L]], run_root = run_root
  )
}

lisa_verify_code_identity_ledger <- function(ledger, identity) {
  if (!inherits(ledger, "lisa_code_identity_ledger") ||
      !inherits(identity, "lisa_post_script_identity") ||
      !file.exists(ledger$path)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-001",
      "the executable-code ledger is absent or malformed."
    )
  }
  current_ledger_hash <- lisa_sha256_file(ledger$path)
  if (!identical(current_ledger_hash, ledger$sha256)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-002",
      "code_identity.tsv changed after its serial pre-execution snapshot.",
      ledger$path
    )
  }
  expected <- ledger$entries[
    ledger$entries$entrypoint == identity$script_name,
    c("relative_path", "sha256"), drop = FALSE
  ]
  observed <- identity$files[, c("relative_path", "sha256"), drop = FALSE]
  rownames(expected) <- NULL
  rownames(observed) <- NULL
  if (!nrow(expected) || !identical(expected, observed)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-003",
      paste0(
        "the current script/helper closure does not match code_identity.tsv for ",
        identity$script_name, "."
      ),
      identity$files$absolute_path
    )
  }
  invisible(TRUE)
}

lisa_post_script_renderer_sha256 <- function(identity, ledger = NULL) {
  renderer_path <- file.path("scripts", "reproduce_lisa_figure.R")
  identity_rows <- identity$files[
    identity$files$relative_path == renderer_path,
    c("relative_path", "sha256"), drop = FALSE
  ]
  if (!nrow(identity_rows)) return(NULL)
  if (nrow(identity_rows) != 1L ||
      !lisa_sha256_all_valid(identity_rows$sha256, 1L)) {
    lisa_code_identity_abort(
      "LISA-CODE-IDENTITY-003",
      "the verified post-processing closure has an ambiguous renderer identity."
    )
  }
  if (is.null(ledger)) return(as.character(identity_rows$sha256)[[1L]])

  ledger_rows <- ledger$entries[
    ledger$entries$entrypoint == identity$script_name &
      ledger$entries$relative_path == renderer_path,
    c("relative_path", "sha256"), drop = FALSE
  ]
  rownames(identity_rows) <- NULL
  rownames(ledger_rows) <- NULL
  if (nrow(ledger_rows) != 1L ||
      !lisa_sha256_all_valid(ledger_rows$sha256, 1L) ||
      !identical(ledger_rows, identity_rows)) {
    lisa_code_identity_abort(
      "LISA-CODE-LEDGER-003",
      paste0(
        "the renderer identity is absent or ambiguous in code_identity.tsv for ",
        identity$script_name, "."
      )
    )
  }
  as.character(ledger_rows$sha256)[[1L]]
}

lisa_assert_child_runtime_identity <- function(
    package_dir, .namespace = NULL,
    .source_environment = environment(lisa_run_post_script)) {
  if (!isTRUE(attr(package_dir, "lisaR.test_injection", exact = TRUE))) {
    return(invisible(TRUE))
  }
  if (identical(
      attr(package_dir, "lisaR.test_injection_kind", exact = TRUE), "pkgload"
  )) {
    lisa_code_identity_abort(
      "LISA-CODE-RUNTIME-001",
      paste0(
        "post-processing subprocesses cannot be launched from a pkgload ",
        "source namespace because external Rscript would load a different ",
        "installed lisaR; install this reviewed tree and run its installed package."
      )
    )
  }
  if (is.null(.namespace)) {
    .namespace <- tryCatch(asNamespace("lisaR"), error = function(error) NULL)
  }
  source_names <- getOption("lisaR.source_runtime_functions", NULL)
  if (!is.environment(.namespace) || !is.character(source_names) ||
      !length(source_names) || anyNA(source_names) || any(!nzchar(source_names))) {
    lisa_code_identity_abort(
      "LISA-CODE-RUNTIME-001",
      paste0(
        "the source-only runner cannot prove that child Rscript processes ",
        "will load the same lisaR implementation; install this reviewed tree first."
      )
    )
  }
  source_names <- sort(unique(source_names))
  mismatched <- source_names[vapply(source_names, function(name) {
    if (!exists(name, envir = .source_environment, inherits = FALSE) ||
        !exists(name, envir = .namespace, inherits = FALSE)) return(TRUE)
    source_function <- get(name, envir = .source_environment, inherits = FALSE)
    child_function <- get(name, envir = .namespace, inherits = FALSE)
    !is.function(source_function) || !is.function(child_function) ||
      !identical(formals(source_function), formals(child_function)) ||
      !identical(
        body(source_function), body(child_function),
        ignore.bytecode = TRUE, ignore.environment = TRUE,
        ignore.srcref = TRUE
      )
  }, logical(1))]
  if (length(mismatched)) {
    lisa_code_identity_abort(
      "LISA-CODE-RUNTIME-001",
      paste0(
        "the installed lisaR namespace differs from the reviewed source-only ",
        "runner implementation for: ", paste(mismatched, collapse = ", "),
        ". Install this reviewed tree before running."
      )
    )
  }
  invisible(TRUE)
}

lisa_run_post_script <- function(package_dir, script_name, args, trusted_run_root = NULL,
                                 code_ledger = NULL) {
  lisa_validate_post_script_args(script_name, args)
  code_identity <- lisa_post_script_identity(package_dir, script_name)
  script_path <- code_identity$files$absolute_path[[1L]]
  cmd <- lisa_rscript_executable()
  project_flag <- match("--project-dir", args)
  if (is.na(project_flag) || project_flag == length(args)) {
    stop("Post-processing scripts require an explicit --project-dir run_root.", call. = FALSE)
  }
  project_root <- lisa_existing_run_root(args[[project_flag + 1L]])
  if (is.null(trusted_run_root)) {
    run_root <- lisa_assert_run_tree_safe(project_root)
  } else {
    run_root <- lisa_existing_run_root(trusted_run_root)
    if (!identical(project_root, run_root)) {
      stop("Trusted post-processing run_root does not match --project-dir.", call. = FALSE)
    }
  }
  # Re-hash immediately before the subprocess boundary. A missing or changed
  # script/helper must fail before any child-process marker can be produced.
  verified_identity <- lisa_verify_post_script_identity(code_identity)
  if (!is.null(code_ledger)) {
    lisa_verify_code_identity_ledger(code_ledger, verified_identity)
  }
  renderer_sha256 <- lisa_post_script_renderer_sha256(
    verified_identity, code_ledger
  )
  lisa_assert_child_runtime_identity(package_dir)
  all_args <- c("--vanilla", script_path, args)
  if (!is.null(renderer_sha256)) {
    # This flag is intentionally appended only after public argument
    # validation. It is not part of lisa_post_script_interfaces(), so callers
    # cannot inject or override the parent's frozen renderer identity.
    all_args <- c(
      all_args, "--lisa-internal-renderer-sha256", renderer_sha256
    )
  }
  package_library <- normalizePath(
    dirname(package_dir), winslash = "/", mustWork = TRUE
  )
  child_libraries <- unique(c(package_library, .libPaths()))
  result <- lisa_run_subprocess(
    cmd, all_args, required = FALSE, stage = script_name, run_root = run_root,
    env = c(R_LIBS = paste(child_libraries, collapse = .Platform$path.sep))
  )
  msg <- trimws(as.character(result$diagnostics %||% ""))
  status <- if (identical(result$exit_code, 0L)) "completed" else "failed"
  if (!identical(status, "completed") && grepl("No KEGG|No canonical KEGG|No KEGG-derived", msg, ignore.case = TRUE)) {
    status <- "skipped_no_kegg_rows"
  }
  if (!identical(result$exit_code, 0L)) {
    msg <- paste0(
      "exit_code=", as.integer(result$exit_code), "; ",
      if (nzchar(msg)) msg else "no stdout/stderr"
    )
  }
  list(
    status = status, output_path = "", message = msg,
    exit_code = as.integer(result$exit_code),
    code_identity = code_identity$files[, c("relative_path", "sha256"), drop = FALSE]
  )
}


# Presentation products of run_lisa() rendered inside the run's own staging tree
# before its one report generator: FULL products (mode = "full") and native KEGG
# maps (run_kegg_maps). Nothing here recomputes DE or enrichment. A failure of a
# FULL product fails the run (it is part of the same transaction); an absent or
# incomplete KEGG snapshot is reported visibly instead of failing the run.
lisa_render_run_presentation <- function(presentation, output_dir, outputs_root,
                                         de_index, contrast_index, registry,
                                         species_contracts, completed_single,
                                         completed_contrasts, package_dir,
                                         workers, add_row) {
  formats <- unlist(presentation$report$formats)
  formats <- vapply(c("png", "svg", "pdf"), function(f) isTRUE(as.logical(formats[[f]])), logical(1))
  names(formats) <- c("png", "svg", "pdf")
  requested_formats <- names(formats)[formats]
  if (identical(presentation$report_mode, "full")) {
    warning("LISA-REPORT-SCALE-001 full mode renders every applicable category figure inside this run and may create thousands of figures and use gigabytes.", call. = FALSE)
    started <- Sys.time()
    lisa_render_run_full_products(output_dir, presentation$report, workers = workers)
    add_row("full_report_products", status = "completed",
      output_path = file.path(output_dir, "artifacts"),
      message = sprintf("FULL category products rendered in %.1f min.",
        as.numeric(difftime(Sys.time(), started, units = "mins"))))
    started <- Sys.time()
    lisa_complete_full_products(output_dir, kegg_maps = FALSE,
      formats = requested_formats,
      families = c("recurrent_genes", "paired_gene_products", "gene_category_network"))
    add_row("full_gene_level_products", status = "completed",
      output_path = file.path(output_dir, "outputs", "gene_level"),
      message = sprintf("Recurrent genes and paired contrast products rendered in %.1f min.",
        as.numeric(difftime(Sys.time(), started, units = "mins"))))
  }
  if (isTRUE(presentation$kegg_maps)) {
    lisa_render_run_kegg_maps(presentation, output_dir, outputs_root, de_index,
      contrast_index, registry, species_contracts, completed_single,
      completed_contrasts, package_dir, add_row)
  }
  lisa_apply_format_policy(file.path(output_dir, "outputs", "gene_level"), formats)
  if (identical(presentation$report_mode, "full")) {
    lisa_write_full_products_inventory(output_dir)
  }
  invisible(TRUE)
}

lisa_kegg_snapshot_problem <- function(cache_root, snapshot_id, organisms) {
  if (!is.character(cache_root) || length(cache_root) != 1L || !nzchar(cache_root) ||
      !dir.exists(cache_root)) {
    return(paste0("the KEGG cache root does not exist: ", cache_root))
  }
  for (organism in unique(organisms)) {
    problem <- tryCatch({
      lisa_read_kegg_snapshot(cache_root, snapshot_id, organism, "pathway_list", "all")
      lisa_read_kegg_snapshot(cache_root, snapshot_id, organism, "pathway_links", "all")
      ""
    }, error = function(error) conditionMessage(error))
    if (nzchar(problem)) return(problem)
  }
  ""
}

lisa_render_run_kegg_maps <- function(presentation, output_dir, outputs_root,
                                      de_index, contrast_index, registry,
                                      species_contracts, completed_single,
                                      completed_contrasts, package_dir, add_row) {
  cache_root <- as.character(presentation$kegg_cache_root)
  snapshot_id <- as.character(presentation$kegg_snapshot_id)
  status_path <- file.path(output_dir, "kegg_maps_status.tsv")
  kegg_collections <- registry$analysis_collection[registry$run_kegg_layers]
  single <- completed_single[completed_single$collection %in% kegg_collections, , drop = FALSE]
  contrasts <- completed_contrasts[completed_contrasts$collection %in% kegg_collections, , drop = FALSE]
  organisms <- vapply(as.character(single$analysis_id), function(id) {
    row <- species_contracts[species_contracts$analysis_id == id, , drop = FALSE]
    lisa_species_contract(row$scientific_name[[1L]])$kegg_code
  }, character(1))
  problem <- lisa_kegg_snapshot_problem(cache_root, snapshot_id, organisms)
  status_rows <- list()
  note <- function(kind, owner, collection, status, message) {
    status_rows[[length(status_rows) + 1L]] <<- data.frame(kind = kind,
      owner = owner, collection = collection, status = status,
      message = gsub("[\r\n\t]+", " | ", message), stringsAsFactors = FALSE)
    # Persist before any error propagates; partial runs retain the original
    # painter log and the reason they failed.
    write_lisa_tsv(do.call(rbind, status_rows), status_path)
  }
  if (nzchar(problem)) {
    message <- paste0("KEGG maps were requested but no usable KEGG snapshot is available (",
      "snapshot '", snapshot_id, "' under ", cache_root, "): ", problem,
      ". The report is complete without native KEGG maps.")
    warning("LISA-KEGG-031 ", message, call. = FALSE)
    note("run", "", "", "unavailable", message)
    for (i in seq_len(nrow(single))) add_row("single_de_kegg_pathway_painter",
      single$analysis_id[[i]], "", single$collection[[i]], "skipped_missing_input", "", message)
    for (i in seq_len(nrow(contrasts))) add_row("contrast_kegg_pathway_painter", "",
      contrasts$contrast_id[[i]], contrasts$collection[[i]], "skipped_missing_input", "", message)
    write_lisa_tsv(do.call(rbind, status_rows), status_path)
    return(invisible(FALSE))
  }
  audit <- file.path(output_dir, "audit")
  if (!dir.exists(audit)) lisa_guarded_dir_create(audit)
  scripts <- c("build_single_de_kegg_pathway_painter.R", "build_contrast_kegg_pathway_painter.R",
    "build_contrast_gene_level_product.R")
  code_ledger <- lisa_write_code_identity_ledger(package_dir, scripts,
    file.path(audit, "kegg_maps_code_identity.tsv"))
  trusted_run_root <- lisa_assert_run_tree_safe(output_dir)
  failed <- function(kind, owner, collection, painter_dir, res) {
    message <- paste0("Native KEGG maps could not be rendered from the declared snapshot '",
      snapshot_id, "': ", res$message)
    note(kind, owner, collection, "failed", message)
    add_row(if (kind == "single_de") "single_de_kegg_pathway_painter" else
      "contrast_kegg_pathway_painter", if (kind == "single_de") owner else "",
      if (kind == "contrast") owner else "", collection, "failed",
      if (is.null(res$output_path)) "" else res$output_path, message)
    stop("LISA-KEGG-032 ", kind, " ", owner, " ", collection, ": ", message,
      call. = FALSE)
  }
  for (i in seq_len(nrow(single))) {
    analysis <- as.character(single$analysis_id[[i]])
    collection <- as.character(single$collection[[i]])
    species <- species_contracts$scientific_name[species_contracts$analysis_id == analysis][[1L]]
    gene_root <- file.path(outputs_root, "gene_level", "single_de", analysis,
      paste0("collection_", collection))
    contributions <- file.path(gene_root, paste0(analysis, "_", collection,
      "_gene_level_gene_category_contributions.tsv"))
    if (!file.exists(contributions)) {
      gsea <- list.files(file.path(outputs_root, "single_de", analysis,
        paste0("collection_", collection), "enrichment"),
        pattern = "_GSEA_.*_annotated[.]tsv$", full.names = TRUE)
      if (length(gsea) != 1L) {
        failed("single_de", analysis, collection, gene_root,
          list(message = "No unique saved GSEA table for a completed analysis."))
      }
      built <- lisa_build_single_gene_level_tables(analysis, collection, outputs_root,
        lisa_extension_file_label_prefix(gsea[[1L]], analysis))
      if (!identical(built$status, "completed")) {
        failed("single_de", analysis, collection, gene_root, built)
      }
    }
    res <- lisa_run_post_script(package_dir, "build_single_de_kegg_pathway_painter.R",
      c("--project-dir", output_dir, "--analysis-id", analysis, "--universe", collection,
        "--species", species, "--kegg-cache-root", cache_root,
        "--kegg-snapshot-id", snapshot_id, "--kegg-access-mode", "cache_only"),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger)
    if (identical(res$status, "completed")) {
      note("single_de", analysis, collection, "completed", "")
      add_row("single_de_kegg_pathway_painter", analysis, "", collection, "completed", res$output_path, res$message)
    } else {
      failed("single_de", analysis, collection, file.path(gene_root, "kegg_painter"), res)
    }
  }
  for (i in seq_len(nrow(contrasts))) {
    contrast_id <- as.character(contrasts$contrast_id[[i]])
    collection <- as.character(contrasts$collection[[i]])
    contrast_def <- contrast_index[contrast_index$contrast_id == contrast_id, , drop = FALSE]
    species <- lisa_contrast_species(contrast_def, species_contracts)
    if (length(species) != 1L) {
      note("contrast", contrast_id, collection, "skipped_species_mismatch",
        "KEGG contrast map omitted: the two analyses use different organisms.")
      add_row("contrast_kegg_pathway_painter", "", contrast_id, collection,
        "skipped_species_mismatch", "", "KEGG contrast map omitted: the two analyses use different organisms.")
      next
    }
    owner <- lisa_contrast_name(contrast_def)
    gene_root <- file.path(outputs_root, "gene_level", "category_contrasts", owner,
      paste0("collection_", collection))
    tables <- list.files(gene_root, pattern = "_contrast_gene_level_paired[.]tsv$|paired_gene_evidence[.]tsv$")
    if (!length(tables)) {
      prep <- lisa_run_post_script(package_dir, "build_contrast_gene_level_product.R",
        c("--project-dir", output_dir, "--contrast-id", contrast_id, "--universe", collection,
          "--tables-only", "true"),
        trusted_run_root = trusted_run_root, code_ledger = code_ledger)
      if (!identical(prep$status, "completed")) {
        failed("contrast", contrast_id, collection, gene_root, prep)
      }
    }
    res <- lisa_run_post_script(package_dir, "build_contrast_kegg_pathway_painter.R",
      c("--project-dir", output_dir, "--contrast-id", contrast_id, "--universe", collection,
        "--species", species[[1L]], "--kegg-cache-root", cache_root,
        "--kegg-snapshot-id", snapshot_id, "--kegg-access-mode", "cache_only"),
      trusted_run_root = trusted_run_root, code_ledger = code_ledger)
    if (identical(res$status, "completed")) {
      note("contrast", contrast_id, collection, "completed", "")
      add_row("contrast_kegg_pathway_painter", "", contrast_id, collection, "completed", res$output_path, res$message)
    } else {
      failed("contrast", contrast_id, collection,
        file.path(gene_root, "contrast_kegg_pathway_painter"), res)
    }
  }
  write_lisa_tsv(if (length(status_rows)) do.call(rbind, status_rows) else
    data.frame(kind = character(), owner = character(), collection = character(),
      status = character(), message = character()), status_path)
  invisible(TRUE)
}
