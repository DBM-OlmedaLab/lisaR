# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

#' Run a LISA analysis
#'
#' `run_lisa()` is the single entry point for a complete LISA analysis, and it
#' can be called in two ways.
#'
#' Supply `config` to run an existing YAML or JSON study configuration exactly
#' as before. Supply the simple arguments instead to run a study directly from
#' one or more differential-expression results, without writing a
#' configuration file and without naming any dictionary, category map or
#' gene-set membership resource: those are resolved internally for the declared
#' species. The simple form assembles an ordinary study, writes it where you
#' can read it, and hands it to the same engine. No threshold, method, ranking
#' orientation or scientific value differs between the two forms.
#'
#' The two forms are mutually exclusive. Mixing `config` with any simple
#' argument is an error rather than a silent preference for one of them, and an
#' unknown argument name is an error rather than an ignored option.
#'
#' @param config Path to a YAML or JSON configuration file, or a named list
#'   containing the same configuration fields. Relative paths in a list are
#'   resolved from the current working directory. Leave it `NULL` to use the
#'   simple arguments below.
#' @param ... Must be empty. It is present so that every simple argument is
#'   supplied by name and a second positional argument can never be taken for
#'   one of them.
#' @param de One already computed differential-expression result, or a named
#'   list of them whose names become the analysis identifiers. Each element may
#'   be a data frame, a path to a `.tsv`, `.txt`, `.tab` or `.csv` table, a
#'   DESeq2 `DESeqResults` object, or an edgeR `TopTags`, `DGEExact` or
#'   `DGELRT` result. Serialized `.rds` inputs are refused here because
#'   deserializing one is a local trust decision; make it explicitly with
#'   [run_lisa_de()]. There is no limit of two analyses.
#' @param output_dir A new directory that will hold the generated `study.json`,
#'   the adapted inputs and receipts under `inputs/`, and the pipeline results
#'   under `results/`. An existing directory is refused, never overwritten.
#' @param species Species represented by the input identifiers, for example
#'   `"Homo sapiens"`. It is never inferred.
#' @param positive_direction One non-empty description of what a positive
#'   effect means, for example `"higher in treated than in control"`, or a
#'   named list with one entry per analysis. lisaR never infers which side of a
#'   comparison is positive.
#' @param method Declared table convention: `"DESeq2"`, `"edgeR"` or
#'   `"limma"`, or a named list per analysis. This is a caller declaration, not
#'   verification of the producing software. It defaults to the convention
#'   implied by a supplied S4 object, otherwise to `"DESeq2"`.
#' @param columns Named list mapping `symbol`, `logfc`, `rank`, `pvalue` and
#'   `padj` to source columns, or a named list of such mappings per analysis.
#'   Passed unchanged to [prepare_lisa_de_input()]. An edgeR result requires an
#'   explicit signed `rank`: `F` and `LR` are unsigned and are never signed for
#'   you.
#' @param contrast Two analysis identifiers, minuend first:
#'   `c("ON", "PRE")` means *ON minus PRE*. The roles may also be named, as
#'   `c(minuend = "ON", subtrahend = "PRE")`. A list of such pairs declares
#'   several contrasts. The orientation is never taken from the order of `de`.
#' @param contrast_label,contrast_positive_direction Optional wording for the
#'   contrast. The default restates the declared orientation and claims no
#'   statistical significance for any category.
#' @param collections Gene-set collections to analyse. Defaults to the
#'   registered analysis collections.
#' @param mode `"standard"` (default) or `"full"`. `"full"` writes one report in
#'   the same `output_dir`: the standard report plus every applicable FULL
#'   figure, rendered from the run's own saved results inside the same run
#'   (under `artifacts/`) without recomputing enrichment, and sealed by the same
#'   run manifest. If a FULL product fails, the run fails. It can create
#'   thousands of files.
#' @param png,svg,pdf Figure formats. In this form PNG and SVG default to
#'   `TRUE` and PDF to `FALSE`. A configuration file keeps its own documented
#'   defaults, which are unchanged.
#' @param tables Whether to export the source-data tables. `TRUE` by default
#'   here.
#' @param recipes Whether to export the figure reproduction recipes. `TRUE` by
#'   default here.
#' @param title Optional report title.
#' @param workers Worker processes for the pipeline.
#' @param dry_run Whether to plan the run instead of computing it. A dry run
#'   writes the output plan, performs no enrichment and promises no report.
#' @param resources Optional resource selection. `"core"` or `"expanded"`
#'   chooses the dictionary tier; a named list with any of `dictionary`,
#'   `term2gene` and `category_map` names registered resources explicitly.
#'   Omitted, the registered defaults for the declared species are resolved
#'   internally and you name nothing.
#' @param pipeline,report Named lists of advanced overrides merged over the
#'   defaults of this form before validation. Explicit values win, including
#'   `report = list(formats = list(svg = FALSE))`. Supplying the same setting
#'   both directly and here is refused, even when the values agree.
#' @details `report.category_nes_variants` selects a nonempty unique array of
#'   `clean`, `percentages`, `direction`, and/or `dispersion`, defaulting to all four. These
#'   views retain the mean NES, optionally add median and direction percentages,
#'   and optionally add type-7 P25--P75 dispersion. They do not rerun GSEA.
#'   Study files must declare `pipeline.schema_version: "1.0.0"`.
#'   Earlier schema versions require explicit migration in a copy of the file;
#'   they are not silently upgraded. At least one variant must be selected,
#'   independently of the requested output formats.
#'
#'   In the simple form every checkable argument is validated and every input
#'   adapted in memory before anything is written, and the whole study is
#'   materialized in a staging directory beside `output_dir` and validated
#'   there. A wrong flag or a missing resource therefore leaves `output_dir`
#'   untouched, and correcting the argument and calling again works without
#'   deleting anything.
#'
#' @return Invisibly, a list describing the completed run, including its output
#'   directory, run identifier and validation gate. The simple form adds
#'   `study_dir`, `config_path`, the adapted `inputs`, the validated `analyses`
#'   and `contrasts`, `effective_config`, and `reports`: a table of the report
#'   pages with the path of each and whether it was actually written. Report
#'   paths are resolved by existence, never constructed.
#' @export
#'
#' @examples
#' project <- lisa_init_project(tempfile("lisa-plan-project-"))
#' config <- file.path(project, "study.yml")
#' plan <- yaml::read_yaml(config)
#' plan$pipeline$output_dir <- "results/example"
#' plan$pipeline$dry_run <- TRUE
#' plan$pipeline$sufficiency_action <- "ignore_with_reason"
#' plan$pipeline$sufficiency_ignore_reason <- "Synthetic sample fixture."
#' yaml::write_yaml(plan, config)
#' result <- run_lisa(config)
#' result$gate
#'
#' # The simple form. Running it needs the installed semantic resources, so the
#' # call itself stays interactive.
#' results <- data.frame(
#'   symbol = c("GENE01", "GENE02", "GENE03"),
#'   log2FoldChange = c(1.2, -0.8, 0.2),
#'   stat = c(4.1, -2.6, 0.5),
#'   pvalue = c(1e-4, 0.01, 0.6),
#'   padj = c(1e-3, 0.04, 0.8)
#' )
#' if (interactive()) {
#'   run <- run_lisa(
#'     de = results,
#'     output_dir = tempfile("lisa-study-"),
#'     species = "Homo sapiens",
#'     positive_direction = "higher in treated than in control"
#'   )
#'   run$reports
#' }
run_lisa <- function(config = NULL, ..., de = NULL, output_dir = NULL,
                     species = NULL, positive_direction = NULL, method = NULL,
                     columns = NULL, contrast = NULL, contrast_label = NULL,
                     contrast_positive_direction = NULL, collections = NULL,
                     mode = NULL, png = NULL, svg = NULL, pdf = NULL,
                     tables = NULL, recipes = NULL, title = NULL,
                     workers = NULL, dry_run = NULL, resources = NULL,
                     pipeline = NULL, report = NULL) {
  extra <- names(list(...))
  if (length(list(...))) {
    stop("LISA-RUN-ARGS-003 unknown argument(s): ",
         paste(if (is.null(extra)) "<unnamed>" else extra, collapse = ", "),
         ". Every run_lisa() option is a named argument; nothing supplied is ",
         "ignored. See ?run_lisa.", call. = FALSE)
  }
  supplied <- Filter(Negate(is.null), list(
    de = de, output_dir = output_dir, species = species,
    positive_direction = positive_direction, method = method,
    columns = columns, contrast = contrast, contrast_label = contrast_label,
    contrast_positive_direction = contrast_positive_direction,
    collections = collections, mode = mode, png = png, svg = svg, pdf = pdf,
    tables = tables, recipes = recipes, title = title, workers = workers,
    dry_run = dry_run, resources = resources, pipeline = pipeline,
    report = report
  ))
  if (is.null(config)) {
    if (!length(supplied)) {
      stop("LISA-RUN-ARGS-001 supply either `config`, a YAML/JSON study ",
           "configuration, or the simple arguments `de`, `output_dir`, ",
           "`species` and `positive_direction`. See ?run_lisa.", call. = FALSE)
    }
    return(lisa_run_from_simple_arguments(
      de = de, output_dir = output_dir, species = species,
      positive_direction = positive_direction, method = method,
      columns = columns, contrast = contrast, contrast_label = contrast_label,
      contrast_positive_direction = contrast_positive_direction,
      collections = collections,
      flags = list(mode = mode, png = png, svg = svg, pdf = pdf,
                   tables = tables, recipes = recipes),
      title = title, workers = workers, dry_run = dry_run,
      resources = resources, pipeline = pipeline, report = report
    ))
  }
  if (length(supplied)) {
    stop("LISA-RUN-ARGS-002 `config` describes the whole study, so it cannot ",
         "be combined with the simple argument(s): ",
         paste(names(supplied), collapse = ", "),
         ". Edit the configuration, or drop `config` and use the simple form. ",
         "Neither is applied silently.", call. = FALSE)
  }
  if (is.character(config) && length(config) == 1L && !is.na(config)) {
    validation <- validate_lisa_config(config, check_files = TRUE)
    lisa_assert_execution_ready(validation)
    return(run_lisa_pipeline_from_config(config))
  }
  if (!is.list(config) || is.null(names(config))) {
    stop("`config` must be a YAML/JSON path or a named configuration list.", call. = FALSE)
  }
  lisa_require_optional("jsonlite", "running LISA from an in-memory configuration")
  config_dir <- normalizePath(getwd(), mustWork = TRUE)
  config_path <- tempfile("lisa-config-", tmpdir = config_dir, fileext = ".json")
  on.exit(if (lisa_path_entry_exists(config_path)) {
    lisa_guarded_delete(config_path, run_root = NULL)
  }, add = TRUE)
  lisa_guarded_write(
    config_path,
    function(path) jsonlite::write_json(
      config, path, auto_unbox = TRUE, pretty = TRUE,
      null = "null", na = "null", digits = 17
    ),
    run_root = NULL
  )
  validation <- validate_lisa_config(config_path, check_files = TRUE)
  lisa_assert_execution_ready(validation)
  run_lisa_pipeline_from_config(config_path)
}

lisa_assert_execution_ready <- function(validation) {
  if (isTRUE(validation$execution_ready)) return(invisible(TRUE))
  blocked <- validation$readiness$stage[
    !is.na(validation$readiness$ready) & !validation$readiness$ready
  ]
  stop(
    "LISA-READINESS-001 configuration is not execution-ready; blocked stage(s): ",
    paste(blocked, collapse = ", "),
    ". Inspect validate_lisa_config(...)$readiness and the corresponding status tables.",
    call. = FALSE
  )
}

#' Validate a LISA study configuration
#'
#' Parses and validates a YAML, JSON or in-memory LISA configuration without
#' running enrichment or writing analysis outputs. Relative paths in a file
#' configuration are resolved from the directory containing that file.
#'
#' @param config Path to a YAML or JSON configuration file, or a named list
#'   containing the equivalent structure.
#' @param check_files Whether to verify referenced input files and resolve all
#'   three registered semantic resources against the active species, registry
#'   and managed cache. `FALSE` performs structural prevalidation only and
#'   cannot establish execution readiness.
#' @param strict Whether unknown configuration keys, including nested analysis,
#'   contrast, duplicate-policy and report keys, are errors. Setting this to
#'   `FALSE` only ignores unknown keys during validation; known fields and
#'   security prohibitions are still checked. Execution is always strict.
#'
#' @return A validation summary containing explicit `structural_valid`,
#'   `inputs_ready`, `dependencies_ready`, `resources_ready` and
#'   `execution_ready` states, a stage-level `readiness` table, detailed input,
#'   dependency and resource status tables, validated indexes and the
#'   normalized configuration. `valid` is a backward-compatible alias for
#'   `execution_ready`. Structural contract violations raise an error; missing
#'   files, dependencies or registered resources return a non-ready summary.
#' @export
#'
#' @examples
#' config <- system.file("examples", "minimal-study.json", package = "lisaR")
#' if (nzchar(config)) {
#'   validation <- validate_lisa_config(config)
#'   validation$readiness
#'   validation$execution_ready
#' }
validate_lisa_config <- function(config, check_files = TRUE, strict = TRUE) {
  check_files <- lisa_config_bool(check_files, "check_files")
  resolved <- lisa_resolve_config(config, strict = strict)
  cfg <- resolved$raw_config
  config_dir <- resolved$config_dir
  contract <- resolved$contract
  pipeline <- lisa_config_get(cfg, "pipeline", list())
  de_index <- resolved$de_index
  input_status <- lisa_config_input_readiness(de_index, check_files)
  source_data_status <- lisa_config_source_data_readiness(
    cfg, config_dir, check_files
  )
  policy_input_status <- lisa_duplicate_policy_input_readiness(
    contract$duplicate_policies, config_dir, check_files
  )
  input_status <- rbind(
    input_status, source_data_status, policy_input_status
  )
  row.names(input_status) <- NULL
  inputs_ready <- if (isTRUE(check_files)) {
    nrow(input_status) > 0L && all(input_status$ready)
  } else {
    NA
  }

  contrast_index <- resolved$contrast_index
  species <- unique(as.character(de_index$species))
  species_ready <- length(species) == 1L && !is.na(species[[1L]]) &&
    nzchar(species[[1L]])
  dependency_status <- if (species_ready) {
    lisa_config_dependency_requirements(
      pipeline, species[[1L]], report = contract$report
    )
  } else {
    data.frame()
  }
  dependencies_ready <- species_ready && nrow(dependency_status) > 0L &&
    all(dependency_status$available)

  resource_spec <- lisa_pipeline_resource_spec()
  resource_selection <- lisa_select_pipeline_resource_references(pipeline)
  effective_cfg <- resolved$cfg
  if (isTRUE(check_files) && species_ready) {
    resource_check <- lisa_resolve_pipeline_resources(
      pipeline, species[[1L]], stop_on_error = FALSE
    )
    resource_status <- resource_check$status
    resources_ready <- resource_check$ready
  } else {
    configured_ids <- vapply(
      resource_spec$config_key,
      function(key) resource_selection$references[[key]]$configured_id,
      character(1)
    )
    effective_ids <- vapply(
      resource_spec$config_key,
      function(key) resource_selection$references[[key]]$resource_id,
      character(1)
    )
    selection_sources <- vapply(
      resource_spec$config_key,
      function(key) resource_selection$references[[key]]$selection_source,
      character(1)
    )
    resource_status <- data.frame(
      config_key = resource_spec$config_key,
      configured_id = configured_ids,
      selection_source = selection_sources,
      resource_id = effective_ids,
      species = if (species_ready) species[[1L]] else NA_character_,
      modality = contract$profile,
      schema = resource_spec$expected_schema,
      compatibility = NA_character_,
      sha256 = NA_character_,
      path = NA_character_,
      checked = FALSE,
      ready = NA,
      status = "not_checked",
      message = if (species_ready) {
        "Resource resolution was skipped because check_files = FALSE."
      } else {
        "Resource resolution requires exactly one declared species."
      },
      stringsAsFactors = FALSE
    )
    resources_ready <- if (isTRUE(check_files)) FALSE else NA
  }

  execution_ready <- isTRUE(inputs_ready) &&
    isTRUE(dependencies_ready) && isTRUE(resources_ready)
  readiness <- data.frame(
    stage = c(
      "structure", "inputs", "dependencies", "resources", "execution"
    ),
    checked = c(
      TRUE, check_files, species_ready, check_files && species_ready, TRUE
    ),
    ready = c(
      TRUE, inputs_ready, dependencies_ready, resources_ready, execution_ready
    ),
    status = c(
      "ready",
      if (isTRUE(check_files)) {
        if (isTRUE(inputs_ready)) "ready" else "not_ready"
      } else "not_checked",
      if (!species_ready) {
        "not_checked"
      } else if (isTRUE(dependencies_ready)) {
        "ready"
      } else {
        "not_ready"
      },
      if (isTRUE(check_files) && species_ready) {
        if (isTRUE(resources_ready)) "ready" else "not_ready"
      } else "not_checked",
      if (isTRUE(execution_ready)) "ready" else "not_ready"
    ),
    stringsAsFactors = FALSE
  )

  structure(list(
    valid = execution_ready,
    structural_valid = TRUE,
    inputs_ready = inputs_ready,
    dependencies_ready = dependencies_ready,
    resources_ready = resources_ready,
    execution_ready = execution_ready,
    readiness = readiness,
    schema_version = contract$schema_version,
    gsea_padj_cutoff = contract$gsea_padj_cutoff,
    config_dir = config_dir,
    analyses = nrow(de_index),
    contrasts = nrow(contrast_index),
    de_index = de_index,
    contrast_index = contrast_index,
    inputs = input_status,
    dependencies = dependency_status,
    resources = resource_status,
    config = effective_cfg,
    configuration_provenance = resolved$provenance,
    resolved_config = resolved
  ), class = "lisa_config_validation")
}

#' Run LISA for one differential-expression analysis
#'
#' This lower-level interface runs gene-set enrichment, semantic category
#' annotation and category summaries for one differential-expression result.
#' Most studies should use [run_lisa()] instead.
#'
#' @param input Differential-expression table path, data frame or supported
#'   differential-expression object.
#' @param output_dir Directory in which to write the analysis.
#' @param comparison_name Portable ASCII technical identifier used in output
#'   paths and file names. Human-readable Unicode belongs in a separate title
#'   or label field.
#' @param species Species represented by the input identifiers.
#' @param symbol_col,rank_col,logfc_col,pvalue_col,padj_col Explicit input
#'   column names. `rank_col` is optional when an effect size is available.
#' @param universes Gene-set collections to analyse. The five options are
#'   `"GOBP-C2"`, `"GOMF"`, `"GOCC"`, `"PATHWAYS"` and `"HALLMARKS"`.
#' @param lisa_dictionary Dictionary tier: `"core"` resolves the registered
#'   `lisa_dictionary_core@1.0.0` resource and `"expanded"` resolves
#'   `lisa_dictionary_expanded@1.0.0`. Both dictionaries are bundled with lisaR;
#'   the scientific TERM2GENE membership resource is obtained separately.
#' @param term2gene_path Optional species-matched technical table of gene-set
#'   memberships used to build the lists passed to GSEA.
#' @param make_heatmaps Reserved compatibility flag. In version 1.0 this
#'   lower-level function warns and does not generate heatmaps. Use a configured
#'   [run_lisa()] study with an expression matrix, then request the `heatmap`
#'   product through [render_lisa_categories()].
#' @param ... Named advanced arguments passed to the single-analysis engine;
#'   see the function-reference vignette. Every argument supplied through
#'   `...` must have an explicit name.
#'   `category_nes_variants` selects a nonempty unique subset of
#'   `c("clean", "percentages", "direction", "dispersion")`, default all four for requested
#'   GSEA category plots. Direction adds median and significant-member direction
#'   percentages; dispersion adds type-7 P25--P75, not confidence intervals.
#'   PNG/PDF/SVG sources, settings and recipes accompany the selected views.
#'   Direct HALLMARKS views do not receive semantic category variants.
#' @param trusted_rds Whether this local R call explicitly trusts a local,
#'   stable and fully trusted serialized `.rds` input. It is `FALSE` by default
#'   and cannot be enabled by a YAML/JSON study configuration. This argument is
#'   after `...` and therefore must always be named explicitly.
#' @param rds_max_bytes Maximum serialized-file and approximate in-memory object
#'   size accepted when `trusted_rds = TRUE`. This is a diagnostic limit, not a
#'   sandbox or a hard memory-safety boundary for hostile RDS files. This
#'   argument is after `...` and therefore must always be named explicitly.
#'
#' @return A list containing `comparison_name`, `output_dir`, ranking QC, GSEA
#'   and ORA tables, category summaries, and trusted-input receipt details when
#'   applicable.
#' @export
#'
#' @examples
#' # The installed sample data make this complete example offline. The
#' # calculation is left for an interactive session because it creates a full
#' # set of scientific tables and figures.
#' if (interactive()) {
#'   de_path <- system.file(
#'     "examples", "quick-start", "data", "de_a.tsv",
#'     package = "lisaR"
#'   )
#'   resource <- function(file) system.file(
#'     "extdata", "quick-start", file, package = "lisaR"
#'   )
#'   de_result <- run_lisa_de(
#'     input = de_path,
#'     output_dir = tempfile("lisa-de-example-"),
#'     comparison_name = "response_a",
#'     species = "Homo sapiens",
#'     symbol_col = "symbol",
#'     rank_col = "stat",
#'     logfc_col = "log2FoldChange",
#'     pvalue_col = "pvalue",
#'     padj_col = "padj",
#'     universes = "GOBP-C2",
#'     lisa_dictionary = "sample",
#'     term2gene_path = resource("example_term2gene_v2_1.tsv"),
#'     lisa_dictionary_path = resource("example_dictionary_v3_0.tsv"),
#'     category_map_path = resource("example_category_map_v1_1.tsv"),
#'     export_formats = "tsv",
#'     verbose = FALSE
#'   )
#'   de_result$output_dir
#'   names(de_result$summaries)
#' }
run_lisa_de <- function(input, output_dir, comparison_name = NULL,
                        species = "Homo sapiens", symbol_col = NULL,
                        rank_col = NULL, logfc_col = NULL, pvalue_col = NULL,
                        padj_col = NULL,
                        universes = lisa_default_collections(),
                        lisa_dictionary = "core", term2gene_path = NULL,
                        make_heatmaps = FALSE, ...,
  trusted_rds = FALSE, rds_max_bytes = 512 * 1024^2) {
  dots_call <- match.call(expand.dots = FALSE)$...
  if (is.null(dots_call) || !length(dots_call)) {
    dot_names <- character()
  } else {
    dot_names <- names(dots_call)
    if (is.null(dot_names) || length(dot_names) != length(dots_call) ||
        anyNA(dot_names) || any(!nzchar(dot_names))) {
      lisa_rds_abort(
        "LISA-RDS-016",
        paste0(
          "every advanced argument supplied through `...` must be named ",
          "explicitly."
        )
      )
    }
  }
  private_formals <- c("input_receipt", ".receipt_object_verified")
  private_matches <- pmatch(
    dot_names %||% character(), private_formals,
    nomatch = 0L, duplicates.ok = TRUE
  )
  private_controls <- dot_names[private_matches > 0L]
  if (length(private_controls)) {
    lisa_rds_abort(
      "LISA-RDS-016",
      paste0(
        "internal trusted-RDS controls cannot be supplied through `...`: ",
        paste(private_controls, collapse = ", "), "."
      )
    )
  }
  run_LISA_DE(
    input = input, output_dir = output_dir,
    comparison_name = comparison_name, species = species,
    symbol_col = symbol_col, rank_col = rank_col, logfc_col = logfc_col,
    pvalue_col = pvalue_col, padj_col = padj_col, universes = universes,
    lisa_dictionary = lisa_dictionary, term2gene_path = term2gene_path,
    trusted_rds = trusted_rds, rds_max_bytes = rds_max_bytes,
    make_heatmaps = make_heatmaps, ...
  )
}

#' Compare two completed LISA analyses
#'
#' This lower-level interface compares category summaries from two compatible
#' LISA analyses. Most studies should declare contrasts in the configuration
#' supplied to [run_lisa()].
#'
#' @param contrast_a,contrast_b Paths to two compatible completed
#'   single-analysis result directories.
#' @param output_dir Directory in which to write the contrast.
#' @param comparison_name Contrast identifier.
#' @param contrast_a_label,contrast_b_label Labels shown in tables and plots.
#' @param universes Gene-set collections to compare. The five options are
#'   `"GOBP-C2"`, `"GOMF"`, `"GOCC"`, `"PATHWAYS"` and `"HALLMARKS"`.
#' @param ... Additional presentation and export arguments passed to the
#'   contrast engine.
#'   `category_nes_variants` selects a nonempty unique subset of
#'   `c("clean", "percentages", "direction", "dispersion")`, default all four. Distributions
#'   are separate for A/B and the existing delta is unchanged. Exact annotated
#'   GSEA member tables are required; legacy summary-only inputs retain their
#'   plots and explicitly record unavailable distribution variants.
#'
#' @return A list containing `comparison_name`, `output_dir`, `universe`, the
#'   contrast `summary`, and plot/input manifests.
#' @export
#'
#' @examples
#' # A contrast consumes two completed single-analysis directories. The
#' # sample project is fully offline; execution is interactive because the
#' # prerequisite analyses create scientific results and figures.
#' if (interactive()) {
#'   project <- lisa_init_project(tempfile("lisa-contrast-project-"))
#'   result <- run_lisa(file.path(project, "study.yml"))
#'   analyses <- utils::read.delim(
#'     file.path(result$output_dir, "single_de_status.tsv"),
#'     check.names = FALSE,
#'     stringsAsFactors = FALSE
#'   )
#'   contrast <- run_lisa_contrast(
#'     contrast_a = analyses$output_dir[analyses$analysis_id == "response_a"],
#'     contrast_b = analyses$output_dir[analyses$analysis_id == "response_b"],
#'     output_dir = tempfile("lisa-contrast-example-"),
#'     comparison_name = "response_a_vs_b",
#'     contrast_a_label = "Response A",
#'     contrast_b_label = "Response B",
#'     universes = "GOBP-C2",
#'     export_formats = "tsv",
#'     verbose = FALSE
#'   )
#'   contrast$output_dir
#'   head(contrast$summary)
#' }
run_lisa_contrast <- function(contrast_a, contrast_b, output_dir,
                              comparison_name = NULL,
                              contrast_a_label = "Contrast A",
                              contrast_b_label = "Contrast B",
                              universes = lisa_default_collections(), ...) {
  run_LISA_contrast(
    contrast_a = contrast_a, contrast_b = contrast_b,
    output_dir = output_dir, comparison_name = comparison_name,
    contrast_a_label = contrast_a_label,
    contrast_b_label = contrast_b_label, universes = universes, ...
  )
}

#' Check a completed LISA run
#'
#' @param run_dir Path to a completed LISA results folder.
#'
#' @return A result containing `gate`, `findings` and the path to the file list
#'   used for verification. An incomplete, changed or planning artifact returns
#'   `gate = "FAIL"` with explanatory findings; an unsafe or unreadable path can
#'   still raise an error before verification.
#' @export
#'
#' @examples
#' project <- lisa_init_project(tempfile("lisa-verify-project-"))
#' config <- file.path(project, "study.yml")
#' plan <- yaml::read_yaml(config)
#' plan$pipeline$output_dir <- "results/example"
#' plan$pipeline$dry_run <- TRUE
#' plan$pipeline$sufficiency_action <- "ignore_with_reason"
#' plan$pipeline$sufficiency_ignore_reason <- "Synthetic sample fixture."
#' yaml::write_yaml(plan, config)
#' planned <- run_lisa(config)
#' check <- verify_lisa_run(planned$plan_dir)
#' check$gate
#' check$findings
verify_lisa_run <- function(run_dir) {
  verify_run(run_dir)
}
