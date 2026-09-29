# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Internal support for the simple-argument form of run_lisa(). Everything here
# is configuration assembly, input adaptation and refusal to guess. No engine,
# threshold, ranking orientation or scientific value is reinterpreted, and no
# new public entry point is introduced: run_lisa() remains the only door.

lisa_run_simple_arguments <- function() {
  c("de", "output_dir", "species", "positive_direction", "method", "columns",
    "contrast", "contrast_label", "contrast_positive_direction", "collections",
    "mode", "png", "svg", "pdf", "tables", "recipes", "title", "workers",
    "dry_run", "resources", "pipeline", "report")
}

lisa_run_text <- function(value, field) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop("LISA-RUN-ARGS-005 `", field, "` must be one non-empty string.",
         call. = FALSE)
  }
  trimws(value)
}

lisa_run_flag <- function(value, field) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("LISA-RUN-ARGS-006 `", field, "` must be TRUE or FALSE, not ",
         paste(utils::capture.output(utils::str(value)), collapse = " "), ".",
         call. = FALSE)
  }
  value
}

lisa_run_count <- function(value, field) {
  if ((!is.numeric(value) && !is.integer(value)) || length(value) != 1L ||
      is.na(value) || value != floor(value) || value < 1) {
    stop("LISA-RUN-ARGS-007 `", field, "` must be one positive whole number.",
         call. = FALSE)
  }
  as.integer(value)
}

# Analysis identifiers travel into directory names, file names and report
# anchors, so they follow the same portable-ASCII rule the study schema already
# applies. Human-readable wording belongs in `label`.
lisa_run_check_id <- function(id, field) {
  id <- lisa_run_text(id, field)
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", id)) {
    stop("LISA-RUN-ARGS-008 `", field, "` must be a portable technical ",
         "identifier: a letter followed by letters, digits or underscores. ",
         "Got: ", id, ". Put human-readable wording in `title` or a label.",
         call. = FALSE)
  }
  id
}

# ---------------------------------------------------------------------------
# Differential-expression inputs
# ---------------------------------------------------------------------------

lisa_run_read_de_table <- function(path, field) {
  path <- lisa_run_text(path, field)
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "rds")) {
    stop(
      "LISA-RUN-DE-001 `", field, "` names a serialized .rds file. ",
      "Deserializing an .rds is a local trust decision that run_lisa() will ",
      "not take for you: make it explicitly with run_lisa_de(input = ",
      shQuote(path), ", trusted_rds = TRUE).",
      call. = FALSE
    )
  }
  if (!file.exists(path) || dir.exists(path)) {
    stop("LISA-RUN-DE-002 differential-expression table does not exist: ",
         path, call. = FALSE)
  }
  data <- if (identical(extension, "csv")) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (extension %in% c("tsv", "txt", "tab")) {
    read_lisa_tsv(path)
  } else {
    stop(
      "LISA-RUN-DE-003 unsupported table extension '.", extension,
      "'. Use .tsv, .txt, .tab or .csv, or pass the data frame itself.",
      call. = FALSE
    )
  }
  if (!is.data.frame(data) || !nrow(data)) {
    stop("LISA-RUN-DE-002 differential-expression table is empty: ", path,
         call. = FALSE)
  }
  list(data = data, source = path)
}

# DESeq2 and edgeR results are ordinary tables wrapped in an S4 container. The
# conversion below is deliberately short and explicit: it moves identifiers out
# of the row names and nothing else. It never invents the sign of a ranking and
# never infers which side of a contrast is positive.
lisa_run_coerce_de_object <- function(value, field, columns = list()) {
  if (is.data.frame(value)) {
    return(list(data = value, source = "", method = NULL, columns = list()))
  }
  if (is.character(value)) {
    table <- lisa_run_read_de_table(value, field)
    return(list(data = table$data, source = table$source, method = NULL,
                columns = list()))
  }
  classes <- class(value)

  if ("DESeqResults" %in% classes) {
    frame <- as.data.frame(value)
    ids <- rownames(value)
    if (is.null(ids) || !length(ids) ||
        identical(as.character(ids), as.character(seq_len(nrow(frame))))) {
      stop(
        "LISA-RUN-DE-005 this DESeqResults object has no usable feature row ",
        "names. Convert it yourself and supply the identifier column: ",
        "de = cbind(symbol = <ids>, as.data.frame(res)).",
        call. = FALSE
      )
    }
    # An LRT `stat` is a likelihood-ratio statistic: non-negative, and with no
    # defensible sign for an unspecified contrast. Ranking on it would fabricate
    # a direction, so it is refused rather than silently signed.
    # Read the column metadata through the slot rather than through
    # S4Vectors::mcols(), so this check adds no package dependency for a
    # conversion that is otherwise pure base R.
    description <- tryCatch({
      metadata <- methods::slot(value, "elementMetadata")
      as.character(metadata[["description"]])
    }, error = function(e) character())
    stat_note <- description[grepl("statistic", description, ignore.case = TRUE)]
    if (length(stat_note) && any(grepl("LRT|likelihood ratio", stat_note,
                                       ignore.case = TRUE)) &&
        (is.null(columns$rank) || identical(columns$rank, "stat"))) {
      stop(
        "LISA-RUN-DE-006 this DESeqResults object comes from a likelihood-ratio ",
        "test; its `stat` column is unsigned and has no defensible direction ",
        "for an unspecified contrast. Supply a separately computed signed ",
        "ranking with columns = list(rank = \"<your signed statistic>\").",
        call. = FALSE
      )
    }
    frame <- cbind(symbol = as.character(ids), frame, stringsAsFactors = FALSE)
    return(list(data = frame, source = "DESeqResults object",
                method = "DESeq2", columns = list(symbol = "symbol")))
  }

  if (any(c("TopTags", "DGEExact", "DGELRT") %in% classes)) {
    frame <- if ("TopTags" %in% classes) {
      as.data.frame(value$table)
    } else {
      as.data.frame(value$table)
    }
    padj_column <- if (is.null(columns$padj)) "FDR" else columns$padj
    if (!padj_column %in% names(frame)) {
      stop("LISA-RUN-DE-008 this edgeR object has no adjusted-p-value column '",
        padj_column, "'. For a complete edgeR table use edgeR::topTags(result, ",
        "n = Inf, sort.by = \"none\") (BH adjustment by default), then add the ",
        "signed ranking for your one-dimensional comparison. Alternatively, ",
        "provide your existing adjusted p-values and name their column with ",
        "columns = list(padj = \"your_q_value\", rank = \"your_signed_stat\"). ",
        "lisaR does not fit a model or choose a p-value adjustment for you.", call. = FALSE)
    }
    ids <- rownames(frame)
    if (is.null(ids) || !length(ids) ||
        identical(as.character(ids), as.character(seq_len(nrow(frame))))) {
      stop(
        "LISA-RUN-DE-005 this edgeR result has no usable feature row names. ",
        "Convert it yourself and supply the identifier column.",
        call. = FALSE
      )
    }
    frame <- cbind(symbol = as.character(ids), frame, stringsAsFactors = FALSE)
    # edgeR publishes F or LR, both unsigned. lisaR will not sign them from
    # logFC, because a multi-coefficient test has no single direction. The
    # caller declares the signed statistic explicitly or the call fails here.
    return(list(data = frame, source = paste(classes[[1L]], "object"),
                method = "edgeR", columns = list(symbol = "symbol"),
                needs_rank = TRUE))
  }

  stop(
    "LISA-RUN-DE-004 `", field, "` must be a data frame, a .tsv/.txt/.tab/.csv ",
    "path, a DESeq2 DESeqResults object, an edgeR result, or a named list of ",
    "those. Got: ", paste(classes, collapse = "/"), ".",
    call. = FALSE
  )
}

# Turn the `de` argument into a named list of declarations, before anything is
# written. `positive_direction`, `method` and `columns` may be given once for a
# single analysis or per analysis with matching names.
lisa_run_collect_de <- function(de, positive_direction, method, columns) {
  if (is.null(de)) {
    stop("LISA-RUN-ARGS-009 supply `de`: one differential-expression result, ",
         "or a named list of them.", call. = FALSE)
  }
  named_list <- is.list(de) && !is.data.frame(de) && !isS4(de) &&
    !any(c("TopTags", "DGEExact", "DGELRT") %in% class(de))
  items <- if (named_list) de else stats::setNames(list(de), "analysis")
  if (!length(items)) {
    stop("LISA-RUN-ARGS-009 `de` is an empty list.", call. = FALSE)
  }
  ids <- names(items)
  if (named_list) {
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids))) {
      stop("LISA-RUN-ARGS-010 name every element of `de`; the names become the ",
           "analysis identifiers used in paths and in the report.",
           call. = FALSE)
    }
    if (anyDuplicated(ids)) {
      stop("LISA-RUN-ARGS-010 `de` names must be unique; duplicated: ",
           paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
    }
  }
  ids <- vapply(seq_along(items), function(i)
    lisa_run_check_id(ids[[i]], paste0("names(de)[", i, "]")), character(1))

  per_analysis <- function(value, field) {
    if (is.null(value)) return(stats::setNames(vector("list", length(ids)), ids))
    if (is.list(value) && !is.null(names(value))) {
      missing_ids <- setdiff(ids, names(value))
      extra <- setdiff(names(value), ids)
      if (length(extra)) {
        stop("LISA-RUN-ARGS-011 `", field, "` names analyses that are not in ",
             "`de`: ", paste(extra, collapse = ", "), ".", call. = FALSE)
      }
      if (length(missing_ids)) {
        stop("LISA-RUN-ARGS-011 `", field, "` is missing: ",
             paste(missing_ids, collapse = ", "), ". Give one value per ",
             "analysis, or a single value that genuinely applies to all.",
             call. = FALSE)
      }
      return(value[ids])
    }
    stats::setNames(rep(list(value), length(ids)), ids)
  }
  directions <- per_analysis(positive_direction, "positive_direction")
  methods <- per_analysis(method, "method")
  # `columns` is the one argument whose shared form is itself a named list, so
  # list(rank = "stat") could mean one mapping for every analysis or a mapping
  # for an analysis called `rank`. Decide by which vocabulary the names belong
  # to, and refuse rather than guess when they belong to both.
  mappings <- if (is.null(columns)) {
    per_analysis(NULL, "columns")
  } else {
    if (!is.list(columns) || is.null(names(columns)) ||
        any(!nzchar(names(columns)))) {
      stop("LISA-RUN-ARGS-011 `columns` must be a named mapping such as ",
           "list(symbol = \"gene\"), or a named list of such mappings, one ",
           "per analysis.", call. = FALSE)
    }
    canonical <- c("symbol", "logfc", "rank", "pvalue", "padj")
    looks_shared <- all(names(columns) %in% canonical)
    looks_per_analysis <- all(names(columns) %in% ids)
    if (looks_shared && looks_per_analysis) {
      stop("LISA-RUN-ARGS-011 `columns` is ambiguous: its names (",
           paste(names(columns), collapse = ", "), ") are both canonical ",
           "column roles and analysis identifiers. Wrap the mapping once per ",
           "analysis, as columns = list(", ids[[1L]],
           " = list(symbol = \"...\")), to say which you mean.", call. = FALSE)
    }
    if (looks_shared) {
      stats::setNames(rep(list(columns), length(ids)), ids)
    } else if (looks_per_analysis) {
      per_analysis(columns, "columns")
    } else {
      stop("LISA-RUN-ARGS-011 `columns` names (",
           paste(names(columns), collapse = ", "),
           ") are neither canonical column roles (",
           paste(canonical, collapse = ", "), ") nor analyses in `de` (",
           paste(ids, collapse = ", "), ").", call. = FALSE)
    }
  }

  declarations <- lapply(seq_along(ids), function(i) {
    id <- ids[[i]]
    mapping <- mappings[[i]]
    if (is.null(mapping)) mapping <- list()
    coerced <- lisa_run_coerce_de_object(items[[i]], paste0("de$", id), mapping)
    direction <- directions[[i]]
    if (is.null(direction)) {
      stop(
        "LISA-RUN-ARGS-012 declare `positive_direction` for analysis ", id,
        ": lisaR never infers which side of the comparison a positive effect ",
        "belongs to. Example: positive_direction = \"higher in treated than ",
        "in control\".",
        call. = FALSE
      )
    }
    direction <- lisa_run_text(direction, paste0("positive_direction$", id))
    if (!is.null(coerced$columns) && length(coerced$columns)) {
      for (key in names(coerced$columns)) {
        if (is.null(mapping[[key]])) mapping[key] <- coerced$columns[key]
      }
    }
    declared_method <- methods[[i]] %||% coerced$method %||% "DESeq2"
    declared_method <- lisa_run_text(declared_method, paste0("method$", id))
    if (isTRUE(coerced$needs_rank) && is.null(mapping$rank)) {
      stop(
        "LISA-RUN-DE-007 analysis ", id, " is an edgeR result. edgeR reports ",
        "an unsigned F or LR statistic, which has no defensible sign for an ",
        "unspecified contrast, so lisaR will not derive a ranking from it. ",
        "Supply the signed statistic you want to rank on, for example ",
        "columns = list(rank = \"signed_stat\"), having computed it yourself.",
        call. = FALSE
      )
    }
    prepared <- prepare_lisa_de_input(
      coerced$data, method = declared_method, columns = mapping,
      positive_direction = direction
    )
    c(prepared, list(analysis_id = id, source = coerced$source,
                     declared_method = declared_method,
                     positive_direction = direction))
  })
  stats::setNames(declarations, ids)
}

# ---------------------------------------------------------------------------
# Contrasts
# ---------------------------------------------------------------------------

lisa_run_collect_contrasts <- function(contrast, ids, contrast_label,
                                       contrast_positive_direction) {
  if (is.null(contrast)) {
    if (!is.null(contrast_label) || !is.null(contrast_positive_direction)) {
      stop("LISA-RUN-ARGS-013 `contrast_label` and ",
           "`contrast_positive_direction` describe a contrast; declare ",
           "`contrast` as well.", call. = FALSE)
    }
    return(list())
  }
  pairs <- if (is.list(contrast) && !is.character(contrast) &&
               length(contrast) && is.list(contrast[[1L]])) {
    contrast
  } else if (is.list(contrast) && all(vapply(contrast, is.character, logical(1))) &&
             length(contrast) && all(lengths(contrast) == 2L)) {
    contrast
  } else {
    list(contrast)
  }
  labels <- if (is.null(contrast_label)) vector("list", length(pairs)) else
    if (length(pairs) == 1L) list(contrast_label) else as.list(contrast_label)
  wordings <- if (is.null(contrast_positive_direction)) vector("list", length(pairs)) else
    if (length(pairs) == 1L) list(contrast_positive_direction) else
      as.list(contrast_positive_direction)
  if (length(labels) != length(pairs) || length(wordings) != length(pairs)) {
    stop("LISA-RUN-ARGS-013 give one `contrast_label` and one ",
         "`contrast_positive_direction` per contrast, or none.", call. = FALSE)
  }

  lapply(seq_along(pairs), function(i) {
    pair <- unlist(pairs[[i]], use.names = TRUE)
    if (!is.character(pair) || length(pair) != 2L || anyNA(pair)) {
      stop("LISA-RUN-ARGS-014 each contrast names exactly two analyses, ",
           "minuend first: contrast = c(\"ON\", \"PRE\") means ON minus PRE.",
           call. = FALSE)
    }
    if (!is.null(names(pair)) && all(nzchar(names(pair)))) {
      roles <- names(pair)
      if (!setequal(roles, c("minuend", "subtrahend"))) {
        stop("LISA-RUN-ARGS-014 a named contrast must use exactly the roles ",
             "`minuend` and `subtrahend`; got: ",
             paste(roles, collapse = ", "), ".", call. = FALSE)
      }
      pair <- c(pair[["minuend"]], pair[["subtrahend"]])
    }
    minuend <- pair[[1L]]
    subtrahend <- pair[[2L]]
    unknown <- setdiff(c(minuend, subtrahend), ids)
    if (length(unknown)) {
      stop("LISA-RUN-ARGS-015 contrast names analyses that are not in `de`: ",
           paste(unknown, collapse = ", "), ". Available: ",
           paste(ids, collapse = ", "), ".", call. = FALSE)
    }
    if (identical(minuend, subtrahend)) {
      stop("LISA-RUN-ARGS-015 a contrast needs two different analyses; got ",
           minuend, " twice.", call. = FALSE)
    }
    list(
      contrast_id = paste0(minuend, "_vs_", subtrahend),
      analysis_a = minuend,
      analysis_b = subtrahend,
      contrast_label = if (is.null(labels[[i]])) {
        paste(minuend, "versus", subtrahend)
      } else {
        lisa_run_text(labels[[i]], "contrast_label")
      },
      comparison = paste(minuend, "minus", subtrahend),
      # A restatement of the declared orientation, not an inference from it,
      # and deliberately descriptive: a contrast delta is not a category
      # significance test.
      positive_direction = if (is.null(wordings[[i]])) {
        paste0("Positive values indicate a larger mean NES in ", minuend,
               " than in ", subtrahend, ".")
      } else {
        lisa_run_text(wordings[[i]], "contrast_positive_direction")
      }
    )
  })
}

# ---------------------------------------------------------------------------
# Configuration assembly
# ---------------------------------------------------------------------------

lisa_run_simple_report_defaults <- function() {
  list(
    mode = "standard",
    category_evidence = TRUE,
    category_nes_variants = c("clean", "percentages", "direction", "dispersion"),
    legacy_gene_products = FALSE,
    evidence_max_sets = 25L,
    evidence_max_genes = 40L,
    source_data = TRUE,
    recipes = TRUE,
    formats = list(png = TRUE, svg = TRUE, pdf = FALSE)
  )
}

# Resource references stay out of the generated configuration unless the caller
# asked for something specific. Omitting them is what makes the resolution
# internal: lisa_select_pipeline_resource_references() then applies the
# registered defaults for the declared species.
lisa_run_resource_block <- function(resources) {
  if (is.null(resources)) return(list())
  if (is.character(resources) && length(resources) == 1L) {
    tier <- tolower(trimws(resources))
    if (!tier %in% c("core", "expanded")) {
      stop("LISA-RUN-ARGS-016 `resources` as a single string selects the ",
           "dictionary tier and must be \"core\" or \"expanded\". For a custom ",
           "resource use resources = list(dictionary = \"id@version\", ...).",
           call. = FALSE)
    }
    return(list(lisa_dictionary = tier))
  }
  if (!is.list(resources) || is.null(names(resources))) {
    stop("LISA-RUN-ARGS-016 `resources` must be \"core\", \"expanded\", or a ",
         "named list with any of dictionary, term2gene and category_map.",
         call. = FALSE)
  }
  known <- c(dictionary = "dictionary_resource",
             term2gene = "term2gene_resource",
             category_map = "category_map_resource")
  extra <- setdiff(names(resources), names(known))
  if (length(extra)) {
    stop("LISA-RUN-ARGS-016 unknown `resources` entries: ",
         paste(extra, collapse = ", "), ". Known: ",
         paste(names(known), collapse = ", "), ".", call. = FALSE)
  }
  block <- list()
  for (key in names(resources)) {
    block[[known[[key]]]] <- lisa_run_text(resources[[key]],
                                           paste0("resources$", key))
  }
  block
}

lisa_run_merge_overrides <- function(base, overrides, field) {
  if (is.null(overrides)) return(base)
  if (!is.list(overrides) || (length(overrides) &&
      (is.null(names(overrides)) || anyNA(names(overrides)) ||
       any(!nzchar(names(overrides))) || anyDuplicated(names(overrides))))) {
    stop("LISA-RUN-ARGS-017 `", field, "` must be a named list of advanced ",
         "overrides.", call. = FALSE)
  }
  for (key in names(overrides)) {
    value <- overrides[[key]]
    if (is.list(value) && !is.null(names(value)) && is.list(base[[key]])) {
      base[[key]] <- lisa_run_merge_overrides(base[[key]], value,
                                              paste0(field, "$", key))
    } else {
      base[key] <- list(value)
    }
  }
  base
}

# A direct flag and the same key inside `report` are two ways of saying the same
# thing. Choosing one silently would ignore an option the caller supplied, so
# the ambiguity is refused instead.
lisa_run_assert_no_format_conflict <- function(flags, report) {
  if (is.null(report)) return(invisible(TRUE))
  lisa_run_merge_overrides(list(), report, "report")
  conflicts <- character()
  formats <- report$formats
  for (key in c("png", "svg", "pdf")) {
    if (!is.null(flags[[key]]) && is.list(formats) && !is.null(formats[[key]])) {
      conflicts <- c(conflicts, paste0(key, " / report$formats$", key))
    }
  }
  if (!is.null(flags$tables) && !is.null(report$source_data)) {
    conflicts <- c(conflicts, "tables / report$source_data")
  }
  if (!is.null(flags$recipes) && !is.null(report$recipes)) {
    conflicts <- c(conflicts, "recipes / report$recipes")
  }
  if (!is.null(flags$mode) && !is.null(report$mode)) {
    conflicts <- c(conflicts, "mode / report$mode")
  }
  if (length(conflicts)) {
    stop("LISA-RUN-ARGS-004 the same setting was supplied twice, once ",
         "directly and once inside `report`: ",
         paste(conflicts, collapse = "; "),
         ". Remove one of the two; lisaR will not choose for you.",
         call. = FALSE)
  }
  invisible(TRUE)
}

lisa_run_assert_no_pipeline_conflict <- function(title, workers, dry_run,
                                                 resource_block, pipeline) {
  if (is.null(pipeline)) return(invisible(TRUE))
  lisa_run_merge_overrides(list(), pipeline, "pipeline")
  direct <- list(report_title = title, workers = workers, dry_run = dry_run)
  conflicts <- names(direct)[vapply(names(direct), function(key) {
    !is.null(direct[[key]]) && key %in% names(pipeline)
  }, logical(1))]
  # Dictionary tier and dictionary identity select the same resource, so they
  # are conflicting aliases even though their configuration keys differ.
  for (key in names(resource_block)) {
    aliases <- if (key %in% c("lisa_dictionary", "dictionary_resource")) {
      c("lisa_dictionary", "dictionary_resource")
    } else key
    conflicts <- c(conflicts, intersect(aliases, names(pipeline)))
  }
  if (length(conflicts)) {
    stop("LISA-RUN-ARGS-004 the same setting was supplied twice, once ",
      "directly and once inside `pipeline`: ",
      paste(unique(conflicts), collapse = ", "),
      ". Remove one of the two, even if their values agree; lisaR will not choose for you.",
      call. = FALSE)
  }
  invisible(TRUE)
}

# Settings are resolved before any differential-expression table is read, so a
# mistyped flag reports the flag rather than failing later inside the input
# adapter. Nothing here touches the filesystem.
lisa_run_resolve_settings <- function(collections, title, workers, dry_run,
                                      flags, resources, pipeline, report) {
  lisa_run_assert_no_format_conflict(flags, report)
  resource_block <- lisa_run_resource_block(resources)
  lisa_run_assert_no_pipeline_conflict(title, workers, dry_run, resource_block, pipeline)
  report_block <- lisa_run_simple_report_defaults()
  if (!is.null(flags$mode)) {
    mode <- lisa_run_text(flags$mode, "mode")
    if (!mode %in% c("standard", "full")) {
      stop("LISA-RUN-ARGS-018 `mode` must be \"standard\" or \"full\". ",
           "\"selected\" is an extension of a completed run; use ",
           "plan_lisa_extension() and render_lisa_categories() for it.",
           call. = FALSE)
    }
    report_block$mode <- mode
  }
  for (key in c("png", "svg", "pdf")) {
    if (!is.null(flags[[key]])) {
      report_block$formats[[key]] <- lisa_run_flag(flags[[key]], key)
    }
  }
  if (!is.null(flags$tables)) {
    report_block$source_data <- lisa_run_flag(flags$tables, "tables")
  }
  if (!is.null(flags$recipes)) {
    report_block$recipes <- lisa_run_flag(flags$recipes, "recipes")
  }
  report_block <- lisa_run_merge_overrides(report_block, report, "report")

  pipeline_block <- c(
    list(
      schema_version = "1.0.0",
      profile = "transcriptomic/genomic",
      evidence_mode = "full_de",
      output_dir = "results",
      dry_run = if (is.null(dry_run)) FALSE else lisa_run_flag(dry_run, "dry_run"),
      workers = if (is.null(workers)) 4L else lisa_run_count(workers, "workers"),
      run_ora = FALSE,
      run_kegg_maps = FALSE,
      run_hallmarks = FALSE,
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    resource_block
  )
  if (!is.null(title)) {
    pipeline_block$report_title <- lisa_run_text(title, "title")
  }
  pipeline_block <- lisa_run_merge_overrides(pipeline_block, pipeline, "pipeline")

  collections <- if (is.null(collections)) {
    lisa_default_collections()
  } else {
    normalize_lisa_collections(collections)
  }
  list(pipeline = pipeline_block, report = report_block,
       collections = collections)
}

lisa_run_build_config <- function(declarations, relative_paths, species,
                                  contrasts, settings) {
  single_de <- lapply(seq_along(declarations), function(i) {
    declaration <- declarations[[i]]
    list(
      analysis_id = declaration$analysis_id,
      label = declaration$analysis_id,
      comparison = declaration$positive_direction,
      positive_direction = declaration$positive_direction,
      model_note = paste0(
        "Differential expression was computed by the caller and declared as a ",
        declaration$declared_method,
        " table; lisaR adapted it without refitting any model."
      ),
      de_path = relative_paths[[i]],
      species = species,
      symbol_col = "symbol",
      rank_col = "rank_value",
      logfc_col = "log2FoldChange",
      pvalue_col = "pvalue",
      padj_col = "padj"
    )
  })

  config <- list(
    pipeline = settings$pipeline,
    report = settings$report,
    collections = as.list(settings$collections),
    single_de = single_de
  )
  if (length(contrasts)) config$contrasts <- contrasts
  config
}

# ---------------------------------------------------------------------------
# Readiness reporting
# ---------------------------------------------------------------------------

# The stock readiness error names the blocked stage only. An investigator who
# has not written the YAML needs the rows that actually failed, which is
# normally a semantic resource that has not been installed yet.
lisa_run_readiness_error <- function(validation, config_path) {
  blocked <- validation$readiness$stage[
    !is.na(validation$readiness$ready) & !validation$readiness$ready
  ]
  describe <- function(table, stage, id_columns, ready_column = "ready",
                       message_columns = "message") {
    if (!is.data.frame(table) || !nrow(table) ||
        !ready_column %in% names(table)) {
      return(character())
    }
    flag <- table[[ready_column]]
    bad <- table[!is.na(flag) & !flag, , drop = FALSE]
    if (!nrow(bad)) return(character())
    id_columns <- intersect(id_columns, names(bad))
    message_columns <- intersect(message_columns, names(bad))
    vapply(seq_len(nrow(bad)), function(i) {
      label <- paste(vapply(id_columns, function(column)
        as.character(bad[[column]][[i]]), character(1)), collapse = " / ")
      message <- paste(vapply(message_columns, function(column)
        as.character(bad[[column]][[i]]), character(1)), collapse = "; ")
      paste0("  - ", stage, ": ", label,
             if (nzchar(message)) paste0(" -- ", message) else "")
    }, character(1))
  }
  details <- c(
    describe(validation$inputs, "input", c("analysis_id", "column", "path")),
    describe(validation$dependencies, "missing package",
             c("package", "manager"), ready_column = "available",
             message_columns = "reason"),
    describe(validation$resources, "resource", c("config_key", "resource_id"))
  )
  hint <- if (any(grepl("msigdb_term2gene", details, fixed = TRUE))) {
    paste0(
      "\nThe gene-set membership resource is not installed yet. Prepare it ",
      "once with prepare_lisa_msigdb_resource(), which asks you to accept the ",
      "MSigDB terms explicitly; lisaR never accepts them on your behalf.\n"
    )
  } else {
    ""
  }
  stop(
    "LISA-RUN-READY-001 the study is not ready to run; blocked stage(s): ",
    paste(blocked, collapse = ", "), ".\n",
    if (length(details)) paste(c(details, ""), collapse = "\n") else "",
    hint,
    "Nothing was written to the requested output directory, so correcting the ",
    "argument and calling run_lisa() again works without deleting anything.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Report resolution
# ---------------------------------------------------------------------------

lisa_run_report_pages <- function() {
  c(overview = "report_index.html",
    analyses = "report_pages/single_de.html",
    contrasts = "report_pages/contrasts.html",
    genes = "report_pages/gene_evidence/index.html",
    methods = "report_pages/downloads.html")
}

# Report artifacts are resolved by existence, never by constructing a plausible
# path. A field named `contrast_report` that points at a file which was never
# written is worse than an honest absence.
lisa_run_resolve_reports <- function(output_dir, has_contrast, dry_run) {
  pages <- lisa_run_report_pages()
  paths <- file.path(output_dir, unname(pages))
  exists <- file.exists(paths) & !dir.exists(paths)
  reports <- data.frame(
    page = names(pages), file = unname(pages), path = paths,
    exists = exists, stringsAsFactors = FALSE
  )
  if (isTRUE(dry_run)) return(reports)
  overview <- reports$exists[reports$page == "overview"]
  if (!isTRUE(overview)) {
    stop("LISA-RUN-REPORT-001 the run completed but the root HTML report is ",
         "absent: ", reports$path[reports$page == "overview"],
         ". Inspect post_lisa_status.tsv in the run directory.", call. = FALSE)
  }
  if (isTRUE(has_contrast) && !isTRUE(reports$exists[reports$page == "contrasts"])) {
    stop("LISA-RUN-REPORT-001 a contrast was configured but its report page ",
         "was not written: ", reports$path[reports$page == "contrasts"],
         ". Inspect post_lisa_status.tsv in the run directory.", call. = FALSE)
  }
  reports
}

# ---------------------------------------------------------------------------
# The simple form of run_lisa()
# ---------------------------------------------------------------------------

lisa_run_from_simple_arguments <- function(de, output_dir, species,
                                           positive_direction, method, columns,
                                           contrast, contrast_label,
                                           contrast_positive_direction,
                                           collections, flags, title, workers,
                                           dry_run, resources, pipeline,
                                           report) {
  lisa_require_optional("jsonlite", "writing the generated study configuration")
  if (is.null(output_dir)) {
    stop("LISA-RUN-ARGS-019 supply `output_dir`: one new directory that will ",
         "hold the generated configuration, the adapted inputs and the ",
         "results.", call. = FALSE)
  }
  if (is.null(species)) {
    stop("LISA-RUN-ARGS-020 declare `species`, for example ",
         "species = \"Homo sapiens\". lisaR resolves its semantic resources ",
         "against the declared species and never infers it from identifiers.",
         call. = FALSE)
  }
  species <- lisa_run_text(species, "species")
  lisa_run_assert_no_format_conflict(flags, report)

  # Step 1: resolve every scalar setting first. A mistyped flag, an unknown
  # resource selection or an impossible mode is reported as itself, before any
  # differential-expression table is read and before the destination exists.
  settings <- lisa_run_resolve_settings(
    collections = collections, title = title, workers = workers,
    dry_run = dry_run, flags = flags, resources = resources,
    pipeline = pipeline, report = report
  )
  dry_run_flag <- isTRUE(settings$pipeline$dry_run)

  # Step 2: adapt every input in memory. A missing column or an unreadable
  # table fails here, still with nothing written.
  declarations <- lisa_run_collect_de(de, positive_direction, method, columns)
  ids <- names(declarations)
  contrasts <- lisa_run_collect_contrasts(contrast, ids, contrast_label,
                                          contrast_positive_direction)
  relative_paths <- as.list(file.path("inputs", paste0(ids, "_de.tsv")))
  config <- lisa_run_build_config(
    declarations = declarations, relative_paths = relative_paths,
    species = species, contrasts = contrasts, settings = settings
  )

  # Step 3: refuse an occupied destination before writing a single byte.
  destination <- lisa_managed_destination(output_dir, create_parent = TRUE)
  study_root <- destination$path
  parent <- destination$parent
  if (lisa_path_entry_exists(study_root)) {
    stop("LISA-RUN-DEST-001 `output_dir` must be a new directory; lisaR will ",
         "not write into an existing one: ", study_root,
         ". Choose another path, or remove that one deliberately.",
         call. = FALSE)
  }

  # Step 4: materialize the whole study in an exclusive staging directory
  # beside the destination, so validation runs against real files without the
  # destination ever existing.
  staging <- file.path(parent, paste0(".", basename(study_root),
                                      ".staging-", Sys.getpid()))
  if (lisa_path_entry_exists(staging)) {
    stop("LISA-RUN-DEST-002 staging destination already exists: ", staging,
         call. = FALSE)
  }
  staging <- lisa_run_root(staging)
  promoted <- FALSE
  on.exit(if (!promoted && lisa_path_entry_exists(staging)) {
    lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  }, add = TRUE)
  lisa_guarded_dir_create(file.path(staging, "inputs"), staging)

  input_rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    id <- ids[[i]]
    declaration <- declarations[[i]]
    de_path <- file.path(staging, relative_paths[[i]])
    write_lisa_tsv(declaration$data, de_path)
    write_lisa_tsv(
      cbind(analysis_id = id, declaration$receipt, stringsAsFactors = FALSE),
      file.path(staging, "inputs", paste0(id, "_receipt.tsv"))
    )
    write_lisa_tsv(declaration$row_audit,
                   file.path(staging, "inputs", paste0(id, "_row_audit.tsv")))
    write_lisa_tsv(declaration$column_mapping,
                   file.path(staging, "inputs", paste0(id, "_column_mapping.tsv")))
    input_rows[[i]] <- data.frame(
      analysis_id = id,
      source = if (nzchar(declaration$source)) declaration$source else
        "in-memory data frame",
      declared_method = declaration$declared_method,
      positive_direction = declaration$positive_direction,
      effective_rows = declaration$receipt$effective_rows,
      stringsAsFactors = FALSE
    )
  }
  staged_config <- file.path(staging, "study.json")
  lisa_guarded_write(
    staged_config,
    function(path) jsonlite::write_json(
      config, path, auto_unbox = TRUE, pretty = TRUE,
      null = "null", na = "null", digits = 17
    ),
    run_root = staging
  )

  # Step 5-6: validate the materialized study. A failure deletes the staging
  # tree through the on.exit handler and leaves the destination untouched.
  validation <- validate_lisa_config(staged_config, check_files = TRUE)
  if (!isTRUE(validation$execution_ready)) {
    lisa_run_readiness_error(validation, staged_config)
  }

  # Step 7: one rename. lisa_promote_managed_directory() refuses an occupied
  # destination itself, which closes the window opened since step 3.
  lisa_promote_managed_directory(staging, study_root)
  promoted <- TRUE
  study_root <- normalizePath(study_root, winslash = "/", mustWork = TRUE)
  config_path <- file.path(study_root, "study.json")

  # Step 8: the ordinary engine, on the promoted study.
  result <- withCallingHandlers(
    run_lisa_pipeline_from_config(config_path),
    error = function(error) {
      message(
        "The generated study is readable at ", study_root,
        ". Rerun it without rebuilding anything: run_lisa(config = ",
        shQuote(config_path), ")."
      )
    }
  )

  reports <- lisa_run_resolve_reports(result$output_dir,
                                      has_contrast = length(contrasts) > 0L,
                                      dry_run = dry_run_flag)
  invisible(c(result, list(
    study_dir = study_root,
    config_path = config_path,
    inputs = do.call(rbind, input_rows),
    analyses = validation$de_index,
    contrasts = validation$contrast_index,
    reports = reports,
    report_index = reports$path[reports$page == "overview"],
    effective_config = validation$configuration_provenance
  )))
}
