# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé. GPL-3 or later.

lisa_inference_one_file <- function(directory, pattern) {
  paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(paths) != 1L) stop("Category inference requires exactly one ", pattern, " in ", directory, call. = FALSE)
  paths[[1L]]
}

lisa_inference_resource_check <- function(path, inventory, role, require_checksum = FALSE) {
  if (!file.exists(path)) stop("Category inference resource is missing: ", path, call. = FALSE)
  sha <- lisa_sha256_file(path)
  expected <- inventory$sha256[inventory$role == role]
  if (isTRUE(require_checksum) && length(expected) != 1L)
    stop("Saved category inference requires one verified resource checksum: ", role, call. = FALSE)
  if (length(expected) && (length(expected) != 1L || !identical(as.character(expected), sha)))
    stop("Category inference resource differs from the saved run: ", role, call. = FALSE)
  sha
}

# Internal writer: called only in an active new run, or by the explicit
# copy-on-augment API below. Never follow archived output_dir absolute paths.
lisa_run_category_inference <- function(run_dir, term2gene, single_status = NULL,
    plot_formats = "png", source_origin = "current_run") {
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  if (is.null(single_status)) single_status <- read_lisa_tsv(file.path(run_dir, "single_de_status.tsv"))
  lisa_require_columns(single_status, c("analysis_id", "collection", "status"), "single DE status")
  de_index <- read_lisa_tsv(file.path(run_dir, "config", "de_index.tsv"))
  registry <- read_lisa_tsv(file.path(run_dir, "lisa_collection_registry.tsv"))
  declared <- expand.grid(analysis_id = as.character(de_index$analysis_id),
    collection = setdiff(as.character(registry$analysis_collection), "HALLMARKS"),
    stringsAsFactors = FALSE)
  single_status <- single_status[single_status$collection != "HALLMARKS", , drop = FALSE]
  key <- function(z) paste(z$analysis_id, z$collection, sep = "\r")
  if (anyDuplicated(key(single_status)) || !setequal(key(single_status), key(declared)))
    stop("Category inference status rows do not match the complete declared analysis/collection family.", call. = FALSE)
  if (!nrow(single_status)) return(invisible(data.frame()))
  if (any(!single_status$status %in% c("completed", "skipped_existing")))
    stop("Category inference requires all declared analysis/collection results; incomplete families are not allowed.", call. = FALSE)
  plot_formats <- intersect(plot_formats, c("png", "svg", "pdf"))
  inventory_path <- file.path(run_dir, "resource_inventory.tsv")
  inventory <- if (file.exists(inventory_path)) read_lisa_tsv(inventory_path) else data.frame(role = character(), sha256 = character())
  archived <- !identical(source_origin, "current_run")
  t2g_sha <- lisa_inference_resource_check(term2gene, inventory, "term2gene_resource", require_checksum = archived)
  t2g <- read_lisa_tsv(term2gene)
  lisa_require_columns(t2g, c("gs_name", "gene_symbol"), "TERM2GENE")
  memberships <- lapply(split(toupper(as.character(t2g$gene_symbol)), t2g$gs_name), function(x) sort(unique(x[!is.na(x) & nzchar(x)])))
  status <- list()
  for (analysis in unique(single_status$analysis_id)) {
    selected <- single_status[single_status$analysis_id == analysis, , drop = FALSE]
    family_dir <- file.path(run_dir, "outputs", "category_inference", analysis)
    if (dir.exists(family_dir)) stop("Category inference already exists; use a fresh destination, not overwrite: ", family_dir, call. = FALSE)
    ledgers <- assignments <- catalogues <- list()
    inputs <- list(); ranking_hash <- NULL
    for (i in seq_len(nrow(selected))) {
      collection <- selected$collection[[i]]
      directory <- file.path(run_dir, "outputs", "single_de", analysis, paste0("collection_", collection))
      ledger_path <- lisa_inference_one_file(file.path(directory, "qc"), "_GSEA_universe_ledger[.]tsv$")
      rank_path <- lisa_inference_one_file(file.path(directory, "inputs"), "_ranked_genes[.]tsv$")
      annotation_path <- lisa_inference_one_file(file.path(directory, "enrichment"), "_GSEA_.*_annotated[.]tsv$")
      ranking_qc_path <- lisa_inference_one_file(file.path(directory, "qc"), "_ranking_qc[.]tsv$")
      rank_qc <- read_lisa_tsv(ranking_qc_path)
      ranks <- read_lisa_tsv(rank_path)
      rank_sha <- lisa_sha256_file(rank_path)
      if (is.null(ranking_hash)) ranking_hash <- rank_sha else if (!identical(ranking_hash, rank_sha))
        stop("Collections within one DE analysis have different ranked universes; cannot share one deduplicated family.", call. = FALSE)
      ledger <- read_lisa_tsv(ledger_path)
      lisa_require_columns(ledger, c("pathway", "eligible_for_gsea", "ranked_gene_count", "pval", "NES"), "saved GSEA ledger")
      if (any(!ledger$pathway %in% names(memberships)))
        stop("Saved GSEA universe contains gene sets absent from the verified TERM2GENE.", call. = FALSE)
      effective <- lapply(memberships[ledger$pathway], intersect, y = as.character(ranks$symbol))
      if (any(lengths(effective) != ledger$ranked_gene_count))
        stop("Effective gene-set membership differs from the saved GSEA universe.", call. = FALSE)
      ledger$collection <- collection
      ledger$effective_members_sha256 <- vapply(effective, lisa_cache_digest, character(1))
      ledger$ranked_universe_sha256 <- rank_sha
      # The lisaR fgsea wrapper uses scoreType=std, gseaParam=1. Versions and
      # original collection/seed remain in the archived runtime receipts.
      ledger$score_type <- "std"
      ledger$gsea_param <- 1
      ledger$hypothesis_id <- vapply(seq_len(nrow(ledger)), function(j)
        lisa_cache_digest(list(analysis, ledger$pathway[[j]], rank_sha,
          ledger$effective_members_sha256[[j]], "std", 1)), character(1))
      ledgers[[i]] <- ledger[c("collection", "pathway", "hypothesis_id", "eligible_for_gsea", "pval", "NES",
        "effective_members_sha256", "ranked_universe_sha256", "score_type", "gsea_param")]
      palette <- read_lisa_tsv(file.path(directory, "qc", "lisa_category_palette.tsv"))
      palette <- palette[!is.na(palette$category_id) & nzchar(palette$category_id) & palette$category_id != "OTHER_UNCLASSIFIED", , drop = FALSE]
      name_col <- intersect(c("category_display_name", "display_name"), names(palette))
      if (!length(name_col)) stop("Category palette has no display names.", call. = FALSE)
      catalogue <- data.frame(category_id = palette$category_id,
        category_display_name = palette[[name_col[[1L]]]], stringsAsFactors = FALSE)
      catalogue$collection <- collection
      catalogues[[i]] <- catalogue
      member_path <- file.path(directory, "qc", "lisa_category_membership.tsv")
      if (file.exists(member_path)) {
        assignment <- read_lisa_tsv(member_path)
        dictionary_sha <- lisa_sha256_file(member_path)
      } else {
        dictionary_path <- as.character(rank_qc$lisa_dictionary_path[[1]])
        dictionary_sha <- lisa_inference_resource_check(dictionary_path, inventory, "dictionary_resource", require_checksum = archived)
        dictionary <- read_lisa_tsv(dictionary_path)
        assignment <- dictionary[dictionary$universe == collection & dictionary$category_id %in% palette$category_id,
          c("gene_set_id", "category_id"), drop = FALSE]
        if ("n_lisa_gene_sets" %in% names(rank_qc) &&
            length(unique(assignment$gene_set_id)) != rank_qc$n_lisa_gene_sets[[1L]])
          stop("Cannot reconstruct the original selected dictionary scope; a saved category-membership table is required.", call. = FALSE)
        # Saved annotations anchor reconstructed mapping to the original run.
        annotation <- read_lisa_tsv(annotation_path)
        # Historical TSVs use an empty field for unclassified gene sets.
        # They remain in the testing family but are not category assignments.
        observed <- unique(annotation[!is.na(annotation$category_id) & nzchar(annotation$category_id) & annotation$category_id != "OTHER_UNCLASSIFIED", c("pathway", "category_id")])
        expected <- unique(assignment[assignment$gene_set_id %in% ledger$pathway, , drop = FALSE])
        if (!setequal(paste(observed$pathway, observed$category_id, sep = "\r"), paste(expected$gene_set_id, expected$category_id, sep = "\r")))
          stop("Saved and reconstructed category memberships disagree.", call. = FALSE)
      }
      assignment <- unique(assignment[c("gene_set_id", "category_id")])
      names(assignment)[[1L]] <- "pathway"
      assignment$collection <- collection
      assignments[[i]] <- assignment
      inputs[[i]] <- data.frame(collection = collection,
        ledger_file = file.path(paste0("collection_", collection), "qc", basename(ledger_path)),
        ledger_sha256 = lisa_sha256_file(ledger_path), annotation_sha256 = lisa_sha256_file(annotation_path),
        ranked_universe_sha256 = rank_sha, dictionary_membership_sha256 = dictionary_sha,
        stringsAsFactors = FALSE)
    }
    result <- lisa_category_inference(do.call(rbind, ledgers), do.call(rbind, assignments), do.call(rbind, catalogues))
    result$categories$analysis_id <- analysis
    result$categories$family_id <- paste0(analysis, ":", substr(lisa_cache_digest(result$family), 1, 16))
    result$categories <- lisa_category_support_presentation(result$categories)
    lisa_guarded_dir_create(family_dir, run_dir)
    lisa_write_inference_tsv(result$raw_results, file.path(family_dir, "raw_results.tsv"))
    lisa_write_inference_tsv(result$family, file.path(family_dir, "hypothesis_family.tsv"))
    lisa_write_inference_tsv(result$categories, file.path(family_dir, "category_results.tsv"))
    write_lisa_tsv(do.call(rbind, inputs), file.path(family_dir, "input_provenance.tsv"))
    metadata <- list(schema_version = "1.0", analysis_id = analysis, source_origin = source_origin,
      alpha = .05, simultaneous_confidence = .95, family_scope = "within_analysis_across_configured_LISA_category_collections; HALLMARKS_excluded",
      collections = as.character(selected$collection), family_size = nrow(result$family),
      method = "robust Hommel closed testing; hommel(simes=FALSE)", hommel_version = as.character(utils::packageVersion("hommel")),
      lisaR_version = as.character(utils::packageVersion("lisaR")), R_version = R.version.string,
      code_sha256 = lisa_cache_digest(list(category_inference = body(lisa_category_inference),
        closed_bounds = body(lisa_inference_closed_bounds), pipeline = body(lisa_run_category_inference),
        figure_source = body(lisa_category_inference_source), plot = body(lisa_plot_category_inference))),
      term2gene_sha256 = t2g_sha, ranked_universe_sha256 = ranking_hash,
      duplicate_rule = "One null per pathway + effective membership + ranked universe + scoreType/std + gseaParam/1; maximum repeated raw P",
      missing_rule = "Eligible missing raw P replaced by 1 before maximum; retained in family. Pretest-ineligible sets excluded, not dropped from coverage.",
      coverage_rule = "Total fixed dictionary members; analysed = eligible members with a finite raw P, regardless of padj or NES availability.",
      interpretation = "At least d member-set nulls are false; not a count of independent mechanisms or a directional category test.",
      comparison_rule = "Each original DE analysis has its own family. Comparisons between DE profiles remain descriptive; no new contrast P values.",
      bounds_rule = "Invert public localtest partial-conjunction P values using <= alpha consistently; exact-zero shortcut. Avoid hommel1.8 discoveries rounding disagreement at exact machine thresholds.",
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    jsonlite::write_json(metadata, file.path(family_dir, "method.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    jsonlite::write_json(lisa_support_presentation_metadata(), file.path(family_dir, "support_presentation.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    for (collection in selected$collection) {
      directory <- file.path(run_dir, "outputs", "single_de", analysis, paste0("collection_", collection))
      table <- result$categories[result$categories$collection == collection, , drop = FALSE]
      stem <- paste0(safe_file_label(analysis), "_LISA_category_inference")
      table_path <- file.path(directory, "lisa_tables", paste0(stem, ".tsv"))
      lisa_write_inference_tsv(table, table_path)
      lisa_write_support_lollipops(directory, table, analysis, collection, plot_formats)
      source <- lisa_category_inference_source(table,
        result$members[result$members$collection == collection, , drop = FALSE],
        paste(lisa_de_label(de_index, analysis), collection, "enrichment support", sep = " | "))
      rendered <- character()
      if (nrow(source)) {
        source_path <- file.path(directory, "plots", paste0(stem, "_source.tsv"))
        lisa_write_inference_tsv(source, source_path)
        for (format in plot_formats) {
          path <- file.path(directory, "plots", paste0(stem, ".", format))
          lisa_category_inference_save_plot(source, path)
          rendered <- c(rendered, file.path("plots", basename(path)))
        }
      }
      row <- lisa_post_status_row("single_de_category_inference", analysis_id = analysis, collection = collection,
        status = "completed", output_path = table_path,
        message = paste("Robust closed testing; family size", nrow(result$family), ";", sum(table$significant), "significant categories;", length(rendered), "figures"))
      status[[length(status) + 1L]] <- row
    }
  }
  result_status <- do.call(rbind, status)
  write_lisa_tsv(result_status, file.path(run_dir, "category_inference_status.tsv"))
  invisible(result_status)
}

#' Add category inference to a new copy of saved LISA results
#'
#' Reuses raw fgsea results without rerunning enrichment. The source run is
#' untouched; `output_dir` must not exist. Rebuild the report in the new copy
#' with the installed `scripts/build_LISA_report.R` after this operation.
#' New `run_lisa()` runs include this layer automatically.
#' @param source_run Existing complete LISA run directory.
#' @param output_dir New destination directory; never an existing run.
#' @param term2gene Exact TERM2GENE resource used by the source run. Its saved
#'   checksum is verified. When `NULL`, use the saved pipeline metadata path.
#' @param plot_formats Formats to generate: `png`, `svg`, `pdf`.
#' @return Invisibly, the new run directory. Existing fgsea provenance is
#'   preserved and the new inference has its own input hashes and method receipt.
#' @export
add_lisa_category_inference <- function(source_run, output_dir, term2gene = NULL,
    plot_formats = c("png", "svg")) {
  source_run <- lisa_assert_run_tree_safe(source_run)
  proposed <- lisa_path_canonical(output_dir)
  if (identical(proposed, source_run) || startsWith(proposed, paste0(source_run, "/")))
    stop("Category-inference destination must be outside the source run.", call. = FALSE)
  if (file.exists(output_dir) || dir.exists(output_dir)) stop("Category-inference destination must be new.", call. = FALSE)
  lisa_managed_destination(output_dir, create_parent = TRUE)
  if (is.null(term2gene)) {
    metadata <- read_lisa_tsv(file.path(source_run, "lisa_pipeline_metadata.tsv"))
    term2gene <- as.character(metadata$value[metadata$key == "term2gene"])
  }
  if (length(term2gene) != 1L || !file.exists(term2gene)) stop("Provide the verified source TERM2GENE resource.", call. = FALSE)
  manifest_path <- file.path(source_run, "run_manifest.tsv")
  if (!file.exists(manifest_path)) stop("Saved-run category inference requires the original complete run manifest.", call. = FALSE)
  source_check <- lisa_validate_built_run_manifest(source_run, read_lisa_tsv(manifest_path))
  if (!identical(source_check$gate, "PASS"))
    stop("Source run failed integrity verification: ", paste(source_check$findings, collapse = "; "), call. = FALSE)
  inventory <- read_lisa_tsv(file.path(source_run, "resource_inventory.tsv"))
  lisa_inference_resource_check(term2gene, inventory, "term2gene_resource", require_checksum = TRUE)
  if (dir.exists(file.path(source_run, "outputs", "category_inference")))
    stop("Source already contains category inference; do not overwrite an existing analysis.", call. = FALSE)
  fs::dir_copy(source_run, output_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  archive <- file.path(output_dir, "source_provenance")
  lisa_guarded_dir_create(archive, output_dir)
  if (!file.rename(file.path(output_dir, "run_manifest.tsv"), file.path(archive, "original_run_manifest.tsv")))
    stop("Cannot archive inherited manifest before deriving a new run.", call. = FALSE)
  lisa_rebase_run_text_paths(output_dir, source_run, output_dir)
  receipt <- data.frame(operation = "add_category_inference_without_rerunning_fgsea",
    status_at_creation = "inference_added_report_rebuild_pending", source_integrity = source_check$gate,
    source_run = basename(source_run), source_manifest_sha256 = if (file.exists(file.path(source_run, "run_manifest.tsv"))) lisa_sha256_file(file.path(source_run, "run_manifest.tsv")) else NA_character_,
    source_pipeline_metadata_sha256 = lisa_sha256_file(file.path(source_run, "lisa_pipeline_metadata.tsv")),
    stringsAsFactors = FALSE)
  write_lisa_tsv(receipt, file.path(output_dir, "category_inference_derivation.tsv"))
  lisa_run_category_inference(output_dir, term2gene, plot_formats = plot_formats,
    source_origin = "saved_fgsea_results_in_new_copy; original DE/GSEA provenance retained")
  invisible(output_dir)
}

# Only the explicit copy-on-augment route needs a new post-render manifest.
# Native run_lisa() remains owned by its existing transaction coordinator.
lisa_finalize_category_inference_derivation <- function(run_dir) {
  receipt_path <- file.path(run_dir, "category_inference_derivation.tsv")
  if (!file.exists(receipt_path)) return(invisible(FALSE))
  receipt <- read_lisa_tsv(receipt_path)
  if (!"operation" %in% names(receipt) || nrow(receipt) != 1L ||
      receipt$operation != "add_category_inference_without_rerunning_fgsea")
    stop("Invalid category inference derivation receipt.", call. = FALSE)
  # Keep the creation receipt immutable: it may already have a portable copy
  # in Downloads. Record the completed rendering separately, without making
  # an exported pending-status table stale by changing its canonical source.
  finished <- data.frame(operation = receipt$operation,
    status = "report_rebuilt_and_resealed", original_DE_GSEA_recomputed = FALSE,
    lisaR_version = as.character(utils::packageVersion("lisaR")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE), stringsAsFactors = FALSE)
  write_lisa_tsv(finished, file.path(run_dir, "audit", "category_inference_report_receipt.tsv"))
  manifest <- lisa_run_manifest(run_dir, workers = 1L)
  write_lisa_tsv(manifest, file.path(run_dir, "run_manifest.tsv"))
  check <- lisa_validate_built_run_manifest(run_dir, manifest)
  if (!identical(check$gate, "PASS")) stop("Derived report failed final manifest verification.", call. = FALSE)
  invisible(TRUE)
}
