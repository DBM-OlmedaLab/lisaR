# Complete presentation products from saved DE/enrichment evidence only.
# Deliberately separate from the scientific pipeline: no DE or enrichment calls.
lisa_write_native_extended_recipe <- function(data, cfg, png, script, type) {
  builder <- system.file("scripts", script, package = "lisaR")
  if (!nzchar(builder)) stop("Native extended builder is not installed.")
  cfg$project_dir <- NULL
  cfg$lisa_internal_renderer_sha256 <- NULL
  if (!nrow(data)) data <- data.frame(empty_reason = "Gene-category contribution table is empty.")
  data$figure_type <- type
  data$figure_id <- tools::file_path_sans_ext(basename(png))
  data$source_row_order <- seq_len(nrow(data))
  data$selected_for_plot <- TRUE
  data$highlighted <- FALSE
  data$labelled <- TRUE
  data$native_builder_sha256 <- lisa_sha256_file(builder)
  data$native_config_json <- jsonlite::toJSON(cfg, auto_unbox = TRUE, null = "null")
  stem <- tools::file_path_sans_ext(png)
  write_lisa_tsv(data, paste0(stem, "_source.tsv"))
  lisa_copy_verified_figure_recipe(system.file("scripts/reproduce_lisa_figure.R", package = "lisaR"),
    paste0(stem, "_recipe.R"), run_root = dirname(png))
}

lisa_complete_full_products <- function(project_dir, kegg_maps = FALSE,
    kegg_cache_root = "", kegg_snapshot_id = "", kegg_access_mode = "cache_only",
    formats = NULL, families = NULL) {
  root <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  # C3: the formats are the run's requested formats, never a fixed png,pdf.
  if (is.null(formats)) formats <- lisa_requested_plot_formats(root)
  formats <- intersect(c("png", "svg", "pdf"), tolower(as.character(formats)))
  requested_formats <- formats
  if (!"png" %in% formats) formats <- c("png", formats)
  format_arg <- paste(formats, collapse = ",")
  package <- lisa_resolve_package_dir()
  scripts <- file.path(package, "scripts")
  if (!dir.exists(scripts)) scripts <- file.path(package, "inst", "scripts")
  read <- function(p) read_lisa_tsv(p)
  de <- read(file.path(root, "config/de_index.tsv"))
  contrasts <- read(file.path(root, "config/contrast_index.tsv"))
  audit <- file.path(root, "audit/full_products")
  dir.create(audit, recursive = TRUE, showWarnings = FALSE)
  pngs <- function(p, pattern = "[.]png$") lisa_plot_files(p, pattern)
  if (kegg_maps && (!nzchar(kegg_cache_root) || !nzchar(kegg_snapshot_id)))
    stop("FULL maps require an explicit KEGG cache and snapshot.", call. = FALSE)
  if (!kegg_access_mode %in% c("cache_only", "external"))
    stop("Unknown KEGG access mode.", call. = FALSE)
  rows <- list()
  record <- function(owner, collection, family, status, files) {
    rows[[length(rows) + 1L]] <<- data.frame(owner = owner, collection = collection,
      family = family, status = status, n_files = length(files), stringsAsFactors = FALSE)
    write_lisa_tsv(do.call(rbind, rows), file.path(audit, "coverage.tsv"))
  }
  run <- function(script, args, key) {
    # Fixed call sites below form the allowlist: only saved-evidence builders.
    log <- file.path(audit, paste0(gsub("[^A-Za-z0-9_.-]", "_", key), ".log"))
    message("FULL missing product: ", key)
    code <- system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(file.path(scripts, script)), shQuote(args)),
      stdout = log, stderr = log)
    if (code != 0L) {
      details <- if (file.exists(log)) utils::tail(readLines(log, warn = FALSE), 40L) else "Log unavailable."
      stop("FULL product failed: ", key, "; inspect ", log, "\n",
        paste(details, collapse = "\n"), call. = FALSE)
    }
  }
  family <- function(owner, collection, name, targets, script, args,
      empty_receipt = NULL, image_pattern = NULL) {
    # An optional allowlist restricts completion to named families, so a caller
    # that already rendered some products elsewhere does not duplicate them.
    if (!is.null(families) && !name %in% families) return(invisible(character()))
    if (is.null(image_pattern)) image_pattern <- if (name == "kegg_maps") "_painted[.]png$" else "[.]png$"
    images <- function() unique(unlist(lapply(targets, pngs, pattern = image_pattern), use.names = FALSE))
    files <- images()
    # Bind completion receipts to saved scientific inputs and the native builder.
    input_roots <- c(file.path(root, "config"), file.path(root, "outputs/single_de"),
      file.path(root, "outputs/category_contrasts"))
    input_files <- sort(unique(unlist(lapply(input_roots, function(p)
      list.files(p, pattern = "[.]tsv$", recursive = TRUE, full.names = TRUE)))))
    fingerprint <- digest::digest(list(inputs = setNames(vapply(input_files, lisa_sha256_file,
      character(1)), substring(input_files, nchar(root) + 2L)),
      builder = lisa_sha256_file(file.path(scripts, script)), args = args), algo = "sha256")
    key <- paste(owner, collection, name, sep = "__")
    marker <- file.path(audit, paste0(key, ".tsv"))
    pending <- paste0(marker, ".pending")
    if (file.exists(pending)) stop("Interrupted FULL family; preserve existing files and inspect log: ", key, call. = FALSE)
    if (file.exists(marker)) {
      ledger <- read(marker)
      if (!"input_fingerprint" %in% names(ledger) ||
          !identical(unique(ledger$input_fingerprint), fingerprint))
        stop("Completed FULL family input or builder changed: ", key, call. = FALSE)
      if (nrow(ledger) && (any(!file.exists(file.path(root, ledger$path))) ||
          any(vapply(file.path(root, ledger$path), lisa_sha256_file, character(1)) != ledger$sha256)))
        stop("Previously completed FULL family changed: ", key, call. = FALSE)
      record(owner, collection, name, "reused_verified", files)
      return(invisible(files))
    }
    if (length(files)) {
      record(owner, collection, name, "existing_not_regenerated", files)
      return(invisible(files))
    }
    writeLines(key, pending)
    run(script, args, key)
    files <- images()
    if (!length(files) && (is.null(empty_receipt) || !file.exists(empty_receipt)))
      stop("FULL builder produced no product or biological empty receipt: ", key, call. = FALSE)
    products <- unique(unlist(lapply(targets, function(p) list.files(p,
      recursive = TRUE, full.names = TRUE)), use.names = FALSE))
    products <- products[!dir.exists(products)]
    # Temporary rendering formats are removed by the final format policy.
    # Do not bind a reusable completion receipt to those transient files.
    ext <- tolower(tools::file_ext(products))
    products <- products[!ext %in% c("png", "svg", "pdf") | ext %in% requested_formats |
      grepl("_kegg_base[.]png$", products)]
    write_lisa_tsv(data.frame(path = substring(products, nchar(root) + 2L),
      sha256 = vapply(products, lisa_sha256_file, character(1)),
      input_fingerprint = fingerprint), marker)
    unlink(pending)
    record(owner, collection, name, if (length(files)) "generated_missing" else "empty_with_evidence", files)
    invisible(files)
  }
  map_args <- function(species) c("--species", species, "--kegg-cache-root", kegg_cache_root,
    "--kegg-snapshot-id", kegg_snapshot_id, "--kegg-access-mode", kegg_access_mode)
  for (analysis in as.character(de$analysis_id)) {
    dirs <- sort(list.dirs(file.path(root, "outputs/single_de", analysis), recursive = FALSE))
    for (dir in dirs) {
      collection <- sub("^collection_", "", basename(dir))
      gene <- file.path(root, "outputs/gene_level/single_de", analysis, basename(dir))
      prefix <- paste(analysis, collection, "gene_level", sep = "_")
      contributions <- file.path(gene, paste0(prefix, "_gene_category_contributions.tsv"))
      if (!file.exists(contributions)) {
        candidates <- list.files(file.path(dir, "enrichment"), pattern = "_GSEA_.*_annotated[.]tsv$", full.names = TRUE)
        if (length(candidates) != 1L) stop("Ambiguous saved GSEA input: ", dir, call. = FALSE)
        label <- sub("_annotated[.]tsv$", "", substring(basename(candidates), nchar(analysis) + 7L))
        built <- lisa_build_single_gene_level_tables(analysis, collection, file.path(root, "outputs"), label)
        if (!identical(built$status, "completed")) stop(built$message, call. = FALSE)
      }
      args <- c("--project-dir", root, "--analysis-id", analysis, "--universe", collection)
      ext <- file.path(root, "artifacts", analysis, collection)
      family(analysis, collection, "prioritized_genes", c(file.path(gene, "category_gene_cards"), file.path(ext, "gene_cards")),
        "build_single_de_category_gene_cards.R", c(args, "--formats", format_arg))
      family(analysis, collection, "volcano", c(file.path(gene, "category_volcano_overlays"), file.path(ext, "volcano")),
        "build_single_de_category_volcano_overlays.R", c(args, "--formats", format_arg))
      family(analysis, collection, "gene_heatmap", c(file.path(gene, "leading_edge_gene_heatmaps"), file.path(ext, "heatmap")),
        "build_single_de_leading_edge_gene_heatmaps.R", c(args, "--formats", format_arg, "--scale", "zscore"))
      if (collection != "HALLMARKS") family(analysis, collection, "recurrent_genes", file.path(gene, "recurrent_gene_screen"),
        "build_single_de_recurrent_gene_screen.R", args)
      if (kegg_maps) family(analysis, collection, "kegg_maps", file.path(gene, "kegg_painter"),
        "build_single_de_kegg_pathway_painter.R", c(args, map_args(as.character(de$species[de$analysis_id == analysis]))),
        file.path(gene, "kegg_painter", paste0(prefix, "_kegg_pathway_painter_index.tsv")))
    }
  }
  for (i in seq_len(nrow(contrasts))) {
    row <- contrasts[i, , drop = FALSE]
    owner <- paste(row$contrast_id, row$output_id, sep = "_")
    dirs <- sort(list.dirs(file.path(root, "outputs/category_contrasts", owner), recursive = FALSE))
    species <- unique(as.character(de$species[de$analysis_id %in% c(row$contrast_a, row$contrast_b)]))
    for (dir in dirs) {
      collection <- sub("^collection_", "", basename(dir))
      gene <- file.path(root, "outputs/gene_level/category_contrasts", owner, basename(dir))
      args <- c("--project-dir", root, "--contrast-id", as.character(row$contrast_id), "--universe", collection)
      family(owner, collection, "paired_gene_products", gene,
        "build_contrast_gene_level_product.R", args,
        image_pattern = "(_contrast_category_card|_paired_gene_heatmap)[.]png$")
      family(owner, collection, "gene_category_network", file.path(gene, "contrast_gene_category_network"),
        "build_contrast_gene_category_network.R", args)
      if (kegg_maps) {
        if (length(species) != 1L) stop("KEGG contrast requires one consistent species.", call. = FALSE)
        family(owner, collection, "kegg_maps", file.path(gene, "contrast_kegg_pathway_painter"),
          "build_contrast_kegg_pathway_painter.R", c(args, map_args(species)))
      }
    }
  }
  invisible(do.call(rbind, rows))
}
