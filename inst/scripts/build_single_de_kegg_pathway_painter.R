#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(system.file("scripts", "kegg_snapshot_helpers.R", package = "lisaR"), local = TRUE)
source(system.file("scripts", "lisa_plot_metadata.R", package = "lisaR"), local = TRUE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing TSV: %s", path), call. = FALSE)
  if (file.info(path)$size == 0L) return(data.frame())
  first_lines <- readLines(path, n = 2L, warn = FALSE)
  if (length(first_lines) == 0L || !any(nzchar(trimws(first_lines)))) return(data.frame())
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

parse_args <- function(args) {
  out <- list(
    project_dir = NA_character_,
    analysis_id = NA_character_,
    universe = "GOBP-C2",
    species = NA_character_,
    kegg_cache_root = NA_character_,
    kegg_snapshot_id = NA_character_,
    kegg_access_mode = "external",
    top_pathways = 20,
    max_abs_log2fc = 0.5,
    color_power = 1.0,
    # --- exact pathway selection (H2) -------------------------------------
    # These three are opt-in. Without them this script behaves exactly as it
    # always has: rank the whole collection, take the top --top-pathways and
    # paint them.
    #
    #   --kegg-id            paint exactly ONE pathway, chosen from the SAME
    #                        collection-wide ranked summary, keeping its own
    #                        rank, title, node values and colour context.
    #   --emit-pathway-index write that ranked summary to a named path so an
    #                        interface can offer the selector without the
    #                        generator having to run again.
    #   --index-only         stop after the summary. Enumerating what could be
    #                        painted must not paint anything.
    kegg_id = NA_character_,
    emit_pathway_index = NA_character_,
    index_only = "false",
    lisa_internal_renderer_sha256 = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("top_pathways")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("max_abs_log2fc", "color_power")) {
      out[[key]] <- as.numeric(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("LISA-KEGG-011 required argument --project-dir is missing. Repair: provide the validated run root for this analysis.", call. = FALSE)
  if (is.na(out$analysis_id) || out$analysis_id == "") stop("Required argument: --analysis-id", call. = FALSE)
  if (is.na(out$species) || out$species == "") stop("Required argument: --species", call. = FALSE)
  if (is.na(out$kegg_cache_root) || out$kegg_cache_root == "") stop("Required argument: --kegg-cache-root", call. = FALSE)
  if (is.na(out$kegg_snapshot_id) || out$kegg_snapshot_id == "") stop("Required argument: --kegg-snapshot-id", call. = FALSE)
  if (!out$kegg_access_mode %in% c("external", "cache_only")) stop("Invalid --kegg-access-mode; use external or cache_only.", call. = FALSE)
  if (!out$index_only %in% c("true", "false")) stop("--index-only must be true or false.", call. = FALSE)
  out$index_only <- identical(out$index_only, "true")
  if (!is.na(out$kegg_id) && out$kegg_id != "" && !grepl("^[a-z]{3}[0-9]{5}$", out$kegg_id)) {
    stop(sprintf("LISA-KEGG-027 --kegg-id must be one KEGG pathway identifier such as hsa04110; received: %s.", out$kegg_id), call. = FALSE)
  }
  if (!is.na(out$kegg_id) && out$kegg_id != "" && out$index_only) {
    stop("LISA-KEGG-028 --kegg-id and --index-only true are mutually exclusive: one paints a map, the other refuses to paint anything.", call. = FALSE)
  }
  contract <- kegg_species_contract(out$species)
  out <- c(out, contract)
  out
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

paths_for <- function(project_dir, analysis_id, universe) {
  collection_dir <- file.path(project_dir, "outputs", "single_de", analysis_id, paste0("collection_", universe))
  gene_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis_id, paste0("collection_", universe))
  prefix <- paste(analysis_id, universe, "gene_level", sep = "_")
  painter_dir <- file.path(gene_dir, "kegg_painter")
  list(
    collection_dir = collection_dir,
    gene_dir = gene_dir,
    de = file.path(collection_dir, "inputs", paste0(analysis_id, "_standardized_DE.tsv")),
    evidence = file.path(gene_dir, paste0(prefix, "_gene_category_contributions.tsv")),
    painter_dir = painter_dir,
    painter_data_dir = file.path(painter_dir, "painter_data"),
    prefix = prefix
  )
}

load_de_table <- function(path) {
  de <- read_tsv(path)
  symbol_col <- if ("symbol" %in% names(de)) "symbol" else if ("hgnc_symbol" %in% names(de)) "hgnc_symbol" else ""
  lfc_col <- if ("log2FoldChange" %in% names(de)) "log2FoldChange" else if ("log2FC" %in% names(de)) "log2FC" else ""
  padj_col <- if ("padj" %in% names(de)) "padj" else if ("FDR" %in% names(de)) "FDR" else ""
  gene_id_col <- if ("gene_id" %in% names(de)) "gene_id" else if ("original_gene_id" %in% names(de)) "original_gene_id" else ""
  if (symbol_col == "" || lfc_col == "") stop(sprintf("Standardized DE table lacks symbol/log2FC columns: %s", path), call. = FALSE)
  out <- data.frame(
    symbol = safe_chr(de[[symbol_col]]),
    gene_id = if (gene_id_col != "") safe_chr(de[[gene_id_col]]) else "",
    log2FC = safe_num(de[[lfc_col]]),
    padj = if (padj_col != "") safe_num(de[[padj_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  out[out$symbol != "" & is.finite(out$log2FC), , drop = FALSE]
}

map_symbols_to_entrez <- function(de_table, cfg) {
  if (!"gene_id" %in% names(de_table)) de_table$gene_id <- ""
  de_table$input_symbol <- safe_chr(de_table$symbol)
  de_table$gene_id <- safe_chr(de_table$gene_id)
  de_table$symbol_upper <- toupper(de_table$input_symbol)

  annotation_db <- kegg_annotation_database(cfg)
  mouse_symbols <- AnnotationDbi::keys(annotation_db, keytype = "SYMBOL")
  symbol_ref <- suppressMessages(AnnotationDbi::select(
    annotation_db,
    keys = mouse_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ))
  symbol_ref <- symbol_ref[!is.na(symbol_ref$ENTREZID), , drop = FALSE]
  symbol_ref$symbol_upper <- toupper(safe_chr(symbol_ref$SYMBOL))

  merged <- merge(symbol_ref, de_table, by = "symbol_upper", all.x = FALSE, all.y = FALSE)
  ensembl_keys <- unique(de_table$gene_id[grepl("^ENS(MUS)?G", de_table$gene_id)])
  if (length(ensembl_keys) > 0) {
    ensembl_ref <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        annotation_db,
        keys = ensembl_keys,
        keytype = "ENSEMBL",
        columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
      )),
      error = function(e) data.frame(SYMBOL = character(), ENTREZID = character(), ENSEMBL = character())
    )
    ensembl_ref <- ensembl_ref[!is.na(ensembl_ref$ENTREZID), , drop = FALSE]
    if (nrow(ensembl_ref) > 0) {
      merged_ensembl <- merge(ensembl_ref, de_table, by.x = "ENSEMBL", by.y = "gene_id", all.x = FALSE, all.y = FALSE)
      merged_ensembl$gene_id <- safe_chr(merged_ensembl$ENSEMBL)
      merged_ensembl$symbol_upper <- toupper(safe_chr(merged_ensembl$SYMBOL))
      merged <- rbind(merged, merged_ensembl[, names(merged), drop = FALSE])
    }
  }
  if (nrow(merged) > 0) {
    merged <- merged[!duplicated(paste(merged$ENTREZID, merged$input_symbol, merged$log2FC, merged$padj, sep = "\r")), , drop = FALSE]
  }
  merged$ENTREZID <- safe_chr(merged$ENTREZID)
  merged$SYMBOL <- safe_chr(merged$input_symbol)
  merged
}

gene_to_pathway_map <- function(cfg) {
  links <- kegg_snapshot_rds(cfg, "pathway_links", "all")
  data.frame(
    kegg_gene = sub("^[a-z]{3}:", "", names(links)),
    kegg_id = sub("^path:", "", unname(links)),
    stringsAsFactors = FALSE
  )
}

pathway_titles <- function(cfg) {
  klist <- kegg_snapshot_rds(cfg, "pathway_list", "all")
  species_suffix <- if (cfg$kegg_code == "mmu") " - Mus musculus \\(house mouse\\)$" else " - Homo sapiens \\(human\\)$"
  data.frame(
    kegg_id = sub("^path:", "", names(klist)),
    kegg_title = sub(species_suffix, "", unname(klist)),
    stringsAsFactors = FALSE
  )
}

build_pathway_summary <- function(evidence, de_entrez, cfg) {
  if (nrow(evidence) == 0 || nrow(de_entrez) == 0) return(list(summary = data.frame(), gene_map = data.frame()))
  evidence$symbol_upper <- toupper(safe_chr(evidence$symbol))
  evidence$gene_contribution_score_num <- safe_num(evidence$gene_contribution_score)
  evidence$padj_num <- safe_num(evidence$padj)
  evidence$log2FC_num <- safe_num(evidence$log2FC)

  mapped <- de_entrez[, c("symbol_upper", "ENTREZID", "SYMBOL", "log2FC", "padj"), drop = FALSE]
  mapped <- mapped[!duplicated(mapped$symbol_upper), , drop = FALSE]
  ev <- merge(evidence, mapped, by = "symbol_upper", all.x = FALSE, all.y = FALSE)
  if (nrow(ev) == 0) return(list(summary = data.frame(), gene_map = data.frame()))

  g2p <- gene_to_pathway_map(cfg)
  evp <- merge(ev, g2p, by.x = "ENTREZID", by.y = "kegg_gene", all.x = FALSE, all.y = FALSE)
  if (nrow(evp) == 0) return(list(summary = data.frame(), gene_map = data.frame()))
  titles <- pathway_titles(cfg)
  evp <- merge(evp, titles, by = "kegg_id", all.x = TRUE)

  split_rows <- split(evp, evp$kegg_id)
  summary <- lapply(split_rows, function(df) {
    data.frame(
      analysis_id = cfg$analysis_id,
      universe = cfg$universe,
      kegg_id = df$kegg_id[[1]],
      kegg_title = df$kegg_title[[1]],
      n_evidence_genes = length(unique(df$SYMBOL)),
      max_gene_contribution_score = suppressWarnings(max(safe_num(df$gene_contribution_score_num), na.rm = TRUE)),
      min_de_padj = suppressWarnings(min(safe_num(df$padj), na.rm = TRUE)),
      mean_abs_log2FC = mean(abs(safe_num(df$log2FC)), na.rm = TRUE),
      genes = paste(sort(unique(df$SYMBOL)), collapse = "/"),
      categories = paste(sort(unique(safe_chr(df$category_id))), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary)
  summary$max_gene_contribution_score[is.infinite(summary$max_gene_contribution_score)] <- NA_real_
  summary$min_de_padj[is.infinite(summary$min_de_padj)] <- NA_real_
  summary$mean_abs_log2FC[is.nan(summary$mean_abs_log2FC)] <- NA_real_
  summary <- summary[order(-summary$n_evidence_genes, -safe_num(summary$max_gene_contribution_score), summary$min_de_padj, summary$kegg_title), , drop = FALSE]
  summary$rank <- seq_len(nrow(summary))
  summary <- summary[, c("rank", setdiff(names(summary), "rank")), drop = FALSE]
  list(summary = summary, gene_map = evp)
}

attr_value <- function(txt, attr) {
  m <- regexec(paste0(attr, "=\"([^\"]*)\""), txt, perl = TRUE)
  out <- regmatches(txt, m)
  vapply(out, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
}

extract_entry_blocks <- function(kgml) {
  # readLines() returns one element per XML line.  Matching only the first
  # element silently produced zero nodes and therefore unpainted KEGG maps.
  kgml <- paste(kgml, collapse = "\n")
  m <- gregexpr("<entry\\b[\\s\\S]*?</entry>", kgml, perl = TRUE)
  blocks <- regmatches(kgml, m)[[1]]
  if (identical(blocks, character(0))) character() else blocks
}

parse_kegg_gene_nodes <- function(kgml) {
  blocks <- extract_entry_blocks(kgml)
  if (length(blocks) == 0) return(data.frame())
  type <- attr_value(blocks, "type")
  blocks <- blocks[type == "gene"]
  if (length(blocks) == 0) return(data.frame())
  entry_name <- attr_value(blocks, "name")
  graphics <- regmatches(blocks, regexpr("<graphics\\b[^>]*/?>", blocks, perl = TRUE))
  graphics[lengths(graphics) == 0] <- NA_character_
  ids <- lapply(strsplit(entry_name, "\\s+"), function(v) gsub("^[a-z]{3}:", "", v))
  data.frame(
    entry_id = attr_value(blocks, "id"),
    kegg_entry_name = entry_name,
    kegg_gene_ids = vapply(ids, paste, collapse = "/", FUN.VALUE = character(1)),
    x = safe_num(attr_value(graphics, "x")),
    y = safe_num(attr_value(graphics, "y")),
    width = safe_num(attr_value(graphics, "width")),
    height = safe_num(attr_value(graphics, "height")),
    kegg_label = attr_value(graphics, "name"),
    stringsAsFactors = FALSE
  )
}

best_hit_for_ids <- function(de_entrez, ids) {
  hits <- de_entrez[de_entrez$ENTREZID %in% ids, , drop = FALSE]
  if (nrow(hits) == 0) return(list(symbol = "", log2FC = NA_real_, padj = NA_real_))
  hits$abs_lfc <- abs(safe_num(hits$log2FC))
  hits <- hits[order(-hits$abs_lfc, safe_num(hits$padj), hits$SYMBOL), , drop = FALSE]
  list(symbol = paste(unique(safe_chr(hits$SYMBOL)), collapse = "/"), log2FC = safe_num(hits$log2FC)[[1]], padj = safe_num(hits$padj)[[1]])
}

node_values_for_pathway <- function(nodes, de_entrez, kegg_id) {
  if (nrow(nodes) == 0) return(data.frame())
  out <- lapply(seq_len(nrow(nodes)), function(i) {
    ids <- unlist(strsplit(nodes$kegg_gene_ids[[i]], "/", fixed = TRUE), use.names = FALSE)
    hit <- best_hit_for_ids(de_entrez, ids)
    data.frame(
      nodes[i, , drop = FALSE],
      symbol = hit$symbol,
      log2FC = hit$log2FC,
      padj = hit$padj,
      kegg_id = kegg_id,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

lfc_to_rgb <- function(x, max_abs, color_power = 1) {
  if (is.na(x)) return(grDevices::col2rgb("#E8E8E8")[, 1] / 255)
  z <- max(-1, min(1, x / max_abs))
  z <- sign(z) * (abs(z) ^ color_power)
  red <- grDevices::col2rgb("#B2182B")[, 1] / 255
  blue <- grDevices::col2rgb("#2166AC")[, 1] / 255
  white <- c(1, 1, 1)
  if (z >= 0) white * (1 - z) + red * z else white * (1 - abs(z)) + blue * abs(z)
}

paint_node_backgrounds <- function(img, node_values, max_abs, color_power) {
  out <- img
  h <- dim(out)[1]
  w <- dim(out)[2]
  for (i in seq_len(nrow(node_values))) {
    if (is.na(node_values$x[[i]]) || is.na(node_values$y[[i]])) next
    xmin <- max(1, floor(node_values$x[[i]] - node_values$width[[i]] / 2))
    xmax <- min(w, ceiling(node_values$x[[i]] + node_values$width[[i]] / 2))
    ymin <- max(1, floor(node_values$y[[i]] - node_values$height[[i]] / 2))
    ymax <- min(h, ceiling(node_values$y[[i]] + node_values$height[[i]] / 2))
    if (xmin > xmax || ymin > ymax) next
    fill_rgb <- lfc_to_rgb(node_values$log2FC[[i]], max_abs, color_power)
    alpha <- if (is.na(node_values$log2FC[[i]])) 0.30 else 0.95
    region <- out[ymin:ymax, xmin:xmax, 1:3, drop = FALSE]
    luminance <- 0.299 * region[, , 1] + 0.587 * region[, , 2] + 0.114 * region[, , 3]
    preserve <- luminance < 0.42
    for (ch in seq_len(3)) {
      channel <- region[, , ch]
      channel[!preserve] <- channel[!preserve] * (1 - alpha) + fill_rgb[[ch]] * alpha
      region[, , ch] <- channel
    }
    out[ymin:ymax, xmin:xmax, 1:3] <- region
  }
  out
}

draw_painted_map <- function(img, node_values, title, subtitle, out_png, out_pdf, max_abs, color_power) {
  h <- dim(img)[1]
  w <- dim(img)[2]
  top_pad <- 42
  bottom_pad <- 28
  total_h <- h + top_pad + bottom_pad
  painted_img <- paint_node_backgrounds(img, node_values, max_abs, color_power)
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = c(0, w), yscale = c(0, total_h)))
    grid::grid.rect(x = w / 2, y = total_h / 2, width = w, height = total_h, default.units = "native", gp = grid::gpar(fill = "white", col = NA))
    grid::grid.raster(painted_img, x = w / 2, y = bottom_pad + h / 2, width = w, height = h, default.units = "native", interpolate = FALSE)
    for (i in seq_len(nrow(node_values))) {
      if (is.na(node_values$x[[i]])) next
      grid::grid.rect(x = node_values$x[[i]], y = bottom_pad + h - node_values$y[[i]], width = node_values$width[[i]], height = node_values$height[[i]], default.units = "native", gp = grid::gpar(fill = NA, col = "#202020", lwd = 0.45))
    }
    grid::grid.text(title, x = 10, y = total_h - 12, just = c("left", "top"), default.units = "native", gp = grid::gpar(fontsize = 10, fontface = "bold", col = "#111111"))
    grid::grid.text(subtitle, x = 10, y = total_h - 28, just = c("left", "top"), default.units = "native", gp = grid::gpar(fontsize = 8, col = "#111111"))
    grid::grid.text(sprintf("Node fill: log2FC capped at +/- %.2g, color power %.2g. Red=up, blue=down, grey=no mapped value. Source pathway diagram: KEGG (www.kegg.jp).", max_abs, color_power), x = 10, y = 12, just = c("left", "bottom"), default.units = "native", gp = grid::gpar(fontsize = 8, col = "#111111"))
    grid::popViewport()
  }
  grDevices::png(out_png, width = w, height = total_h, units = "px", bg = "white")
  draw()
  grDevices::dev.off()
  grDevices::pdf(out_pdf, width = w / 96, height = total_h / 96, bg = "white")
  draw()
  grDevices::dev.off()
  # D4: the painted map as SVG, with the KEGG raster embedded by cairo and the
  # node outlines and text as vectors. Removed later if SVG was not requested.
  tryCatch({
    grDevices::svg(sub("[.]png$", ".svg", out_png), width = w / 96, height = total_h / 96, bg = "white")
    draw()
    grDevices::dev.off()
  }, error = function(e) warning("SVG painted map not written: ", conditionMessage(e), call. = FALSE))
}

safe_title <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

paint_one_pathway <- function(row, de_entrez, cfg, out_dir) {
  kgml <- readLines(kegg_snapshot_resource(cfg, "kgml", row$kegg_id), warn = FALSE)
  img <- png::readPNG(kegg_snapshot_resource(cfg, "image", row$kegg_id))
  nodes <- parse_kegg_gene_nodes(kgml)
  if (!nrow(nodes)) {
    stop(sprintf("LISA-KEGG-017 KGML %s contains zero parsed gene nodes; refuse to emit an unpainted map.", row$kegg_id), call. = FALSE)
  }
  node_values <- node_values_for_pathway(nodes, de_entrez, row$kegg_id)
  if (!any(!is.na(node_values$log2FC))) {
    stop(sprintf("LISA-KEGG-018 KGML %s has zero DE-mapped gene nodes; refuse to emit an unpainted map.", row$kegg_id), call. = FALSE)
  }
  node_values$analysis_id <- rep(cfg$analysis_id, nrow(node_values))
  node_values$universe <- rep(cfg$universe, nrow(node_values))
  node_values$kegg_title <- rep(row$kegg_title, nrow(node_values))
  node_values$node_has_value <- !is.na(node_values$log2FC)
  stem <- sprintf("%02d_%s_%s", row$rank, row$kegg_id, safe_title(row$kegg_title))
  out_png <- file.path(out_dir, paste0(stem, "_painted.png"))
  out_pdf <- file.path(out_dir, paste0(stem, "_painted.pdf"))
  base_png <- file.path(out_dir, paste0(stem, "_kegg_base.png"))
  draw_painted_map(
    img = img,
    node_values = node_values,
    title = sprintf("%s | %s | %s", cfg$analysis_id, row$kegg_id, row$kegg_title),
    subtitle = lisa_plot_subtitle(
      cfg$plot_metadata,
      suffix = sprintf(
        "Collection=%s; source=collection evidence genes mapped to KEGG",
        cfg$universe
      )
    ),
    out_png = out_png,
    out_pdf = out_pdf,
    max_abs = cfg$max_abs_log2fc,
    color_power = cfg$color_power
  )
  if (!file.copy(kegg_snapshot_resource(cfg, "image", row$kegg_id), base_png, overwrite = TRUE)) {
    stop(sprintf("LISA-KEGG-020 failed to preserve the KEGG base image for %s.", row$kegg_id), call. = FALSE)
  }
  contract_rows <- node_values
  contract_rows$figure_id <- paste0("kegg_single__", stem)
  contract_rows$figure_type <- "kegg_single"
  contract_rows$source_row_order <- seq_len(nrow(contract_rows))
  contract_rows$selected_for_plot <- TRUE
  contract_rows$highlighted <- contract_rows$node_has_value
  contract_rows$labelled <- FALSE
  contract_rows$base_image_file <- basename(base_png)
  contract_rows$max_abs_log2fc <- cfg$max_abs_log2fc
  contract_rows$color_power <- cfg$color_power
  contract_rows$figure_title <- sprintf("%s | %s | %s", cfg$analysis_id, row$kegg_id, row$kegg_title)
  contract_rows$figure_subtitle <- lisa_plot_subtitle(
    cfg$plot_metadata,
    suffix = sprintf("Collection=%s; source=collection evidence genes mapped to KEGG", cfg$universe)
  )
  source_path <- file.path(out_dir, paste0(stem, "_painted_source.tsv"))
  recipe_path <- file.path(out_dir, paste0(stem, "_painted_recipe.R"))
  lisaR:::lisa_write_figure_source_tsv(contract_rows, source_path)
  renderer <- file.path(dirname(normalizePath(script_file)), "reproduce_lisa_figure.R")
  if (!file.exists(renderer)) {
    stop("LISA-FIGURE-SOURCE-004 cannot install single-DE KEGG reproduction recipe.", call. = FALSE)
  }
  lisaR:::lisa_copy_verified_figure_recipe(
    renderer, recipe_path, cfg$lisa_internal_renderer_sha256,
    run_root = cfg$project_dir
  )
  node_values$output_png <- rep(out_png, nrow(node_values))
  node_values$output_pdf <- rep(out_pdf, nrow(node_values))
  attr(node_values, "output_png") <- out_png
  attr(node_values, "output_pdf") <- out_pdf
  node_values
}

write_readme <- function(path, cfg, painted_count) {
  lines <- c(
    "Single-DE canonical KEGG pathway painter",
    "",
    sprintf("Analysis: %s", cfg$analysis_id),
    sprintf("Universe: %s", cfg$universe),
    sprintf("Species map: %s", cfg$species),
    sprintf("KEGG access mode: %s", cfg$kegg_access_mode),
    sprintf("Retrieval date: %s", format(Sys.Date(), "%Y-%m-%d")),
    "",
    "Purpose",
    "This directory contains KEGG-derived pathway figures painted over KGML gene-node coordinates.",
    "Source pathway diagrams: KEGG, https://www.kegg.jp/.",
    "Pathways are selected by mapping this collection's post-LISA evidence genes to KEGG pathways.",
    "Node colors come from the full standardized DE table for this analysis.",
    "",
    "Important interpretation limits",
    "- This is a visual overlay, not a rerun of enrichment or causal evidence.",
    "- Grey nodes have no mapped DE value.",
    "",
    sprintf("Painted pathway maps generated: %d", painted_count),
    "",
    "Files",
    "- *_painted.png / *_painted.pdf: KEGG diagram with DE overlay.",
    "- *_kegg_pathway_painter_nodes.tsv: node-level mapping table for audit.",
    "- *_kegg_pathway_painter_index.tsv: pathway-level output index."
  )
  writeLines(lines, path, useBytes = TRUE)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  cfg$plot_metadata <- lisa_plot_metadata(cfg$project_dir, cfg$analysis_id)
  required_pkgs <- c("lisaR", "AnnotationDbi", "png", "grid", cfg$annotation_db)
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop(sprintf("Missing required R packages: %s", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
  paths <- paths_for(cfg$project_dir, cfg$analysis_id, cfg$universe)
  exact_mode <- !is.na(cfg$kegg_id) && nzchar(cfg$kegg_id)
  dir.create(paths$painter_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$painter_data_dir, recursive = TRUE, showWarnings = FALSE)
  # Rankings can change between runs. Remove only this painter's previously
  # generated maps so obsolete ranks cannot leak into a rebuilt report.
  #
  # In exact mode this sweep would be destructive rather than hygienic: the
  # request is for ONE map, and clearing the directory would delete every other
  # map a previous default run legitimately produced. Exact mode therefore
  # removes only the files belonging to the pathway it is about to redraw.
  stale_pattern <- if (exact_mode) {
    sprintf("^[0-9]+_%s_.*_painted\\.(png|pdf)$", cfg$kegg_id)
  } else "_painted\\.(png|pdf)$"
  stale_maps <- list.files(paths$painter_dir, pattern = stale_pattern, full.names = TRUE)
  if (length(stale_maps) && !all(unlink(stale_maps) == 0L)) {
    stop("LISA-KEGG-019 failed to remove stale painted maps before regeneration.", call. = FALSE)
  }
  stale_contract_pattern <- if (exact_mode) {
    sprintf("^[0-9]+_%s_.*_(painted_source[.]tsv|painted_recipe[.]R|kegg_base[.]png)$", cfg$kegg_id)
  } else "_(painted_source[.]tsv|painted_recipe[.]R|kegg_base[.]png)$"
  stale_contracts <- list.files(paths$painter_dir, pattern = stale_contract_pattern, full.names = TRUE)
  if (length(stale_contracts) && !all(unlink(stale_contracts) == 0L)) {
    stop("LISA-KEGG-019A failed to remove stale KEGG figure contracts.", call. = FALSE)
  }

  de <- load_de_table(paths$de)
  de_entrez <- map_symbols_to_entrez(de, cfg)
  evidence <- read_tsv(paths$evidence)
  built <- build_pathway_summary(evidence, de_entrez, cfg)
  summary <- built$summary

  full_summary_path <- file.path(paths$painter_data_dir, paste0(paths$prefix, "_kegg_pathway_painter_pathways.tsv"))
  gene_map_path <- file.path(paths$painter_data_dir, paste0(paths$prefix, "_kegg_pathway_painter_gene_map.tsv"))
  write_tsv(summary, full_summary_path)
  write_tsv(built$gene_map, gene_map_path)

  # The ranked summary is the same object in every mode, so an interface that
  # offers the selector sees exactly the ranking the painter will use.
  if (!is.na(cfg$emit_pathway_index) && nzchar(cfg$emit_pathway_index)) {
    dir.create(dirname(cfg$emit_pathway_index), recursive = TRUE, showWarnings = FALSE)
    write_tsv(summary, cfg$emit_pathway_index)
    message(sprintf("Wrote KEGG pathway index: %s (%d pathways)", cfg$emit_pathway_index, nrow(summary)))
  }
  if (isTRUE(cfg$index_only)) {
    message("Index-only run: no pathway map was painted.")
    return(invisible(NULL))
  }

  selected <- if (exact_mode) {
    # THE exact selector. It is applied to the fully ranked, collection-wide
    # summary -- not to head(summary, top_pathways) -- because a category's map
    # may legitimately rank below the default cut, and refusing it there would
    # make the selector a filter over the default product rather than a real
    # choice. The selected row keeps its own `rank`, so the stem, the title, the
    # node values and the colour context are exactly those the default generator
    # would have produced for this pathway.
    hit <- summary[as.character(summary$kegg_id) == cfg$kegg_id, , drop = FALSE]
    if (nrow(hit) == 0L) {
      stop(sprintf(
        "LISA-KEGG-029 pathway %s is not among the %d pathways this collection's evidence genes map to; it cannot be painted from %s/%s. Choose a pathway from the emitted pathway index.",
        cfg$kegg_id, nrow(summary), cfg$analysis_id, cfg$universe), call. = FALSE)
    }
    if (nrow(hit) != 1L) {
      stop(sprintf("LISA-KEGG-030 pathway %s resolved to %d summary rows; exactly one is required.", cfg$kegg_id, nrow(hit)), call. = FALSE)
    }
    hit
  } else if (cfg$top_pathways > 0 && nrow(summary) > cfg$top_pathways) {
    utils::head(summary, cfg$top_pathways)
  } else summary

  empty_index <- data.frame(
    rank = integer(),
    analysis_id = character(),
    universe = character(),
    kegg_id = character(),
    kegg_title = character(),
    output_png = character(),
    output_pdf = character(),
    kegg_access_mode = character(),
    retrieved_on = character(),
    stringsAsFactors = FALSE
  )
  # Exact mode writes a pathway-scoped inventory so it can never overwrite, or be
  # mistaken for, the index of a default full-collection run in the same
  # directory. The default names are untouched.
  output_stem <- if (exact_mode) cfg$kegg_id else paths$prefix
  index_path <- file.path(paths$painter_dir, paste0(output_stem, "_kegg_pathway_painter_index.tsv"))
  nodes_path <- file.path(paths$painter_dir, paste0(output_stem, "_kegg_pathway_painter_nodes.tsv"))
  readme_path <- file.path(paths$painter_dir,
    if (exact_mode) sprintf("README_kegg_painter_%s.txt", cfg$kegg_id) else "README_kegg_painter.txt")

  if (nrow(selected) == 0) {
    write_tsv(empty_index, index_path)
    write_tsv(data.frame(), nodes_path)
    write_readme(readme_path, cfg, 0)
    message("No KEGG pathways resolved from collection evidence genes; wrote empty painter outputs.")
    return(invisible(NULL))
  }

  all_nodes <- list()
  index <- data.frame()
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    nodes <- tryCatch(
      paint_one_pathway(row, de_entrez, cfg, paths$painter_dir),
      error = function(e) {
        data.frame(kegg_id = row$kegg_id, kegg_title = row$kegg_title, analysis_id = cfg$analysis_id, universe = cfg$universe, error = conditionMessage(e), stringsAsFactors = FALSE)
      }
    )
    all_nodes[[length(all_nodes) + 1]] <- nodes
    png_path <- if (!is.null(attr(nodes, "output_png"))) attr(nodes, "output_png") else if ("output_png" %in% names(nodes) && nrow(nodes) > 0) unique(nodes$output_png)[[1]] else ""
    pdf_path <- if (!is.null(attr(nodes, "output_pdf"))) attr(nodes, "output_pdf") else if ("output_pdf" %in% names(nodes) && nrow(nodes) > 0) unique(nodes$output_pdf)[[1]] else ""
    index <- rbind(
      index,
      data.frame(
        # The pathway's own collection-wide rank, not its position in this
        # loop. In a default run the two are identical because `selected` is the
        # head of the ranked summary; in exact mode only this one records the
        # true ranking context of the map.
        rank = as.integer(row$rank),
        analysis_id = cfg$analysis_id,
        universe = cfg$universe,
        kegg_id = row$kegg_id,
        kegg_title = row$kegg_title,
        n_evidence_genes = row$n_evidence_genes,
        output_png = png_path,
        output_pdf = pdf_path,
        kegg_access_mode = cfg$kegg_access_mode,
        retrieved_on = format(Sys.Date(), "%Y-%m-%d"),
        stringsAsFactors = FALSE
      )
    )
  }
  node_table <- do.call(rbind, all_nodes)
  write_tsv(node_table, nodes_path)
  write_tsv(index, index_path)
  write_readme(readme_path, cfg, sum(index$output_png != ""))
  lisaR:::lisa_assert_selected_kegg_painter_outputs(index, node_table)
  message(sprintf("Wrote single-DE KEGG pathway painter outputs: %s", paths$painter_dir))
  message(sprintf("Painted pathway maps: %d", sum(index$output_png != "")))
}

if (sys.nframe() == 0L) main()
