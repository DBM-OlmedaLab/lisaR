#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(system.file("scripts", "kegg_snapshot_helpers.R", package = "lisaR"), local = TRUE)

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

parse_args <- function(args) {
  out <- list(
    project_dir = NA_character_,
    contrast_id = NA_character_,
    universe = "GOBP-C2",
    species = NA_character_,
    kegg_cache_root = NA_character_,
    kegg_snapshot_id = NA_character_,
    kegg_access_mode = "external",
    top_pathways = 20,
    top_genes = 45,
    max_abs_log2fc = 0.5,
    color_power = 1.0,
    # --- exact pathway selection (H3) --------------------------------------
    # Analogous to the accepted H2 single-DE painter seam, and opt-in in the
    # same way: without these five this script behaves exactly as it always
    # has -- rank the whole contrast, take the top --top-pathways and paint
    # them.
    #
    #   --kegg-id            paint exactly ONE pathway, chosen from the SAME
    #                        complete collection-wide contrast ranking, keeping
    #                        its own `contrast_kegg_rank`, title, A/B node
    #                        halves, thresholds and colour context.
    #   --emit-pathway-index write that complete ranking, with its category
    #                        associations, to a named path so an interface can
    #                        offer the selector without a paint.
    #   --index-only         stop after the ranking. Enumerating what could be
    #                        painted must not paint anything.
    #   --paired-input       read an already prepared paired evidence table and
    #   --summary-input      its ranked category summary instead of the default
    #                        sibling paths.
    kegg_id = NA_character_,
    emit_pathway_index = NA_character_,
    index_only = "false",
    paired_input = NA_character_,
    summary_input = NA_character_,
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
    if (key %in% c("top_pathways", "top_genes")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("max_abs_log2fc", "color_power")) {
      out[[key]] <- as.numeric(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("LISA-KEGG-011 required argument --project-dir is missing. Repair: provide the validated run root for this contrast.", call. = FALSE)
  if (is.na(out$contrast_id) || out$contrast_id == "") stop("Required argument: --contrast-id", call. = FALSE)
  if (is.na(out$species) || out$species == "") stop("Required argument: --species", call. = FALSE)
  if (is.na(out$kegg_cache_root) || out$kegg_cache_root == "") stop("Required argument: --kegg-cache-root", call. = FALSE)
  if (is.na(out$kegg_snapshot_id) || out$kegg_snapshot_id == "") stop("Required argument: --kegg-snapshot-id", call. = FALSE)
  if (!out$kegg_access_mode %in% c("external", "cache_only")) stop("Invalid --kegg-access-mode; use external or cache_only.", call. = FALSE)
  if (!out$index_only %in% c("true", "false")) stop("--index-only must be true or false.", call. = FALSE)
  out$index_only <- identical(out$index_only, "true")
  if (!is.na(out$kegg_id) && out$kegg_id != "" && !grepl("^[a-z]{3}[0-9]{5}$", out$kegg_id)) {
    stop(sprintf("LISA-KEGG-047 --kegg-id must be one KEGG pathway identifier such as hsa04110; received: %s.", out$kegg_id), call. = FALSE)
  }
  if (!is.na(out$kegg_id) && out$kegg_id != "" && out$index_only) {
    stop("LISA-KEGG-048 --kegg-id and --index-only true are mutually exclusive: one paints a map, the other refuses to paint anything.", call. = FALSE)
  }
  supplied <- c(!is.na(out$paired_input) && nzchar(out$paired_input),
                !is.na(out$summary_input) && nzchar(out$summary_input))
  if (any(supplied) && !all(supplied)) {
    stop("LISA-KEGG-049 --paired-input and --summary-input must be given together; a paired table without its ranked summary is not a preparation.", call. = FALSE)
  }
  out <- c(out, kegg_species_contract(out$species))
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

first_nonempty <- function(x) {
  x <- safe_chr(x)
  x <- x[x != ""]
  if (length(x) == 0) "" else x[[1]]
}

contrast_row <- function(project_dir, contrast_id) {
  idx <- read_tsv(file.path(project_dir, "config", "contrast_index.tsv"))
  row <- idx[idx$contrast_id == contrast_id | paste(idx$contrast_id, idx$output_id, sep = "_") == contrast_id, , drop = FALSE]
  if (nrow(row) != 1) stop(sprintf("Contrast id did not resolve to one row: %s", contrast_id), call. = FALSE)
  row
}

paths_for <- function(project_dir, row, universe) {
  contrast_name <- paste(row$contrast_id, row$output_id, sep = "_")
  out_dir <- file.path(project_dir, "outputs", "gene_level", "category_contrasts", contrast_name, paste0("collection_", universe))
  prefix <- paste(row$output_id, universe, "contrast_gene_level", sep = "_")
  painter_dir <- file.path(out_dir, "contrast_kegg_pathway_painter")
  painter_data_dir <- file.path(painter_dir, "painter_data")
  list(
    contrast_name = contrast_name,
    out_dir = out_dir,
    de_a = file.path(project_dir, "outputs", "single_de", row$contrast_a, paste0("collection_", universe), "inputs", paste0(row$contrast_a, "_standardized_DE.tsv")),
    de_b = file.path(project_dir, "outputs", "single_de", row$contrast_b, paste0("collection_", universe), "inputs", paste0(row$contrast_b, "_standardized_DE.tsv")),
    paired = file.path(out_dir, paste0(prefix, "_paired_gene_evidence.tsv")),
    contrast_summary = file.path(out_dir, paste0(prefix, "_contrast_category_gene_summary.tsv")),
    kegg_summary = file.path(painter_data_dir, paste0(prefix, "_top_contrast_kegg_pathways.tsv")),
    kegg_full_summary = file.path(painter_data_dir, paste0(prefix, "_contrast_kegg_pathway_layer.tsv")),
    gene_matrix = file.path(painter_data_dir, paste0(prefix, "_contrast_kegg_gene_matrix.tsv")),
    painter_dir = painter_dir,
    painter_data_dir = painter_data_dir,
    prefix = prefix
  )
}

load_contrast_de_tables <- function(paths) {
  missing <- c(paths$de_a, paths$de_b)[!file.exists(c(paths$de_a, paths$de_b))]
  if (length(missing) > 0) {
    stop(sprintf("Missing standardized DE input(s) for contrast KEGG painter:\n%s", paste(missing, collapse = "\n")), call. = FALSE)
  }
  read_one <- function(path, side) {
    de <- read_tsv(path)
    symbol_col <- if ("symbol" %in% names(de)) "symbol" else if ("hgnc_symbol" %in% names(de)) "hgnc_symbol" else ""
    lfc_col <- if ("log2FoldChange" %in% names(de)) "log2FoldChange" else if ("log2FC" %in% names(de)) "log2FC" else ""
    padj_col <- if ("padj" %in% names(de)) "padj" else if ("FDR" %in% names(de)) "FDR" else ""
    gene_id_col <- if ("gene_id" %in% names(de)) "gene_id" else if ("original_gene_id" %in% names(de)) "original_gene_id" else ""
    if (symbol_col == "" || lfc_col == "") {
      stop(sprintf("Standardized DE table lacks symbol/log2FC columns: %s", path), call. = FALSE)
    }
    data.frame(
      side = side,
      symbol = safe_chr(de[[symbol_col]]),
      gene_id = if (gene_id_col != "") safe_chr(de[[gene_id_col]]) else "",
      log2FC = safe_num(de[[lfc_col]]),
      padj = if (padj_col != "") safe_num(de[[padj_col]]) else NA_real_,
      kegg_id = "",
      stringsAsFactors = FALSE
    )
  }
  out <- rbind(read_one(paths$de_a, "A"), read_one(paths$de_b, "B"))
  out <- out[out$symbol != "" & is.finite(out$log2FC), , drop = FALSE]
  out
}

load_kegg_data_helpers <- function() {
  helper_env <- new.env(parent = baseenv())
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
  script_path <- file.path(dirname(normalizePath(script_file)), "build_contrast_kegg_map_layer.R")
  if (!file.exists(script_path)) stop("Cannot locate build_contrast_kegg_map_layer.R helper script.", call. = FALSE)
  old <- Sys.getenv("KEGG_MAP_LAYER_LIBRARY_ONLY", unset = NA_character_)
  Sys.setenv(KEGG_MAP_LAYER_LIBRARY_ONLY = "1")
  sys.source(script_path, envir = helper_env)
  if (is.na(old)) {
    Sys.unsetenv("KEGG_MAP_LAYER_LIBRARY_ONLY")
  } else {
    Sys.setenv(KEGG_MAP_LAYER_LIBRARY_ONLY = old)
  }
  helper_env
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

build_painter_inputs <- function(paths, cfg, row) {
  # An explicitly supplied prepared pair takes precedence over the default
  # sibling paths, so a contrast map request is satisfied by the shared
  # preparation rather than by making the card/heatmap builder run first.
  paired_path <- if (!is.na(cfg$paired_input) && nzchar(cfg$paired_input)) cfg$paired_input else paths$paired
  summary_path <- if (!is.na(cfg$summary_input) && nzchar(cfg$summary_input)) cfg$summary_input else paths$contrast_summary
  missing <- c(paired_path, summary_path)[!file.exists(c(paired_path, summary_path))]
  if (length(missing) > 0) stop(sprintf("Missing required input files:\n%s", paste(missing, collapse = "\n")), call. = FALSE)
  dir.create(paths$painter_data_dir, recursive = TRUE, showWarnings = FALSE)
  helper <- load_kegg_data_helpers()
  paired <- read_tsv(paired_path)
  kegg_rows <- helper$build_kegg_rows(paired)
  if (nrow(kegg_rows) == 0) {
    fallback <- build_gene_mapped_contrast_summary(paired, cfg, row)
    write_tsv(fallback$summary, paths$kegg_full_summary)
    write_tsv(if (cfg$top_pathways > 0) utils::head(fallback$summary, cfg$top_pathways) else fallback$summary, paths$kegg_summary)
    write_tsv(fallback$gene_matrix, paths$gene_matrix)
    return(invisible(nrow(fallback$summary) > 0))
  }
  pathway_summary <- helper$summarize_pathways(kegg_rows, cfg, row)
  gene_matrix <- helper$build_gene_matrix(pathway_summary, kegg_rows, cfg)
  write_tsv(pathway_summary, paths$kegg_full_summary)
  write_tsv(if (cfg$top_pathways > 0) utils::head(pathway_summary, cfg$top_pathways) else pathway_summary, paths$kegg_summary)
  write_tsv(gene_matrix, paths$gene_matrix)
  invisible(TRUE)
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

build_gene_mapped_contrast_summary <- function(paired, cfg, row) {
  empty <- list(summary = data.frame(), gene_matrix = data.frame())
  if (nrow(paired) == 0 || !"symbol" %in% names(paired)) return(empty)
  side_rows <- lapply(c("A", "B"), function(side) {
    present_col <- paste0("present_", side)
    lfc_col <- paste0("log2FC_", side)
    padj_col <- paste0("padj_", side)
    gene_id_col <- paste0("gene_id_", side)
    score_col <- paste0("gene_contribution_score_", side)
    if (!all(c(present_col, lfc_col) %in% names(paired))) return(data.frame())
    keep <- paired[[present_col]] %in% TRUE | safe_chr(paired[[present_col]]) == "TRUE"
    df <- paired[keep, , drop = FALSE]
    if (nrow(df) == 0) return(data.frame())
    data.frame(
      side = side,
      symbol = safe_chr(df$symbol),
      gene_id = if (gene_id_col %in% names(df)) safe_chr(df[[gene_id_col]]) else "",
      log2FC = safe_num(df[[lfc_col]]),
      padj = if (padj_col %in% names(df)) safe_num(df[[padj_col]]) else NA_real_,
      category_id = safe_chr(df$category_id),
      category_display_name = safe_chr(df$category_display_name),
      contrast_gene_class = if ("contrast_gene_class" %in% names(df)) safe_chr(df$contrast_gene_class) else "",
      gene_contribution_score = if (score_col %in% names(df)) safe_num(df[[score_col]]) else NA_real_,
      kegg_id = "",
      stringsAsFactors = FALSE
    )
  })
  evidence <- do.call(rbind, side_rows)
  evidence <- evidence[evidence$symbol != "" & is.finite(evidence$log2FC), , drop = FALSE]
  if (nrow(evidence) == 0) return(empty)
  entrez <- map_symbols_to_entrez(evidence, cfg)
  if (nrow(entrez) == 0) return(empty)
  if ("kegg_id" %in% names(entrez)) entrez$kegg_id <- NULL
  g2p <- gene_to_pathway_map(cfg)
  evp <- merge(entrez, g2p, by.x = "ENTREZID", by.y = "kegg_gene", all.x = FALSE, all.y = FALSE)
  if (nrow(evp) == 0) return(empty)
  titles <- pathway_titles(cfg)
  evp <- merge(evp, titles, by = "kegg_id", all.x = TRUE)

  split_rows <- split(evp, evp$kegg_id)
  summary <- lapply(split_rows, function(df) {
    class_vals <- safe_chr(df$contrast_gene_class)
    data.frame(
      contrast_id = row$contrast_id,
      output_id = row$output_id,
      universe = cfg$universe,
      pathway = paste0("KEGG_", gsub("[^A-Za-z0-9]+", "_", toupper(df$kegg_title[[1]]))),
      pathway_display_label = df$kegg_title[[1]],
      kegg_id = df$kegg_id[[1]],
      kegg_title = df$kegg_title[[1]],
      contrast_kegg_rank = NA_integer_,
      n_shared_opposite_genes = sum(class_vals == "shared_opposite_direction", na.rm = TRUE),
      pathway_gene_total = length(unique(df$SYMBOL)),
      max_abs_log2FC = max(abs(safe_num(df$log2FC)), na.rm = TRUE),
      max_gene_contribution_score = max(safe_num(df$gene_contribution_score), na.rm = TRUE),
      genes = paste(sort(unique(safe_chr(df$SYMBOL))), collapse = "/"),
      categories = paste(sort(unique(safe_chr(df$category_id))), collapse = ";"),
      evidence_basis = "collection evidence genes mapped to KEGG",
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary)
  summary$max_abs_log2FC[is.infinite(summary$max_abs_log2FC)] <- NA_real_
  summary$max_gene_contribution_score[is.infinite(summary$max_gene_contribution_score)] <- NA_real_
  summary <- summary[order(-summary$n_shared_opposite_genes, -summary$pathway_gene_total, -safe_num(summary$max_gene_contribution_score), summary$kegg_title), , drop = FALSE]
  summary$contrast_kegg_rank <- seq_len(nrow(summary))
  gene_matrix <- evp[, c("side", "SYMBOL", "ENTREZID", "log2FC", "padj", "category_id", "category_display_name", "contrast_gene_class", "kegg_id", "kegg_title"), drop = FALSE]
  gene_matrix <- gene_matrix[order(gene_matrix$kegg_id, gene_matrix$side, gene_matrix$SYMBOL), , drop = FALSE]
  list(summary = summary, gene_matrix = gene_matrix)
}

norm_pathway_name <- function(x) {
  x <- toupper(safe_chr(x))
  x <- gsub("^KEGG_MEDICUS_(REFERENCE|VARIANT)_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("[^A-Z0-9]+", " ", x)
  gsub("\\s+", " ", trimws(x))
}

attr_value <- function(txt, attr) {
  m <- regexec(paste0(attr, "=\"([^\"]*)\""), txt, perl = TRUE)
  out <- regmatches(txt, m)
  vapply(out, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
}

extract_entry_blocks <- function(kgml) {
  # readLines() returns one element per XML line.  Collapse it before the
  # multiline match or every pathway is rendered with zero parsed nodes.
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

  entry_id <- attr_value(blocks, "id")
  entry_name <- attr_value(blocks, "name")
  graphics <- regmatches(blocks, regexpr("<graphics\\b[^>]*/?>", blocks, perl = TRUE))
  graphics[lengths(graphics) == 0] <- NA_character_

  x <- safe_num(attr_value(graphics, "x"))
  y <- safe_num(attr_value(graphics, "y"))
  width <- safe_num(attr_value(graphics, "width"))
  height <- safe_num(attr_value(graphics, "height"))
  label <- attr_value(graphics, "name")

  ids <- lapply(strsplit(entry_name, "\\s+"), function(v) {
    gsub("^[a-z]{3}:", "", v)
  })
  data.frame(
    entry_id = entry_id,
    kegg_entry_name = entry_name,
    kegg_gene_ids = vapply(ids, paste, collapse = "/", FUN.VALUE = character(1)),
    x = x,
    y = y,
    width = width,
    height = height,
    kegg_label = label,
    stringsAsFactors = FALSE
  )
}

resolve_kegg_pathways <- function(pathway_table, cfg) {
  klist <- kegg_snapshot_rds(cfg, "pathway_list", "all")
  kegg_ids <- sub("^path:", "", names(klist))
  species_suffix <- if (cfg$kegg_code == "mmu") " - Mus musculus \\(house mouse\\)$" else " - Homo sapiens \\(human\\)$"
  titles <- sub(species_suffix, "", unname(klist))
  ref <- data.frame(
    kegg_id = kegg_ids,
    kegg_title = titles,
    normalized_title = norm_pathway_name(titles),
    stringsAsFactors = FALSE
  )
  pathway_table$normalized_pathway <- norm_pathway_name(pathway_table$pathway)
  out <- merge(pathway_table, ref, by.x = "normalized_pathway", by.y = "normalized_title", all.x = TRUE)
  out <- out[order(safe_num(out$contrast_kegg_rank)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

map_symbols_to_entrez <- function(de_table, cfg) {
  if (!"symbol" %in% names(de_table) && "matrix.symbol" %in% names(de_table)) de_table$symbol <- de_table[["matrix.symbol"]]
  if (!"gene_id" %in% names(de_table) && "matrix.gene_id" %in% names(de_table)) de_table$gene_id <- de_table[["matrix.gene_id"]]
  if (!"side" %in% names(de_table) && "matrix.side" %in% names(de_table)) de_table$side <- de_table[["matrix.side"]]
  if (!"log2FC" %in% names(de_table) && "matrix.log2FC" %in% names(de_table)) de_table$log2FC <- de_table[["matrix.log2FC"]]
  if (!"padj" %in% names(de_table) && "matrix.padj" %in% names(de_table)) de_table$padj <- de_table[["matrix.padj"]]
  if (!"kegg_id" %in% names(de_table) && "matrix.pathway" %in% names(de_table)) de_table$kegg_id <- de_table[["matrix.pathway"]]
  if (!"gene_id" %in% names(de_table)) de_table$gene_id <- ""
  if (!"kegg_id" %in% names(de_table)) de_table$kegg_id <- ""
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
    merged <- merged[!duplicated(paste(merged$side, merged$ENTREZID, merged$input_symbol, merged$log2FC, merged$padj, sep = "\r")), , drop = FALSE]
  }
  merged$ENTREZID <- safe_chr(merged$ENTREZID)
  merged$SYMBOL <- safe_chr(merged$input_symbol)
  merged
}

best_hit_for_ids <- function(de_entrez, ids, side, kegg_id) {
  hits <- de_entrez[de_entrez$ENTREZID %in% ids & de_entrez$side == side, , drop = FALSE]
  if ("kegg_id" %in% names(hits)) {
    hits <- hits[hits$kegg_id == kegg_id | is.na(hits$kegg_id) | hits$kegg_id == "", , drop = FALSE]
  }
  if (nrow(hits) == 0) return(list(symbol = "", log2FC = NA_real_, padj = NA_real_))
  hits$abs_lfc <- abs(safe_num(hits$log2FC))
  hits <- hits[order(-hits$abs_lfc, safe_num(hits$padj), hits$SYMBOL), , drop = FALSE]
  list(
    symbol = paste(unique(safe_chr(hits$SYMBOL)), collapse = "/"),
    log2FC = safe_num(hits$log2FC)[[1]],
    padj = safe_num(hits$padj)[[1]]
  )
}

node_values_for_pathway <- function(nodes, de_entrez, kegg_id) {
  if (nrow(nodes) == 0) return(data.frame())
  out <- lapply(seq_len(nrow(nodes)), function(i) {
    ids <- unlist(strsplit(nodes$kegg_gene_ids[[i]], "/", fixed = TRUE), use.names = FALSE)
    hit_a <- best_hit_for_ids(de_entrez, ids, "A", kegg_id)
    hit_b <- best_hit_for_ids(de_entrez, ids, "B", kegg_id)
    data.frame(
      nodes[i, , drop = FALSE],
      symbol_A = hit_a$symbol,
      log2FC_A = hit_a$log2FC,
      padj_A = hit_a$padj,
      symbol_B = hit_b$symbol,
      log2FC_B = hit_b$log2FC,
      padj_B = hit_b$padj,
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
  if (z >= 0) {
    white * (1 - z) + red * z
  } else {
    z <- abs(z)
    white * (1 - z) + blue * z
  }
}

paint_half_region <- function(out, xmin, xmax, ymin, ymax, value, max_abs, color_power, alpha_value = 0.95) {
  if (xmin > xmax || ymin > ymax) return(out)
  fill_rgb <- lfc_to_rgb(value, max_abs, color_power)
  alpha <- if (is.na(value)) 0.30 else alpha_value
  region <- out[ymin:ymax, xmin:xmax, 1:3, drop = FALSE]
  luminance <- 0.299 * region[, , 1] + 0.587 * region[, , 2] + 0.114 * region[, , 3]
  preserve <- luminance < 0.42
  for (ch in seq_len(3)) {
    channel <- region[, , ch]
    channel[!preserve] <- channel[!preserve] * (1 - alpha) + fill_rgb[[ch]] * alpha
    region[, , ch] <- channel
  }
  out[ymin:ymax, xmin:xmax, 1:3] <- region
  out
}

paint_node_backgrounds <- function(img, node_values, max_abs, color_power) {
  out <- img
  h <- dim(out)[1]
  w <- dim(out)[2]
  channels <- dim(out)[3]
  if (channels < 3) stop("Expected an RGB/RGBA KEGG image.", call. = FALSE)

  for (i in seq_len(nrow(node_values))) {
    if (is.na(node_values$x[[i]]) || is.na(node_values$y[[i]])) next
    xmin <- max(1, floor(node_values$x[[i]] - node_values$width[[i]] / 2))
    xmax <- min(w, ceiling(node_values$x[[i]] + node_values$width[[i]] / 2))
    ymin <- max(1, floor(node_values$y[[i]] - node_values$height[[i]] / 2))
    ymax <- min(h, ceiling(node_values$y[[i]] + node_values$height[[i]] / 2))
    if (xmin > xmax || ymin > ymax) next
    xmid <- floor((xmin + xmax) / 2)
    out <- paint_half_region(out, xmin, xmid, ymin, ymax, node_values$log2FC_A[[i]], max_abs, color_power)
    out <- paint_half_region(out, xmid + 1, xmax, ymin, ymax, node_values$log2FC_B[[i]], max_abs, color_power)
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
    grid::grid.rect(
      x = w / 2,
      y = total_h / 2,
      width = w,
      height = total_h,
      default.units = "native",
      gp = grid::gpar(fill = "white", col = NA)
    )
    grid::grid.raster(
      painted_img,
      x = w / 2,
      y = bottom_pad + h / 2,
      width = w,
      height = h,
      default.units = "native",
      interpolate = FALSE
    )
    for (i in seq_len(nrow(node_values))) {
      if (is.na(node_values$x[[i]])) next
      grid::grid.rect(
        x = node_values$x[[i]],
        y = bottom_pad + h - node_values$y[[i]],
        width = node_values$width[[i]],
        height = node_values$height[[i]],
        default.units = "native",
        gp = grid::gpar(fill = NA, col = "#202020", lwd = 0.45)
      )
      grid::grid.segments(
        x0 = node_values$x[[i]],
        x1 = node_values$x[[i]],
        y0 = bottom_pad + h - (node_values$y[[i]] - node_values$height[[i]] / 2),
        y1 = bottom_pad + h - (node_values$y[[i]] + node_values$height[[i]] / 2),
        default.units = "native",
        gp = grid::gpar(col = "#202020", lwd = 0.25)
      )
    }
    grid::grid.text(
      title,
      x = 10,
      y = total_h - 12,
      just = c("left", "top"),
      default.units = "native",
      gp = grid::gpar(fontsize = 10, fontface = "bold", col = "#111111")
    )
    grid::grid.text(
      subtitle,
      x = 10,
      y = total_h - 28,
      just = c("left", "top"),
      default.units = "native",
      gp = grid::gpar(fontsize = 8, col = "#111111")
    )
    grid::grid.text(
      sprintf("Split node fill: left=A, right=B; display log2FC capped at +/- %.2g, color power %.2g. Red=up, blue=down, grey=no mapped value. Source pathway diagram: KEGG (www.kegg.jp).", max_abs, color_power),
      x = 10,
      y = 12,
      just = c("left", "bottom"),
      default.units = "native",
      gp = grid::gpar(fontsize = 8, col = "#111111")
    )
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

paint_one_pathway <- function(kegg_id, kegg_title, de_entrez, cfg, row, out_dir, prefix, rank) {
  kgml <- readLines(kegg_snapshot_resource(cfg, "kgml", kegg_id), warn = FALSE)
  img <- png::readPNG(kegg_snapshot_resource(cfg, "image", kegg_id))
  nodes <- parse_kegg_gene_nodes(kgml)
  if (!nrow(nodes)) {
    stop(sprintf("LISA-KEGG-017 KGML %s contains zero parsed gene nodes; refuse to emit an unpainted contrast map.", kegg_id), call. = FALSE)
  }
  node_values <- node_values_for_pathway(nodes, de_entrez, kegg_id)
  if (!any(!is.na(node_values$log2FC_A)) && !any(!is.na(node_values$log2FC_B))) {
    stop(sprintf("LISA-KEGG-018 KGML %s has zero DE-mapped gene nodes; refuse to emit an unpainted contrast map.", kegg_id), call. = FALSE)
  }
  node_values$kegg_id <- rep(kegg_id, nrow(node_values))
  node_values$kegg_title <- rep(kegg_title, nrow(node_values))
  node_values$contrast_id <- rep(row$contrast_id, nrow(node_values))
  node_values$output_id <- rep(row$output_id, nrow(node_values))
  node_values$universe <- rep(cfg$universe, nrow(node_values))
  node_values$node_has_A_value <- !is.na(node_values$log2FC_A)
  node_values$node_has_B_value <- !is.na(node_values$log2FC_B)

  safe_title <- gsub("[^A-Za-z0-9]+", "_", kegg_title)
  safe_title <- gsub("_+", "_", safe_title)
  safe_title <- gsub("^_|_$", "", safe_title)
  stem <- sprintf("%02d_%s_%s", rank, kegg_id, safe_title)
  out_png <- file.path(out_dir, paste0(stem, "_contrast_painted.png"))
  out_pdf <- file.path(out_dir, paste0(stem, "_contrast_painted.pdf"))
  base_png <- file.path(out_dir, paste0(stem, "_kegg_base.png"))
  draw_painted_map(
    img = img,
    node_values = node_values,
    title = sprintf("%s | %s | %s", row$output_id, kegg_id, kegg_title),
    subtitle = sprintf("A=%s; B=%s", row$label_a, row$label_b),
    out_png = out_png,
    out_pdf = out_pdf,
    max_abs = cfg$max_abs_log2fc,
    color_power = cfg$color_power
  )
  if (!file.copy(kegg_snapshot_resource(cfg, "image", kegg_id), base_png, overwrite = TRUE)) {
    stop(sprintf("LISA-KEGG-020 failed to preserve the KEGG base image for %s.", kegg_id), call. = FALSE)
  }
  contract_rows <- node_values
  contract_rows$figure_id <- paste0("kegg_contrast__", stem)
  contract_rows$figure_type <- "kegg_contrast"
  contract_rows$source_row_order <- seq_len(nrow(contract_rows))
  contract_rows$selected_for_plot <- TRUE
  contract_rows$highlighted <- contract_rows$node_has_A_value | contract_rows$node_has_B_value
  contract_rows$labelled <- FALSE
  contract_rows$base_image_file <- basename(base_png)
  contract_rows$max_abs_log2fc <- cfg$max_abs_log2fc
  contract_rows$color_power <- cfg$color_power
  contract_rows$figure_title <- sprintf("%s | %s | %s", row$output_id, kegg_id, kegg_title)
  contract_rows$figure_subtitle <- sprintf("A=%s; B=%s", row$label_a, row$label_b)
  source_path <- file.path(out_dir, paste0(stem, "_contrast_painted_source.tsv"))
  recipe_path <- file.path(out_dir, paste0(stem, "_contrast_painted_recipe.R"))
  lisaR:::lisa_write_figure_source_tsv(contract_rows, source_path)
  renderer <- file.path(dirname(normalizePath(script_file)), "reproduce_lisa_figure.R")
  if (!file.exists(renderer)) {
    stop("LISA-FIGURE-SOURCE-005 cannot install contrast KEGG reproduction recipe.", call. = FALSE)
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

write_readme <- function(path, cfg, row, painted_count) {
  lines <- c(
    "Contrast canonical KEGG pathway painter",
    "",
    sprintf("Contrast: %s_%s", row$contrast_id, row$output_id),
    sprintf("Universe: %s", cfg$universe),
    sprintf("Species map: %s", cfg$species),
    sprintf("KEGG access mode: %s", cfg$kegg_access_mode),
    sprintf("Retrieval date: %s", format(Sys.Date(), "%Y-%m-%d")),
    sprintf("A: %s (%s)", row$contrast_a, row$label_a),
    sprintf("B: %s (%s)", row$contrast_b, row$label_b),
    "",
    "Purpose",
    "This directory contains KEGG-derived pathway figures painted over KGML gene-node coordinates.",
    "Source pathway diagrams: KEGG, https://www.kegg.jp/.",
    "Each gene node is split in two: left half is A log2FC, right half is B log2FC.",
    "Node colors are mapped from the full standardized DE table for each side, not only from the KEGG contrast evidence subset.",
    sprintf("Display color is capped at +/- %.2g log2FC with color power %.2g. This matches the single-DE painter scale by default.", cfg$max_abs_log2fc, cfg$color_power),
    "",
    "Important interpretation limits",
    "- This is a visual overlay, not a rerun of enrichment or causal evidence.",
    sprintf("- Symbol-to-Entrez mapping uses %s inherited from the run species contract.", cfg$annotation_db),
    "- Grey halves are KEGG node halves without a mapped DE value from that side.",
    "",
    sprintf("Painted pathway maps generated: %d", painted_count),
    "",
    "Files",
    "- *_contrast_painted.png / *_contrast_painted.pdf: KEGG diagram with split A/B DE overlay.",
    "- *_contrast_kegg_pathway_painter_nodes.tsv: node-level mapping table for audit.",
    "- *_contrast_kegg_pathway_painter_index.tsv: pathway-level output index."
  )
  writeLines(lines, path, useBytes = TRUE)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  required_pkgs <- c("lisaR", "AnnotationDbi", "png", "grid", cfg$annotation_db)
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop(sprintf("Missing required R packages: %s", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
  row <- contrast_row(cfg$project_dir, cfg$contrast_id)
  paths <- paths_for(cfg$project_dir, row, cfg$universe)
  exact_mode <- !is.na(cfg$kegg_id) && nzchar(cfg$kegg_id)
  index_requested <- !is.na(cfg$emit_pathway_index) && nzchar(cfg$emit_pathway_index)
  dir.create(paths$painter_dir, recursive = TRUE, showWarnings = FALSE)
  # Rankings can change between runs. Remove only this painter's previously
  # generated maps so obsolete ranks cannot leak into a rebuilt report.
  #
  # In exact mode that sweep would be destructive rather than hygienic: the
  # sibling maps it would delete are valid products of other requests. It is
  # therefore narrowed to the one pathway about to be repainted.
  stale_pattern <- if (exact_mode) {
    sprintf("^[0-9]+_%s_.*_contrast_painted\\.(png|pdf)$", cfg$kegg_id)
  } else "_contrast_painted\\.(png|pdf)$"
  stale_maps <- list.files(paths$painter_dir, pattern = stale_pattern, full.names = TRUE)
  if (length(stale_maps) && !all(unlink(stale_maps) == 0L)) {
    stop("LISA-KEGG-019 failed to remove stale contrast painted maps before regeneration.", call. = FALSE)
  }
  stale_contract_pattern <- if (exact_mode) {
    sprintf("^[0-9]+_%s_.*_(contrast_painted_source[.]tsv|contrast_painted_recipe[.]R|kegg_base[.]png)$", cfg$kegg_id)
  } else "_(contrast_painted_source[.]tsv|contrast_painted_recipe[.]R|kegg_base[.]png)$"
  stale_contracts <- list.files(paths$painter_dir, pattern = stale_contract_pattern, full.names = TRUE)
  if (length(stale_contracts) && !all(unlink(stale_contracts) == 0L)) {
    stop("LISA-KEGG-019A failed to remove stale contrast KEGG figure contracts.", call. = FALSE)
  }
  has_kegg_inputs <- build_painter_inputs(paths, cfg, row)
  if (!isTRUE(has_kegg_inputs)) {
    if (exact_mode) {
      stop(sprintf("LISA-KEGG-052 this contrast/collection has no KEGG source pathway at all, so pathway %s cannot be painted.", cfg$kegg_id), call. = FALSE)
    }
    if (index_requested) {
      # Enumeration with nothing to enumerate writes an empty index and touches
      # no painter output. It must not overwrite a default painter run's files.
      dir.create(dirname(cfg$emit_pathway_index), recursive = TRUE, showWarnings = FALSE)
      write_tsv(data.frame(rank = integer(), contrast_id = character(),
        output_id = character(), universe = character(), kegg_id = character(),
        kegg_title = character(), contrast_kegg_rank = numeric(),
        n_shared_opposite_genes = numeric(), pathway_gene_total = numeric(),
        categories = character(), stringsAsFactors = FALSE),
        cfg$emit_pathway_index)
      message("No KEGG source pathways found; wrote an empty contrast pathway index.")
      if (isTRUE(cfg$index_only)) return(invisible(NULL))
    }
    empty_index <- data.frame(
      rank = integer(),
      contrast_id = character(),
      output_id = character(),
      universe = character(),
      source_pathway = character(),
      kegg_id = character(),
      kegg_title = character(),
      output_png = character(),
      output_pdf = character(),
      kegg_access_mode = character(),
      retrieved_on = character(),
      stringsAsFactors = FALSE
    )
    write_tsv(empty_index, file.path(paths$painter_dir, paste0(paths$prefix, "_contrast_kegg_pathway_painter_index.tsv")))
    write_tsv(data.frame(), file.path(paths$painter_dir, paste0(paths$prefix, "_contrast_kegg_pathway_painter_nodes.tsv")))
    write_readme(file.path(paths$painter_dir, "README_contrast_kegg_pathway_painter.txt"), cfg, row, 0)
    message("No KEGG source pathways found; wrote empty pathway painter outputs.")
    return(invisible(NULL))
  }

  # THE complete collection-wide contrast ranking. The default path keeps
  # reading the `top_pathways` head it always read; exact selection and index
  # emission read the FULL summary, because selecting a pathway from a truncated
  # head would silently make a legitimately ranked map unreachable.
  summary <- read_tsv(if (exact_mode || index_requested) paths$kegg_full_summary
                      else paths$kegg_summary)
  if ("kegg_id" %in% names(summary) && any(nzchar(safe_chr(summary$kegg_id)))) {
    resolved <- summary[nzchar(safe_chr(summary$kegg_id)), , drop = FALSE]
    if (!"kegg_title" %in% names(resolved) || any(!nzchar(safe_chr(resolved$kegg_title)))) {
      titles <- pathway_titles(cfg)
      resolved <- merge(resolved, titles, by = "kegg_id", all.x = TRUE)
    }
  } else {
    summary <- summary[grepl("^KEGG_", safe_chr(summary$pathway)) & !grepl("^KEGG_MEDICUS", safe_chr(summary$pathway)), , drop = FALSE]
    if (nrow(summary) == 0) stop("No canonical KEGG rows in contrast KEGG summary.", call. = FALSE)
    resolved <- resolve_kegg_pathways(summary, cfg)
    resolved <- resolved[!is.na(resolved$kegg_id) & resolved$kegg_id != "", , drop = FALSE]
  }
  if (!"n_shared_opposite_genes" %in% names(resolved)) resolved$n_shared_opposite_genes <- 0
  if (!"pathway_gene_total" %in% names(resolved)) resolved$pathway_gene_total <- 0
  if (!"contrast_kegg_rank" %in% names(resolved)) resolved$contrast_kegg_rank <- seq_len(nrow(resolved))
  resolved <- resolved[order(safe_num(resolved$contrast_kegg_rank), -safe_num(resolved$n_shared_opposite_genes), -safe_num(resolved$pathway_gene_total)), , drop = FALSE]
  # The category associations of each pathway, collected BEFORE the per-pathway
  # deduplication below discards the extra (pathway, category) rows. They are a
  # navigation aid only: a contrast map is one contrast/collection-wide figure
  # and is never identified, cropped or duplicated by category.
  pathway_categories <- if ("category_id" %in% names(resolved)) {
    vapply(split(safe_chr(resolved$category_id), safe_chr(resolved$kegg_id)),
      function(ids) paste(sort(unique(ids[nzchar(ids)])), collapse = ";"), character(1))
  } else character()
  resolved <- resolved[!duplicated(resolved$kegg_id), , drop = FALSE]

  # Enumerating what could be painted must not paint anything. The index is the
  # ranking a paint WOULD use, written by the painter itself, so the selector an
  # interface offers cannot drift from the painter's own ordering.
  if (index_requested) {
    emitted <- data.frame(
      rank = seq_len(nrow(resolved)),
      contrast_id = rep(row$contrast_id, nrow(resolved)),
      output_id = rep(row$output_id, nrow(resolved)),
      universe = rep(cfg$universe, nrow(resolved)),
      kegg_id = safe_chr(resolved$kegg_id),
      kegg_title = safe_chr(resolved$kegg_title),
      contrast_kegg_rank = safe_num(resolved$contrast_kegg_rank),
      n_shared_opposite_genes = safe_num(resolved$n_shared_opposite_genes),
      pathway_gene_total = safe_num(resolved$pathway_gene_total),
      categories = unname(pathway_categories[safe_chr(resolved$kegg_id)]),
      stringsAsFactors = FALSE)
    emitted$categories[is.na(emitted$categories)] <- ""
    dir.create(dirname(cfg$emit_pathway_index), recursive = TRUE, showWarnings = FALSE)
    write_tsv(emitted, cfg$emit_pathway_index)
    message(sprintf("Wrote contrast pathway index with %d pathway(s): %s",
      nrow(emitted), cfg$emit_pathway_index))
  }
  if (isTRUE(cfg$index_only)) return(invisible(NULL))

  selected <- if (exact_mode) {
    # THE exact selector. Applied to the fully ranked, collection-wide contrast
    # summary -- not to head(resolved, top_pathways) -- so the map's
    # `contrast_kegg_rank`, title, A/B endpoint node values, thresholds and
    # colour context are exactly the ones the default generator would produce
    # for that pathway.
    hit <- resolved[safe_chr(resolved$kegg_id) == cfg$kegg_id, , drop = FALSE]
    if (nrow(hit) == 0) {
      stop(sprintf("LISA-KEGG-050 pathway %s is not among the %d pathways this contrast's evidence genes map to in %s / %s.",
        cfg$kegg_id, nrow(resolved), cfg$contrast_id, cfg$universe), call. = FALSE)
    }
    if (nrow(hit) != 1) {
      stop(sprintf("LISA-KEGG-051 pathway %s resolved to %d summary rows; exactly one is required.", cfg$kegg_id, nrow(hit)), call. = FALSE)
    }
    hit
  } else if (cfg$top_pathways > 0) utils::head(resolved, cfg$top_pathways) else resolved
  # The rank carried into the file stem is the pathway's own position in the
  # complete ranking, so an exact map is named exactly as the default loop would
  # have named it.
  selected_ranks <- if (exact_mode) {
    match(safe_chr(selected$kegg_id), safe_chr(resolved$kegg_id))
  } else seq_len(nrow(selected))

  if (nrow(selected) == 0) {
    empty_index <- data.frame(
      rank = integer(),
      contrast_id = character(),
      output_id = character(),
      universe = character(),
      source_pathway = character(),
      kegg_id = character(),
      kegg_title = character(),
      output_png = character(),
      output_pdf = character(),
      kegg_access_mode = character(),
      retrieved_on = character(),
      stringsAsFactors = FALSE
    )
    empty_nodes <- data.frame()
    write_tsv(empty_index, file.path(paths$painter_dir, paste0(paths$prefix, "_contrast_kegg_pathway_painter_index.tsv")))
    write_tsv(empty_nodes, file.path(paths$painter_dir, paste0(paths$prefix, "_contrast_kegg_pathway_painter_nodes.tsv")))
    write_readme(file.path(paths$painter_dir, "README_contrast_kegg_pathway_painter.txt"), cfg, row, 0)
    message("No canonical KEGG pathways resolved; wrote empty pathway painter outputs.")
    return(invisible(NULL))
  }

  contrast_de <- load_contrast_de_tables(paths)
  de_entrez <- map_symbols_to_entrez(contrast_de, cfg)

  all_nodes <- list()
  index <- data.frame()
  for (i in seq_len(nrow(selected))) {
    selected_row <- selected[i, , drop = FALSE]
    rank_i <- selected_ranks[[i]]
    nodes <- tryCatch(
      paint_one_pathway(selected_row$kegg_id, selected_row$kegg_title, de_entrez, cfg, row, paths$painter_dir, paths$prefix, rank_i),
      error = function(e) {
        data.frame(
          kegg_id = selected_row$kegg_id,
          kegg_title = selected_row$kegg_title,
          contrast_id = row$contrast_id,
          output_id = row$output_id,
          universe = cfg$universe,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
    all_nodes[[length(all_nodes) + 1]] <- nodes
    png_path <- if (!is.null(attr(nodes, "output_png"))) attr(nodes, "output_png") else if ("output_png" %in% names(nodes) && nrow(nodes) > 0) unique(nodes$output_png)[[1]] else ""
    pdf_path <- if (!is.null(attr(nodes, "output_pdf"))) attr(nodes, "output_pdf") else if ("output_pdf" %in% names(nodes) && nrow(nodes) > 0) unique(nodes$output_pdf)[[1]] else ""
    index <- rbind(
      index,
      data.frame(
        rank = rank_i,
        contrast_id = row$contrast_id,
        output_id = row$output_id,
        universe = cfg$universe,
        source_pathway = selected_row$pathway,
        kegg_id = selected_row$kegg_id,
        kegg_title = selected_row$kegg_title,
        contrast_kegg_rank = selected_row$contrast_kegg_rank,
        n_shared_opposite_genes = selected_row$n_shared_opposite_genes,
        pathway_gene_total = selected_row$pathway_gene_total,
        output_png = png_path,
        output_pdf = pdf_path,
        kegg_access_mode = cfg$kegg_access_mode,
        retrieved_on = format(Sys.Date(), "%Y-%m-%d"),
        stringsAsFactors = FALSE
      )
    )
  }

  node_table <- do.call(rbind, all_nodes)
  # Exact mode names its inventory, node audit and README after the pathway, so
  # it can never overwrite or be confused with a default full-run set.
  output_stem <- if (exact_mode) cfg$kegg_id else paths$prefix
  node_path <- file.path(paths$painter_dir, paste0(output_stem, "_contrast_kegg_pathway_painter_nodes.tsv"))
  index_path <- file.path(paths$painter_dir, paste0(output_stem, "_contrast_kegg_pathway_painter_index.tsv"))
  readme_path <- file.path(paths$painter_dir,
    if (exact_mode) sprintf("README_contrast_kegg_painter_%s.txt", cfg$kegg_id)
    else "README_contrast_kegg_pathway_painter.txt")
  write_tsv(node_table, node_path)
  write_tsv(index, index_path)
  write_readme(readme_path, cfg, row, sum(index$output_png != ""))
  if (nrow(selected) > 0 && !any(index$output_png != "")) {
    stop("LISA-KEGG-016 all selected KEGG pathway renders failed; inspect the node audit error column.", call. = FALSE)
  }

  message(sprintf("Wrote contrast KEGG pathway painter outputs: %s", paths$painter_dir))
  message(sprintf("Painted pathway maps: %d", sum(index$output_png != "")))
}

if (sys.nframe() == 0L) main()
