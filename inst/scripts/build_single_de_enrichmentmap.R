#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_
script_file <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
)
source(file.path(dirname(normalizePath(script_file)), "lisa_plot_metadata.R"))

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

save_plot_svg <- function(path, plot, width, height) {
  grDevices::svg(filename = path, width = width, height = height, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

ensure_package <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  stop(sprintf("Required R package is not available: %s", pkg), call. = FALSE)
}

parse_args <- function(args) {
  out <- list(
    project_dir = default_project_dir,
    analysis_id = NA_character_,
    universe = "GOBP-C2",
    min_jaccard = 0.10,
    max_nodes = 60,
    max_edges = 180,
    top_genes_per_category = 200,
    render_graphs = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("max_nodes", "max_edges", "top_genes_per_category")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("min_jaccard")) {
      out[[key]] <- as.numeric(val)
    } else if (key %in% c("render_graphs")) {
      out[[key]] <- tolower(val) %in% c("true", "t", "1", "yes", "y")
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("Required argument: --project-dir", call. = FALSE)
  if (is.na(out$analysis_id) || out$analysis_id == "") stop("Required argument: --analysis-id", call. = FALSE)
  out
}

paths_for <- function(project_dir, analysis_id, universe) {
  single_dir <- file.path(project_dir, "outputs", "single_de", analysis_id, paste0("collection_", universe))
  prefix <- paste(analysis_id, universe, "enrichmentmap", sep = "_")
  list(
    single_dir = single_dir,
    gsea = file.path(single_dir, "enrichment", paste0(analysis_id, "_GSEA_semantic_annotated.tsv")),
    ora = file.path(single_dir, "enrichment", paste0(analysis_id, "_ORA_semantic_annotated.tsv")),
    out_dir = file.path(project_dir, "outputs", "enrichmentmap", "single_de", analysis_id, paste0("collection_", universe)),
    prefix = prefix
  )
}

safe_num <- function(x, default = NA_real_) {
  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- default
  out
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

classified_rows <- function(df, layer) {
  if (!nrow(df)) return(df)
  if (!"category_id" %in% names(df)) {
    stop(sprintf("%s enrichment rows lack category_id.", layer), call. = FALSE)
  }
  category <- safe_chr(df$category_id)
  keep <- category != "" & category != "OTHER_UNCLASSIFIED"
  if ("classification_status" %in% names(df)) {
    keep <- keep & safe_chr(df$classification_status) == "classified"
  }
  out <- df[keep, , drop = FALSE]
  if (nrow(out)) {
    leaked <- safe_chr(out$category_id) == "OTHER_UNCLASSIFIED"
    if ("classification_status" %in% names(out)) {
      leaked <- leaked | safe_chr(out$classification_status) != "classified"
    }
    if (any(leaked)) {
      stop("Unclassified gene sets reached EnrichmentMap input.", call. = FALSE)
    }
  }
  out
}

split_genes <- function(x) {
  genes <- unique(unlist(strsplit(paste(safe_chr(x), collapse = "/"), "[/;,|[:space:]]+")))
  genes <- toupper(trimws(genes))
  genes[genes != ""]
}

category_gene_sets <- function(df, layer, cfg) {
  if (nrow(df) == 0 || !"category_id" %in% names(df)) {
    return(list())
  }
  gene_col <- if (layer == "GSEA") "leadingEdge" else "genes"
  if (!gene_col %in% names(df)) {
    return(list())
  }
  lapply(split(df, safe_chr(df$category_id)), function(rows) {
    genes <- split_genes(rows[[gene_col]])
    if (length(genes) > cfg$top_genes_per_category) {
      genes <- genes[seq_len(cfg$top_genes_per_category)]
    }
    genes
  })
}

category_meta <- function(df, category_id) {
  rows <- df[safe_chr(df$category_id) == category_id, , drop = FALSE]
  if (nrow(rows) == 0) rows <- df[1, , drop = FALSE]
  first <- rows[1, , drop = FALSE]
  data.frame(
    category_id = category_id,
    label = if ("category_display_name" %in% names(first)) first$category_display_name[[1]] else category_id,
    macrogroup_id = if ("macrogroup_id" %in% names(first)) first$macrogroup_id[[1]] else "",
    macrogroup_name = if ("macrogroup_name" %in% names(first)) first$macrogroup_name[[1]] else "",
    color = if ("color" %in% names(first)) first$color[[1]] else "#7A8790",
    n_terms = nrow(rows),
    min_padj = if ("padj" %in% names(rows)) suppressWarnings(min(safe_num(rows$padj), na.rm = TRUE)) else NA_real_,
    mean_NES = if ("NES" %in% names(rows)) suppressWarnings(mean(safe_num(rows$NES), na.rm = TRUE)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

build_nodes_edges <- function(df, gene_sets, cfg) {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) > 0]
  if (length(gene_sets) == 0) {
    return(list(nodes = data.frame(), edges = data.frame()))
  }
  node_ids <- names(gene_sets)
  node_scores <- vapply(node_ids, function(id) {
    rows <- df[safe_chr(df$category_id) == id, , drop = FALSE]
    nrow(rows) + length(gene_sets[[id]]) / 1000
  }, numeric(1))
  node_ids <- node_ids[order(-node_scores, node_ids)]
  if (length(node_ids) > cfg$max_nodes) {
    node_ids <- node_ids[seq_len(cfg$max_nodes)]
  }
  gene_sets <- gene_sets[node_ids]
  nodes <- do.call(rbind, lapply(node_ids, function(id) {
    cbind(category_meta(df, id), n_genes = length(gene_sets[[id]]))
  }))
  edges <- list()
  if (length(node_ids) >= 2) {
    for (i in seq_len(length(node_ids) - 1L)) {
      for (j in seq.int(i + 1L, length(node_ids))) {
        a <- gene_sets[[node_ids[[i]]]]
        b <- gene_sets[[node_ids[[j]]]]
        inter <- intersect(a, b)
        union <- union(a, b)
        jaccard <- if (length(union) > 0) length(inter) / length(union) else 0
        if (jaccard >= cfg$min_jaccard) {
          edges[[length(edges) + 1L]] <- data.frame(
            source = node_ids[[i]],
            target = node_ids[[j]],
            shared_genes = length(inter),
            jaccard = jaccard,
            overlap_genes = paste(inter, collapse = ";"),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  edges <- if (length(edges) > 0) do.call(rbind, edges) else data.frame(
    source = character(), target = character(), shared_genes = integer(), jaccard = numeric(), overlap_genes = character()
  )
  if (nrow(edges) > cfg$max_edges) {
    edges <- edges[order(-edges$jaccard, -edges$shared_genes, edges$source, edges$target), , drop = FALSE]
    edges <- edges[seq_len(cfg$max_edges), , drop = FALSE]
  }
  list(nodes = nodes, edges = edges)
}

write_gmt <- function(gene_sets, path) {
  lines <- vapply(names(gene_sets), function(id) {
    paste(c(id, id, gene_sets[[id]]), collapse = "\t")
  }, character(1))
  writeLines(lines, path, useBytes = TRUE)
}

circle_layout <- function(nodes) {
  n <- nrow(nodes)
  theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
  nodes$x <- cos(theta)
  nodes$y <- sin(theta)
  nodes
}

plot_enrichmentmap <- function(nodes, edges, cfg, layer, empty_reason = "") {
  if (nrow(nodes) == 0) {
    return(
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_blank() +
        ggplot2::annotate("text", x = 1, y = 1.06, label = "No EnrichmentMap graph", fontface = "bold", size = 5, color = "grey25") +
        ggplot2::annotate("text", x = 1, y = 0.95, label = empty_reason, size = 3.3, color = "grey40") +
        ggplot2::labs(
          title = "LISA category EnrichmentMap",
          subtitle = lisa_plot_subtitle(
            cfg$plot_metadata,
            suffix = sprintf("%s | %s", cfg$universe, layer)
          )
        ) +
        ggplot2::theme_void(base_size = 10) +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
    )
  }
  nodes <- circle_layout(nodes)
  edge_plot <- merge(edges, nodes[, c("category_id", "x", "y")], by.x = "source", by.y = "category_id", all.x = TRUE)
  edge_plot <- merge(edge_plot, nodes[, c("category_id", "x", "y")], by.x = "target", by.y = "category_id", all.x = TRUE, suffixes = c("", "_target"))
  nodes$node_size <- pmax(3, sqrt(safe_num(nodes$n_genes, 1)))
  p <- ggplot2::ggplot() +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = "LISA category EnrichmentMap",
      subtitle = lisa_plot_subtitle(
        cfg$plot_metadata,
        suffix = sprintf("%s | %s | category overlap graph", cfg$universe, layer)
      )
    )
  if (nrow(edge_plot) > 0) {
    p <- p + ggplot2::geom_segment(
      data = edge_plot,
      ggplot2::aes(x = x, y = y, xend = x_target, yend = y_target, linewidth = jaccard),
      color = "#9EADB7",
      alpha = 0.55
    ) +
      ggplot2::scale_linewidth_continuous(range = c(0.25, 2.2))
  }
  p +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, size = node_size),
      fill = nodes$color,
      color = "#1b252d",
      shape = 21,
      stroke = 0.35,
      alpha = 0.95
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x * 1.13, y = y * 1.13, label = label),
      size = 2.6,
      color = "#1b252d",
      check_overlap = TRUE
    )
}

write_layer <- function(df, layer, paths, cfg) {
  gene_sets <- category_gene_sets(df, layer, cfg)
  built <- build_nodes_edges(df, gene_sets, cfg)
  nodes <- built$nodes
  edges <- built$edges
  prefix <- paste0(paths$prefix, "_", layer)
  nodes_path <- file.path(paths$out_dir, paste0(prefix, "_enrichmentmap_nodes.tsv"))
  edges_path <- file.path(paths$out_dir, paste0(prefix, "_enrichmentmap_edges.tsv"))
  gmt_path <- file.path(paths$out_dir, paste0(prefix, "_category_gene_sets.gmt"))
  results_path <- file.path(paths$out_dir, paste0(prefix, "_enrichment_results.tsv"))
  png_path <- file.path(paths$out_dir, paste0(prefix, "_enrichmentmap.png"))
  pdf_path <- file.path(paths$out_dir, paste0(prefix, "_enrichmentmap.pdf"))
  svg_path <- file.path(paths$out_dir, paste0(prefix, "_enrichmentmap.svg"))
  write_tsv(nodes, nodes_path)
  write_tsv(edges, edges_path)
  write_tsv(df, results_path)
  write_gmt(gene_sets, gmt_path)
  status <- if (nrow(nodes) > 0) "completed" else "empty"
  reason <- if (nrow(nodes) > 0) "" else sprintf("No category gene sets available for %s.", layer)
  if (isTRUE(cfg$render_graphs)) {
    plot <- plot_enrichmentmap(nodes, edges, cfg, layer, reason)
    height <- if (nrow(nodes) > 35) 9.5 else 7.2
    ggplot2::ggsave(png_path, plot, width = 9.5, height = height, dpi = 220, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(pdf_path, plot, width = 9.5, height = height, bg = "white", limitsize = FALSE)
    save_plot_svg(svg_path, plot, width = 9.5, height = height)
  } else {
    png_path <- ""
    pdf_path <- ""
    svg_path <- ""
  }
  data.frame(
    analysis_id = cfg$analysis_id,
    universe = cfg$universe,
    layer = layer,
    status = status,
    empty_reason = reason,
    n_nodes = nrow(nodes),
    n_edges = nrow(edges),
    min_jaccard = cfg$min_jaccard,
    nodes_tsv = nodes_path,
    edges_tsv = edges_path,
    category_gene_sets_gmt = gmt_path,
    enrichment_results_tsv = results_path,
    plot_png = png_path,
    plot_pdf = pdf_path,
    plot_svg = svg_path,
    stringsAsFactors = FALSE
  )
}

write_readme <- function(path) {
  writeLines(c(
    "LISA category EnrichmentMap exports",
    "",
    "This layer creates EnrichmentMap-style category networks for LISA categories.",
    "Nodes are LISA categories. Edges connect categories with overlapping supporting genes.",
    "The output includes Cytoscape-friendly nodes/edges tables, category gene-set GMT files,",
    "the enrichment table used for the graph, and default PNG/PDF/SVG preview graphics.",
    "",
    "This is a network-rendering and interoperability layer. It does not replace the",
    "LISA semantic category summaries, gene cards, volcano overlays or recurrent gene screens."
  ), path, useBytes = TRUE)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  cfg$plot_metadata <- lisa_plot_metadata(cfg$project_dir, cfg$analysis_id)
  if (isTRUE(cfg$render_graphs)) ensure_package("ggplot2")
  paths <- paths_for(cfg$project_dir, cfg$analysis_id, cfg$universe)
  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)
  index <- list()
  if (file.exists(paths$gsea)) {
    gsea <- classified_rows(read_tsv(paths$gsea), "GSEA")
    index[[length(index) + 1L]] <- write_layer(gsea, "GSEA", paths, cfg)
  }
  if (file.exists(paths$ora)) {
    ora <- classified_rows(read_tsv(paths$ora), "ORA")
    if (nrow(ora) > 0) {
      index[[length(index) + 1L]] <- write_layer(ora, "ORA", paths, cfg)
    }
  }
  if (length(index) == 0) {
    empty <- data.frame()
    index[[1L]] <- write_layer(empty, "GSEA", paths, cfg)
  }
  index <- do.call(rbind, index)
  index_path <- file.path(paths$out_dir, paste0(paths$prefix, "_enrichmentmap_index.tsv"))
  write_tsv(index, index_path)
  write_readme(file.path(paths$out_dir, "README_enrichmentmap.txt"))
  print(index[, c("analysis_id", "universe", "layer", "status", "n_nodes", "n_edges")])
}

if (sys.nframe() == 0L) main()
