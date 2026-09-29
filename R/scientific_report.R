# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# R-native production-v1 style report renderer.
#
# This renderer packages an existing checksummed report packaging fixture. It does not
# recalculate LISA results. The visual and navigation reference is
# installed build_LISA_report.R / production-v1; report foundation is not a UI reference.

lisa_scientific_report_collection_order <- function() {
  lisa_default_collections()
}

lisa_scientific_report_esc <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

lisa_scientific_report_slug <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", as.character(x)))
  x <- gsub("(^-+|-+$)", "", x)
  ifelse(nzchar(x), x, "section")
}

lisa_scientific_report_short_slug <- function(x, width = 18L) {
  x <- lisa_scientific_report_slug(x)
  substr(x, 1L, width)
}

lisa_scientific_report_source_name <- function(role, table_id) {
  role <- strsplit(ifelse(is.na(role), "", as.character(role)), "\\|", perl = TRUE)[[1]][[1]]
  stem <- tools::file_path_sans_ext(basename(role))
  known <- c(
    "semantic_GSEA_category_summary", "semantic_ORA_category_summary",
    "GSEA_semantic_annotated", "ORA_semantic_annotated", "semantic_integrated_summary",
    "leading_edge_gene_pathways", "gene_category_contributions",
    "category_gene_support_summary", "recurrent_gene_screen", "top_recurrent_genes",
    "paired_gene_evidence", "contrast_category_gene_summary", "network_edges",
    "kegg_pathway_painter_nodes", "kegg_pathway_painter_index",
    "contrast_kegg_pathway_painter_nodes", "contrast_kegg_pathway_painter_index",
    "standardized_DE", "ranked_genes", "category_gene_cards_index",
    "category_volcano_overlays_index", "leading_edge_gene_heatmaps_index"
  )
  hit <- known[vapply(known, function(x) grepl(x, stem, fixed = TRUE), logical(1))]
  label <- if (length(hit)) hit[[length(hit)]] else stem
  label <- gsub("[^A-Za-z0-9_-]+", "_", label)
  label <- substr(label, 1L, 34L)
  paste0(table_id, "_", label, ".tsv")
}

lisa_scientific_report_source_rel <- function(row) {
  owner_kind <- if (nzchar(row$contrast_id[[1]])) "contrasts" else if (nzchar(row$analysis_id[[1]])) "single_de" else "shared"
  owner <- if (owner_kind == "contrasts") row$contrast_id[[1]] else if (owner_kind == "single_de") row$analysis_id[[1]] else "study"
  collection <- if (nzchar(row$collection[[1]])) row$collection[[1]] else "shared"
  file.path("source_data", owner_kind, lisa_scientific_report_short_slug(owner),
            lisa_scientific_report_short_slug(collection, 14L),
            lisa_scientific_report_source_name(row$role[[1]], row$table_id[[1]]))
}

lisa_scientific_report_read_tsv <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
    return(data.frame())
  }
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
                    comment.char = "", check.names = FALSE,
                    stringsAsFactors = FALSE)
}

lisa_scientific_report_sha256 <- lisa_sha256_file

lisa_scientific_report_alias <- function(x, key) {
  x <- ifelse(is.na(x), "", as.character(x))
  vapply(strsplit(x, "\\|", fixed = FALSE), function(z) key %in% z, logical(1))
}

lisa_scientific_report_classify <- function(role, artifact_type, path = "") {
  role <- tolower(as.character(role)); path <- tolower(as.character(path))
  if (tolower(as.character(artifact_type)) == "table") return("table")
  if (grepl("enrichmentmap", role) || grepl("enrichmentmap", path)) return("enrichment_maps")
  if (grepl("contrast.*painted|contrast_kegg_pathway_painter", role) || grepl("contrast.*painted", path)) return("contrast_kegg_painted_maps")
  if (grepl("contrast_dumbbell|dumbbell", role) || grepl("dumbbell", path)) return("category_shifts")
  if (grepl("gene_category_network", role) || grepl("gene_category_network", path)) return("gene_category_networks")
  if (grepl("contrast_category_card", role) || grepl("contrast_category_card", path)) return("contrast_gene_cards")
  if (grepl("paired_gene_heatmap", role) || grepl("paired_gene_heatmap", path)) return("paired_gene_heatmaps")
  if (grepl("gsea_category_pathways", role) || grepl("gsea_category_pathways", path)) return("lisa_category_gene_sets")
  if (grepl("gsea_lollipop|gsea_pathway_dotplot|ora_barplot", role) ||
      grepl("gsea_lollipop|gsea_pathway_dotplot|ora_barplot", path)) return("lisa_summary")
  if (grepl("recurrent", role) || grepl("recurrent", path)) return("recurrent_genes")
  if (grepl("kegg_pathway_painter|_painted", role) || grepl("kegg_pathway_painter|_painted", path)) return("kegg_painted_maps")
  if (grepl("heatmap|leading_edge", role) || grepl("heatmap|leading_edge", path)) return("supporting_gene_heatmaps")
  if (grepl("volcano", role) || grepl("volcano", path)) return("volcano_overlays")
  if (grepl("gene_card|genecard|category_card", role) || grepl("gene_card|genecard|category_card", path)) return("gene_prioritization")
  NA_character_
}

lisa_scientific_report_layer_label <- function(family, contrast = FALSE) {
  single <- c(
    lisa_summary = "LISA summary",
    lisa_category_gene_sets = "Gene sets by LISA category",
    gene_prioritization = "Gene prioritization",
    recurrent_genes = "Cross-category recurrent genes",
    volcano_overlays = "Volcano overlays",
    supporting_gene_heatmaps = "Supporting-gene heatmaps",
    enrichment_maps = "Enrichment maps",
    kegg_painted_maps = "KEGG pathway painted maps",
    category_shifts = "LISA category shifts",
    contrast_gene_cards = "Contrast GeneCards",
    paired_gene_heatmaps = "Paired gene heatmaps",
    gene_category_networks = "Gene-category networks",
    contrast_kegg_painted_maps = "Contrast KEGG pathway painted maps"
  )
  out <- single[[family]]
  if (is.null(out)) stop("LISA-REPORT-SCIENCE-016 unlabelled scientific family: ", family, call. = FALSE)
  unname(out)
}

lisa_scientific_report_layer_note <- function(family) {
  notes <- c(
    gene_prioritization = paste0(
      "Genes are ranked within each LISA category by DE-weighted contribution: ",
      "sqrt(LISA support) x |log2FC| x (1 + -log10(DE FDR)). ",
      "This score prioritizes supporting genes; it is not a statistical test."
    ),
    recurrent_genes = paste0(
      "Genes recurring across multiple LISA categories. Cross-category recurrence ",
      "is shown separately from within-category gene prioritization."
    )
  )
  key <- as.character(family)
  if (!length(key) || !key %in% names(notes)) return("")
  unname(notes[[key]])
}

lisa_scientific_report_pretty <- function(path) {
  x <- tools::file_path_sans_ext(basename(as.character(path)))
  x <- gsub("^[0-9]+_", "", x)
  x <- gsub("[_-]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

lisa_scientific_report_logical_stem <- function(path) {
  tools::file_path_sans_ext(basename(as.character(path)))
}

lisa_scientific_report_expand_figure_aliases <- function(figures) {
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(figures))) {
    roles <- strsplit(ifelse(is.na(figures$role[[i]]), "", figures$role[[i]]), "\\|", perl = TRUE)[[1]]
    if (!length(roles)) roles <- ""
    for (role in roles) {
      k <- k + 1L; row <- figures[i, , drop = FALSE]; row$role <- role
      collection <- sub(".*?/collection_([^/]+).*", "\\1", role)
      if (!identical(collection, role)) row$collection <- collection
      analysis <- sub(".*?/single_de/([^/]+).*", "\\1", role)
      if (!identical(analysis, role)) row$analysis_id <- analysis
      contrast <- sub(".*?/category_contrasts/(C[0-9]+).*", "\\1", role)
      if (!identical(contrast, role)) row$contrast_id <- contrast
      rows[[k]] <- row
    }
  }
  do.call(rbind, rows)
}

lisa_scientific_report_context_match <- function(tables, row) {
  keep <- rep(TRUE, nrow(tables))
  for (nm in c("analysis_id", "contrast_id", "collection")) {
    if (!nm %in% names(tables) || !nm %in% names(row)) next
    key <- ifelse(is.na(row[[nm]][[1]]), "", as.character(row[[nm]][[1]]))
    if (nzchar(key)) keep <- keep & lisa_scientific_report_alias(tables[[nm]], key)
  }
  tables[keep, , drop = FALSE]
}

lisa_scientific_report_source_score <- function(role, family) {
  role <- tolower(ifelse(is.na(role), "", as.character(role)))
  score <- rep(0, length(role))
  patterns <- switch(family,
    lisa_summary = c("semantic_gsea_category_summary" = 140, "semantic_ora_category_summary" = 130,
                     "gsea_semantic_annotated" = 120, "ora_semantic_annotated" = 110),
    lisa_category_gene_sets = c("gsea_semantic_annotated" = 150, "semantic_gsea_category_summary" = 100,
                                "category" = 30),
    gene_prioritization = c("gene_category_contributions" = 320,
                            "category_gene_support_summary" = 220,
                            "leading_edge_gene_pathways" = 120,
                            "category_gene_cards_index" = 80),
    recurrent_genes = c("recurrent_gene_screen" = 180, "top_recurrent_genes" = 170,
                        "gene_category_contributions" = 100),
    category_shifts = c("category_contrast" = 100, "contrast" = 30, "category" = 10),
    contrast_gene_cards = c("paired_gene_evidence" = 150, "contrast_category_gene_summary" = 140),
    paired_gene_heatmaps = c("paired_gene_evidence" = 150, "contrast_category_gene_summary" = 130),
    gene_category_networks = c("network_edges" = 150, "paired_gene_evidence" = 120,
                               "contrast_category_gene_summary" = 110),
    volcano_overlays = c("paired_gene_evidence" = 130, "category_summary" = 120,
                         "gene_evidence" = 110, "gene_card" = 90, "de_" = 30),
    supporting_gene_heatmaps = c("paired_gene_evidence" = 130, "gene_evidence" = 120,
                 "leading_edge" = 110, "recurrent" = 90,
                 "category_summary" = 120, "gene_card" = 70),
    enrichment_maps = c("enrichmentmap" = 100, "network" = 70, "category" = 30),
    kegg_painted_maps = c("kegg_pathway_painter_nodes" = 180, "kegg_pathway_painter_index" = 160,
                          "kegg" = 100, "gene_evidence" = 20),
    contrast_kegg_painted_maps = c("contrast_kegg_pathway_painter_nodes" = 180,
                                   "contrast_kegg_pathway_painter_index" = 160,
                                   "paired_gene_evidence" = 80, "kegg" = 100),
    c("category" = 10)
  )
  for (p in names(patterns)) score[grepl(p, role)] <- pmax(score[grepl(p, role)], patterns[[p]])
  score
}

lisa_scientific_report_source_for <- function(group_row, tables) {
  candidates <- lisa_scientific_report_context_match(tables, group_row)
  if (!nrow(candidates)) return(NA_character_)
  scores <- lisa_scientific_report_source_score(candidates$role, as.character(group_row$family[[1]]))
  figure_role <- as.character(group_row$original_role[[1]])
  figure_dir <- dirname(figure_role)
  figure_tokens <- unique(strsplit(tolower(lisa_scientific_report_pretty(figure_role)), " ", fixed = TRUE)[[1]])
  table_dirs <- dirname(as.character(candidates$role))
  scores[table_dirs == figure_dir] <- scores[table_dirs == figure_dir] + 200
  token_overlap <- vapply(candidates$role, function(x) {
    tokens <- unique(strsplit(tolower(lisa_scientific_report_pretty(x)), " ", fixed = TRUE)[[1]])
    length(intersect(figure_tokens[nchar(figure_tokens) > 2L], tokens[nchar(tokens) > 2L]))
  }, integer(1))
  scores <- scores + token_overlap * 5
  ord <- order(-scores, candidates$table_id)
  candidates$table_id[[ord[[1]]]]
}

lisa_scientific_report_group_figures <- function(figures, tables) {
  if (!nrow(figures)) return(data.frame())
  figures <- lisa_scientific_report_expand_figure_aliases(figures)
  figures$family <- mapply(lisa_scientific_report_classify, figures$role,
                           figures$artifact_type, figures$path,
                           USE.NAMES = FALSE)
  if (anyNA(figures$family)) {
    unknown <- unique(figures$role[is.na(figures$family)])
    stop("LISA-REPORT-SCIENCE-017 unclassified scientific figures; add an explicit family instead of an Other bucket: ",
         paste(head(unknown, 20L), collapse = "; "), call. = FALSE)
  }
  # role retains the canonical source path; path is the short physical fixture
  # identifier and therefore cannot be used to group PNG/SVG/PDF variants.
  figures$logical_stem <- vapply(figures$role, lisa_scientific_report_logical_stem, character(1))
  keys <- paste(ifelse(is.na(figures$analysis_id), "", figures$analysis_id),
                ifelse(is.na(figures$contrast_id), "", figures$contrast_id),
                ifelse(is.na(figures$collection), "", figures$collection),
                figures$family, figures$logical_stem, sep = "\034")
  split_rows <- split(seq_len(nrow(figures)), keys)
  rows <- lapply(seq_along(split_rows), function(i) {
    z <- figures[split_rows[[i]], , drop = FALSE]
    ext <- tolower(tools::file_ext(z$path))
    value <- function(e) {
      hit <- which(ext == e)
      if (length(hit)) as.character(z$output_path[[hit[[1]]]]) else ""
    }
    preview <- value("png")
    if (!nzchar(preview)) preview <- value("svg")
    if (!nzchar(preview)) preview <- value("pdf")
    row <- data.frame(
      card_id = sprintf("C%05d", i),
      title = lisa_scientific_report_pretty(z$role[[1]]),
      original_role = as.character(z$role[[1]]),
      family = z$family[[1]],
      analysis_id = ifelse(is.na(z$analysis_id[[1]]), "", z$analysis_id[[1]]),
      contrast_id = ifelse(is.na(z$contrast_id[[1]]), "", z$contrast_id[[1]]),
      collection = ifelse(is.na(z$collection[[1]]), "", z$collection[[1]]),
      preview_path = preview,
      png_path = value("png"), svg_path = value("svg"), pdf_path = value("pdf"),
      source_table_id = "", stringsAsFactors = FALSE
    )
    row$source_table_id <- lisa_scientific_report_source_for(row, tables)
    row
  })
  do.call(rbind, rows)
}

lisa_scientific_report_root_layout_ok <- function(root) {
  allowed <- c("index.html", "README.md", "PROJECT.md", "assets", "pages",
               "media", "source_data", "metadata")
  all(file.exists(file.path(root, allowed))) &&
    !length(setdiff(list.files(root, all.files = TRUE, no.. = TRUE), allowed))
}

lisa_scientific_report_rel <- function(target, page_depth) {
  paste0(strrep("../", page_depth), target)
}

lisa_scientific_report_assets <- function(output_dir, compact_logo, rectangular_logo) {
  file.copy(compact_logo, file.path(output_dir, "assets", "logo.svg"), overwrite = TRUE)
  file.copy(rectangular_logo, file.path(output_dir, "assets", "lisa-logo-rectangular.svg"), overwrite = TRUE)
  css <- c(
    ":root{--red:#9e2936;--red2:#c85b4c;--nav:#18232c;--ink:#17232b;--muted:#687780;--line:#d9e0e4;--paper:#f4f6f7;--green:#3f876b}",
    "*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;font-family:Inter,Segoe UI,Arial,sans-serif;color:var(--ink);background:var(--paper)}a{color:#8f2d2a;text-decoration:none}a:hover{text-decoration:underline}",
    ".layout{display:grid;grid-template-columns:360px minmax(0,1fr);min-height:100vh}.sidebar{background:var(--nav);color:#fff;padding:22px 18px;position:sticky;top:0;height:100vh;overflow:auto}.brand{display:flex;align-items:center;gap:12px;margin-bottom:18px}.brand img{width:84px;height:84px}.brand strong{display:block;font-size:17px}.brand span{display:block;font-size:12px;color:#b9c0c5;margin-top:2px}",
    ".nav-section{border-top:1px solid #32404a;padding-top:10px;margin-top:12px}.nav-label{margin:0 6px 7px;color:#d3dde4;font-size:11px;text-transform:uppercase;letter-spacing:.08em;font-weight:800}.nav a{display:block;color:#d9e0e4;padding:7px 10px;border-radius:6px;margin:2px 0;font-size:13px;line-height:1.25}.nav a.active,.nav a:hover{background:#2b3942;color:#fff;text-decoration:none}.nav details{margin:4px 0}.nav summary{cursor:pointer;color:#f0d9d4;font-size:12px;padding:5px 8px}.nav-sub{margin-left:10px;padding-left:8px;border-left:1px solid #33424b}.nav-sub a{font-size:12px;padding:5px 8px}",
    ".content{padding:28px 34px 48px;max-width:1660px;width:100%}.topbar{display:flex;justify-content:space-between;gap:20px;align-items:center;margin-bottom:18px}.breadcrumb{font-size:13px;color:var(--muted)}.pill{display:inline-flex;padding:5px 9px;border-radius:999px;background:#e9f3ef;color:#24634d;font-size:12px;font-weight:700}",
    ".hero,.panel{background:#fff;border:1px solid var(--line);border-radius:9px;margin-bottom:20px}.hero{border-left:5px solid var(--red);padding:22px 24px}.hero-logo{display:block;width:min(760px,100%);height:auto;margin:0 0 18px}.hero h1{margin:0 0 7px;font-size:29px}.hero p{margin:0;color:var(--muted);line-height:1.5}.panel{padding:18px 20px}.section-head{display:flex;justify-content:space-between;gap:16px;align-items:flex-start;margin-bottom:14px}.section-head h2{margin:0;font-size:20px}.eyebrow{text-transform:uppercase;letter-spacing:.08em;color:var(--red);font-weight:800;font-size:11px;margin:0 0 5px}",
    ".metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}.metric{background:#fff;border:1px solid var(--line);border-radius:8px;padding:14px}.metric strong{font-size:25px;display:block}.metric span{font-size:12px;color:var(--muted)}.chips{display:flex;flex-wrap:wrap;gap:8px}.chip{border:1px solid var(--line);background:#fff;border-radius:999px;padding:6px 10px;font-size:12px}",
    ".detail-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:12px}.detail-card{display:block;border:1px solid var(--line);border-left:4px solid var(--red);background:#fff;border-radius:8px;padding:15px}.detail-card:hover{background:#fff8f7;text-decoration:none}",
    ".image-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px}.fig-card{border:1px solid var(--line);border-radius:8px;background:#fff;overflow:hidden}.fig-open{display:block;border:0;background:#fff;padding:0;width:100%;cursor:zoom-in}.fig-card img{width:100%;height:250px;object-fit:contain;background:#fafafa;border-bottom:1px solid var(--line)}.fig-caption{padding:11px}.fig-caption strong{display:block;font-size:13px;line-height:1.3;overflow-wrap:anywhere}.fig-meta{color:var(--muted);font-size:11px;margin-top:4px}.fig-actions{display:flex;gap:6px;flex-wrap:wrap;margin-top:9px}.file-link{display:inline-flex;border:1px solid var(--line);background:#fff;border-radius:6px;padding:5px 8px;font-size:12px}.source-link{border-color:#d7b2b6;color:#862331;font-weight:700}",
    "table{border-collapse:collapse;width:100%;font-size:12px}th,td{border-bottom:1px solid var(--line);padding:7px 8px;text-align:left;vertical-align:top}th{background:#f6f7f8;position:sticky;top:0}.table-wrap{overflow:auto;max-height:70vh}.pager{display:flex;gap:6px;flex-wrap:wrap;margin:16px 0}.pager a,.pager strong{padding:5px 9px;border:1px solid var(--line);border-radius:6px;background:#fff}",
    ".modal{position:fixed;inset:0;background:rgba(10,16,20,.84);display:none;align-items:center;justify-content:center;z-index:1000;padding:24px}.modal.open{display:flex}.modal-inner{background:#fff;border-radius:8px;max-width:96vw;max-height:94vh;padding:12px}.modal img{max-width:92vw;max-height:78vh;display:block}.modal-head{display:flex;justify-content:space-between;gap:12px;align-items:center;margin-bottom:8px}.close{border:1px solid var(--line);background:#fff;border-radius:6px;padding:6px 10px;cursor:pointer}",
    ".warning{border-left:4px solid #c47a24;background:#fff8ed;padding:12px 14px}.ok{border-left:4px solid var(--green);background:#edf7f2;padding:12px 14px}pre{white-space:pre-wrap;overflow-wrap:anywhere}.muted{color:var(--muted)}",
    "@media(max-width:900px){.layout{grid-template-columns:1fr}.sidebar{position:relative;height:auto}.content{padding:20px}.topbar{align-items:flex-start;flex-direction:column}.section-head{display:block}.image-grid{grid-template-columns:1fr}}"
  )
  writeLines(paste(css, collapse = ""), file.path(output_dir, "assets", "report.css"), useBytes = TRUE)
  js <- c(
    "function openFigure(src,title,png,svg,pdf){const m=document.getElementById('modal');m.querySelector('img').src=src;m.querySelector('[data-title]').textContent=title||src;[['png',png],['svg',svg],['pdf',pdf]].forEach(x=>{const a=m.querySelector('[data-'+x[0]+']');a.href=x[1]||'#';a.style.display=x[1]?'inline-flex':'none'});m.classList.add('open')}",
    "function closeFigure(){document.getElementById('modal').classList.remove('open')}",
    "document.addEventListener('keydown',e=>{if(e.key==='Escape')closeFigure()})"
  )
  writeLines(js, file.path(output_dir, "assets", "report.js"), useBytes = TRUE)
}

lisa_scientific_report_nav <- function(nav, active, prefix) {
  link <- function(label, href, key) sprintf('<a class="%s" href="%s">%s</a>',
    if (identical(active, key)) "active" else "", lisa_scientific_report_esc(paste0(prefix, href)), lisa_scientific_report_esc(label))
  de <- paste(vapply(nav$de, function(x) {
    children <- paste(vapply(x$groups, function(g) {
      links <- paste(vapply(g$children, function(y) sprintf('<a href="%s">%s</a>',
        lisa_scientific_report_esc(paste0(prefix, y$href)), lisa_scientific_report_esc(y$label)), character(1)), collapse = "")
      sprintf('<details><summary>%s</summary><div class="nav-sub">%s</div></details>', lisa_scientific_report_esc(g$label), links)
    }, character(1)), collapse = "")
    sprintf('<details%s><summary>%s</summary><div class="nav-sub">%s</div></details>',
      if (identical(active, x$key)) " open" else "", lisa_scientific_report_esc(x$label), children)
  }, character(1)), collapse = "")
  cx <- paste(vapply(nav$cx, function(x) {
    children <- paste(vapply(x$groups, function(g) {
      links <- paste(vapply(g$children, function(y) sprintf('<a href="%s">%s</a>',
        lisa_scientific_report_esc(paste0(prefix, y$href)), lisa_scientific_report_esc(y$label)), character(1)), collapse = "")
      sprintf('<details><summary>%s</summary><div class="nav-sub">%s</div></details>', lisa_scientific_report_esc(g$label), links)
    }, character(1)), collapse = "")
    sprintf('<details%s><summary>%s</summary><div class="nav-sub">%s</div></details>',
      if (identical(active, x$key)) " open" else "", lisa_scientific_report_esc(x$label), children)
  }, character(1)), collapse = "")
  paste0(
    '<div class="nav-section"><p class="nav-label">Study</p>', link("Study index", "index.html", "overview"), '</div>',
    '<div class="nav-section"><p class="nav-label">Single DE</p>', de, '</div>',
    '<div class="nav-section"><p class="nav-label">Contrasts</p>', cx, '</div>',
    '<div class="nav-section"><p class="nav-label">Data</p>',
    link("Data downloads", "pages/downloads.html", "downloads"),
    link("Analysis details", "pages/details.html", "details"), '</div>')
}

lisa_scientific_report_shell <- function(title, subtitle, active, body, nav, depth = 1L) {
  prefix <- strrep("../", depth)
  paste0('<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>', lisa_scientific_report_esc(title), ' | LISA</title>',
    '<link rel="stylesheet" href="', prefix, 'assets/report.css"></head><body>',
    '<div class="layout"><aside class="sidebar"><div class="brand"><img src="', prefix,
    'assets/logo.svg" alt="LISA"><div><strong>LISA report</strong><span>LLM-Inferred Semantic Annotation</span><span>',
    lisa_scientific_report_esc(subtitle), '</span></div></div><nav class="nav">',
    lisa_scientific_report_nav(nav, active, prefix), '</nav></aside><main class="content">',
    body, '</main></div>',
    '<div id="modal" class="modal" onclick="if(event.target.id===\'modal\')closeFigure()">',
    '<div class="modal-inner"><div class="modal-head"><strong data-title></strong><div class="fig-actions">',
    '<a class="file-link" data-png download>PNG</a><a class="file-link" data-svg download>SVG</a>',
    '<a class="file-link" data-pdf download>PDF</a><button class="close" onclick="closeFigure()">Close</button>',
    '</div></div><img alt="Expanded figure"></div></div>',
    '<script src="', prefix, 'assets/report.js"></script></body></html>')
}

lisa_scientific_report_card_html <- function(row, depth = 1L) {
  prefix <- strrep("../", depth)
  p <- function(x) ifelse(is.na(x) || !nzchar(x), "", paste0(prefix, x))
  actions <- c()
  for (ext in c("png", "svg", "pdf")) {
    value <- as.character(row[[paste0(ext, "_path")]][[1]])
    if (!is.na(value) && nzchar(value)) actions <- c(actions,
      sprintf('<a class="file-link" download href="%s">%s</a>', lisa_scientific_report_esc(p(value)), toupper(ext)))
  }
  actions <- c(actions, sprintf('<a class="file-link source-link" download href="%s">Source data</a>',
    lisa_scientific_report_esc(p(row$source_table_path[[1]]))))
  preview_is_pdf <- grepl("\\.pdf$", row$preview_path[[1]], ignore.case = TRUE)
  onclick <- sprintf("openFigure('%s','%s','%s','%s','%s')",
    lisa_scientific_report_esc(p(row$preview_path[[1]])),
    lisa_scientific_report_esc(row$title[[1]]),
    lisa_scientific_report_esc(p(row$png_path[[1]])), lisa_scientific_report_esc(p(row$svg_path[[1]])),
    lisa_scientific_report_esc(p(row$pdf_path[[1]])))
  visual <- if (preview_is_pdf) paste0('<object data="', lisa_scientific_report_esc(p(row$preview_path[[1]])),
    '" type="application/pdf" width="100%" height="250"><a class="file-link" href="',
    lisa_scientific_report_esc(p(row$preview_path[[1]])), '">Open PDF preview</a></object>') else paste0(
    '<button class="fig-open" onclick="', onclick, '"><img loading="lazy" src="',
    lisa_scientific_report_esc(p(row$preview_path[[1]])), '" alt="', lisa_scientific_report_esc(row$title[[1]]), '"></button>')
  paste0('<article class="fig-card">', visual, '<div class="fig-caption"><strong>',
    lisa_scientific_report_esc(row$title[[1]]), '</strong><div class="fig-meta">',
    lisa_scientific_report_esc(paste(row$collection[[1]], lisa_scientific_report_layer_label(row$family[[1]], nzchar(row$contrast_id[[1]])), sep = " | ")),
    '</div><div class="fig-actions">', paste(actions, collapse = ""), '</div></div></article>')
}

lisa_scientific_report_table_html <- function(df, max_rows = 200L) {
  if (!nrow(df)) return('<p class="muted">No rows.</p>')
  show <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  head <- paste0("<tr>", paste0("<th>", lisa_scientific_report_esc(names(show)), "</th>", collapse = ""), "</tr>")
  rows <- paste(vapply(seq_len(nrow(show)), function(i) paste0("<tr>",
    paste0("<td>", lisa_scientific_report_esc(show[i, , drop = TRUE]), "</td>", collapse = ""), "</tr>"), character(1)), collapse = "")
  note <- if (nrow(df) > max_rows) sprintf('<p class="muted">Showing %d of %d rows.</p>', max_rows, nrow(df)) else ""
  paste0(note, '<div class="table-wrap"><table>', head, rows, '</table></div>')
}

lisa_scientific_report_write <- function(path, text) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, path, useBytes = TRUE)
}

lisa_scientific_report_verify_links <- function(root) {
  html <- list.files(root, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
  broken <- character()
  for (page in html) {
    txt <- paste(readLines(page, warn = FALSE), collapse = "\n")
    hits <- regmatches(txt, gregexpr('(?:href|src)="[^"]+"', txt, perl = TRUE))[[1]]
    if (!length(hits) || identical(hits, character(0))) next
    refs <- sub('^[^=]+="', "", hits); refs <- sub('"$', "", refs)
    refs <- refs[!grepl("^(#|https?:|mailto:|javascript:|data:)", refs)]
    refs <- sub("#.*$", "", refs); refs <- refs[nzchar(refs)]
    for (ref in refs) if (!file.exists(file.path(dirname(page), utils::URLdecode(ref)))) {
      broken <- c(broken, paste0(sub(paste0("^", root, "/?"), "", page), " -> ", ref))
    }
  }
  unique(broken)
}

lisa_scientific_report_make_nav <- function(cards, de_index, contrast_index, page_index) {
  label_for <- function(df, id, id_col, label_col) {
    if (!nrow(df) || !id_col %in% names(df)) return(id)
    hit <- which(df[[id_col]] == id)
    if (length(hit) && label_col %in% names(df)) as.character(df[[label_col]][hit[[1]]]) else id
  }
  de_ids <- unique(cards$analysis_id[nzchar(cards$analysis_id)])
  cx_ids <- unique(cards$contrast_id[nzchar(cards$contrast_id)])
  grouped <- function(pages) {
    pages <- pages[pages$page == 1L, , drop = FALSE]
    collections <- unique(pages$collection)
    lapply(collections, function(collection) {
      z <- pages[pages$collection == collection, , drop = FALSE]
      family_order <- c("lisa_summary", "lisa_category_gene_sets", "gene_prioritization",
        "recurrent_genes", "volcano_overlays", "supporting_gene_heatmaps", "kegg_painted_maps",
        "enrichment_maps", "category_shifts",
        "contrast_gene_cards", "paired_gene_heatmaps", "gene_category_networks",
        "contrast_kegg_painted_maps")
      z <- z[order(match(z$family, family_order)), , drop = FALSE]
      list(label = collection, children = lapply(seq_len(nrow(z)), function(i) {
        list(label = lisa_scientific_report_layer_label(z$family[[i]], z$kind[[i]] == "contrast"), href = z$path[[i]])
      }))
    })
  }
  de <- lapply(de_ids, function(id) {
    pages <- page_index[page_index$kind == "single" & page_index$owner_id == id, , drop = FALSE]
    list(key = paste0("de-", id), label = label_for(de_index, id, "analysis_id", "label"),
         groups = grouped(pages))
  })
  cx <- lapply(cx_ids, function(id) {
    pages <- page_index[page_index$kind == "contrast" & page_index$owner_id == id, , drop = FALSE]
    list(key = paste0("cx-", id), label = label_for(contrast_index, id, "contrast_id", "contrast_label"),
         groups = grouped(pages))
  })
  list(de = de, cx = cx)
}

#' Build a production-v1 style multipage report from a report packaging fixture
#'
#' @param fixture_path Checksummed report packaging fixture root.
#' @param output_dir Empty destination directory.
#' @param title Report title.
#' @param study_label Short sidebar label.
#' @param cards_per_page Maximum logical figure cards per page.
#' @param include_enrichment_maps Show archived EnrichmentMap products. The
#'   standard scientific report keeps this opt-in and defaults to `FALSE`.
#' @return Report entry point and inventory metrics.
#' @keywords internal
build_lisa_scientific_report <- function(fixture_path, output_dir, title = NULL,
                                    study_label = NULL, cards_per_page = 72L,
                                    include_enrichment_maps = FALSE) {
  fixture_path <- normalizePath(fixture_path, winslash = "/", mustWork = TRUE)
  required <- c("m/registry.tsv", "m/sha256.tsv", "m/de.tsv", "m/cx.tsv", "m/col.tsv", "m/prov.tsv")
  missing <- required[!file.exists(file.path(fixture_path, required))]
  if (length(missing)) stop("LISA-REPORT-SCIENCE-002 incomplete fixture: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.numeric(cards_per_page) || length(cards_per_page) != 1L || cards_per_page < 1L || cards_per_page > 200L) {
    stop("LISA-REPORT-SCIENCE-005 cards_per_page must be between 1 and 200.", call. = FALSE)
  }
  if (is.null(title)) title <- basename(fixture_path)
  if (is.null(study_label)) study_label <- title
  output_dir <- normalizePath(path.expand(output_dir), winslash = "/", mustWork = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("LISA-REPORT-SCIENCE-003 destination must be empty.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (d in c("assets", "pages", "media", "source_data", "metadata")) dir.create(file.path(output_dir, d), showWarnings = FALSE)

  registry <- lisa_scientific_report_read_tsv(file.path(fixture_path, "m/registry.tsv"))
  lisa_require_columns(registry, c("table_id", "artifact_type", "path", "role", "sha256"), "scientific report registry")
  if (!nrow(registry) || !is.character(registry$path) || anyNA(registry$path) ||
      any(!nzchar(registry$path)) || anyDuplicated(registry$path) ||
      !lisa_sha256_all_valid(registry$sha256, nrow(registry))) {
    stop("LISA-REPORT-SCIENCE-004 malformed fixture registry or SHA-256 values.", call. = FALSE)
  }
  for (nm in c("analysis_id", "contrast_id", "collection")) if (!nm %in% names(registry)) registry[[nm]] <- ""
  registry$analysis_id[is.na(registry$analysis_id)] <- ""
  registry$contrast_id[is.na(registry$contrast_id)] <- ""
  registry$collection[is.na(registry$collection)] <- ""
  sums <- lisa_scientific_report_read_tsv(file.path(fixture_path, "m/sha256.tsv"))
  lisa_require_columns(sums, c("path", "sha256"), "scientific report SHA-256 manifest")
  actual <- sort(setdiff(list.files(fixture_path, recursive = TRUE, full.names = FALSE), "m/sha256.tsv"))
  if (!is.character(sums$path) || anyNA(sums$path) || any(!nzchar(sums$path)) ||
      anyDuplicated(sums$path) || !lisa_sha256_all_valid(sums$sha256, nrow(sums)) ||
      !identical(sort(sums$path), actual)) {
    stop("LISA-REPORT-SCIENCE-004 malformed or incomplete fixture SHA-256 manifest.", call. = FALSE)
  }
  for (i in seq_len(nrow(sums))) {
    path <- file.path(fixture_path, sums$path[[i]])
    if (!file.exists(path) || !identical(lisa_scientific_report_sha256(path), as.character(sums$sha256[[i]]))) {
      stop("LISA-REPORT-SCIENCE-004 fixture checksum mismatch: ", sums$path[[i]], call. = FALSE)
    }
  }
  matched <- match(registry$path, sums$path)
  if (anyNA(matched) ||
      !identical(unname(registry$sha256), unname(sums$sha256[matched]))) {
    stop("LISA-REPORT-SCIENCE-004 registry hashes do not match the fixture manifest.", call. = FALSE)
  }
  if (!isTRUE(include_enrichment_maps)) {
    registry <- registry[!grepl("enrichmentmap", registry$role, ignore.case = TRUE), , drop = FALSE]
  }

  registry$output_path <- ""
  for (i in seq_len(nrow(registry))) {
    src <- file.path(fixture_path, registry$path[[i]])
    if (!file.exists(src)) stop("LISA-REPORT-SCIENCE-006 registered artifact missing: ", registry$path[[i]], call. = FALSE)
    ext <- tolower(tools::file_ext(src))
    rel <- if (registry$artifact_type[[i]] == "table") {
      lisa_scientific_report_source_rel(registry[i, , drop = FALSE])
    } else file.path("media", paste0(registry$table_id[[i]], ".", ext))
    dir.create(dirname(file.path(output_dir, rel)), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, file.path(output_dir, rel), overwrite = FALSE)) stop("LISA-REPORT-SCIENCE-007 copy failed: ", rel, call. = FALSE)
    registry$output_path[[i]] <- rel
  }
  tables <- registry[registry$artifact_type == "table", , drop = FALSE]
  figures <- registry[registry$artifact_type == "figure", , drop = FALSE]
  cards <- lisa_scientific_report_group_figures(figures, tables)
  if (!nrow(cards)) stop("LISA-REPORT-SCIENCE-008 no previewable figure groups found.", call. = FALSE)
  missing_preview <- !nzchar(cards$preview_path)
  missing_source <- is.na(cards$source_table_id) | !nzchar(cards$source_table_id)
  if (any(missing_preview)) stop("LISA-REPORT-SCIENCE-009 figure groups lack a previewable PNG/SVG/PDF asset: ", paste(cards$card_id[missing_preview], collapse = ", "), call. = FALSE)
  if (any(missing_source)) stop("LISA-REPORT-SCIENCE-010 figure groups lack explicit source data: ", paste(cards$card_id[missing_source], collapse = ", "), call. = FALSE)
  cards$source_table_path <- tables$output_path[match(cards$source_table_id, tables$table_id)]
  if (anyNA(cards$source_table_path) || any(!nzchar(cards$source_table_path))) {
    stop("LISA-REPORT-SCIENCE-010 source-data path resolution failed.", call. = FALSE)
  }

  de_index <- lisa_scientific_report_read_tsv(file.path(fixture_path, "m/de.tsv"))
  contrast_index <- lisa_scientific_report_read_tsv(file.path(fixture_path, "m/cx.tsv"))
  collection_index <- lisa_scientific_report_read_tsv(file.path(fixture_path, "m/col.tsv"))
  provenance <- lisa_scientific_report_read_tsv(file.path(fixture_path, "m/prov.tsv"))

  # Build deterministic page inventory before rendering so every page gets the same nav.
  rows <- list(); n <- 0L
  add_pages <- function(z, kind, owner_id, collection, family) {
    count <- max(1L, ceiling(nrow(z) / cards_per_page))
    for (page in seq_len(count)) {
      n <<- n + 1L
      rows[[n]] <<- data.frame(kind = kind, owner_id = owner_id, collection = collection,
        family = family, page = page, pages = count,
        content_type = "figures", table_ids = "",
        path = sprintf("pages/%s%03d.html", if (kind == "single") "s" else "c", n),
        nav_label = paste(collection, lisa_scientific_report_layer_label(family, kind == "contrast"), sep = " | "),
        stringsAsFactors = FALSE)
    }
  }
  add_table_page <- function(kind, owner_id, collection, family, table_ids) {
    n <<- n + 1L
    rows[[n]] <<- data.frame(kind = kind, owner_id = owner_id, collection = collection,
      family = family, page = 1L, pages = 1L, content_type = "tables",
      table_ids = paste(table_ids, collapse = "|"),
      path = sprintf("pages/%s%03d.html", if (kind == "single") "s" else "c", n),
      nav_label = paste(collection, lisa_scientific_report_layer_label(family, kind == "contrast"), sep = " | "),
      stringsAsFactors = FALSE)
  }
  single_family_order <- c("lisa_summary", "lisa_category_gene_sets", "gene_prioritization",
                           "recurrent_genes", "volcano_overlays", "supporting_gene_heatmaps",
                           "kegg_painted_maps", "enrichment_maps")
  contrast_family_order <- c("category_shifts", "contrast_gene_cards", "paired_gene_heatmaps",
                             "gene_category_networks", "contrast_kegg_painted_maps")
  order_collections <- function(x) {
    known <- lisa_scientific_report_collection_order(); c(intersect(known, x), sort(setdiff(x, known)))
  }
  for (kind in c("single", "contrast")) {
    ids <- if (kind == "single") unique(cards$analysis_id[nzchar(cards$analysis_id)]) else unique(cards$contrast_id[nzchar(cards$contrast_id)])
    for (id in ids) {
      own <- if (kind == "single") cards[cards$analysis_id == id, , drop = FALSE] else cards[cards$contrast_id == id, , drop = FALSE]
      family_order <- if (kind == "single") single_family_order else contrast_family_order
      for (collection in order_collections(unique(own$collection))) for (family in family_order) {
        z <- own[own$collection == collection & own$family == family, , drop = FALSE]
        if (nrow(z)) add_pages(z, kind, id, collection, family)
      }
    }
  }
  page_index <- if (length(rows)) do.call(rbind, rows) else data.frame()
  nav <- lisa_scientific_report_make_nav(cards, de_index, contrast_index, page_index)

  report_asset <- function(name) {
    installed <- system.file("report_assets", name, package = "lisaR")
    if (nzchar(installed) && file.exists(installed)) return(installed)
    local <- file.path(getwd(), "inst", "report_assets", name)
    if (file.exists(local)) return(local)
    ""
  }
  compact_logo <- report_asset("LISA_logo_C_compact_icon_muted_red_S.svg")
  rectangular_logo <- report_asset("LISA_logo_A1_muted_red_S_automated_annotation_final.svg")
  if (!nzchar(compact_logo) || !nzchar(rectangular_logo)) {
    stop("LISA-REPORT-SCIENCE-011 LISA compact or rectangular logo asset is absent.", call. = FALSE)
  }
  lisa_scientific_report_assets(output_dir, compact_logo, rectangular_logo)

  # Render scientific pages.
  for (i in seq_len(nrow(page_index))) {
    p <- page_index[i, , drop = FALSE]
    own <- if (p$kind == "single") cards[cards$analysis_id == p$owner_id, , drop = FALSE] else cards[cards$contrast_id == p$owner_id, , drop = FALSE]
    z <- own[own$collection == p$collection & own$family == p$family, , drop = FALSE]
    if (p$content_type == "figures") {
      from <- (p$page - 1L) * cards_per_page + 1L; to <- min(p$page * cards_per_page, nrow(z))
      part <- z[from:to, , drop = FALSE]
    } else part <- z[0, , drop = FALSE]
    pager <- if (p$pages > 1L) paste0('<div class="pager">', paste(vapply(seq_len(p$pages), function(j) {
      hit <- page_index[page_index$kind == p$kind & page_index$owner_id == p$owner_id &
        page_index$collection == p$collection & page_index$family == p$family & page_index$page == j, , drop = FALSE]
      if (j == p$page) sprintf("<strong>%d</strong>", j) else sprintf('<a href="%s">%d</a>', basename(hit$path[[1]]), j)
    }, character(1)), collapse = ""), '</div>') else ""
    owner_label <- p$owner_id
    if (p$kind == "single" && "label" %in% names(de_index)) {
      hit <- which(de_index$analysis_id == p$owner_id); if (length(hit)) owner_label <- de_index$label[[hit[[1]]]]
    }
    if (p$kind == "contrast" && "contrast_label" %in% names(contrast_index)) {
      hit <- which(contrast_index$contrast_id == p$owner_id); if (length(hit)) owner_label <- contrast_index$contrast_label[[hit[[1]]]]
    }
    content <- if (p$content_type == "figures") {
      paste0(pager, '<section class="image-grid">', paste(vapply(seq_len(nrow(part)), function(j) lisa_scientific_report_card_html(part[j, , drop = FALSE]), character(1)), collapse = ""), '</section>', pager)
    } else {
      ids <- strsplit(p$table_ids, "\\|", perl = TRUE)[[1]]
      data_tables <- tables[tables$table_id %in% ids, , drop = FALSE]
      paste0('<section class="detail-grid">', paste(vapply(seq_len(nrow(data_tables)), function(j) {
        row <- data_tables[j, , drop = FALSE]
        sprintf('<a class="detail-card" download href="../%s"><strong>%s</strong><span>Download TSV</span></a>',
          lisa_scientific_report_esc(row$output_path[[1]]), lisa_scientific_report_esc(lisa_scientific_report_source_name(row$role[[1]], "")))
      }, character(1)), collapse = ""), '</section>')
    }
    layer_note <- lisa_scientific_report_layer_note(p$family)
    body <- paste0('<section class="hero"><h1>', lisa_scientific_report_esc(lisa_scientific_report_layer_label(p$family, p$kind == "contrast")),
      '</h1><p>', lisa_scientific_report_esc(owner_label), ' | ', lisa_scientific_report_esc(p$collection), '</p>',
      if (nzchar(layer_note)) paste0('<p>', lisa_scientific_report_esc(layer_note), '</p>') else "",
      '</section>', content)
    html <- lisa_scientific_report_shell(lisa_scientific_report_layer_label(p$family, p$kind == "contrast"), study_label,
      paste0(if (p$kind == "single") "de-" else "cx-", p$owner_id), body, nav)
    lisa_scientific_report_write(file.path(output_dir, p$path), html)
  }

  metric <- function(value, label) sprintf('<div class="metric"><strong>%s</strong><span>%s</span></div>', value, label)
  de_links <- paste(vapply(nav$de, function(x) sprintf('<a class="detail-card" href="%s"><strong>%s</strong></a>',
    lisa_scientific_report_esc(x$groups[[1]]$children[[1]]$href), lisa_scientific_report_esc(x$label)), character(1)), collapse = "")
  cx_links <- paste(vapply(nav$cx, function(x) sprintf('<a class="detail-card" href="%s"><strong>%s</strong></a>',
    lisa_scientific_report_esc(x$groups[[1]]$children[[1]]$href), lisa_scientific_report_esc(x$label)), character(1)), collapse = "")
  overview <- paste0('<section class="hero"><img class="hero-logo" src="assets/lisa-logo-rectangular.svg" alt="LISA"><p class="eyebrow">LLM-Inferred Semantic Annotation</p><h1>', lisa_scientific_report_esc(title),
    '</h1><p>', lisa_scientific_report_esc(study_label), '</p></section>',
    '<section class="panel"><div class="section-head"><h2>Single DE analyses</h2></div><div class="detail-grid">', de_links, '</div></section>',
    if (nzchar(cx_links)) paste0('<section class="panel"><div class="section-head"><h2>Contrasts</h2></div><div class="detail-grid">', cx_links, '</div></section>') else "",
    '<section class="panel"><div class="section-head"><h2>Collections</h2></div><div class="chips">',
    paste(sprintf('<span class="chip">%s</span>', lisa_scientific_report_esc(order_collections(unique(cards$collection)))), collapse = ""), '</div></section>')
  lisa_scientific_report_write(file.path(output_dir, "index.html"), lisa_scientific_report_shell(title, study_label, "overview", overview, nav, depth = 0L))

  registry_out <- registry[, unique(c("table_id", "artifact_type", "role", "analysis_id", "contrast_id", "collection", "path", "output_path", "sha256")), drop = FALSE]
  utils::write.table(registry_out, file.path(output_dir, "metadata", "registry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(cards, file.path(output_dir, "metadata", "figure_cards.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(page_index, file.path(output_dir, "metadata", "page_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(collection_index, file.path(output_dir, "metadata", "collections.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(provenance, file.path(output_dir, "metadata", "source_provenance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  source_index <- registry_out[registry_out$artifact_type == "table",
    c("analysis_id", "contrast_id", "collection", "role", "output_path"), drop = FALSE]
  names(source_index) <- c("analysis_id", "contrast_id", "collection", "data_description", "file")
  utils::write.table(source_index, file.path(output_dir, "source_data", "INDEX.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

  # The downloads index lists canonical tables. Figure downloads stay beside
  # each visual card, which is more useful and keeps the HTML safely bounded.
  downloads <- registry_out[registry_out$artifact_type == "table", , drop = FALSE]
  downloads$link <- ifelse(downloads$artifact_type == "table",
    paste0('<a class="file-link" download href="../', downloads$output_path, '">TSV</a>'),
    paste0('<a class="file-link" download href="../', downloads$output_path, '">', toupper(tools::file_ext(downloads$output_path)), '</a>'))
  downloads$data_table <- vapply(downloads$role, function(x) lisa_scientific_report_source_name(x, ""), character(1))
  downloads$data_table <- sub("^_", "", downloads$data_table)
  show <- downloads[, c("analysis_id", "contrast_id", "collection", "data_table", "link"), drop = FALSE]
  head <- paste0("<tr>", paste0("<th>", lisa_scientific_report_esc(names(show)), "</th>", collapse = ""), "</tr>")
  rows_html <- paste(vapply(seq_len(nrow(show)), function(i) paste0("<tr>", paste(vapply(names(show), function(nm) {
    if (nm == "link") paste0("<td>", show[[nm]][[i]], "</td>") else paste0("<td>", lisa_scientific_report_esc(show[[nm]][[i]]), "</td>")
  }, character(1)), collapse = ""), "</tr>"), character(1)), collapse = "")
  downloads_body <- paste0('<section class="hero"><h1>Data downloads</h1><p>Tables used to generate and interpret the figures.</p></section><section class="panel"><div class="table-wrap"><table>', head, rows_html, '</table></div></section>')
  lisa_scientific_report_write(file.path(output_dir, "pages", "downloads.html"), lisa_scientific_report_shell("Data downloads", study_label, "downloads", downloads_body, nav))

  qc <- data.frame(check = c("fixture_checksums", "registered_files", "figure_previews", "figure_source_data", "clean_root", "offline_links", "path_budget", "html_budget"),
    status = c("PASS", "PASS", "PASS", "PASS", "PENDING_FINAL", "PENDING_FINAL", "PENDING_FINAL", "PENDING_FINAL"),
    detail = c(sprintf("%d fixture entries verified", nrow(sums)), sprintf("%d artifacts copied", nrow(registry)),
      sprintf("%d cards have inline PNG/SVG/PDF preview", nrow(cards)), sprintf("%d cards have source TSV", nrow(cards)), rep("validated after final manifest", 4)), stringsAsFactors = FALSE)
  utils::write.table(qc, file.path(output_dir, "metadata", "qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # Only expose fields that help biological interpretation. Paths, fixture
  # identifiers and renderer controls remain available under metadata/.
  scientific_columns <- function(x, wanted) {
    keep <- intersect(wanted, names(x))
    if (!length(keep)) return(x[, 0, drop = FALSE])
    x[, keep, drop = FALSE]
  }
  de_details <- scientific_columns(de_index,
    c("analysis_id", "label", "species", "organism", "analysis_role", "fraction", "comparison"))
  contrast_details <- scientific_columns(contrast_index,
    c("contrast_id", "contrast_label", "analysis_a", "label_a", "analysis_b", "label_b", "interpretation"))
  collection_details <- scientific_columns(collection_index,
    c("analysis_id", "contrast_id", "analysis_collection", "collection", "dictionary_id", "dictionary_label"))
  details_body <- paste0('<section class="hero"><h1>Analysis details</h1><p>Study comparisons and analysis collections.</p></section>',
    '<section class="panel"><h2>Single DE analyses</h2>', lisa_scientific_report_table_html(de_details, 500L), '</section>',
    if (nrow(contrast_details)) paste0('<section class="panel"><h2>Contrasts</h2>', lisa_scientific_report_table_html(contrast_details, 500L), '</section>') else "",
    '<section class="panel"><h2>Collections</h2>', lisa_scientific_report_table_html(collection_details, 500L), '</section>')
  lisa_scientific_report_write(file.path(output_dir, "pages", "details.html"), lisa_scientific_report_shell("Analysis details", study_label, "details", details_body, nav))

  lisa_scientific_report_write(file.path(output_dir, "README.md"), c(
    paste0("# ", title), "", "LISA: LLM-Inferred Semantic Annotation.", "",
    paste0("Updated: ", format(Sys.Date(), "%Y-%m-%d")), "Status: complete", "",
    "## Start here", "", "Open `index.html` in a modern browser.", "",
    "## Contents", "", "- `pages/`: Single-DE and contrast result pages.",
    "- `media/`: downloadable PNG, SVG and PDF figures.",
    "- `source_data/`: tables supporting the figures and gene-prioritization results.",
    "- `PROJECT.md`: analysis design, inputs, parameters and validation record.", "",
    "## Interpretation", "", "Any archived KEGG painted maps included in this report are imported expression overlays, not independent pathway-enrichment tests. The report does not retrieve or paint new KEGG maps.",
    "EnrichmentMap is not included unless explicitly requested."))
  lisa_scientific_report_write(file.path(output_dir, "PROJECT.md"), c(
    paste0("# ", title, " - analysis record"), "", paste0("Updated: ", format(Sys.Date(), "%Y-%m-%d")), "",
    "## Scope", "", "LISA (LLM-Inferred Semantic Annotation) Single-DE analyses and DE-vs-DE contrasts, organized by analysis and biological dictionary.", "",
    "## Analyses", "", paste0("- ", paste(unique(de_index$analysis_id), collapse = "\n- ")), "",
    "## Contrasts", "", if (nrow(contrast_index)) paste0("- ", paste(unique(contrast_index$contrast_id), collapse = "\n- ")) else "None", "",
    "## Collections", "", paste0("- ", paste(order_collections(unique(cards$collection)), collapse = "\n- ")), "",
    "## Scientific outputs", "", "Depending on the imported source, single-DE products include LISA summaries, member gene sets, within-category gene prioritization, recurrent genes, volcano overlays and supporting-gene heatmaps.",
    "Contrast products can include LISA category shifts, Contrast GeneCards, paired gene heatmaps and gene-category networks. Only products present in the source are displayed.", "",
    "## KEGG display", "", "The current pipeline supports figures of enriched KEGG gene sets, not painted native pathway maps. This report can also display archived painted maps when explicitly present in its imported source. Read their recorded scale, orientation and provenance; import does not recreate or revalidate their original painting analysis.", "",
    "## Inputs and provenance", "", "See `source_data/INDEX.tsv` for scientific tables and `metadata/source_provenance.tsv` for source paths and hashes.", "",
    "## Validation", "", "Local links, figure previews, source-data links, path limits and SHA-256 manifests are validated during report construction.",
    "For imported KEGG maps, this verifies the recorded files and hashes, not a new KGML retrieval or node-mapping analysis.", "",
    "## Limits", "", "The report does not infer causality. Contrast interpretation depends on the comparability of the two underlying DE designs."))

  # Final validation and immutable manifest.
  if (!lisa_scientific_report_root_layout_ok(output_dir)) stop("LISA-REPORT-SCIENCE-012 output root is not clean.", call. = FALSE)
  broken <- lisa_scientific_report_verify_links(output_dir)
  if (length(broken)) stop("LISA-REPORT-SCIENCE-013 broken local links: ", paste(head(broken, 20L), collapse = "; "), call. = FALSE)
  rel_files <- sort(list.files(output_dir, recursive = TRUE, full.names = FALSE))
  html <- rel_files[grepl("\\.html$", rel_files)]
  if (any(nchar(rel_files) > 120L) || any(nchar(basename(rel_files)) > 48L)) stop("LISA-REPORT-SCIENCE-014 Windows path budget exceeded.", call. = FALSE)
  if (any(file.info(file.path(output_dir, html))$size > 1024^2)) stop("LISA-REPORT-SCIENCE-015 HTML page budget exceeded.", call. = FALSE)
  qc$status[qc$check %in% c("clean_root", "offline_links", "path_budget", "html_budget")] <- "PASS"
  utils::write.table(qc, file.path(output_dir, "metadata", "qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rel_files <- sort(list.files(output_dir, recursive = TRUE, full.names = FALSE))
  rel_files <- rel_files[!rel_files %in% c("metadata/manifest.tsv", "metadata/SHA256SUMS")]
  manifest <- data.frame(path = rel_files,
    sha256 = vapply(file.path(output_dir, rel_files), lisa_scientific_report_sha256, character(1)),
    stringsAsFactors = FALSE)
  utils::write.table(manifest, file.path(output_dir, "metadata", "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sum_files <- c(rel_files, "metadata/manifest.tsv")
  sum_hashes <- c(manifest$sha256, lisa_scientific_report_sha256(file.path(output_dir, "metadata", "manifest.tsv")))
  writeLines(paste(sum_hashes, sum_files),
    file.path(output_dir, "metadata", "SHA256SUMS"), useBytes = TRUE)
  list(index = file.path(output_dir, "index.html"), files = length(list.files(output_dir, recursive = TRUE)),
    cards = nrow(cards), figures = nrow(figures), tables = nrow(tables), pages = length(html),
    max_path = max(nchar(list.files(output_dir, recursive = TRUE))),
    max_html = max(file.info(file.path(output_dir, html))$size), broken_links = 0L)
}
