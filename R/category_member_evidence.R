# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Individual-category enrichment at the end of an evidence sheet. This module
# renders existing GSEA results only; matrix display caps do not select sets.
lisa_category_member_evidence_source <- function(evidence, category_id) {
  if (!inherits(evidence, "lisa_category_evidence")) stop("Expected lisa_category_evidence object.", call. = FALSE)
  ii <- match(as.character(category_id), as.character(evidence$categories$category_id))
  if (length(ii) != 1L || is.na(ii)) stop("Unknown exact category ID.", call. = FALSE)
  if (identical(as.character(category_id), "OTHER_UNCLASSIFIED")) stop("Member figures require a classified category.", call. = FALSE)
  category <- evidence$categories[ii, , drop = FALSE]
  sets <- evidence$sets[as.character(evidence$sets$category_id) == category_id, , drop = FALSE]
  if (anyDuplicated(sets$pathway)) stop("Duplicate gene sets within category evidence.", call. = FALSE)
  flag <- function(x) as.character(x) %in% c("TRUE", "true")
  cutoff <- lisa_evidence_num(evidence$metadata$gsea_padj_cutoff)
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff < 0 || cutoff > 1) stop("Invalid recorded GSEA cutoff.", call. = FALSE)
  nes <- lisa_evidence_num(sets$NES); fdr <- lisa_evidence_num(sets$gsea_fdr)
  valid <- is.finite(nes) & is.finite(fdr) & fdr >= 0 & fdr <= 1
  evaluable <- flag(sets$evaluable); significant <- flag(sets$significant)
  if (any(evaluable & !valid) || any(significant != (evaluable & valid & fdr <= cutoff)))
    stop("Member evidence significance conflicts with its recorded values/cutoff.", call. = FALSE)
  selected <- which(significant)
  selected <- selected[lisa_evidence_order(fdr[selected], -abs(nes[selected]), as.character(sets$pathway[selected]))]
  # Shared symmetric axis across this evidence scope, not a category-specific
  # rescaling that would make equally sized marks represent different NES.
  all_nes <- lisa_evidence_num(evidence$sets$NES)
  all_sig <- flag(evidence$sets$significant) & flag(evidence$sets$evaluable) & is.finite(all_nes)
  limit <- max(c(1, abs(all_nes[all_sig]))) * 1.08
  columns <- c("row_type", "figure_type", "analysis_id", "collection", "tier", "positive_contrast",
    "category_id", "category_display_name", "gsea_padj_cutoff", "n_mapped_sets", "n_evaluable_sets",
    "n_significant_sets", "selection_rule", "set_sort_label", "plot_x_min", "plot_x_max",
    "pathway", "NES", "gsea_fdr", "evaluable", "significant", "selected_for_plot", "plot_order", "color")
  base <- as.data.frame(stats::setNames(rep(list(""), length(columns)), columns), stringsAsFactors = FALSE)
  base$row_type <- "metadata"; base$figure_type <- "category_member_evidence"
  for (key in intersect(names(evidence$metadata), columns)) base[[key]] <- evidence$metadata[[key]]
  base$category_id <- as.character(category_id)
  base$category_display_name <- as.character(category$category_display_name)
  base$n_mapped_sets <- nrow(sets); base$n_evaluable_sets <- sum(evaluable)
  base$n_significant_sets <- length(selected)
  base$selection_rule <- "All significant evaluable member sets at the recorded GSEA FDR cutoff; no matrix display cap."
  base$set_sort_label <- "GSEA FDR ascending, absolute NES descending, exact pathway ID ascending"
  base$plot_x_min <- -limit; base$plot_x_max <- limit
  rows <- base
  if (nrow(sets)) {
    members <- base[rep(1L, nrow(sets)), , drop = FALSE]
    members$row_type <- "set"; members$pathway <- as.character(sets$pathway)
    members$NES <- nes; members$gsea_fdr <- fdr
    members$evaluable <- evaluable; members$significant <- significant
    members$selected_for_plot <- significant
    members$plot_order <- match(seq_len(nrow(sets)), selected)
    members$color <- ifelse(!is.finite(nes), "#9ca3af", ifelse(nes > 0, "#b2182b", ifelse(nes < 0, "#2166ac", "#64748b")))
    # Put selected rows first but retain excluded candidates and their states.
    members <- members[lisa_evidence_order(members$plot_order, members$pathway), , drop = FALSE]
    # Convert before rbind: a blank metadata NES cell otherwise coerces the
    # numeric member column through default (15-digit) as.character first.
    base[] <- lapply(base, lisa_evidence_text)
    members[] <- lapply(members, lisa_evidence_text)
    rows <- rbind(base, members)
  }
  rows[] <- lapply(rows, lisa_evidence_text)
  rownames(rows) <- NULL
  rows
}

# Self-contained renderer copied verbatim into the standalone base-R recipe.
# Position is original NES; explicit text gives original GSEA FDR. Sign colour
# has no implication of a new test or a combined category significance.
lisa_category_member_evidence_draw <- function(source) {
  num <- function(x) suppressWarnings(as.numeric(x))
  fmt <- function(x) ifelse(is.finite(num(x)), formatC(num(x), format = "g", digits = 3), "NA")
  device_text <- function(x) gsub("[\u2013\u2014\u2212]", "-", as.character(x), perl = TRUE)
  shorten <- function(x, n = 67L) ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 3L), "..."), x)
  meta <- source[source$row_type == "metadata", , drop = FALSE]
  if (nrow(meta) != 1L) stop("Expected exactly one figure metadata row.")
  sets <- source[source$row_type == "set" & source$selected_for_plot %in% c("TRUE", "true"), , drop = FALSE]
  sets <- sets[order(num(sets$plot_order)), , drop = FALSE]
  n <- nrow(sets); limit <- max(abs(num(c(meta$plot_x_min, meta$plot_x_max))))
  if (!is.finite(limit) || limit <= 0) stop("Invalid NES display axis.")
  old <- graphics::par(mar = c(6, 27, 7, 1), xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(-limit, limit * 1.75), ylim = c(.5, max(1, n) + 1.6), xaxs = "i")
  graphics::title(main = device_text(paste0(meta$category_display_name, " [", meta$category_id, "]")), line = 5.1, cex.main = .95)
  graphics::mtext(device_text(paste(meta$analysis_id, meta$collection, paste0("tier: ", meta$tier), sep = " | ")), side = 3, line = 3.9, cex = .74)
  graphics::mtext(device_text(paste0("Positive direction: ", meta$positive_contrast)), side = 3, line = 2.8, cex = .69)
  graphics::mtext(paste0("All ", n, " significant member sets; GSEA FDR <= ", meta$gsea_padj_cutoff,
    " | mapped / evaluable: ", meta$n_mapped_sets, " / ", meta$n_evaluable_sets), side = 3, line = 1.6, cex = .7)
  if (!n) {
    graphics::text(-limit, 1, "No significant evaluable member sets at this cutoff.\nAll mapped set states remain in the source table.", adj = 0, cex = .8)
  } else {
    y <- n + 1L - seq_len(n)
    for (i in seq_len(n)) {
      graphics::rect(-limit, y[[i]] - .44, limit * 1.72, y[[i]] + .44,
        border = NA, col = if (i %% 2L) "#f1f5f9" else "white")
      graphics::segments(0, y[[i]], num(sets$NES[[i]]), y[[i]], col = sets$color[[i]], lwd = 1.1)
      graphics::points(num(sets$NES[[i]]), y[[i]], pch = 16, cex = .85, col = sets$color[[i]])
      graphics::text(-limit * 1.025, y[[i]], device_text(shorten(sets$pathway[[i]])), adj = 1, cex = .62)
      graphics::text(limit * 1.12, y[[i]], fmt(sets$NES[[i]]), adj = 0, cex = .64, col = sets$color[[i]])
      graphics::text(limit * 1.43, y[[i]], fmt(sets$gsea_fdr[[i]]), adj = 0, cex = .61)
    }
    graphics::segments(0, .5, 0, n + .55, col = "#8998a4", lty = 2, lwd = .7)
    graphics::text(limit * c(1.12, 1.43), n + 1, c("NES", "FDR GSEA"), adj = 0, font = 2, cex = .66)
  }
  ticks <- pretty(c(-limit, limit), n = 5L); ticks <- ticks[ticks >= -limit & ticks <= limit]
  graphics::axis(1, at = ticks, labels = fmt(ticks), cex.axis = .7)
  graphics::mtext("Gene-set NES", side = 1, line = 2, at = 0, cex = .8)
  graphics::mtext("Red = positive NES; blue = negative; exact zero is grey. Points have equal size; GSEA FDR is listed separately.", side = 1, line = 3.4, cex = .63)
  graphics::mtext("Order: GSEA FDR, then |NES|, then exact ID. Complete identifiers and values are in the source TSV.", side = 1, line = 4.5, cex = .63)
  invisible(NULL)
}

render_lisa_category_member_evidence <- function(evidence, output_dir, formats = "png") {
  if (!inherits(evidence, "lisa_category_evidence")) stop("Expected lisa_category_evidence object.", call. = FALSE)
  formats <- unique(tolower(formats[nzchar(formats)]))
  if (any(!formats %in% c("png", "pdf", "svg"))) stop("Member figure formats must be png, pdf and/or svg.", call. = FALSE)
  ids <- as.character(evidence$categories$category_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Category IDs must be nonempty and unique.", call. = FALSE)
  if (any(ids == "OTHER_UNCLASSIFIED")) stop("Member figures require classified categories.", call. = FALSE)
  lisa_guarded_dir_create(output_dir)
  figure_dir <- file.path(output_dir, "figures"); lisa_guarded_dir_create(figure_dir)
  recipe_rel <- "figures/reproduce_category_member_evidence.R"
  recipe_path <- file.path(output_dir, recipe_rel)
  recipe <- c("#!/usr/bin/env Rscript", "# Reproduce existing per-category GSEA member evidence without lisaR or GSEA.",
    "# Usage: Rscript --vanilla reproduce_category_member_evidence.R SOURCE.tsv OUTPUT.png",
    paste0("lisa_category_member_evidence_draw <- ", paste(deparse(lisa_category_member_evidence_draw), collapse = "\n")),
    "args <- commandArgs(TRUE)", "if (length(args) != 2L) stop('Expected source TSV and output PNG/PDF/SVG paths.')",
    "source <- utils::read.delim(args[[1]], sep='\\t', quote='', comment.char='', check.names=FALSE, colClasses='character', na.strings=NULL)",
    "n <- sum(source$row_type == 'set' & source$selected_for_plot %in% c('TRUE','true'))",
    "width <- 15; height <- max(5.5, 2.8 + .18 * n)",
    "format <- tolower(tools::file_ext(args[[2]])); if (!format %in% c('png','pdf','svg')) stop('Use PNG, PDF or SVG.')",
    "if (format == 'png') grDevices::png(args[[2]], width=width,height=height,units='in',res=140,bg='white')",
    "if (format == 'pdf') grDevices::pdf(args[[2]], width=width,height=height,useDingbats=FALSE)",
    "if (format == 'svg') grDevices::svg(args[[2]], width=width,height=height,bg='white')",
    "tryCatch(lisa_category_member_evidence_draw(source),finally=grDevices::dev.off())")
  lisa_guarded_write(recipe_path, function(target) writeLines(recipe, target, useBytes = TRUE))
  index <- data.frame(category_id = ids, stem = sprintf("%03d_category_members", seq_along(ids)),
    source_tsv = rep("", length(ids)), recipe_r = rep(recipe_rel, length(ids)), n_significant_sets = integer(length(ids)),
    source_sha256 = rep("", length(ids)), recipe_sha256 = rep(digest::digest(file = recipe_path, algo = "sha256"), length(ids)), stringsAsFactors = FALSE)
  for (i in seq_along(ids)) {
    source <- lisa_category_member_evidence_source(evidence, ids[[i]])
    stem <- file.path(figure_dir, index$stem[[i]])
    source_path <- paste0(stem, "_source.tsv")
    write_lisa_tsv(source, source_path)
    n <- sum(source$row_type == "set" & source$selected_for_plot %in% c("TRUE", "true"))
    for (format in formats) lisa_guarded_write(paste0(stem, ".", format), function(target) {
      width <- 15; height <- max(5.5, 2.8 + .18 * n)
      if (format == "png") grDevices::png(target, width = width, height = height, units = "in", res = 140, bg = "white")
      if (format == "pdf") grDevices::pdf(target, width = width, height = height, useDingbats = FALSE)
      if (format == "svg") grDevices::svg(target, width = width, height = height, bg = "white")
      tryCatch(lisa_category_member_evidence_draw(source), finally = grDevices::dev.off())
    })
    index$source_tsv[[i]] <- paste0("figures/", index$stem[[i]], "_source.tsv")
    index$n_significant_sets[[i]] <- n
    index$source_sha256[[i]] <- digest::digest(file = source_path, algo = "sha256")
  }
  write_lisa_tsv(index, file.path(figure_dir, "category_members_index.tsv"))
  index
}
