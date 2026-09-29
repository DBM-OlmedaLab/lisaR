# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Shared human-readable plot metadata for standalone LISA postprocessors.

lisa_plot_metadata <- function(project_dir, analysis_id) {
  path <- file.path(project_dir, "config", "de_index.tsv")
  fallback <- list(
    label = analysis_id,
    comparison = analysis_id,
    positive_direction = "",
    model_note = ""
  )
  if (!file.exists(path)) return(fallback)
  index <- utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  row <- index[index$analysis_id == analysis_id, , drop = FALSE]
  if (nrow(row) != 1L) return(fallback)
  value <- function(column, default = "") {
    if (!column %in% names(row)) return(default)
    x <- as.character(row[[column]][[1]])
    if (is.na(x) || !nzchar(trimws(x))) default else x
  }
  list(
    label = value("label", analysis_id),
    comparison = value("comparison", analysis_id),
    positive_direction = value("positive_direction", ""),
    model_note = value("model_note", "")
  )
}

lisa_plot_subtitle <- function(metadata, prefix = "", suffix = "") {
  values <- c(
    prefix,
    metadata$comparison,
    metadata$positive_direction,
    suffix
  )
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  lisa_plot_device_text(paste(values, collapse = "\n"))
}

lisa_plot_device_text <- function(x) {
  # Base PDF and some bitmap devices cannot encode typographic dash glyphs
  # with their default font. Keep the source metadata unchanged and normalize
  # only the human-readable text passed to graphics devices.
  gsub("[\u2013\u2014\u2212]", "-", as.character(x), perl = TRUE)
}
