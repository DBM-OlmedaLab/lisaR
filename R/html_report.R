# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

lisa_html_escape <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

build_lisa_html_index <- function(report_rows, output_file, title = "LISA report index") {
  if (is.character(report_rows) && length(report_rows) == 1) {
    report_rows <- read_lisa_tsv(report_rows)
  }
  lisa_require_columns(report_rows, c("label", "href"), "report_rows")

  items <- paste(
    sprintf(
      '<li><a href="%s">%s</a></li>',
      lisa_html_escape(report_rows$href),
      lisa_html_escape(report_rows$label)
    ),
    collapse = "\n"
  )

  html <- paste0(
    "<!doctype html>\n",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>", lisa_html_escape(title), "</title>",
    "<style>body{font-family:system-ui,sans-serif;max-width:980px;margin:40px auto;padding:0 24px;line-height:1.45}code{background:#f4f4f4;padding:2px 4px}</style>",
    "</head><body><main>",
    "<h1>", lisa_html_escape(title), "</h1>",
    "<ul>\n", items, "\n</ul>",
    "</main></body></html>\n"
  )

  lisa_guarded_write(output_file, function(target) writeLines(html, target, useBytes = TRUE))
  invisible(output_file)
}
