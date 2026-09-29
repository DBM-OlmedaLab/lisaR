# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# report foundation static report.  This reader/copying layer deliberately never changes
# a scientific value: it validates and presents caller-supplied canonical TSVs.

lisa_report_foundation_required_files <- function(root) file.path(root, c("m/col.tsv", "m/de.tsv", "m/cx.tsv", "m/prov.tsv", "d/sd.tsv", "d/cx.tsv", "d/ge.tsv", "m/sha256.tsv"))

lisa_report_foundation_sha256 <- lisa_sha256_file

lisa_report_foundation_safe_relative <- function(x) {
  is.character(x) && length(x) == 1L && nzchar(x) && nchar(x) <= 80L &&
    grepl("^[A-Za-z0-9._/-]+$", x) && !grepl("(^/|(^|/)\\.\\.(/|$))", x)
}

lisa_report_foundation_fixture <- function(fixture_path) {
  root <- normalizePath(fixture_path, winslash = "/", mustWork = FALSE)
  required <- lisa_report_foundation_required_files(root); missing <- required[!file.exists(required)]
  if (length(missing)) stop("LISA-REPORT-FOUNDATION-001 external fixture is absent or incomplete: ", paste(missing, collapse = ", "), ". Supply an explicit external audit fixture.", call. = FALSE)
  manifest <- read_lisa_tsv(file.path(root, "m", "sha256.tsv")); lisa_require_columns(manifest, c("path", "sha256"), "report foundation checksum manifest")
  required_rel <- c("m/col.tsv", "m/de.tsv", "m/cx.tsv", "m/prov.tsv", "d/sd.tsv", "d/cx.tsv", "d/ge.tsv")
  if (!is.character(manifest$path) || anyNA(manifest$path) ||
      any(!nzchar(manifest$path)) || anyDuplicated(manifest$path) ||
      any(!required_rel %in% manifest$path) ||
      any(!vapply(as.character(manifest$path), lisa_report_foundation_safe_relative, logical(1))) ||
      !lisa_sha256_all_valid(manifest$sha256, nrow(manifest))) {
    stop("LISA-REPORT-FOUNDATION-002 malformed SHA-256 manifest.", call. = FALSE)
  }
  bad <- vapply(seq_len(nrow(manifest)), function(i) { p <- file.path(root, manifest$path[[i]]); !file.exists(p) || !identical(lisa_report_foundation_sha256(p), as.character(manifest$sha256[[i]])) }, logical(1))
  if (any(bad)) stop("LISA-REPORT-FOUNDATION-003 fixture file missing or mismatched while validating checksums: ", paste(manifest$path[bad], collapse = ", "), call. = FALSE)
  out <- list(root = root, collection = read_lisa_tsv(file.path(root, "m", "col.tsv")), de = read_lisa_tsv(file.path(root, "m", "de.tsv")), contrast_index = read_lisa_tsv(file.path(root, "m", "cx.tsv")), provenance = read_lisa_tsv(file.path(root, "m", "prov.tsv")), single = read_lisa_tsv(file.path(root, "d", "sd.tsv")), contrast = read_lisa_tsv(file.path(root, "d", "cx.tsv")), evidence = read_lisa_tsv(file.path(root, "d", "ge.tsv")), manifest = manifest)
  lisa_require_columns(out$single, c("category_id", "category_display_name", "n_genesets", "mean_NES", "same_direction_pct"), "canonical single.tsv")
  lisa_require_columns(out$contrast, c("category_id", "category_display_name", "mean_NES_A", "mean_NES_B", "delta_mean_NES"), "canonical contrast.tsv")
  lisa_require_columns(out$evidence, c("category_id", "symbol"), "canonical gene evidence")
  out
}

lisa_report_foundation_escape <- lisa_html_escape
lisa_report_foundation_rel <- function(from, to) { from <- strsplit(from, "/", fixed = TRUE)[[1]]; to <- strsplit(to, "/", fixed = TRUE)[[1]]; paste0(paste(rep("..", max(0L, length(from) - 1L)), collapse = "/"), if (length(from) > 1L) "/" else "", paste(to, collapse = "/")) }
lisa_report_foundation_page <- function(title, body, prefix = "") {
  asset <- function(x) paste0(prefix, x)
  paste0("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>", lisa_report_foundation_escape(title), "</title><link rel=\"stylesheet\" href=\"", asset("a/report.css"), "\"></head><body><header><img src=\"", asset("a/logo.svg"), "\" alt=\"LISA project logo\" height=\"45\"><h1>", lisa_report_foundation_escape(title), "</h1><nav><a href=\"", asset("index.html"), "\">Index</a><a href=\"", asset("sd-1.html"), "\">Single-DE</a><a href=\"", asset("cx-1.html"), "\">Contrast</a><a href=\"", asset("co.html"), "\">Collection</a><a href=\"", asset("qc.html"), "\">QC & provenance</a></nav><p><input id=\"q\" oninput=\"lisaSearch(this.value)\" placeholder=\"Search categories, genes, pages\"><span id=\"sr\"></span></p></header><main>", body, "</main><script src=\"", asset("a/search.js"), "\"></script></body></html>")
}

lisa_report_foundation_table <- function(x, columns = names(x)) {
  x <- x[, intersect(columns, names(x)), drop = FALSE]
  rows <- vapply(seq_len(nrow(x)), function(i) paste0("<tr data-s=\"", lisa_report_foundation_escape(paste(x[i, , drop = TRUE], collapse = " ")), "\">", paste0("<td>", lisa_report_foundation_escape(as.character(x[i, ])), "</td>", collapse = ""), "</tr>"), character(1))
  paste0("<table><thead><tr>", paste0("<th>", lisa_report_foundation_escape(names(x)), "</th>", collapse = ""), "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

lisa_report_foundation_bars <- function(x, value, label, title, metric) {
  value_num <- suppressWarnings(as.numeric(x[[value]])); finite <- value_num[is.finite(value_num)]; scale <- if (length(finite) && max(abs(finite)) > 0) 260 / max(abs(finite)) else 1
  cards <- vapply(seq_len(nrow(x)), function(i) paste0("<article class=\"card bar\" data-s=\"", lisa_report_foundation_escape(paste(x[i, , drop = TRUE], collapse = " ")), "\"><b>", lisa_report_foundation_escape(x[[label]][[i]]), "</b><span style=\"width:", max(2, min(260, abs(value_num[[i]]) * scale)), "px\"></span><output>", lisa_report_foundation_escape(x[[value]][[i]]), "</output></article>"), character(1))
  paste0("<section><h2>", lisa_report_foundation_escape(title), "</h2><p class=\"muted\">Metric: ", lisa_report_foundation_escape(metric), ". Every canonical row on this page is shown; pagination preserves all rows.</p>", paste(cards, collapse = ""), "</section>")
}

lisa_report_foundation_contrast_bars <- function(x) {
  a <- suppressWarnings(as.numeric(x$mean_NES_A)); b <- suppressWarnings(as.numeric(x$mean_NES_B)); fin <- c(a[is.finite(a)], b[is.finite(b)]); scale <- if (length(fin) && max(abs(fin)) > 0) 180 / max(abs(fin)) else 1
  cards <- vapply(seq_len(nrow(x)), function(i) paste0("<article class=\"card paired\" data-s=\"", lisa_report_foundation_escape(paste(x[i, , drop = TRUE], collapse = " ")), "\"><b>", lisa_report_foundation_escape(x$category_display_name[[i]]), "</b><div>A <span style=\"width:", max(2, min(180, abs(a[[i]]) * scale)), "px\"></span> ", lisa_report_foundation_escape(x$mean_NES_A[[i]]), "</div><div>B <span style=\"width:", max(2, min(180, abs(b[[i]]) * scale)), "px\"></span> ", lisa_report_foundation_escape(x$mean_NES_B[[i]]), "</div><output>Delta (A - B): ", lisa_report_foundation_escape(x$delta_mean_NES[[i]]), "</output></article>"), character(1))
  paste0("<section><h2>Paired A/B mean NES with delta</h2><p class=\"muted\">This is a paired A/B comparison from canonical contrast semantics, not a same-direction lollipop. Delta is canonical A - B. Pagination preserves all rows.</p>", paste(cards, collapse = ""), "</section>")
}

lisa_report_foundation_write_csv_bom <- function(x, path) { con <- file(path, open = "wb"); on.exit(close(con)); writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con); utils::write.csv(x, con, row.names = FALSE, na = "") }
lisa_report_foundation_copy <- function(from, to) { if (!file.copy(from, to, overwrite = TRUE)) stop("LISA-REPORT-FOUNDATION-010 unable to copy canonical source data.", call. = FALSE) }

lisa_report_foundation_local_links <- function(path) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  hits <- regmatches(html, gregexpr('(?:href|src)="[^"]+"', html, perl = TRUE))[[1]]
  if (!length(hits) || identical(hits, character(0))) return(character(0))
  links <- sub('^[^=]+="', "", hits)
  links <- sub('"$', "", links)
  links[!grepl("^(?:[a-z]+:|#|data:)", links, ignore.case = TRUE)]
}

lisa_report_foundation_verify_links <- function(output_dir) {
  pages <- list.files(output_dir, pattern = "\\.html$", full.names = TRUE)
  broken <- unlist(lapply(pages, function(page) {
    links <- lisa_report_foundation_local_links(page)
    links <- sub("[?#].*$", "", links)
    links <- links[nzchar(links)]
    target <- normalizePath(file.path(dirname(page), links), winslash = "/", mustWork = FALSE)
    missing <- links[!file.exists(target)]
    if (!length(missing)) character(0) else paste0(basename(page), " -> ", missing)
  }), use.names = FALSE)
  if (length(broken)) stop("LISA-REPORT-FOUNDATION-011 broken local report links: ", paste(broken, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Build an offline report foundation report from an explicit external fixture
#' @param fixture_path External audit-fixture directory, never a package fixture.
#' @param output_dir Destination directory for the static report.
#' @return Named paths to generated pages and primary downloads.
#' @keywords internal
build_lisa_foundation_report <- function(fixture_path, output_dir) {
  f <- lisa_report_foundation_fixture(fixture_path); output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE); invisible(vapply(file.path(output_dir, c("a", "d", "m")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
  logo <- system.file("report_assets", "LISA_logo_C_compact_icon_muted_red_S.svg", package = "lisaR")
  if (!nzchar(logo) || !file.exists(logo)) stop("LISA-REPORT-FOUNDATION-004 missing installed project logo asset. No logo was invented.", call. = FALSE)
  lisa_report_foundation_copy(logo, file.path(output_dir, "a/logo.svg"))
  lisa_report_foundation_copy(file.path(f$root, "d/sd.tsv"), file.path(output_dir, "d/s.tsv")); lisa_report_foundation_copy(file.path(f$root, "d/cx.tsv"), file.path(output_dir, "d/c.tsv")); lisa_report_foundation_copy(file.path(f$root, "d/ge.tsv"), file.path(output_dir, "d/g.tsv"))
  lisa_report_foundation_write_csv_bom(f$single, file.path(output_dir, "d/s.csv")); lisa_report_foundation_write_csv_bom(f$contrast, file.path(output_dir, "d/c.csv"))
  xlsx <- file.path(output_dir, "d/report.xlsx"); xlsx_state <- "not_available: optional package openxlsx is not installed"
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    sheets <- list(single_de = f$single, contrast = f$contrast, gene_evidence = f$evidence)
    for (sheet_name in names(sheets)) {
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, sheets[[sheet_name]])
    }
    openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)
    xlsx_state <- "available: primary workbook d/report.xlsx"
  }
  css <- "body{font-family:system-ui,sans-serif;max-width:1160px;margin:2rem auto;padding:0 1rem;color:#1d2630}header{border-bottom:4px solid #9e2936;padding-bottom:1rem}nav a{margin-right:1rem}table{border-collapse:collapse;width:100%;font-size:.86rem}th,td{border-bottom:1px solid #ddd;padding:.35rem;text-align:left}.card{padding:.65rem;margin:.4rem 0;background:#f7f7f7;border-radius:.4rem}.bar span,.paired span{display:inline-block;height:9px;background:#9e2936;margin:.2rem}.paired div{margin:.15rem}.muted{color:#566;font-size:.9rem}input{padding:.45rem;width:20rem;max-width:100%}"
  writeLines(css, file.path(output_dir, "a/report.css"), useBytes = TRUE)
  search <- list(categories = unique(as.character(f$single$category_display_name)), genes = unique(as.character(f$evidence$symbol)), pages = c("index.html", "sd-1.html", "cx-1.html", "co.html", "qc.html"))
  writeLines(jsonlite::toJSON(search, auto_unbox = TRUE), file.path(output_dir, "a/search.json"), useBytes = TRUE)
  # Embed the compact search index in JavaScript so search works under file://
  # on Windows without a local web server or browser fetch permissions.
  js_index <- jsonlite::toJSON(search, auto_unbox = TRUE)
  js <- paste0("const LISA_SEARCH_INDEX=", js_index, ";function lisaSearch(q){q=(q||'').toLowerCase();var local=document.querySelectorAll('[data-s]');local.forEach(x=>x.style.display=x.dataset.s.toLowerCase().includes(q)?'':'none');var r=document.getElementById('sr');if(!q){r.textContent='';return}var x=LISA_SEARCH_INDEX;var n=[...x.categories,...x.genes,...x.pages].filter(z=>z.toLowerCase().includes(q)).slice(0,12);r.textContent=n.length?' Matches: '+n.join(', '):' No index match'}")
  writeLines(js, file.path(output_dir, "a/search.js"), useBytes = TRUE)
  write_page <- function(name, title, body) writeLines(lisa_report_foundation_page(title, body), file.path(output_dir, name), useBytes = TRUE)
  pages_for <- function(x, size) {
    if (!nrow(x)) return(list(integer(0)))
    split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / size))
  }
  # Single-DE renders three cards per canonical row. Sixty rows therefore
  # remain below the approved 200-card budget (3 * 60 = 180).
  sp <- pages_for(f$single, 60L); cp <- pages_for(f$contrast, 100L)
  pager <- function(stem, i, n) {
    links <- c(if (i > 1L) sprintf('<a href="%s-%d.html">Previous</a>', stem, i - 1L),
               if (i < n) sprintf('<a href="%s-%d.html">Next</a>', stem, i + 1L))
    paste0("<p>", paste(links, collapse = " | "), "</p>")
  }
  single_cols <- c("category_id", "category_display_name", "n_genesets", "mean_NES", "min_padj", "same_direction_pct", "mean_NES_direction")
  for (i in seq_along(sp)) { z <- f$single[sp[[i]], , drop = FALSE]; nav <- paste0("<p>Single-DE page ", i, " of ", length(sp), "; rows ", min(sp[[i]]), "-", max(sp[[i]]), " of ", nrow(f$single), ".</p>", pager("sd", i, length(sp))); body <- paste0(nav, lisa_report_foundation_bars(z, "n_genesets", "category_display_name", "Category LISA signal", "n_genesets"), lisa_report_foundation_bars(z, "mean_NES", "category_display_name", "Mean NES", "mean_NES"), lisa_report_foundation_bars(z, "same_direction_pct", "category_display_name", "Same-direction percentage", "same_direction_pct"), lisa_report_foundation_table(z, single_cols), pager("sd", i, length(sp))); write_page(paste0("sd-", i, ".html"), "Single-DE canonical category summary", body) }
  cxcols <- c("category_id", "category_display_name", "mean_NES_A", "mean_NES_B", "delta_mean_NES", "direction_A", "direction_B", "direction_class", "min_padj_A", "min_padj_B")
  for (i in seq_along(cp)) { z <- f$contrast[cp[[i]], , drop = FALSE]; body <- paste0("<p>Contrast page ", i, " of ", length(cp), "; rows ", min(cp[[i]]), "-", max(cp[[i]]), " of ", nrow(f$contrast), ".</p>", pager("cx", i, length(cp)), lisa_report_foundation_contrast_bars(z), lisa_report_foundation_table(z, cxcols), pager("cx", i, length(cp))); write_page(paste0("cx-", i, ".html"), "Contrast canonical category summary", body) }
  gps <- lisa_gps_development_contract(); collection <- paste0("<h2>Collection</h2>", lisa_report_foundation_table(f$collection), "<h2>Gene evidence</h2>", lisa_report_foundation_table(utils::head(f$evidence, 100L), intersect(c("category_id", "symbol", "log2FC_A", "padj_A", "gene_contribution_score_A", "gene_set_id"), names(f$evidence))), "<p>Gene evidence has ", nrow(f$evidence), " canonical rows. The table preview does not alter the copied source data.</p><h2>lisa-gps</h2><p>Name: lisa-gps. Version: ", lisa_report_foundation_escape(gps$development_contract_version), ". Formula: ", lisa_report_foundation_escape(gps$transformations$score), ". Components: effect strength; DE confidence; independent category breadth. No historical score table is shown in this report.</p>")
  write_page("co.html", "Collection and gene evidence", collection)
  downloads <- paste0("<div class=\"card\"><h2>Downloads</h2><p>Primary downloads: ", if (file.exists(xlsx)) "<a href=\"d/report.xlsx\">XLSX workbook</a>" else lisa_report_foundation_escape(xlsx_state), ". UTF-8 BOM CSV exports: <a href=\"d/s.csv\">single-DE CSV</a>, <a href=\"d/c.csv\">contrast CSV</a>.</p><p>Canonical source copies: <a href=\"d/s.tsv\">single canonical data</a>, <a href=\"d/c.tsv\">contrast canonical data</a>, <a href=\"d/g.tsv\">gene evidence canonical data</a>.</p></div>")
  write_page("index.html", "LISA report foundation offline report", paste0("<p>Canonical values are copied and displayed without recalculation. Single-DE pages: ", paste(sprintf("<a href=\"sd-%d.html\">%d</a>", seq_along(sp), seq_along(sp)), collapse = " "), ". Contrast pages: ", paste(sprintf("<a href=\"cx-%d.html\">%d</a>", seq_along(cp), seq_along(cp)), collapse = " "), ".</p>", downloads))
  qc <- paste0("<h2>QC status: PASS</h2><p>Fixture SHA-256 verification: PASS. Report SHA-256 verification: PASS. No scientific computation was performed.</p><h2>Source lineage</h2>", lisa_report_foundation_table(f$provenance), "<h2>Versions</h2><p>lisaR ", utils::packageVersion("lisaR"), "; report contract report foundation; ", lisa_report_foundation_escape(gps$development_contract_version), ".</p><h2>Assets</h2><p>Project logo: PASS (package LISA SVG). Institution logo: PASS_WITH_WARNINGS - no institutional-logo asset exists in this package; none was invented.</p><h2>Workbook</h2><p>", lisa_report_foundation_escape(xlsx_state), ".</p><h2>Checksum design</h2><p>m/SHA256SUMS includes m/manifest.tsv and every other report artifact, but cannot include itself (a cryptographic self-reference); this is the documented non-recursive exception.</p>")
  write_page("qc.html", "QC and provenance", qc)
  manifest_rows <- sort(list.files(output_dir, recursive = TRUE, full.names = FALSE)); manifest_rows <- manifest_rows[!grepl("^m/(manifest.tsv|SHA256SUMS)$", manifest_rows)]
  manifest <- data.frame(path = manifest_rows, sha256 = vapply(file.path(output_dir, manifest_rows), lisa_report_foundation_sha256, character(1)), stringsAsFactors = FALSE); write.table(manifest, file.path(output_dir, "m/manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sums <- sort(c(manifest_rows, "m/manifest.tsv")); sum_hashes <- c(stats::setNames(manifest$sha256, manifest$path), "m/manifest.tsv" = lisa_report_foundation_sha256(file.path(output_dir, "m/manifest.tsv")))[sums]; writeLines(paste(sum_hashes, sums), file.path(output_dir, "m/SHA256SUMS"), useBytes = TRUE)
  pages <- file.path(output_dir, c("index.html", "co.html", "qc.html", paste0("sd-", seq_along(sp), ".html"), paste0("cx-", seq_along(cp), ".html")))
  if (any(file.info(pages)$size > 1024^2)) stop("LISA-REPORT-FOUNDATION-005 generated HTML exceeds 1 MB.", call. = FALSE)
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE); rel <- substring(files, nchar(output_dir) + 2L); if (any(nchar(rel) > 120L) || any(nchar(basename(rel)) > 48L)) stop("LISA-REPORT-FOUNDATION-007 report path budget exceeded.", call. = FALSE)
  lisa_report_foundation_verify_links(output_dir)
  c(index = file.path(output_dir, "index.html"), single_de = file.path(output_dir, "sd-1.html"), contrast = file.path(output_dir, "cx-1.html"), collection = file.path(output_dir, "co.html"), qc = file.path(output_dir, "qc.html"), xlsx = if (file.exists(xlsx)) xlsx else NA_character_, single_csv = file.path(output_dir, "d/s.csv"), contrast_csv = file.path(output_dir, "d/c.csv"))
}

lisa_report_foundation_build_report <- function(fixture_dir, output_dir) build_lisa_foundation_report(fixture_dir, output_dir)
