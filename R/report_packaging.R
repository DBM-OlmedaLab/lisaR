# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Report packaging builder. It only packages validated,
# caller-supplied archived tables; it never runs or changes LISA analysis.

lisa_report_package_sha256 <- lisa_report_foundation_sha256
lisa_report_package_safe_relative <- lisa_report_foundation_safe_relative
lisa_report_package_escape <- lisa_html_escape

lisa_report_package_read_fixture <- function(fixture_path) {
  root <- normalizePath(fixture_path, winslash = "/", mustWork = FALSE)
  need <- file.path(root, c("m/registry.tsv", "m/sha256.tsv", "m/de.tsv", "m/cx.tsv", "m/col.tsv", "m/prov.tsv"))
  if (any(!file.exists(need))) stop("LISA-REPORT-PACKAGE-001 fixture is incomplete. Supply a checksummed external report packaging fixture.", call. = FALSE)
  registry <- read_lisa_tsv(file.path(root, "m/registry.tsv"))
  lisa_require_columns(registry, c("table_id", "path", "role", "sha256"), "report packaging table registry")
  if (!"artifact_type" %in% names(registry)) registry$artifact_type <- "table"
  registry$artifact_type <- tolower(as.character(registry$artifact_type))
  valid_id <- (registry$artifact_type == "table" & grepl("^T[0-9]{4,}$", registry$table_id)) |
    (registry$artifact_type == "figure" & grepl("^F[0-9]{4,}$", registry$table_id))
  if (!nrow(registry) || anyDuplicated(registry$table_id) || anyDuplicated(registry$path) ||
      !is.character(registry$path) || anyNA(registry$path) || any(!nzchar(registry$path)) ||
      any(!registry$artifact_type %in% c("table", "figure")) || any(!valid_id) ||
      any(!vapply(as.character(registry$path), lisa_report_package_safe_relative, logical(1))) ||
      !lisa_sha256_all_valid(registry$sha256, nrow(registry)))
    stop("LISA-REPORT-PACKAGE-002 malformed short registry.", call. = FALSE)
  sums <- read_lisa_tsv(file.path(root, "m/sha256.tsv")); lisa_require_columns(sums, c("path", "sha256"), "report packaging fixture manifest")
  if (!is.character(sums$path) || anyNA(sums$path) || any(!nzchar(sums$path)) ||
      anyDuplicated(sums$path) ||
      any(!vapply(as.character(sums$path), lisa_report_package_safe_relative, logical(1))) ||
      !lisa_sha256_all_valid(sums$sha256, nrow(sums))) {
    stop("LISA-REPORT-PACKAGE-003 fixture checksum manifest is malformed.", call. = FALSE)
  }
  declared <- sort(as.character(sums$path)); actual <- sort(list.files(root, recursive = TRUE, full.names = FALSE))
  actual <- actual[actual != "m/sha256.tsv"]
  if (!identical(declared, actual)) stop("LISA-REPORT-PACKAGE-003 fixture has undeclared or missing files.", call. = FALSE)
  bad <- vapply(seq_len(nrow(sums)), function(i) { p <- file.path(root, sums$path[[i]]); !file.exists(p) || !identical(lisa_report_package_sha256(p), as.character(sums$sha256[[i]])) }, logical(1))
  if (any(bad)) stop("LISA-REPORT-PACKAGE-004 fixture checksum mismatch: ", paste(sums$path[bad], collapse = ", "), call. = FALSE)
  if (any(!file.exists(file.path(root, registry$path)))) stop("LISA-REPORT-PACKAGE-005 registered canonical table is absent.", call. = FALSE)
  rb <- vapply(seq_len(nrow(registry)), function(i) identical(lisa_report_package_sha256(file.path(root, registry$path[[i]])), as.character(registry$sha256[[i]])), logical(1))
  if (any(!rb)) stop("LISA-REPORT-PACKAGE-006 registry checksum mismatch.", call. = FALSE)
  list(root = root, registry = registry, de = read_lisa_tsv(file.path(root, "m/de.tsv")), contrast = read_lisa_tsv(file.path(root, "m/cx.tsv")), collection = read_lisa_tsv(file.path(root, "m/col.tsv")), provenance = read_lisa_tsv(file.path(root, "m/prov.tsv")))
}

lisa_report_package_csv_bom <- function(x, path) lisa_report_foundation_write_csv_bom(x, path)
lisa_report_package_copy <- function(from, to) { if (!file.copy(from, to, overwrite = TRUE)) stop("LISA-REPORT-PACKAGE-010 unable to copy validated canonical table.", call. = FALSE) }
lisa_report_package_table <- function(x, n = 50L) lisa_report_foundation_table(utils::head(x, n))
lisa_report_package_has_alias <- function(x, value) {
  value <- as.character(value)
  unname(vapply(as.character(x), function(z) value %in% strsplit(z, "|", fixed = TRUE)[[1L]], logical(1)))
}

lisa_report_package_option_names <- function() {
  c("gene_cards", "volcano_overlays", "heatmaps", "kegg_maps",
    "enrichment_maps", "recurrent_screens", "other_visuals")
}

lisa_report_package_report_options <- function(publication_profile, report_options = NULL) {
  option_names <- lisa_report_package_option_names()
  if (identical(publication_profile, "minimal")) {
    out <- stats::setNames(as.list(rep(FALSE, length(option_names))), option_names)
  } else if (identical(publication_profile, "annotated")) {
    out <- stats::setNames(as.list(rep(TRUE, length(option_names))), option_names)
  } else {
    if (is.null(report_options) || !is.list(report_options)) {
      stop("LISA-REPORT-PACKAGE-014 custom publication profile requires named report_options.", call. = FALSE)
    }
    unknown <- setdiff(names(report_options), option_names)
    if (length(unknown)) stop("LISA-REPORT-PACKAGE-015 unknown report option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    out <- stats::setNames(as.list(rep(FALSE, length(option_names))), option_names)
    out[names(report_options)] <- lapply(report_options, isTRUE)
  }
  out
}

lisa_report_package_visual_product <- function(role) {
  role <- tolower(as.character(role))
  if (grepl("gene.?card", role)) return("gene_cards")
  if (grepl("volcano", role)) return("volcano_overlays")
  if (grepl("heatmap|leading_edge_ge", role)) return("heatmaps")
  if (grepl("kegg|painted.?pathway", role)) return("kegg_maps")
  if (grepl("enrichment.?map", role)) return("enrichment_maps")
  if (grepl("recurrent", role)) return("recurrent_screens")
  "other_visuals"
}

lisa_report_package_select_registry <- function(registry, options) {
  registry$product <- ifelse(
    registry$artifact_type == "table", "canonical_tables",
    vapply(registry$role, lisa_report_package_visual_product, character(1)))
  requested <- vapply(registry$product, function(product) {
    identical(product, "canonical_tables") || isTRUE(options[[product]])
  }, logical(1))
  registry[requested, , drop = FALSE]
}
lisa_report_package_page <- function(title, body) {
  paste0("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>", lisa_report_package_escape(title), "</title><link rel=\"stylesheet\" href=\"a/r.css\"></head><body><header><img src=\"a/logo.svg\" alt=\"LISA project logo\" height=\"44\"><h1>", lisa_report_package_escape(title), "</h1><nav><a href=\"index.html\">Home</a> <a href=\"dl.html\">Downloads</a> <a href=\"qc.html\">QC</a> <a href=\"methods.html\">Methods</a> <a href=\"prov.html\">Provenance</a></nav><input id=\"q\" placeholder=\"Search report\" oninput=\"lisaSearch(this.value)\"><span id=\"sr\"></span></header><main>", body, "</main><script src=\"a/s.js\"></script></body></html>")
}

lisa_report_package_status <- function(f, options, selected) {
  visual_products <- lisa_report_package_option_names()
  available <- vapply(visual_products, function(product) {
    any(f$registry$artifact_type == "figure" &
        vapply(f$registry$role, lisa_report_package_visual_product, character(1)) == product)
  }, logical(1))
  included <- vapply(visual_products, function(product) {
    any(selected$artifact_type == "figure" & selected$product == product)
  }, logical(1))
  requested <- vapply(visual_products, function(product) isTRUE(options[[product]]), logical(1))
  visual_state <- ifelse(!requested, "disabled", ifelse(!available, "not_available", ifelse(included, "PASS", "FAIL")))
  core <- data.frame(
    product = c("report", "canonical_tables", "xlsx", "csv_zip"),
    requested = TRUE, available = TRUE, included = TRUE, state = "PASS",
    reason = c("offline multipage bundle", "checksummed single-copy registry", "consolidated workbook", "UTF-8-BOM CSV archive"),
    stringsAsFactors = FALSE)
  visuals <- data.frame(
    product = visual_products, requested = requested, available = available,
    included = included, state = visual_state,
    reason = ifelse(!requested, "disabled by effective report profile",
      ifelse(!available, "requested but absent from the canonical archived source", "archived visual artifacts included")),
    stringsAsFactors = FALSE)
  rbind(core, visuals)
}

lisa_report_package_write_xlsx <- function(registry, data_dir, output) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("LISA-REPORT-PACKAGE-011 openxlsx is required for the report packaging primary workbook.", call. = FALSE)
  # Excel hard limits are not sufficient operational limits. Keep this
  # workbook memory-bounded and usable; every excluded table remains complete
  # in TSV/CSV and is explicitly recorded in README.
  selected <- seq_len(min(nrow(registry), 249L)); overflow <- registry$table_id[-selected]
  wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "README")
  dict <- registry; dict$xlsx_state <- ifelse(seq_len(nrow(dict)) %in% selected, "pending", "CSV_only_sheet_limit")
  openxlsx::writeData(wb, "README", dict); openxlsx::freezePane(wb, "README", firstRow = TRUE); openxlsx::addFilter(wb, "README", rows = 1, cols = seq_len(ncol(dict)))
  cumulative_cells <- 0
  for (i in selected) {
    x <- read_lisa_tsv(file.path(data_dir, paste0(registry$table_id[[i]], ".tsv")))
    if (!ncol(x)) {
      dict$xlsx_state[[i]] <- "CSV_only_empty_table"; next
    }
    cells <- as.double(nrow(x)) * as.double(ncol(x))
    if (nrow(x) > 100000L || cells > 1000000 || cumulative_cells + cells > 2000000 ||
        any(nchar(as.matrix(x), type = "chars") > 32767L, na.rm = TRUE)) {
      dict$xlsx_state[[i]] <- "CSV_only_excel_or_workbook_budget"; next
    }
    sn <- registry$table_id[[i]]; openxlsx::addWorksheet(wb, sn); openxlsx::writeData(wb, sn, x, keepNA = TRUE); openxlsx::freezePane(wb, sn, firstRow = TRUE); openxlsx::addFilter(wb, sn, rows = 1, cols = seq_len(ncol(x)))
    dict$xlsx_state[[i]] <- "included"
    cumulative_cells <- cumulative_cells + cells
  }
  openxlsx::writeData(wb, "README", dict, startRow = 1, colNames = TRUE)
  openxlsx::saveWorkbook(wb, output, overwrite = TRUE)
  list(included = sum(dict$xlsx_state == "included"), csv_only = sum(dict$xlsx_state != "included"), overflow = overflow)
}

#' Build the report packaging offline full-report bundle
#'
#' @param fixture_path External checksummed fixture created from a canonical run.
#' @param output_dir Destination bundle directory.
#' @param publication_profile One of `minimal`, `annotated`, or `custom`.
#'   `annotated` is the effective all-true profile for archived visual products.
#' @param report_options Named logical list used only by the `custom` profile.
#' @return Paths and deterministic report inventory.
#' @keywords internal
build_lisa_report_package <- function(fixture_path, output_dir, publication_profile = c("minimal", "annotated", "custom"), report_options = NULL) {
  publication_profile <- match.arg(publication_profile); f <- lisa_report_package_read_fixture(fixture_path)
  report_options <- lisa_report_package_report_options(publication_profile, report_options)
  selected <- lisa_report_package_select_registry(f$registry, report_options)
  tables <- selected[selected$artifact_type == "table", , drop = FALSE]
  figures <- selected[selected$artifact_type == "figure", , drop = FALSE]
  output_dir <- path.expand(output_dir)
  if (!lisa_is_absolute_path(output_dir)) {
    output_dir <- file.path(getwd(), output_dir)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) stop("LISA-REPORT-PACKAGE-012 destination must be empty; refusing to merge report artifacts.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(
    output_dir, winslash = "/", mustWork = TRUE
  )
  invisible(vapply(file.path(output_dir, c("a", "d", "g", "m", "p")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
  logo <- system.file("report_assets", "LISA_logo_C_compact_icon_muted_red_S.svg", package = "lisaR")
  if (!nzchar(logo) || !file.exists(logo)) stop("LISA-REPORT-PACKAGE-007 installed LISA logo asset is absent.", call. = FALSE)
  lisa_report_package_copy(logo, file.path(output_dir, "a", "logo.svg"))
  # Canonical tables occur exactly once in d/Tnnnn.tsv. CSV is a declared
  # companion rendering, never a second metadata/report copy.
  for (i in seq_len(nrow(tables))) {
    id <- tables$table_id[[i]]; lisa_report_package_copy(file.path(f$root, tables$path[[i]]), file.path(output_dir, "d", paste0(id, ".tsv")))
    lisa_report_package_csv_bom(read_lisa_tsv(file.path(output_dir, "d", paste0(id, ".tsv"))), file.path(output_dir, "d", paste0(id, ".csv")))
  }
  if (nrow(figures)) for (i in seq_len(nrow(figures))) {
    ext <- tolower(tools::file_ext(figures$path[[i]])); target <- file.path(output_dir, "g", paste0(figures$table_id[[i]], ".", ext))
    lisa_report_package_copy(file.path(f$root, figures$path[[i]]), target)
  }
  selected$output_path <- ifelse(selected$artifact_type == "table",
    file.path("d", paste0(selected$table_id, ".tsv")),
    file.path("g", paste0(selected$table_id, ".", tolower(tools::file_ext(selected$path)))))
  xlsx_info <- lisa_report_package_write_xlsx(tables, file.path(output_dir, "d"), file.path(output_dir, "d", "report.xlsx"))
  csv_files <- list.files(file.path(output_dir, "d"), pattern = "\\.csv$", full.names = TRUE)
  zip_cmd <- Sys.getenv("R_ZIPCMD", unset = "")
  if (!nzchar(zip_cmd)) zip_cmd <- unname(Sys.which("zip"))
  if (!nzchar(zip_cmd)) stop("LISA-REPORT-PACKAGE-014 a zip executable is required to build the CSV archive.", call. = FALSE)
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(file.path(output_dir, "d")); utils::zip("tables.zip", basename(csv_files), flags = "-q", zip = zip_cmd)
  css <- "body{font-family:system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#18232c}header{border-bottom:4px solid #9e2936;padding-bottom:1rem}nav a{margin-right:1rem}table{border-collapse:collapse;width:100%;font-size:.82rem}th,td{padding:.35rem;border-bottom:1px solid #ddd;text-align:left}.card{margin:.5rem 0;padding:.7rem;background:#f4f6f7;border-radius:.3rem}.warn{color:#8a3b00}.muted{color:#52616b}input{padding:.45rem;width:20rem;max-width:100%}"
  writeLines(css, file.path(output_dir, "a/r.css"), useBytes = TRUE)
  search <- list(artifacts = as.character(selected$table_id), roles = as.character(selected$role), analyses = if ("analysis_id" %in% names(selected)) unique(as.character(selected$analysis_id)) else character(), contrasts = if ("contrast_id" %in% names(selected)) unique(as.character(selected$contrast_id)) else character())
  js <- paste0("const LISA_SEARCH=", jsonlite::toJSON(search, auto_unbox=TRUE), ";function lisaSearch(q){q=(q||'').toLowerCase();document.querySelectorAll('[data-s]').forEach(e=>e.style.display=e.dataset.s.toLowerCase().includes(q)?'':'none');let a=Object.values(LISA_SEARCH).flat().filter(x=>String(x).toLowerCase().includes(q)).slice(0,15);document.getElementById('sr').textContent=q?(a.length?' Matches: '+a.join(', '):' No match'):''}")
  writeLines(js, file.path(output_dir, "a/s.js"), useBytes = TRUE)
  writep <- function(path, title, body) writeLines(lisa_report_package_page(title, body), file.path(output_dir, path), useBytes = TRUE)
  link_rows <- function(z) paste(vapply(seq_len(nrow(z)), function(i) {
    links <- if (z$artifact_type[[i]] == "table") paste0("<a href=\"d/", z$table_id[[i]], ".tsv\">TSV</a> <a href=\"d/", z$table_id[[i]], ".csv\">CSV</a>") else {
      preview <- if (tolower(tools::file_ext(z$output_path[[i]])) %in% c("png", "svg")) paste0("<br><a href=\"", z$output_path[[i]], "\"><img loading=\"lazy\" src=\"", z$output_path[[i]], "\" alt=\"", lisa_report_package_escape(z$table_id[[i]]), "\" height=\"110\"></a>") else ""
      paste0("<a href=\"", z$output_path[[i]], "\">Figure</a>", preview)
    }
    paste0("<tr data-s=\"", lisa_report_package_escape(paste(z[i, , drop=TRUE], collapse=" ")), "\"><td>", lisa_report_package_escape(z$table_id[[i]]), "</td><td>", lisa_report_package_escape(z$role[[i]]), "</td><td>", links, "</td></tr>")
  }, character(1)), collapse="")
  analysis_links <- paste(vapply(seq_len(nrow(f$de)), function(i) sprintf("<li><a href=\"A%02d.html\">%s</a></li>", i, lisa_report_package_escape(if ("label" %in% names(f$de)) f$de$label[[i]] else f$de$analysis_id[[i]])), character(1)), collapse="")
  contrast_links <- paste(vapply(seq_len(nrow(f$contrast)), function(i) sprintf("<li><a href=\"C%02d.html\">%s</a></li>", i, lisa_report_package_escape(if ("contrast_label" %in% names(f$contrast)) f$contrast$contrast_label[[i]] else f$contrast$contrast_id[[i]])), character(1)), collapse="")
  write_registry_pages <- function(prefix, title, intro, z, page_size = 180L) {
    page_count <- max(1L, ceiling(nrow(z) / page_size))
    for (page in seq_len(page_count)) {
      idx <- if (nrow(z)) seq.int((page - 1L) * page_size + 1L, min(page * page_size, nrow(z))) else integer()
      part <- z[idx, , drop = FALSE]
      page_path <- if (page == 1L) paste0(prefix, ".html") else sprintf("%s_%02d.html", prefix, page)
      nav <- if (page_count > 1L) paste0("<p>Pages: ", paste(vapply(seq_len(page_count), function(j) {
        href <- if (j == 1L) paste0(prefix, ".html") else sprintf("%s_%02d.html", prefix, j)
        if (j == page) paste0("<strong>", j, "</strong>") else paste0("<a href=\"", href, "\">", j, "</a>")
      }, character(1)), collapse = " "), "</p>") else ""
      writep(page_path, paste0(title, if (page_count > 1L) paste0(" - page ", page, "/", page_count) else ""), paste0(intro, nav, "<table><tr><th>ID</th><th>Role</th><th>Data</th></tr>", link_rows(part), "</table>", nav))
    }
  }
  writep("index.html", "LISA offline report", paste0("<section class=\"card\"><h2>Project identity and global QC: PASS</h2><p>Registry-driven archived-data bundle. <a href=\"qc.html\">QC</a> | <a href=\"dl.html\">Downloads</a> | <a href=\"methods.html\">Methods</a> | <a href=\"prov.html\">Provenance</a></p></section><section><h2>Analyses</h2><ul>", analysis_links, "</ul><h2>Contrasts</h2><ul>", contrast_links, "</ul></section><p class=\"warn\">Scientific warnings and disabled/not-applicable products are shown on QC; no missing evidence is rendered as zero.</p>"))
  for (i in seq_len(nrow(f$de))) { key <- as.character(f$de$analysis_id[[i]]); z <- if ("analysis_id" %in% names(selected)) selected[lisa_report_package_has_alias(selected$analysis_id, key), , drop=FALSE] else selected; write_registry_pages(sprintf("A%02d", i), paste("Analysis", key), "<h2>Collections, gene sets, gene/protein evidence and requested archived visuals</h2><p>Every listed artifact is checksummed and directly downloadable.</p>", z) }
  for (i in seq_len(nrow(f$contrast))) { key <- as.character(f$contrast$contrast_id[[i]]); z <- if ("contrast_id" %in% names(selected)) selected[lisa_report_package_has_alias(selected$contrast_id, key), , drop=FALSE] else selected; write_registry_pages(sprintf("C%02d", i), paste("Contrast", key), "<h2>A, B, delta, paired evidence and requested archived visuals</h2><p>Archived contrast tables retain A/B/delta semantics and status fields.</p>", z) }
  if ("analysis_collection" %in% names(f$collection)) for (i in seq_len(nrow(f$collection))) { key <- as.character(f$collection$analysis_collection[[i]]); z <- if ("collection" %in% names(selected)) selected[lisa_report_package_has_alias(selected$collection, key), , drop=FALSE] else selected; write_registry_pages(sprintf("K%02d", i), paste("Collection", key), "<h2>Complete collection registry</h2>", z) }
  write_registry_pages("dl", "Downloads", paste0("<p><a href=\"d/report.xlsx\">Consolidated XLSX</a> | <a href=\"d/tables.zip\">UTF-8-BOM CSV ZIP</a></p><p>Workbook tables: ", xlsx_info$included, "; explicit CSV-only tables: ", xlsx_info$csv_only, "; included figures: ", nrow(figures), ".</p>"), selected)
  status <- lisa_report_package_status(f, report_options, selected); write.table(status, file.path(output_dir, "m/qc.tsv"), sep="\t", quote=FALSE,row.names=FALSE); writep("qc.html", "Quality control", paste0("<h2>Global status: PASS</h2><p>File, link, checksum, registry and export validations passed. Requested products absent from the canonical source are explicitly marked not_available.</p>", lisa_report_package_table(status)))
  option_table <- data.frame(option = names(report_options), value = unlist(report_options, use.names = FALSE), source = publication_profile, stringsAsFactors = FALSE)
  write.table(option_table, file.path(output_dir, "m/options.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  methods <- paste0("# Methods\n\nThis report packages existing archived LISA outputs without recalculation. Publication profile: `", publication_profile, "`. Effective report options are recorded in `m/options.tsv`; requested products absent from the canonical archive are marked `not_available` in `m/qc.tsv`. lisa-gps name/version/formula/components are retained from the executed contract.\n")
  writeLines(methods, file.path(output_dir, "m/methods.md"), useBytes=TRUE); writep("methods.html", "Methods", paste0("<pre>", lisa_report_package_escape(methods), "</pre>"))
  prov <- data.frame(key=c("report_contract","publication_profile","fixture_root","lisaR_version"), value=c("report packaging final-user contract",publication_profile,"logical external fixture",as.character(utils::packageVersion("lisaR"))),stringsAsFactors=FALSE); write.table(prov,file.path(output_dir,"m/prov.tsv"),sep="\t",quote=FALSE,row.names=FALSE); writep("prov.html","Provenance",paste0("<p>Host-specific source paths are excluded from this portable report. Canonical input hashes are retained in the fixture and bundle manifests.</p>",lisa_report_package_table(prov),"<h2>Archived provenance</h2>",lisa_report_package_table(f$provenance)))
  write.table(selected, file.path(output_dir, "m/registry.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  files <- sort(list.files(output_dir, recursive=TRUE, full.names=FALSE)); files <- files[!grepl("^m/(manifest.tsv|SHA256SUMS)$", files)]
  manifest <- data.frame(path=files, sha256=vapply(file.path(output_dir,files), lisa_report_package_sha256, character(1)), stringsAsFactors=FALSE); write.table(manifest,file.path(output_dir,"m/manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  sums <- c(files,"m/manifest.tsv"); sum_hashes <- c(manifest$sha256, lisa_report_package_sha256(file.path(output_dir, "m/manifest.tsv"))); writeLines(paste(sum_hashes,sums),file.path(output_dir,"m/SHA256SUMS"),useBytes=TRUE)
  html <- list.files(output_dir, pattern="\\.html$", full.names=TRUE); rel <- list.files(output_dir, recursive=TRUE, full.names=FALSE)
  if (any(file.info(html)$size > 1024^2) || any(nchar(rel)>120L) || any(nchar(basename(rel))>48L)) stop("LISA-REPORT-PACKAGE-013 page or Windows path budget exceeded.",call.=FALSE)
  lisa_report_foundation_verify_links(output_dir)
  list(index=file.path(output_dir,"index.html"), files=length(rel), max_path=max(nchar(rel)), max_html=max(file.info(html)$size), xlsx_tables=xlsx_info$included, csv_only=xlsx_info$csv_only, figures=nrow(figures), publication_profile=publication_profile)
}

lisa_report_package_build_report <- function(fixture_dir, output_dir, publication_profile = "minimal", report_options = NULL) build_lisa_report_package(fixture_dir, output_dir, publication_profile, report_options)
