#!/usr/bin/env Rscript
# Copyright (C) 2026 CSIC. Author: David Olmeda Casadomé.
# GPL-3.0-or-later. Assemble a scoped report from saved results; no science.
# This is a report entry point, not a post-hoc HTML editor. It creates a fresh,
# independent input view and invokes the installed normal report generator.

selected_fail <- function(...) stop(..., call. = FALSE)
selected_read <- function(path) {
  if (!file.exists(path)) selected_fail("Missing saved table: ", path)
  # Canonical no-contrast runs serialize data.frame() as a blank line.
  if (file.info(path)$size <= 16L &&
      !any(nzchar(trimws(readLines(path, warn = FALSE))))) return(data.frame())
  utils::read.delim(path, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE)
}
selected_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
selected_sha <- function(path) digest::digest(file = path, algo = "sha256")
selected_ids <- function(value) {
  if (identical(value, "none")) return(character())
  ids <- strsplit(value, ",", fixed = TRUE)[[1L]]
  if (!length(ids) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", ids)) ||
      anyDuplicated(ids)) selected_fail("Supply unique, exact comma-separated IDs.")
  ids
}
# Named parse_args like every other installed post script, so the registered
# interface in lisa_post_script_interfaces() is validated against this exact
# parser by test-post-script-interfaces.R.
parse_args <- function(argv) {
  args <- list(source_dir = "", output_dir = "", analyses = "", contrasts = "",
    artifacts_dir = "", artifact_manifest = "", title = "Selected LISA FULL report",
    plan_only = "false", complete_missing = "false", assemble = "true", kegg_maps = "false",
    kegg_cache_root = "", kegg_snapshot_id = "", kegg_access_mode = "cache_only")
  if (!length(argv) || length(argv) %% 2L) selected_fail("Arguments must be --name value pairs.")
  seen <- character()
  for (i in seq.int(1L, length(argv), 2L)) {
    key <- gsub("-", "_", sub("^--", "", argv[[i]]))
    if (!startsWith(argv[[i]], "--") || !key %in% names(args) || key %in% seen)
      selected_fail("Unknown or repeated argument: ", argv[[i]])
    args[[key]] <- argv[[i + 1L]]; seen <- c(seen, key)
  }
  if (any(!nzchar(unlist(args[c("source_dir", "output_dir", "analyses", "contrasts")]))))
    selected_fail("Required: --source-dir --output-dir --analyses --contrasts; artifacts are optional.")
  for (key in c("plan_only", "complete_missing", "assemble", "kegg_maps"))
    if (!args[[key]] %in% c("true", "false")) selected_fail(key, " must be true or false.")
  if (xor(nzchar(args$artifacts_dir), nzchar(args$artifact_manifest)))
    selected_fail("Supply both artifact directory and manifest, or neither.")
  args
}
selected_contract <- function(source, analyses, contrasts) {
  de <- selected_read(file.path(source, "config", "de_index.tsv"))
  co <- selected_read(file.path(source, "config", "contrast_index.tsv"))
  if (!nrow(co) && !ncol(co)) co <- data.frame(contrast_id = character(),
    output_id = character(), contrast_a = character(), contrast_b = character())
  if (!all(c("analysis_id") %in% names(de)) || anyDuplicated(de$analysis_id) ||
      !all(c("contrast_id", "output_id", "contrast_a", "contrast_b") %in% names(co)) ||
      anyDuplicated(co$contrast_id)) selected_fail("Ambiguous or incomplete saved indices.")
  if (!all(analyses %in% de$analysis_id) || !all(contrasts %in% co$contrast_id))
    selected_fail("Requested analysis/contrast is absent from saved indices.")
  de <- de[match(analyses, de$analysis_id), , drop = FALSE]
  co <- co[match(contrasts, co$contrast_id), , drop = FALSE]
  if (!all(c(co$contrast_a, co$contrast_b) %in% analyses))
    selected_fail("Every contrast must refer only to selected analyses; no implicit additions.")
  owners <- if (nrow(co)) paste(co$contrast_id, co$output_id, sep = "_") else character()
  list(de = de, contrasts = co, owners = owners)
}
selected_rows <- function(x, contract) {
  keep <- rep(TRUE, nrow(x))
  for (col in intersect(c("analysis_id", "contrast_id", "contrast_name", "output_id"), names(x))) {
    allowed <- switch(col, analysis_id = contract$de$analysis_id,
      contrast_id = c(contract$contrasts$contrast_id, contract$owners),
      contrast_name = contract$owners, output_id = contract$contrasts$output_id)
    value <- as.character(x[[col]])
    keep <- keep & (is.na(value) | !nzchar(value) | value %in% allowed)
  }
  x[keep, , drop = FALSE]
}
selected_tree_files <- function(root) {
  if (!dir.exists(root) || nzchar(Sys.readlink(root))) selected_fail("Missing or linked tree: ", root)
  paths <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
    no.. = TRUE, include.dirs = TRUE)
  if (any(nzchar(Sys.readlink(paths)))) selected_fail("Symbolic links are not input files: ", root)
  paths <- paths[!dir.exists(paths)]
  if (length(paths) && !all(utils::file_test("-f", paths))) selected_fail("Nonregular input file.")
  sort(paths)
}
selected_gene_evidence <- function(source, target, analyses) {
  root <- file.path(source, "report_pages", "gene_evidence")
  if (!dir.exists(root)) return(invisible(NULL))
  tables <- c("analyses", "scopes", "identifiers", "de", "sets", "leading_edges",
              "memberships", "source_symbols", "provenance")
  evidence <- stats::setNames(lapply(tables, function(name)
    selected_read(file.path(root, "tables", paste0(name, ".tsv")))), tables)
  for (name in tables) if ("analysis_id" %in% names(evidence[[name]])) {
    value <- evidence[[name]]$analysis_id
    evidence[[name]] <- evidence[[name]][is.na(value) | !nzchar(value) | value %in% analyses, , drop = FALSE]
  }
  evidence$memberships <- evidence$memberships[evidence$memberships$pathway %in% evidence$sets$pathway, , drop = FALSE]
  symbols <- unique(c(evidence$de$symbol, evidence$memberships$symbol))
  evidence$identifiers <- evidence$identifiers[evidence$identifiers$symbol %in% symbols, , drop = FALSE]
  evidence$source_symbols <- evidence$source_symbols[evidence$source_symbols$symbol %in% symbols, , drop = FALSE]
  # Preserve the saved scientific metadata, not guessed thresholds. The JSON
  # payload begins with metadata followed by analyses in the installed writer.
  html <- readLines(file.path(root, "index.html"), warn = FALSE)
  line <- html[grepl('id="gene-evidence-data"', html, fixed = TRUE)]
  if (length(line) != 1L) selected_fail("Missing saved gene-explorer metadata.")
  prefix <- sub('^.*id="gene-evidence-data">\\{"metadata":', "", line)
  meta <- strsplit(prefix, ',"analyses":', fixed = TRUE)[[1L]][[1L]]
  evidence$metadata <- jsonlite::fromJSON(meta, simplifyVector = FALSE)
  rm(html, line, prefix); invisible(gc())
  class(evidence) <- "lisa_gene_evidence"
  get("render_lisa_gene_evidence", asNamespace("lisaR"))(
    evidence, file.path(target, "report_pages", "gene_evidence"))
}

selected_main <- function(argv = commandArgs(TRUE)) {
  args <- parse_args(argv)
  for (pkg in c("lisaR", "digest", "jsonlite"))
    if (!requireNamespace(pkg, quietly = TRUE)) selected_fail("Missing installed package: ", pkg)
  source <- normalizePath(args$source_dir, winslash = "/", mustWork = TRUE)
  artifact_root <- if (nzchar(args$artifacts_dir)) normalizePath(args$artifacts_dir, winslash = "/", mustWork = TRUE) else ""
  if (nzchar(artifact_root) != nzchar(args$artifact_manifest)) selected_fail("Artifacts require a manifest and vice versa.")
  out <- file.path(normalizePath(dirname(args$output_dir), winslash = "/", mustWork = TRUE), basename(args$output_dir))
  output_link <- Sys.readlink(out)
  if (file.exists(out) || (!is.na(output_link) && nzchar(output_link)) ||
      out == source || startsWith(out, paste0(source, "/")) ||
      startsWith(source, paste0(out, "/")) || (nzchar(artifact_root) && startsWith(out, paste0(artifact_root, "/"))))
    selected_fail("Output must be a new directory outside both source trees.")
  if (!length(selected_ids(args$analyses))) selected_fail("At least one analysis is required.")
  contract <- selected_contract(source, selected_ids(args$analyses), selected_ids(args$contrasts))
  roots <- c(file.path("outputs", "single_de", contract$de$analysis_id),
    if (length(contract$owners)) file.path("outputs", "category_contrasts", contract$owners),
    file.path("report_pages", "evidence", contract$de$analysis_id),
    if (length(contract$owners)) file.path("report_pages", "contrast_evidence", contract$owners),
    file.path("report_pages", "category_navigation", c(contract$de$analysis_id, contract$owners)))
  # Historical projects may not have category navigation; required scientific
  # results and evidence are never treated as an unselected scope when absent.
  roots <- roots[!startsWith(roots, "report_pages/category_navigation/") |
                   dir.exists(file.path(source, roots))]
  gene_roots <- c(file.path("outputs/gene_level/single_de", contract$de$analysis_id),
    if (length(contract$owners)) file.path("outputs/gene_level/category_contrasts", contract$owners))
  roots <- c(roots, gene_roots[dir.exists(file.path(source, gene_roots))])
  sources <- unlist(lapply(file.path(source, roots), selected_tree_files), use.names = FALSE)
  rel <- substring(sources, nchar(source) + 2L)
  aro <- c(contract$de$analysis_id, if (length(contract$owners)) file.path("contrasts", contract$owners))
  if (nzchar(artifact_root)) {
  afiles <- unlist(lapply(file.path(artifact_root, aro), selected_tree_files), use.names = FALSE)
  arel <- paste0("artifacts/", substring(afiles, nchar(artifact_root) + 2L))
  ledger <- selected_read(args$artifact_manifest)
  if (!all(c("path", "bytes", "sha256") %in% names(ledger)) || anyDuplicated(ledger$path))
    selected_fail("Invalid artifact manifest.")
  idx <- match(arel, ledger$path)
  if (anyNA(idx) || !length(afiles)) selected_fail("Selected artifact absent from manifest.")
  if (any(as.numeric(file.info(afiles)$size) != as.numeric(ledger$bytes[idx])) ||
      any(vapply(afiles, selected_sha, character(1L)) != ledger$sha256[idx]))
    selected_fail("Saved extra figures fail manifest identity checks.")
  } else { afiles <- arel <- character() }
  if (anyDuplicated(c(rel, arel))) selected_fail("Colliding source paths.")
  summary <- list(analyses = as.character(contract$de$analysis_id),
    contrasts = as.character(contract$contrasts$contrast_id),
    contrast_a = as.character(contract$contrasts$contrast_a),
    contrast_b = as.character(contract$contrasts$contrast_b),
    extra_files = length(afiles), source_files = length(sources),
    extra_by_format = as.list(table(tolower(tools::file_ext(afiles)))),
    science_run = FALSE, figures_rendered = FALSE)
  cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
  if (args$plan_only == "true") return(invisible(summary))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) selected_fail("Cannot create output directory.")
  # Independent regular copies: no write in this generation can reach a
  # standard/source inode. Original receipts are not rewritten or relabelled.
  all_source <- c(sources, afiles); all_rel <- c(rel, arel)
  before <- vapply(all_source, selected_sha, character(1L))
  saved_manifest <- selected_read(file.path(source, "run_manifest.tsv"))
  smatch <- match(rel, saved_manifest$path)
  if (anyNA(smatch) || any(before[seq_along(sources)] != saved_manifest$sha256[smatch]))
    selected_fail("Saved standard files fail their original manifest.")
  audit <- paste0(out, ".audit")
  dir.create(audit, showWarnings = FALSE)
  selected_write(data.frame(path = all_rel, sha256 = unname(before)),
    file.path(audit, "source_inventory_before.tsv"))
  for (i in seq_along(all_source)) {
    dest <- file.path(out, all_rel[[i]])
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(all_source[[i]], dest, overwrite = FALSE, copy.mode = TRUE))
      selected_fail("Cannot copy saved input: ", all_rel[[i]])
  }
  # The indices and status tables are report views of saved rows, not a new
  # analysis receipt. Their scientific source identities are recorded below.
  selected_write(contract$de, file.path(out, "config", "de_index.tsv"))
  selected_write(contract$contrasts, file.path(out, "config", "contrast_index.tsv"))
  for (name in c("single_de_status.tsv", "contrast_status.tsv", "post_lisa_status.tsv")) {
    src <- file.path(source, name)
    if (file.exists(src)) selected_write(selected_rows(selected_read(src), contract), file.path(out, name))
  }
  for (name in c("lisa_pipeline_metadata.tsv", "contract_manifest.tsv", "report_output_policy.tsv")) {
    src <- file.path(source, name)
    if (!file.exists(src)) selected_fail("Required saved contract is missing: ", name)
    file.copy(src, file.path(out, name), overwrite = FALSE)
  }
  selection <- c(schema_version = "1", report_mode = "full",
    analyses = paste(contract$de$analysis_id, collapse = ";"),
    contrasts = paste(contract$contrasts$contrast_id, collapse = ";"),
    contrast_owners = paste(contract$owners, collapse = ";"),
    assembly = if (args$complete_missing == "true") "saved results; existing figures reused; missing presentation products only" else "saved results only; no DE, GSEA or figure regeneration",
    de_index_source_sha256 = selected_sha(file.path(source, "config", "de_index.tsv")),
    contrast_index_source_sha256 = selected_sha(file.path(source, "config", "contrast_index.tsv")),
    artifact_manifest_sha256 = if (nzchar(args$artifact_manifest)) selected_sha(args$artifact_manifest) else "none")
  selected_write(data.frame(key = names(selection), value = unname(selection)),
    file.path(out, "report_selection.tsv"))
  # Rebuild only the gene-search presentation from saved tables, filtered to
  # the same analyses. Copying the global search would leak excluded analyses.
  selected_gene_evidence(source, out, contract$de$analysis_id)
  generator <- system.file("scripts", "build_LISA_report.R", package = "lisaR")
  if (!nzchar(generator)) selected_fail("Installed report generator not found.")
  if (args$complete_missing == "true") {
    get("lisa_complete_full_products", asNamespace("lisaR"))(out,
      kegg_maps = args$kegg_maps == "true", kegg_cache_root = args$kegg_cache_root,
      kegg_snapshot_id = args$kegg_snapshot_id, kegg_access_mode = args$kegg_access_mode)
  }
  if (args$assemble == "false") {
    jsonlite::write_json(summary, file.path(audit, "prepared.json"), auto_unbox = TRUE, pretty = TRUE)
    return(invisible(summary))
  }
  command <- c("--vanilla", shQuote(generator), "--project-dir", shQuote(out),
    "--title", shQuote(args$title))
  status <- system2(file.path(R.home("bin"), "Rscript"), command)
  if (status != 0L) selected_fail("Report generator failed; unaccepted output retained: ", out)
  nav <- jsonlite::fromJSON(file.path(out, "report_pages", "navigation_inventory.json"), simplifyVector = FALSE)
  kinds <- vapply(nav$contexts, function(x) x$kind, character(1L))
  ids <- vapply(nav$contexts, function(x) x$id, character(1L))
  if (!setequal(ids[kinds == "analysis"], contract$de$analysis_id) ||
      !setequal(ids[kinds == "contrast"], contract$owners)) selected_fail("Rendered navigation differs from exact selection.")
  for (context in nav$contexts) for (collection in context$collections) {
    if (!any(vapply(collection$sections, function(x) identical(x$kind, "full-products"), logical(1L))))
      selected_fail("Visible FULL section missing: ", context$id, "/", collection$id)
  }
  renderer <- new.env(parent = asNamespace("lisaR"))
  sys.source(generator, renderer)
  attached <- renderer$report_category_product_reference_rows(out)
  attached_paths <- unique(as.character(attached$source_path))
  if (!setequal(attached_paths[startsWith(attached_paths, "artifacts/")], arel))
    selected_fail("Saved FULL artifacts are not all attached to visible category products: ",
      length(setdiff(arel, attached_paths)), " missing files.")
  after <- vapply(all_source, selected_sha, character(1L))
  if (!identical(before, after)) selected_fail("Source input changed during assembly.")
  preserved <- startsWith(all_rel, "outputs/") | startsWith(all_rel, "artifacts/") |
    grepl("/figures/", all_rel, fixed = TRUE)
  out_sha <- vapply(file.path(out, all_rel[preserved]), selected_sha, character(1L))
  if (!identical(unname(out_sha), unname(before[preserved]))) selected_fail("Reused scientific files changed.")
  selected_write(data.frame(path = all_rel, sha256 = unname(before),
    disposition = ifelse(preserved, "reused byte-identically", "saved presentation refreshed")),
    file.path(audit, "source_inventory.tsv"))
  summary$status <- "ASSEMBLED_AND_CHECKED; visual review still required"
  summary$generator_sha256 <- selected_sha(generator)
  jsonlite::write_json(summary, file.path(audit, "result.json"), auto_unbox = TRUE, pretty = TRUE)
  cat("Selected FULL report: ", file.path(out, "report_index.html"), "\n", sep = "")
  invisible(summary)
}
if (sys.nframe() == 0L) selected_main()
