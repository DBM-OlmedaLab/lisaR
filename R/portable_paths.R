# Physical names are separate from scientific identifiers and labels.
# Existing runs retain their recorded paths and are never renamed.
lisa_path_utf16_length <- function(path) {
  vapply(enc2utf8(path), function(x) {
    bytes <- iconv(x, from = "UTF-8", to = "UTF-16LE", toRaw = TRUE)[[1L]]
    if (is.null(bytes)) stop("A path cannot be encoded for Windows.", call. = FALSE)
    length(bytes) / 2L
  }, numeric(1), USE.NAMES = FALSE)
}

lisa_assert_portable_path <- function(path, windows = .Platform$OS.type == "windows",
                                      limit = 240L) {
  if (!isTRUE(windows) || !length(path)) return(invisible(path))
  sizes <- lisa_path_utf16_length(gsub("\\\\", "/", path))
  bad <- which(sizes > limit)
  if (length(bad)) {
    i <- bad[[1L]]
    stop("LISA-PATH-001 path requires ", sizes[[i]],
      " UTF-16 units; the portable budget is ", limit, ": ", path[[i]],
      ". Use a shorter project/output directory (for example D:/lisa/riaz), ",
      "or shorter analysis and contrast IDs, then validate again. ",
      "Existing results have not been renamed.", call. = FALSE)
  }
  parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)
  bad_name <- vapply(parts, function(x) any(grepl(
    "^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])([.]|$)|[. ]$",
    x[nzchar(x)], ignore.case = TRUE)), logical(1))
  if (any(bad_name)) stop("LISA-PATH-002 Windows reserved filename: ",
    path[which(bad_name)[[1L]]], ". Choose another name.", call. = FALSE)
  invisible(path)
}

lisa_short_staging_path <- function(output_dir, run_id) {
  # Bind both identities, retain exclusive creation and collision refusal.
  key <- digest::digest(paste(output_dir, run_id, sep = "\n"),
    algo = "sha256", serialize = FALSE)
  file.path(dirname(output_dir), paste0(".ls-", substr(key, 1L, 12L)))
}

lisa_contrast_plot_stem <- function(comparison, prefix, plot_set) {
  key <- digest::digest(paste(comparison, prefix, sep = "\n"),
    algo = "sha256", serialize = FALSE)
  subset <- c(all = "all", same_direction = "same", opposite_direction = "opp")
  stopifnot(plot_set %in% names(subset))
  paste0("cx_", substr(key, 1L, 12L), "_", subset[[plot_set]])
}

lisa_check_run_paths <- function(output_dir, de_index, contrast_index,
    collections, prefix, report, windows = .Platform$OS.type == "windows") {
  if (!isTRUE(windows)) return(invisible(TRUE))
  ledger <- lisa_expected_artifacts(de_index, contrast_index,
    lisa_collection_registry(collections), output_dir,
    file_label_prefix = prefix,
    plot_formats = names(report$formats)[unlist(report$formats)],
    source_data = isTRUE(report$source_data), recipes = isTRUE(report$recipes),
    requested_report_mode = report$mode,
    category_evidence = isTRUE(report$category_evidence),
    category_nes_variants = report$category_nes_variants)
  paths <- unlist(strsplit(ledger$witnesses, "|", fixed = TRUE), use.names = FALSE)
  paths <- sub("^(file|dir):", "", paths[grepl("^(file|dir):", paths)])
  stage <- lisa_short_staging_path(output_dir, "planned")
  destinations <- c(file.path(output_dir, paths), file.path(stage, paths))
  temporary <- file.path(dirname(destinations), ".lisa-write-123456789abc.json")
  lisa_assert_portable_path(c(destinations, temporary, paste0(output_dir, ".lisa.lock")), windows = windows)
  invisible(TRUE)
}
