# Called only on a new, private prepared-project staging directory after the
# original bundle and copied inputs have passed their existing checks.
lisa_compact_example_project <- function(project, example, config_path,
    resource_cache = getOption("lisaR.dictionary_cache_root", NULL)) {
  project <- lisa_existing_run_root(project)
  mapping <- read_lisa_tsv(system.file("examples", "project_paths.tsv", package = "lisaR"))
  mapping <- mapping[mapping$example == example, , drop = FALSE]
  stopifnot(nrow(mapping) > 0L, !anyDuplicated(tolower(mapping$path)))
  old <- file.path(project, mapping$source)
  target <- file.path(project, mapping$path)
  if (!all(file.exists(old)) || any(file.exists(target)))
    stop("Prepared-project inputs are missing or short destinations already exist.", call. = FALSE)
  for (p in c(old, target)) lisa_guarded_path(p, project)
  mapping$sha256 <- vapply(old, lisa_sha256_file, character(1))
  mapping$bytes <- as.numeric(file.info(old)$size)
  cfg <- read_lisa_pipeline_config(config_path)
  rewrite <- function(x) {
    if (is.list(x)) return(lapply(x, rewrite))
    if (!is.character(x)) return(x)
    for (i in seq_along(x)) {
      hit <- match(sub("^[.][.]/", "", x[[i]]), mapping$source)
      if (!is.na(hit)) x[[i]] <- mapping$path[[hit]]
    }
    x
  }
  cfg <- rewrite(cfg)
  if (identical(example, "riaz")) {
    # The tutorial executes PRE and ON only, not the full upstream study.
    ids <- c("responders_vs_pd_pre", "responders_vs_pd_on")
    cfg$single_de <- Filter(function(x) x$analysis_id %in% ids, cfg$single_de)
    cfg$contrasts <- Filter(function(x)
      identical(x$contrast_id, "response_separation_by_time"), cfg$contrasts)
    stopifnot(identical(vapply(cfg$single_de, `[[`, character(1), "analysis_id"), ids),
      length(cfg$contrasts) == 1L,
      identical(cfg$contrasts[[1L]]$analysis_a, ids[[2L]]),
      identical(cfg$contrasts[[1L]]$analysis_b, ids[[1L]]))
    unused <- mapping$path[grepl("/(responders_on_vs_pre|pd_on_vs_pre|differential_longitudinal_response)[.]tsv$", mapping$source)]
    cfg$allowlisted_source_paths <- cfg$allowlisted_source_paths[!cfg$allowlisted_source_paths %in% unused]
    cfg$source_data <- Filter(function(x) !x$path %in% unused, cfg$source_data)
  }
  cfg$pipeline$output_dir <- "out"
  cfg$pipeline$file_label_prefix <- example
  for (i in seq_along(old)) {
    lisa_guarded_dir_create(dirname(target[[i]]), project)
    lisa_guarded_rename(old[[i]], target[[i]], overwrite = FALSE, run_root = project)
    if (!identical(lisa_sha256_file(target[[i]]), mapping$sha256[[i]]))
      stop("Short input copy changed its source bytes.", call. = FALSE)
  }
  # Keep the pre-mapping configuration as provenance, not as a runnable guide.
  original <- file.path(project, "provenance", "source_study.yml")
  lisa_guarded_copy(config_path, original, overwrite = FALSE, run_root = project)
  source_json <- sub("[.]yml$", ".json", config_path)
  if (file.exists(source_json)) lisa_guarded_copy(source_json,
    file.path(project, "provenance", "source_study.json"), FALSE, project)
  old_manifest <- file.path(dirname(config_path), "config_manifest.tsv")
  if (file.exists(old_manifest)) lisa_guarded_rename(old_manifest,
    file.path(project, "provenance", "source_config_manifest.tsv"), FALSE, project)
  for (p in c(config_path, sub("[.]yml$", ".json", config_path))) {
    if (file.exists(p)) lisa_guarded_delete(p, run_root = project)
  }
  new_config <- file.path(project, "study.yml")
  lisa_guarded_write(new_config, function(p) yaml::write_yaml(cfg, p), project, FALSE)
  json_path <- file.path(project, "study.json")
  lisa_guarded_write(json_path, function(p) jsonlite::write_json(cfg, p,
    auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), project, FALSE)
  write_lisa_tsv(mapping, file.path(project, "provenance", "path_map.tsv"))
  if (!is.null(resource_cache)) write_lisa_tsv(
    data.frame(resource_cache = resource_cache),
    file.path(project, "provenance", "resource_cache.tsv"))
  config_manifest <- data.frame(path = c("study.yml", "study.json"),
    sha256 = vapply(c(new_config, json_path), lisa_sha256_file, character(1)))
  write_lisa_tsv(config_manifest, file.path(project, "provenance", "short_config.tsv"))
  lisa_verify_compact_example_project(project)
  invisible(new_config)
}

lisa_verify_compact_example_project <- function(project) {
  project <- lisa_existing_run_root(project)
  cache_record <- file.path(project, "provenance", "resource_cache.tsv")
  if (file.exists(cache_record)) {
    old <- options(lisaR.dictionary_cache_root = read_lisa_tsv(cache_record)$resource_cache[[1L]])
    on.exit(options(old), add = TRUE)
  }
  mapping <- read_lisa_tsv(file.path(project, "provenance", "path_map.tsv"))
  configs <- read_lisa_tsv(file.path(project, "provenance", "short_config.tsv"))
  stopifnot(all(c("path", "sha256", "bytes") %in% names(mapping)),
    nrow(mapping) > 0L, !anyDuplicated(tolower(mapping$path)),
    identical(as.character(configs$path), c("study.yml", "study.json")))
  paths <- file.path(project, c(mapping$path, configs$path))
  for (p in paths) lisa_assert_regular_managed_file(p, project)
  if (!identical(unname(vapply(paths, lisa_sha256_file, character(1))),
      as.character(c(mapping$sha256, configs$sha256))) ||
      !identical(as.numeric(file.info(file.path(project, mapping$path))$size),
        as.numeric(mapping$bytes)))
    stop("Prepared project differs from its short-path manifest.", call. = FALSE)
  validation <- validate_lisa_config(file.path(project, "study.yml"),
    check_files = TRUE, strict = TRUE)
  if (!isTRUE(validation$execution_ready))
    stop("The short-path project is not execution ready; inspect validation.", call. = FALSE)
  invisible(TRUE)
}
