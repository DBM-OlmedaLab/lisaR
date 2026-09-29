skip_unless_installed_lisar <- function() {
  loaded_root <- system.file(package = "lisaR")
  skip_if_not(
    nzchar(loaded_root) &&
      file.exists(file.path(loaded_root, "DESCRIPTION")) &&
      file.exists(file.path(loaded_root, "Meta", "package.rds")),
    "requires an installed lisaR tree, not a pkgload source namespace"
  )
}

expect_lisa_parallel_contract <- function(output_dir, workers_requested) {
  path <- file.path(output_dir, "parallel_execution.tsv")
  expect_true(file.exists(path))
  execution <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE
  )
  required <- c(
    "platform", "backend_effective", "workers_requested", "workers_effective"
  )
  expect_gt(nrow(execution), 0L)
  expect_true(all(required %in% names(execution)))
  runtime_platform <- lisaR:::lisa_runtime_platform()
  requested <- as.integer(execution$workers_requested)
  effective <- as.integer(execution$workers_effective)
  backend <- as.character(execution$backend_effective)
  expect_true(all(execution$platform == runtime_platform))
  expect_true(all(requested == as.integer(workers_requested)))
  expect_true(all(effective >= 1L & effective <= requested))
  if (identical(runtime_platform, "linux")) {
    expected_backend <- ifelse(effective > 1L, "multicore", "sequential")
    expect_identical(backend, unname(expected_backend))
  } else {
    expect_true(all(effective == 1L))
    expect_true(all(backend == "sequential"))
  }
  invisible(execution)
}

test_that("sample-project initializer is additive and preserves existing paths", {
  target <- tempfile("lisaR-quick-start-")
  project <- lisa_init_project(target)
  expect_identical(project, normalizePath(target, winslash = "/"))
  expect_true(file.exists(file.path(project, "study.yml")))
  expect_true(file.exists(file.path(project, "data", "de_a.tsv")))
  expect_true(file.exists(file.path(project, "data", "de_b.tsv")))
  expect_true(file.exists(file.path(project, "data", "expr.tsv")))
  config <- yaml::read_yaml(file.path(project, "study.yml"))
  expect_identical(config$pipeline$schema_version, "1.0.0")
  expect_identical(config$pipeline$dictionary_resource, "lisa_dictionary_quickstart@1.1.0")
  expect_identical(config$report$mode, "standard")
  expect_identical(
    unlist(config$report$category_nes_variants, use.names = FALSE),
    c("clean", "percentages", "direction", "dispersion")
  )
  expect_error(lisa_init_project(target), "will not overwrite")

  occupied <- tempfile("lisaR-quick-start-occupied-")
  dir.create(occupied)
  writeLines("preserve", file.path(occupied, "user-file.txt"))
  expect_error(lisa_init_project(occupied), "will not overwrite")
  expect_identical(readLines(file.path(occupied, "user-file.txt")), "preserve")
})

test_that("sample-project initializer preserves an explicit parent", {
  fake_home <- tempfile("lisaR-quick-home-")
  dir.create(fake_home)
  on.exit(unlink(fake_home, recursive = TRUE, force = TRUE), add = TRUE)
  if (.Platform$OS.type != "windows") Sys.chmod(fake_home, "0755")
  home_mode <- as.integer(file.info(fake_home)$mode)
  project <- lisa_init_project(file.path(fake_home, "sample-project"))
  expect_true(file.exists(file.path(project, "study.yml")))
  expect_identical(dirname(project), normalizePath(fake_home, winslash = "/"))
  expect_identical(as.integer(file.info(fake_home)$mode), home_mode)
})

test_that("a fresh R process resolves a tilde child from startup HOME", {
  skip_unless_installed_lisar()
  fake_home <- tempfile("lisaR-startup-home-")
  dir.create(fake_home)
  on.exit(lisa_test_cleanup_path(fake_home), add = TRUE)
  child_script <- tempfile("lisaR-startup-home-", fileext = ".R")
  on.exit(unlink(child_script, force = TRUE), add = TRUE)
  writeLines(c(
    "library(lisaR)",
    paste0(
      "cat('LISAR_PACKAGE_ROOT=', normalizePath(system.file(package = ",
      "'lisaR'), winslash = '/', mustWork = TRUE), '\\n', sep = '')"
    ),
    "project <- lisa_init_project('~/child')",
    paste0(
      "cat('LISAR_HOME_CHILD=', normalizePath(project, winslash = '/', ",
      "mustWork = TRUE), '\\n', sep = '')"
    )
  ), child_script, useBytes = TRUE)

  library_value <- paste(unique(.libPaths()), collapse = .Platform$path.sep)
  withr::local_envvar(c(
    HOME = fake_home,
    R_USER = fake_home,
    R_LIBS = library_value
  ))
  output <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(),
    c("--vanilla", shQuote(child_script)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(output, collapse = "\n"))

  package_line <- grep("^LISAR_PACKAGE_ROOT=", output, value = TRUE)
  child_line <- grep("^LISAR_HOME_CHILD=", output, value = TRUE)
  expect_length(package_line, 1L)
  expect_length(child_line, 1L)
  expect_identical(
    sub("^LISAR_PACKAGE_ROOT=", "", package_line),
    normalizePath(system.file(package = "lisaR"), winslash = "/",
                  mustWork = TRUE)
  )
  expect_identical(
    sub("^LISAR_HOME_CHILD=", "", child_line),
    normalizePath(file.path(fake_home, "child"), winslash = "/",
                  mustWork = TRUE)
  )
})

test_that("installation smoke stays temporary and has no visible installed example", {
  skip_unless_installed_lisar()
  fixture <- testthat::test_path("fixtures", "installation-smoke")
  project <- tempfile("lisaR-installation-smoke-")
  dir.create(project)
  expect_true(all(file.copy(list.files(fixture, full.names = TRUE), project)))
  config <- file.path(project, "study.yml")
  expect_true(validate_lisa_config(config)$valid)
  result <- suppressWarnings(run_lisa(config))
  expect_identical(verify_lisa_run(result$output_dir)$gate, "PASS")
  expect_false(nzchar(system.file("examples", "installation-smoke", package = "lisaR")))
})

test_that("installed sample project is a real offline full run", {
  # This is deliberately an installed-package integration test. Under
  # pkgload, the parent uses source code but each external Rscript would load
  # whichever lisaR happens to be installed in .libPaths(), so that mixed-code
  # execution must not be treated as evidence. R CMD check and the clean
  # installed-package gate exercise the complete test.
  skip_unless_installed_lisar()
  proxy_values <- c(
    http_proxy = "http://127.0.0.1:9",
    https_proxy = "http://127.0.0.1:9",
    HTTP_PROXY = "http://127.0.0.1:9",
    HTTPS_PROXY = "http://127.0.0.1:9"
  )
  previous_proxy <- Sys.getenv(names(proxy_values), unset = NA_character_)
  on.exit({
    present <- !is.na(previous_proxy)
    if (any(present)) do.call(Sys.setenv, as.list(previous_proxy[present]))
    if (any(!present)) Sys.unsetenv(names(previous_proxy)[!present])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(proxy_values))
  project <- lisa_init_project(tempfile("lisaR-quick-start-run-"))
  config <- file.path(project, "study.yml")
  raw_config <- yaml::read_yaml(config)
  # The distributed Quick Start is standard. Request full explicitly here
  # to retain coverage of the integrated FULL report contract.
  expect_identical(raw_config$report$mode, "standard")
  raw_config$report$mode <- "full"
  raw_config$report$recipes <- TRUE
  yaml::write_yaml(raw_config, config)
  validation <- validate_lisa_config(config)
  expect_true(validation$valid)
  expect_identical(validation$analyses, 2L)
  expect_identical(validation$contrasts, 1L)
  expect_identical(as.integer(validation$config$pipeline$workers), 4L)
  expect_false(isTRUE(validation$config$pipeline$dry_run))
  expect_identical(validation$config$report$mode, "full")
  expect_true(isTRUE(validation$config$report$formats$png))
  expect_true(isTRUE(validation$config$report$source_data))
  expect_false(isTRUE(validation$config$pipeline$run_kegg_maps))
  expect_false(isTRUE(validation$config$pipeline$run_ora))

  # Keep the installed-package smoke inside one process. Forked graphics
  # devices can deadlock under the R CMD check/testthat harness; multicore
  # behavior is covered independently by test-workers.R and the clean-install
  # integration gate. The requested worker count remains part of the receipt.
  previous_ncpus <- Sys.getenv("NCPUS", unset = NA_character_)
  on.exit({
    if (is.na(previous_ncpus)) Sys.unsetenv("NCPUS") else Sys.setenv(NCPUS = previous_ncpus)
  }, add = TRUE)
  Sys.setenv(NCPUS = "1")
  result <- suppressWarnings(run_lisa(config))
  check <- verify_lisa_run(result$output_dir)
  expect_lisa_parallel_contract(result$output_dir, workers_requested = 4L)
  expect_identical(check$gate, "PASS")
  code_ledger <- file.path(result$output_dir, "code_identity.tsv")
  run_manifest <- utils::read.delim(
    file.path(result$output_dir, "run_manifest.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  expect_true(file.exists(code_ledger))
  expect_true("code_identity.tsv" %in% run_manifest$path)
  expect_identical(
    run_manifest$sha256[match("code_identity.tsv", run_manifest$path)],
    lisaR:::lisa_sha256_file(code_ledger)
  )
  expect_true(file.exists(file.path(result$output_dir, "report_index.html")))
  expect_true(file.exists(file.path(result$output_dir, "report_output_policy.tsv")))
  # ONE run, ONE report, ONE output_dir: the FULL products live in this run's
  # artifacts/ and are sealed by the same run manifest (no sibling folders).
  expect_null(result$extension_output_dir)
  expect_false(dir.exists(paste0(result$output_dir, "-full-report")))
  expect_false(dir.exists(paste0(result$output_dir, "-kegg-maps")))
  expect_true(dir.exists(file.path(result$output_dir, "artifacts")))
  metadata <- utils::read.delim(file.path(result$output_dir, "lisa_pipeline_metadata.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(metadata$value[metadata$key == "report_mode"], "full")
  cover <- paste(readLines(file.path(result$output_dir, "report_index.html"), warn = FALSE), collapse = "\n")
  expect_match(cover, "Full LISA report", fixed = TRUE)
  single_page <- paste(readLines(file.path(result$output_dir, "report_pages", "single_de.html"),
    warn = FALSE), collapse = "\n")
  expect_gt(lengths(regmatches(single_page, gregexpr("<h2>FULL figures</h2>", single_page, fixed = TRUE))), 0L)
  post_writer <- utils::read.delim(file.path(result$output_dir, "post_lisa_status.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(sum(post_writer$stage == "root_html_report"), 1L)
  extension_ledger <- file.path(result$output_dir, "artifacts", "code_identity.tsv")
  expect_true(file.exists(extension_ledger))
  expect_true("artifacts/code_identity.tsv" %in% run_manifest$path)
  expect_identical(
    run_manifest$sha256[match("artifacts/code_identity.tsv", run_manifest$path)],
    lisaR:::lisa_sha256_file(extension_ledger)
  )
  expect_true(any(startsWith(run_manifest$path, "artifacts/")))
  code_entries <- utils::read.delim(
    code_ledger, check.names = FALSE, stringsAsFactors = FALSE
  )
  renderer_sha256 <- unique(code_entries$sha256[
    code_entries$relative_path == file.path("scripts", "reproduce_lisa_figure.R")
  ])
  expect_length(renderer_sha256, 1L)
  recipe_files <- c(
    list.files(
      result$output_dir, pattern = "_recipe[.]R$", recursive = TRUE,
      full.names = TRUE
    ),
    list.files(
      file.path(result$output_dir, "artifacts"), pattern = "_recipe[.]R$",
      recursive = TRUE, full.names = TRUE
    )
  )
  recipe_files <- unique(normalizePath(recipe_files, winslash = "/", mustWork = TRUE))
  expect_gt(length(recipe_files), 0L)
  nes_renderer_sha256 <- unique(code_entries$sha256[
    code_entries$relative_path == file.path("scripts", "reproduce_lisa_category_nes.R")
  ])
  expect_length(nes_renderer_sha256, 1L)
  # A real run emits THREE recipe classes, not two. The two generic classes are
  # verified byte copies of an installed renderer; leading-edge heatmaps use a
  # self-contained native recipe on purpose (see test-native-heatmap-recipe.R),
  # so its provenance is established against its generating helper instead of
  # against a renderer hash.
  is_nes_recipe <- grepl("/category_nes/", recipe_files, fixed = TRUE)
  is_native_heatmap_recipe <- grepl(
    "_leading_edge_gene_heatmap_recipe[.]R$", basename(recipe_files)
  )
  expect_true(any(is_nes_recipe))
  expect_true(any(is_native_heatmap_recipe))
  expect_true(any(!is_nes_recipe & !is_native_heatmap_recipe))
  # The three classes partition the emitted recipes: no file is both.
  expect_false(any(is_nes_recipe & is_native_heatmap_recipe))

  # Classes 1 and 2 - verified copies of the installed renderers.
  copied_recipes <- recipe_files[!is_native_heatmap_recipe]
  expected_recipe_sha256 <- ifelse(
    is_nes_recipe[!is_native_heatmap_recipe],
    nes_renderer_sha256[[1L]], renderer_sha256[[1L]]
  )
  expect_identical(unname(vapply(
    copied_recipes, lisaR:::lisa_sha256_file, character(1)
  )), unname(expected_recipe_sha256))

  # Class 3 - native recipes. Provenance comes from the generating helper, whose
  # own identity is proved from THIS run's actual ledger, not from a self-report.
  native_recipes <- recipe_files[is_native_heatmap_recipe]
  native_helper_relative <- file.path(
    "scripts", "lisa_leading_edge_heatmap_native.R"
  )
  # Leading-edge heatmaps are a FULL post script, so their executable closure is
  # recorded in the FULL products ledger (artifacts/), not in the main ledger.
  # That is where the native helper's provenance actually lives.
  extension_code_entries <- utils::read.delim(
    extension_ledger, check.names = FALSE, stringsAsFactors = FALSE
  )
  heatmap_entrypoint <- "build_single_de_leading_edge_gene_heatmaps.R"
  native_ledger_rows <- extension_code_entries[
    extension_code_entries$relative_path == native_helper_relative &
      extension_code_entries$entrypoint == heatmap_entrypoint, , drop = FALSE
  ]
  # Exactly one entry, bound to its own builder as entrypoint.
  expect_equal(nrow(native_ledger_rows), 1L)
  native_helper_sha256 <- unique(native_ledger_rows$sha256)
  expect_length(native_helper_sha256, 1L)
  # The main run ledger must NOT claim the native helper: the two ledgers have
  # distinct responsibilities and this gate previously conflated them.
  expect_false(native_helper_relative %in% code_entries$relative_path)
  native_helper_path <- system.file(
    "scripts", "lisa_leading_edge_heatmap_native.R", package = "lisaR"
  )
  expect_true(nzchar(native_helper_path))
  # The ledger entry is the installed helper that actually ran, byte for byte.
  expect_identical(
    native_helper_sha256[[1L]], lisaR:::lisa_sha256_file(native_helper_path)
  )
  # The helper is a registered closure of its builder, so the pre-launch code
  # ledger covers it; it is deliberately not a generic-renderer dependency.
  expect_true(native_helper_relative %in% file.path(
    "scripts",
    lisaR:::lisa_post_script_dependencies()[[heatmap_entrypoint]]
  ))
  expect_false("reproduce_lisa_figure.R" %in%
    lisaR:::lisa_post_script_dependencies()[[heatmap_entrypoint]])
  # Native recipes must NOT be copies of either generic renderer - that is the
  # exact confusion this gate previously made.
  native_hashes <- unname(vapply(
    native_recipes, lisaR:::lisa_sha256_file, character(1)
  ))
  expect_false(any(native_hashes %in%
    c(renderer_sha256[[1L]], nes_renderer_sha256[[1L]])))

  # Deterministic provenance: every function body embedded in an emitted recipe
  # is byte-identical to the same function deparsed from the ledgered helper.
  helper_env <- new.env(parent = globalenv())
  sys.source(native_helper_path, envir = helper_env)
  native_functions <- c(
    "lisa_heatmap_num", "lisa_heatmap_fdr", "lisa_heatmap_subtitle",
    "lisa_heatmap_device_text", "lisa_heatmap_plot", "lisa_heatmap_size",
    "lisa_heatmap_save_native", "lisa_heatmap_read_bound_source",
    "lisa_heatmap_render_saved"
  )
  expected_definitions <- vapply(native_functions, function(name) {
    paste0(name, " <- ", paste(deparse(
      get(name, envir = helper_env, inherits = FALSE), width.cutoff = 500L
    ), collapse = "\n"))
  }, character(1))
  for (recipe in native_recipes) {
    code <- paste(readLines(recipe, warn = FALSE), collapse = "\n")
    expect_match(code, "Native leading-edge heatmap recipe", fixed = TRUE)
    # Self-contained: no source() call and no lisaR dependency. The negative
    # lookbehind keeps lisa_heatmap_read_bound_source() from matching.
    expect_false(grepl("(^|[^A-Za-z0-9._])source\\(", code))
    expect_false(grepl("lisaR:::", code, fixed = TRUE))
    expect_false(grepl("library(lisaR)", code, fixed = TRUE))
    # Bound to its own saved source and refuses a substituted matrix.
    expect_match(code, "LISA-HEATMAP-RECIPE-002", fixed = TRUE)
    expect_match(code, "lisa_heatmap_saved_context <- ", fixed = TRUE)
    for (definition in expected_definitions) {
      expect_true(grepl(definition, code, fixed = TRUE))
    }
  }

  # Tamper coverage: a single altered byte in a native recipe breaks both its
  # hash and its helper-derived provenance, so this gate cannot pass silently.
  tampered_recipe <- file.path(tempdir(), "tampered_native_heatmap_recipe.R")
  on.exit(unlink(tampered_recipe), add = TRUE)
  original_lines <- readLines(native_recipes[[1L]], warn = FALSE)
  writeLines(sub(
    "lisa_heatmap_num <- ", "lisa_heatmap_num <- ## tampered\n",
    paste(original_lines, collapse = "\n"), fixed = TRUE
  ), tampered_recipe)
  expect_false(identical(
    lisaR:::lisa_sha256_file(tampered_recipe),
    lisaR:::lisa_sha256_file(native_recipes[[1L]])
  ))
  tampered_code <- paste(readLines(tampered_recipe, warn = FALSE), collapse = "\n")
  expect_false(grepl(
    expected_definitions[["lisa_heatmap_num"]], tampered_code, fixed = TRUE
  ))
  member_root <- file.path(result$output_dir, "outputs", "single_de", "response_a",
    "collection_GOBP-C2", "plots")
  member_png <- list.files(member_root, pattern = "[.]png$", recursive = TRUE, full.names = TRUE)
  member_source <- list.files(member_root, pattern = "_source[.]tsv$", recursive = TRUE, full.names = TRUE)
  member_recipe <- list.files(member_root, pattern = "_recipe[.]R$", recursive = TRUE, full.names = TRUE)
  expect_false(any(grepl("GSEA_category_pathways", member_png, fixed = TRUE)))
  expect_false(any(grepl("GSEA_category_pathways", member_source, fixed = TRUE)))
  # The default member plots and their recipes are replaced by evidence.
  # Their optional extension and the new evidence recipes are checked below.
  expect_false(any(grepl("GSEA_category_pathways", member_recipe, fixed = TRUE)))
  extension_files <- list.files(file.path(result$output_dir, "artifacts"),
    recursive = TRUE, full.names = FALSE)
  # GeneCards are a requested FULL product, so the extension carries them.
  # Each card is a complete triplet: figure, its saved source and its recipe.
  expect_true(any(grepl("gene_cards", extension_files, fixed = TRUE)))
  gene_card_files <- extension_files[grepl("gene_cards", extension_files, fixed = TRUE)]
  card_suffixes <- c(
    png = "_category_gene_card[.]png$",
    source = "_category_gene_card_source[.]tsv$",
    recipe = "_category_gene_card_recipe[.]R$"
  )
  card_stems <- lapply(card_suffixes, function(pattern) {
    sub(pattern, "", gene_card_files[grepl(pattern, gene_card_files)])
  })
  expect_gt(length(card_stems$png), 0L)
  expect_setequal(card_stems$source, card_stems$png)
  expect_setequal(card_stems$recipe, card_stems$png)
  # Every gene_cards file belongs to one of the three triplet roles.
  expect_length(gene_card_files, 3L * length(card_stems$png))
  expect_true(any(grepl("member_gene_sets", extension_files, fixed = TRUE)))
  evidence_root <- file.path(result$output_dir, "report_pages", "evidence")
  evidence_pages <- list.files(evidence_root, pattern = "^index[.]html$",
    recursive = TRUE, full.names = TRUE)
  expect_length(evidence_pages, 8L)
  navigation_root <- file.path(result$output_dir, "report_pages", "category_navigation")
  navigation_pages <- list.files(navigation_root,
    pattern = "^index[.]html$", recursive = TRUE, full.names = TRUE)
  # Analyses/Contrasts navigation parity: category evidence exists only per
  # analysis (2 x 4 collections = 8), but navigation also covers the contrast
  # scope (1 x 4 collections = 4), for 12 pages. The extra 4 are the contrast
  # side of the parity, not duplicates of the analysis pages.
  navigation_scopes <- dirname(dirname(
    list.files(navigation_root, pattern = "^index[.]html$", recursive = TRUE)
  ))
  evidence_scopes <- dirname(dirname(
    list.files(evidence_root, pattern = "^index[.]html$", recursive = TRUE)
  ))
  quick_start_collections <- c("GOBP-C2", "GOCC", "GOMF", "PATHWAYS")
  expect_setequal(unique(evidence_scopes), c("response_a", "response_b"))
  expect_setequal(unique(navigation_scopes), c(
    "response_a", "response_b", "response_a_vs_b_response_profiles"
  ))
  expect_length(navigation_pages, 12L)
  expect_length(navigation_pages,
    length(unique(navigation_scopes)) * length(quick_start_collections))
  # Each scope is navigated once per collection - no scope is short-changed.
  for (scope in unique(navigation_scopes)) {
    expect_setequal(
      basename(dirname(list.files(file.path(navigation_root, scope),
        pattern = "^index[.]html$", recursive = TRUE))),
      quick_start_collections
    )
  }
  # The contrast navigation scope carries no category-evidence page of its own.
  expect_false("response_a_vs_b_response_profiles" %in% evidence_scopes)
  expect_true(file.exists(file.path(result$output_dir, "report_pages", "gene_evidence", "index.html")))
  expect_true(any(grepl("GSEA_lollipop_direction_stats.png", member_png, fixed = TRUE)))
  expect_true(any(grepl("GSEA_pathway_dotplot.png", member_png, fixed = TRUE)))
  expect_false(any(grepl("ORA", member_png, fixed = TRUE)))
  for (page in evidence_pages) {
    expect_true(file.exists(file.path(dirname(page), "tables", "leading_edges.tsv")))
    expect_true(file.exists(file.path(dirname(page), "tables", "categories.tsv")))
    expect_true(file.exists(file.path(dirname(page), "reproduce_category_evidence.R")))
  }
  post_status <- utils::read.delim(file.path(result$output_dir, "post_lisa_status.tsv"))
  evidence_status <- post_status[post_status$stage == "single_de_category_evidence", ]
  expect_equal(nrow(evidence_status), 8L)
  expect_true(all(evidence_status$status == "completed"))
  opt_in <- lisaR:::lisa_extension_read_selection(list(mode = "full",
    report = list(legacy_gene_products = TRUE)))
  dummy_products <- data.frame(product = c("gene_cards", "volcano"))
  expect_identical(lisaR:::lisa_extension_filter(dummy_products, opt_in), dummy_products)
  expect_true(any(grepl("volcano", extension_files, fixed = TRUE)))
  expect_true(any(grepl("heatmap", extension_files, fixed = TRUE)))
  expect_true(any(grepl("contrast_profile", extension_files, fixed = TRUE)))
  heatmap_files <- extension_files[grepl("/heatmap/", extension_files, fixed = TRUE)]
  expect_false(any(grepl("gene_card|volcano", heatmap_files)))
  contrast_profiles <- extension_files[
    grepl("/contrast_profile/", extension_files, fixed = TRUE) &
      grepl("[.]png$", extension_files)
  ]
  contrast_status <- utils::read.delim(
    file.path(result$output_dir, "contrast_status.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  contrast_tables <- unlist(lapply(contrast_status$output_dir, function(output_dir) {
    list.files(
      file.path(output_dir, "lisa_tables"),
      pattern = "_GSEA_category_contrast[.]tsv$",
      full.names = TRUE
    )
  }), use.names = FALSE)
  expect_length(contrast_tables, nrow(contrast_status))
  contrast_rows <- do.call(rbind, lapply(contrast_tables, function(path) {
    utils::read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  # Full reports retain every category in map order, including visually blank
  # categories for which neither side crosses the configured GSEA FDR cutoff.
  expect_length(contrast_profiles, nrow(contrast_rows))
  contrast_profile_sources <- list.files(
    file.path(result$output_dir, "artifacts"),
    pattern = "_source[.]tsv$", recursive = TRUE, full.names = TRUE
  )
  contrast_profile_sources <- contrast_profile_sources[
    grepl("/contrast_profile/", contrast_profile_sources, fixed = TRUE)
  ]
  expect_length(contrast_profile_sources, nrow(contrast_rows))
  profile_source <- lisaR:::read_lisa_tsv(contrast_profile_sources[[1L]])
  expect_true(all(c(
    "plot_set", "annotation_variant", "group_by_supracategory",
    "figure_width", "figure_height", "figure_dpi"
  ) %in% names(profile_source)))
  expect_true(all(profile_source$plot_set == "contrast_profile"))
  expect_true(all(profile_source$annotation_variant == "plain"))
  expect_true(all(!profile_source$group_by_supracategory))
  expect_true(all(profile_source$figure_width == 7))
  expect_true(all(profile_source$figure_height == 4.2))
  expect_true(all(profile_source$figure_dpi == 300))
})

test_that("custom category macrogroups receive a deterministic fallback palette", {
  category_map <- data.frame(
    category_id = c("CUSTOM_A", "CUSTOM_B"),
    macrogroup_id = c("CUSTOM_GROUP", "CUSTOM_GROUP"),
    stringsAsFactors = FALSE
  )
  colours <- lisaR:::build_lisa_palette(category_map, "lisa_default")
  expect_length(colours, 2L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colours)))
})
