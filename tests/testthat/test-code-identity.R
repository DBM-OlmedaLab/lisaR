# Generic recipe handoff fixtures use volcano: heatmaps now emit a native recipe.
code_identity_fixture <- function(include_script = TRUE, include_helper = TRUE,
                                  include_renderer = TRUE) {
  root <- tempfile("lisa-code-identity-")
  script_dir <- file.path(root, "inst", "scripts")
  dir.create(script_dir, recursive = TRUE)
  writeLines(
    c(
      "Package: lisaR",
      "Version: 0.0.0",
      "Title: Explicit internal code fixture",
      "Description: Test-only package tree for executable identity checks.",
      "License: MIT"
    ),
    file.path(root, "DESCRIPTION")
  )
  if (isTRUE(include_script)) {
    writeLines(
      c(
        "args <- commandArgs(trailingOnly = TRUE)",
        "project_flag <- match('--project-dir', args)",
        "renderer_flag <- match('--lisa-internal-renderer-sha256', args)",
        "if (!is.na(renderer_flag)) writeLines(args[[renderer_flag + 1L]], file.path(args[[project_flag + 1L]], 'renderer-sha256.txt'))",
        "writeLines('executed', file.path(args[[project_flag + 1L]], 'child-marker.txt'))"
      ),
      file.path(script_dir, "build_single_de_category_volcano_overlays.R")
    )
  }
  if (isTRUE(include_helper)) {
    writeLines("identity_helper <- TRUE", file.path(script_dir, "lisa_plot_metadata.R"))
  }
  if (isTRUE(include_renderer)) {
    writeLines("message('renderer')", file.path(script_dir, "reproduce_lisa_figure.R"))
  }
  lisaR:::lisa_resolve_package_dir(.test_package_dir = root)
}

code_identity_args <- function(run_root) {
  c(
    "--project-dir", run_root,
    "--analysis-id", "analysis-a",
    "--universe", "GOBP-C2"
  )
}

test_that("registered script closures cover every executable postprocessor", {
  closures <- lisaR:::lisa_post_script_dependencies()
  expect_setequal(
    names(closures),
    c(
      "build_LISA_report.R",
      "build_category_evidence.R",
      "build_category_navigation.R",
      "build_gene_evidence.R",
      "build_contrast_evidence.R",
      "build_contrast_navigation.R",
      "build_contrast_gene_category_network.R",
      "build_contrast_gene_level_product.R",
      "build_contrast_kegg_pathway_painter.R",
      "build_contrast_macrogroup_heatmaps.R",
      "build_single_de_category_gene_cards.R",
      "build_single_de_category_volcano_overlays.R",
      "build_single_de_enrichmentmap.R",
      "build_single_de_kegg_pathway_painter.R",
      "build_single_de_leading_edge_gene_heatmaps.R",
      "build_single_de_recurrent_gene_screen.R"
    )
  )
  expect_setequal(
    unique(unlist(closures, use.names = FALSE)),
    c(
      "reproduce_lisa_figure.R", "reproduce_lisa_category_nes.R", "kegg_snapshot_helpers.R",
      "build_contrast_kegg_map_layer.R", "lisa_plot_metadata.R", "lisa_leading_edge_heatmap_native.R"
    )
  )

  script_dir <- system.file("scripts", package = "lisaR")
  renderer_users <- names(closures)[vapply(names(closures), function(script) {
    code <- readLines(file.path(script_dir, script), warn = FALSE)
    any(grepl("reproduce_lisa_figure.R", code, fixed = TRUE))
  }, logical(1))]
  expect_true(length(renderer_users) > 0L)
  expect_true(all(vapply(
    renderer_users,
    function(script) "reproduce_lisa_figure.R" %in% closures[[script]],
    logical(1)
  )))
})

test_that("the report closure binds its installed NES reproduction recipe", {
  root <- code_identity_fixture()
  on.exit(unlink(as.character(root), recursive = TRUE), add = TRUE)
  scripts <- attr(root, "lisaR.script_dir", exact = TRUE)
  writeLines("message('report fixture')", file.path(scripts, "build_LISA_report.R"))
  recipe <- file.path(scripts, "reproduce_lisa_category_nes.R")
  writeLines("message('NES recipe fixture')", recipe)
  assets <- file.path(dirname(scripts), "report_assets")
  dir.create(assets)
  names <- c("lisa_shell.css", "lisa_shell.js", "LISA_logo_C_compact_icon_muted_red_S.svg")
  expect_true(all(file.copy(system.file("report_assets", names, package = "lisaR"), assets)))
  identity <- lisaR:::lisa_post_script_identity(root, "build_LISA_report.R")
  expect_true(file.path("scripts", basename(recipe)) %in% identity$files$relative_path)
  expect_identical(lisaR:::lisa_verify_post_script_identity(identity)$files, identity$files)
  writeLines("message('changed NES recipe')", recipe)
  changed <- tryCatch(lisaR:::lisa_verify_post_script_identity(identity), error = function(e) e)
  expect_s3_class(changed, "lisa_code_identity_error")
  expect_identical(changed$code, "LISA-CODE-IDENTITY-002")
  expect_true(file.remove(recipe))
  missing <- tryCatch(lisaR:::lisa_post_script_identity(root, "build_LISA_report.R"), error = function(e) e)
  expect_s3_class(missing, "lisa_code_identity_error")
  expect_identical(missing$code, "LISA-CODE-IDENTITY-001")
})

test_that("an absent required script or helper aborts before child execution", {
  for (missing in c("script", "helper", "renderer")) {
    code_root <- code_identity_fixture(
      include_script = !identical(missing, "script"),
      include_helper = !identical(missing, "helper"),
      include_renderer = !identical(missing, "renderer")
    )
    run_root <- lisaR:::lisa_run_root(tempfile(paste0("lisa-code-missing-", missing, "-")))
    error <- tryCatch(
      lisaR:::lisa_run_post_script(
        code_root,
        "build_single_de_category_volcano_overlays.R",
        code_identity_args(run_root),
        trusted_run_root = run_root
      ),
      error = identity
    )
    expect_s3_class(error, "lisa_code_identity_error")
    expect_identical(error$code, "LISA-CODE-IDENTITY-001")
    expect_false(file.exists(file.path(run_root, "child-marker.txt")))
  }
})

test_that("symlinked entrypoints, helpers and renderers fail before hashing", {
  targets <- c(
    script = "build_single_de_category_volcano_overlays.R",
    helper = "lisa_plot_metadata.R",
    renderer = "reproduce_lisa_figure.R"
  )
  for (role in names(targets)) {
    code_root <- code_identity_fixture()
    script_dir <- attr(code_root, "lisaR.script_dir", exact = TRUE)
    candidate <- file.path(script_dir, targets[[role]])
    external <- tempfile(paste0("lisa-external-", role, "-"), fileext = ".R")
    writeLines(paste("external", role), external)
    expect_true(file.remove(candidate))
    linked <- suppressWarnings(file.symlink(external, candidate))
    if (!isTRUE(linked)) skip("symbolic links are unavailable on this platform")

    error <- tryCatch(
      lisaR:::lisa_post_script_identity(
        code_root, "build_single_de_category_volcano_overlays.R"
      ),
      error = identity
    )
    expect_s3_class(error, "lisa_code_identity_error")
    expect_identical(error$code, "LISA-CODE-IDENTITY-001")
    expect_identical(readLines(external, warn = FALSE), paste("external", role))
  }
})

test_that("a scripts directory replaced by a symlink fails after resolution", {
  code_root <- code_identity_fixture()
  script_dir <- attr(code_root, "lisaR.script_dir", exact = TRUE)
  external <- tempfile("lisa-external-scripts-")
  dir.create(external)
  expect_true(all(file.copy(
    list.files(script_dir, full.names = TRUE), external
  )))
  unlink(script_dir, recursive = TRUE, force = TRUE)
  linked <- suppressWarnings(file.symlink(external, script_dir))
  if (!isTRUE(linked)) skip("symbolic links are unavailable on this platform")

  error <- tryCatch(
    lisaR:::lisa_post_script_path(
      code_root, "reproduce_lisa_figure.R"
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-ROOT-005")
  expect_true(file.exists(file.path(external, "reproduce_lisa_figure.R")))
})

test_that("a helper mutation between inventory and launch aborts before child execution", {
  code_root <- code_identity_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-code-mutated-"))
  ledger <- lisaR:::lisa_write_code_identity_ledger(
    code_root, "build_single_de_category_volcano_overlays.R",
    file.path(run_root, "code_identity.tsv")
  )
  original_verify <- lisaR:::lisa_verify_post_script_identity
  testthat::local_mocked_bindings(
    lisa_verify_post_script_identity = function(identity) {
      helper <- identity$files$absolute_path[
        basename(identity$files$absolute_path) == "lisa_plot_metadata.R"
      ][[1L]]
      writeLines("identity_helper <- FALSE", helper)
      original_verify(identity)
    },
    .package = "lisaR"
  )
  error <- tryCatch(
    lisaR:::lisa_run_post_script(
      code_root,
      "build_single_de_category_volcano_overlays.R",
      code_identity_args(run_root),
      trusted_run_root = run_root,
      code_ledger = ledger
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-IDENTITY-002")
  expect_false(file.exists(file.path(run_root, "child-marker.txt")))
})

test_that("a renderer mutation between inventory and launch aborts before child execution", {
  code_root <- code_identity_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-renderer-mutated-"))
  ledger <- lisaR:::lisa_write_code_identity_ledger(
    code_root, "build_single_de_category_volcano_overlays.R",
    file.path(run_root, "code_identity.tsv")
  )
  original_verify <- lisaR:::lisa_verify_post_script_identity
  testthat::local_mocked_bindings(
    lisa_verify_post_script_identity = function(identity) {
      renderer <- identity$files$absolute_path[
        basename(identity$files$absolute_path) == "reproduce_lisa_figure.R"
      ][[1L]]
      writeLines("message('mutated renderer')", renderer)
      original_verify(identity)
    },
    .package = "lisaR"
  )
  error <- tryCatch(
    lisaR:::lisa_run_post_script(
      code_root,
      "build_single_de_category_volcano_overlays.R",
      code_identity_args(run_root),
      trusted_run_root = run_root,
      code_ledger = ledger
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-IDENTITY-002")
  expect_false(file.exists(file.path(run_root, "child-marker.txt")))
})

test_that("an unchanged explicit test closure executes and reports its hashes", {
  code_root <- code_identity_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-code-unchanged-"))
  testthat::local_mocked_bindings(
    lisa_assert_child_runtime_identity = function(...) invisible(TRUE),
    .package = "lisaR"
  )
  ledger <- lisaR:::lisa_write_code_identity_ledger(
    code_root, "build_single_de_category_volcano_overlays.R",
    file.path(run_root, "code_identity.tsv")
  )
  result <- lisaR:::lisa_run_post_script(
    code_root,
    "build_single_de_category_volcano_overlays.R",
    code_identity_args(run_root),
    trusted_run_root = run_root,
    code_ledger = ledger
  )
  expect_identical(result$status, "completed")
  expect_true(file.exists(file.path(run_root, "child-marker.txt")))
  renderer_row <- ledger$entries[
    ledger$entries$entrypoint == "build_single_de_category_volcano_overlays.R" &
      ledger$entries$relative_path == file.path("scripts", "reproduce_lisa_figure.R"),
    , drop = FALSE
  ]
  expect_equal(nrow(renderer_row), 1L)
  expect_identical(
    readLines(file.path(run_root, "renderer-sha256.txt"), warn = FALSE),
    renderer_row$sha256
  )
  expect_setequal(
    result$code_identity$relative_path,
    file.path("scripts", c(
      "build_single_de_category_volcano_overlays.R", "lisa_plot_metadata.R",
      "reproduce_lisa_figure.R"
    ))
  )
  expect_true(all(grepl("^[0-9a-f]{64}$", result$code_identity$sha256)))
  expect_true(file.exists(ledger$path))
  expect_false(any(grepl(tempdir(), ledger$entries$package_root, fixed = TRUE)))
  expect_true(all(ledger$entries$package_version == "0.0.0"))
})

test_that("post-processing pins the selected package library in a vanilla child", {
  code_root <- code_identity_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-child-library-"))
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    lisa_assert_child_runtime_identity = function(...) invisible(TRUE),
    lisa_run_subprocess = function(executable, args, required, stage,
                                   run_root, env) {
      captured$executable <- executable
      captured$args <- args
      captured$env <- env
      list(exit_code = 0L, diagnostics = "")
    },
    .package = "lisaR"
  )

  result <- lisaR:::lisa_run_post_script(
    code_root,
    "build_single_de_category_volcano_overlays.R",
    code_identity_args(run_root),
    trusted_run_root = run_root
  )

  expect_identical(result$status, "completed")
  expect_identical(captured$args[[1L]], "--vanilla")
  expect_identical(captured$args[[2L]], file.path(
    attr(code_root, "lisaR.script_dir", exact = TRUE),
    "build_single_de_category_volcano_overlays.R"
  ))
  libraries <- strsplit(
    unname(captured$env[["R_LIBS"]]), .Platform$path.sep, fixed = TRUE
  )[[1L]]
  expect_identical(
    libraries[[1L]],
    normalizePath(dirname(code_root), winslash = "/", mustWork = TRUE)
  )
})

test_that("a rebuilt ledger records a changed helper digest", {
  code_root <- code_identity_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-code-ledger-change-"))
  first <- lisaR:::lisa_write_code_identity_ledger(
    code_root, "build_single_de_category_volcano_overlays.R",
    file.path(run_root, "code_identity-before.tsv")
  )
  helper <- lisaR:::lisa_post_script_path(code_root, "lisa_plot_metadata.R")
  writeLines("identity_helper <- 'changed'", helper)
  second <- lisaR:::lisa_write_code_identity_ledger(
    code_root, "build_single_de_category_volcano_overlays.R",
    file.path(run_root, "code_identity-after.tsv")
  )
  first_hash <- first$entries$sha256[
    basename(first$entries$relative_path) == "lisa_plot_metadata.R"
  ]
  second_hash <- second$entries$sha256[
    basename(second$entries$relative_path) == "lisa_plot_metadata.R"
  ]
  expect_false(identical(first_hash, second_hash))
  expect_false(identical(first$sha256, second$sha256))
})

test_that("extension recipe copy is bound to its serial code ledger", {
  code_root <- code_identity_fixture()
  script_dir <- attr(code_root, "lisaR.script_dir", exact = TRUE)
  recipe <- file.path(script_dir, "reproduce_lisa_figure.R")
  writeLines("message('verified recipe')", recipe)
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-code-recipe-"))
  ledger <- lisaR:::lisa_write_code_identity_ledger(
    code_root, character(), file.path(run_root, "code_identity.tsv"),
    executable_recipes = "reproduce_lisa_figure.R"
  )
  writeLines("message('mutated recipe')", recipe)
  target <- file.path(run_root, "figure")
  error <- tryCatch(
    lisaR:::lisa_extension_write_recipe(target, code_root, ledger),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-LEDGER-003")
  expect_false(file.exists(paste0(target, "_recipe.R")))
})

test_that("a corrupt executable copy is removed", {
  source <- tempfile("lisa-code-copy-source-", fileext = ".R")
  target <- tempfile("lisa-code-copy-target-", fileext = ".R")
  writeLines("message('source')", source)
  expected <- lisaR:::lisa_sha256_file(source)
  corrupt_copy <- function(from, to, overwrite = FALSE) {
    writeLines("message('corrupt')", to)
    TRUE
  }
  error <- tryCatch(
    lisaR:::lisa_copy_verified_code_file(
      source, target, expected, .copy = corrupt_copy
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-COPY-001")
  expect_false(file.exists(target))
})

test_that("a source mutation during executable copy removes the target", {
  source <- tempfile("lisa-code-race-source-", fileext = ".R")
  target <- tempfile("lisa-code-race-target-", fileext = ".R")
  writeLines("message('source')", source)
  expected <- lisaR:::lisa_sha256_file(source)
  mutating_copy <- function(from, to, overwrite = FALSE) {
    copied <- file.copy(from, to, overwrite = overwrite)
    writeLines("message('mutated after copy')", from)
    copied
  }
  error <- tryCatch(
    lisaR:::lisa_copy_verified_code_file(
      source, target, expected, .copy = mutating_copy
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-COPY-001")
  expect_false(file.exists(target))
})

test_that("the renderer handoff fails closed and standalone mode snapshots once", {
  source <- tempfile("lisa-renderer-source-", fileext = ".R")
  writeLines("message('verified renderer')", source)
  expected <- lisaR:::lisa_sha256_file(source)

  malformed_target <- tempfile("lisa-renderer-malformed-", fileext = ".R")
  malformed <- tryCatch(
    lisaR:::lisa_copy_verified_figure_recipe(
      source, malformed_target, "not-a-sha256"
    ),
    error = identity
  )
  expect_s3_class(malformed, "lisa_code_identity_error")
  expect_identical(malformed$code, "LISA-CODE-COPY-001")
  expect_false(file.exists(malformed_target))

  writeLines("message('mutated renderer')", source)
  mutated_target <- tempfile("lisa-renderer-mutated-target-", fileext = ".R")
  mutated <- tryCatch(
    lisaR:::lisa_copy_verified_figure_recipe(
      source, mutated_target, expected
    ),
    error = identity
  )
  expect_s3_class(mutated, "lisa_code_identity_error")
  expect_identical(mutated$code, "LISA-CODE-COPY-001")
  expect_false(file.exists(mutated_target))

  standalone_target <- tempfile("lisa-renderer-standalone-", fileext = ".R")
  expect_no_error(lisaR:::lisa_copy_verified_figure_recipe(
    source, standalone_target
  ))
  expect_identical(
    lisaR:::lisa_sha256_file(standalone_target),
    lisaR:::lisa_sha256_file(source)
  )
})

test_that("source-only launches reject a different installed child runtime", {
  code_root <- code_identity_fixture()
  source_environment <- new.env(parent = emptyenv())
  child_namespace <- new.env(parent = emptyenv())
  source_environment$verified_copy <- function(path) paste0("source:", path)
  child_namespace$verified_copy <- source_environment$verified_copy
  withr::local_options(list(
    lisaR.source_runtime_functions = "verified_copy"
  ))
  expect_no_error(lisaR:::lisa_assert_child_runtime_identity(
    code_root,
    .namespace = child_namespace,
    .source_environment = source_environment
  ))

  child_namespace$verified_copy <- function(path) paste0("other:", path)
  error <- tryCatch(
    lisaR:::lisa_assert_child_runtime_identity(
      code_root,
      .namespace = child_namespace,
      .source_environment = source_environment
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_code_identity_error")
  expect_identical(error$code, "LISA-CODE-RUNTIME-001")
  expect_match(conditionMessage(error), "verified_copy", fixed = TRUE)
})

test_that("all renderer builders use the private verified-copy boundary", {
  script_dir <- system.file("scripts", package = "lisaR")
  if (!nzchar(script_dir) || !dir.exists(script_dir)) {
    script_dir <- testthat::test_path("..", "..", "inst", "scripts")
  }
  renderer_scripts <- names(Filter(
    function(dependencies) "reproduce_lisa_figure.R" %in% dependencies,
    lisaR:::lisa_post_script_dependencies()
  ))
  # This exact set—not a stale aggregate count—is the generic renderer
  # boundary. Leading-edge heatmaps intentionally use their native helper and
  # are tested separately in test-native-heatmap-recipe.R.
  expect_setequal(renderer_scripts, c(
    "build_LISA_report.R",
    "build_contrast_gene_level_product.R",
    "build_contrast_kegg_pathway_painter.R",
    "build_single_de_category_gene_cards.R",
    "build_single_de_category_volcano_overlays.R",
    "build_single_de_kegg_pathway_painter.R"
  ))
  expect_false("build_single_de_leading_edge_gene_heatmaps.R" %in% renderer_scripts)
  for (script_name in renderer_scripts) {
    code <- paste(
      readLines(file.path(script_dir, script_name), warn = FALSE),
      collapse = "\n"
    )
    expect_match(
      code, "lisaR:::lisa_copy_verified_figure_recipe", fixed = TRUE,
      info = script_name
    )
    expect_false(
      grepl("file[.]copy[[:space:]]*[(][[:space:]]*renderer", code),
      info = script_name
    )
    expect_false(
      grepl("lisa_sha256_file[[:space:]]*[(][[:space:]]*renderer", code),
      info = script_name
    )
    expect_match(
      code, "lisa_internal_renderer_sha256", fixed = TRUE,
      info = script_name
    )
  }
  report_script <- paste(
    readLines(file.path(script_dir, "build_LISA_report.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    report_script,
    "allow mutable run data to select code copied into the final report",
    fixed = TRUE
  )
})
