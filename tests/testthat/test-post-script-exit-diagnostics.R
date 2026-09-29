post_script_exit_fixture <- function() {
  root <- tempfile("lisa-post-exit-")
  script_dir <- file.path(root, "inst", "scripts")
  dir.create(script_dir, recursive = TRUE)
  writeLines(
    c(
      "Package: lisaR",
      "Version: 0.0.0",
      "Title: Post-script exit diagnostic fixture",
      "Description: Test-only package tree for subprocess diagnostics.",
      "License: MIT"
    ),
    file.path(root, "DESCRIPTION"), useBytes = TRUE
  )
  writeLines(
    "message('fixture entrypoint')",
    file.path(script_dir, "build_single_de_leading_edge_gene_heatmaps.R"),
    useBytes = TRUE
  )
  writeLines(
    "fixture_helper <- TRUE",
    file.path(script_dir, "lisa_plot_metadata.R"), useBytes = TRUE
  )
  writeLines(
    "message('fixture renderer')",
    file.path(script_dir, "lisa_leading_edge_heatmap_native.R"), useBytes = TRUE
  )
  lisaR:::lisa_resolve_package_dir(.test_package_dir = root)
}

post_script_exit_args <- function(run_root) {
  c(
    "--project-dir", run_root,
    "--analysis-id", "analysis-a",
    "--universe", "GOBP-C2"
  )
}

test_that("post-script failures retain a silent nonzero exit code", {
  code_root <- post_script_exit_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-post-exit-run-"))
  testthat::local_mocked_bindings(
    lisa_assert_child_runtime_identity = function(...) invisible(TRUE),
    lisa_run_subprocess = function(...) {
      list(exit_code = 137L, diagnostics = "")
    },
    .package = "lisaR"
  )

  result <- lisaR:::lisa_run_post_script(
    code_root,
    "build_single_de_leading_edge_gene_heatmaps.R",
    post_script_exit_args(run_root),
    trusted_run_root = run_root
  )

  expect_identical(result$status, "failed")
  expect_identical(result$exit_code, 137L)
  expect_identical(result$message, "exit_code=137; no stdout/stderr")
})

test_that("post-script failures retain both exit code and child diagnostics", {
  code_root <- post_script_exit_fixture()
  run_root <- lisaR:::lisa_run_root(tempfile("lisa-post-detail-run-"))
  testthat::local_mocked_bindings(
    lisa_assert_child_runtime_identity = function(...) invisible(TRUE),
    lisa_run_subprocess = function(...) {
      list(exit_code = 7L, diagnostics = "child diagnostic")
    },
    .package = "lisaR"
  )

  result <- lisaR:::lisa_run_post_script(
    code_root,
    "build_single_de_leading_edge_gene_heatmaps.R",
    post_script_exit_args(run_root),
    trusted_run_root = run_root
  )

  expect_identical(result$status, "failed")
  expect_identical(result$exit_code, 7L)
  expect_identical(result$message, "exit_code=7; child diagnostic")
})
