test_that("build exclusions name every operational runtime directory", {
  source_root <- normalizePath(test_path("..", ".."), mustWork = TRUE)
  buildignore_path <- file.path(source_root, ".Rbuildignore")
  # .Rbuildignore is intentionally excluded from source tarballs. The
  # standalone verify_source_tarball.R test validates the built artifact.
  skip_if_not(file.exists(buildignore_path), ".Rbuildignore is source-only; built-tarball boundary validation is performed by verify_source_tarball.R.")
  buildignore <- readLines(buildignore_path, warn = FALSE)
  runtime_names <- c("runs", "staging", "input_staging", "exports", "backups", "audits", "hpc_runs", "smoke-output")

  for (runtime_name in runtime_names) {
    expect_true(
      any(grepl(paste0("\\^", runtime_name, "\\$"), buildignore)),
      info = paste("missing .Rbuildignore entry for", runtime_name)
    )
  }
})
