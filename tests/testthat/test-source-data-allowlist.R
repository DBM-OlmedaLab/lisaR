test_that("source_data requires a non-empty explicit allowlist before copying", {
  root <- tempfile("lisa-source-allowlist-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  first <- file.path(root, "first.tsv")
  second <- file.path(root, "second.tsv")
  writeLines("value\nfirst", first)
  writeLines("value\nsecond", second)
  output <- lisaR:::lisa_run_root(file.path(root, "result"))
  old_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = output)
  on.exit(options(lisaR.run_root = old_root), add = TRUE)

  absent <- list(source_data = list(list(path = first)))
  absent_error <- tryCatch(
    lisaR:::lisa_copy_config_source_data(absent, output, root),
    lisa_source_data_error = identity
  )
  expect_s3_class(absent_error, "lisa_source_data_error")
  expect_s3_class(absent_error, "lisa_error")
  expect_identical(absent_error$code, "LISA-SOURCE-003")
  expect_match(conditionMessage(absent_error), "explicit non-empty", fixed = TRUE)
  expect_false(dir.exists(file.path(output, "source_data")))

  empty <- list(
    allowlisted_source_paths = character(),
    source_data = list(list(path = first))
  )
  empty_error <- tryCatch(
    lisaR:::lisa_copy_config_source_data(empty, output, root),
    lisa_source_data_error = identity
  )
  expect_identical(empty_error$code, "LISA-SOURCE-003")
  expect_false(dir.exists(file.path(output, "source_data")))

  # All entries are authorized before the first copy. A bad second entry must
  # not leave the valid first file behind.
  partial <- list(source_data = list(
    list(path = first, allowlisted_paths = first),
    list(path = second, allowlisted_paths = character())
  ))
  partial_error <- tryCatch(
    lisaR:::lisa_copy_config_source_data(partial, output, root),
    lisa_source_data_error = identity
  )
  expect_identical(partial_error$code, "LISA-SOURCE-003")
  expect_false(dir.exists(file.path(output, "source_data")))
})

test_that("source_data allowlists require an exact resolved file match", {
  root <- tempfile("lisa-source-exact-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  source_dir <- file.path(root, "inputs")
  dir.create(source_dir)
  source <- file.path(source_dir, "source.tsv")
  other <- file.path(source_dir, "other.tsv")
  writeLines("value\nsource", source)
  writeLines("value\nother", other)

  mismatch <- list(
    allowlisted_source_paths = source_dir,
    source_data = list(list(path = source))
  )
  mismatch_error <- tryCatch(
    lisaR:::lisa_prepare_config_source_data(mismatch, root),
    lisa_source_data_error = identity
  )
  expect_s3_class(mismatch_error, "lisa_source_data_error")
  expect_identical(mismatch_error$code, "LISA-SOURCE-001")
  expect_match(conditionMessage(mismatch_error), "exact allowlist match", fixed = TRUE)

  output <- lisaR:::lisa_run_root(file.path(root, "result"))
  old_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = output)
  on.exit(options(lisaR.run_root = old_root), add = TRUE)
  exact <- list(
    allowlisted_source_paths = file.path("inputs", "source.tsv"),
    source_data = list(list(
      path = file.path("inputs", "source.tsv"),
      role = "audit_source", target_subdir = "registered_sources"
    ))
  )
  manifest_path <- lisaR:::lisa_copy_config_source_data(exact, output, root)
  manifest <- read_lisa_tsv(manifest_path)
  copied <- file.path(output, "source_data", "registered_sources", "source.tsv")
  expect_true(file.exists(copied))
  expect_identical(readLines(copied), readLines(source))
  expect_identical(manifest$source_role, "audit_source")
  expect_identical(manifest$source_path, normalizePath(source, winslash = "/"))
  expect_identical(manifest$source_sha256, manifest$sha256)

  row_override <- list(
    allowlisted_source_paths = other,
    source_data = list(list(path = source, allowlisted_paths = source))
  )
  prepared <- lisaR:::lisa_prepare_config_source_data(row_override, root)
  expect_identical(prepared[[1]]$source_path, normalizePath(source, winslash = "/"))
})

test_that("the installed JSON Schema exposes the fail-closed source_data contract", {
  schema_path <- system.file("schema", "lisa-config.schema.json", package = "lisaR")
  schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
  source_schema <- schema$properties$source_data
  global_allowlist <- schema$properties$allowlisted_source_paths
  row_allowlist <- source_schema$items$properties$allowlisted_paths

  expect_identical(source_schema$items$additionalProperties, FALSE)
  expect_setequal(
    names(source_schema$items$properties),
    lisaR:::lisa_config_allowed_keys("source_data")
  )
  expect_identical(row_allowlist$minItems, 1L)
  expect_identical(row_allowlist$items$minLength, 1L)
  expect_true(length(schema$allOf) >= 1L)
  expect_true("then" %in% names(schema$allOf[[1]]))
  global_branch <- schema$allOf[[1]]$then$anyOf[[1]]
  expect_identical(
    global_branch$properties$allowlisted_source_paths$minItems, 1L
  )
})
