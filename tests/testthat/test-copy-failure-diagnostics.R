test_that("guarded copy folds unique write warnings into one atomic failure", {
  base <- tempfile("lisa-copy-diagnostic-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  root <- lisaR:::lisa_run_root(file.path(base, "run"))
  source <- file.path(base, "source.txt")
  target <- file.path(root, "nested", "target.txt")
  writeLines("copy source", source, useBytes = TRUE)

  testthat::local_mocked_bindings(
    lisa_file_copy = function(...) {
      warning("write error during file append")
      warning("write error during file append")
      FALSE
    },
    .package = "lisaR"
  )

  observed_warnings <- character()
  error <- withCallingHandlers(
    tryCatch(
      lisaR:::lisa_guarded_copy(source, target, run_root = root),
      error = identity
    ),
    warning = function(warning) {
      observed_warnings <<- c(
        observed_warnings, conditionMessage(warning)
      )
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(error, "error")
  expect_match(
    conditionMessage(error),
    "Guarded file copy failed:.*filesystem diagnostic: write error during file append"
  )
  matches <- gregexpr(
    "write error during file append", conditionMessage(error), fixed = TRUE
  )[[1L]]
  expect_identical(sum(matches > 0L), 1L)
  expect_length(observed_warnings, 0L)
  expect_false(lisaR:::lisa_path_entry_exists(target))
  expect_false(any(grepl(
    "^[.]lisa-copy-", list.files(dirname(target), all.files = TRUE)
  )))
})

test_that("guarded copy preserves a warning when the copy succeeds", {
  base <- tempfile("lisa-copy-success-warning-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  root <- lisaR:::lisa_run_root(file.path(base, "run"))
  source <- file.path(base, "source.txt")
  target <- file.path(root, "target.txt")
  writeLines("copy source", source, useBytes = TRUE)

  testthat::local_mocked_bindings(
    lisa_file_copy = function(from, to) {
      copied <- file.copy(
        from, to, overwrite = TRUE, copy.mode = FALSE, copy.date = FALSE
      )
      warning("nonfatal copy diagnostic")
      copied
    },
    .package = "lisaR"
  )

  observed_warnings <- character()
  withCallingHandlers(
    lisaR:::lisa_guarded_copy(source, target, run_root = root),
    warning = function(warning) {
      observed_warnings <<- c(
        observed_warnings, conditionMessage(warning)
      )
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(observed_warnings, "nonfatal copy diagnostic")
  expect_identical(
    readLines(target, warn = FALSE), readLines(source, warn = FALSE)
  )
})
