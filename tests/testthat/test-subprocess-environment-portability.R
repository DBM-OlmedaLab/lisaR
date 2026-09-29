# Regression coverage for the C4 native-platform subprocess defect.
#
# `system2(env = )` is not portable. Per ?system2, "On Windows, 'env' is only
# supported for commands such as 'R' and 'make' which accept environment
# variables on their command line": R appends the NAME=value tokens to the
# child's command line, and `Rscript` takes the first non-option argument as
# the script to execute. Every lisaR post-processing child is launched through
# `lisa_run_subprocess()` with `env = c(R_LIBS = ...)`, so on Windows the child
# ran `R_LIBS=...` instead of the requested script and produced nothing.
#
# These tests pin the replacement contract: the variables reach the child, they
# are restored in the parent afterwards, and the mechanism is the same on every
# platform (so Linux CI is real coverage for it, unlike the path-alias half of
# this repair).

test_that("child environment is exported, then restored exactly", {
  present_name <- "LISAR_PORTABILITY_PRESENT"
  absent_name <- "LISAR_PORTABILITY_ABSENT"
  Sys.setenv(LISAR_PORTABILITY_PRESENT = "original value")
  Sys.unsetenv(absent_name)
  on.exit(Sys.unsetenv(c(present_name, absent_name)), add = TRUE)

  observed <- lisaR:::lisa_with_child_environment(
    c(LISAR_PORTABILITY_PRESENT = "overridden", LISAR_PORTABILITY_ABSENT = "added"),
    c(
      Sys.getenv(present_name, unset = NA_character_),
      Sys.getenv(absent_name, unset = NA_character_)
    )
  )
  expect_identical(observed, c("overridden", "added"))

  # A previously set variable keeps its original value ...
  expect_identical(Sys.getenv(present_name, unset = NA_character_), "original value")
  # ... and a previously unset variable is unset again, not left empty.
  expect_identical(Sys.getenv(absent_name, unset = NA_character_), NA_character_)
})

test_that("the parent environment is restored even when the wrapped call fails", {
  name <- "LISAR_PORTABILITY_ERROR"
  Sys.unsetenv(name)
  on.exit(Sys.unsetenv(name), add = TRUE)
  expect_error(
    lisaR:::lisa_with_child_environment(
      c(LISAR_PORTABILITY_ERROR = "leaked"),
      stop("wrapped failure")
    ),
    "wrapped failure"
  )
  expect_identical(Sys.getenv(name, unset = NA_character_), NA_character_)
})

test_that("an empty environment is a transparent pass-through", {
  expect_identical(lisaR:::lisa_with_child_environment(character(), 42L), 42L)
})

test_that("a subprocess child receives its environment and the parent keeps none of it", {
  root <- lisaR:::lisa_run_root(tempfile("lisa subprocess environment "))
  rscript <- lisaR:::lisa_rscript_executable()
  skip_if_not(file.exists(rscript), "requires the pinned Rscript executable")
  name <- "LISAR_PORTABILITY_CHILD"
  Sys.unsetenv(name)
  on.exit(Sys.unsetenv(name), add = TRUE)
  value <- "child value with spaces ; $ metacharacters"

  result <- lisaR:::lisa_run_subprocess(
    rscript,
    c("--vanilla", "-e", "cat(Sys.getenv('LISAR_PORTABILITY_CHILD'))"),
    stage = "environment portability", run_root = root,
    env = structure(value, names = name)
  )
  expect_identical(result$exit_code, 0L)
  expect_identical(paste(result$stdout, collapse = "\n"), value)
  # The child got it; this process must not have kept it.
  expect_identical(Sys.getenv(name, unset = NA_character_), NA_character_)
})

test_that("a child launched with an environment still runs its script argument", {
  # The exact Windows symptom: with system2(env = ), Rscript treated the
  # leading `R_LIBS=...` token as the file to execute, so the script never ran
  # and no marker was produced. Assert the script itself executed.
  root <- lisaR:::lisa_run_root(tempfile("lisa subprocess script "))
  rscript <- lisaR:::lisa_rscript_executable()
  skip_if_not(file.exists(rscript), "requires the pinned Rscript executable")
  script <- file.path(root, "child.R")
  marker <- file.path(root, "child-marker.txt")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "writeLines(Sys.getenv('LISAR_PORTABILITY_SCRIPT'), args[[1L]])"
  ), script)
  on.exit(Sys.unsetenv("LISAR_PORTABILITY_SCRIPT"), add = TRUE)

  result <- lisaR:::lisa_run_subprocess(
    rscript, c("--vanilla", script, marker),
    stage = "script with environment", run_root = root,
    env = c(LISAR_PORTABILITY_SCRIPT = "executed"),
    required = FALSE
  )
  expect_identical(result$exit_code, 0L)
  expect_true(file.exists(marker))
  expect_identical(readLines(marker, warn = FALSE), "executed")
})

test_that("unsafe subprocess environments are still refused before any child starts", {
  # The portability fix must not have widened the guard: these rejections
  # happen before Sys.setenv() is reached, so nothing is exported either.
  root <- lisaR:::lisa_run_root(tempfile("lisa subprocess guard "))
  rscript <- lisaR:::lisa_rscript_executable()
  skip_if_not(file.exists(rscript), "requires the pinned Rscript executable")
  expect_error(
    lisaR:::lisa_run_subprocess(
      rscript, c("--vanilla", "-e", "quit(status=0)"),
      stage = "bad name", run_root = root,
      env = structure("value", names = "BAD-NAME")
    ),
    "Unsafe subprocess environment"
  )
  expect_identical(Sys.getenv("BAD-NAME", unset = NA_character_), NA_character_)
  expect_error(
    lisaR:::lisa_run_subprocess(
      rscript, c("--vanilla", "-e", "quit(status=0)"),
      stage = "control character", run_root = root,
      env = c(LISAR_PORTABILITY_CONTROL = "value\nwith newline")
    ),
    "Unsafe subprocess environment"
  )
  expect_identical(
    Sys.getenv("LISAR_PORTABILITY_CONTROL", unset = NA_character_), NA_character_
  )
})

test_that("the managed dictionary cache root is canonical and forward-slashed", {
  # R/resources.R used plain normalizePath(), which returns backslashes on
  # Windows and leaves a not-yet-created root in the caller spelling. Paths
  # derived from it then could not be compared with the canonical registry
  # path reported by install_lisa_resource().
  cache <- file.path(tempfile("lisa cache root "), "cache")
  withr::local_options(list(lisaR.dictionary_cache_root = cache))
  # Not created yet: the future path must already be absolute and canonical.
  future <- lisaR:::lisa_dictionary_cache_root()
  expect_false(grepl("\\\\", future))
  expect_identical(future, lisaR:::lisa_path_canonical(cache))

  dir.create(cache, recursive = TRUE)
  existing <- lisaR:::lisa_dictionary_cache_root()
  expect_false(grepl("\\\\", existing))
  expect_identical(
    existing, normalizePath(cache, winslash = "/", mustWork = TRUE)
  )
  # The identity the failing Windows test compared: the default registry path
  # derived from this root matches its own canonical spelling.
  registry <- file.path(existing, "resource_registry.tsv")
  file.create(registry)
  expect_identical(
    registry, normalizePath(registry, winslash = "/", mustWork = TRUE)
  )
})
