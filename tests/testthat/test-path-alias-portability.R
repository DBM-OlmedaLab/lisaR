# Regression coverage for C4 native-platform path portability.
#
# One directory is reachable under several absolute spellings: macOS resolves
# `/var` through a symbolic link to `/private/var`, and Windows still hands out
# 8.3 short names such as `C:/Users/FIXTUR~1`. `normalizePath(mustWork=FALSE)`
# returns a path that does not exist yet unchanged, so a planned destination
# keeps the caller spelling while its existing parent gets canonicalized, and
# spelling-based comparisons then read one directory as two.
#
# A symlinked base directory reproduces the POSIX half of that on any platform.
# The Windows 8.3 half is covered only where it can be expressed portably
# (separator and case handling); the rest is verified by the Windows CI job.

alias_fixture <- function(prefix = "lisa-alias-") {
  base <- tempfile(prefix)
  real <- file.path(base, "real")
  dir.create(real, recursive = TRUE)
  alias <- file.path(base, "alias")
  created <- isTRUE(suppressWarnings(file.symlink(real, alias)))
  list(base = base, real = normalizePath(real, winslash = "/",
                                         mustWork = TRUE),
       alias = alias, created = created)
}

test_that("canonicalization resolves an alias prefix under a future path", {
  fixture <- alias_fixture("lisa-canonical-")
  on.exit(lisa_test_cleanup_path(fixture$base), add = TRUE)
  skip_if_not(fixture$created,
              "this platform did not create a directory symbolic link")

  future <- file.path(fixture$alias, "planned", "deeper", "artifact.png")
  canonical <- lisaR:::lisa_path_canonical(future)

  # The existing prefix is resolved, the planned tail is preserved exactly, and
  # nothing is invented on disk.
  expect_identical(
    canonical,
    paste0(fixture$real, "/planned/deeper/artifact.png")
  )
  expect_false(file.exists(canonical))
  expect_identical(
    canonical,
    lisaR:::lisa_path_canonical(
      file.path(fixture$real, "planned", "deeper", "artifact.png")
    )
  )
  # A doubled separator, as `TMPDIR` with a trailing slash produces on macOS,
  # must not change the result.
  expect_identical(
    canonical,
    lisaR:::lisa_path_canonical(
      paste0(fixture$alias, "//planned/deeper/artifact.png")
    )
  )
})

test_that("a future path keeps its own identity and is not followed", {
  fixture <- alias_fixture("lisa-entry-alias-")
  on.exit(lisa_test_cleanup_path(fixture$base), add = TRUE)
  skip_if_not(fixture$created,
              "this platform did not create a directory symbolic link")

  # The final component of an entry alias must never be resolved, otherwise a
  # symbolic link would silently report the identity of its target and escape
  # the link guards.
  inner <- file.path(fixture$real, "inner")
  dir.create(inner)
  link <- file.path(fixture$real, "inner-link")
  skip_if_not(isTRUE(suppressWarnings(file.symlink(inner, link))),
              "this platform did not create a directory symbolic link")

  aliases <- lisaR:::lisa_path_entry_aliases(link)
  expect_true(all(basename(aliases) == "inner-link"))
  expect_false(any(basename(aliases) == "inner"))
})

test_that("relative presentation links stay short across a path alias", {
  fixture <- alias_fixture("lisa-relative-")
  on.exit(lisa_test_cleanup_path(fixture$base), add = TRUE)
  skip_if_not(fixture$created,
              "this platform did not create a directory symbolic link")

  page_dir <- file.path(fixture$real, "report_pages", "analyses", "A")
  dir.create(page_dir, recursive = TRUE)
  # The route target is planned, not yet built: this is the ordering the
  # report builder relies on, and the case `normalizePath()` cannot resolve.
  planned <- file.path(fixture$alias, "report_pages", "genes", "index.html")

  expect_identical(
    lisaR:::lisa_presentation_relative(planned, page_dir),
    "../../genes/index.html"
  )
  # The alias and the real spelling must produce the same link.
  expect_identical(
    lisaR:::lisa_presentation_relative(planned, page_dir),
    lisaR:::lisa_presentation_relative(
      file.path(fixture$real, "report_pages", "genes", "index.html"), page_dir
    )
  )
  expect_false(grepl("^([.][.]/){4}", lisaR:::lisa_presentation_relative(
    planned, page_dir
  )))
})

test_that("guarded paths accept a run_root alias and still reject escapes", {
  fixture <- alias_fixture("lisa-guarded-")
  on.exit(lisa_test_cleanup_path(fixture$base), add = TRUE)
  skip_if_not(fixture$created,
              "this platform did not create a directory symbolic link")

  run_root <- file.path(fixture$real, "run")
  dir.create(run_root)
  alias_root <- file.path(fixture$alias, "run")

  # Requested through the alias, resolved inside the canonical run root.
  resolved <- lisaR:::lisa_guarded_path(
    file.path(alias_root, "outputs", "planned.tsv"), run_root
  )
  expect_identical(resolved, paste0(run_root, "/outputs/planned.tsv"))
  expect_identical(
    resolved,
    lisaR:::lisa_guarded_path(
      file.path(run_root, "outputs", "planned.tsv"), alias_root
    )
  )

  # Containment is unchanged: a sibling of the run root is still an escape,
  # whichever spelling is used to request it.
  sibling <- file.path(fixture$real, "outside")
  dir.create(sibling)
  expect_error(
    lisaR:::lisa_guarded_path(file.path(sibling, "leak.tsv"), run_root),
    "escapes the authorized run_root"
  )
  expect_error(
    lisaR:::lisa_guarded_path(
      file.path(fixture$alias, "outside", "leak.tsv"), run_root
    ),
    "escapes the authorized run_root"
  )
})

test_that("a link inside the run_root is rejected through an alias spelling", {
  fixture <- alias_fixture("lisa-guarded-link-")
  on.exit(lisa_test_cleanup_path(fixture$base), add = TRUE)
  skip_if_not(fixture$created,
              "this platform did not create a directory symbolic link")

  run_root <- file.path(fixture$real, "run")
  dir.create(run_root)
  outside <- file.path(fixture$real, "outside")
  dir.create(outside)
  link <- file.path(run_root, "escape")
  skip_if_not(isTRUE(suppressWarnings(file.symlink(outside, link))),
              "this platform did not create a directory symbolic link")

  # Before the alias repair this guard only compared one spelling, so a link
  # requested through the alias of the run root slipped past it.
  expect_error(
    lisaR:::lisa_guarded_path(file.path(link, "leak.tsv"), run_root),
    "Symbolic links are not allowed inside the run_root"
  )
  expect_error(
    lisaR:::lisa_guarded_path(
      file.path(fixture$alias, "run", "escape", "leak.tsv"), run_root
    ),
    "Symbolic links are not allowed inside the run_root"
  )
})

test_that("path keys carry Windows separator and case semantics", {
  # Portable on every platform: lisaR compares managed paths with forward
  # slashes, and only Windows folds case.
  expect_identical(
    lisaR:::lisa_path_key("C:\\Users\\fixture-user\\work"),
    lisaR:::lisa_path_key("C:/Users/fixture-user/work")
  )
  expect_identical(lisaR:::lisa_path_key("/a/b/"), lisaR:::lisa_path_key("/a/b"))
  expect_true(lisaR:::lisa_path_within("/a/b/c", "/a/b"))
  expect_false(lisaR:::lisa_path_within("/a/bc", "/a/b"))
  folds <- identical(lisaR:::lisa_path_key("/A/B"),
                     lisaR:::lisa_path_key("/a/b"))
  expect_identical(folds, .Platform$OS.type == "windows")
})

test_that("staged KEGG receipts are rebased from every root spelling", {
  fixture <- alias_fixture("lisa-kegg-alias-")
  on.exit(lisa_test_cleanup_path(fixture$base), add = TRUE)
  skip_if_not(fixture$created,
              "this platform did not create a directory symbolic link")

  work_root <- file.path(fixture$alias, "work")
  dir.create(work_root)
  promoted <- file.path(fixture$real, "promoted")
  # The receipt stores the alias spelling, which is what a writer that used
  # `tempdir()` would have recorded on macOS or Windows.
  status <- data.frame(message = file.path(work_root, "outputs", "map.png"),
                       stringsAsFactors = FALSE)
  lisaR:::write_lisa_tsv(status, file.path(work_root,
                                           "kegg_extension_status.tsv"))

  expect_true(lisaR:::lisa_rebase_kegg_extension_paths(work_root, promoted))
  rebased <- paste(
    lisaR:::read_lisa_tsv(file.path(work_root, "kegg_extension_status.tsv")),
    collapse = "\n"
  )
  promoted_root <- normalizePath(promoted, winslash = "/", mustWork = FALSE)
  requested_root <- gsub("\\\\", "/", work_root)
  canonical_root <- normalizePath(work_root, winslash = "/", mustWork = TRUE)
  expect_true(grepl(promoted_root, rebased, fixed = TRUE))
  expect_false(grepl(requested_root, rebased, fixed = TRUE))
  expect_false(grepl(canonical_root, rebased, fixed = TRUE))
})
