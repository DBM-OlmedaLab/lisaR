test_that("Technical identifiers use one portable ASCII grammar", {
  invalid <- list(
    "non-empty" = c("", "   "),
    "Unsafe" = c(
      "../escape", "/absolute", "a/b", "a\\b", ".", "..", "bad\nname"
    ),
    "ASCII" = c(
      "A:B", "A B", "A*B", "x|y", "\"quoted\"", "<tag>", "what?",
      ".hidden", "_leading", "-leading", "trailing.", "trailing ",
      "caf\u00e9", "nai\u0308ve", "\u03b1", "sample\U0001f9ec"
    ),
    "reserved Windows" = c(
      "CON", "con", "PRN", "Aux", "NUL", "COM1", "com9", "LPT1", "lpt9",
      "CON.txt", "prn.csv", "AUX.rds", "nul.log", "COM1.tsv", "lPt9.data"
    )
  )
  for (error_pattern in names(invalid)) {
    for (value in invalid[[error_pattern]]) {
      expect_error(
        lisaR:::lisa_safe_id(value, "analysis_id"), error_pattern,
        info = value
      )
    }
  }
  expect_error(lisaR:::lisa_safe_id(NULL), "non-empty")
  expect_error(lisaR:::lisa_safe_id(character()), "non-empty")
  expect_error(lisaR:::lisa_safe_id(NA_character_), "non-empty")

  valid <- c(
    "A", "a1", "sample-A_1", "GOBP-C2", "RIAZ_GSE91061",
    "snapshot.2026-08-30", "sample_", "sample-", "COM10", "LPT10",
    "CON_value", "PRN-report"
  )
  observed <- vapply(valid, lisaR:::lisa_safe_id, character(1))
  expect_identical(unname(observed), valid)
  expect_false(identical(
    lisaR:::lisa_safe_id("sample-A"),
    lisaR:::lisa_safe_id("sample_A")
  ))
})

test_that("Scientific resource keys map injectively to portable path components", {
  keys <- c(
    "all", "hsa00010", "hsa:1", "hsa-1", "CON",
    "lisa-key-v1-6873613a31", "se\u00f1al:\u03b1"
  )
  first <- vapply(
    keys, lisaR:::lisa_external_key_path_component, character(1)
  )
  second <- vapply(
    keys, lisaR:::lisa_external_key_path_component, character(1)
  )
  expect_identical(first, second)
  expect_false(anyDuplicated(unname(first)) > 0L)
  expect_identical(unname(first[c("all", "hsa00010")]), c("all", "hsa00010"))
  expect_true(startsWith(unname(first[["hsa:1"]]), "lisa-key-v1-"))
  expect_false(grepl(":", unname(first[["hsa:1"]]), fixed = TRUE))
  expect_silent(vapply(unname(first), lisaR:::lisa_safe_id, character(1)))
})

test_that("Windows CI proves supported file and directory link detection", {
  skip_if(.Platform$OS.type != "windows", "Windows-specific link probe")
  base <- tempfile("lisa-windows-link-probe-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  file_target <- file.path(base, "target.txt")
  directory_target <- file.path(base, "target-directory")
  file_link <- file.path(base, "file-link")
  directory_link <- file.path(base, "directory-link")
  writeLines("outside-original", file_target, useBytes = TRUE)
  dir.create(directory_target)
  writeLines("directory-original", file.path(directory_target, "sentinel.txt"))

  expect_no_error(fs::link_create(file_target, file_link, symbolic = TRUE))
  expect_no_error(
    fs::link_create(directory_target, directory_link, symbolic = TRUE)
  )
  on.exit(for (path in c(file_link, directory_link)) {
    if (lisaR:::lisa_path_is_link(path)) lisaR:::lisa_link_delete(path)
  }, add = TRUE)
  expect_true(lisaR:::lisa_path_is_link(file_link))
  expect_true(lisaR:::lisa_path_is_link(directory_link))
  expect_true(lisaR:::lisa_path_entry_exists(file_link))
  expect_true(lisaR:::lisa_path_entry_exists(directory_link))

  lisaR:::lisa_link_delete(file_link)
  lisaR:::lisa_link_delete(directory_link)
  expect_false(lisaR:::lisa_path_entry_exists(file_link))
  expect_false(lisaR:::lisa_path_entry_exists(directory_link))
  expect_identical(readLines(file_target), "outside-original")
  expect_identical(
    readLines(file.path(directory_target, "sentinel.txt")),
    "directory-original"
  )
})

test_that("Security confines writes, copies, cleanup, and symlinks to a safe run root", {
  root <- tempfile("lisa-contract-root-")
  lisaR:::lisa_run_root(root)
  expect_error(lisaR:::lisa_guarded_path(file.path(root, "..", "escape"), root), "escapes")
  expect_error(lisaR:::lisa_run_root("/"), "Unsafe")
  lisaR:::lisa_guarded_write(file.path(root, "nested", "ok.txt"), function(path) writeLines("ok", path), root)
  expect_true(file.exists(file.path(root, "nested", "ok.txt")))
  src <- tempfile("lisa contract source "); writeLines("copied", src)
  lisaR:::lisa_guarded_copy(src, file.path(root, "nested", "copy.txt"), TRUE, root)
  expect_true(file.exists(file.path(root, "nested", "copy.txt")))
  expect_error(lisaR:::lisa_guarded_copy(src, file.path(root, "nested", "copy.txt"), FALSE, root), "overwrite")
  lisaR:::lisa_guarded_rename(file.path(root, "nested", "copy.txt"), file.path(root, "nested", "renamed file.txt"), run_root = root)
  expect_true(file.exists(file.path(root, "nested", "renamed file.txt")))
  lisaR:::lisa_guarded_delete(file.path(root, "nested", "renamed file.txt"), run_root = root)
  expect_false(file.exists(file.path(root, "nested", "renamed file.txt")))
  expect_error(lisaR:::lisa_guarded_delete(root, recursive = TRUE, run_root = root), "Refusing")
  outside <- tempfile("lisa-contract-outside-"); dir.create(outside)
  link <- file.path(root, "link")
  linked <- suppressWarnings(file.symlink(outside, link))
  if (!isTRUE(linked)) skip("file.symlink() is unavailable in this test environment")
  expect_error(lisaR:::lisa_guarded_write(file.path(root, "link", "escape.txt"), function(path) writeLines("no", path), root), "Symbolic")
})

test_that("Lexical checks reject root, parent, destination, relative, and dangling links", {
  base <- tempfile("lisa-symlink-adversarial-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  outside <- file.path(base, "outside")
  real_root <- file.path(base, "real-root")
  dir.create(outside)
  dir.create(real_root)
  sentinel <- file.path(outside, "sentinel.txt")
  writeLines("outside-original", sentinel, useBytes = TRUE)

  root_link <- file.path(base, "root-link")
  linked <- suppressWarnings(file.symlink(real_root, root_link))
  if (!isTRUE(linked)) skip("file.symlink() is unavailable in this test environment")
  expect_error(lisaR:::lisa_run_root(root_link), "Symbolic")

  parent_link <- file.path(real_root, "linked-parent")
  expect_true(suppressWarnings(file.symlink(outside, parent_link)))
  expect_error(
    lisaR:::lisa_guarded_write(
      file.path(parent_link, "sentinel.txt"),
      function(path) writeLines("escaped", path), real_root
    ),
    "[Ss]ymbolic"
  )

  destination_link <- file.path(real_root, "destination.txt")
  expect_true(suppressWarnings(file.symlink(sentinel, destination_link)))
  source <- file.path(base, "source.txt")
  writeLines("source", source)
  expect_error(
    lisaR:::lisa_guarded_write(
      destination_link, function(path) writeLines("escaped", path), real_root
    ),
    "Symbolic"
  )
  expect_error(
    lisaR:::lisa_guarded_copy(source, destination_link, TRUE, real_root),
    "Symbolic"
  )
  rename_source <- file.path(real_root, "rename-source.txt")
  writeLines("rename-source", rename_source)
  expect_error(
    lisaR:::lisa_guarded_rename(
      rename_source, destination_link, TRUE, real_root
    ),
    "Symbolic"
  )
  expect_true(file.exists(rename_source))
  expect_error(
    lisaR:::lisa_guarded_delete(destination_link, recursive = TRUE,
                                run_root = real_root),
    "Symbolic"
  )

  relative_link <- file.path(real_root, "relative-link")
  expect_true(suppressWarnings(file.symlink("../outside", relative_link)))
  expect_error(
    lisaR:::lisa_guarded_write(
      file.path(relative_link, "..", "safe-looking.txt"),
      function(path) writeLines("escaped", path), real_root
    ),
    "Symbolic"
  )

  dangling_link <- file.path(real_root, "dangling-link")
  expect_true(suppressWarnings(file.symlink("../absent-outside", dangling_link)))
  # Base R may report a dangling Windows reparse entry as existing. lisaR's
  # contract is its link identity and rejection, not file.exists() spelling.
  expect_true(lisaR:::lisa_path_is_link(dangling_link))
  expect_error(
    lisaR:::lisa_guarded_write(
      dangling_link, function(path) writeLines("escaped", path), real_root
    ),
    "Symbolic"
  )
  expect_error(
    lisaR:::lisa_guarded_copy(source, dangling_link, TRUE, real_root),
    "Symbolic"
  )
  dangling_rename_source <- file.path(real_root, "dangling-rename-source.txt")
  writeLines("rename-source", dangling_rename_source)
  expect_error(
    lisaR:::lisa_guarded_rename(
      dangling_rename_source, dangling_link, TRUE, real_root
    ),
    "Symbolic"
  )
  expect_true(file.exists(dangling_rename_source))
  expect_error(
    lisaR:::lisa_guarded_delete(dangling_link, recursive = TRUE,
                                run_root = real_root),
    "Symbolic"
  )
  expect_identical(readLines(sentinel, warn = FALSE), "outside-original")
  expect_false(file.exists(file.path(base, "absent-outside")))
})

test_that("Run-tree scans never follow nested links for manifests or deletion", {
  base <- tempfile("lisa-symlink-tree-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  root <- file.path(base, "root")
  outside <- file.path(base, "outside")
  nested <- file.path(root, "one", "two")
  dir.create(nested, recursive = TRUE)
  dir.create(outside)
  sentinel <- file.path(outside, "sentinel.txt")
  writeLines("outside-original", sentinel, useBytes = TRUE)
  link <- file.path(nested, "relative-outside")
  linked <- suppressWarnings(file.symlink("../../../outside", link))
  if (!isTRUE(linked)) skip("file.symlink() is unavailable in this test environment")

  expect_error(lisaR:::lisa_assert_run_tree_safe(root), "Symbolic")
  expect_error(lisaR:::lisa_run_manifest(root), "Symbolic")
  expect_error(
    lisaR:::lisa_guarded_delete(file.path(root, "one"), recursive = TRUE,
                                run_root = root),
    "Symbolic"
  )
  expect_identical(readLines(sentinel, warn = FALSE), "outside-original")
  expect_true(file.exists(sentinel))
})

test_that("Run-tree scans inspect links before metadata or recursive traversal", {
  base <- tempfile("lisa-symlink-call-order-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  root <- file.path(base, "root")
  outside <- file.path(base, "outside")
  dir.create(root)
  dir.create(outside)
  link <- file.path(root, "outside-link")
  linked <- suppressWarnings(file.symlink(outside, link))
  if (!isTRUE(linked)) skip("file.symlink() is unavailable in this test environment")

  real_file_info <- base::file.info
  real_list_files <- base::list.files
  calls <- new.env(parent = emptyenv())
  calls$file_info_link <- FALSE
  calls$list_files_link <- FALSE
  testthat::local_mocked_bindings(
    file.info = function(path, ...) {
      targets <- vapply(
        as.character(path), lisaR:::lisa_path_is_link, logical(1)
      )
      if (any(targets)) {
        calls$file_info_link <- TRUE
        stop("file.info followed a symbolic link")
      }
      real_file_info(path, ...)
    },
    list.files = function(path = ".", ...) {
      if (lisaR:::lisa_path_is_link(path)) {
        calls$list_files_link <- TRUE
        stop("list.files followed a symbolic link")
      }
      real_list_files(path, ...)
    },
    .package = "base"
  )
  expect_error(lisaR:::lisa_scan_run_tree(root), "Symbolic")
  expect_false(calls$file_info_link)
  expect_false(calls$list_files_link)
})

test_that("Guarded mutations are atomic and clean temporary files after writer failure", {
  root <- lisaR:::lisa_run_root(tempfile("lisa-atomic-write-"))
  destination <- file.path(root, "nested", "artifact.txt")
  lisaR:::lisa_guarded_write(
    destination, function(path) writeLines("original", path), root
  )
  expect_error(
    lisaR:::lisa_guarded_write(destination, function(path) {
      writeLines("partial", path)
      stop("simulated writer failure")
    }, root),
    "simulated writer failure"
  )
  expect_identical(readLines(destination, warn = FALSE), "original")

  absent <- file.path(root, "nested", "absent.txt")
  expect_error(
    lisaR:::lisa_guarded_write(absent, function(path) {
      writeLines("partial", path)
      stop("simulated writer failure")
    }, root),
    "simulated writer failure"
  )
  expect_false(file.exists(absent))

  outside <- tempfile("lisa-writer-symlink-outside-")
  writeLines("outside-original", outside, useBytes = TRUE)
  if (.Platform$OS.type != "windows") Sys.chmod(outside, "0640")
  outside_mode <- as.integer(file.info(outside)$mode)
  probe <- tempfile("lisa-symlink-probe-", tmpdir = root)
  if (!isTRUE(suppressWarnings(file.symlink(outside, probe)))) {
    skip("file.symlink() is unavailable in this test environment")
  }
  lisaR:::lisa_link_delete(probe)
  expect_error(
    lisaR:::lisa_guarded_write(destination, function(path) {
      unlink(path)
      if (!isTRUE(suppressWarnings(file.symlink(outside, path)))) {
        stop("file.symlink() became unavailable")
      }
    }, root),
    "[Ss]ymbolic"
  )
  expect_identical(readLines(destination, warn = FALSE), "original")
  expect_identical(readLines(outside, warn = FALSE), "outside-original")
  expect_identical(as.integer(file.info(outside)$mode), outside_mode)

  png <- file.path(root, "nested", "extension-required.png")
  lisaR:::lisa_guarded_write(png, function(path) {
    if (!identical(tools::file_ext(path), "png")) {
      stop("writer requires a .png extension")
    }
    writeBin(as.raw(1:4), path)
  }, root)
  expect_identical(readBin(png, what = "raw", n = 4L), as.raw(1:4))

  source <- tempfile("lisa-copy-source-")
  writeLines("copied", source)
  copied <- file.path(root, "nested", "copied.txt")
  lisaR:::lisa_guarded_copy(source, copied, run_root = root)
  renamed <- file.path(root, "other", "renamed.txt")
  lisaR:::lisa_guarded_rename(copied, renamed, run_root = root)
  expect_false(file.exists(copied))
  expect_identical(readLines(renamed, warn = FALSE), "copied")

  rollback_source <- file.path(root, "nested", "rollback-source.txt")
  rollback_target <- file.path(root, "other", "rollback-target.txt")
  writeLines("must-survive", rollback_source, useBytes = TRUE)
  real_rename <- base::file.rename
  rename_calls <- 0L
  testthat::local_mocked_bindings(
    file.rename = function(from, to) {
      rename_calls <<- rename_calls + 1L
      if (rename_calls == 2L) return(FALSE)
      real_rename(from, to)
    },
    .package = "base"
  )
  expect_error(
    lisaR:::lisa_guarded_rename(rollback_source, rollback_target,
                                run_root = root),
    "Atomic filesystem promotion failed"
  )
  expect_identical(readLines(rollback_source, warn = FALSE), "must-survive")
  expect_false(file.exists(rollback_target))

  debris <- list.files(root, pattern = "^[.]lisa-(write|copy|rename|backup)-",
                       recursive = TRUE, all.files = TRUE, full.names = TRUE)
  expect_length(debris, 0L)
})

test_that("Failed guarded-rename restoration exposes the staged source", {
  root <- lisaR:::lisa_run_root(tempfile("lisa-rename-recovery-"))
  on.exit(lisa_test_cleanup_path(root), add = TRUE)
  source <- file.path(root, "source.txt")
  destination <- file.path(root, "destination.txt")
  writeLines("must-remain-recoverable", source, useBytes = TRUE)

  real_rename <- base::file.rename
  rename_calls <- 0L
  testthat::local_mocked_bindings(
    file.rename = function(from, to) {
      rename_calls <<- rename_calls + 1L
      if (rename_calls %in% c(2L, 3L)) return(FALSE)
      real_rename(from, to)
    },
    .package = "base"
  )

  error <- tryCatch(
    lisaR:::lisa_guarded_rename(source, destination, run_root = root),
    lisa_filesystem_recovery_error = identity
  )
  expect_s3_class(error, "lisa_filesystem_recovery_error")
  expect_identical(error$source, source)
  expect_identical(error$destination, destination)
  expect_true(nzchar(error$staged))
  expect_true(error$staged %in% error$recovery_paths)
  expect_true(file.exists(error$staged))
  expect_identical(
    readLines(error$staged, warn = FALSE), "must-remain-recoverable"
  )
  expect_false(file.exists(source))
  expect_false(file.exists(destination))
})

test_that("Promotion rechecks the complete staging tree and destination", {
  base <- tempfile("lisa-promotion-symlink-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  outside <- file.path(base, "outside")
  dir.create(outside)
  sentinel <- file.path(outside, "sentinel.txt")
  writeLines("outside-original", sentinel, useBytes = TRUE)
  final <- file.path(base, "final")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "symlinkpromotion"
  )
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  nested <- file.path(tx$staging_dir, "nested")
  dir.create(nested)
  linked <- suppressWarnings(file.symlink(outside, file.path(nested, "escape")))
  if (!isTRUE(linked)) {
    lisaR:::lisa_transaction_abort(tx)
    skip("file.symlink() is unavailable in this test environment")
  }

  expect_error(lisaR:::lisa_transaction_promote(tx), "Symbolic")
  expect_false(dir.exists(final))
  expect_identical(readLines(sentinel, warn = FALSE), "outside-original")
  expect_false(file.exists(file.path(outside, "run_manifest.tsv")))
  lisaR:::lisa_transaction_abort(tx)
})

test_that("Promotion refuses a destination replaced by a symbolic link", {
  base <- tempfile("lisa-promotion-destination-link-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  outside <- file.path(base, "outside")
  dir.create(outside)
  sentinel <- file.path(outside, "sentinel.txt")
  writeLines("outside-original", sentinel, useBytes = TRUE)
  final <- file.path(base, "final")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "destinationlink"
  )
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  linked <- suppressWarnings(file.symlink(outside, final))
  if (!isTRUE(linked)) {
    lisaR:::lisa_transaction_abort(tx)
    skip("file.symlink() is unavailable in this test environment")
  }

  expect_error(lisaR:::lisa_transaction_promote(tx), "Symbolic")
  expect_true(lisaR:::lisa_path_is_link(final))
  expect_identical(readLines(sentinel, warn = FALSE), "outside-original")
  expect_false(file.exists(file.path(outside, "run_manifest.tsv")))
  lisaR:::lisa_transaction_abort(tx)
})

test_that("Security subprocesses preserve arguments and diagnose required failures", {
  root <- lisaR:::lisa_run_root(tempfile("lisa contract subprocess "))
  rscript <- lisaR:::lisa_rscript_executable()
  expect_true(file.exists(rscript))
  r <- lisaR:::lisa_run_subprocess(
    rscript,
    c("-e", "cat('out value'); message('err value')", "spaces ' quotes ; $ metacharacters"),
    stage = "argument safety", run_root = root
  )
  expect_identical(r$exit_code, 0L)
  expect_match(paste(r$stdout, collapse = "\n"), "out value")
  expect_match(paste(r$stderr, collapse = "\n"), "err value")
  optional <- lisaR:::lisa_run_subprocess(rscript, c("-e", "message('bad'); quit(status=7)"),
                                          required = FALSE, stage = "optional", run_root = root)
  expect_identical(optional$exit_code, 7L)
  expect_match(optional$diagnostics, "bad")
  expect_error(lisaR:::lisa_run_subprocess(rscript, c("-e", "quit(status=3)"),
                                            required = TRUE, stage = "required", run_root = root),
               "Required stage failed.*required.*3")

  environment_value <- "child value with spaces ; $ metacharacters"
  inherited <- lisaR:::lisa_run_subprocess(
    rscript,
    c("--vanilla", "-e", "cat(Sys.getenv('LISAR_CHILD_TEST'))"),
    stage = "environment safety", run_root = root,
    env = c(LISAR_CHILD_TEST = environment_value)
  )
  expect_identical(inherited$exit_code, 0L)
  expect_identical(paste(inherited$stdout, collapse = "\n"), environment_value)
  expect_error(
    lisaR:::lisa_run_subprocess(
      rscript, c("--vanilla", "-e", "quit(status=0)"),
      stage = "invalid environment", run_root = root,
      env = structure("value", names = "BAD-NAME")
    ),
    "Unsafe subprocess environment"
  )
  expect_error(
    lisaR:::lisa_run_subprocess(
      rscript, c("--vanilla", "-e", "quit(status=0)"),
      stage = "missing environment name", run_root = root,
      env = structure("value", names = NA_character_)
    ),
    "unique, non-empty names"
  )
  expect_error(
    lisaR:::lisa_run_subprocess(
      rscript, c("--vanilla", "-e", "quit(status=0)"),
      stage = "duplicate environment name", run_root = root,
      env = c(LISAR_DUPLICATE = "first", LISAR_DUPLICATE = "second")
    ),
    "unique, non-empty names"
  )
})

test_that("External parent aliases are allowed but managed-tree aliases are not", {
  base <- tempfile("lisa-external-alias-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  real_parent <- file.path(base, "real-parent")
  alias_parent <- file.path(base, "parent-alias")
  outside <- file.path(base, "outside")
  dir.create(real_parent)
  dir.create(outside)
  if (!isTRUE(suppressWarnings(file.symlink(real_parent, alias_parent)))) {
    skip("file.symlink() is unavailable in this test environment")
  }

  requested_root <- file.path(alias_parent, "run")
  root <- lisaR:::lisa_run_root(requested_root)
  expect_identical(
    root,
    normalizePath(file.path(real_parent, "run"), winslash = "/",
                  mustWork = TRUE)
  )
  lisaR:::lisa_guarded_write(
    file.path(requested_root, "plain.txt"),
    function(path) writeLines("ok", path), requested_root
  )
  expect_true(file.exists(file.path(root, "plain.txt")))

  inside_target <- file.path(root, "inside-target")
  dir.create(inside_target)
  expect_true(suppressWarnings(file.symlink(
    inside_target, file.path(root, "inside-link")
  )))
  expect_true(suppressWarnings(file.symlink(
    outside, file.path(root, "outside-link")
  )))
  expect_error(
    lisaR:::lisa_guarded_write(
      file.path(requested_root, "inside-link", "blocked.txt"),
      function(path) writeLines("blocked", path), requested_root
    ),
    "Symbolic"
  )
  expect_error(
    lisaR:::lisa_guarded_write(
      file.path(requested_root, "outside-link", "blocked.txt"),
      function(path) writeLines("blocked", path), requested_root
    ),
    "Symbolic"
  )

  leaf_alias <- file.path(base, "run-leaf-alias")
  expect_true(suppressWarnings(file.symlink(root, leaf_alias)))
  expect_error(lisaR:::lisa_run_root(leaf_alias), "Symbolic")
})

test_that("Managed boundaries reject dot traversal and preserve parent modes", {
  base <- tempfile("lisa-boundary-modes-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  if (.Platform$OS.type != "windows") Sys.chmod(base, "0755")
  base_mode <- as.integer(file.info(base)$mode)

  expect_error(
    lisaR:::lisa_run_root(file.path(base, "child", "..")),
    "may not contain"
  )
  expect_identical(as.integer(file.info(base)$mode), base_mode)

  existing <- file.path(base, "existing")
  dir.create(existing)
  if (.Platform$OS.type != "windows") Sys.chmod(existing, "0755")
  existing_mode <- as.integer(file.info(existing)$mode)
  lisaR:::lisa_run_root(existing)
  expect_identical(as.integer(file.info(existing)$mode), existing_mode)

  root <- lisaR:::lisa_run_root(file.path(base, "new-root"))
  nested <- lisaR:::lisa_guarded_dir_create(file.path(root, "nested"), root)
  artifact <- file.path(nested, "artifact.txt")
  lisaR:::lisa_guarded_write(
    artifact, function(path) writeLines("private", path), root
  )
  source <- file.path(base, "source.txt")
  writeLines("copy", source)
  if (.Platform$OS.type != "windows") Sys.chmod(source, "0600")
  copied <- file.path(nested, "copied.txt")
  lisaR:::lisa_guarded_copy(source, copied, run_root = root)

  if (.Platform$OS.type != "windows") {
    expect_identical(as.integer(file.info(root)$mode),
                     as.integer(as.octmode("0700")))
    expect_identical(as.integer(file.info(nested)$mode),
                     as.integer(as.octmode("0700")))
    expect_identical(as.integer(file.info(artifact)$mode),
                     as.integer(as.octmode("0600")))
    expect_identical(
      bitwAnd(as.integer(file.info(copied)$mode),
              bitwNot(as.integer(file.info(source)$mode))),
      0L
    )
  }
  expect_identical(as.integer(file.info(base)$mode), base_mode)
})

test_that("Configured managed paths resolve existing parent traversal safely", {
  base <- tempfile("lisa-config-managed-")
  config_dir <- file.path(base, "project", "config")
  dir.create(config_dir, recursive = TRUE)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  # A managed path is returned in the canonical spelling of its existing
  # parent. `tempfile()` returns a platform alias of that directory on macOS
  # (`/var` -> `/private/var`) and Windows (8.3 short names), so the fixture
  # declares the canonical spelling. Every assertion below stays exact.
  base <- normalizePath(base, winslash = "/", mustWork = TRUE)
  config_dir <- file.path(base, "project", "config")

  resolved <- lisaR:::lisa_config_managed_path(
    file.path("..", "results", "study"), config_dir
  )
  expect_identical(
    resolved,
    gsub("\\\\", "/", file.path(base, "project", "results", "study"))
  )
  expect_false(any(strsplit(resolved, "/", fixed = TRUE)[[1L]] %in%
                     c(".", "..")))
  expect_error(
    lisaR:::lisa_managed_destination(
      file.path(config_dir, "child", "..", "study")
    ),
    "may not contain"
  )
  absolute_raw <- lisaR:::lisa_config_managed_path(
    file.path(config_dir, "child", "..", "study"), config_dir
  )
  expect_error(
    lisaR:::lisa_managed_destination(absolute_raw),
    "may not contain"
  )
  expect_error(
    lisaR:::lisa_config_managed_path(
      file.path("missing", "..", "study"), config_dir
    ),
    "may not traverse"
  )
  base_file <- file.path(base, "base-file")
  writeLines("not a directory", base_file)
  expect_error(
    lisaR:::lisa_config_managed_path("study", base_file),
    "base must be one existing directory"
  )
  broken_base <- file.path(base, "broken-base")
  if (isTRUE(suppressWarnings(file.symlink(
    file.path(base, "absent-base-target"), broken_base
  )))) {
    expect_error(
      lisaR:::lisa_config_managed_path("study", broken_base),
      "base must be one existing directory"
    )
  }

  existing_file <- file.path(config_dir, "existing-file")
  writeLines("not a directory", existing_file)
  expect_error(
    lisaR:::lisa_config_managed_path(
      file.path("existing-file", "child"), config_dir
    ),
    "not a directory"
  )
  expect_error(
    lisaR:::lisa_config_managed_path(
      file.path("existing-file", "..", "study"), config_dir
    ),
    "not a directory"
  )

  link_target <- file.path(base, "link-target", "nested")
  dir.create(link_target, recursive = TRUE)
  link <- file.path(config_dir, "linked")
  if (isTRUE(suppressWarnings(file.symlink(link_target, link)))) {
    expect_identical(
      lisaR:::lisa_config_managed_path(
        file.path("linked", "..", "study"), config_dir
      ),
      gsub("\\\\", "/", file.path(dirname(link_target), "study"))
    )
  }

  link_file <- file.path(config_dir, "linked-file")
  if (isTRUE(suppressWarnings(file.symlink(existing_file, link_file)))) {
    expect_error(
      lisaR:::lisa_config_managed_path(
        file.path("linked-file", "..", "study"), config_dir
      ),
      "not a directory"
    )
  }

  dangling <- file.path(config_dir, "dangling")
  if (isTRUE(suppressWarnings(file.symlink(
    file.path(base, "absent-target"), dangling
  )))) {
    # Windows may expose a dangling file.symlink() as a non-directory entry
    # rather than as a link. Both diagnoses must fail before child traversal.
    expect_error(
      lisaR:::lisa_config_managed_path(
        file.path("dangling", "study"), config_dir
      ),
      "broken symbolic link|not a directory"
    )
  }
})

test_that("HOME is an external parent, never the managed run root", {
  fake_home <- tempfile("lisa-fake-home-")
  dir.create(fake_home)
  on.exit(lisa_test_cleanup_path(fake_home), add = TRUE)
  if (.Platform$OS.type != "windows") Sys.chmod(fake_home, "0755")
  home_mode <- as.integer(file.info(fake_home)$mode)

  expect_error(lisaR:::lisa_run_root("~"), "non-home")
  expect_error(lisaR:::lisa_existing_run_root("~"), "non-home")
  expect_error(
    lisaR:::lisa_run_root(file.path(fake_home, "child", "..")),
    "may not contain"
  )
  root <- lisaR:::lisa_run_root(file.path(fake_home, "managed-run"))
  expect_true(dir.exists(root))
  expect_identical(as.integer(file.info(fake_home)$mode), home_mode)
})

test_that("Missing roots are checked without being created", {
  missing <- tempfile("lisa-missing-root-")
  expect_error(lisaR:::lisa_existing_run_root(missing), "does not exist")
  expect_false(lisaR:::lisa_path_entry_exists(missing))
  expect_error(lisaR:::verify_run(missing), "does not exist")
  expect_false(lisaR:::lisa_path_entry_exists(missing))
})

test_that("Boundary deletion does not turn the parent into a managed root", {
  parent <- tempfile("lisa-boundary-delete-")
  dir.create(parent)
  on.exit(lisa_test_cleanup_path(parent), add = TRUE)
  file <- file.path(parent, "one.txt")
  directory <- file.path(parent, "tree")
  writeLines("one", file)
  dir.create(directory)
  writeLines("two", file.path(directory, "two.txt"))

  lisaR:::lisa_guarded_delete(file, run_root = NULL)
  lisaR:::lisa_guarded_delete(directory, recursive = TRUE, run_root = NULL)
  expect_false(file.exists(file))
  expect_false(dir.exists(directory))
  expect_true(dir.exists(parent))
})

test_that("Run-tree scan rechecks links after metadata inspection", {
  base <- tempfile("lisa-scan-swap-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  root <- file.path(base, "root")
  outside <- file.path(base, "outside.txt")
  victim <- file.path(root, "victim.txt")
  dir.create(root)
  writeLines("outside", outside)
  writeLines("victim", victim)
  victim_key <- lisaR:::lisa_path_key(normalizePath(
    victim, winslash = "/", mustWork = TRUE
  ))
  probe <- file.path(base, "probe")
  if (!isTRUE(suppressWarnings(file.symlink(outside, probe)))) {
    skip("file.symlink() is unavailable in this test environment")
  }
  lisaR:::lisa_link_delete(probe)
  real_path_state <- lisaR:::lisa_path_state
  swapped <- FALSE
  victim_inspections <- 0L
  testthat::local_mocked_bindings(
    lisa_path_state = function(path) {
      if (identical(lisaR:::lisa_path_key(path), victim_key)) {
        victim_inspections <<- victim_inspections + 1L
        if (!swapped && identical(victim_inspections, 2L)) {
          swapped <<- TRUE
          unlink(victim)
          file.symlink(outside, victim)
        }
      }
      real_path_state(path)
    },
    .package = "lisaR"
  )
  expect_error(lisaR:::lisa_scan_run_tree(root), "changed")
  expect_true(swapped)
  expect_identical(readLines(outside, warn = FALSE), "outside")
})

test_that("Private tree modes never chmod a path replaced by a link", {
  if (.Platform$OS.type == "windows") {
    skip("Windows does not apply POSIX chmod transitions")
  }
  base <- tempfile("lisa-mode-swap-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  root <- lisaR:::lisa_run_root(file.path(base, "root"))
  victim <- file.path(root, "victim.txt")
  outside <- file.path(base, "outside.txt")
  writeLines("victim", victim)
  writeLines("outside", outside)
  if (.Platform$OS.type != "windows") Sys.chmod(outside, "0640")
  outside_mode <- as.integer(file.info(outside)$mode)
  probe <- file.path(base, "probe")
  if (!isTRUE(suppressWarnings(file.symlink(outside, probe)))) {
    skip("file.symlink() is unavailable in this test environment")
  }
  lisaR:::lisa_link_delete(probe)
  real_scan <- lisaR:::lisa_scan_run_tree
  scans <- 0L
  testthat::local_mocked_bindings(
    lisa_scan_run_tree = function(run_root) {
      scans <<- scans + 1L
      observed <- real_scan(run_root)
      if (scans == 2L) {
        unlink(victim)
        file.symlink(outside, victim)
      }
      observed
    },
    .package = "lisaR"
  )
  expect_error(lisaR:::lisa_private_tree_modes(root), "Symbolic")
  expect_identical(readLines(outside, warn = FALSE), "outside")
  expect_identical(as.integer(file.info(outside)$mode), outside_mode)
})

test_that("Recursive deletion quarantines before no-follow traversal", {
  parent <- tempfile("lisa-delete-quarantine-")
  dir.create(parent)
  on.exit(lisa_test_cleanup_path(parent), add = TRUE)
  target <- file.path(parent, "target")
  dir.create(target)
  writeLines("artifact", file.path(target, "artifact.txt"))
  observed <- new.env(parent = emptyenv())
  observed$path <- ""
  real_unlink <- lisaR:::lisa_unlink_no_follow
  testthat::local_mocked_bindings(
    lisa_unlink_no_follow = function(path) {
      if (!nzchar(observed$path)) {
        observed$path <- path
        expect_false(lisaR:::lisa_path_entry_exists(target))
        expect_match(basename(path), "^[.]lisa-delete-")
      }
      real_unlink(path)
    },
    .package = "lisaR"
  )
  lisaR:::lisa_guarded_delete(target, recursive = TRUE, run_root = NULL)
  expect_false(lisaR:::lisa_path_entry_exists(target))
  expect_true(nzchar(observed$path))
})

test_that("Managed directory promotion is exactly one validated rename", {
  parent <- tempfile("lisa-single-rename-")
  dir.create(parent)
  on.exit(lisa_test_cleanup_path(parent), add = TRUE)
  staging <- lisaR:::lisa_run_root(file.path(parent, "staging"))
  destination <- file.path(parent, "final")
  writeLines("artifact", file.path(staging, "artifact.txt"))
  real_rename <- base::file.rename
  rename_calls <- 0L
  testthat::local_mocked_bindings(
    file.rename = function(from, to) {
      rename_calls <<- rename_calls + 1L
      real_rename(from, to)
    },
    .package = "base"
  )
  lisaR:::lisa_promote_managed_directory(staging, destination)
  expect_identical(rename_calls, 1L)
  expect_true(file.exists(file.path(destination, "artifact.txt")))
  expect_false(lisaR:::lisa_path_entry_exists(staging))
  if (.Platform$OS.type != "windows") {
    expect_identical(
      as.integer(file.info(destination)$mode),
      as.integer(as.octmode("0700"))
    )
    expect_identical(
      as.integer(file.info(file.path(destination, "artifact.txt"))$mode),
      as.integer(as.octmode("0600"))
    )
  }
})

test_that("Failed directory-promotion rollback exposes the surviving copy", {
  parent <- tempfile("lisa-directory-recovery-")
  dir.create(parent)
  on.exit(lisa_test_cleanup_path(parent), add = TRUE)
  staging <- lisaR:::lisa_run_root(file.path(parent, "staging"))
  destination <- file.path(parent, "final")
  canonical_destination <- file.path(
    normalizePath(dirname(destination), winslash = "/", mustWork = TRUE),
    basename(destination)
  )
  writeLines("artifact", file.path(staging, "artifact.txt"))
  real_assert <- lisaR:::lisa_assert_run_tree_safe
  real_rename <- base::file.rename
  rename_calls <- 0L
  testthat::local_mocked_bindings(
    lisa_assert_run_tree_safe = function(path) {
      if (identical(lisaR:::lisa_path_key(path),
                    lisaR:::lisa_path_key(canonical_destination))) {
        stop("injected post-promotion verification failure")
      }
      real_assert(path)
    },
    .package = "lisaR"
  )
  testthat::local_mocked_bindings(
    file.rename = function(from, to) {
      rename_calls <<- rename_calls + 1L
      if (rename_calls == 2L) return(FALSE)
      real_rename(from, to)
    },
    .package = "base"
  )
  error <- tryCatch(
    lisaR:::lisa_promote_managed_directory(staging, destination),
    lisa_filesystem_recovery_error = identity
  )
  expect_s3_class(error, "lisa_filesystem_recovery_error")
  expect_identical(rename_calls, 2L)
  expect_identical(error$destination, canonical_destination)
  expect_true(canonical_destination %in% error$recovery_paths)
  expect_true(file.exists(file.path(canonical_destination, "artifact.txt")))
  expect_false(lisaR:::lisa_path_entry_exists(staging))
})
