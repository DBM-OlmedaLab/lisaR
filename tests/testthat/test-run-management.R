test_that("transaction promotion is atomic and manifests are verified", {
  final <- tempfile("lisa-contract-")
  contract <- lisa_test_run_contract()
  tx <- lisaR:::lisa_transaction_begin(final, contract, run_id = "runmanagementrun")
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  expect_equal(lisaR:::lisa_run_state(tx$staging_dir), "staged")
  lisaR:::lisa_transaction_promote(tx)
  expect_true(dir.exists(final))
  expect_equal(verify_run(final)$gate, "PASS")
  expect_equal(lisaR:::lisa_run_state(final), "completed")
  writeLines("tampered", file.path(final, "artifact.txt"))
  expect_equal(verify_run(final)$gate, "FAIL")
})

test_that("SHA-256 uses portable known vectors without operating-system tools", {
  root <- tempfile("lisa sha vectors-")
  dir.create(root)
  paths <- file.path(root, c("empty file.txt", "ascii file.txt", "unicodé file.txt"))
  file.create(paths[[1]])
  writeBin(charToRaw("abc"), paths[[2]])
  writeBin(charToRaw(enc2utf8("café")), paths[[3]])

  expected <- c(
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
    "850f7dc43910ff890f8879c0ed26fe697c93a067ad93a7d50f466a7028a9bf4e"
  )
  old_path <- Sys.getenv("PATH", unset = NA_character_)
  on.exit(if (is.na(old_path)) Sys.unsetenv("PATH") else Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = "")
  observed <- vapply(paths, lisaR:::lisa_sha256_file, character(1))
  expect_identical(unname(observed), expected)
  expect_true(all(grepl("^[0-9a-f]{64}$", observed)))
})

test_that("SHA-256 rejects absent and non-regular paths with a stable condition", {
  root <- tempfile("lisa-sha-path-")
  dir.create(root)
  missing <- file.path(root, "absent.txt")
  capture <- function(path) {
    tryCatch(lisaR:::lisa_sha256_file(path), lisa_sha256_error = identity)
  }

  absent_error <- capture(missing)
  directory_error <- capture(root)
  expect_s3_class(absent_error, "lisa_sha256_error")
  expect_s3_class(directory_error, "lisa_sha256_error")
  expect_identical(absent_error$code, "LISA-SHA256-002")
  expect_identical(directory_error$code, "LISA-SHA256-002")
  expect_match(conditionMessage(absent_error), "LISA-SHA256-002", fixed = TRUE)
})

test_that("SHA-256 backend output is cardinality-, name-, and format-checked", {
  root <- tempfile("lisa-sha-backend-")
  dir.create(root)
  paths <- file.path(root, c("one.txt", "two.txt"))
  writeLines("one", paths[[1]])
  writeLines("two", paths[[2]])
  valid <- paste(rep("0", 64L), collapse = "")
  backends <- list(
    empty = function(paths) stats::setNames(character(), character()),
    missing = function(paths) stats::setNames(rep(valid, length(paths) - 1L), paths[-length(paths)]),
    na = function(paths) stats::setNames(rep(NA_character_, length(paths)), paths),
    malformed = function(paths) stats::setNames(rep("not-a-sha256", length(paths)), paths),
    misnamed = function(paths) stats::setNames(rep(valid, length(paths)), rev(paths))
  )

  for (backend in backends) {
    error <- tryCatch(
      lisaR:::lisa_sha256_files_batch(paths, .backend = backend),
      lisa_sha256_error = identity
    )
    expect_s3_class(error, "lisa_sha256_error")
    expect_identical(error$code, "LISA-SHA256-004")
  }
  failed <- tryCatch(
    lisaR:::lisa_sha256_file(paths[[1]], .backend = function(paths) stop("simulated")),
    lisa_sha256_error = identity
  )
  expect_identical(failed$code, "LISA-SHA256-003")
})

test_that("blank, NA, malformed, and uppercase run-manifest hashes fail closed", {
  root <- tempfile("lisa-manifest-invalid-sha-")
  dir.create(root)
  writeLines("artifact", file.path(root, "artifact.txt"))
  manifest <- lisaR:::lisa_run_manifest(root)
  invalid <- c("", NA_character_, "not-a-sha256", toupper(manifest$sha256[[1]]))
  for (hash in invalid) {
    candidate <- manifest
    candidate$sha256[[1]] <- hash
    check <- lisaR:::lisa_validate_built_run_manifest(root, candidate)
    expect_identical(check$gate, "FAIL")
    expect_true("manifest_invalid_sha256" %in% check$findings)
  }

  manifest$sha256[[1]] <- ""
  lisaR:::write_lisa_tsv(manifest, file.path(root, "run_manifest.tsv"))
  lisaR:::write_lisa_tsv(
    data.frame(event = "validated", stringsAsFactors = FALSE),
    file.path(root, "run_events.tsv")
  )
  verified <- lisaR:::verify_run(root)
  expect_identical(verified$gate, "FAIL")
  expect_true("manifest_invalid_sha256" %in% verified$findings)
})

test_that("promotion builds once and independently validates one manifest", {
  final <- tempfile("lisa-manifest-once-")
  calls <- new.env(parent = emptyenv())
  calls$manifest <- 0L
  calls$validate <- 0L
  calls$hash <- 0L
  real_manifest <- lisaR:::lisa_run_manifest
  real_validate <- lisaR:::lisa_validate_built_run_manifest
  real_hash <- lisaR:::lisa_hash_files
  testthat::local_mocked_bindings(
    lisa_run_manifest = function(...) {
      calls$manifest <- calls$manifest + 1L
      real_manifest(...)
    },
    lisa_validate_built_run_manifest = function(...) {
      calls$validate <- calls$validate + 1L
      real_validate(...)
    },
    lisa_hash_files = function(...) {
      calls$hash <- calls$hash + 1L
      real_hash(...)
    },
    verify_run = function(...) stop("promotion must not call external verification"),
    .package = "lisaR"
  )
  tx <- lisaR:::lisa_transaction_begin(final, lisa_test_run_contract(), run_id = "manifestonce")
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  lisaR:::lisa_transaction_promote(tx)
  expect_identical(calls$manifest, 1L)
  expect_identical(calls$validate, 1L)
  expect_identical(calls$hash, 2L)
  expect_true(dir.exists(final))
})

test_that("only the serial promotion coordinator removes shared diagnostics", {
  final <- tempfile("lisa-diagnostics-finalize-")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "diagnosticsfinalize"
  )
  diagnostics <- file.path(tx$staging_dir, ".lisa_subprocess")
  dir.create(diagnostics, showWarnings = FALSE)
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  lisaR:::lisa_transaction_promote(tx)
  expect_false(dir.exists(file.path(final, ".lisa_subprocess")))

  blocked <- tempfile("lisa-diagnostics-blocked-")
  tx_blocked <- lisaR:::lisa_transaction_begin(
    blocked, lisa_test_run_contract(), run_id = "diagnosticsblocked"
  )
  active <- file.path(tx_blocked$staging_dir, ".lisa_subprocess", "worker")
  dir.create(active, recursive = TRUE)
  writeLines("still active", file.path(active, "stdout.txt"))
  expect_error(
    lisaR:::lisa_transaction_promote(tx_blocked),
    "diagnostics remained active"
  )
  expect_false(dir.exists(blocked))
})

test_that("external verification recalculates a completed manifest once", {
  final <- tempfile("lisa-verify-once-")
  tx <- lisaR:::lisa_transaction_begin(final, lisa_test_run_contract(), run_id = "verifyonce")
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  lisaR:::lisa_transaction_promote(tx)
  calls <- new.env(parent = emptyenv())
  calls$manifest <- 0L
  real_manifest <- lisaR:::lisa_run_manifest
  testthat::local_mocked_bindings(
    lisa_run_manifest = function(...) {
      calls$manifest <- calls$manifest + 1L
      real_manifest(...)
    },
    .package = "lisaR"
  )
  expect_identical(verify_run(final)$gate, "PASS")
  expect_identical(calls$manifest, 1L)
})

test_that("promotion manifest validation detects same-size content changes", {
  root <- tempfile("lisa-manifest-content-")
  dir.create(root)
  artifact <- file.path(root, "artifact.txt")
  writeLines("alpha", artifact, useBytes = TRUE)
  manifest <- lisaR:::lisa_run_manifest(root)
  writeLines("bravo", artifact, useBytes = TRUE)
  check <- lisaR:::lisa_validate_built_run_manifest(root, manifest)
  expect_identical(check$gate, "FAIL")
  expect_match(paste(check$findings, collapse = " "), "modified:artifact.txt", fixed = TRUE)
})

test_that("promotion failure remains failed without a final directory", {
  final <- tempfile("lisa-contract-")
  tx <- lisaR:::lisa_transaction_begin(final, lisa_test_run_contract(), run_id = "renamefailure")
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  # Fail the final guarded directory promotion without also disabling the
  # atomic file renames now used to write manifests and failure events.
  testthat::local_mocked_bindings(
    lisa_promote_managed_directory = function(...) {
      stop("simulated promotion failure")
    },
    .package = "lisaR"
  )

  expect_error(lisaR:::lisa_transaction_promote(tx), "atomic promotion failed")
  expect_false(dir.exists(final))
  expect_equal(lisaR:::lisa_run_state(tx$staging_dir), "failed")
  expect_false(dir.exists(tx$lock))
})

test_that("aborted staging is preserved and cannot be resumed or reused", {
  final <- tempfile("lisa-contract-")
  contract <- lisa_test_run_contract(environment = "test")
  abandoned <- lisaR:::lisa_transaction_begin(
    final, contract, run_id = "abandonedrun"
  )
  obsolete <- file.path(abandoned$staging_dir, "obsolete.tsv")
  writeLines("preserve for inspection", obsolete)

  expect_error(
    lisaR:::lisa_transaction_begin(final, contract, run_id = "other"),
    "exclusive run lock"
  )
  lisaR:::lisa_transaction_abort(abandoned, "simulated interruption")

  same_id_error <- tryCatch(
    lisaR:::lisa_transaction_begin(
      final, contract, run_id = "abandonedrun", resume = FALSE
    ),
    error = identity
  )
  expect_match(conditionMessage(same_id_error), "LISA-RUN-002", fixed = TRUE)
  expect_true(file.exists(obsolete))
  expect_identical(readLines(obsolete), "preserve for inspection")

  resume_error <- tryCatch(
    lisaR:::lisa_transaction_begin(
      final, contract, run_id = "abandonedrun", resume = TRUE
    ),
    lisa_resume_error = identity
  )
  expect_s3_class(resume_error, "lisa_resume_error")
  expect_s3_class(resume_error, "lisa_error")
  expect_identical(resume_error$code, "LISA-RESUME-001")
  expect_match(conditionMessage(resume_error), "LISA-RESUME-001", fixed = TRUE)
  expect_true(file.exists(obsolete))

  retry <- lisaR:::lisa_transaction_begin(final, contract)
  on.exit(lisaR:::lisa_transaction_abort(retry, "test cleanup"), add = TRUE)
  expect_false(identical(retry$run_id, abandoned$run_id))
  expect_false(identical(retry$staging_dir, abandoned$staging_dir))
  expect_false(file.exists(file.path(retry$staging_dir, "obsolete.tsv")))
  expect_true(file.exists(obsolete))
  expect_identical(readLines(obsolete), "preserve for inspection")
})

test_that("resume false remains valid while resume true fails with a stable condition", {
  config <- list(pipeline = list(
    schema_version = "1.0.0", profile = "targeted",
    evidence_mode = "full_de", resume = FALSE,
    duplicate_policies = list(
      de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error",
      mapped_id_collision_policy = "error"
    )
  ))
  expect_identical(
    lisaR:::lisa_validate_pipeline_config(config)$schema_version,
    "1.0.0"
  )

  config$pipeline$resume <- TRUE
  resume_error <- tryCatch(
    lisaR:::lisa_validate_pipeline_config(config),
    lisa_resume_error = identity
  )
  expect_s3_class(resume_error, "lisa_resume_error")
  expect_identical(resume_error$code, "LISA-RESUME-001")
  expect_match(conditionMessage(resume_error), "start a new run", fixed = TRUE)
})

test_that("an existing final output blocks before the scientific engine", {
  root <- tempfile("lisa-final-collision-")
  dir.create(root)
  final <- file.path(root, "completed-run")
  dir.create(final)
  sentinel <- file.path(final, "sentinel.txt")
  writeLines("immutable", sentinel)
  config_path <- file.path(root, "study.json")
  jsonlite::write_json(list(pipeline = list(
    schema_version = "1.0.0", profile = "targeted",
    evidence_mode = "full_de", output_dir = final, resume = FALSE,
    duplicate_policies = list(
      de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error",
      mapped_id_collision_policy = "error"
    )
  )), config_path, auto_unbox = TRUE)

  engine_called <- FALSE
  testthat::local_mocked_bindings(
    run_lisa_pipeline = function(...) {
      engine_called <<- TRUE
      stop("scientific engine must not be called")
    },
    .package = "lisaR"
  )

  expect_error(
    lisaR:::run_lisa_pipeline_from_config(config_path),
    "LISA-RUN-005"
  )
  expect_false(engine_called)
  expect_identical(readLines(sentinel), "immutable")
  expect_false(dir.exists(paste0(final, ".lisa.lock")))
})

test_that("product planning prevents unsupported products from computing", {
  plan <- lisaR:::lisa_product_plan(lisaR:::lisa_evidence_contract("rank_only"), run_gene_level = TRUE, run_kegg_maps = TRUE)
  expect_true(plan$compute[plan$product == "ranked_enrichment"])
  expect_false(any(plan$compute[plan$product %in% c("gene_level", "kegg_maps")]))
})

test_that("failed execution stages cannot be promoted as PASS", {
  failed <- data.frame(
    analysis_id = "A", collection = "GOBP-C2", status = "failed",
    message = "simulated failure", stringsAsFactors = FALSE
  )
  expect_error(
    lisaR:::lisa_assert_no_failed_stages(failed),
    "LISA-RUN-FAILED.*simulated failure"
  )
  expect_true(lisaR:::lisa_assert_no_failed_stages(
    data.frame(status = "completed", stringsAsFactors = FALSE)
  ))
})

test_that("configured rank_col reaches the single-DE engine", {
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) {
      captured$args <- list(...)
      list(status = "completed")
    },
    .package = "lisaR"
  )
  de <- data.frame(
    analysis_id = "A", de_path = "input.tsv", species = "Homo sapiens",
    rank_col = "custom_score", stringsAsFactors = FALSE
  )
  registry <- data.frame(
    analysis_collection = "GOBP-C2", dictionary_id = "LISA_GOBP_C2",
    run_ora = FALSE, run_category_pathway_plots = FALSE,
    post_lisa_profile = "full", allow_missing_dictionary = FALSE,
    stringsAsFactors = FALSE
  )
  lisaR:::run_lisa_pipeline_single_de(
    de, registry, getwd(), tempfile("single-output-"), "term2gene.tsv",
    "dictionary.tsv", "category-map.tsv", NULL, "core", "semantic",
    "png", "tsv", registered_category_map = TRUE
  )
  expect_identical(captured$args$rank_col, "custom_score")
  expect_true(captured$args$registered_category_map)
})

test_that("run identity validation fails closed on ambiguous identifiers", {
  idx <- data.frame(analysis_id = c("A", "A"), de_path = c("x", "y"), stringsAsFactors = FALSE)
  expect_error(lisaR:::lisa_validate_run_identity(idx), "duplicated|duplicate")
})

test_that("verification detects absent and stale outputs, and low space fails before staging", {
  final <- tempfile("lisa-contract-")
  tx <- lisaR:::lisa_transaction_begin(final, lisa_test_run_contract(), run_id = "integrity")
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  lisaR:::lisa_transaction_promote(tx)
  writeLines("stale", file.path(final, "stale.txt"))
  expect_match(paste(verify_run(final)$findings, collapse = " "), "undeclared:stale.txt", fixed = TRUE)
  unlink(file.path(final, "artifact.txt"))
  expect_match(paste(verify_run(final)$findings, collapse = " "), "absent:artifact.txt", fixed = TRUE)
  old <- getOption("lisaR.staging_space_check")
  options(lisaR.staging_space_check = function(path) FALSE)
  on.exit(options(lisaR.staging_space_check = old), add = TRUE)
  expect_error(lisaR:::lisa_transaction_begin(tempfile("lisa-low-space-"), lisa_test_run_contract()), "insufficient staging disk space")
})

test_that("Invalid staging entries release their transaction locks", {
  base <- tempfile("lisa-invalid-staging-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)
  contract <- lisa_test_run_contract()

  final_file <- file.path(base, "file-final")
  staging_file <- lisaR:::lisa_short_staging_path(final_file, "stagefile")
  writeLines("not a directory", staging_file)
  expect_error(
    lisaR:::lisa_transaction_begin(
      final_file, contract, run_id = "stagefile", resume = TRUE
    ),
    "LISA-RESUME-001"
  )
  expect_false(lisaR:::lisa_path_entry_exists(paste0(final_file, ".lisa.lock")))

  final_link <- file.path(base, "link-final")
  staging_link <- lisaR:::lisa_short_staging_path(final_link, "stagelink")
  if (!isTRUE(suppressWarnings(file.symlink(
    file.path(base, "absent"), staging_link
  )))) {
    skip("file.symlink() is unavailable in this test environment")
  }
  expect_error(
    lisaR:::lisa_transaction_begin(
      final_link, contract, run_id = "stagelink", resume = TRUE
    ),
    "LISA-RESUME-001"
  )
  expect_false(lisaR:::lisa_path_entry_exists(paste0(final_link, ".lisa.lock")))
})

test_that("Transaction abort does not recreate absent staging or output", {
  final <- tempfile("lisa-abort-absent-")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "abortabsent"
  )
  unlink(tx$staging_dir, recursive = TRUE, force = TRUE)
  expect_false(lisaR:::lisa_path_entry_exists(tx$staging_dir))
  lisaR:::lisa_transaction_abort(tx, "already removed")
  expect_false(lisaR:::lisa_path_entry_exists(tx$staging_dir))
  expect_false(lisaR:::lisa_path_entry_exists(final))
  expect_false(lisaR:::lisa_path_entry_exists(tx$lock))
})

test_that("Transaction locks are private directories", {
  final <- tempfile("lisa-private-lock-")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "privatelock"
  )
  on.exit(lisaR:::lisa_transaction_abort(tx), add = TRUE)
  expect_true(dir.exists(tx$lock))
  if (.Platform$OS.type != "windows") {
    expect_identical(as.integer(file.info(tx$lock)$mode),
                     as.integer(as.octmode("0700")))
  }
})

test_that("Promotion recovery conditions survive the transaction boundary", {
  final <- tempfile("lisa-recovery-propagation-")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "recoverypropagation"
  )
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  staged <- tx$staging_dir
  testthat::local_mocked_bindings(
    lisa_promote_managed_directory = function(staging, destination, ...) {
      lisaR:::lisa_filesystem_recovery_abort(
        "simulated directory rollback failure",
        source = staging, destination = destination, staged = staging,
        recovery_paths = staging
      )
    },
    .package = "lisaR"
  )
  error <- tryCatch(
    lisaR:::lisa_transaction_promote(tx),
    lisa_filesystem_recovery_error = identity
  )
  expect_s3_class(error, "lisa_filesystem_recovery_error")
  expect_identical(error$staged, staged)
  expect_true(dir.exists(staged))
  expect_false(lisaR:::lisa_path_entry_exists(tx$lock))
})

test_that("Failed install rollback retains every recovery copy and its paths", {
  root <- lisaR:::lisa_run_root(tempfile("lisa-install-recovery-"))
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  destination <- file.path(root, "artifact.txt")
  writeLines("original", destination, useBytes = TRUE)
  real_hash <- lisaR:::lisa_sha256_file
  real_rename <- base::file.rename
  hash_calls <- 0L
  rename_calls <- 0L
  testthat::local_mocked_bindings(
    lisa_sha256_file = function(path, ...) {
      hash_calls <<- hash_calls + 1L
      if (hash_calls == 2L) return(paste(rep("0", 64L), collapse = ""))
      real_hash(path, ...)
    },
    .package = "lisaR"
  )
  testthat::local_mocked_bindings(
    file.rename = function(from, to) {
      rename_calls <<- rename_calls + 1L
      if (rename_calls == 4L) return(FALSE)
      real_rename(from, to)
    },
    .package = "base"
  )
  error <- tryCatch(
    lisaR:::lisa_guarded_write(
      destination, function(path) writeLines("replacement", path), root
    ),
    lisa_filesystem_recovery_error = identity
  )
  expect_s3_class(error, "lisa_filesystem_recovery_error")
  expect_identical(error$code, "LISA-FS-RECOVERY-001")
  expect_true(all(c("source", "destination", "backup", "staged",
                    "recovery_paths") %in% names(error)))
  expect_gte(length(error$recovery_paths), 2L)
  expect_true(all(vapply(error$recovery_paths,
                         lisaR:::lisa_path_entry_exists, logical(1))))
  surviving_text <- sort(unlist(lapply(
    error$recovery_paths[file.exists(error$recovery_paths) &
                           !dir.exists(error$recovery_paths)],
    readLines, warn = FALSE
  )))
  expect_true(all(c("original", "replacement") %in% surviving_text))
})
