test_that("scientific candidate artifacts are packaged with immutable hashes", {
  root <- system.file("extdata", "dictionaries", package = "lisaR")
  expect_true(dir.exists(root))
  manifest_path <- file.path(root, "RESOURCE_BUNDLE_MANIFEST.tsv")
  manifest <- lisaR:::lisa_read_dictionary_registry(manifest_path)
  expected <- c(
    "lisa_dictionary_core@1.0.0", "lisa_dictionary_expanded@1.0.0",
    "lisa_category_map@1.0.0"
  )
  scientific <- manifest[paste0(manifest$logical_id, "@", manifest$version) %in% expected, ]
  expect_setequal(paste0(scientific$logical_id, "@", scientific$version), expected)
  expect_false(any(grepl("term2gene", scientific$logical_id, ignore.case = TRUE)))
  paths <- file.path(system.file(package = "lisaR"), scientific$artifact)
  expect_true(all(file.exists(paths)))
  expect_identical(
    unname(vapply(paths, lisaR:::lisa_sha256_file, character(1))),
    unname(scientific$sha256)
  )
})

test_that("precomputed validation outputs and implicit semantic routes stay external", {
  expect_identical(
    system.file("extdata", "riaz-gse91061", package = "lisaR"),
    ""
  )
  expect_identical(
    system.file(
      "examples", "riaz-gse91061", "benchmarks",
      "canonical-model-benchmark-v1.tsv", package = "lisaR"
    ),
    ""
  )
  expect_identical(
    system.file(
      "examples", "riaz-gse91061", "benchmarks",
      "canonical-model-effects-v1.tsv.gz", package = "lisaR"
    ),
    ""
  )

  namespace <- asNamespace("lisaR")
  for (retired in c(
    paste0("builtin_", "extra_category_map"),
    paste0("find_pathways_", "category_metadata"),
    paste0("apply_pathways_", "category_hierarchy")
  )) {
    expect_false(exists(retired, envir = namespace, inherits = FALSE))
  }
})

test_that("resource catalogue distinguishes the bundled LISA resources", {
  root <- system.file("extdata", "dictionaries", package = "lisaR")
  manifest_path <- file.path(root, "DICTIONARY_RESOURCE_MANIFEST.tsv")
  expect_true(file.exists(manifest_path))
  manifest <- read.delim(
    manifest_path, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_true(all(c(
    "role", "resource_id", "artifact", "bundled_in_package",
    "distribution_status"
  ) %in% names(manifest)))
  expect_false(any(c("rows", "bytes", "sha256") %in% names(manifest)))

  scientific_ids <- c(
    "lisa_dictionary_core@1.0.0", "lisa_dictionary_expanded@1.0.0",
    "lisa_category_map@1.0.0"
  )
  scientific <- manifest[
    manifest$resource_id %in% scientific_ids, , drop = FALSE
  ]
  expect_setequal(scientific$resource_id, scientific_ids)
  expect_identical(nrow(scientific), 3L)
  bundled <- tolower(as.character(scientific$bundled_in_package))
  expect_true(all(bundled %in% c("true", "yes", "1")))
  expect_true(all(scientific$distribution_status == "source_terms_apply"))
  expect_true(all(file.exists(file.path(root, scientific$artifact))))
})

test_that("technical rebuild metadata distinguishes bundled and external resources", {
  script <- system.file(
    "scripts", "build_compact_dictionary_resources.R", package = "lisaR"
  )
  root <- system.file("extdata", "dictionaries", package = "lisaR")
  contract_path <- file.path(root, "TECHNICAL_REBUILD_CONTRACT.tsv")
  provenance_path <- file.path(root, "TECHNICAL_REBUILD_PROVENANCE.md")
  expect_true(all(file.exists(c(script, contract_path, provenance_path))))

  contract <- read.delim(
    contract_path, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_true(all(c(
    "role", "artifact", "resource_id", "content_location",
    "reconstructability", "distribution_status"
  ) %in% names(contract)))
  expect_false(any(c("rows", "bytes", "sha256") %in% names(contract)))
  expect_setequal(
    contract$resource_id[contract$role == "bundled_runtime_resource"],
    c("lisa_dictionary_core@1.0.0", "lisa_dictionary_expanded@1.0.0",
      "lisa_category_map@1.0.0")
  )
  external_outputs <- contract[contract$role == "external_runtime_output", ]
  expect_setequal(
    external_outputs$resource_id,
    c(
      "lisa_category_map@0.1.0"
    )
  )
  expect_true(all(external_outputs$content_location == "external_archive"))

  code <- paste(readLines(script, warn = FALSE), collapse = "\n")
  # The builder emits the score-free runtime artifacts under the names the
  # active bundle manifest pins; the v0_1 files are retained in the package as
  # history but this builder no longer writes them, so a published identity
  # cannot be overwritten by a rebuild.
  expect_match(code, "lisa_dictionary_core_runtime_v1_0.tsv", fixed = TRUE)
  expect_match(code, "lisa_dictionary_expanded_runtime_v1_0.tsv", fixed = TRUE)
  expect_false(grepl("lisa_core_runtime_v0_1.tsv", code, fixed = TRUE))
  expect_match(code, "SOURCE_FAMILY_INVENTORY.tsv", fixed = TRUE)
  expect_match(code, "retained_in_complete_v0_1_resource", fixed = TRUE)
  expect_false(grepl("restricted_identifier", code, fixed = TRUE))
  expect_false(grepl("RESTRICTED_EXCLUSIONS.tsv", code, fixed = TRUE))
  expect_match(code, "--category-map-sha256", fixed = TRUE)
  expect_false(grepl("[0-9a-f]{64}", code))
  expect_false(grepl("/home/", code, fixed = TRUE))

  active_map_script <- system.file(
    "scripts", "build_active_category_map.R", package = "lisaR"
  )
  expect_true(file.exists(active_map_script))
  active_map_code <- paste(readLines(active_map_script, warn = FALSE),
                           collapse = "\n")
  expect_match(active_map_code, "removed_zero_assignment_category",
               fixed = TRUE)
  expect_match(active_map_code, "lisa_category_map@1.0.0", fixed = TRUE)
  expect_false(grepl("/home/",
                     active_map_code, fixed = TRUE))

  provenance <- paste(
    readLines(provenance_path, warn = FALSE), collapse = "\n"
  )
  expect_match(provenance, "not bundled", ignore.case = TRUE)
})

test_that("active category-map builder removes every zero-assignment row", {
  script <- system.file(
    "scripts", "build_active_category_map.R", package = "lisaR"
  )
  expect_true(file.exists(script))
  root <- tempfile("lisa-active-map-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  write_tsv <- function(x, name) {
    path <- file.path(root, name)
    utils::write.table(
      x, path, sep = "\t", quote = FALSE, row.names = FALSE,
      col.names = TRUE, na = ""
    )
    path
  }
  dictionary <- data.frame(
    universe = "PATHWAYS", gene_set_id = c("GS_A", "GS_B"),
    category_id = c("CAT_A", "CAT_B"), stringsAsFactors = FALSE
  )
  source_map <- data.frame(
    category_id = c("CAT_A", "CAT_B", "CAT_EMPTY", "OTHER_UNCLASSIFIED"),
    display_name = c("A", "B", "Empty", "Other"),
    macrogroup_id = c("SUPER_A", "SUPER_A", "SUPER_EMPTY", "OTHER_UNCLASSIFIED"),
    macrogroup_name = c("A", "A", "Empty", "Other"),
    macrogroup_order = c(1L, 1L, 2L, 3L),
    category_order_within_macrogroup = c(1L, 2L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  core <- write_tsv(dictionary, "core.tsv")
  expanded <- write_tsv(dictionary, "expanded.tsv")
  source <- write_tsv(source_map, "source-map.tsv")
  output <- file.path(root, "active-map.tsv")
  audit <- file.path(root, "audit.tsv")
  args <- c(
    "--vanilla", script,
    "--core", core,
    "--expanded", expanded,
    "--source-map", source,
    "--source-map-sha256", lisaR:::lisa_sha256_file(source),
    "--output", output,
    "--audit", audit
  )
  command <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(), args,
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(command, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(command, collapse = "\n"))
  active <- lisaR:::read_lisa_tsv(output)
  decisions <- lisaR:::read_lisa_tsv(audit)
  expect_setequal(active$category_id, c("CAT_A", "CAT_B"))
  expect_setequal(unique(active$macrogroup_id), "SUPER_A")
  expect_setequal(
    decisions$category_id[!decisions$retained],
    c("CAT_EMPTY", "OTHER_UNCLASSIFIED")
  )
  expect_true(all(
    decisions$disposition[!decisions$retained] ==
      "removed_zero_assignment_category"
  ))
})

test_that("resource provenance preserves bundled and external scientific boundaries", {
  root <- system.file("extdata", "dictionaries", package = "lisaR")
  ledger_path <- file.path(root, "RESOURCE_PROVENANCE.tsv")
  notices_path <- system.file("THIRD_PARTY_NOTICES.md", package = "lisaR")
  expect_true(all(file.exists(c(ledger_path, notices_path))))

  ledger <- read.delim(
    ledger_path, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_true(all(c(
    "component", "source", "source_version", "content_in_package",
    "transformation", "licence_or_terms", "attribution",
    "public_release_status"
  ) %in% names(ledger)))
  synthetic <- grepl(
    "quick start|minimal|custom-resource", ledger$component, ignore.case = TRUE
  )
  scientific <- ledger[!synthetic, , drop = FALSE]
  expect_gt(nrow(scientific), 0L)
  lisa <- scientific[scientific$component == "LISA semantic classification", ]
  expect_identical(nrow(lisa), 1L)
  expect_match(lisa$content_in_package, "Exact bundled", fixed = TRUE)
  expect_match(lisa$public_release_status, "Source terms and attribution apply", fixed = TRUE)
  term2gene <- scientific[scientific$component == "Scientific TERM2GENE membership", ]
  expect_identical(nrow(term2gene), 1L)
  expect_match(term2gene$content_in_package, "None", fixed = TRUE)

  notices <- paste(readLines(notices_path, warn = FALSE), collapse = "\n")
  expect_match(notices, "source-specific attribution and conditions", fixed = TRUE)
  expect_match(notices, "TERM2GENE", ignore.case = TRUE)
})
