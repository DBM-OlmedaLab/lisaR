# C6 runtime score-removal contract.
#
# These tests assert the properties that make removing `LISA_score` from the
# runtime contract safe rather than merely tidy:
#   * the premise that the column was a redundant encoding of `tier`;
#   * that the new resources are an EXACT projection of the old ones;
#   * that historical bytes are preserved and served only for their exact old
#     identities, never substituted with the new score-free bytes;
#   * that both dictionary schemas remain readable and produce identical runtime
#     data (the transitional reader);
#   * that a new custom resource needs no empty score column.
#

c6_read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = character()
  )
}

c6_artifact <- function(relative) {
  file.path(system.file(package = "lisaR"), relative)
}

c6_transform_rows <- function() {
  migration <- lisaR:::lisa_read_resource_migration_manifest()
  migration[migration$change_kind == "supersede_transform", , drop = FALSE]
}

test_that("the runtime dictionary contract no longer contains LISA_score", {
  runtime <- lisaR:::lisa_dictionary_runtime_columns()
  expect_identical(
    runtime,
    c("universe", "gene_set_id", "gene_set_name", "source_id",
      "category_id", "category_display_name", "tier")
  )
  expect_false("LISA_score" %in% runtime)
  # tier survives: it is the column that actually carries the consensus level.
  expect_true("tier" %in% runtime)

  # Both schemas remain declarable; the legacy one is not deleted from history.
  expect_identical(
    lisaR:::lisa_dictionary_schemas(),
    c("lisa_dictionary@1", "lisa_dictionary@2")
  )
  expect_true(all(
    lisaR:::lisa_dictionary_schemas() %in% lisaR:::lisa_resource_schemas()
  ))
})

test_that("the active dictionaries are score-free and registered as such", {
  builtin <- lisaR:::lisa_builtin_resource_registry()
  for (logical_id in c("lisa_dictionary_core", "lisa_dictionary_expanded")) {
    rows <- builtin[builtin$logical_id == logical_id, , drop = FALSE]
    expect_identical(nrow(rows), 1L, info = logical_id)
    expect_identical(rows$version[[1L]], "1.0.0", info = logical_id)
    expect_identical(rows$schema[[1L]], "lisa_dictionary@2", info = logical_id)

    table <- c6_read_tsv(c6_artifact(rows$artifact[[1L]]))
    expect_false("LISA_score" %in% names(table), info = logical_id)
    expect_identical(names(table), lisaR:::lisa_dictionary_runtime_columns(),
                     info = logical_id)
  }
  # The category map never carried a score and must not have moved.
  map <- builtin[builtin$logical_id == "lisa_category_map", , drop = FALSE]
  expect_identical(map$version[[1L]], "1.0.0")
  expect_identical(map$schema[[1L]], "category_map@1")
})

test_that("the new resources are an exact projection of the retired ones", {
  transform <- c6_transform_rows()
  expect_identical(nrow(transform), 2L)

  for (i in seq_len(nrow(transform))) {
    id <- transform$logical_id[[i]]
    old <- c6_read_tsv(c6_artifact(transform$old_artifact[[i]]))
    new <- c6_read_tsv(c6_artifact(transform$new_artifact[[i]]))

    # Exactly one column removed, and it is the expected one.
    expect_true("LISA_score" %in% names(old), info = id)
    expect_false("LISA_score" %in% names(new), info = id)
    expect_identical(names(new), setdiff(names(old), "LISA_score"), info = id)

    # No row added, removed or reordered, and every retained cell identical.
    expect_identical(nrow(new), nrow(old), info = id)
    for (column in names(new)) {
      expect_identical(new[[column]], old[[column]],
                       info = paste(id, column, sep = "/"))
    }

    # The premise the whole migration rests on: the score was a function of
    # tier, so dropping it discarded no information.
    expected <- ifelse(old$tier == "core", "4", "3")
    expect_identical(old$LISA_score, expected, info = id)

    # The projection did not merge rows, and the assignment key stays unique.
    expect_identical(anyDuplicated(new), 0L, info = id)
    expect_identical(
      anyDuplicated(paste(new$universe, new$gene_set_id, new$category_id, sep = "\r")),
      0L, info = id
    )
  }
})

test_that("core remains contained in expanded after the projection", {
  builtin <- lisaR:::lisa_builtin_resource_registry()
  path_for <- function(id) {
    c6_artifact(builtin$artifact[builtin$logical_id == id][[1L]])
  }
  core <- c6_read_tsv(path_for("lisa_dictionary_core"))
  expanded <- c6_read_tsv(path_for("lisa_dictionary_expanded"))

  key <- function(x) do.call(paste, c(x[lisaR:::lisa_dictionary_runtime_columns()], sep = "\r"))
  expect_true(all(key(core) %in% key(expanded)))
  # And the containment is over assignments, not merely over gene sets.
  assignment <- function(x) paste(x$universe, x$gene_set_id, x$category_id, sep = "\r")
  expect_true(all(assignment(core) %in% assignment(expanded)))
})

test_that("historical artifacts are preserved byte for byte", {
  transform <- c6_transform_rows()
  for (i in seq_len(nrow(transform))) {
    old_path <- c6_artifact(transform$old_artifact[[i]])
    expect_true(file.exists(old_path), info = transform$logical_id[[i]])
    expect_identical(
      lisaR:::lisa_sha256_file(old_path), transform$old_sha256[[i]],
      info = transform$logical_id[[i]]
    )
    # The new artifact is a different file, not an overwrite of the old one.
    expect_false(identical(transform$old_artifact[[i]], transform$new_artifact[[i]]))
    expect_identical(
      lisaR:::lisa_sha256_file(c6_artifact(transform$new_artifact[[i]])),
      transform$new_sha256[[i]], info = transform$logical_id[[i]]
    )
    expect_false(identical(transform$old_sha256[[i]], transform$new_sha256[[i]]))
  }
})

test_that("transformed identities stay outside the active/default registry", {
  transform <- c6_transform_rows()
  expect_true(all(transform$content_transformed == "TRUE"))
  expect_true(all(transform$old_state == "retired_transformed"))

  # The generic relabel ledger still excludes transform rows: historical
  # dictionary resolution has its own code-authenticated identity authority.
  superseded <- lisaR:::lisa_superseded_resource_identities()
  expect_false(any(transform$old_resource_id %in% superseded$old_resource_id))

  # Historical versions are not active registry rows, which keeps defaults at
  # 2.0.0 and unversioned selection deterministic.
  builtin <- lisaR:::lisa_builtin_resource_registry()
  registry_key <- paste(builtin$logical_id, builtin$version, sep = "@")
  expect_false(any(transform$old_resource_id %in% registry_key))
})

test_that("exact historical pins resolve their original bytes and schema", {
  cache <- withr::local_tempdir("lisa-c6-historical-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  for (tier in c("core", "expanded")) {
    historical <- lisaR:::lisa_resolve_dictionary_resource(
      paste0("lisa_", tier), "1.0.0", "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
    authority <- lisaR:::lisa_historical_dictionary_definitions()
    authority <- authority[authority$logical_id == paste0("lisa_", tier), , drop = FALSE]
    expect_identical(historical$resource_id, paste0("lisa_", tier, "@1.0.0"),
                     info = tier)
    expect_identical(historical$requested_resource_id, historical$resource_id,
                     info = tier)
    expect_identical(historical$schema, "lisa_dictionary@1", info = tier)
    expect_identical(historical$sha256, authority$sha256[[1L]], info = tier)
    expect_identical(lisaR:::lisa_sha256_file(historical$path), historical$sha256,
                     info = tier)
    expect_true("LISA_score" %in% names(c6_read_tsv(historical$path)), info = tier)
    active <- lisaR:::lisa_resolve_dictionary_resource(
      paste0("lisa_", tier), "2.0.0", "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
    expect_false(identical(historical$path, active$path), info = tier)
    expect_false(identical(historical$sha256, active$sha256), info = tier)
  }
  # The public pre-C5 alias reaches the same exact historical file, not 2.0.0.
  old <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens",
    modality = "transcriptomic/genomic", cache_root = cache
  )
  one <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "1.0.0", "Homo sapiens",
    modality = "transcriptomic/genomic", cache_root = cache
  )
  expect_identical(old$path, one$path)
  expect_identical(old$sha256, one$sha256)
  expect_identical(old$schema, one$schema)
  expect_identical(old$requested_resource_id, "lisa_core@0.1.0")
  expect_identical(old$superseded_resource_id, "lisa_core@0.1.0")
})

test_that("the transitional reader accepts both dictionary schemas", {
  legacy <- g2_dictionary(tier = "custom", legacy_score = TRUE)
  modern <- g2_dictionary(tier = "custom")

  expect_true("LISA_score" %in% names(legacy))
  expect_false("LISA_score" %in% names(modern))

  expect_silent(lisaR:::lisa_validate_custom_dictionary(legacy, "legacy@1.0.0"))
  expect_silent(lisaR:::lisa_validate_custom_dictionary(modern, "modern@1.0.0"))

  # Both normalise to the same runtime table, which is the whole point: an old
  # resource and the equivalent new one cannot produce different science.
  expect_identical(
    lisaR:::lisa_dictionary_runtime_projection(legacy),
    lisaR:::lisa_dictionary_runtime_projection(modern)
  )
  # Projection is idempotent and leaves a modern table untouched.
  expect_identical(lisaR:::lisa_dictionary_runtime_projection(modern), modern)
})

test_that("a standard dictionary is validated on tier, legacy score included", {
  core <- g2_dictionary(tier = "core")
  expect_silent(
    lisaR:::lisa_validate_standard_dictionary(core, "lisa_core@2.0.0", tier = "core")
  )
  # The historical coherence is still enforced when the column is present, so a
  # legacy resource is never checked less strictly than before C6.
  legacy_ok <- g2_dictionary(tier = "core", score = 4, legacy_score = TRUE)
  expect_silent(
    lisaR:::lisa_validate_standard_dictionary(legacy_ok, "lisa_core@1.0.0", tier = "core")
  )
  legacy_bad <- g2_dictionary(tier = "core", score = 3, legacy_score = TRUE)
  expect_error(
    lisaR:::lisa_validate_standard_dictionary(legacy_bad, "lisa_core@1.0.0", tier = "core"),
    "LISA-RESOURCE-024.*tier/score coherence"
  )
  # A reserved tier is still refused for a core resource.
  expect_error(
    lisaR:::lisa_validate_standard_dictionary(
      g2_dictionary(tier = "expanded"), "lisa_core@2.0.0", tier = "core"
    ),
    "LISA-RESOURCE-024.*core tier contract"
  )
})

test_that("a declared schema must match the table it describes", {
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      g2_dictionary(legacy_score = TRUE), "x@1.0.0", schema = "lisa_dictionary@2"
    ),
    "LISA-CUSTOM-001.*carries a LISA_score"
  )
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      g2_dictionary(), "x@1.0.0", schema = "lisa_dictionary@1"
    ),
    "LISA-CUSTOM-001.*omits a LISA_score"
  )
})

test_that("annotation carries every dictionary column except the removed score", {
  dictionary <- g2_dictionary(tier = "core")
  dictionary$source_family <- "GO:BP"
  gsea <- data.frame(
    pathway = c("GS_A", "GS_B"), NES = c(1.5, -2.0), padj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  annotated <- lisaR:::annotate_lisa(gsea, dictionary, by_col = "pathway")

  expect_false("LISA_score" %in% names(annotated))
  for (column in c("universe", "gene_set_name", "source_id", "category_id",
                   "category_display_name", "tier", "source_family")) {
    expect_true(column %in% names(annotated), info = column)
  }
  # Row multiplicity is decided by the gene set join, not by the score: GS_A is
  # assigned to two categories and must still produce two annotated rows.
  expect_identical(sum(annotated$pathway == "GS_A"), 2L)
  expect_identical(sum(annotated$pathway == "GS_B"), 1L)
  # The statistics are carried through untouched.
  expect_identical(sort(unique(annotated$NES)), c(-2.0, 1.5))
})

test_that("the direct HALLMARKS collection needs no score column", {
  term2gene <- data.frame(
    gs_collection = "H", gs_subcollection = "",
    gs_name = c("HALLMARK_APOPTOSIS", "HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA"),
    gs_exact_source = "", gene_symbol = c("CASP3", "BAX", "VEGFA"),
    stringsAsFactors = FALSE
  )
  built <- lisaR:::build_hallmarks_lisa_inputs(term2gene)
  expect_false("LISA_score" %in% names(built$lisa_dict))
  expect_true(all(built$lisa_dict$tier == "msigdb_direct"))
  expect_setequal(
    built$lisa_dict$gene_set_id, c("HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA")
  )
})

test_that("bundled fixtures expose the score-free contract", {
  quick <- c6_read_tsv(system.file(
    "extdata", "quick-start", "lisa_dictionary_quickstart_v1_0.tsv",
    package = "lisaR", mustWork = TRUE
  ))
  expect_false("LISA_score" %in% names(quick))

  custom <- c6_read_tsv(system.file(
    "extdata", "custom-resources", "lisa_dictionary_custom_example_v1_0.tsv",
    package = "lisaR", mustWork = TRUE
  ))
  expect_false("LISA_score" %in% names(custom))
  expect_true(all(custom$tier == "custom"))

  # The legacy fixtures are retained so the old format keeps being exercised.
  legacy <- c6_read_tsv(system.file(
    "extdata", "custom-resources", "example_custom_dictionary.tsv",
    package = "lisaR", mustWork = TRUE
  ))
  expect_true("LISA_score" %in% names(legacy))
})
