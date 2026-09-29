test_that("installed bundled resources resolve directly, immutably and idempotently", {
  cache <- tempfile("lisa-empty-cache-")
  dir.create(cache)
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))

  expected <- data.frame(
    logical_id = c("lisa_core", "lisa_expanded", "lisa_category_map"),
    version = c("2.0.0", "2.0.0", "1.0.0"),
    schema = c("lisa_dictionary@2", "lisa_dictionary@2", "category_map@1"),
    stringsAsFactors = FALSE
  )
  resolved <- lapply(seq_len(nrow(expected)), function(i) {
    lisaR:::lisa_resolve_dictionary_resource(
      expected$logical_id[[i]], expected$version[[i]], "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
  })
  expect_true(all(vapply(resolved, function(x) file.exists(x$path), logical(1))))
  expect_identical(vapply(resolved, `[[`, character(1), "schema"), expected$schema)
  expect_false(file.exists(file.path(cache, "resource_registry.tsv")))
  expect_identical(
    lapply(seq_len(nrow(expected)), function(i) lisaR:::lisa_resolve_dictionary_resource(
      expected$logical_id[[i]], expected$version[[i]], "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )$path),
    lapply(resolved, `[[`, "path")
  )
})

test_that("bundled scientific resources are profile-independent, including global proteomic", {
  cache <- tempfile("lisa-empty-proteomic-cache-")
  dir.create(cache)
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  bundled <- lisaR:::lisa_builtin_resource_registry()
  selected <- bundled[bundled$logical_id %in% c(
    "lisa_dictionary_core", "lisa_dictionary_expanded", "lisa_category_map"
  ), , drop = FALSE]
  expect_equal(nrow(selected), 3L)
  expect_true(all(selected$modality == "all"))
  for (profile in c("transcriptomic/genomic", "global proteomic")) {
    resolved <- lapply(seq_len(nrow(selected)), function(i) {
      lisaR:::lisa_resolve_dictionary_resource(
        selected$logical_id[[i]], selected$version[[i]], "Homo sapiens",
        modality = profile, cache_root = cache
      )
    })
    expect_true(all(vapply(resolved, function(x) x$modality == "all", logical(1))))
    expect_true(all(vapply(resolved, function(x) file.exists(x$path), logical(1))))
  }
})

test_that("external registries cannot replace immutable bundled identities", {
  builtin <- lisaR:::lisa_builtin_resource_registry()
  collision <- builtin[builtin$logical_id == "lisa_dictionary_core" &
                         builtin$version == "1.0.0", ]
  external <- collision[, lisaR:::lisa_dictionary_registry_columns(), drop = FALSE]
  expect_error(
    lisaR:::lisa_active_resource_registry(external),
    "LISA-RESOURCE-020.*cannot be hidden or replaced"
  )
})
