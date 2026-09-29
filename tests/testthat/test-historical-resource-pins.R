# C6 R1: immutable public dictionary pins remain resolvable without entering
# active/default selection or trusting transform-manifest metadata.

c6_historical_variant <- function(cache, version, profile) {
  source <- lisaR:::lisa_builtin_dictionary_path("core")
  lines <- readLines(source, warn = FALSE)
  variant <- file.path(cache, paste0("variant-", gsub("[^a-z]", "-", profile), ".tsv"))
  writeLines(lines[-length(lines)], variant)
  sha <- lisaR:::lisa_sha256_file(variant)
  install_lisa_resource(
    variant, "lisa_core", version, "Homo sapiens", profile,
    "lisa_dictionary@2", "C6 historical-pin adversarial fixture", sha,
    cache_root = cache
  )
  list(
    registry = file.path(cache, "resource_registry.tsv"),
    path = variant,
    sha256 = sha
  )
}

test_that("historical authority is closed and separate from active defaults", {
  active <- lisaR:::lisa_builtin_resource_registry()
  historical <- lisaR:::lisa_historical_dictionary_registry()
  aliases <- lisaR:::lisa_historical_dictionary_aliases()

  expect_setequal(historical$logical_id, c("lisa_core", "lisa_expanded"))
  expect_true(all(historical$version == "1.0.0"))
  expect_true(all(historical$schema == "lisa_dictionary@1"))
  expect_true(all(historical$modality == "all"))
  expect_true(all(historical$.storage == "historical_builtin"))
  expect_identical(
    unname(vapply(historical$.path, lisaR:::lisa_sha256_file, character(1))),
    unname(historical$sha256)
  )
  expect_false(any(paste(historical$logical_id, historical$version) %in%
                     paste(active$logical_id, active$version)))
  expect_setequal(aliases$old_resource_id,
                  c("lisa_core@0.1.0", "lisa_expanded@0.1.0"))
  expect_true(all(aliases$new_version == "1.0.0"))

  # A version-less reference to the retired scientific namespace still
  # resolves, and it resolves forward to the current first-publication content,
  # never back to the score-carrying historical bytes.
  for (tier in c("core", "expanded")) {
    default <- lisaR:::lisa_resolve_dictionary_resource(
      paste0("lisa_", tier), species = "Homo sapiens",
      modality = "transcriptomic/genomic"
    )
    expect_identical(default$logical_id, paste0("lisa_dictionary_", tier),
                     info = tier)
    expect_identical(default$version, "1.0.0", info = tier)
    expect_identical(default$schema, "lisa_dictionary@2", info = tier)
    # The receipt keeps what was asked for alongside what was read.
    expect_identical(default$requested_logical_id, paste0("lisa_", tier),
                     info = tier)
    expect_identical(default$superseded_resource_id, paste0("lisa_", tier),
                     info = tier)
    expect_identical(
      unname(lisaR:::lisa_pipeline_resource_defaults(tier)[["dictionary_resource"]]),
      paste0("lisa_dictionary_", tier, "@1.0.0"), info = tier
    )
  }
})

test_that("transform-manifest metadata cannot steer historical resolution", {
  manifest <- lisaR:::lisa_read_resource_migration_manifest()
  expected <- lisaR:::lisa_historical_dictionary_definitions()
  expected <- expected[expected$logical_id == "lisa_core", , drop = FALSE]

  mutations <- list(
    deleted = manifest[manifest$change_kind != "supersede_transform", , drop = FALSE],
    forged = within(manifest, {
      hit <- change_kind == "supersede_transform"
      old_sha256[hit] <- new_sha256[hit]
      old_artifact[hit] <- new_artifact[hit]
      schema[hit] <- "lisa_dictionary@2"
      content_transformed[hit] <- "FALSE"
    })
  )

  for (tampered in mutations) local({
    forged_manifest <- tampered
    testthat::local_mocked_bindings(
      lisa_read_resource_migration_manifest = function(...) forged_manifest,
      .package = "lisaR"
    )
    one <- lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "1.0.0", "Homo sapiens",
      modality = "transcriptomic/genomic"
    )
    old <- lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Homo sapiens",
      modality = "transcriptomic/genomic"
    )
    expect_identical(one$schema, "lisa_dictionary@1")
    expect_identical(one$sha256, expected$sha256[[1L]])
    expect_identical(basename(one$path), "lisa_core_runtime_v0_1.tsv")
    expect_identical(old$path, one$path)
    expect_identical(old$sha256, one$sha256)
  })
})

test_that("profile overrides and historical fallback retain the C5 contract", {
  cache <- withr::local_tempdir("lisa-c6-old-override-")
  fixture <- c6_historical_variant(cache, "0.1.0", "transcriptomic/genomic")

  own <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens",
    modality = "transcriptomic/genomic", registry = fixture$registry,
    cache_root = cache
  )
  expect_identical(own$resource_id, "lisa_core@0.1.0")
  expect_identical(own$sha256, fixture$sha256)
  expect_identical(own$schema, "lisa_dictionary@2")

  other <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens",
    modality = "global proteomic", registry = fixture$registry,
    cache_root = cache
  )
  authority <- lisaR:::lisa_historical_dictionary_definitions()
  authority <- authority[authority$logical_id == "lisa_core", , drop = FALSE]
  expect_identical(other$resource_id, "lisa_core@1.0.0")
  expect_identical(other$modality, "all")
  expect_identical(other$schema, "lisa_dictionary@1")
  expect_identical(other$sha256, authority$sha256[[1L]])
})

test_that("a different-byte alias target fails with RESOURCE-030", {
  cache <- withr::local_tempdir("lisa-c6-target-shadow-")
  fixture <- c6_historical_variant(cache, "1.0.0", "transcriptomic/genomic")

  direct <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "1.0.0", "Homo sapiens",
    modality = "transcriptomic/genomic", registry = fixture$registry,
    cache_root = cache
  )
  expect_identical(direct$resource_id, "lisa_core@1.0.0")
  expect_identical(direct$sha256, fixture$sha256)

  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Homo sapiens",
      modality = "transcriptomic/genomic", registry = fixture$registry,
      cache_root = cache
    ),
    "LISA-RESOURCE-030"
  )

  fallback <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens", modality = "global proteomic",
    registry = fixture$registry, cache_root = cache
  )
  expect_identical(fallback$schema, "lisa_dictionary@1")
  expect_false(identical(fallback$sha256, fixture$sha256))
})

test_that("historical all-modality identities are reserved before installation writes", {
  definitions <- lisaR:::lisa_historical_dictionary_definitions()
  core <- definitions[definitions$logical_id == "lisa_core", , drop = FALSE]
  expect_error(
    lisaR:::lisa_active_resource_registry(core),
    "LISA-RESOURCE-030.*cannot be re-minted"
  )

  cache <- withr::local_tempdir("lisa-c6-install-reserved-")
  registry <- file.path(cache, "resource_registry.tsv")
  source <- file.path(system.file(package = "lisaR"), core$artifact[[1L]])
  expect_error(
    install_lisa_resource(
      source, "lisa_core", "1.0.0", "Homo sapiens", "all",
      "lisa_dictionary@1", "attempted historical identity replacement",
      core$sha256[[1L]], cache_root = cache, registry_path = registry
    ),
    "LISA-RESOURCE-030"
  )
  expect_false(file.exists(registry))
  expect_false(dir.exists(file.path(cache, "lisa_core", "1.0.0")))
})

test_that("historical aliases keep loud species, ambiguity and unknown-pin guards", {
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Mus musculus",
      modality = "transcriptomic/genomic"
    ),
    "LISA-RESOURCE-006"
  )
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.9.9", "Homo sapiens",
      modality = "transcriptomic/genomic"
    ),
    "LISA-RESOURCE-006"
  )

  cache <- withr::local_tempdir("lisa-c6-ambiguous-")
  c6_historical_variant(cache, "0.1.0", "transcriptomic/genomic")
  c6_historical_variant(cache, "0.1.0", "global proteomic")
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Homo sapiens",
      registry = file.path(cache, "resource_registry.tsv"), cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )
})

test_that("original public Riaz and CPTAC configs retain dictionary readiness", {
  release_root <- Sys.getenv("LISAR_PUBLIC_RELEASE_ROOT", "")
  skip_if(!nzchar(release_root), "set LISAR_PUBLIC_RELEASE_ROOT for original public-config audit")
  skip_if_not_installed("yaml")
  configs <- c(
    file.path(
      release_root, "riaz-gse91061 anonymous-project", "config",
      "riaz-gse91061.yml"
    ),
    file.path(
      release_root, "cptac-ccrcc anonymous-project", "config",
      "cptac-ccrcc-global-proteomic.yml"
    )
  )
  expect_true(all(file.exists(configs)))
  before <- vapply(configs, lisaR:::lisa_sha256_file, character(1))
  for (path in configs) {
    config <- yaml::read_yaml(path)
    expect_identical(config$pipeline$dictionary_resource, "lisa_core@0.1.0")
    expect_identical(config$pipeline$lisa_dictionary, "core")
    reference <- lisaR:::lisa_parse_resource_reference(
      config$pipeline$dictionary_resource, "dictionary_resource",
      require_version = TRUE
    )
    resource <- lisaR:::lisa_resolve_dictionary_resource(
      reference$logical_id, reference$version, "Homo sapiens",
      modality = config$pipeline$profile
    )
    expect_identical(resource$resource_id, "lisa_core@1.0.0")
    expect_identical(resource$schema, "lisa_dictionary@1")
    dictionary <- utils::read.delim(
      resource$path, sep = "\t", check.names = FALSE,
      stringsAsFactors = FALSE
    )
    expect_silent(lisaR:::lisa_validate_standard_dictionary(
      dictionary, resource$resource_id, tier = "core", schema = resource$schema
    ))
    runtime <- lisaR:::lisa_dictionary_runtime_projection(dictionary)
    expect_false("LISA_score" %in% names(runtime))
    expect_identical(names(runtime), lisaR:::lisa_dictionary_runtime_columns())

    # Match the public project's own preflight granularity: TERM2GENE may be
    # absent in an isolated cache, but the dictionary row itself must be ready.
    cache <- withr::local_tempdir("lisa-public-pin-")
    bundle <- lisaR:::lisa_resolve_pipeline_resources(
      config$pipeline, "Homo sapiens", stop_on_error = FALSE,
      cache_root = cache
    )
    dictionary_status <- bundle$status[
      bundle$status$config_key == "dictionary_resource", , drop = FALSE
    ]
    expect_identical(nrow(dictionary_status), 1L)
    expect_true(dictionary_status$ready[[1L]])
    expect_identical(dictionary_status$resource_id[[1L]], "lisa_core@1.0.0")
  }
  expect_identical(
    unname(vapply(configs, lisaR:::lisa_sha256_file, character(1))),
    unname(before)
  )
})
