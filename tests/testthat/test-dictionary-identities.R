# First-publication dictionary identities.
#
# The current contents are the first LISA dictionaries ever published, so they
# carry version 1.0.0 under their own logical namespace. These tests guard the
# two properties that make that safe:
#
#   1. nothing scientific moved - the same bytes, schema and tier semantics are
#      served under the new names; and
#   2. no previously published identity was re-minted - every old pin still
#      reads exactly what it always read, and the new alias can never redirect
#      a historical one.
#
# The negative cases matter as much as the positive ones: an alias that
# silently answered `lisa_core@1.0.0` would be a regression even though every
# "does it resolve?" assertion would still pass.

e1_ids <- function() {
  c(
    core = "lisa_dictionary_core@1.0.0",
    expanded = "lisa_dictionary_expanded@1.0.0",
    quickstart = "lisa_dictionary_quickstart@1.0.0",
    custom_example = "lisa_dictionary_custom_example@1.0.0"
  )
}

# The exact scientific bytes this release publishes. Hard-coded on purpose: the
# point of the whole change is that these digests do NOT move, so reading them
# back out of the manifest under test would prove nothing.
e1_hashes <- function() {
  c(
    core = "f24b5bd8d9ac66b6d64c1a91c56ab567aa64a37f05dee5ce0ed3011b2bed8994",
    expanded = "0914c2d2e40369041a2313e913898f9c8e4065779e3973b076b70ea473045736",
    quickstart = "72f36fe1a82cece837435e09939450e7bf8b90e57d179c218fd907c909235d0d",
    custom_example = "2f1b4161824e700f8252374dd3f878088ae89f54d7f384d16a712000e1de7026",
    category_map = "d61fcb2e1d40d0448d459daaab952975477203bb00b76cb53145aea8f61333b1"
  )
}

e1_empty_cache <- function() {
  cache <- withr::local_tempdir("lisa-e1-", .local_envir = parent.frame())
  withr::local_options(
    list(
      lisaR.dictionary_registry = NULL,
      lisaR.dictionary_cache_root = cache,
      lisaR.shared_dictionary_root = NULL
    ),
    .local_envir = parent.frame()
  )
  cache
}

test_that("the copied sample project and current JSON resolve the score-free v1 dictionary", {
  cache <- e1_empty_cache()
  scratch <- withr::local_tempdir("lisa-e1-sample-")
  project <- lisa_init_project(file.path(scratch, "project"))
  json_path <- system.file("examples", "minimal-study.json", package = "lisaR")
  expect_true(nzchar(json_path))
  configurations <- list(
    copied_yaml = yaml::read_yaml(file.path(project, "study.yml")),
    installed_json = jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  )
  for (name in names(configurations)) {
    pipeline <- configurations[[name]]$pipeline
    expect_identical(pipeline$dictionary_resource,
                     "lisa_dictionary_quickstart@1.1.0", info = name)
    reference <- lisaR:::lisa_parse_resource_reference(
      pipeline$dictionary_resource, "dictionary_resource", require_version = TRUE
    )
    resolved <- lisaR:::lisa_resolve_dictionary_resource(
      reference$logical_id, reference$version, "Homo sapiens",
      modality = pipeline$profile, cache_root = cache
    )
    expect_identical(resolved$schema, "lisa_dictionary@2", info = name)
    expect_identical(resolved$sha256, "09f1b17f89d81e861b533db193e4d2af200084ffffb974e0ec0c97f50ec6db84",
                     info = name)
    expect_identical(lisaR:::lisa_sha256_file(resolved$path),
                     "09f1b17f89d81e861b533db193e4d2af200084ffffb974e0ec0c97f50ec6db84", info = name)
    dictionary <- utils::read.delim(resolved$path, check.names = FALSE)
    expect_identical(names(dictionary), lisaR:::lisa_dictionary_runtime_columns(),
                     info = name)
    expect_false("LISA_score" %in% names(dictionary), info = name)
  }
})

test_that("the active scientific defaults are the first-publication v1 identities", {
  ids <- e1_ids()
  expect_identical(
    unname(lisaR:::lisa_dictionary_tier_resources()),
    unname(ids[c("core", "expanded")])
  )
  for (tier in c("core", "expanded")) {
    defaults <- lisaR:::lisa_pipeline_resource_defaults(tier)
    expect_identical(
      unname(defaults[["dictionary_resource"]]), unname(ids[[tier]]),
      info = tier
    )
    # The category map was already at its first-publication version and must
    # not have been renumbered along with the dictionaries.
    expect_identical(
      unname(defaults[["category_map_resource"]]), "lisa_category_map@1.0.0",
      info = tier
    )
    # Upstream provenance is a different version space entirely.
    expect_identical(
      unname(defaults[["term2gene_resource"]]), "msigdb_term2gene@2026.1",
      info = tier
    )
  }
})

test_that("the new v1 identities serve the exact current scientific bytes", {
  cache <- e1_empty_cache()
  hashes <- e1_hashes()
  expected <- c(core = "core", expanded = "expanded")

  for (tier in names(expected)) {
    resolved <- lisaR:::lisa_resolve_dictionary_resource(
      paste0("lisa_dictionary_", tier), "1.0.0", "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
    expect_identical(resolved$sha256, unname(hashes[[tier]]), info = tier)
    expect_identical(
      lisaR:::lisa_sha256_file(resolved$path), unname(hashes[[tier]]),
      info = tier
    )
    expect_identical(resolved$schema, "lisa_dictionary@2", info = tier)
    expect_identical(resolved$modality, "all", info = tier)
    # A direct request is not an alias: the receipt must not claim one.
    expect_true(is.na(resolved$superseded_resource_id), info = tier)

    # The score column must not reappear anywhere in the active runtime.
    table <- utils::read.delim(
      resolved$path, sep = "\t", header = TRUE, quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
      na.strings = character()
    )
    expect_identical(
      names(table), lisaR:::lisa_dictionary_runtime_columns(), info = tier
    )
    expect_false("LISA_score" %in% names(table), info = tier)
  }

  map <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_category_map", "1.0.0", "Homo sapiens",
    modality = "transcriptomic/genomic", cache_root = cache
  )
  expect_identical(map$sha256, unname(hashes[["category_map"]]))
  expect_identical(map$schema, "category_map@1")
})

test_that("the new standard identities keep standard-tier semantics", {
  # A new logical ID that the tier map did not know would be treated as a
  # custom dictionary: wrong validator on install, wrong tier in a pipeline.
  for (tier in c("core", "expanded")) {
    expect_identical(
      lisaR:::lisa_dictionary_tier_for_resource(
        paste0("lisa_dictionary_", tier, "@1.0.0")
      ),
      tier, info = tier
    )
    # The retired namespace keeps its standard meaning too.
    expect_identical(
      lisaR:::lisa_dictionary_tier_for_resource(paste0("lisa_", tier, "@2.0.0")),
      tier, info = tier
    )
  }
  # The synthetic fixtures are emphatically not standard tiers.
  expect_null(lisaR:::lisa_dictionary_tier_for_resource(
    "lisa_dictionary_quickstart@1.0.0"
  ))
  expect_null(lisaR:::lisa_dictionary_tier_for_resource(
    "lisa_dictionary_custom_example@1.0.0"
  ))

  # A default tier selection resolves to the new identity and reports the tier.
  selection <- lisaR:::lisa_select_pipeline_resource_references(
    list(profile = "transcriptomic/genomic", lisa_dictionary = "expanded")
  )
  expect_identical(selection$dictionary_tier, "expanded")
  expect_identical(
    selection$references$dictionary_resource$resource_id,
    "lisa_dictionary_expanded@1.0.0"
  )
})

test_that("superseded current-content identities resolve to the same bytes", {
  cache <- e1_empty_cache()
  hashes <- e1_hashes()
  cases <- list(
    list(old = "lisa_core", version = "2.0.0",
         new = "lisa_dictionary_core@1.0.0", hash = hashes[["core"]],
         modality = "global proteomic"),
    list(old = "lisa_expanded", version = "2.0.0",
         new = "lisa_dictionary_expanded@1.0.0", hash = hashes[["expanded"]],
         modality = "transcriptomic/genomic"),
    list(old = "lisa_quickstart_dictionary", version = "4.0.0",
         new = "lisa_dictionary_quickstart@1.0.0", hash = hashes[["quickstart"]],
         modality = "transcriptomic/genomic")
  )
  for (case in cases) {
    resolved <- lisaR:::lisa_resolve_dictionary_resource(
      case$old, case$version, "Homo sapiens",
      modality = case$modality, cache_root = cache
    )
    expect_identical(resolved$resource_id, case$new, info = case$old)
    expect_identical(resolved$sha256, unname(case$hash), info = case$old)
    # Both halves of the receipt: what the run asked for and what it read.
    expect_identical(
      resolved$requested_resource_id,
      paste0(case$old, "@", case$version), info = case$old
    )
    expect_identical(resolved$requested_logical_id, case$old, info = case$old)
    expect_identical(
      resolved$superseded_resource_id,
      paste0(case$old, "@", case$version), info = case$old
    )
  }
})

test_that("version-less old scientific names resolve deterministically forward", {
  cache <- e1_empty_cache()
  hashes <- e1_hashes()
  for (tier in c("core", "expanded")) {
    resolved <- lisaR:::lisa_resolve_dictionary_resource(
      paste0("lisa_", tier), species = "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
    # Forward to the current contents, never back to the score-carrying bytes.
    expect_identical(
      resolved$resource_id, paste0("lisa_dictionary_", tier, "@1.0.0"),
      info = tier
    )
    expect_identical(resolved$sha256, unname(hashes[[tier]]), info = tier)
    expect_identical(resolved$requested_logical_id, paste0("lisa_", tier),
                     info = tier)
    expect_identical(resolved$superseded_resource_id, paste0("lisa_", tier),
                     info = tier)
  }
})

test_that("the alias is closed and cannot redirect a historical pin", {
  cache <- e1_empty_cache()
  historical <- c(
    core = "32d5ccafaec18fd6bdc51d3a23d80c21b7e996f5fe061a510c1b24ef3fafc876",
    expanded = "f8dd8bdbf61278d080c6349cebe9fd74993dc40bb07051164412e783358edaf4"
  )
  # The whole risk of a cross-namespace alias is that it swallows an old pin.
  for (tier in c("core", "expanded")) {
    for (version in c("0.1.0", "1.0.0")) {
      resolved <- lisaR:::lisa_resolve_dictionary_resource(
        paste0("lisa_", tier), version, "Homo sapiens",
        modality = "transcriptomic/genomic", cache_root = cache
      )
      info <- paste(tier, version)
      expect_identical(resolved$sha256, unname(historical[[tier]]), info = info)
      expect_identical(resolved$schema, "lisa_dictionary@1", info = info)
      # Still the historical logical namespace, at its own version.
      expect_identical(resolved$logical_id, paste0("lisa_", tier), info = info)
      expect_identical(resolved$version, "1.0.0", info = info)
    }
  }

  aliases <- lisaR:::lisa_current_release_dictionary_aliases()
  # The closed table must not contain a single historical version, in either
  # direction.
  expect_false(any(
    aliases$old_logical_id %in% c("lisa_core", "lisa_expanded") &
      aliases$old_version %in% c("0.1.0", "1.0.0")
  ))
  expect_true(all(aliases$new_version == "1.0.0"))
  expect_true(all(startsWith(aliases$new_logical_id, "lisa_dictionary_")))
  # Every declared digest is a real installed artifact under the new identity.
  builtin <- lisaR:::lisa_builtin_resource_registry()
  bundled <- aliases[aliases$new_logical_id != "lisa_dictionary_custom_example", ]
  for (i in seq_len(nrow(bundled))) {
    row <- builtin[builtin$logical_id == bundled$new_logical_id[[i]] &
                     builtin$version == bundled$new_version[[i]], , drop = FALSE]
    expect_identical(nrow(row), 1L, info = bundled$new_logical_id[[i]])
    expect_identical(row$sha256[[1L]], bundled$sha256[[i]],
                     info = bundled$new_logical_id[[i]])
    expect_identical(row$schema[[1L]], bundled$schema[[i]],
                     info = bundled$new_logical_id[[i]])
  }
})

test_that("an unknown version in either namespace is still refused", {
  cache <- e1_empty_cache()
  for (id in c("lisa_core", "lisa_dictionary_core")) {
    expect_error(
      lisaR:::lisa_resolve_dictionary_resource(
        id, "9.9.9", "Homo sapiens",
        modality = "transcriptomic/genomic", cache_root = cache
      ),
      "LISA-RESOURCE-006", info = id
    )
  }
  # A version-less reference to a name that was never published must not be
  # rescued by the alias table either.
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_not_a_dictionary", species = "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )
})

test_that("superseded current identities stay reserved against re-minting", {
  builtin <- lisaR:::lisa_builtin_resource_registry()
  columns <- lisaR:::lisa_dictionary_registry_columns()
  carrier <- builtin[builtin$logical_id == "lisa_dictionary_core", columns,
                     drop = FALSE]

  reserved <- lisaR:::lisa_reserved_historical_identities()
  expect_true(all(
    c("lisa_core@2.0.0", "lisa_expanded@2.0.0",
      "lisa_quickstart_dictionary@4.0.0") %in% reserved$resource_id
  ))
  # These left the bundle manifest, so the built-in collision guard no longer
  # covers them; without the reservation an external registry could re-mint
  # `lisa_core@2.0.0` over different bytes.
  for (id in c("lisa_core", "lisa_expanded")) {
    row <- carrier
    row$logical_id <- id
    row$version <- "2.0.0"
    expect_error(
      lisaR:::lisa_active_resource_registry(row),
      "LISA-RESOURCE-030.*cannot be re-minted", info = id
    )
  }

  # The synthetic custom-example fixture was never a built-in registry row, so
  # a user may still register it. Reserving it would break the documented
  # custom-resource workflow.
  expect_false(
    "lisa_example_custom_dictionary@2.0.0" %in% reserved$resource_id
  )
})

test_that("corrupt or mismatched alias content fails closed", {
  cache <- e1_empty_cache()
  aliases <- lisaR:::lisa_current_release_dictionary_aliases()
  # A mapping row whose declared digest no longer matches the identity it
  # points at must refuse to resolve rather than quietly hand over other bytes.
  broken <- aliases
  broken$sha256 <- paste0(strrep("0", 63), "1")
  with_mocked_bindings(
    lisa_current_release_dictionary_aliases = function() broken,
    expect_error(
      lisaR:::lisa_resolve_dictionary_resource(
        "lisa_core", "2.0.0", "Homo sapiens",
        modality = "transcriptomic/genomic", cache_root = cache
      ),
      "LISA-RESOURCE-030.*must not be reused for different content"
    ),
    .package = "lisaR"
  )

  # Likewise for a schema the alias does not promise.
  wrong_schema <- aliases
  wrong_schema$schema <- "lisa_dictionary@1"
  with_mocked_bindings(
    lisa_current_release_dictionary_aliases = function() wrong_schema,
    expect_error(
      lisaR:::lisa_resolve_dictionary_resource(
        "lisa_core", "2.0.0", "Homo sapiens",
        modality = "transcriptomic/genomic", cache_root = cache
      ),
      "LISA-RESOURCE-030"
    ),
    .package = "lisaR"
  )
})

test_that("profile, species and custom-override precedence is unchanged", {
  cache <- e1_empty_cache()
  # The scientific dictionaries are profile-independent: every supported
  # profile resolves the same all-modality row, through the new identity and
  # through the alias alike.
  for (profile in c("transcriptomic/genomic", "global proteomic", "targeted")) {
    direct <- lisaR:::lisa_resolve_dictionary_resource(
      "lisa_dictionary_core", "1.0.0", "Homo sapiens",
      modality = profile, cache_root = cache
    )
    aliased <- lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "2.0.0", "Homo sapiens",
      modality = profile, cache_root = cache
    )
    expect_identical(direct$path, aliased$path, info = profile)
    expect_identical(direct$modality, "all", info = profile)
  }

  # An unregistered species is refused, not silently answered by the alias.
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "2.0.0", "Mus musculus",
      modality = "transcriptomic/genomic", cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )

  # A profile-scoped external row for the NEW identity still wins for its own
  # profile, and the all-modality built-in still answers every other profile.
  builtin <- lisaR:::lisa_builtin_resource_registry()
  columns <- lisaR:::lisa_dictionary_registry_columns()
  override <- builtin[builtin$logical_id == "lisa_dictionary_core", columns,
                      drop = FALSE]
  override$modality <- "global proteomic"
  registry <- lisaR:::lisa_active_resource_registry(override)
  rows <- registry[registry$logical_id == "lisa_dictionary_core", , drop = FALSE]
  expect_identical(nrow(rows), 2L)
  expect_setequal(rows$modality, c("all", "global proteomic"))
})

test_that("the migration ledger documents the republication without steering it", {
  migration <- lisaR:::lisa_read_resource_migration_manifest()
  republish <- migration[
    migration$change_kind == "republish_first_version", , drop = FALSE
  ]
  expect_gt(nrow(republish), 0L)
  # A republication moves identity and filename, never bytes.
  expect_identical(republish$old_sha256, republish$new_sha256)
  expect_true(all(republish$content_transformed == "FALSE"))
  expect_true(all(nzchar(republish$rationale)))

  # Ledger and code must agree, but the code is the authority: the resolver
  # reads lisa_current_release_dictionary_aliases(), not this file.
  aliases <- lisaR:::lisa_current_release_dictionary_aliases()
  versioned <- aliases[!is.na(aliases$old_version), , drop = FALSE]
  ledger <- republish[
    republish$old_resource_id != republish$new_resource_id, , drop = FALSE
  ]
  expect_setequal(
    paste0(versioned$old_logical_id, "@", versioned$old_version),
    ledger$old_resource_id
  )
  expect_setequal(
    paste0(versioned$new_logical_id, "@", versioned$new_version),
    ledger$new_resource_id
  )

  # The C5 and C6 history is retained verbatim alongside it.
  expect_true(all(
    c("relabel", "hold", "supersede_transform") %in% migration$change_kind
  ))
  transform <- migration[
    migration$change_kind == "supersede_transform", , drop = FALSE
  ]
  expect_true(all(transform$old_sha256 != transform$new_sha256))
})

test_that("the retired artifacts are still present and byte-identical", {
  # A rename must not have become a move-over-the-top of the historical files.
  retained <- c(
    "extdata/dictionaries/lisa_core_runtime_v0_1.tsv" =
      "32d5ccafaec18fd6bdc51d3a23d80c21b7e996f5fe061a510c1b24ef3fafc876",
    "extdata/dictionaries/lisa_expanded_runtime_v0_1.tsv" =
      "f8dd8bdbf61278d080c6349cebe9fd74993dc40bb07051164412e783358edaf4",
    "extdata/quick-start/example_dictionary.tsv" =
      "4c9cb6f9513b5ee24b7b7bf50081db86558ee5993346c243fb69fc747de8901c",
    "extdata/custom-resources/example_custom_dictionary.tsv" =
      "cf96afe68a10be5f93372e7ce296b7410c4edbbdb7cc783739c9138000457d53"
  )
  for (artifact in names(retained)) {
    path <- system.file(artifact, package = "lisaR")
    expect_true(nzchar(path) && file.exists(path), info = artifact)
    expect_identical(
      lisaR:::lisa_sha256_file(path), unname(retained[[artifact]]),
      info = artifact
    )
  }
})


test_that("sample dictionary label update preserves all scientific fields and old bytes", {
  old <- system.file("extdata/quick-start/lisa_dictionary_quickstart_v1_0.tsv",package="lisaR")
  new <- system.file("extdata/quick-start/lisa_dictionary_quickstart_v1_1.tsv",package="lisaR")
  expect_identical(lisaR:::lisa_sha256_file(old),unname(e1_hashes()[["quickstart"]]))
  a <- read.delim(old,check.names=FALSE); b <- read.delim(new,check.names=FALSE)
  cols <- setdiff(names(a),"source_id")
  expect_identical(a[cols],b[cols])
  expect_true(all(b$source_id=="SYNTHETIC_SAMPLE"))
})
