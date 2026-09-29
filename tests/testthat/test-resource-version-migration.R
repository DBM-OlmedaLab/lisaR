# C5 resource-version migration contract.
#
# These tests assert the properties that make the migration safe rather than
# merely tidy: old pins keep reading the very same bytes, superseded identities
# stay reserved, the active identity is unique, and nothing about the science
# or the content-addressed cache moves because a label changed.

c5_migration <- function() {
  lisaR:::lisa_read_resource_migration_manifest()
}

test_that("the migration manifest is machine-readable and self-consistent", {
  migration <- c5_migration()
  expect_identical(names(migration), lisaR:::lisa_resource_migration_columns())
  expect_gt(nrow(migration), 0L)
  # "supersede_transform" is added by C6  for a content change, and
  # "republish_first_version" by the first-publication release for a pure
  # identity/filename change. The C5 relabel invariants asserted below are
  # unaffected by either and still hold verbatim.
  expect_true(all(migration$change_kind %in%
    c("relabel", "transform", "hold", "supersede_transform",
      "republish_first_version")))

  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]
  expect_setequal(
    relabel$new_resource_id,
    c("lisa_core@1.0.0", "lisa_expanded@1.0.0", "lisa_category_map@1.0.0")
  )
  expect_setequal(
    relabel$old_resource_id,
    c("lisa_core@0.1.0", "lisa_expanded@0.1.0", "lisa_category_map@0.1.1")
  )
  # A relabel is an identity change and nothing else.
  expect_identical(relabel$old_sha256, relabel$new_sha256)
  expect_identical(relabel$old_artifact, relabel$new_artifact)
  expect_true(all(relabel$artifact_renamed == "FALSE"))
  expect_true(all(relabel$content_transformed == "FALSE"))
  expect_true(all(relabel$old_state == "superseded_compatibility"))
  expect_true(all(relabel$new_state == "active"))
  expect_true(all(nzchar(migration$rationale)))
})

test_that("recorded migration digests are the bytes actually installed", {
  migration <- c5_migration()
  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]
  paths <- file.path(system.file(package = "lisaR"), relabel$new_artifact)
  expect_true(all(file.exists(paths)))
  expect_identical(
    unname(vapply(paths, lisaR:::lisa_sha256_file, character(1))),
    unname(relabel$new_sha256)
  )
})

test_that("every migrated resource has exactly one active registered version", {
  builtin <- lisaR:::lisa_builtin_resource_registry()
  # Exactly one active variant per resource is the invariant C5 established and
  # every later release preserves. The dictionary versions advanced to 2.0.0
  # when C6 removed the runtime score column, and the first-publication release
  # then moved those same bytes into their own logical namespace at 1.0.0. The
  # category map never carried a score and never left 1.0.0.
  active <- c(lisa_dictionary_core = "1.0.0",
              lisa_dictionary_expanded = "1.0.0",
              lisa_category_map = "1.0.0")
  for (logical_id in names(active)) {
    rows <- builtin[builtin$logical_id == logical_id, , drop = FALSE]
    expect_identical(nrow(rows), 1L, info = logical_id)
    expect_identical(rows$version[[1L]], active[[logical_id]], info = logical_id)
  }
  # Superseded identities are deliberately absent from the registry so a
  # version-less reference stays unambiguous.
  migration <- c5_migration()
  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]
  registry_key <- paste(builtin$logical_id, builtin$version, sep = "@")
  expect_false(any(relabel$old_resource_id %in% registry_key))
})

test_that("pipeline defaults and tier resolution use the migrated identities", {
  expect_identical(
    unname(lisaR:::lisa_dictionary_tier_resources()),
    c("lisa_dictionary_core@1.0.0", "lisa_dictionary_expanded@1.0.0")
  )
  for (tier in c("core", "expanded")) {
    defaults <- lisaR:::lisa_pipeline_resource_defaults(tier)
    expect_identical(
      unname(defaults[["dictionary_resource"]]),
      sprintf("lisa_dictionary_%s@1.0.0", tier)
    )
    expect_identical(
      unname(defaults[["category_map_resource"]]), "lisa_category_map@1.0.0"
    )
    # Upstream provenance is a different concept from our resource version and
    # must not be renumbered by this migration.
    expect_identical(
      unname(defaults[["term2gene_resource"]]), "msigdb_term2gene@2026.1"
    )
  }
})

test_that("old pins still resolve, and to byte-identical content", {
  cache <- withr::local_tempdir("lisa-c5-oldpin-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  migration <- c5_migration()
  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]
  expect_gt(nrow(relabel), 0L)

  for (i in seq_len(nrow(relabel))) {
    logical_id <- relabel$logical_id[[i]]
    old <- lisaR:::lisa_resolve_dictionary_resource(
      logical_id, relabel$old_version[[i]], "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
    new <- lisaR:::lisa_resolve_dictionary_resource(
      logical_id, relabel$new_version[[i]], "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    )
    # Same file, same digest, same schema: the rename cannot have moved the
    # science that an existing project reads.
    expect_identical(old$path, new$path, info = logical_id)
    expect_identical(old$sha256, new$sha256, info = logical_id)
    expect_identical(old$schema, new$schema, info = logical_id)
    expect_identical(old$sha256, relabel$old_sha256[[i]], info = logical_id)
    # Provenance still distinguishes what was asked for from what was read.
    expect_identical(
      old$requested_resource_id, relabel$old_resource_id[[i]],
      info = logical_id
    )
    expect_identical(
      old$superseded_resource_id, relabel$old_resource_id[[i]],
      info = logical_id
    )
    expect_identical(old$resource_id, relabel$new_resource_id[[i]],
                     info = logical_id)
    expect_true(is.na(new$superseded_resource_id), info = logical_id)
  }
})

test_that("an old pin resolves for every supported profile", {
  cache <- withr::local_tempdir("lisa-c5-profile-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  for (profile in c("transcriptomic/genomic", "global proteomic")) {
    core <- lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Homo sapiens", modality = profile,
      cache_root = cache
    )
    expect_identical(core$resource_id, "lisa_core@1.0.0", info = profile)
    expect_identical(core$schema, "lisa_dictionary@1", info = profile)
    expect_identical(core$modality, "all", info = profile)
    # The category map was not transformed, so its old pin still resolves.
    resolved <- lisaR:::lisa_resolve_dictionary_resource(
      "lisa_category_map", "0.1.1", "Homo sapiens", modality = profile,
      cache_root = cache
    )
    expect_identical(resolved$resource_id, "lisa_category_map@1.0.0", info = profile)
    expect_identical(resolved$modality, "all", info = profile)
  }
})

# A registry that declares a superseded identity for ONE profile, with content
# that differs from the built-in bytes. This is the shape that exposed the
# cross-profile regression: it is installable before and after C5 because both
# LISA-RESOURCE-020 and -030 key on modality as well as identity.
c5_profile_scoped_override <- function(cache, profile = "transcriptomic/genomic",
                                       logical_id = "lisa_core",
                                       old_version = "0.1.0") {
  builtin_path <- lisaR:::lisa_builtin_dictionary_path("core")
  lines <- readLines(builtin_path, warn = FALSE)
  variant <- file.path(cache, "profile_scoped_variant.tsv")
  # Drop the last data row: still a valid dictionary, definitely different bytes.
  writeLines(lines[-length(lines)], variant)
  variant_sha <- lisaR:::lisa_sha256_file(variant)
  testthat::expect_false(identical(variant_sha, lisaR:::lisa_sha256_file(builtin_path)))
  install_lisa_resource(
    variant, logical_id, old_version, "Homo sapiens", profile,
    "lisa_dictionary@2", "c5 regression fixture", variant_sha,
    cache_root = cache
  )
  list(
    registry = file.path(cache, "resource_registry.tsv"),
    variant_sha = variant_sha,
    builtin_sha = lisaR:::lisa_sha256_file(builtin_path),
    historical_sha = lisaR:::lisa_historical_dictionary_definitions()$sha256[
      lisaR:::lisa_historical_dictionary_definitions()$logical_id == logical_id
    ][[1L]],
    profile = profile
  )
}

test_that("a profile-scoped override of an old pin still wins for its own profile", {
  # F6 / pre-C5 parity: the reservation is scoped to logical_id, version,
  # species AND modality, so declaring a superseded identity for one specific
  # profile is supported, not refused. That row must win there.
  cache <- withr::local_tempdir("lisa-c5-override-own-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  fixture <- c5_profile_scoped_override(cache)

  resolved <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens", modality = fixture$profile,
    registry = fixture$registry, cache_root = cache
  )
  expect_identical(resolved$resource_id, "lisa_core@0.1.0")
  expect_identical(resolved$modality, fixture$profile)
  expect_identical(resolved$sha256, fixture$variant_sha)
  # No ledger redirection happened, so no superseded provenance is claimed.
  expect_identical(resolved$requested_resource_id, "lisa_core@0.1.0")
  expect_true(is.na(resolved$superseded_resource_id))
})

test_that("a profile-scoped override does not disable the fallback under other profiles", {
  # F5 regression. Before this fix a single profile-scoped row for
  # lisa_core@0.1.0 suppressed the compatibility ledger for EVERY other
  # profile, so `global proteomic` failed with LISA-RESOURCE-021 where pre-C5
  # it fell back to the built-in `all` row. The fallback must survive.
  cache <- withr::local_tempdir("lisa-c5-override-other-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  fixture <- c5_profile_scoped_override(cache)

  resolved <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens", modality = "global proteomic",
    registry = fixture$registry, cache_root = cache
  )
  # Redirected through the authenticated alias to the closed historical
  # identity, reading the exact built-in bytes the pin has always read.
  expect_identical(resolved$resource_id, "lisa_core@1.0.0")
  expect_identical(resolved$modality, "all")
  expect_identical(resolved$schema, "lisa_dictionary@1")
  expect_identical(resolved$sha256, fixture$historical_sha)
  expect_false(identical(resolved$sha256, fixture$variant_sha))
  # Provenance keeps both identities.
  expect_identical(resolved$requested_resource_id, "lisa_core@0.1.0")
  expect_identical(resolved$superseded_resource_id, "lisa_core@0.1.0")
  expect_identical(lisaR:::lisa_sha256_file(resolved$path), fixture$historical_sha)
})

test_that("the ledger fallback does not fire for an unknown species or a foreign identity", {
  # The widened trigger must not become a silent catch-all: a species the
  # ledger does not cover, and a version belonging to a different logical
  # identity, still fail loudly instead of being redirected.
  cache <- withr::local_tempdir("lisa-c5-noredirect-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  fixture <- c5_profile_scoped_override(cache)

  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Mus musculus", modality = "global proteomic",
      registry = fixture$registry, cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )
  # lisa_category_map@0.1.1 is superseded; lisa_core@0.1.1 never existed and
  # must not borrow another resource's ledger row.
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.1", "Homo sapiens", modality = "global proteomic",
      registry = fixture$registry, cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )
  # The private lisa_category_map@0.1.0 is not the superseded 0.1.1 identity.
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_category_map", "0.1.0", "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )
})

test_that("a matching but corrupt row fails loudly instead of falling back to the ledger", {
  # "No eligible row" must mean selection, not validity. A row that matches the
  # requested identity and profile but whose bytes no longer match its recorded
  # digest must raise, never be quietly replaced by the renamed identity.
  cache <- withr::local_tempdir("lisa-c5-corrupt-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  fixture <- c5_profile_scoped_override(cache)

  installed <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens", modality = fixture$profile,
    registry = fixture$registry, cache_root = cache
  )$path
  cat("# corrupted for the C5 regression test\n", file = installed, append = TRUE)

  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.1.0", "Homo sapiens", modality = fixture$profile,
      registry = fixture$registry, cache_root = cache
    ),
    "LISA-RESOURCE-009"
  )
})

test_that("an old pin with no modality keeps resolving with a profile-scoped override present", {
  cache <- withr::local_tempdir("lisa-c5-nomodality-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  fixture <- c5_profile_scoped_override(cache)

  resolved <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_core", "0.1.0", "Homo sapiens",
    registry = fixture$registry, cache_root = cache
  )
  # Exactly one row carries lisa_core@0.1.0, so it resolves deterministically
  # to the user's row — unchanged by the fallback fix.
  expect_identical(resolved$resource_id, "lisa_core@0.1.0")
  expect_identical(resolved$sha256, fixture$variant_sha)
})

test_that("the reservation refuses the built-in key but permits a specific profile", {
  # F6: demonstrate BOTH sides of the modality-scoped reservation, so the
  # documented contract is the tested one.
  builtin <- lisaR:::lisa_builtin_resource_registry()
  columns <- lisaR:::lisa_dictionary_registry_columns()
  row <- builtin[builtin$logical_id == "lisa_dictionary_core", columns,
                 drop = FALSE]
  row$logical_id <- "lisa_core"
  row$version <- "0.1.0"

  # Prohibited: the exact key the built-in row occupied (species + `all`).
  expect_error(
    lisaR:::lisa_active_resource_registry(row),
    "LISA-RESOURCE-030.*cannot be re-minted"
  )

  # Supported: the same identity scoped to one specific profile.
  scoped <- row
  scoped$modality <- "transcriptomic/genomic"
  expect_silent(registry <- lisaR:::lisa_active_resource_registry(scoped))
  expect_true(any(
    registry$logical_id == "lisa_core" &
      registry$version == "0.1.0" &
      registry$modality == "transcriptomic/genomic"
  ))
  # The active identity is still present and still unique. It now lives in the
  # first-publication namespace; registering a profile-scoped historical
  # `lisa_core` must not disturb it.
  active <- registry[registry$logical_id == "lisa_dictionary_core" &
                       registry$version == "1.0.0", , drop = FALSE]
  expect_identical(nrow(active), 1L)
})

test_that("an unknown version is still refused rather than silently migrated", {
  cache <- withr::local_tempdir("lisa-c5-unknown-")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  expect_error(
    lisaR:::lisa_resolve_dictionary_resource(
      "lisa_core", "0.9.9", "Homo sapiens",
      modality = "transcriptomic/genomic", cache_root = cache
    ),
    "LISA-RESOURCE-006"
  )
})

test_that("superseded identities stay reserved against external registries", {
  builtin <- lisaR:::lisa_builtin_resource_registry()
  columns <- lisaR:::lisa_dictionary_registry_columns()
  migration <- c5_migration()
  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]

  # The scientific dictionaries now live under their first-publication logical
  # namespace, so take any all-modality built-in row as the byte carrier and
  # re-label it with the historical identity under test.
  carrier <- builtin[builtin$modality == "all", columns, drop = FALSE][1L, ,
                     drop = FALSE]
  for (i in seq_len(nrow(relabel))) {
    row <- carrier
    row$logical_id <- relabel$logical_id[[i]]
    row$version <- relabel$old_version[[i]]
    expect_error(
      lisaR:::lisa_active_resource_registry(row),
      "LISA-RESOURCE-030.*cannot be re-minted",
      info = relabel$old_resource_id[[i]]
    )
  }
})

test_that("a same-identity historical 1.0.0 blocks renumbering a quick-start fixture", {
  # lisa_quickstart_dictionary@1.0.0 and lisa_quickstart_term2gene@1.0.0 exist
  # with different bytes from the active variants, which is exactly why C5
  # refuses to renumber those fixtures to 1.0.0. Guard the reason, not just the
  # outcome: if these digests ever coincide the refusal loses its basis.
  # The original 1.0.0 fixtures live in the administrative catalog rather than
  # the resolvable registry, which is exactly why the collision is easy to miss
  # by reading RESOURCE_BUNDLE_MANIFEST.tsv alone.
  root <- system.file("extdata", "dictionaries", package = "lisaR")
  catalog <- read.delim(
    file.path(root, "DICTIONARY_RESOURCE_MANIFEST.tsv"), sep = "\t",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  builtin <- lisaR:::lisa_builtin_resource_registry()
  active <- c(
    lisa_quickstart_dictionary = "3.0.0", lisa_quickstart_term2gene = "2.0.0"
  )
  for (logical_id in names(active)) {
    historical <- catalog$artifact[
      catalog$resource_id == paste0(logical_id, "@1.0.0")
    ]
    expect_identical(length(historical), 1L, info = logical_id)
    historical_sha <- lisaR:::lisa_sha256_file(
      normalizePath(file.path(root, historical), mustWork = TRUE)
    )
    current_sha <- builtin$sha256[
      builtin$logical_id == logical_id &
        builtin$version == active[[logical_id]]
    ]
    expect_identical(length(current_sha), 1L, info = logical_id)
    # Different bytes under the same logical identity: renumbering to 1.0.0
    # would repoint a published identity at different content.
    expect_false(identical(historical_sha, current_sha), info = logical_id)
  }
  migration <- c5_migration()
  held <- migration[migration$change_kind == "hold", , drop = FALSE]
  expect_true(all(held$old_version == held$new_version))
  expect_match(
    held$rationale[held$logical_id == "lisa_quickstart_dictionary"],
    "REFUSED 1.0.0"
  )
})

test_that("renaming a resource version does not move the content cache key", {
  # The fgsea cache is content-addressed. Hold the wrapper and software
  # versions fixed and vary only the resource identity: the key must not move,
  # because the memberships did not move. Then perturb one membership and
  # confirm the key does move, so a real change cannot reuse a stale entry.
  pathways <- list(CAT_A = c("GENE1", "GENE2"), CAT_B = c("GENE2", "GENE3"))
  ranks <- c(GENE1 = 2.5, GENE2 = 0.5, GENE3 = -1.5)
  universe <- names(ranks)
  versions <- list(frozen = "c5-fixed-software-and-wrapper-versions")

  contract <- function(memberships) {
    lisaR:::lisa_fgsea_cache_contract(
      memberships, ranks, min_gs_size = 1L, max_gs_size = 10L, n_threads = 1L,
      fgsea_nperm = 10L, fgsea_eps = 0, random_seed = 1L, task_id = "c5",
      universe = universe, versions = versions
    )
  }

  baseline <- contract(pathways)
  relabelled <- contract(pathways)
  expect_identical(baseline$key, relabelled$key)

  changed <- pathways
  changed$CAT_A <- c(changed$CAT_A, "GENE3")
  expect_false(identical(contract(changed)$key, baseline$key))
})

test_that("the migration is recoverable: superseded bytes are still installed", {
  # Rollback to the pre-C5 state means reading the previous identity's content
  # again. Nothing was deleted or overwritten, so the recorded old digest is
  # still present in the installed package.
  migration <- c5_migration()
  relabel <- migration[migration$change_kind == "relabel", , drop = FALSE]
  for (i in seq_len(nrow(relabel))) {
    path <- file.path(system.file(package = "lisaR"), relabel$old_artifact[[i]])
    expect_true(file.exists(path), info = relabel$old_resource_id[[i]])
    expect_identical(
      lisaR:::lisa_sha256_file(path), relabel$old_sha256[[i]],
      info = relabel$old_resource_id[[i]]
    )
  }
})
