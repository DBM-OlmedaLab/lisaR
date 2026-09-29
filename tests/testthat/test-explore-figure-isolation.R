
g2_repair_workspace <- function(label = "ws", manifest = strrep("a", 64L)) {
  root <- file.path(tempdir(), paste0("lisa-explore-g2-repair-", label))
  unlink(root, recursive = TRUE, force = TRUE)
  dir.create(file.path(root, "prepared"), recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE),
               envir = parent.frame())
  structure(list(root = root, source_run = file.path(root, "source"),
                 source_manifest_hash = manifest,
                 index_path = file.path(root, "index.json"),
                 session_path = file.path(root, "session.json")),
            class = "lisa_explore_workspace")
}

g2_repair_declaration <- function(ws, snapshot = "snapshot-a") {
  declaration <- list(
    schema = "lisa-explore-kegg/2", cache_root = "/synthetic/cache",
    snapshot_id = snapshot, species = "Homo sapiens", organism = "hsa",
    max_abs_log2fc = "0.5", color_power = "1", access_mode = "cache_only")
  lisaR:::lisa_explore_write_json(
    declaration, lisaR:::lisa_explore_kegg_declaration_path(ws),
    run_root = ws$root)
  declaration
}

g2_repair_map_rows <- function(analysis, collection, category = "CAT_A",
                               kegg_id = "hsa04110", available = TRUE) {
  data.frame(
    analysis_id = analysis, collection = collection, category_id = category,
    kegg_id = kegg_id, kegg_title = "Synthetic pathway", rank = 1L,
    n_evidence_genes = 3L, resources_available = available,
    kegg_resource_digest = if (available) strrep("a", 64L) else "",
    stringsAsFactors = FALSE)
}

g2_repair_write_preparation <- function(ws, declaration, units, rows,
                                        record = TRUE) {
  context <- lisaR:::lisa_explore_kegg_index_context(ws, declaration)
  if (!nrow(rows)) rows <- lisaR:::lisa_explore_kegg_empty_index(context)
  if (nrow(rows)) {
    for (field in names(context)) rows[[field]] <- context[[field]]
  }
  target <- lisaR:::lisa_explore_kegg_index_path(ws)
  withr::local_options(list(lisaR.run_root = ws$root))
  lisaR:::write_lisa_tsv(rows, target)
  if (isTRUE(record)) {
    lisaR:::lisa_explore_write_json(
      lisaR:::lisa_explore_kegg_preparation_record(context, units, target),
      lisaR:::lisa_explore_kegg_preparation_path(ws), run_root = ws$root)
  }
  invisible(rows)
}

g2_repair_request <- function(analysis = "A", collection = "C1",
                              category = "CAT_A", entity = NA_character_) {
  list(unit_type = "single_de", analysis_id = analysis,
       contrast_id = NA_character_, collection = collection,
       category_id = category, product = "kegg_pathway_map",
       entity = entity, variant = NA_character_)
}

test_that("G2 reasons distinguish unprepared, prepared-empty, and unavailable resources", {
  ws <- g2_repair_workspace("reasons")
  declaration <- g2_repair_declaration(ws)
  request <- g2_repair_request()

  reason <- lisaR:::lisa_explore_unavailable_reason(ws, request)
  expect_match(reason, "has not been prepared for this workspace")
  expect_false(grepl("genes map to no pathway", reason, fixed = TRUE))

  units <- data.frame(analysis_id = "OTHER", collection = "C2",
                      stringsAsFactors = FALSE)
  g2_repair_write_preparation(
    ws, declaration, units,
    g2_repair_map_rows("OTHER", "C2", "CAT_OTHER"))
  reason <- lisaR:::lisa_explore_unavailable_reason(ws, request)
  expect_match(reason, "has not been prepared for A / C1", fixed = TRUE)
  expect_false(grepl("genes map to no pathway", reason, fixed = TRUE))

  units <- data.frame(analysis_id = "A", collection = "C1",
                      stringsAsFactors = FALSE)
  g2_repair_write_preparation(
    ws, declaration, units, lisaR:::lisa_explore_kegg_empty_index())
  reason <- lisaR:::lisa_explore_unavailable_reason(ws, request)
  expect_match(reason, "was prepared for A / C1", fixed = TRUE)
  expect_match(reason, "found no pathway associations", fixed = TRUE)

  g2_repair_write_preparation(
    ws, declaration, units,
    g2_repair_map_rows("A", "C1", available = FALSE))
  reason <- lisaR:::lisa_explore_unavailable_reason(ws, request)
  expect_match(reason, "known pathway association", fixed = TRUE)
  expect_match(reason, "no validated KGML diagram and base image", fixed = TRUE)
  expect_match(reason, "does not download", fixed = TRUE)
})

test_that("G2 scoped preparation is additive, replaces empty units, and drops stale context", {
  ws <- g2_repair_workspace("additive")
  g2_repair_declaration(ws)
  catalog <- data.frame(
    unit_type = "single_de", product = "volcano",
    analysis_id = c("A", "A", "B"), collection = c("C1", "C2", "C3"),
    gsea_path = c("a", "b", "c"), de_path = c("d", "e", "f"),
    summary_path = c("g", "h", "i"), stringsAsFactors = FALSE)
  synthetic <- new.env(parent = emptyenv())
  synthetic$rows <- list(
    `A\rC1` = g2_repair_map_rows("A", "C1", "CAT_1", "hsa00010"),
    `A\rC2` = g2_repair_map_rows("A", "C2", "CAT_2", "hsa00020"),
    `B\rC3` = g2_repair_map_rows("B", "C3", "CAT_3", "hsa00030"))
  testthat::local_mocked_bindings(
    lisa_extension_discover_catalog = function(...) catalog,
    lisa_resolve_package_dir = function() "/synthetic/package",
    lisa_explore_prepare_one_kegg_index = function(ws, unit, ...) {
      synthetic$rows[[paste(unit$analysis_id[[1L]], unit$collection[[1L]],
                            sep = "\r")]]
    },
    .package = "lisaR")

  first <- lisa_explore_prepare_kegg_index(
    ws, analyses = "A", collections = "C1")
  expect_identical(first$collection, "C1")

  second <- lisa_explore_prepare_kegg_index(
    ws, analyses = "A", collections = "C2")
  expect_setequal(second$collection, c("C1", "C2"))
  status <- lisaR:::lisa_explore_kegg_preparation_status(ws)
  expect_identical(status$state, "current")
  expect_setequal(status$units$collection, c("C1", "C2"))

  synthetic$rows[["A\rC1"]] <- NULL
  replaced_empty <- lisa_explore_prepare_kegg_index(
    ws, analyses = "A", collections = "C1")
  expect_identical(replaced_empty$collection, "C2")
  status <- lisaR:::lisa_explore_kegg_preparation_status(ws)
  expect_setequal(status$units$collection, c("C1", "C2"))
  expect_match(
    lisaR:::lisa_explore_unavailable_reason(ws, g2_repair_request()),
    "was prepared for A / C1", fixed = TRUE)

  stale_ws <- ws
  stale_ws$source_manifest_hash <- strrep("b", 64L)
  stale_result <- lisa_explore_prepare_kegg_index(
    stale_ws, analyses = "B", collections = "C3")
  expect_identical(stale_result$collection, "C3")
  stale_status <- lisaR:::lisa_explore_kegg_preparation_status(stale_ws)
  expect_identical(stale_status$state, "current")
  expect_identical(stale_status$units$collection, "C3")
})

test_that("G2 legacy indexes are explicit and are not merged without provenance", {
  ws <- g2_repair_workspace("legacy")
  declaration <- g2_repair_declaration(ws)
  g2_repair_write_preparation(
    ws, declaration,
    data.frame(analysis_id = "OLD", collection = "OLD_COLLECTION"),
    g2_repair_map_rows("OLD", "OLD_COLLECTION", "OLD_CAT"),
    record = FALSE)
  reason <- lisaR:::lisa_explore_unavailable_reason(ws, g2_repair_request())
  expect_match(reason, "legacy native KEGG pathway index", fixed = TRUE)
  expect_match(reason, "cannot be determined", fixed = TRUE)

  catalog <- data.frame(
    unit_type = "single_de", product = "volcano", analysis_id = "A",
    collection = "C1", gsea_path = "a", de_path = "b", summary_path = "c",
    stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(
    lisa_extension_discover_catalog = function(...) catalog,
    lisa_resolve_package_dir = function() "/synthetic/package",
    lisa_explore_prepare_one_kegg_index = function(...)
      g2_repair_map_rows("A", "C1"),
    .package = "lisaR")
  expect_message(
    result <- lisa_explore_prepare_kegg_index(ws, "A", "C1"),
    "LISA-EXPLORE-069")
  expect_identical(result$analysis_id, "A")
  expect_false("OLD" %in% result$analysis_id)
})

test_that("G2 live and export association sets come from the same exact key", {
  status <- data.frame(
    key = c(rep("fig-shared", 3L), "fig-other"),
    category_id = c("CAT_A", "CAT_B", "CAT_C", "CAT_D"),
    stringsAsFactors = FALSE)
  recorded_submits <- c("CAT_A", "CAT_B")

  # The live shell has one status row per currently valid catalogue association.
  live <- sort(unique(status$category_id[status$key == "fig-shared"]))
  # The exporter calls this same helper for the available artifact instead of
  # limiting itself to the two categories that happened to submit it.
  exported <- lisaR:::lisa_explore_catalogue_associations(
    status, "fig-shared", fallback = recorded_submits)
  expect_identical(exported, live)
  expect_identical(exported, c("CAT_A", "CAT_B", "CAT_C"))
  expect_false(identical(exported, recorded_submits))
  expect_identical(
    lisaR:::lisa_explore_catalogue_associations(status, "fig-other"),
    "CAT_D")
})
