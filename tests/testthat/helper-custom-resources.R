g2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# Builds a dictionary in either schema . `legacy_score = FALSE`
# gives the canonical `lisa_dictionary@2` contract with no score column;
# `legacy_score = TRUE` reinstates the historical `lisa_dictionary@1` column so
# the transitional reader keeps being exercised by the suite.
g2_dictionary <- function(tier = "custom", score = NA_real_,
                          legacy_score = FALSE) {
  out <- data.frame(
    universe = c("GOBP-C2", "GOBP-C2", "PATHWAYS"),
    gene_set_id = c("GS_A", "GS_A", "GS_B"),
    gene_set_name = c("Gene set A", "Gene set A", "Gene set B"),
    source_id = c("TEST", "TEST", "TEST"),
    category_id = c("CAT_A", "CAT_B", "CAT_A"),
    category_display_name = c("Category A", "Category B", "Category A"),
    tier = rep(tier, 3L),
    stringsAsFactors = FALSE
  )
  if (isTRUE(legacy_score)) {
    out$LISA_score <- rep(score, 3L)
    out <- out[, c(setdiff(names(out), c("LISA_score", "tier")),
                   "LISA_score", "tier"), drop = FALSE]
  }
  out
}

g2_expanded_dictionary <- function(legacy_score = FALSE) {
  out <- g2_dictionary(tier = "expanded", score = 3, legacy_score = legacy_score)
  out$tier[[1L]] <- "core"
  if (isTRUE(legacy_score)) out$LISA_score[[1L]] <- 4
  out
}

g2_category_map <- function() {
  data.frame(
    category_id = c("CAT_A", "CAT_B"),
    display_name = c("Category A", "Category B"),
    macrogroup_id = c("SUPER_A", "SUPER_B"),
    macrogroup_name = c("Supercategory A", "Supercategory B"),
    macrogroup_order = c(1L, 2L),
    category_order_within_macrogroup = c(1L, 1L),
    notes = c("Local G2 fixture", "Local G2 fixture"),
    stringsAsFactors = FALSE
  )
}

g2_term2gene <- function() {
  data.frame(
    gs_collection = c("C5", "C5", "C2", "C2"),
    gs_subcollection = c("GO:BP", "GO:BP", "CP", "CP"),
    gs_name = c("GS_A", "GS_A", "GS_B", "GS_B"),
    gs_exact_source = rep("LOCAL_G2_FIXTURE", 4L),
    gene_symbol = c("GENE1", "GENE2", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )
}

g2_resource_definitions <- function() {
  list(
    fixture_lisa_core = list(
      version = "0.1.0", schema = "lisa_dictionary@2",
      artifact = "dictionary.tsv", data = g2_dictionary("custom")
    ),
    fixture_lisa_expanded = list(
      version = "0.1.0", schema = "lisa_dictionary@2",
      artifact = "dictionary.tsv", data = g2_dictionary("custom")
    ),
    fixture_msigdb_term2gene = list(
      version = "2026.1", schema = "term2gene@1",
      artifact = "term2gene.tsv", data = g2_term2gene()
    ),
    fixture_lisa_category_map = list(
      version = "0.1.1", schema = "category_map@1",
      artifact = "category_map.tsv", data = g2_category_map()
    ),
    lab_dictionary = list(
      version = "1.0.0", schema = "lisa_dictionary@2",
      artifact = "dictionary.tsv", data = g2_dictionary()
    ),
    lab_term2gene = list(
      version = "1.0.0", schema = "term2gene@1",
      artifact = "term2gene.tsv", data = g2_term2gene()
    ),
    lab_category_map = list(
      version = "1.0.0", schema = "category_map@1",
      artifact = "category_map.tsv", data = g2_category_map()
    )
  )
}

g2_resource_fixture <- function(
  species = "Homo sapiens",
  modality = "transcriptomic/genomic",
  compatibility = "lisaR>=0.6.0"
) {
  root <- tempfile("lisa-g2-resources-")
  cache <- file.path(root, "cache")
  definitions <- g2_resource_definitions()

  rows <- lapply(names(definitions), function(logical_id) {
    item <- definitions[[logical_id]]
    path <- file.path(cache, logical_id, item$version, item$artifact)
    g2_write_tsv(item$data, path)
    data.frame(
      logical_id = logical_id,
      version = item$version,
      species = species,
      modality = modality,
      schema = item$schema,
      sha256 = lisaR:::lisa_sha256_file(path),
      compatibility = compatibility,
      approved_origin = "Local resource fixture",
      artifact = item$artifact,
      stringsAsFactors = FALSE
    )
  })
  registry <- do.call(rbind, rows)
  rownames(registry) <- NULL

  refs <- stats::setNames(
    paste0(registry$logical_id, "@", registry$version),
    registry$logical_id
  )
  paths <- stats::setNames(
    file.path(cache, registry$logical_id, registry$version, registry$artifact),
    registry$logical_id
  )
  list(
    root = root,
    cache = cache,
    species = species,
    modality = modality,
    registry = registry,
    refs = refs,
    paths = paths
  )
}

# `schema` lets a test swap a fixture to the legacy `lisa_dictionary@1`
# contract. The declared schema and the table must agree , so a
# replacement that reintroduces the score column must say so.
g2_replace_resource_table <- function(fixture, logical_id, table, schema = NULL) {
  g2_write_tsv(table, fixture$paths[[logical_id]])
  selected <- fixture$registry$logical_id == logical_id
  fixture$registry$sha256[selected] <- lisaR:::lisa_sha256_file(
    fixture$paths[[logical_id]]
  )
  if (!is.null(schema)) fixture$registry$schema[selected] <- schema
  fixture
}

g2_resolve_resources <- function(pipeline = list(), fixture,
                                 stop_on_error = TRUE) {
  if (is.null(pipeline$profile)) {
    pipeline$profile <- "transcriptomic/genomic"
  }
  # These are deliberately non-builtin identities.  The C1 bundled resource
  # contract now rejects an external registry row that impersonates a bundled
  # identity, so tests select their synthetic fixtures explicitly.
  for (key in c("dictionary_resource", "term2gene_resource", "category_map_resource")) {
    if (is.null(pipeline[[key]])) {
      fixture_key <- c(
        dictionary_resource = "fixture_lisa_core",
        term2gene_resource = "fixture_msigdb_term2gene",
        category_map_resource = "fixture_lisa_category_map"
      )[[key]]
      pipeline[[key]] <- fixture$refs[[fixture_key]]
    }
  }
  lisaR:::lisa_resolve_pipeline_resources(
    pipeline = pipeline,
    species = fixture$species,
    registry = fixture$registry,
    cache_root = fixture$cache,
    shared_root = NULL,
    stop_on_error = stop_on_error
  )
}

g2_validate_custom_resources <- function(fixture,
                                         species = fixture$species) {
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))
  validate_lisa_custom_resources(
    dictionary_resource = fixture$refs[["lab_dictionary"]],
    category_map_resource = fixture$refs[["lab_category_map"]],
    term2gene_resource = fixture$refs[["lab_term2gene"]],
    species = species
  )
}
