read_audit_vignette <- function(name) {
  read_vignette_source(name)
}

test_that("getting started verifies the returned run before opening its report", {
  text <- read_audit_vignette("getting-started")
  verify_at <- regexpr("verify_lisa_run(result$output_dir)", text, fixed = TRUE)
  gate_at <- regexpr('stopifnot(identical(check$gate, "PASS"))', text, fixed = TRUE)
  browse_at <- regexpr("browseURL(file.path(result$output_dir", text, fixed = TRUE)
  expect_gt(as.integer(verify_at), 0L)
  expect_gt(as.integer(gate_at), as.integer(verify_at))
  expect_gt(as.integer(browse_at), as.integer(gate_at))
  expect_false(grepl('verify_lisa_run("results/example")', text, fixed = TRUE))
})

test_that("Riaz guide reads the prepared study and names its pinned resources", {
  yaml_path <- system.file("examples", "riaz-gse91061", "config",
                           "riaz-gse91061.yml", package = "lisaR", mustWork = TRUE)
  yaml <- paste(readLines(yaml_path, warn = FALSE), collapse = "\n")
  text <- read_audit_vignette("riaz-gse91061-worked-example")
  expect_match(text, 'config <- file.path(project, "study.yml")', fixed = TRUE)
  expect_match(text, "study <- yaml::read_yaml(config)", fixed = TRUE)
  expect_match(text, "11_prepare_prepared_project.R", fixed = TRUE)
  resources <- yaml::yaml.load(yaml)$pipeline
  for (resource in c("dictionary_resource", "term2gene_resource",
                     "category_map_resource")) {
    expect_match(text, resources[[resource]], fixed = TRUE)
  }
  expect_match(yaml, "GSE91061", fixed = TRUE)
  expect_match(yaml, "lisa_category_map@1.0.0", fixed = TRUE)
})

test_that("supercategory membership separates distinct sets from assignments", {
  root <- system.file("extdata", "semantic-catalog", package = "lisaR",
                      mustWork = TRUE)
  counts <- utils::read.delim(file.path(root, "semantic_supercategory_counts.tsv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  categories <- utils::read.delim(file.path(root, "semantic_category_catalog.tsv"),
                                 stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(counts), 42L)
  expect_false(anyDuplicated(paste(counts$universe, counts$macrogroup_id)) > 0L)
  expect_true(all(counts$core_distinct_gene_sets <= counts$core_member_assignments))
  # Independently audited overlapping memberships: summing category counts
  # gives 88/279, not the verified core union 59/252.
  pi3k <- counts[counts$macrogroup_id == "PI3K_AKT_MTOR", ]
  rtk <- counts[counts$macrogroup_id == "RTK_RAS_MAPK", ]
  expect_equal(pi3k$core_distinct_gene_sets, 59L)
  expect_equal(pi3k$core_member_assignments, 88L)
  expect_equal(rtk$core_distinct_gene_sets, 252L)
  expect_equal(rtk$core_member_assignments, 279L)
  text <- read_audit_vignette("semantic-category-guide")
  for (i in seq_len(nrow(counts))) {
    row <- counts[i, ]
    part <- categories[categories$universe == row$universe &
                         categories$macrogroup_id == row$macrogroup_id, ]
    expect_equal(row$category_count, nrow(part))
    expect_equal(row$core_member_assignments, sum(part$core_member_gene_sets))
    sentence <- sprintf(
      "Verified membership: %d %s, %d distinct core gene sets and %d category assignments.",
      row$category_count, if (row$category_count == 1L) "category" else "categories",
      row$core_distinct_gene_sets, row$core_member_assignments
    )
    expect_true(grepl(sentence, text, fixed = TRUE))
  }
  expect_false(grepl("Verified coverage:", text, fixed = TRUE))
})

test_that("documented category deltas remain descriptive even with formal source DE", {
  how <- gsub("[[:space:]]+", " ", read_audit_vignette("how-lisar-works"))
  report <- gsub("[[:space:]]+", " ", read_audit_vignette("scientific-report"))
  expect_match(how, "It is not a category-level test", fixed = TRUE)
  expect_match(how, "not a category false-discovery rate", fixed = TRUE)
  expect_match(report, "always a descriptive comparison", fixed = TRUE)
  expect_match(report, "not a formal interaction test", fixed = TRUE)
  expect_match(report, "remains gene-level evidence", fixed = TRUE)
  expect_false(grepl("unless the (source|corresponding)", report))
})

test_that("resource guide distinguishes cache-only maps from enrichment", {
  text <- gsub("[[:space:]]+", " ", read_audit_vignette("resources-and-provenance"))
  expect_match(text, "immutable `pipeline.kegg_snapshot_id`", fixed = TRUE)
  expect_match(text, "pipeline.kegg_access_mode: cache_only", fixed = TRUE)
  expect_match(text, "not fetch a native map automatically", fixed = TRUE)
  expect_match(text, "route authorised for", fixed = TRUE)
  expect_match(text, "KEGG gene-set enrichment remains a separate product", fixed = TRUE)
})
