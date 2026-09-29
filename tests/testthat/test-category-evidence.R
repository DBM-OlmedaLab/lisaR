evidence_fixture <- function() {
  list(gsea = data.frame(category_id = c(rep("CAT_MIXED", 5L), "CAT_NONE", "OTHER_UNCLASSIFIED"),
    category_display_name = c(rep("Mixed support", 5L), "No significant support", "Residual"),
    pathway = c("SET_A", "SET_B", "SET_C", "SET_D", "SET_EMPTY", "SET_NONE", "UNCLASSIFIED_SET"),
    NES = c(2, -1.5, .6, NA, 1.1, 1, 2), padj = c(.01, .02, .8, NA, .03, .9, .001),
    leadingEdge = c("GeneA;GeneB;GeneB", "GeneB;GeneC", "GeneD", NA, "", "GeneA", "RESIDUAL"),
    eligible_for_gsea = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    gsea_eligibility_status = c("eligible", "eligible", "eligible", "below_min_size", "eligible", "eligible", "eligible")),
    de = data.frame(symbol = c("GeneA", "GeneB", "GeneC", "UNRELATED"),
      gene_id = paste0("id", 1:4), log2FC = c(1, -.8, 2, .1), padj = c(.01, .2, .001, .8)))
}

test_that("evidence retains all states and separates gene and assignment counts", {
  x <- evidence_fixture()
  e <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, analysis_id = "A", collection = "GOBP-C2", tier = "core", max_sets = 2L, max_genes = 2L)
  c <- e$categories[e$categories$category_id == "CAT_MIXED", ]
  expect_equal(c$n_mapped_sets, 5L); expect_equal(c$n_evaluable_sets, 4L)
  expect_equal(c$n_significant_sets, 3L)
  expect_equal(c$n_significant_positive_sets, 2L); expect_equal(c$n_significant_negative_sets, 1L)
  expect_equal(c$n_unique_le_genes, 4L); expect_equal(c$n_le_assignments, 5L)
  expect_equal(c$n_significant_unique_le_genes, 3L); expect_equal(c$n_significant_le_assignments, 4L)
  expect_identical(c$support_state, "mixed_direction_support")
  expect_true(c$sets_truncated); expect_true(c$genes_truncated)
  expect_setequal(e$gene_selection$symbol, c("GeneA", "GeneB"))
  expect_true("GeneC" %in% e$leading_edges$symbol)
  expect_equal(e$sets$non_evaluable_reason[e$sets$pathway == "SET_D"], "below_min_size")
  expect_equal(e$sets$state[e$sets$pathway == "SET_C"], "non_significant")
  expect_equal(e$categories$support_state[e$categories$category_id == "CAT_NONE"], "no_significant_sets")
  expect_false(any(e$sets$category_id == "OTHER_UNCLASSIFIED"))
  expect_false("RESIDUAL" %in% e$genes$symbol)
  expect_true("UNRELATED" %in% e$genes$symbol)
  expect_equal(e$genes$de_state[e$genes$symbol == "GeneD"], "not_in_standardized_DE")
  expect_false(any(grepl("score|prioriti", names(e$genes))))
})

test_that("overlaps are descriptive unique membership with explicit empty and direction states", {
  x <- evidence_fixture(); e <- lisaR:::build_lisa_category_evidence(x$gsea, x$de)
  p <- lisaR:::lisa_category_evidence_overlaps(e, "CAT_MIXED")
  ab <- p[p$pathway_a == "SET_A" & p$pathway_b == "SET_B", ]
  expect_equal(nrow(p), choose(5, 2)); expect_equal(ab$intersection_n, 1L)
  expect_equal(ab$union_n, 3L); expect_equal(ab$jaccard, 1 / 3)
  expect_equal(ab$n_a, 2L); expect_equal(ab$n_b, 2L)
  expect_identical(ab$direction_relation, "opposite")
  empty <- p[p$pathway_b == "SET_EMPTY", ]
  expect_true(all(is.na(empty$jaccard))); expect_true(all(empty$overlap_state == "empty_or_unrecorded_membership"))
  expect_equal(nrow(lisaR:::lisa_category_evidence_overlaps(e, "CAT_NONE")), 0L)
  member <- data.frame(gs_name = c("SET_A", "SET_A", "SET_B", "SET_B"), gene_symbol = c("M1", "M2", "M2", "M3"))
  e2 <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, memberships = member)
  pc <- lisaR:::lisa_category_evidence_overlaps(e2, "CAT_MIXED", "complete")
  expect_equal(pc$jaccard[pc$pathway_a == "SET_A" & pc$pathway_b == "SET_B"], 1 / 3)
  expect_identical(pc$membership_type[[1L]], "complete")
  expect_false(any(e2$leading_edges$symbol %in% member$gene_symbol))
  chunks <- list()
  lisaR:::lisa_category_evidence_overlaps(e, "CAT_MIXED", emit = function(chunk) chunks[[length(chunks) + 1L]] <<- chunk)
  streamed <- do.call(rbind, chunks); rownames(streamed) <- NULL
  expect_identical(streamed, p)
})

test_that("dictionary plus ledger restores missing results and zero-support categories", {
  x <- evidence_fixture()
  dict <- data.frame(category_id = c("CAT_MIXED", "CAT_EMPTY"), gene_set_id = c("MISSING_RESULT", "NO_RESULT"), description = c("Category description", "No support"))
  ledger <- data.frame(pathway = c(x$gsea$pathway, "MISSING_RESULT", "NO_RESULT"),
    eligible_for_gsea = c(x$gsea$eligible_for_gsea, FALSE, FALSE),
    gsea_eligibility_status = c(x$gsea$gsea_eligibility_status, "no_rank_overlap", "above_max_size"))
  e <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, universe_ledger = ledger, category_dictionary = dict)
  expect_true("CAT_EMPTY" %in% e$categories$category_id)
  expect_equal(e$sets$non_evaluable_reason[e$sets$pathway == "MISSING_RESULT"], "no_rank_overlap")
  expect_equal(e$categories$category_description[e$categories$category_id == "CAT_MIXED"], "Category description")
  expect_equal(e$categories$support_state[e$categories$category_id == "CAT_EMPTY"], "no_evaluable_sets")
})

test_that("exact joins, conflicting analyses and invalid caps fail closed", {
  x <- evidence_fixture()
  lower <- x$de; lower$symbol[[1L]] <- "genea"
  e <- lisaR:::build_lisa_category_evidence(x$gsea, lower)
  expect_equal(e$genes$de_state[e$genes$symbol == "GeneA"], "not_in_standardized_DE")
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, rbind(x$de, x$de[1, ])), "duplicate symbols")
  duplicate <- x$gsea[1, ]; duplicate$NES <- 9; duplicate$category_id <- "SECOND_CATEGORY"
  expect_error(lisaR:::build_lisa_category_evidence(rbind(x$gsea, duplicate), x$de), "Conflicting GSEA")
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, max_sets = 0), "positive integers")
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, gsea_padj_cutoff = NA_real_), "FDR cutoffs")
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, gsea_padj_cutoff = 1.1), "FDR cutoffs")
})

test_that("normalised tables and display selection are deterministic under row permutation", {
  x <- evidence_fixture()
  a <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, max_sets = 2L, max_genes = 2L)
  b <- lisaR:::build_lisa_category_evidence(x$gsea[rev(seq_len(nrow(x$gsea))), ], x$de[4:1, ], max_sets = 2L, max_genes = 2L)
  for (name in c("categories", "sets", "genes", "leading_edges", "gene_selection")) expect_identical(a[[name]], b[[name]])
})

test_that("requested evidence tier cannot relabel broader input assignments", {
  x <- evidence_fixture()
  x$gsea$tier <- "core"
  expect_identical(lisaR:::build_lisa_category_evidence(x$gsea, x$de, tier = "core")$metadata$tier, "core")
  for (broader in c("expanded", "broad")) {
    mixed <- x$gsea; mixed$tier[[1L]] <- broader
    expect_error(lisaR:::build_lisa_category_evidence(mixed, x$de, tier = "core"), "tier scope.*core.*input row tiers")
  }
  expanded <- x$gsea; expanded$tier[[1L]] <- "expanded"
  expect_identical(lisaR:::build_lisa_category_evidence(expanded, x$de, tier = "expanded")$metadata$tier, "expanded")
  dictionary <- data.frame(category_id = "CAT_MIXED", tier = "expanded")
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, category_dictionary = dictionary, tier = "core"), "tier scope")
  expect_identical(lisaR:::build_lisa_category_evidence(x$gsea, x$de, category_dictionary = dictionary, tier = "expanded")$metadata$tier, "expanded")
  unscoped <- evidence_fixture()
  custom <- lisaR:::build_lisa_category_evidence(unscoped$gsea, unscoped$de, tier = "core (synthetic)", collection = "SYNTHETIC")
  expect_identical(custom$metadata$tier, "core (synthetic)")
  expect_identical(custom$metadata$collection, "SYNTHETIC")
})

test_that("evidence checks every available analysis collection and universe", {
  x <- evidence_fixture()
  x$gsea$analysis_collection <- "GOBP-C2"; x$gsea$universe <- "GOBP-C2"
  ledger <- data.frame(pathway = x$gsea$pathway, analysis_collection = "GOBP-C2")
  dictionary <- data.frame(category_id = "CAT_MIXED", universe = "GOBP-C2")
  e <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, universe_ledger = ledger, category_dictionary = dictionary, collection = "GOBP-C2")
  expect_identical(e$metadata$collection, "GOBP-C2")
  expect_identical(lisaR:::build_lisa_category_evidence(x$gsea, x$de)$metadata$collection, "GOBP-C2")
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, collection = "GOCC"), "collection scope")
  mismatch <- x$gsea; mismatch$universe[[1L]] <- "GOCC"
  expect_error(lisaR:::build_lisa_category_evidence(mismatch, x$de, collection = "GOBP-C2"), "collection scope")
  ledger$analysis_collection[[1L]] <- "GOMF"
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, universe_ledger = ledger, collection = "GOBP-C2"), "collection scope")
  dictionary$universe <- "PATHWAYS"
  expect_error(lisaR:::build_lisa_category_evidence(x$gsea, x$de, category_dictionary = dictionary, collection = "GOBP-C2"), "collection scope")
})

test_that("duplicate enrichment statistics differing beyond 15 digits still conflict", {
  x <- evidence_fixture()
  for (column in c("NES", "padj")) {
    input <- x$gsea[rep(1L, 2L), , drop = FALSE]
    input$category_id <- c("CAT_ONE", "CAT_TWO")
    input[[column]] <- c(1, 1 + .Machine$double.eps)
    expect_error(lisaR:::build_lisa_category_evidence(input, x$de), paste0("Conflicting GSEA ", column))
    input[[column]][[2L]] <- input[[column]][[1L]]
    expect_equal(nrow(lisaR:::build_lisa_category_evidence(input, x$de)$sets), 2L)
  }
})

test_that("portable render exports full edges, explicit source and offline recipe without DE duplication", {
  skip_if_not(nzchar(system.file("category-evidence", package = "lisaR")))
  x <- evidence_fixture(); x$gsea$category_display_name[[1L]] <- "</script><script>bad()</script>"
  e <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, max_sets = 2L, max_genes = 1L)
  out <- tempfile("category-evidence-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  r <- lisaR:::render_lisa_category_evidence(e, out, formats = c("png", "pdf", "svg"))
  expect_true(file.exists(r$html)); expect_true(file.exists(file.path(out, "reproduce_category_evidence.R")))
  expect_true(all(file.exists(file.path(out, "figures", paste0("001_category_evidence.", c("png", "pdf", "svg"))))))
  expect_identical(r$member_figure_index$category_id, e$categories$category_id)
  expect_true(all(file.exists(file.path(out, "figures", paste0("001_category_members.", c("png", "pdf", "svg"))))))
  expect_true(all(file.exists(file.path(out, r$member_figure_index$source_tsv))))
  expect_true(all(file.exists(file.path(out, r$member_figure_index$recipe_r))))
  html <- paste(readLines(r$html, warn = FALSE), collapse = "\n")
  expect_false(grepl("</script><script>bad()", html, fixed = TRUE))
  expect_false(grepl("fetch(", html, fixed = TRUE)); expect_false(grepl("https://", html, fixed = TRUE))
  expect_match(html, "evidence-data", fixed = TRUE)
  edge <- utils::read.delim(file.path(out, "tables", "leading_edges.tsv"))
  expect_equal(nrow(edge), nrow(e$leading_edges)); expect_true("GeneC" %in% edge$symbol)
  expect_equal(length(list.files(out, "^genes[.]tsv$", recursive = TRUE)), 1L)
  source <- utils::read.delim(file.path(out, "figures", "001_category_evidence_source.tsv"))
  expect_equal(sum(source$row_type == "gene"), 1L)
  expect_equal(sum(source$row_type == "set"), 2L)
  expect_true(all(c("NES", "gsea_fdr", "log2FC", "de_fdr", "positive_contrast") %in% names(source)))
  pairs <- utils::read.delim(file.path(out, "tables", "displayed_leading_edge_overlaps.tsv"))
  expect_equal(nrow(pairs), 1L)
  expect_identical(pairs$pair_scope, "default_displayed_sets")
  expect_false(file.exists(file.path(out, "tables", "all_leading_edge_overlaps.tsv")))
  minimal <- tempfile("category-evidence-min-"); on.exit(unlink(minimal, recursive = TRUE), add = TRUE)
  minimal_result <- lisaR:::render_lisa_category_evidence(e, minimal, formats = character(), source_data = FALSE, recipes = FALSE)
  expect_false(dir.exists(file.path(minimal, "figures")))
  expect_equal(nrow(minimal_result$member_figure_index), 0L)
  expect_false(file.exists(file.path(minimal, "reproduce_category_evidence.R")))
  expect_true(file.exists(file.path(minimal, "tables", "leading_edges.tsv")))
  all_out <- tempfile("category-evidence-all-"); on.exit(unlink(all_out, recursive = TRUE), add = TRUE)
  all_result <- lisaR:::render_lisa_category_evidence(e, all_out, formats = character(), overlap_export = "all")
  # Explicit no-figure mode takes precedence even with default sidecar flags.
  expect_false(dir.exists(file.path(all_out, "figures")))
  expect_equal(nrow(all_result$member_figure_index), 0L)
  expect_false(file.exists(file.path(all_out, "reproduce_category_evidence.R")))
  all_html <- readLines(all_result$html, warn = FALSE)
  data <- all_html[grepl('id="evidence-data"', all_html, fixed = TRUE)]
  payload <- sub('^.*id="evidence-data">', '', data)
  payload <- sub('</script>$', '', payload)
  decoded <- jsonlite::fromJSON(payload)
  expect_length(decoded$formats, 0L)
  expect_length(decoded$member_figure_index, 0L)
  all_pairs <- utils::read.delim(file.path(all_out, "tables", "all_leading_edge_overlaps.tsv"))
  expect_equal(nrow(all_pairs), 10L)
  expect_true(all(all_pairs$pair_scope == "all_mapped_sets"))
})

test_that("category metadata retains empty categories and JSON keeps full cutoff precision", {
  x <- evidence_fixture()
  summary <- data.frame(category_id = c("CAT_MIXED", "CAT_NO_MAPPED"), category_display_name = c("Mixed support", "No mapped sets"))
  threshold <- 0.12345678901234567
  e <- lisaR:::build_lisa_category_evidence(x$gsea, NULL, category_dictionary = summary, gsea_padj_cutoff = threshold)
  expect_equal(e$categories$support_state[e$categories$category_id == "CAT_NO_MAPPED"], "no_mapped_sets")
  expect_true(all(e$genes$de_state == "not_in_standardized_DE"))
  skip_if_not(nzchar(system.file("category-evidence", package = "lisaR")))
  out <- tempfile("category-evidence-precision-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- lisaR:::render_lisa_category_evidence(e, out, formats = character(), overlap_export = "none")
  html <- readLines(result$html, warn = FALSE)
  data <- html[grepl('id="evidence-data"', html, fixed = TRUE)]
  payload <- sub('^.*id="evidence-data">', '', data)
  payload <- sub('</script>$', '', payload)
  decoded <- jsonlite::fromJSON(payload)
  expect_identical(decoded$metadata$gsea_padj_cutoff, threshold)
  expect_false(any(grepl("overlaps.tsv", list.files(file.path(out, "tables")), fixed = TRUE)))
})

test_that("portable evidence uses exact relative assets and passes the actual privacy gate", {
  skip_if_not(nzchar(system.file("category-evidence", package = "lisaR")))
  x <- evidence_fixture()
  e <- lisaR:::build_lisa_category_evidence(x$gsea, x$de, analysis_id = "privacy_fixture", collection = "GOBP-C2")
  out <- tempfile("category-evidence-private-gate-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- lisaR:::render_lisa_category_evidence(e, out, formats = character())
  html <- paste(readLines(result$html, warn = FALSE), collapse = "\n")
  expect_match(html, 'href="assets/category-evidence.css"', fixed = TRUE)
  expect_match(html, 'src="assets/category-evidence.js"', fixed = TRUE)
  # The emitted data are still embedded; no browser fetch of local JSON is used.
  expect_match(html, 'id="evidence-data"', fixed = TRUE)
  for (extension in c("css", "js")) {
    source <- system.file("category-evidence", paste0("viewer.", extension), package = "lisaR")
    copied <- file.path(out, "assets", paste0("category-evidence.", extension))
    expect_true(file.exists(copied))
    expect_identical(digest::digest(file = copied, algo = "sha256"), digest::digest(file = source, algo = "sha256"))
  }
  report <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = report)
  # This is the exact gate that rejected the preserved installed QuickStart,
  # not a weakened or approximate regular expression used only by this test.
  expect_no_error(report$report_assert_no_private_paths(out, result$html))
  e$metadata$positive_contrast <- "/srv/private-study/design.tsv"
  unsafe <- lisaR:::render_lisa_category_evidence(e, out, formats = character())
  expect_error(report$report_assert_no_private_paths(out, unsafe$html), "LISA-REPORT-PRIVACY-004")
})

test_that("category product payloads require exact IDs and safe canonical assets", {
  product <- list(list(category_id = "CAT_ONE", product = "volcano", label = "Volcano overlay",
    assets = list(list(format = "png", href = "assets/category_products/volcano.png"))))
  normalized <- lisaR:::lisa_category_evidence_category_products(product, "CAT_ONE")
  expect_identical(normalized[[1L]]$category_id, "CAT_ONE")
  canonical <- product
  canonical[[1L]]$assets[[1L]]$href <-
    "../../../../outputs/gene_level/A/category%20cards/CAT_ONE.png"
  expect_no_error(lisaR:::lisa_category_evidence_category_products(
    canonical, "CAT_ONE"
  ))
  wrong_id <- product; wrong_id[[1L]]$category_id <- "OTHER"
  expect_error(lisaR:::lisa_category_evidence_category_products(wrong_id, "CAT_ONE"), "exact evidence category")
  outside <- product; outside[[1L]]$assets[[1L]]$href <- "../outside.png"
  expect_error(lisaR:::lisa_category_evidence_category_products(outside, "CAT_ONE"), "safe relative file")
})

test_that("emitted standalone figure recipe runs and preserves exact numeric-looking IDs", {
  skip_if_not(nzchar(system.file("category-evidence", package = "lisaR")))
  gsea <- data.frame(category_id = "001", category_display_name = "Numeric-looking exact identifiers",
    pathway = c("0010", "0020"), NES = c(1.25, -1.1), padj = c(.001, .02),
    leadingEdge = c("0008;0009", "0009"))
  de <- data.frame(symbol = c("0008", "0009"), gene_id = c("00018", "00019"),
    log2FC = c(.75, -1), padj = c(.01, .03))
  e <- lisaR:::build_lisa_category_evidence(gsea, de, analysis_id = "0007", collection = "synthetic",
    tier = "custom", positive_contrast = "synthetic treated minus control")
  out <- tempfile("category-evidence-real-recipe-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  lisaR:::render_lisa_category_evidence(e, out, formats = "png")
  recipe <- file.path(out, "reproduce_category_evidence.R")
  source <- file.path(out, "figures", "001_category_evidence_source.tsv")
  original <- file.path(out, "figures", "001_category_evidence.png")
  reproduced <- file.path(out, "reproduced.png")
  log <- file.path(out, "recipe-cli.log")
  executable <- file.path(R.home("bin"), "Rscript")
  status <- system2(executable, c("--vanilla", shQuote(recipe), shQuote(source), shQuote(reproduced)),
    stdout = log, stderr = log)
  expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  expect_true(file.exists(reproduced))
  if (file.exists(reproduced) && requireNamespace("png", quietly = TRUE)) {
    # Compare decoded image data, not PNG container metadata or file bytes.
    expect_identical(png::readPNG(reproduced), png::readPNG(original))
  }
})


test_that("evidence default ordering uses all significant-set recurrence before caps", {
  g <- data.frame(category_id = "CAT", pathway = c("SET_Z", "SET_A", "SET_B", "SET_NONSIG"),
    NES = c(-3, 2, 1, 4), padj = c(.01, .01, .02, .8),
    leadingEdge = c("GeneZ;GeneA", "GeneB;GeneA", "GeneZ", "GeneA;GeneB;GeneC"))
  e <- lisaR:::build_lisa_category_evidence(g, max_sets = 1L, max_genes = 1L)
  expect_identical(e$sets$pathway[!is.na(e$sets$display_order)], "SET_Z")
  # GeneA/GeneZ each occur twice among ALL significant sets, so exact symbol ties.
  expect_identical(e$gene_selection$symbol, "GeneA")
  expect_equal(e$gene_selection$n_significant_set_connections, 2L)
  expect_equal(e$gene_selection$n_significant_sets_denominator, 3L)
  expect_equal(e$categories$mean_significant_set_NES, 0)
  expect_equal(e$categories$n_available_plot_genes, 2L)
  expect_equal(nrow(e$leading_edges), 8L)
  # Add support outside the capped visible set: recurrence, not visible count,
  # must switch the selected gene without changing the significant set cutoff.
  g$leadingEdge[[2L]] <- "GeneB;GeneZ"
  e2 <- lisaR:::build_lisa_category_evidence(g, max_sets = 1L, max_genes = 1L)
  expect_identical(e2$gene_selection$symbol, "GeneZ")
  expect_equal(e2$gene_selection$n_significant_set_connections, 3L)
  expect_equal(e2$gene_selection$n_significant_sets_denominator, 3L)
})

test_that("evidence alignment verifies the canonical category mean, not the capped mean", {
  g <- data.frame(category_id = "CAT", pathway = c("A", "B"), NES = c(3, -1),
    padj = c(.01, .02), leadingEdge = c("G1", "G2"))
  summary <- data.frame(category_id = "CAT", mean_NES = 1)
  e <- lisaR:::build_lisa_category_evidence(g, category_dictionary = summary, max_sets = 1L)
  expect_equal(e$categories$mean_significant_set_NES, 1)
  expect_identical(e$categories$category_mean_source, "verified_category_summary_mean_NES")
  expect_identical(e$categories$support_state, "mixed_direction_support")
  summary$mean_NES <- 3
  expect_error(lisaR:::build_lisa_category_evidence(g, category_dictionary = summary), "Category mean_NES conflicts")
})


test_that("file-backed evidence preserves numeric-looking exact identifiers", {
  f <- tempfile(fileext = ".tsv"); on.exit(unlink(f), add = TRUE)
  writeLines(c("category_id\tpathway\tNES\tpadj\tleadingEdge", "001\t0002\t2\t0.01\t0003"), f)
  e <- lisaR:::build_lisa_category_evidence(f)
  expect_identical(e$categories$category_id, "001")
  expect_identical(e$sets$pathway, "0002")
  expect_identical(e$genes$symbol, "0003")
  expect_equal(e$sets$NES, 2)
})
