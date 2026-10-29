# Four synthetic sets exercise the saved report-evidence boundary, not fgsea.
# The expectations below come from explicit finite/non-finite test evidence;
# fixture counters are not calculated by the acceptance implementation.
empty_eligibility_members <- function() {
  data.frame(category_id = "CAT_SYNTHETIC",
    pathway = c("GS_ELIGIBLE", "GS_BELOW", "GS_ABOVE", "GS_NO_OVERLAP"),
    classification_status = "classified", analysis_id = "A", collection = "HALLMARKS",
    analysis_collection = "HALLMARKS", NES = c(1.5, NA, NA, NA),
    padj = c(.8, NA, NA, NA), eligible_for_gsea = c(TRUE, FALSE, FALSE, FALSE),
    gsea_eligibility_status = c("eligible", "below_min_size", "above_max_size", "no_rank_overlap"),
    gsea_result_status = c("tested", "below_min_size", "above_max_size", "no_rank_overlap"),
    source_gene_count = c(2L, 1L, 4L, 1L),
    ranked_gene_count = c(2L, 1L, 4L, 0L), min_gs_size = 2L, max_gs_size = 3L)
}

empty_eligibility_run <- function(members = empty_eligibility_members(), mapped = 4L,
    evaluable = 1L, significant = 0L, mean = NA_real_) {
  root <- tempfile("public-empty-eligibility-")
  scope <- file.path(root, "outputs/single_de/A/collection_HALLMARKS")
  s <- file.path(scope, "lisa_tables/A_synthetic_GSEA_category_summary.tsv")
  g <- file.path(scope, "enrichment/A_GSEA_synthetic_annotated.tsv")
  summary <- data.frame(category_id = unique(members$category_id[members$classification_status == "classified"]),
    n_genesets = significant, n_genesets_mapped = mapped, n_genesets_evaluable = evaluable,
    n_genesets_significant = significant, mean_NES = mean, gsea_padj_cutoff = .05,
    analysis_id = "A", collection = "HALLMARKS")
  write <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    lisaR:::write_lisa_tsv(x, path)
  }
  write(summary, s); write(members, g)
  write(data.frame(key = "gsea_padj_cutoff", value = .05), file.path(root, "contract_manifest.tsv"))
  write(data.frame(analysis_id = "A", species = "Homo sapiens"), file.path(root, "config/de_index.tsv"))
  write(data.frame(contrast_id = character(), output_id = character()), file.path(root, "config/contrast_index.tsv"))
  list(root = root, summary = s, gsea = g)
}

empty_eligibility_state <- function(f) lisaR:::lisa_report_selection_state(f$root, "A", "HALLMARKS")

test_that("saved eligible non-significant evidence tolerates only justified pretest exclusions", {
  x <- empty_eligibility_members()
  fixtures <- list(empty_eligibility_run(x))
  for (i in 2:4) fixtures[[length(fixtures) + 1L]] <-
    empty_eligibility_run(x[c(1L, i), ], mapped = 2L)
  for (f in fixtures) {
    before <- tools::md5sum(c(f$gsea, f$summary))
    expect_identical(empty_eligibility_state(f), "completed_empty")
    expect_identical(tools::md5sum(c(f$gsea, f$summary)), before)
    observed <- lisaR:::lisa_category_nes_statistics(lisaR:::read_lisa_tsv(f$gsea), .05)
    expect_equal(observed$n_genesets_evaluable, 1L)
    expect_equal(observed$n_genesets_significant, 0L)
    expect_equal(observed$n_genesets_unavailable, observed$n_genesets_mapped - 1L)
  }
  f <- fixtures[[1L]]
  inventory <- lisaR:::lisa_full_product_inventory(f$root)
  expect_setequal(inventory$family, c("prioritized_genes", "volcano", "gene_heatmap"))
  expect_true(all(inventory$state == "completed_empty"))
})

test_that("actual public ledger completion yields the mixed empty-selection evidence", {
  membership <- list(GS_ELIGIBLE = c("G1", "G2"), GS_BELOW = "G1",
    GS_ABOVE = c("G1", "G2", "G3", "G4"), GS_NO_OVERLAP = "GX")
  term <- do.call(rbind, lapply(names(membership), function(id) data.frame(
    gs_collection = "H", gs_subcollection = "SYNTHETIC", gs_name = id,
    gs_exact_source = "SYNTHETIC", gene_symbol = membership[[id]])))
  dictionary <- data.frame(gene_set_id = names(membership), category_id = "CAT_SYNTHETIC")
  ranks <- stats::setNames(4:1, paste0("G", 1:4))
  ledger <- lisaR:::build_lisa_gsea_ledger(term, dictionary, ranks, 2L, 3L, "HALLMARKS")
  result <- data.frame(pathway = "GS_ELIGIBLE", NES = 1.5, padj = .8)
  completed <- lisaR:::complete_lisa_gsea_results(result, ledger)
  completed$category_id <- "CAT_SYNTHETIC"
  f <- empty_eligibility_run(completed)
  expect_identical(empty_eligibility_state(f), "completed_empty")
  expect_equal(sum(completed$eligible_for_gsea), 1L)
  expect_equal(sum(is.na(completed$NES) & is.na(completed$padj)), 3L)
})

test_that("positive members remain required and wholly excluded scopes are not negative tests", {
  x <- empty_eligibility_members(); x$padj[[1L]] <- .05
  f <- empty_eligibility_run(x, significant = 1L, mean = 1.5)
  expect_identical(empty_eligibility_state(f), "missing")
  expect_error(lisaR:::lisa_full_product_inventory(f$root), "FULL required family absent")
  for (indices in list(2:4, 2L, 3L, 4L)) {
    f <- empty_eligibility_run(empty_eligibility_members()[indices, ],
      mapped = length(indices), evaluable = 0L)
    expect_error(empty_eligibility_state(f), "no eligible tests")
  }
  # This separate presentation boundary must not be broadened by the helper.
  f <- empty_eligibility_run()
  old <- options(lisaR.run_root = f$root); on.exit(options(old), add = TRUE)
  report <- list(mode = "full", formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
    source_data = TRUE, recipes = TRUE)
  expect_error(lisaR:::lisa_render_run_full_products(f$root, report), "LISA-EXTENSION-009")
  expect_false(dir.exists(file.path(f$root, "artifacts")))
})

test_that("tested non-finite statistics are never justified by other excluded members", {
  for (field in c("NES", "padj")) for (value in c(NA_real_, NaN, Inf, -Inf)) {
    x <- empty_eligibility_members(); x[[field]][[1L]] <- value
    f <- empty_eligibility_run(x, evaluable = 0L)
    expect_error(empty_eligibility_state(f), "unavailable statistics")
  }
})

test_that("eligibility cannot be inferred from missing fields or unknown values", {
  fields <- c("eligible_for_gsea", "gsea_eligibility_status", "gsea_result_status")
  for (field in fields) {
    x <- empty_eligibility_members(); x[[field]] <- NULL
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "recorded GSEA eligibility")
  }
  x <- empty_eligibility_members(); x[fields] <- NULL
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "unavailable statistics|eligibility|size bounds")
  for (value in list(NA, "", "false", "FALSE ", "unknown", 0L)) {
    x <- empty_eligibility_members(); x$eligible_for_gsea[[2L]] <- value
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "recorded GSEA eligibility")
  }
  for (field in fields[-1L]) for (value in c(NA_character_, "", "unknown", "not_tested", "BELOW_MIN_SIZE")) {
    x <- empty_eligibility_members(); x[[field]][[2L]] <- value
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "recorded GSEA eligibility")
  }
})

test_that("explicit eligibility contradictions reject even when all statistics are finite", {
  mutations <- list(
    function(x) { x$eligible_for_gsea[[2L]] <- TRUE; x },
    function(x) { x$gsea_result_status[[2L]] <- "tested"; x },
    function(x) { x$gsea_result_status[[2L]] <- "above_max_size"; x },
    function(x) { x$gsea_eligibility_status[[1L]] <- "below_min_size"; x },
    function(x) { x$eligible_for_gsea[[1L]] <- FALSE; x$gsea_eligibility_status[[1L]] <- "below_min_size"; x },
    function(x) { x$NES[[2L]] <- 1; x },
    function(x) { x$padj[[2L]] <- .8; x },
    function(x) { x$NES[[2L]] <- Inf; x })
  for (mutate in mutations) {
    x <- mutate(empty_eligibility_members())
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "contradictory|eligibility")
  }
  x <- empty_eligibility_members()[1L, ]; x$eligible_for_gsea <- FALSE
  x$gsea_eligibility_status <- "below_min_size"
  expect_error(empty_eligibility_state(empty_eligibility_run(x, mapped = 1L)), "contradictory")
})

test_that("structural exclusions require complete coherent saved size evidence", {
  for (field in c("source_gene_count", "ranked_gene_count", "min_gs_size", "max_gs_size")) {
    x <- empty_eligibility_members(); x[[field]] <- NULL
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
    for (value in c(NA_real_, Inf, -1, .5)) {
      x <- empty_eligibility_members(); x[[field]][[2L]] <- value
      expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
    }
  }
  for (change in list(c(2L, 2L), c(3L, 3L), c(4L, 1L), c(1L, 1L))) {
    x <- empty_eligibility_members(); x$ranked_gene_count[[change[[1L]]]] <- change[[2L]]
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
  }
  x <- empty_eligibility_members(); x$min_gs_size[[2L]] <- 0
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
  x <- empty_eligibility_members(); x$max_gs_size[[2L]] <- 1
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
  x <- empty_eligibility_members(); x$source_gene_count[[1L]] <- 1
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
  x <- empty_eligibility_members(); x$min_gs_size[[2L]] <- 3
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
})

test_that("identical repeats and dictionary fan-out preserve support but metadata conflicts reject", {
  x <- empty_eligibility_members()
  f <- empty_eligibility_run(rbind(x, x))
  expect_identical(empty_eligibility_state(f), "completed_empty")
  second <- x; second$category_id <- "CAT_SECOND"
  f <- empty_eligibility_run(rbind(x, second))
  expect_identical(empty_eligibility_state(f), "completed_empty")
  x <- rbind(x, x[1L, ]); x$ranked_gene_count[[5L]] <- 3L; x$source_gene_count[[5L]] <- 3L
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "conflicting duplicate GSEA eligibility")
  x <- empty_eligibility_members(); second$ranked_gene_count[[1L]] <- 3L; second$source_gene_count[[1L]] <- 3L
  expect_error(empty_eligibility_state(empty_eligibility_run(rbind(x, second))), "conflicting duplicate GSEA eligibility")
  x <- rbind(empty_eligibility_members(), empty_eligibility_members()[1L, ])
  x$NES[[5L]] <- 2
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "Conflicting duplicate GSEA member")
})

test_that("historical finite evidence stays supported without inventing eligibility for NA", {
  x <- empty_eligibility_members()[1L, ]
  x[c("eligible_for_gsea", "gsea_eligibility_status", "gsea_result_status",
    "source_gene_count", "ranked_gene_count", "min_gs_size", "max_gs_size")] <- NULL
  expect_identical(empty_eligibility_state(empty_eligibility_run(x, mapped = 1L)), "completed_empty")
  x$NES <- NA_real_
  expect_error(empty_eligibility_state(empty_eligibility_run(x, mapped = 1L, evaluable = 0L)), "unavailable statistics")
  x <- empty_eligibility_members(); names(x)[names(x) == "padj"] <- "gsea_fdr"
  names(x)[names(x) == "pathway"] <- "gene_set_id"
  expect_identical(empty_eligibility_state(empty_eligibility_run(x)), "completed_empty")
})

test_that("mixed exclusions do not weaken cutoff owner schema and summary evidence guards", {
  for (field in c("analysis_id", "analysis_collection", "collection")) {
    x <- empty_eligibility_members(); x[[field]][[2L]] <- "FOREIGN"
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "foreign scope evidence")
  }
  for (field in c("n_genesets", "n_genesets_mapped", "n_genesets_evaluable", "n_genesets_significant", "mean_NES", "gsea_padj_cutoff")) {
    f <- empty_eligibility_run(); s <- lisaR:::read_lisa_tsv(f$summary)
    s[[field]] <- if (field == "gsea_padj_cutoff") .1 else 7
    lisaR:::write_lisa_tsv(s, f$summary)
    expect_error(empty_eligibility_state(f), "LISA-REPORT-EMPTY-001")
  }
  for (field in c("NES", "padj")) {
    x <- empty_eligibility_members(); x[[field]][[2L]] <- "corrupt"
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "Invalid numeric")
  }
  for (value in c(-.1, 1.1)) {
    x <- empty_eligibility_members(); x$padj[[1L]] <- value
    expect_error(empty_eligibility_state(empty_eligibility_run(x)), "\\[0,1\\]")
  }
  for (which in c("gsea", "summary")) {
    f <- empty_eligibility_run(); tab <- lisaR:::read_lisa_tsv(f[[which]])
    duplicate <- if (which == "gsea") "eligible_for_gsea" else "n_genesets_mapped"
    tab[[ncol(tab) + 1L]] <- tab[[duplicate]]; names(tab)[ncol(tab)] <- duplicate
    lisaR:::write_lisa_tsv(tab, f[[which]])
    expect_error(empty_eligibility_state(f), "duplicate evidence columns")
  }
  f <- empty_eligibility_run(); s <- lisaR:::read_lisa_tsv(f$summary)
  s$category_id <- "OTHER_CATEGORY"; lisaR:::write_lisa_tsv(s, f$summary)
  expect_error(empty_eligibility_state(f), "omits mapped members")
  for (which in c("gsea", "summary")) {
    f <- empty_eligibility_run(); unlink(f[[which]])
    expect_error(empty_eligibility_state(f), "missing or ambiguous")
  }
})

test_that("filtering and positive products do not hide contradictory raw eligibility", {
  x <- empty_eligibility_members(); x$padj[[1L]] <- .01
  x$gsea_result_status[[2L]] <- "tested"
  expect_error(empty_eligibility_state(empty_eligibility_run(x, significant = 1L, mean = 1.5)), "contradictory")
  x <- empty_eligibility_members(); hidden <- x[1L, ]
  hidden$category_id <- "OTHER_UNCLASSIFIED"; hidden$classification_status <- "unclassified"
  hidden$pathway <- "GS_BACKGROUND"; hidden$gsea_result_status <- "below_min_size"
  expect_error(empty_eligibility_state(empty_eligibility_run(rbind(x, hidden))), "contradictory")
  x <- empty_eligibility_members(); hidden$gsea_result_status <- "tested"
  f <- empty_eligibility_run(rbind(x, hidden))
  expect_identical(empty_eligibility_state(f), "completed_empty")
  # Finite background cannot make a classified reporting scope evaluable.
  f <- empty_eligibility_run(rbind(x[2:4, ], hidden), mapped = 3L, evaluable = 0L)
  expect_error(empty_eligibility_state(f), "no eligible tests")
  f <- empty_eligibility_run()
  for (which in c("gsea", "summary")) {
    tab <- lisaR:::read_lisa_tsv(f[[which]])
    lisaR:::write_lisa_tsv(tab[FALSE, , drop = FALSE], f[[which]])
  }
  expect_error(empty_eligibility_state(f), "no eligible tests")
})

test_that("recorded bounds are validated without introducing a new scientific minimum", {
  x <- empty_eligibility_members()[1L, ]; x$min_gs_size <- 0L
  expect_identical(empty_eligibility_state(empty_eligibility_run(x, mapped = 1L)), "completed_empty")
  x <- empty_eligibility_members(); x$min_gs_size <- 0L
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "size bounds")
  x$eligible_for_gsea[[4L]] <- TRUE; x$gsea_result_status[[4L]] <- "tested"
  expect_error(empty_eligibility_state(empty_eligibility_run(x)), "contradictory")
})

test_that("NaN-encoded excluded evidence is not a structural absence", {
  f <- empty_eligibility_run()
  x <- lisaR:::read_lisa_tsv(f$gsea)
  x$NES[[2L]] <- NaN
  # Preserve the NaN token, unlike the package writer's intentional NA output.
  utils::write.table(x, f$gsea, sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
  expect_error(empty_eligibility_state(f), "contradictory statistics")
})

test_that("explicit finite ledger evidence cannot omit its saved sizes or eligibility triad", {
  bounds <- c("source_gene_count", "ranked_gene_count", "min_gs_size", "max_gs_size")
  for (drop in c(list(bounds), as.list(bounds))) {
    x <- empty_eligibility_members()[1L, ]; x[drop] <- NULL
    expect_error(empty_eligibility_state(empty_eligibility_run(x, mapped = 1L)), "size bounds")
  }
  x <- empty_eligibility_members()[1L, ]
  x[c("eligible_for_gsea", "gsea_eligibility_status", "gsea_result_status")] <- NULL
  expect_error(empty_eligibility_state(empty_eligibility_run(x, mapped = 1L)), "incomplete recorded GSEA eligibility")
})

test_that("one gene-set result cannot change with dictionary category fan-out", {
  x <- empty_eligibility_members(); second <- x; second$category_id <- "CAT_SECOND"
  for (field in c("NES", "padj")) {
    conflicting <- second
    conflicting[[field]][[1L]] <- if (field == "NES") 2 else .9
    f <- empty_eligibility_run(rbind(x, conflicting))
    expect_error(empty_eligibility_state(f), "conflicting duplicate GSEA member across category assignments")
  }
})

test_that("textual NaN representations cannot excuse excluded results", {
  for (field in c("NES", "padj")) local({
    selected <- field
    f <- empty_eligibility_run()
    read <- lisaR:::read_lisa_tsv
    testthat::local_mocked_bindings(read_lisa_tsv = function(path, ...) {
      tab <- read(path, ...)
      if (identical(path, f$gsea)) {
        tab[[selected]] <- as.character(tab[[selected]])
        tab[[selected]][[2L]] <- "NaN"
      }
      tab
    }, .package = "lisaR")
    expect_error(empty_eligibility_state(f), "contradictory statistics")
  })
})
