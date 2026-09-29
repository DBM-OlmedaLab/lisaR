category_member_fixture <- function() {
  g <- data.frame(category_id = rep("001", 6L), pathway = c("0007", "NEG", "ZERO", "WEAK", "UNTESTED", "TIE"),
    NES = c(1.2345678912345678, -2, 0, 1, NA_real_, 2), padj = c(.01, .02, .03, .8, NA_real_, .02),
    leadingEdge = c("G1", "G2", "G0", "GW", NA_character_, "G3"), stringsAsFactors = FALSE)
  dictionary <- data.frame(category_id = c("001", "EMPTY"), category_display_name = c("Numeric category", "No mapped support"))
  lisaR:::build_lisa_category_evidence(g, category_dictionary = dictionary,
    analysis_id = "A_001", collection = "GOMF", tier = "core", positive_contrast = "Treatment minus baseline",
    max_sets = 1L, max_genes = 1L)
}

test_that("member chart uses all significant unique sets despite matrix caps", {
  e <- category_member_fixture()
  s <- lisaR:::lisa_category_member_evidence_source(e, "001")
  m <- s[s$row_type == "metadata", ]; candidates <- s[s$row_type == "set", ]
  selected <- candidates[candidates$selected_for_plot == "TRUE", ]
  expect_equal(nrow(selected), 4L)
  expect_identical(selected$pathway, c("0007", "NEG", "TIE", "ZERO"))
  expect_identical(selected$plot_order, as.character(1:4))
  expect_equal(nrow(candidates), 6L)
  expect_identical(m$n_mapped_sets, "6")
  expect_identical(m$n_evaluable_sets, "5")
  expect_identical(m$n_significant_sets, "4")
  expect_identical(m$analysis_id, "A_001"); expect_identical(m$category_id, "001")
  expect_equal(as.numeric(selected$NES[1L]), e$sets$NES[e$sets$pathway == "0007"], tolerance = 0)
  expect_identical(selected$color, c("#b2182b", "#2166ac", "#b2182b", "#64748b"))
  expect_true(is.na(candidates$NES[candidates$pathway == "UNTESTED"]))
  expect_false(candidates$selected_for_plot[candidates$pathway == "UNTESTED"] == "TRUE")
  expect_false(candidates$selected_for_plot[candidates$pathway == "WEAK"] == "TRUE")
})

test_that("empty member categories retain scope and honest zero support", {
  e <- category_member_fixture()
  s <- lisaR:::lisa_category_member_evidence_source(e, "EMPTY")
  expect_equal(nrow(s), 1L)
  expect_identical(s$row_type, "metadata")
  expect_identical(s$n_significant_sets, "0")
  expect_identical(s$n_evaluable_sets, "0")
  expect_identical(s$n_mapped_sets, "0")
  expect_identical(s$category_id, "EMPTY")
  # The shared scale does not falsely magnify an empty/weak category.
  nonempty <- lisaR:::lisa_category_member_evidence_source(e, "001")
  expect_identical(s$plot_x_min, nonempty$plot_x_min[[1L]])
  expect_identical(s$plot_x_max, nonempty$plot_x_max[[1L]])
  expect_error(lisaR:::lisa_category_member_evidence_source(e, "unknown"), "Unknown exact category")
})

test_that("member source checks scope and recorded significance without changing it", {
  e <- category_member_fixture(); before <- e
  lisaR:::lisa_category_member_evidence_source(e, "001")
  expect_identical(e, before)
  e$sets$significant[e$sets$pathway == "WEAK"] <- TRUE
  expect_error(lisaR:::lisa_category_member_evidence_source(e, "001"), "significance conflicts")
  e <- before; e$sets <- rbind(e$sets, e$sets[1L, ])
  expect_error(lisaR:::lisa_category_member_evidence_source(e, "001"), "Duplicate gene sets")
  e <- before; e$categories$category_id[[1L]] <- "OTHER_UNCLASSIFIED"
  expect_error(lisaR:::render_lisa_category_member_evidence(e, tempfile()), "classified")
})

test_that("member figure index and standalone source preserve full identifiers", {
  e <- category_member_fixture(); out <- tempfile("member-chart-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  index <- lisaR:::render_lisa_category_member_evidence(e, out, formats = character())
  expect_identical(index$category_id, e$categories$category_id)
  expect_true(all(file.exists(file.path(out, index$source_tsv))))
  expect_true(all(file.exists(file.path(out, index$recipe_r))))
  expect_true(all(nchar(index$source_sha256) == 64L))
  expect_true(all(nchar(index$recipe_sha256) == 64L))
  expect_true(file.exists(file.path(out, "figures", "category_members_index.tsv")))
  source <- utils::read.delim(file.path(out, index$source_tsv[[1L]]), sep = "\t", quote = "",
    comment.char = "", check.names = FALSE, colClasses = "character", na.strings = NULL)
  expect_true("0007" %in% source$pathway)
  expect_identical(unique(source$category_id), "001")
  expect_identical(digest::digest(file = file.path(out, index$source_tsv[[1L]]), algo = "sha256"), index$source_sha256[[1L]])
})

test_that("standalone member recipe reproduces the figure pixels", {
  skip_if_not_installed("png")
  e <- category_member_fixture(); out <- tempfile("member-recipe-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  index <- lisaR:::render_lisa_category_member_evidence(e, out, formats = "png")
  for (i in seq_len(nrow(index))) {
    original <- file.path(out, "figures", paste0(index$stem[[i]], ".png"))
    rendered <- file.path(out, paste0(index$stem[[i]], "_reproduced.png"))
    result <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(file.path(out, index$recipe_r[[i]])),
      shQuote(file.path(out, index$source_tsv[[i]])), shQuote(rendered)), stdout = TRUE, stderr = TRUE)
    expect_null(attr(result, "status"))
    expect_true(file.exists(rendered))
    expect_identical(png::readPNG(rendered), png::readPNG(original))
  }
})

test_that("standard category HTML carries the full member chart contract", {
  e <- category_member_fixture(); out <- tempfile("member-integrated-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  rendered <- lisaR:::render_lisa_category_evidence(e, out, formats = "png")
  expect_identical(rendered$member_figure_index$category_id, e$categories$category_id)
  html <- paste(readLines(rendered$html, warn = FALSE), collapse = "\n")
  expect_match(html, '"member_figure_index":', fixed = TRUE)
  expect_match(html, 'id="member-chart-panel"', fixed = TRUE)
  expect_match(html, 'Gene sets in this category', fixed = TRUE)
  index <- rendered$member_figure_index
  expect_true(all(file.exists(file.path(out, index$source_tsv))))
  expect_true(all(file.exists(file.path(out, index$recipe_r))))
  expect_equal(index$n_significant_sets[index$category_id == "001"], 4L)
  expect_identical(e$metadata$max_sets, 1L)
})

test_that("member and matrix exports obey independent source and recipe flags", {
  e <- category_member_fixture()
  for (sources in c(TRUE, FALSE)) for (recipes in c(TRUE, FALSE)) {
    out <- tempfile("member-export-policy-")
    rendered <- lisaR:::render_lisa_category_evidence(e, out, formats = "png",
      source_data = sources, recipes = recipes)
    html <- paste(readLines(rendered$html, warn = FALSE), collapse = "\n")
    json <- sub('.*<script type="application/json" id="evidence-data">', '', html)
    json <- sub('</script>.*', '', json)
    payload <- jsonlite::fromJSON(json)
    expect_identical(payload$source_data, sources)
    expect_identical(payload$recipes, recipes)
    expect_true(all(nzchar(payload$member_figure_index$source_tsv) == sources))
    expect_true(all(nzchar(payload$member_figure_index$recipe_r) == recipes))
    if (recipes) {
      expect_true(file.exists(file.path(out, 'reproduce_category_evidence.R')))
      expect_true(file.exists(file.path(out, 'figures/001_category_evidence_source.tsv')))
    }
    # Canonical scientific tables and the displayed figure are not removed.
    expect_true(file.exists(file.path(out, 'tables/sets.tsv')))
    expect_true(file.exists(file.path(out, 'figures/001_category_members.png')))
    unlink(out, recursive = TRUE)
  }
})
