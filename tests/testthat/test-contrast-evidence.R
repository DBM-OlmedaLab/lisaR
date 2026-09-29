contrast_evidence_fixture <- function() {
  a <- data.frame(category_id = rep("001", 5), pathway = c("0001", "neg", "zero", "weak", "missing"),
    NES = c(2, -1, 0, 1, NA_real_), padj = c(.01, .03, .04, .8, NA_real_),
    leadingEdge = c("001;G2", "G2", "G0", "GW", NA_character_), stringsAsFactors = FALSE)
  b <- data.frame(category_id = rep("001", 5), pathway = c("0001", "neg", "zero", "b_only", "weak"),
    NES = c(1, -2, 0, 2, NA_real_), padj = c(.02, .5, .04, .02, NA_real_),
    leadingEdge = c("001;G3", "G2", "G0", "G4", NA_character_), stringsAsFactors = FALSE)
  da <- data.frame(symbol = c("001", "G2", "G0", "GW", "UNRELATED"), log2FC = c(.5, -1, .2, 1, .9), padj = c(.01, .2, .8, .3, .001))
  db <- data.frame(symbol = c("001", "G3", "G0", "G4", "UNRELATED"), log2FC = c(-2, 1, .1, 2, -.3), padj = c(.02, .1, .9, .004, .8))
  dictionary <- data.frame(category_id = c("001", "EMPTY"), category_display_name = c("Numeric ID category", "No support"))
  list(a = lisaR:::build_lisa_category_evidence(a, da, category_dictionary = dictionary,
      analysis_id = "analysis A", collection = "GOMF", tier = "core", positive_contrast = "treatment vs baseline"),
    b = lisaR:::build_lisa_category_evidence(b, db, category_dictionary = dictionary,
      analysis_id = "analysis B", collection = "GOMF", tier = "core", positive_contrast = "treatment vs baseline"))
}

test_that("paired evidence uses independent full-support distributions, not visible means", {
  x <- contrast_evidence_fixture()
  e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", max_sets = 1L, max_genes = 1L)
  c <- e$categories[e$categories$category_id == "001", ]
  expect_equal(c$mean_NES_a, 1/3); expect_equal(c$median_NES_a, 0)
  expect_equal(c$p25_NES_a, -.5); expect_equal(c$p75_NES_a, 1)
  expect_equal(c$mean_NES_b, 1); expect_equal(c$delta_mean_NES, -2/3)
  expect_equal(c$n_genesets_significant_a, 3); expect_equal(c$n_genesets_significant_b, 3)
  expect_equal(c$positive_pct_a, 100/3); expect_equal(c$negative_pct_a, 100/3); expect_equal(c$zero_pct_a, 100/3)
  expect_equal(c$same_direction_pct_a, 100/3)
  expect_equal(c$n_significant_both, 2); expect_equal(c$n_significant_A_only, 1); expect_equal(c$n_significant_B_only, 1)
  expect_equal(c$n_mapped_a, 5); expect_equal(c$n_mapped_b, 5)
  empty <- e$categories[e$categories$category_id == "EMPTY", ]
  expect_equal(empty$n_genesets_significant_a, 0)
  expect_true(is.na(empty$mean_NES_a)); expect_true(is.na(empty$positive_pct_a)); expect_true(is.na(empty$delta_mean_NES))
  expect_false(any(grepl("delta.*(padj|fdr)|combined.*fdr|delta.*(p25|p75)", names(e$categories))))
})

test_that("canonical contextual endpoints preserve the existing descriptive delta", {
  x <- contrast_evidence_fixture()
  x$b$sets$gsea_fdr[x$b$sets$evaluable] <- .9
  x$b$sets$significant[] <- FALSE
  x$b$categories$mean_significant_set_NES[] <- NA_real_
  canonical <- data.frame(category_id = "001", mean_NES_A = 1/3, mean_NES_B = NA_real_,
    display_mean_NES_A = 1/3, display_mean_NES_B = .25, delta_mean_NES = 1/3 - .25,
    universe = "GOMF", gsea_padj_cutoff = .25)
  e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", contrast_summary = canonical)
  c <- e$categories[e$categories$category_id == "001", ]
  expect_equal(c$display_mean_NES_b, .25)
  expect_equal(c$delta_mean_NES, canonical$delta_mean_NES)
  expect_equal(c$endpoint_source_b, "contextual_mean_not_significant")
  expect_true(is.na(c$mean_NES_b)); expect_true(is.na(c$median_NES_b))
  expect_true(is.na(c$p25_NES_b)); expect_true(is.na(c$positive_pct_b))
  expect_equal(c$n_genesets_significant_b, 0)
  reversed <- canonical; reversed$delta_mean_NES <- -reversed$delta_mean_NES
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", contrast_summary = reversed), "A minus B")
  wrong <- canonical; wrong$display_mean_NES_B <- 0
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", contrast_summary = wrong), "contextual support")
})

test_that("categories outside the canonical contrast do not acquire new deltas", {
  x <- contrast_evidence_fixture()
  for (side in c("a", "b")) {
    extra_category <- x[[side]]$categories[1L, , drop = FALSE]
    extra_category$category_id <- "outside"
    extra_sets <- x[[side]]$sets
    extra_sets$category_id <- "outside"
    x[[side]]$categories <- rbind(x[[side]]$categories, extra_category)
    x[[side]]$sets <- rbind(x[[side]]$sets, extra_sets)
  }
  canonical <- data.frame(category_id = "001", mean_NES_A = 1/3, mean_NES_B = 1,
    display_mean_NES_A = 1/3, display_mean_NES_B = 1, delta_mean_NES = -2/3,
    universe = "GOMF", gsea_padj_cutoff = .25)
  e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", contrast_summary = canonical)
  inside <- e$categories[e$categories$category_id == "001", ]
  outside <- e$categories[e$categories$category_id == "outside", ]
  expect_identical(inside$delta_state, "canonical")
  expect_equal(inside$delta_mean_NES, canonical$delta_mean_NES)
  expect_equal(outside$mean_NES_a, 1/3)
  expect_equal(outside$mean_NES_b, 1)
  expect_identical(outside$delta_state, "not_in_canonical_contrast")
  expect_true(is.na(outside$delta_mean_NES))
})

test_that("union intersection and evaluable selection align rows without losing full results", {
  x <- contrast_evidence_fixture(); e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01")
  union <- lisaR:::lisa_contrast_evidence_view(e, "001", "union")
  intersection <- lisaR:::lisa_contrast_evidence_view(e, "001", "intersection")
  all <- lisaR:::lisa_contrast_evidence_view(e, "001", "evaluable")
  expect_setequal(union$sets$pathway, c("0001", "neg", "zero", "b_only"))
  expect_setequal(intersection$sets$pathway, c("0001", "zero"))
  expect_setequal(all$sets$pathway, c("0001", "neg", "zero", "weak", "b_only"))
  expect_equal(nrow(e$sets), 6L)
  missing <- e$sets[e$sets$pathway == "missing", ]
  expect_false(missing$present_b); expect_false(missing$evaluable_b)
  expect_true(is.na(missing$NES_b)); expect_true(is.na(missing$gsea_fdr_b))
  weak <- e$sets[e$sets$pathway == "weak", ]
  expect_true(weak$present_b); expect_false(weak$evaluable_b)
  expect_true("UNRELATED" %in% e$genes$symbol)
  expect_equal(e$genes$log2FC_a[e$genes$symbol == "001"], .5)
  expect_equal(e$genes$log2FC_b[e$genes$symbol == "001"], -2)
})

test_that("same-set LE conservation distinguishes absent from disjoint evidence", {
  x <- contrast_evidence_fixture(); e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01")
  first <- e$conservation[e$conservation$pathway == "0001", ]
  expect_equal(first$n_a, 2L); expect_equal(first$n_b, 2L)
  expect_equal(first$intersection_n, 1L); expect_equal(first$union_n, 3L); expect_equal(first$jaccard, 1/3)
  expect_identical(first$overlap_state, "defined")
  no <- e$conservation[e$conservation$pathway == "weak", ]
  expect_identical(no$overlap_state, "not_available"); expect_true(is.na(no$jaccard)); expect_true(is.na(no$intersection_n))
  expect_equal(no$n_a, 1L); expect_true(is.na(no$n_b))
  x$b$leading_edges$symbol[x$b$leading_edges$pathway == "0001"] <- c("different1", "different2")
  disjoint <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01")$conservation
  expect_equal(disjoint$jaccard[disjoint$pathway == "0001"], 0)
  expect_equal(disjoint$overlap_state[disjoint$pathway == "0001"], "defined")
})

test_that("scope or significance conflicts fail instead of relabelling a contrast", {
  x <- contrast_evidence_fixture()
  wrong <- x$b; wrong$metadata$tier <- "expanded"
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, wrong, "C01"), "matching tier")
  wrong <- x$b; wrong$metadata$collection <- "GOCC"
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, wrong, "C01"), "matching collection")
  wrong <- x$b; wrong$metadata$analysis_id <- "analysis A"
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, wrong, "C01"), "different analyses")
  wrong <- x$b; wrong$sets$significant[[1L]] <- FALSE
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, wrong, "C01"), "significance conflicts")
  wrong <- x$b; wrong$categories$mean_significant_set_NES[[1L]] <- 42
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, wrong, "C01"), "Recorded category mean")
  wrong <- x$b; wrong$sets <- rbind(wrong$sets, wrong$sets[1L, ])
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, wrong, "C01"), "Duplicate exact")
  expect_error(lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", max_sets = 0), "positive integers")
})

test_that("paired source retains numeric-looking identifiers and identical row/column order", {
  x <- contrast_evidence_fixture(); e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01", max_sets = 2L, max_genes = 1L)
  s <- lisaR:::lisa_contrast_evidence_figure_source(e, "001")
  a <- s[s$row_type == "set" & s$side == "a", ]; b <- s[s$row_type == "set" & s$side == "b", ]
  expect_identical(a$pathway, b$pathway); expect_identical(a$set_order, b$set_order)
  expect_true("0001" %in% a$pathway); expect_true(all(s$category_id == "001"))
  ga <- s[s$row_type == "gene" & s$side == "a", ]; gb <- s[s$row_type == "gene" & s$side == "b", ]
  expect_identical(ga$symbol, gb$symbol); expect_identical(ga$gene_order, gb$gene_order)
  expect_equal(sum(s$row_type == "set"), 4L)
  expect_equal(unique(s$n_matching_sets), "4")
  reverse <- x; reverse$a$sets <- reverse$a$sets[nrow(reverse$a$sets):1, ]
  reverse$b$sets <- reverse$b$sets[nrow(reverse$b$sets):1, ]
  e2 <- lisaR:::build_lisa_contrast_evidence(reverse$a, reverse$b, "C01", max_sets = 2L, max_genes = 1L)
  expect_identical(e$sets, e2$sets); expect_identical(s, lisaR:::lisa_contrast_evidence_figure_source(e2, "001"))
})

test_that("contrast reader restores normalized tables and binds exact provenance", {
  x <- contrast_evidence_fixture(); root <- tempfile("normalized-evidence-")
  dir.create(file.path(root, "tables"), recursive = TRUE); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (key in c("categories", "sets", "genes", "leading_edges")) lisaR:::write_lisa_tsv(x$a[[key]], file.path(root, "tables", paste0(key, ".tsv")))
  metadata <- data.frame(key = names(x$a$metadata), value = vapply(x$a$metadata, as.character, character(1L)))
  lisaR:::write_lisa_tsv(metadata, file.path(root, "tables", "metadata.tsv"))
  e <- lisaR:::build_lisa_contrast_evidence(root, x$b, "C01")
  expect_equal(nrow(e$provenance), 5L); expect_true(all(nchar(e$provenance$sha256) == 64L))
  expect_true("001" %in% e$genes$symbol); expect_true("0001" %in% e$sets$pathway)
  expect_false(any(grepl(root, e$provenance$table, fixed = TRUE)))
  expect_false(any(grepl(root, rownames(e$provenance), fixed = TRUE)))
  out <- file.path(root, "portable")
  rendered <- lisaR:::render_lisa_contrast_evidence(e, out, formats = character())
  html <- paste(readLines(rendered$html, warn = FALSE), collapse = "\n")
  expect_false(grepl(root, html, fixed = TRUE))
  expect_false(grepl('"_row":', html, fixed = TRUE))
  expect_true(all(vapply(e$provenance$sha256, grepl, logical(1L), x = html, fixed = TRUE)))
})

test_that("offline paired render writes full tables and source-based recipe", {
  skip_if_not(nzchar(system.file("contrast-evidence", package = "lisaR")))
  x <- contrast_evidence_fixture(); e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01")
  out <- tempfile("contrast-evidence-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  r <- lisaR:::render_lisa_contrast_evidence(e, out, formats = "png", variants = "direction")
  expect_true(file.exists(r$html)); expect_true(file.exists(file.path(out, "figures", "001_contrast_evidence.png")))
  expect_true(file.exists(file.path(out, "figures", "001_contrast_evidence_source.tsv")))
  expect_true(file.exists(file.path(out, "reproduce_contrast_evidence.R")))
  html <- paste(readLines(r$html, warn = FALSE), collapse = "\n")
  expect_match(html, 'id="contrast-evidence-main"', fixed = TRUE)
  expect_match(html, 'id="contrast-evidence-data"', fixed = TRUE)
  expect_match(html, 'data-lisa-shell', fixed = TRUE)
  expect_match(html, 'href="#contrast-evidence-main"', fixed = TRUE)
  expect_length(gregexpr('<main\\b', html, perl = TRUE)[[1L]], 1L)
  expect_false(grepl('class="evidence-header"', html, fixed = TRUE))
  expect_false(grepl('aria-current="page"', html, fixed = TRUE)) # honest standalone
  expect_true(file.exists(file.path(out, "lisa-shell", "lisa_shell.js")))
  expect_false(grepl("fetch(", html, fixed = TRUE)); expect_false(grepl("https://", html, fixed = TRUE))
  expect_match(html, 'id="category"', fixed = TRUE)
  expect_match(html, 'id="category-products"', fixed = TRUE)
  expect_false(grepl('id="navigator-section"', html, fixed = TRUE))
  expect_false(grepl('id="sheet" hidden', html, fixed = TRUE))
  sets <- utils::read.delim(file.path(out, "tables", "sets.tsv"), colClasses = "character")
  expect_equal(nrow(sets), nrow(e$sets)); expect_true("0001" %in% sets$pathway)
  expect_error(lisaR:::render_lisa_contrast_evidence(e, out, variants = "unknown"), "Unknown")
})

test_that("contrast category products are exact-ID scoped and remain optional", {
  x <- contrast_evidence_fixture(); e <- lisaR:::build_lisa_contrast_evidence(x$a, x$b, "C01")
  products <- list(list(category_id = "001", product = "gene_cards", label = "Contrast GeneCards",
    assets = list(list(format = "png", href = "saved/category_001_cards.png"))))
  normalized <- lisaR:::lisa_contrast_evidence_category_products(products, e$categories$category_id)
  expect_identical(normalized[[1L]]$category_id, "001")
  expect_identical(normalized[[1L]]$assets[[1L]]$href, "saved/category_001_cards.png")
  canonical <- products
  canonical[[1L]]$assets[[1L]]$href <-
    "../../../../artifacts/contrasts/C01/GO/contrast_heatmap/001.png"
  expect_no_error(lisaR:::lisa_contrast_evidence_category_products(
    canonical, e$categories$category_id
  ))
  bad_id <- products; bad_id[[1L]]$category_id <- "OTHER_UNCLASSIFIED"
  expect_error(lisaR:::lisa_contrast_evidence_category_products(bad_id, e$categories$category_id), "exact contrast category")
  bad_path <- products; bad_path[[1L]]$assets[[1L]]$href <- "../outside.png"
  expect_error(lisaR:::lisa_contrast_evidence_category_products(bad_path, e$categories$category_id), "safe relative file")
})

test_that("static paired renderer uses the same A-circle B-square semantics", {
  # Inspect the renderer actually loaded by this test. An installed package
  # has no source-tree R/contrast_evidence.R; its namespace is authoritative
  # both under load_all() and under R CMD check of the installed tarball.
  renderer <- paste(deparse(body(lisaR:::lisa_contrast_evidence_draw),
                            width.cutoff = 500L), collapse = "\n")
  renderer <- gsub("[[:space:]]+", " ", renderer)
  expect_match(renderer, 'pch = if (side == "a") 16 else 15', fixed = TRUE)
  expect_match(renderer, "A circle; B square", fixed = TRUE)
})

test_that("contrast gene links use query parameters and preserve both exact side scopes", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; browser link regression needs its URL API")
  script <- system.file("contrast-evidence", "viewer.js", package = "lisaR")
  skip_if(!nzchar(script))
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  # Execute the actual link builder with punctuation and numeric-looking IDs,
  # then read it through the same URLSearchParams API as the receiving page.
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "const source=fs.readFileSync(process.argv[2],'utf8');",
    "const start=source.indexOf('  function geneExplorerHref(side) {');",
    "const end=source.indexOf('  function geneDetails()',start);",
    "assert.ok(start>=0 && end>start);",
    "const M={analysis_a:'A + 001',analysis_b:'B & 002',collection:'GOBP-C2',tier:'core',contrast_id:'C: A/B',gene_explorer_href:'../../../gene_evidence/index.html'};",
    "const selectedGene='001+A&B/GENE',active='0007';",
    "const location=new URL('file:///evaluation/report_pages/contrast_evidence/C_output/GOBP-C2/index.html');",
    "const hrefs=vm.runInNewContext(source.slice(start,end)+';[geneExplorerHref(\"a\"),geneExplorerHref(\"b\")]',{M,selectedGene,active,URLSearchParams,location});",
    "hrefs.forEach((href,i)=>{",
    "  assert.ok(href.startsWith(M.gene_explorer_href+'?'));",
    "  const url=new URL(href,location);",
    "  assert.equal(url.hash,'');",
    "  assert.deepEqual(Object.fromEntries(new URLSearchParams(url.search)),{gene:selectedGene,analysis_id:M[i===0?'analysis_a':'analysis_b'],collection:M.collection,tier:M.tier,return_contrast_id:M.contrast_id,return_category_id:active,return_contrast_scope:'C_output',return_analysis_a:M.analysis_a,return_analysis_b:M.analysis_b,return_side:i===0?'a':'b'});",
    "});",
    "assert.notEqual(hrefs[0],hrefs[1]);",
    "process.stdout.write('PASS exact query-scoped A/B gene links\\n');"
  ), probe, useBytes = TRUE)
  result <- suppressWarnings(system2(node, c(shQuote(probe), shQuote(script)), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS exact query-scoped A/B gene links", fixed = TRUE)
})
