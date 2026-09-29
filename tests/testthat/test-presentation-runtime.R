test_that("pipeline forwards a selected NES variant to both scientific engines", {
  root <- tempfile("presentation-runtime-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  calls <- list()
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) { calls$single <<- list(...); NULL },
    run_LISA_contrast = function(...) { calls$contrast <<- list(...); NULL },
    .package = "lisaR")
  registry <- lisaR:::lisa_collection_registry("GOBP-C2")
  de <- data.frame(analysis_id = c("A", "B"), de_path = c("a.tsv", "b.tsv"),
    species = "Homo sapiens", stringsAsFactors = FALSE)
  status <- lisaR:::run_lisa_pipeline_single_de(de, registry, root, root,
    "members.tsv", "dictionary.tsv", "map.tsv", root, "core", "semantic",
    "png", "tsv", workers = 1L, category_nes_variants = "dispersion")
  expect_true(all(status$status == "completed"))
  expect_identical(calls$single$category_nes_variants, "dispersion")
  contrast <- data.frame(contrast_id = "AB", output_id = "profile", analysis_a = "A", analysis_b = "B")
  status <- lisaR:::run_lisa_pipeline_contrasts(de, contrast, registry, root, status,
    "semantic", "png", "tsv", workers = 1L, category_nes_variants = "direction")
  expect_identical(status$status, "completed")
  expect_identical(calls$contrast$category_nes_variants, "direction")
  expect_identical(calls$contrast$comparison_name, "profile")
})

test_that("contrast evidence requires each exact side before launching and uses portable routes", {
  root <- tempfile("contrast-presentation-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  # The unit deliberately launches TWO post scripts: build_contrast_evidence.R
  # and then, only when the evidence stage completed, build_contrast_navigation.R
  # (FULL navigation parity, 2026-09-09). Collect the calls keyed by script name
  # so each stage's arguments are asserted against its own launch instead of
  # against whichever call happened to be last.
  calls <- list()
  testthat::local_mocked_bindings(lisa_run_post_script = function(package_dir, script_name, args, ...) {
    calls[[script_name]] <<- stats::setNames(args[seq.int(2L, length(args), 2L)],
      args[seq.int(1L, length(args), 2L)])
    list(status = "completed", output_path = root, message = "")
  }, .package = "lisaR")
  row <- data.frame(contrast_id = "AB", collection = "GOBP-C2")
  index <- data.frame(contrast_id = "AB", output_id = "profile", analysis_a = "A", analysis_b = "B")
  status <- data.frame(stage = "single_de_category_evidence", analysis_id = c("A", "B"),
    collection = "GOBP-C2", status = "completed")
  run <- function(x = status) lisaR:::lisa_run_contrast_evidence_unit(row, index, x,
    root, "package", root, NULL, variants = "direction")
  result <- run()
  # Both stages are reported, in launch order, for the same contrast identity.
  expect_identical(result$stage,
    c("contrast_category_evidence", "contrast_category_navigation"))
  expect_identical(result$contrast_id, c("AB", "AB"))
  expect_identical(result$status, c("completed", "completed"))
  expect_setequal(names(calls),
    c("build_contrast_evidence.R", "build_contrast_navigation.R"))

  # Stage 1 - contrast evidence: exact sides, selected variant, portable routes.
  evidence <- calls[["build_contrast_evidence.R"]]
  expect_identical(unname(evidence["--analysis-a"]), "A")
  expect_identical(unname(evidence["--analysis-b"]), "B")
  expect_identical(unname(evidence["--variants"]), "direction")
  expect_identical(unname(evidence["--evidence-a"]), file.path(root, "report_pages", "evidence", "A", "GOBP-C2"))
  expect_identical(unname(evidence["--evidence-b"]), file.path(root, "report_pages", "evidence", "B", "GOBP-C2"))
  expect_identical(unname(evidence["--output-dir"]), file.path(root, "report_pages", "contrast_evidence", "AB_profile", "GOBP-C2"))
  expect_identical(unname(evidence["--report-href"]), "../../../contrasts.html")
  expect_identical(unname(evidence["--gene-explorer-href"]), "../../../gene_evidence/index.html")

  # Stage 2 - contrast navigation: same sides and scope, its own output root.
  navigation <- calls[["build_contrast_navigation.R"]]
  expect_identical(unname(navigation["--analysis-a"]), "A")
  expect_identical(unname(navigation["--analysis-b"]), "B")
  expect_identical(unname(navigation["--scope"]), "AB_profile")
  expect_identical(unname(navigation["--universe"]), "GOBP-C2")
  expect_identical(unname(navigation["--label"]), "AB")
  expect_identical(unname(navigation["--output-dir"]), file.path(root, "report_pages", "category_navigation", "AB_profile", "GOBP-C2"))
  # Navigation must not be pointed at the evidence output root.
  expect_false(identical(unname(navigation["--output-dir"]), unname(evidence["--output-dir"])))

  # Missing-side guards: neither script may launch when an exact side is absent.
  failed <- status; failed$status[[2L]] <- "failed"; calls <- list()
  expect_identical(run(failed)$status, "skipped_individual_evidence_not_completed")
  expect_length(calls, 0L)
  failed <- status; failed$collection[[2L]] <- "GOCC"
  expect_identical(run(failed)$status, "skipped_individual_evidence_not_completed")
  expect_length(calls, 0L)
  expect_identical(run(rbind(status, status[1L, ]))$status, "skipped_individual_evidence_not_completed")
  expect_length(calls, 0L)

  # Navigation is gated on the evidence stage: if evidence fails, only stage 1
  # is reported and build_contrast_navigation.R is never launched.
  calls <- list()
  testthat::local_mocked_bindings(lisa_run_post_script = function(package_dir, script_name, args, ...) {
    calls[[script_name]] <<- stats::setNames(args[seq.int(2L, length(args), 2L)],
      args[seq.int(1L, length(args), 2L)])
    list(status = "failed", output_path = root, message = "boom")
  }, .package = "lisaR")
  gated <- run()
  expect_identical(gated$stage, "contrast_category_evidence")
  expect_identical(gated$status, "failed")
  expect_identical(names(calls), "build_contrast_evidence.R")
})

test_that("required contrast evidence artifacts and selected NES files are accounted", {
  root <- tempfile("presentation-witnesses-")
  registry <- lisaR:::lisa_collection_registry(c("GOBP-C2", "HALLMARKS"))
  de <- data.frame(analysis_id = c("A", "B"), species = "Homo sapiens")
  contrast <- data.frame(contrast_id = "AB", output_id = "profile", analysis_a = "A", analysis_b = "B")
  expected <- lisaR:::lisa_expected_artifacts(de, contrast, registry, root,
    category_evidence = TRUE, category_nes_variants = "dispersion", plot_formats = c("png", "svg"))
  single_evidence <- expected[expected$stage == "single_de_category_evidence", ]
  expect_true(all(grepl("figures/category_members_index.tsv", single_evidence$witnesses, fixed = TRUE)))
  expect_true(all(grepl("figures/reproduce_category_member_evidence.R", single_evidence$witnesses, fixed = TRUE)))
  evidence <- expected[expected$stage == "contrast_category_evidence", ]
  expect_identical(evidence$expectation, "required")
  expect_identical(evidence$collection, "GOBP-C2")
  for (file in c("index.html", "tables/categories.tsv", "tables/leading_edges.tsv", "tables/conservation.tsv", "tables/provenance.tsv")) {
    expect_match(evidence$witnesses, file.path("report_pages", "contrast_evidence", "AB_profile", "GOBP-C2", file), fixed = TRUE)
  }
  nes <- expected$witnesses[expected$stage == "single_de_lisa" & expected$collection == "GOBP-C2"]
  expect_true(all(grepl("_dispersion.png", nes, fixed = TRUE)))
  expect_true(all(grepl("_dispersion.svg", nes, fixed = TRUE)))
  expect_false(any(grepl("_clean.png|_direction.png", nes)))
  paired_nes <- expected$witnesses[expected$stage == "category_contrasts" & expected$collection == "GOBP-C2"]
  expect_gt(length(paired_nes), 0L)
  for (subset in c("same_direction", "opposite_direction")) {
    expect_true(all(grepl(paste0("_", subset, "_dispersion.png"), paired_nes, fixed = TRUE)))
    expect_true(all(grepl(paste0("_", subset, "_source.tsv"), paired_nes, fixed = TRUE)))
    expect_true(all(grepl(paste0("_", subset, "_dispersion_settings.json"), paired_nes, fixed = TRUE)))
  }
  expect_false(any(grepl("category_nes", expected$witnesses[expected$collection == "HALLMARKS"])))
  no_image <- lisaR:::lisa_nes_expected_witnesses("scope", "A", "semantic", "clean", character())
  expect_identical(no_image, "file:scope/qc/A_category_nes_variants.tsv")
  no_evidence <- lisaR:::lisa_expected_artifacts(de, contrast, registry, root, category_evidence = FALSE)
  expect_false(any(no_evidence$stage == "contrast_category_evidence"))
})

test_that("planning declares contrast presentation after both prerequisite stages", {
  de <- data.frame(analysis_id = c("A", "B"), de_path = c("a.tsv", "b.tsv"), species = "Homo sapiens")
  contrast <- data.frame(contrast_id = "AB", analysis_a = "A", analysis_b = "B")
  plan <- lisaR:::lisa_pipeline_plan(de, contrast, collections = "GOBP-C2", category_evidence = TRUE)
  expect_lt(match("single_de_category_evidence", plan$stage), match("contrast_category_evidence", plan$stage))
  expect_lt(match("category_contrasts", plan$stage), match("contrast_category_evidence", plan$stage))
  expect_lt(match("contrast_category_evidence", plan$stage), match("root_html_report", plan$stage))
  expect_identical(plan$n_analyses[plan$stage == "contrast_category_evidence"], 0L)
})

test_that("report estimate tracks chosen variants without changing evidence population", {
  cfg <- list(pipeline = list(schema_version = "1.0.0", profile = "targeted", evidence_mode = "full_de",
    duplicate_policies = list(de_table_duplicate_policy = "error", matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")),
    collections = list("GOBP-C2"), single_de = list(
      list(analysis_id = "A", de_path = "a.tsv", species = "Homo sapiens"),
      list(analysis_id = "B", de_path = "b.tsv", species = "Homo sapiens")),
    contrasts = list(list(contrast_id = "AB", analysis_a = "A", analysis_b = "B")))
  all <- plan_lisa_outputs(cfg, category_counts = c(`GOBP-C2` = 5L))
  cfg$report <- list(category_nes_variants = "clean")
  one <- plan_lisa_outputs(cfg, category_counts = c(`GOBP-C2` = 5L))
  # Three additional views for two analyses and the three contrast subsets.
  expect_equal(all$summary$expected_figure_files - one$summary$expected_figure_files, 15)
  expect_identical(all$summary$category_evidence_categories, one$summary$category_evidence_categories)
  expect_identical(all$summary$contrast_evidence_categories, one$summary$contrast_evidence_categories)
})
