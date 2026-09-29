# Reuse existing fixture definitions without sourcing/rerunning their tests.
# This also makes this focused file independent of test-file execution order.
evidence_contract_helpers <- function() {
  helpers <- new.env(parent = environment(evidence_contract_helpers))
  requested <- list(
    "test-expected-artifacts.R" = c("lisa_expected_artifact_fixture", "lisa_materialize_expected_witnesses"),
    "test-code-identity.R" = "code_identity_fixture",
    "test-category-evidence.R" = "evidence_fixture")
  for (file in names(requested)) {
    expressions <- as.list(parse(testthat::test_path(file)))
    for (name in requested[[file]]) {
      hit <- vapply(expressions, function(x) is.call(x) &&
        identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name(name)), logical(1))
      stopifnot(sum(hit) == 1L)
      eval(expressions[[which(hit)]], envir = helpers)
    }
  }
  helpers
}

test_that("production category evidence needs completed status and real HTML witnesses", {
  helpers <- evidence_contract_helpers()
  root <- tempfile("evidence-artifact-contract-")
  fixture <- helpers$lisa_expected_artifact_fixture(root)
  fixture$expected <- lisaR:::lisa_expected_artifacts(
    de_index = data.frame(analysis_id = c("A", "B"), species = "Homo sapiens"),
    contrast_index = data.frame(contrast_id = "A_vs_B", analysis_a = "A", analysis_b = "B"),
    registry = lisaR:::lisa_collection_registry("GOBP-C2"), output_dir = root,
    category_evidence = TRUE)
  evidence_rows <- fixture$expected[fixture$expected$stage == "single_de_category_evidence", , drop = FALSE]
  expect_identical(evidence_rows$analysis_id, c("A", "B"))
  expect_identical(evidence_rows$expectation, rep("required", 2L))
  expect_identical(evidence_rows$status_source, rep("post_lisa_status", 2L))
  expect_identical(evidence_rows$witnesses,
    unname(vapply(c("A", "B"), function(analysis) paste(paste0("file:",
      file.path("report_pages", "evidence", analysis, "GOBP-C2",
        c("index.html", "figures/category_members_index.tsv", "figures/reproduce_category_member_evidence.R"))),
      collapse = "|"), character(1))))

  # Other scientific stages are existing fixture witnesses. The new stage is
  # built and rendered, not replaced by an arbitrary nonempty HTML marker.
  helpers$lisa_materialize_expected_witnesses(
    fixture$expected[fixture$expected$stage != "single_de_category_evidence", , drop = FALSE], root)
  input <- helpers$evidence_fixture()
  for (analysis in c("A", "B")) {
    evidence <- lisaR:::build_lisa_category_evidence(input$gsea, input$de,
      analysis_id = analysis, collection = "GOBP-C2", tier = "core")
    target <- file.path(root, "report_pages", "evidence", analysis, "GOBP-C2")
    rendered <- lisaR:::render_lisa_category_evidence(evidence, target,
      formats = "png", source_data = TRUE, recipes = TRUE)
    expect_identical(rendered$status, "completed")
    html <- paste(readLines(rendered$html, warn = FALSE), collapse = "\n")
    payload_json <- sub('(?s).*<script type="application/json" id="evidence-data">(.*?)</script>.*',
      "\\1", html, perl = TRUE)
    payload <- jsonlite::fromJSON(payload_json, simplifyVector = FALSE)
    expect_identical(payload$metadata$analysis_id, analysis)
    expect_identical(payload$metadata$collection, "GOBP-C2")
    downloads <- unlist(payload$downloads, use.names = FALSE)
    expect_true(length(downloads) > 0L)
    expect_true(all(file.exists(file.path(dirname(rendered$html), downloads))))
    sets <- utils::read.delim(file.path(target, "tables", "sets.tsv"), check.names = FALSE)
    expect_setequal(sets$state, c("significant", "non_significant", "non_evaluable"))
  }
  evidence_status <- data.frame(stage = "single_de_category_evidence",
    analysis_id = c("A", "B"), contrast_id = "", collection = "GOBP-C2",
    status = "completed", stringsAsFactors = FALSE)
  navigation_status <- evidence_status
  navigation_status$stage <- "single_de_category_navigation"
  gene_status <- evidence_status[1L, , drop = FALSE]
  gene_status$stage <- "single_de_gene_evidence"
  gene_status$analysis_id <- gene_status$collection <- ""
  contrast_status <- gene_status
  contrast_status$stage <- "contrast_category_evidence"
  contrast_status$contrast_id <- "A_vs_B"
  contrast_status$collection <- "GOBP-C2"
  fixture$post_lisa_status <- rbind(fixture$post_lisa_status, evidence_status, navigation_status, gene_status, contrast_status)
  new_stages <- fixture$expected[fixture$expected$stage %in% c("single_de_category_navigation", "single_de_gene_evidence"), , drop = FALSE]
  expect_identical(new_stages$contrast_id, rep("", nrow(new_stages)))
  expect_true(all(new_stages$expectation == "required"))
  reconcile <- function(post = fixture$post_lisa_status) lisaR:::lisa_reconcile_expected_artifacts(
    fixture$expected, fixture$single_status, fixture$contrast_status, post, root)
  validation <- reconcile()
  expect_true(lisaR:::lisa_assert_required_artifacts(validation))
  expect_identical(validation$validation_status[validation$stage == "single_de_category_evidence"],
    rep("validated", 2L))

  missing_status <- reconcile(fixture$post_lisa_status[
    !(fixture$post_lisa_status$stage == "single_de_category_evidence" &
      fixture$post_lisa_status$analysis_id == "A"), , drop = FALSE])
  id <- "single_de_category_evidence:A:GOBP-C2"
  expect_identical(missing_status$validation_status[missing_status$artifact_id == id], "missing_status")
  expect_error(lisaR:::lisa_assert_required_artifacts(missing_status),
    "single_de_category_evidence:A:GOBP-C2", class = "lisa_artifact_validation_error")

  path <- file.path(root, "report_pages", "evidence", "A", "GOBP-C2", "index.html")
  saved <- paste0(path, ".saved")
  expect_true(file.rename(path, saved))
  missing_witness <- reconcile()
  expect_identical(missing_witness$validation_status[missing_witness$artifact_id == id], "missing_witness")
  expect_error(lisaR:::lisa_assert_required_artifacts(missing_witness),
    "single_de_category_evidence:A:GOBP-C2", class = "lisa_artifact_validation_error")
  expect_true(file.rename(saved, path))
  expect_true(lisaR:::lisa_assert_required_artifacts(reconcile()))
})

test_that("no-static-format evidence reconciles its real HTML without member artifacts", {
  root <- tempfile("evidence-no-static-contract-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  report <- lisaR:::lisa_validate_report_config(list(report = list(
    formats = list(png = FALSE, svg = FALSE, pdf = FALSE))))
  formats <- names(report$formats)[report$formats]
  expect_length(formats, 0L)
  expected <- lisaR:::lisa_expected_artifacts(
    de_index = data.frame(analysis_id = "A", species = "Homo sapiens"),
    contrast_index = NULL,
    registry = lisaR:::lisa_collection_registry("GOBP-C2"), output_dir = root,
    category_evidence = TRUE, plot_formats = formats)
  expected <- expected[expected$stage == "single_de_category_evidence", , drop = FALSE]
  html_path <- file.path("report_pages", "evidence", "A", "GOBP-C2", "index.html")
  expect_identical(expected$expectation, "required")
  expect_identical(expected$witnesses, paste0("file:", html_path))

  input <- evidence_contract_helpers()$evidence_fixture()
  evidence <- lisaR:::build_lisa_category_evidence(input$gsea, input$de,
    analysis_id = "A", collection = "GOBP-C2", tier = "core")
  target <- file.path(root, dirname(html_path))
  rendered <- lisaR:::render_lisa_category_evidence(evidence, target,
    formats = formats, source_data = TRUE, recipes = TRUE)
  expect_identical(rendered$status, "completed")
  expect_false(dir.exists(file.path(target, "figures")))
  expect_true(file.exists(file.path(target, "tables", "sets.tsv")))
  status <- data.frame(stage = "single_de_category_evidence", analysis_id = "A",
    contrast_id = "", collection = "GOBP-C2", status = "completed")
  reconcile <- function(post = status) lisaR:::lisa_reconcile_expected_artifacts(
    expected, post_lisa_status = post, output_dir = root)
  result <- reconcile()
  expect_identical(result$validation_status, "validated")
  expect_equal(result$witnesses_expected, 1L)
  expect_true(lisaR:::lisa_assert_required_artifacts(result))
  expect_identical(reconcile(status[FALSE, ])$validation_status, "missing_status")

  expect_true(file.rename(rendered$html, paste0(rendered$html, ".saved")))
  missing <- reconcile()
  expect_identical(missing$validation_status, "missing_witness")
  expect_error(lisaR:::lisa_assert_required_artifacts(missing),
    "single_de_category_evidence:A:GOBP-C2", class = "lisa_artifact_validation_error")
})

test_that("each requested static format retains member index and recipe obligations", {
  for (formats in list("png", "pdf", "svg", c("png", "svg"))) {
    expected <- lisaR:::lisa_expected_artifacts(
      de_index = data.frame(analysis_id = "A", species = "Homo sapiens"),
      contrast_index = NULL,
      registry = lisaR:::lisa_collection_registry("GOBP-C2"), output_dir = tempfile(),
      category_evidence = TRUE, plot_formats = formats)
    row <- expected[expected$stage == "single_de_category_evidence", , drop = FALSE]
    expected_paths <- paste0("file:", file.path("report_pages", "evidence", "A", "GOBP-C2",
      c("index.html", "figures/category_members_index.tsv", "figures/reproduce_category_member_evidence.R")))
    expect_setequal(strsplit(row$witnesses, "|", fixed = TRUE)[[1L]], expected_paths)
    expect_identical(row$expectation, "required")
  }
})

test_that("evidence viewer assets are manifest-bound and mutation blocks launch", {
  helpers <- evidence_contract_helpers()
  code_root <- helpers$code_identity_fixture(include_script = FALSE,
    include_helper = FALSE, include_renderer = FALSE)
  script_dir <- attr(code_root, "lisaR.script_dir", exact = TRUE)
  asset_dir <- file.path(dirname(script_dir), "category-evidence")
  dir.create(asset_dir)
  entrypoint <- "build_category_evidence.R"
  assets <- c("viewer.html", "viewer.css", "viewer.js")
  expect_true(file.copy(system.file("scripts", entrypoint, package = "lisaR"), script_dir))
  expect_true(all(file.copy(system.file("category-evidence", assets, package = "lisaR"), asset_dir)))
  shell <- c("lisa_shell.css", "lisa_shell.js", "LISA_logo_C_compact_icon_muted_red_S.svg")
  shell_dir <- file.path(dirname(script_dir), "report_assets"); dir.create(shell_dir)
  expect_true(all(file.copy(system.file("report_assets", shell, package = "lisaR"), shell_dir)))
  identity <- lisaR:::lisa_post_script_identity(code_root, entrypoint)
  expected_paths <- c(file.path("scripts", entrypoint), file.path("category-evidence", assets),
    file.path("report_assets", shell))
  expect_setequal(identity$files$relative_path, expected_paths)
  expect_identical(identity$files$sha256,
    unname(vapply(identity$files$absolute_path, lisaR:::lisa_sha256_file, character(1))))
  expect_identical(lisaR:::lisa_verify_post_script_identity(identity)$files, identity$files)

  run_root <- lisaR:::lisa_run_root(tempfile("evidence-code-ledger-"))
  ledger_path <- file.path(run_root, "code_identity.tsv")
  ledger <- lisaR:::lisa_write_code_identity_ledger(code_root, entrypoint, ledger_path)
  persisted <- lisaR:::read_lisa_tsv(ledger_path)
  expect_setequal(persisted$relative_path, expected_paths)
  expect_identical(persisted$sha256, identity$files$sha256)
  expect_true(lisaR:::lisa_verify_code_identity_ledger(ledger, identity))
  subprocess_called <- FALSE
  testthat::local_mocked_bindings(lisa_run_subprocess = function(...) {
    subprocess_called <<- TRUE
    stop("Unexpected subprocess after an asset identity mismatch.")
  }, .package = "lisaR")

  for (asset in assets) {
    path <- file.path(asset_dir, asset)
    original <- readBin(path, "raw", n = file.info(path)$size)
    writeBin(c(original, charToRaw("\nmutated-viewer-asset\n")), path)
    expect_error(lisaR:::lisa_verify_post_script_identity(identity),
      "LISA-CODE-IDENTITY-002", class = "lisa_code_identity_error")
    fresh <- lisaR:::lisa_post_script_identity(code_root, entrypoint)
    expect_error(lisaR:::lisa_verify_code_identity_ledger(ledger, fresh),
      "LISA-CODE-LEDGER-003", class = "lisa_code_identity_error")
    expect_error(lisaR:::lisa_run_post_script(code_root, entrypoint,
      c("--project-dir", run_root, "--analysis-id", "A", "--universe", "GOBP-C2"),
      trusted_run_root = run_root, code_ledger = ledger),
      "LISA-CODE-LEDGER-003", class = "lisa_code_identity_error")
    expect_false(subprocess_called)
    writeBin(original, path)
    expect_true(lisaR:::lisa_verify_code_identity_ledger(ledger,
      lisaR:::lisa_verify_post_script_identity(identity)))
  }
})
