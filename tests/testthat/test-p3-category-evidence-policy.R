# Installed report regression: exercise the real shared figure renderer for
# both single-analysis and contrast callsites, without rerunning enrichment.
test_that("category evidence policy controls links, never exact saved bitmaps", {
  root <- withr::local_tempdir("p3-evidence-policy-")
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts/build_LISA_report.R", package = "lisaR"), env)
  env$report_project_dir <- root
  env$report_package_dir <- find.package("lisaR")
  page <- file.path(root, "report_pages", "categories.html")
  dir.create(dirname(page))
  file.create(page)
  for (owner in c("single_de", "category_contrasts")) {
    collection <- file.path(root, "outputs", owner, "A", "collection_GOBP-C2")
    dir.create(file.path(collection, "qc"), recursive = TRUE)
    image <- file.path(collection, "saved.png")
    grDevices::png(image, width = 320, height = 240)
    graphics::plot(c(1, 2), c(2, 1), main = "Saved category figure")
    grDevices::dev.off()
    source <- file.path(collection, "source.tsv")
    writeLines("category_id\tmean_NES\nA\t2\nB\t1", source)
    settings <- file.path(collection, "settings.json")
    jsonlite::write_json(list(variant = "clean", width = 320, height = 240,
      source_sha256 = digest::digest(file = source, algo = "sha256")),
      settings, auto_unbox = TRUE)
    writeLines("# saved figure recipe", file.path(collection, "recipe.R"))
    utils::write.table(data.frame(variant = "clean", status = "rendered",
      format = "png", path = "saved.png", source = "source.tsv",
      settings = "settings.json", recipe = "recipe.R"),
      file.path(collection, "qc", "saved_category_nes_variants.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
    evidence <- file.path(root, paste0(owner, "-missing-evidence.html"))
    # Match the real run: export-policy TSV has no category-evidence row.
    policy_file <- file.path(root, "report_output_policy.tsv")
    utils::write.table(data.frame(product = c("png", "source_data", "recipes"),
      requested = c(TRUE, FALSE, FALSE), status = c("requested", "not_requested", "not_requested")),
      policy_file, sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(data.frame(key = "report_category_evidence", value = "FALSE"),
      file.path(root, "contract_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    policy_sha <- digest::digest(file = policy_file, algo = "sha256")
    env$report_output_policy <- env$report_load_output_policy(root)
    expect_false(env$report_requested("category_evidence", TRUE))
    expect_identical(digest::digest(file = policy_file, algo = "sha256"), policy_sha)
    html <- env$report_lisa_categories_section("Saved categories", collection,
      image, page, evidence_index = evidence)
    expect_match(html, '<img class="nes-variant-image"', fixed = TRUE)
    expect_false(grepl("category-figure-hit|missing-evidence[.]html|category_id=", html))
    # Resolve the real relative media URL and verify the copied bitmap bytes.
    href <- sub('.*<img class="nes-variant-image"[^>]*src="([^"]+)".*', '\\1', html)
    expect_identical(digest::digest(file = file.path(dirname(page), href), algo = "sha256"),
      digest::digest(file = image, algo = "sha256"))
    # Explicit TRUE and the default TRUE must still reject required missing
    # evidence; the fix must not merely ignore file.exists() for every run.
    for (explicit in c(TRUE, FALSE)) {
      env$report_output_policy <- data.frame(
        product = c(if (explicit) "category_evidence", "source_data", "recipes"),
        requested = c(if (explicit) TRUE, FALSE, FALSE))
      expect_error(env$report_lisa_categories_section("Saved categories", collection,
        image, page, evidence_index = evidence), "Category evidence target is missing", fixed = TRUE)
    }
  }
})
