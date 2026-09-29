trusted_rds_fixture <- function(object) {
  path <- tempfile(fileext = ".rds")
  saveRDS(object, path)
  path
}

trusted_edger_fixture <- function(class_name, as_matrix = FALSE) {
  table <- data.frame(
    logFC = c(1.5, -0.75),
    logCPM = c(5, 4),
    PValue = c(0.01, 0.2),
    FDR = c(0.02, 0.2),
    row.names = c("GENE1", "GENE2")
  )
  if (isTRUE(as_matrix)) table <- as.matrix(table)
  structure(list(table = table), class = class_name)
}

write_trusted_rds_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  path
}

test_that("trusted RDS parameters are strict", {
  fixture <- trusted_rds_fixture(data.frame(symbol = "G1", score = 1))
  expect_error(
    lisaR:::lisa_prepare_de_input(fixture, "auto", trusted_rds = "true"),
    "LISA-RDS-001"
  )
  for (value in list(0, NA_real_, Inf, 1.5, "1024")) {
    expect_error(
      lisaR:::lisa_prepare_de_input(
        fixture, "auto", trusted_rds = TRUE, rds_max_bytes = value
      ),
      "LISA-RDS-002"
    )
  }
})

test_that("public run_lisa_de cannot inject internal RDS provenance controls", {
  object <- data.frame(symbol = "G1", score = 1)
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) stop("engine must not be reached"),
    .package = "lisaR"
  )
  private_names <- c(
    "input_receipt", ".receipt_object_verified",
    "input_rece", ".receipt_object_ver"
  )
  for (private_name in private_names) {
    dots <- stats::setNames(list(if (startsWith(private_name, "input")) list() else TRUE),
                            private_name)
    error <- tryCatch(
      do.call(
        run_lisa_de,
        c(
          list(
            input = object, output_dir = tempfile("lisa-rds-private-"),
            trusted_rds = TRUE
          ),
          dots
        )
      ),
      error = identity
    )
    expect_s3_class(error, "lisa_rds_error")
    expect_s3_class(error, "lisa_error")
    expect_identical(error$code, "LISA-RDS-016")
    expect_match(conditionMessage(error), private_name, fixed = TRUE)
  }
})

test_that("public run_lisa_de rejects unnamed advanced arguments", {
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) stop("engine must not be reached"),
    .package = "lisaR"
  )
  historical_positional <- list(
    "table.tsv", "unused", NULL, "Homo sapiens",
    NULL, NULL, NULL, NULL, NULL,
    c("GOBP-C2", "GOMF", "GOCC"), "core", NULL, TRUE
  )
  error <- tryCatch(
    do.call(run_lisa_de, c(historical_positional, list("de_table"))),
    error = identity
  )
  expect_s3_class(error, "lisa_rds_error")
  expect_s3_class(error, "lisa_error")
  expect_identical(error$code, "LISA-RDS-016")
  expect_match(conditionMessage(error), "must be named explicitly", fixed = TRUE)
})

test_that("public RDS control inspection does not force advanced arguments", {
  forced <- FALSE
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) invisible(TRUE),
    .package = "lisaR"
  )
  expect_true(run_lisa_de(
    input = data.frame(symbol = "G1", score = 1),
    output_dir = tempfile("lisa-rds-lazy-"),
    delayed_advanced_argument = {
      forced <<- TRUE
      42
    }
  ))
  expect_false(forced)
})

test_that("failure to fingerprint an accepted RDS is a typed lisa error", {
  fixture <- trusted_rds_fixture(data.frame(symbol = "G1", score = 1))
  testthat::local_mocked_bindings(
    lisa_rds_object_sha256 = function(object) stop("simulated digest failure"),
    .package = "lisaR"
  )
  error <- tryCatch(
    lisaR:::lisa_prepare_de_input(
      fixture, "auto", trusted_rds = TRUE, rds_max_bytes = 1024^2
    ),
    error = identity
  )
  expect_s3_class(error, "lisa_rds_error")
  expect_s3_class(error, "lisa_error")
  expect_identical(error$code, "LISA-RDS-015")
})

test_that("untrusted, oversized and non-regular RDS paths are rejected before reading", {
  fixture <- trusted_rds_fixture(data.frame(symbol = "G1", score = 1))
  reads <- 0L
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      readRDS(path)
    },
    .package = "lisaR"
  )
  untrusted_error <- tryCatch(
    lisaR:::lisa_prepare_de_input(fixture, "auto", trusted_rds = FALSE),
    error = identity
  )
  expect_s3_class(untrusted_error, "lisa_rds_error")
  expect_s3_class(untrusted_error, "lisa_error")
  expect_match(conditionMessage(untrusted_error), "LISA-RDS-003")
  expect_identical(reads, 0L)

  expect_error(
    lisaR:::lisa_prepare_de_input(
      fixture, "auto", trusted_rds = TRUE,
      rds_max_bytes = file.info(fixture)$size - 1
    ),
    "LISA-RDS-005"
  )
  expect_identical(reads, 0L)

  directory <- tempfile(fileext = ".rds")
  dir.create(directory)
  expect_error(
    lisaR:::lisa_prepare_de_input(directory, "auto", trusted_rds = TRUE),
    "LISA-RDS-004"
  )
  expect_identical(reads, 0L)
})

test_that("trusted RDS file names with control characters fail before reading", {
  reads <- 0L
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      readRDS(path)
    },
    .package = "lisaR"
  )

  for (control_character in c("\t", "\n")) {
    fixture <- file.path(tempdir(), paste0("trusted", control_character, "input.rds"))
    expect_error(
      lisaR:::lisa_prepare_de_input(
        fixture, "auto", trusted_rds = TRUE, rds_max_bytes = 1024^2
      ),
      "LISA-RDS-014"
    )
  }
  expect_identical(reads, 0L)
})

test_that("an authorized data frame is read exactly once and type checked", {
  object <- data.frame(
    symbol = c("G1", "G2"), score = c(2, -1),
    stringsAsFactors = FALSE
  )
  fixture <- trusted_rds_fixture(object)
  reads <- 0L
  object_digests <- 0L
  real_object_digest <- lisaR:::lisa_rds_object_sha256
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      readRDS(path)
    },
    .package = "lisaR"
  )
  testthat::local_mocked_bindings(
    lisa_rds_object_sha256 = function(object) {
      object_digests <<- object_digests + 1L
      real_object_digest(object)
    },
    .package = "lisaR"
  )

  prepared <- lisaR:::lisa_prepare_de_input(
    fixture, "auto", trusted_rds = TRUE, rds_max_bytes = 1024^2
  )
  expect_identical(reads, 1L)
  expect_identical(object_digests, 1L)
  expect_identical(prepared$input, object)
  expect_identical(prepared$input_type, "de_table")
  expect_identical(prepared$receipt$schema_version, "lisa_rds_input_receipt_v2")
  expect_identical(prepared$receipt$source_file, basename(fixture))
  expect_identical(prepared$receipt$source_sha256, lisaR:::lisa_sha256_file(fixture))
  expect_identical(
    prepared$receipt$object_sha256,
    lisaR:::lisa_rds_object_sha256(object)
  )
  expect_identical(prepared$receipt$input_type, "de_table")
  expect_identical(prepared$receipt$object_class, "data.frame")
  expect_identical(prepared$receipt$n_rows, 2L)
  expect_identical(prepared$receipt$n_columns, 2L)
  expect_identical(prepared$receipt$serialized_bytes, as.numeric(file.info(fixture)$size))
  expect_identical(prepared$receipt$memory_bytes, as.numeric(utils::object.size(object)))
  expect_identical(prepared$receipt$rds_max_bytes, 1024^2)
  expect_true(prepared$receipt$trusted_rds_authorized)
  expect_identical(prepared$receipt$authorization_scope, "explicit_local_R_call")
  expect_identical(prepared$receipt$path_policy, "basename_only")
  expect_identical(prepared$receipt$trust_assumption, "local_stable_fully_trusted_file")

  expect_error(
    lisaR:::lisa_prepare_de_input(
      fixture, "edger", trusted_rds = TRUE, rds_max_bytes = 1024^2
    ),
    "LISA-RDS-007"
  )
  expect_identical(reads, 2L)
})

test_that("only differential-expression edgeR result classes cross the RDS boundary", {
  for (class_name in c("DGELRT", "DGEExact", "TopTags")) {
    object <- trusted_edger_fixture(class_name)
    fixture <- trusted_rds_fixture(object)
    prepared <- lisaR:::lisa_prepare_de_input(
      fixture, "auto", trusted_rds = TRUE, rds_max_bytes = 1024^2
    )
    expect_identical(prepared$input_type, "edger")
    expect_identical(prepared$receipt$object_class, class_name)
    expect_identical(c(prepared$receipt$n_rows, prepared$receipt$n_columns), c(2L, 4L))

    extracted <- lisaR:::extract_de_table(
      prepared$input, prepared$input_type, NULL, NULL
    )
    standardized <- lisaR:::standardize_de_table(
      extracted,
      symbol_col = "rowname", rank_col = NULL, logfc_col = "logFC",
      padj_col = "FDR", pvalue_col = "PValue", gene_id_col = NULL
    )
    expect_identical(standardized$de$symbol, c("GENE1", "GENE2"))
  }

  matrix_object <- trusted_edger_fixture("DGELRT", as_matrix = TRUE)
  matrix_prepared <- lisaR:::lisa_prepare_de_input(
    trusted_rds_fixture(matrix_object), "auto",
    trusted_rds = TRUE, rds_max_bytes = 1024^2
  )
  expect_identical(matrix_prepared$input_type, "edger")

  counts <- structure(
    list(counts = matrix(1:4, nrow = 2)),
    class = "DGEList"
  )
  expect_error(
    lisaR:::lisa_prepare_de_input(
      trusted_rds_fixture(counts), "auto",
      trusted_rds = TRUE, rds_max_bytes = 1024^2
    ),
    "LISA-RDS-007.*contains counts"
  )
})

test_that("spoofed or malformed edgeR result classes fail closed", {
  missing_table <- structure(list(), class = "DGELRT")
  missing_pvalue <- trusted_edger_fixture("DGEExact")
  missing_pvalue$table$PValue <- NULL
  character_logfc <- trusted_edger_fixture("TopTags")
  character_logfc$table$logFC <- as.character(character_logfc$table$logFC)

  for (object in list(missing_table, missing_pvalue, character_logfc)) {
    expect_error(
      lisaR:::lisa_prepare_de_input(
        trusted_rds_fixture(object), "auto",
        trusted_rds = TRUE, rds_max_bytes = 1024^2
      ),
      "LISA-RDS-012"
    )
  }
  expect_error(
    lisaR:::extract_de_table(list(table = data.frame()), "edger", NULL, NULL),
    "LISA-RDS-012"
  )
})

test_that("real edgeR DE results are accepted and DGEList is rejected", {
  skip_if_not_installed("edgeR")
  counts <- matrix(
    c(
      10, 12, 40, 42,
      20, 18, 55, 58,
      30, 28, 70, 75,
      12, 11, 33, 36
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("GENE", 1:4), paste0("S", 1:4))
  )
  group <- factor(c("A", "A", "B", "B"))
  dge <- edgeR::DGEList(counts = counts, group = group)
  exact <- edgeR::exactTest(dge, dispersion = 0.1)
  design <- stats::model.matrix(~ group)
  fit <- edgeR::glmFit(dge, design, dispersion = 0.1)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  top <- edgeR::topTags(exact, n = Inf)

  for (object in list(exact, lrt, top)) {
    prepared <- lisaR:::lisa_prepare_de_input(
      trusted_rds_fixture(object), "auto",
      trusted_rds = TRUE, rds_max_bytes = 1024^2
    )
    expect_identical(prepared$input_type, "edger")
    extracted <- lisaR:::extract_de_table(prepared$input, "edger", NULL, NULL)
    standardized <- lisaR:::standardize_de_table(
      extracted,
      symbol_col = "rowname", rank_col = NULL, logfc_col = "logFC",
      padj_col = if ("FDR" %in% names(extracted)) "FDR" else NULL,
      pvalue_col = "PValue", gene_id_col = NULL
    )
    expect_true(nrow(standardized$de) > 0L)
    expect_true(all(standardized$de$symbol %in% rownames(extracted)))
  }

  expect_error(
    lisaR:::lisa_prepare_de_input(
      trusted_rds_fixture(dge), "auto",
      trusted_rds = TRUE, rds_max_bytes = 1024^2
    ),
    "LISA-RDS-007.*contains counts"
  )
})

test_that("real DESeqResults row names can be selected explicitly as symbols", {
  skip_if_not_installed("DESeq2")
  de_result <- DESeq2::DESeqResults(data.frame(
    baseMean = c(25, 40, 15, 60),
    log2FoldChange = c(1.5, -0.8, 0.4, -1.2),
    lfcSE = c(0.2, 0.3, 0.25, 0.4),
    stat = c(7.5, -2.67, 1.6, -3),
    pvalue = c(0.001, 0.02, 0.1, 0.01),
    padj = c(0.004, 0.04, 0.13, 0.03),
    row.names = paste0("GENE", 1:4)
  ))
  prepared <- lisaR:::lisa_prepare_de_input(
    trusted_rds_fixture(de_result), "auto",
    trusted_rds = TRUE, rds_max_bytes = 1024^2
  )
  expect_identical(prepared$input_type, "deseq2_results")
  extracted <- lisaR:::extract_de_table(
    prepared$input, prepared$input_type, NULL, NULL
  )
  standardized <- lisaR:::standardize_de_table(
    extracted,
    symbol_col = "rowname", rank_col = "stat",
    logfc_col = "log2FoldChange", padj_col = "padj",
    pvalue_col = "pvalue", gene_id_col = NULL
  )
  expect_true(nrow(standardized$de) > 0L)
  expect_true(all(standardized$de$symbol %in% rownames(extracted)))
})

test_that("unsupported classes, empty dimensions and expanded objects fail after one read", {
  unsupported <- trusted_rds_fixture(list(symbol = "G1", score = 1))
  empty <- trusted_rds_fixture(data.frame())
  expanded_object <- data.frame(
    symbol = rep(sprintf("GENE_%06d", seq_len(20000)), each = 1L),
    score = seq_len(20000), stringsAsFactors = FALSE
  )
  expanded <- trusted_rds_fixture(expanded_object)
  expanded_file_size <- as.numeric(file.info(expanded)$size)
  expanded_memory_size <- as.numeric(utils::object.size(expanded_object))
  expect_lt(expanded_file_size, expanded_memory_size)
  expanded_limit <- floor((expanded_file_size + expanded_memory_size) / 2)

  reads <- 0L
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      readRDS(path)
    },
    .package = "lisaR"
  )
  expect_error(
    lisaR:::lisa_prepare_de_input(
      unsupported, "auto", trusted_rds = TRUE, rds_max_bytes = 1024^2
    ),
    "LISA-RDS-007"
  )
  expect_identical(reads, 1L)
  expect_error(
    lisaR:::lisa_prepare_de_input(
      empty, "auto", trusted_rds = TRUE, rds_max_bytes = 1024^2
    ),
    "LISA-RDS-008"
  )
  expect_identical(reads, 2L)
  expect_error(
    lisaR:::lisa_prepare_de_input(
      expanded, "auto", trusted_rds = TRUE, rds_max_bytes = expanded_limit
    ),
    "LISA-RDS-009"
  )
  expect_identical(reads, 3L)
})

test_that("objects already loaded in memory retain the existing path", {
  object <- data.frame(symbol = c("G1", "G2"), score = c(1, -1))
  reads <- 0L
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      stop("must not be called")
    },
    .package = "lisaR"
  )
  prepared <- lisaR:::lisa_prepare_de_input(
    object, "auto", trusted_rds = FALSE, rds_max_bytes = 1024^2
  )
  expect_identical(prepared$input, object)
  expect_identical(prepared$input_type, "auto")
  expect_identical(lisaR:::resolve_input_type(prepared$input, prepared$input_type),
                   "de_table")
  expect_identical(reads, 0L)
})

test_that("a trusted-RDS receipt survives internal branch preparation", {
  object <- data.frame(symbol = c("G1", "G2"), score = c(1, -1))
  prepared <- lisaR:::lisa_prepare_de_input(
    trusted_rds_fixture(object), "auto",
    trusted_rds = TRUE, rds_max_bytes = 1024^2
  )
  carried <- lisaR:::lisa_prepare_de_input(
    prepared$input, prepared$input_type,
    trusted_rds = TRUE, rds_max_bytes = 1024^2,
    input_receipt = prepared$receipt
  )
  expect_identical(carried$source, "trusted_rds_propagated")
  expect_identical(carried$receipt, prepared$receipt)

  malformed_receipts <- list(
    within(prepared$receipt, source_sha256 <- ""),
    within(prepared$receipt, object_sha256 <- ""),
    within(prepared$receipt, source_file <- "trusted\tinput.rds"),
    within(prepared$receipt, source_file <- "trusted\ninput.rds"),
    within(prepared$receipt, object_class <- "matrix"),
    within(prepared$receipt, n_rows <- n_rows + 1L),
    within(prepared$receipt, n_rows <- n_rows + 0.5),
    within(prepared$receipt, n_columns <- n_columns + 0.5),
    within(prepared$receipt, serialized_bytes <- serialized_bytes + 0.5),
    within(prepared$receipt, memory_bytes <- memory_bytes + 1),
    within(prepared$receipt, memory_bytes <- memory_bytes + 0.5),
    within(prepared$receipt, input_type <- "edger"),
    within(prepared$receipt, rds_max_bytes <- rds_max_bytes + 1),
    within(prepared$receipt, rds_max_bytes <- rds_max_bytes + 0.5)
  )
  for (malformed in malformed_receipts) {
    expect_error(
      lisaR:::lisa_prepare_de_input(
        prepared$input, prepared$input_type,
        trusted_rds = TRUE, rds_max_bytes = 1024^2,
        input_receipt = malformed
      ),
      "LISA-RDS-013"
    )
  }
  expect_error(
    lisaR:::lisa_prepare_de_input(
      rbind(prepared$input, data.frame(symbol = "G3", score = 0)),
      prepared$input_type, trusted_rds = TRUE, rds_max_bytes = 1024^2,
      input_receipt = prepared$receipt
    ),
    "LISA-RDS-013"
  )
  same_shape_mutation <- prepared$input
  same_shape_mutation$score[[1L]] <- 999
  expect_identical(
    dim(same_shape_mutation), dim(prepared$input)
  )
  expect_identical(
    as.numeric(utils::object.size(same_shape_mutation)),
    as.numeric(utils::object.size(prepared$input))
  )
  expect_error(
    lisaR:::lisa_prepare_de_input(
      same_shape_mutation, prepared$input_type,
      trusted_rds = TRUE, rds_max_bytes = 1024^2,
      input_receipt = prepared$receipt
    ),
    "LISA-RDS-013"
  )
  expect_error(
    lisaR:::lisa_prepare_de_input(
      prepared$input, prepared$input_type,
      trusted_rds = FALSE, rds_max_bytes = 1024^2,
      input_receipt = prepared$receipt
    ),
    "LISA-RDS-013"
  )
})

test_that("one RDS read and one portable receipt are preserved in every branch", {
  de <- data.frame(
    symbol = c("GENE1", "GENE2", "GENE3", "GENE4"),
    logFC = c(2, -1, 0.5, -0.25),
    PValue = c(0.001, 0.01, 0.2, 0.5),
    stringsAsFactors = FALSE
  )
  input_path <- trusted_rds_fixture(de)
  fixture_dir <- tempfile("lisa-rds-branches-")
  dir.create(fixture_dir)
  dictionary_path <- write_trusted_rds_tsv(
    data.frame(
      universe = c("GOBP-C2", "GOMF"),
      gene_set_id = c("GOBP_TEST_SET", "GOMF_TEST_SET"),
      gene_set_name = c("Test BP set", "Test MF set"),
      source_id = c("GO:TEST:BP", "GO:TEST:MF"),
      category_id = c("TEST_BP", "TEST_MF"),
      category_display_name = c("Test BP", "Test MF"),
      tier = c("fixture", "fixture"),
      stringsAsFactors = FALSE
    ),
    file.path(fixture_dir, "dictionary.tsv")
  )
  category_map_path <- write_trusted_rds_tsv(
    data.frame(
      category_id = c("TEST_BP", "TEST_MF"),
      display_name = c("Test BP", "Test MF"),
      macrogroup_id = c("TEST", "TEST"),
      macrogroup_name = c("Test", "Test"),
      macrogroup_order = c(1, 1),
      category_order_within_macrogroup = c(1, 2),
      notes = c("fixture", "fixture"),
      stringsAsFactors = FALSE
    ),
    file.path(fixture_dir, "category_map.tsv")
  )
  term2gene_path <- write_trusted_rds_tsv(
    data.frame(
      gs_collection = rep(c("GOBP-C2", "GOMF"), each = 2),
      gs_subcollection = "FIXTURE",
      gs_name = rep(c("GOBP_TEST_SET", "GOMF_TEST_SET"), each = 2),
      gs_exact_source = "FIXTURE",
      gene_symbol = rep(c("GENE1", "GENE2"), 2),
      stringsAsFactors = FALSE
    ),
    file.path(fixture_dir, "term2gene.tsv")
  )

  reads <- 0L
  object_digests <- 0L
  real_object_digest <- lisaR:::lisa_rds_object_sha256
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      readRDS(path)
    },
    .package = "lisaR"
  )
  testthat::local_mocked_bindings(
    lisa_rds_object_sha256 = function(object) {
      object_digests <<- object_digests + 1L
      real_object_digest(object)
    },
    .package = "lisaR"
  )
  collection_output <- file.path(fixture_dir, "collection-output")
  collection_results <- lisaR:::run_LISA_DE(
    input = input_path, output_dir = collection_output,
    comparison_name = "branch_fixture", symbol_col = "symbol",
    logfc_col = "logFC", pvalue_col = "PValue",
    lisa_dictionary = "fixture", lisa_dictionary_path = dictionary_path,
    term2gene_path = term2gene_path, category_map_path = category_map_path,
    universes = c("GOBP-C2", "GOMF"), run_gsea = FALSE, run_ora = FALSE,
    outputs = c("tables", "qc"), plots = character(),
    export_formats = "tsv", trusted_rds = TRUE, rds_max_bytes = 1024^2,
    verbose = FALSE
  )
  expect_identical(reads, 1L)
  expect_identical(object_digests, 1L)
  expect_setequal(names(collection_results), c("GOBP-C2", "GOMF"))
  collection_receipts <- lapply(collection_results, function(result) {
    expect_true(file.exists(result$input_receipt_path))
    lisaR:::read_lisa_tsv(result$input_receipt_path)
  })
  expect_identical(collection_receipts[[1L]], collection_receipts[[2L]])
  expect_identical(
    collection_receipts[[1L]]$source_sha256,
    lisaR:::lisa_sha256_file(input_path)
  )
  expect_identical(collection_receipts[[1L]]$source_file, basename(input_path))
  receipt_text <- readLines(collection_results[[1L]]$input_receipt_path, warn = FALSE)
  expect_false(any(grepl(dirname(normalizePath(input_path)), receipt_text, fixed = TRUE)))

  resource_cache <- file.path(fixture_dir, "registered-resources")
  core_path <- write_trusted_rds_tsv(
    transform(
      lisaR:::read_lisa_tsv(dictionary_path)[1, , drop = FALSE],
      gene_set_id = "GOBP_CORE_SET", gene_set_name = "Core fixture set",
      tier = "core"
    ),
    file.path(
      resource_cache, "lisa_dictionary_core", "1.0.0", "dictionary.tsv"
    )
  )
  expanded_path <- write_trusted_rds_tsv(
    transform(
      lisaR:::read_lisa_tsv(dictionary_path)[1, , drop = FALSE],
      gene_set_id = "GOBP_EXPANDED_SET",
      gene_set_name = "Expanded fixture set",
      tier = "expanded"
    ),
    file.path(
      resource_cache, "lisa_dictionary_expanded", "1.0.0", "dictionary.tsv"
    )
  )
  core <- lisaR:::read_lisa_tsv(core_path)
  expanded <- lisaR:::read_lisa_tsv(expanded_path)
  core_set <- core$gene_set_id[core$universe == "GOBP-C2"][[1L]]
  expanded_set <- expanded$gene_set_id[expanded$universe == "GOBP-C2"][[1L]]
  dictionary_term2gene <- write_trusted_rds_tsv(
    data.frame(
      gs_collection = "GOBP-C2", gs_subcollection = "FIXTURE",
      gs_name = rep(unique(c(core_set, expanded_set)), each = 2),
      gs_exact_source = "FIXTURE",
      gene_symbol = rep(c("GENE1", "GENE2"), length(unique(c(core_set, expanded_set)))),
      stringsAsFactors = FALSE
    ),
    file.path(
      resource_cache, "msigdb_term2gene", "2026.1", "term2gene.tsv"
    )
  )
  registered_category_map <- write_trusted_rds_tsv(
    lisaR:::read_lisa_tsv(category_map_path),
    file.path(
      resource_cache, "lisa_category_map", "1.0.0", "category_map.tsv"
    )
  )
  resource_paths <- c(
    core_path, expanded_path, dictionary_term2gene, registered_category_map
  )
  resource_registry <- data.frame(
    # The fixture stands in for the registered scientific defaults, so it must
    # declare the identities the pipeline actually selects.
    logical_id = c(
      "lisa_dictionary_core", "lisa_dictionary_expanded", "msigdb_term2gene",
      "lisa_category_map"
    ),
    version = c("1.0.0", "1.0.0", "2026.1", "1.0.0"),
    species = "Homo sapiens",
    modality = "transcriptomic/genomic",
    schema = c(
      "lisa_dictionary@2", "lisa_dictionary@2", "term2gene@1",
      "category_map@1"
    ),
    sha256 = vapply(resource_paths, lisaR:::lisa_sha256_file, character(1)),
    compatibility = "lisaR>=0.6.0",
    approved_origin = "trusted-RDS test fixture",
    artifact = c(
      "dictionary.tsv", "dictionary.tsv", "term2gene.tsv",
      "category_map.tsv"
    ),
    stringsAsFactors = FALSE
  )
  withr::local_options(list(
    lisaR.dictionary_registry = resource_registry,
    lisaR.dictionary_cache_root = resource_cache,
    lisaR.shared_dictionary_root = NULL
  ))
  reads <- 0L
  object_digests <- 0L
  dictionary_output <- file.path(fixture_dir, "dictionary-output")
  dictionary_results <- lisaR:::run_LISA_DE(
    input = input_path, output_dir = dictionary_output,
    comparison_name = "dictionary_fixture", symbol_col = "symbol",
    logfc_col = "logFC", pvalue_col = "PValue",
    lisa_dictionary = c("core", "expanded"),
    universes = "GOBP-C2",
    run_gsea = FALSE, run_ora = FALSE,
    outputs = c("tables", "qc"), plots = character(),
    export_formats = "tsv", trusted_rds = TRUE, rds_max_bytes = 1024^2,
    verbose = FALSE
  )
  expect_identical(reads, 1L)
  expect_identical(object_digests, 1L)
  expect_setequal(names(dictionary_results), c("core", "expanded"))
  dictionary_receipts <- lapply(dictionary_results, function(result) {
    expect_true(file.exists(result$input_receipt_path))
    lisaR:::read_lisa_tsv(result$input_receipt_path)
  })
  expect_identical(dictionary_receipts[[1L]], dictionary_receipts[[2L]])
  expect_identical(dictionary_receipts[[1L]], collection_receipts[[1L]])
})

test_that("configuration indexes cannot grant trusted-RDS capability", {
  index <- data.frame(
    analysis_id = "A", de_path = "input.rds", trusted_rds = TRUE,
    stringsAsFactors = FALSE
  )
  reads <- 0L
  testthat::local_mocked_bindings(
    lisa_read_rds_once = function(path) {
      reads <<- reads + 1L
      stop("must not be called")
    },
    .package = "lisaR"
  )
  expect_error(
    validate_lisa_de_index(index, require_files = FALSE),
    "LISA-RDS-011"
  )
  expect_identical(reads, 0L)
})

test_that("YAML and JSON study objects cannot grant trusted-RDS capability", {
  config <- list(
    pipeline = list(
      schema_version = "1.0.0",
      profile = "transcriptomic/genomic",
      evidence_mode = "full_de",
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    single_de = list(list(
      analysis_id = "A", de_path = "input.rds", species = "Homo sapiens",
      trusted_rds = TRUE
    ))
  )

  expect_false(
    "trusted_rds" %in% lisaR:::lisa_config_allowed_keys("analysis")
  )
  expect_error(
    lisaR:::lisa_validate_raw_config(config),
    "LISA-CONFIG-UNKNOWN-001.*trusted_rds"
  )

  config$single_de[[1L]]$trusted_rds <- NULL
  config$single_de[[1L]]$input_type <- "deseq2_dds"
  expect_false(
    "input_type" %in% lisaR:::lisa_config_allowed_keys("analysis")
  )
  expect_error(
    lisaR:::lisa_validate_raw_config(config),
    "LISA-CONFIG-UNKNOWN-001.*input_type"
  )
})

test_that("run_lisa_de exposes and forwards the local RDS capability", {
  captured <- NULL
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) {
      captured <<- list(...)
      invisible(captured)
    },
    .package = "lisaR"
  )
  run_lisa_de(
    input = "trusted.rds", output_dir = "unused", trusted_rds = TRUE,
    rds_max_bytes = 4096
  )
  expect_true(captured$trusted_rds)
  expect_identical(captured$rds_max_bytes, 4096)

  formal_names <- names(formals(run_lisa_de))
  expect_lt(match("make_heatmaps", formal_names), match("...", formal_names))
  expect_gt(match("trusted_rds", formal_names), match("...", formal_names))
  expect_gt(match("rds_max_bytes", formal_names), match("...", formal_names))

  positional <- list(
    "table.tsv", "unused", NULL, "Homo sapiens",
    NULL, NULL, NULL, NULL, NULL,
    c("GOBP-C2", "GOMF", "GOCC"), "core", NULL, TRUE
  )
  captured <- NULL
  do.call(run_lisa_de, positional)
  expect_true(captured$make_heatmaps)
  expect_false(captured$trusted_rds)
  expect_identical(captured$rds_max_bytes, 512 * 1024^2)
})
