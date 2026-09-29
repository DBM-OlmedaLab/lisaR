contract_file <- function(path, contents) {
  writeLines(contents, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

contract_input_row <- function(path, role = "analysis.de_path", owner = "a") {
  lisaR:::lisa_contract_file_row(path, role, owner)[,
    c("role", "owner", "basename", "bytes", "sha256"), drop = FALSE
  ]
}

contract_resource_rows <- function(paths) {
  roles <- c(
    "dictionary_resource", "term2gene_resource", "category_map_resource"
  )
  do.call(rbind, lapply(seq_along(paths), function(i) {
    data.frame(
      role = roles[[i]],
      resource_id = paste0("fixture_", i, "@1.0.0"),
      basename = basename(paths[[i]]),
      species = "Homo sapiens",
      modality = "transcriptomic/genomic",
      schema = c("lisa_dictionary@2", "term2gene@1", "category_map@1")[[i]],
      bytes = as.numeric(file.info(paths[[i]])$size),
      sha256 = lisaR:::lisa_sha256_file(paths[[i]]),
      stringsAsFactors = FALSE
    )
  }))
}

contract_code_row <- function(path) {
  data.frame(
    role = "r_code", portable_id = file.path("R", basename(path)),
    basename = basename(path), bytes = as.numeric(file.info(path)$size),
    sha256 = lisaR:::lisa_sha256_file(path), stringsAsFactors = FALSE
  )
}

contract_environment <- function(version = "1.0.0") {
  data.frame(
    component = c("runtime", "package"), name = c("R", "fixture"),
    version = c("R fixture", version), platform = c("fixture-platform", ""),
    requirement = c("R runtime", "active dependency"),
    stringsAsFactors = FALSE
  )
}

test_that("legacy run-contract callers retain their exact public shape", {
  root <- tempfile("lisa-legacy-contract-")
  dir.create(root)
  config <- contract_file(file.path(root, "study.yml"), "pipeline: {}")
  dictionary <- contract_file(file.path(root, "dictionary.tsv"), "term\tcategory")

  contract <- lisaR:::lisa_run_contract(
    config, dictionary, environment = "R fixture"
  )

  expect_identical(
    names(contract),
    c("config", "code", "dictionary", "environment", "contract")
  )
  expect_true(lisaR:::lisa_sha256_all_valid(unname(contract), 5L))
  expect_null(attr(contract, "lisaR.contract_receipts", exact = TRUE))
})

test_that("text and transported contract hashes use frozen canonical bytes", {
  expect_identical(
    lisaR:::lisa_sha256_text(character()),
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  )
  expect_identical(
    lisaR:::lisa_sha256_text("abc"),
    "edeaaff3f1774ad2888673770c6d64097e391bc362d7d6fb34982ddf0efd18cb"
  )
  expect_identical(
    lisaR:::lisa_sha256_text(c("alpha", "\u03b2eta")),
    "239a0e52ce9e1a69c573e6d22538de140c1de7878e8252dd918092e7303ca1dc"
  )
  expect_identical(
    lisaR:::lisa_sha256_text(c("a\nb", "c")),
    "880553fca8fcea94e325ee2cfb48e5a985cc797f39a14cc6d3cedecfeb2ae4d2"
  )

  transported <- c(
    config = "5579d3262c6e66c324eda5ab33024ea8490b388f7800ff5a5eb0e46be5895102",
    inputs = "ca374d41d57ba1799df3fbb679a52a0c757473563cf1cbb6a12453c09df214f3",
    resources = "08b637f743bff66df2cdbdf7eb23bf2de0a62a4aa4c8eda1a461c891b1997935",
    code = "cae41a19cc73ed83128f137ed367457399c0e0847000c0713a1f2bbfef66c3cd",
    environment = "a5c69bfa9841e4e5f516be60741411d409440271357cb5fb83c614c4ca1306c5",
    contract = "87296526c6bcc82ec34bffc74f3fbe1f6ada5c55e93e5e219863f8178bde5f96"
  )
  expect_silent(lisaR:::lisa_validate_transaction_contract(transported))

  root <- tempfile("lisa-transported-contract-")
  dir.create(root)
  legacy_windows <- transported
  legacy_windows[["contract"]] <-
    "7e634e9fb7fde332d1f48eda35ca1d7f2838f4d492b65966463c2ec3d0011901"
  lisaR:::write_lisa_tsv(
    data.frame(key = names(legacy_windows), value = unname(legacy_windows),
               stringsAsFactors = FALSE),
    file.path(root, "run_contract.tsv")
  )
  expect_error(
    lisaR:::lisa_validate_transaction_contract(legacy_windows),
    "aggregate contract digest"
  )
  expect_identical(
    unname(lisaR:::lisa_read_transaction_contract(root)),
    unname(legacy_windows)
  )

  mutated <- legacy_windows
  mutated[["inputs"]] <- paste(rep("0", 64L), collapse = "")
  lisaR:::write_lisa_tsv(
    data.frame(key = names(mutated), value = unname(mutated),
               stringsAsFactors = FALSE),
    file.path(root, "run_contract.tsv")
  )
  expect_error(
    lisaR:::lisa_read_transaction_contract(root),
    "aggregate contract digest"
  )
})

test_that("incomplete or malformed contracts fail before transaction writes", {
  root <- tempfile("lisa-invalid-transaction-contract-")
  dir.create(root)
  final <- file.path(root, "new-parent", "result")

  incomplete <- lisaR:::lisa_run_contract()
  error <- tryCatch(
    lisaR:::lisa_transaction_begin(
      final, incomplete, run_id = "incompletecontract"
    ),
    lisa_contract_error = identity
  )
  expect_s3_class(error, "lisa_contract_error")
  expect_s3_class(error, "lisa_error")
  expect_identical(error$code, "LISA-CONTRACT-010")
  expect_false(dir.exists(dirname(final)))
  expect_false(lisaR:::lisa_path_entry_exists(paste0(final, ".lisa.lock")))

  valid <- lisa_test_run_contract()
  invalid_values <- list(
    blank = "",
    missing = NA_character_,
    malformed = "not-a-sha256",
    uppercase = toupper(valid[["config"]])
  )
  for (value in invalid_values) {
    candidate <- valid
    candidate[["config"]] <- value
    expect_error(
      lisaR:::lisa_transaction_begin(
        file.path(root, paste0("result-", length(value))), candidate,
        run_id = "invaliddigest"
      ),
      class = "lisa_contract_error"
    )
  }

  wrong_aggregate <- valid
  wrong_aggregate[["contract"]] <- paste(rep("0", 64L), collapse = "")
  expect_error(
    lisaR:::lisa_transaction_begin(
      file.path(root, "wrong-aggregate"), wrong_aggregate,
      run_id = "wrongaggregate"
    ),
    "aggregate contract digest"
  )

  wrong_shape <- valid[c("config", "dictionary", "code", "environment", "contract")]
  expect_error(
    lisaR:::lisa_transaction_begin(
      file.path(root, "wrong-shape"), wrong_shape,
      run_id = "wrongshape"
    ),
    "complete canonical shape"
  )
})

test_that("promotion revalidates its in-memory and persisted contracts", {
  root <- tempfile("lisa-promote-contract-revalidation-")
  dir.create(root)

  final_memory <- file.path(root, "memory-result")
  tx_memory <- lisaR:::lisa_transaction_begin(
    final_memory, lisa_test_run_contract(), run_id = "memorycontract"
  )
  on.exit(lisaR:::lisa_transaction_abort(tx_memory, "test cleanup"),
          add = TRUE)
  writeLines("artifact", file.path(tx_memory$staging_dir, "artifact.txt"))
  tx_memory$contract[["code"]] <- ""
  expect_error(
    lisaR:::lisa_transaction_promote(tx_memory),
    class = "lisa_contract_error"
  )
  expect_false(dir.exists(final_memory))

  final_receipt <- file.path(root, "receipt-result")
  tx_receipt <- lisaR:::lisa_transaction_begin(
    final_receipt, lisa_test_run_contract(), run_id = "receiptcontract"
  )
  on.exit(lisaR:::lisa_transaction_abort(tx_receipt, "test cleanup"),
          add = TRUE)
  writeLines("artifact", file.path(tx_receipt$staging_dir, "artifact.txt"))
  replacement <- lisa_test_run_contract(environment = "different fixture")
  lisaR:::write_lisa_tsv(
    data.frame(
      key = names(replacement), value = unname(replacement),
      stringsAsFactors = FALSE
    ),
    file.path(tx_receipt$staging_dir, "run_contract.tsv")
  )
  expect_error(
    lisaR:::lisa_transaction_promote(tx_receipt),
    "no longer matches"
  )
  expect_false(dir.exists(final_receipt))
})

test_that("verification rejects an invalid contract even with a rebuilt manifest", {
  final <- tempfile("lisa-verify-invalid-contract-")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "verifyinvalidcontract"
  )
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  lisaR:::lisa_transaction_promote(tx)

  receipt_path <- file.path(final, "run_contract.tsv")
  receipt <- lisaR:::read_lisa_tsv(receipt_path)
  receipt$value[receipt$key == "config"] <- ""
  lisaR:::write_lisa_tsv(receipt, receipt_path)
  rebuilt <- lisaR:::lisa_run_manifest(final)
  lisaR:::write_lisa_tsv(rebuilt, file.path(final, "run_manifest.tsv"))

  verification <- lisaR:::verify_run(final)
  expect_identical(verification$gate, "FAIL")
  expect_true("run_contract_invalid" %in% verification$findings)
})

test_that("canonical contract is portable, ordered, and sensitive to every bound input", {
  root <- tempfile("lisa-canonical-contract-")
  dir.create(root)
  config <- contract_file(file.path(root, "study.yml"), "pipeline: baseline")
  input <- contract_file(file.path(root, "de.tsv"), "symbol\tvalue\nA\t1")
  resources <- vapply(seq_len(3L), function(i) {
    contract_file(
      file.path(root, paste0("resource-", i, ".tsv")),
      paste0("resource\t", i)
    )
  }, character(1))
  code <- contract_file(file.path(root, "engine.R"), "engine <- 1")

  build <- function(environment_version = "1.0.0") {
    lisaR:::lisa_run_contract(
      config_path = config,
      input_inventory = contract_input_row(input),
      resource_inventory = contract_resource_rows(resources),
      code_inventory = contract_code_row(code),
      effective_environment = contract_environment(environment_version)
    )
  }
  baseline <- build()
  expect_identical(
    names(baseline),
    c("config", "inputs", "resources", "code", "environment", "contract")
  )
  expect_true(lisaR:::lisa_sha256_all_valid(unname(baseline), 6L))

  receipts <- attr(baseline, "lisaR.contract_receipts", exact = TRUE)
  expect_identical(
    names(receipts),
    c(
      "input_inventory.tsv", "resource_inventory.tsv",
      "executable_code_inventory.tsv", "effective_environment.tsv"
    )
  )
  expect_true(all(vapply(receipts, is.data.frame, logical(1))))
  portable_text <- paste(unlist(receipts, use.names = FALSE), collapse = "\n")
  expect_false(grepl(normalizePath(root, winslash = "/"), portable_text,
                     fixed = TRUE))

  reversed_resources <- contract_resource_rows(resources)[3:1, , drop = FALSE]
  reordered <- lisaR:::lisa_run_contract(
    config_path = config,
    input_inventory = contract_input_row(input),
    resource_inventory = reversed_resources,
    code_inventory = contract_code_row(code),
    effective_environment = contract_environment()
  )
  expect_identical(unname(reordered), unname(baseline))

  contract_file(config, "pipeline: changed")
  expect_false(identical(build()[["contract"]], baseline[["contract"]]))
  contract_file(config, "pipeline: baseline")

  contract_file(input, "symbol\tvalue\nA\t2")
  expect_false(identical(build()[["contract"]], baseline[["contract"]]))
  contract_file(input, "symbol\tvalue\nA\t1")

  contract_file(resources[[2L]], "resource\tchanged")
  expect_false(identical(build()[["contract"]], baseline[["contract"]]))
  contract_file(resources[[2L]], "resource\t2")

  changed_modality <- contract_resource_rows(resources)
  changed_modality$modality[[2L]] <- "all"
  modality_contract <- lisaR:::lisa_run_contract(
    config_path = config,
    input_inventory = contract_input_row(input),
    resource_inventory = changed_modality,
    code_inventory = contract_code_row(code),
    effective_environment = contract_environment()
  )
  expect_false(identical(
    modality_contract[["resources"]], baseline[["resources"]]
  ))

  contract_file(code, "engine <- 2")
  expect_false(identical(build()[["contract"]], baseline[["contract"]]))
  contract_file(code, "engine <- 1")

  expect_false(identical(
    build("1.0.1")[["contract"]], baseline[["contract"]]
  ))
})

test_that("effective environment records the recursive installed dependency closure", {
  db <- rbind(
    c(Package = "lisaR", Version = "0.6.0", Depends = "R", Imports = "A", LinkingTo = NA),
    c(Package = "A", Version = "1.0.0", Depends = NA, Imports = "B", LinkingTo = NA),
    c(Package = "B", Version = "2.0.0", Depends = NA, Imports = NA, LinkingTo = NA),
    c(Package = "C", Version = "3.0.0", Depends = NA, Imports = "D", LinkingTo = NA),
    c(Package = "D", Version = "4.0.0", Depends = NA, Imports = NA, LinkingTo = NA)
  )
  dependency_status <- data.frame(
    package = "C", reason = "configured capability", stringsAsFactors = FALSE
  )

  observed <- lisaR:::lisa_effective_environment(
    dependency_status,
    installed_db = db,
    runtime_version = "R test",
    runtime_platform = "test-platform"
  )

  expect_identical(observed$name, c("R", "A", "B", "C", "D", "lisaR"))
  expect_identical(
    observed$version,
    c("R test", "1.0.0", "2.0.0", "3.0.0", "4.0.0", "0.6.0")
  )
  expect_identical(
    observed$requirement[observed$name == "C"], "configured capability"
  )
  expect_true(all(
    observed$requirement[observed$name %in% c("A", "B", "D")] ==
      "dependency closure"
  ))
})

test_that("effective environment excludes transitive build-only dependencies", {
  db <- rbind(
    c(Package = "lisaR", Version = "0.6.0", Depends = "R", Imports = "A", LinkingTo = NA),
    c(Package = "A", Version = "1.0.0", Depends = NA, Imports = "B", LinkingTo = "Toolchain"),
    c(Package = "B", Version = "2.0.0", Depends = NA, Imports = NA, LinkingTo = NA)
  )
  db <- db[, setdiff(colnames(db), "LinkingTo"), drop = FALSE]

  observed <- lisaR:::lisa_effective_environment(
    installed_db = db,
    runtime_version = "R test",
    runtime_platform = "test-platform"
  )

  expect_identical(observed$name, c("R", "A", "B", "lisaR"))
  expect_false("Toolchain" %in% observed$name)
})

test_that("scientific inventory covers external indices and every original input role", {
  root <- tempfile("lisa-scientific-inventory-")
  dir.create(root)
  files <- c(
    de_index = "de-index.tsv", contrast_index = "contrast-index.tsv",
    de = "de.tsv", matrix = "matrix.tsv", mapped = "mapped.tsv",
    source = "source.tsv", mapping = "mapping.tsv"
  )
  paths <- vapply(names(files), function(name) {
    contract_file(file.path(root, files[[name]]), paste0("fixture\t", name))
  }, character(1))
  cfg <- list(
    pipeline = list(
      de_index_path = basename(paths[["de_index"]]),
      contrast_index_path = basename(paths[["contrast_index"]])
    ),
    allowlisted_source_paths = basename(paths[["source"]]),
    source_data = list(list(
      path = basename(paths[["source"]]), role = "sample_metadata"
    ))
  )
  de_index <- data.frame(
    analysis_id = "a", de_path = paths[["de"]],
    matrix_path = paths[["matrix"]], mapped_id_path = paths[["mapped"]],
    stringsAsFactors = FALSE
  )
  policies <- list(de_table_duplicate_policy = list(
    type = "mapping_file", mapping_file = basename(paths[["mapping"]])
  ))

  inventory <- lisaR:::lisa_scientific_input_inventory(
    cfg, root, de_index, data.frame(), policies
  )

  expect_setequal(
    inventory$portable$role,
    c(
      "de_index", "contrast_index", "analysis.de_path",
      "analysis.matrix_path", "analysis.mapped_id_path", "source_data",
      "duplicate_policy.de_table_duplicate_policy.mapping_file"
    )
  )
  expect_setequal(inventory$portable$basename, unname(files))
  expect_false("source_path" %in% names(inventory$portable))
  expect_true(all(startsWith(inventory$private$source_path,
                             normalizePath(root, winslash = "/"))))
  expect_true(lisaR:::lisa_sha256_all_valid(
    inventory$portable$sha256, nrow(inventory$portable)
  ))
})

test_that("transaction writes only portable receipts and revalidates private sources", {
  root <- tempfile("lisa-contract-receipts-")
  dir.create(root)
  config <- contract_file(file.path(root, "study.yml"), "pipeline: baseline")
  input <- contract_file(file.path(root, "de.tsv"), "symbol\tvalue\nA\t1")
  resources <- vapply(seq_len(3L), function(i) {
    contract_file(file.path(root, paste0("resource-", i, ".tsv")),
                  paste0("resource\t", i))
  }, character(1))
  code <- contract_file(file.path(root, "engine.R"), "engine <- 1")
  private_input <- lisaR:::lisa_contract_file_row(
    input, "analysis.de_path", "a"
  )
  contract <- lisaR:::lisa_run_contract(
    config_path = config,
    input_inventory = private_input[, c("role", "owner", "basename", "bytes", "sha256")],
    resource_inventory = contract_resource_rows(resources),
    code_inventory = contract_code_row(code),
    effective_environment = contract_environment()
  )
  attr(contract, "lisaR.private_inputs") <- private_input

  contract_file(config, "pipeline: changed")
  expect_error(
    lisaR:::lisa_transaction_begin(
      file.path(root, "changed-config-result"), contract,
      run_id = "changedconfig"
    ),
    "LISA-CONTRACT-008"
  )
  expect_false(dir.exists(file.path(root, "changed-config-result.lisa.lock")))
  contract_file(config, "pipeline: baseline")

  final <- file.path(root, "result")
  tx <- lisaR:::lisa_transaction_begin(
    final, contract, run_id = "canonicalreceipts"
  )
  on.exit(lisaR:::lisa_transaction_abort(tx, "test cleanup"), add = TRUE)
  receipt_names <- c(
    "run_contract.tsv", "input_inventory.tsv", "resource_inventory.tsv",
    "executable_code_inventory.tsv", "effective_environment.tsv"
  )
  expect_true(all(file.exists(file.path(tx$staging_dir, receipt_names))))
  receipt_text <- paste(vapply(
    file.path(tx$staging_dir, receipt_names),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  ), collapse = "\n")
  expect_false(grepl(normalizePath(root, winslash = "/"), receipt_text,
                     fixed = TRUE))

  lisaR:::lisa_transaction_abort(tx, "first transaction complete")
  identical_target <- contract_file(
    file.path(root, "same-input.tsv"), "symbol\tvalue\nA\t1"
  )
  unlink(input)
  linked <- .Platform$OS.type != "windows" && file.symlink(
    identical_target, input
  )
  if (!isTRUE(linked)) {
    contract_file(input, "symbol\tvalue\nA\t2")
  }
  expect_error(
    lisaR:::lisa_transaction_begin(
      file.path(root, "changed-result"), contract,
      run_id = "changedsource"
    ),
    "LISA-CONTRACT-008"
  )
  expect_false(dir.exists(file.path(root, "changed-result.lisa.lock")))
})

test_that("promotion rejects a resource changed after transaction start", {
  root <- tempfile("lisa-resource-promotion-recheck-")
  dir.create(root)
  config <- contract_file(file.path(root, "study.yml"), "pipeline: baseline")
  input <- contract_file(file.path(root, "de.tsv"), "symbol\tvalue\nA\t1")
  resource_paths <- vapply(seq_len(3L), function(i) {
    contract_file(
      file.path(root, paste0("resource-", i, ".tsv")),
      paste0("resource\t", i)
    )
  }, character(1))
  resources <- lapply(seq_along(resource_paths), function(i) list(
    path = resource_paths[[i]],
    resource_id = paste0("fixture_", i, "@1.0.0"),
    species = "Homo sapiens",
    modality = "transcriptomic/genomic",
    schema = c(
      "lisa_dictionary@2", "term2gene@1", "category_map@1"
    )[[i]],
    sha256 = lisaR:::lisa_sha256_file(resource_paths[[i]])
  ))
  names(resources) <- c(
    "dictionary_resource", "term2gene_resource", "category_map_resource"
  )
  resource_inventory <- lisaR:::lisa_resource_inventory(resources)
  code <- contract_file(file.path(root, "engine.R"), "engine <- 1")
  contract <- lisaR:::lisa_run_contract(
    config_path = config,
    input_inventory = contract_input_row(input),
    resource_inventory = resource_inventory$portable,
    code_inventory = contract_code_row(code),
    effective_environment = contract_environment()
  )
  attr(contract, "lisaR.private_resources") <- resource_inventory$private

  final <- file.path(root, "result")
  tx <- lisaR:::lisa_transaction_begin(
    final, contract, run_id = "changedresourcepromotion"
  )
  on.exit(lisaR:::lisa_transaction_abort(tx, "test cleanup"), add = TRUE)
  writeLines("artifact", file.path(tx$staging_dir, "artifact.txt"))
  contract_file(resource_paths[[2L]], "resource\tchanged")
  expect_error(
    lisaR:::lisa_transaction_promote(tx),
    "LISA-CONTRACT-008 resources changed"
  )
  expect_false(dir.exists(final))
  expect_true(dir.exists(tx$staging_dir))
})

test_that("contract inventories reject lexical symlinks before normalization", {
  root <- tempfile("lisa-contract-symlink-")
  dir.create(root)
  on.exit(lisa_test_cleanup_path(root), add = TRUE)
  target <- contract_file(file.path(root, "target.tsv"), "value\n1")
  link <- file.path(root, "input.tsv")
  skip_if_not(file.symlink(target, link), "symbolic links unavailable")

  expect_error(
    lisaR:::lisa_contract_file_row(link, "analysis.de_path", "a"),
    "LISA-CONTRACT-002"
  )

  resources <- lapply(seq_len(3L), function(i) list(
    path = if (i == 1L) link else target,
    resource_id = paste0("fixture_", i, "@1.0.0"),
    species = "Homo sapiens",
    modality = "transcriptomic/genomic",
    schema = c("lisa_dictionary@2", "term2gene@1", "category_map@1")[[i]],
    sha256 = lisaR:::lisa_sha256_file(target)
  ))
  names(resources) <- c(
    "dictionary_resource", "term2gene_resource", "category_map_resource"
  )
  expect_error(lisaR:::lisa_resource_inventory(resources), "LISA-CONTRACT-003")

  package_root <- file.path(root, "package")
  dir.create(package_root)
  dir.create(file.path(package_root, "R"))
  dir.create(file.path(package_root, "inst", "scripts"), recursive = TRUE)
  contract_file(file.path(package_root, "DESCRIPTION"), "Package: lisaR\nVersion: 0.0.1")
  contract_file(file.path(package_root, "NAMESPACE"), "exportPattern(\"^[[:alpha:]]+\")")
  contract_file(file.path(package_root, "inst", "scripts", "runner.R"), "invisible(TRUE)")
  expect_true(file.symlink(target, file.path(package_root, "R", "linked.R")))
  package_dir <- lisaR:::lisa_resolve_package_dir(
    .test_package_dir = package_root
  )
  expect_error(
    lisaR:::lisa_executable_code_inventory(package_dir),
    "LISA-CONTRACT-005"
  )
})
