lisa_test_run_contract <- function(environment = "R test fixture") {
  root <- tempfile("lisa-test-contract-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  config <- file.path(root, "study.yml")
  dictionary <- file.path(root, "dictionary.tsv")
  writeLines("pipeline: {}", config, useBytes = TRUE)
  writeLines("term\tcategory", dictionary, useBytes = TRUE)
  lisaR:::lisa_run_contract(
    config_path = config,
    dictionary_path = dictionary,
    environment = environment
  )
}

lisa_test_cleanup_path <- function(path) {
  tryCatch({
    if (isTRUE(unname(fs::link_exists(path)))) {
      fs::link_delete(path)
    } else if (isTRUE(unname(fs::dir_exists(path)))) {
      fs::dir_delete(path)
    } else if (isTRUE(unname(fs::file_exists(path)))) {
      fs::file_delete(path)
    }
  }, error = function(error) {
    warning(
      "Test cleanup left an entry in place because no-follow inspection failed: ",
      conditionMessage(error), call. = FALSE
    )
  })
  invisible(path)
}
