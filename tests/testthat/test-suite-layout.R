test_that("risk-based suite and fixture paths are stable", {
  suite_names <- c("unit", "integration", "security", "report-contract", "network-cache")
  source_root <- normalizePath(test_path("..", ".."), mustWork = TRUE)

  expect_true(all(dir.exists(file.path(source_root, "tests", "suites", suite_names))))
  expect_true(all(dir.exists(file.path(source_root, "tests", "fixtures", suite_names))))
})

test_that("network-cache suite has no live-network requirement", {
  expect_false(identical(Sys.getenv("LISAR_ALLOW_NETWORK_TESTS"), "true"))
})
