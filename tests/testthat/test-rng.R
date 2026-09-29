test_that("task seeds are stable and distinct by task identifier", {
  expect_identical(lisa_task_seed(1729L, "analysis-A"), lisa_task_seed(1729L, "analysis-A"))
  expect_false(identical(lisa_task_seed(1729L, "analysis-A"), lisa_task_seed(1729L, "analysis-B")))
})

test_that("RNG helpers preserve the caller RNG state", {
  set.seed(991L)
  before <- .Random.seed
  value <- lisa_with_task_seed(1729L, "state-test", runif(4))

  expect_length(value, 4L)
  expect_identical(.Random.seed, before)
})

test_that("sequential and supported multicore task maps are reproducible", {
  task_ids <- c("A", "B", "C")
  sequential_one <- lisa_map_with_task_seeds(task_ids, function(id) runif(3), backend = "sequential")
  sequential_two <- lisa_map_with_task_seeds(task_ids, function(id) runif(3), backend = "sequential")
  expect_identical(sequential_one, sequential_two)

  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  skip_if_not(parallel::detectCores(logical = FALSE) > 1L, "Multicore verification requires at least two cores.")
  multicore <- lisa_map_with_task_seeds(task_ids, function(id) runif(3), backend = "multicore", workers = 2L)
  expect_identical(multicore, sequential_one)
})

test_that("runtime manifests record RNG, execution, versions, and lock hashes", {
  manifest <- lisa_runtime_manifest(root_seed = 1729L, task_ids = c("A", "B"), backend = "sequential", workers = 1L)

  expect_identical(manifest$root_seed, 1729L)
  expect_equal(manifest$task_seeds$task_id, c("A", "B"))
  expect_match(manifest$rng_kind[[1]], "L'Ecuyer-CMRG", fixed = TRUE)
  expect_identical(manifest$backend, "sequential")
  expect_identical(manifest$workers, 1L)
  expect_identical(manifest$backend_requested, "sequential")
  expect_identical(manifest$backend_effective, "sequential")
  expect_identical(manifest$workers_requested, 1L)
  expect_identical(manifest$workers_effective, 1L)
  expect_identical(manifest$cap_reason, "none")
  expect_identical(manifest$inner_threads, 1L)
  expect_true(all(lisaR:::lisa_inner_thread_variables() %in% names(manifest$thread_environment)))
  expect_true(all(c("R", "lisaR") %in% manifest$versions$package))
  expect_true(all(c("renv.lock", "envs/lisa-core.yml", "envs/lisa-full.yml") %in% manifest$lock_hashes$file))
  expect_true(any(manifest$lock_hashes$file == "envs/LOCK_GENERATION_BLOCKER.md"))
})

test_that("runtime manifests record unavailable optional packages without querying their versions", {
  skip_if(requireNamespace("yaml", quietly = TRUE), "This regression fixture requires yaml to be unavailable.")
  manifest <- lisa_runtime_manifest(task_ids = "unavailable-optional-package")
  yaml_row <- manifest$versions[manifest$versions$package == "yaml", , drop = FALSE]

  expect_equal(nrow(yaml_row), 1L)
  expect_true(is.na(yaml_row$version[[1]]))
})
