test_that("standalone dependency profiles match DESCRIPTION", {
  profile_path <- system.file(
    "standalone",
    "dependency-profiles.tsv",
    package = "lisaR"
  )
  expect_true(file.exists(profile_path))
  profiles <- utils::read.delim(
    profile_path,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  description <- read.dcf(
    system.file("DESCRIPTION", package = "lisaR"),
    fields = c("Imports", "Suggests")
  )
  declared <- paste(description[1L, ], collapse = ",")
  declared <- trimws(unlist(strsplit(gsub("\\s*\\([^)]*\\)", "", declared), ",")))
  declared <- declared[nzchar(declared)]

  expect_true(all(profiles$package %in% declared))
  expect_setequal(
    profiles$package[profiles$profile == "core"],
    c(
      "digest", "fgsea", "fs", "ggplot2", "hommel", "jsonlite", "openxlsx", "readr",
      "yaml"
    )
  )
  installed <- utils::installed.packages()
  system_packages <- rownames(installed)[
    !is.na(installed[, "Priority"]) &
      installed[, "Priority"] %in% c("base", "recommended")
  ]
  runtime_imports <- setdiff(
    trimws(unlist(strsplit(
      gsub("\\s*\\([^)]*\\)", "", description[1L, "Imports"]), ","
    ))),
    system_packages
  )
  expect_true(all(runtime_imports %in%
    profiles$package[profiles$profile == "core"]))
  expect_setequal(
    profiles$package[profiles$profile == "validation"],
    "testthat"
  )
  expect_setequal(
    unique(profiles$profile),
    c("core", "full", "validation")
  )
  expect_false(
    "testthat" %in% profiles$package[
      profiles$profile %in% c("core", "full")
    ]
  )
  expect_setequal(unique(profiles$manager), c("CRAN", "Bioconductor"))
})

test_that("standalone installer is shipped beside its dependency profile", {
  installer <- system.file(
    "standalone",
    "install_lisaR.R",
    package = "lisaR"
  )
  expect_true(file.exists(installer))
  expect_silent(parse(installer))
  installer_text <- paste(readLines(installer, warn = FALSE), collapse = "\n")
  expect_match(installer_text, "core\\|full\\|validation")
  expect_match(
    installer_text,
    'identical\\(profile, "validation"\\)'
  )
  expect_match(installer_text, "package_dependencies")
  expect_match(installer_text, "recursive = TRUE", fixed = TRUE)
  expect_match(installer_text, "missing_closure")
  expect_match(installer_text, "requireNamespace")
})
