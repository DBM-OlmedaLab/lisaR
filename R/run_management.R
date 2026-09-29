# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Run management transactional execution and verified run artifacts.

lisa_resume_abort <- function() {
  code <- "LISA-RESUME-001"
  condition <- structure(
    list(
      message = paste(
        code,
        "resume=TRUE is not supported because lisaR cannot yet prove that an",
        "abandoned staging tree is complete and compatible; start a new run",
        "with resume=false and a new run_id. The abandoned staging tree is",
        "preserved for inspection."
      ),
      call = NULL,
      code = code
    ),
    class = c("lisa_resume_error", "lisa_error", "error", "condition")
  )
  stop(condition)
}

lisa_assert_resume_disabled <- function(resume) {
  if (isTRUE(resume)) lisa_resume_abort()
  invisible(TRUE)
}

lisa_assert_new_final_output <- function(output_dir) {
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop(
      "LISA-RUN-005 final output directory already exists and is immutable.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lisa_sha256_serialized_text <- function(text, line_separator = "\n") {
  if (!is.character(line_separator) || length(line_separator) != 1L ||
      is.na(line_separator) || !line_separator %in% c("\n", "\r\n")) {
    stop("Text SHA-256 serialization requires LF or CRLF separators.",
         call. = FALSE)
  }
  path <- lisa_atomic_temp_path(
    tempdir(), "sha-text", destination = "text.txt", create_file = TRUE
  )
  on.exit(lisa_remove_private_temp(path), add = TRUE)
  values <- enc2utf8(as.character(text))
  if (identical(line_separator, "\r\n")) {
    # Reproduce the exact historical Windows text-connection translation for
    # read compatibility: embedded LF bytes and writeLines separators both
    # became CRLF. Canonical creation never enters this branch.
    values <- gsub("\n", "\r\n", values, fixed = TRUE)
  }
  serialized <- if (length(values)) {
    paste0(paste(values, collapse = line_separator), line_separator)
  } else {
    ""
  }
  connection <- file(path, open = "wb")
  connection_open <- TRUE
  on.exit(if (connection_open) close(connection), add = TRUE)
  writeBin(charToRaw(serialized), connection)
  close(connection)
  connection_open <- FALSE
  lisa_sha256_file(path)
}

lisa_sha256_text <- function(text) {
  # Canonical text hashing is byte-exact UTF-8 with LF between vector elements
  # and one terminal LF for a non-empty vector. A binary connection prevents a
  # native Windows text connection from translating those separators to CRLF.
  lisa_sha256_serialized_text(text, line_separator = "\n")
}

lisa_contract_escape <- function(value) {
  value <- as.character(value)
  value[is.na(value)] <- "<NA>"
  value <- gsub("\\\\", "\\\\\\\\", value)
  value <- gsub("\t", "\\\\t", value, fixed = TRUE)
  value <- gsub("\r", "\\\\r", value, fixed = TRUE)
  gsub("\n", "\\\\n", value, fixed = TRUE)
}

lisa_contract_order_table <- function(table, columns) {
  table <- as.data.frame(table, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(columns, names(table))
  if (length(missing)) {
    stop(
      "LISA-CONTRACT-001 inventory is missing canonical column(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  table <- table[, columns, drop = FALSE]
  escaped <- lapply(table, lisa_contract_escape)
  if (nrow(table)) {
    ordering <- do.call(order, c(escaped, list(method = "radix")))
    table <- table[ordering, , drop = FALSE]
  }
  row.names(table) <- NULL
  table
}

lisa_contract_canonical_table <- function(table, columns) {
  table <- lisa_contract_order_table(table, columns)
  escaped <- lapply(table, lisa_contract_escape)
  names(escaped) <- columns
  header <- paste(columns, collapse = "\t")
  if (!nrow(table)) return(header)
  rows <- vapply(seq_len(nrow(table)), function(i) {
    paste(vapply(escaped, `[[`, character(1), i), collapse = "\t")
  }, character(1))
  paste(c(header, rows), collapse = "\n")
}

lisa_contract_table_hash <- function(table, columns) {
  lisa_sha256_text(lisa_contract_canonical_table(table, columns))
}

lisa_contract_file_identity <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("LISA-CONTRACT-002 contract input paths must be one non-empty string.",
         call. = FALSE)
  }
  lexical_path <- lisa_path_walk_lexically(path)$path
  if (lisa_path_is_link(lexical_path) ||
      !file.exists(lexical_path) || dir.exists(lexical_path) ||
      !isTRUE(utils::file_test("-f", lexical_path))) {
    stop(
      "LISA-CONTRACT-002 contract inputs must be regular non-symlink files: ",
      basename(lexical_path), ".",
      call. = FALSE
    )
  }
  path <- normalizePath(lexical_path, winslash = "/", mustWork = TRUE)
  list(
    path = path,
    basename = basename(lexical_path),
    bytes = as.numeric(file.info(path)$size),
    sha256 = lisa_sha256_file(path)
  )
}

lisa_contract_file_row <- function(path, role, owner) {
  identity <- lisa_contract_file_identity(path)
  data.frame(
    role = as.character(role),
    owner = as.character(owner),
    basename = identity$basename,
    bytes = identity$bytes,
    sha256 = identity$sha256,
    source_path = identity$path,
    stringsAsFactors = FALSE
  )
}

lisa_scientific_input_inventory <- function(cfg, config_dir,
                                            de_index, contrast_index = NULL,
                                            duplicate_policies = list()) {
  rows <- list()
  add <- function(path, role, owner) {
    if (is.null(path) || !length(path)) {
      return(invisible(NULL))
    }
    path <- as.character(path)
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      return(invisible(NULL))
    }
    rows[[length(rows) + 1L]] <<- lisa_contract_file_row(path, role, owner)
    invisible(NULL)
  }

  pipeline <- lisa_config_get(cfg, "pipeline", list())
  external_indices <- list(
    de_index = lisa_config_get(
      pipeline, "de_index_path", lisa_config_get(cfg, "de_index_path", NULL)
    ),
    contrast_index = lisa_config_get(
      pipeline, "contrast_index_path",
      lisa_config_get(cfg, "contrast_index_path", NULL)
    )
  )
  for (role in names(external_indices)) {
    value <- external_indices[[role]]
    if (!is.null(value) && length(value) == 1L && !is.na(value) &&
        nzchar(as.character(value))) {
      add(lisa_config_path(value, config_dir), role, "study")
    }
  }

  add_index_paths <- function(table, role_prefix, owner_columns) {
    if (is.null(table) || !nrow(table)) return(invisible(NULL))
    path_columns <- grep(
      "(^path$|_path$|_file$|_excel$|_rds$)",
      names(table), value = TRUE
    )
    owner_column <- intersect(owner_columns, names(table))
    for (i in seq_len(nrow(table))) {
      owner <- if (length(owner_column)) {
        as.character(table[[owner_column[[1L]]]][[i]])
      } else {
        sprintf("row_%d", i)
      }
      for (column in path_columns) {
        value <- as.character(table[[column]][[i]])
        if (!is.na(value) && nzchar(value)) {
          add(value, paste0(role_prefix, ".", column), owner)
        }
      }
    }
    invisible(NULL)
  }
  add_index_paths(de_index, "analysis", c("analysis_id"))
  add_index_paths(
    contrast_index, "contrast", c("contrast_id", "comparison_id")
  )

  prepared_sources <- lisa_prepare_config_source_data(cfg, config_dir)
  for (i in seq_along(prepared_sources)) {
    source <- prepared_sources[[i]]
    source_role <- as.character(source$source_role %||% "")
    owner <- if (length(source_role) == 1L && !is.na(source_role) &&
                 nzchar(source_role)) {
      source_role
    } else {
      sprintf("source_data[%d]", i)
    }
    add(source$source_path, "source_data", owner)
  }

  for (policy_name in names(duplicate_policies)) {
    policy <- duplicate_policies[[policy_name]]
    if (identical(lisa_policy_type(policy), "mapping_file")) {
      add(
        lisa_config_path(policy$mapping_file, config_dir),
        paste0("duplicate_policy.", policy_name, ".mapping_file"),
        "study"
      )
    }
  }

  private <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      role = character(), owner = character(), basename = character(),
      bytes = numeric(), sha256 = character(), source_path = character(),
      stringsAsFactors = FALSE
    )
  }
  canonical_columns <- c("role", "owner", "basename", "bytes", "sha256")
  if (nrow(private)) {
    key <- do.call(order, c(
      lapply(private[, canonical_columns, drop = FALSE], as.character),
      list(method = "radix")
    ))
    private <- private[key, , drop = FALSE]
  }
  row.names(private) <- NULL
  list(
    portable = private[, canonical_columns, drop = FALSE],
    private = private
  )
}

lisa_resource_inventory <- function(resources) {
  spec <- lisa_pipeline_resource_spec()
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    role <- spec$config_key[[i]]
    resource <- resources[[role]]
    if (is.null(resource) || is.null(resource$path)) {
      stop("LISA-CONTRACT-003 resolved resource is absent: ", role, ".",
           call. = FALSE)
    }
    lexical_path <- lisa_path_walk_lexically(resource$path)$path
    if (lisa_path_is_link(lexical_path) ||
        !file.exists(lexical_path) || dir.exists(lexical_path) ||
        !isTRUE(utils::file_test("-f", lexical_path))) {
      stop(
        "LISA-CONTRACT-003 resolved resource is not one regular non-symlink file: ",
        role, ".", call. = FALSE
      )
    }
    path <- normalizePath(lexical_path, winslash = "/", mustWork = TRUE)
    observed <- lisa_sha256_file(path)
    expected <- as.character(resource$sha256)
    if (length(expected) != 1L || !lisa_sha256_all_valid(expected, 1L) ||
        !identical(unname(observed), unname(expected))) {
      stop(
        "LISA-CONTRACT-004 resolved resource changed after registry verification: ",
        role, ".",
        call. = FALSE
      )
    }
    data.frame(
      role = role,
      resource_id = as.character(resource$resource_id),
      basename = basename(path),
      species = as.character(resource$species),
      modality = as.character(resource$modality),
      schema = as.character(resource$schema),
      bytes = as.numeric(file.info(path)$size),
      sha256 = observed,
      source_path = path,
      stringsAsFactors = FALSE
    )
  })
  private <- do.call(rbind, rows)
  canonical_columns <- c(
    "role", "resource_id", "basename", "species", "modality", "schema",
    "bytes", "sha256"
  )
  private <- private[order(private$role, method = "radix"), , drop = FALSE]
  row.names(private) <- NULL
  list(
    portable = private[, canonical_columns, drop = FALSE],
    private = private
  )
}

lisa_executable_code_inventory <- function(
  package_dir = lisa_resolve_package_dir(.test_package_dir = NULL)
) {
  package_root <- normalizePath(
    as.character(package_dir)[[1L]], winslash = "/", mustWork = TRUE
  )
  script_dir <- lisa_package_script_dir(package_dir)
  roots <- c(
    r_code = file.path(package_root, "R"),
    scripts = script_dir,
    native_code = file.path(package_root, "libs")
  )
  rows <- list()
  resolve_code_file <- function(path) {
    lexical <- gsub("\\\\", "/", as.character(path))
    if (!lisa_path_within(lexical, package_root)) {
      stop(
        "LISA-CONTRACT-005 executable code escaped the package root: ",
        basename(lexical), ".", call. = FALSE
      )
    }
    relative <- substring(lexical, nchar(package_root) + 2L)
    current <- package_root
    for (part in strsplit(relative, "/", fixed = TRUE)[[1L]]) {
      if (!nzchar(part)) next
      current <- paste0(sub("/+$", "", current), "/", part)
      if (lisa_path_is_link(current)) {
        stop(
          "LISA-CONTRACT-005 executable code must not contain symlinks: ",
          basename(lexical), ".", call. = FALSE
        )
      }
    }
    resolved <- normalizePath(lexical, winslash = "/", mustWork = TRUE)
    if (!lisa_path_within(resolved, package_root) || dir.exists(resolved) ||
        !isTRUE(utils::file_test("-f", resolved))) {
      stop(
        "LISA-CONTRACT-005 executable code is not one regular package file: ",
        basename(lexical), ".", call. = FALSE
      )
    }
    resolved
  }
  for (role in names(roots)) {
    root <- roots[[role]]
    if (!dir.exists(root)) next
    files <- sort(list.files(root, recursive = TRUE, full.names = TRUE,
                             all.files = TRUE, no.. = TRUE))
    files <- files[file.exists(files) & !dir.exists(files)]
    for (path in files) {
      path <- resolve_code_file(path)
      relative <- if (startsWith(path, paste0(package_root, "/"))) {
        substring(path, nchar(package_root) + 2L)
      } else {
        file.path("scripts", basename(path))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        role = role,
        portable_id = gsub("\\\\", "/", relative),
        basename = basename(path),
        bytes = as.numeric(file.info(path)$size),
        sha256 = lisa_sha256_file(path),
        source_path = path,
        stringsAsFactors = FALSE
      )
    }
  }
  for (name in c("DESCRIPTION", "NAMESPACE")) {
    path <- file.path(package_root, name)
    if (!file.exists(path)) next
    path <- resolve_code_file(path)
    rows[[length(rows) + 1L]] <- data.frame(
      role = "package_metadata", portable_id = name, basename = name,
      bytes = as.numeric(file.info(path)$size), sha256 = lisa_sha256_file(path),
      source_path = path,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    stop("LISA-CONTRACT-006 no executable lisaR code files were found.",
         call. = FALSE)
  }
  private <- do.call(rbind, rows)
  private <- private[order(private$role, private$portable_id, method = "radix"),
                     , drop = FALSE]
  row.names(private) <- NULL
  canonical_columns <- c(
    "role", "portable_id", "basename", "bytes", "sha256"
  )
  list(
    portable = private[, canonical_columns, drop = FALSE],
    private = private
  )
}

lisa_effective_environment <- function(
  dependency_status = data.frame(),
  installed_db = utils::installed.packages(),
  runtime_version = R.version.string,
  runtime_platform = R.version$platform
) {
  installed_db <- as.matrix(installed_db)
  required_columns <- c("Package", "Version", "Depends", "Imports")
  if (!all(required_columns %in% colnames(installed_db))) {
    stop(
      "LISA-CONTRACT-007 installed package metadata is incomplete; required columns are: ",
      paste(required_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  active_packages <- sort(unique(c(
    "lisaR", as.character(dependency_status$package %||% character())
  )))
  active_packages <- active_packages[
    !is.na(active_packages) & nzchar(active_packages) & active_packages != "R"
  ]
  absent <- setdiff(active_packages, installed_db[, "Package"])
  if (length(absent)) {
    stop(
      "LISA-CONTRACT-007 effective environment is missing required package metadata: ",
      paste(absent, collapse = ", "), ".",
      call. = FALSE
    )
  }
  closure <- tools::package_dependencies(
    active_packages,
    db = installed_db,
    # LinkingTo is a build-time relationship. Binary installations need not
    # retain those toolchain packages, so they are not part of the effective
    # runtime environment.
    which = c("Depends", "Imports"),
    recursive = TRUE
  )
  packages <- sort(unique(c(active_packages, unlist(closure, use.names = FALSE))))
  packages <- packages[!is.na(packages) & nzchar(packages) & packages != "R"]
  absent <- setdiff(packages, installed_db[, "Package"])
  if (length(absent)) {
    stop(
      "LISA-CONTRACT-007 dependency closure is missing installed metadata: ",
      paste(absent, collapse = ", "), ".",
      call. = FALSE
    )
  }
  versions <- installed_db[
    match(packages, installed_db[, "Package"]), "Version"
  ]
  reasons <- stats::setNames(rep("dependency closure", length(packages)), packages)
  reasons[["lisaR"]] <- "executing package"
  if (nrow(dependency_status)) {
    for (i in seq_len(nrow(dependency_status))) {
      reasons[[as.character(dependency_status$package[[i]])]] <-
        as.character(dependency_status$reason[[i]])
    }
  }
  out <- rbind(
    data.frame(
      component = "runtime", name = "R", version = runtime_version,
      platform = runtime_platform, requirement = "R runtime",
      stringsAsFactors = FALSE
    ),
    data.frame(
      component = "package", name = packages, version = unname(versions),
      platform = "", requirement = unname(reasons[packages]),
      stringsAsFactors = FALSE
    )
  )
  row.names(out) <- NULL
  out
}

lisa_configured_run_contract <- function(
  config_path, cfg, config_dir, de_index, contrast_index, duplicate_policies,
  resources, dependency_status,
  package_dir = lisa_resolve_package_dir(.test_package_dir = NULL)
) {
  input_inventory <- lisa_scientific_input_inventory(
    cfg, config_dir, de_index, contrast_index,
    duplicate_policies = duplicate_policies
  )
  resource_inventory <- lisa_resource_inventory(resources)
  code_inventory <- lisa_executable_code_inventory(package_dir)
  effective_environment <- lisa_effective_environment(dependency_status)
  contract <- lisa_run_contract(
    config_path = config_path,
    input_inventory = input_inventory$portable,
    resource_inventory = resource_inventory$portable,
    code_inventory = code_inventory$portable,
    effective_environment = effective_environment
  )
  attr(contract, "lisaR.private_inputs") <- input_inventory$private
  attr(contract, "lisaR.private_resources") <- resource_inventory$private
  attr(contract, "lisaR.private_code") <- code_inventory$private
  contract
}

lisa_run_contract <- function(config_path = NULL, dictionary_path = NULL,
                              environment = R.version.string,
                              input_inventory = NULL,
                              resource_inventory = NULL,
                              code_inventory = NULL,
                              effective_environment = NULL) {
  detailed <- !is.null(input_inventory) || !is.null(resource_inventory) ||
    !is.null(code_inventory) || !is.null(effective_environment)
  if (!detailed) {
    values <- c(
      config = if (!is.null(config_path) && file.exists(config_path)) lisa_sha256_file(config_path) else "",
      code = lisa_sha256_text(paste(utils::packageDescription("lisaR")$Version %||% "development", lisa_contract_versions()$pipeline_schema_version)),
      dictionary = if (!is.null(dictionary_path) && file.exists(dictionary_path)) lisa_sha256_file(dictionary_path) else "",
      environment = lisa_sha256_text(environment)
    )
    return(c(values, contract = lisa_sha256_text(paste(values, collapse = "|"))))
  }

  inventories <- list(
    input_inventory = input_inventory,
    resource_inventory = resource_inventory,
    code_inventory = code_inventory,
    effective_environment = effective_environment
  )
  if (any(vapply(inventories, is.null, logical(1)))) {
    stop(
      "LISA-CONTRACT-001 a detailed run contract requires input, resource, code, and effective-environment inventories.",
      call. = FALSE
    )
  }
  if (is.null(config_path) || !is.character(config_path) ||
      length(config_path) != 1L || is.na(config_path) ||
      !file.exists(config_path) || dir.exists(config_path)) {
    stop(
      "LISA-CONTRACT-001 a detailed run contract requires one existing config file.",
      call. = FALSE
    )
  }
  config_source <- lisa_contract_file_row(config_path, "config", "study")

  input_columns <- c("role", "owner", "basename", "bytes", "sha256")
  resource_columns <- c(
    "role", "resource_id", "basename", "species", "modality", "schema",
    "bytes", "sha256"
  )
  code_columns <- c("role", "portable_id", "basename", "bytes", "sha256")
  environment_columns <- c(
    "component", "name", "version", "platform", "requirement"
  )
  inventory_specs <- list(
    input_inventory = input_columns,
    resource_inventory = resource_columns,
    code_inventory = code_columns
  )
  for (label in names(inventory_specs)) {
    table <- as.data.frame(inventories[[label]], stringsAsFactors = FALSE)
    columns <- inventory_specs[[label]]
    missing <- setdiff(columns, names(table))
    if (length(missing) ||
        !lisa_sha256_all_valid(as.character(table$sha256), nrow(table)) ||
        anyNA(table$bytes) || any(!is.finite(as.numeric(table$bytes))) ||
        any(as.numeric(table$bytes) < 0)) {
      stop("LISA-CONTRACT-001 invalid ", label, ".", call. = FALSE)
    }
  }
  values <- c(
    config = config_source$sha256[[1L]],
    inputs = lisa_contract_table_hash(input_inventory, input_columns),
    resources = lisa_contract_table_hash(resource_inventory, resource_columns),
    code = lisa_contract_table_hash(code_inventory, code_columns),
    environment = lisa_contract_table_hash(
      effective_environment, environment_columns
    ),
    contract = ""
  )
  values[["contract"]] <- lisa_sha256_text(paste(
    names(values)[names(values) != "contract"],
    values[names(values) != "contract"], sep = "=", collapse = "\n"
  ))
  attr(values, "lisaR.contract_receipts") <- list(
    input_inventory.tsv = lisa_contract_order_table(
      input_inventory, input_columns
    ),
    resource_inventory.tsv = lisa_contract_order_table(
      resource_inventory, resource_columns
    ),
    executable_code_inventory.tsv = lisa_contract_order_table(
      code_inventory, code_columns
    ),
    effective_environment.tsv = lisa_contract_order_table(
      effective_environment, environment_columns
    )
  )
  attr(values, "lisaR.private_config") <- config_source
  values
}

lisa_transaction_contract_abort <- function(message) {
  code <- "LISA-CONTRACT-010"
  condition <- structure(
    list(
      message = paste(code, message),
      call = NULL,
      code = code
    ),
    class = c("lisa_contract_error", "lisa_error", "error", "condition")
  )
  stop(condition)
}

lisa_validate_transaction_contract <- function(
  contract, allow_legacy_windows_crlf = FALSE
) {
  legacy_shape <- c(
    "config", "code", "dictionary", "environment", "contract"
  )
  detailed_shape <- c(
    "config", "inputs", "resources", "code", "environment", "contract"
  )
  observed_names <- names(contract)
  known_shape <- is.character(contract) &&
    (identical(observed_names, legacy_shape) ||
       identical(observed_names, detailed_shape))
  if (!known_shape) {
    lisa_transaction_contract_abort(
      "a transactional run contract must use one complete canonical shape."
    )
  }
  if (!lisa_sha256_all_valid(unname(contract), length(contract))) {
    lisa_transaction_contract_abort(
      paste(
        "every transactional run-contract value must be one lowercase",
        "64-hex SHA-256 digest; empty, missing, or malformed values are not",
        "verifiable."
      )
    )
  }

  components <- contract[observed_names != "contract"]
  aggregate_text <- if (identical(observed_names, legacy_shape)) {
    paste(unname(components), collapse = "|")
  } else {
    paste(
      names(components), unname(components), sep = "=", collapse = "\n"
    )
  }
  expected <- lisa_sha256_text(aggregate_text)
  observed <- unname(contract[["contract"]])
  legacy_windows <- isTRUE(allow_legacy_windows_crlf) &&
    identical(observed, lisa_sha256_serialized_text(
      aggregate_text, line_separator = "\r\n"
    ))
  if (!identical(observed, expected) && !legacy_windows) {
    lisa_transaction_contract_abort(
      "the aggregate contract digest does not match its component digests."
    )
  }
  invisible(TRUE)
}

lisa_read_transaction_contract <- function(
  run_dir, allow_legacy_windows_crlf = TRUE
) {
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  path <- file.path(run_dir, "run_contract.tsv")
  if (!file.exists(path) || dir.exists(path) || !utils::file_test("-f", path)) {
    lisa_transaction_contract_abort(
      "run_contract.tsv is absent or is not a regular file."
    )
  }
  receipt <- tryCatch(
    utils::read.delim(
      path,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character"
    ),
    error = function(error) {
      lisa_transaction_contract_abort(
        paste0("run_contract.tsv cannot be read: ", conditionMessage(error))
      )
    }
  )
  if (!is.data.frame(receipt) ||
      !identical(names(receipt), c("key", "value")) ||
      !is.character(receipt$key) || anyNA(receipt$key) ||
      any(!nzchar(receipt$key)) || anyDuplicated(receipt$key)) {
    lisa_transaction_contract_abort(
      "run_contract.tsv does not contain one canonical key/value table."
    )
  }
  contract <- stats::setNames(as.character(receipt$value), receipt$key)
  # Read-only compatibility is limited to the exact aggregate serialization
  # written historically by a native Windows text connection. New in-memory
  # contracts and transaction creation remain canonical-LF only.
  lisa_validate_transaction_contract(
    contract,
    allow_legacy_windows_crlf = allow_legacy_windows_crlf
  )
  contract
}

lisa_verify_contract_sources <- function(contract) {
  groups <- list(
    config = attr(contract, "lisaR.private_config", exact = TRUE),
    inputs = attr(contract, "lisaR.private_inputs", exact = TRUE),
    resources = attr(contract, "lisaR.private_resources", exact = TRUE),
    code = attr(contract, "lisaR.private_code", exact = TRUE)
  )
  for (group in names(groups)) {
    inventory <- groups[[group]]
    if (is.null(inventory) || !nrow(inventory)) next
    observed <- tryCatch(
      lapply(inventory$source_path, lisa_contract_file_identity),
      error = identity
    )
    failed_identity <- inherits(observed, "error")
    if (!failed_identity) {
      observed_sha256 <- vapply(observed, `[[`, character(1), "sha256")
      observed_bytes <- vapply(observed, `[[`, numeric(1), "bytes")
    }
    if (failed_identity ||
        !identical(unname(observed_sha256), as.character(inventory$sha256)) ||
        !identical(unname(observed_bytes), as.numeric(inventory$bytes))) {
      stop(
        "LISA-CONTRACT-008 ", group,
        " changed or ceased to be a regular non-symlink file after the canonical contract was constructed.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

lisa_new_run_id <- function(contract) {
  paste0("run-", format(Sys.time(), "%Y%m%dT%H%M%S"), "-", substr(lisa_sha256_text(paste(contract, stats::runif(1), sep = "|")), 1L, 12L))
}

lisa_transaction_begin <- function(output_dir, contract, run_id = NULL,
                                   resume = FALSE,
                                   artifact_kind = c("scientific_run", "plan")) {
  lisa_assert_resume_disabled(resume)
  artifact_kind <- match.arg(artifact_kind)
  # Validate before resolving or creating the destination parent. An
  # incomplete legacy contract may still be inspected by callers, but it can
  # never acquire a lock or become a verifiable transactional artifact.
  lisa_validate_transaction_contract(contract)
  destination <- lisa_managed_destination(output_dir, create_parent = TRUE)
  output_dir <- destination$path
  parent <- destination$parent
  lisa_assert_new_final_output(output_dir)
  space_check <- getOption("lisaR.staging_space_check", NULL)
  if (is.function(space_check) && !isTRUE(space_check(parent))) {
    stop("LISA-RUN-007 insufficient staging disk space.", call. = FALSE)
  }
  lock <- paste0(output_dir, ".lisa.lock")
  if (lisa_path_entry_exists(lock) ||
      !dir.create(lock, showWarnings = FALSE, mode = "0700")) {
    stop("LISA-RUN-001 an exclusive run lock already exists: ", lock,
         call. = FALSE)
  }
  lisa_set_private_mode(lock, directory = TRUE)
  lock_owned <- TRUE
  on.exit(if (lock_owned && lisa_path_entry_exists(lock)) {
    try(lisa_guarded_delete(lock, recursive = TRUE, run_root = NULL),
        silent = TRUE)
  }, add = TRUE)
  # Re-hash every private source after acquiring the exclusive destination
  # lock. Absolute paths remain process-local; only portable receipts are
  # written into the run.
  lisa_verify_contract_sources(contract)
  run_id <- run_id %||% lisa_new_run_id(contract[["contract"]])
  lisa_safe_id(run_id, "run_id")
  staging <- lisa_short_staging_path(output_dir, run_id)
  if (lisa_path_entry_exists(staging)) {
    stop(
      "LISA-RUN-002 staging already exists and will not be reused; choose a new run_id. The existing staging tree is preserved for inspection.",
      call. = FALSE
    )
  }
  staging <- lisa_run_root(staging)
  old_root <- getOption("lisaR.run_root", NULL); options(lisaR.run_root = staging)
  on.exit(options(lisaR.run_root = old_root), add = TRUE)
  write_lisa_tsv(data.frame(key = names(contract), value = unname(contract), stringsAsFactors = FALSE), file.path(staging, "run_contract.tsv"))
  receipts <- attr(contract, "lisaR.contract_receipts", exact = TRUE)
  if (!is.null(receipts)) {
    expected_receipts <- c(
      "input_inventory.tsv", "resource_inventory.tsv",
      "executable_code_inventory.tsv", "effective_environment.tsv"
    )
    if (!identical(names(receipts), expected_receipts) ||
        any(!vapply(receipts, is.data.frame, logical(1)))) {
      stop("LISA-CONTRACT-009 canonical contract receipts are malformed.",
           call. = FALSE)
    }
    for (receipt_name in expected_receipts) {
      write_lisa_tsv(
        receipts[[receipt_name]], file.path(staging, receipt_name)
      )
    }
  }
  tx <- list(
    run_id = run_id,
    output_dir = output_dir,
    staging_dir = staging,
    lock = lock,
    contract = contract,
    artifact_kind = artifact_kind
  )
  class(tx) <- "lisa_transaction"
  lisa_transaction_event(tx, "started")
  lock_owned <- FALSE
  tx
}

lisa_transaction_event <- function(tx, event, detail = "") {
  event_dir <- if (dir.exists(tx$staging_dir)) {
    lisa_existing_run_root(tx$staging_dir)
  } else if (dir.exists(tx$output_dir)) {
    lisa_existing_run_root(tx$output_dir)
  } else {
    return(invisible(tx))
  }
  old_root <- getOption("lisaR.run_root", NULL); options(lisaR.run_root = event_dir)
  on.exit(options(lisaR.run_root = old_root), add = TRUE)
  path <- file.path(event_dir, "run_events.tsv")
  prior <- if (file.exists(path)) read_lisa_tsv(path) else data.frame(timestamp = character(), event = character(), detail = character(), stringsAsFactors = FALSE)
  write_lisa_tsv(rbind(prior, data.frame(timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), event = event, detail = as.character(detail), stringsAsFactors = FALSE)), path)
  invisible(tx)
}

lisa_transaction_abort <- function(tx, detail = "") {
  # Idempotent: callers can invoke abort after the promote step or the outer
  # pipeline on.exit have already released the lock; in that case this is a no-op.
  if (!lisa_path_entry_exists(tx$lock)) return(invisible(tx))
  try(lisa_transaction_event(tx, "failed", detail), silent = TRUE)
  lisa_guarded_delete(tx$lock, recursive = TRUE, run_root = NULL)
  invisible(tx)
}

lisa_sha256_files_batch <- function(files, run_root = getOption("lisaR.run_root", NULL),
                                    .backend = lisa_sha256_digest_backend) {
  files <- lisa_sha256_regular_files(files, run_root = run_root)
  if (!length(files)) return(stats::setNames(character(), character()))
  # Bound the number of paths materialised by each in-process digest call. This
  # keeps memory predictable for manifests containing tens of thousands of
  # files without exposing paths to an operating-system argument vector.
  batch_size <- getOption("lisaR.sha256_batch_size", 128L)
  if (length(batch_size) != 1L || is.na(batch_size) || !is.numeric(batch_size) || !is.finite(batch_size) ||
      batch_size < 1 || batch_size != as.integer(batch_size)) {
    lisa_sha256_abort(
      "LISA-SHA256-005",
      "lisaR.sha256_batch_size must be one positive integer."
    )
  }
  batches <- parallel::splitIndices(length(files), ceiling(length(files) / as.integer(batch_size)))
  values <- lapply(
    batches,
    function(index) lisa_sha256_files_primitive(files[index], .backend = .backend)
  )
  hashes <- unlist(values, use.names = TRUE)
  lisa_sha256_validate_backend_result(hashes, files)
}

lisa_hash_files <- function(files, workers = 1L, run_root = getOption("lisaR.run_root", NULL)) {
  workers <- lisa_validate_workers(workers)
  files <- lisa_sha256_regular_files(files, run_root = run_root)
  if (!length(files)) return(character())
  plan <- lisa_parallel_plan(
    workers = workers,
    tasks = length(files),
    backend_requested = "auto"
  )
  if (plan$workers_effective == 1L) {
    hashes <- lisa_sha256_files_batch(files, run_root = NULL)
    return(unname(lisa_sha256_validate_backend_result(hashes, files)))
  }
  chunks <- parallel::splitIndices(length(files), plan$workers_effective)
  task_ids <- sprintf("hash:%04d", seq_along(chunks))
  values <- tryCatch(
    lisa_map_task_values(
      task_ids,
      function(task_id) {
        index <- chunks[[match(task_id, task_ids)]]
        lisa_sha256_files_batch(files[index], run_root = NULL)
      },
      plan = lisa_parallel_plan(
        workers = plan$workers_effective,
        tasks = length(chunks),
        platform = plan$platform,
        allocation = NULL,
        backend_requested = plan$backend_effective
      ),
      root_seed = NULL,
      preschedule = TRUE
    ),
    error = function(error) {
      lisa_sha256_abort(
        "LISA-SHA256-006",
        paste0("parallel SHA-256 execution failed: ", conditionMessage(error)),
        files
      )
    }
  )
  hashes <- character(length(files))
  for (i in seq_along(chunks)) {
    expected <- files[chunks[[i]]]
    observed <- lisa_sha256_validate_backend_result(values[[i]], expected)
    hashes[chunks[[i]]] <- unname(observed)
  }
  names(hashes) <- files
  unname(lisa_sha256_validate_backend_result(hashes, files))
}

lisa_run_manifest_snapshot <- function(run_dir) {
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  tree <- lisa_scan_run_tree(run_dir)
  files <- tree$path[!tree$isdir]
  rel <- substring(files, nchar(run_dir) + 2L)
  keep <- !rel %in% c("run_manifest.tsv", "run_events.tsv") &
    !startsWith(rel, ".lisa_subprocess/")
  data.frame(path = rel[keep], bytes = tree$bytes[!tree$isdir][keep],
             stringsAsFactors = FALSE)
}

lisa_finalize_subprocess_diagnostics <- function(run_dir) {
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  diagnostics_root <- file.path(run_dir, ".lisa_subprocess")
  if (!dir.exists(diagnostics_root)) return(invisible(FALSE))
  remaining <- lisa_scan_run_tree(diagnostics_root)
  if (nrow(remaining)) {
    stop("Subprocess diagnostics remained active at manifest finalization.", call. = FALSE)
  }
  lisa_guarded_delete(diagnostics_root, recursive = TRUE, run_root = run_dir)
  invisible(TRUE)
}

lisa_run_manifest <- function(run_dir, workers = 1L) {
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  snapshot <- lisa_run_manifest_snapshot(run_dir)
  files <- file.path(run_dir, snapshot$path)
  data.frame(path = snapshot$path, bytes = snapshot$bytes,
             sha256 = lisa_hash_files(files, workers, run_root = run_dir), stringsAsFactors = FALSE)
}

lisa_validate_built_run_manifest <- function(run_dir, manifest) {
  required <- c("path", "bytes", "sha256")
  findings <- character()
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      !is.character(manifest$path) || anyNA(manifest$path) ||
      any(!nzchar(manifest$path)) || anyDuplicated(manifest$path)) {
    findings <- c(findings, "manifest_incompatible_or_duplicate")
  } else if (!lisa_sha256_all_valid(manifest$sha256, nrow(manifest))) {
    findings <- c(findings, "manifest_invalid_sha256")
  } else {
    # Promotion performs its required independent verification pass here.
    # Keep this separate from lisa_run_manifest(): the manifest is constructed
    # once, then all declared artifact bytes are read and SHA-256 checked once
    # against it before promotion.
    run_dir <- lisa_assert_run_tree_safe(run_dir)
    snapshot <- lisa_run_manifest_snapshot(run_dir)
    missing <- setdiff(manifest$path, snapshot$path)
    undeclared <- setdiff(snapshot$path, manifest$path)
    shared <- intersect(manifest$path, snapshot$path)
    declared_bytes <- as.numeric(manifest$bytes[match(shared, manifest$path)])
    observed_bytes <- snapshot$bytes[match(shared, snapshot$path)]
    resized <- shared[is.na(declared_bytes) | declared_bytes != observed_bytes]
    observed_sha256 <- lisa_hash_files(file.path(run_dir, snapshot$path), run_root = run_dir)
    declared_sha256 <- as.character(manifest$sha256[match(shared, manifest$path)])
    actual_sha256 <- observed_sha256[match(shared, snapshot$path)]
    modified <- shared[is.na(declared_sha256) | declared_sha256 != actual_sha256]
    if (length(missing)) findings <- c(findings, paste0("absent:", missing))
    if (length(undeclared)) findings <- c(findings, paste0("undeclared:", undeclared))
    if (length(resized)) findings <- c(findings, paste0("resized:", resized))
    if (length(modified)) findings <- c(findings, paste0("modified:", modified))
  }
  list(gate = if (length(findings)) "FAIL" else "PASS", findings = unique(findings))
}

#' Verify a completed transactional lisaR run
#'
#' @param run_dir A promoted run directory or an explicit staging directory.
#' @return A list containing one gate status and detailed findings.
#' @keywords internal
verify_run <- function(run_dir, workers = 1L) {
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  manifest_path <- file.path(run_dir, "run_manifest.tsv")
  plan_status_path <- file.path(run_dir, "lisa_plan_status.tsv")
  findings <- character()
  contract_valid <- tryCatch(
    {
      lisa_read_transaction_contract(run_dir)
      TRUE
    },
    lisa_contract_error = function(error) FALSE,
    error = function(error) FALSE
  )
  if (!contract_valid) findings <- c(findings, "run_contract_invalid")
  if (!file.exists(manifest_path)) findings <- c(findings, "manifest_absent") else {
    manifest <- read_lisa_tsv(manifest_path)
    required <- c("path", "bytes", "sha256")
    if (!all(required %in% names(manifest)) || !is.character(manifest$path) ||
        anyNA(manifest$path) || any(!nzchar(manifest$path)) ||
        anyDuplicated(manifest$path)) {
      findings <- c(findings, "manifest_incompatible_or_duplicate")
    } else if (!lisa_sha256_all_valid(manifest$sha256, nrow(manifest))) {
      findings <- c(findings, "manifest_invalid_sha256")
    }
    if (length(findings) == 0) {
      # External verification intentionally rehashes the completed run once.
      # Promotion uses lisa_validate_built_run_manifest() instead, so it does
      # not immediately rebuild this same manifest after hashing it.
      actual <- lisa_run_manifest(run_dir, workers = workers)
      missing <- setdiff(manifest$path, actual$path); undeclared <- setdiff(actual$path, manifest$path)
      shared <- intersect(manifest$path, actual$path)
      declared_index <- match(shared, manifest$path)
      actual_index <- match(shared, actual$path)
      changed <- shared[manifest$sha256[declared_index] != actual$sha256[actual_index]]
      declared_bytes <- suppressWarnings(as.numeric(manifest$bytes[declared_index]))
      observed_bytes <- as.numeric(actual$bytes[actual_index])
      resized <- shared[!is.finite(declared_bytes) | declared_bytes != observed_bytes]
      if (length(missing)) findings <- c(findings, paste0("absent:", missing))
      if (length(undeclared)) findings <- c(findings, paste0("undeclared:", undeclared))
      if (length(changed)) findings <- c(findings, paste0("modified:", changed))
      if (length(resized)) findings <- c(findings, paste0("resized:", resized))
    }
  }
  events <- if (file.exists(file.path(run_dir, "run_events.tsv"))) read_lisa_tsv(file.path(run_dir, "run_events.tsv")) else data.frame()
  plan_events <- "event" %in% names(events) &&
    any(events$event %in% c("plan_validated", "plan_promoted"))
  is_plan <- file.exists(plan_status_path) || plan_events
  if (is_plan) {
    findings <- c(findings, "planning_artifact_not_scientific_run")
  } else if (!("event" %in% names(events) && any(events$event == "validated"))) {
    findings <- c(findings, "validation_incomplete")
  }
  gate <- if (length(findings)) "FAIL" else "PASS"
  list(
    gate = gate,
    findings = unique(findings),
    manifest = manifest_path,
    artifact_type = if (is_plan) "plan" else "scientific_run"
  )
}

lisa_transaction_promote <- function(tx, workers = 1L) {
  artifact_kind <- tx$artifact_kind %||% "scientific_run"
  if (!artifact_kind %in% c("scientific_run", "plan")) {
    stop("LISA-RUN-008 transaction artifact kind is invalid.", call. = FALSE)
  }
  is_plan <- identical(artifact_kind, "plan")
  lisa_validate_transaction_contract(tx$contract)
  receipt_contract <- lisa_read_transaction_contract(
    tx$staging_dir, allow_legacy_windows_crlf = FALSE
  )
  receipt_matches_transaction <-
    identical(names(receipt_contract), names(tx$contract)) &&
    identical(unname(receipt_contract), unname(as.character(tx$contract)))
  if (!receipt_matches_transaction) {
    lisa_transaction_contract_abort(
      "run_contract.tsv no longer matches the contract held by its transaction."
    )
  }
  old_root <- getOption("lisaR.run_root", NULL); options(lisaR.run_root = tx$staging_dir)
  on.exit(options(lisaR.run_root = old_root), add = TRUE)
  # Post-processors legitimately record absolute paths while writing inside
  # the transaction staging tree. Rebase those text references to the final
  # immutable run root before hashing and atomic promotion; otherwise a
  # successful promoted run contains paths to a directory that no longer
  # exists.
  lisa_assert_run_tree_safe(tx$staging_dir)
  lisa_rebase_run_text_paths(tx$staging_dir, tx$staging_dir, tx$output_dir)
  # All scientific workers have joined before promotion.  Only this serial
  # coordinator may remove the shared diagnostics root, eliminating the race
  # caused by child-level root cleanup while also preventing diagnostic debris
  # from entering a promoted run.
  lisa_finalize_subprocess_diagnostics(tx$staging_dir)
  # Resources, inputs and executable code are external to the staging tree.
  # Revalidate their byte identities after every worker has joined so a run
  # cannot be promoted if a cache or input changed during computation.
  lisa_verify_contract_sources(tx$contract)
  manifest <- lisa_run_manifest(tx$staging_dir, workers = workers)
  write_lisa_tsv(manifest, file.path(tx$staging_dir, "run_manifest.tsv"))
  # Independently re-read and SHA-256 verify the newly built manifest once.
  # This is deliberately not verify_run(), which remains the external
  # completed-run verifier and should not be invoked by promotion.
  check <- lisa_validate_built_run_manifest(tx$staging_dir, manifest)
  if (!identical(check$gate, "PASS")) stop("LISA-RUN-004 validation gate failed: ", paste(check$findings, collapse = "; "), call. = FALSE)
  lisa_transaction_event(tx, if (is_plan) "plan_validated" else "validated")
  lisa_assert_run_tree_safe(tx$staging_dir)
  lisa_managed_destination(tx$output_dir, create_parent = FALSE)
  if (lisa_path_entry_exists(tx$output_dir)) stop("LISA-RUN-005 final output directory already exists and is immutable.", call. = FALSE)
  promotion_error <- NULL
  promoted <- tryCatch({
    lisa_promote_managed_directory(tx$staging_dir, tx$output_dir)
    TRUE
  }, error = function(error) {
    promotion_error <<- error
    FALSE
  })
  if (!promoted) {
    if (inherits(promotion_error, "lisa_filesystem_recovery_error")) {
      # Do not write into or remove any recovery copy. The independent lock is
      # the only entry released here; the typed condition tells the operator
      # exactly which paths must be reconciled.
      if (lisa_path_entry_exists(tx$lock)) {
        try(lisa_guarded_delete(tx$lock, recursive = TRUE, run_root = NULL),
            silent = TRUE)
      }
      stop(promotion_error)
    }
    detail <- paste0(
      "atomic promotion failed",
      if (!is.null(promotion_error)) {
        paste0(": ", conditionMessage(promotion_error))
      } else {
        ""
      }
    )
    lisa_transaction_abort(tx, detail)
    stop("LISA-RUN-006 ", detail, ".", call. = FALSE)
  }
  lisa_transaction_event(tx, if (is_plan) "plan_promoted" else "promoted")
  lisa_guarded_delete(tx$lock, recursive = TRUE, run_root = NULL)
  invisible(list(
    run_id = tx$run_id,
    output_dir = tx$output_dir,
    gate = if (is_plan) "PLAN_PASS" else "PASS",
    artifact_type = artifact_kind
  ))
}

lisa_rebase_run_text_paths <- function(run_dir, from_root, to_root) {
  from_raw <- as.character(from_root)
  run_dir <- lisa_assert_run_tree_safe(run_dir)
  from_root <- normalizePath(from_root, winslash = "/", mustWork = FALSE)
  to_root <- normalizePath(to_root, winslash = "/", mustWork = FALSE)
  from_candidates <- unique(c(from_raw, from_root))
  if (all(from_candidates == to_root)) return(invisible(character()))

  tree <- lisa_scan_run_tree(run_dir)
  files <- tree$path[!tree$isdir]
  text_extensions <- c(
    "tsv", "csv", "txt", "md", "html", "htm", "json", "jsonl",
    "yaml", "yml", "xml", "log", "r", "rmd"
  )
  extensions <- tolower(tools::file_ext(files))
  files <- files[extensions %in% text_extensions]
  changed <- character()

  for (path in files) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    has_staging_path <- any(vapply(
      from_candidates,
      function(candidate) any(grepl(candidate, lines, fixed = TRUE)),
      logical(1)
    ))
    if (!has_staging_path) next
    rebased <- lines
    for (candidate in from_candidates) {
      rebased <- gsub(candidate, to_root, rebased, fixed = TRUE)
    }
    writeLines(rebased, path, useBytes = TRUE)
    changed <- c(changed, path)
  }
  invisible(changed)
}

lisa_validate_run_identity <- function(de_index, contrast_index = NULL) {
  validate_lisa_de_index(de_index, require_files = FALSE)
  validate_lisa_contrast_index(contrast_index, de_index)
  if (anyDuplicated(de_index$analysis_id)) stop("LISA-IDENTITY-001 duplicate analysis identity.", call. = FALSE)
  invisible(TRUE)
}

lisa_product_plan <- function(evidence, run_gene_level = FALSE, run_kegg_maps = FALSE) {
  plan <- lisa_product_applicability(evidence, run_gene_level, run_kegg_maps)
  plan$compute <- plan$status == "enabled"
  plan
}

# States are a projection of immutable events plus verified artifacts; callers
# must not infer completion merely from a directory being present.
lisa_run_state <- function(run_dir) {
  events <- if (file.exists(file.path(run_dir, "run_events.tsv"))) read_lisa_tsv(file.path(run_dir, "run_events.tsv")) else data.frame()
  if (!nrow(events)) return("not_started")
  if (any(events$event == "plan_promoted")) return("planned")
  if (any(events$event == "promoted")) return(if (identical(verify_run(run_dir)$gate, "PASS")) "completed" else "failed")
  if (any(events$event == "failed")) return("failed")
  if (any(events$event == "plan_validated")) return("planned")
  if (any(events$event == "plan_computed")) return("planned")
  if (any(events$event == "computed")) return("validating")
  if (any(events$event == "planned")) return("running")
  "staged"
}
