# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Explicit local trust boundary for serialized differential-expression inputs.

lisa_rds_abort <- function(code, message) {
  stop(structure(
    list(message = paste0(code, ": ", message), call = NULL, code = code),
    class = c("lisa_rds_error", "lisa_error", "error", "condition")
  ))
}

lisa_validate_trusted_rds <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    lisa_rds_abort(
      "LISA-RDS-001",
      "`trusted_rds` must be one non-missing logical value supplied by the local R caller."
    )
  }
  isTRUE(value)
}

lisa_validate_rds_max_bytes <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value)) {
    lisa_rds_abort(
      "LISA-RDS-002",
      "`rds_max_bytes` must be one positive, finite whole number of bytes."
    )
  }
  as.numeric(value)
}

lisa_is_rds_path <- function(input) {
  is.character(input) && length(input) == 1L && !is.na(input) &&
    grepl("\\.rds$", input, ignore.case = TRUE)
}

lisa_read_rds_once <- function(path) {
  readRDS(path)
}

lisa_rds_object_sha256 <- function(object) {
  path <- lisa_atomic_temp_path(
    tempdir(), "rds-object-digest", destination = "object.bin",
    create_file = TRUE
  )
  on.exit(lisa_remove_private_temp(path), add = TRUE)
  connection <- file(path, open = "wb")
  tryCatch(
    serialize(object, connection, ascii = FALSE, xdr = TRUE, version = 3L),
    finally = close(connection)
  )
  lisa_sha256_file(path)
}

lisa_detect_rds_input_type <- function(object) {
  if (inherits(object, "DESeqDataSet")) return("deseq2_dds")
  if (inherits(object, "DESeqResults")) return("deseq2_results")
  if (is.data.frame(object) || is.matrix(object)) return("de_table")
  if (inherits(object, c("DGELRT", "DGEExact", "TopTags"))) return("edger")
  NA_character_
}

lisa_edger_result_table <- function(object) {
  accepted <- c("DGELRT", "DGEExact", "TopTags")
  if (!inherits(object, accepted)) {
    lisa_rds_abort(
      "LISA-RDS-012",
      paste0(
        "the edgeR input must be a DGELRT, DGEExact or TopTags differential-expression result; class=",
        paste(class(object), collapse = "/"), "."
      )
    )
  }
  table <- object$table
  if (!(is.data.frame(table) || is.matrix(table))) {
    lisa_rds_abort(
      "LISA-RDS-012",
      paste0(
        "the edgeR result class `", paste(class(object), collapse = "/"),
        "` must contain a tabular `table` component."
      )
    )
  }
  dimensions <- dim(table)
  if (length(dimensions) != 2L || anyNA(dimensions) || any(dimensions < 1L)) {
    lisa_rds_abort(
      "LISA-RDS-012",
      "the edgeR result `table` component must have at least one row and one column."
    )
  }
  lower_names <- tolower(colnames(table))
  required <- c("logfc", "pvalue")
  missing <- required[!required %in% lower_names]
  if (length(missing)) {
    lisa_rds_abort(
      "LISA-RDS-012",
      paste0(
        "the edgeR result `table` component is missing differential-expression column(s): ",
        paste(missing, collapse = ", "), "."
      )
    )
  }
  for (column in required) {
    column_index <- match(column, lower_names)
    value <- table[, column_index, drop = TRUE]
    if (!is.numeric(value)) {
      lisa_rds_abort(
        "LISA-RDS-012",
        paste0("the edgeR result column `", colnames(table)[column_index],
               "` must be numeric.")
      )
    }
  }
  table
}

lisa_rds_dimensions <- function(object, input_type) {
  target <- object
  if (identical(input_type, "edger")) {
    target <- lisa_edger_result_table(object)
  }
  dim(target)
}

lisa_validate_rds_object <- function(object, input_type, rds_max_bytes) {
  if (inherits(object, "DGEList")) {
    lisa_rds_abort(
      "LISA-RDS-007",
      paste0(
        "the serialized object is a `DGEList`, which contains counts rather than ",
        "a differential-expression result. Supply a DGELRT, DGEExact or TopTags result instead."
      )
    )
  }
  detected_type <- lisa_detect_rds_input_type(object)
  if (is.na(detected_type)) {
    lisa_rds_abort(
      "LISA-RDS-007",
      paste0(
        "the serialized object has unsupported class `",
        paste(class(object), collapse = "/"),
        "`; allowed inputs are data.frame, matrix, DESeqDataSet, DESeqResults, DGELRT, DGEExact and TopTags."
      )
    )
  }
  if (!identical(input_type, detected_type)) {
    lisa_rds_abort(
      "LISA-RDS-007",
      sprintf(
        "`input_type = %s` is incompatible with the serialized object class `%s` (detected type: %s).",
        input_type, paste(class(object), collapse = "/"), detected_type
      )
    )
  }

  if (identical(detected_type, "edger")) {
    lisa_edger_result_table(object)
  }

  dimensions <- suppressWarnings(lisa_rds_dimensions(object, input_type))
  valid_dimensions <- is.numeric(dimensions) && length(dimensions) == 2L &&
    all(is.finite(dimensions)) && all(dimensions >= 1) &&
    all(dimensions == floor(dimensions))
  max_cells <- floor(rds_max_bytes / 8)
  too_many_cells <- valid_dimensions &&
    (max_cells < 1 || dimensions[[1L]] > max_cells ||
       dimensions[[2L]] > max_cells ||
       dimensions[[1L]] > max_cells / dimensions[[2L]])
  if (!valid_dimensions || too_many_cells) {
    observed <- if (is.null(dimensions)) "missing" else paste(dimensions, collapse = " x ")
    lisa_rds_abort(
      "LISA-RDS-008",
      paste0(
        "the serialized object must be a non-empty two-dimensional DE input within the configured size budget; dimensions=",
        observed, "."
      )
    )
  }

  memory_bytes <- as.numeric(utils::object.size(object))
  if (!is.finite(memory_bytes) || memory_bytes > rds_max_bytes) {
    lisa_rds_abort(
      "LISA-RDS-009",
      sprintf(
        "the deserialized object occupies %s bytes, exceeding `rds_max_bytes = %s`.",
        format(memory_bytes, scientific = FALSE),
        format(rds_max_bytes, scientific = FALSE)
      )
    )
  }
  invisible(list(input_type = detected_type, dimensions = dimensions,
                 memory_bytes = memory_bytes))
}

lisa_validate_carried_rds_receipt <- function(receipt) {
  required <- c(
    "schema_version", "source_kind", "source_file", "source_sha256",
    "object_sha256", "input_type",
    "object_class", "n_rows", "n_columns", "serialized_bytes",
    "memory_bytes", "rds_max_bytes", "trusted_rds_authorized",
    "authorization_scope", "path_policy", "trust_assumption"
  )
  numeric_fields <- c(
    "n_rows", "n_columns", "serialized_bytes", "memory_bytes", "rds_max_bytes"
  )
  numeric_valid <- is.list(receipt) && all(numeric_fields %in% names(receipt)) &&
    all(vapply(numeric_fields, function(field) {
      value <- receipt[[field]]
      is.numeric(value) && length(value) == 1L && !is.na(value) &&
        is.finite(value) && value >= 0 && value == floor(value)
    }, logical(1)))
  if (!is.list(receipt) || !all(required %in% names(receipt)) ||
      !identical(receipt$schema_version, "lisa_rds_input_receipt_v2") ||
      !identical(receipt$source_kind, "trusted_rds") ||
      !is.character(receipt$input_type) || length(receipt$input_type) != 1L ||
      !receipt$input_type %in% c("de_table", "deseq2_dds", "deseq2_results", "edger") ||
      !is.character(receipt$source_sha256) || length(receipt$source_sha256) != 1L ||
      is.na(receipt$source_sha256) ||
      !grepl("^[0-9a-f]{64}$", receipt$source_sha256) ||
      !is.character(receipt$object_sha256) ||
      length(receipt$object_sha256) != 1L || is.na(receipt$object_sha256) ||
      !grepl("^[0-9a-f]{64}$", receipt$object_sha256) ||
      !is.character(receipt$source_file) || length(receipt$source_file) != 1L ||
      is.na(receipt$source_file) || !nzchar(receipt$source_file) ||
      grepl("[[:cntrl:]]", receipt$source_file, useBytes = TRUE) ||
      !identical(receipt$source_file, basename(receipt$source_file)) ||
      !is.character(receipt$object_class) || length(receipt$object_class) != 1L ||
      is.na(receipt$object_class) || !nzchar(receipt$object_class) ||
      !numeric_valid || receipt$n_rows < 1 || receipt$n_columns < 1 ||
      receipt$rds_max_bytes < 1 || receipt$serialized_bytes > receipt$rds_max_bytes ||
      receipt$memory_bytes > receipt$rds_max_bytes ||
      !isTRUE(receipt$trusted_rds_authorized) ||
      !identical(receipt$authorization_scope, "explicit_local_R_call") ||
      !identical(receipt$path_policy, "basename_only") ||
      !identical(receipt$trust_assumption, "local_stable_fully_trusted_file")) {
    lisa_rds_abort(
      "LISA-RDS-013",
      "an invalid internal trusted-RDS receipt reached a branched analysis."
    )
  }
  receipt
}

lisa_rds_receipt_table <- function(receipt) {
  receipt <- lisa_validate_carried_rds_receipt(receipt)
  data.frame(
    schema_version = receipt$schema_version,
    source_kind = receipt$source_kind,
    source_file = receipt$source_file,
    source_sha256 = receipt$source_sha256,
    object_sha256 = receipt$object_sha256,
    input_type = receipt$input_type,
    object_class = receipt$object_class,
    n_rows = as.integer(receipt$n_rows),
    n_columns = as.integer(receipt$n_columns),
    serialized_bytes = as.numeric(receipt$serialized_bytes),
    memory_bytes = as.numeric(receipt$memory_bytes),
    rds_max_bytes = as.numeric(receipt$rds_max_bytes),
    trusted_rds_authorized = isTRUE(receipt$trusted_rds_authorized),
    authorization_scope = receipt$authorization_scope,
    path_policy = receipt$path_policy,
    trust_assumption = receipt$trust_assumption,
    stringsAsFactors = FALSE
  )
}

lisa_prepare_de_input <- function(input, input_type, trusted_rds = FALSE,
                                  rds_max_bytes = 512 * 1024^2,
                                  input_receipt = NULL,
                                  .receipt_object_verified = FALSE) {
  trusted_rds <- lisa_validate_trusted_rds(trusted_rds)
  rds_max_bytes <- lisa_validate_rds_max_bytes(rds_max_bytes)
  if (!is.logical(.receipt_object_verified) ||
      length(.receipt_object_verified) != 1L ||
      is.na(.receipt_object_verified)) {
    lisa_rds_abort(
      "LISA-RDS-013",
      "the internal trusted-RDS verification state is invalid."
    )
  }
  if (!is.null(input_receipt)) {
    receipt <- lisa_validate_carried_rds_receipt(input_receipt)
    validation <- tryCatch(
      lisa_validate_rds_object(input, input_type, rds_max_bytes),
      error = function(error) NULL
    )
    observed_class <- paste(class(input), collapse = "/")
    observed_object_sha256 <- if (isTRUE(.receipt_object_verified)) {
      receipt$object_sha256
    } else {
      tryCatch(
        lisa_rds_object_sha256(input),
        error = function(error) NA_character_
      )
    }
    if (!trusted_rds || is.null(validation) ||
        !identical(receipt$input_type, input_type) ||
        !identical(receipt$object_class, observed_class) ||
        !identical(as.integer(receipt$n_rows), as.integer(validation$dimensions[[1L]])) ||
        !identical(as.integer(receipt$n_columns), as.integer(validation$dimensions[[2L]])) ||
        !identical(as.numeric(receipt$memory_bytes), as.numeric(validation$memory_bytes)) ||
        !identical(as.numeric(receipt$rds_max_bytes), as.numeric(rds_max_bytes)) ||
        !identical(receipt$object_sha256, observed_object_sha256)) {
      lisa_rds_abort(
        "LISA-RDS-013",
        "the internal trusted-RDS receipt does not match the carried object, type or size limit."
      )
    }
    return(list(
      input = input,
      input_type = input_type,
      source = "trusted_rds_propagated",
      receipt = receipt
    ))
  }
  if (!lisa_is_rds_path(input)) {
    return(list(input = input, input_type = input_type, source = "direct", receipt = NULL))
  }
  if (!trusted_rds) {
    lisa_rds_abort(
      "LISA-RDS-003",
      paste0(
        "RDS input is disabled by default. Only a local R caller may opt in with ",
        "`run_lisa_de(..., trusted_rds = TRUE, rds_max_bytes = ...)`; YAML/JSON configurations cannot grant this capability."
      )
    )
  }

  source_file <- basename(input)
  if (grepl("[[:cntrl:]]", source_file, useBytes = TRUE)) {
    lisa_rds_abort(
      "LISA-RDS-014",
      "the trusted RDS base file name contains a control character and cannot be recorded safely."
    )
  }

  info <- suppressWarnings(file.info(input))
  if (nrow(info) != 1L || is.na(info$size[[1L]]) ||
      !isTRUE(utils::file_test("-f", input))) {
    lisa_rds_abort(
      "LISA-RDS-004",
      paste0("the RDS path must identify an existing regular file: ", input)
    )
  }
  file_bytes <- as.numeric(info$size[[1L]])
  if (file_bytes > rds_max_bytes) {
    lisa_rds_abort(
      "LISA-RDS-005",
      sprintf(
        "the RDS file is %s bytes, exceeding `rds_max_bytes = %s`; it was not deserialized.",
        format(file_bytes, scientific = FALSE),
        format(rds_max_bytes, scientific = FALSE)
      )
    )
  }

  source_sha256 <- lisa_sha256_file(input)

  object <- tryCatch(
    lisa_read_rds_once(input),
    error = function(error) {
      lisa_rds_abort(
        "LISA-RDS-006",
        paste0("trusted RDS deserialization failed: ", conditionMessage(error))
      )
    }
  )
  detected_type <- lisa_detect_rds_input_type(object)
  resolved_type <- if (identical(input_type, "auto")) detected_type else input_type
  validation <- lisa_validate_rds_object(object, resolved_type, rds_max_bytes)
  object_sha256 <- tryCatch(
    lisa_rds_object_sha256(object),
    error = function(error) {
      lisa_rds_abort(
        "LISA-RDS-015",
        paste0(
          "the trusted RDS object could not be fingerprinted after validation: ",
          conditionMessage(error)
        )
      )
    }
  )
  receipt <- list(
    schema_version = "lisa_rds_input_receipt_v2",
    source_kind = "trusted_rds",
    source_file = source_file,
    source_sha256 = source_sha256,
    object_sha256 = object_sha256,
    input_type = resolved_type,
    object_class = paste(class(object), collapse = "/"),
    n_rows = as.integer(validation$dimensions[[1L]]),
    n_columns = as.integer(validation$dimensions[[2L]]),
    serialized_bytes = file_bytes,
    memory_bytes = validation$memory_bytes,
    rds_max_bytes = rds_max_bytes,
    trusted_rds_authorized = TRUE,
    authorization_scope = "explicit_local_R_call",
    path_policy = "basename_only",
    trust_assumption = "local_stable_fully_trusted_file"
  )
  list(
    input = object,
    input_type = resolved_type,
    source = "trusted_rds",
    receipt = receipt
  )
}
