# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

validate_lisa_de_index <- function(de_index, base_dir = getwd(), require_files = TRUE) {
  if (is.character(de_index) && length(de_index) == 1) {
    de_index <- read_lisa_tsv(de_index)
  }

  lisa_require_columns(
    de_index,
    c("analysis_id", "de_path"),
    label = "de_index"
  )
  invisible(vapply(de_index$analysis_id, lisa_safe_id, character(1), field = "analysis_id"))

  if (anyDuplicated(de_index$analysis_id)) {
    duplicated_ids <- unique(de_index$analysis_id[duplicated(de_index$analysis_id)])
    stop("de_index has duplicated analysis_id value(s): ", paste(duplicated_ids, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(tolower(de_index$analysis_id))) {
    stop("analysis_id values must also be distinct on case-insensitive filesystems.", call. = FALSE)
  }

  rds_rows <- which(
    !is.na(de_index$de_path) &
      tolower(tools::file_ext(as.character(de_index$de_path))) == "rds"
  )
  if (length(rds_rows)) {
    lisa_rds_abort(
      "LISA-RDS-011",
      paste0(
        "YAML/JSON and de_index workflows accept text DE tables, not RDS files ",
        "(row(s): ", paste(rds_rows, collapse = ", "), "). Use a TSV/CSV input, ",
        "or make an explicit local call to `run_lisa_de(..., trusted_rds = TRUE, rds_max_bytes = ...)`."
      )
    )
  }

  if (isTRUE(require_files)) {
    missing_paths <- character()
    for (path in de_index$de_path) {
      resolved <- lisa_norm_path(path, base_dir)
      if (!file.exists(resolved)) {
        missing_paths <- c(missing_paths, resolved)
      }
    }
    if (length(missing_paths) > 0) {
      stop("de_index points to missing DE table(s): ", paste(missing_paths, collapse = "; "), call. = FALSE)
    }
  }

  invisible(de_index)
}

lisa_config_input_readiness <- function(de_index, check_files = TRUE) {
  path_columns <- grep(
    "(^path$|_path$|_file$|_excel$|_rds$)",
    names(de_index), value = TRUE
  )
  rows <- lapply(seq_len(nrow(de_index)), function(i) {
    analysis_id <- as.character(de_index$analysis_id[[i]])
    lapply(path_columns, function(column) {
      path <- as.character(de_index[[column]][[i]])
      if (is.na(path) || !nzchar(path)) return(NULL)
      ready <- if (isTRUE(check_files)) {
        file.exists(path) && !isTRUE(file.info(path)$isdir)
      } else {
        NA
      }
      data.frame(
        analysis_id = analysis_id,
        input = column,
        path = path,
        checked = isTRUE(check_files),
        ready = ready,
        status = if (!isTRUE(check_files)) {
          "not_checked"
        } else if (isTRUE(ready)) {
          "ready"
        } else {
          "not_ready"
        },
        message = if (isTRUE(check_files) && !isTRUE(ready)) {
          "Referenced input file is absent or is not a regular file."
        } else {
          ""
        },
        stringsAsFactors = FALSE
      )
    })
  })
  rows <- unlist(rows, recursive = FALSE)
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      analysis_id = character(), input = character(), path = character(),
      checked = logical(), ready = logical(), status = character(),
      message = character(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

lisa_config_source_data_readiness <- function(cfg, config_dir,
                                              check_files = TRUE) {
  rows <- lisa_config_get(cfg, "source_data", list())
  if (!length(rows)) {
    return(data.frame(
      analysis_id = character(), input = character(), path = character(),
      checked = logical(), ready = logical(), status = character(),
      message = character(), stringsAsFactors = FALSE
    ))
  }
  if (!isTRUE(check_files)) {
    return(do.call(rbind, lapply(seq_along(rows), function(index) {
      row <- rows[[index]]
      value <- if (is.list(row)) {
        lisa_config_get(row, "path", lisa_config_get(row, "source_path", ""))
      } else {
        ""
      }
      data.frame(
        analysis_id = "config",
        input = sprintf("source_data[%d]", index),
        path = as.character(value)[[1L]],
        checked = FALSE,
        ready = NA,
        status = "not_checked",
        message = "",
        stringsAsFactors = FALSE
      )
    })))
  }

  prepared <- tryCatch(
    lisa_prepare_config_source_data(cfg, config_dir),
    error = identity
  )
  if (inherits(prepared, "error")) {
    return(data.frame(
      analysis_id = "config",
      input = "source_data",
      path = as.character(prepared$source_path %||% ""),
      checked = TRUE,
      ready = FALSE,
      status = "not_ready",
      message = conditionMessage(prepared),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, lapply(seq_along(prepared), function(index) {
    data.frame(
      analysis_id = "config",
      input = sprintf("source_data[%d]", index),
      path = prepared[[index]]$source_path,
      checked = TRUE,
      ready = TRUE,
      status = "ready",
      message = "",
      stringsAsFactors = FALSE
    )
  }))
  row.names(out) <- NULL
  out
}

lisa_duplicate_policy_input_readiness <- function(policies, config_dir,
                                                  check_files = TRUE) {
  policy_names <- names(policies)[vapply(policies, function(policy) {
    identical(lisa_policy_type(policy), "mapping_file")
  }, logical(1))]
  if (!length(policy_names)) {
    return(data.frame(
      analysis_id = character(), input = character(), path = character(),
      checked = logical(), ready = logical(), status = character(),
      message = character(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, lapply(policy_names, function(policy_name) {
    path <- lisa_config_path(
      policies[[policy_name]]$mapping_file, config_dir
    )
    ready <- if (isTRUE(check_files)) {
      file.exists(path) && !isTRUE(file.info(path)$isdir)
    } else {
      NA
    }
    data.frame(
      analysis_id = "config",
      input = paste0("duplicate_policies.", policy_name, ".mapping_file"),
      path = path,
      checked = isTRUE(check_files),
      ready = ready,
      status = if (!isTRUE(check_files)) {
        "not_checked"
      } else if (isTRUE(ready)) {
        "ready"
      } else {
        "not_ready"
      },
      message = if (isTRUE(check_files) && !isTRUE(ready)) {
        "Referenced duplicate-resolution mapping file is absent or is not a regular file."
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
  }))
  row.names(out) <- NULL
  out
}

validate_lisa_contrast_index <- function(contrast_index, de_index = NULL) {
  if (is.null(contrast_index)) {
    return(invisible(data.frame()))
  }
  if (is.character(contrast_index) && length(contrast_index) == 1) {
    contrast_index <- read_lisa_tsv(contrast_index)
  }
  if (nrow(contrast_index) == 0) {
    return(invisible(contrast_index))
  }

  lisa_require_columns(
    contrast_index,
    c("contrast_id", "analysis_a", "analysis_b"),
    label = "contrast_index"
  )
  invisible(vapply(contrast_index$contrast_id, lisa_safe_id, character(1), field = "contrast_id"))
  if ("output_id" %in% names(contrast_index)) {
    present <- !is.na(contrast_index$output_id) & nzchar(as.character(contrast_index$output_id))
    invisible(vapply(contrast_index$output_id[present], lisa_safe_id, character(1), field = "output_id"))
  }

  if (anyDuplicated(contrast_index$contrast_id)) {
    duplicated_ids <- unique(contrast_index$contrast_id[duplicated(contrast_index$contrast_id)])
    stop("contrast_index has duplicated contrast_id value(s): ", paste(duplicated_ids, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(tolower(contrast_index$contrast_id))) {
    stop("contrast_id values must also be distinct on case-insensitive filesystems.", call. = FALSE)
  }

  if (!is.null(de_index)) {
    known <- de_index$analysis_id
    used <- unique(c(contrast_index$analysis_a, contrast_index$analysis_b))
    missing <- setdiff(used, known)
    if (length(missing) > 0) {
      stop("contrast_index references unknown analysis_id value(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  invisible(contrast_index)
}
