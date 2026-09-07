# General package utilities.

# Helper for list access
`%||%` <- function(a, b) if (!is.null(a)) a else b

validate_conf_level <- function(conf_level, arg = "conf_level") {
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop("`", arg, "` must be a single finite number between 0 and 1.")
  }

  conf_level
}

has_nested_columns <- function(data) {
  any(vapply(data, function(x) is.list(x) || !is.null(dim(x)), logical(1)))
}

# Bind data frames row-wise while filling columns missing from individual inputs.
# This scoped helper is not a drop-in replacement for dplyr::bind_rows().
bind_rows_fill <- function(x) {
  x <- Filter(Negate(is.null), x)

  if (length(x) == 0) {
    return(data.frame())
  }

  x <- lapply(x, as.data.frame, stringsAsFactors = FALSE, optional = TRUE)
  # Keep matrix-valued predictors intact and preserve base list-column filling.
  if (any(vapply(x, has_nested_columns, logical(1)))) {
    all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))
    x <- lapply(x, function(df) {
      for (col in setdiff(all_cols, names(df))) {
        df[[col]] <- rep(NA, nrow(df))
      }
      df[, all_cols, drop = FALSE]
    })
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    return(out)
  }
  as.data.frame(data.table::rbindlist(x, use.names = TRUE, fill = TRUE))
}

# Left join while preserving left-hand order and repeated-key expansion.
# This scoped helper is not a drop-in replacement for dplyr::left_join().
left_join_preserve_order <- function(x, y, by) {
  if (
    length(by) == 0L ||
      has_nested_columns(x) ||
      has_nested_columns(y) ||
      any(setdiff(names(y), by) %in% names(x))
  ) {
    left <- as.data.frame(x)[, by, drop = FALSE]
    right <- as.data.frame(y)[, by, drop = FALSE]
    key_names <- if (length(by)) paste0("key", seq_along(by)) else character()
    names(left) <- names(right) <- key_names
    left$.left <- seq_len(nrow(x))
    right$.right <- seq_len(nrow(y))
    rows <- if (length(by)) {
      merge(left, right, by = key_names, all.x = TRUE, sort = FALSE)
    } else {
      expand.grid(
        .right = if (nrow(y)) seq_len(nrow(y)) else NA_integer_,
        .left = seq_len(nrow(x))
      )
    }
    rows <- rows[order(rows$.left, rows$.right), , drop = FALSE]
    out <- cbind(
      x[rows$.left, , drop = FALSE],
      y[rows$.right, setdiff(names(y), by), drop = FALSE]
    )
    rownames(out) <- NULL
    return(out)
  }
  out <- merge(
    data.table::as.data.table(x),
    data.table::as.data.table(y),
    by = by,
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
  )
  data.table::setcolorder(out, c(names(x), setdiff(names(y), by)))
  as.data.frame(out)
}

matrix_to_long <- function(
  mat,
  id_name = "id",
  time_name = "time",
  value_name = "value"
) {
  ids <- rownames(mat) %||% seq_len(nrow(mat))
  times <- colnames(mat) %||% seq_len(ncol(mat))

  out <- data.frame(
    rep(ids, each = ncol(mat)),
    rep(times, times = nrow(mat)),
    as.vector(t(mat)),
    check.names = FALSE
  )
  names(out) <- c(id_name, time_name, value_name)
  out
}

named_list_to_wide <- function(x, id = seq_along(x), id_name = "boot_id") {
  value_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- data.frame(id, check.names = FALSE)
  names(out) <- id_name

  for (name in value_names) {
    out[[name]] <- vapply(
      x,
      function(item) {
        if (is.null(item) || !name %in% names(item)) {
          return(NA_real_)
        }
        as.numeric(item[[name]][1])
      },
      numeric(1)
    )
  }

  out
}

pivot_state_columns_long <- function(
  data,
  values_to = "probability",
  names_to = "state",
  names_prefix = "state_"
) {
  state_cols <- grep(paste0("^", names_prefix), names(data), value = TRUE)
  id_cols <- setdiff(names(data), state_cols)

  if (length(state_cols) == 0) {
    stop("No state columns found.")
  }

  row_index <- rep(seq_len(nrow(data)), each = length(state_cols))
  out <- data[row_index, id_cols, drop = FALSE]
  rownames(out) <- NULL
  out[[names_to]] <- sub(
    paste0("^", names_prefix),
    "",
    rep(state_cols, times = nrow(data))
  )
  out[[values_to]] <- as.vector(t(as.matrix(data[, state_cols, drop = FALSE])))
  out
}

reorder_columns <- function(data, first) {
  first <- intersect(first, names(data))
  data[, c(first, setdiff(names(data), first)), drop = FALSE]
}

terms_has_offset <- function(terms) {
  offset_terms <- attr(terms, "offset")
  !is.null(offset_terms) && length(offset_terms) > 0
}

model_uses_offset <- function(model) {
  model_chk <- if (inherits(model, "robcov_vglm")) {
    model$vglm_fit
  } else {
    model
  }

  terms_obj <- tryCatch(stats::terms(model_chk), error = function(e) NULL)
  if (!is.null(terms_obj) && terms_has_offset(terms_obj)) {
    return(TRUE)
  }

  offset_obj <- tryCatch(
    {
      if (inherits(model_chk, "vglm")) {
        methods::slot(model_chk, "offset")
      } else {
        model_chk$offset
      }
    },
    error = function(e) NULL
  )
  offset_is_zero <- FALSE
  if (!is.null(offset_obj) && length(offset_obj) > 0) {
    offset_num <- suppressWarnings(as.numeric(offset_obj))
    if (any(!is.na(offset_num) & offset_num != 0)) {
      return(TRUE)
    }
    offset_is_zero <- all(!is.na(offset_num) & offset_num == 0)
  }

  call_obj <- tryCatch(
    {
      if (inherits(model_chk, "vglm")) {
        methods::slot(model_chk, "call")
      } else {
        model_chk$call
      }
    },
    error = function(e) NULL
  )
  if (!is.null(call_obj)) {
    call_args <- as.list(call_obj)
    if (
      "offset" %in%
        names(call_args) &&
        !is.null(call_args[["offset"]]) &&
        !offset_is_zero
    ) {
      return(TRUE)
    }
  }

  FALSE
}

stop_unsupported_offset <- function() {
  stop(
    "Model offsets are not supported by mostr Markov SOP workflows. ",
    "Please remove offset() terms or the offset argument before fitting the model.",
    call. = FALSE
  )
}

with_sop_fallback_notification_scope <- function(code) {
  old <- getOption("mostr.sop_fallback_notification")
  if (is.environment(old)) {
    return(force(code))
  }
  state <- new.env(parent = emptyenv())
  state$shown <- FALSE
  options(mostr.sop_fallback_notification = state)
  on.exit(options(mostr.sop_fallback_notification = old), add = TRUE)
  force(code)
}

notify_sop_reference_fallback <- function(reason = NULL) {
  state <- getOption("mostr.sop_fallback_notification")
  if (is.environment(state) && isTRUE(state$shown)) {
    return(invisible(FALSE))
  }
  if (is.environment(state)) {
    state$shown <- TRUE
  }

  detail <- if (is.null(reason) || !nzchar(reason)) {
    ""
  } else {
    paste0(" Reason: ", reason)
  }
  message(
    "Compiled C++ SOP calculations were not used; ",
    "falling back to the R implementation.",
    detail
  )
  invisible(TRUE)
}
