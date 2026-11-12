# SOP result shaping helpers.

validate_sops_by <- function(by, data, data_arg = "newdata") {
  if (is.null(by)) {
    return(invisible(NULL))
  }
  if (!is.character(by)) {
    stop("'by' must be a character vector of variable names.")
  }
  missing_vars <- setdiff(by, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Variables specified in 'by' not found in ",
      data_arg,
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }
  invisible(NULL)
}

# The arguments of a SOP call that are not settings: the resolved visits and,
# for blrm, the posterior options. The settings are stored separately, in the
# result's `settings` attribute.
sops_call_args <- function(times, ...) {
  c(list(times = times), list(...))
}

# `settings` is the resolved settings object of the call
# (`resolve_markov_settings()`), read back with `markov_result_settings()`.
set_sops_attrs <- function(
  result,
  class_name,
  model,
  call_args,
  settings,
  newdata_orig = NULL,
  avg_args = NULL,
  extra_attrs = list()
) {
  if (!inherits(settings, "markov_settings")) {
    stop("Internal SOP result settings are invalid.")
  }
  result <- as.data.frame(result)
  attr(result, "model") <- model
  attr(result, "call_args") <- call_args
  attr(result, "settings") <- settings
  attr(result, "newdata_orig") <- newdata_orig
  attr(result, "avg_args") <- avg_args

  for (nm in names(extra_attrs)) {
    if (!is.null(extra_attrs[[nm]])) {
      attr(result, nm) <- extra_attrs[[nm]]
    }
  }

  class(result) <- c(class_name, "data.frame")
  result
}

restore_sops_attrs <- function(result, object) {
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  result <- as.data.frame(result)
  class(result) <- class(object)
  result
}

inference_columns <- function() {
  c(
    "estimate",
    "std.error",
    "statistic",
    "p.value",
    "s.value",
    "conf.low",
    "conf.high"
  )
}

order_estimate_columns <- function(x) {
  original_class <- class(x)
  extra_attrs <- attributes(x)
  extra_attrs <- extra_attrs[setdiff(
    names(extra_attrs),
    c("names", "row.names", "class")
  )]
  x <- as.data.frame(x)
  measures <- intersect(inference_columns(), names(x))
  identifiers <- setdiff(names(x), inference_columns())
  out <- x[, c(identifiers, measures), drop = FALSE]
  for (nm in names(extra_attrs)) {
    attr(out, nm) <- extra_attrs[[nm]]
  }
  class(out) <- original_class
  out
}

# TRUE for a result of `comparisons()` or `avg_comparisons()` whose compared
# variable is not among the columns the fitted model uses: those in its formula
# and the previous-state columns. Setting such a column can't change any
# prediction, so the comparison is exactly 0 (or 1 for a ratio) by
# construction. FALSE whenever the model's columns can't be read.
comparison_variable_unused <- function(x) {
  if (!inherits(x, c("markov_comparisons", "markov_avg_comparisons"))) {
    return(FALSE)
  }
  variable <- names(attr(x, "avg_args")$variables)
  if (length(variable) != 1L) {
    return(FALSE)
  }
  used <- markov_model_variables(attr(x, "model"))
  if (is.null(used)) {
    return(FALSE)
  }
  settings <- markov_result_settings(x)
  !variable %in% c(used, settings$p_var, settings$p2_var)
}

# A comparison that is exactly 0 by construction has no sampling uncertainty
# to report. As marginaleffects does, its standard error, test and interval
# are NA, rather than an SE of 0 with a NaN statistic or a zero-width
# interval. Results of other functions, such as predictions under the same
# `variables`, are returned unchanged.
drop_unused_comparison_uncertainty <- function(x) {
  if (!comparison_variable_unused(x)) {
    return(x)
  }
  for (column in intersect(setdiff(inference_columns(), "estimate"), names(x))) {
    x[[column]] <- rep(NA_real_, nrow(x))
  }
  x
}

add_null_test <- function(x, null, warn = TRUE) {
  if (is.null(null)) {
    return(order_estimate_columns(x))
  }
  if (!is.numeric(null) || length(null) != 1L || !is.finite(null)) {
    stop("`null` must be a single finite numeric value.")
  }
  if (!"std.error" %in% names(x)) {
    if (warn) {
      warning(
        "A null hypothesis was requested, but no standard errors are available.",
        call. = FALSE
      )
    }
    return(order_estimate_columns(x))
  }
  x$statistic <- (x$estimate - null) / x$std.error
  x$p.value <- 2 * stats::pnorm(-abs(x$statistic))
  x$s.value <- -log2(x$p.value)
  attr(x, "null") <- null
  order_estimate_columns(x)
}

with_local_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be a single finite number or NULL.")
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(as.integer(seed))
  force(code)
}

# `by` as in the marginaleffects package: `FALSE` or `NULL` keeps one row per
# profile, `TRUE` averages over all profiles, and column names average within
# groups. Grouping by no columns (`character(0)`) is the overall average.
# Returns whether to average and the grouping columns (`NULL` for the overall
# average).
resolve_by_choice <- function(by) {
  if (is.null(by) || isFALSE(by)) {
    return(list(average = FALSE, by = NULL))
  }
  if (isTRUE(by) || identical(by, character(0))) {
    return(list(average = TRUE, by = NULL))
  }
  if (is.character(by) && length(by) > 0L && !anyNA(by) && all(nzchar(by))) {
    return(list(average = TRUE, by = by))
  }
  stop(
    "`by` must be `FALSE`, `TRUE`, or a character vector of column names.",
    call. = FALSE
  )
}

# Grouping by the compared variable of `avg_comparisons()`. `by` naming that
# variable groups the profiles by the value they were observed with. The
# averaged SOPs behind the comparison keep the value that was set in the
# variable's own column, so the observed value is copied, before it is set,
# into a column with this internal name, and the averages are grouped on the
# copy. `restore_observed_by_names()` gives the copy the variable's name once
# the comparison has replaced the set values with contrasts.
observed_by_column <- function(variable) {
  paste0(".mostr_observed_", variable)
}

# `by` with every variable set by `variables` replaced by its observed copy.
observed_by <- function(by, variables) {
  if (is.null(by)) {
    return(NULL)
  }
  set <- by %in% names(variables)
  by[set] <- observed_by_column(by[set])
  by
}

# Adds the observed copies that `by` asks for to profiles that still hold
# their observed values. A copy that is already present is kept, because
# profiles after `create_counterfactual_data()` hold set values.
add_observed_by_columns <- function(data, variables, by) {
  for (variable in names(variables)) {
    column <- observed_by_column(variable)
    if (
      column %in% by &&
        variable %in% names(data) &&
        !column %in% names(data)
    ) {
      data[[column]] <- data[[variable]]
    }
  }
  data
}

# Names observed copies after their variable, in a table that no longer has
# a column of set values.
restore_observed_by_names <- function(data, variables) {
  if (is.null(data)) {
    return(data)
  }
  for (variable in names(variables)) {
    column <- observed_by_column(variable)
    if (column %in% names(data) && !variable %in% names(data)) {
      names(data)[names(data) == column] <- variable
    }
  }
  data
}

# The user's name for a grouping column, for messages.
observed_by_label <- function(by) {
  sub("^\\.mostr_observed_", "", by)
}

# Refuses a `by` column with a missing value among the profiles that are
# averaged. Such a profile would belong to no group and be left out of every
# average without notice, which changes who is averaged. `profiles` holds the
# averaged profiles, one row each, with observed values: `newdata`, or the
# starting profiles stored by `mostm()`. Only these rows count, so a missing
# value at a later visit of the fitted data does not matter. An internal
# observed copy named in `by` is checked on its variable's observed column,
# and the message uses the user's column name. Columns that are absent are
# left to the checks that report them. Every averaging function resolves its
# profiles in `avg_sops()` or `avg_comparison_setup()`, which call this check
# before predicting, so inference replays and bootstrap refits never reach a
# grouping step with a missing `by` value.
stop_missing_by_values <- function(profiles, by, variables = NULL) {
  if (is.null(by) || is.null(profiles)) {
    return(invisible(NULL))
  }
  profiles <- add_observed_by_columns(profiles, variables, by)
  for (column in intersect(by, names(profiles))) {
    n_missing <- sum(is.na(profiles[[column]]))
    if (n_missing == 0L) {
      next
    }
    n_profiles <- nrow(profiles)
    stop(
      "`by` names `",
      observed_by_label(column),
      "`, which is missing for ",
      n_missing,
      " of the ",
      n_profiles,
      " profiles being averaged. A profile with a missing `by` value belongs ",
      "to no group, so it would be left out of every average. Drop or impute ",
      "those patients first, for example by removing them from `newdata` or ",
      "from the data given to `mostm()`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# `avg_sops()` and `avg_state_time()` return the set value in the column of
# each variable in `variables`, so a `by` group of the observed value would
# need the same column name.
stop_by_names_set_variable <- function(by, variables, fn) {
  clash <- intersect(by, names(variables))
  if (length(clash) == 0L) {
    return(invisible(NULL))
  }
  name <- clash[1L]
  stop(
    "`by` names `",
    name,
    "`, which `variables` also sets. In the result of `",
    fn,
    "`, the `",
    name,
    "` column holds the value that was set, so it cannot also hold groups of ",
    "the value each patient was observed with. The averages for each set ",
    "value are already returned without `by`. To average within groups of ",
    "the observed value, copy the column first, in `newdata` or in the data ",
    "given to `mostm()`, for example ",
    "`newdata$observed_",
    name,
    " <- newdata$",
    name,
    "`, and use `by = \"observed_",
    name,
    "\"`.",
    call. = FALSE
  )
}

create_counterfactual_grid <- function(variables) {
  if (length(variables)) {
    do.call(expand.grid, variables)
  } else {
    data.frame(row.names = 1L)
  }
}

#' Create Counterfactual Datasets for G-Computation
#'
#' Creates copies of baseline data with treatment variable set to each level.
#'
#' @param baseline_data Data frame with one row per patient.
#' @param grid Data frame of variable combinations.
#' @param variables Named list of variable values.
#' @param by Optional grouping columns. For each internal observed-value
#'   column that `by` names, the observed value of its variable is copied
#'   before the variable is set.
#'
#' @return Data frame with counterfactual data stacked.
#'
#' @keywords internal
create_counterfactual_data <- function(
  baseline_data,
  grid,
  variables,
  by = NULL
) {
  baseline_data <- add_observed_by_columns(baseline_data, variables, by)
  cf_data_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt_copy <- baseline_data
    for (v in names(grid)) {
      dt_copy[[v]] <- grid[i, v]
    }
    cf_data_list[[i]] <- dt_copy
  }
  do.call(rbind, cf_data_list)
}

# Stacks one copy of the profiles for each combination of `variables`, with
# those columns set, as marginaleffects' predictions() does: `rowid` numbers
# the stacked rows and `rowidcf` records each row's profile.
stack_prediction_profiles <- function(profiles, variables) {
  missing_vars <- setdiff(names(variables), names(profiles))
  if (length(missing_vars) > 0L) {
    stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }
  grid <- create_counterfactual_grid(variables)
  stacked <- create_counterfactual_data(profiles, grid, variables)
  stacked$rowidcf <- rep(profiles$rowid, times = nrow(grid))
  rownames(stacked) <- NULL
  ensure_markov_rowid(stacked)
}

# The distinct patients among prediction rows, one row each: the first copy
# of each profile when `variables` stacked copies of them.
prediction_weight_profiles <- function(data) {
  if (!"rowidcf" %in% names(data)) {
    return(data)
  }
  data[!duplicated(data$rowidcf), , drop = FALSE]
}

validate_sops_weights <- function(weights, n, arg = "weights") {
  if (is.null(weights)) {
    return(NULL)
  }
  if (
    !is.numeric(weights) ||
      length(weights) != n ||
      anyNA(weights) ||
      any(weights < 0)
  ) {
    stop(
      "`",
      arg,
      "` must be a non-negative numeric vector of length ",
      n,
      "."
    )
  }
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("`", arg, "` must have a positive finite sum.")
  }
  weights
}

weighted_sop_mean <- function(x, weights) {
  keep <- !is.na(x)
  if (!any(keep)) {
    return(NA_real_)
  }

  x <- x[keep]
  weights <- weights[keep]
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    return(NA_real_)
  }

  sum(x * weights) / weight_sum
}

aggregate_sops_estimates <- function(result, group_cols, weight_col = NULL) {
  .datatable.aware <- TRUE
  result <- as.data.frame(result)
  missing_vars <- setdiff(group_cols, names(result))
  if (length(missing_vars) > 0) {
    stop("Grouping variables missing: ", paste(missing_vars, collapse = ", "))
  }

  if (is.null(weight_col)) {
    return(aggregate_value(result, "estimate", group_cols, mean))
  }

  if (!weight_col %in% names(result)) {
    stop("Weight column '", weight_col, "' not found in result.")
  }

  work <- data.table::as.data.table(
    result[, c(group_cols, "estimate", weight_col), drop = FALSE]
  )
  work <- work[stats::complete.cases(result[, group_cols, drop = FALSE])]
  out <- as.data.frame(work[,
    list(estimate = weighted_sop_mean(estimate, .SD[[1L]])),
    by = group_cols,
    .SDcols = weight_col
  ])
  rows <- do.call(order, lapply(out[rev(group_cols)], as.factor))
  out <- out[rows, , drop = FALSE]
  rownames(out) <- NULL
  out
}

validate_prediction_weight_ids <- function(prediction_data, id_var) {
  if (is.null(id_var) || !id_var %in% names(prediction_data)) {
    stop(
      "Prediction data must contain unique `id_var` values before bootstrap ",
      "draw weights can be attached."
    )
  }

  ids <- as.character(prediction_data[[id_var]])
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(
      "Prediction data must contain exactly one non-missing row per `id_var` ",
      "before bootstrap draw weights can be attached."
    )
  }

  invisible(ids)
}

validate_draw_weight_column_available <- function(newdata, weight_col) {
  if (is.null(weight_col)) {
    return(invisible(NULL))
  }
  if (weight_col %in% names(newdata)) {
    stop(
      "Prediction data already contains a column named `",
      weight_col,
      "`, which is reserved for bootstrap draw weights. Rename this column ",
      "before calling inferences().",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Marginalize SOPs Array Over Patients
#'
#' Averages individual-level SOPs to get population-average (marginal) SOPs.
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param grid Data frame of variable combinations.
#' @param times Vector of time points.
#' @param y_levels Vector of state levels.
#' @param variables Named list of variables.
#' @param n_cf Number of counterfactual scenarios.
#' @param n_each Number of patients per scenario.
#' @param weights Optional patient-level weights for one set of patients.
#' @param by Optional character vector of variables to aggregate by.
#' @param newdata Counterfactual prediction data, required when `by` is used.
#'
#' @return Data frame with marginalized SOPs.
#'
#' @keywords internal
marginalize_sops_array <- function(
  sops_array,
  grid,
  times,
  y_levels,
  variables,
  n_cf,
  n_each,
  weights = NULL,
  by = NULL,
  newdata = NULL
) {
  n_times <- length(times)
  n_states <- length(y_levels)
  weights <- validate_sops_weights(weights, n_each)

  if (!is.null(by)) {
    if (is.null(newdata) || !is.data.frame(newdata)) {
      stop("`newdata` must be supplied when `by` is used.")
    }
    if (nrow(newdata) != n_cf * n_each) {
      stop("`newdata` is not aligned with the counterfactual SOP array.")
    }
    validate_sops_by(by, newdata)

    expanded_weights <- if (is.null(weights)) {
      NULL
    } else {
      rep(weights, times = n_cf)
    }
    weight_col <- if (is.null(expanded_weights)) NULL else ".mostr_weight"

    result <- array_to_df_individual(
      sops_array = sops_array,
      times = times,
      y_levels = y_levels,
      newdata = newdata,
      weights = expanded_weights,
      weight_col = weight_col
    )
    group_cols <- unique(c("time", "state", names(variables), by))
    return(aggregate_sops_estimates(
      result,
      group_cols,
      weight_col = weight_col
    ))
  }

  # Average within each counterfactual group
  avg_sops_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    start_idx <- (cf_i - 1) * n_each + 1
    end_idx <- cf_i * n_each

    sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
    if (is.null(weights)) {
      avg_sops_mat <- apply(sops_cf, c(2, 3), mean)
    } else {
      weights_norm <- weights / sum(weights)
      sops_cf_mat <- matrix(
        aperm(sops_cf, c(2, 3, 1)),
        nrow = n_times * n_states,
        ncol = n_each
      )
      avg_sops_mat <- matrix(
        as.vector(sops_cf_mat %*% weights_norm),
        nrow = n_times,
        ncol = n_states
      )
    }

    avg_sops_list[[cf_i]] <- avg_sops_mat
  }

  # Format as data frame
  result_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    avg_sops_mat <- avg_sops_list[[cf_i]]

    df <- expand.grid(time = times, state = y_levels)
    df$estimate <- as.vector(avg_sops_mat)

    for (v in names(grid)) {
      df[[v]] <- grid[cf_i, v]
    }

    result_list[[cf_i]] <- df
  }

  bind_rows_fill(result_list)
}


#' Convert SOPs Array to Individual-Level Data Frame
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param times Vector of time points.
#' @param y_levels Vector of state levels.
#' @param newdata Original data with rowid.
#' @param weights Optional patient-level weights aligned to `newdata`.
#' @param weight_col Optional output column name for draw-specific weights.
#'
#' @return Data frame with one row per profile, time, and state.
#'
#' @keywords internal
array_to_df_individual <- function(
  sops_array,
  times,
  y_levels,
  newdata,
  weights = NULL,
  weight_col = NULL
) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]
  weights <- validate_sops_weights(weights, n_pat)

  # Flatten array
  probs_flat <- as.vector(sops_array)

  # Construct indices
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(y_levels, each = n_pat * n_times)

  # Build result by repeating newdata rows
  result <- newdata[idx_pat, , drop = FALSE]

  # Add rowid if it exists
  if ("rowid" %in% names(newdata)) {
    result$rowid <- newdata$rowid[idx_pat]
  } else {
    result$rowid <- idx_pat
  }

  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  if (!is.null(weights)) {
    weight_col <- weight_col %||% ".mostr_weight"
    validate_draw_weight_column_available(newdata, weight_col)
    result[[weight_col]] <- weights[idx_pat]
  }
  rownames(result) <- NULL

  result
}
