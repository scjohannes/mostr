# Average time in states and inference through the average-SOP engines.

#' Calculate Average Total Time in States
#'
#' Estimates the average total time spent in each state over the requested
#' follow-up period, including time accumulated across repeated visits to a
#' state. With `variables = NULL`, averages over patients using their observed
#' starting states and covariates. Supply `variables` to calculate averages
#' under specified counterfactual scenarios.
#'
#' @inheritParams avg_sops
#' @inheritParams avg_comparisons
#' @param variables Optional named list of covariates to set before averaging,
#'   for example `list(tx = c(0, 1))` for average times under each treatment.
#'   The default, `NULL`, uses each patient's observed covariates.
#' @param state_sets States to summarize. `NULL` returns every state separately.
#'   A vector pools its states into one total; a named list returns one total
#'   per named state set, for example
#'   `list(home = "1", hospital = c("2", "3"))`. Names must be unique.
#' @param times Required visits to include, each listed once. Numeric time
#'   variables use numeric values; factor-valued visit indices use fitted visit
#'   levels.
#' @param time_map Optional named numeric vector or data frame giving the
#'   elapsed time of each visit label, for example
#'   `c("1" = 3, "2" = 7, "3" = 14)` for visits on days 3, 7, and 14. Supplying
#'   it measures time on this elapsed-time scale instead of counting visits.
#' @param baseline_time Elapsed time at which patients are in their starting
#'   state, used only with `time_map`. The default `0` means that totals are
#'   counted from time 0, when each patient is in their starting state. Set it
#'   to `NULL` to count from the first mapped visit instead. Supplying a value
#'   without `time_map` is an error.
#' @param target_times Optional numeric vector of elapsed times, requiring
#'   `time_map`. Probabilities are linearly interpolated at these times and
#'   totals are the area under the resulting curve (trapezoidal rule), so the
#'   totals cover the period from the first to the last of these times. When
#'   omitted, the times are `baseline_time` followed by the mapped visit times,
#'   covering the whole period from the starting state to the last visit; with
#'   `baseline_time = NULL`, they are the mapped visit times only. Supply, for
#'   example, `target_times = 7:28` to total only days 7 to 28.
#' @param posterior_summary Summary of total-time posterior draws for `blrm`:
#'   `"mean"` or `"median"`.
#' @param return_draws Logical. For `blrm`, retain posterior time-total draws
#'   for [get_draws()]. Frequentist draws are retained through [inferences()].
#'
#' @return A data frame of class `markov_avg_time`, with `state_set`, any
#'   counterfactual and grouping columns, and `estimate` (average total time).
#'   Rows follow the order of `state_sets`. If supplied, `time_unit` is
#'   included as a label; it does not rescale time. For `blrm`, posterior
#'   `std.error`, `conf.low`, and `conf.high` are returned directly. For
#'   frequentist models, add these columns with [inferences()].
#'
#' @details
#' Without `time_map`, each requested visit counts as one unit of time: the
#' total for a state is the sum of its occupancy probabilities over the
#' requested visits. For example, with `times = 1:30` and daily visits, the
#' total is the expected number of days spent in the state.
#'
#' With `time_map`, visits are placed on the elapsed-time scale. Occupancy
#' probabilities are joined by straight lines between the mapped times, and
#' the total is the area under that curve. By default the curve starts at
#' `baseline_time`, where each patient is in their observed starting state, so
#' the interval before the first visit is included. For example, with visits
#' on days 3, 7, 14, and 28 and `baseline_time = 0`, totals cover days 0 to 28.
#'
#' `inferences(method = "delta")` provides analytical standard errors and
#' confidence intervals of the form estimate plus or minus a normal quantile
#' times the standard error. These intervals are not restricted to the
#' possible range of time and can, for example, extend below zero. Its default
#' `vcov = "unconditional"` accounts for coefficient estimation and for which
#' patients were sampled, including their observed starting states. Use
#' `vcov = "conditional"` to treat the patients' starting states and
#' covariates as given and account only for coefficient estimation. Supplied
#' `newdata` requires conditional variance. When `variables` sets the
#' previous-state variable (a counterfactual starting state), the baseline
#' starts in that state and only conditional variance is available
#' (`method = "delta", vcov = "conditional"` or `method = "mvn"`). The
#' analytical method supports first-order full proportional-odds frequentist
#' models; analytical inference with `by` is not supported. See [inferences()]
#' for additional model restrictions.
#'
#' Extract the analytical covariance matrix with `vcov(result)`. It contains a
#' covariance for every pair of rows, including different state sets and
#' counterfactual scenarios. Its diagonal entries are variances; their square
#' roots are the reported standard errors.
#'
#' The simulation and resampling methods of [inferences()] differ in which
#' sources of uncertainty they include:
#' * `method = "mvn"` draws model coefficients from their estimated
#'   multivariate normal distribution and keeps the patients, and therefore
#'   the distribution of starting states, fixed. It is comparable to
#'   conditional variance, so its intervals are typically narrower than the
#'   default analytical intervals.
#' * `method = "score_bootstrap"`, `method = "bootstrap"`, and `method = "fwb"`
#'   (fractional weighted bootstrap) re-draw or re-weight patients, including
#'   their starting states. They are comparable to unconditional variance.
#'
#' Each of these methods, and posterior summaries for `blrm` models, first
#' converts every draw of occupancy probabilities into time totals and then
#' summarizes the totals, rather than adding up separately summarized
#' probabilities. Use [get_draws()] to inspect the draws when they have been
#' retained with `return_draws = TRUE`.
#'
#' @seealso [avg_sops()], [avg_comparisons()], [time_in_state()], [inferences()]
#' @examples
#' \dontrun{
#' # Average over the fitted patients' observed starting states and covariates.
#' avg_time(fit, times = 1:30) |>
#'   inferences(method = "delta")
#'
#' avg_time(
#'   fit, variables = list(tx = c(0, 1)), times = 1:30,
#'   state_sets = list(home = "1", hospital = c("2", "3"))
#' ) |>
#'   inferences(method = "delta", vcov = "conditional")
#'
#' # Visits on days 3, 7, 14, and 28; totals cover days 0 to 28.
#' avg_time(
#'   fit, times = 1:4, time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
#'   time_unit = "days"
#' ) |>
#'   inferences(method = "fwb", n_draws = 500)
#' }
#' @export
avg_time <- function(
  model,
  newdata = NULL,
  variables = NULL,
  state_sets = NULL,
  by = NULL,
  times,
  y_levels = NULL,
  absorb = NULL,
  time_map = NULL,
  baseline_time = 0,
  target_times = NULL,
  time_unit = NULL,
  id_var = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied to `avg_time()`.")
  }
  if (anyDuplicated(times)) {
    stop(
      "`times` must not contain duplicate visits; each visit is counted once."
    )
  }
  if (!missing(baseline_time) && !is.null(baseline_time) && is.null(time_map)) {
    stop(
      "`baseline_time` requires `time_map`; without `time_map`, totals sum ",
      "occupancy probabilities at the requested visits."
    )
  }
  if (!is.null(target_times) && is.null(time_map)) {
    stop("`target_times` requires `time_map` for real-time integration.")
  }
  if (!is.null(target_times) && !length(target_times)) {
    stop("`target_times` must contain at least one time.")
  }
  if (
    !is.null(time_unit) &&
      (!is.character(time_unit) || length(time_unit) != 1L || is.na(time_unit))
  ) {
    stop("`time_unit` must be a single non-missing character label or NULL.")
  }
  posterior_summary <- match.arg(posterior_summary)
  reserved <- intersect(
    c(names(variables), by),
    c("state_set", if (!is.null(time_unit)) "time_unit")
  )
  if (length(reserved)) {
    stop(
      "Grouping or counterfactual variables use reserved output names: ",
      paste(reserved, collapse = ", "),
      ". Rename these variables."
    )
  }
  if (!is.null(state_sets)) {
    validate_avg_time_state_sets(state_sets)
  }
  extra_args <- list(...)
  avg <- avg_comparison_replay_avg_sops(
    model = model,
    newdata = newdata,
    variables = variables,
    by = by,
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    id_var = id_var,
    p2_var = p2_var,
    gap_var = gap_var,
    time_covariates = time_covariates,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = inherits(model, "blrm"),
    extra_args = extra_args
  )
  state_sets <- normalize_comparison_state_sets(
    state_sets,
    attr(avg, "y_levels")
  )
  args <- list(
    state_sets = state_sets,
    time_map = time_map,
    baseline_time = baseline_time,
    target_times = target_times,
    time_unit = time_unit,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    extra_args = extra_args
  )
  result <- avg_time_from_sops(avg, args, return_draws)
  set_avg_time_attrs(result, avg, args)
}

# Structural `state_sets` checks that do not need the model's state levels, so
# they can run before predictions. Names mirror
# `normalize_comparison_state_sets()`, which also checks the states themselves.
validate_avg_time_state_sets <- function(state_sets) {
  sets <- if (is.list(state_sets)) state_sets else list(state_sets)
  labels <- names(sets)
  if (!is.list(state_sets) || is.null(labels)) {
    labels <- vapply(
      sets,
      function(states) state_set_label(as.character(states)),
      character(1)
    )
  }
  if (
    !length(sets) ||
      any(lengths(sets) == 0L) ||
      anyNA(unlist(sets, use.names = FALSE)) ||
      anyNA(labels) ||
      any(!nzchar(labels)) ||
      anyDuplicated(labels)
  ) {
    stop(
      "`state_sets` must contain nonempty state sets with unique, non-missing names."
    )
  }
  invisible(state_sets)
}

replay_avg_time_sops <- function(object, conf_level, return_draws = FALSE) {
  args <- attr(object, "time_args")
  avg_args <- attr(object, "avg_args")
  if (is.null(args) || is.null(avg_args)) {
    stop("Stored average-time arguments are required for inference.")
  }
  avg_comparison_replay_avg_sops(
    model = attr(object, "model"),
    newdata = if (isTRUE(attr(object, "newdata_supplied"))) {
      attr(object, "newdata_orig")
    } else {
      NULL
    },
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    y_levels = attr(object, "y_levels"),
    absorb = attr(object, "absorb"),
    time_var = attr(object, "time_var"),
    p_var = attr(object, "p_var"),
    id_var = avg_args$id_var,
    p2_var = attr(object, "p2_var"),
    gap_var = attr(object, "gap_var"),
    time_covariates = attr(object, "time_covariates"),
    include_re = args$include_re,
    n_draws = args$n_draws,
    seed = args$seed,
    posterior_summary = args$posterior_summary,
    conf_level = conf_level,
    return_draws = return_draws,
    extra_args = args$extra_args
  )
}

reduce_avg_time_df <- function(
  data,
  state_sets,
  group_cols,
  value_col,
  real_time
) {
  group_cols <- unique(c(intersect("draw_id", names(data)), group_cols))
  pieces <- lapply(names(state_sets), function(label) {
    out <- reduce_state_time_df(
      data,
      state_sets[[label]],
      value_col,
      group_cols,
      real_time
    )
    out$state_set <- label
    out
  })
  reorder_columns(bind_rows_fill(pieces), c("state_set", group_cols, value_col))
}

avg_time_from_sops <- function(avg, args, return_draws) {
  if (!is.null(args$time_map)) {
    avg <- interpolate_sops(
      avg,
      time_map = args$time_map,
      baseline_time = args$baseline_time,
      target_times = comparison_real_time_target_times(
        avg,
        args$time_map,
        args$target_times,
        args$baseline_time
      )
    )
  }
  avg_args <- attr(avg, "avg_args")
  group_cols <- unique(c(names(avg_args$variables), avg_args$by))
  real_time <- !is.null(args$time_map)
  result <- reduce_avg_time_df(
    avg,
    args$state_sets,
    group_cols,
    "estimate",
    real_time
  )
  draws <- attr(avg, "draws")
  if (!is.null(draws)) {
    value_col <- sop_draw_value_col(draws)
    reduced <- reduce_avg_time_df(
      draws,
      args$state_sets,
      group_cols,
      value_col,
      real_time
    )
    names(reduced)[names(reduced) == value_col] <- "estimate"
    keys <- c("state_set", group_cols)
    if (identical(attr(avg, "method"), "posterior")) {
      summary <- summarize_comparison_draws(
        reduced,
        keys,
        args$posterior_summary,
        attr(avg, "conf_level")
      )
      # Keep the row order of the point reduction.
      result <- left_join_preserve_order(
        result[, keys, drop = FALSE],
        summary,
        by = keys
      )
    } else {
      ci <- compute_ci_from_draws(
        reduced,
        keys,
        conf_level = attr(avg, "conf_level"),
        conf_type = attr(avg, "conf_type"),
        point_estimates = result
      )
      result <- left_join_preserve_order(result, ci, by = keys)
    }
    if (!is.null(args$time_unit)) {
      reduced$time_unit <- args$time_unit
    }
    if (isTRUE(return_draws)) attr(result, "draws") <- reduced
  }
  if (!is.null(args$time_unit)) {
    result$time_unit <- args$time_unit
  }
  order_estimate_columns(result)
}

set_avg_time_attrs <- function(result, avg, args) {
  for (name in setdiff(
    names(attributes(avg)),
    c(
      "names",
      "row.names",
      "class",
      "draws",
      "analytical",
      "baseline_anchor_draws"
    )
  )) {
    attr(result, name) <- attr(avg, name)
  }
  attr(result, "time_args") <- args
  class(result) <- c("markov_avg_time", "data.frame")
  result
}

inferences_avg_time <- function(
  object,
  method,
  n_draws,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  update_datadist,
  use_coefstart
) {
  avg <- replay_avg_time_sops(object, conf_level)
  avg <- inferences(
    avg,
    method = method,
    n_draws = n_draws,
    vcov = vcov,
    cluster = cluster,
    workers = workers,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = TRUE,
    update_datadist = update_datadist,
    use_coefstart = use_coefstart
  )
  args <- attr(object, "time_args")
  result <- avg_time_from_sops(avg, args, return_draws)
  set_avg_time_attrs(result, avg, args)
}
