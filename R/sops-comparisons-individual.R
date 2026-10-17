# SOP comparisons for each prediction profile.

#' Comparisons From Markov State Occupancy Probabilities for Each Profile
#'
#' Computes a separate contrast for every prediction profile. Each profile, by
#' default each fitted patient's starting state and covariates, is predicted
#' once with the variable in `variables` set to each requested value. Every
#' later value is then compared with the first value for that same profile.
#' [avg_comparisons()] averages these contrasts over the profiles.
#'
#' @inheritParams avg_comparisons
#' @param newdata Optional data frame of prediction profiles. Each row is one
#'   profile with its own starting state and covariates, and receives its own
#'   comparisons. If `NULL`, uses the starting profiles stored by
#'   [orm_markov()], [blrm_markov()], or [vglm_markov()], one per fitted
#'   patient.
#' @param variables A named list with one variable and at least two values,
#'   for example `list(tx = c(0, 1))`. Each profile is predicted with the
#'   variable set to each value. The first value is the reference level.
#' @param by `FALSE` (the default) returns one comparison per profile. `TRUE`
#'   averages the comparisons over all profiles, and a character vector of
#'   column names averages them within the groups those columns define. Both
#'   return the result of [avg_comparisons()] with the same arguments. For
#'   `comparison = "ratio"`, that average is the ratio of the averaged
#'   summaries, not the average of the profiles' ratios.
#' @param return_draws Logical. For `blrm`, store the posterior draws of each
#'   comparison for extraction with [get_draws()].
#'
#' @return A data frame of class `markov_comparisons` with one row per profile,
#'   comparison level, state set, and (for `estimand = "sop"`) visit:
#'   \item{rowid}{Profile number: the row of `newdata`, or the position of the
#'     patient among the stored starting profiles.}
#'   \item{estimand}{`"sop"`, `"time_in_state"`, or `"time_benefit"`.}
#'   \item{time}{Visit, for `estimand = "sop"`.}
#'   \item{state_set}{Name of the state set. Not used for time benefit.}
#'   \item{term, reference_level, comparison_level}{The compared variable, its
#'     reference value, and the value compared with it.}
#'   \item{contrast, comparison}{A label such as `"1 - 0"`, and whether the
#'     contrast is a `"difference"` or a `"ratio"`.}
#'   \item{(profile columns)}{The profile's covariates and starting state as
#'     observed. The compared variable keeps the profile's observed value, so
#'     results can be summarized, for example, among treated patients.}
#'   \item{estimate}{The comparison for this profile.}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     summaries of the comparison.}
#'
#' @details
#' For one profile, write \eqn{S(A, t \mid z)} for the probability that the
#' profile occupies a state in the state set \eqn{A} at visit \eqn{t} when the
#' compared variable is set to \eqn{z}. With reference value \eqn{z_0} and
#' comparison value \eqn{z_1}, `estimand = "sop"` returns
#' \eqn{S(A, t \mid z_1) - S(A, t \mid z_0)}, or their ratio. Time in state
#' adds these probabilities over the requested visits; with `time_map`, it
#' integrates them over real time from `baseline_time`, starting from the
#' profile's starting state, or from the state set in `variables` when it sets
#' the previous state. Time benefit compares the profile's two state
#' distributions at each visit, as described in [avg_comparisons()].
#'
#' For differences and time benefit, the mean of these comparisons over the
#' profiles equals the [avg_comparisons()] estimate. For ratios it does not,
#' because [avg_comparisons()] divides averaged probabilities.
#'
#' For frequentist models, pass the result to [inferences()]. Because each
#' profile's starting state and covariates are held as given, every method
#' reflects uncertainty in the model coefficients only:
#' \itemize{
#'   \item `method = "delta"` uses conditional variance. For one profile, the
#'     derivative of a difference with respect to the coefficients is the
#'     derivative of the comparison scenario's summary minus that of the
#'     reference scenario, and its variance is \eqn{G V_\beta G^\top}, where
#'     \eqn{G} holds these derivatives and \eqn{V_\beta} is the coefficient
#'     covariance. It supports differences of SOPs and of time in states for
#'     first-order full proportional-odds `orm`, `vglm`, and `robcov_vglm`
#'     fits; [stats::vcov()] then returns the covariance between rows.
#'   \item `method = "mvn"`, `"score_bootstrap"`, or `"fwb"` predict every
#'     profile under each compared value within each coefficient draw or
#'     weighted refit, and compute the comparisons within that draw before
#'     forming the intervals. For the stored profiles, score-bootstrap and FWB
#'     draws also carry each patient's weight in that draw (`score_weight` or
#'     `fwb_weight`); averaging a draw's comparisons with these weights
#'     reproduces the [avg_comparisons()] draw. See [get_draws()].
#' }
#' The ordinary resampling bootstrap (`method = "bootstrap"`) is not
#' available, because a resampled data set can lack outcome states that the
#' fixed profiles need.
#'
#' For `rmsb::blrm()` models, each posterior draw gives a comparison for each
#' profile. The estimate and interval summarize those paired draws, so they
#' can differ from a comparison of [sops()] posterior summaries, especially
#' with `posterior_summary = "median"`.
#'
#' The result has one row per profile, visit (for SOPs), state set, and
#' comparison level, plus a column for each profile variable. Requesting fewer
#' `times`, fewer state sets, or fewer profiles reduces its size. Interval
#' calculations keep every draw of every row, so [inferences()] stops before
#' drawing when `n_draws` times the number of rows exceeds the option
#' `mostr.max_sops_draw_cells` (50 million by default, about 400 MB).
#'
#' @seealso [avg_comparisons()] for averaged comparisons, [sops()] for each
#'   profile's state occupancy probabilities.
#'
#' @examples
#' \dontrun{
#' # One SOP difference per patient, visit, and state, with intervals from
#' # coefficient draws.
#' comparisons(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:30,
#'   absorb = 6
#' ) |>
#'   inferences(method = "mvn", n_draws = 500)
#'
#' # Time in state 1 under treatment versus control for chosen profiles.
#' comparisons(
#'   fit,
#'   newdata = profiles,
#'   variables = list(tx = c(0, 1)),
#'   estimand = "time_in_state",
#'   state_sets = "1",
#'   times = 1:30,
#'   absorb = 6
#' )
#'
#' # The same comparisons averaged within age groups.
#' comparisons(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   estimand = "time_in_state",
#'   by = "age_group",
#'   times = 1:30,
#'   absorb = 6
#' )
#' }
#'
#' @export
comparisons <- function(
  model,
  newdata = NULL,
  variables,
  estimand = c("sop", "time_in_state", "time_benefit"),
  state_sets = NULL,
  comparison = c("difference", "ratio"),
  by = FALSE,
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
  extra_args <- list(...)
  stop_legacy_baseline_args(extra_args)
  estimand <- match.arg(estimand)
  comparison <- match.arg(comparison)
  posterior_summary <- match.arg(posterior_summary)
  by <- resolve_comparisons_by(by)
  validate_comparison_times(times, "comparisons()")
  conf_level <- validate_conf_level(conf_level)
  variables <- validate_avg_comparison_variables(variables, "comparisons()")
  validate_avg_comparison_estimand(estimand, state_sets, comparison)

  if (by$average) {
    return(avg_comparisons(
      model = model,
      newdata = newdata,
      variables = variables,
      estimand = estimand,
      state_sets = state_sets,
      comparison = comparison,
      by = by$by,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_map = time_map,
      baseline_time = baseline_time,
      target_times = target_times,
      time_unit = time_unit,
      id_var = id_var,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      ...
    ))
  }

  setup <- avg_comparison_setup(
    model = model,
    newdata = newdata,
    variables = variables,
    by = NULL,
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    id_var = id_var,
    p2_var = p2_var,
    gap_var = gap_var,
    time_covariates = time_covariates,
    include_re = include_re
  )
  profile_cols <- comparison_profile_columns(
    setup$baseline_data,
    time_var = time_var,
    estimand = estimand
  )

  # Posterior comparisons are summarized from paired draws, so `blrm` SOP
  # draws are always needed here.
  ind <- sops(
    model = model,
    newdata = setup$newdata_pred,
    times = setup$times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    id_var = setup$id_var,
    p2_var = p2_var,
    gap_var = gap_var,
    time_covariates = time_covariates,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = inherits(model, "blrm"),
    ...
  )
  setup <- copy_sops_attrs_to_setup(setup, ind)
  if (estimand != "time_benefit") {
    state_sets <- normalize_comparison_state_sets(
      state_sets,
      attr(ind, "y_levels")
    )
  }

  result <- comparisons_from_individual_sops(
    x = ind,
    setup = setup,
    estimand = estimand,
    state_sets = state_sets,
    comparison = comparison,
    time_map = time_map,
    baseline_time = baseline_time,
    target_times = target_times,
    time_unit = time_unit,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = return_draws
  )
  result <- add_comparison_profile_columns(
    result,
    setup$baseline_data,
    profile_cols
  )

  set_avg_comparison_attrs(
    result = result,
    model = model,
    setup = setup,
    estimand = estimand,
    state_sets = state_sets,
    comparison = comparison,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    time_map = time_map,
    baseline_time = baseline_time,
    target_times = target_times,
    time_unit = time_unit,
    extra_args = extra_args,
    class_name = "markov_comparisons"
  )
}

resolve_comparisons_by <- function(by) {
  if (is.null(by) || isFALSE(by)) {
    return(list(average = FALSE, by = NULL))
  }
  if (isTRUE(by)) {
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

# Profile columns are copied into the output next to the comparison columns,
# so they must not reuse those names. The profile's visit column is dropped:
# it holds the starting visit, not a prediction time. So are columns that
# mostr adds for its own bookkeeping, such as `.markov_source_row`.
comparison_profile_columns <- function(baseline_data, time_var, estimand) {
  profile_cols <- setdiff(names(baseline_data), c("rowid", time_var))
  profile_cols <- profile_cols[!startsWith(profile_cols, ".markov_")]
  reserved <- c(
    "estimand",
    "state_set",
    "term",
    "reference_level",
    "comparison_level",
    "contrast",
    "comparison",
    "time_unit",
    "draw_id",
    "draw",
    inference_columns()
  )
  if (identical(estimand, "sop")) {
    reserved <- c(reserved, "time")
  }
  clash <- intersect(profile_cols, reserved)
  if (length(clash) > 0L) {
    stop(
      "`comparisons()` uses the column name(s) ",
      paste0("`", clash, "`", collapse = ", "),
      " for its output, so the prediction profiles cannot contain them. ",
      "Rename these columns.",
      call. = FALSE
    )
  }
  profile_cols
}

add_comparison_profile_columns <- function(result, baseline_data, profile_cols) {
  index <- result$rowid
  for (col in profile_cols) {
    values <- baseline_data[[col]]
    result[[col]] <- if (is.matrix(values)) {
      values[index, , drop = FALSE]
    } else {
      values[index]
    }
  }
  result
}

comparisons_from_individual_sops <- function(
  x,
  setup,
  estimand,
  state_sets,
  comparison,
  time_map,
  baseline_time,
  target_times,
  time_unit,
  posterior_summary,
  conf_level,
  return_draws
) {
  # As in `avg_comparisons()`, `time_map` applies to time totals only.
  real_time <- !is.null(time_map) && !identical(estimand, "sop")
  if (real_time) {
    target_times <- comparison_real_time_target_times(
      x,
      time_map,
      target_times,
      baseline_time
    )
    x <- interpolate_sops(
      x,
      time_map = time_map,
      target_times = target_times,
      baseline_time = baseline_time
    )
  }
  reduce_setup <- setup
  if (real_time) {
    reduce_setup$times <- sort(unique(x$time))
  }
  reduce_setup$y_levels <- attr(x, "y_levels")

  reduce <- function(data, value_col) {
    individual_comparisons_from_array(
      sops_array = counterfactual_array_from_tidy(
        data = data,
        value_col = value_col,
        setup = reduce_setup
      ),
      setup = reduce_setup,
      estimand = estimand,
      state_sets = state_sets,
      comparison = comparison,
      real_time = real_time,
      value_col = value_col,
      time_unit = time_unit
    )
  }
  point <- reduce(as.data.frame(x), "estimate")

  draws <- attr(x, "draws")
  if (is.null(draws)) {
    return(point)
  }

  # Every draw reduces to the same rows in the same order as the point
  # estimate, so the draws can be summarized column by column.
  value_col <- sop_draw_value_col(draws)
  draw_ids <- unique(draws$draw_id)
  draw_rows <- split(
    seq_len(nrow(draws)),
    factor(draws$draw_id, levels = draw_ids)
  )
  draw_values <- matrix(NA_real_, nrow = length(draw_ids), ncol = nrow(point))
  for (i in seq_along(draw_ids)) {
    reduced <- reduce(draws[draw_rows[[i]], , drop = FALSE], value_col)
    draw_values[i, ] <- reduced[[value_col]]
  }
  cells <- point[, setdiff(names(point), "estimate"), drop = FALSE]
  result <- cbind(
    cells,
    summarize_posterior_draw_matrix(
      draw_values,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    )
  )
  if (isTRUE(return_draws)) {
    attr(result, "draws") <- posterior_draw_matrix_to_df(
      draw_values,
      cells,
      draw_ids
    )
  }
  result
}

# `sops_array` is [scenario, profile, time, state], with the reference level
# of `setup$variables[[1]]` as the first scenario. Returns one row per profile,
# comparison level, state set, and (for SOPs) time, with the profile fastest.
individual_comparisons_from_array <- function(
  sops_array,
  setup,
  estimand,
  state_sets,
  comparison,
  real_time,
  value_col = "estimate",
  time_unit = NULL
) {
  values <- setup$variables[[1]]
  varname <- names(setup$variables)[1]
  times <- setup$times
  n_each <- dim(sops_array)[2]
  n_times <- dim(sops_array)[3]
  level_index <- seq.int(2L, length(values))

  rows <- function(estimate, level_idx, state_set = NULL, by_time = FALSE) {
    out <- data.frame(
      rowid = if (by_time) rep(seq_len(n_each), times = n_times) else
        seq_len(n_each)
    )
    out$estimand <- estimand
    if (by_time) {
      out$time <- rep(times, each = n_each)
    }
    if (!is.null(state_set)) {
      out$state_set <- state_set
    }
    out$term <- varname
    out$reference_level <- values[1]
    out$comparison_level <- values[level_idx]
    out$contrast <- comparison_contrast_label(
      values[level_idx],
      values[1],
      comparison
    )
    out$comparison <- comparison
    if (!is.null(time_unit)) {
      out$time_unit <- time_unit
    }
    out[[value_col]] <- estimate
    out
  }

  pieces <- list()
  if (identical(estimand, "time_benefit")) {
    reference <- comparison_scenario(sops_array, 1L)
    for (level_idx in level_index) {
      benefit <- time_benefit_profile_values(
        reference,
        comparison_scenario(sops_array, level_idx),
        times = times,
        real_time = real_time
      )
      pieces[[length(pieces) + 1L]] <- rows(benefit, level_idx)
    }
    return(bind_rows_fill(pieces))
  }

  state_labels <- as_state_labels(setup$y_levels)
  for (set_i in seq_along(state_sets)) {
    occupancy <- state_set_occupancy(
      sops_array,
      match(as_state_labels(state_sets[[set_i]]), state_labels)
    )
    by_time <- identical(estimand, "sop")
    metric <- if (by_time) {
      occupancy
    } else {
      profile_time_totals(occupancy, times, real_time)
    }
    for (level_idx in level_index) {
      estimate <- comparison_value(
        hi = as.vector(scenario_slice(metric, level_idx)),
        lo = as.vector(scenario_slice(metric, 1L)),
        comparison = comparison
      )
      pieces[[length(pieces) + 1L]] <- rows(
        estimate,
        level_idx,
        state_set = names(state_sets)[set_i],
        by_time = by_time
      )
    }
  }
  bind_rows_fill(pieces)
}

scenario_slice <- function(x, index) {
  if (length(dim(x)) == 3L) {
    x[index, , , drop = FALSE]
  } else {
    x[index, , drop = FALSE]
  }
}

# Sum of the occupancy probabilities in one state set:
# [scenario, profile, time, state] -> [scenario, profile, time].
state_set_occupancy <- function(sops_array, state_index) {
  dims <- dim(sops_array)
  out <- array(0, dim = dims[seq_len(3L)])
  for (state in state_index) {
    out <- out +
      array(sops_array[, , , state, drop = FALSE], dim = dims[seq_len(3L)])
  }
  out
}

# Time in a state set for each scenario and profile:
# [scenario, profile, time] -> [scenario, profile].
profile_time_totals <- function(occupancy, times, real_time) {
  dims <- dim(occupancy)
  flat <- matrix(occupancy, nrow = dims[1L] * dims[2L], ncol = dims[3L])
  totals <- if (isTRUE(real_time)) {
    trapezoid_row_totals(flat, times)
  } else {
    rowSums(flat)
  }
  matrix(totals, nrow = dims[1L], ncol = dims[2L])
}
