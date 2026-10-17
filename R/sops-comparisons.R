# Average SOP comparisons.

#' Average Comparisons From Markov State Occupancy Probabilities
#'
#' Computes average contrasts between counterfactual state occupancy summaries.
#' The first value in `variables[[1]]` is the reference level, and every later
#' value is compared to it. Reverse the order of values to reverse the contrast.
#'
#' @param model A fitted Markov transition model supported by [avg_sops()].
#' @param newdata Optional data frame of standardization profiles. See
#'   [avg_sops()].
#' @param variables The variable to compare: a bare name such as `"tx"`, or a
#'   named list giving its values, such as `list(tx = c(0, 1))`. The first
#'   value is the reference. A bare name compares every observed value with
#'   the first, as in the marginaleffects package: factor levels in level
#'   order, character values in sorted order, `FALSE` then `TRUE`, or `0` then
#'   `1` for a numeric variable coded 0/1. Other numeric variables need their
#'   values listed. The values are read from the data stored with the model,
#'   or from `newdata` for a model without stored data.
#' @param estimand Character scalar. One of `"sop"`, `"time_in_state"`, or
#'   `"time_benefit"`.
#' @param state_sets State selection for `"sop"` and `"time_in_state"`. `NULL`
#'   returns one result per state. An atomic vector is treated as one lumped
#'   state set. A named list returns one result per named state set.
#' @param comparison Character scalar. `"difference"` computes comparison level
#'   minus reference level. `"ratio"` computes comparison level divided by
#'   reference level. `"time_benefit"` only supports `"difference"`.
#' @param by How to average the comparisons of the individual profiles, as in
#'   the marginaleffects package. `TRUE` (the default) averages over all
#'   profiles. A character vector of column names averages within the groups
#'   those columns define. `FALSE` or `NULL` does not average and returns the
#'   [comparisons()] result, one comparison per profile; analytical intervals
#'   for it use conditional variance only, and the ordinary bootstrap is not
#'   available.
#' @param times Required visit-scale time points. See [avg_sops()].
#' @param y_levels A vector of state levels. See [avg_sops()].
#' @param absorb The absorbing state. See [avg_sops()].
#' @param time_map Optional named numeric vector or data frame mapping visit
#'   labels to real elapsed times. Used by `"time_in_state"` and
#'   `"time_benefit"`.
#' @param baseline_time Real time at which patients are in their starting state,
#'   used with `time_map`. The default `0` means that real-time totals start
#'   from the starting state at time 0. Set to `NULL` to start at the first
#'   mapped visit instead. See [interpolate_sops()].
#' @param target_times Optional numeric grid of real times, used with
#'   `time_map`. Probabilities are linearly interpolated at these times, and
#'   totals are calculated with the trapezoidal rule over the grid. When
#'   omitted, the grid is `baseline_time` followed by the mapped visit times, so
#'   totals cover the whole period from the starting state to the last visit.
#'   With `baseline_time = NULL`, the omitted grid is the mapped visit times
#'   only.
#' @param time_unit Optional label stored in output.
#' @param id_var Name of the patient ID variable. See [avg_sops()].
#' @param time_var Name of the time variable in the model. See [avg_sops()].
#' @param p_var Name of the previous state variable in the model. See
#'   [avg_sops()].
#' @param p2_var Optional second previous-state variable. See [avg_sops()].
#' @param gap_var Name of the time gap_var variable, if used. See [avg_sops()].
#' @param time_covariates Optional time-varying covariate lookup table for explicit
#'   precomputed time-basis columns. See [avg_sops()].
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws.
#' @param seed Optional random seed for reproducible posterior draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior comparison draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior comparison draws.
#'   Frequentist uncertainty draws are stored by [inferences()].
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A data frame of class `markov_avg_comparisons`, or, with
#'   `by = FALSE`, of class `markov_comparisons`.
#'
#' @details
#' `estimand = "sop"` and `estimand = "time_in_state"` are computed from marginal
#' SOPs. `estimand = "time_benefit"` is computed before marginalization from
#' paired patient/profile-level counterfactual state distributions, because it
#' is nonlinear in the two state distributions.
#'
#' Frequentist uncertainty follows the existing package workflow:
#' `avg_comparisons(...) |> inferences(...)`. Bayesian `blrm` models use their
#' posterior SOP draws directly and therefore return uncertainty from
#' `avg_comparisons()`.
#'
#' For supported frequentist models, `inferences(method = "delta")` provides
#' analytical intervals for differences in SOPs or time in state. The default
#' `vcov = "unconditional"` accounts for coefficient estimation and sampling
#' the patients used in the averages. Use `vcov = "conditional"` to treat
#' those patients' starting states and covariates as given and account only
#' for coefficient estimation. User-supplied `newdata` requires the conditional
#' choice. See [inferences()] for model restrictions.
#'
#' @seealso [comparisons()] for the comparison of each profile before
#'   averaging, [avg_sops()], [sops()], [inferences()]
#'
#' @examples
#' \dontrun{
#' avg_comparisons(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   estimand = "time_in_state",
#'   state_sets = "1",
#'   times = 1:30,
#'   y_levels = 1:6,
#'   absorb = 6
#' ) |>
#'   inferences(method = "mvn", n_draws = 500)
#' }
#'
#' @export
avg_comparisons <- function(
  model,
  newdata = NULL,
  variables,
  estimand = c("sop", "time_in_state", "time_benefit"),
  state_sets = NULL,
  comparison = c("difference", "ratio"),
  by = TRUE,
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
  by_choice <- resolve_comparisons_by(by)

  validate_comparison_times(times, "avg_comparisons()")
  conf_level <- validate_conf_level(conf_level)
  variables <- resolve_comparison_variables(
    variables,
    model = model,
    newdata = newdata,
    fn = "avg_comparisons()"
  )
  validate_avg_comparison_estimand(estimand, state_sets, comparison)

  if (!by_choice$average) {
    return(comparisons(
      model = model,
      newdata = newdata,
      variables = variables,
      estimand = estimand,
      state_sets = state_sets,
      comparison = comparison,
      by = FALSE,
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
  # From here on, `by` holds the grouping columns, or `NULL` to average over
  # all profiles.
  by <- by_choice$by

  if (estimand == "time_benefit") {
    setup <- avg_comparison_setup(
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
      ...
    )
    tb <- avg_comparison_time_benefit_point(
      model = model,
      setup = setup,
      y_levels = y_levels,
      absorb = absorb,
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
      comparison = comparison,
      time_map = time_map,
      baseline_time = baseline_time,
      target_times = target_times,
      time_unit = time_unit,
      ...
    )
    result <- tb$result
    setup <- tb$setup
  } else {
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
      return_draws = inherits(model, "blrm") || isTRUE(return_draws),
      extra_args = extra_args
    )
    state_sets <- normalize_comparison_state_sets(
      state_sets,
      attr(avg, "y_levels")
    )
    result <- avg_comparison_from_avg_sops(
      avg,
      estimand = estimand,
      state_sets = state_sets,
      comparison = comparison,
      time_map = time_map,
      baseline_time = baseline_time,
      target_times = target_times,
      time_unit = time_unit,
      return_draws = return_draws
    )
    setup <- avg_comparison_setup_from_sops(avg)
  }

  result <- set_avg_comparison_attrs(
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
    extra_args = extra_args
  )

  result
}
