# Time in states for each prediction profile.

#' Time in States for Each Profile
#'
#' Estimates each prediction profile's expected total time in each state set
#' over the requested follow-up, including time accumulated across repeated
#' visits to a state. Profiles are, by default, the fitted patients' starting
#' states and covariates. With `variables`, each profile is predicted once
#' for every requested value. [avg_state_time()] averages these times over
#' the profiles.
#'
#' @inheritParams avg_state_time
#' @param newdata Optional data frame of prediction profiles. Each row is one
#'   profile with its own starting state and covariates. If `NULL`, uses the
#'   starting profiles stored by [orm_markov()], [blrm_markov()], or
#'   [vglm_markov()], one per fitted patient.
#' @param variables Optional covariates to set: a named list of values, such
#'   as `list(tx = c(0, 1))`, or bare names, such as `"tx"`, whose values are
#'   chosen as in [avg_sops()]. Each profile is predicted once for every
#'   combination of the values. The default, `NULL`, uses each profile's own
#'   covariates.
#' @param by `FALSE` (the default) or `NULL` returns one total per profile.
#'   `TRUE` averages over all profiles, and a character vector of column names
#'   averages within the groups those columns define; both return the result
#'   of [avg_state_time()] with the same arguments.
#' @param return_draws Logical. For `blrm`, store the posterior draws of each
#'   total for extraction with [get_draws()].
#'
#' @return A data frame of class `markov_state_time` with one row per state
#'   set, value of `variables`, and profile:
#'   \item{rowid}{Row number of the prediction: the profile, or with
#'     `variables`, the profile and value combination.}
#'   \item{rowidcf}{With `variables`, the profile that the row predicts.}
#'   \item{state_set}{Name of the state set.}
#'   \item{(profile columns)}{The profile's covariates and starting state. The
#'     `variables` columns hold the value used for the row.}
#'   \item{estimate}{The profile's expected total time in the state set.}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     summaries of the total.}
#'   If supplied, `time_unit` is included as a label.
#'
#' @details
#' Totals are calculated as in [avg_state_time()], for one profile at a time:
#' without `time_map`, the sum of the profile's occupancy probabilities over
#' the requested visits; with `time_map`, the area under those probabilities
#' on the elapsed-time scale, starting at `baseline_time` from the profile's
#' starting state, or from the state set in `variables` when it sets the
#' previous state. The mean of the totals over the profiles equals the
#' [avg_state_time()] estimate.
#'
#' For frequentist models, pass the result to [inferences()]. Because each
#' profile's starting state and covariates are held as given, every method
#' reflects uncertainty in the model coefficients only:
#' \itemize{
#'   \item `method = "delta"` uses conditional variance, \eqn{G V_\beta
#'     G^\top}, where \eqn{G} holds the derivatives of the totals with respect
#'     to the model coefficients and \eqn{V_\beta} is their covariance. It
#'     supports first-order full proportional-odds `orm`, `vglm`, and
#'     `robcov_vglm` fits; [stats::vcov()] then returns the covariance between
#'     rows.
#'   \item `method = "mvn"`, `"score_bootstrap"`, or `"fwb"` compute every
#'     profile's totals within each coefficient draw or weighted refit. For
#'     the stored profiles, score-bootstrap and FWB draws also carry each
#'     patient's weight in that draw (`score_weight` or `fwb_weight`);
#'     averaging a draw's totals with these weights reproduces the
#'     [avg_state_time()] draw. See [get_draws()].
#' }
#' The ordinary resampling bootstrap (`method = "bootstrap"`) is not
#' available, because a resampled data set can lack outcome states that the
#' fixed profiles need.
#'
#' For `rmsb::blrm()` models, each posterior draw gives every profile's
#' totals, and the estimate and interval summarize those draws.
#'
#' Interval calculations keep every draw of every row, so [inferences()]
#' stops before drawing when `n_draws` times the number of rows exceeds the
#' option `mostr.max_sops_draw_cells` (50 million by default, about 400 MB).
#'
#' @inheritSection sops Choosing visits
#'
#' @seealso [avg_state_time()] for averaged times, [sops()] for each profile's
#'   state occupancy probabilities, [comparisons()] for each profile's
#'   comparisons.
#'
#' @examples
#' \dontrun{
#' # Each patient's expected days at home under control and under treatment.
#' state_time(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   state_sets = list(home = "1"),
#'   times = 1:30
#' ) |>
#'   inferences(method = "delta")
#' }
#'
#' @export
state_time <- function(
  model,
  newdata = NULL,
  variables = NULL,
  state_sets = NULL,
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
  validate_state_time_arguments(
    times,
    baseline_time = baseline_time,
    baseline_time_supplied = !missing(baseline_time),
    time_map = time_map,
    target_times = target_times,
    time_unit = time_unit,
    fn = "state_time()"
  )
  posterior_summary <- match.arg(posterior_summary)
  conf_level <- validate_conf_level(conf_level)
  extra_args <- list(...)

  by_choice <- resolve_by_choice(by)
  if (by_choice$average) {
    # `baseline_time` is passed on only if supplied, because supplying it
    # without `time_map` is an error.
    return(do.call(
      avg_state_time,
      c(
        list(
          model = model,
          newdata = newdata,
          variables = variables,
          state_sets = state_sets,
          by = by_choice$by %||% TRUE,
          times = times,
          y_levels = y_levels,
          absorb = absorb,
          time_map = time_map,
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
          return_draws = return_draws
        ),
        if (!missing(baseline_time)) list(baseline_time = baseline_time),
        extra_args
      )
    ))
  }

  variables <- resolve_prediction_variables(variables, model, newdata)
  if (!is.null(state_sets)) {
    validate_avg_time_state_sets(state_sets)
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
  profile_cols <- profile_output_columns(
    setup$newdata_pred,
    time_var = time_var,
    reserved = c(
      "state_set",
      if (!is.null(time_unit)) "time_unit",
      "draw_id",
      "draw",
      inference_columns()
    ),
    fn = "state_time()"
  )
  # `sops()` below receives the stacked prediction rows as supplied `newdata`,
  # so check `times` here against the original profiles.
  warn_skipped_sop_visits(
    model,
    setup$newdata_orig,
    times,
    time_var,
    gap_var,
    setup$newdata_supplied
  )

  # Posterior totals are summarized from each draw's totals, so `blrm` SOP
  # draws are always needed here.
  ind <- muffle_skipped_sop_visits(sops(
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
  ))
  setup <- copy_sops_attrs_to_setup(setup, ind)
  state_sets <- normalize_comparison_state_sets(
    state_sets,
    attr(ind, "y_levels")
  )

  reducer <- profile_array_reducer(
    "time",
    list(state_sets = state_sets, time_map = time_map, time_unit = time_unit)
  )
  result <- profile_results_from_individual_sops(
    x = ind,
    setup = setup,
    reducer = reducer,
    time_map = time_map,
    baseline_time = baseline_time,
    target_times = target_times,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = return_draws
  )
  # `rowid` numbers the stacked scenario-profile rows, so it indexes the
  # stacked prediction data, whose `variables` columns hold the values used.
  result <- add_profile_columns(
    result,
    data = setup$newdata_pred,
    index = result$rowid,
    cols = profile_cols
  )

  set_avg_comparison_attrs(
    result = result,
    model = model,
    setup = setup,
    estimand = "time_in_state",
    state_sets = state_sets,
    comparison = NULL,
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
    class_name = "markov_state_time"
  )
}
