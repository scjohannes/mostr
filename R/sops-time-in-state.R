# SOP time-in-state summaries.

# Trapezoidal integral of `value` over `time`, for values in any order:
# missing pairs are dropped and values at the same time are added.
trapezoid_auc <- function(time, value) {
  ok <- !is.na(time) & !is.na(value)
  time <- time[ok]
  value <- value[ok]
  if (length(time) < 2) {
    return(0)
  }
  agg <- stats::aggregate(value ~ time, FUN = sum)
  agg <- agg[order(agg$time), , drop = FALSE]
  if (nrow(agg) < 2) {
    return(0)
  }
  trapezoid_row_totals(matrix(agg$value, nrow = 1L), agg$time)
}

time_in_state_tidy <- function(x, target_states, real_time = FALSE) {
  value_col <- if ("estimate" %in% names(x)) {
    "estimate"
  } else if ("draw" %in% names(x)) {
    "draw"
  } else {
    mostr_abort("Tidy SOP input must contain an `estimate` or `draw` column.")
  }

  keep <- as.character(x$state) %in% as.character(target_states)
  if (!any(keep)) {
    mostr_abort("Target states not found in SOP data frame.")
  }
  x <- x[keep, , drop = FALSE]
  measure_cols <- sop_measure_cols(x)
  weight_cols <- intersect("fwb_weight", names(x))
  # Every other column identifies a group: a profile, scenario, or draw.
  # The tests from `inferences(null = )` are not identifiers; they are
  # dropped, because a null value for probabilities says nothing about a
  # total time.
  group_cols <- setdiff(
    names(x),
    c("time", "state", measure_cols, sop_test_cols(x), weight_cols)
  )
  weight_metadata <- if (length(weight_cols) > 0L) {
    unique(x[, c(group_cols, weight_cols), drop = FALSE])
  } else {
    NULL
  }

  # Profile columns that the model does not use may be missing; such a
  # profile still keeps its own group.
  by_time <- aggregate_value(
    x,
    value_col,
    c(group_cols, "time"),
    sum,
    keep_na_groups = TRUE
  )

  if (!real_time) {
    if (length(group_cols) == 0) {
      out <- data.frame(total_time = sum(by_time[[value_col]], na.rm = TRUE))
      return(out)
    }
    out <- aggregate_value(
      by_time,
      value_col,
      group_cols,
      sum,
      keep_na_groups = TRUE
    )
    names(out)[names(out) == value_col] <- "total_time"
    if (!is.null(weight_metadata)) {
      out <- left_join_preserve_order(out, weight_metadata, by = group_cols)
    }
    return(out)
  }

  groups <- split(
    seq_len(nrow(by_time)),
    split_key(by_time, group_cols, keep_na = TRUE),
    drop = TRUE
  )
  pieces <- lapply(groups, function(idx) {
    group <- by_time[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    meta$total_time <- trapezoid_auc(group$time, group[[value_col]])
    meta
  })
  out <- bind_rows_fill(pieces)
  if (!is.null(weight_metadata)) {
    out <- left_join_preserve_order(out, weight_metadata, by = group_cols)
  }
  out
}

# Whether `x` stores its draws in a table with a `draw_id` column, from which
# the limits of any total of its SOPs can be recomputed.
sop_has_draw_table <- function(x) {
  draws <- attr(x, "draws")
  is.data.frame(draws) && "draw_id" %in% names(draws)
}

# The uncertainty columns that the totals of SOPs with draws carry, in that
# order.
time_in_state_uncertainty_cols <- function(x) {
  intersect(c("conf.low", "conf.high", "std.error"), names(x))
}

# `x` without its uncertainty and test columns, so that its totals come from
# the estimates alone.
sop_point_estimates <- function(x) {
  for (column in setdiff(inference_columns(), "estimate")) {
    x[[column]] <- NULL
  }
  x
}

# Totals of SOPs that have limits but no draws: posterior SOPs from `sops()`
# or `avg_sops()` with `return_draws = FALSE`, or SOPs from
# `inferences(return_draws = FALSE)`. The limits of a sum cannot be computed
# from the limits of its terms, which ignore how the SOPs vary together, so the
# totals get NA limits and a warning instead of losing the columns silently.
# `source` is the SOPs as the user passed them, before any interval columns
# were dropped for interpolation.
time_in_state_without_draws <- function(result, source) {
  columns <- time_in_state_uncertainty_cols(source)
  if (length(columns) == 0L) {
    return(result)
  }
  rerun <- if (identical(attr(source, "method"), "posterior")) {
    paste0(
      "rerun `sops()` or `avg_sops()` with `return_draws = TRUE`, or use ",
      "`state_time()` or `avg_state_time()`."
    )
  } else {
    "rerun `inferences()` with `return_draws = TRUE`."
  }
  warning(
    "`time_in_state()` found uncertainty columns but no stored draws, and ",
    "the limits of a total cannot be computed from the limits of its SOPs. ",
    "The totals have NA `",
    paste(columns, collapse = "`, `"),
    "`. To compute them, ",
    rerun,
    call. = FALSE
  )
  for (column in columns) {
    result[[column]] <- rep(NA_real_, nrow(result))
  }
  result
}

time_in_state_tidy_inference <- function(x, target_states, real_time, source = x) {
  result <- time_in_state_tidy(x, target_states, real_time = real_time)
  if (!sop_has_draw_table(x)) {
    return(time_in_state_without_draws(result, source))
  }
  draws <- attr(x, "draws")
  draws <- draws[draws$time %in% x$time, , drop = FALSE]
  draw_result <- time_in_state_tidy(
    draws,
    target_states,
    real_time = real_time
  )
  weight_cols <- intersect(
    "fwb_weight",
    names(draw_result)
  )
  group_cols <- setdiff(
    names(draw_result),
    c("draw_id", "total_time", weight_cols)
  )
  draws_for_ci <- draw_result
  draws_for_ci$estimate <- draws_for_ci$total_time
  point_estimates <- result
  point_estimates$estimate <- point_estimates$total_time
  ci <- compute_ci_from_draws(
    draws_df = draws_for_ci,
    group_cols = group_cols,
    conf_level = attr(x, "conf_level") %||% 0.95,
    conf_type = attr(x, "conf_type") %||% "perc",
    point_estimates = point_estimates,
    keep_na_groups = TRUE
  )

  result <- left_join_preserve_order(result, ci, by = group_cols)
  attr(result, "draws") <- draw_result
  result
}

# The recorded method alone decides, since `inferences()` replaces it on every
# run. Delta SOPs that lost their derivatives stop in `time_in_state_delta()`.
time_in_state_is_delta <- function(x) {
  identical(attr(x, "method"), "delta")
}

# Delta-method uncertainty of the totals of SOPs from
# `inferences(method = "delta")`. Each total is a weighted sum of SOP rows,
# plus a starting-state term when real time starts at `baseline_time`:
#   total_g = sum_r w_gr * p_r + b * a_g,
# where w_gr is the time weight of row r (1 per visit without `time_map`, the
# trapezoid weights of the interpolation grid with it) when row r is a target
# state of group g and 0 otherwise, b is the weight of `baseline_time`, and a_g
# is the probability that group g starts in a target state. The weights are
# applied to what the SOPs store: their derivatives with respect to the
# coefficients (conditional variance, w J V J' w') or their per-patient
# influence values (unconditional variance), so the covariance of every pair of
# SOP rows is never formed. Under unconditional variance, a_g is the share of
# sampled patients starting in a target state, whose influence is added.
time_in_state_delta <- function(
  x,
  target_states,
  time_map,
  baseline_time,
  target_times
) {
  analytical <- attr(x, "analytical", exact = TRUE)
  if (is.null(analytical)) {
    mostr_abort(
      "These delta-method SOPs no longer carry the derivatives needed for ",
      "standard errors of the total time. Run `inferences(method = \"delta\")` ",
      "again and pass its result to `time_in_state()` unchanged."
    )
  }
  # The totals come from the estimates alone. Without the interval columns,
  # `interpolate_sops()` does not warn that it interpolates interval endpoints,
  # which are replaced below.
  point <- sop_point_estimates(x)
  result <- if (is.null(time_map)) {
    time_in_state_tidy(point, target_states, real_time = FALSE)
  } else {
    grid <- comparison_real_time_target_times(
      point,
      time_map,
      target_times,
      baseline_time
    )
    interpolated <- interpolate_sops(
      point,
      time_map = time_map,
      target_times = grid,
      baseline_time = baseline_time
    )
    time_in_state_tidy(interpolated, target_states, real_time = TRUE)
  }

  source <- as.data.frame(x)
  weight_cols <- intersect("fwb_weight", names(source))
  group_cols <- setdiff(
    names(source),
    c(
      "time",
      "state",
      sop_measure_cols(source),
      sop_test_cols(source),
      weight_cols
    )
  )
  source_row <- match(
    delta_row_key(delta_row_key_frame(x)),
    analytical$row_key
  )
  if (anyNA(source_row) || !all(group_cols %in% names(result))) {
    mostr_abort(
      "The SOP rows do not align with their stored delta-method derivatives. ",
      "Pass the result of `inferences(method = \"delta\")` to ",
      "`time_in_state()` unchanged."
    )
  }
  group <- if (length(group_cols) == 0L) {
    rep(1L, nrow(source))
  } else {
    match(
      delta_row_key(source[, group_cols, drop = FALSE]),
      delta_row_key(as.data.frame(result)[, group_cols, drop = FALSE])
    )
  }

  source_time <- as.character(source$time)
  integration <- if (is.null(time_map)) {
    list(
      visit = stats::setNames(
        rep(1, length(unique(source_time))),
        unique(source_time)
      ),
      baseline = 0
    )
  } else {
    delta_real_time_weights(
      x,
      list(
        time_map = time_map,
        baseline_time = baseline_time,
        target_times = target_times
      )
    )
  }
  time_weight <- unname(integration$visit[match(
    source_time,
    names(integration$visit)
  )])
  targets <- as.character(target_states)
  keep <- which(
    as.character(source$state) %in% targets & time_weight != 0
  )
  delta_assert_bytes(
    as.double(nrow(source)) * 64 + as.double(length(keep)) * 12,
    "The analytical time-in-state operator"
  )
  if (anyNA(group[keep]) || anyNA(time_weight[keep])) {
    mostr_abort(
      "The SOP rows do not align with the time-in-state groups."
    )
  }
  rows_by_group <- split(
    keep,
    factor(group[keep], levels = seq_len(nrow(result)))
  )
  operator <- lapply(rows_by_group, function(rows) {
    list(index = source_row[rows], weight = time_weight[rows])
  })

  # Probability that each group starts in a target state.
  p_var <- markov_result_settings(x)$p_var
  start_share <- numeric(nrow(result))
  start <- NULL
  if (integration$baseline != 0) {
    if (p_var %in% group_cols) {
      # Each profile, or each scenario that sets the starting state, starts
      # in its own state.
      start_share <- as.numeric(as.character(result[[p_var]]) %in% targets)
    } else if (inherits(x, "markov_avg_sops")) {
      baseline <- baseline_rows_for_anchor(x)
      if (!p_var %in% names(baseline)) {
        mostr_abort(
          "Previous-state variable `",
          p_var,
          "` not found in baseline data."
        )
      }
      start <- as.character(baseline[[p_var]])
      start_share <- rep(mean(start %in% targets), nrow(result))
    } else {
      mostr_abort(
        "Previous-state variable `",
        p_var,
        "` not found in SOP output."
      )
    }
  }

  estimate <- rep(0, length(analytical$row_key))
  estimate[source_row] <- source$estimate
  calculated <- drop(delta_apply_comparison_operator(operator, estimate)) +
    integration$baseline * start_share
  # With `time_map`, the weights assume one shared visit grid and probabilities
  # that sum to one over the states (so that interpolation's renormalization
  # changes nothing). Rows filtered out by the user can break either; the
  # standard errors would then belong to a different total, so stop instead.
  if (
    any(!is.finite(calculated)) ||
      estimates_differ(calculated, result$total_time)
  ) {
    mostr_abort(
      "With `time_map`, delta-method standard errors of the total time need ",
      "the same visits and all states in every group. Pass the result of ",
      "`inferences(method = \"delta\")` to `time_in_state()` unchanged."
    )
  }

  propagated <- delta_propagate_comparison_state(x, operator)
  standard_error <- propagated$std.error
  if (
    identical(propagated$analytical$representation, "influence") &&
      !is.null(start)
  ) {
    profile_ids <- as.character(propagated$analytical$profile_ids)
    ids <- delta_profile_ids(baseline, markov_result_settings(x)$id_var)
    index <- match(profile_ids, ids)
    if (anyNA(index) || length(ids) != length(profile_ids)) {
      mostr_abort(
        "Baseline-anchor patients do not align with the analytical patient influences."
      )
    }
    indicator <- as.numeric(start[index] %in% targets)
    influence <- propagated$analytical$influence +
      integration$baseline * (indicator - mean(indicator))
    standard_error <- sqrt(
      apply(influence, 2L, stats::var) / length(profile_ids)
    )
  }

  bounds <- delta_interval_bounds(
    estimate = result$total_time,
    standard_error = standard_error,
    conf_level = attr(x, "conf_level") %||% 0.95,
    conf_type = "wald"
  )
  result$conf.low <- bounds$conf.low
  result$conf.high <- bounds$conf.high
  result$std.error <- standard_error
  result
}


#' Compute Total Time in Target State(s)
#'
#' Calculates the expected total time spent in specified target state(s).
#' By default this sums state occupancy probabilities over visit-scale unit
#' steps. When `time_map` is supplied, it maps/interpolates to
#' real time and uses trapezoidal AUC. Already interpolated SOP output from
#' [interpolate_sops()] is integrated on its current real-time grid.
#' It automatically adapts to the input format:
#' \itemize{
#'   \item **Patient-Level:** If input is an array from \code{soprob_markov}, it returns time-in-state for each patient.
#'   \item **Tidy SOP:** If input is a `markov_sops` or `markov_avg_sops`
#'     object, it returns grouped total time. Stored simulation, bootstrap, or
#'     posterior draws are reduced identically and used to recompute confidence
#'     intervals and standard errors. For SOPs from
#'     `inferences(method = "delta")`, standard errors and confidence limits of
#'     the totals are computed with the delta method; see Details. SOPs with
#'     confidence limits but no stored draws give totals with `NA` limits and
#'     a warning; see Details.
#' }
#'
#' @param x Input object.
#'   \itemize{
#'     \item **Array:** `[Patients x Time x States]` (Frequentist) or `[Draws x Patients x Time x States]` (Bayes).
#'     \item **Tidy SOP:** A `markov_sops`, `markov_avg_sops`, or
#'       `markov_interpolated_sops` object.
#'   }
#' @param target_states Vector of target state(s) to include in the time calculation.
#'   Can be integer indices or character names (e.g., \code{1} or \code{c("Home", "Rehab")}).
#' @param time_map Optional named numeric vector or data frame mapping visit
#'   labels to real elapsed times. Supplying this switches tidy SOP and array
#'   inputs to trapezoidal real-time AUC.
#' @param baseline_time Real time of the observed baseline state for tidy SOP
#'   outputs. The default `0` supplies a fixed interpolation anchor to
#'   [interpolate_sops()]. Set to `NULL` to disable baseline anchoring.
#' @param target_times Optional numeric real-time grid for tidy SOP outputs when
#'   `time_map` is supplied. This controls the interpolation grid used for AUC,
#'   for example `target_times = 1:28` with `baseline_time = 0` uses day 0 as an
#'   anchor but starts the AUC at day 1. If `NULL`, integration runs from
#'   `baseline_time` through the mapped follow-up times, matching [avg_state_time()]
#'   and [avg_comparisons()]. With `baseline_time = NULL`, it starts at the
#'   first mapped follow-up time.
#'
#' @return
#' \itemize{
#'   \item **Input = Array (Frequentist):** A named numeric vector of length \code{n_patients}.
#'   \item **Input = Array (Bayesian):** A matrix of dimension \code{[n_draws x n_patients]}.
#'   \item **Input = Tidy SOP:** A data frame with `total_time` and grouping
#'     columns. When stored draws are available, draw-based `conf.low`,
#'     `conf.high`, and `std.error` columns are recomputed and the reduced draws
#'     are retained in a `draws` attribute. Any draw-specific baseline anchors
#'     have already been propagated by [interpolate_sops()]. For SOPs from
#'     `inferences(method = "delta")`, `conf.low`, `conf.high`, and
#'     `std.error` are delta-method results for the totals. For SOPs with
#'     limits but neither draws nor delta-method results, these columns are
#'     `NA`. Tests from
#'     `inferences(null = )` are dropped, because a null value for
#'     probabilities does not apply to a total time.
#' }
#'
#' @details
#' ## SOPs with limits but no draws
#'
#' The confidence limits of a total cannot be computed from the limits of the
#' SOPs it adds up, because those limits ignore how the SOPs vary together.
#' `time_in_state()` therefore needs the draws behind them. Posterior SOPs from
#' [sops()] or [avg_sops()] keep their draws only with `return_draws = TRUE`,
#' and SOPs from [inferences()] keep them unless `return_draws = FALSE`.
#' Without the draws, the totals are still computed from the estimates, but
#' `conf.low`, `conf.high`, and `std.error` are `NA`, and a warning says which
#' call to rerun. [state_time()] and [avg_state_time()] compute posterior
#' totals with their limits directly, without `return_draws = TRUE`.
#'
#' ## Delta-method uncertainty
#'
#' A total time is a weighted sum of SOPs: each visit counts once without
#' `time_map`, and with `time_map` each visit gets the weight it carries in
#' the area under the interpolated curve. With `baseline_time`, the curve
#' starts at that time from the starting state, so the total also depends on
#' the share of profiles that start in a target state. For SOPs from
#' `inferences(method = "delta")`,
#' `time_in_state()` applies these weights to how each SOP changes with each
#' model coefficient, and so obtains the standard error of the total. It keeps
#' the variance chosen with `vcov` in [inferences()]:
#' \itemize{
#'   \item **Conditional variance** accounts for coefficient estimation and
#'     treats the patients' starting states and covariates as given. It is the
#'     only choice for [sops()] results, one row per profile.
#'   \item **Unconditional variance** (the default for [avg_sops()]) also
#'     accounts for which patients were sampled. With `time_map` and
#'     `baseline_time`, this includes the share of sampled patients who start
#'     in a target state.
#' }
#' Confidence limits are the total plus and minus the normal quantile for the
#' stored `conf_level` times the standard error, as in [avg_state_time()] and
#' [state_time()] with `method = "delta"`, which give the same results. The
#' calculation uses only the weighted sums, never the covariance of every pair
#' of SOPs, within the memory limit `mostr.delta_max_bytes` (see
#' [mostr_options]).
#'
#' Pass the result of [inferences()] unchanged. `time_in_state()` stops with
#' an error, instead of returning totals without uncertainty, for delta
#' results already passed through [interpolate_sops()] (use `time_map` here
#' instead), and, with `time_map`, for delta results whose rows were filtered
#' so that the groups no longer share the same visits or some states are
#' missing.
#'
#' @examples
#' \dontrun{
#' # --- Scenario 1: Patient-Level Estimates ---
#' sops_arr <- soprob_markov(model, data, times = 1:30)
#' days_home <- time_in_state(sops_arr, target_states = 1)
#'
#' # Real-time AUC after fitting on factor visit indices
#' real_days_home <- time_in_state(
#'   avg,
#'   target_states = 1,
#'   time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
#'   baseline_time = 0,
#'   target_times = 1:28
#' )
#'
#' # --- Scenario 2: Bootstrap Inference ---
#' avg_boot <- avg_sops(
#'   model,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:30
#' ) |>
#'   inferences(method = "bootstrap", n_draws = 100)
#' boot_effects <- time_in_state(avg_boot, target_states = 1)
#' }
#'
#' @keywords time-in-state auc bootstrap
#' @seealso [avg_state_time()] for model-level average time estimates with analytical,
#'   simulation, bootstrap, or posterior uncertainty.
#' @export
time_in_state <- function(
  x,
  target_states = 1,
  time_map = NULL,
  baseline_time = 0,
  target_times = NULL
) {
  sops <- x
  baseline_time_supplied <- !missing(baseline_time)
  use_real_time <- !is.null(time_map)

  if (inherits(sops, "markov_interpolated_sops")) {
    # The stored derivatives belong to the visit SOPs, not to the interpolated
    # rows, and `interpolate_sops()` interpolated only the interval endpoints.
    if (time_in_state_is_delta(sops)) {
      mostr_abort(
        "`time_in_state()` cannot compute delta-method standard errors from ",
        "SOPs already interpolated with `interpolate_sops()`. Pass the SOPs ",
        "from `inferences(method = \"delta\")` with `time_map` instead."
      )
    }
    if (!is.null(time_map) || baseline_time_supplied) {
      mostr_abort(
        "`markov_interpolated_sops` is already on the real-time scale; ",
        "omit `time_map` and `baseline_time`."
      )
    }
    if (!is.null(target_times)) {
      target_times <- validate_sop_xout(
        target_times,
        min(sops$time),
        max(sops$time)
      )
      missing_xout <- setdiff(target_times, unique(sops$time))
      if (length(missing_xout) > 0) {
        mostr_abort(
          "`target_times` for already interpolated SOPs must be contained in ",
          "the existing `time` values."
        )
      }
      sops <- sops[sops$time %in% target_times, , drop = FALSE]
    }
    return(time_in_state_tidy_inference(
      sops,
      target_states,
      real_time = TRUE
    ))
  }

  if (
    baseline_time_supplied &&
      !is.null(baseline_time) &&
      is.null(time_map)
  ) {
    mostr_abort(
      "`baseline_time` requires `time_map` for real-time interpolation."
    )
  }

  if (inherits(sops, c("markov_sops", "markov_avg_sops"))) {
    if (!is.null(target_times) && !use_real_time) {
      mostr_abort(
        "`target_times` requires `time_map` or an already interpolated SOP object."
      )
    }
    if (time_in_state_is_delta(sops)) {
      return(time_in_state_delta(
        sops,
        target_states,
        time_map = time_map,
        baseline_time = baseline_time,
        target_times = target_times
      ))
    }
    if (use_real_time) {
      target_times <- comparison_real_time_target_times(
        sops,
        time_map,
        target_times,
        baseline_time
      )
      # Without draws the totals get NA limits, so the interval columns are
      # dropped first: `interpolate_sops()` would otherwise warn that it
      # interpolates interval endpoints, which would then be discarded.
      interpolated <- interpolate_sops(
        if (sop_has_draw_table(sops)) sops else sop_point_estimates(sops),
        time_map = time_map,
        target_times = target_times,
        baseline_time = baseline_time
      )
      return(time_in_state_tidy_inference(
        interpolated,
        target_states,
        real_time = TRUE,
        source = sops
      ))
    }
    return(time_in_state_tidy_inference(
      sops,
      target_states,
      real_time = FALSE
    ))
  }

  if (inherits(sops, "data.frame")) {
    mostr_abort(
      "A data frame given to `time_in_state()` must be SOPs from `sops()`, ",
      "`avg_sops()`, or `interpolate_sops()`."
    )
  }

  # Array input.
  if (!is.null(target_times)) {
    mostr_abort("`target_times` is only supported for tidy SOP outputs.")
  }

  # --- 1. Detect Dimensions ---
  dims <- dim(sops)
  if (is.null(dims)) {
    mostr_abort("Input 'sops' must be SOPs or an array.")
  }

  ndim <- length(dims)
  dnames <- dimnames(sops)

  if (ndim == 3) {
    # Frequentist: [Patients, Time, States]
    state_dim <- 3
  } else if (ndim == 4) {
    # Bayesian: [Draws, Patients, Time, States]
    state_dim <- 4
  } else {
    mostr_abort("Input array must be 3D or 4D.")
  }

  # --- 2. Resolve Target States ---
  state_names <- dnames[[state_dim]]
  if (is.null(state_names)) {
    state_names <- as.character(1:dims[state_dim])
  }

  target_chars <- as.character(target_states)
  missing_states <- setdiff(target_chars, state_names)
  if (length(missing_states) > 0) {
    mostr_abort(paste(
      "Target states not found in input array:",
      paste(missing_states, collapse = ", ")
    ))
  }

  # --- 3. Filter & Aggregate States ---
  subset_args <- rep(list(TRUE), ndim)
  subset_args[[state_dim]] <- match(target_chars, state_names)

  sops_slice <- do.call("[", c(list(sops), subset_args, drop = FALSE))
  prob_in_target <- apply(sops_slice, (1:ndim)[-state_dim], sum)

  # --- 4. Integrate Over Time ---
  if (use_real_time) {
    if (baseline_time_supplied && !is.null(baseline_time)) {
      mostr_abort(
        "Empirical `baseline_time` anchoring is only supported for tidy SOP ",
        "outputs."
      )
    }
    time_names <- dnames[[state_dim - 1L]]
    if (is.null(time_names)) {
      time_names <- seq_len(dims[state_dim - 1L])
    }
    time_map <- standardize_time_map(time_map)
    real_time <- map_sop_time_values(time_names, time_map)

    if (ndim == 3) {
      return(apply(prob_in_target, 1, function(z) trapezoid_auc(real_time, z)))
    }
    return(apply(prob_in_target, c(1, 2), function(z) {
      trapezoid_auc(real_time, z)
    }))
  }

  if (ndim == 3) {
    # Freq Input -> [Pat, Time] -> Apply over rows
    res <- rowSums(prob_in_target)
  } else {
    # Bayes Input -> [Draws, Pat, Time] -> Apply over Draws & Pat
    res <- apply(prob_in_target, c(1, 2), sum)
  }

  return(res)
}
