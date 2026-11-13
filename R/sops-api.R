# SOP public APIs and BLRM summaries.

#' Calculate Individual State Occupation Probabilities
#'
#' Computes state occupation probabilities (SOPs) for each prediction profile,
#' by default each fitted patient's starting state and covariates. With
#' `variables`, each profile is predicted once for every requested value, as
#' in marginaleffects' `predictions()`. [avg_sops()] averages these SOPs over
#' the profiles.
#'
#' For supported frequentist models, pass the result to
#' `inferences(method = "delta", vcov = "conditional")` for analytical
#' confidence intervals. These account for coefficient estimation while
#' treating the starting states and covariates used for prediction as given.
#' `sops()` does not support unconditional analytical variance. See [inferences()]
#' for model restrictions and other inference methods.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, or `blrm`). For
#'   `vglm` models, the family must be cumulative-logit with `reverse = TRUE`.
#'   For `orm` models, the family must be logistic.
#' @param newdata Optional. A data frame of prediction profiles. When supplied,
#'   every row supplies a separate starting state and set of covariates. The
#'   internal `rowid` column is regenerated. If `NULL`, uses the data stored by
#'   [mostm()] and requires exactly one complete designated starting profile
#'   per fitted patient. Every column the model uses needs a value in every
#'   row, except the time and `gap_var` columns, which are set at each visit,
#'   and columns that `variables` sets; other columns may have missing values.
#'   The previous-state column holds states of the fitted model, as a factor
#'   (with any levels), numbers, or text.
#' @param times Required. The consecutive visits at which to calculate SOPs,
#'   in increasing order, for example `times = 1:30` for daily visits up to day
#'   30. Each value is one transition, so leaving out visits does not skip
#'   time; see the section "Choosing visits". Numeric time variables use
#'   numeric values; factor-valued visit indices use fitted visit levels.
#' @param variables Optional covariates to set before predicting: a named list
#'   of values, such as `list(tx = c(0, 1))`, or bare names, such as `"tx"`.
#'   Each profile is then predicted once for every combination of the values,
#'   with those columns set to the combination. For the values a bare name
#'   selects, see [avg_sops()]. The default, `NULL`, predicts each profile with
#'   its own covariates.
#' @param by How to average the profiles' SOPs, as in the marginaleffects
#'   package. `FALSE` (the default) or `NULL` returns one row per profile,
#'   time, and state. `TRUE` averages over all profiles, and a character
#'   vector of column names averages within the groups those columns define;
#'   both return the result of [avg_sops()] with the same arguments. A `by`
#'   column must have no missing values among the averaged profiles; drop or
#'   impute those patients first.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs. Only
#'   `blrm` fits have random effects, so `TRUE` gives an error for other fits.
#'   `newdata` must then contain the patient ID column (`id_var` of [mostm()]),
#'   so that each row gets its patient's estimated random effect.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible random draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior SOP draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior SOP draws for
#'   extraction with `get_draws()`. Large draw objects are guarded to protect
#'   memory.
#' @param ... Additional arguments (currently unused). The column names and
#'   states recorded by [mostm()] can't be given here; see the section
#'   "Settings taken from the fit".
#'
#' @return A data frame of class `markov_sops`, or, with `by = TRUE` or
#'   grouping columns, the `markov_avg_sops` result of [avg_sops()]. The
#'   per-profile result contains:
#'   \item{rowid}{Row number of the prediction: the profile, or with
#'     `variables`, the profile and value combination.}
#'   \item{rowidcf}{With `variables`, the profile that the row predicts.}
#'   \item{time}{Time point}
#'   \item{state}{State name}
#'   \item{estimate}{Probability that the profile is in the state.}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     uncertainty summaries computed directly from SOP draws}
#'   Plus all profile columns, with the `variables` columns holding the value
#'   used for the row.
#'
#' @details
#' This function wraps `soprob_markov()` and converts its array output to a tidy
#' data frame. The output contains one row per profile (and value of
#' `variables`), time, and state.
#'
#' Automatic prediction from a model fitted with [mostm()] uses the cohort-wide
#' designated `first_followup_time` retained before response-driven model-frame
#' omission. The transition response is not required on that row, but ID,
#' modeled predictors, time information, and previous state must be complete.
#' Each included patient must contribute at least one usable fitted transition
#' somewhere; later rows are not substituted for a missing starting profile.
#'
#' For `rmsb::blrm()` models, SOPs are computed on sampled posterior draws and
#' then summarized. The point estimate is the requested posterior summary of the
#' SOP draws, not a plug-in calculation from summarized model parameters.
#' State-wise medians and interval bounds are not constrained to sum to one
#' across states; use `posterior_summary = "mean"` when the displayed estimates
#' themselves need to preserve total probability. Draw-level probabilities
#' stored with `return_draws = TRUE` remain normalized within each draw.
#'
#' **Model Requirements:**
#'
#' For `vglm` models, only the cumulative-logit family with `reverse = TRUE` is
#' supported. For `orm` models, only the logistic family is supported. This is
#' because the package's Markov simulation logic expects
#' higher-numbered states to represent worse outcomes, and uses reverse cumulative
#' probabilities to model the probability of being in state k or worse.
#'
#' **Averages:**
#'
#' The average of the profiles' SOPs over all profiles (`by = TRUE`), or within
#' groups such as ECOG levels (`by = "ecog"`), is the [avg_sops()] result, so
#' `sops()` returns it directly. With `variables`, the average is the
#' standardized (G-computation) probability under each value.
#'
#' @section Choosing visits:
#' SOPs are calculated one visit at a time. The first value of `times` is the
#' first transition from the starting state, and each later value is one more
#' transition. The values set the time variable, and so every time term of the
#' model formula, such as `rms::rcs(time, 4)` or `log(time)`, for each
#' transition; they do not say how much time passes between transitions.
#' Leaving out visits therefore does not skip time. With daily visits,
#' `times = c(1, 7, 14)` gives the probabilities after 1, 2, and 3
#' transitions, labelled 1, 7, and 14, instead of the probabilities on days 7
#' and 14. To get SOPs at selected visits, request every visit, such as
#' `times = 1:14`, and keep the rows you need.
#'
#' A model fitted with `gap_var` (see [mostm()]) includes the time since the
#' previous visit. Each transition then gets the difference between its
#' requested time and the one before it, and the first transition gets the
#' first requested time itself, the time since time 0. With
#' `times = c(1, 2, 4, 7)` the gaps are 1, 1, 2 and 3, so the visits may be
#' spaced unevenly. When `newdata` starts at a later visit, the first gap is
#' still the first requested time, not the time since that visit.
#'
#' When `newdata` is omitted, each patient starts from their state before the
#' first follow-up visit (`first_followup_time` of [mostm()]), so
#' `times` must begin at that visit. Each row of a supplied `newdata` starts
#' from its previous state at the first requested visit, which may be a later
#' visit.
#'
#' A warning is given when `times` leaves out visits that occur in the data the
#' model was fitted to, is not in increasing order, or, when `newdata` is
#' omitted, does not begin at the first follow-up visit. Times after the last
#' fitted visit are allowed. There is no warning for a model fitted with
#' `gap_var`, because each transition then uses its own elapsed time, so the
#' visits may be spaced unevenly.
#'
#' @section Settings taken from the fit:
#' The patient ID column, the time column, the previous-state column, the
#' column with the state two visits back, the column with the time since the
#' previous visit, the states and the absorbing states are those of the
#' fitted model: the `id_var`, `time_var`, `p_var`, `p2_var`, `gap_var` and
#' `absorb` given to [mostm()], and the states of the response in fitted
#' order. A model fitted with `p2_var` is therefore always predicted with the
#' second-order recursion. The `state` column of the result holds the fitted
#' state labels. To predict with other columns or other absorbing states,
#' refit the model with [mostm()].
#'
#' A patient who is in an absorbing state stays there, so every profile must
#' start outside the absorbing states: a profile whose previous state (the
#' `p_var` column) is absorbing, or, for a second-order model, whose state two
#' visits back (the `p2_var` column) is, gives an error. Remove those patients
#' from `newdata` first.
#'
#' @seealso [avg_sops()] for averaged SOPs, [soprob_markov()] for the
#'   underlying computation.
#'
#' @examples
#' \dontrun{
#' # SOPs for each profile
#' sops_ind <- sops(fit, newdata = baseline_data, times = 1:30)
#'
#' # Each profile's SOPs under control and under treatment
#' sops_cf <- sops(fit, newdata = baseline_data, variables = list(tx = c(0, 1)),
#'                 times = 1:30)
#'
#' # Averaged within ECOG levels; the same as avg_sops(by = "ecog")
#' sops_ecog <- sops(fit, newdata = baseline_data, times = 1:30, by = "ecog")
#'
#' # Compatible with inferences()
#' sops_cf |> inferences(method = "mvn", n_draws = 1000)
#' }
#'
#' @importFrom stats model.frame
#' @export
sops <- function(
  model,
  newdata = NULL,
  times,
  variables = NULL,
  by = FALSE,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  # --- 1. Setup & Defaults ---
  stop_removed_fit_arguments(...)
  # Validate model compatibility
  validate_markov_model(model)
  stop_include_re_without_blrm(model, include_re)
  if (missing(times) || is.null(times)) {
    mostr_abort("`times` must be supplied to `sops()`.")
  }
  conf_level <- validate_conf_level(conf_level)
  settings <- resolve_markov_settings(model)
  by_choice <- resolve_by_choice(by)
  if (by_choice$average) {
    return(avg_sops_impl(
      model = model,
      newdata = newdata,
      variables = variables,
      by = by_choice$by,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      ...
    ))
  }
  sops_impl(
    model = model,
    newdata = newdata,
    variables = variables,
    times = times,
    settings = settings,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = return_draws,
    ...
  )
}

# `sops()` once its settings are resolved: one row per profile, time, and
# state. Averaging is left to `avg_sops_impl()`.
sops_impl <- function(
  model,
  newdata,
  variables,
  times,
  settings,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  variables <- resolve_prediction_variables(
    variables,
    model,
    newdata,
    settings
  )

  data_res <- resolve_markov_source_data(model, newdata, settings)
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied
  id_var <- markov_require_id_var(model, settings$id_var)
  if (inherits(model, "blrm") && isTRUE(include_re)) {
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  }
  if (!is.null(refit_data)) {
    validate_markov_id_var(
      id_var,
      refit_data,
      "the refit data stored on the fitted model"
    )
  }
  if (!newdata_supplied) {
    validate_markov_id_var(id_var, newdata_orig, "stored model data")
  }

  # Validate and coerce supplied times before extracting one prediction row per ID.
  time_res <- resolve_sop_times(
    model,
    newdata_orig,
    times,
    settings$time_var,
    default = "unique"
  )
  times <- time_res$times
  warn_skipped_sop_visits(
    model,
    newdata_orig,
    times,
    settings,
    newdata_supplied
  )

  newdata <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = settings$time_var
    )
  }
  # Every profile column is copied into the result, and the prediction rows
  # are stored for `inferences()`, so bookkeeping columns are removed here
  # once for the estimates, every replay, and the draws.
  newdata <- drop_markov_internal_columns(newdata)
  newdata <- ensure_markov_rowid(newdata)
  if (length(variables) > 0L) {
    newdata <- stack_prediction_profiles(newdata, variables)
  }

  if (is.null(settings$y_levels)) {
    stop_no_fitted_states()
  }

  if (inherits(model, "blrm")) {
    posterior_summary <- match.arg(posterior_summary)
    return(sops_blrm(
      model = model,
      newdata = newdata,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      newdata_orig = newdata_orig,
      refit_data = refit_data,
      newdata_supplied = newdata_supplied
    ))
  }

  # --- 2. Compute SOPs (Vectorized) ---
  sops_array <- predict_markov_sops(
    model = model,
    newdata = newdata,
    times = times,
    settings = settings
  )

  # --- 3. Tidy the Output ---
  result <- array_to_df_individual(
    sops_array,
    times,
    settings$y_levels,
    newdata
  )

  # --- 4. Store Attributes for Downstream Use ---
  set_sops_attrs(
    result,
    class_name = "markov_sops",
    model = model,
    call_args = sops_call_args(times),
    settings = settings,
    newdata_orig = newdata_orig,
    extra_attrs = list(
      newdata_pred = newdata,
      refit_data = refit_data,
      newdata_supplied = newdata_supplied
    )
  )
}

sops_blrm <- function(
  model,
  newdata,
  times,
  settings,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  newdata_orig,
  refit_data,
  newdata_supplied
) {
  y_levels <- settings$y_levels
  draw_indices <- select_posterior_draws(model, n_draws, seed)
  n_pat <- nrow(newdata)
  n_times <- length(times)
  n_states <- length(y_levels)
  n_cells <- n_pat * n_times * n_states
  n_requested <- length(draw_indices)
  # `sops()` keeps individual-level draw summaries in memory. The default
  # guard is 50 million numeric cells, about 400 MB before data-frame overhead.
  max_draw_cells <- getOption("mostr.max_sops_draw_cells", 5e7)

  if (n_cells * n_requested > max_draw_cells) {
    mostr_abort(
      "`sops()` for `blrm` would require ",
      format(n_cells * n_requested, scientific = FALSE),
      " posterior draw cells, above the current limit of ",
      format(max_draw_cells, scientific = FALSE),
      ". This limit protects memory because individual-level `sops()` stores ",
      "draw x patient x time x state values. Lower `n_draws`, use ",
      "`avg_sops()` for marginal summaries, or increase option ",
      "`mostr.max_sops_draw_cells` if this allocation is intentional."
    )
  }

  draw_values <- matrix(NA_real_, nrow = n_requested, ncol = n_cells)
  rownames(draw_values) <- draw_indices
  chunks <- split_draw_indices(
    draw_indices,
    chunk_size = getOption("mostr.blrm_chunk_size") %||%
      posterior_chunk_size(n_cells, n_requested)
  )
  gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  prediction_cache <- if (is.null(settings$p2_var)) {
    cache <- new.env(parent = emptyenv())
    cache$cursor <- 0L
    cache$X <- list()
    cache$Z <- list()
    cache
  } else {
    NULL
  }
  row_cursor <- 1L
  for (chunk in chunks) {
    if (!is.null(prediction_cache)) {
      prediction_cache$cursor <- 0L
    }
    arr <- predict_markov_sops(
      model = model,
      newdata = newdata,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = NULL,
      .draw_indices = chunk,
      .n_draws_total = length(draw_indices),
      .gamma_draws = subset_cached_draw_matrix(
        gamma_draws,
        draw_indices,
        chunk
      ),
      .prediction_cache = prediction_cache
    )
    for (i in seq_along(chunk)) {
      draw_values[row_cursor, ] <- as.vector(arr[i, , , ])
      row_cursor <- row_cursor + 1L
    }
  }

  stats <- summarize_posterior_draw_matrix(
    draw_values,
    posterior_summary = posterior_summary,
    conf_level = conf_level
  )
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  result <- newdata[idx_pat, , drop = FALSE]
  result$time <- rep(rep(times, each = n_pat), times = n_states)
  result$state <- rep(y_levels, each = n_pat * n_times)
  result$estimate <- stats$estimate
  result$conf.low <- stats$conf.low
  result$conf.high <- stats$conf.high
  result$std.error <- stats$std.error
  rownames(result) <- NULL

  draws_df <- if (return_draws) {
    sops_draw_matrix_to_df(draw_values, result, draw_indices)
  } else {
    NULL
  }

  set_sops_attrs(
    result,
    class_name = "markov_sops",
    model = model,
    call_args = sops_call_args(
      times,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    ),
    settings = settings,
    newdata_orig = newdata_orig,
    extra_attrs = list(
      method = "posterior",
      engine = "posterior",
      n_draws = n_requested,
      draw_ids = draw_indices,
      conf_level = conf_level,
      posterior_summary = posterior_summary,
      draws = draws_df,
      newdata_pred = newdata,
      refit_data = refit_data,
      newdata_supplied = newdata_supplied
    )
  )
}

posterior_group_setup <- function(data, group_cols) {
  rows <- which(stats::complete.cases(data[, group_cols, drop = FALSE]))
  if (!length(rows)) {
    mostr_abort(
      "No complete prediction rows remain after omitting missing grouping values."
    )
  }
  complete_data <- data[rows, , drop = FALSE]
  keys <- interaction(
    complete_data[, group_cols, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- as.integer(keys)
  n_groups <- nlevels(keys)
  first <- vapply(
    seq_len(n_groups),
    function(i) which(groups == i)[1L],
    integer(1)
  )
  metadata <- complete_data[first, group_cols, drop = FALSE]
  rownames(metadata) <- NULL
  list(
    rows = rows,
    groups = groups,
    n_groups = n_groups,
    metadata = metadata
  )
}

posterior_group_cells <- function(metadata, times, y_levels) {
  pieces <- vector("list", nrow(metadata))
  for (i in seq_len(nrow(metadata))) {
    piece <- expand.grid(
      time = times,
      state = y_levels,
      KEEP.OUT.ATTRS = FALSE
    )
    for (nm in names(metadata)) {
      piece[[nm]] <- metadata[[nm]][i]
    }
    pieces[[i]] <- piece
  }
  bind_rows_fill(pieces)
}

posterior_draw_matrix_to_df <- function(draw_values, cells, draw_indices) {
  n_draws <- nrow(draw_values)
  n_cells <- ncol(draw_values)
  out <- cells[rep(seq_len(n_cells), times = n_draws), , drop = FALSE]
  out$draw_id <- rep(draw_indices, each = n_cells)
  out$estimate <- as.vector(t(draw_values))
  rownames(out) <- NULL
  out
}

split_draw_indices <- function(draw_indices, chunk_size = NULL) {
  if (is.null(chunk_size)) {
    # Fallback only: package callers pass a chunk size (the option, or one sized
    # by mostr.posterior_working_bytes). Nov 2026, sops() with 200 profiles x 15
    # visits and 400 draws: 1, 10 and 50 draws per chunk took 18.5, 12.7 and
    # 11.7 s, so 10 is close to the best.
    chunk_size <- getOption("mostr.blrm_chunk_size", 10L)
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  split(draw_indices, ceiling(seq_along(draw_indices) / chunk_size))
}

summarize_posterior_draw_matrix <- function(
  draw_values,
  posterior_summary,
  conf_level
) {
  alpha <- 1 - conf_level
  estimate <- switch(
    posterior_summary,
    mean = colMeans(draw_values, na.rm = TRUE),
    median = apply(draw_values, 2, stats::median, na.rm = TRUE)
  )
  conf.low <- apply(
    draw_values,
    2,
    stats::quantile,
    probs = alpha / 2,
    na.rm = TRUE
  )
  conf.high <- apply(
    draw_values,
    2,
    stats::quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  std.error <- apply(draw_values, 2, stats::sd, na.rm = TRUE)
  data.frame(
    estimate = as.numeric(estimate),
    conf.low = as.numeric(conf.low),
    conf.high = as.numeric(conf.high),
    std.error = as.numeric(std.error)
  )
}

sops_draw_matrix_to_df <- function(draw_values, result, draw_indices) {
  meta_cols <- setdiff(
    names(result),
    c("estimate", "conf.low", "conf.high", "std.error")
  )
  n_draws <- nrow(draw_values)
  n_cells <- ncol(draw_values)
  out <- result[rep(seq_len(n_cells), times = n_draws), meta_cols, drop = FALSE]
  out$draw_id <- rep(draw_indices, each = n_cells)
  out$estimate <- as.vector(t(draw_values))
  rownames(out) <- NULL
  out
}


#' Calculate Averaged State Occupation Probabilities (Marginal Effects)
#'
#' Averages state occupation probabilities over patients' starting states and
#' covariates. With `variables = NULL`, uses each patient's observed covariates.
#' Supply `variables` to set the specified covariates to common values before
#' averaging, for example to estimate probabilities under each treatment.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, or `blrm`). For
#'   `vglm` models, the family must be cumulative-logit with `reverse = TRUE`.
#'   For `orm` models, the family must be logistic.
#' @param newdata Optional data frame of standardization profiles. When
#'   supplied, every row supplies a separate starting state and set of
#'   covariates. The internal `rowid` column is regenerated. If `NULL`, uses
#'   data stored by [mostm()] and requires exactly one complete designated
#'   starting profile per fitted patient. Every column the model uses needs a
#'   value in every row, except the time and `gap_var` columns, which are set
#'   at each visit, and columns that `variables` sets; other columns may have
#'   missing values.
#'   The previous-state column holds states of the fitted model, as a factor
#'   (with any levels), numbers, or text.
#' @param variables Optional covariates to set before averaging: a named list
#'   of values, such as `list(tx = c(0, 1))` for probabilities if every patient
#'   received control or treatment, or bare names, such as `"tx"`. As in
#'   marginaleffects' `predictions()`, a bare name uses every observed value:
#'   factor levels in level order, character values in sorted order, `FALSE`
#'   and `TRUE`, or `0` and `1` for a 0/1 numeric variable. Other numeric
#'   variables use their five-number summary (minimum, lower quartile,
#'   median, upper quartile, maximum). Values are read from the data stored
#'   with the model, or from `newdata` for a model without stored data. The
#'   default, `NULL`, averages with each patient's observed covariates,
#'   including their observed treatment.
#' @param by How to average, as in the marginaleffects package. `TRUE` (the
#'   default) averages over all profiles. A character vector of column names
#'   averages within the groups those columns define after standardization.
#'   A `by` column must have no missing values among the averaged profiles,
#'   because a profile without a group would be left out of every average;
#'   drop or impute those patients first, in `newdata` or in the data given to
#'   [mostm()]. Without `variables`, groups use each profile's observed
#'   values, so `by = "tx"` averages separately among the untreated and the
#'   treated. `by` cannot name a variable that `variables` also sets: the
#'   result's column for that variable holds the value that was set, so it
#'   cannot also hold groups of the observed value, and this is an error. The
#'   averages for each set value are already returned without `by`. To
#'   average under each set treatment within groups of the observed
#'   treatment, copy the column first, for example
#'   `newdata$observed_tx <- newdata$tx`, and use `by = "observed_tx"`.
#'   [avg_comparisons()] with `by = "tx"` does this for treatment effects.
#'   `FALSE` or `NULL` does not average and returns the [sops()] result, one
#'   row per profile, value of `variables`, time, and state.
#' @inheritParams sops
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible random draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior SOP draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior SOP draws for
#'   extraction with `get_draws()`.
#' @param ... Additional arguments reserved for future extensions. The column
#'   names and states recorded by [mostm()] can't be given here; see the
#'   section "Settings taken from the fit".
#'
#' @return With `by = FALSE`, the `markov_sops` result of [sops()]. Otherwise
#'   a data frame of class `markov_avg_sops` with columns:
#'   \item{time}{Time point}
#'   \item{state}{State level}
#'   \item{(variables)}{Values requested in `variables`, when supplied}
#'   \item{estimate}{Average probability across individuals}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     uncertainty summaries computed directly from marginalized SOP draws}
#'
#' @details
#' With `variables = NULL`, this function averages patients' SOPs using their
#' observed covariates, returning one estimate per time and state. Supply `by`
#' to average separately within groups.
#'
#' Supplying `variables` implements G-computation (standardization):
#'
#' 1. **Counterfactual Creation**: For each value in `variables`, creates a
#'    copy of `newdata` with that variable set to the specified value.
#'
#' 2. **Individual Prediction**: Calls `sops()` to compute individual-level
#'    SOPs for each patient under each counterfactual scenario.
#'
#' 3. **Marginalization**: Averages individual SOPs across patients within
#'    each time-state-treatment combination, yielding population-average
#'    (marginal) probabilities.
#'
#' With `variables = list(tx = c(0, 1))`, the result represents the expected SOP
#' if the entire population received treatment or control, averaged over the
#' observed covariate distribution.
#' When `newdata` is omitted, the average uses the starting profiles stored by
#' [mostm()]: exactly the fitted patients with a complete designated starting
#' profile and
#' at least one usable likelihood transition. A missing response at the
#' designated row alone does not exclude a patient who contributes a later
#' transition. Patients with no usable fitted transition are excluded, and a
#' fitted patient without a complete designated profile is an error.
#'
#' The average gives every profile the same weight. Case weights given to
#' [mostm()] with `weights` are used only to fit the model, not in this
#' average, so the result is not the same as repeating each patient according
#' to their weight.
#'
#' For `rmsb::blrm()` models, patient-level posterior arrays are chunked and
#' marginalized by draw before summarizing to reduce memory use.
#' State-wise medians and interval bounds are not constrained to sum to one
#' across states; use `posterior_summary = "mean"` when the displayed estimates
#' themselves need to preserve total probability. Draw-level probabilities
#' stored with `return_draws = TRUE` remain normalized within each draw.
#'
#' For supported frequentist models, `inferences(method = "delta")` adds
#' analytical confidence intervals. Its default `vcov = "unconditional"`
#' accounts for coefficient estimation and sampling the patients used in the
#' average. Use `vcov = "conditional"` to treat those patients' starting states
#' and covariates as given and account only for coefficient estimation.
#' User-supplied `newdata`, and models fitted with case weights, require the
#' conditional choice. See [inferences()] for model restrictions.
#'
#' @inheritSection sops Choosing visits
#' @inheritSection sops Settings taken from the fit
#'
#' @seealso [sops()] for individual-level SOPs, [inferences()] for analytical,
#'   simulation, or bootstrap
#'   uncertainty and [avg_comparisons()] for specialized estimands.
#'
#' @examples
#' \dontrun{
#' fit <- mostm(
#'   ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
#'   data = data,
#'   id_var = "id",
#'   absorb = 6,
#'   engine = "vglm"
#' )
#'
#' # Models fitted with mostm() can use their stored starting profiles.
#' observed <- avg_sops(fit, times = 1:60)
#'
#' # Set treatment to each requested value before averaging.
#' result <- avg_sops(
#'   model = fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60
#' ) |> inferences(method = "mvn", n_draws = 500)
#'
#' # Refit bootstrap resamples the patients stored on the fitted model.
#' result_boot <- avg_sops(
#'   model = fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60
#' ) |> inferences(method = "fwb", n_draws = 500)
#' }
#'
#' @export
avg_sops <- function(
  model,
  newdata = NULL,
  variables = NULL,
  by = TRUE,
  times,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  # --- 1. Input Validation ---
  stop_removed_fit_arguments(...)
  # Validate model compatibility
  validate_markov_model(model)
  stop_include_re_without_blrm(model, include_re)
  if (missing(times) || is.null(times)) {
    mostr_abort("`times` must be supplied to `avg_sops()`.")
  }
  conf_level <- validate_conf_level(conf_level)
  settings <- resolve_markov_settings(model)

  by_choice <- resolve_by_choice(by)
  if (!by_choice$average) {
    return(sops_impl(
      model = model,
      newdata = newdata,
      variables = variables,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      ...
    ))
  }
  avg_sops_impl(
    model = model,
    newdata = newdata,
    variables = variables,
    by = by_choice$by,
    times = times,
    settings = settings,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = return_draws,
    ...
  )
}

# `avg_sops()` once its settings are resolved. `by` holds the grouping
# columns, or `NULL` to average over all profiles. Inference replays call this
# with the settings stored on a result.
avg_sops_impl <- function(
  model,
  newdata,
  variables,
  by,
  times,
  settings,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  var_list <- resolve_prediction_variables(
    variables,
    model,
    newdata,
    settings
  )
  stop_by_names_set_variable(by, var_list, "avg_sops()")

  data_res <- resolve_markov_source_data(model, newdata, settings)
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied

  id_var <- markov_require_id_var(model, settings$id_var)
  if (inherits(model, "blrm") && isTRUE(include_re)) {
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  }

  if (!is.null(refit_data)) {
    validate_markov_id_var(
      id_var,
      refit_data,
      "the refit data stored on the fitted model"
    )
  }

  if (!newdata_supplied) {
    validate_markov_id_var(id_var, newdata_orig, "stored model data")
  }

  # Validate variables exist in data
  missing_vars <- setdiff(names(var_list), names(newdata_orig))
  if (length(missing_vars) > 0) {
    mostr_abort("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- 2. Extract Baseline Data (One Row Per Patient) ---
  # For standardization, we need unique patient baseline covariates
  # Average the counterfactual SOP predictions over the target population.
  baseline_data <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = settings$time_var
    )
  }
  stop_missing_by_values(baseline_data, by, var_list)
  baseline_data <- ensure_markov_rowid(baseline_data)

  # --- 3. Create Counterfactual Datasets ---
  # For each combination in variables, create a copy of baseline_data with
  # the variable(s) set to that value

  grid <- create_counterfactual_grid(var_list)

  newdata_expanded <- create_counterfactual_data(
    baseline_data,
    grid,
    var_list,
    by = by
  )
  warn_skipped_sop_visits(
    model,
    newdata_orig,
    times,
    settings,
    newdata_supplied
  )

  if (is.null(settings$y_levels)) {
    stop_no_fitted_states()
  }

  if (inherits(model, "blrm")) {
    posterior_summary <- match.arg(posterior_summary)
    result <- avg_sops_blrm(
      model = model,
      newdata_expanded = newdata_expanded,
      grid = grid,
      variables = var_list,
      by = by,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      ...
    )
    attr(result, "newdata_orig") <- newdata_orig
    attr(result, "newdata_pred") <- newdata_expanded
    attr(result, "refit_data") <- refit_data
    attr(result, "newdata_supplied") <- newdata_supplied
    return(result)
  }

  # --- 4. Compute and marginalize without materializing individual rows ---
  time_res <- resolve_sop_times(
    model,
    newdata_expanded,
    times,
    settings$time_var,
    default = "unique"
  )
  resolved_times <- time_res$times
  sops_array <- predict_markov_sops(
    model = model,
    newdata = newdata_expanded,
    times = resolved_times,
    settings = settings,
    ...
  )
  result <- marginalize_sops_array(
    sops_array = sops_array,
    grid = grid,
    times = resolved_times,
    y_levels = settings$y_levels,
    variables = var_list,
    n_cf = nrow(grid),
    n_each = nrow(baseline_data),
    by = by,
    newdata = newdata_expanded
  )

  set_sops_attrs(
    result,
    class_name = "markov_avg_sops",
    model = model,
    call_args = sops_call_args(resolved_times),
    settings = settings,
    newdata_orig = newdata_orig,
    avg_args = list(
      variables = var_list,
      by = by,
      times = resolved_times
    ),
    extra_attrs = list(
      newdata_pred = newdata_expanded,
      refit_data = refit_data,
      newdata_supplied = newdata_supplied
    )
  )
}

avg_sops_blrm <- function(
  model,
  newdata_expanded,
  grid,
  variables,
  by,
  times,
  settings,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  ...
) {
  y_levels <- settings$y_levels

  time_res <- resolve_sop_times(
    model,
    newdata_expanded,
    times,
    settings$time_var,
    default = "unique"
  )
  times <- time_res$times

  validate_sops_by(by, newdata_expanded, data_arg = "expanded prediction data")

  n_cf <- nrow(grid)
  n_each <- nrow(newdata_expanded) / n_cf
  group_cols <- unique(c(names(variables), by))
  group_setup <- if (is.null(by)) {
    list(
      rows = seq_len(nrow(newdata_expanded)),
      groups = rep(seq_len(n_cf), each = n_each),
      n_groups = n_cf,
      metadata = grid[, names(variables), drop = FALSE]
    )
  } else {
    posterior_group_setup(newdata_expanded, group_cols)
  }
  prediction_data <- newdata_expanded[group_setup$rows, , drop = FALSE]
  draw_indices <- select_posterior_draws(model, n_draws, seed)
  chunks <- split_draw_indices(
    draw_indices,
    chunk_size = getOption("mostr.blrm_avg_chunk_size") %||%
      posterior_chunk_size(
        nrow(prediction_data) * length(times) * length(y_levels),
        length(draw_indices)
      )
  )
  gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  n_cells <- group_setup$n_groups * length(times) * length(y_levels)
  draw_values <- matrix(
    NA_real_,
    nrow = length(draw_indices),
    ncol = n_cells,
    dimnames = list(draw_indices, NULL)
  )
  prediction_cache <- if (is.null(settings$p2_var)) {
    cache <- new.env(parent = emptyenv())
    cache$cursor <- 0L
    cache$X <- list()
    cache$Z <- list()
    cache
  } else {
    NULL
  }
  row_cursor <- 1L

  for (chunk in chunks) {
    if (!is.null(prediction_cache)) {
      prediction_cache$cursor <- 0L
    }
    arr <- predict_markov_sops(
      model = model,
      newdata = prediction_data,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = NULL,
      .draw_indices = chunk,
      .n_draws_total = length(draw_indices),
      .gamma_draws = subset_cached_draw_matrix(
        gamma_draws,
        draw_indices,
        chunk
      ),
      .prediction_cache = prediction_cache
    )
    reduced <- reduce_sops_draw_array(
      arr,
      group_setup$groups,
      group_setup$n_groups
    )
    reduced <- aperm(reduced, c(1L, 3L, 4L, 2L))
    rows <- row_cursor:(row_cursor + length(chunk) - 1L)
    draw_values[rows, ] <- matrix(as.numeric(reduced), nrow = length(chunk))
    row_cursor <- row_cursor + length(chunk)
  }

  cells <- posterior_group_cells(group_setup$metadata, times, y_levels)
  stats <- summarize_posterior_draw_matrix(
    draw_values,
    posterior_summary = posterior_summary,
    conf_level = conf_level
  )
  result <- cbind(cells, stats)
  draws_df <- if (return_draws) {
    posterior_draw_matrix_to_df(draw_values, cells, draw_indices)
  } else {
    NULL
  }

  set_sops_attrs(
    result,
    class_name = "markov_avg_sops",
    model = model,
    call_args = sops_call_args(
      times,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    ),
    settings = settings,
    avg_args = list(
      variables = variables,
      by = by,
      times = times
    ),
    extra_attrs = list(
      method = "posterior",
      engine = "posterior",
      n_draws = length(draw_indices),
      draw_ids = draw_indices,
      conf_level = conf_level,
      posterior_summary = posterior_summary,
      draws = draws_df
    )
  )
}

# Draws per chunk = budget / (prediction rows x visits x states x 8 x 3); most
# analyses fit all draws in one chunk. 256 MiB was a choice. Nov 2026 (fake
# blrm, 2 x 1000 profiles x 29 visits, 100 draws): chunks of 8, 32 and 100 draws
# took 14.2, 14.5 and 14.9 s with 0.17, 0.42 and 0.80 GB, so larger chunks cost
# memory and save no time.
posterior_chunk_size <- function(cells_per_draw, n_draws) {
  budget <- getOption("mostr.posterior_working_bytes", 256 * 1024^2)
  bytes_per_draw <- max(1, cells_per_draw * 8 * 3)
  max(1L, min(as.integer(n_draws), floor(budget / bytes_per_draw)))
}
