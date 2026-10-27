# Fitting entry point for Markov ordinal state transition models.

#' Fit a Markov ordinal state transition model
#'
#' `mostm()` fits a model for the transitions between ordinal states, such as
#' daily health states of hospitalized patients, and stores the data that
#' [sops()], [avg_sops()], [state_time()], [avg_comparisons()], [inferences()],
#' and the plotting functions need. Each row of `data` is one transition: the
#' response is a patient's state at one visit, and the predictors usually
#' include the state at the previous visit (`yprev`), time, and treatment.
#'
#' The model is a cumulative logit (proportional odds) model. `engine`
#' chooses the package that fits it:
#'
#' * `"orm"` (the default) fits by maximum likelihood with [rms::orm()].
#' * `"vglm"` fits by maximum likelihood with [VGAM::vglm()]. Use it for a
#'   partial proportional odds model, in which some predictors have a separate
#'   effect at each threshold between states: supply
#'   `family = VGAM::cumulative(reverse = TRUE, parallel = FALSE)` together with
#'   `constraints`. Without `family`, this engine fits the same proportional
#'   odds model as `"orm"`, so changing `engine` is enough to compare the two.
#' * `"blrm"` fits a Bayesian model with [rmsb::blrm()], which samples from the
#'   posterior distribution with Stan. Partial proportional odds effects are
#'   set with `ppo` and `cppo`. This engine needs the rmsb package.
#'
#' @param formula A model formula with the current state as the response, for
#'   example `ordered(y) ~ yprev + rms::rcs(time, 4) + tx`. Offsets are not
#'   supported.
#' @param data A data frame with one row per patient and visit.
#' @param ... Arguments for the engine, given by name. Each must be an
#'   argument of the chosen engine; an argument that belongs to a different
#'   engine is an error that names the engine it belongs to.
#'   * `"orm"`: arguments of [rms::orm()] and [rms::orm.fit()], such as
#'     `penalty`, `maxit`, or `na.action`.
#'   * `"vglm"`: arguments of [VGAM::vglm()], such as `family`, `constraints`,
#'     `coefstart`, or `control`, and of [VGAM::vglm.control()], such as
#'     `maxit` or `trace`. `family` defaults to
#'     `VGAM::cumulative(reverse = TRUE, parallel = TRUE)`.
#'   * `"blrm"`: arguments of [rmsb::blrm()], such as `ppo`, `cppo`, `iter`,
#'     `warmup`, `chains`, the prior settings `conc` and `psigma`, or
#'     `backend`, which chooses between rstan and cmdstan. Settings of the
#'     Stan sampler that [rmsb::blrm()] has no argument for go in
#'     `sampling.args`: for example, `sampling.args = list(seed = 123)` makes
#'     the posterior draws reproducible, and `cores = 4` in the same list runs
#'     four chains in parallel. Without `cores` or `options(mc.cores = ...)`,
#'     the chains run one after another. Other arguments are passed on by
#'     [rmsb::blrm()] to the Stan optimizer only, which finds the posterior
#'     mode before sampling, so a `seed` given directly to `mostm()` does not
#'     make the draws reproducible.
#'
#'   The design matrix and response are always kept on the fit, because the
#'   state occupancy and inference functions need them, so the engine
#'   arguments that control this (`x` and `y` for rms and rmsb, `x.arg` and
#'   `y.arg` for VGAM) are not accepted.
#' @param engine The package that fits the model: `"orm"`, `"vglm"`, or
#'   `"blrm"`. See the description above.
#' @param id_var Name of the column in `data` that identifies patients. With
#'   it, `mostm()` stores each patient's starting row for automatic
#'   predictions, and the `"orm"` and `"vglm"` engines compute standard errors
#'   that allow for correlation between the transitions of one patient. Without
#'   it, a warning is given, and these features, as well as the fractional
#'   weighted bootstrap, are not available from the fit.
#' @param time_var Name of the time column in `data`. It is used to find each
#'   patient's starting row and, with `missing_states` or an [rms::Ocens()]
#'   response, to order each patient's visits. Time must then be numeric or a
#'   factor whose levels are in time order: character time is an error,
#'   because it would sort as text, putting `"day 10"` before `"day 2"`.
#' @param p_var Name of the previous-state column in `data`, as in [sops()].
#'   It is used with `missing_states` or an [rms::Ocens()] response, which
#'   need to know which predictor is the previous state, and to describe
#'   rows that are not used because of missing values.
#' @param first_followup_time Time of the first visit after baseline. Each
#'   patient's row at this time is their starting row for automatic
#'   predictions. For numeric time, `NULL` (the default) uses 1, time may not
#'   be below 1, and a value you supply must be the earliest time in the data.
#'   Factor and character time need an explicit value, such as `"day 1"`.
#' @param missing_states The states a patient may be in at a visit whose state
#'   was not recorded (`NA` in the response), for example `1:7` when state 8 is
#'   death and deaths are always recorded. With it, visits with a missing state
#'   are used: the likelihood adds up the probabilities of these states
#'   instead of dropping the transitions into and out of the visit. `NULL`
#'   (the default) drops rows with a missing response, with a warning. See
#'   "Missing and partially known states" below.
#' @param weights Optional non-negative row weights: a column of `data` or a
#'   vector with one value per row of `data`. Not available for
#'   `engine = "blrm"`.
#' @param subset Optional selection of the rows of `data` to fit, as a logical
#'   expression evaluated in `data`, for example `time <= 28`.
#' @param type Small-sample correction of the patient-clustered standard
#'   errors of the `"orm"` and `"vglm"` engines. `"HC0"` (the default) applies
#'   none. `"HC1"` multiplies the variance by \eqn{(m - 1) / (m - p)}, where
#'   \eqn{m} is the number of fitted rows with positive weight and \eqn{p} is
#'   the number of model coefficients.
#' @param cadjust Logical. Multiply the patient-clustered variance of the
#'   `"orm"` and `"vglm"` engines by \eqn{n / (n - 1)}, where \eqn{n} is the
#'   number of patients? `NULL` (the default) applies this correction when
#'   `id_var` is supplied. It is independent of `type`, so `type = "HC1"` with
#'   `cadjust = TRUE` applies both factors.
#'
#' @details
#' ## Starting rows
#'
#' Functions such as [sops()] and [avg_sops()] predict, by default, from each
#' fitted patient's row at `first_followup_time`. `mostm()` keeps these rows
#' before rows with a missing response are dropped, so a patient whose first
#' outcome is missing still has a starting row, as long as the patient has at
#' least one fitted transition. The starting row needs the patient ID, every
#' predictor of the model, the time, and the previous state; the response may
#' be missing. A fitted patient without a complete starting row is an error
#' when predictions are made: a later row is never used instead. Patients
#' without any fitted transition are not part of the fitted cohort.
#'
#' ## Weights and subsets
#'
#' `weights` and `subset` are evaluated once, first in `data` and then where
#' `mostm()` was called, so they can refer to vectors outside `data`. The
#' weights must be finite and non-negative. The data stored on the fit contain
#' only the rows selected by `subset`, and each row keeps its weight, so
#' bootstrap refits resample rows together with their weights.
#'
#' ## Standard errors
#'
#' For the `"orm"` and `"vglm"` engines, `id_var` gives standard errors that
#' remain valid when the transitions of one patient are correlated. They come
#' from a sandwich estimator that adds up the contributions of each patient
#' (cluster-robust standard errors); `type` and `cadjust` set its small-sample
#' corrections. These standard errors are used by `inferences(method =
#' "delta")` with conditional variance and by `inferences(method = "mvn")`.
#' Unconditional variance also accounts for which patients were sampled; see
#' [inferences()]. For the `"blrm"` engine, uncertainty comes from the
#' posterior draws.
#'
#' ## Missing and partially known states
#'
#' By default, `mostm()` fits the rows without missing values. A visit whose
#' state is missing (`NA`) then removes two transitions: the one into the
#' visit (missing response) and the one out of it (missing previous state).
#' `mostm()` warns when rows are dropped. Dropping these transitions can bias
#' the estimates. For example, when deaths are always recorded but other
#' states are often missing, the remaining transitions overrepresent deaths.
#'
#' Give `missing_states` to use these visits instead. Each patient's rows,
#' ordered by `time_var`, then form one chain, so time must be numeric or a
#' factor with its levels in time order; convert character time such as
#' `"day 1"` to such a factor or to numbers. For a visit with a missing
#' state, the likelihood adds up the probabilities of all paths through the
#' states in `missing_states`. Each patient contributes the probability of
#' everything observed for that patient, given the first previous state and
#' the covariates, and `mostm()` maximises the sum of their logarithms. When a
#' state is known only up to a range, such as "in hospital, level unknown",
#' code the response as `Ocens(lower, upper)` with [rms::Ocens()], for example
#' `Ocens(y_lo, y_hi) ~ yprev + rms::rcs(time, 4) + tx`. A row's range then
#' replaces `missing_states`, `Ocens(b, b)` is an observed state `b`, and a
#' row with both endpoints missing is a missing state. A range allows the
#' states in it that are observed exactly somewhere in the data, so a finite
#' endpoint may lie beyond the lowest or highest such state, but each range
#' must contain at least one of them. An `Ocens()` response
#' alone, without `missing_states`, uses the ranges and drops rows with a
#' missing state. Ranges and missing states work with all three engines,
#' because `mostm()` evaluates the response itself: the engine sees only
#' exactly observed states.
#'
#' With nothing missing or range-coded, `missing_states` changes nothing: the
#' fit is the usual fit of the engine. Otherwise:
#'
#' * Keep one row per scheduled visit until death or the end of follow-up,
#'   including the rows whose state is missing. A visit without a row is
#'   treated as not scheduled. Only the first row's previous state is used
#'   (it must be observed); afterwards, the previous state comes from the
#'   chain, so `p_var` may be missing or wrong after a missing state. Where
#'   the previous visit's state is observed, `p_var` must agree with it.
#' * A missing state cannot be a state that the data never show being left,
#'   such as death (`absorb` in [sops()]), unless the visit is the patient's
#'   last. If a patient may have died at a visit without a recorded state,
#'   the patient is lost to follow-up from then on: end the patient's rows at
#'   the last visit with a known state.
#' * The model must be first order: apart from `p_var`, every covariate must
#'   be fixed per patient (such as treatment) or a function of time only. A
#'   covariate computed from the state history, such as the state two visits
#'   earlier, would be wrong after a missing state. Rows with missing
#'   covariates are dropped, with a warning, and the patient's chain restarts
#'   at the next row whose previous state is observed.
#' * The estimates are valid when states are missing at random: whether a
#'   state is missing may depend on what was observed, such as earlier states
#'   or treatment, but not on the unobserved state itself.
#' * `id_var` is required, and `weights` cannot be used. Random effects
#'   (`cluster()` in `formula`) are available only with `engine = "blrm"`.
#'
#' For the `"orm"` and `"vglm"` engines, `mostm()` starts from the fit to the
#' transitions with both states observed and maximises the likelihood with
#' the EM algorithm followed by Newton steps. Spline knots in `rms::rcs()`
#' terms are set from these transitions. Standard errors use the exact
#' observed information (by Louis' identity) and the patients' scores, in the
#' same patient-clustered sandwich as above; `type = "HC1"` counts the rows
#' whose observation rules out at least one state. The returned object is the
#' engine's fit to the expected numbers of transitions at the estimates, which
#' predicts like any other fit of that engine. [logLik()], [AIC()], and [BIC()]
#' use the maximised log-likelihood of the observed data, and [nobs()]
#' counts the rows that rule out at least one state. Printing the fit shows
#' this log-likelihood, how many states were missing, and the coefficients
#' with their robust standard errors. All [inferences()] methods except the
#' fractional weighted bootstrap are available; bootstrap refits resample
#' whole patients, including their missing rows.
#'
#' For `engine = "blrm"`, `mostm()` samples the posterior of the
#' [rmsb::blrm()] model with the same priors (`conc`, `pcontrast`,
#' `npcontrast`, `psigma`, `rsdmean`, and `rsdsd`), using its own Stan program
#' with this likelihood. It needs the rstan package and a C++ toolchain; the
#' program is compiled on first use (which takes a few minutes) and saved in
#' `tools::R_user_dir("mostr", "cache")`. Set `options(mostr.stan_cache_dir =
#' FALSE)` to keep it only for the current session. Constrained partial
#' proportional odds (`ppo` with `cppo`) and a random intercept per patient
#' (`cluster()` with the patient ID column) are available; `sampling.args`
#' goes to [rstan::sampling()], for example `list(seed = 1, cores = 4)`.
#' There is no Stan optimizer run, so arguments that [rmsb::blrm()] would
#' pass to it, such as `seed`, are errors.
#' Check the convergence diagnostics in `$diagnostics` (`Rhat` below 1.01 and
#' effective sample sizes of at least 100 per chain).
#'
#' The time for a fit with missing states grows with the number of rows next
#' to a missing or range-coded state and with the square of the number of
#' states; fully observed transitions add little. As a guide, on a recent
#' laptop with states missing completely at random, an `"orm"` fit with
#' standard errors for 5,000 patients followed for 30 visits with 100 states
#' took about 15 seconds with 0.5% of the states missing and about three
#' minutes with 10%. Before fitting large data with many states and many
#' missing states, time a fit on a random subset of patients.
#'
#' ## Refits
#'
#' Bootstrap inference (`inferences(method = "bootstrap")` and
#' `method = "fwb"`) and [bootstrap_model_coefs()] refit the model to
#' resampled patients by calling `mostm()` again with the same engine and
#' arguments. The refits use the argument values from the original fit, so
#' changing a variable after fitting does not change the refits.
#'
#' The state occupancy, inference, and diagnostic functions accept only fits
#' made by `mostm()`; a model fitted directly with [rms::orm()],
#' [VGAM::vglm()], or [rmsb::blrm()] lacks the stored data they need.
#'
#' @return The fitted model, with the data needed by the state occupancy and
#'   inference functions attached. The returned object's `call` is the
#'   `mostm()` call. Its class depends on `engine`:
#'   * `"orm"`: an `orm` object. With `id_var`, its covariance matrix
#'     (`vcov(fit)`) holds the patient-clustered variance.
#'   * `"vglm"`: with `id_var`, a `robcov_vglm` object that holds the
#'     patient-clustered variance and the fitted `vglm` object in `$vglm_fit`
#'     (see [robcov_vglm()]); without `id_var`, the S4 `vglm` object.
#'   * `"blrm"`: a `blrm` object with the posterior draws.
#'
#'   A fit that uses missing or range-coded states also has class `mostm_ml`
#'   and a component `ml`, a list that includes:
#'   * `loglik`: the maximised log-likelihood of the observed data.
#'   * `missing_info_fraction`: the largest share of information lost because
#'     states are missing or known only up to a range, from 0 to 1.
#'   * `n_missing`, `n_range`: the numbers of rows with a missing state and
#'     with a range of states.
#'   * `patient_scores`, `information`, `var_model`: one row per patient of
#'     the derivatives of that patient's log-likelihood, the observed
#'     information, and its inverse, the model-based covariance. The robust
#'     covariance is `var_robust`.
#'   * `iterations`, `newton_steps`, `max_abs_score`: EM iterations, Newton
#'     steps, and the largest absolute derivative of the log-likelihood at the
#'     estimates (close to 0 at a maximum).
#'
#'   For `engine = "blrm"`, `ml` holds the counts of missing and range-coded
#'   rows, `diagnostics` the convergence diagnostics, and `rstan` the
#'   [rstan::stanfit-class] object.
#'
#' @examples
#' # Daily states of 200 simulated patients over 14 days: 1 is the best state
#' # and 8 is death
#' trajectories <- sim_markov(
#'   "actt2",
#'   n_patients = 200,
#'   follow_up_time = 14,
#'   seed = 1
#' )
#' model_data <- prepare_markov_data(trajectories, absorbing_state = 8)
#'
#' fit <- mostm(
#'   ordered(y) ~ yprev + rms::rcs(time, 3) + tx,
#'   data = model_data,
#'   id_var = "id"
#' )
#'
#' # The same model fitted with VGAM
#' fit_vglm <- mostm(
#'   ordered(y) ~ yprev + rms::rcs(time, 3) + tx,
#'   data = model_data,
#'   id_var = "id",
#'   engine = "vglm"
#' )
#' # Both engines reach the same maximum of the log-likelihood
#' c(orm = stats::logLik(fit), vglm = stats::logLik(fit_vglm$vglm_fit))
#'
#' # Probability of each state on day 14 with and without treatment. SOPs are
#' # calculated one day at a time, so request every day up to day 14.
#' sops_by_arm <- avg_sops(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:14,
#'   absorb = "8"
#' )
#' sops_by_arm[sops_by_arm$time == 14, ]
#'
#' # Remove the recorded state of 10% of the visits of living patients. The
#' # previous state after a missing visit is then unknown as well.
#' set.seed(2)
#' with_missing <- model_data
#' with_missing$y <- as.numeric(as.character(with_missing$y))
#' lost <- with_missing$y != 8 & stats::runif(nrow(with_missing)) < 0.1
#' with_missing$y[lost] <- NA
#' after <- c(FALSE, lost[-length(lost)]) &
#'   c(FALSE, with_missing$id[-1] == with_missing$id[-nrow(with_missing)])
#' with_missing$yprev[after] <- NA
#'
#' # A missing state is any state but death
#' fit_missing <- mostm(
#'   ordered(y) ~ yprev + rms::rcs(time, 3) + tx,
#'   data = with_missing,
#'   id_var = "id",
#'   missing_states = 1:7
#' )
#' fit_missing
#'
#' \donttest{
#' if (requireNamespace("rmsb", quietly = TRUE)) {
#'   # A Bayesian model with a separate time effect at each state threshold,
#'   # fitted to the first 60 patients to keep the example short
#'   fit_blrm <- mostm(
#'     ordered(y) ~ yprev + time + tx,
#'     data = model_data,
#'     id_var = "id",
#'     engine = "blrm",
#'     subset = id <= 60,
#'     ppo = ~time,
#'     iter = 1000,
#'     chains = 2,
#'     sampling.args = list(seed = 1)
#'   )
#' }
#' }
#'
#' @seealso [sops()] and [avg_sops()] for state occupancy probabilities,
#'   [inferences()] for confidence intervals, and [prepare_markov_data()] for
#'   preparing the data.
#' @export
mostm <- function(
  formula,
  data,
  ...,
  engine = c("orm", "vglm", "blrm"),
  id_var = NULL,
  time_var = "time",
  p_var = "yprev",
  first_followup_time = NULL,
  missing_states = NULL,
  weights = NULL,
  subset = NULL,
  type = c("HC0", "HC1"),
  cadjust = NULL
) {
  refit_args <- markov_capture_refit_args(
    match.call(),
    environment(),
    names(formals(sys.function())),
    ...
  )
  engine <- match.arg(engine)
  model_call <- match.call()
  model_call[[1L]] <- quote(mostm)
  model_call$engine <- engine
  dots <- match.call(expand.dots = FALSE)$...
  mostm_check_engine_args(engine, dots, refit_args)

  if (engine == "blrm") {
    if (!missing(type) || !is.null(cadjust)) {
      stop(
        "`type` and `cadjust` set the patient-clustered standard errors of the ",
        "\"orm\" and \"vglm\" engines. `engine = \"blrm\"` takes its ",
        "uncertainty from the posterior draws.",
        call. = FALSE
      )
    }
    if (!is.null(model_call$weights)) {
      stop(
        "`weights` is not available for `engine = \"blrm\"`, because ",
        "`rmsb::blrm()` does not accept row weights.",
        call. = FALSE
      )
    }
    if (!requireNamespace("rmsb", quietly = TRUE)) {
      stop(
        "Package 'rmsb' is required for `engine = \"blrm\"`.",
        call. = FALSE
      )
    }
  }
  type <- match.arg(type)
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be supplied as a data frame.", call. = FALSE)
  }
  if (terms_has_offset(stats::terms(formula, specials = "offset", data = data))) {
    stop_unsupported_offset()
  }

  id_var <- validate_markov_id_var(id_var, data, "data")
  if (is.null(id_var)) {
    warn_missing_markov_id_var()
  } else {
    warn_duplicate_markov_id_time(data, id_var, time_var = time_var)
  }

  # Row weights and `subset` are evaluated once here, in the caller's frame,
  # so they may refer to variables that exist only there. The weights are
  # stored as a column of `data`, which keeps them aligned with resampled
  # rows during bootstrap refits.
  eval_env <- parent.frame()
  prior_weight_column <- NULL
  if (!is.null(model_call$weights)) {
    prior_weights <- markov_eval_prior_weights(
      model_call$weights,
      data,
      eval_env
    )
    if (!is.null(prior_weights)) {
      prior_weight_column <- markov_prior_weight_column()
      data[[prior_weight_column]] <- prior_weights
    }
  }
  subset_rows <- if (!is.null(model_call$subset)) {
    eval(model_call$subset, data, eval_env)
  }

  setup <- list(
    formula = formula,
    data = data,
    dots = dots,
    subset_rows = subset_rows,
    prior_weight_column = prior_weight_column,
    eval_env = eval_env,
    data_name = if (is.language(model_call$data)) {
      as.character(model_call$data)
    } else {
      character()
    }
  )
  # With `missing_states` or an Ocens() response, rows with missing or
  # range-coded states are kept and the observed-data likelihood is
  # maximised (markov-ml-*.R). Otherwise the engine fits the complete rows.
  ml_problem <- mostm_ml_problem(
    setup,
    engine = engine,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    missing_states = missing_states,
    engine_args = list(...)
  )
  if (is.null(ml_problem)) {
    fitted <- switch(
      engine,
      orm = mostm_fit_orm(setup, ...),
      vglm = mostm_fit_vglm(setup, ...),
      blrm = mostm_fit_blrm(setup, ...)
    )
    # An explicit `na.action` already says how to handle missing values.
    if (is.null(dots$na.action)) {
      markov_warn_dropped_complete_rows(
        markov_subset_source_data(data, subset_rows, eval_env)$data,
        fitted$fit_data,
        formula = formula,
        p_var = p_var,
        missing_states_given = !is.null(missing_states)
      )
    }
  } else {
    dropped <- ml_problem$dropped
    markov_warn_dropped_rows(
      nrow(ml_problem$selected),
      response = dropped$response,
      previous = dropped$previous,
      covariate = dropped$covariate,
      p_var = p_var,
      missing_states_given = !is.null(missing_states)
    )
    fitted <- if (engine == "blrm") {
      mostm_fit_ml_blrm(setup, ml_problem, engine_args = list(...))
    } else {
      mostm_fit_ml(
        setup,
        engine = engine,
        problem = ml_problem,
        type = type,
        cadjust = cadjust %||% TRUE,
        control = markov_ml_control(variance = !markov_ml_coefficients_only()),
        engine_args = list(...)
      )
    }
    fitted$fit_data <- markov_ml_fit_data(ml_problem, ml_problem$selected)
  }

  fit_data <- fitted$fit_data
  stored <- markov_prepare_stored_data(
    data = data,
    formula = formula,
    subset = subset_rows,
    eval_env = eval_env,
    fit_data = fit_data,
    id_var = id_var,
    time_var = time_var,
    first_followup_time = first_followup_time
  )
  attach_data <- function(model) {
    model <- markov_attach_model_data(
      model,
      data = fit_data,
      id_var = id_var,
      refit_data = stored$refit_data,
      starting_profile_data = stored$starting_profile_data,
      starting_profile_metadata = stored$starting_profile_metadata,
      refit_args = refit_args,
      prior_weight_column = prior_weight_column
    )
    if (!is.null(ml_problem)) {
      model <- markov_ml_mark_fit(model)
    }
    markov_set_engine(model, engine)
  }

  fit <- attach_data(markov_set_model_call(fitted$fit, model_call))
  if (!is.null(ml_problem)) {
    if (engine == "blrm") {
      return(fit)
    }
    fit <- markov_ml_finish(
      fit,
      fitted$ml,
      ml_problem,
      engine = engine,
      id_var = id_var,
      fitted_ids = stored$starting_profile_metadata$fitted_ids
    )
    return(attach_data(fit))
  }
  # Bootstrap refits of fits with missing states keep `id_var` for their
  # chains; a resample with nothing missing needs no covariance either.
  if (is.null(id_var) || engine == "blrm" || markov_ml_coefficients_only()) {
    return(fit)
  }

  robust <- switch(
    engine,
    orm = robcov_orm(
      fit,
      cluster = fit_data[[id_var]],
      type = type,
      cadjust = cadjust
    ),
    vglm = robcov_vglm(
      fit,
      cluster = fit_data[[id_var]],
      type = type,
      cadjust = cadjust
    )
  )
  attach_data(robust)
}

mostm_fit_orm <- function(setup, ...) {
  formula <- setup$formula
  data <- setup$data
  subset_rows <- setup$subset_rows
  fit_call <- quote(rms::orm(
    formula = formula,
    data = data,
    ...,
    x = TRUE,
    y = TRUE
  ))
  if (!is.null(setup$prior_weight_column)) {
    fit_call$weights <- as.name(setup$prior_weight_column)
  }
  if (!is.null(subset_rows)) {
    fit_call$subset <- quote(subset_rows)
  }
  fit <- eval(fit_call)
  # When its fitting algorithm fails, rms::orm() returns an object that holds
  # only `fail = TRUE` instead of stopping.
  if (isTRUE(fit$fail)) {
    stop(
      "`rms::orm()` could not fit the model: its fitting algorithm stopped ",
      "without finding the coefficient estimates. You can try another ",
      "algorithm with `opt_method = \"LM\"`, scale the covariates with ",
      "`scale = TRUE`, allow more iterations with a larger `maxit`, or ",
      "simplify the model, for example by removing covariates that are ",
      "nearly collinear.",
      call. = FALSE
    )
  }

  fit_frame <- markov_rms_model_frame(
    formula,
    data,
    setup$dots,
    subset = subset_rows,
    weight_column = setup$prior_weight_column,
    eval_env = setup$eval_env
  )
  list(fit = fit, fit_data = markov_align_model_data(data, fit_frame))
}

mostm_fit_blrm <- function(setup, ...) {
  formula <- setup$formula
  data <- setup$data
  subset_rows <- setup$subset_rows
  fit_call <- quote(rmsb::blrm(
    formula = formula,
    data = data,
    ...,
    x = TRUE,
    y = TRUE
  ))
  if (!is.null(subset_rows)) {
    fit_call$subset <- quote(subset_rows)
  }
  fit <- eval(fit_call)

  fit_frame <- markov_rms_model_frame(
    formula,
    data,
    setup$dots,
    subset = subset_rows,
    eval_env = setup$eval_env
  )
  list(fit = fit, fit_data = markov_align_model_data(data, fit_frame))
}

# Arguments that mostm() handles itself or that every engine fit must keep.
# They cannot be passed to an engine through `...`.
mostm_fixed_args <- function() {
  c(
    names(formals(mostm)),
    "x",
    "y",
    "x.arg",
    "y.arg",
    "offset",
    "setup"
  )
}

# Names of the arguments each engine accepts through `...`. `NULL` means the
# engine's package is not installed.
mostm_engine_args <- function(engine) {
  args <- switch(
    engine,
    orm = c(names(formals(rms::orm)), names(formals(rms::orm.fit))),
    vglm = c(
      names(formals(mostm_fit_vglm)),
      names(formals(VGAM::vglm.control)),
      names(formals(vglm_categorical_control()))
    ),
    blrm = if (requireNamespace("rmsb", quietly = TRUE)) {
      names(formals(rmsb::blrm))
    }
  )
  if (is.null(args)) {
    return(NULL)
  }
  setdiff(unique(args), c(mostm_fixed_args(), "..."))
}

# VGAM adds the control settings of categorical families, such as
# VGAM::cumulative(), to those of VGAM::vglm.control().
vglm_categorical_control <- function() {
  tryCatch(
    utils::getFromNamespace("vglm.VGAMcategorical.control", "VGAM"),
    error = function(e) function() NULL
  )
}

mostm_check_engine_args <- function(engine, dots, values) {
  supplied <- names(dots)
  if (length(dots) && (is.null(supplied) || any(supplied == ""))) {
    stop(
      "Arguments passed to the engine through `...` must be named.",
      call. = FALSE
    )
  }
  markov_reject_legacy_args(dots)
  if ("offset" %in% supplied) {
    stop_unsupported_offset()
  }
  kept <- intersect(supplied, c("x", "y", "x.arg", "y.arg"))
  if (length(kept)) {
    stop(
      "`mostm()` always keeps the design matrix and response on the fit; ",
      "remove ",
      paste0("`", kept, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # rmsb::blrm() passes other arguments on to the Stan optimizer, so only
  # arguments that belong to another engine are rejected there. The other
  # engines' arguments are looked up only for names the engine does not
  # accept: looking up those of "blrm" loads rmsb and rstan, which takes
  # seconds in every orm or vglm fit and bootstrap refit.
  others <- setdiff(markov_engines(), engine)
  rejected <- setdiff(supplied, mostm_engine_args(engine))
  if (length(rejected)) {
    other_args <- lapply(others, mostm_engine_args)
    if (engine == "blrm") {
      rejected <- intersect(rejected, unlist(other_args))
    }
  }
  if (length(rejected)) {
    described <- vapply(
      rejected,
      function(arg) {
        owners <- others[vapply(other_args, function(a) arg %in% a, TRUE)]
        if (!length(owners)) {
          return(paste0("`", arg, "`"))
        }
        paste0(
          "`",
          arg,
          "` (an argument of ",
          paste0("`engine = \"", owners, "\"`", collapse = " and "),
          ")"
        )
      },
      character(1)
    )
    stop(
      "Not an argument of `engine = \"",
      engine,
      "\"`: ",
      paste(described, collapse = ", "),
      ". See `?mostm` for the arguments each engine accepts.",
      call. = FALSE
    )
  }

  if (
    engine == "orm" &&
      !is.null(values$family) &&
      !is.character(values$family)
  ) {
    stop(
      "For `engine = \"orm\"`, `family` names the link function of ",
      "`rms::orm()`, such as \"logistic\". VGAM family objects apply to ",
      "`engine = \"vglm\"`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

markov_set_model_call <- function(model, call) {
  if (isS4(model)) {
    methods::slot(model, "call") <- call
  } else {
    model$call <- call
  }
  model
}
