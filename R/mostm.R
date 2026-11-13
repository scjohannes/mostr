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
#' In a partial proportional odds model the thresholds (cut-points) between
#' states move with the predictors that have a separate effect at each
#' threshold, so for some covariate values, often outside the range of the
#' data, two thresholds can cross. The model then gives a negative probability
#' for a state, and state occupancy probabilities aren't defined. [sops()],
#' [avg_sops()], [state_time()], [comparisons()] and the other functions that
#' use the fit stop with an error that says for how many profiles and at which
#' time this happens; for `"blrm"` fits it also says in how many posterior
#' draws, and [inferences()] says in how many simulation draws or bootstrap
#' samples. They never set the negative probability to zero. Check the
#' covariate values of those profiles, or simplify the terms with a separate
#' effect at each threshold. A negative value between \eqn{-10^{-10}} and zero
#' comes from rounding, not from crossed thresholds, and is set to zero.
#'
#' @param formula A model formula with the current state as the response, for
#'   example `ordered(y) ~ yprev + rms::rcs(time, 4) + tx`. Offsets are not
#'   supported. For the `"orm"` and `"blrm"` engines, categorical predictors,
#'   including time, must be factor or character columns of `data` rather
#'   than converted with `factor()` inside the formula; for factor time, also
#'   set `first_followup_time`. Random-effect terms are available only for the
#'   `"blrm"` engine; see "Random effects" below.
#'
#'   Time terms are written inline, for example `rms::rcs(time, 4)`,
#'   `log(time)` or `rms::pol(time, 2)`; the prediction functions evaluate
#'   them at each requested time. For the `"orm"` and `"blrm"` engines, use
#'   `rms::rcs()` and `rms::pol()` rather than `splines::ns()`,
#'   `splines::bs()` or `stats::poly()`, for time and for any other variable:
#'   rms evaluates the latter again on the prediction rows, which gives a
#'   different basis, so `mostm()` refuses them.
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
#' @param id_var Name of the column in `data` that identifies patients, for
#'   example `"id"`. It is required, because the rows of one patient belong
#'   together: `mostm()` uses it to store each patient's starting row for
#'   automatic predictions, the `"orm"` and `"vglm"` engines use it for
#'   standard errors that allow for correlation between the transitions of one
#'   patient, and the bootstrap resamples or reweights whole patients by it.
#' @param time_var Name of the time column in `data`. It is used to find each
#'   patient's starting row, and [sops()] and the other prediction functions
#'   set this column to each requested visit.
#' @param p_var Name of the column in `data` that holds each patient's state at
#'   the previous visit, `"yprev"` by default. The prediction functions replace
#'   it with each possible previous state as they move from visit to visit.
#' @param p2_var For a second-order model, the name of the column that holds
#'   each patient's state two visits back, for example `"ypprev"`. `NULL` (the
#'   default) fits a first-order model, in which the next state depends only on
#'   the previous one. The prediction functions take it from the fit, as they
#'   do `id_var`, `time_var` and `p_var`, and then use second-order
#'   calculations. Without it, a column with the older state in the formula
#'   would be held at each patient's starting value.
#' @param absorb The absorbing states, such as death: states that a patient
#'   never leaves, given as state labels, for example `absorb = "8"` or
#'   `absorb = c("home", "dead")`. The prediction functions keep a patient who
#'   reaches one of them there, and they never predict a transition out of
#'   it. `NULL` (the default) means that no state is absorbing. See
#'   "Absorbing states" below.
#' @param gap_var Name of the column that holds the time since each
#'   patient's previous visit, if the model includes it, for example `"gap"`.
#'   The prediction functions set this column to the difference between
#'   consecutive requested times, so the visits in `times` may be spaced
#'   unevenly. Time must be numeric. See "Time since the previous visit"
#'   below.
#' @param first_followup_time Time of the first visit after baseline. Each
#'   patient's row at this time is their starting row for automatic
#'   predictions. For numeric time, `NULL` (the default) uses 1, time may not
#'   be below 1, and a value you supply must be the earliest time in the data.
#'   Factor and character time need an explicit value, such as `"day 1"`.
#' @param weights Optional non-negative row weights: a column of `data` or a
#'   vector with one value per row of `data`. Not available for
#'   `engine = "blrm"`. The weights are used only to fit the model:
#'   [avg_sops()], [avg_state_time()], and [avg_comparisons()] give every
#'   profile the same weight, and [inferences()] offers only conditional
#'   variance for a weighted fit.
#' @param subset Optional selection of the rows of `data` to fit, as a logical
#'   expression evaluated in `data`, for example `time <= 28`.
#' @param type Small-sample correction of the patient-clustered standard
#'   errors of the `"orm"` and `"vglm"` engines. `"HC0"` (the default) applies
#'   none. `"HC1"` multiplies the variance by \eqn{(m - 1) / (m - p)}, where
#'   \eqn{m} is the number of fitted rows with positive weight and \eqn{p} is
#'   the number of model coefficients.
#' @param cadjust Logical. Multiply the patient-clustered variance of the
#'   `"orm"` and `"vglm"` engines by \eqn{n / (n - 1)}, where \eqn{n} is the
#'   number of patients? `NULL` (the default) applies this correction. It is
#'   independent of `type`, so `type = "HC1"` with `cadjust = TRUE` applies
#'   both factors.
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
#' ## Absorbing states
#'
#' A patient who reaches an absorbing state, such as death, stays there. The
#' data therefore contain transitions into it but none out of it, and the
#' model can't say what happens after it. Name such states in `absorb`: the
#' prediction functions then keep the probability of being in them, and
#' predict transitions only from the other states. `mostm()` never works out
#' the absorbing states from the data, because a state that is seen only as
#' an outcome, at a patient's last visit, need not be absorbing.
#'
#' `mostm()` checks `absorb` against the fitted rows of `data`:
#'
#' * Each label must be a state of the response.
#' * No fitted row may have an absorbing state as its previous state
#'   (`p_var`) or, for a second-order model, as the state two visits back
#'   (`p2_var`): a patient in an absorbing state has no further transitions.
#'   Remove such rows first; [prepare_markov_data()] does this for the
#'   previous state.
#' * A state that occurs as an outcome but never as a previous state, and is
#'   not named in `absorb`, gives a warning: it may be absorbing. Without
#'   `absorb`, a model with a numeric previous state, such as
#'   `rms::rcs(yprev, 4)`, predicts transitions out of such a state from the
#'   other states' transitions, so dead patients can come back to life. With
#'   a factor previous state, the prediction functions stop instead, because
#'   the model has no coefficient for that previous state.
#'
#' Bootstrap refits keep the absorbing states of the fit. A resample in which
#' no patient reaches one of them has no transitions into it, so its refit
#' leaves that state out, and the state's probability is 0 in that resample.
#'
#' ## Time since the previous visit
#'
#' When visits are spaced unevenly, the model can include the time since each
#' patient's previous visit as a predictor, in a column named by `gap_var`.
#' The prediction functions then set that column for every transition: to the
#' difference between the requested time and the one before it, and, for the
#' first requested time, to that time itself, that is, the time since time 0.
#' So with numeric time and `times = c(1, 2, 4, 7)`, the gaps are 1, 1, 2 and
#' 3. A patient's first visit after baseline in the fitted data should
#' therefore have a gap equal to its time.
#'
#' `mostm()` checks that the gap column is numeric, that time is numeric, and
#' that on every row after a patient's first the gap equals the time since
#' the patient's previous row, so that both columns use the same unit. Factor
#' or character time can't be used with `gap_var`, because visit labels don't
#' say how much time passes between visits.
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
#' For the `"orm"` and `"vglm"` engines, the covariance matrix of the
#' coefficients stored on the fit gives standard errors that remain valid when
#' the transitions of one patient are correlated. It comes from a sandwich
#' estimator that adds up the contributions of each patient, as identified by
#' `id_var` (cluster-robust standard errors); `type` and `cadjust` set its
#' small-sample corrections. These standard errors are used by
#' `inferences(method = "delta")` with conditional variance and by
#' `inferences(method = "mvn")`. Unconditional variance also accounts for which
#' patients were sampled; see [inferences()]. For the `"blrm"` engine,
#' uncertainty comes from the posterior draws.
#'
#' `mostm()` stops when this covariance cannot be computed, for example when
#' the fitting algorithm has not converged or some coefficients grow without
#' bound because a transition never occurs in the data. The `"orm"` engine
#' also needs the logit link (`family = "logistic"`, the default), which the
#' state occupancy functions require as well.
#'
#' ## Random effects
#'
#' A random effect is a shift on the log-odds scale that belongs to one patient
#' and applies at all of that patient's visits. It describes patients who do
#' better or worse throughout than their previous state and covariates predict.
#'
#' Only `engine = "blrm"` supports random effects. For a random intercept for
#' each patient, add `cluster(id)` to the formula, where `id` is the column that
#' identifies patients, and attach rms and rmsb with `library()` so that
#' `cluster()` can be found. Write the term without a package prefix:
#' [rmsb::blrm()] would fit `rms::cluster(id)` as an ordinary predictor, a
#' straight-line effect of the ID number, so `mostm()` stops with an error. To
#' include each patient's estimated random effect in predictions, use
#' `include_re = TRUE` in [sops()] or [avg_sops()]; see
#' [the Bayesian article](https://scjohannes.github.io/mostr/articles/bayesian-rmsb-sops.html).
#'
#' With `engine = "orm"` or `"vglm"`, `cluster()` and `mix_re()` terms are
#' errors, although [rms::orm()] can fit them from rms version 8.2 on. State
#' occupancy probabilities from such a fit would describe only a patient whose
#' random effect is 0, which is not the average over patients. With these
#' engines, the patient-clustered standard errors described under "Standard
#' errors" already allow for correlation between the transitions of one
#' patient. `mix_re()`, which lets the
#' size of the random effect differ between the first follow-up and later visits,
#' is not available for any engine.
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
#'   * `"orm"`: an `orm` object whose covariance matrix (`vcov(fit)`) holds
#'     the patient-clustered variance.
#'   * `"vglm"`: a `robcov_vglm` object that holds the patient-clustered
#'     variance and the fitted `vglm` object in `$vglm_fit` (see
#'     [robcov_vglm()]).
#'   * `"blrm"`: a `blrm` object with the posterior draws.
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
#' # State 8, death, is absorbing: patients who reach it stay there.
#' fit <- mostm(
#'   ordered(y) ~ yprev + rms::rcs(time, 3) + tx,
#'   data = model_data,
#'   id_var = "id",
#'   absorb = "8"
#' )
#'
#' # The same model fitted with VGAM
#' fit_vglm <- mostm(
#'   ordered(y) ~ yprev + rms::rcs(time, 3) + tx,
#'   data = model_data,
#'   id_var = "id",
#'   absorb = "8",
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
#'   times = 1:14
#' )
#' sops_by_arm[sops_by_arm$time == 14, ]
#'
#' \donttest{
#' if (requireNamespace("rmsb", quietly = TRUE)) {
#'   # A Bayesian model with a separate time effect at each state threshold,
#'   # fitted to the first 60 patients to keep the example short
#'   fit_blrm <- mostm(
#'     ordered(y) ~ yprev + time + tx,
#'     data = model_data,
#'     id_var = "id",
#'     absorb = "8",
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
  id_var,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  absorb = NULL,
  gap_var = NULL,
  first_followup_time = NULL,
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
      mostr_abort(
        "`type` and `cadjust` set the patient-clustered standard errors of the ",
        "\"orm\" and \"vglm\" engines. `engine = \"blrm\"` takes its ",
        "uncertainty from the posterior draws."
      )
    }
    if (!is.null(model_call$weights)) {
      mostr_abort(
        "`weights` is not available for `engine = \"blrm\"`, because ",
        "`rmsb::blrm()` does not accept row weights."
      )
    }
    if (!requireNamespace("rmsb", quietly = TRUE)) {
      mostr_abort(
        "Package 'rmsb' is required for `engine = \"blrm\"`."
      )
    }
  }
  type <- match.arg(type)
  if (missing(data) || !is.data.frame(data)) {
    mostr_abort("`data` must be supplied as a data frame.")
  }
  if (terms_has_offset(stats::terms(formula, specials = "offset", data = data))) {
    stop_unsupported_offset()
  }
  if (engine %in% c("orm", "blrm")) {
    mostm_check_inline_factors(formula, engine, time_var)
    mostm_check_data_dependent_bases(formula, engine)
  }
  mostm_check_random_effects(formula, engine)

  if (missing(id_var) || is.null(id_var)) {
    stop_missing_mostm_id_var()
  }
  id_var <- validate_markov_id_var(id_var, data, "data")
  if (
    !is.character(time_var) ||
      length(time_var) != 1L ||
      is.na(time_var) ||
      !time_var %in% names(data)
  ) {
    mostr_abort(
      "The time column `",
      format(time_var),
      "` is not in `data`. Set `time_var` to the column that holds the ",
      "visit time; `mostm()` uses it to find each patient's starting row."
    )
  }
  validate_mostm_previous_state_vars(p_var, p2_var, data)
  if (!missing(p_var) && !p_var %in% names(data)) {
    mostr_abort(
      "The column `p_var = \"",
      p_var,
      "\"` is not in `data`. Set `p_var` to the column that holds each ",
      "patient's previous state."
    )
  }
  # Bootstrap refits use only the coefficients, so they skip the duplicate
  # check, the starting rows, and the patient-clustered covariance.
  coefficients_only <- markov_refit_coefficients_only()
  if (!coefficients_only) {
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

  # Bootstrap refits reuse the checked absorbing states and gap column; a
  # resample without a death would otherwise fail the check that `absorb`
  # names states of the response.
  if (!coefficients_only) {
    absorb <- validate_mostm_absorb(
      absorb,
      formula = formula,
      data = data,
      rows = mostm_subset_row_index(subset_rows, nrow(data)),
      p_var = p_var,
      p2_var = p2_var
    )
    validate_mostm_gap_var(gap_var, data, id_var, time_var)
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
  fitted <- switch(
    engine,
    orm = mostm_fit_orm(setup, ...),
    vglm = mostm_fit_vglm(setup, ...),
    blrm = mostm_fit_blrm(setup, ...)
  )

  fit_data <- fitted$fit_data
  stored <- markov_prepare_stored_data(
    data = data,
    formula = formula,
    subset = subset_rows,
    eval_env = eval_env,
    fit_data = fit_data,
    id_var = id_var,
    time_var = time_var,
    first_followup_time = first_followup_time,
    starting_profiles = !coefficients_only
  )
  # The column names and model settings that every later call on the fit
  # reads, built once here. The engines record the states in fitted order
  # (VGAM as the response's column names, rms and rmsb as its levels); an
  # object without them can't be used to predict (`stop_no_fitted_states()`).
  y_levels <- markov_model_ylevels_or_null(fitted$fit)
  settings <- markov_settings(
    engine = engine,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    first_followup_time = stored$starting_profile_metadata$first_followup_time,
    from_fit = setdiff(
      markov_settings_fit_fields(),
      if (is.null(y_levels)) "y_levels"
    )
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
    model <- markov_set_settings(model, settings)
    markov_set_engine(model, engine)
  }

  fit <- attach_data(markov_set_model_call(fitted$fit, model_call))
  if (coefficients_only || engine == "blrm") {
    return(fit)
  }
  # The patient-clustered covariance of orm fits needs the logit link, which
  # every state occupancy function also requires.
  if (engine == "orm" && !identical(fit$family, "logistic")) {
    stop_unsupported_markov_link(fit$family)
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
  # With weights, rms::orm() warns that its own bootstrap ignores them; see
  # suppress_orm_weights_warning() for why mostr hides that one warning.
  fit <- suppress_orm_weights_warning(eval(fit_call))
  # When its fitting algorithm fails, rms::orm() returns an object that holds
  # only `fail = TRUE` instead of stopping.
  if (isTRUE(fit$fail)) {
    stop_failed_orm_fit()
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
    mostr_abort(
      "Arguments passed to the engine through `...` must be named."
    )
  }
  markov_reject_legacy_args(dots)
  if ("offset" %in% supplied) {
    stop_unsupported_offset()
  }
  kept <- intersect(supplied, c("x", "y", "x.arg", "y.arg"))
  if (length(kept)) {
    mostr_abort(
      "`mostm()` always keeps the design matrix and response on the fit; ",
      "remove ",
      paste0("`", kept, "`", collapse = ", "),
      "."
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
    mostr_abort(
      "Not an argument of `engine = \"",
      engine,
      "\"`: ",
      paste(described, collapse = ", "),
      ". See `?mostm` for the arguments each engine accepts."
    )
  }

  if (
    engine == "orm" &&
      !is.null(values$family) &&
      !is.character(values$family)
  ) {
    mostr_abort(
      "For `engine = \"orm\"`, `family` names the link function of ",
      "`rms::orm()`, such as \"logistic\". VGAM family objects apply to ",
      "`engine = \"vglm\"`."
    )
  }
  invisible(NULL)
}

# rms::orm() and rmsb::blrm() name the design matrix columns of a factor after
# the variable (`x=2`), but stats::model.matrix() names those of a variable
# converted inside the formula after the call (`factor(x)2`). rms then fails
# with "subscript out of bounds", so mostm() rejects such terms before fitting.
mostm_check_inline_factors <- function(formula, engine, time_var) {
  calls <- formula_inline_factor_calls(formula)
  if (!length(calls)) {
    return(invisible(NULL))
  }
  labels <- unique(vapply(calls, deparse1, character(1)))
  several <- length(labels) > 1L
  fitter <- switch(engine, orm = "`rms::orm()`", blrm = "`rmsb::blrm()`")

  # A call whose first argument is a variable name gives a concrete example.
  bare <- Filter(function(call) length(call) >= 2L && is.name(call[[2L]]), calls)
  example <- NULL
  if (length(bare)) {
    converted <- bare[[1L]]
    column <- call("$", as.name("data"), converted[[2L]])
    converted[[2L]] <- column
    example <- deparse1(call("<-", column, converted))
  }

  converts_time <- is.character(time_var) &&
    any(vapply(calls, function(call) any(time_var %in% all.vars(call)), TRUE))

  mostr_abort(
    fitter,
    " cannot fit ",
    if (several) {
      "variables that are converted to factors"
    } else {
      "a variable that is converted to a factor"
    },
    " inside the formula: ",
    format_sop_visit_list(paste0("`", labels, "`")),
    ". Create ",
    if (several) "each factor" else "the factor",
    " as a column of `data` and use that column in the formula instead",
    if (!is.null(example)) paste0(", for example `", example, "`"),
    ".",
    if (converts_time) {
      paste0(
        " To use time as a factor, store it as a factor column of `data` ",
        "and set `first_followup_time` to the level of the first visit ",
        "after baseline, for example `first_followup_time = \"1\"`. See the ",
        "article on time as a factor: ",
        "https://scjohannes.github.io/mostr/articles/factor-time-orm-sops.html"
      )
    }
  )
}

# The rows of `data` that `subset` selects, as row numbers; all rows without
# `subset`. A missing value in a logical `subset` drops the row, as in the
# engines.
mostm_subset_row_index <- function(subset_rows, n) {
  if (is.null(subset_rows)) {
    return(seq_len(n))
  }
  if (is.logical(subset_rows)) {
    return(which(rep_len(subset_rows, n)))
  }
  seq_len(n)[subset_rows]
}

# `splines::ns()`, `splines::bs()` and `stats::poly()` choose their knots or
# orthogonal polynomials from the data they are evaluated on. `rms::orm()` and
# `rmsb::blrm()` keep only the formula, so mostr's orm design matrices and
# `rms::predictrms()` evaluate these terms again on the prediction rows and
# get another basis, without any error. Measured on 2026-11-13
# (benchmarks/local/settings-refactor/spline-covariates.R): with `ns(age, 3)`
# for a covariate, the design matrix of 40 profiles differed from the fitted
# one by up to 0.42 and their average SOPs by up to 0.042; `bs()` and `poly()`
# were as wrong, and one profile could even fail. `rms::rcs()` and `rms::pol()`
# store their knots or powers and agreed exactly. VGAM keeps the fitted basis,
# so `engine = "vglm"` is not checked. `poly(raw = TRUE)` has no data-dependent
# part and is allowed.
mostm_check_data_dependent_bases <- function(formula, engine) {
  calls <- formula_rhs_calls(formula, c("ns", "bs", "poly"))
  calls <- Filter(
    function(call) {
      fun <- call[[1L]]
      name <- as.character(if (is.call(fun)) fun[[3L]] else fun)
      if (!identical(name, "poly")) {
        return(TRUE)
      }
      matched <- tryCatch(
        match.call(stats::poly, call),
        error = function(e) call
      )
      !isTRUE(matched$raw)
    },
    calls
  )
  if (!length(calls)) {
    return(invisible(NULL))
  }
  labels <- unique(vapply(calls, deparse1, character(1)))
  fitter <- switch(engine, orm = "`rms::orm()`", blrm = "`rmsb::blrm()`")
  first <- calls[[1L]]
  variable <- if (length(first) >= 2L && is.name(first[[2L]])) {
    as.character(first[[2L]])
  } else {
    "x"
  }
  mostr_abort(
    "With ",
    fitter,
    ", ",
    format_sop_visit_list(paste0("`", labels, "`")),
    if (length(labels) == 1L) " would be" else " would each be",
    " evaluated again on the prediction rows, which gives a different basis ",
    "without any error, so the state occupancy probabilities would be ",
    "wrong. Use `rms::rcs()` for a spline, for example `rcs(",
    variable,
    ", 4)`, or `rms::pol()` for a polynomial, for example `pol(",
    variable,
    ", 2)`: both keep the knots or terms of the fit. `engine = \"vglm\"` ",
    "keeps the fitted basis of these terms."
  )
}

# rms 8.2 and later fit `cluster()` and `mix_re()` terms as random effects, and
# state occupancy probabilities from such a fit would be wrong without any
# message. VGAM::vglm() and rmsb::blrm() would fit these terms as ordinary
# predictors, except that rmsb::blrm() recognizes a bare `cluster(id)` as a
# random intercept. mostm() therefore rejects them before fitting.
mostm_check_random_effects <- function(formula, engine) {
  calls <- formula_rhs_calls(formula, c("cluster", "mix_re"))
  if (engine == "blrm") {
    calls <- Filter(
      function(call) !identical(call[[1L]], as.name("cluster")),
      calls
    )
  }
  if (!length(calls)) {
    return(invisible(NULL))
  }
  labels <- unique(vapply(calls, deparse1, character(1)))
  several <- length(labels) > 1L
  formatted <- format_sop_visit_list(paste0("`", labels, "`"))

  if (engine == "blrm") {
    mostr_abort(
      "`rmsb::blrm()` would fit ",
      formatted,
      if (several) {
        " as ordinary predictors rather than random effects. "
      } else {
        " as an ordinary predictor rather than a random effect. "
      },
      "For a random intercept for each patient, write `cluster(id)` without ",
      "a package prefix and attach rms and rmsb with `library()`. Other ",
      "random-effect terms, such as `mix_re()`, are not available."
    )
  }

  mostr_abort(
    "`engine = \"",
    engine,
    "\"` does not support random effects, so remove ",
    formatted,
    " from the formula. ",
    if (engine == "orm") {
      paste0(
        "State occupancy probabilities from such a fit would describe only a ",
        "patient whose random effect is 0, which is not the average over ",
        "patients. "
      )
    } else if (several) {
      "`VGAM::vglm()` would fit such terms as ordinary predictors. "
    } else {
      "`VGAM::vglm()` would fit such a term as an ordinary predictor. "
    },
    "`id_var` already gives standard errors that allow for correlation ",
    "between the transitions of one patient. For a random intercept for each ",
    "patient, use `engine = \"blrm\"` with `cluster(id)` in the formula."
  )
}

# Calls on the right-hand side of `formula` that convert a variable to a
# factor, such as `factor(time)` or `base::as.ordered(tx)`, wherever they
# appear, including inside interactions and other calls. The response, often
# `ordered(y)`, is not searched.
formula_inline_factor_calls <- function(formula) {
  formula_rhs_calls(formula, c("factor", "as.factor", "ordered", "as.ordered"))
}

# Calls to `functions`, with or without a package prefix, on the right-hand side
# of `formula`, wherever they appear. A matching call is not searched further.
formula_rhs_calls <- function(formula, functions) {
  find_calls <- function(expr) {
    fun <- expr[[1L]]
    if (
      is.call(fun) &&
        length(fun) == 3L &&
        is.name(fun[[1L]]) &&
        as.character(fun[[1L]]) %in% c("::", ":::")
    ) {
      fun <- fun[[3L]]
    }
    if (is.name(fun) && as.character(fun) %in% functions) {
      return(list(expr))
    }
    args <- Filter(is.call, as.list(expr)[-1L])
    unlist(lapply(args, find_calls), recursive = FALSE)
  }
  rhs <- formula[[length(formula)]]
  if (!is.call(rhs)) {
    return(list())
  }
  find_calls(rhs)
}

markov_set_model_call <- function(model, call) {
  if (isS4(model)) {
    methods::slot(model, "call") <- call
  } else {
    model$call <- call
  }
  model
}
