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
#'   patient's starting row.
#' @param first_followup_time Time of the first visit after baseline. Each
#'   patient's row at this time is their starting row for automatic
#'   predictions. For numeric time, `NULL` (the default) uses 1, time may not
#'   be below 1, and a value you supply must be the earliest time in the data.
#'   Factor and character time need an explicit value, such as `"day 1"`.
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
    markov_set_engine(model, engine)
  }

  fit <- attach_data(markov_set_model_call(fitted$fit, model_call))
  if (is.null(id_var) || engine == "blrm") {
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
