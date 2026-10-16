# Observed-data maximum likelihood for first-order rms Markov models.

#' Fit an rms Markov Model by Maximum Likelihood with Missing States
#'
#' `orm_markov_ml()` fits the same first-order proportional-odds Markov
#' transition model as [orm_markov()], but uses every visit of every patient
#' when some outcome states are missing or known only up to a range of states.
#' [orm_markov()] drops every transition that touches such a visit;
#' `orm_markov_ml()` instead sums the likelihood over all states the patient
#' could have been in. When no state is missing, both functions give the same
#' fit. The result works with [sops()], [avg_sops()], [avg_time()],
#' [avg_comparisons()], and [inferences()] like an [orm_markov()] fit.
#'
#' @section The model:
#' Patient \eqn{i} has states \eqn{Y_{i0}} (the previous state on the first
#' row, `p_var`), \eqn{Y_{i1}, \dots, Y_{iT_i}} on an ordinal scale
#' \eqn{1, \dots, K}. For a visit \eqn{t \ge 1} and a non-absorbing previous
#' state \eqn{a}, the model is
#' \deqn{P(Y_{it} \ge j \mid Y_{i,t-1} = a) =
#'   \mathrm{logit}^{-1}(\alpha_j + x_{it}(a)^\top \beta), \quad j = 2, \dots, K,}
#' where \eqn{x_{it}(a)} is the design row built by `formula` from the
#' covariates on that row with `p_var` set to \eqn{a}. Absorbing states are
#' never left. Together these give a \eqn{K \times K} transition matrix
#' \eqn{P^{(t)}} for each row.
#'
#' Each row also records which states the observation allows: one state when
#' it is observed, a range of states for an [rms::Ocens()] response, or all
#' states allowed by `missing_states` when it is missing. Writing \eqn{D_{it}}
#' for the diagonal 0/1 matrix of allowed states, patient \eqn{i} contributes
#' \deqn{L_i = e_{Y_{i0}}^\top \Big[\prod_{t=1}^{T_i} P^{(t)} D_{it}\Big]
#'   \mathbf{1},}
#' the probability of everything observed for that patient given the starting
#' state and covariates. When two consecutive states are observed, the factor
#' for that transition is exactly the row that [orm_markov()] uses.
#' `orm_markov_ml()` maximises \eqn{\sum_i \log L_i}.
#'
#' @section Missing states:
#' Use `NA` in the response for a visit whose state is unknown, keeping the row
#' (for example from [prepare_markov_data()], which keeps rows with a missing
#' response). `missing_states` says what a missing response means:
#'
#' * `"any"`: nothing is known about the visit, so the patient may also have
#'   died (or entered another absorbing state) at or before it. Choose this
#'   when vital status was not recorded at visits with a missing state.
#' * `"non_absorbing"`: the patient was known to be alive, in some
#'   non-absorbing state. Choose this when vital status is known at every
#'   visit, so a missing state never hides a death. A missing visit then still
#'   carries the information that the patient survived to it.
#'
#' The chain follows the rows of each patient in `time_var` order. A visit
#' without a row is treated as not scheduled, not as a missing state, so keep
#' one row for every scheduled visit until an observed absorbing state or the
#' end of follow-up. In particular, do not end a patient's rows at a missing
#' state: if a death may have gone unrecorded, the later visits (showing the
#' death, or missing) must stay in the data, or the death is lost.
#'
#' Rows after a missing state usually have a missing or wrong `p_var`. Only
#' the first row's `p_var` (the starting state, which must be observed) is
#' used; afterwards the previous state is taken from the response. Where the
#' previous row's state is observed, `p_var` must agree with it.
#'
#' @section States known up to a range:
#' The response can be [rms::Ocens()]`(lower, upper)` instead of a single
#' state column, as in [rms::orm()]: `Ocens(y_lo, y_hi) ~ rcs(time, 4) + tx +
#' yprev`. With states 1 to 8 and death as state 8, `Ocens(1, 7)` means
#' "alive, level unknown" and `Ocens(b, b)` is an observed state `b`. Two
#' finite endpoints include both endpoints. A lower endpoint of `-Inf` means
#' "below the upper endpoint", and an upper endpoint of `Inf` means "above the
#' lower endpoint". Both endpoints `NA` is a missing state and follows
#' `missing_states`. Rows with observed, range, and missing states can be
#' mixed freely. Whether a state is reported exactly may not depend on the
#' unreported exact state.
#'
#' @section Assumptions:
#' The fit is valid when the states are *missing at random*: whether a state
#' is missing may depend on what was observed (for example earlier states or
#' treatment), but not on the unobserved state itself. If, say, patients in
#' worse health are more likely to miss a visit than observed patients in the
#' same previous state, the estimates are biased. Sensitivity analyses for
#' such missingness are outside the scope of this function.
#'
#' The model must be first order: apart from `p_var`, every covariate must be
#' fixed per patient (such as treatment) or a function of time only. A
#' covariate computed from the state history, such as a lagged state or days
#' spent in hospital so far, would be wrong after a missing state and is not
#' allowed. The starting state and all covariates must be observed.
#'
#' @section Algorithm:
#' The fit starts from the [rms::orm()] fit on transitions with both states
#' observed. It then alternates between computing, for every row, the
#' probability of each possible transition given everything observed for the
#' patient (a forward-backward pass), and refitting the proportional-odds
#' model with these probabilities as weights (EM). A few Newton steps then
#' finish the maximisation. The returned object is an [rms::orm()] fit
#' evaluated at the maximum likelihood estimates, on one row per distinct
#' combination of covariates, previous state, and state, weighted by the
#' expected number of such transitions. Its `x` and `y` components therefore
#' hold these combinations, not the original rows. For the same reason, the
#' fit statistics that [rms::orm()] stores (`stats`, `freq`, and `deviance`)
#' describe these expected counts, not the patients' data. [logLik()],
#' [AIC()], and [BIC()] instead use the maximised log-likelihood
#' \eqn{\sum_i \log L_i} (also in `$ml$loglik`), with one degree of freedom
#' per coefficient. Their number of observations, returned by [nobs()], is
#' the number of rows whose observation rules out at least one state; when no
#' state is missing, this is the number of transitions that [orm_markov()]
#' uses, so both functions give the same values. Printing the fit shows the
#' log-likelihood, the numbers of patients and rows, and the coefficients
#' with robust standard errors, Wald Z statistics, and p-values.
#'
#' @section Standard errors:
#' Standard errors use the observed information \eqn{J}, the negative second
#' derivative of \eqn{\sum_i \log L_i}, computed from exact patient scores
#' \eqn{U_i}. The model-based covariance is \eqn{J^{-1}}. The covariance
#' returned by [vcov()] and used by [inferences()] is the patient-cluster
#' robust covariance \eqn{J^{-1} (\sum_i U_i U_i^\top) J^{-1}}, multiplied by
#' the same corrections as in [orm_markov()]: `type = "HC1"` uses
#' \eqn{(m-1)/(m-p)} with \eqn{m} the number of rows whose observation
#' excludes at least one state and \eqn{p} the number of coefficients, and
#' `cadjust` uses \eqn{n/(n-1)} with \eqn{n} the number of patients with at
#' least one such row. Unconditional variance in [inferences()] and the score
#' bootstrap use the stored patient scores, so every patient with an observed
#' starting state contributes a starting profile, including patients whose
#' later states are all missing.
#'
#' Fractional weighted bootstrap (`method = "fwb"`) is not available for these
#' fits. Ordinary bootstrap refits resample whole patients, including their
#' missing rows.
#'
#' @param formula Model formula. The response is a state column (with `NA`
#'   for missing states) or `Ocens(lower, upper)`. The right-hand side must
#'   contain `p_var`. Offsets are not supported.
#' @param data Long data with one row per scheduled visit and columns for the
#'   patient ID, time, response, previous state, and covariates, keeping rows
#'   whose response is missing.
#' @param ... Further arguments passed to [rms::orm()], such as `scale` or
#'   `opt_method`. An `opt_method` is used by the starting fit, every EM
#'   step, and the final fit. Without it, Newton-Raphson is used and, if it
#'   fails (for example when some transitions are nearly deterministic), the
#'   step is repeated with Levenberg-Marquardt (`opt_method = "LM"`).
#'   `weights`, `subset`, `offset`, `penalty`, `cluster`, `family`, and
#'   `p2_var` are not supported.
#' @param id_var Character scalar naming the patient ID column. Required: each
#'   patient's rows form one chain.
#' @param absorbing_state State value(s) that are never left, such as death.
#'   Required.
#' @param p_var Name of the previous-state column used in `formula`, as in
#'   [sops()].
#' @param missing_states What a missing response means: `"any"` (any state,
#'   including absorbing ones) or `"non_absorbing"` (a non-absorbing state).
#'   See the Missing states section.
#' @param type,cadjust Corrections for the robust covariance, as in
#'   [orm_markov()].
#' @param time_var Character scalar naming the time column. Rows are ordered by
#'   this column within each patient.
#' @param first_followup_time First scheduled post-baseline time, used to pick
#'   each patient's starting profile as in [orm_markov()].
#' @param control Settings from `orm_markov_ml_control()`.
#' @param x,y Passed to [rms::orm()] for the returned fit.
#'
#' @return An object of class `orm_markov_ml`, inheriting from `orm`. Its
#'   coefficients are the maximum likelihood estimates and `$var` holds the
#'   patient-cluster robust covariance. `$ml` is a list with:
#'   * `loglik`: maximised log-likelihood \eqn{\sum_i \log L_i}.
#'   * `trace`: log-likelihood after each EM iteration (never decreasing).
#'   * `iterations`, `converged`, `newton_steps`: EM iterations, whether EM
#'     met its convergence criterion, and Newton steps taken.
#'   * `max_abs_score`: largest absolute derivative of the log-likelihood at
#'     the estimate; close to zero at a maximum.
#'   * `newton_decrement`: \eqn{g^\top J^{-1} g} for the score \eqn{g} at the
#'     estimate, roughly how much the log-likelihood could still rise; a
#'     warning is given when it exceeds `1e-6`.
#'   * `information`, `bread`: observed information \eqn{J} and its inverse.
#'   * `var_model`, `var_robust`: model-based and robust covariances.
#'   * `patient_scores`: one row per patient, the derivative of that patient's
#'     \eqn{\log L_i}; zeros for patients whose later states are all missing
#'     under `missing_states = "any"`.
#'   * `missing_info_fraction`: largest share of information lost because
#'     states are missing or known only up to a range, from 0 (none lost) to 1.
#'   * `start_coefficients`: estimates from transitions with both states
#'     observed (the [orm_markov()] fit on those rows).
#'   * `n_missing`, `n_range`, `n_gaps`, `n_trailing`: counts of rows with a
#'     missing state, rows with a range of states, runs of such rows that end
#'     in an observed state, and runs at the end of a patient's rows.
#'
#' @examples
#' \dontrun{
#' data <- prepare_markov_data(trajectories_with_missing, absorbing_state = 6)
#' fit <- orm_markov_ml(
#'   y ~ rms::rcs(time, 4) + tx + yprev,
#'   data = data,
#'   id_var = "id",
#'   absorbing_state = 6,
#'   missing_states = "any"
#' )
#' fit$ml$missing_info_fraction
#' avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:28) |>
#'   inferences(method = "delta")
#' }
#'
#' @seealso [orm_markov()], [avg_sops()], [inferences()]
#' @export
orm_markov_ml <- function(
  formula,
  data,
  ...,
  id_var,
  absorbing_state,
  p_var = "yprev",
  missing_states = c("any", "non_absorbing"),
  type = c("HC0", "HC1"),
  cadjust = NULL,
  time_var = "time",
  first_followup_time = NULL,
  control = orm_markov_ml_control(),
  x = TRUE,
  y = TRUE
) {
  refit_args <- markov_capture_refit_args(
    match.call(),
    environment(),
    names(formals(sys.function())),
    ...
  )
  type <- match.arg(type)
  missing_states <- match.arg(missing_states)
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be supplied as a data frame.", call. = FALSE)
  }
  if (missing(id_var) || is.null(id_var)) {
    stop(
      "`id_var` is required: each patient's rows form one Markov chain.",
      call. = FALSE
    )
  }
  id_var <- validate_markov_id_var(id_var, data, "data")
  if (
    missing(absorbing_state) ||
      is.null(absorbing_state) ||
      !length(absorbing_state) ||
      anyNA(absorbing_state)
  ) {
    stop(
      "`absorbing_state` is required: give the state(s) that are never left, ",
      "such as death.",
      call. = FALSE
    )
  }
  validate_markov_id_var(time_var, data, "data")
  if (
    !is.character(p_var) ||
      length(p_var) != 1L ||
      is.na(p_var) ||
      !p_var %in% names(data)
  ) {
    stop("`p_var` must name the previous-state column of `data`.", call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.", call. = FALSE)
  }
  mt <- stats::terms(formula, specials = "offset", data = data)
  if (terms_has_offset(mt)) {
    stop_unsupported_offset()
  }
  if (!p_var %in% all.vars(stats::delete.response(mt))) {
    stop(
      "The previous-state variable `",
      p_var,
      "` must appear on the right-hand side of `formula`.",
      call. = FALSE
    )
  }
  control <- markov_ml_validate_control(control)
  if (is.null(cadjust)) {
    cadjust <- TRUE
  }
  if (!is.logical(cadjust) || length(cadjust) != 1L || is.na(cadjust)) {
    stop("`cadjust` must be TRUE or FALSE.", call. = FALSE)
  }

  ml_call <- match.call(expand.dots = FALSE)
  dots <- ml_call$...
  markov_reject_legacy_wrapper_args(dots, "orm_markov_ml()")
  markov_ml_reject_dots(names(dots))
  dots_fit <- list(...)

  warn_duplicate_markov_id_time(
    data,
    id_var,
    time_var = time_var,
    wrapper = "orm_markov_ml()"
  )

  problem <- markov_ml_setup(
    formula = formula,
    data = data,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    absorbing_state = absorbing_state,
    missing_states = missing_states,
    dots_fit = dots_fit
  )
  theta_start <- stats::coef(problem$start_fit)

  em <- markov_ml_em(problem, theta_start, control)
  if (!em$converged) {
    warning(
      "EM did not converge within `maxit = ",
      control$maxit,
      "` iterations; the fit is marked as not converged.",
      call. = FALSE
    )
  }
  information <- markov_ml_information(problem, em$theta, control$info_step)
  newton <- markov_ml_newton(
    problem,
    em$theta,
    em$state$loglik,
    information,
    control
  )
  theta_hat <- newton$theta
  state <- markov_ml_estep(problem, theta_hat)
  patient_scores <- markov_ml_patient_scores(problem, theta_hat, state)
  score <- colSums(patient_scores)
  max_abs_score <- max(abs(score))
  # Scale-free check that the estimate is a maximum: the Newton decrement
  # g' J^{-1} g approximates how much the log-likelihood could still rise.
  newton_decrement <- tryCatch(
    sum(score * solve(information, score)),
    error = function(e) NA_real_
  )
  if (!is.finite(newton_decrement) || newton_decrement > 1e-6) {
    warning(
      "The fit may not be at the maximum of the likelihood (Newton ",
      "decrement ",
      signif(newton_decrement, 3),
      "). Increase `maxit` or `newton_steps` in `orm_markov_ml_control()`.",
      call. = FALSE
    )
  }

  var_model <- NULL
  var_robust <- NULL
  information_pd <- NA
  if (isTRUE(control$variance)) {
    if (newton$steps > 0L) {
      information <- markov_ml_information(
        problem,
        theta_hat,
        control$info_step
      )
    }
    var_model <- markov_ml_invert_information(information)
    information_pd <- !is.null(var_model)
    if (!information_pd) {
      warning(
        "The observed information is not positive definite; covariances are ",
        "set to NA. The model may not be identified from the observed data.",
        call. = FALSE
      )
      var_model <- matrix(NA_real_, nrow(information), ncol(information))
      dimnames(var_model) <- dimnames(information)
    }
  }

  fit <- markov_ml_final_fit(
    problem,
    state,
    theta_hat,
    score = score,
    dots_fit = dots_fit,
    x = x,
    y = y
  )
  complete_information <- fit$ml_complete_information
  fit$ml_complete_information <- NULL

  stored <- markov_prepare_stored_data(
    data = data,
    formula = formula,
    subset = NULL,
    eval_env = parent.frame(),
    fit_data = data,
    id_var = id_var,
    time_var = time_var,
    first_followup_time = first_followup_time
  )
  fitted_ids <- stored$starting_profile_metadata$fitted_ids
  patient_order <- match(fitted_ids, problem$patient_ids)
  patient_scores <- patient_scores[patient_order, , drop = FALSE]
  contributing <- as.vector(
    rowsum(as.numeric(problem$informative), problem$row_patient, reorder = TRUE)
  )[patient_order] >
    0
  names(contributing) <- fitted_ids
  n_informative_rows <- sum(problem$informative)

  robust <- NULL
  if (isTRUE(control$variance) && isTRUE(information_pd)) {
    robust <- markov_ml_sandwich(
      patient_scores,
      bread = var_model,
      n_rows = n_informative_rows,
      contributing = contributing,
      type = type,
      cadjust = cadjust
    )
    var_robust <- robust$covariance
  } else if (isTRUE(control$variance)) {
    var_robust <- var_model
  }

  counts <- markov_ml_missing_counts(problem)
  fit$ml <- c(
    list(
      loglik = state$loglik,
      trace = em$trace,
      iterations = em$iterations,
      converged = em$converged,
      newton_steps = newton$steps,
      max_abs_score = max_abs_score,
      newton_decrement = newton_decrement,
      information = information,
      information_positive_definite = information_pd,
      bread = var_model,
      var_model = var_model,
      var_robust = var_robust,
      patient_scores = patient_scores,
      contributing = contributing,
      n_informative_rows = n_informative_rows,
      missing_info_fraction = if (isTRUE(information_pd)) {
        markov_ml_missing_info_fraction(complete_information, information)
      } else {
        NA_real_
      },
      complete_information = complete_information,
      start_coefficients = theta_start,
      states = problem$labels,
      absorbing_state = problem$labels[problem$absorbing],
      missing_states = missing_states,
      p_var = p_var,
      id_var = id_var,
      type = type,
      cadjust = cadjust,
      adjustment_factor = robust$adjustment_factor %||% NA_real_
    ),
    counts
  )

  if (isTRUE(control$variance)) {
    fit$var <- var_robust
    fit$orig.var <- var_model
    fit$clusterInfo <- list(name = id_var, n = sum(contributing))
    attr(fit, "markov_robust_covariance") <- list(
      implementation = "mostr::orm_markov_ml",
      backend = "orm",
      cluster = data[[id_var]],
      type = type,
      cadjust = cadjust,
      adjustment_factor = fit$ml$adjustment_factor,
      n_rows = n_informative_rows,
      n_represented_rows = n_informative_rows,
      n_parameters = length(theta_hat),
      n_clusters = sum(contributing),
      bread_convention = "observed_data_information",
      weighted = FALSE,
      covariance_identity = var_robust
    )
  }

  fit$call <- match.call()
  fit$call[[1L]] <- quote(mostr::orm_markov_ml)
  fit <- markov_attach_model_data(
    fit,
    data = data,
    id_var = id_var,
    refit_data = stored$refit_data,
    starting_profile_data = stored$starting_profile_data,
    starting_profile_metadata = stored$starting_profile_metadata,
    refit_args = refit_args
  )
  fit <- markov_set_fit_wrapper(fit, "orm_markov_ml")
  class(fit) <- unique(c("orm_markov_ml", class(fit)))
  fit
}

#' Control Settings for `orm_markov_ml()`
#'
#' @param maxit Maximum number of EM iterations. A fit that reaches it is
#'   marked as not converged.
#' @param reltol EM stops when the relative change in log-likelihood falls
#'   below this value, and Newton steps finish the maximisation. When many
#'   states are missing, EM converges slowly, while Newton steps converge
#'   quickly once close to the maximum.
#' @param coef_tol EM also stops when no coefficient changes by more than this
#'   value.
#' @param newton_steps Maximum number of Newton steps after EM. They stop
#'   early once the Newton decrement (roughly how much the log-likelihood could
#'   still rise) is below `1e-20`.
#' @param info_step Relative step size for the numerical derivative of the
#'   score that gives the observed information.
#' @param trace Print the log-likelihood after each EM iteration?
#' @param variance Compute the observed information and covariances? Bootstrap
#'   refits set this to `FALSE` because they need only the coefficients.
#'
#' @return A list of settings.
#' @export
orm_markov_ml_control <- function(
  maxit = 200L,
  reltol = 1e-6,
  coef_tol = 1e-7,
  newton_steps = 10L,
  info_step = 1e-4,
  trace = FALSE,
  variance = TRUE
) {
  markov_ml_validate_control(list(
    maxit = maxit,
    reltol = reltol,
    coef_tol = coef_tol,
    newton_steps = newton_steps,
    info_step = info_step,
    trace = trace,
    variance = variance
  ))
}

markov_ml_validate_control <- function(control) {
  if (!is.list(control)) {
    stop("`control` must come from `orm_markov_ml_control()`.", call. = FALSE)
  }
  defaults <- list(
    maxit = 200L,
    reltol = 1e-6,
    coef_tol = 1e-7,
    newton_steps = 10L,
    info_step = 1e-4,
    trace = FALSE,
    variance = TRUE
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown `control` setting(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  control <- utils::modifyList(defaults, control)
  positive_number <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }
  count <- function(x, min) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= min && x == round(x)
  }
  flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!count(control$maxit, 1)) stop("`maxit` must be a positive integer.")
  if (!count(control$newton_steps, 0)) {
    stop("`newton_steps` must be a non-negative integer.")
  }
  for (nm in c("reltol", "coef_tol", "info_step")) {
    if (!positive_number(control[[nm]])) {
      stop("`", nm, "` must be a positive number.")
    }
  }
  if (!flag(control$trace) || !flag(control$variance)) {
    stop("`trace` and `variance` must be TRUE or FALSE.")
  }
  control$maxit <- as.integer(control$maxit)
  control$newton_steps <- as.integer(control$newton_steps)
  control
}

markov_ml_reject_dots <- function(dot_names) {
  rejected <- c(
    "weights",
    "subset",
    "offset",
    "penalty",
    "penalty.matrix",
    "var.penalty",
    "cluster",
    "p2_var",
    "family",
    "initial",
    "method"
  )
  bad <- intersect(dot_names, rejected)
  if ("p2_var" %in% bad) {
    stop(
      "`orm_markov_ml()` fits first-order models only and has no `p2_var`; ",
      "the transition may depend on the past only through `p_var`.",
      call. = FALSE
    )
  }
  if (length(bad)) {
    stop(
      "`orm_markov_ml()` does not support the argument(s) ",
      paste0("`", bad, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Replace rcs() knot counts with the knots chosen by the starting fit so a
# refit on other data uses the same design.
markov_ml_lock_knots <- function(formula, design) {
  parms <- design$parms
  is_rcs <- function(fn) {
    identical(fn, as.name("rcs")) ||
      (is.call(fn) &&
        length(fn) == 3L &&
        (identical(fn[[1L]], as.name("::")) ||
          identical(fn[[1L]], as.name(":::"))) &&
        identical(fn[[3L]], as.name("rcs")))
  }
  rewrite <- function(e) {
    if (!is.call(e)) {
      return(e)
    }
    if (is_rcs(e[[1L]])) {
      variable <- paste(deparse(e[[2L]]), collapse = "")
      knots <- parms[[variable]]
      if (!is.numeric(knots) || length(knots) < 3L) {
        stop(
          "Could not find the spline knots chosen for `",
          variable,
          "`; give the knots explicitly in `rcs()`.",
          call. = FALSE
        )
      }
      return(as.call(list(e[[1L]], e[[2L]], as.numeric(knots))))
    }
    for (i in seq_along(e)[-1L]) {
      if (!is.null(e[[i]])) {
        e[[i]] <- rewrite(e[[i]])
      }
    }
    e
  }
  formula[[3L]] <- rewrite(formula[[3L]])
  formula
}

# Weighted rms::orm() evaluated at the ML estimates on the expected transition
# counts, returning a full orm object for prediction. It does not
# iterate (`maxit = 1` returns the fit at `initial`). By Fisher's identity its
# gradient, the posterior-weighted complete-data score, equals the
# observed-data score at `theta`; comparing the two checks that the refit
# describes the same model.
markov_ml_final_fit <- function(problem, state, theta, score, dots_fit, x, y) {
  # One row per (unique design row, destination state) with the expected
  # number of such transitions as its weight. The gradient and information of
  # a weighted orm fit depend only on these sums, so this matches a fit on one
  # row per patient visit and possible transition, with far fewer rows.
  counts <- markov_ml_aggregate_weights(problem, state)
  expanded <- problem$design_data[counts$design, , drop = FALSE]
  expanded[[markov_ml_response_column()]] <- markov_ml_state_values(
    counts$state,
    problem$scale,
    problem$response_type
  )
  expanded[[markov_ml_weight_column()]] <- counts$weight
  rownames(expanded) <- NULL

  formula <- markov_ml_lock_knots(
    problem$working_formula,
    problem$start_fit$Design
  )
  # Scaling and optimiser settings only matter when iterating.
  refit_args <- dots_fit[setdiff(names(dots_fit), c("maxit", "scale", "opt_method"))]
  refit_args$maxit <- 1L
  fit <- markov_ml_fit_orm(
    formula,
    expanded,
    refit_args,
    weights = markov_ml_weight_column(),
    initial = theta
  )
  if (
    isTRUE(fit$fail) ||
      !identical(names(stats::coef(fit)), problem$coef_names) ||
      max(abs(stats::coef(fit) - theta)) > 1e-10 * max(1, max(abs(theta)))
  ) {
    stop("Internal error: the weighted refit is not evaluated at the estimates.")
  }
  gradient_gap <- max(abs(fit$u - score))
  if (
    length(fit$u) != length(score) ||
      !is.finite(gradient_gap) ||
      gradient_gap > 1e-6 * max(1, max(abs(fit$u)), max(abs(score)))
  ) {
    stop(
      "Internal error: the weighted refit's gradient differs from the ",
      "observed-data score (largest difference ",
      signif(gradient_gap, 3),
      ").",
      call. = FALSE
    )
  }
  X_final <- orm_model_matrix(fit, problem$design_data)
  if (
    !identical(dim(X_final), dim(problem$X)) ||
      max(abs(X_final - problem$X)) > 1e-10 * max(1, max(abs(problem$X)))
  ) {
    stop(
      "Internal error: the weighted refit uses a different design from the ",
      "starting fit. Give spline knots and factor levels explicitly.",
      call. = FALSE
    )
  }
  complete_information <- tryCatch(
    solve(stats::vcov(fit, intercepts = "all")),
    error = function(e) NULL
  )
  if (!is.null(complete_information)) {
    complete_information <- (complete_information +
      t(complete_information)) /
      2
    dimnames(complete_information) <- list(problem$coef_names, problem$coef_names)
  }
  fit$coefficients <- theta
  fit$info.matrix <- NULL
  fit$var <- NULL
  if (!isTRUE(x)) {
    fit$x <- NULL
  }
  if (!isTRUE(y)) {
    fit$y <- NULL
  }
  fit$ml_complete_information <- complete_information
  fit
}

markov_ml_missing_counts <- function(problem) {
  exact <- problem$scale$exact
  missing <- problem$scale$missing
  n_range <- sum(!exact & !missing)
  # Runs of rows without an exact state; a run is a gap when an exact row
  # follows it within the same patient and trailing otherwise.
  next_exact <- c(exact[-1L], FALSE)
  next_same <- c(problem$row_patient[-1L], NA) == problem$row_patient
  next_same[is.na(next_same)] <- FALSE
  run_end <- !exact & (!next_same | next_exact)
  closes_gap <- run_end & next_same & next_exact
  list(
    n_missing = sum(missing),
    n_range = n_range,
    n_gaps = sum(closes_gap),
    n_trailing = sum(run_end & !next_same),
    n_patients = problem$N,
    n_rows = problem$n
  )
}

# The rms::orm() print method is not used: its statistics (frequencies, LR
# chi-square, R2, rank correlations) would describe the weighted refit on
# expected transition counts from markov_ml_final_fit(), not the data.
#' @export
print.orm_markov_ml <- function(x, digits = 4, ...) {
  ml <- x$ml
  cat("Observed-data maximum likelihood Markov fit (orm_markov_ml)\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat(
    "  Log-likelihood: ",
    format(ml$loglik, nsmall = 4),
    "\n  EM iterations: ",
    ml$iterations,
    if (isTRUE(ml$converged)) " (converged)" else " (NOT converged)",
    "; Newton steps: ",
    ml$newton_steps,
    "; max |score|: ",
    format(signif(ml$max_abs_score, 3)),
    "\n  Patients: ",
    ml$n_patients,
    "; clusters for robust standard errors: ",
    sum(ml$contributing),
    "\n  Rows: ",
    ml$n_rows,
    " (",
    ml$n_missing,
    " missing, ",
    ml$n_range,
    " with a range of states); gaps: ",
    ml$n_gaps,
    "; trailing runs: ",
    ml$n_trailing,
    "\n  Missing states mean: ",
    ml$missing_states,
    "; missing-information fraction: ",
    format(signif(ml$missing_info_fraction, 3)),
    "\n\n",
    sep = ""
  )
  theta <- stats::coef(x)
  if (is.null(x$var)) {
    cat("Coefficients (standard errors not computed: `variance = FALSE`):\n")
    print(cbind(Coef = round(theta, digits)))
    return(invisible(x))
  }
  se <- sqrt(diag(stats::vcov(x, intercepts = "all")))
  z <- theta / se
  table <- cbind(
    Coef = theta,
    `Robust S.E.` = se,
    `Wald Z` = z,
    `Pr(>|Z|)` = 2 * stats::pnorm(-abs(z))
  )
  cat(
    "Coefficients (patient-cluster robust standard errors, ",
    ml$type,
    if (isTRUE(ml$cadjust)) " with n/(n-1)",
    "):\n",
    sep = ""
  )
  if (!isTRUE(ml$information_positive_definite)) {
    cat(
      "  (unavailable: the observed information is not positive definite)\n"
    )
  }
  stats::printCoefmat(
    table,
    digits = digits,
    signif.stars = FALSE,
    P.values = TRUE,
    has.Pvalue = TRUE,
    eps.Pvalue = 1e-4,
    na.print = "NA"
  )
  invisible(x)
}

# rms::logLik.rms() would use the deviance of the weighted refit on expected
# transition counts; AIC.rms() and stats::BIC() call logLik().
#' @export
logLik.orm_markov_ml <- function(object, ...) {
  structure(
    object$ml$loglik,
    nobs = stats::nobs(object),
    df = length(stats::coef(object)),
    class = "logLik"
  )
}

# Rows whose observation rules out at least one state: the transitions
# orm_markov() counts when nothing is missing, and `m` of the HC1 correction.
#' @export
nobs.orm_markov_ml <- function(object, ...) {
  object$ml$n_informative_rows
}
