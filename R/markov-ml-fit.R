# Maximum likelihood fits of the orm and vglm engines with missing states.
#
# mostm() calls mostm_fit_ml() when `missing_states` is given or the response
# is an Ocens() range and some rows need the forward-backward recursion. The
# engine fits only exactly observed states: first the transitions with both
# states observed (starting values and the design), and at the end the
# expected transition counts at the estimates, which gives an ordinary engine
# object for prediction. Everything in between is computed from the linear
# predictors (markov-ml-lp.R), the same way for both engines.

# Replace the number of knots in rcs() terms with the knots that rms chooses
# for `data`, so the start fit, the design rows of the possible previous
# states, and the final fit on expected counts share one spline basis.
markov_ml_lock_knots <- function(formula, data) {
  env <- environment(formula) %||% globalenv()
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
      basis_call <- e
      basis_call[[1L]] <- quote(rms::rcs)
      basis <- eval(basis_call, data, env)
      knots <- attr(basis, "parms")
      if (!is.numeric(knots) || length(knots) < 3L) {
        stop(
          "Could not find the spline knots of `",
          paste(deparse(e), collapse = ""),
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

# Engine arguments that do not apply to the maximum likelihood fit.
markov_ml_reject_engine_args <- function(engine, dot_names) {
  rejected <- switch(
    engine,
    orm = c(
      "penalty",
      "penalty.matrix",
      "var.penalty",
      "initial",
      "na.action"
    ),
    vglm = c("na.action", "etastart", "mustart", "form2"),
    blrm = c(
      "na.action",
      "iprior",
      "ascale",
      "priorsdppo",
      "backend",
      "loo",
      "ppairs",
      "method",
      "inito",
      "inits",
      "standata",
      "file",
      "progress",
      "debug",
      "showopt"
    )
  )
  bad <- intersect(dot_names, rejected)
  if (length(bad)) {
    stop(
      "With `missing_states` or an `Ocens()` response, `engine = \"",
      engine,
      "\"` does not accept ",
      paste0("`", bad, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Run an engine fit quietly: rms warns that orm ignores weights in its own
# bootstrap, and orm.fit() reports non-convergence as a message; callers
# check the result.
markov_ml_quiet <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (
        grepl(
          "weights are ignored in model validation",
          conditionMessage(w),
          fixed = TRUE
        )
      ) {
        invokeRestart("muffleWarning")
      }
    },
    message = function(m) {
      if (
        grepl(
          "without improving log likelihood|iterations without convergence",
          conditionMessage(m)
        )
      ) {
        invokeRestart("muffleMessage")
      }
    }
  )
}

# ---------------------------------------------------------------------------
# Engine adapters
# ---------------------------------------------------------------------------

# rms::orm() on `data` with optional weight column and starting values.
markov_ml_orm <- function(formula, data, engine_args, weights = NULL, initial = NULL) {
  fit_call <- as.call(c(
    list(quote(rms::orm)),
    list(formula = formula, data = quote(data), x = TRUE, y = TRUE),
    engine_args
  ))
  if (!is.null(weights)) {
    fit_call$weights <- as.name(weights)
  }
  if (!is.null(initial)) {
    fit_call$initial <- unname(initial)
  }
  fit <- markov_ml_quiet(eval(fit_call))
  if (isTRUE(fit$fail) && is.null(engine_args$opt_method)) {
    # Retry near-separated data with Levenberg-Marquardt.
    fit_call$opt_method <- "LM"
    fit <- markov_ml_quiet(eval(fit_call))
  }
  fit
}

# Starting fit on the transitions with both states observed and the design
# rows of the possible previous states.
markov_ml_start <- function(problem, engine, setup, engine_args) {
  start_data <- markov_ml_with_state(
    problem,
    problem$start_data,
    problem$complete_part$state
  )
  if (engine == "vglm") {
    start_args <- engine_args
    start_setup <- list(
      formula = problem$working_formula,
      data = start_data,
      dots = setup$dots,
      subset_rows = NULL,
      prior_weight_column = NULL,
      eval_env = setup$eval_env,
      data_name = character()
    )
    fit <- do.call(mostm_fit_vglm, c(list(start_setup), start_args))$fit
    fit <- markov_set_engine(fit, "vglm")
    validate_markov_model(fit)
    constraints <- VGAM::constraints(fit)
    X <- markov_ml_vglm_design(fit, problem$design_data, names(constraints))
    theta <- stats::coef(fit)
  } else {
    # The blrm engine uses the rms design for its own sampler.
    orm_args <- if (engine == "orm") engine_args else list()
    fit <- markov_ml_orm(problem$working_formula, start_data, orm_args)
    if (isTRUE(fit$fail)) {
      stop(
        "`rms::orm()` could not fit the transitions with both states ",
        "observed, which give the starting values.",
        call. = FALSE
      )
    }
    fit <- markov_set_engine(fit, "orm")
    validate_markov_model(fit)
    M <- problem$K - 1L
    X <- cbind(
      "(Intercept)" = 1,
      orm_model_matrix(fit, problem$design_data)
    )
    constraints <- c(
      list(diag(M)),
      rep(list(matrix(1, nrow = M, ncol = 1L)), ncol(X) - 1L)
    )
    theta <- stats::coef(fit)
  }
  if (length(theta) != sum(vapply(constraints, ncol, integer(1)))) {
    stop("Internal error: the starting fit has unexpected coefficients.")
  }
  if (any(!is.finite(X))) {
    stop("The design matrix for possible previous states is not finite.")
  }
  K <- problem$K
  n_levels <- if (engine == "vglm") VGAM::npred(fit) + 1L else length(fit$yunique)
  if (n_levels != K) {
    stop("Internal error: the starting fit does not use every state.")
  }
  list(
    fit = fit,
    theta = theta,
    lp = markov_lp(X, constraints, names(theta))
  )
}

# Small design matrix of a vglm fit for new rows, with the columns of its
# constraint matrices.
markov_ml_vglm_design <- function(fit, data, columns) {
  tt <- stats::delete.response(stats::terms(fit))
  X <- stats::model.matrix(
    tt,
    data = data,
    contrasts.arg = if (length(fit@contrasts)) fit@contrasts else NULL,
    xlev = fit@xlevels
  )
  absent <- setdiff(columns, colnames(X))
  if (length(absent)) {
    stop(
      "Internal error: design columns ",
      paste(absent, collapse = ", "),
      " are missing for the possible previous states."
    )
  }
  X[, columns, drop = FALSE]
}

# One row per (design row, destination state) with the expected number of
# such transitions at the estimates as its weight.
markov_ml_expected_rows <- function(problem, counts) {
  expanded <- markov_ml_with_state(
    problem,
    problem$design_data[counts$design, , drop = FALSE],
    counts$state
  )
  expanded[[markov_ml_weight_column()]] <- counts$weight
  rownames(expanded) <- NULL
  expanded
}

# Engine fit on the expected transition counts at the estimates `theta`. By
# Fisher's identity its score, the posterior-weighted complete-data score,
# equals the observed-data score, so `theta` maximises it as well; the
# returned object holds exactly `theta` and is used for prediction.
markov_ml_final_fit <- function(
  problem,
  engine,
  setup,
  counts,
  theta,
  engine_args,
  converged = TRUE
) {
  expanded <- markov_ml_expected_rows(problem, counts)
  if (engine == "orm") {
    # Evaluate the fit at `theta` without iterating (`maxit = 1`); the
    # optimizer settings only matter when iterating.
    args <- engine_args[setdiff(names(engine_args), c("maxit", "scale", "opt_method"))]
    args$maxit <- 1L
    fit <- markov_ml_orm(
      problem$working_formula,
      expanded,
      args,
      weights = markov_ml_weight_column(),
      initial = theta
    )
    if (
      isTRUE(fit$fail) ||
        !identical(names(stats::coef(fit)), problem$lp$coef_names) ||
        max(abs(stats::coef(fit) - theta)) > 1e-10 * max(1, max(abs(theta)))
    ) {
      stop("Internal error: the orm fit is not evaluated at the estimates.")
    }
    X_final <- cbind(1, orm_model_matrix(fit, problem$design_data))
    fit$coefficients <- theta
    fit$info.matrix <- NULL
  } else {
    args <- engine_args
    args$coefstart <- unname(theta)
    final_setup <- list(
      formula = problem$working_formula,
      data = expanded,
      dots = setup$dots,
      subset_rows = NULL,
      prior_weight_column = markov_ml_weight_column(),
      eval_env = setup$eval_env,
      data_name = character()
    )
    # Starting at the estimates, VGAM can only take tiny steps and warns about
    # them; the result is checked against `theta` below.
    fit <- suppressWarnings(
      do.call(mostm_fit_vglm, c(list(final_setup), args))$fit
    )
    # Only a maximum of the likelihood also maximises the fit to the expected
    # counts; a fit that did not converge has been warned about.
    gap <- max(abs(stats::coef(fit) - theta))
    if (
      !identical(names(stats::coef(fit)), problem$lp$coef_names) ||
        (converged && (!is.finite(gap) || gap > 1e-4 * max(1, max(abs(theta)))))
    ) {
      stop(
        "Internal error: the vglm fit to the expected transition counts does ",
        "not reproduce the estimates (largest difference ",
        signif(gap, 3),
        ").",
        call. = FALSE
      )
    }
    fit <- set_coef(fit, theta)
    X_final <- markov_ml_vglm_design(
      fit,
      problem$design_data,
      names(VGAM::constraints(fit))
    )
  }
  X <- problem$lp$X
  if (
    !identical(dim(unname(X_final)), dim(unname(X))) ||
      max(abs(X_final - X)) > 1e-10 * max(1, max(abs(X)))
  ) {
    stop(
      "Internal error: the fit to the expected transition counts uses a ",
      "different design from the starting fit. Give spline knots and ",
      "factor levels explicitly.",
      call. = FALSE
    )
  }
  fit
}

# ---------------------------------------------------------------------------
# Maximum likelihood fit
# ---------------------------------------------------------------------------

# Maximise the observed-data likelihood and compute the components of the
# robust covariance. Returns the engine fit at the estimates and a list `ml`
# with the likelihood, information, and patient scores.
mostm_fit_ml <- function(
  setup,
  engine,
  problem,
  type,
  cadjust,
  control,
  engine_args
) {
  start <- markov_ml_start(problem, engine, setup, engine_args)
  problem$lp <- start$lp
  theta_start <- start$theta

  state <- NULL
  if (problem$chain$n) {
    em <- markov_ml_em(problem, theta_start, control)
    if (!em$converged) {
      warning(
        "The EM algorithm did not converge within ",
        control$maxit,
        " iterations; Newton steps continue from its last estimates.",
        call. = FALSE
      )
    }
    information <- markov_ml_information(problem, em$theta, em$state)
    newton <- markov_ml_newton(
      problem,
      em$theta,
      em$state$loglik,
      information,
      control
    )
    theta_hat <- newton$theta
    state <- markov_ml_estep(problem, theta_hat)
    loglik <- state$loglik
  } else {
    # Only exact and range-free transitions: the start fit maximises the
    # likelihood of these transitions.
    theta_hat <- theta_start
    loglik <- markov_ml_loglik(problem, theta_hat)
    em <- list(trace = loglik, iterations = 0L, converged = TRUE)
    newton <- list(theta = theta_hat, loglik = loglik, steps = 0L)
    information <- markov_ml_complete_information(problem, theta_hat)
  }
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
      ").",
      call. = FALSE
    )
  }

  var_model <- NULL
  information_pd <- NA
  if (isTRUE(control$variance)) {
    if (newton$steps > 0L) {
      information <- markov_ml_information(problem, theta_hat, state)
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

  counts <- if (problem$chain$n) {
    markov_ml_aggregate_weights(problem, state)
  } else {
    comp <- problem$complete_part
    list(design = comp$agg_design, state = comp$agg_state, weight = comp$agg_count)
  }
  fit <- markov_ml_final_fit(
    problem,
    engine,
    setup,
    counts,
    theta_hat,
    engine_args,
    converged = is.finite(newton_decrement) && newton_decrement <= 1e-6
  )
  complete_information <- if (problem$chain$n) {
    eta <- markov_lp_eta(problem$lp, theta_hat)
    markov_lp_outcome_information(
      problem$lp,
      markov_lp_outcome(problem$lp, eta, counts$design, counts$state),
      counts$weight
    )
  } else {
    information
  }

  n_informative_rows <- sum(problem$informative)
  contributing <- as.vector(
    rowsum(as.numeric(problem$informative), problem$row_patient, reorder = TRUE)
  ) > 0
  names(contributing) <- problem$patient_ids

  var_robust <- NULL
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

  ml <- c(
    list(
      loglik = loglik,
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
      missing_states = problem$labels[problem$missing_index],
      p_var = problem$p_var,
      type = type,
      cadjust = cadjust,
      adjustment_factor = robust$adjustment_factor %||% NA_real_
    ),
    markov_ml_missing_counts(problem)
  )
  list(fit = fit, ml = ml)
}

# Rows of `data` used by the fit, in their original order.
markov_ml_fit_data <- function(problem, data) {
  data[sort(problem$source_rows), , drop = FALSE]
}
