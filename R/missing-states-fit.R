# Fitting when some states are not observed: the maximum likelihood fit for
# `engine = "orm"` and the pieces shared with `engine = "blrm"`.
#
# The model structure (design matrix columns, spline knots, factor levels,
# state labels) and the starting values come from the engine's own fit to the
# transitions whose state and previous state are both observed. mostr then
# maximises the likelihood that also includes the stretches of unobserved
# visits (`markov_path_plan()`), and puts its estimates into that fit, so that
# everything downstream reads an ordinary `orm` object.

# Arguments of `mostm()` that the path likelihood has not been worked out for.
mostm_check_path_fit_args <- function(engine, dots, p2_var, formula) {
  refuse <- function(what, why) {
    mostr_abort(
      what,
      " can't be combined with unobserved states (a missing state, a missing ",
      "previous state, or a range of states given with `Ocens()`): ",
      why,
      " Use `missing_states = \"drop\"` to fit only the visits whose state ",
      "and previous state are both observed."
    )
  }
  if (!is.null(p2_var)) {
    refuse(
      "A second-order model (`p2_var`)",
      "the sum over the states at unobserved visits is implemented for ",
      "models in which the next state depends on the previous state only."
    )
  }
  if (engine == "orm") {
    for (arg in intersect(names(dots), c("penalty", "penalty.matrix"))) {
      refuse(
        paste0("`", arg, "`"),
        "the penalised likelihood is not implemented for them."
      )
    }
    family <- dots$family
    if (!is.null(family) && !identical(family, "logistic")) {
      refuse(
        "A link other than the logit",
        "only the logit link is implemented for them."
      )
    }
  }
  if (engine == "blrm") {
    for (arg in intersect(names(dots), c("ppo", "cppo"))) {
      refuse(
        paste0("A partial proportional odds model (`", arg, "`)"),
        "the Stan program for them fits the proportional odds model only."
      )
    }
    if (length(formula_rhs_calls(formula, "cluster"))) {
      refuse(
        "A random intercept (`cluster()`)",
        "the Stan program for them has no random effects."
      )
    }
    if (isTRUE(dots$iprior != 0) && !is.null(dots$iprior)) {
      refuse(
        "`iprior`",
        "the Stan program for them uses the Dirichlet prior on the ",
        "intercepts, the default of `rmsb::blrm()`."
      )
    }
  }
  invisible(NULL)
}

# The formula and data for the engine's fit to the fully observed transitions.
# With `Ocens()` the engine must not see the ranges, so the response becomes a
# column that holds the exactly observed states and `NA` elsewhere.
markov_observed_transition_setup <- function(setup, ranges) {
  if (!ranges$ocens) {
    return(setup)
  }
  exact <- !is.na(ranges$lo) & ranges$lo == ranges$hi
  state <- rep(NA_character_, nrow(setup$data))
  state[exact] <- ranges$labels[ranges$lo[exact]]
  setup$data$.markov_state <- factor(
    state,
    levels = ranges$labels,
    ordered = TRUE
  )
  setup$formula[[2L]] <- quote(.markov_state)
  setup
}

markov_check_observed_states <- function(fitted_labels, labels) {
  if (identical(as.character(fitted_labels), labels)) {
    return(invisible(NULL))
  }
  mostr_abort(
    "The state",
    if (length(setdiff(labels, fitted_labels)) > 1L) "s " else " ",
    paste(setdiff(labels, fitted_labels), collapse = ", "),
    " never occur",
    if (length(setdiff(labels, fitted_labels)) == 1L) "s",
    " at a visit whose state and previous state are both observed, so the ",
    "model can't be set up. Merge rarely observed states, or check the ",
    "bounds given to `Ocens()`."
  )
}

# The likelihood, its score and its second derivative at `theta` (intercepts,
# then slopes). `observed` holds the fully observed transitions.
markov_path_loglik <- function(theta, observed, plan, order = 2L) {
  n_intercepts <- length(plan$labels) - 1L
  cpp_path_loglik(
    alpha = theta[seq_len(n_intercepts)],
    beta = theta[-seq_len(n_intercepts)],
    x_obs = observed$x,
    y_obs = observed$y,
    w_obs = observed$weight,
    group_obs = observed$group,
    n_groups = as.integer(observed$n_groups),
    x_gap = plan$x_gap,
    absorb = plan$absorb,
    stretch_first_row = plan$stretch_first_row,
    stretch_first_step = plan$stretch_first_step,
    stretch_n_steps = plan$stretch_n_steps,
    stretch_weight = plan$stretch_weight,
    step_lo = plan$step_lo,
    step_hi = plan$step_hi,
    step_rows = plan$step_rows,
    order = as.integer(order)
  )
}

# Score of each patient: the observed transitions plus the patient's
# stretches, each stretch's score counted with the patient's weight.
markov_path_patient_scores <- function(result, plan, patient_ids) {
  scores <- result$obs_group_scores
  if (nrow(plan$instances)) {
    stretch_scores <- result$stretch_scores[plan$instances$stretch, , drop = FALSE] *
      plan$instances$weight
    patient <- match(plan$instances$id, patient_ids)
    summed <- rowsum(stretch_scores, patient)
    rows <- as.integer(rownames(summed))
    scores[rows, ] <- scores[rows, , drop = FALSE] + summed
  }
  scores
}

# Newton-Raphson with step halving. The second derivative is exact, so the
# iteration converges quadratically from the fit to the observed transitions,
# usually in three to five steps. It stops when a full step would raise the
# log-likelihood by less than `tolerance` (half the Newton decrement
# score' * step): at 1e-10 the estimates are settled to far below their
# standard errors.
markov_path_maximise <- function(
  theta,
  observed,
  plan,
  maxit = 50L,
  tolerance = 1e-10
) {
  evaluate <- function(theta, order) {
    markov_path_loglik(theta, observed, plan, order)
  }
  total_score <- function(result) {
    colSums(result$obs_group_scores) +
      as.numeric(crossprod(result$stretch_scores, plan$stretch_weight))
  }
  current <- evaluate(theta, 2L)
  if (!is.finite(current$loglik)) {
    mostr_abort(
      "The likelihood with unobserved states is not finite at the starting ",
      "values, the estimates from the fully observed transitions."
    )
  }
  converged <- FALSE
  for (iteration in seq_len(maxit)) {
    score <- total_score(current)
    step <- tryCatch(
      solve(-current$hessian, score),
      error = function(e) NULL
    )
    if (is.null(step) || any(!is.finite(step))) {
      break
    }
    gain <- sum(score * step) / 2
    if (is.finite(gain) && abs(gain) < tolerance) {
      converged <- TRUE
      break
    }
    # A Newton step can overshoot far from the maximum; halve it until the
    # log-likelihood goes up.
    accepted <- FALSE
    for (halving in 0:30) {
      candidate <- theta + step / 2^halving
      trial <- evaluate(candidate, 0L)
      if (is.finite(trial$loglik) && trial$loglik >= current$loglik) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      break
    }
    theta <- candidate
    current <- evaluate(theta, 2L)
  }
  if (!converged) {
    mostr_abort(
      "The fit with unobserved states did not converge. Some coefficients ",
      "may have no finite estimate, for example because a transition never ",
      "occurs in the data: simplify the model or merge rarely observed ",
      "states."
    )
  }
  list(
    theta = theta,
    result = current,
    score = total_score(current),
    iterations = iteration - 1L
  )
}

# What a path fit records about itself, read by the covariance code and
# printed by nothing.
markov_path_fit_info <- function(model) {
  attr(model, "markov_path_fit", exact = TRUE)
}

markov_is_path_fit <- function(model) {
  !is.null(markov_path_fit_info(model))
}

# The fully observed transitions of an rms or rmsb fit as the kernel needs
# them, with the patients numbered in the order of `patient_ids`.
markov_path_observed <- function(x, y, labels, weights, ids, patient_ids) {
  x <- unname(as.matrix(x))
  storage.mode(x) <- "double"
  list(
    x = x,
    y = match(as.character(y), labels),
    weight = if (is.null(weights) || !length(weights)) {
      rep(1, nrow(x))
    } else {
      as.numeric(weights)
    },
    group = match(as.character(ids), patient_ids),
    n_groups = length(patient_ids)
  )
}

mostm_fit_orm_paths <- function(setup, ranges, path, ...) {
  observed_setup <- markov_observed_transition_setup(setup, ranges)
  fitted <- mostm_fit_orm(observed_setup, ...)
  start <- fitted$fit
  if (!identical(start$family, "logistic")) {
    stop_unsupported_markov_link(start$family)
  }
  labels <- ranges$labels
  markov_check_observed_states(start$yunique, labels)

  rows <- mostm_subset_row_index(setup$subset_rows, nrow(setup$data))
  data <- setup$data[rows, , drop = FALSE]
  plan <- markov_path_plan(
    data = data,
    ranges = list(
      labels = labels,
      lo = ranges$lo[rows],
      hi = ranges$hi[rows]
    ),
    absorb = match(as_state_labels(path$absorb), labels),
    id_var = path$id_var,
    time_var = path$time_var,
    p_var = path$p_var,
    weights = if (!is.null(setup$prior_weight_column)) {
      data[[setup$prior_weight_column]]
    },
    design = function(rows) orm_model_matrix(start, rows)
  )
  if (is.null(plan$x_gap)) {
    plan$x_gap <- start$x[0L, , drop = FALSE]
  }
  plan$x_gap <- unname(plan$x_gap)
  storage.mode(plan$x_gap) <- "double"

  fit_ids <- as.character(fitted$fit_data[[path$id_var]])
  patient_ids <- unique(c(fit_ids, plan$instances$id))
  observed <- markov_path_observed(
    start$x,
    start$y,
    labels,
    start$weights,
    fit_ids,
    patient_ids
  )
  theta <- stats::coef(start)
  maximum <- markov_path_maximise(unname(theta), observed, plan)

  names_theta <- names(theta)
  information <- -maximum$result$hessian
  dimnames(information) <- list(names_theta, names_theta)
  scores <- markov_path_patient_scores(maximum$result, plan, patient_ids)
  dimnames(scores) <- list(patient_ids, names_theta)

  fit <- start
  fit$coefficients <- stats::setNames(maximum$theta, names_theta)
  fit$u <- stats::setNames(maximum$score, names_theta)
  fit$iter <- maximum$iterations
  n_intercepts <- length(labels) - 1L
  fit$linear.predictors <- stats::setNames(
    fit$coefficients[[fit$interceptRef]] +
      as.numeric(observed$x %*% maximum$theta[-seq_len(n_intercepts)]),
    names(start$linear.predictors)
  )
  # rms reads `info.matrix` in `vcov()` and `print()`. Its list form holds the
  # intercept block, the slope block and the block between them; the
  # intercept block is given as a sparse symmetric matrix by its upper
  # triangle, because with stretches it is no longer tridiagonal.
  intercepts <- seq_len(n_intercepts)
  upper <- which(upper.tri(diag(n_intercepts), diag = TRUE), arr.ind = TRUE)
  fit$info.matrix <- list(
    a = list(
      row = upper[, 1L],
      col = upper[, 2L],
      a = information[intercepts, intercepts, drop = FALSE][upper]
    ),
    b = unname(information[-intercepts, -intercepts, drop = FALSE]),
    ab = unname(information[intercepts, -intercepts, drop = FALSE]),
    iname = names_theta[intercepts],
    xname = names_theta[-intercepts]
  )
  fit$var <- NULL
  fit$deviance <- c(NA_real_, -2 * maximum$result$loglik)
  names(fit$deviance) <- names(start$deviance)
  # These statistics of `rms::orm()` describe its fit to the fully observed
  # transitions (likelihood ratio and score tests against the model without
  # predictors, rank measures of the fitted rows). They are not computed for
  # the likelihood with unobserved states, and stale values would mislead.
  stale <- setdiff(names(fit$stats), c("Obs", "ESS", "Distinct Y", "Median Y", "d.f."))
  fit$stats[stale] <- NA_real_
  fit$stats[["Max Deriv"]] <- max(abs(maximum$score))

  attr(fit, "markov_path_fit") <- list(
    ids = patient_ids,
    scores = scores,
    information = information,
    loglik = maximum$result$loglik,
    iterations = maximum$iterations,
    n_observed_transitions = nrow(observed$x),
    n_stretches = nrow(plan$instances),
    n_distinct_stretches = length(plan$stretch_n_steps),
    n_unobserved_visits = plan$n_steps_total,
    n_likelihood_terms = sum(observed$weight > 0) +
      sum(plan$instances$weight > 0)
  )
  list(
    fit = fit,
    fit_data = fitted$fit_data,
    fitted_ids = patient_ids
  )
}

# The patient-clustered covariance of a path fit, stored on the fit in the
# form `robcov_orm()` stores it. `robcov_orm()` itself can't be used: it
# rebuilds one score per fitted row from `x` and `y`, and a stretch of
# unobserved visits has one score for all its rows together.
robcov_orm_paths <- function(fit, cluster, type, cadjust) {
  info <- markov_path_fit_info(fit)
  coefficients <- stats::coef(fit)
  information <- info$information
  bread <- tryCatch(solve(information), error = function(e) NULL)
  if (is.null(bread) || !markov_correlation_positive_definite(bread)) {
    mostr_abort(
      "The information matrix of the fit with unobserved states can't be ",
      "inverted, so standard errors are not available. Simplify the model, ",
      "for example by removing predictors that are nearly collinear or by ",
      "merging rarely observed states."
    )
  }
  bread <- (bread + t(bread)) / 2
  dimnames(bread) <- list(names(coefficients), names(coefficients))
  if (is.null(cadjust)) {
    cadjust <- TRUE
  }
  n_clusters <- nrow(info$scores)
  if (n_clusters < 2L) {
    mostr_abort(
      "Cluster-robust covariance requires at least two patients."
    )
  }
  assembled <- markov_assemble_sandwich(
    sandwich_scores = info$scores,
    bread = bread,
    n = info$n_likelihood_terms,
    p = length(coefficients),
    n_clusters = n_clusters,
    type = type,
    cadjust = cadjust
  )
  dimnames(assembled$covariance) <- dimnames(bread)

  fit$orig.var <- bread
  fit$var <- assembled$covariance
  fit$clusterInfo <- list(name = "cluster", n = n_clusters)
  attr(fit, "markov_robust_covariance") <- list(
    implementation = "mostr::robcov_orm",
    backend = "orm",
    cluster = cluster,
    type = type,
    cadjust = cadjust,
    adjustment_factor = assembled$adjustment_factor,
    n_rows = nrow(fit$x),
    n_represented_rows = info$n_likelihood_terms,
    n_parameters = length(coefficients),
    n_clusters = n_clusters,
    bread_convention = "inverse observed information of the path likelihood",
    weighted = !is.null(fit$weights) && any(fit$weights != 1),
    covariance_identity = assembled$covariance
  )
  fit
}
