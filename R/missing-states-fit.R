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
  refuse <- function(what, ...) {
    mostr_abort(
      what,
      " can't be combined with unobserved states (a missing state, a missing ",
      "previous state, or a range of states given with `Ocens()`): ",
      ...,
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

# The fully observed transitions: rows whose state is observed exactly and
# whose previous state is recorded. These, and only these, are fitted as
# ordinary transitions; every other row with a state is a step of a stretch
# (or comes before the patient's first known state and is not used), so a
# row counted here must never also end a stretch.
markov_observed_transition_rows <- function(ranges, previous) {
  !is.na(ranges$lo) & ranges$lo == ranges$hi & !is.na(previous)
}

# The formula and data for the engine's fit to the fully observed transitions.
# With `Ocens()` the engine must not see the ranges, so the response becomes a
# column that holds the exactly observed states and `NA` elsewhere.
#
# With `p_var`, the fit is also restricted (through `subset`) to rows whose
# previous state is recorded. The engine drops a row with a missing `p_var`
# by itself only when `p_var` is in the formula; with a formula without it
# (or with the previous state only in a column derived from it), a row with
# an observed state after an unobserved visit was fitted by the engine and
# was also the last step of a stretch, so it entered the likelihood twice.
# Without stretches the engine's own fit is used and `p_var` is not given:
# such rows are then fitted as before.
markov_observed_transition_setup <- function(setup, ranges, p_var = NULL) {
  if (ranges$ocens) {
    exact <- !is.na(ranges$lo) & ranges$lo == ranges$hi
    state <- rep(NA_character_, nrow(setup$data))
    state[exact] <- ranges$labels[ranges$lo[exact]]
    # Only states that are observed exactly become levels: rms stops with an
    # obscure error on a level without observations.
    setup$data$.markov_state <- factor(
      state,
      levels = ranges$labels[ranges$labels %in% state],
      ordered = TRUE
    )
    setup$formula[[2L]] <- quote(.markov_state)
  }
  if (!is.null(p_var)) {
    n <- nrow(setup$data)
    keep <- logical(n)
    keep[mostm_subset_row_index(setup$subset_rows, n)] <- TRUE
    setup$subset_rows <- keep &
      markov_observed_transition_rows(ranges, setup$data[[p_var]])
  }
  setup
}

# Stop before the engine's fit when a state never occurs at a visit whose
# state and previous state are both observed: `rms::orm()` would stop first,
# with "missing value where TRUE/FALSE needed", on the empty factor level.
mostm_check_observed_states <- function(setup, ranges, path) {
  rows <- mostm_subset_row_index(setup$subset_rows, nrow(setup$data))
  observed <- markov_observed_transition_rows(
    list(lo = ranges$lo[rows], hi = ranges$hi[rows]),
    setup$data[[path$p_var]][rows]
  )
  present <- ranges$labels[sort(unique(ranges$lo[rows][observed]))]
  markov_check_observed_states(present, ranges$labels)
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
  beta <- theta[-seq_len(n_intercepts)]
  cpp_path_loglik(
    alpha = theta[seq_len(n_intercepts)],
    beta = beta,
    x_obs = observed$x,
    y_obs = observed$y,
    w_obs = observed$weight,
    group_obs = observed$group,
    n_groups = as.integer(observed$n_groups),
    x_gap = plan$x_gap,
    x_prev = plan$x_prev %||% matrix(0, 0L, length(beta)),
    absorb = plan$absorb,
    stretch_first_row = plan$stretch_first_row,
    stretch_first_prev = plan$stretch_first_prev %||% integer(),
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
    # A patient of weight 0 adds nothing, also when the stretch is impossible
    # and its score NaN.
    stretch_scores[plan$instances$weight == 0, ] <- 0
    patient <- match(plan$instances$id, patient_ids)
    summed <- rowsum(stretch_scores, patient)
    rows <- as.integer(rownames(summed))
    scores[rows, ] <- scores[rows, , drop = FALSE] + summed
  }
  scores
}

# The total score: the observed transitions plus each stretch's score times
# its weight. A stretch of weight 0 adds exactly nothing, also when it is
# impossible at `theta` and its score is NaN (NaN * 0 would be NaN).
markov_path_total_score <- function(result, plan) {
  weight <- plan$stretch_weight
  used <- weight != 0
  stretch_scores <- result$stretch_scores
  if (!all(used)) {
    stretch_scores <- stretch_scores[used, , drop = FALSE]
    weight <- weight[used]
  }
  colSums(result$obs_group_scores) +
    as.numeric(crossprod(stretch_scores, weight))
}

# Maximise the likelihood with unobserved states from `theta`, the estimates
# from the fully observed transitions; see `markov_newton_maximise()`.
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
  current <- evaluate(theta, 2L)
  if (!is.finite(current$loglik)) {
    mostr_abort(
      "The likelihood with unobserved states is not finite at the starting ",
      "values, the estimates from the fully observed transitions."
    )
  }
  maximum <- markov_newton_maximise(
    theta,
    evaluate,
    function(result) markov_path_total_score(result, plan),
    current = current,
    maxit = maxit,
    tolerance = tolerance
  )
  if (!maximum$converged) {
    mostr_abort(
      "The fit with unobserved states did not converge. Some coefficients ",
      "may have no finite estimate, for example because a transition never ",
      "occurs in the data: simplify the model or merge rarely observed ",
      "states."
    )
  }
  maximum$converged <- NULL
  maximum
}

# Newton-Raphson with a line search. `evaluate(theta, order)` returns the
# log-likelihood (`loglik`) and, for `order = 2`, the Hessian (`hessian`);
# `score_of()` turns its result into the score. The second derivative is
# exact, so near the maximum the iteration converges quadratically, usually
# in three to five steps from the fit to the observed transitions.
#
# It stops when the Newton decrement (half of score' * step, the gain a full
# step promises) is below `tolerance`, or below 1e-13 of the log-likelihood
# when that is larger, at a point where the log-likelihood is concave. The
# relative part exists because the score cannot be computed more precisely
# than its sum of rounded terms allows: with `rms::rcs(time, 4)`, 2,000
# patients and a log-likelihood near -2e4, the largest score entries stayed
# near 2e-3 and the decrement near 1.6e-10 (8e-15 of the log-likelihood)
# however often the step was repeated, so an absolute 1e-10 was never met and
# the fit stopped after 50 iterations as "not converged". Up to a
# log-likelihood of -1000 the absolute 1e-10 applies, as before; at it the
# estimates are settled to far below their standard errors.
#
# Returns the maximum, the evaluation there (with the exact Hessian), the
# score, the number of iterations and whether it converged.
markov_newton_maximise <- function(
  theta,
  evaluate,
  score_of,
  current = evaluate(theta, 2L),
  maxit = 50L,
  tolerance = 1e-10
) {
  converged <- FALSE
  steps <- 0L
  for (iteration in seq_len(maxit)) {
    score <- score_of(current)
    direction <- markov_ascent_direction(-current$hessian, score)
    if (is.null(direction)) {
      break
    }
    step <- direction$step
    gain <- sum(score * step) / 2
    scale <- max(1, abs(current$loglik))
    # Differences of log-likelihoods within this are rounding error: the
    # log-likelihood is a sum of many rounded terms.
    rounding <- 1e-13 * scale
    # Only a point where the log-likelihood is concave can be the maximum.
    concave <- !direction$modified
    if (concave && gain < max(tolerance, rounding)) {
      converged <- TRUE
      break
    }
    # The full step is evaluated with its derivatives straight away, because
    # it is nearly always accepted. A Newton step can overshoot far from the
    # maximum; it is halved until the log-likelihood rises by at least 1e-4 of
    # what the slope promises for that step (Armijo's rule), less the
    # rounding error of the log-likelihood.
    accepted <- NULL
    for (halving in 0:30) {
      candidate <- theta + step / 2^halving
      trial <- evaluate(candidate, if (halving == 0L) 2L else 0L)
      needed <- 1e-4 * 2 * gain / 2^halving - rounding
      if (is.finite(trial$loglik) && trial$loglik >= current$loglik + needed) {
        accepted <- trial
        break
      }
    }
    if (is.null(accepted)) {
      break
    }
    change <- accepted$loglik - current$loglik
    steps <- steps + 1L
    theta <- candidate
    current <- if (halving == 0L) accepted else evaluate(theta, 2L)
    # A step that leaves the log-likelihood unchanged within rounding, while
    # the decrement is small (below 1e-10 of the log-likelihood, about 1e-5
    # standard errors), only moves through the noise of the score: the
    # estimates are settled, and further steps would repeat up to `maxit`.
    if (concave && abs(change) <= rounding && gain < 1e-10 * scale) {
      converged <- TRUE
      break
    }
  }
  list(
    theta = theta,
    result = current,
    score = score_of(current),
    iterations = steps,
    converged = converged
  )
}

# A step that raises the log-likelihood, from the information matrix (the
# negative Hessian) and the score. The log-likelihood with unobserved states
# need not be concave away from its maximum. Where the information has a
# negative eigenvalue, the Newton step solve(information, score) can point
# downhill, and then no step length raises the log-likelihood: at the
# starting values of 30 patients x 20 visits, 8 states and 70% of the states
# missing, the information had eigenvalues from -12.7 to 3.0e3 and the
# Newton step lowered the log-likelihood however much it was shortened.
#
# Where the information is positive definite (its Cholesky factor exists)
# the step is the Newton step. Elsewhere each eigenvalue is replaced by its
# absolute value, at least 1e-8 of the largest: the result is positive
# definite, so the step is an ascent direction (score' step > 0), and along
# a direction of negative curvature it moves as far as the curvature
# suggests. Adding a multiple of the identity until the matrix is positive
# definite, the other common choice, overshoots the negative eigenvalue by
# up to the factor it grows by and then crawls: on f = -(x^2 - 1)^2 from
# x = 0.1 it needed more than 40 steps. Returns the step and whether the
# information was modified (`modified`), or NULL when no step can be formed:
# a score or information that is not finite, or an information that is 0.
markov_ascent_direction <- function(information, score) {
  if (!all(is.finite(score)) || !all(is.finite(information))) {
    return(NULL)
  }
  factor <- tryCatch(chol(information), error = function(e) NULL)
  if (!is.null(factor)) {
    step <- backsolve(factor, backsolve(factor, score, transpose = TRUE))
    if (all(is.finite(step))) {
      return(list(step = as.numeric(step), modified = FALSE))
    }
  }
  decomposition <- eigen(information, symmetric = TRUE)
  values <- abs(decomposition$values)
  smallest <- 1e-8 * max(values)
  if (!(smallest > 0)) {
    return(NULL)
  }
  values <- pmax(values, smallest)
  vectors <- decomposition$vectors
  step <- as.numeric(vectors %*% (crossprod(vectors, score) / values))
  list(step = step, modified = TRUE)
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

# The layout of the stretches of the rows `subset` keeps
# (`markov_path_layout()`): whether there are any, and which rows they hold.
mostm_path_layout <- function(setup, ranges, path) {
  rows <- mostm_subset_row_index(setup$subset_rows, nrow(setup$data))
  data <- setup$data[rows, , drop = FALSE]
  markov_path_layout(
    data = data,
    ranges = list(
      labels = ranges$labels,
      lo = ranges$lo[rows],
      hi = ranges$hi[rows]
    ),
    absorb = match(as_state_labels(path$absorb), ranges$labels),
    id_var = path$id_var,
    time_var = path$time_var,
    p_var = path$p_var,
    weights = if (!is.null(setup$prior_weight_column)) {
      data[[setup$prior_weight_column]]
    }
  )
}

mostm_fit_orm_paths <- function(setup, ranges, path, ...) {
  layout <- mostm_path_layout(setup, ranges, path)
  # Unobserved states that say nothing about the model leave no stretch: a
  # patient's first row without a previous state, or unobserved visits after
  # the patient's last recorded one. The engine's fit to the observed
  # transitions is then the maximum likelihood fit.
  if (!layout$n_steps_total) {
    return(mostm_fit_orm(markov_observed_transition_setup(setup, ranges), ...))
  }
  mostm_check_observed_states(setup, ranges, path)
  fitted <- mostm_fit_orm(
    markov_observed_transition_setup(setup, ranges, path$p_var),
    ...
  )
  start <- fitted$fit
  if (!identical(start$family, "logistic")) {
    stop_unsupported_markov_link(start$family)
  }
  labels <- ranges$labels
  markov_check_observed_states(start$yunique, labels)

  plan <- markov_path_plan_design(
    layout,
    function(rows) orm_model_matrix(start, rows)
  )
  if (is.null(plan$x_gap)) {
    plan$x_gap <- start$x[0L, , drop = FALSE]
  }
  plan$x_gap <- unname(plan$x_gap)
  storage.mode(plan$x_gap) <- "double"
  if (!is.null(plan$x_prev)) {
    storage.mode(plan$x_prev) <- "double"
  }

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
