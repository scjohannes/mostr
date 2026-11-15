# Fitting when some states are not observed: the maximum likelihood fit for
# `engine = "orm"` and the pieces shared with `engine = "blrm"`.
#
# The model structure (design matrix columns, spline knots, factor levels,
# state labels) and the starting values come from the engine's own fit to the
# set-up rows: every visit with an exactly observed state, with an unknown
# previous state filled in by the patient's last known state
# (`markov_setup_rows()`). mostr then maximises the likelihood of the fully
# observed transitions and the stretches of unobserved visits
# (`markov_path_plan()`), in which no filled-in state appears, and puts its
# estimates into the engine's fit, restricted to the fully observed
# transitions, so that everything downstream reads an ordinary `orm` object.

# Arguments of `mostm()` that the path likelihood has not been worked out for.
mostm_check_path_fit_args <- function(engine, dots, p2_var) {
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

# The formula and data for an engine fit to exactly observed states. With
# `Ocens()` the engine must not see the ranges, so the response becomes a
# column that holds the exactly observed states and `NA` elsewhere.
markov_observed_transition_setup <- function(setup, ranges) {
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
  setup
}

# The rows the engine's set-up fit uses when there are stretches, and the
# previous states it reads.
#
# The engine's fit gives the model structure (design columns, spline knots,
# factor levels, state labels) and the starting values. Fitted to the fully
# observed transitions alone, it fails when a state, or a previous state,
# occurs only next to unobserved visits: with deaths always recorded and
# other states often missing, every death can follow an unobserved visit, so
# no fully observed transition ends in death (1 of 500 simulated trials with
# 25% of the states missing, not at random), and the factor level of a
# previous state that is always followed by an unobserved visit has no rows,
# so `rms::orm()` finds a singular information matrix. Few rows also make
# that fit fail to converge more often.
#
# So the set-up fit uses every row of `rows` (the rows `subset` keeps) whose
# state is observed exactly and whose previous state is recorded or can be
# filled in with the patient's last exactly observed state. That is the last
# recorded previous state of the patient up to this row (the rows of a
# patient are consecutive visits, which `markov_path_layout()` checks): a
# visit after an unobserved one has no recorded previous state, and the
# visit before it carries the last known state as its own previous state, or
# carries it forward in turn. Rows before the patient's first known state
# get none and are not used; they are not in the likelihood either. The
# filled-in states are a device for the set-up only: they shift the starting
# values (a transition over several visits is fitted as one), which the
# maximisation then corrects, and they never enter the likelihood
# (`recorded` says which fitted rows are fully observed transitions).
#
# Returns `previous` (the previous-state column, filled in), `recorded`
# (whether the previous state is recorded) and `use` (whether the row is
# fitted), each with one entry per row of `data`.
markov_setup_rows <- function(data, ranges, rows, id_var, time_var, p_var) {
  n <- nrow(data)
  previous <- data[[p_var]]
  recorded <- !is.na(previous)
  if (length(rows)) {
    # Visits in time order within patient, as in `markov_path_layout()`.
    time <- data[[time_var]][rows]
    id <- as.character(data[[id_var]][rows])
    ord <- if (is.numeric(time) || is.factor(time)) {
      order(match(id, unique(id)), time)
    } else {
      order(match(id, unique(id)))
    }
    sorted <- rows[ord]
    m <- length(sorted)
    s_id <- id[ord]
    patient_start <- c(TRUE, s_id[-1L] != s_id[-m])
    # Position of the last row with a recorded previous state up to each
    # row, and of the patient's first row: the previous state is carried
    # forward only within the patient.
    last_known <- cummax(ifelse(recorded[sorted], seq_len(m), 0L))
    first_of_patient <- cummax(ifelse(patient_start, seq_len(m), 0L))
    fill <- !recorded[sorted] & last_known >= first_of_patient & last_known > 0L
    previous[sorted[fill]] <- previous[sorted[last_known[fill]]]
  }
  in_rows <- logical(n)
  in_rows[rows] <- TRUE
  exact <- !is.na(ranges$lo) & ranges$lo == ranges$hi
  list(
    previous = previous,
    recorded = recorded,
    use = in_rows & exact & !is.na(previous)
  )
}

# The setup of the engine's set-up fit (`markov_setup_rows()`), and for each
# row of `data` whether its previous state is recorded, whether the set-up
# fit uses it and its previous state as the set-up fit reads it.
# `subset_index` maps the rows of the layout (the rows `subset` keeps) to
# rows of `data`.
markov_path_setup <- function(setup, ranges, path) {
  subset_index <- mostm_subset_row_index(setup$subset_rows, nrow(setup$data))
  rows <- markov_setup_rows(
    setup$data,
    ranges,
    subset_index,
    path$id_var,
    path$time_var,
    path$p_var
  )
  out <- markov_observed_transition_setup(setup, ranges)
  out$data[[path$p_var]] <- rows$previous
  out$subset_rows <- rows$use
  list(
    setup = out,
    recorded = rows$recorded,
    use = rows$use,
    previous = rows$previous,
    subset_index = subset_index
  )
}

# Stop before the engine's fit when a state is never observed exactly in the
# set-up rows: `rms::orm()` would stop first, with an obscure error, on the
# empty factor level. Such a state has no observation in the likelihood
# either, other than as part of a range.
#
# A previous state given as a factor needs set-up rows too: a level without
# rows makes `rms::orm()` stop with a singular information matrix, and its
# effect would be determined by ranges or stretches alone. That happens when
# a state is observed only at a patient's last recorded visit and the
# likelihood still needs it as a previous state (a range after it, or a
# state inside a stretch).
mostm_check_setup_states <- function(ranges, set_up, layout, p_var) {
  use <- set_up$use
  present <- ranges$labels[sort(unique(ranges$lo[use]))]
  markov_check_observed_states(present, ranges$labels)
  in_formula <- p_var %in% all.vars(set_up$setup$formula[[3L]])
  if (!in_formula || (!is.factor(set_up$previous) && !is.character(set_up$previous))) {
    return(invisible(NULL))
  }
  later <- which(!layout$is_first)
  candidates <- sequence(
    layout$step_hi[later - 1L] - layout$step_lo[later - 1L] + 1L,
    from = layout$step_lo[later - 1L]
  )
  needed <- union(
    as.character(set_up$previous[
      set_up$subset_index[layout$step_data_rows[layout$is_first]]
    ]),
    layout$labels[setdiff(candidates, layout$absorb)]
  )
  absent <- setdiff(needed, as.character(set_up$previous[use]))
  if (length(absent)) {
    mostr_abort(
      "The state",
      if (length(absent) > 1L) "s " else " ",
      paste(absent, collapse = ", "),
      if (length(absent) > 1L) " are" else " is",
      " never the previous state (`",
      p_var,
      "`) of a visit with an exactly observed state, although the model ",
      "needs ",
      if (length(absent) > 1L) "their effects" else "its effect",
      " for visits with unobserved states, so the model can't be set up. ",
      "Merge rarely observed states."
    )
  }
  invisible(NULL)
}

markov_check_observed_states <- function(fitted_labels, labels) {
  if (identical(as.character(fitted_labels), labels)) {
    return(invisible(NULL))
  }
  absent <- setdiff(labels, fitted_labels)
  mostr_abort(
    "The state",
    if (length(absent) > 1L) "s " else " ",
    paste(absent, collapse = ", "),
    if (length(absent) > 1L) " are" else " is",
    " never observed exactly at a visit after the patient's first known ",
    "state, so the model can't be set up. Every state must be observed ",
    "exactly at least once; a range given with `Ocens()` does not count. ",
    "Merge rarely observed states, or check the bounds given to `Ocens()`."
  )
}

# The rows of the engine's set-up fit that are fully observed transitions.
# `fit_data` holds the fitted rows of the set-up data (as returned by
# `mostm_fit_orm()`); their row names identify the rows of `data`. A fitted
# row whose previous state was filled in ends a stretch and enters the
# likelihood there, so it must not be counted as a transition as well.
markov_setup_recorded_rows <- function(fit_data, data, recorded, n_fitted) {
  row <- match(rownames(fit_data), rownames(data))
  if (anyNA(row) || length(row) != n_fitted) {
    mostr_abort(
      "Internal error: the rows of the set-up fit for unobserved states ",
      "can't be matched to the rows of `data`."
    )
  }
  recorded[row]
}

# The engine's fit restricted to the fully observed transitions `keep`: the
# fields that hold one entry per fitted row, and the counts rms computes from
# them. The design, the states and the coefficients stay; downstream code
# reads `x`, `y` and `weights` as the fully observed transitions.
markov_orm_keep_rows <- function(fit, keep) {
  n <- length(keep)
  fit$x <- fit$x[keep, , drop = FALSE]
  fit$y <- markov_keep_response_rows(fit$y, keep)
  if (length(fit$weights) == n) {
    fit$weights <- fit$weights[keep]
  }
  if (length(fit$linear.predictors) == n) {
    fit$linear.predictors <- fit$linear.predictors[keep]
  }
  states <- as.character(fit$yunique)
  observed <- markov_blrm_observed_states(list(y = fit$y, ylevels = states))
  fit$freq <- stats::setNames(
    as.integer(table(factor(observed, levels = states))),
    states
  )
  if (!is.null(fit$stats)) {
    fit$stats[["Obs"]] <- sum(keep)
    fit$stats[["Distinct Y"]] <- length(unique(observed))
    # rms computes these from the fitted rows with formulas of its own; they
    # are not recomputed for the fully observed transitions.
    fit$stats[intersect(names(fit$stats), c("ESS", "Median Y"))] <- NA_real_
  }
  fit
}

# Rows `keep` of a response stored as a vector (orm) or as a two-column
# matrix (`Ocens` in rmsb), keeping its class and attributes.
markov_keep_response_rows <- function(y, keep) {
  if (is.null(dim(y))) {
    return(y[keep])
  }
  attributes_kept <- attributes(y)
  attributes_kept[c("dim", "dimnames")] <- NULL
  out <- unclass(y)[keep, , drop = FALSE]
  for (name in names(attributes_kept)) {
    attr(out, name) <- attributes_kept[[name]]
  }
  out
}

# The engine's set-up fit with `fit()`. When rms's default optimiser fails on
# it, the fit is repeated once with `opt_method = "LM", maxit = 100`, which
# the error message recommends to users: the set-up fit only gives the
# structure and the starting values, so the user should not have to choose
# its optimiser. In a simulation with 25% of the states missing not at
# random, 22 of 500 fits to the fully observed transitions needed this.
markov_setup_fit_with_retry <- function(fit, dots) {
  tryCatch(
    fit(list()),
    mostr_error_orm_failed = function(e) {
      if (any(c("opt_method", "maxit") %in% names(dots))) {
        mostr_abort(e)
      }
      fit(list(opt_method = "LM", maxit = 100))
    }
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
    # The orm fit has no random intercepts (`cluster()` is refused for orm).
    cluster_shift = double(),
    obs_cluster = integer(),
    stretch_cluster = integer(),
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
  set_up <- markov_path_setup(setup, ranges, path)
  mostm_check_setup_states(ranges, set_up, layout, path$p_var)
  dots <- list(...)
  fitted <- markov_setup_fit_with_retry(
    function(extra) {
      do.call(mostm_fit_orm, c(list(set_up$setup), utils::modifyList(dots, extra)))
    },
    dots
  )
  start <- fitted$fit
  if (!identical(start$family, "logistic")) {
    stop_unsupported_markov_link(start$family)
  }
  labels <- ranges$labels
  markov_check_observed_states(start$yunique, labels)
  keep <- markov_setup_recorded_rows(
    fitted$fit_data,
    set_up$setup$data,
    set_up$recorded,
    nrow(start$x)
  )
  fit_data <- fitted$fit_data[keep, , drop = FALSE]

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

  theta <- stats::coef(start)
  fit <- markov_orm_keep_rows(start, keep)
  fit_ids <- as.character(fit_data[[path$id_var]])
  patient_ids <- unique(c(fit_ids, plan$instances$id))
  observed <- markov_path_observed(
    fit$x,
    fit$y,
    labels,
    fit$weights,
    fit_ids,
    patient_ids
  )
  maximum <- markov_path_maximise(unname(theta), observed, plan)

  names_theta <- names(theta)
  information <- -maximum$result$hessian
  dimnames(information) <- list(names_theta, names_theta)
  scores <- markov_path_patient_scores(maximum$result, plan, patient_ids)
  dimnames(scores) <- list(patient_ids, names_theta)

  fit$coefficients <- stats::setNames(maximum$theta, names_theta)
  fit$u <- stats::setNames(maximum$score, names_theta)
  fit$iter <- maximum$iterations
  n_intercepts <- length(labels) - 1L
  fit$linear.predictors <- stats::setNames(
    fit$coefficients[[fit$interceptRef]] +
      as.numeric(observed$x %*% maximum$theta[-seq_len(n_intercepts)]),
    names(fit$linear.predictors)
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
  # These statistics of `rms::orm()` describe its set-up fit (likelihood
  # ratio and score tests against the model without predictors, rank
  # measures of the fitted rows). They are not computed for the likelihood
  # with unobserved states, and stale values would mislead.
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
    fit_data = fit_data,
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
