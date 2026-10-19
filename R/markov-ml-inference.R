# Scores, observed information, and covariance for `orm_markov_ml()`.

# Proportional-odds terms of the observed state of each row.
#
# For a design row x with linear predictor eta, F_j = plogis(alpha_j + eta)
# is P(Y >= j) for j = 2..K, with F_1 = 1 and F_{K+1} = 0, and
# f_j = F_j (1 - F_j). The probability of state b, p_b = F_b - F_{b+1},
# depends only on the two thresholds around b: alpha_b (coefficient b - 1,
# present when b >= 2, `upper`) and alpha_{b+1} (coefficient b, present when
# b <= K - 1, `lower`). Only these two are computed, so the cost does not grow
# with the number of states.
po_outcome_terms <- function(x, y, coefficients, n_intercepts) {
  x <- as.matrix(x)
  K <- n_intercepts + 1L
  alpha <- coefficients[seq_len(n_intercepts)]
  beta <- coefficients[-seq_len(n_intercepts)]
  if (length(beta) != ncol(x)) {
    stop("`x` must have one column per slope coefficient.")
  }
  y <- as.integer(y)
  if (anyNA(y) || any(y < 1L | y > K)) {
    stop("`y` must contain state indices between 1 and ", K, ".")
  }
  eta <- drop(x %*% beta)
  upper <- y >= 2L
  lower <- y <= K - 1L
  F_upper <- rep(1, length(y))
  F_upper[upper] <- stats::plogis(eta[upper] + alpha[y[upper] - 1L])
  F_lower <- rep(0, length(y))
  F_lower[lower] <- stats::plogis(eta[lower] + alpha[y[lower]])
  list(
    x = x,
    y = y,
    upper = upper,
    lower = lower,
    F_upper = F_upper,
    F_lower = F_lower,
    f_upper = F_upper * (1 - F_upper),
    f_lower = F_lower * (1 - F_lower),
    prob = F_upper - F_lower
  )
}

# Sums of `values` (a vector or a matrix with one row per element of `index`)
# by `index` in 1..size, as a size-row matrix.
po_sum_by <- function(values, index, size) {
  values <- as.matrix(values)
  out <- matrix(0, nrow = size, ncol = ncol(values))
  if (length(index)) {
    out[sort(unique(index)), ] <- rowsum(values, index, reorder = TRUE)
  }
  out
}

# Proportional-odds log-likelihood scores for individual rows.
#
# The score of log p_b is (1[j = b] f_b - 1[j = b + 1] f_{b+1}) / p_b for
# alpha_j and x (f_b - f_{b+1}) / p_b for the slopes. Returns an n x p matrix
# in coefficient order (thresholds, then slopes).
po_row_scores <- function(x, y, coefficients, n_intercepts) {
  terms <- po_outcome_terms(x, y, coefficients, n_intercepts)
  y <- terms$y
  n <- length(y)
  rows <- seq_len(n)
  upper <- terms$upper
  lower <- terms$lower
  intercept <- matrix(0, nrow = n, ncol = n_intercepts)
  intercept[cbind(rows[upper], y[upper] - 1L)] <-
    terms$f_upper[upper] / terms$prob[upper]
  intercept[cbind(rows[lower], y[lower])] <-
    intercept[cbind(rows[lower], y[lower])] -
    terms$f_lower[lower] / terms$prob[lower]
  slope <- terms$x * ((terms$f_upper - terms$f_lower) / terms$prob)
  out <- cbind(intercept, slope)
  colnames(out) <- names(coefficients)
  out
}

# Weighted sums of the row scores within groups (for example patients), one
# row per group 1..n_groups. Each row's threshold scores have at most two
# nonzero entries, so no row-by-coefficient matrix is built.
po_group_scores <- function(
  x,
  y,
  coefficients,
  n_intercepts,
  group = rep(1L, NROW(x)),
  n_groups = 1L,
  weights = 1
) {
  terms <- po_outcome_terms(x, y, coefficients, n_intercepts)
  y <- terms$y
  upper <- terms$upper
  lower <- terms$lower
  w <- rep_len(weights, length(y)) / terms$prob
  # Threshold j of group g sits at position (j - 1) * n_groups + g of the
  # n_groups x n_intercepts block.
  position <- c(
    (y[upper] - 2) * n_groups + group[upper],
    (y[lower] - 1) * n_groups + group[lower]
  )
  intercept <- po_sum_by(
    c(w[upper] * terms$f_upper[upper], -w[lower] * terms$f_lower[lower]),
    position,
    n_groups * n_intercepts
  )
  slope <- po_sum_by(
    terms$x * (w * (terms$f_upper - terms$f_lower)),
    group,
    n_groups
  )
  out <- cbind(matrix(intercept, nrow = n_groups), slope)
  colnames(out) <- names(coefficients)
  out
}

# Exact observed information (negative Hessian) of the weighted
# proportional-odds log-likelihood sum(weights * log p_b).
#
# Write u and v for the linear predictors of F_b and F_{b+1}, and
# f'(u) = f(u) (1 - 2 F(u)). The second derivatives of log p_b are
# f'(u) / p_b - f(u)^2 / p_b^2 in u, -f'(v) / p_b - f(v)^2 / p_b^2 in v, and
# f(u) f(v) / p_b^2 in u and v. A row touches only the thresholds of u and v,
# so the threshold block is tridiagonal; the slope entries multiply these by
# x (once for a threshold and slope, x x' for two slopes). Returns a dense
# p x p matrix in coefficient order.
po_information <- function(x, y, coefficients, n_intercepts, weights = 1) {
  terms <- po_outcome_terms(x, y, coefficients, n_intercepts)
  y <- terms$y
  upper <- terms$upper
  lower <- terms$lower
  w <- rep_len(weights, length(y))
  prob <- terms$prob
  fu <- terms$f_upper
  fl <- terms$f_lower
  i_uu <- w * (fu^2 / prob^2 - fu * (1 - 2 * terms$F_upper) / prob)
  i_vv <- w * (fl^2 / prob^2 + fl * (1 - 2 * terms$F_lower) / prob)
  i_uv <- -w * fu * fl / prob^2

  x <- terms$x
  p <- length(coefficients)
  slopes <- n_intercepts + seq_len(ncol(x))
  thresholds <- seq_len(n_intercepts)
  J <- matrix(0, nrow = p, ncol = p)
  J[cbind(thresholds, thresholds)] <- drop(
    po_sum_by(i_uu[upper], y[upper] - 1L, n_intercepts) +
      po_sum_by(i_vv[lower], y[lower], n_intercepts)
  )
  if (n_intercepts > 1L) {
    # Rows with both thresholds link coefficients b - 1 and b.
    both <- upper & lower
    off <- drop(po_sum_by(i_uv[both], y[both] - 1L, n_intercepts - 1L))
    j <- seq_len(n_intercepts - 1L)
    J[cbind(j, j + 1L)] <- off
    J[cbind(j + 1L, j)] <- off
  }
  if (ncol(x)) {
    cross <- po_sum_by(
      x[upper, , drop = FALSE] * (i_uu + i_uv)[upper],
      y[upper] - 1L,
      n_intercepts
    ) +
      po_sum_by(
        x[lower, , drop = FALSE] * (i_vv + i_uv)[lower],
        y[lower],
        n_intercepts
      )
    J[thresholds, slopes] <- cross
    J[slopes, thresholds] <- t(cross)
    J[slopes, slopes] <- crossprod(x, x * (i_uu + 2 * i_uv + i_vv))
  }
  dimnames(J) <- list(names(coefficients), names(coefficients))
  J
}

# Exact per-patient observed-data scores by Fisher's identity: the
# complete-data score averaged over the posterior transition weights. With
# xi(a, b) / p_b = FP[a] Q[b], the threshold score for alpha_j is
# FP[a] f_j (Q_j - Q_{j-1}) and the linear-predictor score is its sum over j.
#
# Rows whose own and previous states are observed have point-mass posteriors,
# so their scores are ordinary proportional-odds row scores.
markov_ml_chain_scores <- function(problem, state) {
  chain <- problem$chain
  K <- problem$K
  p <- length(problem$coef_names)
  if (!chain$n) {
    return(matrix(0, nrow = 0L, ncol = p))
  }
  dQ <- state$Q[, -1L, drop = FALSE] - state$Q[, -K, drop = FALSE]
  intercept <- matrix(0, nrow = chain$n, ncol = problem$n_int)
  slope <- matrix(0, nrow = chain$n, ncol = ncol(problem$X))
  for (ai in seq_len(problem$n_origin)) {
    a <- problem$origins[[ai]]
    d <- chain$DI[, ai]
    D <- state$fwd$FP[, a] * state$probs$dens[d, , drop = FALSE] * dQ
    intercept <- intercept + D
    slope <- slope + problem$X[d, , drop = FALSE] * rowSums(D)
  }
  cbind(intercept, slope)
}

# `state` is the forward-backward pass at `theta`; it is needed only when some
# rows need the recursion.
markov_ml_patient_scores <- function(problem, theta, state = NULL) {
  p <- length(problem$coef_names)
  U <- matrix(0, nrow = problem$N, ncol = p)
  comp <- problem$complete_part
  if (length(comp$rows)) {
    U <- U +
      po_group_scores(
        problem$X[comp$design, , drop = FALSE],
        comp$state,
        theta,
        problem$n_int,
        group = comp$patient,
        n_groups = problem$N
      )
  }
  if (problem$chain$n) {
    if (is.null(state)) {
      state <- markov_ml_estep(problem, theta)
    }
    chain <- rowsum(
      markov_ml_chain_scores(problem, state),
      problem$chain$patient,
      reorder = TRUE
    )
    rows <- as.integer(rownames(chain))
    U[rows, ] <- U[rows, ] + chain
  }
  dimnames(U) <- list(problem$patient_ids, problem$coef_names)
  U
}

# Score and observed information of the fully observed transitions, exact
# from the proportional-odds formulas on their aggregated rows.
markov_ml_complete_score <- function(problem, theta) {
  comp <- problem$complete_part
  drop(po_group_scores(
    problem$X[comp$agg_design, , drop = FALSE],
    comp$agg_state,
    theta,
    problem$n_int,
    weights = comp$agg_count
  ))
}

markov_ml_complete_information <- function(problem, theta) {
  comp <- problem$complete_part
  po_information(
    problem$X[comp$agg_design, , drop = FALSE],
    comp$agg_state,
    theta,
    problem$n_int,
    weights = comp$agg_count
  )
}

# Log-likelihood of the fully observed transitions.
markov_ml_complete_loglik <- function(problem, theta) {
  comp <- problem$complete_part
  prob <- po_outcome_terms(
    problem$X[comp$agg_design, , drop = FALSE],
    comp$agg_state,
    theta,
    problem$n_int
  )$prob
  if (any(!(prob > 0))) {
    p_rows <- po_outcome_terms(
      problem$X[comp$design, , drop = FALSE],
      comp$state,
      theta,
      problem$n_int
    )$prob
    markov_ml_stop_impossible(problem$ids[comp$rows[!(p_rows > 0)]])
  }
  sum(comp$agg_count * log(prob))
}

# The rows that need the forward-backward recursion as a problem of their own:
# only the design rows they reference, with `chain$DI` pointing into these, and
# no fully observed transitions. Its score is the recursion rows' part of the
# total score, at a cost that does not grow with the other rows.
markov_ml_chain_problem <- function(problem) {
  chain <- problem$chain
  used <- sort(unique(as.vector(chain$DI)))
  chain$DI[] <- match(chain$DI, used)
  problem$X <- problem$X[used, , drop = FALSE]
  problem$chain <- chain
  problem$complete_part <- list(
    rows = integer(),
    origin = integer(),
    design = integer(),
    state = integer(),
    patient = integer(),
    W = NULL,
    agg_design = integer(),
    agg_state = integer(),
    agg_count = numeric()
  )
  problem$design_data <- NULL
  problem
}

markov_ml_chain_score <- function(chain_problem, theta) {
  state <- markov_ml_estep(chain_problem, theta)
  colSums(markov_ml_chain_scores(chain_problem, state))
}

markov_ml_total_score <- function(problem, theta) {
  out <- markov_ml_complete_score(problem, theta)
  if (problem$chain$n) {
    out <- out +
      markov_ml_chain_score(markov_ml_chain_problem(problem), theta)
  }
  names(out) <- problem$coef_names
  out
}

# Observed information: exact for the fully observed transitions, plus
# central differences of the score of the rows that need the recursion.
markov_ml_information <- function(problem, theta, step = 1e-4) {
  J <- markov_ml_complete_information(problem, theta)
  if (problem$chain$n) {
    chain_problem <- markov_ml_chain_problem(problem)
    for (k in seq_along(theta)) {
      h <- step * max(1, abs(theta[[k]]))
      up <- theta
      down <- theta
      up[[k]] <- up[[k]] + h
      down[[k]] <- down[[k]] - h
      J[, k] <- J[, k] -
        (markov_ml_chain_score(chain_problem, up) -
          markov_ml_chain_score(chain_problem, down)) /
          (2 * h)
    }
    J <- (J + t(J)) / 2
  }
  dimnames(J) <- list(problem$coef_names, problem$coef_names)
  J
}

# Inverse of the observed information, or NULL when it is not positive
# definite.
markov_ml_invert_information <- function(J) {
  R <- tryCatch(chol(J), error = function(e) NULL)
  if (is.null(R)) {
    return(NULL)
  }
  V <- chol2inv(R)
  dimnames(V) <- dimnames(J)
  (V + t(V)) / 2
}

# Patient-cluster sandwich with the conventions of `robcov_orm()`. Only
# patients with at least one informative row are counted as clusters, and HC1
# counts informative rows as observations.
markov_ml_sandwich <- function(
  patient_scores,
  bread,
  n_rows,
  contributing,
  type = "HC0",
  cadjust = TRUE
) {
  scores <- patient_scores[contributing, , drop = FALSE]
  assembled <- markov_assemble_sandwich(
    sandwich_scores = scores,
    bread = bread,
    n = n_rows,
    p = ncol(patient_scores),
    n_clusters = nrow(scores),
    type = type,
    cadjust = cadjust
  )
  dimnames(assembled$covariance) <- dimnames(bread)
  dimnames(assembled$meat) <- dimnames(bread)
  assembled
}

# Largest eigenvalue of I_c^{-1} (I_c - J): the fraction of information lost
# to missing states in the worst direction. When nothing is missing, both
# matrices are the same object and nothing is lost.
markov_ml_missing_info_fraction <- function(complete_information, J) {
  if (identical(complete_information, J)) {
    return(0)
  }
  R <- tryCatch(chol(complete_information), error = function(e) NULL)
  if (is.null(R)) {
    return(NA_real_)
  }
  Rinv <- backsolve(R, diag(nrow(R)))
  M <- t(Rinv) %*% (complete_information - J) %*% Rinv
  max(eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values)
}

# ---------------------------------------------------------------------------
# Patient influence components
# ---------------------------------------------------------------------------

markov_is_ml_fit <- function(model) {
  identical(markov_model_fit_wrapper(model), "orm_markov_ml")
}

markov_ml_reject_cluster <- function(cluster) {
  if (!is.null(cluster)) {
    stop(
      "`cluster` cannot be supplied for an `orm_markov_ml()` fit: its ",
      "patient clusters are the chains defined by `id_var` when fitting.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

markov_ml_stored <- function(model) {
  ml <- model$ml
  if (is.null(ml) || is.null(ml$patient_scores) || is.null(ml$bread)) {
    stop(
      "This `orm_markov_ml()` fit has no stored scores and covariance. ",
      "Refit it with `orm_markov_ml()` and the default ",
      "`orm_markov_ml_control(variance = TRUE)`.",
      call. = FALSE
    )
  }
  if (anyNA(ml$bread)) {
    stop(
      "The observed information of this `orm_markov_ml()` fit is not ",
      "positive definite, so no covariance is available.",
      call. = FALSE
    )
  }
  ml
}

# Coefficients, bread, score contributions, and their cluster ids. For
# `orm_markov_ml()` fits these are the observed-data components stored at fit
# time: one score row per patient, already summed over the patient's chain.
# For other orm fits they are the likelihood-row scores and model-based bread.
markov_influence_components <- function(model, cluster = NULL) {
  if (markov_is_ml_fit(model)) {
    markov_ml_reject_cluster(cluster)
    ml <- markov_ml_stored(model)
    return(list(
      coefficients = stats::coef(model),
      bread = ml$bread,
      bread_source = "observed_data_information",
      scores = ml$patient_scores,
      cluster = rownames(ml$patient_scores),
      level = "patient"
    ))
  }

  scores <- compute_scores_orm(model)
  if (!is.null(cluster)) {
    cluster <- align_cluster_orm(model, cluster, nrow(scores))
  }
  bread <- orm_model_bread(model)
  list(
    coefficients = stats::coef(model),
    bread = bread$bread,
    bread_source = bread$source,
    scores = scores,
    cluster = cluster,
    level = "row"
  )
}

# Robust covariance of an `orm_markov_ml()` fit, recomputed from the stored
# patient scores when a different correction is requested.
markov_ml_vcov <- function(model, cluster = NULL, type = NULL, cadjust = NULL) {
  markov_ml_reject_cluster(cluster)
  ml <- markov_ml_stored(model)
  type <- type %||% ml$type
  cadjust <- cadjust %||% ml$cadjust
  if (identical(type, ml$type) && identical(cadjust, ml$cadjust)) {
    return(ml$var_robust)
  }
  markov_ml_sandwich(
    ml$patient_scores,
    ml$bread,
    n_rows = ml$n_informative_rows,
    contributing = ml$contributing,
    type = type,
    cadjust = cadjust
  )$covariance
}
