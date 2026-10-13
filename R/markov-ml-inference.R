# Scores, observed information, and covariance for `orm_markov_ml()`.

# Proportional-odds log-likelihood scores for individual rows.
#
# For a design row x with linear predictor eta, F_j = plogis(alpha_j + eta)
# and f_j = F_j (1 - F_j) for j = 2..K, with f_1 = f_{K+1} = 0. The score of
# log p_b is (1[j = b] f_b - 1[j = b + 1] f_{b+1}) / p_b for alpha_j and
# x (f_b - f_{b+1}) / p_b for the slopes. Returns an n x p matrix in
# coefficient order (thresholds, then slopes).
po_row_scores <- function(x, y, coefficients, n_intercepts) {
  x <- as.matrix(x)
  n <- nrow(x)
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
  cum <- matrix(stats::plogis(outer(eta, alpha, "+")), nrow = n)
  P <- cbind(1, cum) - cbind(cum, 0)
  dens <- cum * (1 - cum)
  fpad <- cbind(0, dens, 0)
  rows <- seq_len(n)
  pb <- P[cbind(rows, y)]
  fb <- fpad[cbind(rows, y)]
  fb1 <- fpad[cbind(rows, y + 1L)]

  intercept <- matrix(0, nrow = n, ncol = n_intercepts)
  upper <- y >= 2L
  intercept[cbind(rows[upper], y[upper] - 1L)] <- fb[upper] / pb[upper]
  lower <- y <= K - 1L
  intercept[cbind(rows[lower], y[lower])] <-
    intercept[cbind(rows[lower], y[lower])] - fb1[lower] / pb[lower]
  slope <- x * ((fb - fb1) / pb)
  out <- cbind(intercept, slope)
  colnames(out) <- names(coefficients)
  out
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

markov_ml_patient_scores <- function(problem, theta, state = NULL) {
  if (is.null(state)) {
    state <- markov_ml_estep(problem, theta)
  }
  p <- length(problem$coef_names)
  U <- matrix(0, nrow = problem$N, ncol = p)
  comp <- problem$complete_part
  if (length(comp$rows)) {
    complete <- rowsum(
      po_row_scores(
        problem$X[comp$design, , drop = FALSE],
        comp$state,
        theta,
        problem$n_int
      ),
      comp$patient,
      reorder = TRUE
    )
    rows <- as.integer(rownames(complete))
    U[rows, ] <- U[rows, ] + complete
  }
  if (problem$chain$n) {
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

markov_ml_total_score <- function(problem, theta) {
  state <- markov_ml_estep(problem, theta)
  comp <- problem$complete_part
  complete <- colSums(
    po_row_scores(
      problem$X[comp$agg_design, , drop = FALSE],
      comp$agg_state,
      theta,
      problem$n_int
    ) *
      comp$agg_count
  )
  out <- complete + colSums(markov_ml_chain_scores(problem, state))
  names(out) <- problem$coef_names
  out
}

# Observed information by central differences of the exact total score.
markov_ml_information <- function(problem, theta, step = 1e-4) {
  p <- length(theta)
  J <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(p)) {
    h <- step * max(1, abs(theta[[k]]))
    up <- theta
    down <- theta
    up[[k]] <- up[[k]] + h
    down[[k]] <- down[[k]] - h
    J[, k] <- -(markov_ml_total_score(problem, up) -
      markov_ml_total_score(problem, down)) /
      (2 * h)
  }
  J <- (J + t(J)) / 2
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
# to missing states in the worst direction.
markov_ml_missing_info_fraction <- function(complete_information, J) {
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
