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

# Exact observed information: the proportional-odds information of the fully
# observed transitions plus that of the segment rows
# (markov_ml_chain_information()). `state` is the forward-backward pass at
# `theta`; it is computed when NULL and some rows need the recursion.
markov_ml_information <- function(problem, theta, state = NULL) {
  J <- markov_ml_complete_information(problem, theta)
  if (problem$chain$n) {
    if (is.null(state)) {
      state <- markov_ml_estep(problem, theta)
    }
    if (!identical(unname(state$theta), unname(theta))) {
      stop("Internal error: the forward-backward pass is not at `theta`.")
    }
    J <- J + markov_ml_chain_information(problem, state)
  }
  dimnames(J) <- list(problem$coef_names, problem$coef_names)
  J
}

# Observed information of the segment rows by Louis' identity.
#
# In a segment, the state Z_0 before the first row is observed and Z_1, ...,
# Z_T are not (the last may be). Row t contributes log P_t(Z_{t-1}, Z_t) to
# the complete-data log-likelihood, with score s_t and negative Hessian h_t,
# and the observed information is
#   E[sum_t h_t | O] - Var(sum_t s_t | O).
# Segments are independent given the observed states, so their terms add up.
# With dP and d2P the first and second derivatives of P_t(a, b) and the
# posterior xi_t(a, b) = FP_t(a) P_t(a, b) Q_t(b) of the forward-backward pass:
# * Within rows, h - s s' = -d2P / P, so E[h_t - s_t s_t' | O] =
#   -sum_ab FP_t(a) Q_t(b) d2P(a, b).
# * Between rows, the rows before t depend on the rows from t on only through
#   Z_{t-1}, so sum_{u < t} E[s_u s_t' | O] = sum_a W_{t-1}(a) G_t(a)', with
#   G_t(a) = sum_b Q_t(b) dP(a, b) and W_{t-1}(a) = FP_t(a) times the expected
#   score of the earlier rows of the segment given Z_{t-1} = a and the
#   observations up to row t - 1. This term enters twice (u < t and u > t).
# * The squared expectation adds E[S | O] E[S | O]' with E[S | O] the
#   segment's score (Fisher's identity).
# W follows the forward pass, with c_t = C_t / s_t (allowed states over the
# scaling constant): W_0 = 0 and
#   W_t(b) = c_t(b) [sum_a P_t(a, b) W_{t-1}(a) + sum_a FP_t(a) dP(a, b)].
# No term divides by a transition probability. W_1 starts from the observed
# state and is sparse (markov_ml_first_row_terms()), so segments of one or two
# rows cost about one recursion pass; only longer segments carry a dense K x p
# matrix per row.
markov_ml_chain_information <- function(problem, state) {
  chain <- problem$chain
  K <- problem$K
  n_int <- problem$n_int
  X <- problem$X
  thresholds <- seq_len(n_int)
  slopes <- n_int + seq_len(ncol(X))
  p <- n_int + ncol(X)
  origins <- problem$origins
  FP <- state$fwd$FP
  dQ <- state$Q[, -1L, drop = FALSE] - state$Q[, -K, drop = FALSE]

  # Within rows. Here F_j = plogis(alpha_j + eta) = P(Y > j) for coefficient
  # j = 1, ..., K - 1 (`cum`), f_j = F_j (1 - F_j) (`dens`), and
  # P(a, b) = F_{b-1} - F_b. The second derivative of P(a, b) in alpha_j is
  # f'_j = f_j (1 - 2 F_j) for b = j + 1 and -f'_j for b = j; the entries for
  # alpha_j and a slope are x times these, and those for two slopes x x'.
  # So sum_b Q(b) d2P(a, b) needs only f'_j (Q(j + 1) - Q(j)), and
  # FP_t(a) (Q_t(j + 1) - Q_t(j)) is summed by design row first.
  R <- matrix(0, nrow = nrow(X), ncol = n_int)
  for (ai in seq_along(origins)) {
    agg <- rowsum(FP[, origins[[ai]]] * dQ, chain$DI[, ai], reorder = TRUE)
    rows <- as.integer(rownames(agg))
    R[rows, ] <- R[rows, ] + agg
  }
  M <- state$probs$dens * (1 - 2 * state$probs$cum) * R
  J <- matrix(0, nrow = p, ncol = p)
  J[cbind(thresholds, thresholds)] <- -colSums(M)
  cross <- -crossprod(M, X)
  J[thresholds, slopes] <- cross
  J[slopes, thresholds] <- t(cross)
  J[slopes, slopes] <- -crossprod(X, X * rowSums(M))

  segment_scores <- rowsum(
    markov_ml_chain_scores(problem, state),
    chain$segment,
    reorder = TRUE
  )
  between <- markov_ml_chain_between(problem, state, dQ)
  J + crossprod(segment_scores) - between - t(between)
}

# W_1 of the first rows `r` of segments (see markov_ml_chain_information()).
# With b0 the observed state before the segment, f the densities f_j of its
# thresholds, x0 its design row, and c = C_1 / s_1, W_1(b) = c_b dP(b0, b) has
# threshold entries c_b f_{b-1} (alpha_{b-1}) and -c_b f_b (alpha_b) and slope
# entries c_b (f_{b-1} - f_b) x0. A segment starting in an absorbing state has
# f = 0: nothing in it depends on the coefficients.
markov_ml_first_row_terms <- function(problem, state, r) {
  chain <- problem$chain
  ai0 <- match(chain$baseline[chain$segment[r]], problem$origins)
  live <- !is.na(ai0)
  d0 <- rep(1L, length(r))
  d0[live] <- chain$DI[cbind(r[live], ai0[live])]
  f <- state$probs$dens[d0, , drop = FALSE]
  f[!live, ] <- 0
  list(
    c = chain$C[r, , drop = FALSE] / state$fwd$s[r],
    f = f,
    x0 = problem$X[d0, , drop = FALSE]
  )
}

# Between-row terms of markov_ml_chain_information(): the sum over segments
# of sum_t sum_a W_{t-1}(a) G_t(a)' (p x p).
markov_ml_chain_between <- function(problem, state, dQ) {
  chain <- problem$chain
  K <- problem$K
  n_int <- problem$n_int
  X <- problem$X
  thresholds <- seq_len(n_int)
  slopes <- n_int + seq_len(ncol(X))
  p <- n_int + ncol(X)
  origins <- problem$origins
  absorbing <- problem$absorbing
  dens <- state$probs$dens
  out <- matrix(0, nrow = p, ncol = p)
  if (length(chain$rows_by_pos) < 2L) {
    return(out)
  }

  # Second rows, all segments at once using the sparse W_1. A segment's rows
  # are consecutive, so its first row precedes its second. Origin a adds
  # c_a f_{a-1} G(a) to the row of alpha_{a-1}, -c_a f_a G(a) to the row of
  # alpha_a, and c_a (f_{a-1} - f_a) G(a) times x0 to the slope rows.
  r2 <- chain$rows_by_pos[[2L]]
  first <- markov_ml_first_row_terms(problem, state, r2 - 1L)
  f <- first$f
  weight <- first$c * (cbind(0, f) - cbind(f, 0))
  dQ2 <- dQ[r2, , drop = FALSE]
  slope_terms <- matrix(0, nrow = length(r2), ncol = p)
  for (ai in seq_along(origins)) {
    a <- origins[[ai]]
    d <- chain$DI[r2, ai]
    G <- dens[d, , drop = FALSE] * dQ2
    G <- cbind(G, X[d, , drop = FALSE] * rowSums(G))
    if (a >= 2L) {
      out[a - 1L, ] <- out[a - 1L, ] +
        drop(crossprod(first$c[, a] * f[, a - 1L], G))
    }
    if (a <= n_int) {
      out[a, ] <- out[a, ] - drop(crossprod(first$c[, a] * f[, a], G))
    }
    slope_terms <- slope_terms + weight[, a] * G
  }
  out[slopes, ] <- out[slopes, ] + crossprod(first$x0, slope_terms)

  # Third and later rows: segment by segment with a dense W.
  length_of <- tabulate(chain$segment, chain$N)
  long <- which(length_of >= 3L)
  if (!length(long)) {
    return(out)
  }
  start <- integer(chain$N)
  start[chain$pats_by_pos[[1L]]] <- chain$rows_by_pos[[1L]]
  first <- markov_ml_first_row_terms(problem, state, start[long])
  above <- cbind(thresholds + 1L, thresholds)
  diagonal <- cbind(thresholds, thresholds)
  for (k in seq_along(long)) {
    c1 <- first$c[k, ]
    f1 <- first$f[k, ]
    W <- matrix(0, nrow = K, ncol = p)
    W[above] <- c1[-1L] * f1
    W[diagonal] <- -c1[-K] * f1
    W[, slopes] <- outer(c1 * (c(0, f1) - c(f1, 0)), first$x0[k, ])
    rows <- start[long[[k]]] + seq_len(length_of[long[[k]]]) - 1L
    for (t in 2:length(rows)) {
      r <- rows[[t]]
      d <- chain$DI[r, ]
      D <- dens[d, , drop = FALSE]
      if (t >= 3L) {
        G <- D * rep(dQ[r, ], each = length(d))
        G <- cbind(G, X[d, , drop = FALSE] * rowSums(G))
        out <- out + crossprod(W[origins, , drop = FALSE], G)
      }
      if (t < length(rows)) {
        fp <- state$fwd$FP[r, origins]
        V <- crossprod(
          state$probs$P[d, , drop = FALSE],
          W[origins, , drop = FALSE]
        )
        if (length(absorbing)) {
          V[absorbing, ] <- V[absorbing, ] + W[absorbing, ]
        }
        fd <- drop(fp %*% D)
        V[above] <- V[above] + fd
        V[diagonal] <- V[diagonal] - fd
        V[, slopes] <- V[, slopes] +
          crossprod(fp * (cbind(0, D) - cbind(D, 0)), X[d, , drop = FALSE])
        W <- V * (chain$C[r, ] / state$fwd$s[r])
      }
    }
  }
  out
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
