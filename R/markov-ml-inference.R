# Scores, observed information, and covariance with missing states.

# Exact per-row observed-data scores of the rows that need the recursion, by
# Fisher's identity: the complete-data score averaged over the posterior
# transition weights. With xi(a, b) / P(a, b) = FP[a] Q[b], the score for
# eta_j of the transition from origin a is FP[a] f_j (Q_{j+1} - Q_j).
markov_ml_chain_scores <- function(problem, state) {
  chain <- problem$chain
  K <- problem$K
  lp <- problem$lp
  if (!chain$n) {
    return(matrix(0, nrow = 0L, ncol = lp$p))
  }
  dQ <- state$Q[, -1L, drop = FALSE] - state$Q[, -K, drop = FALSE]
  out <- matrix(0, nrow = chain$n, ncol = lp$p)
  for (ai in seq_len(problem$n_origin)) {
    a <- problem$origins[[ai]]
    d <- chain$DI[, ai]
    D <- state$fwd$FP[, a] * state$probs$dens[d, , drop = FALSE] * dQ
    out <- out + markov_lp_map(lp, D, d)
  }
  out
}

# Per-patient observed-data scores (one row per patient of
# `problem$patient_ids`). `state` is the forward-backward pass at `theta`; it
# is needed only when some rows need the recursion.
markov_ml_patient_scores <- function(problem, theta, state = NULL) {
  lp <- problem$lp
  comp <- problem$complete_part
  eta <- if (is.null(state)) markov_lp_eta(lp, theta) else state$probs$eta
  U <- markov_lp_outcome_scores(
    lp,
    markov_lp_outcome(lp, eta, comp$design, comp$state),
    group = comp$patient,
    n_groups = problem$N
  )
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
  dimnames(U) <- list(problem$patient_ids, lp$coef_names)
  U
}

# Score and information of the fully observed transitions, on their
# aggregated rows.
markov_ml_complete_score <- function(problem, theta, eta = NULL) {
  lp <- problem$lp
  eta <- eta %||% markov_lp_eta(lp, theta)
  drop(markov_lp_outcome_scores(
    lp,
    markov_ml_complete_terms(problem, eta),
    weights = problem$complete_part$agg_count
  ))
}

markov_ml_complete_information <- function(problem, theta, eta = NULL) {
  lp <- problem$lp
  eta <- eta %||% markov_lp_eta(lp, theta)
  markov_lp_outcome_information(
    lp,
    markov_ml_complete_terms(problem, eta),
    weights = problem$complete_part$agg_count
  )
}

markov_ml_total_score <- function(problem, theta) {
  if (!problem$chain$n) {
    out <- markov_ml_complete_score(problem, theta)
  } else {
    state <- markov_ml_estep(problem, theta)
    out <- markov_ml_complete_score(problem, theta, state$probs$eta) +
      colSums(markov_ml_chain_scores(problem, state))
  }
  names(out) <- problem$lp$coef_names
  out
}

# Exact observed information: that of the fully observed transitions plus
# that of the segment rows (markov_ml_chain_information()). `state` is the
# forward-backward pass at `theta`; it is computed when NULL and some rows
# need the recursion.
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
  dimnames(J) <- list(problem$lp$coef_names, problem$lp$coef_names)
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
# state and is handled for all segments at once, so segments of one or two
# rows cost about one recursion pass; only longer segments carry a dense
# K x p matrix per row.
markov_ml_chain_information <- function(problem, state) {
  chain <- problem$chain
  K <- problem$K
  lp <- problem$lp
  FP <- state$fwd$FP
  dQ <- state$Q[, -1L, drop = FALSE] - state$Q[, -K, drop = FALSE]

  # Within rows. P(a, b) = F_{b-1} - F_b with F_j = plogis(eta_j), so the
  # second derivative of P(a, b) in eta_j is f'_j = f_j (1 - 2 F_j) for
  # b = j + 1 and -f'_j for b = j, and none links two thresholds. Hence
  # sum_b Q(b) d2P(a, b) needs only f'_j (Q(j + 1) - Q(j)), and
  # FP_t(a) (Q_t(j + 1) - Q_t(j)) is summed by design row first.
  R <- matrix(0, nrow = nrow(lp$X), ncol = lp$M)
  for (ai in seq_len(problem$n_origin)) {
    agg <- rowsum(
      FP[, problem$origins[[ai]]] * dQ,
      chain$DI[, ai],
      reorder = TRUE
    )
    rows <- as.integer(rownames(agg))
    R[rows, ] <- R[rows, ] + agg
  }
  used <- which(rowSums(R != 0) > 0)
  Mw <- (state$probs$dens * (1 - 2 * state$probs$cum) * R)[used, , drop = FALSE]
  J <- -markov_lp_info_diag(lp, Mw, used)

  segment_scores <- rowsum(
    markov_ml_chain_scores(problem, state),
    chain$segment,
    reorder = TRUE
  )
  between <- markov_ml_chain_between(problem, state, dQ)
  J + crossprod(segment_scores) - between - t(between)
}

# W_1 of the first rows `r` of segments (see markov_ml_chain_information()).
# With b0 the observed state before the segment, d0 the design row of the
# transition from it, f the densities f_j of its thresholds, and
# c = C_1 / s_1, W_1(b) = c_b dP(b0, b) = c_b (f_{b-1} g_{b-1} - f_b g_b),
# where g_j is the derivative of eta_j of d0 with respect to the
# coefficients.
markov_ml_first_row_terms <- function(problem, state, r) {
  chain <- problem$chain
  ai0 <- match(chain$baseline[chain$segment[r]], problem$origins)
  d0 <- chain$DI[cbind(r, ai0)]
  list(
    c = chain$C[r, , drop = FALSE] / state$fwd$s[r],
    f = state$probs$dens[d0, , drop = FALSE],
    d0 = d0
  )
}

# Between-row terms of markov_ml_chain_information(): the sum over segments
# of sum_t sum_a W_{t-1}(a) G_t(a)' (p x p).
markov_ml_chain_between <- function(problem, state, dQ) {
  chain <- problem$chain
  K <- problem$K
  lp <- problem$lp
  M <- lp$M
  p <- lp$p
  origins <- problem$origins
  dens <- state$probs$dens
  out <- matrix(0, nrow = p, ncol = p)
  if (length(chain$rows_by_pos) < 2L) {
    return(out)
  }

  # Second rows, all segments at once. A segment's rows are consecutive, so
  # its first row precedes its second. Origin a adds c_a (f_{a-1} g_{a-1} -
  # f_a g_a) G(a)'. For design columns with one coefficient shared by all
  # thresholds (g_j the same for every j), this is c_a (f_{a-1} - f_a) x0
  # G(a)'; for the others, the rows of their constraint matrices pick the
  # coefficients of eta_{a-1} and eta_a.
  r2 <- chain$rows_by_pos[[2L]]
  first <- markov_ml_first_row_terms(problem, state, r2 - 1L)
  f <- first$f
  cc <- first$c
  X0 <- lp$X[first$d0, , drop = FALSE]
  P <- lp$parallel_cols
  Pc <- lp$parallel_coef
  parallel_weight <- cc * (cbind(0, f) - cbind(f, 0))
  parallel_terms <- matrix(0, nrow = length(r2), ncol = p)
  dQ2 <- dQ[r2, , drop = FALSE]
  for (ai in seq_along(origins)) {
    a <- origins[[ai]]
    d <- chain$DI[r2, ai]
    G <- markov_lp_map(lp, dens[d, , drop = FALSE] * dQ2, d)
    if (length(P)) {
      parallel_terms <- parallel_terms + parallel_weight[, a] * G
    }
    for (k in lp$other) {
      ik <- lp$index[[k]]
      Ck <- lp$constraints[[k]]
      xk <- X0[, k]
      if (a >= 2L) {
        out[ik, ] <- out[ik, ] +
          outer(Ck[a - 1L, ], colSums((cc[, a] * f[, a - 1L] * xk) * G))
      }
      if (a <= M) {
        out[ik, ] <- out[ik, ] -
          outer(Ck[a, ], colSums((cc[, a] * f[, a] * xk) * G))
      }
    }
  }
  if (length(P)) {
    out[Pc, ] <- out[Pc, ] +
      crossprod(X0[, P, drop = FALSE], parallel_terms)
  }

  # Third and later rows: segment by segment with a dense W.
  length_of <- tabulate(chain$segment, chain$N)
  long <- which(length_of >= 3L)
  if (!length(long)) {
    return(out)
  }
  start <- integer(chain$N)
  start[chain$pats_by_pos[[1L]]] <- chain$rows_by_pos[[1L]]
  first <- markov_ml_first_row_terms(problem, state, start[long])
  # Rows b of the derivative of P(b0, b): f_{b-1} at eta_{b-1} and -f_b at
  # eta_b.
  step_rows <- cbind(seq_len(M) + 1L, seq_len(M))
  diagonal <- cbind(seq_len(M), seq_len(M))
  for (k in seq_along(long)) {
    E1 <- matrix(0, nrow = K, ncol = M)
    E1[step_rows] <- first$f[k, ]
    E1[diagonal] <- -first$f[k, ]
    W <- first$c[k, ] * markov_lp_map(lp, E1, rep(first$d0[[k]], K))
    rows <- start[long[[k]]] + seq_len(length_of[long[[k]]]) - 1L
    for (t in 2:length(rows)) {
      r <- rows[[t]]
      d <- chain$DI[r, ]
      D <- dens[d, , drop = FALSE]
      if (t >= 3L) {
        G <- markov_lp_map(lp, D * rep(dQ[r, ], each = length(d)), d)
        out <- out + crossprod(W[origins, , drop = FALSE], G)
      }
      if (t < length(rows)) {
        fp <- state$fwd$FP[r, origins]
        V <- crossprod(
          state$probs$P[d, , drop = FALSE],
          W[origins, , drop = FALSE]
        )
        # sum_a FP(a) dP(a, b): row b gains R_{b-1} - R_b, where R_j is
        # sum_a FP(a) f_j(a) g_j(a).
        Rj <- markov_lp_lift_rows(
          lp,
          crossprod(fp * D, lp$X[d, , drop = FALSE])
        )
        V[-1L, ] <- V[-1L, ] + Rj
        V[-K, ] <- V[-K, ] - Rj
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
# Patient influence components of fitted models
# ---------------------------------------------------------------------------

markov_ml_reject_cluster <- function(cluster) {
  if (!is.null(cluster)) {
    stop(
      "`cluster` cannot be supplied for a fit with `missing_states` or an ",
      "`Ocens()` response: its clusters are the patients given by `id_var`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Stored components of a fit with missing states, checked for inference.
markov_ml_stored <- function(model) {
  ml <- markov_ml_components(model)
  if (is.null(ml$patient_scores) || is.null(ml$bread)) {
    stop(
      "This fit has no stored patient scores and covariance, so analytical ",
      "inference is not available.",
      call. = FALSE
    )
  }
  if (anyNA(ml$bread)) {
    stop(
      "The observed information of this fit is not positive definite, so no ",
      "covariance is available.",
      call. = FALSE
    )
  }
  ml
}

# Coefficients, inverse observed information, and per-patient scores of a fit
# with missing states: one score row per patient with fitted rows, in the
# order of the stored starting profiles.
markov_ml_influence_components <- function(model, cluster = NULL) {
  markov_ml_reject_cluster(cluster)
  ml <- markov_ml_stored(model)
  list(
    coefficients = get_coef(model),
    bread = ml$bread,
    scores = ml$patient_scores,
    cluster = rownames(ml$patient_scores)
  )
}

# Robust covariance of a fit with missing states, recomputed from the stored
# patient scores when another correction is requested.
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
