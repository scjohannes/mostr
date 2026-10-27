# Linear predictors of cumulative-logit transition models with missing states.
#
# Both frequentist engines describe the M = K - 1 linear predictors of a design
# row x (whose first column is the intercept) as eta = Gamma x, where column k
# of Gamma is C_k beta_k: C_k is the M x r_k constraint matrix of design column
# k and beta_k its r_k coefficients, as in VGAM::constraints() and
# get_effective_coefs(). An rms::orm() fit is the special case with C = I_M
# for the intercept and C = 1_M for every slope. F_j = plogis(eta_j) is
# P(Y > j) for j = 1, ..., M, and the probability of state b is
# F_{b-1} - F_b with F_0 = 1 and F_K = 0.
#
# Derivatives are computed for eta and mapped to the coefficients with the
# constraint matrices. Design columns whose constraint is 1_M ("parallel")
# share one computation; the others (the intercept and terms with
# threshold-specific effects) are handled one at a time, so a proportional
# odds model costs no more than with formulas written for it alone.

# `X` has one row per design row, starting with the intercept column, and one
# column per element of `constraints` (in the same order). The coefficients
# are the r_k coefficients of each column in turn, named `coef_names`.
markov_lp <- function(X, constraints, coef_names) {
  X <- as.matrix(X)
  constraints <- lapply(constraints, as.matrix)
  if (length(constraints) != ncol(X)) {
    stop("Internal error: one constraint matrix is needed per design column.")
  }
  M <- nrow(constraints[[1L]])
  sizes <- vapply(constraints, ncol, integer(1))
  if (sum(sizes) != length(coef_names)) {
    stop("Internal error: the constraints do not match the coefficients.")
  }
  ends <- cumsum(sizes)
  index <- Map(
    function(end, size) seq.int(end - size + 1L, length.out = size),
    ends,
    sizes
  )
  parallel <- vapply(
    constraints,
    function(C) ncol(C) == 1L && all(C == 1),
    logical(1)
  )
  identity <- vapply(
    constraints,
    function(C) ncol(C) == M && all(C == diag(M)),
    logical(1)
  )
  list(
    X = X,
    M = M,
    p = sum(sizes),
    constraints = constraints,
    index = index,
    coef_names = coef_names,
    parallel_cols = which(parallel),
    parallel_coef = unlist(index[parallel], use.names = FALSE),
    other = which(!parallel),
    identity = identity
  )
}

# Linear predictors of the design rows `rows` (all when NULL), one column per
# threshold.
markov_lp_eta <- function(lp, theta, rows = NULL) {
  X <- if (is.null(rows)) lp$X else lp$X[rows, , drop = FALSE]
  eta <- matrix(0, nrow = nrow(X), ncol = lp$M)
  if (length(lp$parallel_cols)) {
    eta <- eta +
      drop(X[, lp$parallel_cols, drop = FALSE] %*% theta[lp$parallel_coef])
  }
  for (k in lp$other) {
    gamma <- drop(lp$constraints[[k]] %*% theta[lp$index[[k]]])
    eta <- eta + outer(X[, k], gamma)
  }
  eta
}

# Map derivatives with respect to eta (one row per element of `rows`, one
# column per threshold) to derivatives with respect to the coefficients.
markov_lp_map <- function(lp, E, rows) {
  E <- matrix(E, ncol = lp$M)
  X <- lp$X[rows, , drop = FALSE]
  out <- matrix(0, nrow = nrow(E), ncol = lp$p)
  if (length(lp$parallel_cols)) {
    out[, lp$parallel_coef] <- X[, lp$parallel_cols, drop = FALSE] *
      rowSums(E)
  }
  for (k in lp$other) {
    EC <- if (lp$identity[[k]]) E else E %*% lp$constraints[[k]]
    out[, lp$index[[k]]] <- EC * X[, k]
  }
  out
}

# Rows j = 1..M of the derivatives of eta_j with respect to the coefficients,
# each weighted by the design columns: `E` is M x q, and row j of the result
# has entries C_k[j, ] * E[j, k] for the coefficients of design column k.
markov_lp_lift_rows <- function(lp, E) {
  out <- matrix(0, nrow = lp$M, ncol = lp$p)
  if (length(lp$parallel_cols)) {
    out[, lp$parallel_coef] <- E[, lp$parallel_cols, drop = FALSE]
  }
  for (k in lp$other) {
    out[, lp$index[[k]]] <- lp$constraints[[k]] * E[, k]
  }
  out
}

# sum_d sum_j w[d, j] g_dj g_dj', where g_dj is the derivative of eta_j of
# design row rows[d] with respect to the coefficients.
markov_lp_info_diag <- function(lp, w, rows) {
  w <- matrix(w, ncol = lp$M)
  X <- lp$X[rows, , drop = FALSE]
  P <- lp$parallel_cols
  Pc <- lp$parallel_coef
  J <- matrix(0, nrow = lp$p, ncol = lp$p)
  if (length(P)) {
    XP <- X[, P, drop = FALSE]
    J[Pc, Pc] <- crossprod(XP, XP * rowSums(w))
  }
  for (k in lp$other) {
    # T[j, l] = sum_d w[d, j] x_dk x_dl
    T <- crossprod(w * X[, k], X)
    Ck <- lp$constraints[[k]]
    ik <- lp$index[[k]]
    if (length(P)) {
      J[ik, Pc] <- crossprod(Ck, T[, P, drop = FALSE])
      J[Pc, ik] <- t(J[ik, Pc, drop = FALSE])
    }
    for (l in lp$other) {
      J[ik, lp$index[[l]]] <- crossprod(Ck, T[, l] * lp$constraints[[l]])
    }
  }
  J
}

# Proportional-odds terms of the outcome state `b` of design rows `rows`.
# `eta` holds the linear predictors of all design rows. u = b - 1 and v = b
# are the thresholds around b: F_u = P(Y >= b) (1 for b = 1) and
# F_v = P(Y > b) (0 for b = K), so the probability of b is F_u - F_v.
markov_lp_outcome <- function(lp, eta, rows, b) {
  M <- lp$M
  n <- length(b)
  upper <- b >= 2L
  lower <- b <= M
  F_u <- rep(1, n)
  F_u[upper] <- stats::plogis(eta[cbind(rows[upper], b[upper] - 1L)])
  F_v <- rep(0, n)
  F_v[lower] <- stats::plogis(eta[cbind(rows[lower], b[lower])])
  list(
    rows = rows,
    b = b,
    upper = upper,
    lower = lower,
    F_u = F_u,
    F_v = F_v,
    f_u = F_u * (1 - F_u),
    f_v = F_v * (1 - F_v),
    prob = F_u - F_v
  )
}

# Sums of `values` (a vector or a matrix with one row per element of `index`)
# by `index` in 1..size, as a size-row matrix.
markov_sum_by <- function(values, index, size) {
  values <- as.matrix(values)
  out <- matrix(0, nrow = size, ncol = ncol(values))
  if (length(index)) {
    out[sort(unique(index)), ] <- rowsum(values, index, reorder = TRUE)
  }
  out
}

# Weighted sums, within groups 1..n_groups, of the scores of log P(Y = b) for
# outcome terms `terms` (markov_lp_outcome()). The score for eta_u is
# f_u / p and for eta_v is -f_v / p; every other eta has none, so no
# row-by-threshold matrix is built for terms with one coefficient per
# threshold.
markov_lp_outcome_scores <- function(
  lp,
  terms,
  weights = 1,
  group = rep(1L, length(terms$b)),
  n_groups = 1L
) {
  M <- lp$M
  b <- terms$b
  upper <- terms$upper
  lower <- terms$lower
  w <- rep_len(weights, length(b)) / terms$prob
  g_u <- ifelse(upper, w * terms$f_u, 0)
  g_v <- ifelse(lower, -w * terms$f_v, 0)
  X <- lp$X[terms$rows, , drop = FALSE]
  out <- matrix(0, nrow = n_groups, ncol = lp$p)
  if (length(lp$parallel_cols)) {
    out[, lp$parallel_coef] <- markov_sum_by(
      X[, lp$parallel_cols, drop = FALSE] * (g_u + g_v),
      group,
      n_groups
    )
  }
  for (k in lp$other) {
    xk <- X[, k]
    if (lp$identity[[k]]) {
      # Coefficient j of group g sits at position (j - 1) * n_groups + g.
      position <- c(
        (b[upper] - 2L) * n_groups + group[upper],
        (b[lower] - 1L) * n_groups + group[lower]
      )
      sums <- markov_sum_by(
        c((g_u * xk)[upper], (g_v * xk)[lower]),
        position,
        n_groups * M
      )
      out[, lp$index[[k]]] <- matrix(sums, nrow = n_groups)
    } else {
      Ck <- lp$constraints[[k]]
      E <- matrix(0, nrow = length(b), ncol = ncol(Ck))
      E[upper, ] <- Ck[b[upper] - 1L, , drop = FALSE] * (g_u * xk)[upper]
      E[lower, ] <- E[lower, , drop = FALSE] +
        Ck[b[lower], , drop = FALSE] * (g_v * xk)[lower]
      out[, lp$index[[k]]] <- markov_sum_by(E, group, n_groups)
    }
  }
  colnames(out) <- lp$coef_names
  out
}

# Observed information (negative Hessian) of sum(weights * log P(Y = b)) for
# outcome terms `terms`.
#
# With f'(u) = f(u) (1 - 2 F(u)), the second derivatives of log p_b are
# f'(u) / p - f(u)^2 / p^2 in eta_u, -f'(v) / p - f(v)^2 / p^2 in eta_v, and
# f(u) f(v) / p^2 in eta_u and eta_v. The information is therefore
# tridiagonal in eta; it is mapped to the coefficients with the constraint
# matrices.
markov_lp_outcome_information <- function(lp, terms, weights = 1) {
  M <- lp$M
  b <- terms$b
  upper <- terms$upper
  lower <- terms$lower
  both <- upper & lower
  w <- rep_len(weights, length(b))
  prob <- terms$prob
  fu <- terms$f_u
  fv <- terms$f_v
  w_uu <- ifelse(
    upper,
    w * (fu^2 / prob^2 - fu * (1 - 2 * terms$F_u) / prob),
    0
  )
  w_vv <- ifelse(
    lower,
    w * (fv^2 / prob^2 + fv * (1 - 2 * terms$F_v) / prob),
    0
  )
  w_uv <- ifelse(both, -w * fu * fv / prob^2, 0)

  X <- lp$X[terms$rows, , drop = FALSE]
  P <- lp$parallel_cols
  Pc <- lp$parallel_coef
  J <- matrix(0, nrow = lp$p, ncol = lp$p)
  if (length(P)) {
    XP <- X[, P, drop = FALSE]
    J[Pc, Pc] <- crossprod(XP, XP * (w_uu + w_vv + 2 * w_uv))
  }
  for (k in lp$other) {
    xk <- X[, k]
    # Diagonal entries of the eta information by threshold j, and the
    # entries linking thresholds j and j + 1, each times x_k x_l.
    T_diag <- markov_sum_by(
      X[upper, , drop = FALSE] * (w_uu * xk)[upper],
      b[upper] - 1L,
      M
    ) +
      markov_sum_by(
        X[lower, , drop = FALSE] * (w_vv * xk)[lower],
        b[lower],
        M
      )
    T_off <- if (M > 1L) {
      markov_sum_by(
        X[both, , drop = FALSE] * (w_uv * xk)[both],
        b[both] - 1L,
        M - 1L
      )
    }
    Ck <- lp$constraints[[k]]
    ik <- lp$index[[k]]
    block <- function(l, Cl) {
      out <- crossprod(Ck, T_diag[, l] * Cl)
      if (M > 1L) {
        out <- out +
          crossprod(Ck[-M, , drop = FALSE], T_off[, l] * Cl[-1L, , drop = FALSE]) +
          crossprod(Ck[-1L, , drop = FALSE], T_off[, l] * Cl[-M, , drop = FALSE])
      }
      out
    }
    if (length(P)) {
      ones <- matrix(1, nrow = M, ncol = 1L)
      cross <- matrix(0, nrow = length(ik), ncol = length(P))
      for (i in seq_along(P)) {
        cross[, i] <- block(P[[i]], ones)
      }
      J[ik, Pc] <- cross
      J[Pc, ik] <- t(cross)
    }
    for (l in lp$other) {
      J[ik, lp$index[[l]]] <- block(l, lp$constraints[[l]])
    }
  }
  dimnames(J) <- list(lp$coef_names, lp$coef_names)
  J
}
