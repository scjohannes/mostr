# Observed-data likelihood and its maximisation with missing states.
#
# A scaled forward pass over the rows that need the recursion gives the
# log-likelihood; the matching backward pass gives the posterior transition
# weights used by the EM M-step and by the exact (Fisher-identity) score. The
# engines differ only in the linear predictors (`problem$lp`, see
# markov-ml-lp.R).

# ---------------------------------------------------------------------------
# Transition probabilities
# ---------------------------------------------------------------------------

# Probabilities of every design row. Returns NULL when `strict = FALSE` and
# the coefficients give invalid probabilities (used by step halving).
markov_ml_transition_probs <- function(problem, theta, strict = TRUE) {
  eta <- markov_lp_eta(problem$lp, theta)
  cum <- stats::plogis(eta)
  P <- cbind(1, cum) - cbind(cum, 0)
  if (any(P < -1e-10) || any(abs(rowSums(P) - 1) > 1e-12)) {
    if (!strict) {
      return(NULL)
    }
    stop(
      "Transition probabilities are invalid: some state has a negative ",
      "probability. With threshold-specific effects, the coefficients may ",
      "give crossing thresholds.",
      call. = FALSE
    )
  }
  P[P < 0] <- 0
  list(P = P, eta = eta, cum = cum, dens = cum * (1 - cum), theta = theta)
}

# ---------------------------------------------------------------------------
# Forward-backward recursion
# ---------------------------------------------------------------------------

# Scaled forward pass over the segments that need the recursion. `FP` holds
# the forward vector before each segment row (the filtered distribution of the
# previous state) and `s` the scaling constants; sum(log(s)) is the segments'
# log-likelihood. Only states the model can leave (`origins`) occur before
# the last row of a segment (checked by markov_ml_prepare()).
markov_ml_forward <- function(problem, probs) {
  chain <- problem$chain
  K <- problem$K
  P <- probs$P
  origins <- problem$origins

  cur <- matrix(0, nrow = chain$N, ncol = K)
  cur[cbind(seq_len(chain$N), chain$baseline)] <- 1
  FP <- matrix(0, nrow = chain$n, ncol = K)
  s <- numeric(chain$n)

  for (k in seq_along(chain$rows_by_pos)) {
    r <- chain$rows_by_pos[[k]]
    i <- chain$pats_by_pos[[k]]
    fp <- cur[i, , drop = FALSE]
    FP[r, ] <- fp
    new <- matrix(0, nrow = length(r), ncol = K)
    for (ai in seq_along(origins)) {
      new <- new + fp[, origins[[ai]]] * P[chain$DI[r, ai], , drop = FALSE]
    }
    new <- new * chain$C[r, , drop = FALSE]
    sk <- rowSums(new)
    s[r] <- sk
    cur[i, ] <- new / sk
  }

  list(FP = FP, s = s, loglik = sum(log(s)), valid = all(s > 0))
}

markov_ml_stop_impossible <- function(ids) {
  stop(
    "Observed states have zero probability under the transition model for ",
    "patient(s) ",
    markov_format_patient_ids(ids),
    ".",
    call. = FALSE
  )
}

# Scaled backward pass. Returns Q = c * g / s for each segment row, so the
# posterior of a transition is xi(a, b) = FP[a] * P[a, b] * Q[b].
markov_ml_backward <- function(problem, probs, fwd) {
  chain <- problem$chain
  K <- problem$K
  P <- probs$P
  origins <- problem$origins

  Q <- matrix(0, nrow = chain$n, ncol = K)
  gcur <- matrix(1, nrow = chain$N, ncol = K)
  for (k in rev(seq_along(chain$rows_by_pos))) {
    r <- chain$rows_by_pos[[k]]
    i <- chain$pats_by_pos[[k]]
    q <- chain$C[r, , drop = FALSE] * gcur[i, , drop = FALSE] / fwd$s[r]
    Q[r, ] <- q
    gp <- matrix(0, nrow = length(r), ncol = K)
    for (ai in seq_along(origins)) {
      gp[, origins[[ai]]] <- rowSums(P[chain$DI[r, ai], , drop = FALSE] * q)
    }
    gcur[i, ] <- gp
  }
  Q
}

# Posterior transition weights of the segment rows for origin index `ai`.
markov_ml_chain_xi <- function(problem, state, ai) {
  a <- problem$origins[[ai]]
  state$fwd$FP[, a] *
    state$probs$P[problem$chain$DI[, ai], , drop = FALSE] *
    state$Q
}

# Log-likelihood terms of the fully observed transitions, by aggregated
# (design row, state).
markov_ml_complete_terms <- function(problem, eta) {
  comp <- problem$complete_part
  markov_lp_outcome(problem$lp, eta, comp$agg_design, comp$agg_state)
}

# Forward-backward pass at `theta`. With `strict = FALSE` an invalid
# coefficient vector or a zero-probability observation returns NULL instead of
# an error; callers use this only for step control.
markov_ml_estep <- function(problem, theta, strict = TRUE, backward = TRUE) {
  probs <- markov_ml_transition_probs(problem, theta, strict = strict)
  if (is.null(probs)) {
    return(NULL)
  }
  comp <- problem$complete_part
  p_complete <- markov_ml_complete_terms(problem, probs$eta)$prob
  if (any(!(p_complete > 0))) {
    if (!strict) {
      return(NULL)
    }
    p_rows <- probs$P[cbind(comp$design, comp$state)]
    markov_ml_stop_impossible(problem$ids[comp$rows[!(p_rows > 0)]])
  }
  loglik <- sum(comp$agg_count * log(p_complete))
  out <- list(theta = theta, probs = probs, loglik = loglik)
  if (problem$chain$n) {
    fwd <- markov_ml_forward(problem, probs)
    if (!fwd$valid) {
      if (!strict) {
        return(NULL)
      }
      markov_ml_stop_impossible(problem$chain$ids[!(fwd$s > 0)])
    }
    out$fwd <- fwd
    out$loglik <- loglik + fwd$loglik
    if (backward) {
      out$Q <- markov_ml_backward(problem, probs, fwd)
    }
  }
  out
}

markov_ml_loglik <- function(problem, theta, strict = TRUE) {
  state <- markov_ml_estep(problem, theta, strict = strict, backward = FALSE)
  if (is.null(state)) {
    return(NA_real_)
  }
  state$loglik
}

# Expected complete-data counts aggregated by (design row, destination
# state): the observed counts of the fully observed transitions plus the
# posterior weights of the rows that need the recursion.
markov_ml_aggregate_weights <- function(problem, state) {
  comp <- problem$complete_part
  K <- problem$K
  n_design <- nrow(problem$lp$X)
  W <- matrix(0, nrow = n_design, ncol = K)
  W[cbind(comp$agg_design, comp$agg_state)] <- comp$agg_count
  if (problem$chain$n) {
    for (ai in seq_len(problem$n_origin)) {
      xi <- markov_ml_chain_xi(problem, state, ai)
      agg <- rowsum(xi, problem$chain$DI[, ai], reorder = TRUE)
      rows <- as.integer(rownames(agg))
      W[rows, ] <- W[rows, ] + agg
    }
  }
  idx <- which(W > 0, arr.ind = TRUE)
  list(design = idx[, 1L], state = idx[, 2L], weight = W[idx])
}

# ---------------------------------------------------------------------------
# Maximisation
# ---------------------------------------------------------------------------

# Newton direction for gradient `g` and information `H`, falling back to a
# Levenberg-Marquardt step when `H` is not positive definite.
markov_ml_direction <- function(g, H) {
  R <- tryCatch(chol(H), error = function(e) NULL)
  lambda <- 1e-8 * max(1, max(abs(diag(H))))
  while (is.null(R) && lambda < 1e8) {
    R <- tryCatch(chol(H + diag(lambda, nrow(H))), error = function(e) NULL)
    lambda <- lambda * 100
  }
  if (is.null(R)) {
    return(NULL)
  }
  drop(backsolve(R, forwardsolve(t(R), g)))
}

# M-step: maximise the expected complete-data log-likelihood
# sum(weight * log P(Y = state)) over the aggregated rows `counts` by Newton
# steps with step halving, starting from `theta`.
markov_ml_mstep <- function(problem, counts, theta, maxit = 25L) {
  lp <- problem$lp
  objective <- function(theta) {
    eta <- markov_lp_eta(lp, theta)
    terms <- markov_lp_outcome(lp, eta, counts$design, counts$state)
    if (any(!(terms$prob > 0))) {
      return(NULL)
    }
    list(terms = terms, value = sum(counts$weight * log(terms$prob)))
  }
  current <- objective(theta)
  if (is.null(current)) {
    stop("Internal error: the M-step starts from invalid coefficients.")
  }
  for (iter in seq_len(maxit)) {
    g <- drop(markov_lp_outcome_scores(lp, current$terms, counts$weight))
    H <- markov_lp_outcome_information(lp, current$terms, counts$weight)
    direction <- markov_ml_direction(g, H)
    if (is.null(direction) || any(!is.finite(direction))) {
      break
    }
    if (sum(g * direction) < 1e-12) {
      break
    }
    t <- 1
    accepted <- FALSE
    while (t >= 1 / 1024) {
      candidate <- theta + t * direction
      next_value <- objective(candidate)
      if (
        !is.null(next_value) &&
          next_value$value >= current$value - 1e-12 * max(1, abs(current$value))
      ) {
        accepted <- TRUE
        break
      }
      t <- t / 2
    }
    if (!accepted) {
      break
    }
    theta <- candidate
    current <- next_value
  }
  names(theta) <- lp$coef_names
  theta
}

markov_ml_em <- function(problem, theta, control) {
  state <- markov_ml_estep(problem, theta)
  trace <- state$loglik
  converged <- FALSE
  iterations <- 0L
  for (iter in seq_len(control$maxit)) {
    iterations <- iter
    counts <- markov_ml_aggregate_weights(problem, state)
    theta_new <- markov_ml_mstep(problem, counts, theta)
    state_new <- markov_ml_estep(problem, theta_new)
    trace <- c(trace, state_new$loglik)
    dll <- state_new$loglik - state$loglik
    dtheta <- max(abs(theta_new - theta))
    theta <- theta_new
    state <- state_new
    if (
      abs(dll) <= control$reltol * (abs(state$loglik) + control$reltol) ||
        dtheta < control$coef_tol
    ) {
      converged <- TRUE
      break
    }
  }
  list(
    theta = theta,
    state = state,
    trace = trace,
    iterations = iterations,
    converged = converged
  )
}

markov_ml_newton <- function(problem, theta, loglik, information, control) {
  steps <- 0L
  for (step in seq_len(control$newton_steps)) {
    score <- markov_ml_total_score(problem, theta)
    direction <- tryCatch(
      drop(solve(information, score)),
      error = function(e) NULL
    )
    if (is.null(direction) || any(!is.finite(direction))) {
      break
    }
    # Newton decrement: roughly how much the log-likelihood could still rise.
    if (sum(score * direction) < 1e-20) {
      break
    }
    t <- 1
    accepted <- FALSE
    while (t >= 1 / 1024) {
      candidate <- theta + t * direction
      ll <- markov_ml_loglik(problem, candidate, strict = FALSE)
      if (is.finite(ll) && ll >= loglik - 1e-12 * max(1, abs(loglik))) {
        accepted <- TRUE
        break
      }
      t <- t / 2
    }
    if (!accepted) {
      break
    }
    theta <- candidate
    loglik <- ll
    steps <- steps + 1L
  }
  list(theta = theta, loglik = loglik, steps = steps)
}

# Settings of the maximisation. `variance = FALSE` (bootstrap refits) skips
# the covariances.
markov_ml_control <- function(variance = TRUE) {
  list(
    maxit = 200L,
    reltol = 1e-6,
    coef_tol = 1e-7,
    newton_steps = 10L,
    variance = variance
  )
}
