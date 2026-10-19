# Brute-force oracle for orm_markov_ml(): enumerate every completion of the
# missing and range-coarsened states allowed by each row, independently of the
# forward-backward recursion.

# Small K = 4 data set (state 4 absorbing) with fully observed patients plus
# hand-built patterns: single gap, two-visit gap, missing first visit, trailing
# missing visits, a gap ending in death, all post-baseline states missing, and
# Ocens() ranges inside and at the end of chains.
ml_oracle_data <- function(n_complete = 40L, seed = 11L) {
  raw <- make_markov_trajectories(
    n_patients = n_complete,
    follow_up_time = 6,
    n_states = 4,
    seed = seed
  )
  base <- prepare_markov_data(raw, absorbing_state = 4, ordered_response = FALSE)
  base$y <- as.numeric(as.character(base$y))
  base <- data.frame(
    id = base$id,
    time = base$time,
    y_lo = base$y,
    y_hi = base$y,
    yprev = as.numeric(as.character(base$yprev)),
    tx = base$tx
  )

  pattern <- function(id, baseline, tx, lo, hi = lo) {
    n <- length(lo)
    exact <- !is.na(lo) & !is.na(hi) & lo == hi
    yprev <- c(baseline, ifelse(exact[-n], lo[-n], NA))
    data.frame(
      id = id,
      time = seq_len(n),
      y_lo = lo,
      y_hi = hi,
      yprev = yprev,
      tx = tx
    )
  }
  patterns <- rbind(
    pattern(101, 2, 0, c(2, NA, 3, 2, 1)),
    pattern(102, 3, 1, c(3, NA, NA, 2, 2)),
    pattern(103, 1, 0, c(NA, 2, 2, 3, 3)),
    pattern(104, 2, 1, c(2, 3, NA, NA)),
    pattern(105, 3, 0, c(3, NA, 4)),
    pattern(106, 2, 1, c(NA, NA, NA, NA)),
    pattern(107, 2, 0, c(2, 1, 3, 2, -Inf), c(2, 3, 3, 2, 3)),
    pattern(108, 1, 1, c(1, 2, 1), c(1, 2, 3))
  )
  out <- rbind(base, patterns)
  out$yprev <- factor(out$yprev, levels = 1:3)
  out$y <- ifelse(out$y_lo == out$y_hi, out$y_lo, NA)
  out
}

# ml_oracle_data() without the patients who reach state 4: states 1 to 3, none
# of them absorbing, with the same missing-state and range patterns.
ml_oracle_data_no_absorbing <- function() {
  data <- ml_oracle_data()
  reaches_4 <- (!is.na(data$y_lo) & data$y_lo >= 4) |
    (!is.na(data$y_hi) & is.finite(data$y_hi) & data$y_hi >= 4)
  data[!data$id %in% data$id[reaches_4], , drop = FALSE]
}

ml_oracle_allowed <- function(lo, hi, K, absorbing, missing_states) {
  states <- seq_len(K)
  if (is.na(lo) && is.na(hi)) {
    if (identical(missing_states, "non_absorbing")) {
      return(setdiff(states, absorbing))
    }
    return(states)
  }
  if (is.infinite(lo)) {
    return(states[states < hi])
  }
  if (is.infinite(hi)) {
    return(states[states > lo])
  }
  states[states >= lo & states <= hi]
}

# Transition matrix for one data row, from the starting fit's design.
ml_oracle_transition <- function(start_fit, row, theta, K, absorbing) {
  n_int <- K - 1L
  origins <- setdiff(seq_len(K), absorbing)
  nd <- row[rep(1L, length(origins)), , drop = FALSE]
  nd$yprev <- factor(origins, levels = levels(row$yprev))
  X <- orm_model_matrix(start_fit, nd)
  eta <- drop(X %*% theta[-seq_len(n_int)])
  P <- diag(K)
  for (j in seq_along(origins)) {
    cum <- stats::plogis(theta[seq_len(n_int)] + eta[[j]])
    P[origins[[j]], ] <- c(1, cum) - c(cum, 0)
  }
  P
}

# Log-likelihood and pairwise posteriors (keyed "id_time") by enumeration.
ml_oracle <- function(
  data,
  start_fit,
  theta,
  K = 4L,
  absorbing = 4L,
  missing_states = "any"
) {
  data <- data[order(data$id, data$time), , drop = FALSE]
  loglik <- 0
  xi <- list()
  for (id in unique(data$id)) {
    rows <- data[data$id == id, , drop = FALSE]
    n <- nrow(rows)
    baseline <- as.integer(as.character(rows$yprev[[1L]]))
    allowed <- lapply(seq_len(n), function(t) {
      ml_oracle_allowed(rows$y_lo[[t]], rows$y_hi[[t]], K, absorbing, missing_states)
    })
    Ps <- lapply(seq_len(n), function(t) {
      ml_oracle_transition(start_fit, rows[t, , drop = FALSE], theta, K, absorbing)
    })
    paths <- as.matrix(expand.grid(allowed))
    probs <- apply(paths, 1L, function(path) {
      prev <- c(baseline, path)
      prod(vapply(seq_len(n), function(t) Ps[[t]][prev[[t]], prev[[t + 1L]]], 0))
    })
    L <- sum(probs)
    loglik <- loglik + log(L)
    for (t in seq_len(n)) {
      prev <- if (t == 1L) rep(baseline, nrow(paths)) else paths[, t - 1L]
      m <- matrix(0, K, K)
      for (k in seq_len(nrow(paths))) {
        m[prev[[k]], paths[k, t]] <- m[prev[[k]], paths[k, t]] + probs[[k]] / L
      }
      xi[[paste(id, rows$time[[t]], sep = "_")]] <- m
    }
  }
  list(loglik = loglik, xi = xi)
}

ml_oracle_problem <- function(
  data,
  formula = rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
  missing_states = "any"
) {
  markov_ml_setup(
    formula = formula,
    data = data,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    absorbing_state = 4,
    missing_states = missing_states,
    dots_fit = list()
  )
}

# Implementation posteriors as a list of K x K matrices keyed "id_time".
ml_problem_xi <- function(problem, state) {
  K <- problem$K
  xi_by_origin <- lapply(seq_len(problem$n_origin), function(ai) {
    markov_ml_xi(problem, state, ai)
  })
  # Mass already absorbed before a segment row stays in its absorbing state.
  absorbed <- matrix(0, nrow = problem$n, ncol = K)
  if (problem$chain$n) {
    absorbed[problem$chain$rows, ] <- state$fwd$FP * state$Q
  }
  out <- vector("list", problem$n)
  for (r in seq_len(problem$n)) {
    m <- matrix(0, K, K)
    for (ai in seq_len(problem$n_origin)) {
      m[problem$origins[[ai]], ] <- xi_by_origin[[ai]][r, ]
    }
    for (d in problem$absorbing) {
      m[d, d] <- absorbed[r, d]
    }
    out[[r]] <- m
  }
  names(out) <- paste(problem$data$id, problem$data$time, sep = "_")
  out
}

# Random coefficient vectors near the starting estimate with decreasing
# thresholds.
ml_random_theta <- function(theta, n_int, seed) {
  set.seed(seed)
  out <- theta + stats::rnorm(length(theta), sd = 0.25)
  out[seq_len(n_int)] <- sort(out[seq_len(n_int)], decreasing = TRUE)
  out
}

# Simulated K = 6 data with states missing completely at random.
ml_sim_missing <- function(
  n_patients = 120,
  follow_up_time = 12,
  missing = 0.2,
  seed = 1,
  treatment_effect = 0
) {
  data <- make_test_data(
    n_patients = n_patients,
    follow_up_time = follow_up_time,
    seed = seed,
    treatment_effect = treatment_effect
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  set.seed(seed + 1000)
  drop <- stats::runif(nrow(data)) < missing
  data$y[drop] <- NA
  after <- c(FALSE, drop[-nrow(data)]) & c(FALSE, data$id[-1] == data$id[-nrow(data)])
  data$yprev[after] <- NA
  data
}
