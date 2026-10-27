# Brute-force oracle for fits with missing states: enumerate every completion
# of the missing and range-coded states allowed by each row, independently of
# the forward-backward recursion. Transition probabilities come from the
# engine fit's effective coefficients (get_effective_coefs()), as in SOP
# prediction.

# One patient's rows with states `lo` to `hi` (NA for a missing state) and
# `yprev` set from the previous row where that state is exact.
ml_oracle_pattern <- function(id, baseline, tx, lo, hi = lo) {
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

# Small K = 4 data set (state 4 is death) with fully observed patients plus
# hand-built patterns: single gap, two-visit gap, missing first visit,
# trailing missing visits, a gap ending in death, all post-baseline states
# missing, and Ocens() ranges inside and at the end of chains.
# `long_runs = TRUE` adds runs of three to five missing states.
ml_oracle_data <- function(n_complete = 60L, seed = 11L, long_runs = FALSE) {
  raw <- make_markov_trajectories(
    n_patients = n_complete,
    follow_up_time = 6,
    n_states = 4,
    seed = seed
  )
  base <- prepare_markov_data(raw, absorbing_state = 4)
  base <- data.frame(
    id = base$id,
    time = base$time,
    y_lo = as.numeric(as.character(base$y)),
    y_hi = as.numeric(as.character(base$y)),
    yprev = as.numeric(as.character(base$yprev)),
    tx = base$tx
  )

  pattern <- ml_oracle_pattern
  patterns <- rbind(
    pattern(101, 2, 0, c(2, NA, 3, 2, 1)),
    pattern(102, 3, 1, c(3, NA, NA, 2, 2)),
    pattern(103, 1, 0, c(NA, 2, 2, 3, 3)),
    pattern(104, 2, 1, c(2, 3, NA, NA)),
    pattern(105, 3, 0, c(3, NA, 4)),
    pattern(106, 2, 1, c(NA, NA, NA, NA)),
    pattern(107, 2, 0, c(2, 1, 3, 2, -Inf), c(2, 3, 3, 2, 3)),
    pattern(108, 1, 1, c(1, 2, 1), c(1, 2, 4))
  )
  if (long_runs) {
    patterns <- rbind(
      patterns,
      pattern(109, 2, 0, c(2, NA, NA, NA, NA, 3, 2)),
      pattern(110, 3, 1, c(3, NA, NA, NA, 4)),
      pattern(111, 1, 1, c(NA, NA, NA, NA, NA)),
      pattern(112, 2, 0, c(NA, 1, NA, NA, 2, NA), c(NA, 2, NA, NA, 2, NA))
    )
  }
  out <- rbind(base, patterns)
  out$yprev <- factor(out$yprev, levels = 1:3)
  out$y <- ifelse(out$y_lo == out$y_hi, out$y_lo, NA)
  out
}

ml_oracle_allowed <- function(lo, hi, K, missing_states) {
  states <- seq_len(K)
  if (is.na(lo) && is.na(hi)) {
    return(missing_states)
  }
  if (is.infinite(lo)) {
    return(states[states < hi])
  }
  if (is.infinite(hi)) {
    return(states[states > lo])
  }
  states[states >= lo & states <= hi]
}

# Set up a problem as mostm() does and return it with its starting fit.
ml_test_problem <- function(
  data,
  formula = Ocens(y_lo, y_hi) ~ yprev + time + tx,
  engine = "orm",
  missing_states = 1:3,
  engine_args = list()
) {
  setup <- list(
    formula = formula,
    data = data,
    dots = NULL,
    subset_rows = NULL,
    prior_weight_column = NULL,
    eval_env = environment(),
    data_name = "data"
  )
  problem <- mostm_ml_problem(
    setup,
    engine = engine,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    missing_states = missing_states,
    engine_args = engine_args
  )
  start <- markov_ml_start(problem, engine, setup, engine_args)
  problem$lp <- start$lp
  problem$start_fit <- start$fit
  problem$theta <- start$theta
  problem
}

# K x K transition matrix of `row` for the coefficients `theta`, from the
# engine fit's effective coefficients and its design for each origin.
ml_oracle_transition <- function(problem, row, theta) {
  K <- problem$K
  fit <- problem$start_fit
  origins <- problem$origins
  nd <- row[rep(1L, length(origins)), , drop = FALSE]
  nd$yprev <- factor(origins, levels = levels(row$yprev))
  Gamma <- get_effective_coefs(fit, beta = theta)
  X <- if (inherits(fit, "vglm")) {
    markov_ml_vglm_design(fit, nd, colnames(Gamma))
  } else {
    orm_model_matrix(fit, nd, include_intercept = TRUE)[, colnames(Gamma)]
  }
  cum <- stats::plogis(X %*% t(Gamma))
  P <- matrix(0, K, K)
  P[origins, ] <- cbind(1, cum) - cbind(cum, 0)
  P
}

# Log-likelihood by enumeration of every completion of each patient's rows.
ml_oracle_loglik <- function(data, problem, theta, missing_states = 1:3) {
  data <- data[order(data$id, data$time), , drop = FALSE]
  K <- problem$K
  loglik <- 0
  for (id in unique(data$id)) {
    rows <- data[data$id == id, , drop = FALSE]
    n <- nrow(rows)
    baseline <- as.integer(as.character(rows$yprev[[1L]]))
    allowed <- lapply(seq_len(n), function(t) {
      ml_oracle_allowed(rows$y_lo[[t]], rows$y_hi[[t]], K, missing_states)
    })
    Ps <- lapply(seq_len(n), function(t) {
      ml_oracle_transition(problem, rows[t, , drop = FALSE], theta)
    })
    paths <- as.matrix(expand.grid(allowed))
    probs <- apply(paths, 1L, function(path) {
      prev <- c(baseline, path)
      prod(vapply(seq_len(n), function(t) Ps[[t]][prev[[t]], prev[[t + 1L]]], 0))
    })
    loglik <- loglik + log(sum(probs))
  }
  loglik
}

# Central-difference gradient of `f` at `x`.
ml_numeric_gradient <- function(f, x, h = 1e-5) {
  vapply(
    seq_along(x),
    function(j) {
      e <- replace(numeric(length(x)), j, h)
      (f(x + e) - f(x - e)) / (2 * h)
    },
    numeric(1)
  )
}

# Central-difference Jacobian of the vector-valued `f` at `x`.
ml_numeric_jacobian <- function(f, x, h = 1e-5) {
  out <- vapply(
    seq_along(x),
    function(j) {
      e <- replace(numeric(length(x)), j, h)
      (f(x + e) - f(x - e)) / (2 * h)
    },
    numeric(length(f(x)))
  )
  matrix(out, ncol = length(x))
}

# Random coefficient vectors near `theta` with valid probabilities.
ml_random_theta <- function(problem, theta, seed, sd = 0.1) {
  set.seed(seed)
  for (attempt in 1:50) {
    out <- theta + stats::rnorm(length(theta), sd = sd)
    if (!is.null(markov_ml_transition_probs(problem, out, strict = FALSE))) {
      return(out)
    }
  }
  stop("No valid random coefficients found.")
}

# Simulated K = 6 data (state 6 is death) with living states missing
# completely at random; the previous state after a missing visit is unknown.
ml_sim_missing <- function(
  n_patients = 120,
  follow_up_time = 12,
  missing = 0.2,
  seed = 1
) {
  data <- make_test_data(
    n_patients = n_patients,
    follow_up_time = follow_up_time,
    seed = seed
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$y <- as.numeric(as.character(data$y))
  set.seed(seed + 1000)
  drop <- data$y != 6 & stats::runif(nrow(data)) < missing
  data$y[drop] <- NA
  after <- c(FALSE, drop[-nrow(data)]) &
    c(FALSE, data$id[-1] == data$id[-nrow(data)])
  data$yprev[after] <- NA
  data
}
