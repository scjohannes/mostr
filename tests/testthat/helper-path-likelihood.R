# A readable reference for cpp_path_loglik(), the log-likelihood of a
# proportional-odds Markov model when some visits are only known to lie in a
# range of states. It calls no package code and sums over every state path by
# enumeration, so it can be checked by eye against the formulas below.
#
# Model, for K ordered states and a design row x (previous state, time and
# covariates already encoded in x):
#   c_k = plogis(alpha_k + x'beta) = P(Y >= k + 1), k = 1..K-1,
#   P(Y = y | x) = c_{y-1} - c_y, with c_0 = 1 and c_K = 0.
# theta = c(alpha, beta).
#
# Log-likelihood = sum_i w_i log P(y_i | x_i)                 (observed rows)
#                + sum_s weight_s log f_s                      (stretches)
# where a stretch s has steps t = 1..T with the state z_t known to lie in
# lo_t..hi_t, and
#   f_s = sum over paths (z_1, ..., z_T), z_t in lo_t..hi_t, of
#         P(z_1 | first row) * prod_{t >= 2} P(z_t | z_{t-1}),
# with P(z_t | z_{t-1}) = 1(z_t == z_{t-1}) when z_{t-1} is absorbing, and
# otherwise the cell probability for the design row step_rows[u, z_{t-1}] of
# step u (the step's position in the step arrays).

path_cell_probability <- function(alpha, beta, x, y) {
  K <- length(alpha) + 1
  eta <- sum(x * beta)
  cumulative <- c(1, stats::plogis(alpha + eta), 0) # c_0, c_1..c_{K-1}, c_K
  cumulative[y] - cumulative[y + 1]
}

path_loglik_bruteforce <- function(theta, K, plan) {
  alpha <- theta[seq_len(K - 1)]
  beta <- theta[-seq_len(K - 1)]

  obs_group_loglik <- numeric(plan$n_groups)
  for (i in seq_len(nrow(plan$x_obs))) {
    prob <- path_cell_probability(alpha, beta, plan$x_obs[i, ], plan$y_obs[i])
    g <- plan$group_obs[i]
    obs_group_loglik[g] <- obs_group_loglik[g] + plan$w_obs[i] * log(prob)
  }

  # gap_prob[r, y] = P(Y = y | x_gap[r, ]), computed once per row because
  # many paths share a row.
  gap_prob <- matrix(0, nrow(plan$x_gap), K)
  for (r in seq_len(nrow(plan$x_gap))) {
    for (y in seq_len(K)) {
      gap_prob[r, y] <- path_cell_probability(alpha, beta, plan$x_gap[r, ], y)
    }
  }

  m <- length(plan$stretch_first_row)
  stretch_loglik <- numeric(m)
  for (s in seq_len(m)) {
    steps <- plan$stretch_first_step[s] + seq_len(plan$stretch_n_steps[s]) - 1
    ranges <- lapply(steps, function(u) plan$step_lo[u]:plan$step_hi[u])
    paths <- as.matrix(expand.grid(ranges))
    f <- 0
    for (r in seq_len(nrow(paths))) {
      z <- paths[r, ]
      prob <- gap_prob[plan$stretch_first_row[s], z[1]]
      for (t in seq_along(steps)[-1]) {
        previous <- z[t - 1]
        if (previous %in% plan$absorb) {
          prob <- prob * (z[t] == previous)
        } else {
          row <- plan$step_rows[steps[t], previous]
          prob <- prob * gap_prob[row, z[t]]
        }
      }
      f <- f + prob
    }
    stretch_loglik[s] <- log(f)
  }

  list(
    loglik = sum(obs_group_loglik) + sum(plan$stretch_weight * stretch_loglik),
    stretch_loglik = stretch_loglik,
    obs_group_loglik = obs_group_loglik
  )
}

# Calls the kernel with the pieces of a plan.
path_kernel <- function(theta, K, plan, order = 2L) {
  cpp_path_loglik(
    alpha = theta[seq_len(K - 1)],
    beta = theta[-seq_len(K - 1)],
    x_obs = plan$x_obs,
    y_obs = plan$y_obs,
    w_obs = plan$w_obs,
    group_obs = plan$group_obs,
    n_groups = plan$n_groups,
    x_gap = plan$x_gap,
    absorb = plan$absorb,
    stretch_first_row = plan$stretch_first_row,
    stretch_first_step = plan$stretch_first_step,
    stretch_n_steps = plan$stretch_n_steps,
    stretch_weight = plan$stretch_weight,
    step_lo = plan$step_lo,
    step_hi = plan$step_hi,
    step_rows = plan$step_rows,
    order = order
  )
}

# Deterministic test plans. The design row for previous state `prev` at time
# `time` with covariate `z` is (1(prev == 2), ..., 1(prev == K), time / 4, z),
# so p = K + 1. Rows depend only on (prev, time, z), so stretches with the
# same start time and covariate share their later-step rows in x_gap.
path_design_row <- function(K, prev, time, z) {
  c(as.numeric(prev == 2:K), time / 4, z)
}

path_theta <- function(K) {
  alpha <- seq(1.2, -1.4, length.out = K - 1)
  beta <- c(seq(0.4, 1.1, length.out = K - 1), -0.3, 0.6)
  c(alpha, beta)
}

# `obs`: data frame with columns prev, time, z, y, w, group (may have 0 rows).
# `stretches`: list of list(prev, time, z, lo, hi, weight); step t of a
# stretch is at time `time + t - 1`, and step 1's previous state is `prev`.
build_path_plan <- function(K, absorb, obs, stretches, n_groups = 3L) {
  p <- K + 1
  x_obs <- matrix(0, nrow(obs), p)
  for (i in seq_len(nrow(obs))) {
    x_obs[i, ] <- path_design_row(K, obs$prev[i], obs$time[i], obs$z[i])
  }

  rows <- list()
  row_index <- function(x) {
    key <- paste(format(x, digits = 17), collapse = "|")
    if (is.null(rows[[key]])) {
      rows[[key]] <<- list(index = length(rows) + 1L, x = x)
    }
    rows[[key]]$index
  }
  first_row <- first_step <- n_steps <- integer(0)
  weight <- numeric(0)
  step_lo <- step_hi <- integer(0)
  step_rows <- matrix(0L, 0, K)
  for (st in stretches) {
    T_s <- length(st$lo)
    first_row <- c(first_row, row_index(path_design_row(K, st$prev, st$time, st$z)))
    first_step <- c(first_step, length(step_lo) + 1L)
    n_steps <- c(n_steps, T_s)
    weight <- c(weight, st$weight)
    for (t in seq_len(T_s)) {
      step_lo <- c(step_lo, as.integer(st$lo[t]))
      step_hi <- c(step_hi, as.integer(st$hi[t]))
      used <- integer(K)
      if (t >= 2) {
        for (j in setdiff(seq_len(K), absorb)) {
          used[j] <- row_index(path_design_row(K, j, st$time + t - 1, st$z))
        }
      }
      step_rows <- rbind(step_rows, used)
    }
  }
  x_gap <- matrix(0, length(rows), p)
  for (r in rows) x_gap[r$index, ] <- r$x
  dimnames(step_rows) <- NULL

  list(
    x_obs = x_obs,
    y_obs = as.integer(obs$y),
    w_obs = as.numeric(obs$w),
    group_obs = as.integer(obs$group),
    n_groups = as.integer(n_groups),
    x_gap = x_gap,
    absorb = as.integer(absorb),
    stretch_first_row = first_row,
    stretch_first_step = first_step,
    stretch_n_steps = n_steps,
    stretch_weight = weight,
    step_lo = step_lo,
    step_hi = step_hi,
    step_rows = step_rows
  )
}

path_observed_rows <- function(K) {
  data.frame(
    prev = c(1, 2, K - 1, 2, 1, K - 1),
    time = c(1, 2, 3, 4, 2, 5),
    z = c(0.5, -1, 0, 1.5, 0.2, -0.4),
    y = c(1, 2, K, K - 1, 2, 1),
    w = c(1, 2.5, 0.5, 1, 3, 0.7),
    group = c(1, 1, 2, 2, 3, 3)
  )
}

# Stretches covering: closed T = 2, 3, 4; closed ending in the highest state;
# open T = 1, 2, 3; interior ranges; a range at step 1 then an exact state; a
# single-step range (an interval-censored row); ranges that include the
# highest state at an interior step; a final full range [1, K]; all steps
# exact; two stretches sharing x_gap rows. `top` is the highest state allowed
# at an unobserved visit: K - 1 when K is absorbing, else K.
path_stretches <- function(K, top) {
  s <- function(prev, time, z, lo, hi, weight = 1) {
    list(prev = prev, time = time, z = z, lo = lo, hi = hi, weight = weight)
  }
  list(
    s(2, 1, 0.3, c(1, 3), c(top, 3)),
    s(2, 1, 0.3, c(1, 2), c(top, 2), weight = 2.5), # shares rows with the first
    s(1, 2, -0.5, c(1, 1, K), c(top, top, K), weight = 0.5),
    s(K - 1, 1, 1, c(1, 1, 1, 2), c(top, top, top, 2), weight = 3),
    s(1, 3, 0, 1, top),
    s(2, 2, -1, c(1, 1), c(top, top), weight = 1.5),
    s(K - 1, 1, 0.8, c(1, 1, 1), c(top, top, top)),
    s(2, 2, 0.1, c(2, K - 1), c(3, K - 1)),
    s(1, 1, -0.2, c(2, 1), c(3, top), weight = 2),
    s(2, 3, 0.4, 2, 3, weight = 0.8),
    s(1, 2, 0.6, c(3, K), c(K, K)),
    s(2, 1, -0.7, c(3, K - 1), c(K, K)),
    s(1, 1, 0.9, c(1, 1), c(top, K)),
    s(K - 1, 2, -0.3, c(2, 1, 2), c(2, 1, 2), weight = 1.2)
  )
}

path_plans <- function() {
  plans <- list()
  for (K in c(3L, 5L)) {
    for (absorbing in c(TRUE, FALSE)) {
      absorb <- if (absorbing) K else integer(0)
      top <- if (absorbing) K - 1L else K
      name <- sprintf("K = %d, %s", K, if (absorbing) "absorbing K" else "no absorbing state")
      plans[[name]] <- list(
        K = K,
        plan = build_path_plan(K, absorb, path_observed_rows(K), path_stretches(K, top))
      )
    }
  }
  plans
}
