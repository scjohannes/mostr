# The scale on which the Stan program for blrm fits with unobserved states
# samples the intercepts and coefficients (`markov_path_stan_scale()`), and
# the arguments given to the sampler, without compiling or sampling. The
# sampled density is checked in test-slow-missing-states-blrm.R.

test_that("the log-odds of the states and their derivatives are right", {
  # u[j] = log(pi[j] / pi[k]) from intercepts alpha[j] = logit(P(Y >= j + 1)),
  # against the state probabilities written out directly, also with a last
  # state of probability 1e-7, as at the centred predictors of a model whose
  # worst state is rare.
  for (pi in list(c(0.4, 0.3, 0.2, 0.1), c(0.6, 0.25, 0.15 - 1e-7, 1e-7))) {
    alpha <- stats::qlogis(rev(cumsum(rev(pi)))[-1])
    log_odds <- mostr:::markov_path_log_odds(alpha)
    expect_equal(log_odds$u, log(pi[-4] / pi[4]), tolerance = 1e-10)
    # The derivatives against central differences.
    numeric_jacobian <- central_difference_jacobian(
      function(a) mostr:::markov_path_log_odds(a)$u,
      alpha
    )
    expect_equal(log_odds$jacobian, unname(numeric_jacobian), tolerance = 1e-6)
  }
  # Intercepts that do not decrease give no probabilities.
  expect_false(all(is.finite(mostr:::markov_path_log_odds(c(0, 1))$u)))
})

test_that("the sampler gets adapt_delta 0.9 with random intercepts unless the user sets it", {
  build <- function(backend, sampling_args = NULL, random_intercepts = TRUE) {
    mostr:::markov_path_sampler_args(
      backend,
      stan_data = list(k = 3L),
      pars = c("alpha", "beta"),
      iter = 200,
      warmup = 100,
      chains = 2,
      refresh = 0,
      sampling_args = sampling_args,
      random_intercepts = random_intercepts
    )
  }
  rstan <- build("rstan")
  expect_identical(rstan$control, list(adapt_delta = 0.9))
  expect_identical(rstan$iter, 200)
  expect_identical(rstan$warmup, 100)
  expect_identical(rstan$pars, c("alpha", "beta"))
  cmdstan <- build("cmdstan")
  expect_identical(cmdstan$adapt_delta, 0.9)
  expect_identical(cmdstan$iter_sampling, 100)
  expect_identical(cmdstan$iter_warmup, 100)
  expect_null(cmdstan$control)

  # The user's value wins, and the rest of the user's `control` is kept.
  user <- build("rstan", list(control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 1))
  expect_identical(user$control, list(adapt_delta = 0.95, max_treedepth = 12))
  expect_identical(user$seed, 1)
  other <- build("rstan", list(control = list(max_treedepth = 12)))
  expect_identical(other$control, list(max_treedepth = 12, adapt_delta = 0.9))
  expect_identical(build("cmdstan", list(adapt_delta = 0.8))$adapt_delta, 0.8)

  # Without random intercepts Stan's default stays.
  expect_null(build("rstan", random_intercepts = FALSE)$control)
  expect_null(build("cmdstan", random_intercepts = FALSE)$adapt_delta)
  kept <- build("rstan", list(control = list(max_treedepth = 12)), random_intercepts = FALSE)
  expect_identical(kept$control, list(max_treedepth = 12))
})

# The Stan data of the blrm fit for the gap data, built from the orm path fit
# to the same data: its fully observed transitions and a plan of the
# stretches with the same design (yprev = 2, yprev = 3, time, tx). The
# predictors are transformed as in a blrm fit (centred, orthogonalised,
# scaled), so the maximum likelihood fit on the Stan data's scale must be the
# orm fit moved to that scale.
scale_inputs <- function() {
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  labels <- as.character(1:4)
  ranges <- ms_row_ranges(data)
  ranges$lo[is.na(data$y)] <- NA
  ranges$hi[is.na(data$y)] <- NA
  plan <- mostr:::markov_path_plan(
    data = data,
    ranges = list(labels = labels, lo = ranges$lo, hi = ranges$hi),
    absorb = 4L,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    weights = NULL,
    design = function(rows) {
      previous <- as.integer(as.character(rows$yprev))
      cbind(previous == 2, previous == 3, rows$time, rows$tx) * 1
    }
  )
  assembled <- mostr:::markov_blrm_path_stan_data(
    x_setup = fit$x,
    keep = rep(TRUE, nrow(fit$x)),
    observed_state = match(as.character(fit$y), labels),
    plan = plan,
    n_states = 4L,
    conc = 0.5
  )
  list(fit = fit, assembled = assembled)
}

# (alpha, beta) on the original scale -> (alpha, beta) on the Stan data's
# scale: intercepts at the centred predictors, beta_s = rinv^-1 beta. A
# linear map, returned as a matrix.
scale_map <- function(xbar, rinv) {
  n_intercepts <- 3L
  p <- length(xbar)
  out <- matrix(0, n_intercepts + p, n_intercepts + p)
  out[1:n_intercepts, 1:n_intercepts] <- diag(n_intercepts)
  out[1:n_intercepts, n_intercepts + 1:p] <- matrix(xbar, n_intercepts, p, byrow = TRUE)
  out[n_intercepts + 1:p, n_intercepts + 1:p] <- solve(rinv)
  out
}

test_that("the sampling scale is the maximum likelihood fit and its covariance", {
  skip_on_cran() # Slow (about 2 s): loads rmsb for `selectedQr()`.
  skip_if_not_installed("rmsb")
  inputs <- scale_inputs()
  fit <- inputs$fit
  data <- inputs$assembled$data
  map <- scale_map(inputs$assembled$transform$xbar, inputs$assembled$rinv)
  # Newton's method from the state frequencies, a start that needs no fit.
  scale <- mostr:::markov_path_stan_scale(
    data,
    starts = list(state_frequencies = mostr:::markov_path_stan_marginal_start(data))
  )
  expect_true(scale$maximum_likelihood)
  expect_identical(scale$start, "state_frequencies")
  expect_gt(scale$iterations, 0L)

  theta <- as.numeric(map %*% stats::coef(fit))
  log_odds <- mostr:::markov_path_log_odds(theta[1:3])
  expect_equal(scale$theta_hat, c(log_odds$u, theta[-(1:3)]), tolerance = 1e-6)
  # The covariance of (u, beta): the inverse information of the orm fit
  # moved to the Stan data's scale and then to the log-odds.
  covariance <- map %*% solve(attr(fit, "markov_path_fit")$information) %*% t(map)
  jacobian <- diag(7)
  jacobian[1:3, 1:3] <- log_odds$jacobian
  expected <- jacobian %*% covariance %*% t(jacobian)
  expect_equal(scale$theta_L %*% t(scale$theta_L), expected, tolerance = 1e-6)
  expect_true(all(scale$theta_L[upper.tri(scale$theta_L)] == 0))

  # The structure fit's draws give a start on the same scale: draws equal to
  # the orm estimates are the maximum, which Newton's method then does not
  # move.
  draws <- matrix(stats::coef(fit), 3, 7, byrow = TRUE)
  start <- mostr:::markov_path_stan_start_from_draws(
    draws, 4L, inputs$assembled$rinv, inputs$assembled$transform$xbar
  )
  expect_equal(c(start$alpha, start$beta), theta, tolerance = 1e-10)
  from_draws <- mostr:::markov_path_stan_scale(data, starts = list(structure_fit = start))
  expect_identical(from_draws$start, "structure_fit")
  expect_equal(from_draws$theta_hat, scale$theta_hat, tolerance = 1e-6)
  # Draws with other columns give no start.
  expect_null(mostr:::markov_path_stan_start_from_draws(
    draws[, -1], 4L, inputs$assembled$rinv, inputs$assembled$transform$xbar
  ))
})

test_that("without a maximum likelihood fit the program gets a rough but valid scale", {
  skip_on_cran() # Slow (about 2 s): loads rmsb for `selectedQr()`.
  skip_if_not_installed("rmsb")
  data <- scale_inputs()$assembled$data
  start <- mostr:::markov_path_stan_marginal_start(data)
  # No Newton iteration allowed: the maximisation fails.
  rough <- mostr:::markov_path_stan_scale(
    data,
    starts = list(state_frequencies = start),
    maxit = 0L
  )
  expect_false(rough$maximum_likelihood)
  expect_identical(rough$theta_L, diag(7))
  expect_equal(
    rough$theta_hat,
    c(mostr:::markov_path_log_odds(start$alpha)$u, start$beta)
  )
  # A start whose intercepts do not decrease can't be used, neither for the
  # maximisation nor for the rough scale; the next start is.
  bad <- list(alpha = c(0, 1, 2), beta = numeric(4))
  skipped <- mostr:::markov_path_stan_scale(
    data,
    starts = list(structure_fit = bad, state_frequencies = start)
  )
  expect_true(skipped$maximum_likelihood)
  expect_identical(skipped$start, "state_frequencies")
  rough_skipped <- mostr:::markov_path_stan_scale(
    data,
    starts = list(structure_fit = bad, state_frequencies = start),
    maxit = 0L
  )
  expect_identical(rough_skipped$start, "state_frequencies")
  # The Stan data the program reads: theta_hat of length k - 1 + p, finite,
  # and a lower triangular theta_L of that size.
  for (scale in list(rough, skipped)) {
    expect_length(scale$theta_hat, data$k - 1L + data$p)
    expect_true(all(is.finite(scale$theta_hat)))
    expect_identical(dim(scale$theta_L), c(7L, 7L))
    expect_true(all(scale$theta_L[upper.tri(scale$theta_L)] == 0))
    expect_true(all(diag(scale$theta_L) > 0))
  }
})

test_that("the start from the state frequencies is finite for a state never observed", {
  # State 4 occurs only after unobserved visits, so no fully observed
  # transition ends in it.
  data <- list(k = 4L, p = 2L, y_obs = c(1L, 2L, 3L), w_obs = c(5, 3, 2))
  start <- mostr:::markov_path_stan_marginal_start(data)
  expect_true(all(is.finite(start$alpha)))
  expect_true(all(diff(start$alpha) < 0))
  expect_identical(start$beta, c(0, 0))
  expect_equal(stats::plogis(start$alpha[3]), 0.5 / 12)
})
