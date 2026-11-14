# The maximisation of the likelihood with unobserved states
# (`markov_newton_maximise()`, `markov_ascent_direction()`): steps where the
# log-likelihood is not concave, the stopping rule, and stretches of weight 0.

# A one-coefficient objective for `markov_newton_maximise()`: `f`, its
# derivative `d1` and second derivative `d2` as functions of theta.
ms_objective <- function(f, d1, d2) {
  list(
    evaluate = function(theta, order) {
      out <- list(loglik = f(theta), score = d1(theta))
      if (order >= 2L) out$hessian <- matrix(d2(theta), 1L, 1L)
      out
    },
    score_of = function(result) result$score
  )
}

test_that("an indefinite information matrix still gives an ascent direction", {
  # The Newton step solve(information, score) = (0.5, -1) points downhill
  # here: score' step = -0.5.
  information <- diag(c(2, -1))
  score <- c(1, 1)
  expect_lt(sum(score * solve(information, score)), 0)
  direction <- markov_ascent_direction(information, score)
  expect_true(direction$modified)
  expect_gt(sum(score * direction$step), 0)
  # A positive definite information gives the Newton step itself.
  information <- matrix(c(4, 1, 1, 3), 2L, 2L)
  direction <- markov_ascent_direction(information, score)
  expect_false(direction$modified)
  expect_equal(direction$step, solve(information, score), tolerance = 1e-14)
  # Nothing can be formed from a score that is not finite.
  expect_null(markov_ascent_direction(information, c(1, NaN)))
})

test_that("Newton iterations climb out of a region that is not concave", {
  # f = -(theta^2 - 1)^2 has its maxima at -1 and 1 and f'' > 0 at 0.1, where
  # the plain Newton step goes downhill, towards the minimum at 0.
  objective <- ms_objective(
    function(x) -(x^2 - 1)^2,
    function(x) -4 * x * (x^2 - 1),
    function(x) -(12 * x^2 - 4)
  )
  maximum <- markov_newton_maximise(0.1, objective$evaluate, objective$score_of)
  expect_true(maximum$converged)
  expect_equal(maximum$theta, 1, tolerance = 1e-5)
  # The Hessian returned is the exact, unmodified one at the maximum.
  expect_identical(
    maximum$result$hessian[1, 1],
    -(12 * maximum$theta^2 - 4)
  )
})

test_that("the stopping rule allows for rounding in large log-likelihoods", {
  # The score carries an error of 1e-3 that pushes away from the maximum at
  # 1, as rounding in a sum of many terms can, so the iterates alternate
  # between 1 - 1e-3 / h and 1 + 1e-3 / h, where the log-likelihood is the
  # same, and the Newton decrement stays at (2e-3)^2 / (2 * h) = 1.6e-10 (the
  # value measured with 2,000 patients and a log-likelihood near -2e4).
  h <- 1.25e4
  objective <- function(level) {
    ms_objective(
      function(x) level - h / 2 * (x - 1)^2,
      function(x) -h * (x - 1) + if (x > 1) -1e-3 else 1e-3,
      function(x) -h
    )
  }
  # At a log-likelihood of -2e4 that is rounding error: the fit stops.
  large <- objective(-2e4)
  maximum <- markov_newton_maximise(0, large$evaluate, large$score_of)
  expect_true(maximum$converged)
  expect_lte(maximum$iterations, 3L)
  expect_equal(maximum$theta, 1, tolerance = 1e-6)
  # At a log-likelihood of -1 a decrement of 1.6e-10 is not converged.
  small <- objective(-1)
  maximum <- markov_newton_maximise(0, small$evaluate, small$score_of, maxit = 10L)
  expect_false(maximum$converged)
})

test_that("a likelihood without a finite maximum does not converge", {
  # f = theta rises without end and its Hessian is 0: no step can be formed.
  objective <- ms_objective(function(x) x, function(x) 1, function(x) 0)
  maximum <- markov_newton_maximise(0, objective$evaluate, objective$score_of)
  expect_false(maximum$converged)
  # f = theta + theta^3 / 3 rises without end and is convex for theta > 0:
  # every step needs the modified information, so the iteration never stops
  # as converged.
  objective <- ms_objective(
    function(x) x + x^3 / 3,
    function(x) 1 + x^2,
    function(x) 2 * x
  )
  maximum <- markov_newton_maximise(1, objective$evaluate, objective$score_of)
  expect_false(maximum$converged)
})

test_that("a stretch of weight 0 adds nothing, even when it is impossible", {
  K <- 5L
  stretches <- path_stretches(K, K - 1L)
  # From the absorbing state K the patient can't move to state 1.
  impossible <- list(prev = 1, time = 1, z = 0, lo = c(K, 1), hi = c(K, 1),
                     weight = 0)
  with_zero <- build_path_plan(K, K, path_observed_rows(K),
                               c(stretches, list(impossible)))
  without <- build_path_plan(K, K, path_observed_rows(K), stretches)
  with_zero$labels <- without$labels <- as.character(seq_len(K))
  theta <- path_theta(K)
  result <- path_kernel(theta, K, with_zero)
  expect_true(all(is.nan(result$stretch_scores[length(stretches) + 1L, ])))
  expect_identical(
    markov_path_total_score(result, with_zero),
    markov_path_total_score(path_kernel(theta, K, without), without)
  )
  # Newton steps from there are the same with and without it.
  observed <- list(
    x = with_zero$x_obs, y = with_zero$y_obs, weight = with_zero$w_obs,
    group = with_zero$group_obs, n_groups = with_zero$n_groups
  )
  steps <- function(plan) {
    markov_newton_maximise(
      theta,
      function(theta, order) markov_path_loglik(theta, observed, plan, order),
      function(result) markov_path_total_score(result, plan),
      maxit = 3L
    )
  }
  fit_with <- steps(with_zero)
  fit_without <- steps(without)
  expect_identical(fit_with$iterations, 3L)
  expect_identical(fit_with$theta, fit_without$theta)
  # The patient scores of weight-0 patients are 0, not NaN.
  with_zero$instances <- data.frame(
    stretch = seq_along(with_zero$stretch_weight),
    id = c(rep("1", length(stretches)), "2"),
    weight = with_zero$stretch_weight
  )
  scores <- markov_path_patient_scores(fit_with$result, with_zero, c("1", "2", "3"))
  expect_true(all(is.finite(scores)))
})

test_that("a start where the likelihood is not concave reaches the maximum", {
  skip_on_cran() # The pass depends on the seed.
  # 30 patients, 8 states, 70% of states missing: at the starting values the
  # information has an eigenvalue of -12.7, and the plain Newton step lowered
  # the log-likelihood however much it was shortened.
  data <- ms_simulate_many_states(30, 20, K = 8, miss = 0.7, seed = 4)
  fit <- mostm(
    ordered(y) ~ yprev + time + tx + age,
    data = data,
    id_var = "id",
    absorb = 8
  )
  # A maximum: score 0 and a positive definite information matrix (BFGS from
  # the starting values reached the same point within 2e-5).
  info <- attr(fit, "markov_path_fit")
  expect_lt(max(abs(fit$u)), 1e-5)
  expect_gt(min(eigen(info$information, symmetric = TRUE)$values), 0.01)
  expect_lt(max(abs(coef(fit))), 10)
})
