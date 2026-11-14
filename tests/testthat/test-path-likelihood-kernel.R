# cpp_path_loglik() against the path enumeration in helper-path-likelihood.R
# and against central differences. No random numbers are used.

total_score <- function(result, plan) {
  colSums(result$obs_group_scores) +
    colSums(plan$stretch_weight * result$stretch_scores)
}

test_that("the kernel log-likelihood equals the path enumeration", {
  plans <- path_plans()
  for (name in names(plans)) {
    case <- plans[[name]]
    theta <- path_theta(case$K)
    kernel <- path_kernel(theta, case$K, case$plan)
    reference <- path_loglik_bruteforce(theta, case$K, case$plan)
    expect_true(is.finite(kernel$loglik), info = name)
    expect_equal(kernel$loglik, reference$loglik, tolerance = 1e-12, info = name)
    expect_equal(kernel$stretch_loglik, reference$stretch_loglik,
                 tolerance = 1e-12, info = name)
  }
})

test_that("scores equal numerical gradients of the enumeration", {
  plans <- path_plans()
  for (name in names(plans)) {
    case <- plans[[name]]
    K <- case$K
    plan <- case$plan
    theta <- path_theta(K)
    kernel <- path_kernel(theta, K, plan)
    q <- length(theta)
    expect_equal(dim(kernel$obs_group_scores), c(plan$n_groups, q))
    expect_equal(dim(kernel$stretch_scores), c(length(plan$stretch_weight), q))
    m <- length(plan$stretch_weight)
    numeric <- unname(central_difference_jacobian(
      function(th) {
        reference <- path_loglik_bruteforce(th, K, plan)
        c(reference$stretch_loglik, reference$obs_group_loglik)
      },
      theta
    ))
    expect_equal(kernel$stretch_scores, numeric[seq_len(m), , drop = FALSE],
                 tolerance = 1e-7, info = name)
    expect_equal(kernel$obs_group_scores, numeric[-seq_len(m), , drop = FALSE],
                 tolerance = 1e-7, info = name)
  }
})

test_that("the Hessian is the derivative of the analytical score", {
  plans <- path_plans()
  for (name in names(plans)) {
    case <- plans[[name]]
    K <- case$K
    plan <- case$plan
    theta <- path_theta(K)
    kernel <- path_kernel(theta, K, plan)
    numeric_hessian <- central_difference_jacobian(
      function(th) total_score(path_kernel(th, K, plan, order = 1L), plan), theta
    )
    expect_equal(kernel$hessian, t(kernel$hessian), tolerance = 1e-14, info = name)
    expect_equal(kernel$hessian, unname(numeric_hessian),
                 tolerance = 1e-6, info = name)
  }
})

test_that("long stretches with tiny probabilities keep exact derivatives", {
  # The kernel rescales the forward vector at every step. Without that, f
  # below would be about 1e-40 in the first case and would underflow to 0
  # (log f < log(.Machine$double.xmin)) in the second.
  K <- 4L
  T_s <- 60L
  narrow <- list(prev = 1, time = 1, z = 0.2,
                 lo = rep(c(1, 3), length.out = T_s), hi = rep(c(1, 4), length.out = T_s),
                 weight = 1)
  narrow$lo[T_s] <- narrow$hi[T_s] <- 3
  alternating <- list(prev = 1, time = 1, z = 0,
                      lo = rep(c(4, 1), length.out = T_s), hi = rep(c(4, 1), length.out = T_s),
                      weight = 2)
  cases <- list(
    list(stretch = narrow, theta = path_theta(K), below = log(1e-40)),
    # P(4 | 1) = plogis(-10) and P(1 | 4) = plogis(-18): about 1e-364 in all.
    list(stretch = alternating, theta = c(10, 0, -10, 0, 0, 8, 0, 0),
         below = log(.Machine$double.xmin))
  )
  for (case in cases) {
    plan <- build_path_plan(K, integer(0), path_observed_rows(K)[0, ],
                            list(case$stretch), n_groups = 1L)
    theta <- case$theta
    kernel <- path_kernel(theta, K, plan)
    expect_true(is.finite(kernel$loglik))
    expect_lt(kernel$stretch_loglik, case$below)
    numeric_score <- central_difference_jacobian(
      function(th) path_kernel(th, K, plan, order = 0L)$stretch_loglik, theta
    )
    expect_equal(kernel$stretch_scores, unname(numeric_score), tolerance = 1e-7)
    numeric_hessian <- central_difference_jacobian(
      function(th) total_score(path_kernel(th, K, plan, order = 1L), plan), theta
    )
    expect_equal(kernel$hessian, unname(numeric_hessian), tolerance = 1e-6)
  }
})

test_that("orders 0, 1 and 2 agree on what they share", {
  case <- path_plans()[["K = 5, absorbing K"]]
  theta <- path_theta(case$K)
  r0 <- path_kernel(theta, case$K, case$plan, order = 0L)
  r1 <- path_kernel(theta, case$K, case$plan, order = 1L)
  r2 <- path_kernel(theta, case$K, case$plan, order = 2L)
  expect_named(r0, c("loglik", "stretch_loglik"))
  expect_named(r1, c("loglik", "stretch_loglik", "obs_group_scores", "stretch_scores"))
  expect_named(r2, c(names(r1), "hessian"))
  expect_identical(r0$loglik, r1$loglik)
  expect_identical(r0$loglik, r2$loglik)
  expect_identical(r0$stretch_loglik, r2$stretch_loglik)
  expect_identical(r1$obs_group_scores, r2$obs_group_scores)
  expect_identical(r1$stretch_scores, r2$stretch_scores)
})

test_that("an exact single-step stretch equals the same observed transition", {
  K <- 5L
  obs <- data.frame(prev = 2, time = 3, z = 0.4, y = 4, w = 1, group = 1)
  as_row <- build_path_plan(K, K, obs, list(), n_groups = 1L)
  as_stretch <- build_path_plan(
    K, K, obs[0, ],
    list(list(prev = 2, time = 3, z = 0.4, lo = 4, hi = 4, weight = 1)),
    n_groups = 1L
  )
  theta <- path_theta(K)
  row <- path_kernel(theta, K, as_row)
  stretch <- path_kernel(theta, K, as_stretch)
  expect_equal(stretch$loglik, row$loglik, tolerance = 1e-14)
  expect_equal(stretch$stretch_scores, row$obs_group_scores, tolerance = 1e-14)
  expect_equal(stretch$hessian, row$hessian, tolerance = 1e-14)
  expect_equal(row$loglik, path_loglik_bruteforce(theta, K, as_row)$loglik,
               tolerance = 1e-14)
})

test_that("a stretch of exact states equals the product of observed transitions", {
  K <- 5L
  stretch <- list(prev = 3, time = 2, z = -0.3, lo = c(2, 3, 1), hi = c(2, 3, 1),
                  weight = 1)
  obs <- data.frame(prev = c(3, 2, 3), time = 2:4, z = -0.3, y = c(2, 3, 1),
                    w = 1, group = 1)
  theta <- path_theta(K)
  as_stretch <- path_kernel(theta, K, build_path_plan(K, K, obs[0, ], list(stretch), 1L))
  as_rows <- path_kernel(theta, K, build_path_plan(K, K, obs, list(), 1L))
  expect_equal(as_stretch$loglik, as_rows$loglik, tolerance = 1e-13)
  expect_equal(as_stretch$stretch_scores, as_rows$obs_group_scores, tolerance = 1e-13)
  expect_equal(as_stretch$hessian, as_rows$hessian, tolerance = 1e-13)
})

test_that("a single-step range is the log of the summed cell probabilities", {
  K <- 5L
  plan <- build_path_plan(
    K, integer(0), path_observed_rows(K)[0, ],
    list(list(prev = 2, time = 1, z = 0.5, lo = 2, hi = 3, weight = 1)),
    n_groups = 1L
  )
  theta <- path_theta(K)
  x <- path_design_row(K, 2, 1, 0.5)
  expected <- log(
    path_cell_probability(theta[1:4], theta[-(1:4)], x, 2) +
      path_cell_probability(theta[1:4], theta[-(1:4)], x, 3)
  )
  expect_equal(path_kernel(theta, K, plan)$loglik, expected, tolerance = 1e-14)
})

test_that("a full range [1, K] without absorbing states adds nothing", {
  K <- 5L
  s <- function(lo, hi) list(prev = 2, time = 1, z = 0.3, lo = lo, hi = hi, weight = 2)
  theta <- path_theta(K)
  empty_obs <- path_observed_rows(K)[0, ]
  open_one <- path_kernel(theta, K, build_path_plan(K, integer(0), empty_obs, list(s(1, K)), 1L))
  expect_equal(open_one$loglik, 0, tolerance = 1e-15)
  expect_equal(open_one$stretch_scores, matrix(0, 1, length(theta)), tolerance = 1e-15)
  expect_equal(open_one$hessian, matrix(0, length(theta), length(theta)), tolerance = 1e-15)

  with_last <- path_kernel(theta, K, build_path_plan(K, integer(0), empty_obs, list(s(c(2, 1), c(3, K))), 1L))
  without <- path_kernel(theta, K, build_path_plan(K, integer(0), empty_obs, list(s(2, 3)), 1L))
  expect_equal(with_last$loglik, without$loglik, tolerance = 1e-14)
  expect_equal(with_last$stretch_scores, without$stretch_scores, tolerance = 1e-13)
  expect_equal(with_last$hessian, without$hessian, tolerance = 1e-13)
})

test_that("unobserved visits exclude the absorbing state", {
  K <- 5L
  s <- function(top) {
    list(prev = 1, time = 2, z = -0.5, lo = c(1, 1, K), hi = c(top, top, K), weight = 1)
  }
  theta <- path_theta(K)
  empty_obs <- path_observed_rows(K)[0, ]
  restricted <- build_path_plan(K, K, empty_obs, list(s(K - 1)), 1L)
  # Allowing state K at the interior visits adds the paths that die early and
  # then stay dead by the identity transition.
  unrestricted <- build_path_plan(K, K, empty_obs, list(s(K)), 1L)
  f_restricted <- path_kernel(theta, K, restricted)$loglik
  f_unrestricted <- path_kernel(theta, K, unrestricted)$loglik
  expect_equal(f_restricted, path_loglik_bruteforce(theta, K, restricted)$loglik,
               tolerance = 1e-12)
  expect_equal(f_unrestricted, path_loglik_bruteforce(theta, K, unrestricted)$loglik,
               tolerance = 1e-12)
  expect_lt(f_restricted, f_unrestricted)
})

test_that("plans without observed rows or without stretches work", {
  K <- 3L
  theta <- path_theta(K)
  no_obs <- build_path_plan(K, K, path_observed_rows(K)[0, ], path_stretches(K, K - 1L))
  no_stretch <- build_path_plan(K, K, path_observed_rows(K), list())
  for (plan in list(no_obs, no_stretch)) {
    kernel <- path_kernel(theta, K, plan)
    expect_equal(kernel$loglik, path_loglik_bruteforce(theta, K, plan)$loglik,
                 tolerance = 1e-12)
    expect_true(isSymmetric(kernel$hessian))
  }
  expect_equal(dim(path_kernel(theta, K, no_stretch)$stretch_scores), c(0L, length(theta)))
  expect_equal(nrow(no_stretch$x_gap), 0L)
})

test_that("an impossible observation gives -Inf, not NaN", {
  K <- 5L
  plan <- build_path_plan(
    K, K, path_observed_rows(K)[0, ],
    list(list(prev = 1, time = 1, z = 0, lo = c(K, 1), hi = c(K, 1), weight = 1)),
    n_groups = 1L
  )
  result <- path_kernel(path_theta(K), K, plan)
  expect_identical(result$loglik, -Inf)
  expect_identical(result$stretch_loglik, -Inf)
})

test_that("malformed plans are refused", {
  K <- 5L
  plan <- path_plans()[["K = 5, absorbing K"]]$plan
  theta <- path_theta(K)
  bad <- plan
  bad$step_rows[2, 1] <- 0L # needed: step 2 of stretch 1 starts in 1..4
  expect_error(path_kernel(theta, K, bad), "step_rows")
  bad <- plan
  bad$step_hi[1] <- K + 1L
  expect_error(path_kernel(theta, K, bad), "step_lo")
  bad <- plan
  bad$group_obs[1] <- 4L
  expect_error(path_kernel(theta, K, bad), "group_obs")
})
