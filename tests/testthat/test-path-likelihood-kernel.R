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
    theta <- case$theta
    kernel <- path_kernel(theta, case$K, case$plan)
    reference <- path_loglik_bruteforce(theta, case$K, case$plan)
    expect_true(is.finite(kernel$loglik), info = name)
    expect_equal(kernel$loglik, reference$loglik, tolerance = 1e-12, info = name)
    expect_equal(kernel$stretch_loglik, reference$stretch_loglik,
                 tolerance = 1e-12, info = name)
  }
})

test_that("scores equal numerical gradients of the enumeration", {
  skip_on_cran() # Slow (about 2 s); the long-stretch test below compares scores with numerical gradients on CRAN.
  # The interaction layout runs the same code as the general one; its
  # log-likelihood is checked above, and numerical gradients are slow.
  plans <- path_plans(c("general", "additive"))
  for (name in names(plans)) {
    case <- plans[[name]]
    K <- case$K
    plan <- case$plan
    theta <- case$theta
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
  plans <- path_plans(c("general", "additive"))
  for (name in names(plans)) {
    case <- plans[[name]]
    K <- case$K
    plan <- case$plan
    theta <- case$theta
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

test_that("a 2000-step stretch keeps its log-likelihood", {
  # The kernel takes one log per run of forward scales instead of one per
  # step; the run's product must not underflow however long the stretch.
  # Exact states alternate between 4 and 1, so f is the product of
  # P(4 | 1) = plogis(-10) and P(1 | 4) = plogis(-18), 1000 of each
  # (log f about -28000).
  # The plan is written out directly: build_path_plan() makes one design row
  # per step and took about 2 s for 2000 steps. The time and z columns have
  # coefficient 0, so two rows serve every step: x_gap row 1 for previous
  # state 1, row 2 for previous state 4.
  K <- 4L
  T_s <- 2000L
  states <- as.integer(rep(c(4, 1), length.out = T_s))
  step_rows <- matrix(0L, T_s, K)
  step_rows[, 1] <- 1L
  step_rows[, 4] <- 2L
  plan <- list(
    x_obs = matrix(0, 0, 5), y_obs = integer(0), w_obs = numeric(0),
    group_obs = integer(0), n_groups = 1L,
    x_gap = rbind(c(0, 0, 0, 0, 0), c(0, 0, 1, 0, 0)),
    absorb = integer(0),
    stretch_first_row = 1L, stretch_first_prev = 1L,
    stretch_first_step = 1L, stretch_n_steps = T_s, stretch_weight = 1,
    step_lo = states, step_hi = states, step_rows = step_rows
  )
  theta <- c(10, 0, -10, 0, 0, 8, 0, 0)
  kernel <- path_kernel(theta, K, plan)
  expected <- 1000 * (stats::plogis(-10, log.p = TRUE) +
                        stats::plogis(-18, log.p = TRUE))
  expect_equal(kernel$stretch_loglik, expected, tolerance = 1e-13)
  expect_true(all(is.finite(kernel$stretch_scores)))
  expect_true(all(is.finite(kernel$hessian)))
})

test_that("tiny probabilities in a row keep the log-likelihood finite", {
  # The kernel multiplies the probabilities of weight-1 observed rows, and
  # the forward scales of a stretch, and takes one log per run. Here they
  # are about 1e-40, 1e-150 and 1e-260: multiplied naively the product
  # underflows to 0 (1e-450) and the log-likelihood becomes -Inf.
  # Design row (1(prev == 2), time / 4, z) with beta = (0, 1012, -161) and
  # z = 1: eta = 253 time - 161 = 92, 345, 598 at times 1, 2, 3.
  K <- 2L
  theta <- c(0, 0, 1012, -161)
  eta <- 253 * (1:3) - 161
  expected <- sum(stats::plogis(-eta, log.p = TRUE))
  obs <- data.frame(prev = 1, time = 1:3, z = 1, y = 1, w = 1, group = 1)
  as_rows <- path_kernel(theta, K, build_path_plan(K, integer(0), obs, list(), 1L))
  stretch <- list(prev = 1, time = 1, z = 1, lo = c(1, 1, 1), hi = c(1, 1, 1),
                  weight = 1)
  as_stretch <- path_kernel(
    theta, K, build_path_plan(K, integer(0), obs[0, ], list(stretch), 1L)
  )
  expect_equal(as_rows$loglik, expected, tolerance = 1e-14)
  expect_equal(as_stretch$loglik, expected, tolerance = 1e-14)
  expect_equal(as_stretch$stretch_scores, as_rows$obs_group_scores, tolerance = 1e-12)
})

test_that("orders 0, 1 and 2 agree on what they share", {
  case <- path_plans()[["K = 5, absorbing K, additive layout"]]
  theta <- case$theta
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
  plan <- path_plans()[["K = 5, absorbing K, general layout"]]$plan
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

test_that("the additive and the general layout give the same likelihood", {
  # The same model and data, with the previous state's columns either in
  # every design row (general) or added from x_prev (additive). The two
  # differ only in the order of floating-point operations.
  plans <- path_plans(c("general", "additive"))
  for (K in c(3L, 5L)) {
    for (absorbing in c("absorbing K", "no absorbing state")) {
      general <- plans[[sprintf("K = %d, %s, general layout", K, absorbing)]]
      additive <- plans[[sprintf("K = %d, %s, additive layout", K, absorbing)]]
      expect_lt(nrow(additive$plan$x_gap), nrow(general$plan$x_gap))
      theta <- general$theta
      a <- path_kernel(theta, K, general$plan)
      b <- path_kernel(theta, K, additive$plan)
      for (part in names(a)) {
        expect_equal(b[[part]], a[[part]], tolerance = 1e-13, info = part)
      }
    }
  }
})

test_that("an additive plan needs a valid previous state for each stretch", {
  K <- 5L
  plan <- path_plans("additive")[["K = 5, absorbing K, additive layout"]]$plan
  plan$stretch_first_prev[2] <- 0L
  expect_error(path_kernel(path_theta(K), K, plan), "stretch_first_prev")
})
