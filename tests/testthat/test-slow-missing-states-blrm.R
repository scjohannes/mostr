# engine = "blrm" with unobserved states samples from mostr's own Stan
# program. Compiling it takes minutes the first time (it is then cached in
# tools::R_user_dir("mostr", "cache")), so this runs only with
# MOSTR_RUN_SLOW=true. The test uses that default cache on purpose: a
# temporary directory would force the compile on every run. It writes there
# what any first blrm fit with unobserved states writes, and the cache code
# removes the files of earlier versions of the program.
#
# The program samples z, with theta = theta_hat + theta_L z, theta the
# log-odds u = log(pi[-k] / pi[k]) of the state probabilities pi at the
# centred predictors followed by the coefficients beta of the transformed
# predictors (`markov_path_stan_scale()`, stored on the fit as
# `attr(fit, "markov_path_fit")$stan_scale`). `ms_stan_z()` gives the z of a
# point (beta, pi). Because the program samples u and not pi, it adds
# sum(log(pi)), the log Jacobian of u -> pi, to rmsb's log density.
ms_stan_z <- function(scale, beta, pi) {
  k <- length(pi)
  theta <- c(log(pi[-k] / pi[k]), beta)
  as.numeric(forwardsolve(scale$theta_L, theta - scale$theta_hat))
}

test_that("blrm with unobserved states agrees with the orm path fit", {
  skip_on_cran() # The pass depends on the seed.
  skip_if_not_installed("rmsb")
  skip_if_not_installed("rstan")
  skip_if_not(
    identical(Sys.getenv("MOSTR_RUN_SLOW"), "true"),
    "Set MOSTR_RUN_SLOW=true to fit blrm models."
  )
  withr::local_options(mc.cores = 2)

  data <- ms_fixture("gaps")
  data <- data[data$id != attr(data, "special")$unknown_start, ]
  fit_orm <- ms_fit(data)
  fit <- suppressMessages(suppressWarnings(ms_fit(
    data,
    engine = "blrm",
    iter = 2000,
    warmup = 1000,
    chains = 2,
    refresh = 0,
    sampling.args = list(seed = 2026)
  )))

  expect_s3_class(fit, "blrm")
  info <- attr(fit, "markov_path_fit")
  expect_false(is.null(info))
  # The sampling scale came from the maximum likelihood fit.
  expect_true(info$stan_scale$maximum_likelihood)
  expect_identical(info$n_stretches, attr(fit_orm, "markov_path_fit")$n_stretches)
  expect_identical(colnames(fit$draws), names(coef(fit_orm)))
  expect_identical(nrow(fit$draws), 2000L)

  # With the default priors, nearly flat for these coefficients, the
  # posterior means are close to the maximum likelihood estimates: within
  # half a posterior SD (the Monte Carlo error of a mean of 2000 draws is
  # about 0.05 SD). A likelihood that ignored the stretches would move the
  # treatment effect by about 0.3 SD here.
  posterior_mean <- colMeans(fit$draws)
  posterior_sd <- apply(fit$draws, 2, stats::sd)
  expect_lt(max(abs(posterior_mean - coef(fit_orm)) / posterior_sd), 0.5)
  expect_equal(unname(posterior_sd), unname(sqrt(diag(fit_orm$orig.var))), tolerance = 0.3)

  average <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:3, n_draws = 50, seed = 1)
  totals <- stats::aggregate(estimate ~ tx + time, average, sum)
  expect_equal(totals$estimate, rep(1, nrow(totals)), tolerance = 1e-8)
  expect_equal(
    average$estimate,
    avg_sops(fit_orm, variables = list(tx = c(0, 1)), times = 1:3)$estimate,
    tolerance = 0.05
  )
})

test_that("without the maximum likelihood scale the fit still samples the posterior", {
  skip_on_cran() # The pass depends on the seed.
  skip_if_not_installed("rmsb")
  skip_if_not_installed("rstan")
  skip_if_not(
    identical(Sys.getenv("MOSTR_RUN_SLOW"), "true"),
    "Set MOSTR_RUN_SLOW=true to fit blrm models."
  )
  withr::local_options(mc.cores = 2)
  # Newton's method gets no iteration, so the pre-fit fails and the program
  # samples on the rough scale (theta_L the identity): slower, but the same
  # posterior as in the first test.
  original <- mostr:::markov_path_stan_scale
  testthat::local_mocked_bindings(
    markov_path_stan_scale = function(data, starts, maxit = 50L) {
      original(data, starts, maxit = 0L)
    },
    .package = "mostr"
  )
  data <- ms_fixture("downstream")
  fit_orm <- ms_fixture("downstream_fit")
  fit <- suppressMessages(suppressWarnings(ms_fit(
    data,
    engine = "blrm",
    iter = 2000,
    warmup = 1000,
    chains = 2,
    refresh = 0,
    sampling.args = list(seed = 2026)
  )))
  scale <- attr(fit, "markov_path_fit")$stan_scale
  expect_false(scale$maximum_likelihood)
  expect_identical(scale$theta_L, diag(7))
  expect_identical(dim(fit$draws), c(2000L, 7L))
  posterior_mean <- colMeans(fit$draws)
  posterior_sd <- apply(fit$draws, 2, stats::sd)
  expect_lt(max(abs(posterior_mean - coef(fit_orm)) / posterior_sd), 0.5)
})

test_that("the Stan program's density and gradient are the path likelihood's", {
  skip_on_cran() # Compiles a Stan program.
  skip_if_not_installed("rmsb")
  skip_if_not_installed("rstan")
  skip_if_not(
    identical(Sys.getenv("MOSTR_RUN_SLOW"), "true"),
    "Set MOSTR_RUN_SLOW=true to fit blrm models."
  )
  # The program computes the likelihood and its gradient in C++
  # (inst/include/mostr). Its log density is compared with the reference in
  # helper-missing-states.R, which lists every path of the data frame, plus
  # the Dirichlet prior and the Jacobian term sum(log(pi)); differences
  # between two points cancel the constants that `log_prob()` leaves out. The
  # gradient is compared with central differences of the log density. The
  # points are fixed, not draws.
  data <- ms_fixture("gaps")
  data <- data[data$id != attr(data, "special")$unknown_start, ]
  fit <- suppressMessages(suppressWarnings(ms_fit(
    data, engine = "blrm", iter = 10, warmup = 5, chains = 1, refresh = 0,
    sampling.args = list(seed = 1)
  )))
  stanfit <- fit$rstan
  reference <- ms_reference(data)
  # The sampler's coefficients are those of the centred, orthogonalised
  # predictors scaled to SD 1 (see `mostm_fit_blrm_paths()`), computed from
  # the set-up rows, which the fit does not store.
  rinv <- attr(fit, "markov_path_fit")$rinv
  expect_identical(dim(rinv), c(4L, 4L))
  k <- 4L
  conc <- fit$conc %||% (1 / (0.8 + 0.35 * k))
  scale <- attr(fit, "markov_path_fit")$stan_scale
  log_posterior <- function(u) {
    par <- rstan::constrain_pars(stanfit, u)
    beta <- as.numeric(rinv %*% par$beta)
    alpha <- par$alpha - sum(beta * fit$xbar)
    reference(c(alpha, beta)) + sum((conc - 1) * log(par$pi)) + sum(log(par$pi))
  }
  point <- function(beta, pi) {
    rstan::unconstrain_pars(
      stanfit,
      list(z = ms_stan_z(scale, beta, pi), gamma_raw = numeric(), sigmag = numeric())
    )
  }
  points <- list(
    point(c(0.1, -0.2, 0.05, 0.3), c(0.4, 0.3, 0.2, 0.1)),
    point(c(-0.3, 0.1, 0.2, -0.1), c(0.25, 0.25, 0.3, 0.2)),
    point(c(0, 0, 0, 0), c(0.1, 0.2, 0.3, 0.4))
  )
  # The point is the one asked for.
  expect_equal(as.numeric(rstan::constrain_pars(stanfit, points[[1]])$pi), c(0.4, 0.3, 0.2, 0.1))
  density <- vapply(points, function(u) {
    rstan::log_prob(stanfit, u, adjust_transform = FALSE)
  }, numeric(1))
  expected <- vapply(points, log_posterior, numeric(1))
  expect_equal(diff(density), diff(expected), tolerance = 1e-10)
  for (u in points) {
    numeric_gradient <- central_difference_jacobian(
      function(v) rstan::log_prob(stanfit, v), u
    )[1, ]
    expect_equal(
      as.numeric(rstan::grad_log_prob(stanfit, u)),
      unname(numeric_gradient),
      tolerance = 1e-6
    )
  }
})

test_that("blrm fits a state recorded only after unobserved visits", {
  skip_on_cran() # The pass depends on the seed.
  skip_if_not_installed("rmsb")
  skip_if_not_installed("rstan")
  skip_if_not(
    identical(Sys.getenv("MOSTR_RUN_SLOW"), "true"),
    "Set MOSTR_RUN_SLOW=true to fit blrm models."
  )
  withr::local_options(mc.cores = 2)
  # No transition with both states observed ends in death; the orm fit to
  # the same data is the reference (see the first test).
  data <- ms_deaths_after_gaps()
  fit_orm <- ms_fit(data)
  fit <- suppressMessages(suppressWarnings(ms_fit(
    data,
    engine = "blrm",
    iter = 2000,
    warmup = 1000,
    chains = 2,
    refresh = 0,
    sampling.args = list(seed = 2026)
  )))
  expect_identical(colnames(fit$draws), names(coef(fit_orm)))
  observed <- !is.na(data$y) & !is.na(data$yprev)
  expect_identical(nrow(fit$x), sum(observed))
  expect_identical(fit$ylevels, c("1", "2", "3", "4"))
  posterior_mean <- colMeans(fit$draws)
  posterior_sd <- apply(fit$draws, 2, stats::sd)
  expect_lt(max(abs(posterior_mean - coef(fit_orm)) / posterior_sd), 0.5)
})

# Random intercepts (`cluster(id)`). The priors are those of rmsb's
# `lrmconppo`: a Dirichlet prior on pi (the state probabilities at the mean
# predictors), gamma_raw ~ N(0, 1) with gamma = sigmag * gamma_raw, and
# sigmag ~ student_t(4, rsdmean, rsdsd) (rmsb's default, psigma = 1,
# rsdmean = 0, rsdsd = 1). rmsb's parameters are beta, pi, gamma_raw and
# sigmag; ours are z (see `ms_stan_z()`), gamma_raw and sigmag.
# `re_unconstrain()` gives the unconstrained point of a point on the
# original scale (intercepts, coefficients of the design columns, gamma_raw,
# sigmag): for rmsb's program without `scale`, for ours with the fit's
# sampling scale. `re_jacobian()` is the term sum(log(pi)) our program adds.
re_state_probabilities <- function(alpha) -diff(c(1, stats::plogis(alpha), 0))
re_unconstrain <- function(stanfit, point, rinv, xbar, tau = FALSE, scale = NULL) {
  centred <- point$alpha + sum(point$beta * xbar)
  beta <- as.numeric(solve(rinv, point$beta))
  pi <- re_state_probabilities(centred)
  pars <- if (is.null(scale)) {
    list(beta = beta, pi = pi)
  } else {
    list(z = ms_stan_z(scale, beta, pi))
  }
  pars <- c(pars, list(gamma_raw = point$gamma_raw, sigmag = array(point$sigmag, 1)))
  # rmsb's program has partial proportional odds coefficients `tau`, empty
  # here, between beta and pi.
  if (tau) pars <- c(pars[1], list(tau = numeric()), pars[-1])
  rstan::unconstrain_pars(stanfit, pars)
}
re_jacobian <- function(point, xbar) {
  sum(log(re_state_probabilities(point$alpha + sum(point$beta * xbar))))
}
re_points <- function(n_clusters, seed) {
  withr::with_seed(seed, lapply(1:3, function(i) list(
    alpha = c(1.4, -0.9, -3.1) + stats::rnorm(3, 0, 0.2),
    beta = c(1.1, 2.4, -0.12, -0.4) + stats::rnorm(4, 0, 0.2),
    gamma_raw = stats::rnorm(n_clusters),
    sigmag = stats::runif(1, 0.3, 1.5)
  )))
}
re_dirichlet <- function(point, xbar, conc) {
  sum((conc - 1) * log(re_state_probabilities(point$alpha + sum(point$beta * xbar))))
}
# The gap fixture without the patient whose first previous state is
# unknown, and with one more patient (the first of `alive_long` without
# gaps) whose visits 1, 3 and 5 are unobserved, so that every transition of
# that patient touches an unobserved visit.
re_data <- function() {
  data <- ms_fixture("downstream")
  special <- attr(data, "special")
  hidden <- attr(data, "hidden")$id
  only_gaps <- setdiff(special$alive_long, c(hidden, special$unknown_start))[1]
  data <- ms_hide(data, rep(only_gaps, 3), c(1, 3, 5))
  attr(data, "only_gaps") <- only_gaps
  data
}
re_skip <- function() {
  skip_if_not_installed("rmsb")
  skip_if_not_installed("rstan")
  skip_if_not(
    identical(Sys.getenv("MOSTR_RUN_SLOW"), "true"),
    "Set MOSTR_RUN_SLOW=true to fit blrm models."
  )
}

test_that("with random intercepts the Stan density is the path likelihood's", {
  skip_on_cran() # Compiles a Stan program.
  re_skip()
  # rms is attached so that `cluster()` is found; a warning that it was
  # built under a newer R says nothing about the test.
  suppressWarnings(withr::local_package("rms"))
  data <- re_data()
  fit <- suppressMessages(suppressWarnings(ms_fit(
    data,
    formula = ordered(y) ~ yprev + time + tx + cluster(id),
    engine = "blrm", iter = 10, warmup = 5, chains = 1, refresh = 0,
    sampling.args = list(seed = 1)
  )))
  stanfit <- fit$rstan
  rinv <- attr(fit, "markov_path_fit")$rinv
  scale <- attr(fit, "markov_path_fit")$stan_scale
  expect_true(scale$maximum_likelihood)
  clusters <- levels(fit$clusterInfo$cluster)
  only_gaps <- as.character(attr(data, "only_gaps"))
  # Every fitted patient has a random intercept, also the one whose every
  # transition touches an unobserved visit and who has no row of `x`.
  expect_true(only_gaps %in% clusters)
  expect_false(only_gaps %in% as.character(fit$clusterInfo$cluster))
  expect_setequal(clusters, attr(fit, "markov_path_fit")$ids)

  # (1) Against the reference that lists every path of the data frame, with
  # each patient's gamma added to every transition, plus the priors and the
  # Jacobian term sum(log(pi)). The difference is a constant, which
  # `log_prob()` leaves out of the priors.
  reference <- ms_reference(data)
  points <- re_points(length(clusters), seed = 11)
  expected <- vapply(points, function(point) {
    gamma <- stats::setNames(point$sigmag * point$gamma_raw, clusters)
    reference(c(point$alpha, point$beta), gamma = gamma) +
      re_dirichlet(point, fit$xbar, fit$conc) + re_jacobian(point, fit$xbar) -
      sum(point$gamma_raw^2) / 2 - 2.5 * log1p(point$sigmag^2 / 4)
  }, numeric(1))
  unconstrained <- lapply(
    points,
    re_unconstrain,
    stanfit = stanfit,
    rinv = rinv,
    xbar = fit$xbar,
    scale = scale
  )
  density <- vapply(unconstrained, function(u) {
    rstan::log_prob(stanfit, u, adjust_transform = FALSE)
  }, numeric(1))
  expect_equal(diff(density), diff(expected), tolerance = 1e-10)
  # The gradient, including those of gamma_raw and sigmag, against central
  # differences of the log density.
  u <- unconstrained[[1]]
  numeric_gradient <- central_difference_jacobian(
    function(v) rstan::log_prob(stanfit, v), u
  )[1, ]
  expect_equal(
    as.numeric(rstan::grad_log_prob(stanfit, u)),
    unname(numeric_gradient),
    tolerance = 1e-6
  )

  # (2) Against rmsb's own program on the fully observed transitions, with
  # clusters numbered as in the fit, plus the stretches computed in R (the
  # reference above less the observed transitions): the difference, less
  # our Jacobian term sum(log(pi)), is the same at every point. rmsb's
  # predictors are centred at the mean of its own rows, so the Dirichlet
  # terms differ and are removed from both.
  observed <- data[!is.na(data$y) & !is.na(data$yprev), ]
  observed$cl <- factor(observed$id, levels = clusters)
  rmsb_fit <- suppressWarnings(rmsb::blrm(
    ordered(y) ~ yprev + time + tx + cluster(cl),
    data = observed, iter = 10, warmup = 5, chains = 1, refresh = 0,
    x = TRUE, y = TRUE
  ))
  # The last cluster has observed transitions, so rmsb has as many as the fit.
  expect_identical(rmsb_fit$clusterInfo$n, length(clusters))
  rmsb_transform <- rmsb::selectedQr(rmsb_fit$x, center = TRUE)
  previous <- as.integer(as.character(observed$yprev))
  observed_x <- cbind(previous == 2, previous == 3, observed$time, observed$tx)
  difference <- vapply(points, function(point) {
    gamma <- stats::setNames(point$sigmag * point$gamma_raw, clusters)
    eta <- as.numeric(observed_x %*% point$beta) + gamma[as.character(observed$id)]
    cumulative <- cbind(1, stats::plogis(outer(eta, point$alpha, "+")), 0)
    rows <- seq_along(eta)
    observed_loglik <- sum(log(
      cumulative[cbind(rows, observed$y)] - cumulative[cbind(rows, observed$y + 1L)]
    ))
    stretches <- reference(c(point$alpha, point$beta), gamma = gamma) - observed_loglik
    ours <- rstan::log_prob(
      stanfit,
      re_unconstrain(stanfit, point, rinv, fit$xbar, scale = scale),
      adjust_transform = FALSE
    ) - re_dirichlet(point, fit$xbar, fit$conc) - re_jacobian(point, fit$xbar)
    theirs <- rstan::log_prob(
      rmsb_fit$rstan,
      re_unconstrain(rmsb_fit$rstan, point, rmsb_transform$Rinv, rmsb_transform$xbar, tau = TRUE),
      adjust_transform = FALSE
    ) - re_dirichlet(point, rmsb_transform$xbar, rmsb_fit$conc)
    ours - theirs - stretches
  }, numeric(1))
  expect_equal(diff(difference), c(0, 0), tolerance = 1e-8)
})

test_that("without unobserved visits the Stan density is rmsb's", {
  skip_on_cran() # Compiles a Stan program.
  re_skip()
  # rms is attached so that `cluster()` is found; a warning that it was
  # built under a newer R says nothing about the test.
  suppressWarnings(withr::local_package("rms"))
  # Only observed transitions enter: our program's data, built from rmsb's
  # own fit, with no stretch. With the same design rows, centring and
  # clusters, our log density minus rmsb's minus sum(log(pi)), the Jacobian
  # term our program adds because it samples the log-odds of pi, is the same
  # at every point (gamma_raw and sigmag included).
  complete <- ms_fixture("complete")
  rmsb_fit <- suppressWarnings(rmsb::blrm(
    ordered(y) ~ yprev + time + tx + cluster(id),
    data = complete, iter = 10, warmup = 5, chains = 1, refresh = 0,
    x = TRUE, y = TRUE
  ))
  codes <- as.integer(as.factor(rmsb_fit$clusterInfo$cluster))
  plan <- mostr:::markov_empty_path_plan(rmsb_fit$ylevels, 4L)
  plan$stretch_cluster <- integer()
  assembled <- mostr:::markov_blrm_path_stan_data(
    x_setup = rmsb_fit$x,
    keep = rep(TRUE, nrow(rmsb_fit$x)),
    observed_state = as.integer(rmsb_fit$y[, 1L]),
    plan = plan,
    n_states = 4L,
    conc = rmsb_fit$conc,
    clusters = list(n = max(codes), obs = codes)
  )
  expect_identical(assembled$data$Nc, rmsb_fit$clusterInfo$n)
  scale <- mostr:::markov_path_stan_scale(
    assembled$data,
    starts = list(
      state_frequencies = mostr:::markov_path_stan_marginal_start(assembled$data)
    )
  )
  expect_true(scale$maximum_likelihood)
  assembled$data$theta_hat <- scale$theta_hat
  assembled$data$theta_L <- scale$theta_L
  ours <- suppressWarnings(rstan::sampling(
    mostr:::markov_path_stan_model("rstan"),
    data = assembled$data, chains = 1, iter = 2, warmup = 1, refresh = 0
  ))
  rmsb_transform <- rmsb::selectedQr(rmsb_fit$x, center = TRUE)
  difference <- vapply(re_points(max(codes), seed = 12), function(point) {
    rstan::log_prob(
      ours,
      re_unconstrain(ours, point, assembled$rinv, rmsb_transform$xbar, scale = scale),
      adjust_transform = FALSE
    ) - rstan::log_prob(
      rmsb_fit$rstan,
      re_unconstrain(rmsb_fit$rstan, point, rmsb_transform$Rinv, rmsb_transform$xbar, tau = TRUE),
      adjust_transform = FALSE
    ) - re_jacobian(point, rmsb_transform$xbar)
  }, numeric(1))
  expect_equal(diff(difference), c(0, 0), tolerance = 1e-8)
})

test_that("a blrm fit with random intercepts and unobserved states is an rmsb fit", {
  skip_on_cran() # The pass depends on the seed.
  re_skip()
  # rms is attached so that `cluster()` is found; a warning that it was
  # built under a newer R says nothing about the test.
  suppressWarnings(withr::local_package("rms"))
  withr::local_options(mc.cores = 2)
  data <- re_data()
  fit <- suppressMessages(suppressWarnings(ms_fit(
    data,
    formula = ordered(y) ~ yprev + time + tx + cluster(id),
    engine = "blrm", iter = 400, warmup = 200, chains = 2, refresh = 0,
    sampling.args = list(seed = 2026)
  )))
  n_clusters <- fit$clusterInfo$n
  expect_identical(n_clusters, length(attr(fit, "markov_path_fit")$ids))
  # rmsb's fields for random effects.
  expect_identical(dim(fit$omega), c(400L, 1L))
  expect_identical(colnames(fit$omega), "sigmag")
  expect_true(all(fit$omega[, 1] > 0))
  expect_identical(names(fit$gammas), sprintf("gamma[%d]", seq_len(n_clusters)))
  expect_identical(fit$clusterInfo$name, "id")
  expect_length(fit$clusterInfo$cluster, nrow(fit$x))
  expect_identical(
    colnames(fit$draws),
    c("y>=2", "y>=3", "y>=4", "yprev=2", "yprev=3", "time", "tx")
  )
  expect_true("sigmag[1]" %in% rownames(fit$diagnostics))
  gamma <- mostr:::get_blrm_gamma_draws(fit, seq_len(400))
  expect_identical(dim(gamma), c(400L, n_clusters))
  expect_equal(unname(apply(gamma, 2, stats::median)), unname(fit$gammas))

  # Predictions with and without each patient's random effect, including the
  # patient seen only around unobserved visits.
  average <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:3, n_draws = 40, seed = 1)
  totals <- stats::aggregate(estimate ~ tx + time, average, sum)
  expect_equal(totals$estimate, rep(1, nrow(totals)), tolerance = 1e-8)
  with_re <- avg_sops(
    fit, variables = list(tx = c(0, 1)), times = 1:3, n_draws = 40, seed = 1,
    include_re = TRUE
  )
  expect_identical(dim(with_re), dim(average))
  expect_false(isTRUE(all.equal(with_re$estimate, average$estimate)))
  only_gaps <- sops(
    fit,
    newdata = data[data$id == attr(data, "only_gaps"), ][1, ],
    times = 1:2, n_draws = 20, seed = 1, include_re = TRUE
  )
  expect_true(all(is.finite(only_gaps$estimate)))
})

test_that("with random intercepts the posterior is that of rmsb's parametrisation", {
  skip_on_cran() # The pass depends on the seed.
  re_skip()
  # rms is attached so that `cluster()` is found; a warning that it was
  # built under a newer R says nothing about the test.
  suppressWarnings(withr::local_package("rms"))
  withr::local_options(mc.cores = 2)
  # Reference: the program before it sampled on the maximum likelihood scale
  # (rmsb's parameters beta, pi, gamma_raw, sigmag), 4 chains of 5000 warm-up
  # and 5000 sampling iterations on the same data (Nov 2026,
  # benchmarks/local/intermittent-missingness/random-intercepts/ref_small.R):
  # posterior means and SDs. Their Monte Carlo errors are at most 0.02
  # posterior SDs. With 2 x 1000 draws here, the Monte Carlo error of a mean
  # is about 0.05 posterior SD for the coefficients and 0.1 for sigmag
  # (bulk ESS about 100 to 400), so 0.35 SD is more than 3 errors.
  reference_mean <- c(1.9357, -1.0913, -3.7455, 1.1677, 3.1667, -0.071848, -0.63384, 0.62588)
  reference_sd <- c(0.4342, 0.40264, 0.48632, 0.37591, 0.44811, 0.090483, 0.31341, 0.31202)
  fit <- suppressMessages(suppressWarnings(ms_fit(
    re_data(),
    formula = ordered(y) ~ yprev + time + tx + cluster(id),
    engine = "blrm", iter = 2000, warmup = 1000, chains = 2, refresh = 0,
    sampling.args = list(seed = 2026)
  )))
  expect_true(attr(fit, "markov_path_fit")$stan_scale$maximum_likelihood)
  # The package's default target acceptance rate with random intercepts.
  expect_identical(fit$rstan@stan_args[[1]]$control$adapt_delta, 0.9)
  draws <- cbind(fit$draws, sigmag = fit$omega[, 1])
  expect_lt(max(abs(colMeans(draws) - reference_mean) / reference_sd), 0.35)
  expect_equal(unname(apply(draws, 2, stats::sd)), reference_sd, tolerance = 0.2)
})

test_that("without unobserved states a fit with random intercepts is rmsb's own", {
  skip_on_cran() # Fits a blrm model.
  re_skip()
  # rms is attached so that `cluster()` is found; a warning that it was
  # built under a newer R says nothing about the test.
  suppressWarnings(withr::local_package("rms"))
  fit <- suppressWarnings(ms_fit(
    ms_fixture("complete"),
    formula = ordered(y) ~ yprev + time + tx + cluster(id),
    engine = "blrm", iter = 20, warmup = 10, chains = 1, refresh = 0
  ))
  expect_null(attr(fit, "markov_path_fit"))
  expect_false(grepl("path_loglik", fit$stancode))
  expect_identical(length(fit$gammas), fit$clusterInfo$n)
  expect_identical(fit$clusterInfo$n, length(unique(ms_fixture("complete")$id)))
})
