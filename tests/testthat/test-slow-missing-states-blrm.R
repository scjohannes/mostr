# engine = "blrm" with unobserved states samples from mostr's own Stan
# program. Compiling it takes minutes the first time (it is then cached in
# tools::R_user_dir("mostr", "cache")), so this runs only with
# MOSTR_RUN_SLOW=true. The test uses that default cache on purpose: a
# temporary directory would force the compile on every run. It writes there
# what any first blrm fit with unobserved states writes, and the cache code
# removes the files of earlier versions of the program.

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
