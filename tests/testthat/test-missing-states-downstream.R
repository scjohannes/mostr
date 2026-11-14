# The functions that read a fit work on a fit with unobserved states.

test_that("predictions, state times and comparisons work on a path fit", {
  skip_on_cran() # Slow: predictions; the variance tests run delta SEs off CRAN too.
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  expect_false(is.null(attr(fit, "markov_path_fit")))

  individual <- sops(fit, times = 1:3)
  expect_true(all(is.finite(individual$estimate)))
  by_row <- stats::aggregate(estimate ~ rowid + time, individual, sum)
  expect_equal(by_row$estimate, rep(1, nrow(by_row)), tolerance = 1e-12)

  average <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:6)
  totals <- stats::aggregate(estimate ~ tx + time, average, sum)
  expect_equal(totals$estimate, rep(1, nrow(totals)), tolerance = 1e-12)

  times_in_state <- avg_state_time(fit, variables = list(tx = c(0, 1)), times = 1:6)
  totals <- stats::aggregate(estimate ~ tx, times_in_state, sum)
  expect_equal(totals$estimate, c(6, 6), tolerance = 1e-10)
  individual_time <- state_time(fit, times = 1:6)
  expect_true(all(is.finite(individual_time$estimate)))

  contrast <- avg_comparisons(fit, variables = list(tx = c(0, 1)), times = 1:6)
  expect_true(all(is.finite(contrast$estimate)))
  expect_equal(sum(contrast$estimate[contrast$time == 6]), 0, tolerance = 1e-12)

  printed <- utils::capture.output(print(fit))
  expect_true(any(grepl("Logistic (Proportional Odds)", printed, fixed = TRUE)))
})

test_that("MVN draws use the path fit's covariance", {
  skip_on_cran() # The pass depends on the seed.
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  contrast <- avg_comparisons(fit, variables = list(tx = c(0, 1)), times = 1:3)
  mvn <- inferences(contrast, method = "mvn", n_draws = 200, seed = 1)
  delta <- inferences(contrast, method = "delta", vcov = "conditional")
  expect_true(all(is.finite(mvn$std.error)))
  # MVN draws from the same covariance as the conditional delta method.
  expect_equal(mvn$std.error, delta$std.error, tolerance = 0.25)
})

test_that("bootstrap and FWB refits succeed on a path fit", {
  skip_on_cran() # Slow: refits.
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  out <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:3)
  for (method in c("bootstrap", "fwb")) {
    warnings <- character()
    result <- withCallingHandlers(
      inferences(out, method = method, n_draws = 3, seed = 2),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl("fitting failed", warnings)), info = method)
    expect_identical(attr(result, "n_successful"), 3L, info = method)
    expect_true(all(is.finite(result$std.error)), info = method)
  }
})

test_that("bootstrap refits maximise the likelihood with unobserved states", {
  skip_on_cran() # Slow: four fits.
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  refit_data <- attr(fit, "markov_refit_data")
  # The stored refit data keep the unobserved visits.
  expect_identical(sum(is.na(refit_data$y)), sum(is.na(data$y)))

  # Refitting the original sample reproduces the path fit, not the fit that
  # drops the unobserved visits.
  refit <- mostr:::bootstrap_analysis_wrapper(
    boot_data = refit_data,
    model = fit,
    factor_cols = c("y", "yprev"),
    original_data = refit_data
  )$model
  expect_equal(coef(refit), coef(fit), tolerance = 1e-10)
  expect_false(is.null(attr(refit, "markov_path_fit")))
  fit_drop <- ms_fit(data, missing_states = "drop")
  expect_gt(max(abs(coef(refit) - coef(fit_drop))), 1e-3)

  # A resample in which one patient appears twice: the refit equals a fit to
  # data in which that patient's rows are copied under a new ID.
  twice <- attr(data, "special")$alive_long[3]
  copy <- refit_data[refit_data$id == twice, ]
  copy$id <- 1000L
  doubled <- rbind(refit_data, copy)
  refit_doubled <- mostr:::bootstrap_analysis_wrapper(
    boot_data = doubled,
    model = fit,
    factor_cols = c("y", "yprev"),
    original_data = refit_data
  )$model
  expect_equal(
    coef(refit_doubled),
    coef(ms_fit(doubled)),
    tolerance = 1e-10
  )
  expect_equal(
    attr(refit_doubled, "markov_path_fit")$loglik,
    ms_reference(doubled)(coef(refit_doubled)),
    tolerance = 1e-8
  )
})

test_that("a patient with only unobserved transitions stays in the cohort", {
  skip_on_cran() # Slow: a fit and predictions.
  data <- ms_fixture("downstream")
  # Patient 1000: state 2 before day 1, unobserved on day 1, state 3 on day 2.
  # Its only likelihood term is a stretch.
  extra <- data.frame(
    id = 1000L,
    time = 1:2,
    tx = 1,
    yprev = factor(c(2, NA), levels = 1:3),
    y = c(NA, 3L)
  )
  with_extra <- rbind(data, extra)
  fit <- expect_no_warning(ms_fit(with_extra))
  info <- attr(fit, "markov_path_fit")
  expect_true("1000" %in% info$ids)
  expect_true("1000" %in% as.character(attr(fit, "markov_starting_profile_data")$id))
  individual <- sops(fit, times = 1:2)
  expect_true(1000 %in% individual$id)
  reference <- ms_reference(with_extra)
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  # Its score is the gradient of its stretch alone.
  expect_equal(
    unname(info$scores["1000", ]),
    unname(central_difference_jacobian(
      function(theta) reference(theta, per_patient = TRUE)[["1000"]],
      coef(fit)
    )[1, ]),
    tolerance = 1e-6
  )
})

test_that("a patient without a first previous state has no starting row", {
  skip_on_cran() # Slow: two fits.
  # The patient's later transitions are fitted, but predictions start from
  # each patient's first row, whose previous state is unknown here. The
  # same holds with `missing_states = "drop"`.
  data <- ms_fixture("gaps")
  first <- attr(data, "special")$unknown_start
  for (missing_states in c("average", "drop")) {
    fit <- ms_fit(data, missing_states = missing_states)
    expect_error(
      avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:2),
      paste0("incomplete starting-profile predictors: ", first),
      class = "mostr_error"
    )
  }
})
