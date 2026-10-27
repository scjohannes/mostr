# Fits with missing and range-coded states (`missing_states`, Ocens()).

ppo_tx <- function() {
  list(family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx))
}

test_that("the likelihood equals enumeration over completions", {
  data <- ml_oracle_data(long_runs = TRUE)
  cases <- list(
    orm = ml_test_problem(data),
    vglm_ppo = ml_test_problem(data, engine = "vglm", engine_args = ppo_tx())
  )
  for (problem in cases) {
    for (seed in 1:2) {
      theta <- ml_random_theta(problem, problem$theta, seed)
      expect_equal(
        markov_ml_loglik(problem, theta),
        ml_oracle_loglik(data, problem, theta),
        tolerance = 1e-10
      )
    }
  }
})

test_that("scores and observed information are exact derivatives", {
  data <- ml_oracle_data(long_runs = TRUE)
  cases <- list(
    orm = ml_test_problem(data),
    vglm_ppo = ml_test_problem(data, engine = "vglm", engine_args = ppo_tx()),
    vglm_constrained = ml_test_problem(
      data,
      engine = "vglm",
      engine_args = list(
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time),
        constraints = list(
          "(Intercept)" = diag(3),
          yprev = matrix(1, 3, 1),
          time = cbind(1, 1:3),
          tx = matrix(1, 3, 1)
        )
      )
    )
  )
  for (problem in cases) {
    theta <- ml_random_theta(problem, problem$theta, 3)
    score <- markov_ml_total_score(problem, theta)
    expect_equal(
      unname(score),
      ml_numeric_gradient(function(x) markov_ml_loglik(problem, x), theta),
      tolerance = 1e-6
    )
    expect_equal(
      unname(colSums(markov_ml_patient_scores(problem, theta))),
      unname(score),
      tolerance = 1e-10
    )
    information <- markov_ml_information(problem, theta)
    hessian <- ml_numeric_jacobian(
      function(x) markov_ml_total_score(problem, x),
      theta
    )
    expect_equal(unname(information), -hessian, tolerance = 1e-6)
  }
})

test_that("the fit maximises the likelihood of the observed data", {
  data <- ml_oracle_data(long_runs = TRUE)
  fit <- mostm(
    Ocens(y_lo, y_hi) ~ yprev + time + tx,
    data = data,
    id_var = "id",
    missing_states = 1:3
  )
  expect_s3_class(fit, c("mostm_ml", "orm"))
  expect_lt(fit$ml$max_abs_score, 1e-6)
  problem <- ml_test_problem(data)
  expect_equal(
    fit$ml$loglik,
    ml_oracle_loglik(data, problem, stats::coef(fit)),
    tolerance = 1e-10
  )
  # The model-based covariance is the inverse observed information.
  expect_equal(
    fit$ml$var_model,
    solve(markov_ml_information(problem, stats::coef(fit))),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_gt(fit$ml$missing_info_fraction, 0)
  expect_lt(fit$ml$missing_info_fraction, 1)
})

test_that("orm and vglm give the same proportional-odds fit", {
  data <- ml_sim_missing()
  formula <- ordered(y) ~ yprev + rms::rcs(time, 3) + tx
  fit_orm <- mostm(formula, data = data, id_var = "id", missing_states = 1:5)
  fit_vglm <- mostm(
    formula,
    data = data,
    id_var = "id",
    missing_states = 1:5,
    engine = "vglm"
  )
  expect_s3_class(fit_vglm, c("mostm_ml", "robcov_vglm"))
  expect_equal(unname(stats::coef(fit_orm)), unname(get_coef(fit_vglm)), tolerance = 1e-7)
  expect_equal(fit_orm$ml$loglik, fit_vglm$ml$loglik, tolerance = 1e-10)
  expect_equal(
    unname(get_vcov_robust(fit_orm)),
    unname(get_vcov_robust(fit_vglm)),
    tolerance = 1e-6
  )
})

test_that("with nothing missing, missing_states gives the engine's usual fit", {
  data <- make_test_data(n_patients = 60, follow_up_time = 10)
  formula <- ordered(y) ~ yprev + time + tx
  for (engine in c("orm", "vglm")) {
    usual <- mostm(formula, data = data, id_var = "id", engine = engine)
    with_states <- mostm(
      formula,
      data = data,
      id_var = "id",
      engine = engine,
      missing_states = 1:5
    )
    expect_false(inherits(with_states, "mostm_ml"))
    expect_identical(get_coef(with_states), get_coef(usual))
    expect_identical(get_vcov_robust(with_states), get_vcov_robust(usual))
  }
})

test_that("without missing_states, rows with missing values are dropped with a warning", {
  data <- ml_sim_missing()
  formula <- ordered(y) ~ yprev + time + tx
  expect_warning(
    fit <- mostm(formula, data = data, id_var = "id"),
    class = "mostr_dropped_rows_warning"
  )
  expect_warning(
    mostm(formula, data = data, id_var = "id"),
    "missing response.*missing previous state.*missing_states"
  )
  complete <- data[stats::complete.cases(data[c("y", "yprev")]), ]
  expect_no_warning(
    reference <- mostm(formula, data = complete, id_var = "id")
  )
  expect_equal(stats::coef(fit), stats::coef(reference))
})

test_that("Ocens() without missing_states drops missing rows and restarts chains", {
  data <- ml_sim_missing()
  data$y_lo <- data$y
  data$y_hi <- data$y
  expect_warning(
    fit <- mostm(
      Ocens(y_lo, y_hi) ~ yprev + time + tx,
      data = data,
      id_var = "id"
    ),
    class = "mostr_dropped_rows_warning"
  )
  expect_s3_class(fit, "mostm_ml")
  expect_equal(fit$ml$n_missing, 0L)
  expect_equal(fit$ml$missing_info_fraction, 0)
  # Every range is a single state, so the fit is the complete-row fit.
  reference <- suppressWarnings(
    mostm(ordered(y) ~ yprev + time + tx, data = data, id_var = "id")
  )
  expect_equal(unname(stats::coef(fit)), unname(stats::coef(reference)), tolerance = 1e-6)
})

test_that("rows with missing covariates are dropped and chains restart", {
  data <- ml_sim_missing()
  data$tx[data$id == data$id[[1L]] & data$time == 3] <- NA
  expect_warning(
    fit <- mostm(
      ordered(y) ~ yprev + time + tx,
      data = data,
      id_var = "id",
      missing_states = 1:5
    ),
    "missing covariates"
  )
  expect_s3_class(fit, "mostm_ml")
  expect_lt(fit$ml$max_abs_score, 1e-6)
})

test_that("arguments of fits with missing states are checked", {
  data <- ml_sim_missing()
  formula <- ordered(y) ~ yprev + time + tx
  expect_error(
    suppressWarnings(mostm(formula, data = data, missing_states = 1:5)),
    "`id_var` is required"
  )
  expect_error(
    mostm(formula, data = data, id_var = "id", missing_states = c(1, NA)),
    "without missing values"
  )
  expect_error(
    mostm(formula, data = data, id_var = "id", missing_states = 1:9),
    "not among the states"
  )
  expect_error(
    mostm(
      formula,
      data = data,
      id_var = "id",
      missing_states = 1:5,
      weights = rep(1, nrow(data))
    ),
    "`weights` cannot be combined"
  )
  expect_error(
    mostm(
      ordered(y) ~ time + tx,
      data = data,
      id_var = "id",
      missing_states = 1:5
    ),
    "must appear on the right-hand side"
  )
  expect_error(
    mostm(
      formula,
      data = data,
      id_var = "id",
      missing_states = 1:5,
      penalty = 1
    ),
    "does not accept `penalty`"
  )
})

test_that("a missing state before later rows cannot be a state never left", {
  data <- ml_oracle_data()
  expect_error(
    mostm(
      Ocens(y_lo, y_hi) ~ yprev + time + tx,
      data = data,
      id_var = "id",
      missing_states = 1:4
    ),
    "never occur as the previous state"
  )
})

test_that("inference works and agrees between engines", {
  data <- ml_sim_missing()
  formula <- ordered(y) ~ yprev + rms::rcs(time, 3) + tx
  fits <- list(
    orm = mostm(formula, data = data, id_var = "id", missing_states = 1:5),
    vglm = mostm(
      formula,
      data = data,
      id_var = "id",
      missing_states = 1:5,
      engine = "vglm"
    )
  )
  results <- lapply(fits, function(fit) {
    sops <- avg_sops(
      fit,
      variables = list(tx = c(0, 1)),
      times = 1:6,
      absorb = "6"
    )
    list(
      delta = inferences(sops, method = "delta"),
      unconditional = inferences(sops, method = "delta", vcov = "unconditional"),
      mvn = inferences(sops, method = "mvn", n_draws = 50, seed = 1),
      score = inferences(sops, method = "score_bootstrap", n_draws = 50, seed = 1)
    )
  })
  for (method in names(results$orm)) {
    expect_equal(
      results$orm[[method]]$conf.low,
      results$vglm[[method]]$conf.low,
      tolerance = 1e-5
    )
  }

  sops <- avg_sops(fits$orm, variables = list(tx = c(0, 1)), times = 1:6, absorb = "6")
  expect_error(
    inferences(sops, method = "fwb", n_draws = 5),
    "fractional weighted bootstrap is not available"
  )
  boot <- inferences(sops, method = "bootstrap", n_draws = 3, seed = 1)
  expect_true(all(is.finite(boot$conf.low)))
  expect_error(
    avg_sops(fits$orm, times = 1:6, absorb = "6", p2_var = "yprev2"),
    "first-order model"
  )
})

test_that("bootstrap refits keep each resampled patient's chain", {
  data <- ml_sim_missing()
  fit <- mostm(
    ordered(y) ~ yprev + time + tx,
    data = data,
    id_var = "id",
    missing_states = 1:5
  )
  boot_data <- markov_model_refit_data(fit)
  boot_data$new_id <- paste0(boot_data$id, "_1")
  boot_data$id <- 1L
  refit <- update_bootstrap_model(fit, boot_data)
  expect_s3_class(refit, "mostm_ml")
  expect_equal(stats::coef(refit), stats::coef(fit), tolerance = 1e-6)
  expect_null(refit$ml$var_robust)
})

test_that("refits of resamples without some state renumber missing_states", {
  data <- ml_sim_missing()
  formula <- ordered(y) ~ yprev + time + tx
  # Patients who are ever in state 3 are left out of the resample.
  in_3 <- unique(data$id[
    (!is.na(data$y) & data$y == 3) | (!is.na(data$yprev) & data$yprev == 3)
  ])
  resample <- data[!data$id %in% in_3, ]
  resample$y <- factor(resample$y, levels = 1:6)
  releveled <- relevel_factors_consecutive(
    resample,
    factor_cols = c("y", "yprev"),
    y_levels = 1:6
  )
  expect_identical(releveled$missing_levels, "3")
  refit_args <- list(formula = formula, missing_states = 1:5)
  expect_equal(
    markov_ml_refit_missing_states(refit_args, releveled$state_mapping),
    1:4
  )
  ocens_args <- list(formula = Ocens(y_lo, y_hi) ~ yprev + time + tx)
  expect_error(
    markov_ml_refit_missing_states(ocens_args, releveled$state_mapping),
    "cannot be renumbered"
  )
})

test_that("logLik(), nobs(), and print() describe the observed data", {
  data <- ml_sim_missing()
  fit <- mostm(
    ordered(y) ~ yprev + time + tx,
    data = data,
    id_var = "id",
    missing_states = 1:5
  )
  expect_equal(as.numeric(stats::logLik(fit)), fit$ml$loglik)
  expect_equal(attr(stats::logLik(fit), "df"), length(stats::coef(fit)))
  expect_equal(stats::nobs(fit), fit$ml$n_informative_rows)
  expect_output(print(fit), "Robust S.E.")
  expect_output(print(fit), "A missing state may be: 1, 2, 3, 4, 5")
})
