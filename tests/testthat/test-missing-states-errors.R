# The data contract for unobserved states and the arguments that can't be
# combined with them.

test_that("`yprev` must be missing exactly after an unobserved visit", {
  data <- ms_fixture("gaps")
  hidden <- attr(data, "hidden")
  # Carry the last observed state forward into `yprev`.
  filled <- data
  row <- which(filled$id == hidden$id[1] & filled$time == hidden$time[1] + 1)
  filled$yprev[row] <- "2"
  expect_error(
    ms_fit(filled),
    "`yprev` is recorded on 1 row whose previous visit has no exactly observed state",
    class = "mostr_error",
    fixed = TRUE
  )
  # A missing `yprev` after an observed visit.
  dangling <- data
  # A patient followed for six visits without unobserved ones.
  complete_id <- attr(data, "special")$alive_long[10]
  dangling$yprev[dangling$id == complete_id & dangling$time == 3] <- NA
  expect_error(
    ms_fit(dangling),
    paste0(
      "`yprev` is missing on 1 row although the patient's previous row has ",
      "an observed state (first: `id` = ",
      complete_id,
      ")"
    ),
    class = "mostr_error",
    fixed = TRUE
  )
  # Both are fine with `missing_states = "drop"`, which uses neither row.
  expect_no_error(ms_fit(filled, missing_states = "drop"))
  expect_no_error(ms_fit(dangling, missing_states = "drop"))
})

test_that("a stretch needs every predictor", {
  data <- ms_fixture("gaps")
  hidden <- attr(data, "hidden")
  row <- which(data$id == hidden$id[1] & data$time == hidden$time[1])
  data$tx <- as.numeric(data$tx)
  data$tx[row] <- NA
  expect_error(
    ms_fit(data),
    "A visit with an unobserved state needs every predictor",
    class = "mostr_error"
  )
})

test_that("row weights must be constant within a stretch", {
  data <- ms_fixture("gaps")
  per_row <- seq(0.5, 1.5, length.out = nrow(data))
  expect_error(
    ms_fit(data, weights = per_row),
    "Row weights must be the same on all rows of a run of visits",
    class = "mostr_error"
  )
})

test_that("weights constant per patient enter the likelihood per patient", {
  skip_on_cran() # Slow: two weighted fits.
  data <- ms_fixture("gaps")
  patient_weight <- 1 + (data$id %% 3) / 2
  fit <- ms_fit(data, weights = patient_weight)
  reference <- ms_reference(data, weights = patient_weight)
  expect_equal(
    attr(fit, "markov_path_fit")$loglik,
    reference(coef(fit)),
    tolerance = 1e-8
  )
  ms_expect_maximum(fit, reference)

  # Doubling every weight doubles the log-likelihood and the information:
  # the estimates stay, the model-based covariance halves, and the
  # patient-clustered covariance stays.
  doubled <- 2 * patient_weight
  fit_doubled <- ms_fit(data, weights = doubled)
  expect_equal(coef(fit_doubled), coef(fit), tolerance = 1e-8)
  expect_equal(fit_doubled$orig.var, fit$orig.var / 2, tolerance = 1e-8)
  expect_equal(fit_doubled$var, fit$var, tolerance = 1e-8)
})

test_that("unobserved states refuse a second-order model", {
  data <- ms_fixture("gaps")
  data$yprev2 <- data$yprev
  expect_error(
    mostm(
      ordered(y) ~ yprev + yprev2 + time + tx,
      data = data,
      id_var = "id",
      absorb = 4,
      p2_var = "yprev2"
    ),
    "A second-order model (`p2_var`) can't be combined with unobserved states",
    class = "mostr_error",
    fixed = TRUE
  )
})

test_that("orm with unobserved states refuses a penalty and other links", {
  data <- ms_fixture("gaps")
  expect_error(
    ms_fit(data, penalty = 1),
    "`penalty` can't be combined with unobserved states",
    class = "mostr_error",
    fixed = TRUE
  )
  expect_error(
    ms_fit(data, family = "probit"),
    "A link other than the logit can't be combined with unobserved states",
    class = "mostr_error",
    fixed = TRUE
  )
})

test_that("vglm needs `missing_states = \"drop\"` for unobserved states", {
  data <- ms_fixture("gaps")
  expect_error(
    ms_fit(data, engine = "vglm"),
    "`engine = \"vglm\"` can only leave these rows out",
    class = "mostr_error",
    fixed = TRUE
  )
  expect_error(
    ms_fit(data, engine = "vglm", missing_states = "average"),
    "`missing_states = \"average\"` is available only for",
    class = "mostr_error",
    fixed = TRUE
  )
  skip_on_cran() # Slow: a vglm fit and an orm fit.
  fit <- ms_fit(data, engine = "vglm", missing_states = "drop")
  expect_s3_class(fit, "robcov_vglm")
  # vglm and orm reach the same maximum of the likelihood on the same rows.
  fit_orm_drop <- ms_fit(data, missing_states = "drop")
  expect_equal(
    as.numeric(stats::logLik(fit$vglm_fit)),
    ms_loglik(fit_orm_drop),
    tolerance = 1e-6
  )
})

test_that("`missing_states` takes \"average\" or \"drop\"", {
  data <- ms_fixture("complete")
  for (value in list("avg", c("average", "drop"), NA, 1, TRUE)) {
    expect_error(
      ms_fit(data, missing_states = value),
      "`missing_states` must be \"average\" or \"drop\".",
      class = "mostr_error",
      fixed = TRUE
    )
  }
})

test_that("blrm refuses arguments without a path likelihood before fitting", {
  skip_on_cran() # Slow: loads rmsb.
  skip_if_not_installed("rmsb")
  data <- ms_fixture("gaps")
  expect_error(
    ms_fit(data, engine = "blrm", ppo = ~time),
    "A partial proportional odds model (`ppo`) can't be combined with unobserved states",
    class = "mostr_error",
    fixed = TRUE
  )
  expect_error(
    ms_fit(data, engine = "blrm", formula = ordered(y) ~ yprev + time + tx + cluster(id)),
    "A random intercept (`cluster()`) can't be combined with unobserved states",
    class = "mostr_error",
    fixed = TRUE
  )
  expect_error(
    ms_fit(data, engine = "blrm", iprior = 1),
    "`iprior` can't be combined with unobserved states",
    class = "mostr_error",
    fixed = TRUE
  )
})
