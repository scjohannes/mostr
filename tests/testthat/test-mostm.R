test_that("mostm() records its call and engine for each frequentist engine", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6101)

  fit_orm <- mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id")
  expect_s3_class(fit_orm, "orm")
  expect_identical(mostr:::markov_model_engine(fit_orm), "orm")
  expect_identical(fit_orm$call[[1L]], quote(mostm))
  expect_identical(fit_orm$call$engine, "orm")
  expect_false(is.null(attr(fit_orm, "markov_robust_covariance")))

  fit_vglm <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    engine = "vglm"
  )
  expect_s3_class(fit_vglm, "robcov_vglm")
  expect_identical(mostr:::markov_model_engine(fit_vglm), "vglm")
  expect_identical(
    mostr:::markov_model_engine(fit_vglm$vglm_fit),
    "vglm"
  )
  expect_identical(fit_vglm$call[[1L]], quote(mostm))
  expect_identical(
    methods::slot(fit_vglm$vglm_fit, "call")$engine,
    "vglm"
  )
})

test_that("the vglm engine fits the orm proportional odds model by default", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 60, follow_up_time = 8, seed = 6102)
  fit_orm <- mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id")
  fit_default <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    engine = "vglm"
  )
  fit_explicit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    engine = "vglm",
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE)
  )

  expect_equal(stats::coef(fit_default), stats::coef(fit_explicit))
  expect_equal(
    as.numeric(stats::logLik(fit_default$vglm_fit)),
    as.numeric(stats::logLik(fit_orm)),
    tolerance = 1e-6
  )
  expect_equal(
    unname(stats::coef(fit_default)[c("tx", "time")]),
    unname(stats::coef(fit_orm)[c("tx", "time")]),
    tolerance = 1e-4
  )
})

test_that("mostm() names the engine an argument belongs to", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rmsb")
  local_reproducible_output(width = 80)

  data <- make_test_data(n_patients = 20, follow_up_time = 6, seed = 6103)
  expect_snapshot(
    mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", ppo = ~time),
    error = TRUE
  )
  expect_snapshot(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "vglm",
      iter = 100,
      maxitt = 10
    ),
    error = TRUE
  )
  expect_snapshot(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "blrm",
      constraints = list()
    ),
    error = TRUE
  )
  expect_snapshot(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      backend = "cmdstan"
    ),
    error = TRUE
  )
  expect_snapshot(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      family = VGAM::cumulative(reverse = TRUE)
    ),
    error = TRUE
  )
})

test_that("mostm() accepts the arguments of the chosen engine", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rmsb")

  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 6104)
  expect_no_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      maxit = 30,
      family = "logistic"
    )
  )
  expect_no_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "vglm",
      maxit = 30,
      trace = FALSE
    )
  )
  # rmsb::blrm() has its own `backend`, which mostm() passes on.
  expect_no_error(
    mostr:::mostm_check_engine_args(
      "blrm",
      list(backend = "cmdstan", iter = 100, ppo = quote(~time)),
      list()
    )
  )
})

test_that("mostm() rejects arguments that every fit must keep", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 20, follow_up_time = 6, seed = 6105)
  expect_error(
    mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", x = FALSE),
    "always keeps the design matrix and response",
    fixed = TRUE
  )
  expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "vglm",
      x.arg = FALSE,
      y.arg = FALSE
    ),
    "remove `x.arg`, `y.arg`",
    fixed = TRUE
  )
  expect_error(
    mostm(ordered(y) ~ time + tx + yprev, data, "vglm"),
    "must be named",
    fixed = TRUE
  )
  expect_error(
    mostm(ordered(y) ~ time + tx + yprev, data = as.list(data)),
    "`data` must be supplied as a data frame.",
    fixed = TRUE
  )
  expect_error(
    mostm(ordered(y) ~ time + tx + yprev, data = data, engine = "rms"),
    "should be one of"
  )
})

test_that("the blrm engine rejects robust-covariance settings and weights", {
  skip_if_not_installed("rms")
  skip_if_not_installed("rmsb")

  data <- make_test_data(n_patients = 20, follow_up_time = 6, seed = 6106)
  data$w <- 1
  expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "blrm",
      type = "HC1"
    ),
    "posterior draws",
    fixed = TRUE
  )
  expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "blrm",
      cadjust = FALSE
    ),
    "posterior draws",
    fixed = TRUE
  )
  expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "blrm",
      weights = w
    ),
    "does not accept row weights",
    fixed = TRUE
  )

  fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    engine = "blrm",
    standata = TRUE
  )
  expect_identical(mostr:::markov_model_engine(fit), "blrm")
})

test_that("bootstrap refits call mostm() with the original engine", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6107)
  fits <- list(
    orm = mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id"),
    vglm = mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "vglm"
    )
  )

  for (engine in names(fits)) {
    fit <- fits[[engine]]
    boot <- mostr:::bootstrap_analysis_wrapper(
      boot_data = mostr:::markov_model_refit_data(fit),
      model = fit,
      factor_cols = character(0),
      original_data = data,
      update_datadist = FALSE
    )
    refit_call <- mostr:::bootstrap_model_call(boot$model)

    expect_identical(refit_call[[1L]], quote(mostm))
    expect_identical(refit_call$engine, engine)
    # Refits need only the coefficients, so they skip the patient-clustered
    # covariance.
    expect_null(refit_call$id_var)
    expect_identical(mostr:::markov_model_engine(boot$model), engine)
    expect_equal(
      stats::coef(boot$model),
      stats::coef(fit),
      tolerance = 1e-6
    )
  }
  expect_s4_class(
    mostr:::bootstrap_analysis_wrapper(
      boot_data = mostr:::markov_model_refit_data(fits$vglm),
      model = fits$vglm,
      factor_cols = character(0),
      original_data = data,
      update_datadist = FALSE
    )$model,
    "vglm"
  )
})

test_that("fits without mostm() provenance must be refitted", {
  skip_if_not_installed("rms")
  local_reproducible_output(width = 80)

  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 6108)
  fit <- mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id")
  old <- fit
  attr(old, "markov_engine") <- NULL
  attr(old, "markov_fit_wrapper") <- "orm_markov"

  expect_snapshot(
    sops(old, times = 1:2, y_levels = 1:6, absorb = 6),
    error = TRUE
  )
})
