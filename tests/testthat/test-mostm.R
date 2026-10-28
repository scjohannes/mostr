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

test_that("orm and vglm fits do not load the blrm engine's packages", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("rmsb")
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  # rmsb and rstan take seconds to load, so an orm or vglm fit, and every
  # bootstrap worker that refits one, must not load them. A fresh R process
  # shows which namespaces a fit loads.
  dev_path <- if (
    requireNamespace("pkgload", quietly = TRUE) &&
      pkgload::is_dev_package("mostr")
  ) {
    pkgload::pkg_path(testthat::test_path())
  }
  loaded <- callr::r(
    function(dev_path) {
      if (is.null(dev_path)) {
        library(mostr)
      } else {
        pkgload::load_all(dev_path, quiet = TRUE)
      }
      trajectories <- sim_markov(
        "actt2",
        n_patients = 40,
        follow_up_time = 6,
        seed = 1
      )
      data <- prepare_markov_data(trajectories, absorbing_state = 8)
      mostm(ordered(y) ~ yprev + time + tx, data = data, id_var = "id")
      mostm(
        ordered(y) ~ yprev + time + tx,
        data = data,
        id_var = "id",
        engine = "vglm"
      )
      c("rmsb", "rstan") %in% loadedNamespaces()
    },
    args = list(dev_path = dev_path)
  )

  expect_identical(loaded, c(FALSE, FALSE))
})

test_that("mostm() reports an orm fit that rms could not fit", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6201)
  # `maxit = 0` stops rms::orm() before it converges.
  for (id_var in list("id", NULL)) {
    expect_error(
      suppressMessages(suppressWarnings(capture.output(
        mostm(
          ordered(y) ~ time + tx + yprev,
          data = data,
          id_var = id_var,
          maxit = 0
        )
      ))),
      "`rms::orm()` could not fit the model",
      fixed = TRUE
    )
  }
})

test_that("mostm() accepts the engine arguments used in the articles", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  expect_contains(
    mostr:::mostm_engine_args("orm"),
    c(
      "opt_method",
      "scale",
      "maxit",
      "penalty",
      "var.penalty",
      "na.action",
      "family"
    )
  )
  expect_contains(
    mostr:::mostm_engine_args("vglm"),
    c(
      "family",
      "constraints",
      "etastart",
      "mustart",
      "coefstart",
      "control",
      "maxit",
      "trace"
    )
  )
  skip_if_not_installed("rmsb")
  expect_contains(
    mostr:::mostm_engine_args("blrm"),
    c("ppo", "cppo", "iter", "warmup", "chains", "refresh", "method", "backend")
  )
})

test_that("every engine warns without id_var, even without an id column", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 6112)
  data$id <- NULL
  expect_warning(
    mostm(ordered(y) ~ time + tx + yprev, data = data),
    "`id_var` was not supplied to mostm()",
    fixed = TRUE
  )
  expect_warning(
    mostm(ordered(y) ~ time + tx + yprev, data = data, engine = "vglm"),
    "`id_var` was not supplied to mostm()",
    fixed = TRUE
  )
  skip_if_not_installed("rmsb")
  expect_warning(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      engine = "blrm",
      standata = TRUE
    ),
    "`id_var` was not supplied to mostm()",
    fixed = TRUE
  )
})

test_that("orm refits keep the fit's na.action and penalty", {
  skip_if_not_installed("rms")
  skip_if_not_installed("Hmisc")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6111)
  # Missing treatment on later visits only, so every patient keeps a complete
  # starting row.
  later <- which(data$time > 2)
  data$tx[later[c(3, 17, 40)]] <- NA

  # rms::orm() accepts Hmisc::na.delete() but fails with stats::na.omit()
  # when rows have missing values.
  fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    na.action = Hmisc::na.delete,
    penalty = 2,
    var.penalty = "sandwich"
  )
  unpenalized <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    na.action = Hmisc::na.delete
  )
  expect_gt(
    max(abs(stats::coef(fit) - stats::coef(unpenalized))),
    1e-4
  )

  # Refitting the stored data reproduces the penalized fit.
  boot <- mostr:::bootstrap_analysis_wrapper(
    boot_data = mostr:::markov_model_refit_data(fit),
    model = fit,
    factor_cols = character(0),
    original_data = data,
    update_datadist = FALSE
  )
  expect_identical(mostr:::bootstrap_model_call(boot$model)$penalty, 2)
  expect_equal(stats::coef(boot$model), stats::coef(fit), tolerance = 1e-6)

  out <- avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    y_levels = 1:6,
    absorb = 6
  )
  withr::local_seed(6113)
  inferred <- inferences(out, method = "bootstrap", n_draws = 2)
  expect_equal(attr(inferred, "n_successful"), 2L)
})

test_that("orm fits reject factor() written inside the formula", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6301)

  # rms::orm() cannot build the design matrix for such terms and stops with
  # "subscript out of bounds", so mostm() stops before fitting.
  err <- expect_error(
    mostm(ordered(y) ~ factor(time) + tx + yprev, data = data, id_var = "id"),
    "`factor(time)`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "`rms::orm()`", fixed = TRUE)
  expect_match(conditionMessage(err), "column of `data`", fixed = TRUE)
  expect_match(conditionMessage(err), "first_followup_time", fixed = TRUE)
  expect_match(
    conditionMessage(err),
    "vignette(\"factor-time-orm-sops\")",
    fixed = TRUE
  )

  # The time advice follows `time_var`.
  renamed <- data
  names(renamed)[names(renamed) == "time"] <- "day"
  err <- expect_error(
    mostm(
      ordered(y) ~ factor(day) + tx + yprev,
      data = renamed,
      id_var = "id",
      time_var = "day"
    ),
    "`factor(day)`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "first_followup_time", fixed = TRUE)

  formulas <- list(
    ordered(y) ~ time + as.factor(tx) + yprev,
    ordered(y) ~ time + base::factor(tx) + yprev,
    ordered(y) ~ time * ordered(tx) + yprev,
    ordered(y) ~ time + yprev + as.ordered(tx)
  )
  for (formula in formulas) {
    err <- expect_error(
      mostm(formula, data = data, id_var = "id"),
      "column of `data`",
      fixed = TRUE
    )
    expect_no_match(conditionMessage(err), "first_followup_time", fixed = TRUE)
  }

  # A factor column in `data` fits, and `ordered()` on the response is allowed.
  data$tx_factor <- factor(data$tx)
  fit <- mostm(ordered(y) ~ time + tx_factor + yprev, data = data, id_var = "id")
  expect_s3_class(fit, "orm")
})

test_that("blrm fits reject factor() written inside the formula", {
  skip_if_not_installed("rms")
  skip_if_not_installed("rmsb")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6302)
  err <- expect_error(
    mostm(
      ordered(y) ~ factor(time) + tx + yprev,
      data = data,
      id_var = "id",
      engine = "blrm"
    ),
    "`factor(time)`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "`rmsb::blrm()`", fixed = TRUE)
  expect_match(conditionMessage(err), "first_followup_time", fixed = TRUE)
})

test_that("vglm fits accept factor() written inside the formula", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6303)
  fit <- mostm(
    ordered(y) ~ factor(time) + tx + yprev,
    data = data,
    id_var = "id",
    engine = "vglm"
  )
  expect_s3_class(fit, "robcov_vglm")
})
