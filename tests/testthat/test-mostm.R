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
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = as.list(data),
      id_var = "id"
    ),
    "`data` must be supplied as a data frame.",
    fixed = TRUE
  )
  expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "rms"
    ),
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
      original_data = data
    )
    refit_call <- mostr:::bootstrap_model_call(boot$model)

    expect_identical(refit_call[[1L]], quote(mostm))
    expect_identical(refit_call$engine, engine)
    expect_identical(refit_call$id_var, "id")
    # Refits need only the coefficients, so they skip the patient-clustered
    # covariance and the starting rows.
    expect_null(attr(boot$model, "markov_robust_covariance"))
    expect_null(attr(boot$model, "markov_starting_profile_data"))
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
      original_data = data
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
  expect_error(
    suppressMessages(suppressWarnings(capture.output(
      mostm(
        ordered(y) ~ time + tx + yprev,
        data = data,
        id_var = "id",
        maxit = 0
      )
    ))),
    "`rms::orm()` could not fit the model",
    fixed = TRUE
  )
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

test_that("every engine refuses a fit without id_var", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  # A column called `id` must not stand in for the missing argument.
  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 6112)
  engines <- list(
    orm = list(),
    vglm = list(
      engine = "vglm",
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE)
    )
  )
  if (requireNamespace("rmsb", quietly = TRUE)) {
    engines$blrm <- list(engine = "blrm", standata = TRUE)
  }
  for (engine in names(engines)) {
    base_args <- c(
      list(ordered(y) ~ time + tx + yprev, data = data),
      engines[[engine]]
    )
    for (id_args in list(list(), list(id_var = NULL))) {
      expect_error(
        do.call(mostm, c(base_args, id_args)),
        "`mostm()` needs `id_var`",
        fixed = TRUE,
        info = paste(engine, length(id_args))
      )
    }
  }
})

test_that("a patient column with another name gives the same results as `id`", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6113)
  renamed <- data
  names(renamed)[names(renamed) == "id"] <- "pid"
  fit_model <- function(data, id_var) {
    mostm(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = id_var,
      engine = "vglm"
    )
  }
  fit_id <- fit_model(data, "id")
  fit_pid <- fit_model(renamed, "pid")
  expect_equal(stats::vcov(fit_pid), stats::vcov(fit_id))

  args <- list(times = 1:3, absorb = 6, variables = list(tx = c(0, 1)))
  profiles_id <- data[data$time == 1, , drop = FALSE]
  profiles_pid <- renamed[renamed$time == 1, , drop = FALSE]
  for (fun in c("avg_sops", "avg_comparisons")) {
    expected <- do.call(fun, c(list(fit_id, newdata = profiles_id), args))
    out <- do.call(fun, c(list(fit_pid, newdata = profiles_pid), args))
    expect_equal(out$estimate, expected$estimate, info = fun)

    # Conditional variance and the bootstrap both use patient clusters.
    for (method in c("delta", "bootstrap")) {
      extra <- if (method == "bootstrap") {
        list(n_draws = 5, seed = 7)
      } else {
        list(vcov = "conditional")
      }
      expected_ci <- do.call(inferences, c(list(expected, method = method), extra))
      out_ci <- do.call(inferences, c(list(out, method = method), extra))
      expect_equal(out_ci$std.error, expected_ci$std.error, info = paste(fun, method))
    }
  }
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
    original_data = data
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

test_that("orm fits reject random-effect terms before fitting", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6304)
  data$mre <- as.numeric(data$time > 1)

  # Stopping before the fit keeps this the same for every rms version: rms 8.2
  # and later fit these terms, and earlier versions do not have the functions.
  err <- expect_error(
    mostm(ordered(y) ~ time + tx + yprev + cluster(id), data = data, id_var = "id"),
    "`cluster(id)`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "does not support random effects", fixed = TRUE)
  expect_match(conditionMessage(err), "random effect is 0", fixed = TRUE)
  expect_match(conditionMessage(err), "`id_var`", fixed = TRUE)
  expect_match(conditionMessage(err), "`engine = \"blrm\"`", fixed = TRUE)

  # Found with a package prefix, inside other calls, and next to each other.
  err <- expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev + cluster(id) + mix_re(mre),
      data = data,
      id_var = "id"
    ),
    "`cluster(id)` and `mix_re(mre)`",
    fixed = TRUE
  )
  labels <- list(
    "`rms::cluster(id)`" = ordered(y) ~ time + tx + yprev + rms::cluster(id),
    "`rms::mix_re(mre)`" = ordered(y) ~ time + tx + yprev + rms::mix_re(mre),
    "`cluster(id)`" = ordered(y) ~ time * cluster(id) + yprev
  )
  for (label in names(labels)) {
    expect_error(
      mostm(labels[[label]], data = data, id_var = "id"),
      label,
      fixed = TRUE
    )
  }

  # Only calls count: a column named `cluster` is an ordinary predictor.
  data$cluster <- data$tx
  fit <- mostm(ordered(y) ~ time + cluster + yprev, data = data, id_var = "id")
  expect_s3_class(fit, "orm")
})

test_that("vglm fits reject random-effect terms before fitting", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6305)
  # VGAM would otherwise fit the patient ID number as a straight-line predictor.
  cluster <- function(x) x
  err <- expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev + cluster(id),
      data = data,
      id_var = "id",
      engine = "vglm"
    ),
    "`cluster(id)`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "does not support random effects", fixed = TRUE)
  expect_match(conditionMessage(err), "ordinary predictor", fixed = TRUE)
  expect_match(conditionMessage(err), "`engine = \"blrm\"`", fixed = TRUE)
})

test_that("blrm fits accept only an unprefixed cluster() random intercept", {
  skip_if_not_installed("rms")
  skip_if_not_installed("rmsb")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6306)
  data$mre <- as.numeric(data$time > 1)

  # rmsb::blrm() recognizes only `cluster(id)` as a random intercept. It would
  # fit the other forms as ordinary predictors, such as the patient ID number.
  err <- expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev + rms::cluster(id),
      data = data,
      id_var = "id",
      engine = "blrm"
    ),
    "`rms::cluster(id)`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "ordinary predictor", fixed = TRUE)
  expect_match(conditionMessage(err), "without a package prefix", fixed = TRUE)
  err <- expect_error(
    mostm(
      ordered(y) ~ time + tx + yprev + cluster(id) + mix_re(mre),
      data = data,
      id_var = "id",
      engine = "blrm"
    ),
    "`mix_re(mre)`",
    fixed = TRUE
  )
  expect_no_match(conditionMessage(err), "`cluster(id)` and", fixed = TRUE)

  # `standata = TRUE` prepares the data for Stan without sampling. `cluster()`
  # comes from rms from version 8.2 on and from rmsb before that, and rms 8.1
  # does not look for it in the formula's environment, so it is global here.
  assign("cluster", function(x) x, envir = globalenv())
  withr::defer(rm("cluster", envir = globalenv()))
  fit <- mostm(
    ordered(y) ~ time + tx + yprev + cluster(id),
    data = data,
    id_var = "id",
    engine = "blrm",
    standata = TRUE
  )
  expect_identical(mostr:::markov_model_engine(fit), "blrm")
  expect_equal(fit$Nc, length(unique(data$id)))
})

test_that("model-based functions reject orm fits that hold random effects", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6307)
  fit <- mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id")
  expect_null(mostr:::validate_markov_model(fit))

  # A fit saved before mostm() rejected these terms can still hold them.
  for (field in c("sigma", "sigma1", "sigma2", "gamma")) {
    with_random_effects <- fit
    with_random_effects[[field]] <- 1
    expect_error(
      mostr:::validate_markov_model(with_random_effects),
      "random effects",
      fixed = TRUE
    )
  }
  with_random_effects <- fit
  with_random_effects$sigma <- 1
  expect_error(
    sops(with_random_effects, times = 1:2, y_levels = 1:6, absorb = 6),
    "random effects",
    fixed = TRUE
  )
})
