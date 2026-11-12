test_that("mostm() records its call and engine, and bootstrap refits call it with that engine", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6101)
  fits <- list(
    orm = mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id"),
    vglm = mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "vglm"
    )
  )

  # What each engine records on the fit itself.
  expect_s3_class(fits$orm, "orm")
  expect_identical(fits$orm$call$engine, "orm")
  expect_false(is.null(attr(fits$orm, "markov_robust_covariance")))
  expect_s3_class(fits$vglm, "robcov_vglm")
  expect_identical(mostr:::markov_model_engine(fits$vglm$vglm_fit), "vglm")
  expect_identical(
    methods::slot(fits$vglm$vglm_fit, "call")$engine,
    "vglm"
  )

  for (engine in names(fits)) {
    fit <- fits[[engine]]
    expect_identical(mostr:::markov_model_engine(fit), engine)
    expect_identical(fit$call[[1L]], quote(mostm), info = engine)

    boot <- mostr:::bootstrap_analysis_wrapper(
      boot_data = mostr:::markov_model_refit_data(fit),
      model = fit,
      factor_cols = character(0),
      original_data = data
    )
    refit_call <- mostr:::bootstrap_model_call(boot$model)

    expect_identical(refit_call[[1L]], quote(mostm), info = engine)
    expect_identical(refit_call$engine, engine, info = engine)
    expect_identical(refit_call$id_var, "id", info = engine)
    # Refits need only the coefficients, so they skip the patient-clustered
    # covariance and the starting rows.
    expect_null(attr(boot$model, "markov_robust_covariance"))
    expect_null(attr(boot$model, "markov_starting_profile_data"))
    expect_identical(mostr:::markov_model_engine(boot$model), engine)
    expect_equal(
      stats::coef(boot$model),
      stats::coef(fit),
      tolerance = 1e-6,
      info = engine
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

test_that("the vglm engine fits the orm proportional odds model by default", {
  skip_on_cran()
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
  skip_on_cran()
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

test_that("mostm() accepts the engine arguments used in the articles and passes them on", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  # The argument names each engine accepts.
  accepted <- list(
    orm = c(
      "opt_method",
      "scale",
      "maxit",
      "penalty",
      "var.penalty",
      "na.action",
      "family"
    ),
    vglm = c(
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
  has_rmsb <- requireNamespace("rmsb", quietly = TRUE)
  if (has_rmsb) {
    accepted$blrm <- c("ppo", "cppo", "iter", "warmup", "chains", "refresh", "method", "backend")
  }
  for (engine in names(accepted)) {
    expect_contains(mostr:::mostm_engine_args(engine), accepted[[engine]])
  }

  # A real fit passes the arguments on to the engine.
  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 6104)
  fit_orm <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    maxit = 30,
    family = "logistic"
  )
  expect_identical(mostr:::bootstrap_model_call(fit_orm)$maxit, 30)
  expect_identical(mostr:::bootstrap_model_call(fit_orm)$family, "logistic")
  fit_vglm <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    engine = "vglm",
    maxit = 30,
    trace = FALSE
  )
  expect_identical(mostr:::bootstrap_model_call(fit_vglm)$maxit, 30)
  expect_equal(fit_vglm$vglm_fit@control$maxit, 30)
  expect_identical(fit_vglm$vglm_fit@control$trace, FALSE)

  # rmsb::blrm() has its own `backend`, which mostm() passes on.
  if (has_rmsb) {
    expect_no_error(
      mostr:::mostm_check_engine_args(
        "blrm",
        list(backend = "cmdstan", iter = 100, ppo = quote(~time)),
        list()
      )
    )
  }
})

test_that("mostm() rejects arguments that every fit must keep", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 20, follow_up_time = 6, seed = 6105)
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
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      engine = "rms"
    ),
    "should be one of"
  )
})

test_that("the blrm engine rejects `cadjust`", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("rmsb")

  # `type` and case weights with blrm are rows F02 and F01 of test-refusals.R;
  # `cadjust` alone is refused by the same check.
  data <- make_test_data(n_patients = 20, follow_up_time = 6, seed = 6106)
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
})

test_that("fits without mostm() provenance must be refitted", {
  skip_on_cran()
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

test_that("mostm() refuses `id_var = NULL`", {
  skip_if_not_installed("rms")

  # A column called `id` must not stand in for the missing argument. A missing
  # `id_var` is row F22 of test-refusals.R; the check runs before the engine
  # is chosen, so one engine is enough.
  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 6112)
  expect_error(
    mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = NULL),
    "`mostm()` needs `id_var`",
    fixed = TRUE
  )
})

test_that("a patient column with another name gives the same results as `id`", {
  skip_on_cran()
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

# Fits each case of a formula table and checks the error message: it must
# contain every phrase in `must` and none in `must_not`. Cases with
# `engine = "blrm"` run only when rmsb is installed.
expect_formula_refusals <- function(cases, data) {
  has_rmsb <- requireNamespace("rmsb", quietly = TRUE)
  for (case in cases) {
    if (identical(case$engine, "blrm") && !has_rmsb) {
      next
    }
    args <- list(
      case$formula,
      data = if (is.null(case$data)) data else case$data,
      id_var = "id",
      engine = case$engine
    )
    args$time_var <- case$time_var
    err <- expect_error(
      do.call(mostm, args),
      case$must[[1L]],
      fixed = TRUE,
      info = case$label
    )
    for (phrase in case$must[-1L]) {
      expect_match(conditionMessage(err), phrase, fixed = TRUE, info = case$label)
    }
    for (phrase in case$must_not) {
      expect_no_match(conditionMessage(err), phrase, fixed = TRUE)
    }
  }
}

test_that("mostm() refuses factor() written inside the formula unless the engine fits it", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6301)
  renamed <- data
  names(renamed)[names(renamed) == "time"] <- "day"

  # rms::orm() cannot build the design matrix for such terms and stops with
  # "subscript out of bounds", so mostm() stops before fitting. The time
  # advice follows `time_var` and is left out for other variables.
  other_variables <- list(
    ordered(y) ~ time + as.factor(tx) + yprev,
    ordered(y) ~ time + base::factor(tx) + yprev,
    ordered(y) ~ time * ordered(tx) + yprev,
    ordered(y) ~ time + yprev + as.ordered(tx)
  )
  cases <- c(
    list(
      list(
        label = "orm, factor(time)",
        engine = "orm",
        formula = ordered(y) ~ factor(time) + tx + yprev,
        must = c(
          "`factor(time)`",
          "`rms::orm()`",
          "column of `data`",
          "first_followup_time",
          "articles/factor-time-orm-sops.html"
        )
      ),
      list(
        label = "orm, factor() of a renamed time variable",
        engine = "orm",
        formula = ordered(y) ~ factor(day) + tx + yprev,
        data = renamed,
        time_var = "day",
        must = c("`factor(day)`", "first_followup_time")
      ),
      list(
        label = "blrm, factor(time)",
        engine = "blrm",
        formula = ordered(y) ~ factor(time) + tx + yprev,
        must = c("`factor(time)`", "`rmsb::blrm()`", "first_followup_time")
      )
    ),
    lapply(other_variables, function(formula) {
      list(
        label = paste("orm,", deparse(formula)),
        engine = "orm",
        formula = formula,
        must = "column of `data`",
        must_not = "first_followup_time"
      )
    })
  )
  expect_formula_refusals(cases, data)

  # A factor column in `data` fits, and `ordered()` on the response is allowed.
  data$tx_factor <- factor(data$tx)
  fit <- mostm(ordered(y) ~ time + tx_factor + yprev, data = data, id_var = "id")
  expect_s3_class(fit, "orm")

  # VGAM fits factor() written inside the formula: one coefficient for each
  # later visit, and SOPs of the stored patients that add up to 1.
  fit <- mostm(
    ordered(y) ~ factor(time) + tx + yprev,
    data = data,
    id_var = "id",
    engine = "vglm"
  )
  expect_s3_class(fit, "robcov_vglm")
  expect_contains(names(stats::coef(fit)), paste0("factor(time)", 2:6))
  out <- sops(fit, times = 1:2, y_levels = 1:6, absorb = 6)
  totals <- tapply(out$estimate, list(out$rowid, out$time), sum)
  expect_equal(dim(totals), c(length(unique(data$id)), 2L))
  expect_equal(as.vector(totals), rep(1, length(totals)), tolerance = 1e-10)
})

test_that("mostm() refuses random-effect terms that its engine cannot fit", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6304)
  data$mre <- as.numeric(data$time > 1)
  # VGAM would otherwise fit the patient ID number as a straight-line
  # predictor. Stopping before the fit keeps orm the same for every rms
  # version: rms 8.2 and later fit these terms, and earlier versions do not
  # have the functions. rmsb::blrm() recognizes only `cluster(id)` as a random
  # intercept and would fit the other forms as ordinary predictors.
  cluster <- function(x) x
  cases <- list(
    list(
      label = "orm, cluster()",
      engine = "orm",
      formula = ordered(y) ~ time + tx + yprev + cluster(id),
      must = c(
        "`cluster(id)`",
        "does not support random effects",
        "random effect is 0",
        "`id_var`",
        "`engine = \"blrm\"`"
      )
    ),
    # Found next to each other, with a package prefix, and inside other calls.
    list(
      label = "orm, two terms",
      engine = "orm",
      formula = ordered(y) ~ time + tx + yprev + cluster(id) + mix_re(mre),
      must = "`cluster(id)` and `mix_re(mre)`"
    ),
    list(
      label = "orm, rms::cluster()",
      engine = "orm",
      formula = ordered(y) ~ time + tx + yprev + rms::cluster(id),
      must = "`rms::cluster(id)`"
    ),
    list(
      label = "orm, rms::mix_re()",
      engine = "orm",
      formula = ordered(y) ~ time + tx + yprev + rms::mix_re(mre),
      must = "`rms::mix_re(mre)`"
    ),
    list(
      label = "orm, cluster() in an interaction",
      engine = "orm",
      formula = ordered(y) ~ time * cluster(id) + yprev,
      must = "`cluster(id)`"
    ),
    list(
      label = "vglm, cluster()",
      engine = "vglm",
      formula = ordered(y) ~ time + tx + yprev + cluster(id),
      must = c(
        "`cluster(id)`",
        "does not support random effects",
        "ordinary predictor",
        "`engine = \"blrm\"`"
      )
    ),
    list(
      label = "blrm, rms::cluster()",
      engine = "blrm",
      formula = ordered(y) ~ time + tx + yprev + rms::cluster(id),
      must = c("`rms::cluster(id)`", "ordinary predictor", "without a package prefix")
    ),
    list(
      label = "blrm, cluster() and mix_re()",
      engine = "blrm",
      formula = ordered(y) ~ time + tx + yprev + cluster(id) + mix_re(mre),
      must = "`mix_re(mre)`",
      must_not = "`cluster(id)` and"
    )
  )
  expect_formula_refusals(cases, data)

  # Only calls count: a column named `cluster` is an ordinary predictor.
  with_column <- data
  with_column$cluster <- with_column$tx
  fit <- mostm(ordered(y) ~ time + cluster + yprev, data = with_column, id_var = "id")
  expect_s3_class(fit, "orm")

  skip_if_not_installed("rmsb")
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
})

test_that("a weighted orm fit does not show rms's warning that its bootstrap ignores weights", {
  # rms warns that orm's own validation and bootstrap ignore weights. mostr
  # never uses those, and its own bootstrap keeps the weights, so the warning
  # would mislead.
  # grid_fit_to() hides all warnings, so fit directly.
  data <- grid_cell_data(grid_cell("orm", weighted = TRUE))
  expect_no_warning(
    mostm(y ~ time + tx + age + yprev, data = data, id_var = "id", weights = w),
    message = "weights are ignored"
  )
})
