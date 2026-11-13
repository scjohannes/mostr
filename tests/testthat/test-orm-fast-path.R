add_test_age <- function(data) {
  data$age <- 45 + (data$id %% 25)
  data
}

# SOPs from the streamed design builder, markov_msm_build() followed by
# markov_msm_run(), which MVN draws and analytical inference use. `args` are
# the arguments of sops(); the settings, absorbing states included, come from
# the fit as in sops().
streamed_builder_sops <- function(fit, profiles, args) {
  settings <- mostr:::resolve_markov_settings(fit)
  components <- mostr:::markov_msm_build(
    model = fit,
    newdata = profiles,
    times = args$times,
    settings = settings
  )
  mostr:::markov_msm_run(
    components,
    mostr::get_effective_coefs(fit),
    times = args$times,
    absorb = settings$absorb
  )
}

# Public sops() uses the batched design builder and is checked against the
# readable loop on every grid set-up in test-sops-native-slow-loop.R, so
# agreeing with sops() checks the streamed builder against that loop too.
test_that("the streamed design builder matches sops() on every model set-up", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  # The streamed builder handles first-order recursion only, so the
  # second-order cells (orm and vglm with linear, rcs and factor time, and vglm
  # with partial odds) are left out.
  grid_cells <- Filter(
    function(cell) cell$order == 1L,
    grid_consistency_cells()
  )
  cases <- lapply(grid_cells, function(cell) {
    args <- grid_public_args(cell)
    list(
      label = grid_label(cell),
      fit = grid_fit(cell),
      profiles = grid_profiles(cell),
      args = args
    )
  })

  # orm with a spline-by-treatment interaction.
  orm_data <- add_test_age(
    make_test_data(n_patients = 60, follow_up_time = 7, seed = 6101)
  )
  orm_fit <- suppressWarnings(mostm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = orm_data,
    id_var = "id",
    absorb = "6"
  ))
  cases[[length(cases) + 1L]] <- list(
    label = "orm, rcs(time) * tx",
    fit = orm_fit,
    profiles = orm_data[!duplicated(orm_data$id), ][1:12, , drop = FALSE],
    args = list(times = 1:5)
  )

  # vglm with a spline of a numeric previous state.
  numeric_data <- make_markov_trajectories(
    n_patients = 100,
    follow_up_time = 7,
    n_states = 10,
    thresholds = seq(-3, 3, length.out = 9),
    allowed_start_state = as.integer(2:9),
    absorbing_state = 10,
    seed = 4242
  ) |>
    prepare_markov_data(absorbing_state = 10, factor_previous = FALSE)
  expect_type(numeric_data$yprev, "double")
  cases[[length(cases) + 1L]] <- list(
    label = "vglm, rcs(yprev) of a numeric previous state",
    fit = suppressWarnings(mostm(
      ordered(y) ~ rms::rcs(time, 4) + tx + rms::rcs(yprev, 6),
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = numeric_data,
      id_var = "id",
      absorb = "10",
      engine = "vglm"
    )),
    profiles = numeric_data[numeric_data$time == 1, , drop = FALSE][1:8, ],
    args = list(times = 1:7)
  )

  # vglm with an inline rcs(time) spline.
  spline_data <- make_markov_trajectories(
    n_patients = 70,
    follow_up_time = 8,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 4041,
    treatment_effect = 0.18
  ) |>
    prepare_markov_data()
  spline_data$age <- 45 + ((spline_data$id * 11) %% 31)
  spline_profiles <- spline_data[spline_data$time == 1, , drop = FALSE][1:12, ]
  cases[[length(cases) + 1L]] <- list(
    label = "vglm, inline rcs(time) * tx",
    fit = suppressWarnings(mostm(
      ordered(y) ~ rms::rcs(time, 4) * tx + yprev + age,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = spline_data,
      id_var = "id",
      absorb = "6",
      engine = "vglm"
    )),
    profiles = spline_profiles,
    args = list(times = 1:8)
  )

  # vglm with factor visit time, no absorbing state, and numeric `times` that
  # name the visit levels.
  factor_data <- withr::with_seed(2026, {
    visits <- as.character(1:4)
    data <- expand.grid(id = 1:80, time = visits, KEEP.OUT.ATTRS = FALSE)
    data <- data[order(data$id, data$time), , drop = FALSE]
    data$time <- factor(data$time, levels = visits)
    data$tx <- stats::rbinom(80, 1, 0.5)[data$id]
    data$yprev <- factor(
      sample(1:4, nrow(data), replace = TRUE),
      levels = as.character(1:4)
    )
    eta <- -0.35 * data$tx + 0.25 * as.integer(data$time) +
      0.2 * as.integer(as.character(data$yprev))
    worse <- stats::plogis(eta - mean(eta))
    data$y <- ordered(
      pmin(4L, pmax(1L, 1L + stats::rbinom(nrow(data), 3L, worse))),
      levels = as.character(1:4)
    )
    data
  })
  factor_fit <- suppressWarnings(mostm(
    y ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = factor_data,
    id_var = "id",
    first_followup_time = "1",
    engine = "vglm"
  ))
  factor_baseline <- factor_data[!duplicated(factor_data$id), , drop = FALSE]
  factor_args <- list(times = 1:4)
  cases[[length(cases) + 1L]] <- list(
    label = "vglm, factor visit time",
    fit = factor_fit,
    profiles = factor_baseline,
    args = factor_args
  )

  for (case in cases) {
    result <- suppressWarnings(do.call(
      sops,
      c(list(case$fit, newdata = case$profiles), case$args)
    ))
    streamed <- streamed_builder_sops(case$fit, case$profiles, case$args)

    n_states <- length(mostr:::markov_model_settings(case$fit)$y_levels)
    expect_probability_array(
      streamed,
      expected_dim = c(nrow(case$profiles), length(case$args$times), n_states)
    )
    matrices <- lapply(
      seq_len(nrow(case$profiles)),
      function(i) matrix(streamed[i, , ], ncol = n_states)
    )
    expect_identical(
      nrow(result),
      nrow(case$profiles) * length(case$args$times) * n_states,
      label = case$label
    )
    # Both builders feed the same recursion, so only rounding separates them.
    expect_lt(
      grid_largest_gap(result, matrices, result$rowid, case$args$times),
      1e-12,
      label = case$label
    )
  }

  # Numeric `times` on a factor-time fit come back as the visit labels.
  direct <- do.call(
    soprob_markov,
    c(list(factor_fit, newdata = factor_baseline), factor_args)
  )
  expect_equal(dim(direct), c(nrow(factor_baseline), 4L, 4L))
  expect_equal(dimnames(direct)[[2]], as.character(1:4))
  individual <- do.call(
    sops,
    c(list(factor_fit, newdata = factor_baseline), factor_args)
  )
  expect_equal(sort(unique(as.character(individual$time))), as.character(1:4))
  averaged <- do.call(
    avg_sops,
    c(
      list(
        factor_fit,
        newdata = factor_baseline,
        variables = list(tx = c(0, 1))
      ),
      factor_args
    )
  )
  expect_equal(
    as.character(attr(averaged, "avg_args")$times),
    as.character(1:4)
  )
  expect_equal(
    unname(streamed_builder_sops(
      factor_fit,
      factor_baseline,
      utils::modifyList(
        factor_args,
        list(times = attr(averaged, "avg_args")$times)
      )
    )),
    unname(direct),
    tolerance = 1e-10
  )
})

test_that("orm SOPs report a state count that differs from the fit and a missing absorbing state", {
  skip_on_cran()
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 6111)
  data <- add_test_age(data)
  # Death, state 6, is left out of `absorb` on purpose, so `mostm()` warns
  # and prediction then has no transitions out of it to work with.
  expect_warning(
    fit <- mostm(y ~ time + tx + yprev + age, data = data, id_var = "id"),
    class = "mostr_warning_outcome_only_states"
  )
  baseline <- data[!duplicated(data$id), ][1:3, , drop = FALSE]

  # The fit always supplies its own states, so a state count that differs
  # from the fit can reach the recursion only through internal settings.
  fewer_states <- mostr:::markov_settings_update(
    mostr:::resolve_markov_settings(fit),
    y_levels = fit$yunique[-length(fit$yunique)]
  )
  expect_snapshot(
    mostr:::predict_markov_sops(
      fit,
      baseline,
      times = 1,
      settings = fewer_states
    ),
    error = TRUE
  )

  expect_snapshot(
    mostr::soprob_markov(
      fit,
      baseline[1L, , drop = FALSE],
      times = 1
    ),
    error = TRUE
  )
})

test_that("execution-plan budgets bound all retained first-order designs", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6112)
  data <- add_test_age(data)
  fit <- suppressWarnings(mostm(
    y ~ time + tx + yprev + age,
    data = data,
    id_var = "id",
    absorb = "6"
  ))
  baseline <- data[!duplicated(data$id), ][1:4, , drop = FALSE]
  args <- list(
    model = fit,
    newdata = baseline,
    times = 1:5
  )
  # The compiled result, before the budget forces the R fallback.
  expect_no_message(expected <- do.call(mostr::soprob_markov, args))

  withr::local_options(mostr.execution_plan_max_bytes = 1)
  for (builder in c("batched", "streamed")) {
    condition <- tryCatch(
      mostr:::compile_sop_execution_plan(
        fit,
        baseline,
        times = 1:5,
        settings = mostr:::resolve_markov_settings(fit),
        builder = builder
      ),
      error = function(e) e
    )
    expect_s3_class(condition, "mostr_execution_plan_too_large")
    expect_gt(condition$required_bytes, condition$max_bytes)
  }
  expect_message(
    actual <- do.call(mostr::soprob_markov, args),
    "Compiled C\\+\\+ SOP calculations were not used"
  )
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("fallback notifications are emitted once per inference call", {
  messages <- testthat::capture_messages({
    mostr:::with_sop_fallback_notification_scope({
      mostr:::notify_sop_r_engine_fallback("first draw failed.")
      mostr:::with_sop_fallback_notification_scope({
        mostr:::notify_sop_r_engine_fallback("second draw failed.")
      })
    })
  })
  expect_length(messages, 1L)
  expect_match(messages, "Compiled C\\+\\+ SOP calculations were not used")
})

test_that("orm fast path agrees with predict.orm for numeric previous state", {
  skip_if_not_installed("rms")

  data <- make_markov_trajectories(
    n_patients = 70,
    follow_up_time = 7,
    seed = 6102,
    absorbing_state = 6
  ) |>
    mostr::prepare_markov_data(
      absorbing_state = 6,
      factor_previous = FALSE
    )

  fit <- suppressWarnings(
    rms::orm(
      y ~ rms::rcs(time, 3) * tx + rms::rcs(yprev, 3),
      data = data,
      x = TRUE,
      y = TRUE
    )
  )
  newdata <- data[1:15, , drop = FALSE]

  direct <- mostr:::predict_orm_response_markov(fit, newdata)
  predicted <- stats::predict(fit, newdata = newdata, type = "fitted.ind")

  expect_equal(unname(direct), unname(predicted), tolerance = 1e-12)
})

test_that("orm scale=TRUE agrees with default scale in SOP workflows", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 6105)
  data <- add_test_age(data)

  fit_default <- mostm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    id_var = "id",
    absorb = "6"
  )
  fit_scaled <- mostm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    id_var = "id",
    absorb = "6",
    scale = TRUE
  )
  baseline <- data[!duplicated(data$id), ]

  predict_scaled <- stats::predict(
    fit_scaled,
    newdata = baseline,
    type = "fitted.ind"
  )
  fast_scaled <- mostr:::predict_orm_response_markov(
    fit_scaled,
    baseline
  )
  expect_equal(unname(fast_scaled), unname(predict_scaled), tolerance = 1e-12)

  sops_default <- mostr::sops(
    fit_default,
    baseline,
    times = 1:4
  )
  sops_scaled <- mostr::sops(
    fit_scaled,
    baseline,
    times = 1:4
  )
  expect_equal(sops_scaled$estimate, sops_default$estimate, tolerance = 1e-10)

  robust_default <- fit_default
  robust_scaled <- fit_scaled
  avg_default <- mostr::avg_sops(
    robust_default,
    variables = "tx",
    times = 1:4
  )
  avg_scaled <- mostr::avg_sops(
    robust_scaled,
    variables = "tx",
    times = 1:4
  )
  expect_equal(avg_scaled$estimate, avg_default$estimate, tolerance = 1e-10)

  withr::local_seed(61051)
  mvn_default <- mostr::inferences(
    avg_default,
    method = "mvn",
    n_draws = 2,
    return_draws = TRUE
  )
  withr::local_seed(61051)
  mvn_scaled <- mostr::inferences(
    avg_scaled,
    method = "mvn",
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(mvn_scaled$conf.low, mvn_default$conf.low, tolerance = 1e-8)
  expect_equal(mvn_scaled$conf.high, mvn_default$conf.high, tolerance = 1e-8)

  withr::local_seed(61052)
  score_default <- mostr:::inferences_score_bootstrap(
    avg_default,
    cluster = data$id,
    n_draws = 2,
    return_draws = TRUE
  )
  withr::local_seed(61052)
  score_scaled <- mostr:::inferences_score_bootstrap(
    avg_scaled,
    cluster = data$id,
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(
    score_scaled$conf.low,
    score_default$conf.low,
    tolerance = 1e-8
  )
  expect_equal(
    score_scaled$conf.high,
    score_default$conf.high,
    tolerance = 1e-8
  )

  bootstrap_full <- function(fit) {
    full <- mostr::avg_sops(
      fit,
      newdata = data,
      variables = "tx",
      times = 1:4
    )
    withr::with_seed(61053, mostr::inferences(
      full,
      method = "bootstrap",
      n_draws = 1,
      return_draws = TRUE
    ))
  }
  boot_default <- bootstrap_full(fit_default)
  boot_scaled <- bootstrap_full(fit_scaled)
  expect_equal(attr(boot_scaled, "method"), "bootstrap")
  expect_equal(attr(boot_scaled, "n_successful"), 1L)
  expect_false(anyNA(boot_scaled$conf.low))
  expect_equal(boot_scaled$conf.low, boot_default$conf.low, tolerance = 1e-6)
  expect_equal(boot_scaled$conf.high, boot_default$conf.high, tolerance = 1e-6)
})

test_that("orm checks the size of an MVN covariance matrix and supports score-bootstrap inference", {
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  data <- make_test_data(n_patients = 55, follow_up_time = 6, seed = 6103)
  data <- add_test_age(data)

  fit <- mostm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    id_var = "id",
    absorb = "6"
  )
  robust_fit <- fit
  baseline <- data[!duplicated(data$id), ]

  avg <- mostr::avg_sops(
    robust_fit,
    variables = "tx",
    times = 1:4
  )

  expect_error(
    mostr::inferences(
      avg,
      method = "mvn",
      vcov = stats::vcov(fit),
      n_draws = 2
    ),
    "Dimension mismatch"
  )

  withr::local_seed(61032)
  score <- mostr:::inferences_score_bootstrap(
    avg,
    cluster = data$id,
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(attr(score, "method"), "score_bootstrap")
  expect_equal(attr(score, "n_successful"), 2L)
  expect_true(any(score$std.error > 0))
})


test_that("orm SOPs and design matrices do not require the outcome column in newdata", {
  skip_on_cran()
  skip_if_not_installed("rms")

  data <- sim_markov("actt2", n_patients = 200, seed = 1)
  first_data <- prepare_markov_data(
    data,
    absorbing_state = 8,
    ordered_response = FALSE
  )
  first_fit <- mostm(
    y ~ tx + yprev + rms::rcs(time, 4),
    data = first_data,
    id_var = "id",
    absorb = 8
  )
  data$ypprev <- stats::ave(data$yprev, data$id, FUN = function(value) {
    c(value[1L], value[-length(value)])
  })
  second_data <- prepare_markov_data(
    data,
    absorbing_state = 8,
    ordered_response = FALSE
  )
  second_data$ypprev <- factor(second_data$ypprev)
  second_fit <- mostm(
    y ~ tx + yprev + ypprev + rms::rcs(time, 4),
    data = second_data,
    id_var = "id",
    p2_var = "ypprev",
    absorb = 8
  )
  states <- function(values, data, column = "yprev") {
    factor(values, levels = levels(data[[column]]))
  }

  run_plan <- function(newdata, builder) {
    plan <- mostr:::compile_sop_execution_plan(
      first_fit,
      newdata,
      times = 1:3,
      settings = mostr:::markov_settings(y_levels = 1:8, absorb = 8),
      builder = builder
    )
    mostr:::run_sop_execution_plan(plan, mostr::get_effective_coefs(first_fit))
  }

  cases <- list(
    list(
      label = "first-order sops()",
      profile = data.frame(
        id = 1,
        tx = 0,
        time = 1,
        yprev = states(5, first_data)
      ),
      outputs = function(newdata) {
        list(sops = sops(
          first_fit,
          newdata = newdata,
          times = 1:3
        )$estimate)
      }
    ),
    list(
      label = "design-matrix callers",
      profile = data.frame(
        id = 1:2,
        tx = 0:1,
        time = 1,
        yprev = states(c(5, 4), first_data)
      ),
      outputs = function(newdata) {
        list(
          predict_orm_response_markov =
            mostr:::predict_orm_response_markov(first_fit, newdata),
          markov_linear_predictor_matrix =
            mostr:::markov_linear_predictor_matrix(first_fit, newdata),
          batched = run_plan(newdata, "batched"),
          streamed = run_plan(newdata, "streamed")
        )
      }
    ),
    list(
      label = "second-order sops()",
      profile = data.frame(
        id = 1,
        tx = 0,
        time = 1,
        yprev = states(5, second_data),
        ypprev = states(5, second_data, "ypprev")
      ),
      outputs = function(newdata) {
        list(sops = sops(
          second_fit,
          newdata = newdata,
          times = 1:3
        )$estimate)
      }
    )
  )

  for (case in cases) {
    profile_with_y <- case$profile
    profile_with_y$y <- NA_real_
    without_y <- case$outputs(case$profile)
    with_y <- case$outputs(profile_with_y)
    for (output in names(without_y)) {
      expect_equal(
        without_y[[output]],
        with_y[[output]],
        tolerance = 1e-12,
        info = paste(case$label, output)
      )
    }
  }
})
