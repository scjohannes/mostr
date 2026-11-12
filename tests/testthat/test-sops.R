test_that("soprob_markov matches Hmisc::soprobMarkovOrdm for a proportional-odds VGAM model", {
  skip_if_not_installed("Hmisc")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 18L
  data <- make_test_data(
    n_patients = 60,
    seed = 1234,
    follow_up_time = follow_up
  )
  model <- suppressWarnings(
    mostm(
      ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      engine = "vglm"
    )
  )$vglm_fit

  expected <- Hmisc::soprobMarkovOrdm(
    model,
    data = data[1, ],
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6"
  )
  actual <- soprob_markov(
    model,
    newdata = data[1, ],
    times = seq_len(follow_up),
    absorb = "6"
  )

  expect_probability_array(actual, expected_dim = c(1L, follow_up, 6L))
  expect_equal(actual[1, , ], expected, tolerance = 1e-10)
})

test_that("soprob_markov follows custom partial proportional-odds constraints and treatment interactions", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 18L
  data <- make_test_data(
    n_patients = 70,
    seed = 789,
    follow_up_time = follow_up
  )
  time_spline <- rms::rcs(data$time, 4)
  data$time_lin <- as.vector(time_spline[, 1])
  data$time_nlin_1 <- as.vector(time_spline[, 2])
  data$time_nlin_2 <- as.vector(time_spline[, 3])
  time_covariates <- make_time_covariates(
    data,
    "time",
    "time_lin",
    "time_nlin_1",
    "time_nlin_2"
  )

  # The linear time term gets its own effect on each cut-point through a
  # constraint column 0:4, rather than the identity matrix that
  # `parallel = FALSE` would give.
  n_thresholds <- 5L
  main_constraints <- list(
    "(Intercept)" = diag(n_thresholds),
    "yprev" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = 1),
    "time_lin" = cbind(
      time_lin_1 = 1,
      time_lin_5 = 0:(n_thresholds - 1)
    ),
    "time_nlin_1" = cbind(PO_effect = 1),
    "time_nlin_2" = cbind(PO_effect = 1)
  )
  cases <- list(
    list(
      label = "main effects",
      formula = ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + yprev,
      constraints = main_constraints
    ),
    list(
      label = "treatment-by-time and previous-state-by-time interactions",
      formula = ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) *
        tx +
        yprev +
        yprev:time_lin,
      constraints = c(
        main_constraints,
        list(
          "time_lin:tx" = cbind(PO_effect = 1),
          "time_nlin_1:tx" = cbind(PO_effect = 1),
          "time_nlin_2:tx" = cbind(PO_effect = 1),
          "time_lin:yprev" = cbind(PO_effect = 1)
        )
      )
    )
  )

  baseline <- data[data$time == 1, , drop = FALSE][1, , drop = FALSE]
  # The readable loop sets the time covariates of each visit itself.
  visit_columns <- time_covariates[, c("time_lin", "time_nlin_1", "time_nlin_2")]

  for (case in cases) {
    model <- suppressWarnings(
      mostm(
        case$formula,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
        data = data,
        id_var = "id",
        constraints = case$constraints,
        engine = "vglm"
      )
    )$vglm_fit
    model@call$constraints <- case$constraints

    # Treated and control SOPs both follow the model, so the treatment and
    # interaction coefficients enter every transition.
    for (tx in c(0, 1)) {
      profile <- baseline
      profile$tx <- tx
      label <- paste0(case$label, ", tx = ", tx)
      result <- soprob_markov(
        model,
        newdata = profile,
        times = seq_len(follow_up),
        absorb = "6",
        time_covariates = time_covariates
      )
      expect_probability_array(result, expected_dim = c(1L, follow_up, 6L))
      expect_absorbing_probability_monotone(result, absorbing_state = 6)
      loop <- slow_first_order_sops(
        slow_predictor(model),
        profile,
        times = seq_len(follow_up),
        states = 1:6,
        absorb = 6,
        visit_columns = visit_columns
      )
      expect_equal(unname(result[1, , ]), loop, tolerance = 1e-10, info = label)
    }
  }
})

test_that("soprob_markov normalizes transitions when no state is absorbing", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  follow_up <- 10L
  data <- as.data.frame(make_markov_trajectories(
    n_patients = 80,
    follow_up_time = follow_up,
    treatment_prob = 0.5,
    absorbing_state = NULL,
    allowed_start_state = as.integer(1:6),
    thresholds = c(-2, -1, 0, 1, 2),
    seed = 1000
  ))
  data$yprev <- factor(data$yprev, levels = 1:6)
  baseline <- data[data$time == 1, , drop = FALSE][seq_len(5), , drop = FALSE]
  model <- suppressWarnings(
    mostm(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      engine = "vglm"
    )
  )

  result <- soprob_markov(
    model,
    newdata = baseline,
    times = seq_len(follow_up),
    absorb = NULL
  )

  expect_probability_array(result, expected_dim = c(5L, follow_up, 6L))
})

test_that("avg_sops() equals manual G-computation over individual SOPs", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 6L
  data <- make_test_data(
    n_patients = 40,
    follow_up_time = follow_up,
    seed = 123
  )
  time_covariates <- make_time_covariates(
    data,
    "time",
    "time_lin",
    "time_nlin_1"
  )
  robust_model <- mostm(
    ordered(y) ~ (time_lin + time_nlin_1) * tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    engine = "vglm"
  )
  baseline <- data[data$time == 1, , drop = FALSE]

  baseline_tx1 <- baseline
  baseline_tx1$tx <- 1
  manual_tx1 <- sops(
    robust_model,
    newdata = baseline_tx1,
    times = seq_len(follow_up),
    absorb = "6",
    time_covariates = time_covariates
  )
  manual_tx1 <- aggregate(
    estimate ~ time + state,
    data = manual_tx1,
    FUN = mean
  )
  manual_tx1$tx <- 1

  baseline_tx0 <- baseline
  baseline_tx0$tx <- 0
  manual_tx0 <- sops(
    robust_model,
    newdata = baseline_tx0,
    times = seq_len(follow_up),
    absorb = "6",
    time_covariates = time_covariates
  )
  manual_tx0 <- aggregate(
    estimate ~ time + state,
    data = manual_tx0,
    FUN = mean
  )
  manual_tx0$tx <- 0

  manual <- rbind(manual_tx0, manual_tx1)
  manual <- manual[order(manual$tx, manual$state, manual$time), ]
  rownames(manual) <- NULL
  automated <- avg_sops(
    robust_model,
    newdata = baseline,
    variables = "tx",
    times = seq_len(follow_up),
    absorb = "6",
    time_covariates = time_covariates
  )
  automated <- as.data.frame(automated)
  automated <- automated[order(automated$tx, automated$state, automated$time), ]
  rownames(automated) <- NULL

  expect_equal(manual$estimate, automated$estimate, tolerance = 1e-15)
  expect_equal(manual$tx, automated$tx)
  expect_equal(as.character(manual$state), as.character(automated$state))
  expect_equal(manual$time, automated$time)
})
