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
      absorb = "6",
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
    times = seq_len(follow_up)
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
  # The linear time term gets its own effect on each cut-point through a
  # constraint column 0:4, rather than the identity matrix that
  # `parallel = FALSE` would give. Time enters as `time` and `log(time)`,
  # inline terms that prediction evaluates at each visit.
  n_thresholds <- 5L
  main_constraints <- list(
    "(Intercept)" = diag(n_thresholds),
    "yprev" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = 1),
    "time" = cbind(
      time_1 = 1,
      time_5 = 0:(n_thresholds - 1)
    ),
    "log(time)" = cbind(PO_effect = 1)
  )
  cases <- list(
    list(
      label = "main effects",
      formula = ordered(y) ~ time + log(time) + tx + yprev,
      constraints = main_constraints
    ),
    list(
      label = "treatment-by-time and previous-state-by-time interactions",
      formula = ordered(y) ~ (time + log(time)) *
        tx +
        yprev +
        yprev:time,
      constraints = c(
        main_constraints,
        list(
          "time:tx" = cbind(PO_effect = 1),
          "log(time):tx" = cbind(PO_effect = 1),
          "time:yprev" = cbind(PO_effect = 1)
        )
      )
    )
  )

  baseline <- data[data$time == 1, , drop = FALSE][1, , drop = FALSE]

  for (case in cases) {
    model <- suppressWarnings(
      mostm(
        case$formula,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time),
        data = data,
        id_var = "id",
        absorb = "6",
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
        times = seq_len(follow_up)
      )
      expect_probability_array(result, expected_dim = c(1L, follow_up, 6L))
      expect_absorbing_probability_monotone(result, absorbing_state = 6)
      loop <- slow_first_order_sops(
        slow_predictor(model),
        profile,
        times = seq_len(follow_up),
        states = 1:6,
        absorb = 6
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
    times = seq_len(follow_up)
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
  robust_model <- mostm(
    ordered(y) ~ rms::rcs(time, 3) * tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    absorb = "6",
    engine = "vglm"
  )
  baseline <- data[data$time == 1, , drop = FALSE]

  baseline_tx1 <- baseline
  baseline_tx1$tx <- 1
  manual_tx1 <- sops(
    robust_model,
    newdata = baseline_tx1,
    times = seq_len(follow_up)
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
    times = seq_len(follow_up)
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
    times = seq_len(follow_up)
  )
  automated <- as.data.frame(automated)
  automated <- automated[order(automated$tx, automated$state, automated$time), ]
  rownames(automated) <- NULL

  expect_equal(manual$estimate, automated$estimate, tolerance = 1e-15)
  expect_equal(manual$tx, automated$tx)
  expect_equal(as.character(manual$state), as.character(automated$state))
  expect_equal(manual$time, automated$time)
})
