# Exact control-arm state occupancy probabilities of an ACTT preset, from the
# product of its daily transition matrices. Test-only reference for the
# preset calibrations.
exact_actt_control_occupancy <- function(parameters, days) {
  occupancy <- c(0, 0, 0, parameters$baseline_probabilities, 0)
  for (day in seq_len(days)) {
    transition <- matrix(0, nrow = 8, ncol = 8)
    transition[8, 8] <- 1
    for (previous in 1:7) {
      eta <- mostr:::actt2_markov_lp(
        yprev = factor(previous, levels = 1:8),
        t = day,
        tx = 0,
        treatment_effect = 0,
        extra_params = parameters$extra_params
      )
      # Pr(Y >= 2), ..., Pr(Y >= 8)
      cumulative <- stats::plogis(parameters$orm_intercepts + eta)
      transition[previous, ] <- c(1, cumulative) - c(cumulative, 0)
    }
    occupancy <- drop(occupancy %*% transition)
  }
  occupancy
}

test_that("sim_markov requires a valid preset", {
  expect_error(
    sim_markov(n_patients = 10),
    "`preset` is required",
    fixed = TRUE
  )
  expect_error(
    sim_markov("actt3", n_patients = 10),
    '`preset` must be one of "actt1", "actt2", "violet", or NULL',
    fixed = TRUE
  )
  expect_error(
    sim_markov(c("actt1", "actt2")),
    "`preset` must be one of",
    fixed = TRUE
  )
})

test_that("sim_markov validates preset inputs", {
  expect_snapshot(error = TRUE, sim_markov("actt2", n_patients = 0))
  expect_snapshot(error = TRUE, sim_markov("actt2", treatment_prob = 2))
  expect_snapshot(error = TRUE, sim_markov("actt2", treatment_effect = Inf))
  expect_snapshot(
    error = TRUE,
    sim_markov("actt2", treatment_effect_decay = -1)
  )
  expect_snapshot(error = TRUE, sim_markov("actt2", follow_up_time = 0))
})

test_that("sim_markov rejects patients from two sources", {
  baseline <- data.frame(id = 1:2, yprev = c(5L, 6L), tx = c(0L, 1L))

  expect_error(
    sim_markov("actt2", n_patients = 2, baseline_data = baseline),
    "Supply either `baseline_data` or `n_patients` and `treatment_prob`",
    fixed = TRUE
  )
  expect_error(
    sim_markov("actt2", treatment_prob = 1, baseline_data = baseline),
    "Supply either `baseline_data` or `n_patients` and `treatment_prob`",
    fixed = TRUE
  )
})

test_that("sim_markov requires a complete model without a preset", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1), tx = c(0, 1))
  lp <- function(yprev, t, tx, treatment_effect, extra_params) 0

  expect_error(
    sim_markov(NULL, intercepts = 0, lp_function = lp),
    "With `preset = NULL`, supply `baseline_data`.",
    fixed = TRUE
  )
  expect_error(
    sim_markov(NULL, baseline_data = baseline),
    "With `preset = NULL`, supply `intercepts`, `lp_function`.",
    fixed = TRUE
  )
  expect_error(
    sim_markov(
      NULL,
      n_patients = 2,
      baseline_data = baseline,
      intercepts = 0,
      lp_function = lp
    ),
    "Supply either `baseline_data` or `n_patients`",
    fixed = TRUE
  )
  expect_error(
    sim_markov(NULL, n_patients = 2, intercepts = 0, lp_function = lp),
    "`n_patients` and `treatment_prob` sample patients from a preset",
    fixed = TRUE
  )
})

test_that("ACTT-1 preset returns reproducible trajectories", {
  traj <- sim_markov("actt1", n_patients = 80, follow_up_time = 6, seed = 1987)
  repeated <- sim_markov(
    "actt1",
    n_patients = 80,
    follow_up_time = 6,
    treatment_effect = -0.0730,
    seed = 1987
  )

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "time", "y", "yprev", "tx"),
    n_patients = 80,
    follow_up_time = 6,
    states = 1:8
  )
  expect_identical(traj, repeated)
  expect_setequal(unique(traj$yprev[traj$time == 1]), 4:7)
  expect_absorbing_state_sticky(traj, absorbing_state = 8)
})

test_that("ACTT-1 preset is the ACTT-2 model with ACTT-1 starting states", {
  actt1 <- mostr:::actt1_markov_parameters()
  actt2 <- mostr:::actt2_markov_parameters()

  expect_equal(actt1$baseline_states, 4:7)
  expect_equal(actt1$baseline_probabilities, c(138, 435, 193, 285) / 1051)
  expect_identical(actt1$extra_params, actt2$extra_params)
  expect_identical(actt1$orm_intercepts[-7], actt2$orm_intercepts[-7])
  expect_equal(actt1$orm_intercepts[[7]], -21.4071)

  traj <- sim_markov("actt1", n_patients = 10000, follow_up_time = 1, seed = 4815)
  initial <- traj[traj$time == 1, , drop = FALSE]
  observed <- prop.table(table(factor(initial$yprev, levels = 4:7)))
  expect_equal(
    as.numeric(observed),
    actt1$baseline_probabilities,
    tolerance = 0.02
  )
})

test_that("ACTT-1 preset kills 18% of the control arm by day 28", {
  actt1 <- mostr:::actt1_markov_parameters()
  occupancy <- exact_actt_control_occupancy(actt1, days = 28)

  expect_equal(sum(occupancy), 1)
  expect_equal(occupancy[[8]], 0.18, tolerance = 1e-4)

  # The unchanged ACTT-2 death threshold gives about 11% with these patients.
  unchanged <- actt1
  unchanged$orm_intercepts <- mostr:::actt2_markov_parameters()$orm_intercepts
  expect_equal(
    exact_actt_control_occupancy(unchanged, days = 28)[[8]],
    0.1107,
    tolerance = 1e-3
  )
})

test_that("ACTT-2 preset places 9.3% of the control arm in state 2 on day 14", {
  occupancy <- exact_actt_control_occupancy(
    mostr:::actt2_markov_parameters(),
    days = 14
  )

  expect_equal(occupancy[[2]], 0.093, tolerance = 1e-2)
})

test_that("ACTT-2 preset returns reproducible trajectories", {
  traj <- sim_markov("actt2", n_patients = 80, follow_up_time = 6, seed = 1987)
  repeated <- sim_markov(
    "actt2",
    n_patients = 80,
    follow_up_time = 6,
    treatment_effect = -0.0730,
    seed = 1987
  )

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "time", "y", "yprev", "tx"),
    n_patients = 80,
    follow_up_time = 6,
    states = 1:8
  )
  expect_identical(traj, repeated)
  expect_setequal(unique(traj$yprev[traj$time == 1]), 4:7)
  expect_absorbing_state_sticky(traj, absorbing_state = 8)
})

test_that("ACTT-2 preset samples its starting-state distribution", {
  traj <- sim_markov("actt2", n_patients = 10000, follow_up_time = 1, seed = 4815)
  initial <- traj[traj$time == 1, , drop = FALSE]
  observed <- prop.table(table(factor(initial$yprev, levels = 4:7)))
  expected <- mostr:::actt2_markov_parameters()

  expect_equal(
    as.numeric(observed),
    expected$baseline_probabilities,
    tolerance = 0.02
  )
  expect_equal(mean(initial$tx), 0.5, tolerance = 0.03)
})

test_that("ACTT presets code treatment as one", {
  for (preset in c("actt1", "actt2")) {
    control <- sim_markov(
      preset,
      n_patients = 20,
      treatment_prob = 0,
      follow_up_time = 1,
      seed = 91
    )
    treated <- sim_markov(
      preset,
      n_patients = 20,
      treatment_prob = 1,
      follow_up_time = 1,
      seed = 91
    )

    expect_equal(unique(control$tx), 0)
    expect_equal(unique(treated$tx), 1)
  }
})

test_that("ACTT restricted cubic spline reproduces the supplied basis", {
  days <- c(1, 3, 7, 22, 28)
  expected <- rbind(
    c(1, 0, 0, 0, 0),
    c(3, 0.0016, 0, 0, 0),
    c(7, 0.2, 0.0016, 0, 0),
    c(22, 12.7933333333, 6.548, 2.1253333333, 0.3426666667),
    c(28, 25.84, 15.12, 6.4, 1.76)
  )

  expect_equal(
    unname(mostr:::actt_markov_rcs_basis(days)),
    unname(expected),
    tolerance = 1e-9
  )
})

test_that("ACTT spline basis agrees with Hmisc", {
  skip_if_not_installed("Hmisc")
  actual <- mostr:::actt_markov_rcs_basis(1:28)
  expected <- Hmisc::rcspline.eval(
    1:28,
    knots = c(2, 6, 11, 16, 21, 27),
    inclx = TRUE
  )

  expect_equal(dim(actual), dim(expected))
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
})

test_that("ACTT-2 linear predictor preserves shared orm slopes", {
  parameters <- mostr:::actt2_markov_parameters()
  lp <- function(yprev, day, tx, treatment_effect = -0.0730) {
    mostr:::actt2_markov_lp(
      yprev = factor(yprev, levels = 1:8),
      t = day,
      tx = tx,
      treatment_effect = treatment_effect,
      extra_params = parameters$extra_params
    )
  }

  expect_equal(lp(4, 7, 1) - lp(4, 7, 0), -0.0730)
  expect_equal(lp(4, 7, 1, treatment_effect = 0), lp(4, 7, 0))
  expect_equal(lp(5, 7, 1) - lp(4, 7, 1), 9.9809 - 7.1412)

  basis <- mostr:::actt_markov_rcs_basis(7)
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  expect_equal(
    unname(lp(4, 7, 0) - parameters$extra_params["yprev=4"]),
    drop(basis %*% parameters$extra_params[day_names])
  )
})

test_that("ACTT-2 reversed thresholds reproduce orm state probabilities", {
  parameters <- mostr:::actt2_markov_parameters()
  expect_equal(
    parameters$orm_intercepts,
    c(-0.6465, -4.9314, -5.2230, -7.4169, -11.5321, -15.2262, -22.0693)
  )

  eta <- mostr:::actt2_markov_lp(
    yprev = factor(5, levels = 1:8),
    t = 10,
    tx = 0,
    treatment_effect = -0.0730,
    extra_params = parameters$extra_params
  )
  orm_cumulative <- stats::plogis(parameters$orm_intercepts + eta)
  orm_probabilities <- c(
    1 - orm_cumulative[1],
    -diff(orm_cumulative),
    orm_cumulative[length(orm_cumulative)]
  )
  simulator_cumulative <- stats::plogis(rev(parameters$orm_intercepts) + eta)
  simulator_columns <- c(
    simulator_cumulative[1],
    diff(simulator_cumulative),
    1 - simulator_cumulative[length(simulator_cumulative)]
  )
  simulator_probabilities <- rev(simulator_columns)

  expect_equal(simulator_probabilities, orm_probabilities, tolerance = 1e-15)
  expect_gte(min(simulator_probabilities), 0)
  expect_equal(sum(simulator_probabilities), 1)
})

test_that("ACTT-2 transitions from states 1 and 2 are identical", {
  parameters <- mostr:::actt2_markov_parameters()
  lp <- function(yprev, day) {
    mostr:::actt2_markov_lp(yprev, day, 1, -0.0730, parameters$extra_params)
  }

  for (day in c(1, 14, 28)) {
    expect_equal(lp(1, day), lp(2, day))
  }
})

test_that("ACTT presets match a per-patient call of their linear predictor", {
  # Without the attribute, the engine calls the function once per patient.
  per_patient_lp <- function(yprev, t, tx, treatment_effect, extra_params) {
    mostr:::actt2_markov_lp(yprev, t, tx, treatment_effect, extra_params)
  }
  expect_true(attr(mostr:::actt2_markov_lp, "mostr_vectorized_lp"))
  expect_null(attr(per_patient_lp, "mostr_vectorized_lp"))

  baseline <- data.frame(id = 1:16, yprev = rep(1:8, 2), tx = rep(0:1, 8))
  for (preset in c("actt1", "actt2")) {
    simulate <- function(lp_function = NULL) {
      sim_markov(
        preset,
        baseline_data = baseline,
        follow_up_time = 30,
        treatment_effect = log(0.8),
        treatment_effect_decay = 0.1,
        lp_function = lp_function,
        seed = 2026
      )
    }

    expect_identical(simulate(), simulate(per_patient_lp))
  }
})

test_that("presets simulate exactly the patients in baseline_data", {
  baseline <- data.frame(id = c(10L, 20L), yprev = c(5L, 8L), tx = c(0L, 1L))

  for (preset in c("actt1", "actt2")) {
    traj <- sim_markov(
      preset,
      baseline_data = baseline,
      follow_up_time = 3,
      seed = 1
    )

    expect_equal(nrow(traj), 6)
    expect_equal(unique(traj$id), c(10L, 20L))
    expect_equal(traj$tx, rep(c(0L, 1L), each = 3))
    expect_equal(traj$y[traj$id == 20L], rep(8, 3))
  }
})

test_that("a supplied model argument replaces only that part of a preset", {
  baseline <- data.frame(id = 1:30, yprev = rep(4:7, length.out = 30), tx = 0L)
  actt1 <- mostr:::actt1_markov_parameters()

  expect_identical(
    sim_markov(
      "actt2",
      baseline_data = baseline,
      intercepts = rev(actt1$orm_intercepts),
      follow_up_time = 10,
      seed = 7
    ),
    sim_markov("actt1", baseline_data = baseline, follow_up_time = 10, seed = 7)
  )
})

test_that("sim_markov passes the decayed treatment effect to lp_function", {
  seen <- new.env()
  record_lp <- function(yprev, t, tx, treatment_effect, extra_params) {
    assign(as.character(t), treatment_effect, envir = seen)
    0
  }
  baseline <- data.frame(id = 1:2, yprev = c(1, 1), tx = c(0, 1))

  sim_markov(
    NULL,
    baseline_data = baseline,
    follow_up_time = 5,
    intercepts = 0,
    lp_function = record_lp,
    treatment_effect = 0.5,
    treatment_effect_decay = 0.2,
    absorbing_states = integer(),
    seed = 1
  )
  recorded <- vapply(as.character(1:5), get, numeric(1), envir = seen)
  expect_equal(unname(recorded), 0.5 * exp(-0.2 * (0:4)))

  sim_markov(
    NULL,
    baseline_data = baseline,
    follow_up_time = 2,
    intercepts = 0,
    lp_function = record_lp,
    absorbing_states = integer(),
    seed = 1
  )
  expect_equal(get("1", envir = seen), 0)
})

test_that("sim_markov validates a user model and warns about absent covariates", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1))
  lp <- function(yprev, t, treatment_effect, extra_params) 0

  expect_error(
    sim_markov(NULL, baseline_data = list(), intercepts = 0, lp_function = lp),
    "baseline_data must be a data frame",
    fixed = TRUE
  )
  expect_error(
    sim_markov(
      NULL,
      baseline_data = data.frame(id = 1),
      intercepts = 0,
      lp_function = lp
    ),
    "baseline_data must contain columns: yprev",
    fixed = TRUE
  )
  expect_error(
    sim_markov(
      NULL,
      baseline_data = baseline,
      intercepts = character(),
      lp_function = lp
    ),
    "intercepts must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    sim_markov(
      NULL,
      baseline_data = baseline,
      intercepts = 0,
      lp_function = "nope"
    ),
    "lp_function must be a function",
    fixed = TRUE
  )

  expect_warning(
    result <- sim_markov(
      NULL,
      baseline_data = baseline,
      follow_up_time = 2,
      intercepts = 0,
      lp_function = lp,
      covariate_names = "missing",
      seed = 1
    ),
    "None of the specified covariate_names found"
  )

  expect_equal(nrow(result), 4)
  expect_equal(names(result), c("id", "time", "y", "yprev"))
})

test_that("sim_markov carries forward once all patients are absorbing", {
  baseline <- data.frame(id = 1:3, yprev = c(2, 2, 2), tx = 0)

  result <- sim_markov(
    NULL,
    baseline_data = baseline,
    follow_up_time = 4,
    intercepts = -Inf,
    lp_function = function(yprev, t, tx, treatment_effect, extra_params) 0,
    absorbing_states = 2,
    seed = 1,
    covariate_names = "tx"
  )

  expect_equal(unique(result$y), 2)
  expect_equal(nrow(result), 12)
})

test_that("sim_markov passes tx even when not requested as a covariate", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1), tx = c(0, 1))

  expect_no_warning(
    result <- sim_markov(
      NULL,
      baseline_data = baseline,
      follow_up_time = 2,
      intercepts = 0,
      lp_function = function(yprev, t, tx, treatment_effect, extra_params) {
        ifelse(tx == 1, 10, -10)
      },
      covariate_names = character(0),
      seed = 1
    )
  )

  expect_equal(nrow(result), 4)
  expect_true(all(result$tx %in% c(0, 1)))
})

test_that("sim_markov preserves character patient IDs", {
  baseline <- data.frame(id = c("b", "a"), yprev = c(1, 1), tx = c(0, 1))

  result <- sim_markov(
    NULL,
    baseline_data = baseline,
    follow_up_time = 2,
    intercepts = 0,
    lp_function = function(yprev, t, tx, treatment_effect, extra_params) 0,
    covariate_names = "tx",
    seed = 1
  )

  expect_identical(result$id, rep(c("b", "a"), each = 2))
  expect_equal(result$time, rep(1:2, 2))
  expect_equal(result$tx, rep(c(0, 1), each = 2))
})

test_that("sim_markov supports threshold-specific linear predictors", {
  baseline <- data.frame(id = 1:4, yprev = c(1, 1, 1, 1), tx = c(0, 1, 0, 1))
  simulate <- function(lp_function) {
    sim_markov(
      NULL,
      baseline_data = baseline,
      follow_up_time = 3,
      intercepts = c(-1, 1),
      lp_function = lp_function,
      treatment_effect = 0.5,
      absorbing_states = integer(),
      covariate_names = "tx",
      seed = 10
    )
  }

  scalar <- simulate(function(yprev, t, tx, treatment_effect, extra_params) {
    tx * treatment_effect
  })
  scalar_as_vector <- simulate(
    function(yprev, t, tx, treatment_effect, extra_params) {
      rep(tx * treatment_effect, 2)
    }
  )
  expect_identical(scalar, scalar_as_vector)

  ppo <- simulate(function(yprev, t, tx, treatment_effect, extra_params) {
    tx * treatment_effect * c(0, 1)
  })
  expect_named(ppo, c("id", "time", "y", "yprev", "tx"))
  expect_equal(nrow(ppo), 12)
  expect_true(all(ppo$y %in% 1:3))
})

test_that("sim_markov validates threshold-specific linear predictors", {
  baseline <- data.frame(id = 1, yprev = 1, tx = 1)
  simulate <- function(intercepts, lp_function) {
    sim_markov(
      NULL,
      baseline_data = baseline,
      follow_up_time = 1,
      intercepts = intercepts,
      lp_function = lp_function,
      absorbing_states = integer(),
      covariate_names = "tx",
      seed = 1
    )
  }

  expect_error(
    simulate(
      c(-1, 0, 1),
      function(yprev, t, tx, treatment_effect, extra_params) c(0, 0)
    ),
    "lp_function must return either length 1 or length\\(intercepts\\)",
    fixed = FALSE
  )
  expect_error(
    simulate(
      c(-1, 0, 1),
      function(yprev, t, tx, treatment_effect, extra_params) NA_real_
    ),
    "lp_function must return finite numeric values",
    fixed = TRUE
  )
  expect_error(
    simulate(
      c(0, 1),
      function(yprev, t, tx, treatment_effect, extra_params) c(5, -5)
    ),
    "create decreasing cumulative probabilities",
    fixed = TRUE
  )
})

test_that("ACTT-2 Markov trajectories support 60-day follow-up", {
  data <- sim_markov("actt2", n_patients = 20, follow_up_time = 60, seed = 123)
  expect_equal(nrow(data), 1200)
  expect_equal(sort(unique(data$time)), 1:60)
  expect_equal(data$y[data$yprev == 8], rep(8L, sum(data$yprev == 8)))
  expect_identical(
    data,
    sim_markov("actt2", n_patients = 20, follow_up_time = 60, seed = 123)
  )
})
