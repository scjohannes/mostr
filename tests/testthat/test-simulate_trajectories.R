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

# The calls that sim_markov() refuses are rows S01 to S16 of test-refusals.R.

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

test_that("presets return reproducible trajectories with their default effect", {
  actt_columns <- c("id", "time", "y", "yprev", "tx")
  # `follow_up_time = NULL` leaves the argument out, so the default of 28
  # days applies.
  presets <- list(
    list(
      label = "ACTT-1",
      preset = "actt1",
      n_patients = 80,
      follow_up_time = 6,
      seed = 1987,
      default_effect = -0.0730,
      columns = actt_columns,
      states = 1:8,
      start_states = 4:7,
      absorbing = 8
    ),
    list(
      label = "ACTT-2",
      preset = "actt2",
      n_patients = 80,
      follow_up_time = 6,
      seed = 1987,
      default_effect = -0.0730,
      columns = actt_columns,
      states = 1:8,
      start_states = 4:7,
      absorbing = 8
    ),
    list(
      label = "ACTT-2, 60 days",
      preset = "actt2",
      n_patients = 20,
      follow_up_time = 60,
      seed = 123,
      default_effect = -0.0730,
      columns = actt_columns,
      states = 1:8,
      start_states = NULL,
      absorbing = 8
    ),
    list(
      label = "VIOLET",
      preset = "violet",
      n_patients = 30,
      follow_up_time = NULL,
      seed = 1904,
      default_effect = log(0.75),
      columns = c("id", "time", "y", "yprev", "age", "sofa", "tx"),
      states = 1:4,
      start_states = NULL,
      absorbing = NULL
    )
  )

  for (case in presets) {
    simulate <- function(...) {
      args <- list(case$preset, n_patients = case$n_patients, seed = case$seed, ...)
      args$follow_up_time <- case$follow_up_time
      do.call(sim_markov, args)
    }
    follow_up_time <- if (is.null(case$follow_up_time)) 28 else case$follow_up_time
    traj <- simulate()

    expect_trajectory_contract(
      traj,
      expected_cols = case$columns,
      n_patients = case$n_patients,
      follow_up_time = follow_up_time,
      states = case$states
    )
    expect_identical(traj, simulate(), info = case$label)
    expect_identical(
      traj,
      simulate(treatment_effect = case$default_effect),
      info = case$label
    )
    if (!is.null(case$start_states)) {
      expect_setequal(unique(traj$yprev[traj$time == 1]), case$start_states)
    }
    if (!is.null(case$absorbing)) {
      expect_absorbing_state_sticky(traj, absorbing_state = case$absorbing)
      absorbed <- traj$yprev == case$absorbing
      expect_equal(
        traj$y[absorbed],
        rep(as.integer(case$absorbing), sum(absorbed)),
        info = case$label
      )
    }
  }
  expect_equal(sum(violet_baseline$n), 1351)
})

test_that("ACTT-1 and ACTT-2 presets share their slopes and sample their own starting states", {
  actt1 <- mostr:::actt1_markov_parameters()
  actt2 <- mostr:::actt2_markov_parameters()

  # ACTT-1 is the ACTT-2 model with ACTT-1 starting states and death threshold.
  expect_equal(actt1$baseline_states, 4:7)
  expect_equal(actt1$baseline_probabilities, c(138, 435, 193, 285) / 1051)
  expect_identical(actt1$extra_params, actt2$extra_params)
  expect_identical(actt1$orm_intercepts[-7], actt2$orm_intercepts[-7])
  expect_equal(actt1$orm_intercepts[[7]], -21.4071)

  parameters <- list(actt1 = actt1, actt2 = actt2)
  for (preset in names(parameters)) {
    traj <- sim_markov(preset, n_patients = 10000, follow_up_time = 1, seed = 4815)
    initial <- traj[traj$time == 1, , drop = FALSE]
    observed <- prop.table(table(factor(initial$yprev, levels = 4:7)))
    expect_equal(
      as.numeric(observed),
      parameters[[preset]]$baseline_probabilities,
      tolerance = 0.02,
      info = preset
    )
    expect_equal(mean(initial$tx), 0.5, tolerance = 0.03, info = preset)
  }
})

test_that("ACTT presets match their control-arm calibration targets", {
  actt1 <- mostr:::actt1_markov_parameters()
  actt2 <- mostr:::actt2_markov_parameters()
  targets <- list(
    list(
      label = "ACTT-1 kills 18% of the control arm by day 28",
      parameters = actt1,
      days = 28,
      state = 8,
      target = 0.18,
      tolerance = 1e-4
    ),
    list(
      label = "ACTT-2 places 9.3% of the control arm in state 2 on day 14",
      parameters = actt2,
      days = 14,
      state = 2,
      target = 0.093,
      tolerance = 1e-2
    )
  )

  for (case in targets) {
    occupancy <- exact_actt_control_occupancy(case$parameters, days = case$days)
    expect_equal(sum(occupancy), 1, info = case$label)
    expect_equal(
      occupancy[[case$state]],
      case$target,
      tolerance = case$tolerance,
      info = case$label
    )
  }

  # The unchanged ACTT-2 death threshold gives about 11% with ACTT-1 patients.
  unchanged <- actt1
  unchanged$orm_intercepts <- actt2$orm_intercepts
  expect_equal(
    exact_actt_control_occupancy(unchanged, days = 28)[[8]],
    0.1107,
    tolerance = 1e-3
  )
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

test_that("ACTT-2 linear predictor preserves shared orm slopes and merges states 1 and 2", {
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

  # Transitions from states 1 and 2 are identical, also for a numeric `yprev`.
  numeric_lp <- function(yprev, day) {
    mostr:::actt2_markov_lp(yprev, day, 1, -0.0730, parameters$extra_params)
  }
  for (day in c(1, 14, 28)) {
    expect_equal(numeric_lp(1, day), numeric_lp(2, day), info = paste("day", day))
  }
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

test_that("sim_markov warns about absent covariates and still simulates", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1))
  lp <- function(yprev, t, treatment_effect, extra_params) 0

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
  # A linear predictor of +10 or -10 all but fixes the next state, so the
  # states show that lp_function saw each patient's own tx.
  expect_equal(result$y[result$tx == 1], c(2, 2))
  expect_equal(result$y[result$tx == 0], c(1, 1))
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
