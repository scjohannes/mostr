test_that("sim_actt1_markov returns reproducible ACTT-1 trajectories", {
  traj <- sim_actt1_markov(
    n_patients = 80,
    follow_up_time = 6,
    seed = 1987
  )
  repeated <- sim_actt1_markov(
    n_patients = 80,
    follow_up_time = 6,
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

  absorbing <- sim_actt1_markov(
    n_patients = 1,
    treatment_prob = 1,
    follow_up_time = 4,
    baseline_data = data.frame(id = 1L, yprev = 8L, tx = 1L),
    seed = 1
  )
  expect_equal(absorbing$y, rep(8L, 4))
  expect_equal(absorbing$yprev, rep(8L, 4))
})

test_that("sim_actt1_markov samples the reported baseline counts", {
  traj <- sim_actt1_markov(
    n_patients = 10000,
    follow_up_time = 1,
    seed = 4815
  )
  initial <- traj[traj$time == 1, , drop = FALSE]
  observed <- prop.table(table(factor(initial$yprev, levels = 4:7)))
  parameters <- mostr:::actt1_markov_parameters()

  expect_equal(
    parameters$baseline_probabilities,
    c(138, 435, 193, 285) / 1051
  )
  expect_equal(parameters$baseline_states, 4:7)
  expect_equal(sum(parameters$baseline_probabilities), 1)
  expect_equal(
    as.numeric(observed),
    parameters$baseline_probabilities,
    tolerance = 0.02
  )
})

test_that("sim_actt1_markov codes one as Remdesivir", {
  placebo <- sim_actt1_markov(
    n_patients = 20,
    treatment_prob = 0,
    follow_up_time = 1,
    seed = 91
  )
  remdesivir <- sim_actt1_markov(
    n_patients = 20,
    treatment_prob = 1,
    follow_up_time = 1,
    seed = 91
  )

  expect_equal(unique(placebo$tx), 0)
  expect_equal(unique(remdesivir$tx), 1)
})

test_that("ACTT-1 parameters preserve selected posterior means", {
  parameters <- mostr:::actt1_markov_parameters()

  expect_equal(
    parameters$orm_intercepts,
    c(-3.3535, -3.5896, -3.8395, -6.0379, -9.7465, -12.3361, -18.8334)
  )
  expect_equal(
    unname(parameters$extra_params[c(
      "yprev=2",
      "yprev=3",
      "yprev=4",
      "yprev=5",
      "yprev=6",
      "yprev=7"
    )]),
    c(-0.0402, 3.9042, 5.0754, 8.0171, 11.5073, 16.1963)
  )
  expect_equal(
    unname(parameters$extra_params[c(
      "day",
      "day'",
      "day''",
      "day'''",
      "day''''"
    )]),
    c(-0.1110, 1.6444, -3.1537, 1.2387, 1.1910)
  )
  expect_identical(
    names(parameters$extra_params),
    c(
      "day",
      "day'",
      "day''",
      "day'''",
      "day''''",
      paste0("yprev=", 1:7),
      "day x f(y)"
    )
  )
  expect_equal(parameters$extra_params[["day x f(y)"]], -0.0188)
})

test_that("ACTT-1 linear predictor combines time state and treatment", {
  parameters <- mostr:::actt1_markov_parameters()
  extra_params <- c(parameters$extra_params, treatment_effect_decay = 0)
  lp <- function(yprev, day, tx, parameter = -0.0657, decay = 0) {
    extra_params[["treatment_effect_decay"]] <- decay
    mostr:::actt1_markov_lp(
      yprev = factor(yprev, levels = 1:8),
      t = day,
      tx = tx,
      parameter = parameter,
      extra_params = extra_params
    )
  }
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  basis <- mostr:::actt_markov_rcs_basis(7)

  expect_equal(lp(4, 7, 1) - lp(4, 7, 0), rep(-0.0657, 7))
  expect_equal(lp(4, 1, 1) - lp(4, 1, 0), rep(-0.0657, 7))
  expect_equal(
    lp(4, 15, 1, decay = 0.2) - lp(4, 15, 0, decay = 0.2),
    rep(-0.0657 * exp(-0.2 * 14), 7)
  )
  expect_equal(
    lp(5, 7, 0) - lp(4, 7, 0),
    rep(8.0171 - 5.0754, 7)
  )
  expect_equal(lp(4, 7, 1, parameter = 0), lp(4, 7, 0))
  expect_equal(
    unname(lp(4, 7, 0) - parameters$extra_params["yprev=4"]),
    drop(basis %*% parameters$extra_params[day_names]) +
      parameters$extra_params[["day x f(y)"]] * 7 * rev(2:8)
  )
  expect_equal(
    mostr:::actt1_markov_lp(
      yprev = factor(4, levels = 1:8),
      t = 7,
      tx = 1,
      parameter = -0.0657,
      extra_params = parameters$extra_params
    ),
    lp(4, 7, 1)
  )
})

test_that("ACTT-1 reversed thresholds reproduce Bayesian model probabilities", {
  parameters <- mostr:::actt1_markov_parameters()
  eta <- mostr:::actt1_markov_lp(
    yprev = factor(5, levels = 1:8),
    t = 10,
    tx = 1,
    parameter = -0.0657,
    extra_params = c(parameters$extra_params, treatment_effect_decay = 0)
  )
  model_cumulative <- stats::plogis(parameters$orm_intercepts + rev(eta))
  model_probabilities <- c(
    1 - model_cumulative[1],
    -diff(model_cumulative),
    model_cumulative[length(model_cumulative)]
  )
  simulator_cumulative <- stats::plogis(rev(parameters$orm_intercepts) + eta)
  simulator_probabilities <- rev(c(
    simulator_cumulative[1],
    diff(simulator_cumulative),
    1 - simulator_cumulative[length(simulator_cumulative)]
  ))

  expect_equal(simulator_probabilities, model_probabilities, tolerance = 1e-15)
  expect_gte(min(simulator_probabilities), 0)
  expect_equal(sum(simulator_probabilities), 1)
})

test_that("ACTT-1 partial PO time deviation calibrates mortality", {
  parameters <- mostr:::actt1_markov_parameters()
  occupancy <- c(0, 0, 0, parameters$baseline_probabilities, 0)

  for (day in 1:28) {
    transition <- matrix(0, nrow = 8, ncol = 8)
    for (previous_state in 1:7) {
      eta <- mostr:::actt1_markov_lp(
        yprev = factor(previous_state, levels = 1:8),
        t = day,
        tx = 0,
        parameter = 0,
        extra_params = c(
          parameters$extra_params,
          treatment_effect_decay = 0
        )
      )
      cumulative <- stats::plogis(rev(parameters$orm_intercepts) + eta)
      transition[previous_state, ] <- rev(c(
        cumulative[1],
        diff(cumulative),
        1 - cumulative[length(cumulative)]
      ))
    }
    transition[8, 8] <- 1
    occupancy <- drop(occupancy %*% transition)
  }

  expect_equal(occupancy[8], 0.0953, tolerance = 0.001)
})

test_that("sim_actt1_markov validates wrapper inputs", {
  expect_snapshot(error = TRUE, sim_actt1_markov(n_patients = 0))
  expect_snapshot(error = TRUE, sim_actt1_markov(treatment_prob = 2))
  expect_snapshot(error = TRUE, sim_actt1_markov(treatment_effect = Inf))
  expect_snapshot(error = TRUE, sim_actt1_markov(treatment_effect_decay = -1))
})

test_that("sim_actt2_markov returns reproducible ACTT-2 trajectories", {
  traj <- sim_actt2_markov(
    n_patients = 80,
    follow_up_time = 6,
    seed = 1987
  )
  repeated <- sim_actt2_markov(
    n_patients = 80,
    follow_up_time = 6,
    treatment_effect = -0.0657,
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

test_that("sim_actt2_markov samples the supplied baseline distribution", {
  traj <- sim_actt2_markov(
    n_patients = 10000,
    follow_up_time = 1,
    seed = 4815
  )
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

test_that("sim_actt2_markov uses one for the baricitinib arm", {
  placebo <- sim_actt2_markov(
    n_patients = 20,
    treatment_prob = 0,
    follow_up_time = 1,
    seed = 91
  )
  baricitinib <- sim_actt2_markov(
    n_patients = 20,
    treatment_prob = 1,
    follow_up_time = 1,
    seed = 91
  )

  expect_equal(unique(placebo$tx), 0)
  expect_equal(unique(baricitinib$tx), 1)
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
  extra_params <- c(parameters$extra_params, treatment_effect_decay = 0)
  lp <- function(yprev, day, tx, decay = 0) {
    extra_params[["treatment_effect_decay"]] <- decay
    mostr:::actt2_markov_lp(
      yprev = factor(yprev, levels = 1:8),
      t = day,
      tx = tx,
      parameter = -0.0657,
      extra_params = extra_params
    )
  }

  expect_equal(lp(4, 7, 1) - lp(4, 7, 0), -0.0657)
  expect_equal(
    lp(4, 7, 1, decay = 0.2) - lp(4, 7, 0, decay = 0.2),
    -0.0657 * exp(-0.2 * 6)
  )
  expect_equal(lp(5, 7, 1) - lp(4, 7, 1), 7.4148 - 4.5199)

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
    c(-1.7503, -2.4174, -2.6911, -4.8406, -8.9716, -12.6833, -19.5542)
  )

  eta <- mostr:::actt2_markov_lp(
    yprev = factor(5, levels = 1:8),
    t = 10,
    tx = 0,
    parameter = -0.0657,
    extra_params = c(parameters$extra_params, treatment_effect_decay = 0)
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

test_that("sim_actt2_markov validates wrapper inputs", {
  expect_snapshot(error = TRUE, sim_actt2_markov(n_patients = 0))
  expect_snapshot(error = TRUE, sim_actt2_markov(treatment_prob = 2))
  expect_snapshot(error = TRUE, sim_actt2_markov(treatment_effect = Inf))
  expect_snapshot(error = TRUE, sim_actt2_markov(treatment_effect_decay = -1))
})

test_that("sim_trajectories_markov validates inputs and absent covariates", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1))

  expect_error(
    sim_trajectories_markov(baseline_data = list()),
    "baseline_data must be a data frame",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_markov(baseline_data = data.frame(id = 1)),
    "baseline_data must contain columns: yprev",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_markov(baseline_data = baseline, intercepts = character()),
    "intercepts must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_markov(baseline_data = baseline, lp_function = "nope"),
    "lp_function must be a function",
    fixed = TRUE
  )

  expect_warning(
    result <- sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 2,
      intercepts = 0,
      lp_function = function(yprev, t, parameter, extra_params) 0,
      covariate_names = "missing",
      seed = 1
    ),
    "None of the specified covariate_names found"
  )

  expect_equal(nrow(result), 4)
  expect_equal(names(result), c("id", "time", "y", "yprev"))
})

test_that("sim_trajectories_markov carries forward once all patients are absorbing", {
  baseline <- data.frame(id = 1:3, yprev = c(2, 2, 2), tx = 0)

  result <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 4,
    intercepts = -Inf,
    lp_function = function(yprev, t, tx, parameter, extra_params) 0,
    absorbing_states = 2,
    seed = 1,
    covariate_names = "tx"
  )

  expect_equal(unique(result$y), 2)
  expect_equal(nrow(result), 12)
})

test_that("sim_trajectories_markov passes tx even when not requested as a covariate", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1), tx = c(0, 1))

  expect_warning(
    result <- sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 2,
      intercepts = 0,
      lp_function = function(yprev, t, tx, parameter, extra_params) {
        ifelse(tx == 1, 10, -10)
      },
      covariate_names = character(0),
      seed = 1
    ),
    "None of the specified covariate_names found",
    fixed = TRUE
  )

  expect_equal(nrow(result), 4)
  expect_true(all(result$tx %in% c(0, 1)))
})

test_that("sim_trajectories_markov preserves character patient IDs", {
  baseline <- data.frame(id = c("b", "a"), yprev = c(1, 1), tx = c(0, 1))

  result <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 2,
    intercepts = 0,
    lp_function = function(yprev, t, tx, parameter, extra_params) 0,
    covariate_names = "tx",
    seed = 1
  )

  expect_identical(result$id, rep(c("b", "a"), each = 2))
  expect_equal(result$time, rep(1:2, 2))
  expect_equal(result$tx, rep(c(0, 1), each = 2))
})

test_that("sim_trajectories_markov supports threshold-specific linear predictors", {
  baseline <- data.frame(id = 1:4, yprev = c(1, 1, 1, 1), tx = c(0, 1, 0, 1))

  scalar <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 3,
    intercepts = c(-1, 1),
    lp_function = function(yprev, t, tx, parameter, extra_params) {
      tx * parameter
    },
    parameter = 0.5,
    absorbing_states = integer(),
    covariate_names = "tx",
    seed = 10
  )

  scalar_as_vector <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 3,
    intercepts = c(-1, 1),
    lp_function = function(yprev, t, tx, parameter, extra_params) {
      rep(tx * parameter, 2)
    },
    parameter = 0.5,
    absorbing_states = integer(),
    covariate_names = "tx",
    seed = 10
  )

  expect_identical(scalar, scalar_as_vector)

  ppo <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 3,
    intercepts = c(-1, 1),
    lp_function = function(yprev, t, tx, parameter, extra_params) {
      tx * parameter * c(0, 1)
    },
    parameter = 0.5,
    absorbing_states = integer(),
    covariate_names = "tx",
    seed = 10
  )

  expect_named(ppo, c("id", "time", "y", "yprev", "tx"))
  expect_equal(nrow(ppo), 12)
  expect_true(all(ppo$y %in% 1:3))
})

test_that("sim_trajectories_markov validates threshold-specific linear predictors", {
  baseline <- data.frame(id = 1, yprev = 1, tx = 1)

  expect_error(
    sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 1,
      intercepts = c(-1, 0, 1),
      lp_function = function(yprev, t, tx, parameter, extra_params) c(0, 0),
      absorbing_states = integer(),
      covariate_names = "tx",
      seed = 1
    ),
    "lp_function must return either length 1 or length\\(intercepts\\)",
    fixed = FALSE
  )

  expect_error(
    sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 1,
      intercepts = c(-1, 0, 1),
      lp_function = function(yprev, t, tx, parameter, extra_params) NA_real_,
      absorbing_states = integer(),
      covariate_names = "tx",
      seed = 1
    ),
    "lp_function must return finite numeric values",
    fixed = TRUE
  )

  expect_error(
    sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 1,
      intercepts = c(0, 1),
      lp_function = function(yprev, t, tx, parameter, extra_params) c(5, -5),
      absorbing_states = integer(),
      covariate_names = "tx",
      seed = 1
    ),
    "create decreasing cumulative probabilities",
    fixed = TRUE
  )
})


test_that("ACTT-2 Markov trajectories support 60-day follow-up", {
  data <- sim_actt2_markov(n_patients = 20, follow_up_time = 60, seed = 123)
  expect_equal(nrow(data), 1200)
  expect_equal(sort(unique(data$time)), 1:60)
  expect_equal(data$y[data$yprev == 8], rep(8L, sum(data$yprev == 8)))
  expect_identical(
    data,
    sim_actt2_markov(n_patients = 20, follow_up_time = 60, seed = 123)
  )
})
