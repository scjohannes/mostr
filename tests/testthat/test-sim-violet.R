test_that("VIOLET preset returns reproducible trajectories", {
  traj <- sim_markov("violet", n_patients = 30, seed = 1904)

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "time", "y", "yprev", "age", "sofa", "tx"),
    n_patients = 30,
    follow_up_time = 28,
    states = 1:4
  )
  expect_identical(traj, sim_markov("violet", n_patients = 30, seed = 1904))
  expect_identical(
    traj,
    sim_markov(
      "violet",
      n_patients = 30,
      treatment_effect = log(0.75),
      seed = 1904
    )
  )
  expect_equal(sum(violet_baseline$n), 1351)
})
