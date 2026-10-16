test_that("sim_violet_markov returns reproducible VIOLET trajectories", {
  traj <- sim_violet_markov(n_patients = 30, seed = 1904)

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "time", "y", "yprev", "age", "sofa", "tx"),
    n_patients = 30,
    follow_up_time = 28,
    states = 1:4
  )
  expect_identical(traj, sim_violet_markov(n_patients = 30, seed = 1904))
  expect_equal(sum(violet_baseline$n), 1351)
})
