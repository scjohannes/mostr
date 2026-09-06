test_that("the package exposes Markov workflows without removed APIs", {
  excluded <- c(
    "assess_operating_characteristics",
    "format_competing_risks",
    "jackknife_mcse",
    "plot_operchar",
    "recurr_event",
    "run_power_iteration",
    "sample_from_arrow",
    "sim_actt2_brownian",
    "sim_actt2_markov_60day",
    "sim_trajectories_brownian",
    "sim_trajectories_brownian_gap",
    "sim_trajectories_deterministic",
    "sim_trajectories_tte",
    "states_to_drs",
    "states_to_hce",
    "states_to_tte",
    "states_to_tte_v2",
    "states_to_ttest",
    "summarize_oc_results",
    "summarize_power_results",
    "tidy_bootstrap_coefs",
    "tidy_po",
    "calc_time_in_state_diff"
  )
  expect_equal(
    intersect(excluded, ls(asNamespace("mostr"), all.names = TRUE)),
    character()
  )
  retained <- c(
    "sim_trajectories_markov",
    "sim_actt1_markov",
    "sim_actt2_markov",
    "orm_markov",
    "vglm_markov",
    "blrm_markov",
    "sops",
    "avg_sops",
    "avg_comparisons",
    "inferences",
    "time_in_state",
    "apply_to_bootstrap"
  )
  expect_equal(setdiff(retained, getNamespaceExports("mostr")), character())
})
