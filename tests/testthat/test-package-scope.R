# Refit bootstraps always resample the data stored on a mostm() fit, so the
# public API accepts no separate refit data and the low-level bootstrap helpers
# stay internal.
test_that("the package exposes Markov workflows without removed APIs", {
  # `load_all()` exports every function, so read the declared exports.
  namespace_lines <- readLines(system.file("NAMESPACE", package = "mostr"))
  declared_exports <- sub(
    "^export\\((.*)\\)$",
    "\\1",
    grep("^export\\(", namespace_lines, value = TRUE)
  )
  argument_names <- function(funs) {
    unique(unlist(lapply(funs, function(fun) names(formals(fun)))))
  }

  cases <- list(
    list(
      label = "removed functions",
      names = c(
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
        "orm_markov",
        "vglm_markov",
        "blrm_markov",
        "tidy_po",
        "calc_time_in_state_diff",
        "lp_violet",
        "sim_trajectories_markov",
        "sim_actt1_markov",
        "sim_actt2_markov",
        "sim_violet_markov"
      ),
      within = ls(asNamespace("mostr"), all.names = TRUE),
      present = FALSE
    ),
    list(
      label = "main workflow exports",
      names = c(
        "sim_markov",
        "mostm",
        "sops",
        "avg_sops",
        "avg_state_time",
        "avg_comparisons",
        "inferences",
        "time_in_state",
        "bootstrap_model_coefs"
      ),
      within = getNamespaceExports("mostr"),
      present = TRUE
    ),
    list(
      label = "separate refit data in SOP functions",
      names = "refit_data",
      within = argument_names(list(sops, avg_sops, avg_comparisons, avg_state_time)),
      present = FALSE
    ),
    list(
      # The fit records the absorbing states and the gap column, and
      # `time_covariates` was replaced by time terms in the formula.
      label = "fit settings and time_covariates in prediction and plot functions",
      names = c("absorb", "gap_var", "time_covariates"),
      within = argument_names(list(
        sops, avg_sops, state_time, avg_state_time, comparisons,
        avg_comparisons, soprob_markov, plot_lp_difference, plot_transitions,
        plot_correlation, plot_variogram
      )),
      present = FALSE
    ),
    list(
      label = "fit settings recorded by mostm()",
      names = c("absorb", "gap_var"),
      within = argument_names(list(mostm)),
      present = TRUE
    ),
    list(
      label = "data argument of bootstrap_model_coefs()",
      names = "data",
      within = argument_names(list(bootstrap_model_coefs)),
      present = FALSE
    ),
    list(
      label = "low-level bootstrap helpers in the declared exports",
      names = c(
        "bootstrap_analysis_wrapper",
        "apply_to_bootstrap",
        "fast_group_bootstrap",
        "materialize_bootstrap_sample_indexed",
        "drop_absent_states"
      ),
      within = declared_exports,
      present = FALSE
    ),
    list(
      label = "inference functions in the declared exports",
      names = c("inferences", "bootstrap_model_coefs"),
      within = declared_exports,
      present = TRUE
    )
  )

  for (case in cases) {
    unexpected <- if (case$present) {
      setdiff(case$names, case$within)
    } else {
      intersect(case$names, case$within)
    }
    expect_equal(unexpected, character(), info = case$label)
  }
})
