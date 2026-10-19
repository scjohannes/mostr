# mostr (development version)

## Breaking changes

* `sim_markov()` replaces `sim_trajectories_markov()`, `sim_actt1_markov()`,
  and `sim_actt2_markov()`. Choose a model with `preset`: `"actt1"`,
  `"actt2"`, `"violet"`, or `NULL` to supply your own `baseline_data`,
  `intercepts`, and `lp_function`. `preset` has no default. Any model argument
  you supply, such as `intercepts` or `baseline_data`, replaces that part of
  the preset and keeps the rest. Each preset's model is described in its own
  help topic: `?actt1_preset`, `?actt2_preset`, and `?violet_preset`.
  `sim_markov("actt2", ...)` gives the same trajectories as
  `sim_actt2_markov(...)` did with the same arguments and seed.

* `sim_markov()` simulates exactly the patients in a supplied `baseline_data`
  and stops if you also supply `n_patients` or `treatment_prob`, which only
  sample patients from a preset. Previously the trial wrappers recycled the
  columns of a supplied `baseline_data` over `n_patients` generated patients,
  so a 2-patient `baseline_data` with the default `n_patients = 1000` gave
  1,000 patients.

* `preset = "actt1"` uses the ACTT-2 transition model with the starting states
  reported for ACTT-1 (13.1%, 41.4%, 18.4%, and 27.1% in states 4 to 7) and a
  threshold for death raised from -22.0693 to -21.4071, so that 18% of control
  patients have died by day 28. The ACTT-1 model based on Rohde et al. (2024)
  that `sim_actt1_markov()` used has been removed, and
  `vignette("partial-po-vglm-sops")` now simulates VIOLET data.

* An `lp_function` now receives `treatment_effect`, the treatment effect on
  the current day after `treatment_effect_decay`, instead of `parameter`, so
  `treatment_effect_decay` also works for your own models. With
  `preset = NULL`, `treatment_effect` defaults to 0, `extra_params` to `NULL`,
  `absorbing_states` to the worst state (`length(intercepts) + 1`),
  `covariate_names` to every column of `baseline_data` except `id` and
  `yprev`, and `follow_up_time` to 28 days (previously 60). `lp_violet()` has
  been removed.

## New features

* New `preset = "violet"` simulates 4-state daily trajectories (home, in
  hospital, on a ventilator or with ARDS, dead) from the Markov model that
  Frank Harrell fitted to the VIOLET trial of vitamin D3 to create his
  simulated trial `simlongord`. Patients are sampled from the 1,351 VIOLET
  patients' starting states, ages, and SOFA scores. The default treatment odds
  ratio is 0.75, as in `simlongord`, and the daily state proportions agree
  with `simlongord` within simulation error.

* `violet_baseline` now holds the 1,351 VIOLET patients as 858 distinct
  combinations of starting state, age, and SOFA score with their counts `n`
  (previously 250,000 resampled patients in a modified 6-state coding). Its
  states follow `sim_markov("violet")`.

## Minor improvements

* `sim_markov("actt1")` and `sim_markov("actt2")` are much faster for large
  samples, because they now calculate each day's linear predictor for all
  patients in one step instead of one patient at a time, as
  `sim_markov("violet")` already did. The trajectories for a given seed are
  unchanged. A model that you supply through `lp_function` is still called
  once per patient.

# mostr 0.1.1

* `sim_actt2_markov()` uses a model refitted with states 1 and 2 merged,
  because state 2 was mostly recorded only at study visits. The old fit
  rehospitalized about 29% of control patients (2% observed); the new one
  about 2.7%. The default `treatment_effect` of `sim_actt2_markov()` and
  `sim_actt1_markov()` is now -0.0730.

# mostr 0.1.0

* `sim_actt1_markov()` now adds the average age and sex contribution of an
  ACTT-1 patient (0.2253 on the log-odds scale) to every threshold. The
  reported thresholds describe a woman aged 0, so simulated patients
  previously moved to better states too quickly; placebo mortality by day 28
  rises from 9.5% to 13.4%, closer to the reported 15.2% by day 29.

* Bootstrap refits (`inferences(method = "bootstrap")` and `method = "fwb"`)
  now always resample the data stored by `orm_markov()` or `vglm_markov()`:
  the rows used for fitting, after any `subset`, with their
  fitting weights. The `refit_data` argument of `sops()`, `avg_sops()`,
  `avg_comparisons()`, `avg_time()`, `plot_transitions()`,
  `plot_correlation()`, and `plot_variogram()` has been removed; drop it from
  existing calls. `newdata` still sets the patient profiles used for
  prediction and can contain any profiles. Refit bootstraps of models fitted
  without these wrappers now error.

* `bootstrap_model_coefs()` no longer has a `data` argument. It resamples the
  patients stored on a model fitted with `orm_markov()` or `vglm_markov()`, and
  `id_var` now defaults to the patient ID column given when fitting.

* `fast_group_bootstrap()`, `materialize_bootstrap_sample()`,
  `apply_to_bootstrap()`, `bootstrap_analysis_wrapper()`, and
  `relevel_factors_consecutive()` are no longer exported. Use `inferences()` or
  `bootstrap_model_coefs()` for bootstrap calculations.

* Bootstrap refits use starting values (from `use_coefstart = TRUE` or a
  `coefstart` passed to `vglm_markov()`) only when the resampled data contain
  every outcome state. When a state is missing, the refitted model has fewer
  intercepts, so it is fitted without starting values. A refit with starting
  values that fails, for example because a covariate level is absent from the
  resample, is retried once without them.

* `avg_sops()` now accepts `variables = NULL` to average over patients using
  their observed starting states and covariates, with the existing analytical,
  posterior, simulation, and bootstrap inference methods.

* `avg_time()` estimates average total time in separate or pooled states, using
  visit sums or mapped elapsed-time integration. It supports conditional and
  unconditional analytical variance, posterior intervals, and all existing
  simulation and bootstrap methods. Use `vcov()` for analytical covariance and
  `get_draws()` for retained time-total draws.

* With `time_map` and no `target_times`, `avg_comparisons(estimand =
  "time_in_state")` (and `"time_benefit"`) and `time_in_state()` for SOP data
  frames now integrate from `baseline_time`
  (default `0`) to the last mapped visit, including the interval before the
  first visit. Previously, the omitted grid started at the first mapped visit.
  Set `baseline_time = NULL` or supply `target_times` to choose the period.
  `avg_time()` uses the same default.

* When `variables` sets the starting state (for example
  `variables = list(yprev = c("1", "2"))`), real-time summaries from
  `avg_time()`, `avg_comparisons(estimand = "time_in_state")`, and
  `interpolate_sops()` now start each scenario in the set state at
  `baseline_time`. Previously they started from the patients' observed starting
  states. These results support only conditional variance:
  `inferences(method = "delta")` requires `vcov = "conditional"` (or a
  coefficient covariance matrix), and `method = "mvn"` remains available;
  unconditional delta, score-bootstrap, bootstrap, and FWB inference now error.

* `plot_time()` plots `avg_time()` results. `type = "bar"` draws one
  horizontal stacked bar per scenario, split into the average time in each
  state ("Grotta bars"); `scale = "proportion"` shows shares of follow-up.
  `type = "pointrange"` draws one point per state and scenario with its
  vertical confidence interval. Scenarios are drawn in black and successively
  lighter shades of gray.

* `inferences()` now tests `avg_comparisons()` results against "no effect" by
  default: 0 for differences and 1 for ratios. The result gains `statistic`,
  `p.value`, and `s.value` columns. The new default is `null = "auto"`; pass
  a number to test against another value, or `null = NULL` to skip the test.
  SOP and average-time results are still not tested unless you supply `null`.

* `plot_correlation()` and `plot_variogram()` now default to Spearman rank
  correlation for observed data and fitted models, accounting for tied states.
  Use `method = "pearson"` to retain the previous calculation.

* `plot_correlation()` now uses a 0-1 color scale when all available correlations
  are nonnegative, retaining -1 to 1 when negative correlations are present.
  Set `fill_limits` explicitly to keep a common scale across plots.

* Speed up row assembly, grouped summaries, patient bootstrap materialization,
  and baseline selection using data.table internally. Results remain data frames.

* Initial release of Markov Ordinal State Transition modeling tools, migrated
  from markov.misc with its original authorship and GPL license preserved.
* Retain Markov data generation, first- and second-order transition models,
  state occupancy probabilities, treatment comparisons, time in state,
  diagnostics, and analytical, posterior, MVN, and bootstrap inference.
* Remove non-Markov generators, simulation-study and power workflows, Arrow
  sampling, endpoint conversions, competing-risk analyses, and coefficient
  tidiers. The simulation-truth helper `calc_time_in_state_diff()` is also
  excluded; model-based summaries remain available through `time_in_state()`.
* Use Markov-generated data throughout examples, tests, and all nine vignettes.
* `sim_actt2_markov()` supports custom follow-up, including 60 days; the separate
  Brownian-calibrated 60-day generator is excluded.
