# mostr (development version)

* New experimental `orm_markov_ml()` fits the first-order proportional-odds
  Markov model of `orm_markov()` by maximum likelihood when some states are
  missing (`NA`) or known only up to a range (`rms::Ocens()`, for example
  "alive, level unknown"). Instead of dropping every transition that touches
  such a visit, it sums over the states the patient could have been in.
  `missing_states` says whether a missing state may hide a death. With nothing
  missing it gives the same fit as `orm_markov()`. The fit works with
  `sops()`, `avg_sops()`, `avg_time()`, `avg_comparisons()`, and
  `inferences()` (delta, MVN, score bootstrap, and refit bootstrap; not
  fractional weighted bootstrap). Standard errors use the observed
  information of the observed-data likelihood. The estimates are valid when
  states are missing at random.

* New experimental `blrm_markov_ml()` is the Bayesian counterpart of
  `orm_markov_ml()`: it samples the posterior of the `blrm_markov()` model,
  with the priors of `rmsb::blrm()` (Dirichlet prior for the intercepts,
  `pcontrast` and `npcontrast`), using the same observed-data likelihood for
  missing and range-coded states. Constrained partial proportional odds
  (`ppo` with `cppo`) and a random intercept per patient (`cluster(id)`,
  with the priors of `psigma`, `rsdmean`, and `rsdsd`) are supported; the
  SOP functions use the patients' random intercepts with
  `include_re = TRUE`. With nothing missing it samples the same posterior as
  `blrm_markov()`. It runs its own
  Stan program, which needs rstan and a C++ toolchain; the program is compiled
  on first use and cached in `tools::R_user_dir("mostr", "cache")`. The fit
  works with `sops()`, `avg_sops()`, `avg_time()`, and `avg_comparisons()`
  like a `blrm_markov()` fit.

* `orm_markov_ml()` now rejects random effects (`cluster()`) in the formula
  with an informative error.

* `orm_markov_ml()` and `blrm_markov_ml()` no longer require
  `absorbing_state`. Without it (the new default `NULL`), no state is
  absorbing: every state can be left, as on a quality-of-life scale without
  death, and a missing state may be any state under either `missing_states`
  rule. Data with a state that is never left, such as death, should still
  give it in `absorbing_state`; if that state is never a previous state in
  the data, the fit stops with a message saying so.

* `orm_markov_ml()` without missing or range-coded states is now about as fast
  as `orm_markov()`, also with many states: it uses the `rms::orm()` fit and
  its exact information directly instead of running EM, Newton steps, and
  numerical derivatives over every row (with 500 states this took about an
  hour). The missing-information fraction is then exactly 0. With missing
  states, transitions with both states observed now add their exact
  information, so numerical derivatives cover only the rows next to a missing
  or range-coded state. A new "Computing time" section in `?orm_markov_ml`
  explains how the time grows with missing states and the number of states.

* Printing an `orm_markov_ml()` fit now shows the maximised log-likelihood,
  the numbers of patients and rows, and a table of coefficients with
  patient-cluster robust standard errors, Wald Z statistics, and p-values. It
  no longer shows the `rms::orm()` summary (observation count, response
  frequencies, likelihood-ratio test, R2, rank correlations), which described
  the internal weighted fit on expected transition counts rather than the
  data. Printing a fit made with `orm_markov_ml_control(variance = FALSE)`
  no longer fails.

* `logLik()`, `AIC()`, and `BIC()` on an `orm_markov_ml()` fit now use the
  maximised observed-data log-likelihood instead of the internal weighted fit
  on expected transition counts. `nobs()` returns the number of rows whose
  observation rules out at least one state. With nothing missing, all four
  give the same values as for `orm_markov()`.

* The `stats`, `freq`, and `deviance` components of an `orm_markov_ml()` fit
  now describe the data rather than the internal weighted fit: `stats["Obs"]`
  equals `nobs()`, `freq` counts the exactly observed states, and the model
  deviance is -2 times the maximised log-likelihood. Statistics that would
  need an intercept-only fit to the incomplete data (likelihood-ratio and
  score tests, R2) or an observed state on every row (rank correlations such
  as Dxy and rho) are `NA`.

* New vignette `vignette("missing-qol-recorded-deaths")` shows why dropping
  transitions next to missing quality-of-life questionnaires overestimates
  mortality when deaths are always recorded, and how `orm_markov_ml()` with
  `missing_states = "non_absorbing"` avoids this. It also shows how to code
  partially known states with `rms::Ocens()` and how to analyze patients lost
  to follow-up.

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
