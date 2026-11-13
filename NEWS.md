# mostr (development version)

## Breaking changes

* A second-order model is declared once, when it is fitted:
  `mostm(..., p2_var = "ypprev")` names the column with the state two visits
  back. `mostm()` also takes `p_var`, the previous-state column (`"yprev"` by
  default). The fit records these with `id_var` and `time_var`, and its
  states in fitted order, and every model-based function uses them. Before,
  leaving out `p2_var` for a second-order fit, or `time_var` with `newdata`
  for a fit whose time column is not `"time"`, silently held that column at
  each profile's starting value and gave wrong probabilities.

  `sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
  `avg_comparisons()`, `soprob_markov()`, and `plot_lp_difference()` no
  longer have the arguments `id_var`, `time_var`, `p_var`, `p2_var`, and
  `y_levels`; passing one is an error. Drop them from these calls, and give
  `p2_var`, `p_var`, and `time_var` to `mostm()` instead when they differ from
  the defaults. The `state` column of the results now always holds the fitted
  state labels (for example `"1"`, `"2"`, ...); before, it took the type of a
  supplied `y_levels`, such as integers for `y_levels = 1:6`.

  `plot_transitions()`, `plot_correlation()`, and `plot_variogram()` keep these
  arguments, because they also plot observed data, where nothing else names
  the columns. For a fitted model they default to the fit's values, and a
  value that differs from the fit's is an error, now also a `y_levels` with a
  different number of states.

* `mostm()` replaces `orm_markov()`, `vglm_markov()`, and `blrm_markov()`
  (#11). Choose the package that fits the model with `engine`: `"orm"` (the
  default) uses `rms::orm()`, `"vglm"` uses `VGAM::vglm()`, and `"blrm"` uses
  `rmsb::blrm()`. The shared arguments `id_var`, `time_var`,
  `first_followup_time`, `weights`, `subset`, `type`, and `cadjust` work as
  before, except that `id_var` is now required (see below). Translate existing
  calls as follows:

  * `orm_markov(formula, data, ...)` becomes `mostm(formula, data, ...)`.
  * `vglm_markov(formula, family, data, ...)` becomes
    `mostm(formula, data, ..., engine = "vglm", family = family)`.
  * `blrm_markov(formula, ppo, cppo, data, ...)` becomes
    `mostm(formula, data, ..., engine = "blrm", ppo = ppo, cppo = cppo)`.

  `sops()` and the other model-based functions do not accept fits made by the
  removed functions; refit them with `mostm()`.

* Arguments for the engine go through `...` and must be named. An argument
  that belongs to a different engine is an error that names that engine, and
  so is an argument that the `"orm"` or `"vglm"` engine does not have. With
  `engine = "blrm"`, other arguments are passed on to `rmsb::blrm()`, including
  its own `backend`, which chooses between rstan and cmdstan.

* `x`, `y`, `x.arg`, and `y.arg` are no longer accepted: `mostm()` always keeps
  the design matrix and response, which the state occupancy and inference
  functions need.

* `engine = "vglm"` uses `family = VGAM::cumulative(reverse = TRUE,
  parallel = TRUE)` when `family` is not given. This is the proportional-odds
  model that `engine = "orm"` fits; `vglm_markov()` required `family`.

* `data` must be a data frame for every engine. `vglm_markov()` also looked up
  variables in the formula's environment when `data` was missing.

* `mostm()` requires `id_var`, the name of the column that identifies
  patients, for every engine; a fit without it is an error. The removed
  functions only warned, or, for `vglm_markov()`, warned only when `data` had
  an `id` column. Such fits did not store each patient's starting row and used
  the model-based covariance instead of patient-clustered standard errors.
  Later steps then failed or quietly used a column called `id`: for example,
  `avg_sops(newdata = ...)` failed with "ID variable 'id' not found" when the
  patient column had another name, and the bootstrap resampled patients by an
  `id` column the fit did not know about. Because `mostm()` now always finds
  each patient's starting row, it also requires the time column named by
  `time_var`. Every `"orm"` and `"vglm"` fit now computes the
  patient-clustered covariance, so a fit for which it cannot be computed, such
  as one whose coefficients grow without bound, or an `"orm"` fit with a link
  other than logit, stops in `mostm()` instead of later.

* `bootstrap_model_coefs()` no longer has an `id_var` argument. It resamples
  patients by the `id_var` given to `mostm()`.

* `get_effective_coefs()` accepts the `robcov_vglm` object that
  `mostm(engine = "vglm")` returns.

* `weights`, `type`, and `cadjust` are errors with `engine = "blrm"`:
  `rmsb::blrm()` has no row weights, and Bayesian fits take their uncertainty
  from the posterior draws instead of patient-clustered standard errors.

* Bootstrap refits (`inferences(method = "bootstrap")`, `method = "fwb"`, and
  `bootstrap_model_coefs()`) call `mostm()` again with the original engine and
  arguments. They no longer compute patient-clustered standard errors for each
  refit, which only the original fit uses.

* `avg_sops()`, `avg_comparisons()`, and `avg_state_time()`, and `sops()`,
  `comparisons()`, and `state_time()` when they average with `by`, now stop
  when a `by` column has a missing value among the averaged profiles, and say
  how many profiles are affected. Previously those profiles were silently
  left out of every group, which changed who was averaged. Drop or impute
  those patients first, in `newdata` or in the data given to `mostm()`. Only
  the averaged profiles count: with stored data, a missing value at a later
  visit, which is not a starting row, does not matter. Summaries such as
  `time_in_state()` still keep missing values as their own group.

* A patient can no longer start in an absorbing state. With `absorb`,
  `sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
  `avg_comparisons()`, and `soprob_markov()` stop when a profile's previous
  state (`p_var`, usually `yprev`) is absorbing, or, with `p2_var`, when its
  state two visits back is. The error names the column and says how many
  profiles are affected; remove those patients from `newdata` first.
  Previously, first-order models predicted a move out of the absorbing state
  from a previous state never seen in the fitted data (a numeric `yprev` gave
  a probability of staying dead below 1, a factor `yprev` failed with an
  unrelated error), while second-order models kept such patients in the
  absorbing state with probability 1.

* `inferences()` no longer offers `method = "score_bootstrap"`; asking for it
  is an error. Use `method = "mvn"`, `"delta"`, `"bootstrap"`, or `"fwb"`
  instead. The score bootstrap has been removed from the package entirely.

* `inferences(method = "delta")` refuses unconditional variance for every
  model fitted with case weights (`mostm(weights = )`), whatever the engine,
  because the calculation has not been checked with weights. Previously
  weighted `"orm"` fits were refused but weighted `"vglm"` fits were not.
  Unconditional variance is the default for `avg_sops()`,
  `avg_state_time()`, and `avg_comparisons()` results, so for a weighted fit
  pass `vcov = "conditional"`.

* `include_re = TRUE` is an error for fits other than `rmsb::blrm()` in
  `sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
  `avg_comparisons()`, and `soprob_markov()`. Only `blrm` fits have random
  effects; previously the argument was silently ignored for `"orm"` and
  `"vglm"` fits.

* `use_coefstart = TRUE` in `inferences()` and `bootstrap_model_coefs()` is an
  error unless the model was fitted with `engine = "vglm"`. It is meant to
  start each bootstrap refit from the original coefficients, but only `vglm`
  refits take starting coefficients, so for `"orm"` fits the argument was
  silently ignored. `inferences()` also refuses it with `method = "mvn"` or
  `"delta"`, which do not refit the model. Leave `use_coefstart` at its
  default, `FALSE`; the results are the same as before.

* State occupancy probabilities from partial proportional-odds models now stop
  with an error where the model's thresholds (cut-points) cross. For some
  covariate values, often outside the range of the data, a predictor with a
  separate effect at each threshold can push two thresholds past each other,
  and the model then gives a negative probability for a state. Previously the
  compiled engine set such a probability to zero and rescaled the others, the
  reference engine for `"vglm"` fits kept the negative probability, and
  `blrm` fits clipped and rescaled each posterior draw, so the two
  `"vglm"` engines could disagree by as much as 0.5. Now `sops()`,
  `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
  `avg_comparisons()`, and `soprob_markov()` stop and say for how many
  profiles and at which time the thresholds cross (and in how many posterior
  draws for `blrm` fits). `inferences()` with `method = "mvn"`,
  `"bootstrap"`, or `"fwb"` likewise stops, saying in how many simulation draws
  or bootstrap samples the thresholds cross, instead of dropping those draws
  with a warning. A negative value between -1e-10 and zero comes from rounding
  and is still set to zero.

* `inferences()` no longer has an `update_datadist` argument; passing it is an
  error. mostr does not need an rms `datadist`, and bootstrap refits
  (`inferences()` with `method = "bootstrap"` or `"fwb"`, and
  `bootstrap_model_coefs()`) no longer create a `dd` object in the global
  environment or set `options(datadist = "dd")` while they run. Instead, each
  refit runs with the `datadist` option unset and restores it afterwards, so a
  `datadist` you have set is ignored by the refits and left as it is.
  Previously, `update_datadist = FALSE` together with a `datadist` option
  naming an object that does not exist made every refit fail. Remove
  `update_datadist` from your calls; the estimates and standard errors are the
  same as before.

## Minor improvements

* When `rms::orm()` cannot fit a model, the error now suggests
  `opt_method = "LM", maxit = 100` together, and merging rarely observed
  states. Since rms 8.2-0, the `"LM"` algorithm often needs more than the
  default 30 iterations.

* Every error that mostr raises now has class `mostr_error`, so you can catch
  them with `tryCatch(..., mostr_error = function(e) ...)`. The "Error in" line
  names the mostr function you called, such as `sops()`, instead of an
  internal helper. The messages themselves are unchanged.

## Documentation

* The help pages of `avg_sops()`, `avg_state_time()`, `avg_comparisons()`,
  and `mostm()` now say that the averages give every profile the same weight:
  case weights given to `mostm(weights = )` are used only to fit the model,
  so the result is not the same as repeating each patient according to their
  weight.

## Bug fixes

* Parallel bootstrap and simulation inference (`workers` greater than 1) now
  restore every level of a nested future plan when they finish. Before, a
  plan such as `future::plan(list(multisession, multisession))` came back with
  only its outer level.

* `time_in_state()` on SOPs from `inferences(method = "delta")` now returns
  `std.error`, `conf.low`, and `conf.high` for the totals, with conditional or
  unconditional variance as chosen in `inferences()`. They match
  `avg_state_time()` and `state_time()` with delta inference. Before, the
  uncertainty was silently dropped. Delta SOPs that already went through
  `interpolate_sops()` now give an error; pass `time_map` to
  `time_in_state()` instead.

* `plot_lp_difference()` no longer needs `absorb` for a fit made by `mostm()`:
  by default it shows only the previous states that occur in the fitted data.
  Before, it failed on vglm fits ("replacement has 72 rows, data has 96") and
  drew empty panels for the absorbing state on orm fits. Asking for a state
  that never occurs as a previous state through `previous_states` is now an
  error that lists the valid states.

* `sops()` on the patients stored by `mostm()` no longer returns the internal
  column `.markov_source_row`, either in its estimates or in the results and
  draws of `inferences()`.

* `mostm()` no longer shows rms's warning "currently weights are ignored in
  model validation and bootstrapping orm fits" for a weighted `"orm"` fit. The
  warning concerns rms's own `validate()` and `bootcov()`, which mostr does not
  use; mostr's bootstrap and fractional weighted bootstrap keep the weights.

* The ordinary bootstrap (`inferences(method = "bootstrap")`) gave wrong
  results when a resample lacked a state other than the absorbing one. The
  refit renumbered the remaining states, for example 1, 3, 4 as 1, 2, 3, but
  the profiles supplied in `newdata` kept their labels. A profile in the
  missing state 2 was then predicted as if it started in state 3, without
  any message, and a profile in state 3 was read as the absorbing state, so
  the draw was dropped with a misleading warning. Averages within `by`
  groups reported the probability of state 3 under state 2. `orm` refits of
  such resamples always failed, and a numeric `yprev` was renumbered too,
  which changed linear or spline terms such as `rcs(yprev, 3)`. Refits now
  keep the original state labels and numbers and only drop the missing
  state. When some resample lacks a state that a supplied profile starts in,
  `inferences()` now stops before refitting, naming the state and the number
  of resamples that lack it; use `method = "fwb"`, which keeps every patient
  in every refit, or remove those profiles. A state that no profile starts
  in may still be missing from a resample, with probability 0 in that draw.

* Bootstrap and FWB draws whose refit or prediction fails are still left out
  with one warning each, and now also with one warning that says how many of
  the draws failed and why the first one did.

* `variables` that set the state two visits back in a second-order model,
  for example `variables = list(ypprev = c("1", "2"))`, now support only
  conditional variance, as `variables` that set `yprev` do: the bootstrap,
  FWB, and unconditional delta inference are refused. Previously they ran,
  and the bootstrap could predict from the wrong state.

* A comparison of a variable that the model doesn't use, for example
  `avg_comparisons(fit, variables = list(note = c(0, 1)))` where `note` is a
  column of the data but not of the model formula, is 0 (1 for a ratio)
  whatever the coefficients are. `inferences()` now reports its standard
  error, test statistic, p-value, s-value, and confidence interval as `NA`
  with every method, and so do the posterior intervals of `blrm` fits, as the
  marginaleffects package does. Previously the standard error was 0, the
  `statistic` `NaN`, and the interval of zero width. Such calls still run without a message, and
  predictions under such a variable from `sops()`, `avg_sops()`,
  `state_time()`, and `avg_state_time()` keep their standard errors. A column
  counts as used when it appears anywhere in the model formula, including
  inside terms such as `rms::rcs(age, 3)`, or is the previous state.

* `avg_comparisons()` and `comparisons()` with `by` naming the compared
  variable, for example `variables = list(tx = c(0, 1))` with `by = "tx"`,
  now group patients by the treatment they were observed with for every
  estimand, giving the average effect among the untreated and among the
  treated. Each group's estimate equals the same call with `newdata`
  restricted to that group, and intervals from `inferences()` with
  `method = "mvn"`, `"bootstrap"`, or `"fwb"` follow the
  same groups. Previously, SOP and time-in-state comparisons silently ignored
  such a `by` and returned the overall average without a `tx` column, and
  with `time_map` they mixed the starting states of one group with the
  averages of all patients. Time benefit was already grouped correctly.
  Analytical intervals (`method = "delta"`) remain unavailable with `by`.

* `avg_sops()` and `avg_state_time()` now stop with an explanation when `by`
  names a variable that `variables` also sets. Their result's column for that
  variable holds the value that was set, so `by` was silently ignored. To
  average within groups of the observed value, copy the column first, for
  example `newdata$observed_tx <- newdata$tx`, and use `by = "observed_tx"`.
  `by` without `variables`, or on other columns, works as before.

* With `engine = "vglm"`, starting values given through `etastart` or
  `mustart`, with one row per row of `data`, now drop the rows that `subset`
  or missing values remove. Previously such a fit stopped with a dimension
  mismatch.

* When `rms::orm()` cannot fit the model, for example because it does not
  converge, `mostm()` stops with an error that says so and suggests settings to
  try. Previously a fit with `id_var` stopped with an unrelated error about
  the robust covariance, and a fit without `id_var` returned the failed fit.

* For models fitted with `mostm(engine = "orm")` (previously `orm_markov()`),
  `newdata` no longer needs an outcome column. Prediction profiles without
  one, such as a single patient's starting state and covariates, previously
  failed with "object 'y' not found" in `sops()`, `avg_sops()`,
  `state_time()`, `comparisons()`, `avg_comparisons()`, and
  `plot_lp_difference()`, and you had to add a placeholder column such as
  `y = NA`. Results are unchanged (#6).

* With `engine = "orm"` or `"blrm"`, `mostm()` stops before fitting when the
  formula converts a predictor with `factor()`, `as.factor()`, `ordered()`, or
  `as.ordered()`, for example `ordered(y) ~ factor(time) + tx + yprev`. rms
  cannot fit such terms, and the fit previously stopped with "subscript out of
  bounds". The new error asks for a factor column in `data` instead. For time,
  store `time` as a factor column and set `first_followup_time`, as shown in
  `vignette("factor-time-orm-sops")`. `engine = "vglm"` still fits these
  formulas.

* `mostm()` no longer accepts random-effect terms that mostr cannot use. From
  rms 8.2 on, `rms::orm()` fits `cluster(id)` and `mix_re()` terms as random
  effects, and `mostm(engine = "orm")` accepted such a fit without a message.
  `sops()`, `avg_sops()`, and `inferences()` then described only a patient whose
  random effect is 0, which is not the average over patients, and their
  standard errors did not account for the random effects. With
  `engine = "orm"` or `"vglm"`, these terms are now errors, with or without a
  package prefix. `id_var` already gives standard errors that allow for
  correlation between the transitions of one patient. With `engine = "blrm"`,
  `cluster(id)` still adds a random intercept for each patient, but
  `rms::cluster(id)` and `mix_re()` are errors: `rmsb::blrm()` fitted them as
  ordinary predictors, so the first one fitted a straight-line effect of the
  patient ID number.

* The state occupancy and inference functions and `get_vcov_robust()` also
  reject an `rms::orm()` fit that holds random effects, such as a fit saved
  before this change. `get_vcov_robust()` now says so when `rms::orm()` could
  not fit the model, which rms 8.2 and later report by returning a fit without
  coefficients. Previously it stopped with an error about unsupported model
  families.

* Bootstrap refits (`inferences()` with `method = "bootstrap"` or `"fwb"`,
  and `bootstrap_model_coefs()`) now treat state labels as labels. When a
  resample lacked a state, the refit renumbered the remaining states and
  converted the labels to numbers, so text labels on factor columns, such as
  `"home"`, `"ward"`, `"dead"`, gave an "NAs introduced by coercion" warning
  for every refit and a missing state went unnoticed. States are now matched
  by label, in the model's state order, for numbers and text alike, and
  refits keep the labels (see the first bug fix above).

* Bootstrap refits now notice a missing state in the response column
  whatever it is called, for example `state` in
  `ordered(state) ~ yprev + time`. Previously only a column called `y` was
  checked, so with another name a resample that lacked a state was not
  noticed and that refit's state occupancy probabilities failed and were
  dropped. `bootstrap_model_coefs()` now treats only the response and `yprev`
  as states; previously it also renumbered every factor covariate, such as a
  factor `tx`, together with the states, which gave wrong labels or failed
  refits when a resample lacked a state.

# mostr 0.2.0

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

* `by` now works as in the marginaleffects package in every function that
  has it. Per-profile functions (`sops()`, `state_time()`, `comparisons()`)
  default to `by = FALSE`, and averaging functions (`avg_sops()`,
  `avg_state_time()`, `avg_comparisons()`) default to `by = TRUE`. In all of
  them, `TRUE` averages over all profiles, a character vector averages within
  groups, and `FALSE` or `NULL` returns one row per profile. Each function
  hands the other case to its partner. This changes three results:
  * `avg_sops()`, `avg_state_time()`, and `avg_comparisons()` with
    `by = NULL` now return per-profile rows instead of the overall average;
    drop `by` or use `by = TRUE` for the average.
  * `sops(by = "x")` now returns `avg_sops(by = "x")`, a `markov_avg_sops`
    object, instead of a grouped `markov_sops` object. The estimates are the
    same.
  * `sops(by = TRUE)` now returns the overall average instead of an error.

* Bare variable names, such as `variables = "tx"`, choose values as
  marginaleffects does: factor levels in level order, character values in
  sorted order, `FALSE` then `TRUE`, or `0` then `1` for a 0/1 numeric
  variable. Other numeric variables use their five-number summary in
  `sops()`, `avg_sops()`, `state_time()`, and `avg_state_time()`, and must
  be listed with values in `comparisons()` and `avg_comparisons()`, which now
  also accept bare names. Values come from the data stored with the model.
  Previously `avg_sops()` used every unique value of the prediction profiles
  in order of appearance.

* `inferences(method = "bootstrap")` on averages within `by` groups now
  stops before refitting when a resample leaves out every patient of a group,
  naming the group and how many resamples lack it. Previously that group's
  interval silently rested on fewer resamples. Use `method = "fwb"`, which
  keeps every patient, or combine small groups.

* `avg_time()` is renamed `avg_state_time()`, and its result class
  `markov_avg_time` is now `markov_avg_state_time`. The old name has been
  removed, so calls to `avg_time()` must use `avg_state_time()`.

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

* New `comparisons()` returns a separate comparison for each prediction
  profile (by default, each fitted patient): SOP differences or ratios,
  time in states, or time benefit, with the compared variable set to each
  requested value for that same profile. Averaging these per-profile
  differences gives `avg_comparisons()`, and `comparisons(by = ...)` returns
  `avg_comparisons(by = ...)`. `rmsb::blrm()` fits get posterior intervals
  from paired draws. For frequentist fits, `inferences()` adds intervals with
  `method = "mvn"`, `"score_bootstrap"`, or `"fwb"`, computing each
  comparison within every draw; for the stored patients, score-bootstrap and
  FWB draws carry each patient's weight (`score_weight` or `fwb_weight`), and
  weighting a draw's comparisons by it reproduces the `avg_comparisons()`
  draw. Intervals are also available with `method = "delta"` (conditional
  variance) for differences in SOPs or time in states; `vcov()` returns the
  covariance between rows. `plot_comparisons()` draws the comparisons of one
  profile and explains how to plot or average several.
  `vignette("analytical-confidence-intervals")` has a new section on
  comparisons for each patient.

* New `state_time()` returns each prediction profile's expected time in each
  state set, under each value of `variables` if given; its mean over the
  profiles is `avg_state_time()`. `inferences()` adds intervals with the
  same methods as for `comparisons()`.

* `sops()` gains `variables`, as marginaleffects' `predictions()`: each
  profile is predicted once for every value, `rowid` numbers the predictions,
  and the new `rowidcf` column gives the profile. Score-bootstrap and FWB
  draws give every copy of a patient's profile that patient's weight.

## Minor improvements

* `sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
  and `avg_comparisons()` now warn when `times` skips fitted visits. SOPs are
  calculated one transition per value of `times`, so with daily visits
  `times = c(1, 7, 14)` gives the probabilities after 1, 2, and 3 transitions,
  labelled 1, 7, and 14, not the probabilities on days 7 and 14. To get SOPs
  at selected visits, request every visit, such as `times = 1:14`, and keep
  the rows you need. The functions also warn when `times` is not in
  increasing order, and, for the starting profiles stored on a model fitted
  by `orm_markov()`, `blrm_markov()`, or `vglm_markov()`, when `times` does
  not begin at the first follow-up visit. There is no warning when `gap_var`
  is supplied, because each transition then uses its own elapsed time, or for
  times after the last fitted visit. The help for `times` now explains that
  it lists consecutive visits.

* `sim_markov("actt1")` and `sim_markov("actt2")` are much faster for large
  samples, because they now calculate each day's linear predictor for all
  patients in one step instead of one patient at a time, as
  `sim_markov("violet")` already did. The trajectories for a given seed are
  unchanged. A model that you supply through `lp_function` is still called
  once per patient.

## Bug fixes

* Wald intervals from simulation or bootstrap draws (`conf_type = "wald"` in
  `inferences()`) are now centered on the point estimate for SOP and
  time-in-state comparisons from `avg_comparisons()` and for SOPs from
  `interpolate_sops()`. These intervals were previously centered on the mean
  of the draws, which differs from the point estimate, so the reported
  estimate was not at the middle of its interval. Standard errors and
  percentile intervals are unchanged.

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
