# mostr 0.2.1

## Breaking changes

* Exported arguments now use one name for one concept (#12). The old names
  are gone, with no deprecation period, because no version has been
  released yet:

  * `prepare_markov_data()`: `absorbing_state` is now `absorb`, as in
    `mostm()`, and it is required. The old default, `6`, silently removed
    the rows after state 6 even in data where state 6 is not death. Write
    `absorb = 6` to keep the old behaviour, or `absorb = NULL` to keep every
    row. `prepare_markov_data()` also gains `y_var` and `p_var` to name the
    state columns, which were fixed to `y` and `yprev`.
  * `sim_markov()`: `absorbing_states` is now `absorb`, and `follow_up_time`
    is replaced by `times`, the days to simulate, as in `sops()`. Write
    `times = 1:28` for `follow_up_time = 28`; this gives identical
    trajectories. Each value of `times` is one transition, so
    `times = c(1, 7, 14)` simulates three transitions.
  * `plot_lp_difference()`: `x` is now `model` and `profile` is now
    `newdata`, as in the prediction functions.
  * `robcov_vglm()`: `fit` is now `model`.
  * `bootstrap_model_coefs()`: `n_boot` is now `n_draws`, as in
    `inferences()`.
  * `time_in_state()`: `target_states` is now `state_sets`, as in
    `state_time()`. It still takes one set of states.
  * `plot_time()`: `type` is now `geom`, as in `plot_sops()` and
    `plot_comparisons()`.

* `compare_se_orm_vglm()` is no longer exported.

* The `violet_baseline` dataset is no longer exported or documented; it is
  now internal data that `sim_markov("violet")` draws its patients from.
  `sim_markov("violet")` draws exactly the same patients as before. Code that
  used `violet_baseline` or `data(violet_baseline)` directly can simulate
  patients with `sim_markov("violet")` instead.

* A second-order model is declared once, when it is fitted:
  `mostm(..., p2_var = "ypprev")` names the column with the state two visits
  back. `mostm()` also takes `p_var`, the previous-state column (`"yprev"` by
  default). The fit records these with `id_var` and `time_var`, and its
  states in fitted order, and every model-based function uses them. Before,
  leaving out `p2_var` for a second-order fit, or `time_var` with `newdata`
  for a fit whose time column is not `"time"`, silently held that column at
  each profile's starting value and gave wrong probabilities.

  The absorbing states and the column with the time since the previous visit
  are declared when the model is fitted too: `mostm(..., absorb = "8")` and
  `mostm(..., gap_var = "gap")`. `mostm()` checks them against the data:

  * the `p_var` column must be in `data`, also when it is the default
    `"yprev"`, so that `absorb` is checked against the real previous states;
  * every `absorb` label must be a state of the response;
  * no fitted row may have an absorbing state as its previous state, or, for
    a second-order model, as the state two visits back (remove such rows
    first; `prepare_markov_data(absorb = )` does this for the
    previous state);
  * a state that occurs as an outcome but never as a previous state, and is
    not named in `absorb`, gives a warning: it may be absorbing, such as
    death. Before, a model with a numeric previous state (for example
    `rms::rcs(yprev, 4)`) fitted without `absorb` let dead patients move back
    to other states without any message;
  * the `gap_var` column must be numeric and, on every row after a patient's
    first, equal the time since the patient's previous row, so that a gap in
    another unit than the time column is caught. `gap_var` needs numeric
    time.

  As before, the first predicted transition gets the first requested time as
  its gap, the time since time 0. Absorbing states are never worked out from
  the data, and they can't be changed at prediction time.

  `sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
  `avg_comparisons()`, `soprob_markov()`, and `plot_lp_difference()` no
  longer have the arguments `id_var`, `time_var`, `p_var`, `p2_var`,
  `y_levels`, `absorb`, and `gap_var`; passing one is an error. Drop them from
  these calls, and give `absorb` and `gap_var`, and `p2_var`, `p_var`, and
  `time_var` when they differ from the defaults, to `mostm()` instead. The
  `state` column of the results now always holds the fitted state labels (for
  example `"1"`, `"2"`, ...); before, it took the type of a supplied
  `y_levels`, such as integers for `y_levels = 1:6`.

  `time_covariates` is removed from these functions too. Write time terms
  inline in the model formula instead, for example `rms::rcs(time, 4)`,
  `log(time)`, `rms::pol(time, 2)` or `as.numeric(time > 7)`; the prediction
  functions evaluate them at each requested time. They give the same results
  as the precomputed columns that `time_covariates` supplied. For a model
  with factor time there is now no way to give each visit's gap, so
  `gap_var` with factor time is an error.

  `plot_transitions()`, `plot_correlation()`, and `plot_variogram()` keep
  `id_var`, `time_var`, `p_var`, `p2_var`, and `y_levels`, because they also
  plot observed data, where nothing else names the columns. For a fitted
  model they default to the fit's values, and a value that differs from the
  fit's is an error, now also a `y_levels` with a different number of states.
  They no longer have `absorb`, `gap_var`, and `time_covariates`, which only
  model-based plots used; these come from the fit.

* With `engine = "orm"` or `"blrm"`, `mostm()` refuses `splines::ns()`,
  `splines::bs()` and `stats::poly()` (other than `raw = TRUE`) in the
  formula, for time and for any other variable. rms evaluates these terms
  again on the prediction rows, which gives a different basis without any
  error: with `ns(age, 3)` the average probabilities of 40 profiles were off
  by up to 0.04. Use `rms::rcs()` and `rms::pol()`, which keep the fitted
  knots and terms. `engine = "vglm"` keeps the fitted basis and accepts them.

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

* `time_in_state()` no longer accepts a data frame with `boot_id`, `tx` and
  `state_*` columns, which no function of the package has produced for some
  time. Pass SOPs from `sops()` or `avg_sops()`, after `inferences()` if you
  want intervals; any other data frame is refused.

* The deprecated arguments `adjust` of `robcov_vglm()` and `get_vcov_robust()`
  and `parallel` of `bootstrap_model_coefs()` are removed. Use `cadjust` and
  `workers`. The value of `robcov_vglm()` no longer has an `adjust` element;
  read `cadjust`.

* `mostm()`, `comparisons()` and `avg_comparisons()` no longer name the
  replacement when they are given `start_time`, `origin_time`, `origin` or
  `time_map` from earlier versions.

* Results of `inferences()` no longer carry an `n_boot` attribute; `n_draws`
  holds the same number.

* The test that `inferences(null = )` adds now agrees with the interval it is
  reported next to (#57). Before, it was always a Wald test on the scale of
  the estimate, `(estimate - null) / std.error`, so p could be below 0.05
  while the 95% interval contained the null, or the other way round. For a
  ratio, the result also depended on which level was the reference.
  * Percentile intervals, the default for `method = "mvn"`, `"bootstrap"`
    and `"fwb"`, now come without `statistic`, `p.value` and `s.value`, as
    in the marginaleffects package. So `avg_comparisons()` and
    `comparisons()` results from these methods no longer have p-values by
    default. Use `conf_type = "wald"` to get a Wald interval and a test that
    agrees with it. With `return_draws = TRUE`, `get_draws()` gives the draws
    for any other summary. A numeric `null` given with a percentile interval
    now gives a warning instead of a test.
  * SOPs from `method = "delta"` with the default logit interval are tested
    on the logit scale, the scale that interval is built on. `null` must
    then be strictly between 0 and 1. `conf_type = "wald"` keeps the test on
    the probability scale.
  * In simulations of 1344 ratio comparisons, the old test disagreed with the
    percentile interval in 50 of 2688 tests. For delta SOPs tested against
    `null = 0.05`, it disagreed with the logit interval in 87 of 4480.

## Minor improvements

* When `rms::orm()` cannot fit a model, the error now suggests
  `opt_method = "LM", maxit = 100` together, and merging rarely observed
  states. Since rms 8.2-0, the `"LM"` algorithm often needs more than the
  default 30 iterations.

* Every error that mostr raises now has class `mostr_error`, so you can catch
  them with `tryCatch(..., mostr_error = function(e) ...)`. The "Error in" line
  names the mostr function you called, such as `sops()`, instead of an
  internal helper. The messages themselves are unchanged.

* Large predictions use less memory and are faster: the options
  `mostr.max_batched_design_rows` (now `5e4` rows per visit, was `2e6`) and
  `mostr.second_order_working_bytes` (now 32 MiB, was 256 MiB) have smaller
  defaults. They change only how the work is split, not the results.

## Documentation

* Every exported function now has help-page examples that run (#26). Each
  example simulates a small trial with `sim_markov()`, prepares it with
  `prepare_markov_data()`, fits it with `mostm()` and then shows the
  function, so it can be copied and run as it is. Only slow parts, the
  bootstraps that refit the model and the `blrm` fit in `?mostm`, are inside
  `\donttest{}`. `soprob_markov()` has an example for the first time.

* `?inferences` now explains that a bootstrap with `newdata`, even when it
  holds the fitted patients themselves, keeps the profiles fixed: its
  intervals match conditional variance and are narrower than those of the
  same call without `newdata`, which match unconditional variance.

* The new help page `?mostr_options` describes the options that limit memory
  use, their defaults, and when to change them. It replaces
  `?mostr.delta_max_bytes`, which now opens the same page.

* The help pages of `avg_sops()`, `avg_state_time()`, `avg_comparisons()`,
  and `mostm()` now say that the averages give every profile the same weight:
  case weights given to `mostm(weights = )` are used only to fit the model,
  so the result is not the same as repeating each patient according to their
  weight.

## Bug fixes

The fixes below change results that earlier versions returned, mostly without
a warning; recompute analyses that match one of them. The linked issues and
the git history give the details.

* `blrm` fits:
  * With `include_re = TRUE`, random intercepts are matched to patients by ID.
    Whole-number IDs could get other patients' intercepts, and
    `cluster(site)` looked up patient IDs among the sites (#39).
  * With `posterior_summary = "median"`, `time_in_state()` and
    `interpolate_sops()` take the median of each draw's total or interpolated
    value instead of adding or interpolating medians (#40).
  * With `ppo` and `cppo`, predictions use the same `cppo` scores as the fit
    (#46).

* `vglm` fits honour a `parallel` formula that names an inline spline, such
  as `rms::rcs(time, 4)` (#42).

* Variance and intervals:
  * The robust variance of penalized `orm` fits is no longer too large (#49).
  * In the patient-clustered and unconditional variance of `orm` fits, rows
    whose observed state the model makes almost impossible no longer
    contribute 0 (#48).
  * `vcov()` on delta-method results follows the rows of a filtered or sorted
    result (#43).
  * Unconditional variance takes each row's patient from the data stored on
    the fit, so reordering the data frame after the fit no longer changes it.
  * `time_in_state()` on delta-method SOPs keeps their standard errors and
    intervals, and a repeated `inferences()` call fully replaces the earlier
    one.

* Bootstrap and FWB:
  * Refits of resamples that lack a state keep the state labels. They were
    renumbered, so profiles were predicted from the wrong state.
  * Each resampled patient starts from their row at `first_followup_time`,
    also for character or factor time (#47).
  * Resamples in which patients start in a state that nobody enters are kept
    instead of failing, which made the lower bound for that state too high
    (#54).
  * `bootstrap_model_coefs()` gives `NA` for the intercepts and
    previous-state coefficients of a resample without the reference previous
    state, whose refit changed their meaning (#54).

* Real-time results count a starting state that never occurs as an outcome
  (#50).

* `plot_transitions()`, `plot_correlation()` and `plot_variogram()` follow
  the fit's visits instead of daily visits from 1 (#51), and score integer
  states correctly when one is missing (#52).

* Comparisons:
  * `by` naming the compared variable groups patients by their observed
    value. SOP and time-in-state comparisons ignored it.
  * A comparison of a variable the model doesn't use reports its uncertainty
    as `NA` instead of 0.
  * A ratio draw with a zero denominator counts as `Inf` instead of being
    dropped (#57).

* `sim_markov()` reads a factor `yprev` in `baseline_data` by its labels
  (#41).

Inputs that gave wrong results are now errors that name the problem: the VGAM
arguments `xij` and `form2` (#45); `strat()` terms (#59); random-effect terms
that mostr cannot use; a delta-method `cluster` that does not group the rows
by patient (#53); the time column in `variables` (#29); `inferences()` on
interpolated SOPs (#58); missing values in prediction columns, which `vglm`
fits dropped; and bootstrap resamples that a second-order model cannot predict
from (#54, #55). Other fixes remove errors from valid calls, for example with
factor predictors of `vglm` fits (#44), matrix columns in `newdata` (#56), or
a covariate in very small units in the robust and unconditional variance.

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
