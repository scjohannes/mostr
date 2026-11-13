# Decision record

What was decided in `mostr`, why, and which options were rejected and why.
Newest entries first. Each entry has a date, the decision, the reason, and
the options decided against. Where no reason was written down, the entry says
"reason not recorded" instead of guessing one.

This file is tracked by git and left out of the built package
(`.Rbuildignore`). Scripts and raw results may live in the gitignored folders
`design/`, `benchmarks/local/` and `scripts/`, but those stay on one machine, so
the conclusions and the numbers they rest on belong here. Commit hashes refer
to the branch where the work was done; issue and PR numbers refer to the
project's GitHub repository.

## 2026-11-13: Patients without a fitted row are named at the fit and left out of resamples

**Decided.** `mostm()` warns when a patient loses every row to missing values
in model variables and names those patients. The stored refit data holds only
patients with at least one fitted row, so the bootstrap and the fractional
weighted bootstrap resample the cohort the estimate uses.

**Why.** The estimate already averaged over fitted patients only, while each
resample brought the other patients back; prediction then refused their rows
and every draw failed (80 patients, 3 without `age`: 20 of 20 draws). The
maintainer wants the user told when the model is fitted, not at the bootstrap.

**Details.** A missing response or previous state on every row counts as well
as a missing predictor. Patients who lose only some rows are not named.
Patients excluded by `subset` are not named, because the user left them out.

## 2026-11-13: Two review findings left as they are

**Decided.** Two findings from the review of `feature/bartisan-lessons` are
known and not fixed, because neither is serious enough to act on now.

1. **Failed simulation draws in comparison results are not summed up.** With
   `method = "mvn"`, `avg_comparisons()` (time benefit) and `comparisons()`
   drop a draw whose prediction fails, with one warning per draw that names
   the internal `soprob_markov` and no summary
   (`inferences_avg_comparisons_time_benefit_mvn()` in
   `R/sops-comparisons-inference.R`, and the matching loop in
   `R/sops-comparisons-individual-inference.R`). The bootstrap paths count
   failed draws in one warning (`keep_successful_draws()`). No realistic
   trigger is known: crossed cut-points, the one common per-draw failure, stop
   the call and are counted separately, and a finite covariance matrix gives
   finite coefficient draws. If a trigger turns up, send these loops through
   `keep_successful_draws()`.

2. **Second-order models check crossed cut-points for different histories on
   the two prediction paths.** The R recursion
   (`soprob_markov_second_order_run()` in `R/sops-engine.R`) predicts and
   checks only pairs of (state two visits back, previous state) that some
   profile can be in at that visit. The compiled path
   (`R/sops-execution-plan.R`) and the diagnostic trace predict and check
   every pair. They differ only at the second visit: there, a pair whose
   older state no profile starts in has probability 0, while from the third
   visit on every pair of non-absorbing states has some probability. So a
   call with two time points can stop on the compiled path and return SOPs on
   the R path, and a longer call stops on both but names different visits.
   Only models with non-proportional terms are affected. If this is taken up,
   make the compiled path skip pairs without probability, which also saves
   their computation. Read from the code on 2026-11-13, not reproduced by a
   run.

## 2026-11-13: Keep the new regression tests on CRAN at the edge of the time budget

**Decided.** The tests added on 2026-11-12 and 2026-11-13 that run on CRAN
stay there, although the CRAN-mode suite is now at its 60-second budget.

**Why.** The maintainer chose to keep them. Measured on a quiet machine
(Ryzen AI 9 HX 370): 57 s and 65 s with the AGENTS.md command, and 56 s in the
tests step of `R CMD check --as-cran` on a built tarball, whose whole check
took about 6.5 minutes. The 30 tests added since 12 November that run on CRAN
cost 6.1 s together; the slowest (`time_in_state()` on delta and blrm SOPs, a
second-order prediction check and a real-fit `plot_lp_difference()` test) take
0.4 to 1 s each.

**Decided against.** Skipping the slowest of them on CRAN, which would have
saved about 4 s but removed regression tests for bugs fixed this week from
CRAN's runs.

## 2026-11-13: Check prediction rows where they enter the engines

**Decided.** Prediction rows are checked and put in the form the model reads
them in one helper, `prepare_prediction_rows()` (`R/sops-dispatch.R`), called
by `predict_markov_sops()`, which every SOP route passes through, and by
`compile_sop_execution_plan()`, which the analytical inference paths call
directly. It refuses missing values in the columns the fitted formula uses
(`stop_missing_predictors()`), and turns the previous-state columns into the
model's fitted states and type (`fitted_state_values()`), refusing values
that are not fitted states. `plot_lp_difference()` makes the same check on
its `profile`. vglm designs are built with `na.action = na.pass` and refused
when they hold a missing value (`vglm_prediction_design()`), as a second line
of defence.

**Why.** A missing covariate made vglm drop the row and return another
profile's probabilities, and a profile factor with only its own states made
the recursion's later previous states NA on every engine. The rows are
checked after `variables` has set its columns, because a column that
`variables` sets may be missing in `newdata` without harm. The previous-state
columns are converted only in the rows handed to the engines, so the
profile columns that results copy keep the user's types.

**Decided against.**
- Checking in `resolve_markov_source_data()`, where absorbing starts are
  refused: it sees `newdata` before `variables` sets its columns, so it would
  refuse `avg_comparisons(variables = "tx")` on profiles with a missing `tx`.
  The cost is that the error counts prediction rows, which `variables`
  multiplies, rather than the user's profiles.
- Taking the type and levels from the model at each of the dozen places the
  recursion builds a previous-state column: converting the starting rows once
  fixes them all, and new engines can't miss it.
- Converting the columns in the user's `newdata`: results copy the profile
  columns, so a character `yprev` would come back as a factor.

## 2026-11-13: Remove the score bootstrap

**Decided.** The score bootstrap, hidden from `inferences()` since 2026-11-09,
is deleted with its tests: the internal `inferences_score_bootstrap()`,
`generate_score_bootstrap_draws()`, and the `score_weight` draw column. MVN
draws no longer carry patient weights. `inferences(method = "score_bootstrap")`
still stops, with an error that says the method was removed and names the
remaining ones. This carries out the retirement decided on 2026-11-12 in
"Leave the user's parallel plan alone; retire the score bootstrap".

**Why.** The maintainer was never sure that the score bootstrap has a sound
theoretical basis, or where it would fail and where it would not, so it was
never planned for release.

**Decided against.** Keeping the code as an unexported function, as it was
from 2026-11-09. Reason not recorded beyond the doubts above.

## 2026-11-13: Move `absorb` and `gap_var` into `mostm()`, and remove `time_covariates`

**Decided.** The absorbing states and the column with the time since the
previous visit are given once, to `mostm(absorb = , gap_var = )`, recorded in
the fit's settings and checked when fitting. The prediction functions
(`sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`, `comparisons()`,
`avg_comparisons()`, `soprob_markov()`, `plot_lp_difference()`) and the model
plots (`plot_transitions()`, `plot_correlation()`, `plot_variogram()`) lose
`absorb`, `gap_var` and `time_covariates`; the plots only used them for
fitted models. `time_covariates` is removed: time terms are written inline in
the formula. The maintainer decided:

- Absorbing states are never worked out from the data, and there is no
  per-call `absorb`: first-passage probabilities (treating a non-absorbing
  state as absorbing at prediction time) are not wanted.
- The first transition keeps `times[1]` as its gap, the time since time 0,
  also when supplied `newdata` starts at a later visit. The rule is accepted
  and documented in `?mostm` and in "Choosing visits" of `?sops`.
- With factor time there is no longer a way to give each visit's gap, so
  `gap_var` with factor time is refused.

`mostm()` refuses `absorb` labels that are not states of the response and
fitted rows whose previous state, or for a second-order fit the state two
visits back, is absorbing. It warns (class
`mostr_warning_outcome_only_states`) about a state that occurs as an outcome
but never as a previous state and is not in `absorb`; it does not change the
fit. It refuses a `gap_var` column that is not numeric, or that differs from
the time since the patient's previous row on any row after a patient's first
(relative tolerance `sqrt(.Machine$double.eps)`). With `engine = "orm"` or
`"blrm"` it refuses `splines::ns()`, `splines::bs()` and `stats::poly()`
(other than `raw = TRUE`) for any variable.

**Why.** Two silent wrong answers:

- With a numeric previous state (for example `rcs(yprev, 4)`) and `absorb`
  left out for death, dead patients moved back to state 7: in the simulated
  trial measured on 2026-11-12, the day-14 probability of death was 0.012
  instead of 0.057. A factor previous state stops at prediction because the
  model has no coefficient for death as a previous state; a numeric one
  extrapolates without a message. Recording `absorb` on the fit lets
  `mostm()` warn about it once, where the data are at hand.
- rms evaluates `ns()`, `bs()` and orthogonal `poly()` again on the
  prediction rows and gets another basis, without an error. Measured on
  2026-11-13 (`benchmarks/local/settings-refactor/spline-covariates.R`, 300
  simulated patients, `ns(age, 3)` etc. for a covariate, 40 supplied
  profiles): the orm design matrix differed from the fitted one by up to
  0.42 (`ns`), 0.38 (`bs`) and 0.50 (`poly`), and the average SOPs by up to
  0.042, 0.033 and 0.038; one profile with `poly()` failed with an error.
  `rms::predictrms()`, which the blrm route uses, gave the same differences.
  `rms::rcs()` and `rms::pol()` agreed exactly. With `engine = "vglm"`, which
  keeps the fitted `predvars`, the SOPs agreed to 3e-16
  (`spline-covariates-vglm.R`), so vglm is not refused.

Inline time terms replace `time_covariates`. Measured on 2026-11-12 (SOPs and
delta standard errors compared with precomputed columns passed through
`time_covariates`): identical for `log(time)`, `rms::pol(time, 2)`,
`as.numeric(time > 7)`, `rcs()`, `lsp()`, and vglm with `splines::ns(time, 3)`.
orm with `splines::ns(time, 3)` fitted but prediction failed with "all
interior knots match left boundary knot", because rms evaluates `ns()` again
on the new data (now refused at fit time); `I(time > 7)` makes `rms::orm()`
itself fail. A `gap_var` that is checked against the time column catches a
gap recorded in another unit than time, which would otherwise be fitted and
predicted on different scales.

The 97 saved outputs of
`benchmarks/local/settings-refactor/baseline-before-move.R` (on `40480f3`)
and of its translated copy `baseline-after-move.R`, which passes `absorb`
and `gap_var` to `mostm()` and replaces the `time_covariates` spline columns
by `rms::rcs(time, 3)` inline, are identical with `identical()`
(`compare-move.R`), including the SOPs, draws and covariances of the
inline-term fit. The only difference in the full objects is the coefficient
names stored with that fit's delta-method result (`time`, `time'` instead of
`time_lin`, `time_nlin`). The gap fit in both scripts uses data whose gap
equals the time since the previous row; the grid data's own `gap` column did
not, so the test grid's gap set-up now spaces its visits at times 1, 2, 4, 5,
7, 8 with matching gaps.

**Decided against.** Inferring absorbing states from the data (a state seen
only as an outcome): a rare state seen only at a last visit would be
declared absorbing without the user knowing, so the check only warns. A
per-call `absorb` for first-passage probabilities: not wanted, and it would
predict from transitions the model was fitted on as if they could not
happen. Keeping `time_covariates` for factor-time gaps: the only remaining
use, and it kept a second way to set visit-specific columns that the
execution plans, delta method, bootstrap replays and diagnostic traces all
had to carry. Refusing `ns()`, `bs()` and `poly()` only on the time variable:
the measurement showed the same wrong answer for any covariate. Changing the
fit when the warning fires: the user decides whether the state is absorbing.

## 2026-11-13: One positive-definiteness rule for the orm and vglm robust variance

**Decided.** `orm_model_bread()` and the vglm checks of `vcov(fit)`, the
chosen bread and the observed information use one helper,
`markov_correlation_positive_definite()` in `R/robcov_orm.R`. It refuses a
matrix with a non-positive diagonal entry and otherwise checks the
correlation matrix: the smallest eigenvalue must exceed p * 1e-10 * max(1,
largest), the orm rule. `robcov_vglm()` now inverts `vcov(fit)` with
`chol2inv(chol())` instead of `solve()`.

**Why.** vglm checked the raw matrix with tolerance 100 * eps * p * max(1,
largest), so the two engines disagreed (measured with 3 coefficients):

- a covariance with correlation eigenvalues 2.07 and 1.3e-11 (built from
  eigenvalues 2, 1 - 1e-11 and 1e-11): orm refused, vglm accepted;
- `diag(1e-16, 3)`: orm accepted, vglm refused, although its correlation
  matrix is the identity.

On a real fit (40 simulated patients, `ordered(y) ~ time + tx + yprev + bmi`,
12 coefficients) the smallest correlation eigenvalue of `vcov(fit)` is 0.0039,
5e5 times above the cut-off, and it stays there when `bmi` is multiplied by
1e-6 or 1e6; the old raw rule refused both. With `bmi` multiplied by 1e-6 the
vglm fit now works and its conditional delta standard errors match the unit
scale to a relative 1e-10. With `bmi` multiplied by 1e6, `solve()` refused
`vcov(fit)` (reciprocal condition number 2e-18 on the raw scale), hence the
Cholesky inverse; that fit still stops, now at the observed information,
because the finite-difference step of the observed bread,
eps^(1/3) * max(1, |coefficient|), is 6e-6 for a coefficient of 1e-7 and moves
the linear predictor by up to 157. That step is a separate problem, not changed
here.

**Decided against.** Keeping the raw rule for vglm: it makes the result depend
on the units of the covariates. Two separate copies of the correlation rule:
the two had already drifted apart once.

## 2026-11-12: Keep the blrm options as separate arguments

**Decided.** `n_draws`, `seed`, `include_re`, `return_draws` and `conf_level`
stay separate arguments of the prediction functions instead of moving into
one control list.

**Why.** marginaleffects, whose interface mostr follows, keeps such options
as separate arguments, and a list would hide them from autocompletion and
from the help page's argument list.

**Decided against.** A `blrm_control()` list: it would shorten the argument
lists but make the options harder to find, for no gain in what users can do.

## 2026-11-12: Remove the column-name and state arguments from the prediction functions

**Decided.** `sops()`, `avg_sops()`, `state_time()`, `avg_state_time()`,
`comparisons()`, `avg_comparisons()`, `soprob_markov()` and
`plot_lp_difference()` no longer have `id_var`, `time_var`, `p_var`, `p2_var`
and `y_levels`. They take them from the settings that `mostm()` records on the
fit, and passing one is an error that says so. `plot_transitions()`,
`plot_correlation()` and `plot_variogram()` keep them, because they also plot
observed data, where nothing else names the columns; for a fit they still
default to the fit's values and refuse a contradiction, now also a `y_levels`
with a different number of states.

**Why.** Since the settings object, these arguments could only repeat what
the fit records: every model-based function accepts only `mostm()` fits, and a
different value was refused. Keeping them meant documentation, contradiction
checks and argument threading for no choice the user could make. Removing them
before the first CRAN release needs only a NEWS entry; after it, a deprecation
cycle. The test fixtures that are not `mostm()` fits now record settings as
`mostm()` does (`with_fit_settings()` in `tests/testthat/helper-mocks.R`).

The 97 saved outputs of `benchmarks/local/settings-refactor/baseline.R`
(before, on `6c4f0c4`) and of its translated copy
`baseline-after-remove-args.R` are identical with `identical()`, except one:
the call that passed `y_levels = 1:4` now returns the fit's state labels, a
factor of `"1"` to `"4"`, in `state` instead of integers; its estimates and
all other columns are identical (`compare-remove-args.R`).

**Decided against.** Deprecating the arguments first: the package has no
released users. Keeping `y_levels` for the type of the `state` column: the
type can be changed after the call, and keeping one argument would keep the
contradiction checks for it. Removing the arguments from the three data plots
too: for observed data they are the only way to name the columns.

## 2026-11-12: Size limits measured; two defaults open

**Decided.** The size limits stay as they are for now, with comments that
state what each bounds and what was measured. Whether to change two of them is
open (below).

**Why.** None had a recorded reason; all came over from markov.misc without a
measurement. Measured Nov 2026 on an installed build (Ryzen AI 9 HX 370, 32 GB;
scripts and tables in `benchmarks/local/thresholds/`):

- `mostr.max_batched_design_rows = 2e6`: the batched design build is faster
  only below about 25k rows per visit; from 50k it is 1.3 to 1.9 times slower
  with up to 1.8 times the peak memory (it subsets one data frame with
  duplicated rows, which makes unique row names). Proposed: 5e4. Open.
- `mostr.second_order_working_bytes = 256 MiB`: the chunking counts only the
  transition probabilities, so peak use is about 6 times the budget, and
  smaller chunks cost no time. Proposed: 32 MiB, or counting each pair's design.
  Open.
- `mostr.execution_plan_max_bytes`, `mostr.delta_max_bytes` and
  `mostr.posterior_working_bytes` (256 MiB each): reasonable. The plan budget
  is conservative: above it MVN inference ran 3.6 times slower in the R engine,
  and delta inference stops.
- `mostr.blrm_chunk_size = 10`: a fallback no package caller reaches; 10 is
  within about 8% of the fastest chunk size.

**Decided against.** Changing defaults inside a comment-only change: a changed
default needs its own check that results are identical.

**Update 2026-11-13.** The maintainer decided to change both:
`mostr.max_batched_design_rows` is now `5e4` and
`mostr.second_order_working_bytes` 32 MiB. All the size options are now
documented on one help page, `?mostr_options` (before, only
`mostr.delta_max_bytes` had a page); the blrm chunk sizes stay undocumented
because they only tune speed. Check (`benchmarks/local/thresholds/identical-defaults.R`,
saved before and after the change): a first-order orm call with 1e5 design rows
per visit (batched before, one previous state at a time after), a 1e4-row
control, and delta and MVN inference on the large call gave `identical()`
outputs. A second-order orm call forced to the R engine with 20,000 profiles
(36 state pairs per chunk before, 11 after) differed by at most 5.6e-17, about
one rounding step of a double, because the sums over state pairs are added in
another order.

## 2026-11-12: Recommend `"LM"` only with `maxit = 100`

**Decided.** When `rms::orm()` fails, the error suggests
`opt_method = "LM", maxit = 100`, `scale = TRUE`, removing nearly collinear
covariates, or merging rare states. The README example passes `maxit = 100`.

**Why.** rms 8.2-0 changed LM to relax its damping by 1.5 instead of 10 after
an accepted step (upstream: the old rule could cycle on an ill-conditioned
data set). LM then often needs more than the default 30 iterations. On
simulated ACTT-2 data (25 seeds per setting), LM with `scale = TRUE` failed in
up to 16 of 25 small data sets and 4 of 100 README-sized ones under rms 8.2-1,
and never under 8.1-1. With `maxit = 100` it converged everywhere, in at most
98 iterations. Under 8.1-1 the default Newton-Raphson algorithm failed on small
data sets with no finite estimates, where LM succeeded, so LM stays the
suggested alternative. Scripts: `benchmarks/local/lm-convergence/`.

**Decided against.** Dropping the LM suggestion: it still rescues fits the
default algorithm cannot. Recommending LM without `maxit`: it fails often
under rms 8.2. An upstream report to rms is drafted but not filed.

## 2026-11-12: Cluster-corrected robust variance by default

**Decided.** Keep `cadjust = TRUE` (multiply by G / (G - 1), G = number of
patients) as the default for patient-clustered variances, with `type = "HC0"`.

**Why.** It is the default of `sandwich::vcovCL()` for models other than `lm`;
mostr's four type and cadjust combinations agree with `vcovCL()` on a binary
logistic model to glm's convergence tolerance. It also makes the conditional
variance equal `cov(influence) / G`, the convention of the unconditional
variance, which the identity test relies on. `rms::robcov()` applies no factor,
so mostr's standard errors are sqrt(G / (G - 1)) times larger than rms's.

**Decided against.** Following `rms::robcov()` (no factor). Reason recorded
only now; the default dates from the first commit.

## 2026-11-12: Leave the user's parallel plan alone; retire the score bootstrap

**Decided.** The bootstrap currently replaces the user's future plan with
`callr` for the run and restores it afterwards. It should not touch the user's
plan; this is fixed together with the switch to mirai (issue #27), not
separately. The score bootstrap, already hidden from `inferences()`, will be
retired, so its probability clipping gets no comment.

**Why.** Changing global state inside a function, even temporarily, surprises
users and interacts badly with their own parallel set-up; the backend switch
rewrites the same code anyway.

**Why `callr`.** The maintainer chose `future.callr::callr` over other future
backends because with them parallel runs ran out of memory and workers did not
shut down quickly enough; with `callr` this happened less. Not measured: this
rests on the maintainer's experience.

**Restoring the plan is CRAN-safe until then (checked 2026-11-12).** All three
places that switch the plan (`apply_to_bootstrap()`, `apply_to_fwb_bootstrap()`
and `apply_sop_simulation_draws()`) register `on.exit()` before running, so
the plan is restored after an error inside a worker too. Tried with a
sequential, a `multisession` (2 workers) and a `callr` (3 workers) plan, each
with an error in one draw. One gap was found and fixed: they saved the plan
with `plan()`, which returns only the outer level, so a nested plan such as
`list(multisession, multisession)` came back with its inner level dropped.
They now keep the value of `plan(callr, ...)`, which is the whole previous
plan. Tests use at most 2 workers, within CRAN's limit.

## 2026-11-12: Ship one vignette; the other articles are website-only

**Decided.** The package ships only `vignettes/mostr.Rmd`, the get-started
guide. The nine longer analyses are pkgdown articles in `vignettes/articles/`,
which `.Rbuildignore` excludes, so CRAN never builds them. Help pages, the
README and messages link to their website URLs.

**Why.** One `R CMD check --as-cran` with vignettes took 13.7 minutes on the
development machine (busy with other R processes), against CRAN's limit of
about 10. Rebuilding the vignettes took 7.1 minutes of that; four of them took
70 to 92 seconds each. Tests in CRAN mode were already cut to about 50 seconds,
so vignettes were the remaining cost. Many packages, the tidyverse among them,
ship few vignettes and keep the rest on the website.

**Decided against.** Shipping saved results that the vignettes load instead of
recomputing (bartisan's approach): it keeps the vignettes installed, but the
saved results must be regenerated whenever the output changes. Gating the heavy
chunks behind an environment variable, as the Bayesian article does: an
installed vignette would then show code without its results.

## 2026-11-12: Keep a tracked record of decisions

**Decided.** Record decisions in this file and state, before a simulation or
benchmark runs, the claim, the measurement and what each outcome would mean
(rules in `AGENTS.md`).

**Why.** Earlier reasoning lived in gitignored folders, out of reach of
co-authors and cloud agents. A comparison with the bartisan package found
its tracked development records (for example its log of each CRAN rejection
and fix) useful evidence that mostr lacked.

**Decided against.** Keeping notes in `design/` only, because nobody else can
read them. A task log as long as bartisan's (about 10,000 lines): copy the
habit at a much smaller size.

## 2026-11-12: Rename the R engine and remove dead C++ routines

**Decided.** `soprob_markov_reference()` is now `soprob_markov_r_engine()`.
`cpp_markov_propagate()` and `cpp_normalize_probability_array()` are deleted.
Test-only delta helpers moved from `R/` into the test oracle helper. The Hmisc
comparison tests skip when Hmisc is not installed (`2d6a2a6`).

**Why.** The R engine runs in production (for `blrm` fits and designs too large
for the compiled path) and shares code with the compiled path, so the word
"reference" wrongly suggested an independent check. The two C++ routines had
no callers, and no planted bug targeted them.

## 2026-11-11: Validate the unconditional variance with a study, not with tests

**Decided.**
- No package tests that compare bootstrap, FWB or score-bootstrap standard
  errors with unconditional variance, and no coverage tests.
- Instead, a fast test of the exact identity: unconditional variance =
  conditional variance + the variance between patients' predictions + twice
  their covariance with the coefficient influence, all divided by n. It needs
  no random draws and runs on CRAN
  (`test-sops-delta-unconditional-identity.R`).
- A slow test (`MOSTR_RUN_SLOW=true`) of MVN draws against conditional delta
  standard errors, both given the same small covariance
  (`vcov = 1e-3 * fit$var`), with a 2% tolerance (`test-slow-mvn-delta.R`).
- A one-off validation study, rerun only when the unconditional variance code
  changes (`R/sops-delta-unconditional.R`, `R/sops-delta-inference.R`,
  `R/sops-delta-comparisons.R`, `R/sops-delta-avg-time.R`). Its script and
  results are local, in `benchmarks/local/validation-unconditional/`.

Committed in `5d1bbb0`.

**Why.** Comparisons with bootstrap and coverage simulations are too slow and
too noisy to be tests. MVN standard deviations do not converge to the delta
standard error as draws increase: MVN passes draws through the nonlinear
recursion, while the delta method uses a linear approximation. With a small
covariance the two must agree, and the Monte Carlo error of a standard
deviation from B draws is about 1/sqrt(2B) relative (0.5% for B = 20,000).

**Results of the study** (500 patients, B = 5000 bootstrap draws, 1000
repetitions, orm and vglm; acceptance bands written down before each run):
- Every band passed. Bootstrap and FWB standard errors divided by the
  unconditional delta standard error had medians 0.99 to 1.01, and every cell
  was within 0.978 to 1.035.
- Unconditional 95% intervals covered the population value 94.4% to 95.0% of
  the time on average (cells 92.9% to 96.4%).
- With starting state and treatment only, unconditional variance was within
  0.3% of conditional variance, so a second scenario with a patient covariate
  was added, with its bands written down before it ran. There, conditional
  intervals covered the population value only 91.4% of the time for SOPs and
  93.1% for contrasts (83.8% in one cell), while unconditional ones stayed at
  95%. This is why unconditional variance is the default for averages.
- The identity holds to a relative error of about 1e-14 for orm and 1e-12 for
  vglm. vglm's summed scores are 0.013, not 0, because VGAM stops iterating
  earlier; the test therefore asks for 1e-8. A planted 1% error fails it.

**Decided against.** The test-type comparisons listed above. A second-order
set-up for the MVN test, because delta inference supports first-order models
only.

## 2026-11-11: Identify states by label, never by position

**Decided.**
- Bootstrap refits keep the original state labels and numbers and only drop
  empty factor levels (`1c7a164`).
- `inferences(method = "bootstrap")` stops before refitting when some resample
  lacks a state that a supplied profile starts in, and recommends
  `method = "fwb"`.
- `variables` that set the state two visits back (`p2_var`) allow only
  conditional variance, like `variables` that set `yprev`.
- Values that `variables` gives the previous state mean the state whose label
  they match, whatever their type (`dc5fd49`).

**Why.** A resample without state 2 was refitted with states 1, 3, 4 renamed
1, 2, 3, while supplied profiles kept their labels. Profiles were then
predicted from the wrong state without any message, and `by` averages put
state 3's probability under state 2. For matching `variables` by label, the
reason is not recorded beyond the commit message.

**Decided against.** Dropping the affected bootstrap draws instead of
stopping: the maintainer chose to stop. The reason, as for per-profile
`comparisons()` (2026-10-17 entry): a supplied profile that starts in a state
missing from a resample has no prediction in that draw, and dropping the draw
would silently rest the interval on the resamples that happen to contain the
state. Checked 2026-11-13 (200 patients, a state 5 patients were ever in):
averages over the stored patients ran and gave that state probability 0 in
resamples without it, while supplied `newdata` with a profile starting there
stopped because 2 of 200 resamples lacked the state. A resample lacks a state
that m patients were ever in with probability about e^-m, so with m = 5 and
1,000 draws the call almost always stops; with m = 10 it stops about 4% of
the time. (For `by` groups,
the 2026-10-17 entry rejected dropping because the interval would silently rest
on fewer resamples.)

## 2026-11-11: Consolidate tests only when no planted bug is lost

**Decided.** Tests were merged into tables by file group, from 699 to 502 tests
(`bca6836`). The rule: run the planted-bug checks before and after, and no
planted bug may go from caught to missed. All 27 planted bugs were still
caught. These must not be merged away:
- each C++ routine tested against its own pure-R reference
  (`test-sops-native-reference.R`), and the Hmisc comparisons, which are the
  only tests that catch several planted bugs;
- finite-difference and patient-weight-refit checks, the fresh-process test of
  `mostm()`, random-number-stream tests.

**Why.** Merging can silently lose protection. Line coverage cannot show this:
it stayed at 89.3% while five more planted bugs became caught.

**Decided against.** Judging consolidation by the number of tests or by line
coverage.

## 2026-11-11: Bootstrap refits leave the rms `datadist` alone

**Decided.** `inferences()` no longer has `update_datadist`. Refits run with
the `datadist` option unset and restore it afterwards (`50ba8b2`).

**Why.** mostr never needs a `datadist`. Refits used to create a `dd` object in
the global environment, and `update_datadist = FALSE` together with a
`datadist` option naming a missing object made every refit fail.

## 2026-11-10: Choices made while building the refusal table

**Decided** (committed in `481ee0d`):
- A comparison of a variable the model does not use is 0 (1 for a ratio), and
  its standard error, test and interval are `NA` with every method, as in the
  marginaleffects package. Previously the standard error was 0 and the test
  statistic `NaN`.
- `use_coefstart = TRUE` is an error unless the fit used `engine = "vglm"`, and
  with `method = "mvn"` or `"delta"`. Only vglm refits take starting
  coefficients, and mvn and delta do not refit, so the argument was silently
  ignored.
- A numeric `yprev` whose absorbing state is left out of `absorb` is the
  user's error and is not refused. Reason not recorded.
- No check with 30 states: dropped as unnecessary.

**Deferred.** Refusing time, gap and time-covariate columns in `variables`
(issue #29, "not now": the maintainer trusts users not to put time columns
in `variables` for now). Slow `interpolate_sops()` (issue
#28).

## 2026-11-10: A CRAN release is planned

**Decided.** A first CRAN release is planned; its timing and steps are not
decided. Until then:
- a test runs on CRAN only if it passes for any seed; heavy grids and refits
  skip on CRAN or go into `test-slow-*.R` files;
- new NEWS entries go under `# mostr (development version)`, above 0.2.0;
- the examples in `\dontrun{}` must run before a first submission (issue #26).

**Why.** Adopted from bartisan, whose second CRAN submission was rejected for
running over 10 minutes on Windows. mostr's Windows check on GitHub already
takes 7.5 minutes with vignettes off, so the time budget is tight. The reason
for the any-seed rule itself is not recorded.

## 2026-11-10: Where the bartisan lessons are done

**Decided.** The test overhaul and bartisan lessons 2, 4, 5, 7, 8, 10 and 11
are done on `feature/bartisan-lessons` and get no GitHub issues. Lessons 3, 6
and 9 became issues #26 (runnable examples), #25 (bartisan as a fitting engine)
and #27 (let the user choose the parallel backend). On 2026-11-08 it was also
decided that bugs found by the overhaul are fixed on the same branch.

**Decided against.** A separate `fix/` branch per bug, as the test overhaul
plan had proposed: the maintainer's preference, which they see no need to
justify.

## 2026-11-09: Behaviour for the problems found by the test overhaul

**Decided** (committed in `c7cfac1`):
- **Absorbing starting states are refused.** First-order models predicted a
  move out of the absorbing state for such patients (0.199 at visit 1 with a
  numeric `yprev`; a factor `yprev` failed with an unrelated error), while
  second-order models kept them there with probability 1.
- **Unconditional variance is refused for every fit with case weights.** The
  calculation has not been checked with weights, and engines disagreed: orm
  refused, vglm ran.
- **Averages give every profile the same weight.** Case weights are used only
  to fit the model; this is now documented. Weighting the averages by case
  weights was not chosen, to mirror marginaleffects.
- **`include_re = TRUE` is an error except for blrm fits.** It was silently
  ignored for orm and vglm.
- **Crossed thresholds stop with an error.** In a partial proportional-odds
  model, thresholds can cross for some covariate values and the model then
  gives a negative probability, so SOPs are not defined. The compiled engine
  set such values to zero and rescaled, the vglm R engine kept the negative
  value, and blrm clipped each draw; the two vglm engines could differ by 0.5.
  Every route now stops, including any MVN, bootstrap or posterior draw.
  Values between -1e-10 and 0 come from rounding and are still set to zero.
  Clipping and rescaling, and dropping draws with a warning, were rejected.
- **The score bootstrap is removed from `inferences()`.** Its code stays as the
  internal `inferences_score_bootstrap()`. The maintainer was never sure the
  score bootstrap has a sound theoretical basis, or where it would fail and
  where not, so it was never planned for release; it is being removed
  entirely (2026-11-13).

## 2026-11-09: `by` groups

**Decided.**
- When `by` names the compared variable (`variables = list(tx = c(0, 1))`,
  `by = "tx"`), patients are grouped by their observed value for every
  estimand (`d1ed9b8`). SOPs and time in state had silently ignored the
  grouping, while time benefit already followed the `comparisons()` help.
  `avg_sops()` and `avg_state_time()` refuse `by` on a variable they set.
- A missing value in a `by` column among the averaged profiles is an error
  (`8c601b9`). Those profiles had been left out of every group, which changed
  who was averaged. Restricting `by` to model variables was rejected; reason
  not recorded.
- `inferences(method = "delta")` keeps refusing `by` for now. An investigation
  found both variances feasible (about 2 to 3 days): conditional variance from
  each group's averaged derivatives, which matched `newdata` restricted to the
  group exactly; unconditional variance from influence terms centred within
  the group, with the variance taken over all n patients, which agreed with
  score and patient bootstraps. Treating a group as if it were the whole cohort
  overstates standard errors by 30% to 90%. If built: ship both variances
  together and refuse unconditional variance for groups of fewer than two
  patients. Whether to build it is not decided.

## 2026-11-09: `mostm()` requires `id_var`

**Decided.** `id_var` is required for every engine, and code no longer falls
back to a column called `id` (`3e7868f`). orm and vglm fits always carry the
patient-clustered covariance.

**Why.** Fits without `id_var` only warned. Averages then failed when the
patient column had another name, and both bootstraps silently clustered on any
column called `id`.

**Decided against.** Keeping the warning.

## 2026-11-08: Fit orm examples with the default optimizer

**Decided.** Examples fit orm with the default optimizer, not
`opt_method = "LM"` (`4f59e46`).

**Why.** Under rms 8.2.1, `"LM"` stops after 30 iterations without converging
on the 40-patient ACTT-2 example, so R CMD check failed on every platform. The
default optimizer converges under rms 8.1.1 and 8.2.1 with the same
coefficients.

## 2026-11-07: Test plan: independent checks first

**Decided.**
- Every model set-up of a shared grid (`helper-grid.R`) is checked first
  against independent references: a readable loop that calls no package code
  and takes transition probabilities only from the fitting package's
  `predict()` (`helper-slow-loop.R`), and numerical derivatives. Checks that
  compare two routes through the package come second (`2ff2d7f`).
- Fewer tests come mainly from merging similar tests into tables, not from the
  grid.
- Planted bugs, run before and after every step, are the safety net
  (`benchmarks/local/mutation/`).
- Real blrm fits run only with `MOSTR_RUN_SLOW=true`.

**Why.** Measured on 15 planted bugs: checks that compare two routes caught 5;
adding independent checks raised that to 12. A comparison of two routes cannot
see a bug that affects both equally. The readable loop also lets a maintainer
check the SOP recursion without reading the C++. The new tests caught all 20
planted bugs, against 15 before.

**Decided against.** A grid of route comparisons only. Line coverage as the
measure of protection.

**Not true, so not to be tested:** percentile intervals need not contain the
estimate; the unconditional standard error is not always at least the
conditional one; MVN standard deviations do not converge to the delta standard
error; posterior medians and limits need not sum to one; weighted-draw
identities do not hold for ratios.

## 2026-11-07: Reject random effects; support for them paused

**Decided.** `mostm()` rejects `cluster()` and `mix_re()` terms for orm and
vglm, and `mix_re()` and prefixed `cluster()` for blrm (PR #18). Support for
rms 8.2 random effects is paused until it becomes relevant.

**Why.** Without the check, random-intercept orm fits were accepted silently:
SOPs were computed at a random effect of zero, which overstated treatment
effects by 4.5 to 17 percentage points, standard errors were wrong, and FWB was
biased. In simulations (2026-10-20), a fixed-effects fit with
patient-clustered standard errors was about as efficient for SOP treatment
effects as the marginal random-intercept fit, with 1 to 2 points of bias, so it
stays the recommended approach. Still to check: a simulation where the starting
state is correlated with the random effect.

## 2026-11-07: bartisan as a fitting engine: explored, not built

**Decided.** Nothing implemented; tracked in issue #25.

**Why.** Reading bartisan 0.1.0's source (not run): complete-data fits are
possible with `bcf(family = ordinal(), propensity = FALSE)`; the default
propensity model would use `yprev`, which is measured after randomisation.
SOPs would need a new route from draws of the linear predictor and cut-points.
Missing states cannot be handled through `custom_family()`; they need sampling
of hidden paths in a fork of bartisan's C++, or an R loop that imputes paths and
refits.

## 2026-10-28: Missing-state fits: start fit and optimizers (paused)

This work is on `feature/mostm-missing-states`, not on main.

**Findings.** All optimizers reached the same maximum from a start without an
engine fit, so starting values do not matter. SQUAREM cut EM steps about
fourfold when 95% of the information was missing, but gave no gain (1.1 times,
within noise) for 10,000 patients, 10 states, 28 days and 10% missing, where EM
needs 2 steps and the orm start fit (about 6 seconds) dominates. A "bridged
pairs" start fit, which would remove the need for transitions with both states
observed, worked in an in-memory patch but is not implemented.
`rms::orm()` fits `Ocens(-Inf, b)` with one state too few (confirmed on rms
8.2-1); `mostm()` is unaffected because it reads `Ocens()` itself. Not reported
upstream yet: the maintainer wants to investigate it in more detail before
opening an issue on the rms GitHub repository.

**Decided against, for now.** SQUAREM, because it gave no measurable gain at
realistic amounts of missing data.

## 2026-10-27: Missing states are part of `mostm()`, not a separate engine

On `feature/mostm-missing-states`.

**Decided.**
- `mostm(missing_states = 1:7)` or an `Ocens()` response switches orm, vglm
  (including partial proportional odds) and blrm to fits on the observed data.
  The default, `NULL`, uses complete rows and warns about dropped rows.
- `Ocens()` overrides `missing_states` row by row.
- A missing state can never be absorbing, so there is no `absorbing_state`
  argument; if the concept is needed, call it `absorb` as in `sops()`.
- With nothing missing, results must equal the plain fit.

**Why.** One implementation per engine. A patient lost to follow-up should
stop having rows, not enter an absorbing state.

**Decided against.** Separate `orm_markov_ml()` and `blrm_markov_ml()`
functions.

## 2026-10-27: Seeds and `datadist` with `mostm()`

**Decided.** The blrm sampler's seed goes in `sampling.args` (`8ad4de3`). The
README and vignettes no longer set up a `datadist` (`156ac4e`).

**Why.** `rmsb::blrm()` passes other arguments only to its optimizer, so a
`seed` given to `mostm()` never reached Stan. mostr's fits, SOPs, inference
and plots never use a `datadist`; only rms's own summaries do.

## 2026-10-26: One fitting function, `mostm(engine = )`

**Decided.** `mostm()` replaces `orm_markov()`, `vglm_markov()` and
`blrm_markov()` (issue #11). Shared arguments are named and checked once;
engine arguments go through `...`. With `engine = "vglm"` and no `family`, the
proportional-odds model that orm fits is used. `weights`, `type` and
`cadjust` are errors for blrm. Bootstrap refits do not compute
patient-clustered standard errors. The version became 0.2.0, as planned for
when comparisons landed.

**Why.** One place for shared arguments and their checks. `rmsb::blrm()` has
no row weights, and Bayesian fits take their uncertainty from the posterior.
Only the original fit uses the clustered standard errors.

**Decided against.** Keeping the three wrappers.

## 2026-10-23: Wald intervals from draws are centred on the point estimate

**Decided.** Draw-based Wald intervals use the point estimate as their centre
(`348901e`).

**Why.** Some were centred on the mean of the draws, which differs from the
point estimate, so the estimate was not in the middle of its interval.

## 2026-10-20: Warn when `times` skips fitted visits

**Decided.** `sops()` and related functions warn when `times` leaves out fitted
visits, is not increasing, or (for stored profiles) does not start at the
first follow-up visit (`e405750`).

**Why.** The recursion takes one transition per value of `times`, so
`times = c(1, 7, 14)` with daily visits silently gave the probabilities after
1, 2 and 3 transitions, labelled 1, 7 and 14.

**Decided against.** An error instead of a warning; reason not recorded.

## 2026-10-19: Simulators become `sim_markov(preset = )`

**Decided.** `sim_markov()` replaces the four trial simulators. The `"actt1"`
preset is the ACTT-2 transition model with ACTT-1's reported starting states
and the death threshold raised from -22.0693 to -21.4071, which gives 18%
control-arm deaths by day 28 by exact calculation (`44a87fb`).

**Decided against.** Keeping the Rohde et al. ACTT-1 model; reason not
recorded.

## 2026-10-17: `by` follows marginaleffects

**Decided.**
- Every summary comes in a pair (`sops()`/`avg_sops()`,
  `state_time()`/`avg_state_time()`, `comparisons()`/`avg_comparisons()`)
  that reads `by` as marginaleffects does (`66566f6`, `035c099`).
- The ordinary bootstrap stops when a resample lacks every patient of a `by`
  group, and points to FWB (`5b62643`).

**Why.** The commits name matching marginaleffects, which inspired mostr's
interface, as the goal; nothing more is recorded. A missing group had no
estimate in that resample, so its interval silently rested on fewer
resamples.

**Decided against.** Dropping such resamples. The ordinary bootstrap for
`comparisons()` per profile was not offered. The maintainer's reason
(2026-11-13): an average over patients still works when a resample lacks a
rare state, because that state then simply has probability 0 and the other
patients are still predicted. A single profile that starts in the missing
state has no prediction at all in that resample, and nothing can stand in for
it, so the per-profile interval is not defined. The maintainer wants to look
at this in more detail. Checked 2026-11-13: per-profile `comparisons()`
refuses `method = "boot"` and points to FWB, which ran for a profile starting
in a state that 5 of 200 patients ever had.

## 2026-10-12: ACTT-2 model refitted with states 1 and 2 merged

**Decided.** The ACTT-2 simulation model was refitted with states 1 and 2
merged (`cf18c8e`, version 0.1.1).

**Why.** State 2 was mostly recorded only at study visits. The old fit
rehospitalized about 29% of control patients (2% observed); the new one about
2.7%.

## 2026-10-01: Waning-effect vignette kept out of the repository

**Decided.** `vignettes/waning-treatment-effect-sops.Rmd` is local only, "for
now" (`5f46e20`), because it is not ready.

## 2026-09-23: Spearman correlation by default in diagnostics

**Decided.** `plot_correlation()` and `plot_variogram()` use Spearman rank
correlation by default (`67eef4b`).

**Why.** It accounts for tied states. `method = "pearson"` gives the earlier
calculation.
