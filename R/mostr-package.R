#' mostr: Markov Ordinal State Transition Models
#'
#' @description
#' mostr analyses outcomes recorded repeatedly on an ordered scale, such as a
#' patient's daily health state from "at home" to "dead". It fits a
#' *transition model*, an ordinal regression of each day's state on the
#' previous day's state, time, treatment and other covariates, and turns the
#' fitted model into quantities that are easy to report:
#'
#' * **state occupancy probabilities (SOPs)**: the probability of being in each
#'   state on each day;
#' * **expected time in a state**: for example, the expected number of days at
#'   home during follow-up;
#' * **treatment comparisons** of either, as differences or ratios, with
#'   confidence intervals.
#'
#' The sections below list every function by what you want to do. For a worked
#' example, start with `vignette("mostr")`.
#'
#' @section Prepare data and fit a model:
#' The data have one row per patient and visit, with the patient ID (`id`),
#' the visit time (`time`), the current state (`y`), the state at the previous
#' visit (`yprev`) and covariates such as treatment (`tx`).
#'
#' * [prepare_markov_data()] removes the rows after a patient has entered an
#'   absorbing state such as death, and converts `y` and `yprev` to factors.
#' * [mostm()] fits the transition model. `engine = "orm"` (the default) uses
#'   [rms::orm()], `engine = "vglm"` uses [VGAM::vglm()], which also fits
#'   partial proportional-odds models, and `engine = "blrm"` fits a Bayesian
#'   model with `rmsb::blrm()`. The functions for probabilities, times,
#'   comparisons, confidence intervals and model checks need a fit made by
#'   [mostm()], because it stores the patient data they use.
#'
#' @section Estimate state occupancy probabilities:
#' All of these take `times`, the consecutive visits to predict (for example
#' `times = 1:28`), and `absorb`, the absorbing states such as death.
#'
#' * [avg_sops()] averages the SOPs over the patients. With
#'   `variables = list(tx = c(0, 1))` it gives the SOPs if every patient had
#'   received control, and if every patient had received treatment.
#' * [sops()] gives the SOPs of each patient, or of each row of `newdata`.
#' * [interpolate_sops()] moves SOPs estimated at scheduled visits, such as
#'   days 3, 7 and 14, onto a scale of elapsed days.
#' * [soprob_markov()] is the low-level calculation behind [sops()]: it
#'   returns the SOPs as an array with one dimension each for patient, time
#'   and state.
#'
#' @section Estimate time spent in states:
#' * [avg_state_time()] gives the expected total time in each state, or in a
#'   set of states given by `state_sets`, averaged over patients.
#' * [state_time()] gives the same for each patient.
#' * [time_in_state()] adds up the time in chosen states from results you have
#'   already calculated with [sops()], [avg_sops()] or [soprob_markov()].
#'
#' @section Compare treatments:
#' * [avg_comparisons()] compares two or more values of one variable, usually
#'   treatment, averaged over patients. `estimand = "sop"` compares SOPs on each
#'   day, `estimand = "time_in_state"` compares the expected time in states,
#'   and `estimand = "time_benefit"` estimates how much longer patients spend in
#'   better states under one value than under the other. `comparison` chooses a
#'   difference or a ratio.
#' * [comparisons()] gives the same comparisons for each patient.
#'
#' @section Add confidence intervals:
#' [inferences()] adds standard errors, confidence intervals and, for
#' comparisons, tests to any result above. Its `method` argument chooses how:
#'
#' * `method = "delta"` calculates the intervals from formulas, without random
#'   draws, so it is fast and gives the same result every time. It works for
#'   first-order, full proportional-odds models fitted with the `"orm"` or
#'   `"vglm"` engine. Its `vcov` argument chooses the variance. **Conditional
#'   variance** (`vcov = "conditional"`) accounts for the uncertainty in the
#'   estimated model coefficients and treats the patients' starting states and
#'   covariates as given. **Unconditional variance** (`vcov = "unconditional"`,
#'   the default for averages) also accounts for which patients happened to
#'   be sampled and for their role in estimating the coefficients. See
#'   `vignette("analytical-confidence-intervals")`.
#' * `method = "mvn"` (the default) draws coefficient vectors from a
#'   multivariate normal distribution and recalculates the result for each
#'   draw. It works for every frequentist model, including partial
#'   proportional-odds and second-order models.
#' * `method = "bootstrap"` refits the model to patients resampled with
#'   replacement, and `method = "fwb"` (fractional weighted bootstrap) refits
#'   it with a random positive weight for every patient. See
#'   `vignette("mvn-vs-bootstrap-orm-sops")` and
#'   `vignette("fractional-weighted-bootstrap-sops")`.
#'
#' Results from `engine = "blrm"` already contain posterior intervals and need
#' no [inferences()] step; see `vignette("bayesian-rmsb-sops")`.
#'
#' After [inferences()]:
#'
#' * [get_draws()] returns the individual simulation or bootstrap draws, for
#'   summaries of your own.
#' * [vcov()][vcov.markov_sops] returns the covariance matrix of the estimates
#'   after `method = "delta"`. Its `rows` argument selects the estimates, which
#'   saves memory; the option [mostr.delta_max_bytes] caps the memory that
#'   these calculations may use.
#'
#' @section Plot results:
#' * [plot_sops()] plots SOPs over time, as lines or stacked bars, from
#'   observed trajectories or from [sops()] or [avg_sops()] results.
#' * [plot_comparisons()] plots [avg_comparisons()] results, or the
#'   [comparisons()] of one patient.
#' * [plot_time()] plots [avg_state_time()] results as stacked bars of the
#'   time spent in each state, or as points with confidence intervals.
#' * [plot_lp_difference()] shows, for one patient profile, how the variable
#'   in `variables` (usually treatment) shifts the model's linear predictor
#'   (the log odds of a worse state) at each time, separately for each
#'   previous state.
#'
#' @section Check the model against the data:
#' Each of these functions accepts either the observed trajectories or a fit
#' from [mostm()], so you can draw the same summary from both and compare them.
#'
#' * [plot_transitions()] shows the proportion of patients moving from each
#'   state to each other state on selected days.
#' * [plot_correlation()] shows the rank correlation between patients' states
#'   on every pair of days.
#' * [plot_variogram()] plots these correlations against the number of days
#'   between the two days.
#'
#' @section Simulate data:
#' * [sim_markov()] simulates daily ordinal trajectories from a transition
#'   model, either from a model of a published trial (the presets
#'   [actt1_preset], [actt2_preset] and [violet_preset]) or from a model you
#'   describe.
#' * [violet_baseline] holds the starting states, ages and SOFA scores of the
#'   VIOLET trial patients, from which `sim_markov("violet")` samples.
#'
#' @section Work with model coefficients and standard errors:
#' These functions are used inside [inferences()] and are exported for users
#' who want to check or extend its calculations.
#'
#' * [bootstrap_model_coefs()] refits the model to resampled patients and
#'   returns the coefficients of each refit.
#' * [get_vcov_robust()] returns a covariance matrix of the coefficients whose
#'   standard errors allow for correlation between the transitions of one
#'   patient.
#' * [robcov_vglm()] computes such a covariance matrix for a [VGAM::vglm()]
#'   fit; `engine = "vglm"` returns its result.
#' * [compare_se_orm_vglm()] lists the standard errors of an [rms::orm()] fit
#'   and a [VGAM::vglm()] fit of the same model side by side.
#' * [set_coef()] replaces the coefficients of a fitted model.
#' * [get_effective_coefs()] returns the coefficients as a matrix with one row
#'   per threshold between states.
#'
#' @useDynLib mostr, .registration = TRUE
"_PACKAGE"
