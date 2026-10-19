#' Simulate ordinal patient trajectories from a Markov model
#'
#' @description
#' Generates daily ordinal states for each patient from a first-order Markov
#' model: the state on day `t` depends on the state on day `t - 1`, the
#' patient's covariates, and the day, through a proportional-odds or partial
#' proportional-odds model. Choose a ready-made model of a published trial with
#' `preset`, or describe your own model with `preset = NULL`.
#'
#' The presets are:
#'
#' - `"actt1"`: 8 states of the NIAID ordinal scale for a population like that
#'   of the ACTT-1 trial of Remdesivir in COVID-19, built from the ACTT-2 model
#'   (see [actt1_preset]).
#' - `"actt2"`: 8 states of the NIAID ordinal scale, from a model of the ACTT-2
#'   trial of Baricitinib plus Remdesivir in COVID-19 (see [actt2_preset]).
#' - `"violet"`: 4 states (home, hospital, ventilator or ARDS, dead), from
#'   Frank Harrell's model of the VIOLET trial of vitamin D3 in critically ill
#'   adults (see [violet_preset]).
#'
#' With a preset, `n_patients` sets the sample size and `treatment_effect` the
#' treatment effect. Any model argument you supply, such as `intercepts` or
#' `baseline_data`, replaces that part of the preset and keeps the rest.
#'
#' @param preset The model to simulate from: `"actt1"`, `"actt2"`, or
#'   `"violet"`, or `NULL` to supply your own model through `baseline_data`,
#'   `intercepts`, and `lp_function`. There is no default.
#' @param n_patients Number of patients to sample from the preset's population.
#'   Must be a positive integer (default: 1000). Cannot be combined with
#'   `baseline_data`.
#' @param treatment_prob Probability that a sampled patient is assigned to
#'   treatment (`tx = 1`); the other patients are controls (`tx = 0`). Must be
#'   between 0 and 1 (default: 0.5). Cannot be combined with `baseline_data`.
#' @param treatment_effect Log odds ratio of being in a worse state on day 1
#'   for treated (`tx = 1`) versus control (`tx = 0`) patients. Values below 0
#'   mean that treatment helps. `NULL` (the default) uses the preset's value:
#'   -0.0730 for `"actt1"` and `"actt2"`, and `log(0.75)` for `"violet"`. With
#'   `preset = NULL` it means 0, no effect.
#' @param treatment_effect_decay Nonnegative exponential decay rate of the
#'   treatment effect per day. The effect on day `t` is `treatment_effect *
#'   exp(-treatment_effect_decay * (t - 1))`, so `-log(0.01) / 14` shrinks it to
#'   1% of its day-1 size by day 15. The default 0 gives a constant effect.
#' @param follow_up_time Number of days to simulate (default: 28).
#' @param baseline_data A data frame with one row per patient, replacing the
#'   patients that a preset would sample. It must contain `id`, `yprev` (the
#'   state on the day before the first simulated day), and every covariate
#'   that the model uses: `tx` for `"actt1"` and `"actt2"`, and `age`, `sofa`,
#'   and `tx` for `"violet"`. Required when `preset = NULL`.
#' @param intercepts Numeric vector of thresholds, one fewer than the number of
#'   states, in increasing order. With `K` states numbered from 1 (best) to `K`
#'   (worst), `intercepts[j]` plus the linear predictor is the log odds of
#'   being in state `K + 1 - j` or worse: `intercepts[1]` belongs to the worst
#'   state alone and `intercepts[K - 1]` to states 2 and worse. Required when
#'   `preset = NULL`; otherwise `NULL` uses the preset's thresholds.
#' @param lp_function A function that calculates the linear predictor, which
#'   is added to the thresholds. It is called with these named arguments:
#'   - `yprev`: previous state, a factor with levels `1` to `K`
#'   - `t`: current day
#'   - one argument for each column named in `covariate_names`, and `tx` when
#'     `baseline_data` has a `tx` column
#'   - `treatment_effect`: the treatment effect on day `t`, already reduced by
#'     `treatment_effect_decay`
#'   - `extra_params`: the value of `extra_params`
#'
#'   Each argument holds one patient's value. The function must return either
#'   one number, added to every threshold (a proportional-odds model), or a
#'   numeric vector with `length(intercepts)` elements, added to the
#'   corresponding thresholds (a partial proportional-odds model). Larger
#'   values make worse states more likely. Required when `preset = NULL`.
#' @param extra_params Optional object, usually a named numeric vector of
#'   coefficients, passed unchanged to `lp_function`. `NULL` uses the preset's
#'   coefficients.
#' @param absorbing_states Integer vector of states that patients cannot leave
#'   once they enter them. `NULL` (the default) means the worst state,
#'   `length(intercepts) + 1`, which is death in every preset. Use `integer()`
#'   when no state is absorbing.
#' @param covariate_names Character vector of `baseline_data` columns passed to
#'   `lp_function` as arguments with the same names. A `tx` column is always
#'   passed. `NULL` (the default) uses the preset's covariates, or, with
#'   `preset = NULL`, every column except `id` and `yprev`.
#' @param seed Integer. Random seed for reproducibility (default: `NULL`, no
#'   seed set). The seed is set once, before patients and treatment are
#'   sampled, and the trajectories continue from the same random-number
#'   stream.
#'
#' @return A data frame with one row per patient and day and columns:
#'   - `id`: patient identifier
#'   - `time`: day, from 1 to `follow_up_time`
#'   - `y`: state on this day
#'   - `yprev`: state on the previous day
#'   - all other columns of `baseline_data`, such as `tx`
#'
#' @details
#' On each day the function:
#' 1. calculates the linear predictor for each patient not in an absorbing
#'    state,
#' 2. adds it to the thresholds to get the probability of each state or worse,
#' 3. converts these to the probability of each state,
#' 4. samples the next state for each patient, and
#' 5. keeps patients in absorbing states where they are.
#'
#' States must be consecutive integers from 1 (best) to `K` (worst). If
#' `lp_function` returns threshold-specific linear predictors, the probability
#' of being in a state or worse must not increase from one state to the next
#' worse state; otherwise the function stops because the implied state
#' probabilities would be negative.
#'
#' Each preset's thresholds follow the order described under `intercepts`. To
#' change one of them, start from the preset's help page, which lists them.
#'
#' @examples
#' # 50 patients from the ACTT-2 model
#' actt2 <- sim_markov("actt2", n_patients = 50, seed = 123)
#' head(actt2)
#'
#' # A treatment effect that fades to 1% of its day-1 size by day 15
#' decaying <- sim_markov(
#'   "actt2",
#'   n_patients = 50,
#'   treatment_effect = log(0.8),
#'   treatment_effect_decay = -log(0.01) / 14,
#'   seed = 123
#' )
#'
#' # Your own patients with the ACTT-1 model
#' patients <- data.frame(id = 1:4, yprev = c(5, 5, 7, 7), tx = c(0, 1, 0, 1))
#' sim_markov("actt1", baseline_data = patients, follow_up_time = 7, seed = 1)
#'
#' # Your own model with three states: 1 = home, 2 = in hospital, 3 = dead
#' baseline <- data.frame(id = 1:20, yprev = 2, tx = rep(0:1, 10))
#' hospital_lp <- function(yprev, t, tx, treatment_effect, extra_params) {
#'   # Hospitalized patients are more likely to stay in hospital or die
#'   5 * (yprev == "2") - 0.1 * t + treatment_effect * tx
#' }
#' trajectories <- sim_markov(
#'   NULL,
#'   baseline_data = baseline,
#'   follow_up_time = 10,
#'   intercepts = c(-9, -3),
#'   lp_function = hospital_lp,
#'   treatment_effect = log(0.8),
#'   seed = 1
#' )
#' head(trajectories)
#'
#' # Partial proportional odds: treatment affects only the odds of death
#' death_only_lp <- function(yprev, t, tx, treatment_effect, extra_params) {
#'   hospital_lp(yprev, t, tx = 0, treatment_effect, extra_params) +
#'     c(treatment_effect * tx, 0)
#' }
#' trajectories_ppo <- sim_markov(
#'   NULL,
#'   baseline_data = baseline,
#'   follow_up_time = 10,
#'   intercepts = c(-9, -3),
#'   lp_function = death_only_lp,
#'   treatment_effect = log(0.8),
#'   seed = 1
#' )
#'
#' @seealso [actt1_preset], [actt2_preset], and [violet_preset] for the
#'   presets' models.
#'
#' @importFrom stats plogis
#'
#' @export
sim_markov <- function(
  preset,
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = NULL,
  treatment_effect_decay = 0,
  follow_up_time = 28,
  baseline_data = NULL,
  intercepts = NULL,
  lp_function = NULL,
  extra_params = NULL,
  absorbing_states = NULL,
  covariate_names = NULL,
  seed = NULL
) {
  if (missing(preset)) {
    stop(
      "`preset` is required: use one of ",
      paste0('"', markov_preset_names, '"', collapse = ", "),
      ", or NULL to supply your own model.",
      call. = FALSE
    )
  }
  model <- markov_preset(preset)

  sampling_supplied <- !missing(n_patients) || !missing(treatment_prob)
  if (!is.null(baseline_data) && sampling_supplied) {
    stop(
      "Supply either `baseline_data` or `n_patients` and `treatment_prob`, ",
      "not both: `baseline_data` already fixes the patients and their ",
      "treatment.",
      call. = FALSE
    )
  }
  if (is.null(model)) {
    if (sampling_supplied) {
      stop(
        "`n_patients` and `treatment_prob` sample patients from a preset. ",
        "With `preset = NULL`, supply the patients in `baseline_data`.",
        call. = FALSE
      )
    }
    missing_model <- c("baseline_data", "intercepts", "lp_function")[
      c(is.null(baseline_data), is.null(intercepts), is.null(lp_function))
    ]
    if (length(missing_model) > 0) {
      stop(
        "With `preset = NULL`, supply ",
        paste0("`", missing_model, "`", collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  } else if (is.null(baseline_data)) {
    validate_markov_wrapper_inputs(n_patients, treatment_prob)
  }

  treatment_effect <- treatment_effect %||% model$treatment_effect %||% 0
  validate_markov_treatment_effect(treatment_effect, treatment_effect_decay)
  validate_markov_follow_up_time(follow_up_time)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(baseline_data)) {
    baseline_data <- model$sample_patients(
      as.integer(n_patients),
      treatment_prob
    )
  }
  if (is.null(covariate_names)) {
    covariate_names <- if (is.null(model)) {
      setdiff(names(baseline_data), c("id", "yprev"))
    } else {
      model$covariate_names
    }
  }
  intercepts <- intercepts %||% model$intercepts

  simulate_markov_trajectories(
    baseline_data = baseline_data,
    follow_up_time = follow_up_time,
    intercepts = intercepts,
    lp_function = lp_function %||% model$lp_function,
    extra_params = extra_params %||% model$extra_params,
    treatment_effect = treatment_effect,
    treatment_effect_decay = treatment_effect_decay,
    absorbing_states = absorbing_states %||% (length(intercepts) + 1L),
    covariate_names = covariate_names
  )
}


#' ACTT-2 preset for sim_markov()
#'
#' @description
#' `sim_markov("actt2")` generates synthetic daily health states from a
#' first-order Markov proportional-odds model based on data from the
#' [Adaptive COVID-19 Treatment Trial 2 (ACTT-2)](https://doi.org/10.1056/NEJMoa2031994),
#' which compared Baricitinib plus Remdesivir with placebo plus Remdesivir in
#' adults hospitalized with COVID-19 (Kalil et al., 2021). States follow the
#' trial's
#' [NIAID eight-point ordinal scale](https://www.fda.gov/media/144473/download#page=10):
#'
#' 1. outside hospital, with usual activities unrestricted
#' 2. outside hospital, with restricted activities, oxygen at home, or both
#' 3. in hospital, without oxygen and no further medical care needed
#' 4. in hospital, without oxygen but still needing medical care
#' 5. in hospital, needing supplemental oxygen
#' 6. in hospital, needing high-flow oxygen or noninvasive ventilation
#' 7. in hospital, needing invasive ventilation or extracorporeal membrane
#'    oxygenation (ECMO)
#' 8. dead
#'
#' The output has 28 rows per patient by default, and dead patients stay in
#' state 8.
#'
#' @details
#' **Patients.** Each patient's state on the day before follow-up (`yprev`) is
#' drawn from states 4 to 7 with probabilities 0.1374637, 0.5459826,
#' 0.2090997, and 0.1074540. Treatment is assigned independently with
#' probability `treatment_prob`: `tx = 1` stands for Baricitinib plus
#' Remdesivir and `tx = 0` for placebo plus Remdesivir.
#'
#' **Transition model.** Write \eqn{Y_t} for the state on day \eqn{t} and
#' \eqn{j} for the state on the previous day. The log odds of being in state
#' \eqn{k} or worse on day \eqn{t} is
#' \deqn{\mathrm{logit}\,\Pr(Y_t \ge k \mid Y_{t-1} = j) = \alpha_k +
#'   \gamma_j + f(t) + \beta(t)\,\mathrm{tx}}{logit Pr(Y_t >= k | Y_{t-1} = j)
#'   = alpha_k + gamma_j + f(t) + beta(t) tx}
#' for the cutoffs \eqn{k = 2, \ldots, 8}. Here \eqn{\beta(t)} is the
#' treatment effect on day \eqn{t} described under `treatment_effect_decay` in
#' [sim_markov()]. Every term shifts all cutoffs equally, which makes this a
#' proportional-odds model. The thresholds are
#'
#' | Cutoff \eqn{k} | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
#' |---|---|---|---|---|---|---|---|
#' | \eqn{\alpha_k} | -0.6465 | -4.9314 | -5.2230 | -7.4169 | -11.5321 | -15.2262 | -22.0693 |
#'
#' and the previous-state effects are
#'
#' | Previous state \eqn{j} | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
#' |---|---|---|---|---|---|---|---|
#' | \eqn{\gamma_j} | 0 | 0 | 6.2218 | 7.1412 | 9.9809 | 14.1645 | 19.3517 |
#'
#' The large previous-state effects mean that patients usually stay in, or
#' move close to, their previous state. The day enters through a restricted
#' cubic spline, \eqn{f(t) = -0.3084\,t + 2.5186\,t' - 5.0276\,t'' +
#' 3.1714\,t''' - 1.0672\,t''''}{f(t) = -0.3084 t + 2.5186 t' - 5.0276 t'' +
#' 3.1714 t''' - 1.0672 t''''}, where \eqn{t'} to \eqn{t''''} are the
#' nonlinear basis columns that `rms::rcs()` creates with knots at days 2, 6,
#' 11, 16, 21, and 27. Days beyond 28, for example with
#' `follow_up_time = 60`, extrapolate this spline.
#'
#' In ACTT-2, state 2 was mostly recorded only at the day 14, 21, and 28
#' visits, so the transition model was fitted with states 1 and 2 merged. The
#' added lowest threshold, \eqn{\alpha_2}, splits them again, placing the
#' observed 9.3% of the control arm in state 2 on day 14. Both states share
#' the same coefficients, so the split does not change the occupancy of
#' states 3 to 8.
#'
#' **Treatment effect.** The default `treatment_effect` of -0.0730 is the
#' coefficient fitted to the ACTT-2 data.
#'
#' **Changing the model.** [sim_markov()] expects the thresholds from the
#' worst state down, so `intercepts = rev(alpha)` with `alpha` the
#' \eqn{\alpha_k} above in the order of the table reproduces the preset. Change
#' an element of `alpha` before reversing it to move one threshold.
#'
#' @references
#' Kalil AC, Patterson TF, Mehta AK, et al. (2021). "Baricitinib plus
#' Remdesivir for Hospitalized Adults with Covid-19." *New England Journal of
#' Medicine*, 384(9), 795-807. \doi{10.1056/NEJMoa2031994}
#'
#' @examples
#' actt2 <- sim_markov("actt2", n_patients = 20, seed = 123)
#' head(actt2)
#'
#' # The same model with a higher threshold for death
#' alpha <- c(-0.6465, -4.9314, -5.2230, -7.4169, -11.5321, -15.2262, -21.5)
#' sicker <- sim_markov(
#'   "actt2",
#'   n_patients = 20,
#'   intercepts = rev(alpha),
#'   seed = 123
#' )
#'
#' @seealso [sim_markov()] to simulate, and [actt1_preset] for the same model
#'   in a sicker population.
#'
#' @name actt2_preset
NULL


#' ACTT-1 preset for sim_markov()
#'
#' @description
#' `sim_markov("actt1")` generates synthetic daily health states for a
#' population like that of the
#' [Adaptive COVID-19 Treatment Trial 1 (ACTT-1)](https://doi.org/10.1056/NEJMoa2007764),
#' which compared Remdesivir with placebo in adults hospitalized with COVID-19
#' (Beigel et al., 2020). It uses the ACTT-2 transition model described in
#' [actt2_preset], with the same eight states, and changes two things:
#'
#' - the patients' starting states follow the distribution reported for
#'   ACTT-1, which had more patients on ventilation than ACTT-2; and
#' - the threshold for death is raised so that 18% of control patients have
#'   died by day 28.
#'
#' The output has 28 rows per patient by default, and dead patients stay in
#' state 8.
#'
#' @details
#' **Patients.** Each patient's state on the day before follow-up (`yprev`) is
#' drawn from states 4 to 7 in proportion to the 138, 435, 193, and 285
#' ACTT-1 patients reported in these states at enrollment (13.1%, 41.4%,
#' 18.4%, and 27.1% of the 1,051 patients with a reported state). In ACTT-2
#' these proportions were 13.7%, 54.6%, 20.9%, and 10.7%. Treatment is
#' assigned independently with probability `treatment_prob`: `tx = 1` for
#' treated and `tx = 0` for control patients.
#'
#' **Death threshold.** With these starting states and the unchanged ACTT-2
#' model, 11.1% of control patients would die by day 28. The threshold for
#' death, \eqn{\alpha_8}, is therefore raised from -22.0693 to -21.4071. On
#' each day, this multiplies every patient's odds of dying by
#' \eqn{\exp(0.6622) = 1.94} and leaves their probability of being in state 7
#' or worse unchanged, so the additional deaths come from patients who would
#' otherwise have needed invasive ventilation (state 7). The value was found by calculating the
#' control arm's state probabilities exactly, multiplying the 28 daily
#' transition matrices, and solving for 18% in state 8 on day 28. All other
#' thresholds and coefficients are those of [actt2_preset].
#'
#' The resulting control-arm state probabilities, in percent, are
#'
#' | Day | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
#' |---|---|---|---|---|---|---|---|---|
#' | 7 | 18.81 | 7.22 | 1.11 | 10.80 | 19.76 | 10.83 | 22.79 | 8.68 |
#' | 14 | 42.78 | 8.03 | 0.53 | 5.26 | 10.58 | 5.93 | 14.53 | 12.36 |
#' | 28 | 55.70 | 9.79 | 0.21 | 2.02 | 4.44 | 2.68 | 7.17 | 18.00 |
#'
#' **Treatment effect.** The default `treatment_effect` of -0.0730 is the
#' coefficient fitted to the ACTT-2 data for Baricitinib, not an estimate from
#' ACTT-1. Choose the effect you want to study with `treatment_effect`; the
#' calibration to 18% applies to the control arm only.
#'
#' **Changing the model.** The thresholds, from the best cutoff to death, are
#' -0.6465, -4.9314, -5.2230, -7.4169, -11.5321, -15.2262, and -21.4071.
#' Supply them reversed as `intercepts` in [sim_markov()] after changing one,
#' as shown in [actt2_preset].
#'
#' @references
#' Beigel JH, Tomashek KM, Dodd LE, et al. (2020). "Remdesivir for the
#' Treatment of Covid-19 - Final Report." *New England Journal of Medicine*,
#' 383(19), 1813-1826. \doi{10.1056/NEJMoa2007764}
#'
#' @examples
#' actt1 <- sim_markov("actt1", n_patients = 20, seed = 123)
#' head(actt1)
#'
#' @seealso [sim_markov()] to simulate, and [actt2_preset] for the transition
#'   model.
#'
#' @name actt1_preset
NULL


#' VIOLET preset for sim_markov()
#'
#' @description
#' `sim_markov("violet")` simulates daily 4-state trajectories resembling the
#' [VIOLET trial](https://doi.org/10.1056/NEJMoa1911124) of early high-dose
#' vitamin D3 versus placebo in critically ill adults with vitamin D deficiency.
#' The transition model is the partial proportional-odds Markov model that
#' Frank Harrell fitted to the VIOLET data to simulate his large `simlongord`
#' trial (Harrell, 2022), and patients are drawn from the same 1,351 VIOLET
#' patients. Its daily state proportions agree with `simlongord` within
#' simulation error. The states are:
#'
#' 1. at home
#' 2. in hospital or another health care facility
#' 3. on a ventilator or with acute respiratory distress syndrome (ARDS)
#' 4. dead
#'
#' The output has columns `id`, `time`, `y`, `yprev`, `age`, `sofa`, and `tx`,
#' 28 rows per patient by default, and dead patients stay in state 4.
#'
#' @details
#' **Patients.** Each simulated patient's starting state (`yprev`), age, and
#' SOFA (Sequential Organ Failure Assessment) score are drawn together, with
#' replacement, from the 1,351 VIOLET patients summarized in
#' [violet_baseline]. About 65% of patients start in hospital (state 2) and
#' 35% on a ventilator or with ARDS (state 3). Treatment (`tx = 1`) is
#' assigned independently of these characteristics with probability
#' `treatment_prob`; control patients have `tx = 0`. A `baseline_data` that
#' you supply must contain the columns `id`, `yprev`, `age`, `sofa`, and `tx`.
#'
#' **Transition model.** Write \eqn{Y_t} for the state on day \eqn{t} and
#' \eqn{k} for the state on the previous day. The log odds of being in state
#' \eqn{j} or worse on day \eqn{t} is
#' \deqn{\mathrm{logit}\,\Pr(Y_t \ge j \mid Y_{t-1} = k) = \alpha_j +
#'   \gamma_k + \delta_j t + \delta'_j \max(t - 2, 0) +
#'   0.0105\,\mathrm{age} + 0.046\,\mathrm{SOFA} +
#'   \beta(t)\,\mathrm{tx}}{logit Pr(Y_t >= j | Y_{t-1} = k) = alpha_j +
#'   gamma_k + delta_j t + delta'_j max(t - 2, 0) + 0.0105 age + 0.046 SOFA +
#'   beta(t) tx}
#' for the cutoffs \eqn{j = 2, 3, 4}. Here \eqn{\beta(t)} is the treatment
#' effect on day \eqn{t} described under `treatment_effect_decay` in
#' [sim_markov()]. Positive values make worse states more likely: for example,
#' each additional SOFA point multiplies the odds of being in a worse state by
#' \eqn{\exp(0.046) = 1.05}.
#'
#' The thresholds \eqn{\alpha_j} and time slopes \eqn{\delta_j} and
#' \eqn{\delta'_j} differ between the cutoffs, so the mix of states can change
#' over follow-up (a partial proportional odds model for time):
#'
#' | Cutoff \eqn{j} | States counted | \eqn{\alpha_j} | \eqn{\delta_j} | \eqn{\delta'_j} |
#' |---|---|---|---|---|
#' | 4 | dead | -4.713 | 0.412 | -0.410 |
#' | 3 | ventilator/ARDS or dead | 1.741 | -0.481 | 0.525 |
#' | 2 | hospital or worse | 10.034 | -1.315 | 1.384 |
#'
#' Time enters as a straight line whose slope changes after day 2, because
#' very few patients went home on day 1. The previous state shifts all three
#' cutoffs equally: \eqn{\gamma_k} is -15.797 after a day at home (\eqn{k = 1}),
#' -6.364 after a day in hospital (\eqn{k = 2}), and 0 after a day on a
#' ventilator or with ARDS (\eqn{k = 3}). The large negative values mean that
#' patients usually stay in, or improve from, their previous state.
#'
#' **Treatment effect.** The default `treatment_effect` of `log(0.75)` is the
#' odds ratio of 0.75 that Harrell built into `simlongord`, the same for every
#' previous state. In VIOLET itself, vitamin D had no detectable effect
#' (fitted log odds ratio -0.028, standard error 0.050).
#'
#' **Relation to Harrell's model.** Harrell numbered the states from dead
#' (lowest) to home (highest), modeled the probability of a better state, and
#' published the coefficients rounded to three decimals. This preset numbers
#' the states from best to worst, so each intercept and coefficient above is
#' his printed value with the opposite sign; the transition probabilities are
#' the same. The one exception is age. Harrell printed 0.011 (with his sign,
#' -0.011), but `simlongord` was evidently simulated with a value close to
#' 0.0105: refitting his model to all 250,000 simulated patients gives 0.0105
#' (standard error 0.00015), while every other coefficient agrees with its
#' printed value. A difference of 0.0005 per year of age moves a typical
#' patient's log odds by about 0.03 every day, so with 0.011 about 1 percentage
#' point fewer patients are at home by day 28 than in `simlongord`. With 0.0105
#' the daily state proportions agree with `simlongord` within simulation
#' error.
#'
#' He fitted the model to days 1 to 27, because ventilator use on day 28
#' looked anomalous, and simulated 28 days, which is the default here. Later
#' days extend the fitted time trends. Harrell assigned treatment once to each
#' of the 1,351 VIOLET patients before resampling them, so in `simlongord` the
#' arms differ slightly in their starting states and covariates (for example,
#' 33.9% versus 36.0% start on a ventilator). This preset assigns treatment
#' independently to every simulated patient, so the arms differ only by chance.
#'
#' @references
#' National Heart, Lung, and Blood Institute PETAL Clinical Trials Network
#' (2019). "Early High-Dose Vitamin D3 for Critically Ill, Vitamin
#' D-Deficient Patients." *New England Journal of Medicine*, 381(26),
#' 2529-2540. \doi{10.1056/NEJMoa1911124}
#'
#' Harrell FE Jr (2022). "Simulating a Large Clinical Trial From VIOLET."
#' <https://hbiostat.org/data/repo/simlongord>
#'
#' @examples
#' violet <- sim_markov("violet", n_patients = 20, seed = 123)
#' head(violet)
#'
#' # A trial in which treatment has no effect
#' null_trial <- sim_markov(
#'   "violet",
#'   n_patients = 20,
#'   treatment_effect = 0,
#'   seed = 123
#' )
#'
#' @seealso [sim_markov()] to simulate, and [violet_baseline] for the VIOLET
#'   patients.
#'
#' @name violet_preset
NULL


simulate_markov_trajectories <- function(
  baseline_data,
  follow_up_time,
  intercepts,
  lp_function,
  extra_params,
  treatment_effect,
  treatment_effect_decay,
  absorbing_states,
  covariate_names
) {
  # Input validation
  if (!is.data.frame(baseline_data)) {
    stop("baseline_data must be a data frame")
  }

  required_cols <- c("id", "yprev")
  missing_cols <- setdiff(required_cols, names(baseline_data))
  if (length(missing_cols) > 0) {
    stop(
      "baseline_data must contain columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.numeric(intercepts) || length(intercepts) < 1) {
    stop("intercepts must be a numeric vector with at least one element")
  }

  if (!is.function(lp_function)) {
    stop("lp_function must be a function")
  }

  # Setup
  N <- nrow(baseline_data)
  times <- 1:follow_up_time
  n_states <- length(intercepts) + 1

  # Initialize state matrix
  # Rows = patients, Columns = time points (0 to follow_up_time)
  state_matrix <- matrix(
    as.numeric(baseline_data$yprev),
    nrow = N,
    ncol = length(times) + 1,
    dimnames = list(baseline_data$id, 0:max(times))
  )

  # Detect which covariates from baseline_data need to be passed to lp_function
  available_covariates <- intersect(covariate_names, names(baseline_data))
  if (length(covariate_names) > 0 && length(available_covariates) == 0) {
    warning("None of the specified covariate_names found in baseline_data")
  }

  # Main simulation loop
  for (t in times) {
    # Current states
    yprev <- state_matrix[, t]

    # Identify patients not in absorbing states
    active_idx <- which(!yprev %in% absorbing_states)

    # If all patients are in absorbing states, fill remaining time and exit
    if (length(active_idx) == 0) {
      if (t <= max(times)) {
        for (i in (t + 1):ncol(state_matrix)) {
          state_matrix[, i] <- state_matrix[, t]
        }
      }
      break
    }

    # Extract data for active patients
    yprev_active <- yprev[active_idx]
    X_active <- baseline_data[active_idx, , drop = FALSE]

    # Build argument list for lp_function
    # This allows dynamic passing of covariates
    # We pass yprev as a factor to ensure compatibility with lp functions that expect it
    lp_args <- list(
      yprev = factor(yprev_active, levels = 1:n_states)
    )

    # MoreArgs for scalar values that don't vary across patients
    more_args <- list(
      t = t,
      treatment_effect = treatment_effect *
        exp(-treatment_effect_decay * (t - 1)),
      extra_params = extra_params
    )

    # Add covariates dynamically
    for (cov_name in available_covariates) {
      lp_args[[cov_name]] <- X_active[[cov_name]]
    }

    # Special handling for treatment variable (might be tx_val in function)
    if ("tx" %in% names(X_active) && !"tx" %in% names(lp_args)) {
      lp_args[["tx"]] <- X_active[["tx"]]
    }

    n_thresholds <- length(intercepts)
    if (isTRUE(attr(lp_function, "mostr_vectorized_lp"))) {
      # Package models compute all active patients in one call.
      lp_matrix <- normalize_vectorized_markov_lp(
        do.call(lp_function, c(lp_args, more_args)),
        n_active = length(active_idx),
        n_thresholds = n_thresholds
      )
    } else {
      # Arbitrary user functions retain scalar invocation and validation.
      lp_list <- do.call(
        mapply,
        c(
          list(FUN = lp_function, SIMPLIFY = FALSE, USE.NAMES = FALSE),
          lp_args,
          list(MoreArgs = more_args)
        )
      )
      lp_values <- vapply(
        seq_along(lp_list),
        function(i) {
          normalize_markov_lp(lp_list[[i]], n_thresholds, active_idx[i], t)
        },
        numeric(n_thresholds)
      )
      lp_matrix <- matrix(
        as.numeric(lp_values),
        nrow = length(active_idx),
        ncol = n_thresholds,
        byrow = TRUE
      )
    }

    # Calculate transition probabilities using scalar PO or threshold-specific
    # PPO linear predictors.
    cum_probs <- plogis(
      matrix(
        intercepts,
        nrow = length(active_idx),
        ncol = n_thresholds,
        byrow = TRUE
      ) +
        lp_matrix
    )

    if (
      n_thresholds > 1 &&
        any(
          cum_probs[, -1L, drop = FALSE] -
            cum_probs[, -n_thresholds, drop = FALSE] <
            -1e-12
        )
    ) {
      stop(
        "lp_function produced threshold-specific linear predictors that ",
        "create decreasing cumulative probabilities. This implies negative ",
        "state probabilities; use threshold-specific effects that preserve ",
        "nondecreasing cumulative probabilities across intercepts."
      )
    }

    # Convert cumulative probabilities to individual state probabilities
    prob_matrix <- cbind(cum_probs, 1) - cbind(0, cum_probs)

    # Probability columns run from the highest to the lowest state.
    y_new_active <- n_states + 1L - sample_categorical_rows(prob_matrix)

    # Update state matrix
    # First, carry forward previous state for all patients
    state_matrix[, t + 1] <- state_matrix[, t]
    # Then update only active patients with new states
    state_matrix[active_idx, t + 1] <- y_new_active
  }

  # Construct the patient-major long result directly.
  patient_rows <- rep(seq_len(N), each = follow_up_time)
  result <- data.frame(
    id = baseline_data$id[patient_rows],
    time = rep(seq_len(follow_up_time), times = N),
    y = as.vector(t(state_matrix[, -1L, drop = FALSE])),
    check.names = FALSE
  )
  baseline_columns <- setdiff(names(baseline_data), "id")
  result <- cbind(
    result,
    baseline_data[patient_rows, baseline_columns, drop = FALSE]
  )
  result$yprev <- as.vector(t(
    state_matrix[, -(follow_up_time + 1L), drop = FALSE]
  ))
  rownames(result) <- NULL

  return(result)
}


markov_preset_names <- c("actt1", "actt2", "violet")


# The model behind each sim_markov() preset: how to sample patients, the
# thresholds in the order sim_markov() uses (worst state first), the linear
# predictor and its coefficients, the covariates passed to it, and the
# default treatment effect. NULL means the user supplies the model.
markov_preset <- function(preset) {
  if (is.null(preset)) {
    return(NULL)
  }
  if (
    !is.character(preset) ||
      length(preset) != 1L ||
      !preset %in% markov_preset_names
  ) {
    stop(
      "`preset` must be one of ",
      paste0('"', markov_preset_names, '"', collapse = ", "),
      ", or NULL to supply your own model.",
      call. = FALSE
    )
  }

  switch(
    preset,
    actt1 = actt_markov_preset(actt1_markov_parameters()),
    actt2 = actt_markov_preset(actt2_markov_parameters()),
    violet = violet_markov_preset()
  )
}


actt_markov_preset <- function(parameters) {
  list(
    sample_patients = function(n_patients, treatment_prob) {
      data.frame(
        id = seq_len(n_patients),
        yprev = sample(
          parameters$baseline_states,
          size = n_patients,
          replace = TRUE,
          prob = parameters$baseline_probabilities
        ),
        tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
      )
    },
    intercepts = rev(parameters$orm_intercepts),
    lp_function = actt2_markov_lp,
    extra_params = parameters$extra_params,
    covariate_names = "tx",
    treatment_effect = -0.0730
  )
}


violet_markov_preset <- function() {
  violet <- violet_markov_parameters()
  list(
    sample_patients = function(n_patients, treatment_prob) {
      patients <- violet$baseline[
        sample.int(
          nrow(violet$baseline),
          size = n_patients,
          replace = TRUE,
          prob = violet$baseline$n
        ),
        ,
        drop = FALSE
      ]
      data.frame(
        id = seq_len(n_patients),
        yprev = patients$yprev,
        age = patients$age,
        sofa = patients$sofa,
        tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
      )
    },
    intercepts = violet$intercepts,
    lp_function = violet_markov_lp,
    extra_params = violet$extra_params,
    covariate_names = c("age", "sofa", "tx"),
    treatment_effect = log(0.75)
  )
}


actt2_markov_parameters <- function() {
  list(
    baseline_states = 4:7,
    baseline_probabilities = c(0.1374637, 0.5459826, 0.2090997, 0.1074540),
    # Fitted with states 1 and 2 merged; the first threshold splits them again
    # so that 9.3% of the control arm is in state 2 on day 14.
    orm_intercepts = c(
      -0.6465,
      -4.9314,
      -5.2230,
      -7.4169,
      -11.5321,
      -15.2262,
      -22.0693
    ),
    extra_params = c(
      "day" = -0.3084,
      "day'" = 2.5186,
      "day''" = -5.0276,
      "day'''" = 3.1714,
      "day''''" = -1.0672,
      "yprev=1" = 0,
      "yprev=2" = 0,
      "yprev=3" = 6.2218,
      "yprev=4" = 7.1412,
      "yprev=5" = 9.9809,
      "yprev=6" = 14.1645,
      "yprev=7" = 19.3517
    )
  )
}


# The ACTT-2 model with ACTT-1's reported starting states (Beigel et al.,
# 2020). The death threshold is raised from -22.0693 so that, by exact
# calculation, 18% of the control arm has died by day 28.
actt1_markov_parameters <- function() {
  parameters <- actt2_markov_parameters()
  baseline_counts <- c(138, 435, 193, 285)
  parameters$baseline_probabilities <- baseline_counts / sum(baseline_counts)
  parameters$orm_intercepts[7] <- -21.4071
  parameters
}


actt2_markov_lp <- function(yprev, t, tx, treatment_effect, extra_params) {
  basis <- actt_markov_rcs_basis(t)
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  time_effect <- drop(basis %*% extra_params[day_names])
  previous_state_effect <- extra_params[
    paste0("yprev=", as.character(yprev))
  ]

  as.numeric(time_effect + previous_state_effect + treatment_effect * tx)
}


violet_markov_parameters <- function() {
  # Harrell (2022) fitted VGAM::vglm(ordered(status) ~ pstatus + lsp(day, 2) +
  # age + sofa + tx, cumulative(reverse = TRUE, parallel = FALSE ~ lsp(day, 2)))
  # with the states ordered from dead (lowest) to home (highest), and published
  # the coefficients rounded to three decimals:
  # https://hbiostat.org/data/repo/simlongord
  # mostr numbers the states from home (1) to dead (4), so every value below is
  # Harrell's printed value with the opposite sign. Thresholds are ordered for
  # Pr(Y >= 4), Pr(Y >= 3), and Pr(Y >= 2), which are one minus Harrell's
  # cumulative probabilities 1, 2, and 3.
  # Age is the exception: Harrell printed -0.011, but refitting his model to
  # all of `simlongord` gives -0.01047 (SE 0.00015), while every other
  # coefficient agrees with its printed value. 0.0105 is the value within the
  # printed value's rounding interval nearest that fit; with 0.011, simulated
  # patients reach home about 1 percentage point less often than in
  # `simlongord` by day 28.
  list(
    # Lazy data is not visible from the namespace unless mostr is attached.
    baseline = mostr::violet_baseline,
    intercepts = c(-4.713, 1.741, 10.034),
    extra_params = c(
      "yprev=1" = -15.797,
      "yprev=2" = -6.364,
      "yprev=3" = 0,
      "day:1" = 0.412,
      "day:2" = -0.481,
      "day:3" = -1.315,
      "day':1" = -0.410,
      "day':2" = 0.525,
      "day':3" = 1.384,
      "age" = 0.0105,
      "sofa" = 0.046
    )
  )
}


violet_markov_lp <- function(
  yprev,
  t,
  age,
  sofa,
  tx,
  treatment_effect,
  extra_params
) {
  previous_state_effect <- extra_params[paste0("yprev=", as.character(yprev))]
  shared_effect <- previous_state_effect +
    extra_params[["age"]] * age +
    extra_params[["sofa"]] * sofa +
    treatment_effect * tx
  # Linear spline in time with a knot at day 2; its slopes differ by threshold.
  time_effect <- t * extra_params[c("day:1", "day:2", "day:3")] +
    max(t - 2, 0) * extra_params[c("day':1", "day':2", "day':3")]

  outer(unname(shared_effect), unname(time_effect), `+`)
}
attr(violet_markov_lp, "mostr_vectorized_lp") <- TRUE


actt_markov_rcs_basis <- function(day) {
  markov_rcs_basis(day, c(2, 6, 11, 16, 21, 27))
}


markov_rcs_basis <- function(day, knots) {
  n_knots <- length(knots)
  last_knot <- knots[n_knots]
  penultimate_knot <- knots[n_knots - 1L]
  scale <- (last_knot - knots[1L])^(2 / 3)
  nonlinear <- matrix(0, nrow = length(day), ncol = n_knots - 2L)

  for (column in seq_len(n_knots - 2L)) {
    nonlinear[, column] <-
      pmax((day - knots[column]) / scale, 0)^3 +
      ((penultimate_knot - knots[column]) *
        pmax((day - last_knot) / scale, 0)^3 -
        (last_knot - knots[column]) *
          pmax((day - penultimate_knot) / scale, 0)^3) /
        (last_knot - penultimate_knot)
  }

  cbind(day, nonlinear)
}


validate_markov_wrapper_inputs <- function(n_patients, treatment_prob) {
  if (
    length(n_patients) != 1L ||
      !is.numeric(n_patients) ||
      !is.finite(n_patients) ||
      n_patients <= 0 ||
      n_patients > .Machine$integer.max ||
      n_patients != as.integer(n_patients)
  ) {
    stop(simpleError(
      "n_patients must be a positive integer",
      call = sys.call(-1)
    ))
  }

  if (
    length(treatment_prob) != 1L ||
      !is.numeric(treatment_prob) ||
      !is.finite(treatment_prob) ||
      treatment_prob < 0 ||
      treatment_prob > 1
  ) {
    stop(simpleError(
      "treatment_prob must be a numeric scalar between 0 and 1",
      call = sys.call(-1)
    ))
  }
}


validate_markov_treatment_effect <- function(
  treatment_effect,
  treatment_effect_decay
) {
  if (
    length(treatment_effect) != 1L ||
      !is.numeric(treatment_effect) ||
      !is.finite(treatment_effect)
  ) {
    stop(simpleError(
      "treatment_effect must be a finite numeric scalar",
      call = sys.call(-1)
    ))
  }

  if (
    length(treatment_effect_decay) != 1L ||
      !is.numeric(treatment_effect_decay) ||
      !is.finite(treatment_effect_decay) ||
      treatment_effect_decay < 0
  ) {
    stop(simpleError(
      "treatment_effect_decay must be a nonnegative finite numeric scalar",
      call = sys.call(-1)
    ))
  }
}


validate_markov_follow_up_time <- function(follow_up_time) {
  if (
    length(follow_up_time) != 1L ||
      !is.numeric(follow_up_time) ||
      !is.finite(follow_up_time) ||
      follow_up_time < 1 ||
      follow_up_time != as.integer(follow_up_time)
  ) {
    stop(simpleError(
      "follow_up_time must be a positive integer",
      call = sys.call(-1)
    ))
  }
}


normalize_markov_lp <- function(lp, n_thresholds, active_index, time) {
  if (!is.numeric(lp) || is.matrix(lp) || is.array(lp)) {
    stop("lp_function must return a numeric scalar or numeric vector")
  }

  lp <- as.numeric(lp)
  if (!(length(lp) %in% c(1, n_thresholds))) {
    stop(
      "lp_function must return either length 1 or length(intercepts) (",
      n_thresholds,
      ") for each patient; got length ",
      length(lp),
      " for active patient index ",
      active_index,
      " at time ",
      time
    )
  }

  if (any(!is.finite(lp))) {
    stop("lp_function must return finite numeric values")
  }

  if (length(lp) == 1) {
    rep(lp, n_thresholds)
  } else {
    lp
  }
}


# Linear predictors from a package model that handles all active patients in
# one call: one value per patient, or one row per patient and one column per
# threshold.
normalize_vectorized_markov_lp <- function(lp, n_active, n_thresholds) {
  if (is.null(dim(lp)) && length(lp) == n_active) {
    lp <- matrix(lp, nrow = n_active, ncol = n_thresholds)
  }
  if (
    !is.numeric(lp) ||
      !identical(as.integer(dim(lp)), as.integer(c(n_active, n_thresholds)))
  ) {
    stop(
      "A vectorized lp_function must return one value per patient or a ",
      "matrix with one row per patient and one column per threshold"
    )
  }
  if (any(!is.finite(lp))) {
    stop("lp_function must return finite numeric values")
  }

  lp
}
