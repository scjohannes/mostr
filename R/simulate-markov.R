#' Simulate Individual Patient Trajectories Using Markov Model
#'
#' Generates daily (or other discrete-time) ordinal states for each patient
#' from a first-order Markov model: the state at time `t` depends on the state
#' at time `t - 1`, the patient's baseline covariates, and time, through a
#' proportional odds or partial proportional odds model that you specify.
#' For ready-made models of published trials, use [sim_actt1_markov()],
#' [sim_actt2_markov()], or [sim_violet_markov()].
#'
#' @param baseline_data A data frame with one row per patient. Must include
#'   `id` and `yprev`, the state immediately before the first simulated time,
#'   and every covariate that `lp_function` uses (for example `tx`).
#' @param follow_up_time Integer. Number of time periods to simulate (default: 60).
#' @param intercepts Numeric vector of thresholds, one fewer than the number of
#'   states, in increasing order. With `K` states numbered from 1 (best) to `K`
#'   (worst), `intercepts[j]` plus the linear predictor is the log odds of
#'   being in state `K + 1 - j` or worse: `intercepts[1]` belongs to the worst
#'   state alone and `intercepts[K - 1]` to states 2 and worse.
#' @param lp_function A function that calculates the linear predictor for one
#'   patient at one time point. It is called with these named arguments:
#'   - `yprev`: previous state, a factor with levels `1` to `K`
#'   - `t`: current time
#'   - one argument for each column named in `covariate_names`, and `tx` when
#'     `baseline_data` has a `tx` column
#'   - `parameter`: the value of `parameter`
#'   - `extra_params`: the value of `extra_params`
#'
#'   It must return either one number, added to every threshold (a
#'   proportional odds model), or a numeric vector with `length(intercepts)`
#'   elements, added to the corresponding thresholds (a partial proportional
#'   odds model). Larger values make worse states more likely.
#' @param extra_params Optional object, usually a named numeric vector of
#'   coefficients, passed unchanged to `lp_function` (default: `NULL`).
#' @param parameter Numeric. Passed to `lp_function`; typically the treatment
#'   effect on the log odds ratio scale (default: 0, meaning odds ratio = 1).
#' @param absorbing_states Integer vector of states that are absorbing (once
#'   entered, patients cannot leave). Default is the worst state,
#'   `length(intercepts) + 1`, which is death in all models in this package.
#'   Use `integer()` when no state is absorbing.
#' @param seed Integer. Random seed for reproducibility (default: NULL, no seed set).
#' @param covariate_names Character vector of `baseline_data` columns passed to
#'   `lp_function` as arguments with the same names. Default is every column
#'   except `id` and `yprev`. A `tx` column is always passed.
#'
#' @return A data frame with one row per patient and time point and columns:
#'   - id: patient identifier
#'   - time: time point (1 to follow_up_time)
#'   - y: state at this time
#'   - yprev: state at the previous time point
#'   - all other columns from baseline_data
#'
#' @details
#' At each time step the function:
#' 1. Calculates the linear predictor for each patient not in an absorbing
#'    state
#' 2. Adds it to the intercepts to get the probabilities of each state or
#'    worse
#' 3. Converts these to the probability of each state
#' 4. Samples the next state for each patient
#' 5. Keeps patients in absorbing states where they are
#'
#' States must be consecutive integers from 1 (best) to `K` (worst). If
#' `lp_function` returns threshold-specific linear predictors, the probability
#' of being in a state or worse must not increase from one state to the next
#' worse state; otherwise the function stops because the implied state
#' probabilities would be negative.
#'
#' @examples
#' # Three states: 1 = home, 2 = in hospital, 3 = dead (absorbing by default)
#' baseline <- data.frame(id = 1:20, yprev = 2, tx = rep(0:1, 10))
#' hospital_lp <- function(yprev, t, tx, parameter, extra_params) {
#'   # Hospitalized patients are more likely to stay in hospital or die
#'   5 * (yprev == "2") - 0.1 * t + parameter * tx
#' }
#' trajectories <- sim_trajectories_markov(
#'   baseline_data = baseline,
#'   follow_up_time = 10,
#'   intercepts = c(-9, -3),
#'   lp_function = hospital_lp,
#'   parameter = log(0.8),
#'   seed = 1
#' )
#' head(trajectories)
#'
#' # Partial proportional odds: treatment affects only the odds of death
#' death_only_lp <- function(yprev, t, tx, parameter, extra_params) {
#'   hospital_lp(yprev, t, tx = 0, parameter, extra_params) +
#'     c(parameter * tx, 0)
#' }
#' trajectories_ppo <- sim_trajectories_markov(
#'   baseline_data = baseline,
#'   follow_up_time = 10,
#'   intercepts = c(-9, -3),
#'   lp_function = death_only_lp,
#'   parameter = log(0.8),
#'   seed = 1
#' )
#'
#' @seealso [sim_violet_markov()], [sim_actt1_markov()], and
#'   [sim_actt2_markov()] for models of published trials.
#'
#' @importFrom stats plogis
#'
#' @export
sim_trajectories_markov <- function(
  baseline_data,
  follow_up_time = 60,
  intercepts,
  lp_function,
  extra_params = NULL,
  parameter = 0,
  absorbing_states = length(intercepts) + 1L,
  seed = NULL,
  covariate_names = setdiff(names(baseline_data), c("id", "yprev"))
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

  if (!is.null(seed)) {
    set.seed(seed)
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
      parameter = parameter,
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

#' Simulate ACTT-2 Markov Trajectories
#'
#' Generates synthetic daily health states using a proportional-odds model
#' based on data from the
#' [Adaptive COVID-19 Treatment Trial 2 (ACTT-2)](https://doi.org/10.1056/NEJMoa2031994),
#' a clinical trial in adults hospitalized with moderate to severe COVID-19.
#' States follow the trial's
#' [NIAID eight-point ordinal scale](https://www.fda.gov/media/144473/download#page=10),
#' from unrestricted activities outside hospital (1) to death (8).
#' Use `n_patients` to choose the sample size and `treatment_effect` to choose
#' the treatment effect in the simulated data.
#'
#' @param n_patients Integer. Number of patients to simulate (default: 1000).
#' @param treatment_prob Numeric. Probability of assignment to Baricitinib plus
#'   Remdesivir, represented by `tx = 1` (default: 0.5). Placebo plus Remdesivir
#'   is represented by `tx = 0`.
#' @param treatment_effect Numeric. Day-1 log-odds coefficient for Baricitinib
#'   plus Remdesivir (`tx = 1`) relative to Placebo plus Remdesivir (`tx = 0`).
#'   The default -0.0730 is the coefficient fitted to the ACTT-2 data.
#' @param treatment_effect_decay Nonnegative exponential decay rate per day.
#'   The effect on day `t` is `treatment_effect *
#'   exp(-treatment_effect_decay * (t - 1))`. The default 0 gives a constant
#'   effect.
#' @param seed Integer. Random seed for reproducibility (default: `NULL`).
#' @param ... Additional arguments passed to [sim_trajectories_markov()]. These
#'   can override ACTT-2 defaults such as `follow_up_time`. Overrides are applied
#'   after the wrapper generates its baseline data; supplying `baseline_data`
#'   replaces those generated baseline rows.
#'
#' @return A data frame with the same columns returned by
#'   [sim_trajectories_markov()]. The default output contains 28 rows per
#'   patient and states numbered 1 through 8.
#'
#' @details
#' Initial previous states are sampled from states 4 through 7 with
#' probabilities 0.1374637, 0.5459826, 0.2090997, and 0.1074540. They represent
#' the states immediately before the first generated day-1 outcome. State 8 is
#' absorbing by default. Set `follow_up_time = 60` for 60 days; days beyond
#' 28 extrapolate the fitted time effect.
#'
#' In ACTT-2, state 2 was mostly recorded only at the day 14, 21, and 28
#' visits, so the transition model was fitted with states 1 and 2 merged. An
#' added lowest threshold splits them again, placing the observed 9.3% of the
#' control arm in state 2 on day 14. Both states share the same coefficients,
#' so the split does not change the occupancy of states 3 to 8. The treatment
#' effect is supplied through `treatment_effect`.
#'
#' @examples
#' actt2 <- sim_actt2_markov(n_patients = 20, seed = 123)
#' head(actt2)
#' decaying <- sim_actt2_markov(
#'   n_patients = 20,
#'   treatment_effect = log(0.8),
#'   treatment_effect_decay = -log(0.01) / 14,
#'   seed = 123
#' )
#'
#' @export
sim_actt2_markov <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = -0.0730,
  treatment_effect_decay = 0,
  seed = NULL,
  ...
) {
  validate_markov_wrapper_inputs(n_patients, treatment_prob)

  validate_markov_treatment_effect(treatment_effect, treatment_effect_decay)

  n_patients <- as.integer(n_patients)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  actt2 <- actt2_markov_parameters()
  baseline_data <- data.frame(
    id = seq_len(n_patients),
    yprev = sample(
      actt2$baseline_states,
      size = n_patients,
      replace = TRUE,
      prob = actt2$baseline_probabilities
    ),
    tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
  )

  args <- list(
    baseline_data = baseline_data,
    follow_up_time = 28L,
    intercepts = rev(actt2$orm_intercepts),
    lp_function = actt2_markov_lp,
    extra_params = c(
      actt2$extra_params,
      treatment_effect_decay = treatment_effect_decay
    ),
    parameter = treatment_effect,
    absorbing_states = 8L,
    seed = NULL,
    covariate_names = "tx"
  )
  args <- override_markov_wrapper_args(args, list(...))

  do.call(sim_trajectories_markov, args)
}


#' Simulate ACTT-1 ordinal patient trajectories
#'
#' @description
#' Simulates individual 8-state trajectories resembling the
#' [Adaptive COVID-19 Treatment Trial 1 (ACTT-1)](https://doi.org/10.1056/NEJMoa2007764)
#' of Remdesivir versus placebo in adults hospitalized with COVID-19 (Beigel et
#' al., 2020). The transition model is a simplified version of the Bayesian
#' Markov partial proportional-odds model that Rohde et al. (2024) fitted to
#' the ACTT-1 data. The wrapper uses that model's posterior mean coefficients,
#' the baseline-state counts reported by Beigel et al. (2020), and the ACTT-2
#' restricted cubic spline knots because the original ACTT-1 knot locations
#' are unavailable.
#'
#' @param n_patients Number of patients to simulate. Must be a positive integer.
#' @param treatment_prob Probability of assignment to Remdesivir. Must be a
#'   numeric scalar between 0 and 1.
#' @param treatment_effect Day-1 proportional-odds coefficient applied to
#'   Remdesivir (`tx = 1`) relative to placebo (`tx = 0`). The default -0.0730
#'   matches [sim_actt2_markov()].
#' @param treatment_effect_decay Nonnegative exponential decay rate per day.
#'   The effect on day `t` is `treatment_effect *
#'   exp(-treatment_effect_decay * (t - 1))`. The default 0 gives a constant
#'   effect.
#' @param seed Optional random-number seed. The seed is set once before baseline
#'   states and treatment assignments are sampled, and trajectory generation
#'   continues from the same random-number stream.
#' @param ... Advanced overrides passed to [sim_trajectories_markov()]. These
#'   can override defaults such as `follow_up_time`. Overrides are applied after
#'   the wrapper generates its baseline data; supplying `baseline_data` replaces
#'   those generated baseline rows.
#'
#' @return A data frame with the same columns returned by
#'   [sim_trajectories_markov()]. The default output contains 28 rows per
#'   patient and states numbered 1 through 8.
#'
#' @details
#' Initial previous states are sampled from states 4 through 7 according to the
#' reported counts 138, 435, 193, and 285. These counts are normalized over the
#' 1,051 patients with a reported state and represent states immediately before
#' the first generated day-1 outcome. State 8 is absorbing by default.
#'
#' Treatment is coded as `tx = 0` for placebo and `tx = 1` for Remdesivir. The
#' user-specified effect replaces the fitted ACTT-1 treatment main effect and
#' treatment-by-day spline interactions.
#'
#' Age and sex are not simulated because their joint distribution with the
#' baseline state is not reported. The fitted thresholds describe a woman aged
#' 0 years, the only patient whose age and sex terms are all zero, so every
#' threshold is shifted by the average contribution of an ACTT-1 patient,
#' 0.2253 on the log-odds scale. For the spline age effect
#' \eqn{f(a) = 0.0017a + 0.0172a' - 0.0586a''}{f(a) = 0.0017 a + 0.0172 a' - 0.0586 a''},
#' ages are assumed normal with the reported mean 58.9 and standard deviation
#' 15.0 years, truncated at 18, and the unreported knots are placed at the 5th,
#' 35th, 65th, and 95th percentiles of that normal distribution, as
#' `rms::rcs()` would. The average age contribution,
#' \eqn{\int_{18}^{\infty} f(a)\phi(a)\,da / (1 - \Phi(18))}{the integral of f(a) phi(a) / (1 - Phi(18)) over a >= 18},
#' is 0.2215, where \eqn{\phi}{phi} and \eqn{\Phi}{Phi} are the normal density
#' and distribution function; it was computed with [stats::integrate()]. Men add
#' \eqn{0.0059 \times 684/1062 = 0.0038}{0.0059 * 684/1062 = 0.0038}.
#' Simulating individual ages instead, independent of baseline state, changes
#' state occupancy probabilities by less than 0.3 percentage points.
#'
#' The fitted constrained partial proportional-odds time deviation is retained.
#' For the cumulative cutoff `y >= k`, it adds `-0.0188 * day * k` to the linear
#' predictor. This threshold-specific term is essential for reproducing ACTT-1
#' mortality over follow-up even though its printed coefficient is small.
#'
#' @references
#' Beigel JH, Tomashek KM, Dodd LE, et al. (2020). "Remdesivir for the
#' Treatment of Covid-19 - Final Report." *New England Journal of Medicine*,
#' 383(19), 1813-1826. \doi{10.1056/NEJMoa2007764}
#'
#' Rohde MD, French B, Stewart TG, Harrell FE Jr (2024). "Bayesian Transition
#' Models for Ordinal Longitudinal Outcomes." *Statistics in Medicine*, 43(18),
#' 3539-3561. \doi{10.1002/sim.10133}
#'
#' @examples
#' actt1 <- sim_actt1_markov(n_patients = 20, seed = 123)
#' head(actt1)
#' decaying <- sim_actt1_markov(
#'   n_patients = 20,
#'   treatment_effect = log(0.8),
#'   treatment_effect_decay = -log(0.01) / 14,
#'   seed = 123
#' )
#'
#' @export
sim_actt1_markov <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = -0.0730,
  treatment_effect_decay = 0,
  seed = NULL,
  ...
) {
  validate_markov_wrapper_inputs(n_patients, treatment_prob)
  validate_markov_treatment_effect(treatment_effect, treatment_effect_decay)

  n_patients <- as.integer(n_patients)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  actt1 <- actt1_markov_parameters()
  baseline_data <- data.frame(
    id = seq_len(n_patients),
    yprev = sample(
      actt1$baseline_states,
      size = n_patients,
      replace = TRUE,
      prob = actt1$baseline_probabilities
    ),
    tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
  )

  args <- list(
    baseline_data = baseline_data,
    follow_up_time = 28L,
    intercepts = rev(actt1$orm_intercepts + actt1$average_age_sex_effect),
    lp_function = actt1_markov_lp,
    extra_params = c(
      actt1$extra_params,
      treatment_effect_decay = treatment_effect_decay
    ),
    parameter = treatment_effect,
    absorbing_states = 8L,
    seed = NULL,
    covariate_names = "tx"
  )
  args <- override_markov_wrapper_args(args, list(...))

  do.call(sim_trajectories_markov, args)
}


#' Simulate VIOLET ordinal patient trajectories
#'
#' @description
#' Simulates daily 4-state trajectories resembling the
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
#' Use `n_patients` to choose the sample size and `treatment_effect` to choose
#' the treatment effect in the simulated data.
#'
#' @param n_patients Number of patients to simulate. Must be a positive integer.
#' @param treatment_prob Probability of assignment to treatment (`tx = 1`). Must
#'   be a numeric scalar between 0 and 1. Control patients have `tx = 0`.
#' @param treatment_effect Log odds ratio of being in a worse state on the next
#'   day for `tx = 1` relative to `tx = 0`, the same for every previous state.
#'   Values below 0 mean that treatment helps. The default `log(0.75)` is the
#'   odds ratio of 0.75 that Harrell built into `simlongord`. In VIOLET itself,
#'   vitamin D had no detectable effect (fitted log odds ratio -0.028, standard
#'   error 0.050).
#' @param treatment_effect_decay Nonnegative exponential decay rate per day.
#'   The effect on day `t` is `treatment_effect *
#'   exp(-treatment_effect_decay * (t - 1))`. The default 0 gives a constant
#'   effect, as in Harrell's model.
#' @param seed Optional random-number seed. The seed is set once before
#'   patients and treatment assignments are sampled, and trajectory generation
#'   continues from the same random-number stream.
#' @param ... Advanced overrides passed to [sim_trajectories_markov()], such as
#'   `follow_up_time`. A supplied `baseline_data` replaces the sampled patients
#'   and must contain the columns `id`, `yprev`, `age`, `sofa`, and `tx`.
#'
#' @return A data frame with columns `id`, `time`, `y`, `yprev`, `age`, `sofa`,
#'   and `tx`, as described in [sim_trajectories_markov()]. The default output
#'   contains 28 rows per patient, and dead patients stay in state 4.
#'
#' @details
#' **Patients.** Each simulated patient's starting state (`yprev`), age, and
#' SOFA (Sequential Organ Failure Assessment) score are drawn together, with
#' replacement, from the 1,351 VIOLET patients summarized in
#' [violet_baseline]. About 65% of patients start in hospital (state 2) and
#' 35% on a ventilator or with ARDS (state 3). Treatment is assigned
#' independently of these characteristics with probability `treatment_prob`.
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
#' effect on day \eqn{t} described under `treatment_effect_decay`. Positive
#' values make worse states more likely: for example, each additional SOFA
#' point multiplies the odds of being in a worse state by
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
#' **Relation to Harrell's model.** Harrell numbered the states from dead
#' (lowest) to home (highest), modeled the probability of a better state, and
#' published the coefficients rounded to three decimals. This function numbers
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
#' 33.9% versus 36.0% start on a ventilator). This function assigns treatment
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
#' violet <- sim_violet_markov(n_patients = 20, seed = 123)
#' head(violet)
#'
#' # A trial in which treatment has no effect
#' null_trial <- sim_violet_markov(
#'   n_patients = 20,
#'   treatment_effect = 0,
#'   seed = 123
#' )
#'
#' @seealso [violet_baseline] for the VIOLET patients, and
#'   [sim_trajectories_markov()] to simulate from your own model.
#'
#' @export
sim_violet_markov <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = log(0.75),
  treatment_effect_decay = 0,
  seed = NULL,
  ...
) {
  validate_markov_wrapper_inputs(n_patients, treatment_prob)
  validate_markov_treatment_effect(treatment_effect, treatment_effect_decay)

  n_patients <- as.integer(n_patients)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  violet <- violet_markov_parameters()
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
  baseline_data <- data.frame(
    id = seq_len(n_patients),
    yprev = patients$yprev,
    age = patients$age,
    sofa = patients$sofa,
    tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
  )

  args <- list(
    baseline_data = baseline_data,
    follow_up_time = 28L,
    intercepts = violet$intercepts,
    lp_function = violet_markov_lp,
    extra_params = c(
      violet$extra_params,
      treatment_effect_decay = treatment_effect_decay
    ),
    parameter = treatment_effect,
    absorbing_states = 4L,
    seed = NULL,
    covariate_names = c("age", "sofa", "tx")
  )
  args <- override_markov_wrapper_args(args, list(...))

  do.call(sim_trajectories_markov, args)
}


actt1_markov_parameters <- function() {
  baseline_counts <- c(138, 435, 193, 285)

  list(
    baseline_states = 4:7,
    baseline_probabilities = baseline_counts / sum(baseline_counts),
    orm_intercepts = c(
      -3.3535,
      -3.5896,
      -3.8395,
      -6.0379,
      -9.7465,
      -12.3361,
      -18.8334
    ),
    # Average of 0.0017 * age + 0.0172 * age' - 0.0586 * age'' for ages
    # N(58.9, 15) truncated at 18 (knots at the 5th, 35th, 65th, and 95th
    # percentiles), plus 0.0059 * 684 / 1062 for male sex.
    average_age_sex_effect = 0.2253,
    extra_params = c(
      "day" = -0.1110,
      "day'" = 1.6444,
      "day''" = -3.1537,
      "day'''" = 1.2387,
      "day''''" = 1.1910,
      "yprev=1" = 0,
      "yprev=2" = -0.0402,
      "yprev=3" = 3.9042,
      "yprev=4" = 5.0754,
      "yprev=5" = 8.0171,
      "yprev=6" = 11.5073,
      "yprev=7" = 16.1963,
      "day x f(y)" = -0.0188
    )
  )
}


actt1_markov_lp <- function(yprev, t, tx, parameter, extra_params) {
  basis <- actt_markov_rcs_basis(t)
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  time_effect <- drop(basis %*% extra_params[day_names])
  previous_state_effect <- extra_params[
    paste0("yprev=", as.character(yprev))
  ]
  treatment_effect <- decayed_markov_treatment_effect(
    parameter,
    t,
    extra_params
  )
  shared_effect <- time_effect + previous_state_effect + treatment_effect * tx
  cutoff_scores <- rev(2:8)
  partial_po_effect <- extra_params[["day x f(y)"]] * t * cutoff_scores

  as.numeric(shared_effect + partial_po_effect)
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


actt2_markov_lp <- function(yprev, t, tx, parameter, extra_params) {
  basis <- actt_markov_rcs_basis(t)
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  time_effect <- drop(basis %*% extra_params[day_names])
  previous_state_effect <- extra_params[
    paste0("yprev=", as.character(yprev))
  ]

  treatment_effect <- decayed_markov_treatment_effect(
    parameter,
    t,
    extra_params
  )

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


violet_markov_lp <- function(yprev, t, age, sofa, tx, parameter, extra_params) {
  previous_state_effect <- extra_params[paste0("yprev=", as.character(yprev))]
  treatment_effect <- decayed_markov_treatment_effect(
    parameter,
    t,
    extra_params
  )
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


# Replace whole arguments: utils::modifyList() would recurse into a supplied
# baseline_data data frame and merge it column by column into the generated one.
override_markov_wrapper_args <- function(args, overrides) {
  override_names <- names(overrides)
  if (
    length(overrides) > 0 &&
      (is.null(override_names) || !all(nzchar(override_names)))
  ) {
    stop("All arguments passed through `...` must be named.", call. = FALSE)
  }

  args[override_names] <- overrides
  args
}


decayed_markov_treatment_effect <- function(parameter, time, extra_params) {
  decay <- if ("treatment_effect_decay" %in% names(extra_params)) {
    extra_params[["treatment_effect_decay"]]
  } else {
    0
  }

  parameter * exp(-decay * (time - 1))
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
