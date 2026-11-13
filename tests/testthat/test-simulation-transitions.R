# Do trajectories from sim_markov() follow the transition model they are
# simulated from?
#
# sim_markov() is the package's only simulator: it draws first-order Markov
# trajectories from given coefficients, either a preset's or a user's
# `lp_function`. Nothing in the package simulates from a fitted model or from a
# second-order model, so there is nothing of that kind to check here.
#
# Three kinds of check:
# - a degenerate model that puts all probability on one next state, so that
#   the whole path is known in advance and must be reproduced exactly;
# - properties that hold for any random draw (valid states, `yprev` equal to
#   the previous `y`, absorbing states never left);
# - with many patients, transition frequencies against the model's
#   probabilities, computed here with plogis() of the cumulative logits.
#
# The first two hold for any seed and run on CRAN.

# Degenerate model -------------------------------------------------------------

degenerate_states <- 5L

# The next state that the degenerate model makes certain: state 5 on day `x`,
# otherwise a cycle through states 1 to 4 that depends on the previous state,
# the day, `x` and `tx`, so that a wrong day, a wrong previous state or another
# patient's covariates all give a different path.
degenerate_next_state <- function(yprev, t, x, tx) {
  ifelse(t == x, 5, (yprev + t + x + 2 * tx) %% 4 + 1)
}

# sim_markov() adds the linear predictor to its thresholds, and threshold j is
# the cutoff Pr(Y >= K + 1 - j). Both models below put the cumulative logit of
# every cutoff at or below the target state at +50 or more and every other at
# -50 or less. plogis(-50) is about 2e-22, far below the smallest uniform draw
# R makes (about 1e-10), so no seed can draw any state but the target.
#
# Proportional odds: thresholds 100 apart, so one shift of all of them places
# the target's cutoff at +50 and the next worse cutoff at -50.
degenerate_po <- list(
  intercepts = 100 * seq_len(degenerate_states - 1L),
  shift = function(state) -100 * (degenerate_states + 1 - state) + 50
)
# Partial proportional odds: one shift per threshold.
degenerate_ppo <- list(
  intercepts = c(-0.3, -0.1, 0.1, 0.3),
  shift = function(state) {
    cutoff <- degenerate_states + 1 - seq_len(degenerate_states - 1L)
    outer(state, cutoff, function(s, k) ifelse(k <= s, 50, -50))
  }
)

# `vectorized = TRUE` marks the function the way the presets' functions are
# marked, so that sim_markov() calls it once for all patients instead of once
# per patient.
degenerate_lp <- function(model, vectorized) {
  lp <- function(yprev, t, x, tx, treatment_effect, extra_params) {
    state <- degenerate_next_state(as.numeric(as.character(yprev)), t, x, tx)
    shift <- model$shift(state)
    if (vectorized) shift else drop(shift)
  }
  if (vectorized) {
    attr(lp, "mostr_vectorized_lp") <- TRUE
  }
  lp
}

# The path that the degenerate model fixes, in plain loops that call no
# package code: patients in an absorbing state stay there, everyone else moves
# to degenerate_next_state().
degenerate_expected_path <- function(baseline, follow_up_time, absorbing) {
  rows <- lapply(seq_len(nrow(baseline)), function(i) {
    state <- baseline$yprev[i]
    y <- numeric(follow_up_time)
    yprev <- numeric(follow_up_time)
    for (t in seq_len(follow_up_time)) {
      yprev[t] <- state
      if (!state %in% absorbing) {
        state <- degenerate_next_state(state, t, baseline$x[i], baseline$tx[i])
      }
      y[t] <- state
    }
    data.frame(
      id = baseline$id[i],
      time = seq_len(follow_up_time),
      y = y,
      yprev = yprev
    )
  })
  do.call(rbind, rows)
}

expect_degenerate_path <- function(baseline, follow_up_time, model, vectorized,
                                   absorbing, label) {
  simulated <- sim_markov(
    NULL,
    baseline_data = baseline,
    follow_up_time = follow_up_time,
    intercepts = model$intercepts,
    lp_function = degenerate_lp(model, vectorized),
    absorbing_states = absorbing,
    seed = 1
  )
  # `absorbing_states = NULL` means the worst state.
  expected <- degenerate_expected_path(
    baseline,
    follow_up_time,
    absorbing %||% degenerate_states
  )

  # Covariates stay with their patient on every day. One expectation per
  # case, because each expectation costs several milliseconds and this file
  # runs on CRAN.
  expected$x <- rep(baseline$x, each = follow_up_time)
  expected$tx <- rep(baseline$tx, each = follow_up_time)
  expect_equal(as.list(simulated), as.list(expected), info = label)
}

test_that("a degenerate transition model reproduces its one possible path", {
  # `x` is the day on which the patient moves to state 5. Patients 5 and 10
  # start in state 5, and patients die on days 1 to 6, so absorbed and active
  # patients are mixed on most days.
  baseline <- data.frame(
    id = 1:10,
    yprev = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
    x = c(99, 99, 99, 99, 99, 2, 3, 4, 6, 1),
    tx = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
  )
  # Both kinds of model, called per patient and for all patients at once, with
  # the worst state absorbing (the default); then no absorbing state, and a
  # better state absorbing as well.
  cases <- list(
    list(model = degenerate_po, vectorized = FALSE, absorbing = NULL),
    list(model = degenerate_po, vectorized = TRUE, absorbing = NULL),
    list(model = degenerate_ppo, vectorized = FALSE, absorbing = NULL),
    list(model = degenerate_ppo, vectorized = TRUE, absorbing = NULL),
    list(model = degenerate_ppo, vectorized = FALSE, absorbing = integer()),
    list(model = degenerate_po, vectorized = TRUE, absorbing = c(2L, 5L))
  )
  for (case in cases) {
    expect_degenerate_path(
      baseline,
      follow_up_time = 8,
      model = case$model,
      vectorized = case$vectorized,
      absorbing = case$absorbing,
      label = paste(
        if (identical(case$model, degenerate_po)) "PO" else "PPO",
        if (case$vectorized) "vectorized" else "per patient",
        "absorbing:",
        paste(case$absorbing %||% "default", collapse = ", ")
      )
    )
  }
})

test_that("absorbed patients keep their rows after everyone is absorbed", {
  # Everyone reaches state 5 by day 2, so days 3 to 6 take the branch that
  # copies the last states forward without simulating.
  baseline <- data.frame(
    id = 1:3,
    yprev = c(1, 3, 4),
    x = c(1, 2, 2),
    tx = c(0, 1, 0)
  )
  for (vectorized in c(FALSE, TRUE)) {
    expect_degenerate_path(
      baseline,
      follow_up_time = 6,
      model = degenerate_po,
      vectorized = vectorized,
      absorbing = NULL,
      label = if (vectorized) "vectorized" else "per patient"
    )
  }
})

# A random model ---------------------------------------------------------------

# Four states, a partial proportional-odds model in which the previous state,
# the day and treatment shift all cutoffs, and the day also has a slope of its
# own for the best and worst cutoffs. In the usual notation,
#   logit Pr(Y_t >= k | Y_{t-1} = j) =
#     alpha_k + 1.2 (j - 2) - 0.35 t + delta_k t + beta tx,
# with alpha = (1, -1, -3) and delta = (-0.15, 0, 0.25) for k = 2, 3, 4.
random_model_probability <- function(yprev, t, tx, beta) {
  alpha <- c(1, -1, -3)
  delta <- c(-0.15, 0, 0.25)
  cumulative <- stats::plogis(
    alpha + 1.2 * (yprev - 2) - 0.35 * t + delta * t + beta * tx
  )
  c(1, cumulative) - c(cumulative, 0)
}
# The same model as sim_markov() takes it: thresholds and day slopes for the
# cutoffs Pr(Y >= 4), Pr(Y >= 3) and Pr(Y >= 2), in that order.
random_model_intercepts <- c(-3, -1, 1)
random_model_lp <- function(yprev, t, tx, treatment_effect, extra_params) {
  shared <- 1.2 * (as.numeric(as.character(yprev)) - 2) - 0.35 * t +
    treatment_effect * tx
  outer(shared, c(0.25, 0, -0.15) * t, `+`)
}
attr(random_model_lp, "mostr_vectorized_lp") <- TRUE

test_that("random trajectories are valid and consistent for any seed", {
  n_patients <- 200
  # Longer follow-up would make the day slopes cross the thresholds.
  follow_up_time <- 6
  baseline <- data.frame(
    id = seq_len(n_patients),
    yprev = rep(1:4, length.out = n_patients),
    tx = rep(0:1, each = n_patients / 2)
  )

  for (absorbing in list(4L, c(1L, 4L), integer())) {
    simulated <- sim_markov(
      NULL,
      baseline_data = baseline,
      follow_up_time = follow_up_time,
      intercepts = random_model_intercepts,
      lp_function = random_model_lp,
      treatment_effect = 0.6,
      absorbing_states = absorbing
    )
    # Columns are patients and rows are days.
    y <- matrix(simulated$y, nrow = follow_up_time)
    yprev <- matrix(simulated$yprev, nrow = follow_up_time)
    stays <- simulated$yprev %in% absorbing

    # One expectation per case, because each costs several milliseconds and
    # this test runs on CRAN.
    expect_equal(
      list(
        id = simulated$id,
        time = simulated$time,
        valid_states = all(simulated$y %in% 1:4),
        first_yprev = yprev[1, ],
        later_yprev = yprev[-1, ],
        y_after_absorbing = simulated$y[stays]
      ),
      list(
        # One row per patient and day, patient by patient.
        id = rep(baseline$id, each = follow_up_time),
        time = rep(seq_len(follow_up_time), n_patients),
        valid_states = TRUE,
        first_yprev = baseline$yprev,
        later_yprev = y[-follow_up_time, ],
        y_after_absorbing = simulated$yprev[stays]
      ),
      info = paste("absorbing:", paste(absorbing, collapse = ", "))
    )
  }
})

# Transition frequencies -------------------------------------------------------

# Compare, for every previous state, day and covariate value, the number of
# patients moving to each state with the number the model predicts. A count
# passes when it is within 4 binomial standard errors of its expectation, plus
# one patient, so that a single draw of a transition with a probability near 0
# does not fail. Rows starting in an absorbing state are left out; the checks
# above cover them. Returns the number of cells compared.
expect_transition_frequencies <- function(simulated, covariates, n_states,
                                          absorbing, probability, label) {
  by <- c("yprev", "time", covariates)
  simulated <- simulated[!simulated$yprev %in% absorbing, , drop = FALSE]
  y <- simulated$y
  # An integer code for each combination of the `by` columns; interaction()
  # and table() take most of this test's time on 36,000 rows.
  code <- 0
  for (column in by) {
    values <- match(simulated[[column]], sort(unique(simulated[[column]])))
    code <- code * max(values) + values
  }
  cells <- sort(unique(code))
  cell <- match(code, cells)
  observed <- matrix(
    tabulate((cell - 1L) * n_states + y, nbins = length(cells) * n_states),
    ncol = n_states,
    byrow = TRUE
  )
  first_rows <- simulated[match(cells, code), , drop = FALSE]
  p <- t(vapply(
    seq_len(nrow(first_rows)),
    function(i) probability(first_rows[i, , drop = FALSE]),
    numeric(n_states)
  ))
  expected <- rowSums(observed) * p
  excess <- abs(observed - expected) - (4 * sqrt(expected * (1 - p)) + 1)

  worst <- which.max(apply(excess, 1, max))
  expect_lte(max(excess), 0, label = paste0(
    label, ", worst cell ",
    paste(by, unlist(first_rows[worst, by]), sep = " = ", collapse = ", "),
    ": observed ",
    paste(observed[worst, ], collapse = " "), ", expected ",
    paste(round(expected[worst, ], 1), collapse = " "),
    "; largest excess over the allowed difference"
  ))
  nrow(observed)
}

test_that("transition frequencies match the model's transition probabilities", {
  skip_on_cran() # The pass depends on the seed.

  # A user's model with partial proportional odds.
  baseline <- data.frame(
    id = 1:9000,
    yprev = rep(1:3, 3000),
    tx = rep(0:1, each = 4500)
  )
  simulated <- sim_markov(
    NULL,
    baseline_data = baseline,
    follow_up_time = 3,
    intercepts = random_model_intercepts,
    lp_function = random_model_lp,
    treatment_effect = 0.6,
    seed = 20261113
  )
  cells <- expect_transition_frequencies(
    simulated,
    covariates = "tx",
    n_states = 4,
    absorbing = 4,
    probability = function(cell) {
      random_model_probability(cell$yprev, cell$time, cell$tx, beta = 0.6)
    },
    label = "user's model"
  )
  expect_equal(cells, 3 * 3 * 2)

  # The VIOLET preset, against the formula and table on its help page, with
  # its default treatment effect log(0.75). Age and SOFA come in two fixed
  # pairs, so age alone identifies a patient's covariates.
  baseline <- expand.grid(
    rep = 1:1000,
    yprev = 2:3,
    profile = 1:2,
    tx = 0:1
  )
  baseline <- data.frame(
    id = seq_len(nrow(baseline)),
    yprev = baseline$yprev,
    age = c(45, 75)[baseline$profile],
    sofa = c(4, 12)[baseline$profile],
    tx = baseline$tx
  )
  violet_probability <- function(cell) {
    # Cutoffs Pr(Y >= 2), Pr(Y >= 3) and Pr(Y >= 4).
    alpha <- c(10.034, 1.741, -4.713)
    delta <- c(-1.315, -0.481, 0.412)
    delta_after_day_2 <- c(1.384, 0.525, -0.410)
    gamma <- c(-15.797, -6.364, 0)[cell$yprev]
    t <- cell$time
    cumulative <- stats::plogis(
      alpha + gamma + delta * t + delta_after_day_2 * max(t - 2, 0) +
        0.0105 * cell$age + 0.046 * cell$sofa + log(0.75) * cell$tx
    )
    c(1, cumulative) - c(cumulative, 0)
  }
  simulated <- sim_markov(
    "violet",
    baseline_data = baseline,
    follow_up_time = 3,
    seed = 20261114
  )
  expect_transition_frequencies(
    simulated,
    covariates = c("age", "tx"),
    n_states = 4,
    absorbing = 4,
    probability = violet_probability,
    label = "VIOLET"
  )

  # The ACTT-2 preset, against the tables on its help page. On days 1 and 2
  # the nonlinear spline columns are 0, because the first knot is at day 2,
  # so the day enters as -0.3084 t.
  baseline <- expand.grid(rep = 1:1000, yprev = 4:7, tx = 0:1)
  baseline <- data.frame(
    id = seq_len(nrow(baseline)),
    yprev = baseline$yprev,
    tx = baseline$tx
  )
  actt2_probability <- function(cell) {
    # Cutoffs Pr(Y >= 2) to Pr(Y >= 8).
    alpha <- c(-0.6465, -4.9314, -5.2230, -7.4169, -11.5321, -15.2262, -22.0693)
    gamma <- c(0, 0, 6.2218, 7.1412, 9.9809, 14.1645, 19.3517)[cell$yprev]
    cumulative <- stats::plogis(
      alpha + gamma - 0.3084 * cell$time + 0.7 * cell$tx
    )
    c(1, cumulative) - c(cumulative, 0)
  }
  simulated <- sim_markov(
    "actt2",
    baseline_data = baseline,
    follow_up_time = 2,
    treatment_effect = 0.7,
    seed = 20261115
  )
  expect_transition_frequencies(
    simulated,
    covariates = "tx",
    n_states = 8,
    absorbing = 8,
    probability = actt2_probability,
    label = "ACTT-2"
  )
})
