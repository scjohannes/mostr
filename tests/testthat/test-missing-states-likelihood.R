# Fits with unobserved states against the plain-R reference in
# helper-missing-states.R, which reads the data frame itself and lists every
# path of states. These compare two routes to the same numbers, so they pass
# for any simulated data on which the fit converges.

test_that("the fitted log-likelihood equals the path enumeration of the data", {
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  info <- attr(fit, "markov_path_fit")
  expect_false(is.null(info))
  reference <- ms_reference(data)
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  # The deviance stored on the fit is the same likelihood.
  expect_equal(fit$deviance[[2]], -2 * reference(coef(fit)), tolerance = 1e-8)
  # One stretch per run of hidden visits followed by an observed visit, of
  # the run's length plus that visit. A run after the patient's last
  # observed visit says nothing and leaves no stretch. The patient whose
  # first rows have no previous state has none.
  hidden <- attr(data, "hidden")
  n_stretches <- 0L
  n_steps <- 0L
  for (id in unique(hidden$id)) {
    times <- sort(hidden$time[hidden$id == id])
    runs <- split(times, cumsum(c(1, diff(times) != 1)))
    for (run in runs) {
      followed <- any(data$id == id & data$time == max(run) + 1)
      if (!followed) next
      n_stretches <- n_stretches + 1L
      n_steps <- n_steps + length(run) + 1L
    }
  }
  expect_identical(info$n_stretches, n_stretches)
  expect_identical(info$n_unobserved_visits, n_steps)
})

test_that("the fit maximises the reference log-likelihood", {
  skip_on_cran() # Slow: numerical gradients; the log-likelihood check runs on CRAN.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  reference <- ms_reference(data)
  ms_expect_maximum(fit, reference)
  # The score stored on the fit is the reference's gradient.
  expect_equal(
    unname(fit$u),
    unname(central_difference_jacobian(reference, coef(fit))[1, ]),
    tolerance = 1e-4
  )
  expect_lt(max(abs(fit$u)), 1e-3)
})

test_that("a general optimiser does not improve on the fit", {
  skip_on_cran() # Slow: the optimiser evaluates the reference many times.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  reference <- ms_reference(data)
  start <- coef(fit)
  # Trial steps that cross the intercepts give negative probabilities and
  # NaN log-likelihoods; BFGS steps back from them.
  better <- suppressWarnings(stats::optim(
    start,
    reference,
    method = "BFGS",
    control = list(fnscale = -1, reltol = 1e-14, maxit = 200)
  ))
  expect_lt(better$value - reference(start), 1e-6)
  # From a distant start the optimiser reaches the same maximum.
  far <- suppressWarnings(stats::optim(
    start + 0.3,
    reference,
    method = "BFGS",
    control = list(fnscale = -1, reltol = 1e-14, maxit = 500)
  ))
  expect_lt(abs(far$value - reference(start)), 1e-6)
  expect_equal(unname(far$par), unname(start), tolerance = 1e-3)
})

test_that("the information is the negative Hessian of the reference", {
  skip_on_cran() # Slow: about 100 evaluations of the reference.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  reference <- ms_reference(data)
  hessian <- ms_numeric_hessian(reference, coef(fit))
  information <- attr(fit, "markov_path_fit")$information
  expect_equal(unname(information), -hessian, tolerance = 1e-4)
  # `orig.var` is its inverse.
  expect_equal(
    unname(fit$orig.var),
    solve(unname(information)),
    tolerance = 1e-8
  )
})

test_that("a gap before death excludes death during the gap", {
  skip_on_cran() # Slow: the reference per patient; the full fit is checked on CRAN.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  dead <- attr(data, "special")$dead_late[1]
  rows <- data$id == dead
  # The visit before death is unobserved and death follows it.
  expect_true(is.na(data$y[rows][sum(rows) - 1]))
  expect_equal(as.integer(as.character(data$y[rows][sum(rows)])), 4L)

  one <- data[rows, ]
  theta <- coef(fit)
  per_patient <- ms_reference(data)(theta, per_patient = TRUE)
  # The reference written out for this patient: the gap's state is 1, 2 or 3.
  alpha <- theta[1:3]
  beta <- theta[4:7]
  transition <- function(previous, time, tx) {
    x <- c(previous == 2, previous == 3, time, tx)
    cumulative <- c(1, stats::plogis(alpha + sum(x * beta)), 0)
    cumulative[1:4] - cumulative[2:5]
  }
  n <- nrow(one)
  states <- c(as.integer(as.character(one$yprev[1])), as.integer(as.character(one$y)))
  known_part <- 0
  for (t in seq_len(n - 2)) {
    known_part <- known_part +
      log(transition(states[t], one$time[t], one$tx[t])[states[t + 1]])
  }
  gap <- sum(vapply(1:3, function(z) {
    transition(states[n - 1], one$time[n - 1], one$tx[n - 1])[z] *
      transition(z, one$time[n], one$tx[n])[4]
  }, numeric(1)))
  expect_equal(
    per_patient[[as.character(dead)]],
    unname(known_part + log(gap)),
    tolerance = 1e-12
  )

  # Allowing death at the unobserved visit gives a different likelihood.
  with_death <- ms_reference_loglik_death_allowed(theta, data)
  expect_gt(abs(with_death - ms_reference(data)(theta)), 1e-3)
  expect_gt(abs(with_death - attr(fit, "markov_path_fit")$loglik), 1e-3)
})

test_that("missing states after the last observed visit say nothing", {
  skip_on_cran() # Slow: three fits; the plan test of the rule runs on CRAN.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  trailing <- attr(data, "special")$trailing
  # The patient may have been lost to follow-up with the rows left in the
  # data, so the rows say nothing, also with `absorb`: the fit is the one
  # without them.
  shortened <- data[!(data$id == trailing & is.na(data$y)), ]
  fit_short <- ms_fit(shortened)
  expect_equal(coef(fit), coef(fit_short), tolerance = 1e-10)
  expect_equal(fit$var, fit_short$var, tolerance = 1e-10)
  expect_equal(
    attr(fit, "markov_path_fit")$loglik,
    attr(fit_short, "markov_path_fit")$loglik,
    tolerance = 1e-10
  )
  # The reference sums them over every state, which adds a factor of 1.
  per_patient <- ms_reference(data)(coef(fit), per_patient = TRUE)
  per_patient_short <- ms_reference(shortened)(coef(fit), per_patient = TRUE)
  expect_equal(per_patient, per_patient_short, tolerance = 1e-12)

  # Known to be alive at those visits: `Ocens(1, 3)` adds log P(not dead at
  # both visits) to the patient's contribution and changes the fit.
  alive <- data
  alive$lower <- alive$upper <- as.integer(as.character(alive$y))
  rows <- alive$id == trailing & is.na(alive$y)
  alive$lower[rows] <- 1
  alive$upper[rows] <- 3
  fit_alive <- ms_fit(alive, formula = Ocens(lower, upper) ~ yprev + time + tx)
  expect_gt(max(abs(coef(fit_alive) - coef(fit))), 1e-6)
  reference_alive <- ms_reference(
    alive,
    ranges = ms_row_ranges(alive, lower = alive$lower, upper = alive$upper)
  )
  expect_equal(
    attr(fit_alive, "markov_path_fit")$loglik,
    reference_alive(coef(fit_alive)),
    tolerance = 1e-8
  )
  ms_expect_maximum(fit_alive, reference_alive)
  per_patient_alive <- reference_alive(coef(fit), per_patient = TRUE)
  id <- as.character(trailing)
  expect_lt(per_patient_alive[[id]], per_patient[[id]])
})

test_that("a patient without any recorded state is not fitted", {
  skip_on_cran() # Slow: two fits.
  data <- ms_fixture("downstream")
  unseen <- data[data$id == attr(data, "special")$alive_long[7], ]
  unseen$id <- max(data$id) + 1
  unseen$y <- NA
  unseen$yprev[-1] <- NA
  with_unseen <- rbind(data, unseen)
  expect_warning(fit <- ms_fit(with_unseen), "1 patient has")
  fit_without <- ms_fit(data)
  expect_equal(coef(fit), coef(fit_without), tolerance = 1e-10)
  expect_false(as.character(unseen$id[1]) %in% attr(fit, "markov_path_fit")$ids)
})

test_that("without `absorb` a gap at the end of follow-up contributes nothing", {
  skip_on_cran() # Slow: two fits.
  data <- ms_fixture("complete")
  # Without `absorb` a missing state can be any of the four, including 4.
  # Only trailing visits are hidden, so no design row with previous state 4
  # (which `yprev` has no level for) is needed.
  last_time <- tapply(data$time, data$id, max)
  ids <- as.integer(names(last_time)[last_time == 6])[1:3]
  hidden <- data
  for (j in ids) {
    hidden$y[hidden$id == j & hidden$time >= 5] <- NA
    hidden$yprev[hidden$id == j & hidden$time == 6] <- NA
  }
  # State 4 is never a previous state, and without `absorb` mostm() warns.
  fit_quietly <- function(data) {
    suppressWarnings(
      mostm(ordered(y) ~ yprev + time + tx, data = data, id_var = "id")
    )
  }
  fit <- fit_quietly(hidden)
  # With no stretch left, the fit to the observed transitions is the maximum.
  expect_null(attr(fit, "markov_path_fit"))

  removed <- hidden[!is.na(hidden$y), ]
  fit_removed <- fit_quietly(removed)
  expect_identical(ms_numbers(fit), ms_numbers(fit_removed))
  # And it maximises the reference, in which the trailing visits sum to 1.
  reference <- ms_reference(hidden, absorb = integer())
  expect_equal(
    reference(coef(fit)),
    ms_reference(removed, absorb = integer())(coef(fit)),
    tolerance = 1e-12
  )
  ms_expect_maximum(fit, reference)
})

test_that("each kind of stretch alone matches the reference", {
  skip_on_cran() # Slow: seven fits; the fit with all kinds runs on CRAN.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  special <- attr(data, "special")
  theta <- coef(fit)
  reference <- ms_reference(data)(theta, per_patient = TRUE)
  # One fit per kind of stretch, on the patients without unobserved visits
  # plus the one patient with that kind, so an error in one kind can't be
  # hidden by the others.
  kinds <- list(
    "gap of length 1" = special$alive_long[1],
    "gap of length 2" = special$alive_long[2],
    "gap of length 3" = special$alive_long[3],
    "two gaps" = special$alive_long[4],
    "gap to the last row" = special$trailing,
    "gap before death" = special$dead_late[1],
    "unknown first states" = special$unknown_start
  )
  complete_ids <- setdiff(unique(data$id), unlist(attr(data, "hidden")$id))
  complete_ids <- setdiff(complete_ids, special$unknown_start)
  for (kind in names(kinds)) {
    subset_data <- data[data$id %in% c(kinds[[kind]], complete_ids), ]
    sub_fit <- ms_fit(subset_data)
    sub_reference <- ms_reference(subset_data)
    # The patient with unknown first states leaves no stretch, so that fit
    # is an ordinary one; its deviance is the same likelihood.
    expect_equal(
      ms_loglik(sub_fit),
      sub_reference(coef(sub_fit)),
      tolerance = 1e-8,
      info = kind
    )
    ms_expect_maximum(sub_fit, sub_reference, info = kind)
  }
  # The patient with unknown first states contributes only from the first
  # observed state on: their first two rows add nothing.
  first <- special$unknown_start
  rows <- data[data$id == first, ]
  expect_true(is.na(rows$y[1]) && is.na(rows$yprev[1]) && is.na(rows$yprev[2]))
  from_third <- rows[-(1:2), ]
  alone <- ms_reference(from_third)(theta)
  expect_equal(reference[[as.character(first)]], alone, tolerance = 1e-12)
})

test_that("a formula without the previous state counts each row once", {
  skip_on_cran() # Slow: a fit; the two-state test of the same fix runs on CRAN.
  # A row with an observed state after an unobserved visit ends a stretch; it
  # must not also be fitted as an observed transition, which happened when
  # the formula did not contain `yprev` and so the engine kept the row.
  data <- ms_fixture("gaps")
  fit <- ms_fit(data, formula = ordered(y) ~ time + tx)
  info <- attr(fit, "markov_path_fit")
  reference <- ms_reference(data, previous_effect = FALSE)
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  ms_expect_maximum(fit, reference)
  # Every visit with an observed state and a recorded previous state is an
  # observed transition, every other one in a stretch is a step.
  observed <- !is.na(data$y) & !is.na(data$yprev)
  expect_identical(info$n_observed_transitions, sum(observed))
})

# Two states, the second absorbing: an unobserved visit before a later
# recorded one can only have been in state 1. 40 patients over 5 days, death
# on a fixed day for some of them.
ms_two_state_data <- function() {
  rows <- lapply(1:40, function(i) {
    death <- (i * 7) %% 9 + 1
    times <- seq_len(min(5, death))
    data.frame(
      id = i,
      time = times,
      tx = i %% 2,
      yprev = 1,
      y = ifelse(times == death, 2L, 1L)
    )
  })
  data <- do.call(rbind, rows)
  data$yprev <- factor(data$yprev, levels = 1)
  data
}

test_that("an unobserved visit with one possible state is summed over", {
  skip_on_cran() # Slow: two fits; the plan and kernel tests run on CRAN.
  complete <- ms_two_state_data()
  hidden <- complete
  # Interior visits of every third patient, followed by a recorded visit.
  hide <- hidden$id %% 3 == 0 & hidden$time %in% 2:3 &
    hidden$time < stats::ave(hidden$time, hidden$id, FUN = max)
  hidden$y[hide] <- NA
  hidden$yprev[which(hide) + 1L] <- NA
  expect_gt(sum(hide), 10L)
  fit_two <- function(data) {
    mostm(ordered(y) ~ time + tx, data = data, id_var = "id", absorb = 2)
  }
  fit <- fit_two(hidden)
  info <- attr(fit, "markov_path_fit")
  expect_false(is.null(info))
  reference <- ms_reference(hidden, K = 2, absorb = 2, previous_effect = FALSE)
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  ms_expect_maximum(fit, reference)
  # A hidden visit that can only be state 1 is state 1: the fit is the one to
  # the complete data (up to the two fits' convergence tolerances).
  fit_complete <- fit_two(complete)
  expect_equal(coef(fit), coef(fit_complete), tolerance = 1e-6)
  expect_equal(info$loglik, ms_loglik(fit_complete), tolerance = 1e-10)
})

test_that("a stretch with a range and missing states after it matches the reference", {
  skip_on_cran() # Slow: a fit and the reference.
  data <- ms_fixture("complete")
  data$lower <- data$upper <- as.integer(as.character(data$y))
  # A patient with observed states, an unobserved visit, a range and then
  # unobserved visits to the end.
  patient <- attr(ms_fixture("gaps"), "special")$alive_long[7]
  rows <- which(data$id == patient)
  expect_length(rows, 6L)
  data$lower[rows[3]] <- data$upper[rows[3]] <- NA
  data$lower[rows[4]] <- 2
  data$upper[rows[4]] <- 3
  data$lower[rows[5:6]] <- data$upper[rows[5:6]] <- NA
  data$yprev[rows[4:6]] <- NA
  fit <- ms_fit(data, formula = Ocens(lower, upper) ~ yprev + time + tx)
  info <- attr(fit, "markov_path_fit")
  # One stretch: the unobserved visit (states 1 to 3, because a range
  # follows) and the range; the visits after it are trimmed.
  expect_identical(info$n_stretches, 1L)
  expect_identical(info$n_unobserved_visits, 2L)
  reference <- ms_reference(
    data,
    ranges = ms_row_ranges(data, lower = data$lower, upper = data$upper)
  )
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  ms_expect_maximum(fit, reference)
})
