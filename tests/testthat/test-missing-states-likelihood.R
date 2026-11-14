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
  # One stretch per run of hidden visits, of the run's length plus the
  # observed visit after it if there is one. The patient whose first rows
  # have no previous state has none.
  hidden <- attr(data, "hidden")
  n_stretches <- 0L
  n_steps <- 0L
  for (id in unique(hidden$id)) {
    times <- sort(hidden$time[hidden$id == id])
    runs <- split(times, cumsum(c(1, diff(times) != 1)))
    for (run in runs) {
      n_stretches <- n_stretches + 1L
      followed <- any(data$id == id & data$time == max(run) + 1)
      n_steps <- n_steps + length(run) + followed
    }
  }
  expect_identical(info$n_stretches, n_stretches)
  expect_identical(info$n_unobserved_visits, n_steps)
})

test_that("the fit maximises the reference log-likelihood", {
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

test_that("a gap at the end of follow-up counts as survival with `absorb`", {
  skip_on_cran() # Slow: refits; the reference check of the full fit runs on CRAN.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  trailing <- attr(data, "special")$trailing
  # Removing the trailing unobserved rows changes the fit: with death
  # absorbing and always recorded, they say the patient was alive.
  shortened <- data[!(data$id == trailing & is.na(data$y)), ]
  fit_short <- ms_fit(shortened)
  expect_gt(max(abs(coef(fit) - coef(fit_short))), 1e-6)
  reference_short <- ms_reference(shortened)
  expect_equal(
    attr(fit_short, "markov_path_fit")$loglik,
    reference_short(coef(fit_short)),
    tolerance = 1e-8
  )
  ms_expect_maximum(fit_short, reference_short)
  # The trailing rows add log P(not dead at both visits) to the patient's
  # contribution.
  per_patient <- ms_reference(data)(coef(fit), per_patient = TRUE)
  per_patient_short <- reference_short(coef(fit), per_patient = TRUE)
  expect_lt(per_patient[[as.character(trailing)]], per_patient_short[[as.character(trailing)]])
  others <- setdiff(names(per_patient), as.character(trailing))
  expect_equal(per_patient[others], per_patient_short[others], tolerance = 1e-12)
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
