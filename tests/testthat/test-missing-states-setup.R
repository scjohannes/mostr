# The model structure (design columns, knots, factor levels, states) and the
# starting values of a fit with unobserved states come from every visit with
# an exactly observed state, with an unknown previous state filled in for the
# set-up only. A state therefore has to be observed exactly at some visit, not
# at a visit whose previous state is also observed. These tests compare two
# routes to the same numbers (the fit and the reference that lists every path
# of the data frame), so they pass for any simulated data on which the fit
# converges.

test_that("the set-up carries the last known state forward within a patient", {
  # Rows out of order; patient 2's first previous state is unknown.
  data <- data.frame(
    id = c(1, 1, 1, 1, 1, 2, 2, 2),
    time = c(2, 1, 3, 4, 5, 1, 2, 3),
    yprev = factor(c("2", "1", NA, NA, "3", NA, NA, "2"), levels = 1:3),
    y = c(NA, 2L, NA, 3L, 3L, NA, 2L, 1L)
  )
  ranges <- markov_state_ranges(ordered(y) ~ 1, data, environment())
  setup <- markov_setup_rows(data, ranges, seq_len(nrow(data)), "id", "time", "yprev")
  # Patient 1: day 1 from state 1 to 2, day 2 unobserved, day 3 unobserved,
  # day 4 state 3 (last known state 2, from day 1), day 5 from 3 to 3.
  expect_identical(
    as.character(setup$previous),
    c("2", "1", "2", "2", "3", NA, NA, "2")
  )
  expect_identical(setup$recorded, !is.na(data$yprev))
  # Rows with an exact state and a known or carried-forward previous state.
  expect_identical(setup$use, c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
  # `rows` (the rows `subset` keeps) limits the rows used.
  setup_subset <- markov_setup_rows(data, ranges, c(1:5), "id", "time", "yprev")
  expect_identical(setup_subset$use, c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("a state recorded only after unobserved visits is fitted", {
  data <- ms_deaths_after_gaps()
  # No transition with both states observed ends in death.
  expect_false(any(data$y == 4 & !is.na(data$yprev), na.rm = TRUE))
  fit <- ms_fit(data)
  info <- attr(fit, "markov_path_fit")
  expect_false(is.null(info))
  reference <- ms_reference(data)
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  expect_identical(fit$yunique, c("1", "2", "3", "4"))

  # The stored rows are still the fully observed transitions only.
  observed <- !is.na(data$y) & !is.na(data$yprev)
  expect_identical(nrow(fit$x), sum(observed))
  expect_identical(info$n_observed_transitions, sum(observed))
  expect_identical(as.character(fit$y), as.character(data$y[observed]))
  expect_identical(length(fit$linear.predictors), sum(observed))
  expect_equal(fit$stats[["Obs"]], sum(observed))
  expect_identical(nrow(attr(fit, "markov_data")), sum(observed))
  expect_identical(
    unname(fit$linear.predictors),
    unname(coef(fit)[[1]] + as.numeric(fit$x %*% coef(fit)[-(1:3)]))
  )

  skip_on_cran() # Slow: numerical gradients and predictions.
  ms_expect_maximum(fit, reference)
  average <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:6)
  totals <- stats::aggregate(estimate ~ tx + time, average, sum)
  expect_equal(totals$estimate, rep(1, nrow(totals)), tolerance = 1e-10)
  expect_true(all(average$estimate[average$state == 4 & average$time == 6] > 0))
  delta <- inferences(average, method = "delta")
  expect_true(all(is.finite(delta$std.error)))
})

test_that("a previous state absent from fully observed transitions is fitted", {
  data <- ms_previous_state_before_gaps()
  observed <- !is.na(data$y) & !is.na(data$yprev)
  expect_false(any(data$yprev[observed] == "3"))
  expect_true(any(data$y == 3, na.rm = TRUE))
  # Patients who died on their first row from state 3 have no recorded state
  # left and are named as not part of the fit.
  fit <- suppressWarnings(
    ms_fit(data),
    classes = "mostr_warning_dropped_patients"
  )
  reference <- ms_reference(data)
  expect_equal(
    attr(fit, "markov_path_fit")$loglik,
    reference(coef(fit)),
    tolerance = 1e-8
  )
  expect_true("yprev=3" %in% names(coef(fit)))
  expect_identical(nrow(fit$x), sum(observed))
  skip_on_cran() # Slow: numerical gradients.
  ms_expect_maximum(fit, reference)
})

test_that("a state seen only once, after a gap, is fitted (integer and ordered response)", {
  data <- ms_state_only_after_gap()
  reference <- ms_reference(data)
  for (formula in list(y ~ yprev + time + tx, ordered(y) ~ yprev + time + tx)) {
    fit <- ms_fit(data, formula = formula)
    expect_identical(as.character(fit$yunique), c("1", "2", "3", "4"))
    expect_false("1" %in% as.character(fit$y))
    expect_equal(
      attr(fit, "markov_path_fit")$loglik,
      reference(coef(fit)),
      tolerance = 1e-8
    )
  }
})

test_that("spline knots come from the visits with an exactly observed state", {
  data <- ms_fixture("gaps")
  fit <- ms_fit(data, formula = ordered(y) ~ yprev + rms::rcs(time, 3) + tx)
  ranges <- markov_state_ranges(ordered(y) ~ 1, data, environment())
  setup <- markov_setup_rows(data, ranges, seq_len(nrow(data)), "id", "time", "yprev")
  expect_equal(
    fit$Design$parms$time,
    Hmisc::rcspline.eval(data$time[setup$use], nk = 3, knots.only = TRUE)
  )
  # Without missing states the knots are those of the ordinary fit.
  complete <- ms_fit(
    ms_fixture("complete"),
    formula = ordered(y) ~ yprev + rms::rcs(time, 3) + tx
  )
  expect_null(attr(complete, "markov_path_fit"))
})

test_that("a state that is never observed exactly is refused", {
  data <- ms_fixture("complete")
  data$y[data$y == 1] <- 2L
  data$yprev[!is.na(data$yprev) & data$yprev == "1"] <- "2"
  data$yprev <- droplevels(data$yprev)
  data$lower <- data$y
  data$upper <- data$y
  # A last visit known to be in state 1 or 2.
  data$lower[data$id == 1 & data$time == 6] <- 1L
  expect_error(
    ms_fit(data, formula = Ocens(lower, upper) ~ yprev + time + tx),
    paste0(
      "The state 1 is never observed exactly at a visit after the patient's ",
      "first known state, so the model can't be set up."
    ),
    class = "mostr_error",
    fixed = TRUE
  )
})

test_that("a previous state seen only before a range is refused", {
  # State 3 occurs once, at a patient's last recorded visit, followed by a
  # visit known only to be in states 1 to 3. The likelihood needs the effect
  # of state 3 as previous state, and no set-up row has it.
  data <- ms_fixture("complete")
  data$y[data$y == 3] <- 2L
  data$yprev[!is.na(data$yprev) & data$yprev == "3"] <- "2"
  alive <- setdiff(unique(data$id), unique(data$id[data$y == 4]))
  last <- max(which(data$id == alive[1]))
  data$y[last] <- 3L
  extra <- data[last, ]
  extra$time <- extra$time + 1
  extra$yprev <- factor("3", levels = 1:3)
  extra$y <- NA
  data <- rbind(data, extra)
  data$lower <- ifelse(is.na(data$y), 1L, data$y)
  data$upper <- ifelse(is.na(data$y), 3L, data$y)
  expect_error(
    suppressWarnings(
      ms_fit(data, formula = Ocens(lower, upper) ~ yprev + time + tx)
    ),
    paste0(
      "The state 3 is never the previous state (`yprev`) of a visit with an ",
      "exactly observed state, although the model needs its effect"
    ),
    class = "mostr_error",
    fixed = TRUE
  )
  # Without the previous state in the model nothing is missing.
  fit <- suppressWarnings(
    ms_fit(data, formula = Ocens(lower, upper) ~ time + tx)
  )
  expect_false(is.null(attr(fit, "markov_path_fit")))
})
