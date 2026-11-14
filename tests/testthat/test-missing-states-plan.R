# `markov_path_plan()`, which lays out the stretches of unobserved visits, and
# the invariance of the fit to the order of the rows.

# A plan for a small data frame with states 1..4 (4 absorbing unless
# `absorb` says otherwise) and the design (yprev = 2, yprev = 3, time, tx).
ms_plan <- function(data, absorb = 4L, lower = NULL, upper = NULL, weights = NULL) {
  ranges <- ms_row_ranges(data, absorb = absorb, lower = lower, upper = upper)
  # Rows without any record of the state are `NA` for the plan.
  y_missing <- if (is.null(lower)) is.na(data$y) else is.na(lower) & is.na(upper)
  ranges$lo[y_missing] <- NA
  ranges$hi[y_missing] <- NA
  mostr:::markov_path_plan(
    data = data,
    ranges = list(labels = as.character(1:4), lo = ranges$lo, hi = ranges$hi),
    absorb = absorb,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    weights = weights,
    design = function(rows) {
      previous <- as.integer(as.character(rows$yprev))
      cbind(previous == 2, previous == 3, rows$time, rows$tx) * 1
    }
  )
}

ms_two_patients <- function() {
  data.frame(
    id = rep(c(1L, 2L), each = 4),
    time = rep(1:4, 2),
    tx = 0,
    yprev = factor(c(2, 2, NA, NA, 2, 2, NA, NA), levels = 1:3),
    y = c(2L, NA, NA, 3L, 2L, NA, NA, 3L)
  )
}

test_that("identical stretches of two patients are stored once", {
  plan <- ms_plan(ms_two_patients())
  expect_identical(plan$stretch_n_steps, 3L)
  expect_identical(plan$stretch_weight, 2)
  expect_identical(plan$instances$id, c("1", "2"))
  expect_identical(plan$instances$stretch, c(1L, 1L))
  expect_identical(plan$n_steps_total, 6L)
  expect_identical(plan$step_lo, c(1L, 1L, 3L))
  expect_identical(plan$step_hi, c(3L, 3L, 3L))
  # Design rows: the first step's own row, and one row per possible previous
  # state (1, 2, 3; never the absorbing 4) at each later step.
  expect_identical(plan$step_rows[1, ], c(0L, 0L, 0L, 0L))
  expect_true(all(plan$step_rows[2:3, 1:3] > 0))
  expect_identical(plan$step_rows[2:3, 4], c(0L, 0L))
  x <- plan$x_gap
  expect_equal(x[plan$stretch_first_row, ], c(1, 0, 2, 0))
  expect_equal(x[plan$step_rows[3, 3], ], c(0, 1, 4, 0))
  # Different covariates make different stretches.
  different <- ms_two_patients()
  different$tx[different$id == 2] <- 1
  plan <- ms_plan(different)
  expect_identical(plan$stretch_n_steps, c(3L, 3L))
  expect_identical(plan$stretch_weight, c(1, 1))
  # Weights add up over identical stretches.
  plan <- ms_plan(ms_two_patients(), weights = rep(c(0.5, 2), each = 4))
  expect_identical(plan$stretch_weight, 2.5)
})

test_that("steps at the end of a stretch that allow every state are trimmed", {
  data <- ms_two_patients()
  data <- data[data$id == 1, ]
  # Without absorbing states the last two visits say nothing.
  trailing <- data
  trailing$y[4] <- NA
  plan <- ms_plan(trailing, absorb = integer())
  expect_identical(plan$n_steps_total, 0L)
  expect_identical(nrow(plan$instances), 0L)
  # With death absorbing they say nothing either: the patient may have been
  # lost to follow-up.
  plan <- ms_plan(trailing, absorb = 4L)
  expect_identical(plan$n_steps_total, 0L)
  # Known to be alive (`Ocens(1, 3)`), they say "not dead".
  plan <- ms_plan(
    trailing,
    absorb = 4L,
    lower = c(2, NA, 1, 1),
    upper = c(2, NA, 3, 3)
  )
  expect_identical(plan$n_steps_total, 3L)
  expect_identical(plan$step_hi, c(3L, 3L, 3L))
  # A missing state between recorded visits is alive, one after them is not
  # known: the stretch ends with the last range.
  plan <- ms_plan(
    trailing,
    absorb = 4L,
    lower = c(2, NA, 2, NA),
    upper = c(2, NA, 3, NA)
  )
  expect_identical(plan$n_steps_total, 2L)
  expect_identical(plan$step_lo, c(1L, 2L))
  expect_identical(plan$step_hi, c(3L, 3L))

  # A range before the trailing visits is kept, and only it.
  ranged <- trailing
  lower <- c(2, 2, NA, NA)
  upper <- c(2, 3, NA, NA)
  plan <- ms_plan(ranged, absorb = integer(), lower = lower, upper = upper)
  expect_identical(plan$n_steps_total, 1L)
  expect_identical(c(plan$step_lo, plan$step_hi), c(2L, 3L))
  # A full range written out is trimmed as well.
  plan <- ms_plan(
    ranged,
    absorb = integer(),
    lower = c(2, 2, 1, 1),
    upper = c(2, 3, 4, 4)
  )
  expect_identical(plan$n_steps_total, 1L)
})

test_that("rows of a patient are put in time order", {
  data <- ms_two_patients()
  shuffled <- data[c(4, 2, 7, 1, 8, 3, 6, 5), ]
  plan <- ms_plan(data)
  plan_shuffled <- ms_plan(shuffled)
  same <- c("x_gap", "stretch_first_row", "stretch_n_steps", "stretch_weight",
            "step_lo", "step_hi", "step_rows")
  expect_identical(plan_shuffled[same], plan[same])
  # The steps point at the same visits.
  expect_identical(
    shuffled$time[plan_shuffled$step_data_rows],
    data$time[plan$step_data_rows]
  )
})

test_that("factor time is ordered by its levels", {
  data <- ms_two_patients()
  data$numeric_time <- data$time
  # Level order differs from alphabetical order ("day 10" < "day 2").
  levels <- c("day 2", "day 3", "day 10", "day 11")
  data$day <- factor(levels[data$time], levels = levels)
  shuffled <- data[c(4, 2, 7, 1, 8, 3, 6, 5), ]
  plan_factor <- mostr:::markov_path_plan(
    data = shuffled,
    ranges = list(
      labels = as.character(1:4),
      lo = shuffled$y,
      hi = shuffled$y
    ),
    absorb = 4L,
    id_var = "id",
    time_var = "day",
    p_var = "yprev",
    weights = NULL,
    design = function(rows) {
      previous <- as.integer(as.character(rows$yprev))
      cbind(previous == 2, previous == 3, rows$numeric_time, rows$tx) * 1
    }
  )
  plan <- ms_plan(data)
  expect_identical(plan_factor$step_lo, plan$step_lo)
  expect_identical(plan_factor$x_gap, plan$x_gap)
  expect_identical(
    as.character(shuffled$day[plan_factor$step_data_rows]),
    rep(c("day 3", "day 10", "day 11"), 2)
  )
})

test_that("the fit does not depend on the order of rows or patients", {
  skip_on_cran() # Slow: two fits.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  # Patients in reverse order, rows interleaved between patients and in
  # reverse time order within each patient.
  withr::with_seed(7, shuffled <- data[sample(nrow(data)), ])
  fit_shuffled <- ms_fit(shuffled)
  expect_equal(coef(fit_shuffled), coef(fit), tolerance = 1e-10)
  expect_equal(
    attr(fit_shuffled, "markov_path_fit")$loglik,
    attr(fit, "markov_path_fit")$loglik,
    tolerance = 1e-10
  )
  reversed <- data[order(-data$id, -data$time), ]
  fit_reversed <- ms_fit(reversed)
  expect_equal(coef(fit_reversed), coef(fit), tolerance = 1e-10)
  # The patient-clustered covariance too; the patients are listed in another
  # order, so compare their scores by ID.
  expect_equal(fit_reversed$var, fit$var, tolerance = 1e-10)
  scores <- attr(fit, "markov_path_fit")$scores
  scores_reversed <- attr(fit_reversed, "markov_path_fit")$scores
  expect_equal(scores_reversed[rownames(scores), ], scores, tolerance = 1e-8)
})
