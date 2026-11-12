# `variables` may name a column that is in the data but that the fitted model
# doesn't use, as in marginaleffects. Setting such a column can't change any
# prediction, so a comparison of it is exactly 0. Its uncertainty is not
# reported: standard error, test statistic, p-value and interval limits are NA,
# as marginaleffects reports them, instead of an SE of 0 and a NaN statistic.
# Predictions under such a variable (avg_sops()) are ordinary predictions and
# keep their standard errors. A column that is not in the data at all is still
# an error.

unused_profiles <- function(case) {
  profiles <- case$profiles
  profiles$note <- profiles$id %% 2
  profiles
}

uncertainty_columns <- c("std.error", "statistic", "p.value", "conf.low", "conf.high")

expect_unreported_uncertainty <- function(x, info) {
  for (column in intersect(uncertainty_columns, names(x))) {
    expect_true(all(is.na(x[[column]])), info = paste(info, column))
  }
  expect_true("std.error" %in% names(x), info = info)
}

test_that("a comparison of a variable the model doesn't use is 0, without a message", {
  case <- make_comparisons_case()
  profiles <- unused_profiles(case)
  base <- list(
    case$model,
    newdata = profiles,
    variables = list(note = c(0, 1)),
    times = 1:4,
    absorb = 3
  )
  # The averaged comparisons are 0 in every grid set-up; this checks the
  # comparison of each profile.
  per_profile <- do.call(comparisons, c(base, list(state_sets = "1")))
  expect_equal(per_profile$estimate, rep(0, nrow(per_profile)))
  expect_no_warning(do.call(avg_comparisons, c(base, list(state_sets = "1"))))
  expect_no_message(do.call(avg_comparisons, c(base, list(state_sets = "1"))))
})

test_that("its standard error, test and interval are NA for every method", {
  skip_on_cran()
  case <- make_comparisons_case()
  profiles <- unused_profiles(case)
  contrast <- avg_comparisons(
    case$model,
    newdata = profiles,
    variables = list(note = c(0, 1)),
    times = 1:4,
    absorb = 3,
    state_sets = "1"
  )
  expect_unreported_uncertainty(
    inferences(contrast, method = "delta", vcov = "conditional"),
    "delta"
  )
  expect_unreported_uncertainty(
    inferences(contrast, method = "mvn", n_draws = 20, seed = 1),
    "mvn"
  )
  per_profile <- comparisons(
    case$model,
    newdata = profiles,
    variables = list(note = c(0, 1)),
    times = 1:4,
    absorb = 3,
    state_sets = "1"
  )
  expect_unreported_uncertainty(
    inferences(per_profile, method = "delta", vcov = "conditional"),
    "per-profile delta"
  )
  time_benefit <- avg_comparisons(
    case$model,
    variables = list(marker = c(0, 1)),
    times = 1:4,
    absorb = 3,
    estimand = "time_benefit"
  )
  # `marker` is in this model, so this is an ordinary comparison with an SE.
  expect_true(all(is.finite(
    inferences(time_benefit, method = "mvn", n_draws = 20, seed = 1)$std.error
  )))
  skip_on_cran()
  stored <- avg_comparisons(
    case$model,
    variables = list(note = c(0, 1)),
    newdata = profiles,
    times = 1:4,
    absorb = 3,
    state_sets = "1"
  )
  expect_unreported_uncertainty(
    inferences(stored, method = "bootstrap", n_draws = 3, seed = 1),
    "bootstrap"
  )
})

test_that("predictions under a variable the model doesn't use keep their standard errors", {
  case <- make_comparisons_case()
  profiles <- unused_profiles(case)
  average <- avg_sops(
    case$model,
    newdata = profiles,
    variables = list(note = c(0, 1)),
    times = 1:4,
    absorb = 3
  )
  # That the estimates are equal for every value is checked in every grid
  # set-up; this checks that they keep finite standard errors.
  inferred <- inferences(average, method = "delta", vcov = "conditional")
  expect_true(all(is.finite(inferred$std.error)))
})
