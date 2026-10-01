make_null_test_model <- function() {
  skip_if_not_installed("VGAM")
  withr::local_seed(4417)
  data <- make_markov_trajectories(
    n_patients = 80,
    follow_up_time = 4,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = 1:2,
    absorbing_state = 3,
    treatment_effect = 0.3,
    seed = 4417
  ) |>
    prepare_markov_data(absorbing_state = 3)
  vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id"
  )
}

make_null_test_comparison <- function(model, comparison) {
  avg_comparisons(
    model,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    state_sets = "1",
    comparison = comparison,
    times = 1:3,
    absorb = 3
  )
}

test_that("differences are tested against 0 by default", {
  model <- make_null_test_model()
  result <- make_null_test_comparison(model, "difference") |>
    inferences(method = "delta")

  expect_identical(attr(result, "null"), 0)
  expect_equal(result$statistic, result$estimate / result$std.error)
  expect_equal(
    result$p.value,
    2 * stats::pnorm(-abs(result$statistic))
  )
  expect_equal(result$s.value, -log2(result$p.value))
})

test_that("ratios are tested against 1 by default and reject a zero null", {
  model <- make_null_test_model()
  ratio <- make_null_test_comparison(model, "ratio")
  result <- inferences(ratio, method = "mvn", n_draws = 20, seed = 31)

  expect_identical(attr(result, "null"), 1)
  expect_equal(result$statistic, (result$estimate - 1) / result$std.error)
  expect_error(
    inferences(ratio, method = "mvn", n_draws = 20, seed = 31, null = 0),
    "Ratio comparisons cannot be tested against `null = 0`.",
    fixed = TRUE
  )
})

test_that("an explicit null overrides the default and NULL skips the test", {
  model <- make_null_test_model()
  difference <- make_null_test_comparison(model, "difference")

  shifted <- inferences(difference, method = "delta", null = 0.5)
  expect_identical(attr(shifted, "null"), 0.5)
  expect_equal(shifted$statistic, (shifted$estimate - 0.5) / shifted$std.error)

  untested <- inferences(difference, method = "delta", null = NULL)
  expect_null(attr(untested, "null"))
  expect_length(
    intersect(c("statistic", "p.value", "s.value"), names(untested)),
    0L
  )
})

test_that("SOPs and average times are not tested by default", {
  model <- make_null_test_model()
  avg <- avg_sops(model, times = 1:3, absorb = 3) |>
    inferences(method = "delta")
  total <- avg_time(model, times = 1:3, absorb = 3) |>
    inferences(method = "delta")

  for (result in list(avg, total)) {
    expect_null(attr(result, "null"))
    expect_length(
      intersect(c("statistic", "p.value", "s.value"), names(result)),
      0L
    )
  }
})

test_that("character nulls other than \"auto\" are rejected", {
  model <- make_null_test_model()
  difference <- make_null_test_comparison(model, "difference")

  expect_error(
    inferences(difference, method = "delta", null = "zero"),
    "`null` must be \"auto\", NULL, or a single finite numeric value.",
    fixed = TRUE
  )
})
