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
  mostm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
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

test_that("default nulls are 0 for differences, 1 for ratios and none for SOPs and average times; an explicit null overrides them", {
  model <- make_null_test_model()
  difference <- make_null_test_comparison(model, "difference")

  # `null` is the expected null value; NULL means no test columns.
  cases <- list(
    list(
      label = "difference, default null",
      result = inferences(difference, method = "delta"),
      null = 0
    ),
    list(
      label = "ratio, default null",
      result = make_null_test_comparison(model, "ratio") |>
        inferences(method = "mvn", n_draws = 20, seed = 31),
      null = 1
    ),
    list(
      label = "difference, explicit null",
      result = inferences(difference, method = "delta", null = 0.5),
      null = 0.5
    ),
    list(
      label = "difference, null = NULL",
      result = inferences(difference, method = "delta", null = NULL),
      null = NULL
    ),
    list(
      label = "avg_sops(), default null",
      result = avg_sops(model, times = 1:3, absorb = 3) |>
        inferences(method = "delta"),
      null = NULL
    ),
    list(
      label = "avg_state_time(), default null",
      result = avg_state_time(model, times = 1:3, absorb = 3) |>
        inferences(method = "delta"),
      null = NULL
    )
  )

  for (case in cases) {
    result <- case$result
    if (is.null(case$null)) {
      expect_null(attr(result, "null"), label = case$label)
      expect_equal(
        intersect(c("statistic", "p.value", "s.value"), names(result)),
        character(),
        info = case$label
      )
      next
    }
    expect_identical(attr(result, "null"), case$null, label = case$label)
    expect_equal(
      result$statistic,
      (result$estimate - case$null) / result$std.error,
      info = case$label
    )
    expect_equal(
      result$p.value,
      2 * stats::pnorm(-abs(result$statistic)),
      info = case$label
    )
    expect_equal(result$s.value, -log2(result$p.value), info = case$label)
  }
})
