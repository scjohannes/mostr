make_null_test_model <- function() {
  skip_if_not_installed("VGAM")
  withr::local_seed(4417)
  data <- make_markov_trajectories(
    n_patients = 80,
    times = 1:4,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = 1:2,
    absorb = 3,
    treatment_effect = 0.3,
    seed = 4417
  ) |>
    prepare_markov_data(absorb = 3)
  mostm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    absorb = 3,
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
    times = 1:3
  )
}

test_that("default nulls are 0 for differences, 1 for ratios and none for SOPs and average times; an explicit null overrides them", {
  skip_on_cran()
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
      label = "ratio, default null, Wald interval",
      result = make_null_test_comparison(model, "ratio") |>
        inferences(method = "mvn", n_draws = 20, seed = 31, conf_type = "wald"),
      null = 1
    ),
    list(
      label = "ratio, default null, percentile interval",
      result = make_null_test_comparison(model, "ratio") |>
        inferences(method = "mvn", n_draws = 20, seed = 31),
      null = NULL
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
      result = avg_sops(model, times = 1:3) |>
        inferences(method = "delta"),
      null = NULL
    ),
    list(
      label = "avg_state_time(), default null",
      result = avg_state_time(model, times = 1:3) |>
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

test_that("s-values stay finite when the p-value underflows to 0", {
  # 2 * pnorm(-39) is below the smallest double, so it is 0 and
  # -log2(p) would be Inf. The s-value is computed from log(p) instead.
  x <- data.frame(estimate = c(1, 3.9, 40), std.error = 0.1)
  result <- add_null_test(x, null = 0)
  z <- c(10, 39, 400)
  expected <- -(log(2) + stats::pnorm(-z, log.p = TRUE)) / log(2)
  expect_equal(result$s.value, expected)
  expect_true(all(is.finite(result$s.value)))
  expect_equal(result$s.value[1], -log2(2 * stats::pnorm(-10)))
})

test_that("percentile intervals come without a test; Wald intervals keep one that agrees with them", {
  skip_on_cran()
  model <- make_null_test_model()
  difference <- make_null_test_comparison(model, "difference")

  # As in marginaleffects, a percentile interval is reported without a test:
  # a Wald test built from the draws' standard deviation can disagree with it.
  percentile <- expect_no_warning(
    inferences(difference, method = "mvn", n_draws = 20, seed = 31)
  )
  expect_null(attr(percentile, "null"))
  expect_equal(sop_test_cols(percentile), character())

  # A null the user asks for explicitly is not silently ignored.
  expect_warning(
    explicit <- inferences(
      difference,
      method = "mvn",
      n_draws = 20,
      seed = 31,
      null = 0.5
    ),
    "percentile"
  )
  expect_equal(sop_test_cols(explicit), character())

  # With `conf_type = "wald"`, the test and the interval rest on the same
  # standard error, so p < 0.05 exactly when the interval excludes the null.
  wald <- inferences(
    difference,
    method = "mvn",
    n_draws = 20,
    seed = 31,
    conf_type = "wald"
  )
  expect_identical(attr(wald, "null"), 0)
  excluded <- wald$conf.low > 0 | wald$conf.high < 0
  expect_identical(wald$p.value < 0.05, excluded)
})

test_that("logit intervals are tested on the logit scale", {
  # The estimate's logit has standard error SE / (p (1 - p)) by the delta
  # method; the logit interval is built from the same quantity, so a test on
  # that scale excludes the null exactly when the interval does.
  x <- data.frame(
    estimate = c(0.0051, 0.03, 0.2, 0.6),
    std.error = c(0.0013, 0.007, 0.04, 0.05)
  )
  attr(x, "conf_type") <- "logit"
  bounds <- suppressWarnings(
    delta_interval_bounds(x$estimate, x$std.error, 0.95, "logit")
  )
  x$conf.low <- bounds$conf.low
  x$conf.high <- bounds$conf.high

  for (null in c(0.008, 0.02, 0.05, 0.3, 0.5)) {
    result <- add_null_test(x, null)
    logit_se <- x$std.error / (x$estimate * (1 - x$estimate))
    expect_equal(
      result$statistic,
      (stats::qlogis(x$estimate) - stats::qlogis(null)) / logit_se,
      info = paste("null =", null)
    )
    expect_equal(result$p.value, 2 * stats::pnorm(-abs(result$statistic)))
    excluded <- x$conf.low > null | x$conf.high < null
    expect_identical(result$p.value < 0.05, excluded, info = paste("null =", null))
  }
  # P(state 8 at day 2) = 0.0051, interval 0.0031 to 0.0084: a null of 0.008
  # lies inside the interval, so it must not be rejected. On the probability
  # scale, (0.0051 - 0.008) / 0.0013 = -2.23 rejected it.
  expect_gt(add_null_test(x[1, ], 0.008)$p.value, 0.05)

  # The logit of 0 or 1 is infinite, so such a null can't be tested there.
  expect_error(add_null_test(x, 0), class = "mostr_error")
  expect_error(add_null_test(x, 1), class = "mostr_error")

  # An estimate on the boundary has no logit interval and no test.
  boundary <- x[1, ]
  boundary$estimate <- 0
  expect_true(is.na(add_null_test(boundary, 0.5)$statistic))
})

test_that("delta SOPs with the logit interval get a test that agrees with it", {
  skip_on_cran()
  model <- make_null_test_model()
  result <- avg_sops(model, times = 1:3) |>
    inferences(method = "delta", null = 0.2)
  expect_identical(attr(result, "conf_type"), "logit")
  tested <- !is.na(result$p.value)
  expect_true(any(tested))
  excluded <- result$conf.low > 0.2 | result$conf.high < 0.2
  expect_identical(result$p.value[tested] < 0.05, excluded[tested])
})

test_that("ratio draws keep x / 0 as Inf and report 0 / 0 draws they leave out", {
  withr::local_seed(1)
  B <- 200
  lo <- stats::rlnorm(B, log(0.01), 0.5)
  hi <- stats::rlnorm(B, log(0.02), 0.5)
  lo[1:20] <- 0 # 10% of draws have a reference probability of 0 ...
  hi[1:10] <- 0 # ... half of them a comparison probability of 0 too.
  draws <- rbind(
    data.frame(draw_id = seq_len(B), time = 1, tx = 0, estimate = lo),
    data.frame(draw_id = seq_len(B), time = 1, tx = 1, estimate = hi)
  )
  ratio <- compare_counterfactual_levels(
    draws,
    "estimate",
    c("draw_id", "time", "tx"),
    list(tx = c(0, 1)),
    "ratio"
  )
  # x / 0 with x > 0 is +Inf, and 0 / 0 is undefined (NaN).
  expect_equal(sum(is.infinite(ratio$estimate)), 10)
  expect_equal(sum(is.nan(ratio$estimate)), 10)

  # 10 of the remaining 190 draws (5.3%) are Inf, so the 97.5% quantile is
  # Inf; dropping them would report a finite upper limit. The 0 / 0 draws
  # can't be ordered, so they are left out, with a warning giving their count.
  expect_warning(
    ci <- compute_ci_from_draws(ratio, c("time", "contrast")),
    "10 of 200"
  )
  expect_identical(ci$conf.high, Inf)
  expect_true(is.finite(ci$conf.low))
  # A standard deviation that includes infinite draws is not finite.
  expect_true(is.na(ci$std.error))
})
