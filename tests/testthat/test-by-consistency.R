# `by` follows the marginaleffects package everywhere: `FALSE` or `NULL` gives
# one row per profile, `TRUE` averages over all profiles, and column names
# average within groups. Each per-profile function hands averaging to its
# `avg_*()` partner, and each `avg_*()` function hands `by = FALSE` back.

test_that("sops() returns each profile's SOPs under each value of `variables`", {
  case <- make_comparisons_case()
  n <- nrow(case$profiles)
  cf <- sops(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3
  )
  expect_s3_class(cf, "markov_sops")
  expect_setequal(cf$rowid, seq_len(2L * n))
  expect_setequal(cf$rowidcf, seq_len(n))
  for (value in c(0, 1)) {
    expected <- sops_with_value(case, "tx", value, times = 1:3)
    rows <- cf[cf$tx == value, , drop = FALSE]
    expect_equal(rows$rowidcf, expected$rowid)
    expect_equal(rows$estimate, expected$estimate)
  }

  averaged <- avg_sops(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3
  )
  means <- stats::aggregate(estimate ~ time + state + tx, cf, mean)
  index <- match(
    paste(averaged$time, averaged$state, averaged$tx),
    paste(means$time, means$state, means$tx)
  )
  expect_equal(averaged$estimate, means$estimate[index], tolerance = 1e-12)

  # Without `variables`, the rows are the profiles themselves.
  plain <- sops(case$model, times = 1:3, absorb = 3)
  expect_false("rowidcf" %in% names(plain))
})

test_that("sops() and avg_sops() hand `by` to each other", {
  case <- make_comparisons_case()
  args <- list(case$model, times = 1:3, absorb = 3)
  cf_args <- c(args, list(variables = list(tx = c(0, 1))))

  averaged <- do.call(avg_sops, args)
  expect_s3_class(averaged, "markov_avg_sops")
  expect_equal(do.call(avg_sops, c(args, list(by = TRUE))), averaged)
  expect_equal(do.call(sops, c(args, list(by = TRUE))), averaged)
  expect_equal(
    do.call(sops, c(args, list(by = "marker"))),
    do.call(avg_sops, c(args, list(by = "marker")))
  )
  expect_equal(
    do.call(sops, c(cf_args, list(by = TRUE))),
    do.call(avg_sops, cf_args)
  )

  per_profile <- do.call(sops, args)
  expect_equal(do.call(avg_sops, c(args, list(by = FALSE))), per_profile)
  expect_equal(do.call(avg_sops, c(args, list(by = NULL))), per_profile)
  expect_equal(do.call(sops, c(args, list(by = NULL))), per_profile)
  expect_equal(
    do.call(avg_sops, c(cf_args, list(by = FALSE))),
    do.call(sops, cf_args)
  )

  expect_error(do.call(sops, c(args, list(by = 1))), "`by` must be")
  expect_error(do.call(avg_sops, c(args, list(by = NA))), "`by` must be")
})

test_that("avg_state_time() and state_time() hand `by` to each other", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:3,
    absorb = 3
  )
  averaged <- do.call(avg_state_time, args)
  expect_s3_class(averaged, "markov_avg_state_time")
  expect_equal(do.call(avg_state_time, c(args, list(by = TRUE))), averaged)
  expect_equal(do.call(state_time, c(args, list(by = TRUE))), averaged)
  expect_equal(
    do.call(state_time, c(args, list(by = "marker"))),
    do.call(avg_state_time, c(args, list(by = "marker")))
  )

  per_profile <- do.call(state_time, args)
  expect_s3_class(per_profile, "markov_state_time")
  expect_equal(do.call(avg_state_time, c(args, list(by = FALSE))), per_profile)
  expect_equal(do.call(avg_state_time, c(args, list(by = NULL))), per_profile)
  expect_error(do.call(avg_state_time, c(args, list(by = 2))), "`by` must be")
})

test_that("bare variable names choose values like marginaleffects' predictions", {
  case <- make_comparisons_case()
  args <- list(case$model, times = 1:3, absorb = 3)
  expect_equal(
    do.call(avg_sops, c(args, list(variables = "tx"))),
    do.call(avg_sops, c(args, list(variables = list(tx = c(0, 1)))))
  )
  expect_equal(
    do.call(sops, c(args, list(variables = "tx"))),
    do.call(sops, c(args, list(variables = list(tx = c(0, 1)))))
  )
  expect_equal(
    do.call(avg_sops, c(args, list(variables = c("tx", "marker")))),
    do.call(
      avg_sops,
      c(args, list(variables = list(tx = c(0, 1), marker = c(0, 1))))
    )
  )
  expect_equal(
    do.call(avg_state_time, c(args, list(variables = "tx", state_sets = "1"))),
    do.call(
      avg_state_time,
      c(args, list(variables = list(tx = c(0, 1)), state_sets = "1"))
    )
  )
  expect_error(
    do.call(avg_sops, c(args, list(variables = "missing_column"))),
    "Variables not in data: missing_column"
  )

  # Continuous numeric variables use the five-number summary; one observed
  # value is enough for predictions.
  expect_identical(
    prediction_variable_values(c(10, 1, 4, 2, 3, NA), "age"),
    c(1, 2, 3, 4, 10)
  )
  expect_identical(prediction_variable_values(c(5, 5, 7), "dose"), c(5, 6, 7))
  expect_identical(prediction_variable_values(c(1, 0, 1), "tx"), c(0, 1))
  expect_identical(prediction_variable_values(c(1, 1), "tx"), 1)
  expect_identical(
    prediction_variable_values(factor(c("b", "a"), c("a", "b", "c")), "f"),
    factor(c("a", "b"), levels = c("a", "b", "c"))
  )
  expect_error(prediction_variable_values(c(NA, NA), "x"), "at least one")
})

test_that("draw weights let users average sops() draws themselves", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    absorb = 3
  )
  for (setting in list(
    list(method = "mvn", n_draws = 4L, seed = 501, weight = NULL),
    list(method = "score_bootstrap", n_draws = 4L, seed = 502, weight = "score_weight"),
    list(method = "fwb", n_draws = 2L, seed = 503, weight = "fwb_weight")
  )) {
    unit <- inferences(
      do.call(sops, args),
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    )
    average <- inferences(
      do.call(avg_sops, args),
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    )
    draws <- get_draws(unit)
    weights <- if (is.null(setting$weight)) {
      rep(1, nrow(draws))
    } else {
      expect_true(attr(unit, "draw_weights_attached"))
      # Both copies of a patient share that patient's weight in each draw.
      copies <- unique(draws[, c("draw_id", "rowidcf", setting$weight)])
      expect_identical(
        nrow(copies),
        setting$n_draws * length(unique(draws$rowidcf))
      )
      draws[[setting$weight]]
    }
    key <- paste(draws$draw_id, draws$time, draws$state, draws$tx)
    manual <- c(tapply(draws$draw * weights, key, sum)) /
      c(tapply(weights, key, sum))
    average_draws <- attr(average, "draws")
    average_key <- paste(
      average_draws$draw_id,
      average_draws$time,
      average_draws$state,
      average_draws$tx
    )
    expect_equal(
      unname(manual[average_key]),
      average_draws$estimate,
      tolerance = 1e-10
    )
  }
})
