test_that("row binding preserves column types, order and input data", {
  first <- data.frame(
    state = ordered(c("ill", "well"), levels = c("ill", "well")),
    day = as.Date(c("2026-01-01", "2026-01-02"))
  )
  second <- data.frame(
    value = 3,
    state = ordered("well", levels = c("ill", "well"))
  )
  inputs <- list(first, NULL, second)
  before <- serialize(inputs, NULL)
  out <- bind_rows_fill(inputs)
  expect_identical(class(out), "data.frame")
  expect_named(out, c("state", "day", "value"))
  expect_identical(
    out$state,
    ordered(c("ill", "well", "well"), levels = c("ill", "well"))
  )
  expect_identical(out$day, c(first$day, as.Date(NA)))
  expect_identical(out$value, c(NA_real_, NA_real_, 3))
  expect_identical(serialize(inputs, NULL), before)
  expect_identical(bind_rows_fill(list(NULL)), data.frame())
  expect_identical(bind_rows_fill(list(first[FALSE, ], first)), first)
})

test_that("table helpers preserve matrix-valued model predictors", {
  x <- data.frame(id = c(2L, 1L))
  x$basis <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(NULL, c("linear", "curved"))
  )
  y <- data.frame(id = c(1L, 2L, 2L), value = 1:3)
  bound <- bind_rows_fill(list(x, x))
  expect_identical(bound$basis, rbind(x$basis, x$basis))
  joined <- left_join_preserve_order(x, y, "id")
  expect_identical(joined$basis, x$basis[c(1, 1, 2), ])
  expect_identical(joined$value, c(2L, 3L, 1L))
  ids <- data.frame(original_id = c(2L, 2L), new_id = c("a", "b"), boot_id = 1L)
  sample <- materialize_bootstrap_sample_indexed(
    ids,
    x,
    "id",
    bootstrap_row_plan(x, "id")
  )
  expect_identical(sample$basis, x$basis[c(1, 1), ])
})

test_that("left joins keep the left row order, repeated and missing matches, typed keys and overlapping names", {
  ids <- data.frame(id = 1:2)
  keyed_x <- data.frame(key = c("b", "a", "b", NA, "NA", "absent"), left = 1:6)
  keyed_y <- data.frame(key = c("b", "b", "a", NA, "NA"), right = 11:15)
  control_x <- data.frame(a = c("a\rb", "a"), b = c("c", "b\rc"))
  control_y <- control_x
  control_y$value <- 1:2

  cases <- list(
    list(
      label = "no keys and no right rows",
      x = ids,
      y = data.frame(value = integer()),
      by = character(),
      check = function(out, label) {
        expect_identical(
          out,
          data.frame(id = 1:2, value = c(NA_integer_, NA_integer_)),
          info = label
        )
      }
    ),
    list(
      label = "no keys pairs every left row with every right row",
      x = ids,
      y = data.frame(value = 3:4),
      by = character(),
      check = function(out, label) {
        expect_identical(
          out,
          data.frame(id = c(1L, 1L, 2L, 2L), value = c(3L, 4L, 3L, 4L)),
          info = label
        )
      }
    ),
    list(
      label = "a non-key column on both sides",
      x = data.frame(id = 1:2, value = c(1, 2)),
      y = data.frame(id = 2:1, value = c(3, 4)),
      by = "id",
      check = function(out, label) {
        expect_named(out, c("id", "value", "value"))
        expect_identical(out[[3]], c(4, 3), info = label)
      }
    ),
    list(
      label = "repeated keys, NA and the string \"NA\"",
      x = keyed_x,
      y = keyed_y,
      by = "key",
      check = function(out, label) {
        expect_identical(class(out), "data.frame", info = label)
        expect_identical(
          out$left,
          c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 6L),
          info = label
        )
        expect_identical(
          out$right,
          c(11L, 12L, 13L, 11L, 12L, 14L, 15L, NA_integer_),
          info = label
        )
      }
    ),
    list(
      label = "no left rows",
      x = keyed_x[FALSE, ],
      y = keyed_y,
      by = "key",
      check = function(out, label) {
        expect_equal(nrow(out), 0L, info = label)
      }
    ),
    list(
      label = "no right rows",
      x = keyed_x,
      y = keyed_y[FALSE, ],
      by = "key",
      check = function(out, label) {
        expect_identical(out$right, rep(NA_integer_, nrow(keyed_x)), info = label)
      }
    ),
    list(
      label = "keys with control characters",
      x = control_x,
      y = control_y,
      by = c("a", "b"),
      check = function(out, label) {
        expect_identical(out$value, 1:2, info = label)
      }
    )
  )

  for (case in cases) {
    before <- serialize(list(case$x, case$y), NULL)
    out <- left_join_preserve_order(case$x, case$y, case$by)
    case$check(out, case$label)
    expect_identical(serialize(list(case$x, case$y), NULL), before, info = case$label)
  }
})

test_that("table summaries accept copied data.table inputs", {
  x <- data.table::data.table(group = c("a", "a"), estimate = c(1, 3))
  before <- data.table::copy(x)
  expect_equal(
    aggregate_value(x, "estimate", "group", mean),
    data.frame(group = "a", estimate = 2)
  )
  expect_equal(compute_ci_from_draws(x, "group")$std.error, sqrt(2))
  expect_identical(x, before)
})

test_that("grouped reductions match aggregate ordering and missing-value rules", {
  data <- data.frame(
    state = factor(
      c("z", "a", "z", "a", "z", "missing", NA),
      levels = c("z", "a", "missing", "unused")
    ),
    time = c(10, 2, 2, 10, 2, 2, 2),
    estimate = c(1, 2, 3, 4, NA, NA, 7)
  )
  before <- serialize(data, NULL)
  for (fun in list(mean, sum)) {
    expected <- stats::aggregate(
      estimate ~ state + time,
      data,
      fun,
      na.rm = TRUE
    )
    expect_equal(
      aggregate_value(data, "estimate", c("state", "time"), fun),
      expected
    )
  }
  expect_identical(serialize(data, NULL), before)
  expect_equal(
    aggregate_value(data, "estimate", character(), sum),
    data.frame(estimate = 17)
  )
})

test_that("weighted reductions retain missing and zero-weight groups", {
  data <- data.frame(
    group = factor(
      c("z", "z", "a", "a", "b", "b", NA),
      levels = c("z", "a", "b", "unused")
    ),
    estimate = c(1, 3, NA, NA, 4, 5, 6),
    weight = c(1, 3, 1, 1, 0, 0, 1)
  )
  out <- aggregate_sops_estimates(data, "group", "weight")
  expect_identical(out$group, data$group[c(1, 3, 5)])
  expect_identical(out$estimate, c(2.5, NA_real_, NA_real_))
})

test_that("draw summaries preserve missing-value rules and group order", {
  draws <- data.frame(
    state = factor(
      c("z", "z", "a", "a", "a", "empty", NA),
      levels = c("z", "a", "empty", "unused")
    ),
    estimate = c(1, 3, 2, 6, NA, NA, 9)
  )
  before <- serialize(draws, NULL)
  out <- compute_ci_from_draws(draws, "state")
  expect_identical(out$state, draws$state[c(1, 3)])
  expect_equal(out$conf.low, c(1.05, 2.1))
  expect_equal(out$conf.high, c(2.95, 5.9))
  expect_equal(out$std.error, c(sqrt(2), sqrt(8)))
  centers <- data.frame(
    state = factor(c("a", "z"), levels = levels(draws$state)),
    estimate = c(10, 20)
  )
  wald <- compute_ci_from_draws(
    draws,
    "state",
    conf_type = "wald",
    point_estimates = centers
  )
  expected_center <- centers$estimate[match(wald$state, centers$state)]
  expect_equal((wald$conf.low + wald$conf.high) / 2, expected_center)
  expect_identical(serialize(draws, NULL), before)
})
