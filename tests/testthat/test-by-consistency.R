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

test_that("per-profile functions and their avg_*() partners hand `by` to each other", {
  skip_on_cran()
  case <- make_comparisons_case()
  model <- case$model
  tx <- list(tx = c(0, 1))

  # `unit` gives one row per profile and `average` averages over profiles.
  # `bad_by` lists the calls that must refuse an invalid `by`.
  cases <- list(
    list(
      label = "sops()",
      unit = sops,
      average = avg_sops,
      unit_class = "markov_sops",
      average_class = "markov_avg_sops",
      args = list(model, times = 1:3, absorb = 3),
      bad_by = list(list(sops, 1), list(avg_sops, NA))
    ),
    list(
      label = "sops() with variables",
      unit = sops,
      average = avg_sops,
      unit_class = "markov_sops",
      average_class = "markov_avg_sops",
      args = list(model, variables = tx, times = 1:3, absorb = 3)
    ),
    list(
      label = "state_time()",
      unit = state_time,
      average = avg_state_time,
      unit_class = "markov_state_time",
      average_class = "markov_avg_state_time",
      args = list(
        model,
        variables = tx,
        state_sets = "1",
        times = 1:3,
        absorb = 3
      ),
      bad_by = list(list(avg_state_time, 2))
    ),
    list(
      label = "comparisons()",
      unit = comparisons,
      average = avg_comparisons,
      unit_class = "markov_comparisons",
      average_class = "markov_avg_comparisons",
      args = list(
        model,
        variables = tx,
        estimand = "time_in_state",
        comparison = "ratio",
        times = 1:4,
        absorb = 3
      ),
      average_rows = 3L,
      bad_by = list(
        list(comparisons, 1),
        list(comparisons, NA_character_),
        list(avg_comparisons, 1)
      )
    )
  )

  for (row in cases) {
    call <- function(fun, ...) do.call(fun, c(row$args, list(...)))
    averaged <- call(row$average)
    expect_true(inherits(averaged, row$average_class), info = row$label)
    if (!is.null(row$average_rows)) {
      expect_identical(nrow(averaged), row$average_rows, info = row$label)
    }
    expect_equal(call(row$average, by = TRUE), averaged, info = row$label)
    expect_equal(call(row$unit, by = TRUE), averaged, info = row$label)
    expect_equal(
      call(row$unit, by = "marker"),
      call(row$average, by = "marker"),
      info = row$label
    )

    per_profile <- call(row$unit)
    expect_true(inherits(per_profile, row$unit_class), info = row$label)
    expect_equal(call(row$average, by = FALSE), per_profile, info = row$label)
    expect_equal(call(row$average, by = NULL), per_profile, info = row$label)
    expect_equal(call(row$unit, by = NULL), per_profile, info = row$label)

    for (bad in row$bad_by) {
      expect_error(call(bad[[1]], by = bad[[2]]), "`by` must be", info = row$label)
    }
  }

  # Per-patient results have no average over sampled patients.
  expect_error(
    inferences(
      avg_comparisons(
        model,
        variables = tx,
        estimand = "time_in_state",
        comparison = "ratio",
        times = 1:4,
        absorb = 3,
        by = FALSE
      ),
      method = "delta",
      vcov = "unconditional"
    ),
    "avg_comparisons(by = FALSE)",
    fixed = TRUE
  )
})

test_that("bare variable names choose values like marginaleffects", {
  skip_on_cran()
  case <- make_comparisons_case()
  sop_args <- list(case$model, times = 1:3, absorb = 3)
  comparison_args <- list(case$model, state_sets = "1", times = 1:4, absorb = 3)
  starts <- factor(c("1", "2"), levels = levels(case$data$yprev))

  # Each bare name gives the same result as the values written out. Factors
  # and binary variables use their observed values; comparisons compare each
  # later level with the first.
  calls <- list(
    list(
      label = "avg_sops()",
      fun = avg_sops,
      args = sop_args,
      bare = "tx",
      explicit = list(tx = c(0, 1))
    ),
    list(
      label = "sops()",
      fun = sops,
      args = sop_args,
      bare = "tx",
      explicit = list(tx = c(0, 1))
    ),
    list(
      label = "avg_sops(), two variables",
      fun = avg_sops,
      args = sop_args,
      bare = c("tx", "marker"),
      explicit = list(tx = c(0, 1), marker = c(0, 1))
    ),
    list(
      label = "avg_state_time()",
      fun = avg_state_time,
      args = c(sop_args, list(state_sets = "1")),
      bare = "tx",
      explicit = list(tx = c(0, 1))
    ),
    list(
      label = "comparisons()",
      fun = comparisons,
      args = comparison_args,
      bare = "tx",
      explicit = list(tx = c(0, 1))
    ),
    list(
      label = "avg_comparisons()",
      fun = avg_comparisons,
      args = comparison_args,
      bare = "tx",
      explicit = list(tx = c(0, 1))
    ),
    list(
      label = "avg_comparisons(), starting state",
      fun = avg_comparisons,
      args = comparison_args,
      bare = "yprev",
      explicit = list(yprev = starts),
      contrast = "2 - 1"
    )
  )
  for (call in calls) {
    bare <- do.call(call$fun, c(call$args, list(variables = call$bare)))
    explicit <- do.call(call$fun, c(call$args, list(variables = call$explicit)))
    expect_equal(bare, explicit, info = call$label)
    if (!is.null(call$contrast)) {
      expect_identical(unique(bare$contrast), call$contrast, info = call$label)
    }
  }

  refusals <- list(
    list(
      label = "avg_sops(), missing column",
      fun = avg_sops,
      args = sop_args,
      variables = "missing_column",
      phrase = "Variables not in data: missing_column"
    ),
    list(
      label = "comparisons(), two variables",
      fun = comparisons,
      args = comparison_args,
      variables = c("tx", "marker"),
      phrase = "exactly one named variable"
    ),
    list(
      label = "comparisons(), missing column",
      fun = comparisons,
      args = comparison_args,
      variables = "missing_column",
      phrase = "Variables not in data: missing_column"
    )
  )
  for (call in refusals) {
    expect_error(
      do.call(call$fun, c(call$args, list(variables = call$variables))),
      call$phrase,
      fixed = TRUE,
      info = call$label
    )
  }

  # The value rules. Predictions summarize a continuous numeric variable by
  # its five-number summary and need one observed value; comparisons need two
  # observed values and refuse to choose values of a continuous variable.
  # `values` is the expected result, `error` the expected refusal.
  value_rule <- function(kind, x, name, values = NULL, error = NULL) {
    list(kind = kind, x = x, name = name, values = values, error = error)
  }
  rules <- list(
    value_rule("prediction", c(10, 1, 4, 2, 3, NA), "age", values = c(1, 2, 3, 4, 10)),
    value_rule("prediction", c(5, 5, 7), "dose", values = c(5, 6, 7)),
    value_rule("prediction", c(1, 0, 1), "tx", values = c(0, 1)),
    value_rule("prediction", c(1, 1), "tx", values = 1),
    value_rule(
      "prediction",
      factor(c("b", "a"), c("a", "b", "c")),
      "f",
      values = factor(c("a", "b"), levels = c("a", "b", "c"))
    ),
    value_rule("prediction", c(NA, NA), "x", error = "at least one"),
    value_rule(
      "comparison",
      factor(c("b", "a", "b"), c("b", "a", "c")),
      "x",
      values = factor(c("b", "a"), levels = c("b", "a", "c"))
    ),
    value_rule("comparison", c(TRUE, FALSE, NA), "x", values = c(FALSE, TRUE)),
    value_rule("comparison", c("x", "b", "x"), "x", values = c("b", "x")),
    value_rule("comparison", c(1, 0, 1), "x", values = c(0, 1)),
    value_rule("comparison", c(1L, 0L), "x", values = c(0L, 1L)),
    value_rule("comparison", c(40, 52.5, 61), "age", error = "list(age = "),
    value_rule("comparison", c(2, 3), "dose", error = "list(dose = "),
    value_rule("comparison", c(1, 1, NA), "x", error = "at least two"),
    value_rule("comparison", as.Date("2026-01-01"), "x", error = "values")
  )
  for (rule in rules) {
    fun <- switch(
      rule$kind,
      prediction = prediction_variable_values,
      comparison = comparison_variable_values
    )
    label <- paste0(rule$kind, " values of ", paste(format(rule$x), collapse = ", "))
    if (!is.null(rule$error)) {
      expect_error(fun(rule$x, rule$name), rule$error, fixed = TRUE, info = label)
    } else {
      expect_identical(fun(rule$x, rule$name), rule$values, info = label)
    }
  }
})

test_that("averaging per-profile results reproduces the avg_*() result", {
  skip_on_cran()
  case <- make_comparisons_case()
  starts <- factor(c("1", "2"), levels = levels(case$profiles$yprev))
  three <- list(tx = c(0, 1, 2))
  tx <- list(tx = c(0, 1))
  comparison_keys <- c("time", "state_set", "contrast")
  state_time_keys <- function(variables) c("state_set", names(variables))

  # `keys` are the columns that identify an averaged row.
  cases <- list(
    list(
      label = "comparisons(), SOPs of lumped states",
      unit = comparisons,
      average = avg_comparisons,
      args = list(
        variables = three,
        estimand = "sop",
        state_sets = list(alive = c("1", "2"), dead = "3")
      ),
      keys = comparison_keys
    ),
    list(
      label = "comparisons(), time in state",
      unit = comparisons,
      average = avg_comparisons,
      args = list(variables = three, estimand = "time_in_state", state_sets = NULL),
      keys = comparison_keys
    ),
    list(
      label = "comparisons(), real-time time in state",
      unit = comparisons,
      average = avg_comparisons,
      args = list(
        variables = three,
        estimand = "time_in_state",
        state_sets = "1",
        time_map = case$time_map
      ),
      keys = comparison_keys
    ),
    list(
      label = "comparisons(), time benefit",
      unit = comparisons,
      average = avg_comparisons,
      args = list(variables = three, estimand = "time_benefit", state_sets = NULL),
      keys = comparison_keys
    ),
    list(
      label = "comparisons(), real-time time benefit",
      unit = comparisons,
      average = avg_comparisons,
      args = list(
        variables = three,
        estimand = "time_benefit",
        state_sets = NULL,
        time_map = case$time_map
      ),
      keys = comparison_keys
    ),
    list(
      label = "state_time(), observed profiles",
      unit = state_time,
      average = avg_state_time,
      args = list(variables = NULL, time_map = NULL, state_sets = NULL),
      keys = state_time_keys(NULL)
    ),
    list(
      label = "state_time(), treatment",
      unit = state_time,
      average = avg_state_time,
      args = list(variables = tx, time_map = NULL, state_sets = "1"),
      keys = state_time_keys(tx)
    ),
    list(
      label = "state_time(), real-time lumped states",
      unit = state_time,
      average = avg_state_time,
      args = list(
        variables = tx,
        time_map = case$time_map,
        state_sets = list(alive = c("1", "2"), dead = "3")
      ),
      keys = state_time_keys(tx)
    ),
    list(
      label = "state_time(), starting state",
      unit = state_time,
      average = avg_state_time,
      args = list(
        variables = list(yprev = starts),
        time_map = case$time_map,
        state_sets = "1"
      ),
      keys = state_time_keys(list(yprev = starts))
    )
  )

  for (row in cases) {
    args <- c(list(case$model, times = 1:4, absorb = 3), row$args)
    unit <- do.call(row$unit, args)
    average <- do.call(row$average, args)
    keys <- intersect(row$keys, names(average))
    means <- stats::aggregate(unit$estimate, unit[, keys, drop = FALSE], mean)
    index <- match(
      do.call(paste, average[, keys, drop = FALSE]),
      do.call(paste, means[, keys, drop = FALSE])
    )
    expect_false(anyNA(index), info = row$label)
    expect_equal(means$x[index], average$estimate, tolerance = 1e-10, info = row$label)
    # One row per profile for each averaged row.
    profile <- if ("rowidcf" %in% names(unit)) unit$rowidcf else unit$rowid
    expect_identical(
      nrow(unit),
      nrow(average) * length(unique(profile)),
      info = row$label
    )
  }
})

test_that("draw weights let users average per-profile draws themselves", {
  skip_on_cran()
  comparisons_case <- make_comparisons_case()
  observed_case <- make_observed_average_case()
  model <- comparisons_case$model
  tx <- list(tx = c(0, 1))

  sop_args <- list(model, variables = tx, times = 1:2, absorb = 3)
  state_time_args <- list(
    model,
    variables = tx,
    state_sets = "1",
    times = 1:3,
    absorb = 3,
    time_map = c("1" = 1, "2" = 2, "3" = 4)
  )
  observed_args <- list(
    observed_case$model,
    times = observed_case$times,
    absorb = observed_case$absorb
  )
  real_time_args <- function(estimand) {
    list(
      model,
      variables = tx,
      estimand = estimand,
      state_sets = if (estimand == "time_in_state") "1",
      times = 1:4,
      absorb = 3,
      time_map = comparisons_case$time_map
    )
  }
  comparison_args <- list(
    model,
    variables = tx,
    state_sets = "1",
    times = 1:4,
    absorb = 3
  )

  # A setting draws `n_draws` with `method` and `seed` for the per-profile
  # result of `unit` and for the average of `average`. Without `weight`, each
  # profile counts once in each draw; otherwise its weight in each draw is in
  # column `weight`. `keys` identify a row of averaged draws, and `profile`
  # identifies a patient among the per-profile draws.
  setting <- function(label, unit, average, args, method, n_draws, seed,
                      weight = NULL, keys, profile) {
    list(
      label = label,
      unit = unit,
      average = average,
      args = args,
      method = method,
      n_draws = n_draws,
      seed = seed,
      weight = weight,
      keys = keys,
      profile = profile
    )
  }
  sop_keys <- c("draw_id", "time", "state", "tx")
  cases <- list(
    setting("sops(), mvn", sops, avg_sops, sop_args, "mvn", 4L, 501,
      keys = sop_keys, profile = "rowidcf"),
    setting("sops(), score bootstrap", sops, avg_sops, sop_args, "score_bootstrap",
      4L, 502, "score_weight", sop_keys, "rowidcf"),
    setting("sops(), fwb", sops, avg_sops, sop_args, "fwb", 2L, 503, "fwb_weight",
      sop_keys, "rowidcf"),
    setting("sops() on stored profiles, mvn", sops, avg_sops, observed_args, "mvn",
      4L, 4211, keys = c("draw_id", "time", "state"), profile = "rowid"),
    setting("sops() on stored profiles, score bootstrap", sops, avg_sops,
      observed_args, "score_bootstrap", 4L, 4211, "score_weight",
      c("draw_id", "time", "state"), "rowid"),
    setting("state_time(), mvn", state_time, avg_state_time, state_time_args,
      "mvn", 4L, 601, keys = c("draw_id", "tx"), profile = "rowidcf"),
    setting("state_time(), score bootstrap", state_time, avg_state_time,
      state_time_args, "score_bootstrap", 4L, 602, "score_weight",
      c("draw_id", "tx"), "rowidcf"),
    setting("state_time(), fwb", state_time, avg_state_time, state_time_args,
      "fwb", 2L, 603, "fwb_weight", c("draw_id", "tx"), "rowidcf"),
    setting("comparisons(), real-time time in state, mvn", comparisons,
      avg_comparisons, real_time_args("time_in_state"), "mvn", 6L, 202,
      keys = "draw_id", profile = "rowid"),
    setting("comparisons(), real-time time benefit, mvn", comparisons,
      avg_comparisons, real_time_args("time_benefit"), "mvn", 6L, 202,
      keys = "draw_id", profile = "rowid"),
    setting("comparisons(), score bootstrap", comparisons, avg_comparisons,
      comparison_args, "score_bootstrap", 4L, 301, "score_weight",
      c("draw_id", "time"), "rowid"),
    setting("comparisons(), fwb", comparisons, avg_comparisons, comparison_args,
      "fwb", 3L, 302, "fwb_weight", c("draw_id", "time"), "rowid")
  )

  for (row in cases) {
    infer <- function(fun) {
      run_inferences_method(
        do.call(fun, row$args),
        method = row$method,
        n_draws = row$n_draws,
        seed = row$seed
      )
    }
    unit <- infer(row$unit)
    average <- infer(row$average)
    expect_identical(attr(unit, "method"), row$method, info = row$label)
    if (identical(row$unit, state_time)) {
      expect_true(inherits(unit, "markov_state_time"), info = row$label)
      expect_null(attr(unit, "null"), info = row$label)
    }
    draws <- get_draws(unit)
    weights <- if (is.null(row$weight)) {
      rep(1, nrow(draws))
    } else {
      expect_true(attr(unit, "draw_weights_attached"), info = row$label)
      expect_identical(attr(unit, "draw_weight_col"), row$weight, info = row$label)
      expect_true(row$weight %in% names(attr(unit, "draws")), info = row$label)
      # All rows of a patient share that patient's weight in each draw.
      copies <- unique(draws[, c("draw_id", row$profile, row$weight)])
      expect_identical(
        nrow(copies),
        row$n_draws * length(unique(draws[[row$profile]])),
        info = row$label
      )
      draws[[row$weight]]
    }
    key <- do.call(paste, unname(as.list(draws[row$keys])))
    manual <- c(tapply(draws$draw * weights, key, sum)) /
      c(tapply(weights, key, sum))
    average_draws <- attr(average, "draws")
    average_key <- do.call(paste, unname(as.list(average_draws[row$keys])))
    expect_false(anyNA(manual[average_key]), info = row$label)
    expect_equal(
      unname(manual[average_key]),
      average_draws$estimate,
      tolerance = 1e-12,
      info = row$label
    )
    expect_identical(
      nrow(get_draws(average)),
      nrow(average) * row$n_draws,
      info = row$label
    )
  }
})
