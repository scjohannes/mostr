test_that("comparisons() returns paired SOP differences for each profile", {
  skip_on_cran()
  case <- make_comparisons_case()
  out <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    times = 1:4,
    absorb = 3
  )
  lo <- sops_with_value(case, "tx", 0)
  hi <- sops_with_value(case, "tx", 1)
  n <- nrow(case$profiles)

  expect_s3_class(out, "markov_comparisons")
  expect_false(inherits(out, "markov_sops"))
  expect_false(inherits(out, "markov_avg_comparisons"))
  expect_identical(nrow(out), n * 4L * 3L)
  for (state in as.character(1:3)) {
    rows <- out[out$state_set == state, , drop = FALSE]
    keep <- as.character(hi$state) == state
    expect_equal(rows$rowid, hi$rowid[keep])
    expect_equal(rows$time, hi$time[keep])
    expect_equal(rows$estimate, hi$estimate[keep] - lo$estimate[keep])
  }
  expect_true(all(out$estimand == "sop"))
  expect_true(all(out$term == "tx"))
  expect_true(all(out$contrast == "1 - 0"))
  expect_true(all(out$comparison == "difference"))
  # Profile columns hold each patient's observed values.
  expect_equal(out$tx, case$profiles$tx[out$rowid])
  expect_equal(out$id, case$profiles$id[out$rowid])
  expect_equal(out$marker, case$profiles$marker[out$rowid])
  expect_identical(names(out)[ncol(out)], "estimate")
})

test_that("comparisons() lumps state sets and forms ratios for every level", {
  case <- make_comparisons_case()
  out <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1, 2)),
    state_sets = list(alive = c("1", "2")),
    comparison = "ratio",
    times = 1:4,
    absorb = 3
  )
  alive <- function(x) {
    keep <- as.character(x$state) %in% c("1", "2")
    stats::aggregate(estimate ~ rowid + time, x[keep, ], sum)
  }
  lo <- alive(sops_with_value(case, "tx", 0))
  expect_identical(unique(out$state_set), "alive")
  expect_identical(unique(out$contrast), c("1 / 0", "2 / 0"))
  for (level in c(1, 2)) {
    hi <- alive(sops_with_value(case, "tx", level))
    rows <- out[out$comparison_level == level, , drop = FALSE]
    expected <- merge(hi, lo, by = c("rowid", "time"))
    expected <- expected[match(
      paste(rows$rowid, rows$time),
      paste(expected$rowid, expected$time)
    ), ]
    expect_equal(rows$estimate, expected$estimate.x / expected$estimate.y)
  }
})

test_that("comparisons() starts real-time totals from a set starting state", {
  case <- make_comparisons_case()
  levels <- levels(case$profiles$yprev)
  starts <- factor(c("1", "2"), levels = levels)
  args <- list(
    case$model,
    variables = list(yprev = starts),
    estimand = "time_in_state",
    state_sets = "1",
    times = 1:4,
    absorb = 3,
    time_map = case$time_map,
    target_times = seq(0, 7, by = 0.5)
  )
  unit <- do.call(comparisons, args)
  average <- do.call(avg_comparisons, args)
  expect_equal(mean(unit$estimate), average$estimate, tolerance = 1e-10)
  stored <- attr(unit, "comparison_baseline_data")
  expect_equal(unit$yprev, stored$yprev[unit$rowid])
  # Stored starting profiles carry internal bookkeeping columns that are not
  # patient variables.
  expect_true(".markov_source_row" %in% names(stored))
  expect_false(any(startsWith(names(unit), ".markov_")))
})

test_that("per-profile and average summaries validate arguments before predicting", {
  case <- make_comparisons_case()
  model <- case$model
  local_mocked_bindings(
    sops = function(...) {
      stop("Predictions should not run before argument validation.")
    },
    avg_comparison_replay_avg_sops = function(...) {
      stop("Predictions should not run before argument validation.")
    }
  )
  comparison_clash <- case$profiles
  comparison_clash$contrast <- "a"
  comparison_clash$estimate <- 1
  state_set_clash <- case$profiles
  state_set_clash$state_set <- "x"
  tx <- list(tx = c(0, 1))

  cases <- list(
    list(
      label = "comparisons(), profile columns named like output columns",
      fun = comparisons,
      args = list(newdata = comparison_clash, variables = tx, times = 1:4),
      phrase = "`contrast`, `estimate`"
    ),
    list(
      label = "comparisons(), no times",
      fun = comparisons,
      args = list(variables = tx),
      phrase = "`times` must be supplied to `comparisons()`"
    ),
    list(
      label = "comparisons(), one value",
      fun = comparisons,
      args = list(variables = list(tx = 1), times = 1:4),
      phrase = "at least two values"
    ),
    list(
      label = "comparisons(), two variables",
      fun = comparisons,
      args = list(variables = list(tx = c(0, 1), marker = c(0, 1)), times = 1:4),
      phrase = "exactly one named variable"
    ),
    list(
      label = "comparisons(), duplicate times",
      fun = comparisons,
      args = list(variables = tx, times = c(1, 1, 2)),
      phrase = "duplicate visits"
    ),
    list(
      label = "comparisons(), time-benefit ratio",
      fun = comparisons,
      args = list(
        variables = tx,
        estimand = "time_benefit",
        comparison = "ratio",
        times = 1:4
      ),
      phrase = "only supports `comparison = \"difference\"`"
    ),
    list(
      label = "state_time(), profile column named like an output column",
      fun = state_time,
      args = list(newdata = state_set_clash, times = 1:4),
      phrase = "`state_set`"
    ),
    list(
      label = "state_time(), no times",
      fun = state_time,
      args = list(),
      phrase = "`times` must be supplied"
    ),
    list(
      label = "state_time(), baseline time without time map",
      fun = state_time,
      args = list(times = 1:4, baseline_time = 1),
      phrase = "`baseline_time` requires `time_map`"
    ),
    list(
      label = "avg_state_time(), duplicate times",
      fun = avg_state_time,
      args = list(times = c(1, 1, 2), absorb = 3),
      phrase = "`times` must not contain duplicate"
    ),
    list(
      label = "avg_state_time(), empty state sets",
      fun = avg_state_time,
      args = list(times = 1:2, absorb = 3, state_sets = list()),
      phrase = "`state_sets`"
    ),
    list(
      label = "avg_state_time(), duplicate state-set names",
      fun = avg_state_time,
      args = list(times = 1:2, absorb = 3, state_sets = list(a = "1", a = "2")),
      phrase = "`state_sets`"
    ),
    list(
      label = "avg_state_time(), by a reserved output name",
      fun = avg_state_time,
      args = list(times = 1:2, absorb = 3, by = "state_set"),
      phrase = "reserved output names"
    ),
    list(
      label = "avg_state_time(), baseline time without time map",
      fun = avg_state_time,
      args = list(times = 1:2, absorb = 3, baseline_time = 1),
      phrase = "`baseline_time` requires `time_map`"
    )
  )

  for (row in cases) {
    expect_error(
      do.call(row$fun, c(list(model), row$args)),
      row$phrase,
      fixed = TRUE,
      info = row$label
    )
  }
})

test_that("posterior comparisons summarize paired draws for each profile", {
  skip_on_cran()
  coefficients <- rbind(
    c(0.5, -0.5, 0.4, 1.2),
    c(2, 1, -0.8, -1.1),
    c(-0.4, -2, 1.5, 0.2)
  )
  colnames(coefficients) <- c("y>=2", "y>=3", "tx", "time")
  model <- fake_blrm(coefficients)
  profiles <- data.frame(
    id = c("a", "b"),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1
  )
  local_mocked_bindings(blrm_design_matrix = fake_blrm_design(c("tx", "time")))
  posterior_sops <- function(value) {
    newdata <- profiles
    newdata$tx <- value
    attr(
      sops(
        model,
        newdata = newdata,
        times = 1:3,
        n_draws = 3,
        seed = 9325,
        return_draws = TRUE
      ),
      "draws"
    )
  }
  out <- comparisons(
    model,
    newdata = profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "median",
    return_draws = TRUE
  )
  lo <- posterior_sops(0)
  hi <- posterior_sops(1)
  keep <- as.character(lo$state) == "1"
  paired <- data.frame(
    rowid = lo$rowid[keep],
    time = lo$time[keep],
    draw_id = lo$draw_id[keep],
    estimate = hi$estimate[keep] - lo$estimate[keep]
  )
  for (i in seq_len(nrow(out))) {
    values <- paired$estimate[
      paired$rowid == out$rowid[i] & paired$time == out$time[i]
    ]
    expect_equal(out$estimate[i], stats::median(values))
    expect_equal(out$std.error[i], stats::sd(values))
    expect_equal(out$conf.low[i], unname(stats::quantile(values, 0.025)))
    expect_equal(out$conf.high[i], unname(stats::quantile(values, 0.975)))
  }
  expect_equal(out$tx, profiles$tx[out$rowid])
  expect_identical(attr(out, "method"), "posterior")
  draws <- get_draws(out)
  expect_identical(nrow(draws), 3L * nrow(out))
  expect_true(all(c("draw", "draw_id", "contrast", "tx") %in% names(draws)))
  expect_identical(inferences(out, method = "mvn", n_draws = 2), out)

  compact <- comparisons(
    model,
    newdata = profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:3,
    n_draws = 3,
    seed = 9325
  )
  expect_null(attr(compact, "draws"))
  means <- stats::aggregate(estimate ~ rowid + time, paired, mean)
  index <- match(
    paste(compact$rowid, compact$time),
    paste(means$rowid, means$time)
  )
  expect_equal(compact$estimate, means$estimate[index])
})

# Draws of paired SOP differences in state `state` for each profile, from one
# `sops()` call on the profiles stacked under `tx = 0` and `tx = 1`, so both
# scenarios share each coefficient draw or refit. `...` goes to `sops()`.
paired_sop_draws <- function(case, method, n_draws, seed, state = "1", ...) {
  lo <- case$profiles
  lo$tx <- 0
  hi <- case$profiles
  hi$tx <- 1
  # User-supplied profiles get no draw weights; that warning is expected.
  inferred <- suppressWarnings(run_inferences_method(
    sops(case$model, newdata = rbind(lo, hi), times = 1:4, absorb = 3, ...),
    method = method,
    n_draws = n_draws,
    seed = seed
  ))
  draws <- attr(inferred, "draws")
  draws <- draws[
    as.character(draws$state) == state,
    c("rowid", "time", "draw_id", "estimate")
  ]
  n <- nrow(case$profiles)
  lo_draws <- draws[draws$rowid <= n, , drop = FALSE]
  hi_draws <- draws[draws$rowid > n, , drop = FALSE]
  hi_draws$rowid <- hi_draws$rowid - n
  paired <- merge(hi_draws, lo_draws, by = c("rowid", "time", "draw_id"))
  data.frame(
    rowid = paired$rowid,
    time = paired$time,
    draw_id = paired$draw_id,
    estimate = paired$estimate.x - paired$estimate.y
  )
}

test_that("draw-based intervals reduce each draw to paired profile comparisons", {
  skip_on_cran()
  case <- make_comparisons_case()
  cmp <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:4,
    absorb = 3
  )
  for (setting in list(
    list(method = "mvn", n_draws = 8L, seed = 101),
    list(method = "score_bootstrap", n_draws = 6L, seed = 102),
    list(method = "fwb", n_draws = 3L, seed = 103)
  )) {
    run <- function() {
      run_inferences_method(
        cmp,
        method = setting$method,
        n_draws = setting$n_draws,
        seed = setting$seed
      )
    }
    # Supplied profiles cannot be matched to the bootstrap patients, so the
    # patient-weighted methods warn that no draw weights are attached.
    if (identical(setting$method, "mvn")) {
      inferred <- run()
    } else {
      expect_warning(inferred <- run(), "draw weights were not attached")
      expect_false(attr(inferred, "draw_weights_attached"))
      expect_identical(
        attr(inferred, "draw_weight_omission_reason"),
        "user_supplied_newdata"
      )
    }
    expect_s3_class(inferred, "markov_comparisons")
    expect_identical(attr(inferred, "method"), setting$method)
    expect_equal(inferred$estimate, cmp$estimate)
    draws <- attr(inferred, "draws")
    expect_identical(nrow(draws), setting$n_draws * nrow(cmp))
    expect_false(any(c("score_weight", "fwb_weight") %in% names(draws)))
    expected <- paired_sop_draws(
      case,
      setting$method,
      setting$n_draws,
      setting$seed
    )
    index <- match(
      paste(draws$rowid, draws$time, draws$draw_id),
      paste(expected$rowid, expected$time, expected$draw_id)
    )
    expect_false(anyNA(index))
    expect_equal(draws$estimate, expected$estimate[index], tolerance = 1e-10)

    for (i in c(1L, nrow(inferred))) {
      values <- draws$estimate[
        draws$rowid == inferred$rowid[i] & draws$time == inferred$time[i]
      ]
      expect_equal(inferred$std.error[i], stats::sd(values))
      expect_equal(
        inferred$conf.low[i],
        unname(stats::quantile(values, 0.025))
      )
      expect_equal(
        inferred$conf.high[i],
        unname(stats::quantile(values, 0.975))
      )
    }
  }
})

test_that("comparison intervals add null tests and replace earlier inference", {
  skip_on_cran()
  case <- make_comparisons_case()
  cmp <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    times = 1:4,
    absorb = 3
  )
  first <- inferences(cmp, n_draws = 5, seed = 11)
  expect_identical(attr(first, "null"), 0)
  expect_equal(first$statistic, first$estimate / first$std.error)
  expect_true(all(c("p.value", "s.value") %in% names(first)))
  draws <- get_draws(first)
  expect_true(all(c("draw", "draw_id", "tx", "id") %in% names(draws)))

  again <- inferences(first, n_draws = 5, seed = 11)
  expect_identical(names(again), names(first))
  expect_equal(again$std.error, first$std.error)

  compact <- inferences(cmp, n_draws = 5, seed = 11, return_draws = FALSE)
  expect_null(attr(compact, "draws"))
  expect_equal(compact$std.error, first$std.error)

  ratio <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    comparison = "ratio",
    times = 1:4,
    absorb = 3
  )
  expect_identical(attr(inferences(ratio, n_draws = 5, seed = 11), "null"), 1)
})

# Checks one model setup end to end: point estimates against paired `sops()`
# differences, every draw method against paired `sops()` draws, and
# real-time totals against `avg_comparisons()`. `...` goes to every call, and
# `info` names the setup in failures.
expect_comparisons_match_paired_sops <- function(case, ..., info = NULL) {
  cmp <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:4,
    absorb = 3,
    ...
  )
  lo <- sops_with_value(case, "tx", 0, ...)
  hi <- sops_with_value(case, "tx", 1, ...)
  keep <- as.character(hi$state) == "1"
  expect_identical(as.character(cmp$time), as.character(hi$time[keep]), info = info)
  expect_equal(cmp$estimate, hi$estimate[keep] - lo$estimate[keep], info = info)

  for (setting in list(
    list(method = "mvn", n_draws = 5L, seed = 401),
    list(method = "score_bootstrap", n_draws = 4L, seed = 402),
    list(method = "fwb", n_draws = 2L, seed = 403)
  )) {
    inferred <- suppressWarnings(run_inferences_method(
      cmp,
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    ))
    draws <- attr(inferred, "draws")
    expected <- paired_sop_draws(
      case,
      setting$method,
      setting$n_draws,
      setting$seed,
      ...
    )
    index <- match(
      paste(draws$rowid, draws$time, draws$draw_id),
      paste(expected$rowid, expected$time, expected$draw_id)
    )
    expect_false(anyNA(index), info = info)
    expect_equal(
      draws$estimate,
      expected$estimate[index],
      tolerance = 1e-10,
      info = info
    )
  }

  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    state_sets = "1",
    times = 1:4,
    absorb = 3,
    time_map = case$time_map,
    ...
  )
  unit <- inferences(
    do.call(comparisons, args),
    method = "mvn",
    n_draws = 4,
    seed = 404
  )
  average <- inferences(
    do.call(avg_comparisons, args),
    method = "mvn",
    n_draws = 4,
    seed = 404
  )
  expect_equal(
    mean(unit$estimate),
    average$estimate,
    tolerance = 1e-10,
    info = info
  )
  unit_means <- stats::aggregate(estimate ~ draw_id, attr(unit, "draws"), mean)
  average_draws <- attr(average, "draws")
  expect_equal(
    unit_means$estimate,
    average_draws$estimate[match(unit_means$draw_id, average_draws$draw_id)],
    tolerance = 1e-10,
    info = info
  )
  invisible(cmp)
}

test_that("comparisons() matches paired sops() with factor visit times and second-order models", {
  skip_on_cran()
  cases <- list(
    list(label = "factor visit times", case = make_factor_comparisons_case()),
    list(
      label = "second-order model",
      case = make_second_order_comparisons_case(),
      p2_var = "ypprev"
    )
  )
  for (row in cases) {
    extra <- if (!is.null(row$p2_var)) list(p2_var = row$p2_var)
    cmp <- do.call(
      expect_comparisons_match_paired_sops,
      c(list(row$case, info = row$label), extra)
    )
    expect_identical(
      as.character(unique(cmp$time)),
      as.character(1:4),
      info = row$label
    )
    if (!is.null(row$p2_var)) {
      expect_identical(attr(cmp, "p2_var"), row$p2_var, info = row$label)
      expect_error(
        inferences(cmp, method = "delta"),
        "second-order",
        info = row$label
      )
    }
  }
})
