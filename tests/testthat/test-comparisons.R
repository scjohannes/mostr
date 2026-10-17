test_that("comparisons() returns paired SOP differences for each profile", {
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

test_that("averaging comparisons() reproduces avg_comparisons()", {
  case <- make_comparisons_case()
  variables <- list(tx = c(0, 1, 2))
  settings <- list(
    list(estimand = "sop", state_sets = list(alive = c("1", "2"), dead = "3")),
    list(estimand = "time_in_state", state_sets = NULL),
    list(
      estimand = "time_in_state",
      state_sets = "1",
      time_map = case$time_map
    ),
    list(estimand = "time_benefit", state_sets = NULL),
    list(estimand = "time_benefit", state_sets = NULL, time_map = case$time_map)
  )
  for (setting in settings) {
    args <- c(
      list(
        case$model,
        variables = variables,
        times = 1:4,
        absorb = 3
      ),
      setting
    )
    unit <- do.call(comparisons, args)
    average <- do.call(avg_comparisons, args)
    keys <- intersect(c("time", "state_set", "contrast"), names(average))
    expected <- stats::aggregate(
      unit$estimate,
      unit[, keys, drop = FALSE],
      mean
    )
    index <- match(
      do.call(paste, average[, keys, drop = FALSE]),
      do.call(paste, expected[, keys, drop = FALSE])
    )
    expect_false(anyNA(index))
    expect_equal(expected$x[index], average$estimate, tolerance = 1e-10)
    expect_identical(
      nrow(unit),
      nrow(average) * length(unique(unit$rowid))
    )
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

test_that("comparisons(by = ) hands off to avg_comparisons()", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    comparison = "ratio",
    times = 1:4,
    absorb = 3
  )
  expect_equal(
    do.call(comparisons, c(args, list(by = TRUE))),
    do.call(avg_comparisons, args)
  )
  expect_equal(
    do.call(comparisons, c(args, list(by = "marker"))),
    do.call(avg_comparisons, c(args, list(by = "marker")))
  )
  expect_s3_class(
    do.call(comparisons, c(args, list(by = NULL))),
    "markov_comparisons"
  )
  expect_error(
    do.call(comparisons, c(args, list(by = 1))),
    "`by` must be"
  )
  expect_error(
    do.call(comparisons, c(args, list(by = NA_character_))),
    "`by` must be"
  )
})

test_that("comparisons() validates arguments before predicting", {
  case <- make_comparisons_case()
  local_mocked_bindings(
    sops = function(...) {
      stop("Predictions should not run before argument validation.")
    }
  )
  clash <- case$profiles
  clash$contrast <- "a"
  clash$estimate <- 1
  expect_error(
    comparisons(
      case$model,
      newdata = clash,
      variables = list(tx = c(0, 1)),
      times = 1:4
    ),
    "`contrast`, `estimate`"
  )
  expect_error(
    comparisons(case$model, variables = list(tx = c(0, 1))),
    "`times` must be supplied to `comparisons\\(\\)`"
  )
  expect_error(
    comparisons(case$model, variables = list(tx = 1), times = 1:4),
    "at least two values"
  )
  expect_error(
    comparisons(
      case$model,
      variables = list(tx = c(0, 1), marker = c(0, 1)),
      times = 1:4
    ),
    "exactly one named variable"
  )
  expect_error(
    comparisons(
      case$model,
      variables = list(tx = c(0, 1)),
      times = c(1, 1, 2)
    ),
    "duplicate visits"
  )
  expect_error(
    comparisons(
      case$model,
      variables = list(tx = c(0, 1)),
      estimand = "time_benefit",
      comparison = "ratio",
      times = 1:4
    ),
    "only supports `comparison = \"difference\"`"
  )
})

test_that("posterior comparisons summarize paired draws for each profile", {
  coefficients <- rbind(
    c(0.5, -0.5, 0.4, 1.2),
    c(2, 1, -0.8, -1.1),
    c(-0.4, -2, 1.5, 0.2)
  )
  colnames(coefficients) <- c("y>=2", "y>=3", "tx", "time")
  model <- structure(
    list(
      draws = coefficients,
      non.slopes = 2L,
      pppo = 0L,
      tauInfo = data.frame(name = character()),
      ylevels = 1:3,
      yname = "y",
      clusterInfo = list(name = "id", cluster = c("a", "b"))
    ),
    class = c("blrm", "orm"),
    markov_fit_wrapper = "blrm_markov"
  )
  profiles <- data.frame(
    id = c("a", "b"),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1
  )
  local_mocked_bindings(
    blrm_design_matrix = function(model, newdata, second = FALSE) {
      cbind(tx = newdata$tx, time = newdata$time)
    }
  )
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
  inferred <- suppressWarnings(inferences(
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
      inferences(
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

test_that("averaging real-time comparison draws reproduces avg_comparisons()", {
  case <- make_comparisons_case()
  for (estimand in c("time_in_state", "time_benefit")) {
    args <- list(
      case$model,
      variables = list(tx = c(0, 1)),
      estimand = estimand,
      state_sets = if (estimand == "time_in_state") "1",
      times = 1:4,
      absorb = 3,
      time_map = case$time_map
    )
    unit <- inferences(
      do.call(comparisons, args),
      method = "mvn",
      n_draws = 6,
      seed = 202
    )
    average <- inferences(
      do.call(avg_comparisons, args),
      method = "mvn",
      n_draws = 6,
      seed = 202
    )
    unit_means <- stats::aggregate(estimate ~ draw_id, attr(unit, "draws"), mean)
    average_draws <- attr(average, "draws")
    expect_equal(
      unit_means$estimate,
      average_draws$estimate[match(unit_means$draw_id, average_draws$draw_id)],
      tolerance = 1e-10
    )
  }
})

test_that("weighted averages of profile draws reproduce avg_comparisons()", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:4,
    absorb = 3
  )
  for (setting in list(
    list(method = "score_bootstrap", n_draws = 4L, seed = 301),
    list(method = "fwb", n_draws = 3L, seed = 302)
  )) {
    weight_col <- if (setting$method == "fwb") "fwb_weight" else "score_weight"
    unit <- inferences(
      do.call(comparisons, args),
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    )
    average <- inferences(
      do.call(avg_comparisons, args),
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    )
    expect_true(attr(unit, "draw_weights_attached"))
    expect_identical(attr(unit, "draw_weight_col"), weight_col)
    draws <- attr(unit, "draws")
    expect_true(weight_col %in% names(draws))
    # One weight per patient and draw, shared by all of that patient's rows.
    n_profiles <- length(unique(draws$rowid))
    expect_identical(
      nrow(unique(draws[, c("draw_id", "rowid", weight_col)])),
      setting$n_draws * n_profiles
    )
    expect_true(weight_col %in% names(get_draws(unit)))

    key <- paste(draws$draw_id, draws$time)
    weights <- draws[[weight_col]]
    manual <- c(tapply(draws$estimate * weights, key, sum)) /
      c(tapply(weights, key, sum))
    average_draws <- attr(average, "draws")
    expect_equal(
      unname(manual[paste(average_draws$draw_id, average_draws$time)]),
      average_draws$estimate,
      tolerance = 1e-10
    )
  }
})

test_that("comparison intervals add null tests and replace earlier inference", {
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

test_that("comparison inference rejects unsupported methods and large draws", {
  case <- make_comparisons_case()
  cmp <- comparisons(
    case$model,
    variables = list(tx = c(0, 1)),
    times = 1:4,
    absorb = 3
  )
  expect_error(
    inferences(cmp, method = "bootstrap", n_draws = 2),
    "`method = \"fwb\"`"
  )
  withr::local_options(mostr.max_sops_draw_cells = 10)
  expect_error(
    inferences(cmp, n_draws = 5),
    "mostr.max_sops_draw_cells"
  )
})

# Checks one model setup end to end: point estimates against paired `sops()`
# differences, every draw method against paired `sops()` draws, and
# real-time totals against `avg_comparisons()`. `...` goes to every call.
expect_comparisons_match_paired_sops <- function(case, ...) {
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
  expect_identical(as.character(cmp$time), as.character(hi$time[keep]))
  expect_equal(cmp$estimate, hi$estimate[keep] - lo$estimate[keep])

  for (setting in list(
    list(method = "mvn", n_draws = 5L, seed = 401),
    list(method = "score_bootstrap", n_draws = 4L, seed = 402),
    list(method = "fwb", n_draws = 2L, seed = 403)
  )) {
    inferred <- suppressWarnings(inferences(
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
    expect_false(anyNA(index))
    expect_equal(draws$estimate, expected$estimate[index], tolerance = 1e-10)
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
  expect_equal(mean(unit$estimate), average$estimate, tolerance = 1e-10)
  unit_means <- stats::aggregate(estimate ~ draw_id, attr(unit, "draws"), mean)
  average_draws <- attr(average, "draws")
  expect_equal(
    unit_means$estimate,
    average_draws$estimate[match(unit_means$draw_id, average_draws$draw_id)],
    tolerance = 1e-10
  )
  invisible(cmp)
}

test_that("comparisons() works with factor visit times", {
  case <- make_factor_comparisons_case()
  cmp <- expect_comparisons_match_paired_sops(case)
  expect_identical(as.character(unique(cmp$time)), as.character(1:4))
})

test_that("comparisons() works with second-order models", {
  case <- make_second_order_comparisons_case()
  cmp <- expect_comparisons_match_paired_sops(case, p2_var = "ypprev")
  expect_identical(attr(cmp, "p2_var"), "ypprev")
  expect_error(inferences(cmp, method = "delta"), "second-order")
})
