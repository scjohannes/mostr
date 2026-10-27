# Each profile's SOPs for a state set, summed over visits, from `sops()` on
# the profiles with `var` set to `value`. Returned in profile order.
profile_visit_totals <- function(case, var, value, states, ...) {
  ind <- sops_with_value(case, var, value, ...)
  keep <- as.character(ind$state) %in% states
  totals <- stats::aggregate(estimate ~ rowid, ind[keep, ], sum)
  totals$estimate[order(totals$rowid)]
}

test_that("state_time() totals each profile's SOPs under each value", {
  case <- make_comparisons_case()
  n <- nrow(case$profiles)
  out <- state_time(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    state_sets = list(alive = c("1", "2")),
    times = 1:4,
    absorb = 3
  )
  expect_s3_class(out, "markov_state_time")
  expect_identical(nrow(out), 2L * n)
  expect_identical(out$rowid, seq_len(2L * n))
  expect_identical(out$rowidcf, rep(seq_len(n), 2L))
  expect_true(all(out$state_set == "alive"))
  for (value in c(0, 1)) {
    rows <- out[out$tx == value, , drop = FALSE]
    expect_equal(
      rows$estimate,
      profile_visit_totals(case, "tx", value, c("1", "2"))
    )
    expect_equal(rows$marker, case$profiles$marker)
  }
  expect_false(any(startsWith(names(out), ".markov_")))
  expect_false("time" %in% names(out))
  expect_identical(names(out)[ncol(out)], "estimate")

  plain <- state_time(case$model, times = 1:4, absorb = 3)
  expect_false("rowidcf" %in% names(plain))
  expect_setequal(plain$state_set, as.character(1:3))
})

test_that("averaging state_time() reproduces avg_state_time()", {
  case <- make_comparisons_case()
  starts <- factor(c("1", "2"), levels = levels(case$profiles$yprev))
  settings <- list(
    list(variables = NULL, time_map = NULL, state_sets = NULL),
    list(variables = list(tx = c(0, 1)), time_map = NULL, state_sets = "1"),
    list(
      variables = list(tx = c(0, 1)),
      time_map = case$time_map,
      state_sets = list(alive = c("1", "2"), dead = "3")
    ),
    list(variables = list(yprev = starts), time_map = case$time_map, state_sets = "1")
  )
  for (setting in settings) {
    args <- c(
      list(case$model, times = 1:4, absorb = 3),
      setting
    )
    unit <- do.call(state_time, args)
    average <- do.call(avg_state_time, args)
    keys <- intersect(c("state_set", names(setting$variables)), names(average))
    means <- stats::aggregate(unit$estimate, unit[, keys, drop = FALSE], mean)
    index <- match(
      do.call(paste, average[, keys, drop = FALSE]),
      do.call(paste, means[, keys, drop = FALSE])
    )
    expect_false(anyNA(index))
    expect_equal(means$x[index], average$estimate, tolerance = 1e-10)
  }
})

test_that("state_time() validates its arguments before predicting", {
  case <- make_comparisons_case()
  local_mocked_bindings(
    sops = function(...) stop("Predictions should not run.")
  )
  clash <- case$profiles
  clash$state_set <- "x"
  expect_error(
    state_time(case$model, newdata = clash, times = 1:4),
    "`state_set`"
  )
  expect_error(state_time(case$model), "`times` must be supplied")
  expect_error(
    state_time(case$model, times = 1:4, baseline_time = 1),
    "`baseline_time` requires `time_map`"
  )
})

test_that("state_time() draws reduce to avg_state_time() draws", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:3,
    absorb = 3,
    time_map = c("1" = 1, "2" = 2, "3" = 4)
  )
  for (setting in list(
    list(method = "mvn", n_draws = 4L, seed = 601, weight = NULL),
    list(method = "score_bootstrap", n_draws = 4L, seed = 602, weight = "score_weight"),
    list(method = "fwb", n_draws = 2L, seed = 603, weight = "fwb_weight")
  )) {
    unit <- inferences(
      do.call(state_time, args),
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    )
    expect_s3_class(unit, "markov_state_time")
    expect_identical(attr(unit, "method"), setting$method)
    expect_null(attr(unit, "null"))
    average <- inferences(
      do.call(avg_state_time, args),
      method = setting$method,
      n_draws = setting$n_draws,
      seed = setting$seed
    )
    draws <- get_draws(unit)
    weights <- if (is.null(setting$weight)) {
      rep(1, nrow(draws))
    } else {
      draws[[setting$weight]]
    }
    key <- paste(draws$draw_id, draws$tx)
    manual <- c(tapply(draws$draw * weights, key, sum)) /
      c(tapply(weights, key, sum))
    average_draws <- attr(average, "draws")
    expect_equal(
      unname(manual[paste(average_draws$draw_id, average_draws$tx)]),
      average_draws$estimate,
      tolerance = 1e-10
    )
  }
  expect_error(
    inferences(do.call(state_time, args), method = "bootstrap", n_draws = 2),
    "`method = \"fwb\"`"
  )
})

test_that("posterior state_time() summarizes each profile's total per draw", {
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
    markov_engine = "blrm"
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
  out <- state_time(
    model,
    newdata = profiles,
    state_sets = "1",
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "median",
    return_draws = TRUE
  )
  ind <- sops(
    model,
    newdata = profiles,
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    return_draws = TRUE
  )
  draws <- attr(ind, "draws")
  totals <- stats::aggregate(
    estimate ~ rowid + draw_id,
    draws[as.character(draws$state) == "1", ],
    sum
  )
  for (i in seq_len(nrow(out))) {
    values <- totals$estimate[totals$rowid == out$rowid[i]]
    expect_equal(out$estimate[i], stats::median(values))
    expect_equal(out$std.error[i], stats::sd(values))
  }
  expect_identical(nrow(get_draws(out)), 3L * nrow(out))
})
