make_avg_time_case <- function(factor_time = FALSE) {
  skip_if_not_installed("VGAM")
  withr::local_seed(9321)
  data <- make_markov_trajectories(
    n_patients = 100,
    follow_up_time = 5,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = 1:2,
    absorbing_state = 3,
    treatment_effect = 0.3,
    seed = 9321
  ) |>
    prepare_markov_data(absorbing_state = 3)
  data$marker <- data$id %% 2
  data$subgroup <- ifelse(data$id %% 3 == 0, "a", "b")
  if (factor_time) {
    data$time <- factor(data$time, levels = 1:5)
  }
  model <- vglm_markov(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id"
  )
  list(model = model, data = data)
}

test_that("analytical inference clears bootstrap and previous null metadata", {
  case <- make_avg_time_case()
  boot <- avg_time(case$model, times = 1:3, absorb = 3) |>
    inferences(method = "fwb", n_draws = 3, seed = 7211, null = 0)
  result <- inferences(boot, method = "delta")
  expect_identical(attr(result, "method"), "delta")
  expect_null(attr(result, "fwb_weight_type"))
  expect_null(attr(result, "fwb_weight_scale"))
  expect_null(attr(result, "null"))
  expect_null(attr(result, "draws"))
  expect_length(
    intersect(c("statistic", "p.value", "s.value"), names(result)),
    0L
  )
  tested <- inferences(result, method = "delta", null = 1)
  expect_identical(attr(tested, "null"), 1)
  expect_equal(tested$statistic, (tested$estimate - 1) / tested$std.error)
})

test_that("avg_time returns observed-cohort totals and named state sets", {
  case <- make_avg_time_case()
  avg <- avg_sops(case$model, absorb = 3, times = 1:4)
  result <- avg_time(case$model, absorb = 3, times = 1:4)
  expect_s3_class(result, "markov_avg_time")
  expect_identical(inherits(result, "markov_avg_sops"), FALSE)
  expect_setequal(names(result), c("state_set", "estimate"))
  expect_setequal(result$state_set, as.character(1:3))
  for (state in as.character(1:3)) {
    expect_equal(
      result$estimate[result$state_set == state],
      time_in_state(avg, target_states = state)$total_time
    )
  }
  sets <- avg_time(
    case$model,
    absorb = 3,
    times = 1:4,
    state_sets = list(home = "1", alive = c("1", "2"), all = 1:3)
  )
  expect_equal(sets$estimate[sets$state_set == "all"], 4)
  expect_equal(
    sets$estimate[sets$state_set == "alive"],
    time_in_state(avg, target_states = 1:2)$total_time
  )
  pooled <- avg_time(case$model, absorb = 3, times = 1:4, state_sets = 1:2)
  expect_equal(pooled$estimate, sets$estimate[sets$state_set == "alive"])
})

test_that("avg_time preserves multiple scenarios and subgroup averaging", {
  case <- make_avg_time_case()
  profiles <- case$data[!duplicated(case$data$id), , drop = FALSE][1:15, ]
  variables <- list(tx = c(0, 1), marker = c(0, 1))
  avg <- avg_sops(
    case$model,
    absorb = 3,
    newdata = profiles,
    variables = variables,
    by = "subgroup",
    times = 1:3
  )
  result <- avg_time(
    case$model,
    absorb = 3,
    newdata = profiles,
    variables = variables,
    by = "subgroup",
    times = 1:3,
    state_sets = 1:2
  )
  expected <- time_in_state(avg, target_states = 1:2)
  joined <- merge(result, expected, by = c("tx", "marker", "subgroup"))
  expect_equal(nrow(result), 8L)
  expect_equal(nrow(joined), nrow(result))
  expect_equal(joined$estimate, joined$total_time)
  observed <- avg_time(
    case$model,
    absorb = 3,
    newdata = profiles,
    by = "subgroup",
    times = 1:3,
    state_sets = 1:2
  )
  observed_sops <- avg_sops(
    case$model,
    absorb = 3,
    newdata = profiles,
    by = "subgroup",
    times = 1:3
  )
  expected <- time_in_state(observed_sops, target_states = 1:2)
  joined <- merge(observed, expected, by = "subgroup")
  expect_equal(nrow(observed), 2L)
  expect_equal(joined$estimate, joined$total_time)
})

test_that("avg_time integrates irregular factor visits with optional baseline", {
  case <- make_avg_time_case(factor_time = TRUE)
  time_map <- c("1" = 2, "2" = 7, "3" = 9)
  avg <- avg_sops(case$model, absorb = 3, times = 1:3)
  # The omitted grid starts at `baseline_time = 0`, like `c(0, 2, 7, 9)`.
  for (target in list(NULL, c(0, 1, 2, 4, 7, 9), c(2, 7, 9))) {
    result <- avg_time(
      case$model,
      absorb = 3,
      times = 1:3,
      state_sets = list(alive = 1:2, all = 1:3),
      time_map = time_map,
      target_times = target,
      time_unit = "days"
    )
    expected <- time_in_state(
      avg,
      target_states = 1:2,
      time_map = time_map,
      target_times = target %||% c(0, 2, 7, 9)
    )
    expect_equal(
      result$estimate[result$state_set == "alive"],
      expected$total_time
    )
    expect_equal(
      result$estimate[result$state_set == "all"],
      if (identical(target, c(2, 7, 9))) 7 else 9
    )
    expect_equal(result$time_unit, rep("days", nrow(result)))
  }
  no_baseline <- avg_time(
    case$model,
    absorb = 3,
    times = 1:3,
    state_sets = 1:2,
    time_map = time_map,
    baseline_time = NULL
  )
  expect_equal(
    no_baseline$estimate,
    time_in_state(
      avg,
      1:2,
      time_map = time_map,
      baseline_time = NULL
    )$total_time
  )
  single <- avg_time(
    case$model,
    absorb = 3,
    times = 1,
    time_map = c("1" = 2),
    baseline_time = NULL
  )
  expect_equal(single$estimate, rep(0, 3))
  single_anchored <- avg_time(
    case$model,
    absorb = 3,
    times = 1,
    time_map = c("1" = 2)
  )
  expect_equal(sum(single_anchored$estimate), 2)
})

test_that("avg_time rejects grouping names that would overwrite output columns", {
  case <- make_avg_time_case()
  profiles <- case$data[!duplicated(case$data$id), , drop = FALSE][1:8, ]
  profiles$state_set <- profiles$subgroup
  profiles$time_unit <- profiles$subgroup
  expect_snapshot(error = TRUE, {
    avg_time(
      case$model,
      newdata = profiles,
      absorb = 3,
      times = 1:2,
      by = "state_set"
    )
  })
  expect_snapshot(error = TRUE, {
    avg_time(
      case$model,
      newdata = profiles,
      absorb = 3,
      times = 1:2,
      variables = list(state_set = c("a", "b"))
    )
  })
  expect_snapshot(error = TRUE, {
    avg_time(
      case$model,
      newdata = profiles,
      absorb = 3,
      times = 1:2,
      by = "time_unit",
      time_unit = "days"
    )
  })
  expect_snapshot(error = TRUE, {
    avg_time(
      case$model,
      newdata = profiles,
      absorb = 3,
      times = 1:2,
      variables = list(time_unit = c("a", "b")),
      time_unit = "days"
    )
  })
  unlabeled <- avg_time(
    case$model,
    newdata = profiles,
    absorb = 3,
    times = 1:2,
    by = "time_unit"
  )
  expect_setequal(unlabeled$time_unit, c("a", "b"))
})

test_that("all draw engines reduce matched SOP draws before calculating intervals", {
  case <- make_avg_time_case()
  avg <- avg_sops(
    case$model,
    absorb = 3,
    variables = list(tx = 0:1),
    times = 1:3
  )
  total <- avg_time(
    case$model,
    absorb = 3,
    variables = list(tx = 0:1),
    times = 1:3,
    state_sets = list(alive = 1:2, all = 1:3)
  )
  for (method in c("mvn", "score_bootstrap", "bootstrap", "fwb")) {
    inferred <- inferences(avg, method = method, n_draws = 6, seed = 9322)
    result <- inferences(total, method = method, n_draws = 6, seed = 9322)
    expected <- time_in_state(inferred, target_states = 1:2)
    joined <- merge(result[result$state_set == "alive", ], expected, by = "tx")
    expect_equal(joined$estimate, joined$total_time, info = method)
    expect_equal(joined$std.error.x, joined$std.error.y, info = method)
    expect_equal(joined$conf.low.x, joined$conf.low.y, info = method)
    expect_equal(joined$conf.high.x, joined$conf.high.y, info = method)
    draws <- get_draws(result)
    reduced <- attr(expected, "draws")
    actual <- draws[draws$state_set == "alive", ]
    actual <- actual[order(actual$draw_id, actual$tx), ]
    reduced <- reduced[order(reduced$draw_id, reduced$tx), ]
    expect_equal(actual$draw, reduced$total_time, info = method)
    expect_equal(
      draws$draw[draws$state_set == "all"],
      rep(3, sum(draws$state_set == "all")),
      tolerance = 1e-12
    )
    expect_equal(attr(result, "method"), method)
    expect_equal(attr(result, "n_successful"), length(unique(draws$draw_id)))
    expect_identical(
      any(c("time", "state") %in% names(attr(result, "draws"))),
      FALSE
    )
  }
})

test_that("mapped resampling uses the baseline anchor of each draw", {
  case <- make_avg_time_case()
  avg <- avg_sops(case$model, absorb = 3, times = 1:2)
  total <- avg_time(
    case$model,
    absorb = 3,
    times = 1:2,
    state_sets = "1",
    time_map = c("1" = 4, "2" = 7),
    target_times = c(0, 4, 7)
  )
  for (method in c("score_bootstrap", "bootstrap", "fwb")) {
    inferred <- inferences(avg, method = method, n_draws = 6, seed = 9323)
    result <- inferences(total, method = method, n_draws = 6, seed = 9323)
    anchors <- attr(inferred, "baseline_anchor_draws")
    anchors <- anchors[as.character(anchors$state) == "1", ]
    sop_draws <- attr(inferred, "draws")
    sop_draws <- sop_draws[as.character(sop_draws$state) == "1", ]
    first <- sop_draws[sop_draws$time == 1, ]
    last <- sop_draws[sop_draws$time == 2, ]
    first <- first[order(first$draw_id), ]
    last <- last[order(last$draw_id), ]
    anchors <- anchors[order(anchors$draw_id), ]
    expected <- 2 *
      anchors$estimate +
      3.5 * first$estimate +
      1.5 * last$estimate
    draws <- get_draws(result)
    draws <- draws[order(draws$draw_id), ]
    expect_equal(draws$draw, expected, info = method)
    expect_gt(stats::sd(anchors$estimate), 0)
  }
})

test_that("avg_time replaces inference metadata and preserves seeds and Wald centers", {
  case <- make_avg_time_case()
  total <- avg_time(case$model, absorb = 3, times = 1:3, state_sets = "1")
  withr::local_seed(9324)
  seed_before <- .Random.seed
  result <- inferences(
    total,
    n_draws = 8,
    seed = 1,
    conf_type = "wald",
    null = 0
  )
  expect_identical(.Random.seed, seed_before)
  expect_equal(result$estimate, total$estimate)
  expect_equal((result$conf.low + result$conf.high) / 2, total$estimate)
  expect_equal(
    result$p.value,
    2 * stats::pnorm(-abs(result$estimate / result$std.error))
  )
  repeated <- inferences(
    result,
    n_draws = 8,
    seed = 1,
    conf_type = "wald",
    null = 0
  )
  expect_equal(get_draws(repeated), get_draws(result))
  dropped <- inferences(result, n_draws = 5, seed = 2, return_draws = FALSE)
  expect_null(attr(dropped, "draws"))
  expect_equal(attr(dropped, "n_draws"), 5)
  expect_identical(
    any(c("statistic", "p.value", "null") %in% names(dropped)),
    FALSE
  )
})

test_that("posterior medians and intervals are calculated after time reduction", {
  coefficients <- rbind(
    c(0.5, -0.5, 0, 1.2),
    c(2, 1, 0, -1.1),
    c(-0.4, -2, 0, 0.2)
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
  avg <- avg_sops(
    model,
    newdata = profiles,
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "median",
    return_draws = TRUE
  )
  result <- avg_time(
    model,
    newdata = profiles,
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "median",
    return_draws = TRUE
  )
  sop_draws <- attr(avg, "draws")
  reduced <- stats::aggregate(estimate ~ draw_id + state, sop_draws, sum)
  for (state in as.character(1:3)) {
    values <- reduced$estimate[as.character(reduced$state) == state]
    row <- result[result$state_set == state, ]
    expect_equal(row$estimate, stats::median(values))
    expect_equal(row$std.error, stats::sd(values))
    expect_equal(row$conf.low, unname(stats::quantile(values, 0.025)))
    expect_equal(row$conf.high, unname(stats::quantile(values, 0.975)))
  }
  naive <- stats::aggregate(estimate ~ state, avg, sum)
  joined <- merge(result, naive, by.x = "state_set", by.y = "state")
  expect_gt(max(abs(joined$estimate.x - joined$estimate.y)), 1e-5)
  expect_equal(nrow(get_draws(result)), 9L)
  expect_identical(inferences(result, method = "mvn", n_draws = 2), result)
  compact <- avg_time(
    model,
    newdata = profiles,
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "mean",
    return_draws = FALSE
  )
  expect_null(attr(compact, "draws"))
  means <- stats::aggregate(estimate ~ state, reduced, mean)
  joined <- merge(compact, means, by.x = "state_set", by.y = "state")
  expect_equal(joined$estimate.x, joined$estimate.y)
})

test_that("avg_time validates arguments before predicting", {
  case <- make_avg_time_case()
  local_mocked_bindings(
    avg_comparison_replay_avg_sops = function(...) {
      stop("Predictions should not run before argument validation.")
    }
  )
  expect_error(
    avg_time(case$model, times = c(1, 1, 2), absorb = 3),
    "`times` must not contain duplicate"
  )
  expect_error(
    avg_time(case$model, times = 1:2, absorb = 3, state_sets = list()),
    "`state_sets`"
  )
  expect_error(
    avg_time(
      case$model,
      times = 1:2,
      absorb = 3,
      state_sets = list(a = "1", a = "2")
    ),
    "`state_sets`"
  )
  expect_error(
    avg_time(case$model, times = 1:2, absorb = 3, by = "state_set"),
    "reserved output names"
  )
  expect_error(
    avg_time(case$model, times = 1:2, absorb = 3, baseline_time = 1),
    "`baseline_time` requires `time_map`"
  )
})

test_that("avg_comparisons rejects duplicate times", {
  case <- make_avg_time_case()
  expect_error(
    avg_comparisons(
      case$model,
      variables = list(tx = c(0, 1)),
      estimand = "time_in_state",
      times = c(1, 1, 2),
      absorb = 3
    ),
    "`times` must not contain duplicate"
  )
})

test_that("posterior average times keep the requested state-set order", {
  coefficients <- rbind(
    c(0.5, -0.5, 0, 1.2),
    c(2, 1, 0, -1.1),
    c(-0.4, -2, 0, 0.2)
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
  result <- avg_time(
    model,
    newdata = profiles,
    variables = list(tx = c(1, 0)),
    times = 1:3,
    n_draws = 3,
    seed = 9326,
    state_sets = list(zeta = "3", alpha = "1")
  )
  expect_identical(result$state_set, c("zeta", "zeta", "alpha", "alpha"))
  expect_identical(result$tx, c(0, 1, 0, 1))
  comparison <- avg_comparisons(
    model,
    newdata = profiles,
    variables = list(tx = c(1, 0)),
    estimand = "time_in_state",
    times = 1:3,
    n_draws = 3,
    seed = 9326,
    state_sets = list(zeta = "3", alpha = "1")
  )
  expect_identical(comparison$state_set, c("zeta", "alpha"))
})
