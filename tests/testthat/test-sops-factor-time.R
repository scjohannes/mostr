make_factor_visit_case <- function(n_patients = 80, n_visits = 4, seed = 2026) {
  withr::with_seed(seed, {
    visits <- as.character(seq_len(n_visits))
    ids <- seq_len(n_patients)
    tx <- stats::rbinom(n_patients, 1, 0.5)

    data <- expand.grid(
      id = ids,
      time = visits,
      KEEP.OUT.ATTRS = FALSE
    )
    data <- data[order(data$id, data$time), , drop = FALSE]
    data$time <- factor(data$time, levels = visits)
    data$tx <- tx[data$id]
    data$yprev <- factor(
      sample(seq_len(4), nrow(data), replace = TRUE),
      levels = as.character(seq_len(4))
    )
    eta <- -0.35 *
      data$tx +
      0.25 * as.integer(data$time) +
      0.2 * as.integer(as.character(data$yprev))
    p_worse <- stats::plogis(eta - mean(eta))
    y_num <- pmin(
      4L,
      pmax(1L, 1L + stats::rbinom(nrow(data), 3L, p_worse))
    )
    data$y <- ordered(y_num, levels = as.character(seq_len(4)))
    data
  })
}

make_interpolation_case <- function() {
  x <- expand.grid(
    time = factor(c("v1", "v2"), levels = c("v1", "v2")),
    state = factor(c("1", "2"), levels = c("1", "2")),
    tx = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  x$estimate <- mapply(
    function(time, state, tx) {
      if (tx == 0 && time == "v1" && state == "1") {
        return(0.8)
      }
      if (tx == 0 && time == "v1" && state == "2") {
        return(0.2)
      }
      if (tx == 0 && time == "v2" && state == "1") {
        return(0.4)
      }
      if (tx == 0 && time == "v2" && state == "2") {
        return(0.6)
      }
      if (tx == 1 && time == "v1" && state == "1") {
        return(0.7)
      }
      if (tx == 1 && time == "v1" && state == "2") {
        return(0.3)
      }
      if (tx == 1 && time == "v2" && state == "1") {
        return(0.5)
      }
      0.5
    },
    as.character(x$time),
    as.character(x$state),
    x$tx
  )

  class(x) <- c("markov_avg_sops", class(x))
  attr(x, "settings") <- markov_settings(
    p_var = "yprev",
    y_levels = factor(1:2)
  )
  attr(x, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = NULL,
    times = factor(c("v1", "v2"), levels = c("v1", "v2")),
    id_var = "id"
  )
  attr(x, "newdata_orig") <- data.frame(
    id = seq_len(4),
    yprev = factor(c("1", "2", "2", "2"), levels = c("1", "2")),
    tx = c(0, 1, 0, 1)
  )
  x
}

make_absorbing_factor_visit_data <- function(n_patients = 80, seed = 20260531) {
  visit_days <- stats::setNames(c(3, 7, 14, 21), as.character(1:4))
  raw_data <- sim_markov(
    "actt2",
    n_patients = n_patients,
    treatment_effect = -0.03,
    seed = seed
  )

  baseline_state <- raw_data[raw_data$time == 1, c("id", "yprev")]
  names(baseline_state)[2] <- "baseline_yprev"

  visit_data <- raw_data[
    raw_data$time %in% unname(visit_days),
    c("id", "tx", "time", "y", "yprev")
  ]
  visit_data <- merge(
    visit_data,
    baseline_state,
    by = "id",
    all.x = TRUE,
    sort = FALSE
  )
  visit_data <- visit_data[order(visit_data$id, visit_data$time), ]
  visit_data$visit_day <- visit_data$time
  visit_data$time <- factor(
    match(visit_data$visit_day, unname(visit_days)),
    levels = seq_along(visit_days),
    labels = names(visit_days)
  )
  visit_data$yprev <- ave(
    visit_data$y,
    visit_data$id,
    FUN = function(y) c(NA_integer_, y[-length(y)])
  )
  first_visit <- !duplicated(visit_data$id)
  visit_data$yprev[first_visit] <- visit_data$baseline_yprev[first_visit]
  visit_data$baseline_yprev <- NULL

  list(
    data = prepare_markov_data(
      visit_data,
      absorbing_state = 8L,
      factor_previous = TRUE
    ),
    visit_days = visit_days,
    absorb = "8"
  )
}

# An orm fit to `make_absorbing_factor_visit_data(n_patients = 80, seed = 2031)`,
# made once per test run.
absorbing_factor_visit_fit <- function() {
  grid_fixture("absorbing factor visit fit", function() {
    case <- make_absorbing_factor_visit_data(n_patients = 80, seed = 2031)
    case$fit <- mostm(
      y ~ tx + time + yprev,
      data = case$data,
      id_var = "id",
      first_followup_time = "1",
      opt_method = "LM",
      scale = TRUE,
      penalty = 0.1,
      maxit = 100
    )
    case
  })
}

add_interpolation_draws <- function(x) {
  draws <- x[rep(seq_len(nrow(x)), each = 3L), , drop = FALSE]
  draws$draw_id <- rep(seq_len(3), times = nrow(x))
  offset <- c(-0.1, 0, 0.1)[draws$draw_id]
  draws$estimate <- draws$estimate + ifelse(draws$state == "1", offset, -offset)
  draws$estimate <- pmin(pmax(draws$estimate, 0), 1)
  rownames(draws) <- NULL

  attr(x, "draws") <- draws
  attr(x, "conf_level") <- 0.8
  attr(x, "conf_type") <- "perc"
  x
}

add_interpolation_intervals <- function(x) {
  x$conf.low <- pmax(x$estimate - 0.05, 0)
  x$conf.high <- pmin(x$estimate + 0.05, 1)
  x$std.error <- 0.025
  x
}

test_that("factor visit gap requires explicit numeric time_covariates", {
  skip_if_not_installed("VGAM")

  data <- make_factor_visit_case(n_patients = 40)
  fit <- suppressWarnings(
    mostm(
      y ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      first_followup_time = "1",
      engine = "vglm"
    )
  )
  baseline <- data[!duplicated(data$id), , drop = FALSE]

  expect_error(
    soprob_markov(
      fit,
      newdata = baseline,
      times = 1:4,
      y_levels = factor(1:4),
      p_var = "yprev",
      gap_var = "gap"
    ),
    "Factor visit time with `gap_var` requires numeric gap_var values"
  )

  out <- soprob_markov(
    fit,
    newdata = baseline,
    times = 1:4,
    y_levels = factor(1:4),
    p_var = "yprev",
    gap_var = "gap",
    time_covariates = data.frame(gap = c(3, 4, 7, 14))
  )
  expect_equal(dim(out), c(nrow(baseline), 4L, 4L))
})

test_that("factor visit orm SOP interpolation uses absorbing state", {
  skip_if_not_installed("rms")

  case <- absorbing_factor_visit_fit()
  fit <- case$fit

  sop_visit <- avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = names(case$visit_days),
    y_levels = fit$yunique,
    absorb = case$absorb
  )
  expect_equal(
    as.character(attr(sop_visit, "avg_args")$times),
    names(case$visit_days)
  )
  expect_false(anyNA(sop_visit$estimate))

  sop_days <- interpolate_sops(
    sop_visit,
    time_map = case$visit_days,
    target_times = 0:max(case$visit_days),
    baseline_time = 0
  )
  expect_false(anyNA(sop_days$estimate))

  sums <- stats::aggregate(estimate ~ tx + time, sop_days, sum)
  expect_equal(sums$estimate, rep(1, nrow(sums)), tolerance = 1e-10)
})

test_that("interpolate_sops maps visits, anchors baseline, and normalizes states for each time grid", {
  x <- make_interpolation_case()
  time_map <- c(v1 = 3, v2 = 7)

  # `times` is the expected time grid and `baseline_time` the stored baseline
  # time. `values` holds hand-computed estimates; `same_as` names a case whose
  # result this one must equal.
  cases <- list(
    list(
      label = "explicit target times from baseline",
      args = list(
        time_map = time_map,
        target_times = c(0, 1, 3, 5, 7),
        baseline_time = 0
      ),
      times = c(0, 1, 3, 5, 7),
      baseline_time = 0,
      # Linear interpolation between the baseline anchor 0.25 on day 0 and
      # 0.8 on day 3, and between 0.8 on day 3 and 0.4 on day 7.
      values = data.frame(
        tx = c(0, 0),
        time = c(1, 5),
        state = c("1", "1"),
        estimate = c(0.25 + (0.8 - 0.25) / 3, 0.6)
      )
    ),
    list(
      label = "default target times",
      args = list(time_map = time_map),
      times = c(0, 3, 7),
      baseline_time = 0
    ),
    list(
      label = "no baseline",
      args = list(time_map = time_map, baseline_time = NULL),
      times = c(3, 7),
      baseline_time = NULL
    ),
    list(
      label = "unused time-map entry",
      args = list(time_map = c(unused = -10, time_map)),
      times = c(0, 3, 7),
      baseline_time = 0,
      same_as = "default target times"
    )
  )

  results <- list()
  for (case in cases) {
    out <- do.call(interpolate_sops, c(list(x), case$args))
    results[[case$label]] <- out
    expect_equal(sort(unique(out$time)), case$times, info = case$label)
    expect_equal(attr(out, "baseline_time"), case$baseline_time, info = case$label)

    for (i in seq_len(NROW(case$values))) {
      value <- case$values[i, ]
      row <- out[
        out$tx == value$tx & out$time == value$time & out$state == value$state,
      ]
      expect_equal(row$estimate, value$estimate, tolerance = 1e-10, info = case$label)
    }

    # The baseline anchor is the starting-state distribution of the stored
    # profiles (one of four starts in state 1), the same in every scenario.
    if (0 %in% case$times) {
      baseline <- out[out$time == 0, c("tx", "state", "estimate")]
      expect_equal(
        baseline[baseline$tx == 0, c("state", "estimate")],
        baseline[baseline$tx == 1, c("state", "estimate")],
        ignore_attr = TRUE,
        info = case$label
      )
      expect_equal(
        baseline$estimate[baseline$state == "1" & baseline$tx == 0],
        0.25,
        info = case$label
      )
      expect_equal(
        baseline$estimate[baseline$state == "2" & baseline$tx == 0],
        0.75,
        info = case$label
      )
    }

    sums <- stats::aggregate(estimate ~ tx + time, out, sum)
    expect_equal(
      sums$estimate,
      rep(1, nrow(sums)),
      tolerance = 1e-10,
      info = case$label
    )

    if (!is.null(case$same_as)) {
      expect_equal(out, results[[case$same_as]], ignore_attr = TRUE, info = case$label)
    }
  }
})

test_that("interpolate_sops gives individual profiles one-hot baselines", {
  x <- expand.grid(
    id = 1:2,
    time = factor(c("v1", "v2"), levels = c("v1", "v2")),
    state = factor(c("1", "2"), levels = c("1", "2")),
    KEEP.OUT.ATTRS = FALSE
  )
  x$yprev <- factor(ifelse(x$id == 1, "1", "2"), levels = c("1", "2"))
  x$estimate <- 0.5
  class(x) <- c("markov_sops", class(x))
  attr(x, "settings") <- markov_settings(
    p_var = "yprev",
    y_levels = factor(1:2)
  )

  out <- interpolate_sops(x, time_map = c(v1 = 3, v2 = 7))
  baseline <- out[out$time == 0, , drop = FALSE]

  expect_equal(
    baseline$estimate[baseline$id == 1 & baseline$state == "1"],
    1
  )
  expect_equal(
    baseline$estimate[baseline$id == 1 & baseline$state == "2"],
    0
  )
  expect_equal(
    baseline$estimate[baseline$id == 2 & baseline$state == "1"],
    0
  )
  expect_equal(
    baseline$estimate[baseline$id == 2 & baseline$state == "2"],
    1
  )

  # Without its starting state, a profile has no baseline anchor.
  no_start <- x
  no_start$yprev <- NULL
  expect_error(
    interpolate_sops(no_start, time_map = c(v1 = 3, v2 = 7)),
    "Previous-state variable `yprev` not found in SOP output.",
    fixed = TRUE
  )
  unanchored <- interpolate_sops(
    no_start,
    time_map = c(v1 = 3, v2 = 7),
    baseline_time = NULL
  )
  expect_equal(sort(unique(unanchored$time)), c(3, 7))
})

test_that("interpolate_sops takes percentile intervals from interpolated stored draws", {
  simulated <- add_interpolation_draws(make_interpolation_case())
  # blrm results store posterior draws next to their intervals.
  posterior <- add_interpolation_intervals(make_interpolation_case())
  attr(posterior, "draws") <- attr(simulated, "draws")
  attr(posterior, "method") <- "posterior"
  attr(posterior, "conf_level") <- 0.8

  cases <- list(
    list(label = "simulation draws", x = simulated),
    list(label = "posterior draws", x = posterior)
  )
  for (case in cases) {
    expect_warning(
      out <- interpolate_sops(
        case$x,
        time_map = c(v1 = 3, v2 = 7),
        target_times = 0:7,
        baseline_time = 0
      ),
      NA,
      info = case$label
    )
    draws <- attr(out, "draws")
    expect_true(is.data.frame(draws), info = case$label)
    expect_equal(sort(unique(draws$time)), 0:7, info = case$label)

    # Every draw shares the fixed baseline anchor, so the interval is a point.
    baseline <- out[out$time == 0 & out$tx == 0 & out$state == "1", ]
    expect_equal(baseline$estimate, 0.25, tolerance = 1e-10, info = case$label)
    expect_equal(baseline$conf.low, 0.25, tolerance = 1e-10, info = case$label)
    expect_equal(baseline$conf.high, 0.25, tolerance = 1e-10, info = case$label)
    expect_equal(baseline$std.error, 0, tolerance = 1e-10, info = case$label)

    day_1 <- out[out$time == 1 & out$tx == 0 & out$state == "1", ]
    expect_false(is.na(day_1$conf.low), info = case$label)
    expect_false(is.na(day_1$conf.high), info = case$label)
    expect_gt(day_1$conf.high, day_1$conf.low)
    day_1_draws <- draws[draws$time == 1 & draws$tx == 0 & draws$state == "1", ]
    expect_equal(
      day_1$conf.low,
      unname(stats::quantile(day_1_draws$estimate, probs = 0.1)),
      tolerance = 1e-10,
      info = case$label
    )
    expect_equal(
      day_1$conf.high,
      unname(stats::quantile(day_1_draws$estimate, probs = 0.9)),
      tolerance = 1e-10,
      info = case$label
    )
  }
})

test_that("interpolate_sops centers Wald intervals on interpolated estimates", {
  x <- make_interpolation_case()
  draws <- x[rep(seq_len(nrow(x)), each = 3L), , drop = FALSE]
  draws$draw_id <- rep(seq_len(3), times = nrow(x))
  # Unequal offsets put the draw mean above each state-1 estimate, so a
  # draw-mean center would move the interval.
  offset <- c(0, 0.1, 0.2)[draws$draw_id]
  draws$estimate <- draws$estimate + ifelse(draws$state == "1", offset, -offset)
  rownames(draws) <- NULL
  attr(x, "draws") <- draws
  attr(x, "conf_level") <- 0.95
  attr(x, "conf_type") <- "wald"

  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    target_times = 0:7,
    baseline_time = 0
  )

  expect_equal(
    (out$conf.low + out$conf.high) / 2,
    out$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    out$conf.high - out$estimate,
    stats::qnorm(0.975) * out$std.error,
    tolerance = 1e-10
  )
  expect_gt(max(out$std.error), 0)
})

test_that("interpolate_sops warns when interval-bearing objects do not store draws", {
  x <- add_interpolation_intervals(make_interpolation_case())
  expect_warning(
    out <- interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      target_times = 1:7,
      baseline_time = 0
    ),
    "no stored draws"
  )

  day_1 <- out[out$time == 1 & out$tx == 0 & out$state == "1", ]
  expect_false(is.na(day_1$conf.low))
  expect_false(is.na(day_1$conf.high))

  attr(x, "method") <- "posterior"
  expect_warning(
    interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      target_times = 1:7,
      baseline_time = 0
    ),
    "posterior draws"
  )
})

test_that("plot_sops bars work with derived facet labels on interpolated draws", {
  x <- add_interpolation_draws(make_interpolation_case())
  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    target_times = 1:7,
    baseline_time = 0
  )
  label_tx <- function(x) {
    x <- as.character(x)
    factor(
      ifelse(x %in% c("1", "Treatment"), "Treatment", "Control"),
      levels = c("Control", "Treatment")
    )
  }
  out$tx_label <- label_tx(out$tx)
  draws <- attr(out, "draws")
  draws$tx_label <- label_tx(draws$tx)
  attr(out, "draws") <- draws

  expect_equal(anyDuplicated(out[c("time", "state", "tx_label")]), 0L)
  expect_s3_class(
    plot_sops(out, geom = "bar", facet_var = "tx_label"),
    "ggplot"
  )
})

test_that("compiled interpolation preserves local support", {
  source_times <- c(2, 4, 7)
  target_times <- c(1, 2, 3, 7, 8)
  values <- c(0.2, 0.5, 0.9)
  plan <- mostr:::compile_linear_interpolation_plan(
    source_times,
    target_times
  )
  actual <- drop(mostr:::apply_linear_interpolation_plan(
    matrix(values, nrow = 1L),
    plan
  ))
  expected <- stats::approx(
    source_times,
    values,
    xout = target_times,
    rule = 1
  )$y

  expect_equal(actual, expected)
})

test_that("interpolate_sops does not extrapolate estimates or draws", {
  x <- data.frame(
    time = c("v2", "v3", "v1", "v2", "v3"),
    state = c("1", "1", "2", "2", "2"),
    estimate = c(0.4, 0.8, 0.7, 0.6, 0.2)
  )
  class(x) <- c("markov_sops", "data.frame")
  attr(x, "settings") <- markov_settings(y_levels = c("1", "2"))
  attr(x, "conf_level") <- 0.8
  attr(x, "draws") <- do.call(
    rbind,
    lapply(1:3, function(draw_id) {
      out <- x
      out$draw_id <- draw_id
      out$estimate <- out$estimate + (draw_id - 2) * 0.01
      out
    })
  )

  out <- interpolate_sops(
    x,
    time_map = c(v1 = 1, v2 = 2, v3 = 3),
    target_times = c(1, 1.5, 2, 2.5, 3),
    baseline_time = NULL,
    normalize = FALSE
  )
  draws <- attr(out, "draws")
  early <- out[out$state == "1" & out$time < 2, , drop = FALSE]
  middle <- out[out$state == "1" & out$time == 2.5, , drop = FALSE]
  early_draws <- draws[
    draws$state == "1" & draws$time < 2,
    ,
    drop = FALSE
  ]

  expect_equal(early$estimate, rep(NA_real_, nrow(early)))
  expect_equal(early$conf.low, rep(NA_real_, nrow(early)))
  expect_equal(early$conf.high, rep(NA_real_, nrow(early)))
  expect_equal(
    early_draws$estimate,
    rep(NA_real_, nrow(early_draws))
  )
  expect_equal(middle$estimate, 0.6, tolerance = 1e-12)
})

test_that("individual bootstrap weights propagate through baseline interpolation", {
  profiles <- data.frame(
    id = 1:2,
    rowid = 1:2,
    yprev = factor(c("1", "2"), levels = c("1", "2")),
    tx = c(0, 1)
  )
  x <- merge(
    expand.grid(
      rowid = 1:2,
      time = factor(c("v1", "v2"), levels = c("v1", "v2")),
      state = factor(c("1", "2"), levels = c("1", "2"))
    ),
    profiles,
    by = "rowid",
    sort = FALSE
  )
  x$estimate <- ifelse(x$state == "1", 0.6, 0.4)
  class(x) <- c("markov_sops", class(x))
  attr(x, "settings") <- markov_settings(
    p_var = "yprev",
    y_levels = factor(c("1", "2"))
  )
  attr(x, "newdata_orig") <- profiles
  attr(x, "newdata_pred") <- profiles

  draws <- do.call(
    rbind,
    lapply(1:2, function(draw_id) {
      draw <- x
      draw$draw_id <- draw_id
      draw$score_weight <- c(0.25, 1.75)[draw$rowid]
      draw
    })
  )
  attr(x, "draws") <- draws

  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    target_times = 0:7
  )
  out_draws <- attr(out, "draws")
  early <- out_draws[out_draws$time < 3, , drop = FALSE]

  expect_false(anyNA(early$estimate))
  expect_false(anyNA(early$score_weight))
  expect_equal(
    early$score_weight,
    c(0.25, 1.75)[early$rowid]
  )
  baseline <- early[early$time == 0, , drop = FALSE]
  expect_equal(
    baseline$estimate,
    as.numeric(as.character(baseline$state) == as.character(baseline$yprev))
  )
})

test_that("averaged bootstrap draws use draw-specific baseline distributions", {
  x <- add_interpolation_draws(make_interpolation_case())
  baseline <- attr(x, "newdata_orig")
  anchors <- lapply(
    list(c(3, 1, 1, 1), c(1, 1, 1, 3), c(1, 1, 1, 1)),
    function(weights) {
      mostr:::empirical_baseline_anchor_from_data(
        x,
        baseline = baseline,
        baseline_time = 0,
        weights = weights
      )
    }
  )
  attr(x, "baseline_anchor_draws") <-
    mostr:::combine_baseline_anchor_draws(anchors, 1:3)

  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    target_times = 0:7
  )
  draws <- attr(out, "draws")
  baseline_draws <- draws[
    draws$time == 0 & draws$tx == 0 & draws$state == "1",
    ,
    drop = FALSE
  ]

  expect_equal(baseline_draws$estimate, c(0.5, 1 / 6, 0.25))
  baseline_result <- out[out$time == 0 & out$tx == 0 & out$state == "1", ]
  expect_gt(baseline_result$std.error, 0)

  auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7),
    target_times = 0:7
  )
  expect_gt(auc$std.error[auc$tx == 0], 0)
})

test_that("interpolate_sops guardrails are clear", {
  x <- make_interpolation_case()
  expect_error(
    interpolate_sops(x, time_map = c(v1 = 3)),
    "`time_map` is missing entries"
  )
  expect_error(
    interpolate_sops(x, time_map = c(v1 = 3, v2 = 7), target_times = 8),
    "`target_times` must stay within the supported time range"
  )
  expect_error(
    interpolate_sops(x, time_map = c(v1 = 3, v2 = 7), baseline_time = Inf),
    "`baseline_time` must be a single finite numeric value or `NULL`"
  )
  expect_error(
    interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      target_times = -1
    ),
    "`target_times` must stay within the supported time range"
  )
  expect_error(
    interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      target_times = 1:7,
      baseline_time = NULL
    ),
    "`target_times` must stay within the supported time range"
  )
})

test_that("time_in_state restricts stored draws to the requested time range", {
  x <- add_interpolation_draws(make_interpolation_case())
  interpolated <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    target_times = 0:7
  )

  result <- time_in_state(interpolated, target_states = "1", target_times = 3:5)
  result <- result[order(result$tx), , drop = FALSE]
  draws <- attr(result, "draws")
  draws <- draws[order(draws$tx, draws$draw_id), , drop = FALSE]

  expect_equal(result$total_time, c(1.4, 1.3))
  expect_equal(draws$total_time, c(1.2, 1.4, 1.6, 1.1, 1.3, 1.5))
  expect_equal(result$conf.low, c(1.24, 1.14))
  expect_equal(result$conf.high, c(1.56, 1.46))
  expect_equal(result$std.error, c(0.2, 0.2))
})

test_that("time_in_state uses trapezoidal AUC on mapped real time", {
  x <- make_interpolation_case()

  auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7),
    baseline_time = NULL
  )
  auc <- auc[order(auc$tx), , drop = FALSE]

  expect_equal(auc$total_time[auc$tx == 0], 2.4, tolerance = 1e-10)
  expect_equal(auc$total_time[auc$tx == 1], 2.4, tolerance = 1e-10)

  baseline_auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7),
    baseline_time = 0,
    target_times = 0:7
  )
  baseline_auc <- baseline_auc[order(baseline_auc$tx), , drop = FALSE]

  expect_equal(
    baseline_auc$total_time[baseline_auc$tx == 0],
    3.975,
    tolerance = 1e-10
  )
  expect_equal(
    baseline_auc$total_time[baseline_auc$tx == 1],
    3.825,
    tolerance = 1e-10
  )
  default_auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7)
  )
  default_auc <- default_auc[order(default_auc$tx), , drop = FALSE]
  expect_equal(default_auc$total_time, c(3.975, 3.825), tolerance = 1e-10)

  day_1_auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7),
    baseline_time = 0,
    target_times = 1:7
  )
  day_1_auc <- day_1_auc[order(day_1_auc$tx), , drop = FALSE]

  expect_equal(
    day_1_auc$total_time[day_1_auc$tx == 0],
    3.633333,
    tolerance = 1e-6
  )
  expect_equal(day_1_auc$total_time[day_1_auc$tx == 1], 3.5, tolerance = 1e-10)

  interpolated <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    target_times = 1:7,
    baseline_time = 0
  )
  interpolated_auc <- time_in_state(interpolated, target_states = "1")
  interpolated_auc <- interpolated_auc[
    order(interpolated_auc$tx),
    ,
    drop = FALSE
  ]

  expect_equal(interpolated_auc, day_1_auc)
})
