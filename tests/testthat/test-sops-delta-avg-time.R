delta_avg_time_case <- local({
  cases <- list()
  function(factor_time = FALSE) {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    key <- as.character(factor_time)
    if (is.null(cases[[key]])) {
      data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 9461)
      if (factor_time) {
        data$time <- factor(data$time)
      }
      fit <- suppressWarnings(vglm_markov(
        ordered(y) ~ time + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        first_followup_time = if (factor_time) levels(data$time)[1L] else 1
      ))
      cases[[key]] <<- list(
        fit = fit,
        times = if (factor_time) levels(data$time)[1:3] else 1:3,
        time_map = c("1" = 3, "2" = 7, "3" = 14)
      )
    }
    cases[[key]]
  }
})

avg_time_numerical_jacobian <- function(model, point) {
  beta <- get_coef(model)
  out <- matrix(
    NA_real_,
    nrow = nrow(point(model)),
    ncol = length(beta),
    dimnames = list(NULL, names(beta))
  )
  for (j in seq_along(beta)) {
    step <- .Machine$double.eps^(1 / 3) * max(1, abs(beta[j]))
    upper <- lower <- beta
    upper[j] <- upper[j] + step
    lower[j] <- lower[j] - step
    out[, j] <- (point(set_coef(model, upper))$estimate -
      point(set_coef(model, lower))$estimate) /
      (2 * step)
  }
  out
}

test_that("average-time conditional variance agrees with numerical derivatives", {
  case <- delta_avg_time_case(factor_time = TRUE)
  for (real_time in c(FALSE, TRUE)) {
    point <- function(model) {
      avg_time(
        model,
        variables = list(tx = c(0, 1)),
        times = case$times,
        state_sets = list(lower = c("1", "2"), all = as.character(1:6)),
        absorb = 6,
        time_map = if (real_time) case$time_map else NULL,
        target_times = if (real_time) c(0, 1, 4, 9, 14) else NULL
      )
    }
    raw <- point(case$fit)
    inferred <- inferences(raw, method = "delta", vcov = case$fit$var)
    numerical <- avg_time_numerical_jacobian(case$fit, point)
    covariance <- case$fit$var[colnames(numerical), colnames(numerical)]
    expected <- numerical %*% covariance %*% t(numerical)

    expect_equal(inferred$estimate, raw$estimate, tolerance = 1e-12)
    expect_equal(
      unname(get_jacobian(inferred)),
      unname(numerical),
      tolerance = 2e-6
    )
    expect_equal(
      unname(stats::vcov(inferred)),
      unname(expected),
      tolerance = 2e-6
    )
    expect_equal(
      inferred$std.error,
      sqrt(pmax(diag(expected), 0)),
      tolerance = 2e-6
    )
    all <- inferred$state_set == "all"
    expect_equal(
      inferred$estimate[all],
      rep(if (real_time) 14 else 3, sum(all))
    )
    expect_lt(max(inferred$std.error[all]), 1e-10)
    expect_identical(attr(inferred, "conf_type"), "wald")
    expect_null(attr(inferred, "draws"))
  }
})

test_that("unconditional average time includes baseline sampling and joint covariance", {
  case <- delta_avg_time_case()
  states <- list(
    lower = c("1", "2"),
    upper = as.character(3:6),
    all = as.character(1:6)
  )
  variables <- list(tx = c(0, 1))
  target_times <- c(0, 1, 4, 9, 14)
  point <- avg_time(
    case$fit,
    variables = variables,
    times = case$times,
    state_sets = states,
    absorb = 6,
    time_map = case$time_map,
    target_times = target_times
  )
  inferred <- inferences(point, method = "delta")
  avg <- avg_sops(
    case$fit,
    variables = variables,
    times = case$times,
    absorb = 6
  ) |>
    inferences(method = "delta", conf_type = "wald")
  source <- attr(avg, "analytical")
  baseline <- attr(avg, "newdata_pred")[
    seq_along(source$profile_ids),
    ,
    drop = FALSE
  ]

  nodes <- c(0, unname(case$time_map))
  weights <- vapply(
    seq_along(nodes),
    function(i) {
      values <- stats::approx(
        nodes,
        as.numeric(seq_along(nodes) == i),
        xout = target_times
      )$y
      sum(diff(target_times) * (head(values, -1L) + tail(values, -1L)) / 2)
    },
    numeric(1)
  )
  operator <- matrix(0, nrow(inferred), nrow(avg))
  anchor <- matrix(0, nrow(baseline), nrow(inferred))
  for (i in seq_len(nrow(inferred))) {
    selected <- as.character(avg$state) %in%
      states[[inferred$state_set[i]]] &
      avg$tx == inferred$tx[i]
    operator[i, ] <- selected * weights[-1L][match(avg$time, case$times)]
    indicator <- as.numeric(
      as.character(baseline$yprev) %in% states[[inferred$state_set[i]]]
    )
    anchor[, i] <- weights[1L] * (indicator - mean(indicator))
  }
  without_anchor <- source$influence %*% t(operator)
  influence <- without_anchor + anchor
  covariance <- stats::cov(influence) / nrow(influence)

  expect_equal(
    unname(attr(inferred, "analytical")$influence),
    unname(influence),
    tolerance = 1e-12
  )
  expect_equal(
    unname(stats::vcov(inferred)),
    unname(covariance),
    tolerance = 1e-12
  )
  expect_equal(
    inferred$std.error^2,
    unname(diag(covariance)),
    tolerance = 1e-12
  )
  expect_gt(
    max(abs(covariance - stats::cov(without_anchor) / nrow(influence))),
    1e-4
  )
  expect_lt(max(inferred$std.error[inferred$state_set == "all"]), 1e-10)
  selected_rows <- c(1L, 3L)
  expect_equal(
    unname(stats::vcov(inferred, rows = selected_rows)),
    unname(covariance[selected_rows, selected_rows]),
    tolerance = 1e-12
  )
})

test_that("average-time contrasts agree with analytical average comparisons", {
  case <- delta_avg_time_case()
  args <- list(
    model = case$fit,
    variables = list(tx = c(0, 1)),
    times = case$times,
    state_sets = list(lower = c("1", "2")),
    absorb = 6,
    time_map = case$time_map,
    target_times = 0:14
  )
  point <- do.call(avg_time, args)
  comparison <- do.call(
    avg_comparisons,
    c(args, list(estimand = "time_in_state"))
  )
  for (variance in c("conditional", "unconditional")) {
    inferred <- inferences(point, method = "delta", vcov = variance)
    contrasted <- inferences(comparison, method = "delta", vcov = variance)
    contrast <- as.numeric(inferred$tx == 1) - as.numeric(inferred$tx == 0)
    expect_equal(
      drop(contrast %*% inferred$estimate),
      contrasted$estimate,
      tolerance = 1e-12
    )
    expect_equal(
      drop(contrast %*% stats::vcov(inferred) %*% contrast),
      contrasted$std.error^2,
      tolerance = 1e-12
    )
  }
})

test_that("real-time covariance respects omitted anchors and zero-length intervals", {
  case <- delta_avg_time_case()
  # A grid of mapped visit times only never reaches the starting-state time, so
  # the result matches omitting the starting state altogether.
  anchored <- avg_time(
    case$fit,
    times = case$times,
    time_map = case$time_map,
    target_times = unname(case$time_map),
    absorb = 6
  ) |>
    inferences(method = "delta")
  unanchored <- avg_time(
    case$fit,
    times = case$times,
    time_map = case$time_map,
    baseline_time = NULL,
    absorb = 6
  ) |>
    inferences(method = "delta")
  expect_equal(anchored$estimate, unanchored$estimate, tolerance = 1e-12)
  expect_equal(
    stats::vcov(anchored),
    stats::vcov(unanchored),
    tolerance = 1e-12
  )

  point <- avg_time(
    case$fit,
    times = case$times,
    time_map = case$time_map,
    target_times = 7,
    absorb = 6
  )
  for (variance in c("conditional", "unconditional")) {
    inferred <- inferences(point, method = "delta", vcov = variance)
    expect_equal(inferred$estimate, rep(0, nrow(inferred)))
    expect_equal(inferred$std.error, rep(0, nrow(inferred)))
  }
})

test_that("unconditional observed-cohort times agree with patient-weight refits", {
  case <- delta_avg_time_case()
  point <- avg_time(
    case$fit,
    times = case$times,
    state_sets = list(lower = c("1", "2")),
    absorb = 6,
    time_map = case$time_map,
    target_times = 0:14
  )
  inferred <- inferences(point, method = "delta")
  analytical <- attr(inferred, "analytical")
  fit_data <- attr(case$fit, "markov_data")
  baseline <- attr(point, "newdata_pred")
  n <- nrow(baseline)
  patient_ids <- as.character(baseline$id[c(
    1L,
    which(as.character(baseline$yprev) %in% c("1", "2"))[1L]
  )])
  step <- 1e-2
  for (patient_id in unique(patient_ids)) {
    estimates <- vapply(
      c(-1, 1),
      function(direction) {
        fit_weights <- 1 +
          direction * step * (as.character(fit_data$id) == patient_id)
        fitted <- VGAM::vglm(
          ordered(y) ~ time + tx + yprev,
          family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
          data = fit_data,
          weights = fit_weights,
          coefstart = get_coef(case$fit),
          control = VGAM::vglm.control(epsilon = 1e-10, maxit = 100)
        )
        individual <- sops(
          set_coef(case$fit, stats::coef(fitted)),
          newdata = baseline,
          times = case$times,
          absorb = 6
        )
        total <- time_in_state(
          individual,
          target_states = c("1", "2"),
          time_map = case$time_map,
          target_times = 0:14
        )
        profile_weights <- 1 +
          direction * step * (as.character(total$id) == patient_id)
        stats::weighted.mean(total$total_time, profile_weights)
      },
      numeric(1)
    )
    numerical <- n * diff(estimates) / (2 * step)
    expected <- analytical$influence[
      match(patient_id, analytical$profile_ids),
      1L
    ]
    expect_equal(unname(numerical), unname(expected), tolerance = 2e-3)
  }
})

test_that("average-time analytical inference preserves restrictions and row alignment", {
  case <- delta_avg_time_case()
  point <- avg_time(case$fit, times = case$times, absorb = 6)
  inferred <- inferences(point, method = "delta")
  reordered <- point[nrow(point):1L, , drop = FALSE]
  result <- inferences(reordered, method = "delta")
  expect_equal(result$std.error, rev(inferred$std.error), tolerance = 1e-12)

  supplied <- avg_time(
    case$fit,
    newdata = attr(point, "newdata_pred"),
    times = case$times,
    absorb = 6
  )
  expect_snapshot(inferences(supplied, method = "delta"), error = TRUE)
  expect_no_error(inferences(supplied, method = "delta", vcov = "conditional"))
  expect_snapshot(
    inferences(point, method = "delta", conf_type = "logit"),
    error = TRUE
  )
  grouped <- avg_time(case$fit, times = case$times, absorb = 6, by = "tx")
  expect_snapshot(inferences(grouped, method = "delta"), error = TRUE)
})

starting_state_manual_time <- function(case, state, target_states, target_times) {
  baseline <- attr(
    avg_sops(case$fit, times = case$times, absorb = 6),
    "newdata_pred"
  )
  baseline$yprev[] <- as.character(state)
  individual <- sops(
    case$fit,
    newdata = baseline,
    times = case$times,
    absorb = 6
  )
  mean(time_in_state(
    individual,
    target_states = target_states,
    time_map = case$time_map,
    target_times = target_times
  )$total_time)
}

test_that("counterfactual starting states anchor the baseline at the set state", {
  case <- delta_avg_time_case()
  levels <- levels(attr(
    avg_sops(case$fit, times = case$times, absorb = 6),
    "newdata_pred"
  )$yprev)
  starts <- factor(c("1", "2"), levels = levels)
  result <- avg_time(
    case$fit,
    variables = list(yprev = starts),
    times = case$times,
    state_sets = "1",
    absorb = 6,
    time_map = case$time_map,
    target_times = 0:14
  )
  manual <- vapply(
    c("1", "2"),
    function(state) starting_state_manual_time(case, state, "1", 0:14),
    numeric(1)
  )
  expect_equal(
    result$estimate[match(c("1", "2"), as.character(result$yprev))],
    unname(manual),
    tolerance = 1e-10
  )

  comparison <- avg_comparisons(
    case$fit,
    variables = list(yprev = starts),
    times = case$times,
    state_sets = "1",
    absorb = 6,
    estimand = "time_in_state",
    time_map = case$time_map,
    target_times = 0:14
  )
  expect_equal(
    comparison$estimate,
    unname(manual[["2"]] - manual[["1"]]),
    tolerance = 1e-10
  )
})

test_that("counterfactual starting states allow only conditional variance", {
  case <- delta_avg_time_case()
  levels <- levels(attr(
    avg_sops(case$fit, times = case$times, absorb = 6),
    "newdata_pred"
  )$yprev)
  variables <- list(yprev = factor(c("1", "2"), levels = levels))
  point_fun <- function(model) {
    avg_time(
      model,
      variables = variables,
      times = case$times,
      state_sets = list(lower = c("1", "2")),
      absorb = 6,
      time_map = case$time_map,
      target_times = 0:14
    )
  }
  point <- point_fun(case$fit)
  inferred <- inferences(point, method = "delta", vcov = "conditional")
  numerical <- avg_time_numerical_jacobian(case$fit, point_fun)
  covariance <- stats::vcov(case$fit)[colnames(numerical), colnames(numerical)]
  expect_equal(inferred$estimate, point$estimate, tolerance = 1e-12)
  expect_equal(
    unname(stats::vcov(inferred)),
    unname(numerical %*% covariance %*% t(numerical)),
    tolerance = 2e-6
  )
  expect_no_error(inferences(point, method = "mvn", n_draws = 3, seed = 1))

  unconditional <- "starting state.*conditional"
  expect_error(inferences(point, method = "delta"), unconditional)
  expect_error(
    inferences(point, method = "delta", vcov = "unconditional"),
    unconditional
  )
  for (method in c("score_bootstrap", "bootstrap", "fwb")) {
    expect_error(
      inferences(point, method = method, n_draws = 3, seed = 1),
      unconditional,
      info = method
    )
  }
  sop_average <- avg_sops(
    case$fit,
    variables = variables,
    times = case$times,
    absorb = 6
  )
  expect_error(inferences(sop_average, method = "delta"), unconditional)
  expect_no_error(
    inferences(sop_average, method = "delta", vcov = "conditional")
  )
  comparison <- avg_comparisons(
    case$fit,
    variables = variables,
    times = case$times,
    state_sets = "1",
    absorb = 6,
    estimand = "time_in_state",
    time_map = case$time_map,
    target_times = 0:14
  )
  expect_error(inferences(comparison, method = "delta"), unconditional)
  expect_no_error(
    inferences(comparison, method = "delta", vcov = "conditional")
  )
})

test_that("interpolated SOP estimates and draws start at the set starting state", {
  case <- delta_avg_time_case()
  levels <- levels(attr(
    avg_sops(case$fit, times = case$times, absorb = 6),
    "newdata_pred"
  )$yprev)
  avg <- avg_sops(
    case$fit,
    variables = list(yprev = factor(c("1", "3"), levels = levels)),
    times = case$times,
    absorb = 6
  ) |>
    inferences(method = "mvn", n_draws = 4, seed = 1)
  interpolated <- interpolate_sops(avg, time_map = case$time_map)
  start <- interpolated[interpolated$time == 0, , drop = FALSE]
  expect_equal(
    start$estimate,
    as.numeric(as.character(start$state) == as.character(start$yprev))
  )
  draws <- attr(interpolated, "draws")
  start_draws <- draws[draws$time == 0, , drop = FALSE]
  expect_equal(nrow(start_draws), 4L * nrow(start))
  value <- if ("draw" %in% names(start_draws)) "draw" else "estimate"
  expect_equal(
    start_draws[[value]],
    as.numeric(
      as.character(start_draws$state) == as.character(start_draws$yprev)
    )
  )
})

test_that("mapped totals include the baseline interval by default", {
  case <- delta_avg_time_case()
  default <- avg_time(
    case$fit,
    variables = list(tx = c(0, 1)),
    times = case$times,
    state_sets = list(lower = c("1", "2"), all = as.character(1:6)),
    absorb = 6,
    time_map = case$time_map
  )
  explicit <- avg_time(
    case$fit,
    variables = list(tx = c(0, 1)),
    times = case$times,
    state_sets = list(lower = c("1", "2"), all = as.character(1:6)),
    absorb = 6,
    time_map = case$time_map,
    target_times = c(0, 3, 7, 14)
  )
  expect_equal(default$estimate, explicit$estimate, tolerance = 1e-12)
  expect_equal(default$estimate[default$state_set == "all"], c(14, 14))
  expect_equal(
    inferences(default, method = "delta")$std.error,
    inferences(explicit, method = "delta")$std.error,
    tolerance = 1e-12
  )
  comparison_args <- list(
    model = case$fit,
    variables = list(tx = c(0, 1)),
    times = case$times,
    state_sets = "1",
    absorb = 6,
    estimand = "time_in_state",
    time_map = case$time_map
  )
  expect_equal(
    do.call(avg_comparisons, comparison_args)$estimate,
    do.call(
      avg_comparisons,
      c(comparison_args, list(target_times = c(0, 3, 7, 14)))
    )$estimate,
    tolerance = 1e-12
  )
  unanchored <- avg_time(
    case$fit,
    times = case$times,
    state_sets = list(all = as.character(1:6)),
    absorb = 6,
    time_map = case$time_map,
    baseline_time = NULL
  )
  expect_equal(unanchored$estimate, 11)
})
