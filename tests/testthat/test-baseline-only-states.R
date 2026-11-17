# A starting state that never occurs as an outcome, such as an admission-only
# screening state, is not one of the fitted states. The real-time summaries
# still start from it at `baseline_time`, where it holds its share of the
# patients, and give it probability 0 at every visit (issue #50).

baseline_only_case <- local({
  case <- NULL
  function() {
    if (!is.null(case)) {
      return(case)
    }
    visit_days <- c(3, 7, 14, 28)
    trajectories <- sim_markov(
      "actt2",
      n_patients = 200,
      times = visit_days,
      seed = 1
    )
    # State 3 occurs only as a starting state: outcome 3 becomes 4 at every
    # visit, previous state 3 becomes 4 after the first visit, and every
    # fourth patient starts in state 3.
    trajectories$y[trajectories$y == 3] <- 4
    later <- trajectories$time != visit_days[1L]
    trajectories$yprev[later & trajectories$yprev == 3] <- 4
    first <- trajectories$time == visit_days[1L]
    trajectories$yprev[first & trajectories$id %% 4 == 0] <- 3
    trajectories$visit <- factor(match(trajectories$time, visit_days))
    model_data <- prepare_markov_data(trajectories, absorb = 8)
    fit <- mostm(
      ordered(y) ~ yprev + visit + tx,
      data = model_data,
      id_var = "id",
      time_var = "visit",
      first_followup_time = "1",
      absorb = "8"
    )
    start <- as.character(model_data$yprev[model_data$visit == "1"])
    case <<- list(
      fit = fit,
      time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
      share = table(factor(start, levels = as.character(1:8))) / length(start)
    )
    case
  }
})

# Occupancy curves of avg_sops() output joined by straight lines from the
# starting shares at day 0, and the area under them on `grid` (trapezoids).
baseline_only_hand_totals <- function(avg, share, time_map, grid) {
  states <- names(share)
  knots <- c(0, unname(time_map))
  occupancy <- matrix(0, length(knots), length(states))
  colnames(occupancy) <- states
  occupancy[1L, ] <- as.numeric(share)
  for (v in seq_along(time_map)) {
    rows <- as.character(avg$time) == names(time_map)[v]
    occupancy[v + 1L, as.character(avg$state[rows])] <- avg$estimate[rows]
  }
  on_grid <- apply(occupancy, 2L, function(p) {
    stats::approx(knots, p, xout = grid)$y
  })
  colSums(diff(grid) * (on_grid[-1L, ] + on_grid[-nrow(on_grid), ]) / 2)
}

test_that("mostm() fits a starting state that never occurs as an outcome", {
  case <- baseline_only_case()
  expect_false("3" %in% as_state_labels(resolve_markov_settings(case$fit)$y_levels))
  expect_equal(unname(case$share[["3"]]), 0.25)
})

test_that("interpolated SOPs start from every starting state's share", {
  case <- baseline_only_case()
  avg <- avg_sops(case$fit, times = factor(1:4))
  out <- interpolate_sops(
    avg,
    time_map = case$time_map,
    target_times = c(0, 1, 2, 3, 7)
  )
  sums <- tapply(out$estimate, out$time, sum)
  expect_equal(as.numeric(sums), rep(1, 5))

  at_start <- out[out$time == 0, ]
  expected <- case$share[as.character(at_start$state)]
  expect_equal(at_start$estimate, as.numeric(expected))
  expect_equal(at_start$estimate[as.character(at_start$state) == "3"], 0.25)

  # Nobody is in state 3 after the first transition.
  expect_equal(out$estimate[as.character(out$state) == "3" & out$time == 3], 0)
  # Halfway to day 3 is still a third of the way through the starting share.
  expect_equal(
    out$estimate[as.character(out$state) == "3" & out$time == 1],
    0.25 * 2 / 3
  )

  # Without normalising, the probabilities already sum to one.
  raw <- interpolate_sops(
    avg,
    time_map = case$time_map,
    target_times = c(0, 1, 2, 3, 7),
    normalize = FALSE
  )
  expect_equal(as.numeric(tapply(raw$estimate, raw$time, sum)), rep(1, 5))
})

test_that("a patient's own SOPs start in a starting state that is not fitted", {
  case <- baseline_only_case()
  profile_sops <- sops(case$fit, times = factor(1:4))
  starts_in_3 <- unique(profile_sops$id[as.character(profile_sops$yprev) == "3"])
  one <- profile_sops[profile_sops$id == starts_in_3[1L], ]
  out <- interpolate_sops(one, time_map = case$time_map, target_times = c(0, 1, 3))
  expect_equal(as.numeric(tapply(out$estimate, out$time, sum)), rep(1, 3))
  expect_equal(out$estimate[as.character(out$state) == "3" & out$time == 0], 1)
  expect_equal(out$estimate[as.character(out$state) == "3" & out$time == 3], 0)
})

test_that("average time in states counts the time in a starting state that is not fitted", {
  case <- baseline_only_case()
  avg <- avg_sops(case$fit, times = factor(1:4))
  totals <- avg_state_time(case$fit, times = factor(1:4), time_map = case$time_map)
  hand <- baseline_only_hand_totals(avg, case$share, case$time_map, c(0, 3, 7, 14, 28))
  hand <- hand[hand > 0 | names(hand) %in% as.character(totals$state_set)]
  expect_setequal(as.character(totals$state_set), names(hand))
  expect_equal(
    totals$estimate,
    unname(hand[as.character(totals$state_set)]),
    tolerance = 1e-10
  )
  # A quarter of the patients spend half of the first three days in state 3.
  expect_equal(totals$estimate[as.character(totals$state_set) == "3"], 0.375)
  expect_equal(sum(totals$estimate), 28)

  # state_sets may name the starting state.
  pooled <- avg_state_time(
    case$fit,
    state_sets = list(early = c("3", "4")),
    times = factor(1:4),
    time_map = case$time_map
  )
  expect_equal(pooled$estimate, unname(hand[["3"]] + hand[["4"]]), tolerance = 1e-10)

  # Each patient's totals cover the whole follow-up.
  per_patient <- state_time(case$fit, times = factor(1:4), time_map = case$time_map)
  expect_equal(
    as.numeric(tapply(per_patient$estimate, per_patient$id, sum)),
    rep(28, length(unique(per_patient$id)))
  )
})

test_that("delta inference works for real-time totals with a starting state that is not fitted", {
  case <- baseline_only_case()
  totals <- avg_state_time(case$fit, times = factor(1:4), time_map = case$time_map)
  delta <- inferences(totals, method = "delta")
  expect_equal(delta$estimate, totals$estimate)
  expect_true(all(is.finite(delta$std.error)))
  conditional <- inferences(totals, method = "delta", vcov = "conditional")
  expect_true(all(is.finite(conditional$std.error)))

  # Target times between the starting state and the first visit.
  grid <- 0:28
  contrast <- avg_comparisons(
    case$fit,
    variables = list(tx = c(0, 1)),
    times = factor(1:4),
    estimand = "time_in_state",
    time_map = case$time_map,
    target_times = grid
  )
  arms <- avg_sops(case$fit, variables = list(tx = c(0, 1)), times = factor(1:4))
  hand <- baseline_only_hand_totals(arms[arms$tx == 1, ], case$share, case$time_map, grid) -
    baseline_only_hand_totals(arms[arms$tx == 0, ], case$share, case$time_map, grid)
  expect_equal(
    contrast$estimate,
    unname(hand[as.character(contrast$state_set)]),
    tolerance = 1e-10
  )
  contrast_delta <- inferences(contrast, method = "delta")
  expect_equal(contrast_delta$estimate, contrast$estimate)
  expect_true(all(is.finite(contrast_delta$std.error)))
})
