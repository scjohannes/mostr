# Posterior summaries of totals and interpolations (issue #40). With
# `posterior_summary = "median"`, a total or an interpolation must be computed
# in every posterior draw before the median is taken: the median of a sum is
# not the sum of the medians, and medians of probabilities do not sum to one.
# With the default "mean" both orders give the same result.

median_fixture_model <- function() {
  # Five draws with decreasing intercepts, so that no cut-points cross.
  coefficients <- rbind(
    c(0.5, -0.5, 0, 1.2),
    c(2, 1, 0, -1.1),
    c(-0.4, -2, 0, 0.2),
    c(1.1, 0.3, 0, -0.6),
    c(-0.2, -1, 0, 0.9)
  )
  colnames(coefficients) <- c("y>=2", "y>=3", "tx", "time")
  fake_blrm(coefficients)
}

median_fixture_profiles <- function() {
  data.frame(
    id = c("a", "b"),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1
  )
}

median_fixture_args <- function(summary) {
  list(
    median_fixture_model(),
    newdata = median_fixture_profiles(),
    times = 1:3,
    n_draws = 5,
    seed = 1,
    posterior_summary = summary
  )
}

for (summary in c("median", "mean")) {
  test_that(paste0(
    "time_in_state() of posterior SOPs agrees with avg_state_time() and ",
    "state_time() (posterior_summary = \"", summary, "\")"
  ), {
    local_mocked_bindings(
      blrm_design_matrix = fake_blrm_design(c("tx", "time"))
    )
    args <- median_fixture_args(summary)
    avg <- do.call(avg_sops, c(args, return_draws = TRUE))
    individual <- do.call(sops, c(args, return_draws = TRUE))
    average_time <- do.call(avg_state_time, args)
    profile_time <- do.call(state_time, args)
    columns <- c("estimate", "conf.low", "conf.high", "std.error")

    for (state in as.character(1:3)) {
      total <- time_in_state(avg, state_sets = state)
      expected <- average_time[average_time$state_set == state, ]
      expect_equal(total$total_time, expected$estimate, info = state)
      expect_equal(
        unlist(total[, c("total_time", columns[-1])]),
        unlist(expected[, columns]),
        ignore_attr = TRUE,
        info = state
      )

      totals <- time_in_state(individual, state_sets = state)
      expected <- profile_time[profile_time$state_set == state, ]
      expected <- expected[match(totals$id, expected$id), ]
      expect_equal(totals$total_time, expected$estimate, info = state)
      expect_equal(totals$conf.low, expected$conf.low, info = state)
      expect_equal(totals$conf.high, expected$conf.high, info = state)
    }
  })

  test_that(paste0(
    "time_in_state(time_map = ) of posterior SOPs agrees with ",
    "avg_state_time(time_map = ) (posterior_summary = \"", summary, "\")"
  ), {
    local_mocked_bindings(
      blrm_design_matrix = fake_blrm_design(c("tx", "time"))
    )
    args <- median_fixture_args(summary)
    day_of_visit <- c("1" = 2, "2" = 5, "3" = 9)
    avg <- do.call(avg_sops, c(args, return_draws = TRUE))
    average_time <- do.call(
      avg_state_time,
      c(args, list(time_map = day_of_visit))
    )
    for (state in as.character(1:3)) {
      total <- time_in_state(avg, state_sets = state, time_map = day_of_visit)
      expected <- average_time[average_time$state_set == state, ]
      expect_equal(total$total_time, expected$estimate, info = state)
      expect_equal(total$conf.low, expected$conf.low, info = state)
      expect_equal(total$conf.high, expected$conf.high, info = state)
    }
  })

  test_that(paste0(
    "interpolate_sops() with an identity time_map returns the posterior ",
    "estimates (posterior_summary = \"", summary, "\")"
  ), {
    local_mocked_bindings(
      blrm_design_matrix = fake_blrm_design(c("tx", "time"))
    )
    avg <- do.call(
      avg_sops,
      c(median_fixture_args(summary), return_draws = TRUE)
    )
    identity_map <- c("1" = 1, "2" = 2, "3" = 3)
    for (normalize in c(TRUE, FALSE)) {
      interpolated <- interpolate_sops(
        avg,
        time_map = identity_map,
        baseline_time = NULL,
        normalize = normalize
      )
      key <- paste(avg$time, avg$state)
      row <- match(key, paste(interpolated$time, interpolated$state))
      for (column in c("estimate", "conf.low", "conf.high", "std.error")) {
        expect_equal(
          interpolated[[column]][row],
          avg[[column]],
          info = paste(column, normalize)
        )
      }
    }
  })
}

test_that("interpolated posterior medians are medians of the interpolated draws", {
  # Between visits, the median of interpolated draws is not the interpolated
  # median, so each draw is interpolated first.
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time"))
  )
  avg <- do.call(
    avg_sops,
    c(median_fixture_args("median"), return_draws = TRUE)
  )
  day_of_visit <- c("1" = 2, "2" = 5, "3" = 9)
  for (normalize in c(TRUE, FALSE)) {
    interpolated <- interpolate_sops(
      avg,
      time_map = day_of_visit,
      target_times = c(0, 1, 3.5, 7, 9),
      normalize = normalize
    )
    draws <- attr(interpolated, "draws")
    medians <- stats::aggregate(estimate ~ time + state, draws, stats::median)
    row <- match(
      paste(interpolated$time, interpolated$state),
      paste(medians$time, medians$state)
    )
    expect_equal(
      interpolated$estimate,
      medians$estimate[row],
      info = paste("normalize", normalize)
    )
    expect_true(all(
      interpolated$estimate >= interpolated$conf.low - 1e-12 &
        interpolated$estimate <= interpolated$conf.high + 1e-12
    ))
  }
})

test_that("a posterior-median total lies inside its own interval", {
  # Three visits and three draws: in draw d, every patient is in state 1 at
  # visit d only. Each draw spends exactly one visit in state 1, so the total
  # is 1 in every draw. The median probability of state 1 is 0 at every visit,
  # so adding up the visit-wise medians would give a total of 0, outside the
  # interval [1, 1].
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time"))
  )
  args <- median_fixture_args("median")
  args$n_draws <- 3
  avg <- do.call(avg_sops, c(args, return_draws = TRUE))

  draws <- attr(avg, "draws")
  in_state_one <- as.numeric(draws$time == draws$draw_id)
  draws$estimate <- ifelse(
    draws$state == 1,
    in_state_one,
    ifelse(draws$state == 2, 1 - in_state_one, 0)
  )
  attr(avg, "draws") <- draws
  medians <- stats::aggregate(estimate ~ time + state, draws, stats::median)
  avg$estimate <- medians$estimate[match(
    paste(avg$time, avg$state),
    paste(medians$time, medians$state)
  )]

  total <- time_in_state(avg, state_sets = 1)
  expect_equal(total$total_time, 1)
  expect_equal(total$conf.low, 1)
  expect_equal(total$conf.high, 1)
  expect_gte(total$total_time, total$conf.low)
  expect_lte(total$total_time, total$conf.high)
})

test_that("plot_time(scale = \"proportion\") divides by the follow-up length", {
  # With posterior medians, the state totals need not add up to the
  # follow-up, so their sum is the wrong denominator for a share of it.
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time"))
  )
  args <- median_fixture_args("median")
  average_time <- do.call(avg_state_time, args)
  # The fixture's medians add up to less than the 3 visits.
  expect_false(isTRUE(all.equal(sum(average_time$estimate), 3)))
  plotted <- ggplot2::layer_data(plot_time(average_time, scale = "proportion"))
  built <- plot_time(average_time, scale = "proportion")$data
  expect_equal(
    built$estimate[match(average_time$state_set, built$state_set)],
    average_time$estimate / 3
  )
  expect_equal(
    built$conf.low[match(average_time$state_set, built$state_set)],
    average_time$conf.low / 3
  )
  expect_true(nrow(plotted) > 0L)

  # With time_map, follow-up runs from baseline_time = 0 to day 9.
  day_of_visit <- c("1" = 2, "2" = 5, "3" = 9)
  real_time <- do.call(avg_state_time, c(args, list(time_map = day_of_visit)))
  built <- plot_time(real_time, scale = "proportion")$data
  expect_equal(
    built$estimate[match(real_time$state_set, built$state_set)],
    real_time$estimate / 9
  )

  # Posterior means add up to the follow-up, so the plot is unchanged.
  mean_time <- do.call(avg_state_time, median_fixture_args("mean"))
  expect_equal(sum(mean_time$estimate), 3)
  built <- plot_time(mean_time, scale = "proportion")$data
  expect_equal(
    built$estimate[match(mean_time$state_set, built$state_set)],
    mean_time$estimate / 3
  )
})
