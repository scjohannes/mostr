make_plot_time_model <- function() {
  skip_if_not_installed("VGAM")
  grid_fixture("plot_time_vglm", function() {
    withr::local_seed(5521)
    data <- make_markov_trajectories(
      n_patients = 80,
      follow_up_time = 4,
      n_states = 3,
      thresholds = c(-1, 1),
      allowed_start_state = 1:2,
      absorbing_state = 3,
      treatment_effect = 0.3,
      seed = 5521
    ) |>
      prepare_markov_data(absorbing_state = 3)
    mostm(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      first_followup_time = 1,
      id_var = "id",
      engine = "vglm"
    )
  })
}

test_that("plot_time stacks every state into bars spanning follow-up", {
  skip_on_cran()
  model <- make_plot_time_model()
  by_arm <- avg_state_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3
  )
  by_arm_days <- avg_state_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3,
    time_unit = "days"
  )
  cohort <- avg_state_time(model, times = 1:3, absorb = 3)
  by_sex_estimate <- c(1, 3, 2, 2, 3, 5, 4, 4)
  by_sex <- structure(
    data.frame(
      state_set = rep(c("1", "2"), 4),
      tx = rep(c(0, 0, 1, 1), 2),
      sex = rep(c("f", "m"), each = 4),
      estimate = by_sex_estimate
    ),
    class = c("markov_avg_state_time", "data.frame"),
    settings = markov_settings(y_levels = c("1", "2")),
    time_args = list(state_sets = list("1" = "1", "2" = "2"))
  )

  # Three visits, each spent in exactly one state, so time-scale bars end at 3
  # and proportion bars at 1. The first treatment level is drawn at the top.
  cases <- list(
    list(
      label = "two arms on the time scale",
      x = by_arm,
      args = list(),
      n_layers = 2L,
      n_bars = 2L,
      bar_end = 3,
      segments = by_arm$estimate,
      x_label = "Average time in state",
      y_label = "tx",
      groups = c("1", "0"),
      states = c("1", "2", "3")
    ),
    list(
      label = "two arms on the time scale with a time unit",
      x = by_arm_days,
      args = list(),
      n_layers = 2L,
      n_bars = 2L,
      bar_end = 3,
      segments = by_arm_days$estimate,
      x_label = "Average time in state (days)",
      y_label = "tx",
      groups = c("1", "0"),
      states = c("1", "2", "3")
    ),
    list(
      label = "two arms as proportions without labels",
      x = by_arm_days,
      args = list(scale = "proportion", show_labels = FALSE),
      n_layers = 1L,
      n_bars = 2L,
      bar_end = 1,
      segments = by_arm_days$estimate / 3,
      x_label = "Proportion of follow-up",
      y_label = "tx",
      groups = c("1", "0"),
      states = c("1", "2", "3")
    ),
    list(
      label = "one bar per arm within each facet",
      x = by_sex,
      args = list(facet_var = "sex", scale = "proportion"),
      n_layers = 2L,
      n_bars = 4L,
      bar_end = 1,
      segments = by_sex_estimate / rep(c(4, 8), each = 4),
      x_label = "Proportion of follow-up",
      y_label = "tx",
      groups = c("1", "0"),
      states = c("1", "2")
    ),
    list(
      label = "single observed-cohort bar",
      x = cohort,
      args = list(),
      n_layers = 2L,
      n_bars = 1L,
      bar_end = 3,
      segments = cohort$estimate,
      x_label = "Average time in state",
      y_label = NULL,
      groups = "All patients",
      states = c("1", "2", "3")
    )
  )

  for (case in cases) {
    label <- case$label
    plot <- do.call(plot_time, c(list(case$x), case$args))
    bars <- ggplot2::layer_data(plot, 1L)
    bar_ends <- tapply(bars$xmax, list(bars$PANEL, bars$y), max)
    bar_ends <- as.vector(bar_ends[!is.na(bar_ends)])

    expect_length(plot$layers, case$n_layers)
    expect_length(bar_ends, case$n_bars)
    expect_equal(bar_ends, rep(case$bar_end, case$n_bars), info = label)
    expect_equal(
      sort(bars$x - bars$xmin),
      sort(case$segments),
      info = label
    )
    expect_identical(plot$labels$x, case$x_label, info = label)
    expect_identical(plot$labels$y, case$y_label, info = label)
    expect_identical(levels(plot$data$.group), case$groups, info = label)
    expect_identical(levels(plot$data$state_set), case$states, info = label)
  }
})

test_that("plot_time pointrange draws confidence intervals per state", {
  skip_on_cran()
  model <- make_plot_time_model()
  total <- avg_state_time(model, variables = list(tx = c(0, 1)), times = 1:3, absorb = 3) |>
    inferences(method = "delta")

  plot <- plot_time(total, type = "pointrange")
  built <- ggplot2::ggplot_build(plot)
  intervals <- built$data[[1]]

  expect_s3_class(plot$layers[[1]]$geom, "GeomLinerange")
  expect_s3_class(plot$layers[[2]]$geom, "GeomPoint")
  expect_equal(sort(intervals$ymin), sort(total$conf.low))
  expect_equal(sort(intervals$ymax), sort(total$conf.high))
  expect_identical(plot$labels$x, "State")
  expect_identical(plot$labels$y, "Average time in state")
  expect_identical(plot$labels$colour, "tx")
  expect_null(plot$mapping$shape)

  hidden <- plot_time(total, type = "pointrange", show_uncertainty = FALSE)
  expect_length(hidden$layers, 1L)
})

test_that("plot_time pointrange draws scenarios from black to lighter grays", {
  total <- structure(
    data.frame(
      state_set = "1",
      arm = factor(paste("arm", 1:8), levels = paste("arm", 1:8)),
      estimate = 1:8
    ),
    class = c("markov_avg_state_time", "data.frame"),
    settings = markov_settings(y_levels = "1"),
    time_args = list(state_sets = list("1" = "1"))
  )

  points <- ggplot2::ggplot_build(plot_time(total, type = "pointrange"))$data[[1]]
  points <- points[order(points$group), ]
  rgb <- grDevices::col2rgb(points$colour)
  gray <- colSums(rgb * c(0.299, 0.587, 0.114))

  expect_length(unique(points$colour), 8L)
  expect_identical(points$colour[1], "#000000")
  expect_true(all(rgb[1, ] == rgb[2, ] & rgb[2, ] == rgb[3, ]))
  expect_true(all(diff(gray) > 0))
})

test_that("plot_time skips labels on narrow segments", {
  total <- structure(
    data.frame(
      state_set = c("1", "2", "1", "2"),
      tx = c(0, 0, 1, 1),
      estimate = c(0.1, 9.9, 5, 5)
    ),
    class = c("markov_avg_state_time", "data.frame"),
    settings = markov_settings(y_levels = c("1", "2")),
    time_args = list(state_sets = list("1" = "1", "2" = "2"))
  )

  labels <- ggplot2::ggplot_build(plot_time(total, digits = 0))$data[[2]]

  expect_setequal(labels$label, c("", "10", "5"))
  expect_identical(sum(labels$label == ""), 1L)
})

test_that("plot_time accepts or rejects its inputs", {
  skip_on_cran()
  model <- make_plot_time_model()
  total <- avg_state_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3
  )
  overlapping <- avg_state_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3,
    state_sets = list(alive = c("1", "2"), well = "1")
  )
  partial <- avg_state_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3,
    state_sets = list(well = "1", ill = "2")
  )

  # `error = NULL` marks a call that must draw a plot.
  cases <- list(
    list(
      label = "bars that would count a state twice",
      call = function() plot_time(overlapping),
      error = "do not overlap"
    ),
    list(
      label = "overlapping state sets as point ranges",
      call = function() plot_time(overlapping, type = "pointrange"),
      error = NULL
    ),
    list(
      label = "bars that leave out a state, on the time scale",
      call = function() plot_time(partial),
      error = NULL
    ),
    list(
      label = "proportions that leave out a state",
      call = function() plot_time(partial, scale = "proportion"),
      error = "missing state(s): 3"
    ),
    list(
      label = "not a state-time summary",
      call = function() plot_time(data.frame(state_set = "1", estimate = 1)),
      error = "`x` must be a `markov_avg_state_time` object"
    ),
    list(
      label = "unknown group_var",
      call = function() plot_time(total, group_var = "arm"),
      error = "`group_var` column not found in `data`: arm"
    ),
    list(
      label = "missing show_labels",
      call = function() plot_time(total, show_labels = NA),
      error = "`show_labels` must be TRUE or FALSE."
    )
  )

  for (case in cases) {
    if (is.null(case$error)) {
      expect_true(ggplot2::is_ggplot(case$call()), info = case$label)
    } else {
      expect_error(case$call(), case$error, fixed = TRUE, info = case$label)
    }
  }
})
