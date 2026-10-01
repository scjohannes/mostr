make_plot_time_model <- function() {
  skip_if_not_installed("VGAM")
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
  vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id"
  )
}

test_that("plot_time stacks every state into bars spanning follow-up", {
  model <- make_plot_time_model()
  total <- avg_time(model, variables = list(tx = c(0, 1)), times = 1:3, absorb = 3)

  plot <- plot_time(total)
  built <- ggplot2::ggplot_build(plot)
  bars <- built$data[[1]]

  expect_s3_class(plot, "ggplot")
  expect_length(built$data, 2L)
  # Three visits, each spent in exactly one state.
  expect_equal(as.vector(tapply(bars$xmax, bars$y, max)), c(3, 3))
  expect_equal(sort(bars$x - bars$xmin), sort(total$estimate))
  expect_identical(plot$labels$y, "tx")
  expect_identical(plot$labels$x, "Average time in state")
  expect_identical(
    levels(plot$data$state_set),
    c("1", "2", "3")
  )
  # The first treatment level is drawn at the top.
  expect_identical(levels(plot$data$.group), c("1", "0"))
})

test_that("plot_time proportions divide each bar by its follow-up", {
  model <- make_plot_time_model()
  total <- avg_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3,
    time_unit = "days"
  )

  plot <- plot_time(total, scale = "proportion", show_labels = FALSE)
  bars <- ggplot2::ggplot_build(plot)$data[[1]]

  expect_length(plot$layers, 1L)
  expect_equal(as.vector(tapply(bars$xmax, bars$y, max)), c(1, 1))
  expect_identical(plot$labels$x, "Proportion of follow-up")
  expect_identical(
    plot_time(total)$labels$x,
    "Average time in state (days)"
  )
})

test_that("plot_time pointrange draws confidence intervals per state", {
  model <- make_plot_time_model()
  total <- avg_time(model, variables = list(tx = c(0, 1)), times = 1:3, absorb = 3) |>
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
    class = c("markov_avg_time", "data.frame"),
    y_levels = "1",
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
    class = c("markov_avg_time", "data.frame"),
    y_levels = c("1", "2"),
    time_args = list(state_sets = list("1" = "1", "2" = "2"))
  )

  labels <- ggplot2::ggplot_build(plot_time(total, digits = 0))$data[[2]]

  expect_setequal(labels$label, c("", "10", "5"))
  expect_identical(sum(labels$label == ""), 1L)
})

test_that("plot_time draws one bar per scenario within each facet", {
  total <- structure(
    data.frame(
      state_set = rep(c("1", "2"), 4),
      tx = rep(c(0, 0, 1, 1), 2),
      sex = rep(c("f", "m"), each = 4),
      estimate = c(1, 3, 2, 2, 3, 5, 4, 4)
    ),
    class = c("markov_avg_time", "data.frame"),
    y_levels = c("1", "2"),
    time_args = list(state_sets = list("1" = "1", "2" = "2"))
  )

  plot <- plot_time(total, facet_var = "sex", scale = "proportion")
  bars <- ggplot2::ggplot_build(plot)$data[[1]]

  expect_identical(plot$labels$y, "tx")
  expect_equal(as.vector(tapply(bars$xmax, list(bars$PANEL, bars$y), max)), rep(1, 4))
  expect_equal(sort(bars$x - bars$xmin), sort(c(1, 3, 2, 2, 3, 5, 4, 4) / rep(c(4, 8), each = 4)))
})

test_that("plot_time labels a single observed-cohort bar", {
  model <- make_plot_time_model()
  total <- avg_time(model, times = 1:3, absorb = 3)

  plot <- plot_time(total)

  expect_identical(levels(plot$data$.group), "All patients")
  expect_null(plot$labels$y)
})

test_that("plot_time rejects bars that would count a state twice", {
  model <- make_plot_time_model()
  overlapping <- avg_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3,
    state_sets = list(alive = c("1", "2"), well = "1")
  )
  partial <- avg_time(
    model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = 3,
    state_sets = list(well = "1", ill = "2")
  )

  expect_error(plot_time(overlapping), "do not overlap", fixed = TRUE)
  expect_s3_class(plot_time(overlapping, type = "pointrange"), "ggplot")
  expect_s3_class(plot_time(partial), "ggplot")
  expect_error(
    plot_time(partial, scale = "proportion"),
    "missing state(s): 3",
    fixed = TRUE
  )
})

test_that("plot_time validates its inputs", {
  expect_error(
    plot_time(data.frame(state_set = "1", estimate = 1)),
    "`x` must be a `markov_avg_time` object",
    fixed = TRUE
  )
  model <- make_plot_time_model()
  total <- avg_time(model, variables = list(tx = c(0, 1)), times = 1:3, absorb = 3)
  expect_error(
    plot_time(total, group_var = "arm"),
    "`group_var` column not found in `data`: arm",
    fixed = TRUE
  )
  expect_error(
    plot_time(total, show_labels = NA),
    "`show_labels` must be TRUE or FALSE.",
    fixed = TRUE
  )
})
