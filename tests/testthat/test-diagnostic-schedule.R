# The model-based diagnostic plots run the SOP recursion one transition per
# visit, filling in fitted visits that `times` leaves out. They used to fill in
# every whole number from 1, whatever the fit's visits were, so a fit with
# `gap_var`, weekly visits, landmark `newdata` or character time was plotted
# from a recursion the model was not fitted for (issue #51).

# The times each prediction of `expr` was made at, in order.
diagnostic_predicted_times <- function(expr) {
  seen <- list()
  predict <- predict_orm_response_markov
  with_mocked_bindings(
    predict_orm_response_markov = function(object, newdata) {
      seen[[length(seen) + 1L]] <<- unique(as.character(newdata$time))
      predict(object, newdata)
    },
    expr
  )
  unlist(seen)
}

# The recursion schedule of a diagnostic plot and the SOPs it produces at the
# plotted times, averaged over patients.
diagnostic_plot_sops <- function(fit, times, newdata = NULL) {
  setup <- plot_transition_model_setup(
    fit,
    newdata = newdata,
    variables = NULL,
    times = times,
    settings = plot_markov_settings(fit)
  )
  trace <- markov_transition_trace(
    fit,
    setup$data,
    setup$times,
    setup$settings
  )
  average <- apply(trace$sops, c(2, 3), mean)
  list(
    schedule = as.character(setup$times),
    sops = average[as.character(setup$plot_times), , drop = FALSE]
  )
}

# `avg_sops()` as a time-by-state matrix.
diagnostic_avg_sops <- function(fit, times, newdata = NULL) {
  result <- as.data.frame(avg_sops(fit, newdata = newdata, times = times))
  states <- unique(as.character(result$state))
  out <- t(vapply(
    as.character(times),
    function(time) {
      rows <- as.character(result$time) == time
      result$estimate[rows][match(states, as.character(result$state[rows]))]
    },
    numeric(length(states))
  ))
  dimnames(out) <- list(as.character(times), states)
  out
}

# Compares the plotted SOPs with `avg_sops()` over every visit of `schedule`,
# the visits the model is fitted for, at the plotted times.
expect_plot_sops_match <- function(plotted, fit, schedule, newdata = NULL) {
  reference <- diagnostic_avg_sops(fit, schedule, newdata)
  shown <- rownames(plotted$sops)
  expect_equal(
    unname(plotted$sops[, colnames(reference), drop = FALSE]),
    unname(reference[shown, , drop = FALSE]),
    tolerance = 1e-10
  )
}

test_that("diagnostic plots of a `gap_var` fit run the requested visits (#51)", {
  skip_if_not_installed("rms")
  fit <- grid_fit(grid_cell("orm", visit = "gap"))
  # Fitted visits 1, 2, 4, 5, 7, 8; the gaps of `times` are 1, 1, 2, 1.
  times <- c(1, 2, 4, 5)

  plotted <- diagnostic_plot_sops(fit, times)
  expect_equal(plotted$schedule, as.character(times))
  expect_plot_sops_match(plotted, fit, times)
})

test_that("all three diagnostic plots run the schedule of the set-up (#51)", {
  # Slow, and repeats the schedule check above through the public plots.
  skip_on_cran()
  skip_if_not_installed("rms")
  fit <- grid_fit(grid_cell("orm", visit = "gap"))
  times <- c(1, 2, 4, 5)
  for (plot in list(plot_transitions, plot_correlation, plot_variogram)) {
    expect_equal(
      diagnostic_predicted_times(plot(fit, times = times)),
      as.character(times)
    )
  }
})

test_that("diagnostic plots of a weekly fit run the weekly visits (#51)", {
  skip_if_not_installed("rms")
  data <- grid_data()
  data$time <- data$time * 7
  fit <- mostm(
    y ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorb = 4,
    first_followup_time = 7
  )

  plotted <- diagnostic_plot_sops(fit, c(7, 14, 21, 28))
  expect_equal(plotted$schedule, as.character(c(7, 14, 21, 28)))
  expect_plot_sops_match(plotted, fit, c(7, 14, 21, 28))

  # Sparse weekly times fill in the weekly visits between them.
  sparse <- diagnostic_plot_sops(fit, c(7, 28))
  expect_equal(sparse$schedule, as.character(c(7, 14, 21, 28)))
  expect_plot_sops_match(sparse, fit, c(7, 14, 21, 28))

  # A time between weekly visits would add a transition the model does not
  # have, and the stored profiles hold the state before the first follow-up
  # visit, so neither can be plotted.
  expect_error(
    plot_transitions(fit, times = c(7, 10, 14)),
    "fitted visit",
    class = "mostr_error"
  )
  expect_error(
    plot_transitions(fit, times = c(1, 7)),
    "first follow-up visit",
    class = "mostr_error"
  )
})

test_that("diagnostic plots of landmark `newdata` start at its first visit (#51)", {
  skip_if_not_installed("rms")
  fit <- grid_fit(grid_cell("orm"))
  landmark <- grid_data()
  landmark <- landmark[landmark$time == 3, , drop = FALSE]

  plotted <- diagnostic_plot_sops(fit, 3:6, newdata = landmark)
  expect_equal(plotted$schedule, as.character(3:6))
  expect_plot_sops_match(plotted, fit, 3:6, newdata = landmark)

  sparse <- diagnostic_plot_sops(fit, c(3, 6), newdata = landmark)
  expect_equal(sparse$schedule, as.character(3:6))
  expect_plot_sops_match(sparse, fit, 3:6, newdata = landmark)

  # Daily stored profiles still fill in the visits from the first follow-up.
  daily <- diagnostic_plot_sops(fit, c(2, 4))
  expect_equal(daily$schedule, as.character(1:4))
  expect_plot_sops_match(daily, fit, 1:4)

  # The same for factor time: the recursion starts at the first requested
  # level of supplied rows, not at the first level.
  factor_fit <- grid_fit(grid_cell("orm", time = "factor"))
  factor_landmark <- grid_cell_data(grid_cell("orm", time = "factor"))
  factor_landmark <- factor_landmark[factor_landmark$time == 3, , drop = FALSE]
  factor_landmark$time <- factor(
    factor_landmark$time,
    levels = levels(factor(grid_data()$time))
  )
  factor_plotted <- diagnostic_plot_sops(
    factor_fit,
    c("3", "5"),
    newdata = factor_landmark
  )
  expect_equal(factor_plotted$schedule, c("3", "4", "5"))
})

test_that("diagnostic plots of character time keep the requested order (#51)", {
  skip_if_not_installed("rms")
  labels <- c("one", "two", "three", "four", "five", "six")
  data <- grid_data()
  data$time <- labels[data$time]
  fit <- mostm(
    y ~ tx + yprev,
    data = data,
    id_var = "id",
    absorb = 4,
    first_followup_time = "one"
  )

  # Alphabetical order would put "six" fourth.
  plotted <- diagnostic_plot_sops(fit, labels)
  expect_equal(plotted$schedule, labels)
  expect_plot_sops_match(plotted, fit, labels)
})

# A character time column that is in the formula is still character time: the
# recursion runs the visits in the requested order, as for character time
# outside the formula. The fit's levels of such a column are in alphabetical
# order, so treating it as factor time ran "six" fourth, refused "one" as
# coming after the first follow-up visit "five", and warned that `times` was
# not in increasing order (the #44 change to get_time_info() undid #51 here).
test_that("character time in the formula keeps the requested order (#44, #51)", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  labels <- c("one", "two", "three", "four", "five", "six")
  data <- grid_data()
  data$time <- labels[data$time]
  factor_data <- data
  factor_data$time <- factor(factor_data$time, levels = labels)

  for (engine in c("orm", "vglm")) {
    fit_args <- list(
      formula = y ~ tx + yprev + time,
      id_var = "id",
      absorb = 4,
      first_followup_time = "one",
      engine = engine
    )
    fit <- do.call(mostm, c(list(data = data), fit_args))
    factor_fit <- do.call(mostm, c(list(data = factor_data), fit_args))

    # The same model as with factor time in visit order, so the same SOPs,
    # and no warning about the order of `times`.
    expect_no_warning(sops_chr <- diagnostic_avg_sops(fit, labels))
    expect_equal(
      sops_chr,
      diagnostic_avg_sops(factor_fit, labels),
      tolerance = 1e-6,
      label = paste(engine, "character-time SOPs")
    )

    plotted <- diagnostic_plot_sops(fit, labels)
    expect_equal(plotted$schedule, labels, label = paste(engine, "schedule"))
    expect_plot_sops_match(plotted, fit, labels)
  }
})
