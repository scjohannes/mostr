# The column names and model settings that `mostm()` records on a fit. The SOP
# functions take them from the fit and have no arguments for them, so a call
# can't predict a different model than the one fitted. The diagnostic plots,
# which also take observed data, accept them but refuse a value that
# contradicts the fit.

settings_second_order_fit <- function() {
  grid_fixture("settings second-order orm", function() {
    suppressWarnings(mostm(
      y ~ time + tx + age + yprev + ypprev,
      data = grid_data(),
      id_var = "id",
      p2_var = "ypprev",
      absorb = 4
    ))
  })
}

settings_day_fit <- function() {
  grid_fixture("settings orm with time column day", function() {
    data <- grid_data()
    names(data)[names(data) == "time"] <- "day"
    suppressWarnings(mostm(
      y ~ day + tx + age + yprev,
      data = data,
      id_var = "id",
      time_var = "day",
      absorb = 4
    ))
  })
}

settings_time_fit <- function() {
  grid_fixture("settings orm with time column time", function() {
    suppressWarnings(mostm(
      y ~ time + tx + age + yprev,
      data = grid_data(),
      id_var = "id",
      absorb = 4
    ))
  })
}

# The same fit with its settings changed to first order: `ypprev` is then an
# ordinary covariate, held at each profile's starting value.
settings_as_first_order <- function(fit) {
  markov_set_settings(
    fit,
    markov_settings_update(markov_model_settings(fit), p2_var = NULL)
  )
}

test_that("a second-order fit predicts with the second-order recursion", {
  fit <- settings_second_order_fit()
  held <- settings_as_first_order(fit)
  profiles <- grid_profiles(grid_cell("orm", order = 2L))
  # Every prediction function takes `p2_var` from the fit, so each gives
  # results that differ from holding `ypprev` fixed. The grid tests check that
  # the second-order results themselves are right.
  cases <- list(
    list(fun = sops, args = list(newdata = profiles)),
    list(fun = avg_sops, args = list(variables = list(tx = c(0, 1)))),
    list(
      fun = comparisons,
      args = list(newdata = profiles, variables = list(tx = c(0, 1)))
    ),
    list(
      fun = avg_comparisons,
      args = list(variables = list(tx = c(0, 1)))
    ),
    list(fun = state_time, args = list(newdata = profiles)),
    list(fun = avg_state_time, args = list(variables = list(tx = c(0, 1)))),
    list(fun = soprob_markov, args = list(newdata = profiles))
  )
  for (case in cases) {
    args <- c(list(times = 1:4), case$args)
    second <- suppressWarnings(do.call(case$fun, c(list(fit), args)))
    first <- suppressWarnings(do.call(case$fun, c(list(held), args)))
    if (is.data.frame(second)) {
      second <- second$estimate
      first <- first$estimate
    }
    expect_false(isTRUE(all.equal(second, first)))
  }
})

test_that("predictions use the fit's time column", {
  day <- settings_day_fit()
  time <- settings_time_fit()
  profiles <- grid_data()[grid_data()$time == 1, , drop = FALSE][1:5, ]
  day_profiles <- profiles
  names(day_profiles)[names(day_profiles) == "time"] <- "day"
  # The two fits differ only in the name of the time column.
  expect_equal(
    suppressWarnings(
      sops(day, newdata = day_profiles, times = 1:4)
    )$estimate,
    suppressWarnings(
      sops(time, newdata = profiles, times = 1:4)
    )$estimate
  )
  expect_equal(
    suppressWarnings(avg_sops(day, times = 1:4))$estimate,
    suppressWarnings(avg_sops(time, times = 1:4))$estimate
  )
})

test_that("the prediction functions refuse the settings that the fit records", {
  first <- grid_fit(grid_cell("orm"))
  profiles <- grid_profiles(grid_cell("orm"))
  expect_snapshot(error = TRUE, {
    sops(first, times = 1:4, p2_var = "ypprev")
    avg_sops(first, times = 1:4, p_var = "yprev")
    state_time(first, times = 1:4, time_var = "time")
    comparisons(first, variables = "tx", times = 1:4, id_var = "id")
    avg_comparisons(first, variables = "tx", times = 1:4, y_levels = 1:4)
    sops(first, times = 1:4, absorb = 4)
    avg_sops(first, times = 1:4, gap_var = "gap")
    avg_state_time(first, times = 1:4, time_covariates = data.frame(time = 1:4))
  })
  # Each function refuses each argument, also when the value repeats the fit.
  # `absorb` passed to a prediction function would otherwise be swallowed by
  # `...` and let patients leave a state the fit does not record as absorbing.
  calls <- list(
    sops = list(fun = sops, args = list()),
    avg_sops = list(fun = avg_sops, args = list()),
    state_time = list(fun = state_time, args = list()),
    avg_state_time = list(fun = avg_state_time, args = list()),
    comparisons = list(fun = comparisons, args = list(variables = "tx")),
    avg_comparisons = list(fun = avg_comparisons, args = list(variables = "tx")),
    soprob_markov = list(fun = soprob_markov, args = list(newdata = profiles))
  )
  fitted <- list(
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    p2_var = "ypprev",
    y_levels = as.character(1:4),
    absorb = 4,
    gap_var = "gap"
  )
  removed <- list(time_covariates = data.frame(time = 1:4))
  for (name in names(calls)) {
    for (arg in c(names(fitted), names(removed))) {
      message <- if (arg %in% names(removed)) {
        paste0("`", arg, "` was removed")
      } else {
        paste0("`", arg, "` is not an argument of this function")
      }
      expect_error(
        do.call(
          calls[[name]]$fun,
          c(
            list(first, times = 1:4),
            calls[[name]]$args,
            c(fitted, removed)[arg]
          )
        ),
        message,
        fixed = TRUE,
        class = "mostr_error",
        info = paste(name, arg)
      )
    }
  }
  # The plots have no `...`, so R itself refuses the arguments they lost.
  for (arg in c("absorb", "gap_var", "time_covariates")) {
    value <- c(fitted, removed)[arg]
    expect_error(
      do.call(
        plot_lp_difference,
        c(
          list(
            first,
            profile = profiles[1, , drop = FALSE],
            variables = list(tx = c(0, 1)),
            times = 1:4
          ),
          value
        )
      ),
      "unused argument",
      info = paste("plot_lp_difference", arg)
    )
    for (plot in list(plot_transitions, plot_correlation, plot_variogram)) {
      expect_error(
        do.call(plot, c(list(first, times = 1:4), value)),
        "unused argument",
        info = arg
      )
    }
  }
})

test_that("diagnostic plots refuse a setting that contradicts the fit", {
  second <- settings_second_order_fit()
  first <- grid_fit(grid_cell("orm"))
  day <- settings_day_fit()
  expect_snapshot(error = TRUE, {
    plot_transitions(second, times = 1:4, p2_var = "yprev")
    plot_transitions(first, times = 1:4, p2_var = "ypprev")
    plot_correlation(first, times = 1:4, p_var = "ypprev")
    plot_variogram(day, times = 1:4, time_var = "time")
    plot_transitions(first, times = 1:4, id_var = "age")
    plot_correlation(first, times = 1:4, y_levels = c(1, 2, 4, 3))
  })
})

test_that("`mostm()` records the settings and checks the previous-state columns", {
  fit <- settings_second_order_fit()
  settings <- markov_model_settings(fit)
  expect_s3_class(settings, "markov_settings")
  expect_identical(settings$id_var, "id")
  expect_identical(settings$time_var, "time")
  expect_identical(settings$p_var, "yprev")
  expect_identical(settings$p2_var, "ypprev")
  expect_identical(settings$order, 2L)
  expect_identical(as_state_labels(settings$y_levels), as.character(1:4))

  vglm_fit <- grid_fit(grid_cell("vglm", order = 2L))
  expect_identical(markov_model_settings(vglm_fit)$p2_var, "ypprev")
  expect_identical(markov_model_settings(vglm_fit$vglm_fit)$p2_var, "ypprev")

  data <- grid_data()
  expect_snapshot(error = TRUE, {
    mostm(y ~ time + yprev, data = data, id_var = "id", p2_var = "missing")
    mostm(y ~ time + yprev, data = data, id_var = "id", p_var = "yprev", p2_var = "yprev")
    mostm(y ~ time + yprev, data = data, id_var = "id", p_var = c("a", "b"))
  })
})

test_that("bootstrap refits keep the fit's settings", {
  fit <- settings_second_order_fit()
  refit <- update_bootstrap_model(fit, markov_model_refit_data(fit))
  expect_identical(markov_model_settings(refit)$p2_var, "ypprev")
  expect_identical(markov_model_settings(refit)$time_var, "time")
  expect_identical(markov_model_settings(refit)$absorb, 4)
})

test_that("`mostm()` records `absorb` and `gap_var`, and the predictions use them", {
  fit <- settings_time_fit()
  settings <- markov_model_settings(fit)
  expect_identical(settings$absorb, 4)
  expect_null(settings$gap_var)
  expect_contains(settings$from_fit, c("absorb", "gap_var"))
  gap_fit <- grid_fit(grid_cell("orm", visit = "gap"))
  expect_identical(markov_model_settings(gap_fit)$gap_var, "gap")

  result <- sops(fit, newdata = grid_profiles(grid_cell("orm")), times = 1:4)
  expect_identical(markov_result_settings(result)$absorb, 4)
  # State 4 keeps its probability: a patient who reaches it stays there.
  expect_lt(grid_largest_fall(result, grid_n_states, "rowid", 1:4), 1e-12)
})

test_that("`mostm()` refuses `absorb` labels that are not states of the response", {
  data <- grid_data()
  expect_snapshot(error = TRUE, {
    mostm(y ~ time + yprev, data = data, id_var = "id", absorb = 5)
    mostm(y ~ time + yprev, data = data, id_var = "id", absorb = c("4", "dead"))
    mostm(y ~ time + yprev, data = data, id_var = "id", absorb = c(4, NA))
  })
  skip_on_cran() # Slow: four fits.
  # The states are compared by label, so a number, a string and a factor
  # name the same state.
  for (absorb in list(4, 4L, "4", factor("4"))) {
    expect_no_error(mostm(
      y ~ time + yprev,
      data = data,
      id_var = "id",
      absorb = absorb
    ))
  }
})

test_that("`mostm()` refuses fitted rows that start in an absorbing state", {
  data <- grid_data()
  # State 2 is the previous state of many fitted rows, so it can't be
  # absorbing in these data.
  expect_snapshot(error = TRUE, {
    mostm(y ~ time + yprev, data = data, id_var = "id", absorb = c(2, 4))
  })
  # Without the rows that start in state 2, the state two visits back is
  # still 2 in some rows: a history that an absorbing state 2 makes
  # impossible.
  no_start_in_2 <- data[data$yprev != "2", , drop = FALSE]
  no_start_in_2$yprev <- droplevels(no_start_in_2$yprev)
  expect_snapshot(error = TRUE, {
    mostm(
      y ~ time + yprev + ypprev,
      data = no_start_in_2,
      id_var = "id",
      p2_var = "ypprev",
      absorb = c(2, 4)
    )
  })
  # The refusals above stop before fitting; the fits below are slower.
  skip_on_cran()
  # A first-order fit does not look at the older state.
  expect_no_error(mostm(
    y ~ time + yprev,
    data = no_start_in_2,
    id_var = "id",
    absorb = c(2, 4)
  ))
  # Rows left out by `subset` are not fitted, so they are not checked. A
  # numeric previous state, because rms can't fit an empty factor level.
  numeric <- data
  numeric$yprev <- as.numeric(as.character(numeric$yprev))
  expect_no_error(mostm(
    y ~ time + yprev,
    data = numeric,
    id_var = "id",
    absorb = c(2, 4),
    subset = yprev != 2
  ))
})

test_that("`mostm()` warns about a state seen only as an outcome that `absorb` leaves out", {
  data <- grid_data()
  expect_snapshot(fit <- mostm(y ~ time + yprev, data = data, id_var = "id"))
  # The warning names every such state and suggests the full `absorb`.
  no_start_in_2 <- data[data$yprev != "2", , drop = FALSE]
  no_start_in_2$yprev <- droplevels(no_start_in_2$yprev)
  expect_warning(
    mostm(y ~ time + yprev, data = no_start_in_2, id_var = "id"),
    "absorb = c(\"2\", \"4\")",
    fixed = TRUE,
    class = "mostr_warning_outcome_only_states"
  )
  skip_on_cran() # Slow: three more fits.
  expect_no_warning(
    mostm(y ~ time + yprev, data = data, id_var = "id", absorb = 4),
    class = "mostr_warning_outcome_only_states"
  )
  expect_warning(
    mostm(y ~ time + yprev, data = no_start_in_2, id_var = "id", absorb = 4),
    "absorb = c(\"4\", \"2\")",
    fixed = TRUE,
    class = "mostr_warning_outcome_only_states"
  )
  # Bootstrap refits don't warn again.
  fit <- settings_time_fit()
  expect_no_warning(
    update_bootstrap_model(fit, markov_model_refit_data(fit)),
    class = "mostr_warning_outcome_only_states"
  )
})

test_that("a numeric previous state without `absorb` warns at fit time (regression)", {
  # With a numeric previous state, a model without `absorb` predicts
  # transitions out of death from the other states' transitions, so dead
  # patients come back to life. Fitting such a model gave no message before.
  data <- grid_data()
  data$yprev <- as.numeric(as.character(data$yprev))
  expect_warning(
    without <- mostm(y ~ time + tx + yprev, data = data, id_var = "id"),
    class = "mostr_warning_outcome_only_states"
  )
  # Slow, and the size of the leak below depends on the seed.
  skip_on_cran()
  with <- mostm(y ~ time + tx + yprev, data = data, id_var = "id", absorb = 4)
  profiles <- grid_profiles(grid_cell("orm", previous = "numeric"))
  absorbed <- sops(with, newdata = profiles, times = 1:4)
  expect_lt(grid_largest_fall(absorbed, grid_n_states, "rowid", 1:4), 1e-12)
  leaking <- sops(without, newdata = profiles, times = 1:4)
  dead <- function(x) x$estimate[x$state == "4" & x$time == 4]
  expect_gt(max(abs(dead(leaking) - dead(absorbed))), 1e-3)
})

test_that("`mostm()` checks that `gap_var` is the time since the previous row", {
  cell <- grid_cell("orm", visit = "gap")
  data <- grid_cell_data(cell)
  formula <- y ~ time + gap + yprev
  fit_gap <- function(data, ...) {
    mostm(formula, data = data, id_var = "id", absorb = 4, gap_var = "gap", ...)
  }
  days <- data
  days$gap <- days$gap * 7
  text <- data
  text$gap <- as.character(text$gap)
  factor_time <- data
  factor_time$time <- factor(factor_time$time)
  expect_snapshot(error = TRUE, {
    fit_gap(days)
    fit_gap(text)
    mostm(formula, data = data, id_var = "id", absorb = 4, gap_var = "interval")
    mostm(formula, data = data, id_var = "id", absorb = 4, gap_var = c("gap", "time"))
    mostm(
      formula,
      data = factor_time,
      id_var = "id",
      absorb = 4,
      gap_var = "gap",
      first_followup_time = "1"
    )
  })

  # The refusals above stop before fitting; the fits below are slower.
  skip_on_cran()
  expect_no_error(fit_gap(data))
  # Gaps computed in another way may differ from the time differences by
  # rounding only.
  rounded <- data
  rounded$gap <- rounded$gap * (1 + 1e-12)
  expect_no_error(fit_gap(rounded))
  # A patient's first row is not checked: its gap is the time since time 0.
  first_rows <- data
  first_rows$gap[!duplicated(first_rows$id)] <- 3
  expect_no_error(fit_gap(first_rows))
  # The rows' order does not matter.
  expect_no_error(fit_gap(data[rev(seq_len(nrow(data))), , drop = FALSE]))
})

# The default `p_var = "yprev"` is recorded on the fit like a supplied one, so
# it must name a column of `data` too. Otherwise the fit succeeds, `absorb` is
# never checked against the real previous states, and the mistake shows only
# at the first prediction.
test_that("`mostm()` refuses a default `p_var` that is not a column of `data`", {
  skip_if_not_installed("rms")
  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 9461)
  names(data)[names(data) == "yprev"] <- "prev_state"
  expect_error(
    mostm(y ~ time + tx + prev_state, data = data, id_var = "id", absorb = 6),
    "p_var",
    class = "mostr_error"
  )
  expect_no_error(
    mostm(
      y ~ time + tx + prev_state,
      data = data,
      id_var = "id",
      p_var = "prev_state",
      absorb = 6
    )
  )
})
