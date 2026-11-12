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
      p2_var = "ypprev"
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
      time_var = "day"
    ))
  })
}

settings_time_fit <- function() {
  grid_fixture("settings orm with time column time", function() {
    suppressWarnings(mostm(
      y ~ time + tx + age + yprev,
      data = grid_data(),
      id_var = "id"
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
    args <- c(list(times = 1:4, absorb = 4), case$args)
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
      sops(day, newdata = day_profiles, times = 1:4, absorb = 4)
    )$estimate,
    suppressWarnings(
      sops(time, newdata = profiles, times = 1:4, absorb = 4)
    )$estimate
  )
  expect_equal(
    suppressWarnings(avg_sops(day, times = 1:4, absorb = 4))$estimate,
    suppressWarnings(avg_sops(time, times = 1:4, absorb = 4))$estimate
  )
})

test_that("the prediction functions refuse the settings that the fit records", {
  first <- grid_fit(grid_cell("orm"))
  profiles <- grid_profiles(grid_cell("orm"))
  expect_snapshot(error = TRUE, {
    sops(first, times = 1:4, absorb = 4, p2_var = "ypprev")
    avg_sops(first, times = 1:4, absorb = 4, p_var = "yprev")
    state_time(first, times = 1:4, absorb = 4, time_var = "time")
    comparisons(first, variables = "tx", times = 1:4, absorb = 4, id_var = "id")
    avg_comparisons(first, variables = "tx", times = 1:4, absorb = 4, y_levels = 1:4)
  })
  # Each function refuses each argument, also when the value repeats the fit.
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
    y_levels = as.character(1:4)
  )
  for (name in names(calls)) {
    for (arg in names(fitted)) {
      expect_error(
        do.call(
          calls[[name]]$fun,
          c(
            list(first, times = 1:4, absorb = 4),
            calls[[name]]$args,
            fitted[arg]
          )
        ),
        paste0("`", arg, "` is not an argument of this function"),
        fixed = TRUE,
        info = paste(name, arg)
      )
    }
  }
  expect_error(
    plot_lp_difference(
      first,
      profile = profiles[1, , drop = FALSE],
      variables = list(tx = c(0, 1)),
      times = 1:4,
      absorb = 4,
      p_var = "yprev"
    ),
    "unused argument"
  )
})

test_that("diagnostic plots refuse a setting that contradicts the fit", {
  second <- settings_second_order_fit()
  first <- grid_fit(grid_cell("orm"))
  day <- settings_day_fit()
  expect_snapshot(error = TRUE, {
    plot_transitions(second, times = 1:4, absorb = 4, p2_var = "yprev")
    plot_transitions(first, times = 1:4, absorb = 4, p2_var = "ypprev")
    plot_correlation(first, times = 1:4, absorb = 4, p_var = "ypprev")
    plot_variogram(day, times = 1:4, absorb = 4, time_var = "time")
    plot_transitions(first, times = 1:4, absorb = 4, id_var = "age")
    plot_correlation(first, times = 1:4, absorb = 4, y_levels = c(1, 2, 4, 3))
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
  expect_identical(settings$engine, "orm")
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
})
