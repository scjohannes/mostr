# The column names and model settings that `mostm()` records on a fit. The SOP
# functions take them from the fit, so leaving one out can't silently predict
# a different model, and a value that contradicts the fit is an error.

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

test_that("a second-order fit predicts second order without `p2_var`", {
  fit <- settings_second_order_fit()
  profiles <- grid_profiles(grid_cell("orm", order = 2L))
  # Each call without `p2_var` must give the second-order results that naming
  # the column gives, which differ from holding `ypprev` fixed.
  cases <- list(
    list(fun = sops, args = list(newdata = profiles)),
    list(fun = avg_sops, args = list(variables = list(tx = c(0, 1)))),
    list(
      fun = comparisons,
      args = list(newdata = profiles, variables = list(tx = c(0, 1)))
    ),
    list(fun = avg_state_time, args = list(variables = list(tx = c(0, 1)))),
    list(
      fun = soprob_markov,
      args = list(newdata = profiles, y_levels = 1:4)
    )
  )
  for (case in cases) {
    common <- c(list(fit, times = 1:4, absorb = 4), case$args)
    implicit <- suppressWarnings(do.call(case$fun, common))
    explicit <- suppressWarnings(
      do.call(case$fun, c(common, list(p2_var = "ypprev")))
    )
    if (is.data.frame(implicit)) {
      expect_identical(implicit$estimate, explicit$estimate)
    } else {
      expect_identical(implicit, explicit)
    }
  }
  first_order <- suppressWarnings(
    avg_sops(fit, times = 1:4, absorb = 4)
  )
  expect_identical(
    first_order$estimate,
    suppressWarnings(avg_sops(fit, times = 1:4, absorb = 4, p2_var = "ypprev"))$estimate
  )
})

test_that("predictions use the fit's time column without `time_var`", {
  fit <- settings_day_fit()
  profiles <- grid_data()[grid_data()$time == 1, , drop = FALSE][1:5, ]
  names(profiles)[names(profiles) == "time"] <- "day"
  implicit <- suppressWarnings(
    sops(fit, newdata = profiles, times = 1:4, absorb = 4)
  )
  explicit <- suppressWarnings(
    sops(fit, newdata = profiles, times = 1:4, absorb = 4, time_var = "day")
  )
  expect_identical(implicit$estimate, explicit$estimate)
  stored <- suppressWarnings(avg_sops(fit, times = 1:4, absorb = 4))
  expect_identical(
    stored$estimate,
    suppressWarnings(avg_sops(fit, times = 1:4, absorb = 4, time_var = "day"))$estimate
  )
})

