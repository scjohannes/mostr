# Checks of the prediction rows before they reach an engine
# (`prepare_prediction_rows()`): missing predictor values are refused, and the
# previous-state columns take the model's fitted states. test-known-bugs.R
# has the first tests of both, through sops() and avg_sops().

# First- and second-order orm and vglm fits of one data set, shared by the
# tests below. `bmi` is constant within each patient; `ypprev` is the state
# two visits back.
prediction_rows_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6101)
      data$bmi <- 20 + data$id %% 7
      data <- data[order(data$id, data$time), , drop = FALSE]
      data$ypprev <- stats::ave(
        as.integer(as.character(data$yprev)),
        data$id,
        FUN = function(x) c(x[1L], utils::head(x, -1L))
      )
      data$ypprev <- factor(data$ypprev, levels = levels(data$yprev))
      first <- ordered(y) ~ time + tx + yprev + bmi
      second <- ordered(y) ~ time + tx + yprev + ypprev + bmi
      cache <<- list(
        profiles = data[data$time == 1, , drop = FALSE],
        orm = mostm(first, data = data, id_var = "id", absorb = 6),
        vglm = mostm(first, data = data, id_var = "id", absorb = 6, engine = "vglm"),
        orm2 = mostm(second, data = data, id_var = "id", absorb = 6, p2_var = "ypprev"),
        vglm2 = mostm(
          second,
          data = data,
          id_var = "id",
          absorb = 6,
          p2_var = "ypprev",
          engine = "vglm"
        )
      )
    }
    cache
  }
})

test_that("missing values in columns the model does not use, or that `variables` sets, are allowed", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- prediction_rows_fits()
  profiles <- fits$profiles[1:4, ]
  unused <- profiles
  unused$note <- NA
  set_by_variables <- profiles
  set_by_variables$tx <- NA

  for (engine in c("orm", "vglm")) {
    expect_equal(
      sops(fits[[engine]], newdata = unused, times = 1:3)$estimate,
      sops(fits[[engine]], newdata = profiles, times = 1:3)$estimate,
      info = engine
    )
    expect_equal(
      avg_comparisons(
        fits[[engine]],
        newdata = set_by_variables,
        variables = "tx",
        times = 1:3
      )$estimate,
      avg_comparisons(
        fits[[engine]],
        newdata = profiles,
        variables = "tx",
        times = 1:3
      )$estimate,
      info = engine
    )
  }
})

test_that("soprob_markov() and plot_lp_difference() refuse a missing predictor, naming it", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- prediction_rows_fits()
  profiles <- fits$profiles[1:3, ]
  profiles$bmi[2] <- NA

  for (engine in c("orm", "vglm")) {
    expect_error(
      soprob_markov(fits[[engine]], newdata = profiles, times = 1:3),
      "`bmi` in 1 of the 3 prediction rows",
      class = "mostr_error",
      info = engine
    )
    expect_error(
      plot_lp_difference(
        fits[[engine]],
        profile = profiles[2, ],
        variables = list(tx = c(0, 1)),
        times = 1:3
      ),
      "`bmi`",
      class = "mostr_error",
      info = engine
    )
  }
})

test_that("a vglm design with a missing value is refused rather than shortened", {
  # stats::model.matrix() drops such a row by default, which shifted the
  # probabilities of the other rows onto the wrong profiles.
  skip_on_cran()
  skip_if_not_installed("VGAM")
  fit <- prediction_rows_fits()$vglm
  profiles <- prediction_rows_fits()$profiles[1:2, ]
  profiles$bmi[1] <- NA

  expect_error(
    vglm_prediction_design(fit$vglm_fit %||% fit, profiles),
    "missing values in 1 of the 2 prediction rows, in columns `bmi`",
    class = "mostr_error"
  )
})

test_that("second-order SOPs take the fitted states whatever the levels or type of both previous-state columns", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- prediction_rows_fits()
  profiles <- fits$profiles[1:4, ]
  own_levels <- profiles
  own_levels$yprev <- factor(as.character(profiles$yprev))
  own_levels$ypprev <- factor(as.character(profiles$ypprev))
  integers <- profiles
  integers$yprev <- as.integer(as.character(profiles$yprev))
  integers$ypprev <- as.integer(as.character(profiles$ypprev))
  not_a_state <- profiles
  not_a_state$ypprev <- as.character(not_a_state$ypprev)
  not_a_state$ypprev[1] <- "9"

  for (engine in c("orm2", "vglm2")) {
    expected <- sops(fits[[engine]], newdata = profiles, times = 1:3)$estimate
    expect_equal(
      sops(fits[[engine]], newdata = own_levels, times = 1:3)$estimate,
      expected,
      info = engine
    )
    expect_equal(
      sops(fits[[engine]], newdata = integers, times = 1:3)$estimate,
      expected,
      info = engine
    )
    expect_error(
      sops(fits[[engine]], newdata = not_a_state, times = 1:3),
      "`ypprev` holds 9",
      class = "mostr_error",
      info = engine
    )
  }
})
