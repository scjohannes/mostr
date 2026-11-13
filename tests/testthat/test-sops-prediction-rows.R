# Checks of the prediction rows before they reach an engine
# (`prepare_prediction_rows()`): missing predictor values are refused, and the
# previous-state columns take the model's fitted states. The tests at the end
# of this file, written before the fix, check both through sops() and
# avg_sops().

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

test_that("plot_lp_difference() plots a profile whose factor yprev has only its own level", {
  # vglm used to drop the rows it could not predict, so this plot stopped with
  # "The model returned linear predictors for 8 of the 24 plotted rows".
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- prediction_rows_fits()
  full <- fits$profiles[1, , drop = FALSE]
  own_level <- full
  own_level$yprev <- factor(as.character(own_level$yprev))

  for (engine in c("orm", "vglm")) {
    plot_data <- function(profile) {
      plot_lp_difference(
        fits[[engine]],
        profile = profile,
        variables = list(tx = c(0, 1)),
        times = 1:3
      )$data
    }
    expect_equal(plot_data(own_level), plot_data(full), info = engine)
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

# Profiles with their own factor levels, integer or unfitted previous states, and missing covariates: tests written before the fix, when each was a known bug.

# One orm and one vglm fit of the same first-order model, shared by the tests
# below. `bmi` is a numeric covariate, constant within each patient.
known_bug_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6101)
      data$bmi <- 20 + data$id %% 7
      formula <- ordered(y) ~ time + tx + yprev + bmi
      cache <<- list(
        yprev_levels = levels(data$yprev),
        orm = mostm(formula, data = data, id_var = "id", absorb = 6),
        vglm = mostm(formula, data = data, id_var = "id", absorb = 6, engine = "vglm")
      )
    }
    cache
  }
})

# The condition that `expr` raises, or NULL when it raises none. Checking the
# condition afterwards, rather than with `expect_error(class = )`, lets every
# engine in a loop report its own failure: an error of another class would
# otherwise end the test at the first engine.
known_bug_condition <- function(expr) {
  tryCatch(
    {
      expr
      NULL
    },
    error = function(e) e
  )
}

# Two hand-built starting profiles whose factor `yprev` has the fitted levels.
known_bug_profiles <- function() {
  data.frame(
    id = 1:2,
    tx = c(0, 1),
    time = 1,
    bmi = c(22, 25),
    yprev = factor(c(2, 3), levels = known_bug_fits()$yprev_levels)
  )
}

test_that("a profile whose factor `yprev` has only its own levels gives the SOPs of the fitted levels", {
  # The recursion builds the later previous states 1, 4 and 5 with the levels
  # of the profile's own `yprev` column, so they become NA. orm then stops
  # with an error about state levels missing from the fit, and vglm with
  # "PO predictor dimensions do not match the SOP state plan." A character
  # `yprev` already works; the same states with fewer factor levels should
  # too.
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- known_bug_fits()
  full <- known_bug_profiles()
  own_levels <- full
  own_levels$yprev <- factor(c(2, 3))

  for (engine in c("orm", "vglm")) {
    expected <- sops(fits[[engine]], newdata = full, times = 1:3)
    result <- tryCatch(
      sops(fits[[engine]], newdata = own_levels, times = 1:3),
      error = function(e) e
    )
    expect_false(
      inherits(result, "error"),
      info = paste(engine, if (inherits(result, "error")) conditionMessage(result))
    )
    if (!inherits(result, "error")) {
      expect_equal(result$estimate, expected$estimate, info = engine)
    }
  }
})

test_that("a profile whose `yprev` holds integers gives the SOPs of the fitted factor levels for every engine", {
  # The model reads `yprev` as a factor. orm turns the integers into its
  # fitted levels; vglm stops in stats::model.matrix() with "contrasts apply
  # only to factors".
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- known_bug_fits()
  full <- known_bug_profiles()
  integers <- full
  integers$yprev <- c(2L, 3L)

  for (engine in c("orm", "vglm")) {
    expected <- sops(fits[[engine]], newdata = full, times = 1:3)
    result <- tryCatch(
      sops(fits[[engine]], newdata = integers, times = 1:3),
      error = function(e) e
    )
    expect_false(
      inherits(result, "error"),
      info = paste(engine, if (inherits(result, "error")) conditionMessage(result))
    )
    if (!inherits(result, "error")) {
      expect_equal(result$estimate, expected$estimate, info = engine)
    }
  }
})

# A fake blrm fit (no Stan) whose design reads `tx`, `time` and `yprev` as
# numbers, with two hand-built profiles whose factor `yprev` has the fitted
# states as levels. `sformula` is the formula rmsb stores, from which mostr
# reads the columns the model uses.
known_bug_blrm <- function() {
  coefficients <- rbind(
    c(0.5, -0.5, 0, 1.2, 0.1),
    c(2, 1, 0, -1.1, 0.2),
    c(-0.4, -2, 0, 0.2, 0.3)
  )
  colnames(coefficients) <- c("y>=2", "y>=3", "tx", "time", "yprev")
  model <- fake_blrm(coefficients)
  model$sformula <- y ~ tx + time + yprev
  list(
    model = model,
    profiles = data.frame(
      id = c("a", "b"),
      tx = c(0, 1),
      yprev = factor(c(1, 2), levels = 1:3),
      time = 1
    )
  )
}

test_that("blrm SOPs of a profile whose factor `yprev` has only its own levels are those of the fitted states", {
  # As for orm and vglm, the later previous state 3 becomes NA, and blrm
  # stops with "BLRM prediction produced invalid probability mass".
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time", "yprev"))
  )
  fake <- known_bug_blrm()
  own_levels <- fake$profiles
  own_levels$yprev <- factor(c(1, 2))

  expected <- sops(fake$model, newdata = fake$profiles, times = 1:3, n_draws = 3, seed = 4)
  result <- tryCatch(
    sops(fake$model, newdata = own_levels, times = 1:3, n_draws = 3, seed = 4),
    error = function(e) e
  )
  expect_false(
    inherits(result, "error"),
    info = if (inherits(result, "error")) conditionMessage(result)
  )
  if (!inherits(result, "error")) {
    expect_equal(result$estimate, expected$estimate)
  }
})

test_that("blrm refuses a profile `yprev` that is not a fitted state, and a missing covariate, naming the column", {
  # A `yprev` of 9 is predicted as if it were a state, and a missing `tx`
  # stops with "BLRM prediction produced invalid probability mass", which
  # names no column.
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time", "yprev"))
  )
  fake <- known_bug_blrm()
  not_a_state <- fake$profiles
  not_a_state$yprev <- c("1", "9")
  missing_tx <- fake$profiles
  missing_tx$tx <- c(NA, 1)

  for (case in list(list(not_a_state, "yprev"), list(missing_tx, "tx"))) {
    for (fn in c("sops", "avg_sops")) {
      condition <- known_bug_condition(
        get(fn)(fake$model, newdata = case[[1]], times = 1:3, n_draws = 3, seed = 4)
      )
      expect_true(inherits(condition, "mostr_error"), info = paste(case[[2]], fn))
      expect_match(
        conditionMessage(condition %||% simpleCondition("")),
        case[[2]],
        info = paste(case[[2]], fn)
      )
    }
  }
})

test_that("a profile value that is not a fitted level gives a mostr_error naming the column for every engine", {
  # orm refuses through mostr_abort(); vglm stops in stats::model.frame() with
  # "factor yprev has new levels 9", a plain simpleError.
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- known_bug_fits()
  profiles <- known_bug_profiles()
  profiles$yprev <- c("2", "9")

  for (engine in c("orm", "vglm")) {
    condition <- known_bug_condition(
      sops(fits[[engine]], newdata = profiles, times = 1:3)
    )
    expect_true(inherits(condition, "mostr_error"), info = engine)
    expect_match(
      conditionMessage(condition %||% simpleCondition("")),
      "yprev",
      info = engine
    )
  }
})

test_that("a profile with a missing covariate gives a mostr_error naming the column for every engine", {
  # vglm builds the design with stats::model.matrix(), which drops the row
  # with NA. With two profiles, the probabilities of the other profile are
  # then spread over both rows, out of step with the states, and returned
  # without a warning; with one or three profiles, base R errors such as
  # "replacement has length zero" stop the call. orm
  # keeps the NA row and stops later, with an error about missing transition
  # probabilities and state levels that does not name the column.
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  fits <- known_bug_fits()
  profiles <- known_bug_profiles()
  profiles$bmi <- c(NA, 25)

  for (engine in c("orm", "vglm")) {
    for (fn in c("sops", "avg_sops")) {
      condition <- known_bug_condition(
        get(fn)(fits[[engine]], newdata = profiles, times = 1:3)
      )
      expect_true(
        inherits(condition, "mostr_error"),
        info = paste(engine, fn)
      )
      expect_match(
        conditionMessage(condition %||% simpleCondition("")),
        "bmi",
        info = paste(engine, fn)
      )
    }
  }
})
