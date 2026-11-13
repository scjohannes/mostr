# Known bugs, each pinned by a test that fails until it is fixed. Each test
# names the bug and states the behaviour expected once it is fixed. They live
# in one file so that fixes on parallel branches do not touch the same lines;
# move a test next to the code it checks once its bug is fixed.

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

test_that("known bug: a profile whose factor `yprev` has only its own levels gives the SOPs of the fitted levels", {
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

test_that("known bug: a profile whose `yprev` holds integers gives the SOPs of the fitted factor levels for every engine", {
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

test_that("known bug: blrm SOPs of a profile whose factor `yprev` has only its own levels are those of the fitted states", {
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

test_that("known bug: blrm refuses a profile `yprev` that is not a fitted state, and a missing covariate, naming the column", {
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

test_that("known bug: a profile value that is not a fitted level gives a mostr_error naming the column for every engine", {
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

test_that("known bug: a profile with a missing covariate gives a mostr_error naming the column for every engine", {
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

test_that("known bug: orm and vglm accept or refuse a covariance by the same rule, the orm rule", {
  # robcov_orm() (orm_model_bread()) checks the correlation matrix: it refuses
  # when the smallest eigenvalue is at most p * 1e-10 times the largest, the
  # reason given in its comment. The vglm check (validate_vglm_matrix()) uses
  # the raw matrix and a tolerance of 100 * machine epsilon * p * max(1, largest
  # eigenvalue). So vglm accepts a correlation nearly as singular as 1e-11,
  # and refuses a well-conditioned covariance whose coefficients are on a
  # small scale (variances 1e-16), which orm accepts.
  coefficients <- c("y>=2" = 0.4, "y>=3" = -0.6, x = 0.3)
  coefficient_names <- names(coefficients)
  accepts <- function(engine, covariance) {
    dimnames(covariance) <- list(coefficient_names, coefficient_names)
    tryCatch(
      {
        if (engine == "orm") {
          fit <- structure(
            list(
              coefficients = coefficients,
              orig.var = covariance,
              non.slopes = 2L
            ),
            class = "orm"
          )
          orm_model_bread(fit)
        } else {
          validate_vglm_matrix(
            covariance,
            coefficient_names,
            "bread",
            positive_definite = TRUE
          )
        }
        TRUE
      },
      mostr_error = function(e) FALSE
    )
  }

  # Unit variances; eigenvalues 2, 1 - 1e-11 and 1e-11.
  rotation <- qr.Q(qr(matrix(c(1, 2, 3, 4, 5, 7, 2, 1, 9), 3)))
  nearly_singular <- rotation %*% diag(c(2, 1 - 1e-11, 1e-11)) %*% t(rotation)
  nearly_singular <- (nearly_singular + t(nearly_singular)) / 2
  small_scale <- diag(1e-16, 3)

  expect_false(accepts("orm", nearly_singular))
  expect_true(accepts("orm", small_scale))
  expect_identical(
    accepts("vglm", nearly_singular),
    accepts("orm", nearly_singular)
  )
  expect_identical(accepts("vglm", small_scale), accepts("orm", small_scale))
})
