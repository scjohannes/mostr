test_that("the vglm engine supplies rms formula helpers and keeps %ia% semantics", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2468, follow_up_time = 12)
  baseline <- data[data$time == 1, , drop = FALSE]
  # `%ia%` of a linear term gives one interaction column per non-reference
  # previous state; of a spline basis, more, including nonlinear (time') ones.
  n_linear <- nlevels(data$yprev) - 1

  cases <- list(
    list(
      label = "%ia% and rcs() without rms in the formula environment",
      formula = "ordered(y) ~ rcs(time, 4) + tx + yprev + time %ia% yprev",
      pattern = "%ia%",
      spline = FALSE,
      check = function(fit, ia_terms, label) {
        expect_equal(
          grep("time * yprev=", ia_terms, fixed = TRUE),
          seq_along(ia_terms),
          info = label
        )
        out <- mostr::soprob_markov(
          fit,
          newdata = baseline[1:3, , drop = FALSE],
          times = 1:4
        )
        expect_equal(dim(out), c(3L, 4L, 6L), info = label)
        expect_equal(
          unname(rowSums(out[, 4, ])),
          rep(1, 3),
          tolerance = 1e-10,
          info = label
        )
      }
    ),
    list(
      label = "an explicit rms::%ia% call",
      formula = paste(
        "ordered(y) ~ rms::rcs(time, 4) + tx + yprev +",
        "rms::`%ia%`(time, yprev)"
      ),
      pattern = "rms::`%ia%`",
      spline = FALSE
    ),
    list(
      label = "rcs(time, 4) %ia% yprev",
      formula = "ordered(y) ~ rcs(time, 4) + tx + yprev + rcs(time, 4) %ia% yprev",
      pattern = "%ia%",
      spline = TRUE
    )
  )

  for (case in cases) {
    form <- stats::as.formula(case$formula, env = baseenv())
    expect_false(
      exists("%ia%", envir = environment(form), inherits = TRUE),
      info = case$label
    )
    expect_false(
      exists("rcs", envir = environment(form), inherits = TRUE),
      info = case$label
    )

    fit <- suppressWarnings(
      mostr::mostm(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        absorb = 6,
        engine = "vglm"
      )
    )
    ia_terms <- grep(
      case$pattern,
      names(fit$vglm_fit@constraints),
      value = TRUE,
      fixed = TRUE
    )
    spline_terms <- grep("time'", ia_terms, fixed = TRUE)
    if (case$spline) {
      expect_gt(length(ia_terms), n_linear, label = case$label)
      expect_gt(length(spline_terms), 0, label = case$label)
    } else {
      expect_length(ia_terms, n_linear)
      expect_equal(spline_terms, integer(0), info = case$label)
    }
    if (!is.null(case$check)) {
      case$check(fit, ia_terms, case$label)
    }
  }
})

test_that("the vglm engine accepts family names and requires a data frame", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3),
    w = c(1, 2, 1, 2, 1, 2, 1, 2)
  )
  env <- list2env(data, parent = baseenv())
  form <- stats::as.formula("y ~ x", env = env)

  fit_from_name <- suppressWarnings(
    mostm(
      form,
      family = "binomialff",
      data = data,
      id_var = "id",
      engine = "vglm"
    )
  )
  expect_error(
    mostm(form, family = "binomialff", id_var = "id", engine = "vglm"),
    "`data` must be supplied as a data frame.",
    fixed = TRUE
  )
  weighted <- suppressWarnings(
    mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      weights = w,
      smart = FALSE,
      engine = "vglm"
    )
  )

  expect_s4_class(fit_from_name$vglm_fit, "vglm")
  weighted <- weighted$vglm_fit
  expect_s4_class(weighted, "vglm")
  expect_false(isTRUE(weighted@smart.prediction$smart.arg))
  expect_equal(nrow(weighted@prior.weights), nrow(data))
  expect_equal(as.vector(weighted@prior.weights), data$w)
})

test_that("the vglm engine rejects invalid families, weights, and form2 subsets", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(0, 1, 0, 1),
    x = c(0, 0, 1, 1),
    bad_w = c(1, 1, -1, 1)
  )

  expect_error(
    suppressWarnings(mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      weights = bad_w,
      engine = "vglm"
    )),
    "negative weights",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      mostm(
        y ~ x,
        family = list(),
        data = data,
        id_var = "id",
        engine = "vglm"
      )
    ),
    "is not a VGAM family function",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      form2 = ~x,
      subset = data$x == 1,
      engine = "vglm"
    )),
    "subset",
    fixed = TRUE
  )
})

test_that("the vglm engine records NA actions after dropping incomplete rows", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    id = c(1, 1, 2, 2, 3, 3),
    time = c(1, 2, 1, 2, 1, 2),
    y = c(0, 1, 0, 1, NA, 1),
    x = c(0, 0, 1, 1, 2, 2)
  )
  fit_binomial <- function(formula) {
    mostm(
      formula,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      engine = "vglm"
    )
  }

  expect_error(suppressWarnings(fit_binomial(y ~ 0)), "model matrix")
  with_na <- suppressWarnings(fit_binomial(y ~ x))$vglm_fit

  expect_s4_class(with_na, "vglm")
  expect_length(with_na@na.action, 1)
})

test_that("the vglm engine aligns etastart with the rows it fits", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6110)
  full <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorb = 6,
    engine = "vglm"
  )
  # Starting values have one row per row of `data`.
  start <- full$vglm_fit@predictors
  expect_equal(nrow(start), nrow(data))

  # Rows removed by `subset` lose their starting values too, so the fit
  # equals one on the kept rows with their starting values.
  kept <- data$time <= 4
  subset_fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorb = 6,
    engine = "vglm",
    subset = time <= 4,
    etastart = start
  )
  kept_rows <- data[kept, , drop = FALSE]
  kept_start <- start[kept, , drop = FALSE]
  kept_fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = kept_rows,
    id_var = "id",
    absorb = 6,
    engine = "vglm",
    etastart = kept_start
  )
  expect_equal(stats::coef(subset_fit), stats::coef(kept_fit), tolerance = 1e-8)

  # So do rows removed because of a missing predictor.
  incomplete <- data
  incomplete$tx[which(incomplete$time > 2)[c(2, 9)]] <- NA
  complete <- !is.na(incomplete$tx)
  incomplete_fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = incomplete,
    id_var = "id",
    absorb = 6,
    engine = "vglm",
    etastart = start
  )
  complete_rows <- incomplete[complete, , drop = FALSE]
  complete_start <- start[complete, , drop = FALSE]
  complete_fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = complete_rows,
    id_var = "id",
    absorb = 6,
    engine = "vglm",
    etastart = complete_start
  )
  expect_equal(
    stats::coef(incomplete_fit),
    stats::coef(complete_fit),
    tolerance = 1e-8
  )
})

test_that("the vglm engine stores form2 design slots", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3),
    z = c(1, 2, 1, 2, 1, 2, 1, 2)
  )

  fit <- suppressWarnings(
    mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      form2 = y ~ z,
      qr.arg = TRUE,
      engine = "vglm"
    )
  )$vglm_fit

  expect_s4_class(fit, "vglm")
  expect_true(length(fit@Xm2) > 0)
  expect_true(length(fit@Ym2) > 0)
  expect_equal(deparse(fit@misc$form2), deparse(y ~ z))
  expect_true(length(fit@callXm2) > 0)
})

test_that("the vglm engine validates form2 row alignment", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3)
  )
  env <- list2env(
    list(z = 1:3, z_y = c(0, 1, 0)),
    parent = baseenv()
  )

  expect_error(
    suppressWarnings(mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      form2 = stats::as.formula("z_y ~ z", env = env),
      engine = "vglm"
    )),
    "number of rows of 'y' and 'Ym2' are unequal",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      id_var = "id",
      form2 = stats::as.formula("~ z", env = env),
      engine = "vglm"
    )),
    "number of rows of 'x' and 'Xm2' are unequal",
    fixed = TRUE
  )
})
