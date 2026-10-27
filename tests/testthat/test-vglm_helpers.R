test_that("the vglm engine supplies rms formula helpers for %ia%", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2468, follow_up_time = 12)
  baseline <- data[data$time == 1, , drop = FALSE]
  stripped_env <- baseenv()

  form <- stats::as.formula(
    "ordered(y) ~ rcs(time, 4) + tx + yprev + time %ia% yprev",
    env = stripped_env
  )
  expect_equal(
    exists("%ia%", envir = environment(form), inherits = TRUE),
    FALSE
  )
  expect_equal(exists("rcs", envir = environment(form), inherits = TRUE), FALSE)

  fit <- suppressWarnings(
    mostr::mostm(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      engine = "vglm"
    )
  )

  ia_terms <- grep("%ia%", names(fit@constraints), value = TRUE, fixed = TRUE)
  expect_length(ia_terms, nlevels(data$yprev) - 1)
  expect_equal(grep("time'", ia_terms, fixed = TRUE), integer(0))
  expect_equal(
    grep("time * yprev=", ia_terms, fixed = TRUE),
    seq_along(ia_terms)
  )

  out <- mostr::soprob_markov(
    fit,
    newdata = baseline[1:3, , drop = FALSE],
    times = 1:4,
    y_levels = 1:6,
    absorb = "6"
  )

  expect_equal(dim(out), c(3L, 4L, 6L))
  expect_equal(unname(rowSums(out[, 4, ])), rep(1, 3), tolerance = 1e-10)
})

test_that("the vglm engine keeps explicit rms::%ia% calls working", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2469, follow_up_time = 12)
  form <- stats::as.formula(
    paste(
      "ordered(y) ~ rms::rcs(time, 4) + tx + yprev +",
      "rms::`%ia%`(time, yprev)"
    ),
    env = baseenv()
  )

  fit <- suppressWarnings(
    mostr::mostm(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      engine = "vglm"
    )
  )

  ia_terms <- grep(
    "rms::`%ia%`",
    names(fit@constraints),
    value = TRUE,
    fixed = TRUE
  )
  expect_length(ia_terms, nlevels(data$yprev) - 1)
  expect_equal(grep("time'", ia_terms, fixed = TRUE), integer(0))
})

test_that("the vglm engine preserves rms %ia% semantics for spline bases", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2470, follow_up_time = 12)
  form <- stats::as.formula(
    "ordered(y) ~ rcs(time, 4) + tx + yprev + rcs(time, 4) %ia% yprev",
    env = baseenv()
  )

  fit <- suppressWarnings(
    mostr::mostm(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      engine = "vglm"
    )
  )

  ia_terms <- grep("%ia%", names(fit@constraints), value = TRUE, fixed = TRUE)
  expect_gt(length(ia_terms), nlevels(data$yprev) - 1)
  expect_gt(length(grep("time'", ia_terms, fixed = TRUE)), 0)
})

test_that("the vglm engine accepts family names and requires a data frame", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3),
    w = c(1, 2, 1, 2, 1, 2, 1, 2)
  )
  env <- list2env(data, parent = baseenv())
  form <- stats::as.formula("y ~ x", env = env)

  fit_from_name <- suppressWarnings(
    mostm(form, family = "binomialff", data = data, engine = "vglm")
  )
  expect_error(
    mostm(form, family = "binomialff", engine = "vglm"),
    "`data` must be supplied as a data frame.",
    fixed = TRUE
  )
  weighted <- suppressWarnings(
    mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      weights = w,
      smart = FALSE,
      engine = "vglm"
    )
  )

  expect_s4_class(fit_from_name, "vglm")
  expect_s4_class(weighted, "vglm")
  expect_false(isTRUE(weighted@smart.prediction$smart.arg))
  expect_equal(nrow(weighted@prior.weights), nrow(data))
  expect_equal(as.vector(weighted@prior.weights), data$w)
})

test_that("the vglm engine rejects invalid families, weights, and form2 subsets", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1),
    x = c(0, 0, 1, 1),
    off = c(0.1, 0.2, 0.3, 0.4),
    bad_w = c(1, 1, -1, 1)
  )

  expect_error(
    suppressWarnings(mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      weights = bad_w,
      engine = "vglm"
    )),
    "negative weights",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      mostm(y ~ x, family = list(), data = data, engine = "vglm")
    ),
    "is not a VGAM family function",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      form2 = ~x,
      subset = data$x == 1,
      engine = "vglm"
    )),
    "subset",
    fixed = TRUE
  )
  expect_error(
    mostm(
      y ~ x + offset(off),
      family = VGAM::binomialff,
      data = data,
      engine = "vglm"
    ),
    "Model offsets are not supported",
    fixed = TRUE
  )
  expect_error(
    mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      offset = off,
      engine = "vglm"
    ),
    "Model offsets are not supported",
    fixed = TRUE
  )
})

test_that("the vglm engine records NA actions after dropping incomplete rows", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1, NA, 1),
    x = c(0, 0, 1, 1, 2, 2)
  )

  expect_error(
    suppressWarnings(
      mostm(y ~ 0, family = VGAM::binomialff, data = data, engine = "vglm")
    ),
    "model matrix"
  )
  with_na <- suppressWarnings(
    mostm(y ~ x, family = VGAM::binomialff, data = data, engine = "vglm")
  )

  expect_s4_class(with_na, "vglm")
  expect_length(with_na@na.action, 1)
})

test_that("the vglm engine aligns etastart with the rows it fits", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6110)
  full <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
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
    engine = "vglm",
    etastart = start
  )
  complete_rows <- incomplete[complete, , drop = FALSE]
  complete_start <- start[complete, , drop = FALSE]
  complete_fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = complete_rows,
    id_var = "id",
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
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3),
    z = c(1, 2, 1, 2, 1, 2, 1, 2)
  )

  fit <- suppressWarnings(
    mostm(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      form2 = y ~ z,
      qr.arg = TRUE,
      engine = "vglm"
    )
  )

  expect_s4_class(fit, "vglm")
  expect_true(length(fit@Xm2) > 0)
  expect_true(length(fit@Ym2) > 0)
  expect_equal(deparse(fit@misc$form2), deparse(y ~ z))
  expect_true(length(fit@callXm2) > 0)
})

test_that("the vglm engine validates form2 row alignment", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
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
      form2 = stats::as.formula("~ z", env = env),
      engine = "vglm"
    )),
    "number of rows of 'x' and 'Xm2' are unequal",
    fixed = TRUE
  )
})

test_that("add_rms_formula_helpers() returns non-formulas unchanged", {
  expect_equal(mostr:::add_rms_formula_helpers(1), 1)

  form <- y ~ x
  environment(form) <- NULL
  out <- mostr:::add_rms_formula_helpers(form)

  expect_s3_class(out, "formula")
  expect_equal(exists("rcs", envir = environment(out), inherits = TRUE), TRUE)
  expect_equal(exists("lsp", envir = environment(out), inherits = TRUE), TRUE)
})
