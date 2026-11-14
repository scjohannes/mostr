# rms codes an ordered factor with numeric levels as "scored": a column with
# the level's number and indicators from the third level on (`yo`, `yo=3`).
# A plain factor with the same levels spans the same model with treatment
# indicators (`yprev=2`, `yprev=3`), so the two fits must give the same
# likelihood and the same predictions. Each test compares two routes to one
# number on the same data, so it passes for any seed.

# The standard simulated data with the previous state also as an ordered
# factor, `yo`, and a patient-level ordered covariate `sev` (levels 1..3)
# next to the same values as a plain factor, `sevf`.
ordered_factor_data <- function(data) {
  data$yo <- factor(as.character(data$yprev), levels = 1:3, ordered = TRUE)
  sev <- (data$id %% 3) + 1
  data$sev <- factor(sev, levels = 1:3, ordered = TRUE)
  data$sevf <- factor(sev, levels = 1:3)
  data
}

test_that("orm design rows of an ordered factor are rms's own", {
  data <- ordered_factor_data(ms_fixture("complete"))
  fit <- ms_fit(data, formula = ordered(y) ~ yo + time + tx + sev, p_var = "yo")
  # The scored coding, not treatment indicators.
  expect_identical(fit$Design$assume.code[1:2], c(7L, 1L))
  rows <- data[!is.na(data$yo), ]
  expect_equal(
    unname(mostr:::orm_model_matrix(fit, rows)),
    unname(fit$x),
    tolerance = 0
  )
  # One row, and rows in which some levels do not occur. The blrm route
  # goes through rms::predictrms(), which needs every level in its rows; it
  # reads only the Design of the fit, so the orm fit stands in for a blrm
  # fit here.
  subsets <- list(1L, which(rows$yo == "3"), which(rows$yo != "1"))
  for (i in subsets) {
    expect_equal(
      unname(mostr:::orm_model_matrix(fit, rows[i, ])),
      unname(fit$x[i, , drop = FALSE]),
      tolerance = 0
    )
    expect_equal(
      unname(mostr:::blrm_design_matrix(fit, rows[i, ])),
      unname(fit$x[i, , drop = FALSE]),
      tolerance = 0
    )
  }
  # The same columns as the rows rms fitted, for every level.
  one <- rows[1, ]
  for (k in 1:3) {
    one$yo <- factor(k, levels = 1:3, ordered = TRUE)
    x <- mostr:::orm_model_matrix(fit, one)
    expect_equal(unname(x[1, c("yo", "yo=3")]), c(k, k == 3))
  }
})

test_that("an ordered previous state gives the same SOPs as a plain factor", {
  skip_on_cran() # Slow (about 1.3 s); the design rows it rests on are checked on CRAN above.
  data <- ordered_factor_data(ms_fixture("complete"))
  plain <- ms_fit(data)
  ordered <- ms_fit(data, formula = ordered(y) ~ yo + time + tx, p_var = "yo")
  expect_equal(ms_loglik(ordered), ms_loglik(plain), tolerance = 1e-10)
  start <- data[data$time == 1, ]
  expect_equal(
    avg_sops(ordered, newdata = start, times = 1:4)$estimate,
    avg_sops(plain, newdata = start, times = 1:4)$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    sops(ordered, newdata = start[1:3, ], times = 1:4)$estimate,
    sops(plain, newdata = start[1:3, ], times = 1:4)$estimate,
    tolerance = 1e-8
  )
})

test_that("an ordered covariate gives the same SOPs as a plain factor", {
  data <- ordered_factor_data(ms_fixture("complete"))
  plain <- ms_fit(data, formula = ordered(y) ~ yprev + time + tx + sevf)
  ordered <- ms_fit(data, formula = ordered(y) ~ yprev + time + tx + sev)
  expect_equal(ms_loglik(ordered), ms_loglik(plain), tolerance = 1e-10)
  start <- data[data$time == 1, ]
  expect_equal(
    avg_sops(ordered, newdata = start, times = 1:4)$estimate,
    avg_sops(plain, newdata = start, times = 1:4)$estimate,
    tolerance = 1e-8
  )
})

test_that("a path fit with an ordered previous state maximises the same likelihood", {
  data <- ordered_factor_data(ms_fixture("gaps"))
  plain <- ms_fixture("gap_fit")
  ordered <- ms_fit(data, formula = ordered(y) ~ yo + time + tx, p_var = "yo")
  reference <- ms_reference(data)
  expect_equal(ms_loglik(ordered), reference(coef(plain)), tolerance = 1e-8)
  expect_equal(
    attr(ordered, "markov_path_fit")$loglik,
    attr(plain, "markov_path_fit")$loglik,
    tolerance = 1e-8
  )
  # The slopes of the other predictors do not depend on the coding.
  expect_equal(
    coef(ordered)[c("time", "tx")],
    coef(plain)[c("time", "tx")],
    tolerance = 1e-5
  )
})
