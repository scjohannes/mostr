# A response `Ocens(lower, upper)` gives a range of possible states per row.
# The fits are compared with the `y = NA` fit and with the plain-R reference
# in helper-missing-states.R extended to ranges; these comparisons hold for any
# data.

ms_ocens_columns <- function(data) {
  y <- as.integer(as.character(data$y))
  data$lower <- y
  data$upper <- y
  data
}

ms_ocens_fit <- function(data, ...) {
  ms_fit(data, formula = Ocens(lower, upper) ~ yprev + time + tx, ...)
}

test_that("the full non-absorbing range equals a missing state", {
  data <- ms_ocens_columns(ms_fixture("gaps"))
  fit_na <- ms_fixture("gap_fit")
  expect_same_fit <- function(fit) {
    expect_equal(coef(fit), coef(fit_na), tolerance = 1e-8)
    expect_equal(
      vcov(fit, intercepts = "all"),
      vcov(fit_na, intercepts = "all"),
      tolerance = 1e-8
    )
    expect_equal(
      attr(fit, "markov_path_fit")$loglik,
      attr(fit_na, "markov_path_fit")$loglik,
      tolerance = 1e-10
    )
  }
  # Every unknown visit coded as states 1 to 3.
  ranged <- data
  ranged$lower[is.na(ranged$y)] <- 1
  ranged$upper[is.na(ranged$y)] <- 3
  expect_same_fit(ms_ocens_fit(ranged))
  skip_on_cran() # Slow: one more fit.
  # Both bounds missing.
  expect_same_fit(ms_ocens_fit(data))
})

test_that("infinite bounds equal the lowest and highest state", {
  skip_on_cran() # Slow: two fits; the full-range check runs on CRAN.
  data <- ms_ocens_columns(ms_fixture("gaps"))
  unknown <- which(is.na(data$y))
  # Interior unknown visits are known to be in states 1-2 or in 2-4 (4 is
  # death: a range that includes it at a visit followed by a live state
  # admits only the live states).
  low <- unknown[c(TRUE, FALSE)]
  high <- unknown[c(FALSE, TRUE)]
  finite <- data
  finite$lower[low] <- 1
  finite$upper[low] <- 2
  finite$lower[high] <- 2
  finite$upper[high] <- 4
  infinite <- finite
  infinite$lower[low] <- -Inf
  infinite$upper[high] <- Inf

  fit_finite <- ms_ocens_fit(finite)
  fit_infinite <- ms_ocens_fit(infinite)
  expect_identical(coef(fit_infinite), coef(fit_finite))
  expect_identical(fit_infinite$var, fit_finite$var)

  reference <- ms_reference(
    finite,
    ranges = ms_row_ranges(finite, lower = infinite$lower, upper = infinite$upper)
  )
  expect_equal(
    attr(fit_finite, "markov_path_fit")$loglik,
    reference(coef(fit_finite)),
    tolerance = 1e-8
  )
  ms_expect_maximum(fit_finite, reference)
})

test_that("a narrower range changes the fit and matches the reference", {
  skip_on_cran() # Slow: a fit and the reference.
  data <- ms_ocens_columns(ms_fixture("gaps"))
  fit_na <- ms_fixture("gap_fit")
  narrow <- data
  unknown <- which(is.na(data$y))
  narrow$lower[unknown] <- 2
  narrow$upper[unknown] <- 3
  fit_narrow <- ms_ocens_fit(narrow)
  expect_gt(max(abs(coef(fit_narrow) - coef(fit_na))), 1e-4)

  ranges <- ms_row_ranges(narrow, lower = narrow$lower, upper = narrow$upper)
  expect_identical(ranges$lo[unknown], rep(2L, length(unknown)))
  reference <- ms_reference(narrow, ranges = ranges)
  expect_equal(
    attr(fit_narrow, "markov_path_fit")$loglik,
    reference(coef(fit_narrow)),
    tolerance = 1e-8
  )
  ms_expect_maximum(fit_narrow, reference)
  # The patient whose first two rows have no previous state: a range on the
  # first row is not used either, because nothing precedes it.
  first <- attr(data, "special")$unknown_start
  per_patient <- reference(coef(fit_narrow), per_patient = TRUE)
  rows <- narrow[narrow$id == first, ]
  expect_equal(
    per_patient[[as.character(first)]],
    ms_reference(rows[-(1:2), ])(coef(fit_narrow)),
    tolerance = 1e-12
  )
})

test_that("exact ranges everywhere give the ordinary fit", {
  skip_on_cran() # Slow: two fits and predictions.
  data <- ms_ocens_columns(ms_fixture("complete"))
  fit_ocens <- ms_ocens_fit(data)
  fit_plain <- ms_fit(data)
  expect_null(attr(fit_ocens, "markov_path_fit"))
  expect_equal(coef(fit_ocens), coef(fit_plain), tolerance = 1e-12)
  expect_equal(
    vcov(fit_ocens, intercepts = "all"),
    vcov(fit_plain, intercepts = "all"),
    tolerance = 1e-12
  )
  expect_equal(
    as.data.frame(avg_sops(fit_ocens, variables = list(tx = c(0, 1)), times = 1:4))$estimate,
    as.data.frame(avg_sops(fit_plain, variables = list(tx = c(0, 1)), times = 1:4))$estimate,
    tolerance = 1e-12
  )
})

test_that("a single interval-censored row with a known previous state works", {
  skip_on_cran() # Slow: a fit and the reference.
  data <- ms_ocens_columns(ms_fixture("complete"))
  special <- attr(ms_gap_data(data), "special")
  # At the patient's last row: a stretch of one step.
  last_patient <- special$alive_long[1]
  last_row <- max(which(data$id == last_patient))
  one <- data
  one$lower[last_row] <- 2
  one$upper[last_row] <- 3
  # At an interior row: the next row's previous state is unknown.
  inner_patient <- special$alive_long[2]
  inner_row <- which(data$id == inner_patient & data$time == 3)
  one$lower[inner_row] <- 1
  one$upper[inner_row] <- 2
  one$yprev[inner_row + 1] <- NA

  fit <- ms_ocens_fit(one)
  info <- attr(fit, "markov_path_fit")
  expect_identical(info$n_stretches, 2L)
  expect_identical(info$n_unobserved_visits, 3L)
  reference <- ms_reference(
    one,
    ranges = ms_row_ranges(one, lower = one$lower, upper = one$upper)
  )
  expect_equal(info$loglik, reference(coef(fit)), tolerance = 1e-8)
  ms_expect_maximum(fit, reference)
})

test_that("`Ocens()` refuses incomplete or reversed bounds and vglm", {
  data <- ms_ocens_columns(ms_fixture("gaps"))
  one_sided <- data
  one_sided$upper[which(is.na(data$y))[1]] <- 3
  expect_error(
    ms_ocens_fit(one_sided),
    "missing bound on 1 row where the other bound is given",
    class = "mostr_error"
  )
  reversed <- data
  reversed$lower[1] <- 3
  reversed$upper[1] <- 2
  expect_error(
    ms_ocens_fit(reversed),
    "needs `lower <= upper`",
    class = "mostr_error"
  )
  expect_error(
    ms_ocens_fit(data, engine = "vglm"),
    "is available only for `engine = \"orm\"` and `engine = \"blrm\"`",
    class = "mostr_error",
    fixed = TRUE
  )
  expect_error(
    ms_fit(data, formula = Ocens(lower) ~ yprev + time + tx),
    "with the lowest and the highest state",
    class = "mostr_error"
  )
})
