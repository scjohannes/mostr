# With every state observed, `missing_states` changes nothing: mostm() takes
# the ordinary code path whatever its value. No random numbers are drawn after
# the data are simulated, and the comparisons hold for any data.

ms_avg_sops_numbers <- function(fit) {
  out <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:6)
  out <- as.data.frame(out)
  attributes(out) <- attributes(out)[c("names", "row.names", "class")]
  out
}

test_that("complete data give the same fit for every `missing_states`", {
  data <- ms_fixture("complete")
  fit_default <- ms_fit(data)
  fit_average <- ms_fit(data, missing_states = "average")
  fit_drop <- ms_fit(data, missing_states = "drop")

  expect_null(attr(fit_default, "markov_path_fit"))
  expect_null(attr(fit_average, "markov_path_fit"))
  expect_null(attr(fit_drop, "markov_path_fit"))

  expect_identical(ms_numbers(fit_average), ms_numbers(fit_default))
  expect_identical(ms_numbers(fit_drop), ms_numbers(fit_default))
  skip_on_cran() # Slow: six sets of predictions.
  expect_identical(
    ms_avg_sops_numbers(fit_average),
    ms_avg_sops_numbers(fit_default)
  )
  expect_identical(
    ms_avg_sops_numbers(fit_drop),
    ms_avg_sops_numbers(fit_default)
  )
})

test_that("complete data give the same vglm fit with `missing_states = \"drop\"`", {
  skip_on_cran() # Slow: two vglm fits.
  data <- ms_fixture("complete")
  fit <- ms_fit(data, engine = "vglm")
  fit_drop <- ms_fit(data, engine = "vglm", missing_states = "drop")
  expect_identical(stats::coef(fit_drop), stats::coef(fit))
  expect_identical(stats::vcov(fit_drop), stats::vcov(fit))
})

test_that("`missing_states = \"drop\"` fits only the fully observed rows", {
  skip_on_cran() # Slow: two fits.
  data <- ms_fixture("gaps")
  fit_drop <- ms_fit(data, missing_states = "drop")
  expect_null(attr(fit_drop, "markov_path_fit"))
  observed <- data[!is.na(data$y) & !is.na(data$yprev), ]
  fit_observed <- ms_fit(observed)
  expect_identical(stats::coef(fit_drop), stats::coef(fit_observed))
  expect_identical(fit_drop$var, fit_observed$var)
})
