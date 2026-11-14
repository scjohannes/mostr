# Patients whose rows are all left out of the fit (a predictor is missing on
# every row) are not part of the fitted cohort. `mostm()` says so when it fits,
# and the bootstrap leaves them out of each resample's cohort, as the point
# estimate does.

dropped_patients_data <- function() {
  data <- make_test_data(n_patients = 80, follow_up_time = 6, seed = 9461)
  ids <- unique(data$id)
  set.seed(4127)
  ages <- stats::setNames(round(stats::rnorm(length(ids), 60, 10)), ids)
  data$age <- unname(ages[as.character(data$id)])
  data
}

fit_dropped_patients <- function(data) {
  mostm(
    y ~ time + tx + age + yprev,
    data = data,
    id_var = "id",
    absorb = 6
  )
}

test_that("mostm() warns about patients left out of the fit and names them", {
  skip_if_not_installed("rms")
  data <- dropped_patients_data()
  ids <- unique(data$id)
  dropped <- ids[c(2, 5, 9)]
  data$age[data$id %in% dropped] <- NA
  # A patient who loses one later row stays in the fit and is not named.
  partial <- ids[12]
  partial_rows <- which(data$id == partial)
  expect_gt(length(partial_rows), 1L)
  data$age[partial_rows[2L]] <- NA

  warning <- expect_warning(
    fit_dropped_patients(data),
    class = "mostr_warning_dropped_patients"
  )
  for (id in dropped) {
    expect_match(conditionMessage(warning), as.character(id), fixed = TRUE)
  }
  expect_match(conditionMessage(warning), "3 patients", fixed = TRUE)
  expect_setequal(as.character(warning$ids), as.character(dropped))
})

test_that("mostm() gives no dropped-patient warning when every patient is fitted", {
  skip_if_not_installed("rms")
  data <- dropped_patients_data()
  expect_no_warning(
    fit_dropped_patients(data),
    class = "mostr_warning_dropped_patients"
  )
})

test_that("the bootstrap leaves out patients the fit left out", {
  skip_if_not_installed("rms")
  skip_on_cran()
  data <- dropped_patients_data()
  dropped <- unique(data$id)[c(2, 5, 9)]
  data$age[data$id %in% dropped] <- NA
  fit <- suppressWarnings(
    fit_dropped_patients(data),
    classes = "mostr_warning_dropped_patients"
  )
  point <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:3)

  # Before the fix every draw failed: the resampled cohort held the dropped
  # patients, whose `age` is missing, and prediction refused their rows.
  for (method in c("bootstrap", "fwb")) {
    result <- expect_no_warning(
      inferences(point, method = method, n_draws = 6, seed = 311)
    )
    expect_equal(result$estimate, point$estimate)
    expect_false(anyNA(result$conf.low))
    expect_length(unique(attr(result, "draws")$draw_id), 6L)
  }
})
