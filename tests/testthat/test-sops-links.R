test_that("Markov workflows reject links other than the logit", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  cases <- list(
    list(
      label = "vglm with a probit link",
      link = "probitlink",
      message = paste0(
        "Only cumulative-logit models are supported; ",
        "refit the model with a logit link."
      ),
      error = function() {
        data <- suppressWarnings(
          make_test_data(n_patients = 80, seed = 841, follow_up_time = 5)
        )
        fit <- suppressWarnings(
          mostm(
            ordered(y) ~ time + tx + yprev,
            family = VGAM::cumulative(
              reverse = TRUE,
              parallel = TRUE,
              link = "probitlink"
            ),
            data = data,
            id_var = "id",
            absorb = 6,
            engine = "vglm"
          )
        )
        tryCatch(mostr:::validate_markov_model(fit), error = identity)
      }
    ),
    list(
      label = "orm with the probit family",
      link = "probit",
      error = function() {
        data <- suppressWarnings(
          make_test_data(n_patients = 80, seed = 842, follow_up_time = 5)
        )
        # `mostm()` stops before computing the patient-clustered covariance.
        tryCatch(
          suppressWarnings(
            mostm(
              ordered(y) ~ time + tx + yprev,
              family = "probit",
              data = data,
              id_var = "id",
              absorb = 6
            )
          ),
          error = identity
        )
      }
    )
  )

  for (case in cases) {
    error <- case$error()
    expect_true(inherits(error, "mostr_unsupported_link"), info = case$label)
    expect_equal(error$link, case$link, info = case$label)
    if (!is.null(case$message)) {
      expect_equal(conditionMessage(error), case$message, info = case$label)
    }
  }
})

test_that("Markov workflows reject raw and manually robust backend fits", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- suppressWarnings(
    make_test_data(n_patients = 40, seed = 844, follow_up_time = 5)
  )
  raw_fit <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  robust_fit <- robcov_vglm(raw_fit, cluster = data$id)
  raw_orm_fit <- rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    x = TRUE,
    y = TRUE
  )

  expect_error(
    mostr:::validate_markov_model(raw_fit),
    "Refit this model with `mostm(engine = \"vglm\")`",
    fixed = TRUE
  )
  expect_error(
    mostr:::validate_markov_model(robust_fit),
    "Refit this model with `mostm(engine = \"vglm\")`",
    fixed = TRUE
  )
  expect_error(
    mostr:::validate_markov_model(raw_orm_fit),
    "Refit this model with `mostm(engine = \"orm\")`",
    fixed = TRUE
  )
})
