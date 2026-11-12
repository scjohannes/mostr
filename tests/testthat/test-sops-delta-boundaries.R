test_that("logit delta intervals distinguish structural boundaries", {
  skip_on_cran()
  local_reproducible_output(width = 80)
  # The warning for an estimate of exactly 0 or 1 that is not a structural
  # boundary.
  expect_snapshot(
    nonstructural <- delta_interval_bounds(
      estimate = c(0, 1),
      standard_error = c(0.1, 0.1),
      conf_level = 0.95,
      conf_type = "logit"
    )
  )

  # Zero or tiny uncertainty does not certify a structural boundary either:
  # logit limits stay undefined, while Wald limits are computed as usual.
  z <- stats::qnorm(0.975)
  tiny_estimate <- rep(c(0, 1), 3)
  tiny_se <- rep(c(0, 1e-20, 1e-7), each = 2)
  cases <- list(
    list(
      label = "logit, standard error 0.1",
      estimate = c(0, 1),
      standard_error = c(0.1, 0.1),
      conf_type = "logit",
      warns = TRUE,
      conf.low = rep(NA_real_, 2),
      conf.high = rep(NA_real_, 2)
    ),
    list(
      label = "logit, zero or tiny standard error",
      estimate = tiny_estimate,
      standard_error = tiny_se,
      conf_type = "logit",
      warns = TRUE,
      conf.low = rep(NA_real_, 6),
      conf.high = rep(NA_real_, 6)
    ),
    list(
      label = "Wald, zero or tiny standard error",
      estimate = tiny_estimate,
      standard_error = tiny_se,
      conf_type = "wald",
      warns = FALSE,
      conf.low = tiny_estimate - z * tiny_se,
      conf.high = tiny_estimate + z * tiny_se
    )
  )
  for (case in cases) {
    warning <- NULL
    bounds <- withCallingHandlers(
      delta_interval_bounds(
        case$estimate,
        case$standard_error,
        0.95,
        case$conf_type
      ),
      warning = function(condition) {
        warning <<- conditionMessage(condition)
        invokeRestart("muffleWarning")
      }
    )
    if (case$warns) {
      expect_match(
        warning,
        "nonstructural boundary",
        fixed = TRUE,
        info = case$label
      )
      expect_identical(bounds$conf.low, case$conf.low, info = case$label)
      expect_identical(bounds$conf.high, case$conf.high, info = case$label)
    } else {
      expect_null(warning, info = case$label)
      expect_equal(bounds$conf.low, case$conf.low, info = case$label)
      expect_equal(bounds$conf.high, case$conf.high, info = case$label)
    }
  }
})

test_that("absorption during prediction does not certify numerical saturation", {
  skip_on_cran()
  data <- make_test_data(n_patients = 60, follow_up_time = 8, seed = 8118)
  model <- suppressWarnings(mostm(
    ordered(y) ~ tx,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    engine = "vglm"
  ))
  states <- seq_len(nrow(get_effective_coefs(model)) + 1L)
  baseline <- data[!is.na(data$y) & data$yprev != max(states), , drop = FALSE]
  baseline <- baseline[!duplicated(baseline$id), , drop = FALSE][1:4, ]
  individual <- sops(
    model,
    newdata = baseline,
    times = 1:3,
    absorb = max(states)
  )
  expect_gt(min(individual$estimate), 0)
  expect_lt(max(individual$estimate), 1)
  absorbing <- individual[individual$state == max(states), ]
  expect_equal(
    all(vapply(
      split(absorbing, absorbing$rowid),
      function(profile) {
        all(diff(profile$estimate[order(profile$time)]) > 0)
      },
      logical(1)
    )),
    TRUE
  )

  coefficients <- get_coef(model)
  coefficients[grepl("Intercept", names(coefficients))] <-
    coefficients[grepl("Intercept", names(coefficients))] + 1000
  saturated_model <- set_coef(model, coefficients)
  points <- list(
    sops(
      saturated_model,
      newdata = baseline,
      times = 1:3,
      absorb = max(states)
    ),
    avg_sops(
      saturated_model,
      newdata = baseline,
      times = 1:3,
      variables = list(tx = c(0, 1)),
      absorb = max(states)
    )
  )
  for (point in points) {
    expect_equal(all(point$estimate %in% c(0, 1)), TRUE)
    warning <- NULL
    result <- withCallingHandlers(
      inferences(point, method = "delta", vcov = stats::vcov(model)),
      warning = function(condition) {
        warning <<- conditionMessage(condition)
        invokeRestart("muffleWarning")
      }
    )
    expect_match(warning, "nonstructural boundary", fixed = TRUE)
    expect_equal(result$std.error, rep(0, nrow(result)))
    expect_identical(result$conf.low, rep(NA_real_, nrow(result)))
    expect_identical(result$conf.high, rep(NA_real_, nrow(result)))
  }
})
