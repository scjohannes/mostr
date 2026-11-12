test_that("soprobMarkovOrdm and soprob_markov yield identical results - single patient", {
  skip_on_cran()
  skip_if_not_installed("Hmisc")
  withr::local_seed(1234)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 100,
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()
  time_covariates <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  baseline_data <- data[data$time == 1, ]

  absorb <- 6
  y_levels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up

  fit.b <- VGAM::vglm(
    y ~ rms::rcs(time, 3) *
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  fit.a <- suppressWarnings(mostm(
    y ~
      (time_lin +
        time_nlin_1) *
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    engine = "vglm"
  ))$vglm_fit

  # Test single patient (first patient)
  a <- Hmisc::soprobMarkovOrdm(
    fit.b,
    baseline_data[1, ],
    times = times,
    ylevels = y_levels,
    absorb = absorb
  )

  b <- soprob_markov(
    model = fit.a,
    newdata = baseline_data[1, ],
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    p_var = "yprev",
    time_covariates = time_covariates
  )[1, , ]

  expect_equal(a, b, tolerance = 1e-15)
})

test_that("soprobMarkovOrdm and soprob_markov yield identical results - single patient; partial PO", {
  skip_on_cran()
  skip_if_not_installed("Hmisc")
  withr::local_seed(1234)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 100,
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  time_covariates <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  baseline_data <- data[data$time == 1, ]

  absorb <- 6
  y_levels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up

  fit.b <- VGAM::vglm(
    y ~ rms::rcs(time, 3) +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data
  )

  fit.a <- suppressWarnings(mostm(
    y ~
      time_lin +
      time_nlin_1 +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data,
    id_var = "id",
    engine = "vglm"
  ))$vglm_fit

  # Test single patient (first patient)
  a <- Hmisc::soprobMarkovOrdm(
    fit.b,
    baseline_data[1, ],
    times = times,
    ylevels = y_levels,
    absorb = absorb
  )

  b <- soprob_markov(
    model = fit.a,
    newdata = baseline_data[1, ],
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    p_var = "yprev",
    time_covariates = time_covariates
  )[1, , ]

  expect_equal(a, b, tolerance = 1e-15)
})

test_that("soprobMarkovOrdm and soprob_markov yield identical results - all patients", {
  skip_on_cran()
  skip_if_not_installed("Hmisc")
  withr::local_seed(7890)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 50, # Use smaller sample for full comparison
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  time_covariates <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  baseline_data <- data[data$time == 1, ]
  y_levels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up
  absorb <- 6

  fit.b <- VGAM::vglm(
    y ~ rms::rcs(time, 3) +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data
  )

  fit.a <- suppressWarnings(mostm(
    y ~
      time_lin +
      time_nlin_1 +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data,
    id_var = "id",
    engine = "vglm"
  ))$vglm_fit

  # Test ALL patients
  # Get results from soprobMarkovOrdm (one patient at a time)
  results_hmisc <- array(
    dim = c(nrow(baseline_data), length(times), length(y_levels)),
    dimnames = list(NULL, NULL, levels(y_levels))
  )

  for (i in seq_len(nrow(baseline_data))) {
    results_hmisc[i, , ] <- Hmisc::soprobMarkovOrdm(
      fit.b,
      baseline_data[i, ],
      times = times,
      ylevels = y_levels,
      absorb = absorb
    )
  }

  # Get results from soprob_markov (all at once)
  results_new <- soprob_markov(
    model = fit.a,
    newdata = baseline_data,
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    p_var = "yprev",
    time_covariates = time_covariates
  )

  # Strip dimnames for fair comparison (soprob_markov adds them, manual loop doesn't)
  dimnames(results_new)[1:2] <- list(NULL, NULL)

  expect_equal(results_hmisc, results_new, tolerance = 1e-15)
})


test_that("compiled plan matches markov_msm_build() and the R reference for PPO, with and without constraints", {
  skip_on_cran()
  follow_up <- 30
  times <- 1:follow_up
  absorb <- 6

  # User constraint matrices: a linear trend in the time effect across
  # thresholds and a separate treatment effect at the last threshold.
  constraints <- list(
    "(Intercept)" = diag(5),
    "time_lin" = cbind(time_lin_1 = 1, time_lin_5 = 0:4),
    "time_nlin_1" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = rep(1, 5), Deviation_5th = c(rep(0, 4), 1)),
    "yprev" = cbind(PO_effect = 1),
    "age" = cbind(PO_effect = 1),
    "time_lin:tx" = cbind(PO_effect = 1),
    "time_nlin_1:tx" = cbind(PO_effect = 1),
    "time_lin:yprev" = cbind(PO_effect = 1),
    "tx:age" = cbind(PO_effect = 1)
  )
  cases <- list(
    list(
      label = "constrained PPO",
      n_patients = 100,
      age_seed = 12345,
      constraints = constraints,
      manual_tolerance = 1e-15
    ),
    list(
      label = "PPO",
      n_patients = 200,
      age_seed = 1234567,
      constraints = NULL,
      manual_tolerance = 1e-10
    )
  )

  for (case in cases) {
    data <- withr::with_seed(case$age_seed, {
      make_test_data(
        n_patients = case$n_patients,
        follow_up_time = follow_up,
        seed = 12345
      ) |>
        dplyr::group_by(id) |>
        dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
        dplyr::ungroup()
    })
    time_covariates <- data |>
      dplyr::select(time_lin, time_nlin_1) |>
      dplyr::distinct() |>
      as.data.frame()
    y_levels <- factor(sort(unique(data$y)), ordered = FALSE)

    # `constraints = NULL` keeps VGAM's constraints from `parallel`.
    fit <- suppressWarnings(mostm(
      y ~ time_lin +
        time_nlin_1 +
        tx +
        yprev +
        age +
        time_lin:tx +
        time_nlin_1:tx +
        yprev:time_lin +
        tx:age,
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx + time_lin),
      data = data,
      id_var = "id",
      constraints = case$constraints,
      engine = "vglm"
    ))$vglm_fit

    newdata <- data[data$time == 1, ] # Baseline rows

    # Both paths run the same native recursion in markov_msm_run(); they differ
    # only in how the design matrices and effective coefficients are built.
    # --- COMPILED PLAN (soprob_markov) ---
    res_plan <- soprob_markov(
      model = fit,
      newdata = newdata,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_covariates = time_covariates
    )

    # --- MANUAL BUILD (markov_msm_build + markov_msm_run) ---
    components <- markov_msm_build(
      model = fit,
      newdata = newdata,
      settings = markov_settings(
        y_levels = y_levels,
        time_covariates = time_covariates
      )
    )
    Gamma_hat <- compute_Gamma(get_coef(fit), VGAM::constraints(fit))
    res_manual <- markov_msm_run(
      components = components,
      Gamma = Gamma_hat,
      times = times,
      absorb = absorb
    )

    # --- R REFERENCE (pure-R recursion, coefficients read from VGAM) ---
    plan <- compile_sop_execution_plan(
      model = fit,
      newdata = newdata,
      times = times,
      settings = markov_settings(
        y_levels = y_levels,
        absorb = absorb,
        time_var = "time_lin",
        time_covariates = time_covariates
      )
    )
    res_reference <- reference_first_order_plan_sops(
      plan,
      reference_backend_gamma(fit, plan$components$col_names),
      category_probabilities = reference_checked_probabilities
    )

    # Fix dimnames for comparison
    dimnames(res_plan)[1:2] <- list(NULL)

    expect_equal(
      res_plan,
      res_manual,
      tolerance = case$manual_tolerance,
      info = case$label
    )
    expect_equal(
      unname(res_plan),
      res_reference,
      tolerance = 1e-12,
      info = case$label
    )
  }
})
