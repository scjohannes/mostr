local_orm_datadist <- function(data, env = parent.frame()) {
  dd <- rms::datadist(data)
  old_dd_exists <- exists("dd", envir = globalenv(), inherits = FALSE)
  old_dd <- if (old_dd_exists) get("dd", envir = globalenv()) else NULL
  assign("dd", dd, envir = globalenv())
  withr::local_options(datadist = "dd", .local_envir = env)
  withr::defer(
    {
      if (old_dd_exists) {
        assign("dd", old_dd, envir = globalenv())
      } else if (exists("dd", envir = globalenv(), inherits = FALSE)) {
        remove(list = "dd", envir = globalenv())
      }
    },
    envir = env
  )
}

# The wrapper arguments below are expressions over variables that exist only
# in the fitting function's frame, so a refit must not re-evaluate them.
fit_vglm_in_closure <- function(data, parallel_flag, first_time) {
  vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = parallel_flag),
    data = data,
    id_var = "id",
    first_followup_time = if (first_time > 0) first_time else 1
  )
}

fit_orm_in_closure <- function(data, use_logistic, first_time) {
  orm_markov(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    family = if (use_logistic) "logistic" else "probit",
    first_followup_time = if (first_time > 0) first_time else 1
  )
}

expect_refits_succeed <- function(fit, method) {
  out <- avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    y_levels = 1:6,
    absorb = 6
  )
  warnings <- character()
  inferred <- withCallingHandlers(
    inferences(out, method = method, n_draws = 2, use_coefstart = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("Bootstrap model fitting failed", warnings)))
  expect_equal(attr(inferred, "n_successful"), 2L)
}

test_that("vglm_markov bootstrap refits use the fit-time argument values", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5101)
  fit <- fit_vglm_in_closure(data, parallel_flag = TRUE, first_time = 1)

  withr::local_seed(5102)
  expect_refits_succeed(fit, "bootstrap")
  withr::local_seed(5103)
  expect_refits_succeed(fit, "fwb")
})

test_that("orm_markov bootstrap refits use the fit-time argument values", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5201)
  local_orm_datadist(data)
  fit <- fit_orm_in_closure(data, use_logistic = TRUE, first_time = 1)

  withr::local_seed(5202)
  expect_refits_succeed(fit, "bootstrap")
  withr::local_seed(5203)
  expect_refits_succeed(fit, "fwb")
})

test_that("refits ignore later changes to variables used by the original call", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5301)
  parallel_flag <- TRUE
  fit <- vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = parallel_flag),
    data = data,
    id_var = "id"
  )
  # A non-parallel refit would produce a different number of coefficients.
  parallel_flag <- FALSE

  boot <- mostr:::bootstrap_analysis_wrapper(
    boot_data = data,
    model = fit,
    factor_cols = character(0),
    original_data = data,
    update_datadist = FALSE
  )

  expect_equal(length(stats::coef(boot$model)), length(stats::coef(fit)))
})

# Row weights supplied as a vector that is not a column of `data` must follow
# the resampled rows. The vector below lives only in the fitting function's
# frame, and the data carry a copy (`w_copy`) so the tests can check which
# weight each resampled row should receive.
fit_weighted_in_closure <- function(data, backend) {
  row_weights <- data$w_copy
  switch(
    backend,
    vglm = vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      weights = row_weights,
      id_var = "id"
    ),
    orm = orm_markov(
      ordered(y) ~ time + tx + yprev,
      data = data,
      weights = row_weights,
      id_var = "id"
    )
  )
}

make_weighted_refit_data <- function(seed) {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = seed)
  withr::with_seed(seed, data$w_copy <- stats::runif(nrow(data), 0.5, 1.5))
  data
}

resample_patients <- function(data, ids) {
  boot_ids <- data.frame(
    original_id = ids,
    new_id = paste0(ids, "_", seq_along(ids)),
    boot_id = 1L
  )
  boot_data <- mostr:::materialize_bootstrap_sample(boot_ids, data, "id")
  boot_data$id <- boot_data$new_id
  boot_data
}

expect_refit_weights_follow_rows <- function(backend) {
  data <- make_weighted_refit_data(5501)
  if (backend == "orm") {
    skip_if_not_installed("rms")
    local_orm_datadist(data, env = parent.frame())
  }
  fit <- suppressWarnings(fit_weighted_in_closure(data, backend))
  ids <- unique(data$id)
  # Bootstrap inference resamples the refit data stored on the model, which
  # carries each row's weight. Fewer, repeated patients: the resampled data
  # has a different row count.
  boot_data <- resample_patients(
    mostr:::markov_model_refit_data(fit),
    c(ids[1:25], ids[1:5])
  )

  boot <- suppressWarnings(mostr:::bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = fit,
    factor_cols = character(0),
    original_data = data,
    update_datadist = FALSE
  ))

  expect_false(is.null(boot$model))
  expect_equal(
    mostr:::bootstrap_stored_fit_weights(boot$model, boot$data),
    boot$data$w_copy
  )

  # FWB multiplies the original row weights by the random weights.
  fwb <- suppressWarnings(mostr:::bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = fit,
    factor_cols = character(0),
    original_data = data,
    update_datadist = FALSE,
    fit_weights = rep(2, nrow(boot_data))
  ))
  expect_equal(fwb$data$.mostr_fit_weight, 2 * boot_data$w_copy)
}

test_that("vglm_markov refits resample weights given as an outside vector", {
  expect_refit_weights_follow_rows("vglm")
})

test_that("orm_markov refits resample weights given as an outside vector", {
  expect_refit_weights_follow_rows("orm")
})

test_that("bootstrap and FWB succeed with weights given as an outside vector", {
  data <- make_weighted_refit_data(5601)
  fit <- suppressWarnings(fit_weighted_in_closure(data, "vglm"))

  withr::local_seed(5602)
  expect_refits_succeed(fit, "bootstrap")
  withr::local_seed(5603)
  expect_refits_succeed(fit, "fwb")
})

test_that("FWB never falls back to unweighted refits for wrapper fits", {
  data <- make_weighted_refit_data(5701)
  fit <- suppressWarnings(fit_weighted_in_closure(data, "vglm"))
  out <- avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    y_levels = 1:6,
    absorb = 6
  )

  withr::local_seed(5702)
  warnings <- character()
  withCallingHandlers(
    inferences(out, method = "fwb", n_draws = 2, use_coefstart = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("Could not recover original model weights", warnings)))
})

test_that("stored row weights stay out of starting profiles and SOP output", {
  data <- make_weighted_refit_data(5801)
  fit <- suppressWarnings(fit_weighted_in_closure(data, "vglm"))
  weight_column <- mostr:::markov_prior_weight_column()

  expect_equal(
    mostr:::markov_model_refit_data(fit)[[weight_column]],
    data$w_copy
  )
  expect_false(
    weight_column %in% names(attr(fit, "markov_starting_profile_data"))
  )
  out <- sops(fit, times = 1:2, y_levels = 1:6, absorb = 6)
  expect_false(weight_column %in% names(out))
})

# A `subset` that refers to a variable local to the fitting function. The
# stored refit data already exclude the removed rows, so refits must not
# evaluate the expression again.
fit_subset_in_closure <- function(data, backend) {
  keep_rows <- data$id != data$id[[1]]
  switch(
    backend,
    vglm = vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      subset = keep_rows,
      id_var = "id"
    ),
    orm = orm_markov(
      ordered(y) ~ time + tx + yprev,
      data = data,
      subset = keep_rows,
      id_var = "id"
    )
  )
}

test_that("refits of wrapper fits with an outside subset vector succeed", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5901)
  local_orm_datadist(data)
  for (backend in c("vglm", "orm")) {
    fit <- fit_subset_in_closure(data, backend)
    refit_data <- mostr:::markov_model_refit_data(fit)
    expect_false(data$id[[1]] %in% refit_data$id)
    expect_equal(nrow(refit_data), sum(data$id != data$id[[1]]))

    withr::local_seed(5902)
    expect_refits_succeed(fit, "bootstrap")
    withr::local_seed(5903)
    expect_refits_succeed(fit, "fwb")
  }
})

test_that("refits reuse a fit-time coefstart given as a local variable", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5801)
  fit_with_local_coefstart <- function(data) {
    pilot <- vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id"
    )
    start_values <- stats::coef(pilot)
    vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      coefstart = start_values
    )
  }
  fit <- fit_with_local_coefstart(data)

  # Without `use_coefstart`, the refit keeps the user's own starting values.
  boot <- mostr:::bootstrap_analysis_wrapper(
    boot_data = mostr:::markov_model_refit_data(fit),
    model = fit,
    factor_cols = character(0),
    original_data = data,
    update_datadist = FALSE
  )
  expect_false(is.null(boot$model))

  withr::local_seed(5802)
  expect_refits_succeed(fit, "bootstrap")
})

fit_vglm_with_local_coefstart <- function(data) {
  pilot <- vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id"
  )
  start_values <- stats::coef(pilot)
  vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    coefstart = start_values
  )
}

test_that("refits on samples with every state reuse the fit-time coefstart", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5901)
  fit <- fit_vglm_with_local_coefstart(data)
  refit_data <- mostr:::markov_model_refit_data(fit)

  boot <- mostr:::bootstrap_analysis_wrapper(
    boot_data = refit_data,
    model = fit,
    factor_cols = intersect(c("y", "yprev"), names(refit_data)),
    original_data = refit_data,
    y_levels = 1:6,
    absorb = 6,
    update_datadist = FALSE
  )

  expect_length(boot$missing_states, 0L)
  expect_equal(
    methods::slot(boot$model$vglm_fit, "call")$coefstart,
    mostr:::markov_model_refit_args(fit)$coefstart
  )
})

test_that("refits on samples missing a state are fit without coefstart", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5902)
  fit <- fit_vglm_with_local_coefstart(data)
  refit_data <- mostr:::markov_model_refit_data(fit)
  absorbed_ids <- unique(refit_data$id[refit_data$y == 6])
  expect_gt(length(absorbed_ids), 0L)
  kept_ids <- setdiff(unique(refit_data$id), absorbed_ids)
  # No sampled patient reaches the absorbing state, so the refit has fewer
  # intercepts than the fit-time starting coefficients.
  boot_data <- resample_patients(refit_data, kept_ids)

  for (use_coefstart in c(FALSE, TRUE)) {
    # The sample deliberately lacks the absorbing state and has few patients,
    # which triggers expected warnings about both.
    boot <- suppressWarnings(mostr:::bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = fit,
      factor_cols = intersect(c("y", "yprev"), names(boot_data)),
      original_data = refit_data,
      y_levels = 1:6,
      absorb = 6,
      update_datadist = FALSE,
      use_coefstart = use_coefstart
    ))

    expect_true(6 %in% boot$missing_states)
    expect_false(is.null(boot$model))
    expect_null(methods::slot(boot$model$vglm_fit, "call")$coefstart)
  }
})
