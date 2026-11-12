make_delta_average_fixture <- function() {
  out <- expand.grid(
    time = 1:3,
    state = c("1", "2"),
    tx = c(0, 1, 2),
    KEEP.OUT.ATTRS = FALSE
  )
  state_one <- c(
    "0" = 0.2,
    "1" = 0.5,
    "2" = 0.35
  )
  out$estimate <- state_one[as.character(out$tx)] + 0.05 * (out$time - 1)
  out$estimate[out$state == "2"] <- 1 - out$estimate[out$state == "2"]
  class(out) <- c("markov_avg_sops", "data.frame")
  attr(out, "avg_args") <- list(
    variables = list(tx = c(0, 1, 2)),
    by = NULL,
    times = 1:3,
    id_var = "id"
  )
  attr(out, "settings") <- markov_settings(id_var = "id", y_levels = c("1", "2"))
  attr(out, "newdata_orig") <- data.frame(
    id = 1:4,
    tx = c(0, 0, 1, 1),
    yprev = c("1", "1", "2", "2")
  )
  attr(out, "newdata_supplied") <- TRUE
  out
}

test_that("get_jacobian() stays unexported", {
  expect_equal("get_jacobian" %in% getNamespaceExports("mostr"), FALSE)
  expect_type(getFromNamespace("get_jacobian", "mostr"), "closure")
})

test_that("vcov selection preserves established positional arguments", {
  expect_identical(
    names(formals(inferences)),
    c(
      "x",
      "method",
      "n_draws",
      "vcov",
      "cluster",
      "workers",
      "seed",
      "conf_level",
      "conf_type",
      "null",
      "return_draws",
      "use_coefstart"
    )
  )
})

delta_public_vglm_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 60,
        follow_up_time = 8,
        seed = 8101
      )
      model <- make_test_model(data)
      baseline <- data[!duplicated(data$id), , drop = FALSE][
        1:12,
        ,
        drop = FALSE
      ]
      value <<- list(
        model = model,
        baseline = baseline,
        covariance = stats::vcov(model),
        y_levels = seq_len(nrow(get_effective_coefs(model)) + 1L)
      )
    }
    value
  }
})

delta_public_orm_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 55,
        follow_up_time = 7,
        seed = 8102
      )
      model <- suppressWarnings(mostm(
        y ~ time + tx + yprev,
        data = data,
        id_var = "id"
      ))
      baseline <- data[!duplicated(data$id), , drop = FALSE][
        1:8,
        ,
        drop = FALSE
      ]
      covariance <- orm_model_bread(model)$bread
      value <<- list(
        model = model,
        baseline = baseline,
        covariance = covariance,
        y_levels = model$yunique
      )
    }
    value
  }
})

delta_public_unconditional_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 32,
        follow_up_time = 6,
        seed = 8103
      )
      model <- mostm(
        ordered(y) ~ time_lin + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        engine = "vglm"
      )
      value <<- list(
        data = data,
        model = model,
        y_levels = seq_len(nrow(get_effective_coefs(model$vglm_fit)) + 1L)
      )
    }
    value
  }
})

replay_public_sops <- function(object, model) {
  args <- attr(object, "call_args")
  settings <- markov_result_settings(object)
  sops(
    model = model,
    newdata = attr(object, "newdata_pred"),
    times = args$times,
    absorb = settings$absorb,
    gap_var = settings$gap_var,
    time_covariates = settings$time_covariates
  )
}

central_public_sops_jacobian <- function(object, model) {
  object_key <- sop_draw_cell_key(object, c("rowid", "time", "state"))
  numeric_jacobian(model, function(model) {
    result <- replay_public_sops(object, model)
    index <- match(
      object_key,
      sop_draw_cell_key(result, c("rowid", "time", "state"))
    )
    if (anyNA(index)) {
      stop("Finite-difference SOP replay did not return every public cell.")
    }
    result$estimate[index]
  })
}

test_that("analytical comparison operators reproduce existing reductions", {
  avg <- make_delta_average_fixture()
  avg_args <- attr(avg, "avg_args")
  state_sets <- list(low = "1", all = c("1", "2"))

  sop <- avg_comparison_from_avg_sops(
    avg,
    estimand = "sop",
    state_sets = state_sets,
    comparison = "difference",
    time_map = NULL,
    baseline_time = 0,
    target_times = NULL,
    time_unit = NULL,
    return_draws = FALSE
  )
  sop_args <- list(
    estimand = "sop",
    state_sets = state_sets,
    comparison = "difference"
  )
  sop_operator <- delta_comparison_operator(sop, avg, sop_args, avg_args)
  expect_equal(
    drop(delta_apply_comparison_operator(sop_operator, avg$estimate)),
    sop$estimate
  )

  real_time <- avg_comparison_from_avg_sops(
    avg,
    estimand = "time_in_state",
    state_sets = state_sets,
    comparison = "difference",
    time_map = c("1" = 3, "2" = 7, "3" = 14),
    baseline_time = 0,
    target_times = 0:14,
    time_unit = "days",
    return_draws = FALSE
  )
  real_time_args <- list(
    estimand = "time_in_state",
    state_sets = state_sets,
    comparison = "difference",
    time_map = c("1" = 3, "2" = 7, "3" = 14),
    baseline_time = 0,
    target_times = 0:14
  )
  real_time_operator <- delta_comparison_operator(
    real_time,
    avg,
    real_time_args,
    avg_args
  )
  expect_equal(
    drop(delta_apply_comparison_operator(real_time_operator, avg$estimate)),
    real_time$estimate,
    tolerance = 1e-12
  )
})

test_that("low-rank analytical accessors select rows without dense storage", {
  object <- data.frame(
    time = 1:3,
    state = "1",
    estimate = c(0.2, 0.4, 0.6)
  )
  class(object) <- c("markov_avg_sops", "data.frame")
  jacobian <- matrix(
    c(1, 0, 1, 1, 0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("beta_1", "beta_2"))
  )
  covariance <- matrix(
    c(0.04, 0.01, 0.01, 0.09),
    nrow = 2,
    dimnames = list(colnames(jacobian), colnames(jacobian))
  )
  variance <- delta_jacobian_variance(jacobian, covariance)
  result <- delta_finalize_result(
    object,
    standard_error = sqrt(variance),
    conf_level = 0.95,
    conf_type = "logit",
    target = "empirical",
    analytical = list(
      representation = "coefficient",
      jacobian = jacobian,
      coefficient_vcov = covariance,
      coefficient_names = colnames(jacobian),
      covariance_metadata = list(source = "test")
    )
  )

  selected <- get_jacobian(result, rows = c(3, 1))
  expect_equal(unname(selected), unname(jacobian[c(3, 1), , drop = FALSE]))
  expected <- jacobian[c(3, 1), , drop = FALSE] %*%
    covariance %*%
    t(jacobian[c(3, 1), , drop = FALSE])
  expect_equal(
    unname(vcov.markov_avg_sops(result, rows = c(3, 1))),
    unname(expected)
  )
  expect_equal(
    unname(stats::vcov(result, rows = c(3, 1))),
    unname(expected)
  )
  expect_null(attr(result, "draws"))
  expect_identical(attr(result, "method"), "delta")
  expect_identical(attr(result, "conf_type"), "logit")
})

test_that("delta replay rejects stale estimates even when cell keys align", {
  skip_on_cran()
  case <- delta_public_vglm_case()
  point <- avg_sops(
    case$model,
    newdata = case$baseline,
    variables = list(tx = c(0, 1)),
    times = 1,
    absorb = max(case$y_levels)
  )
  point$estimate[1] <- point$estimate[1] + 1e-8

  condition <- tryCatch(
    inferences(point, method = "delta", vcov = case$covariance),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(
    conditionMessage(condition),
    "did not reproduce the stored point estimates",
    fixed = TRUE
  )
})

test_that("public delta workflow returns bounded SOP and Wald comparison intervals", {
  skip_on_cran()
  case <- delta_public_vglm_case()
  avg <- avg_sops(
    case$model,
    newdata = case$baseline,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = max(case$y_levels)
  )
  inferred <- inferences(
    avg,
    method = "delta",
    vcov = case$covariance
  )

  expect_identical(attr(inferred, "method"), "delta")
  expect_identical(attr(inferred, "target"), "empirical")
  expect_identical(attr(inferred, "conf_type"), "logit")
  expect_equal(nrow(get_jacobian(inferred)), nrow(inferred))
  expect_equal(
    unname(diag(vcov.markov_avg_sops(inferred, rows = 1:4))),
    inferred$std.error[1:4]^2,
    tolerance = 1e-12
  )
  expect_gte(min(inferred$conf.low), 0)
  expect_lte(max(inferred$conf.high), 1)

  state_rows <- which(
    inferred$time == 1 & as.character(inferred$tx) == "0"
  )
  expect_length(state_rows, length(case$y_levels))
  state_covariance <- stats::vcov(inferred, rows = state_rows)
  expect_identical(state_covariance, t(state_covariance))
  eigenvalues <- eigen(
    state_covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  eigen_tolerance <-
    sqrt(.Machine$double.eps) *
    max(abs(eigenvalues), .Machine$double.eps) *
    length(eigenvalues)
  expect_gte(min(eigenvalues), -eigen_tolerance)
  expect_lt(
    max(abs(state_covariance %*% rep(1, length(state_rows)))),
    1e-10
  )
  expect_lte(
    sum(abs(eigenvalues) > eigen_tolerance),
    length(state_rows) - 1L
  )

  comparison <- avg_comparisons(
    case$model,
    newdata = case$baseline,
    variables = list(tx = c(0, 1)),
    estimand = "sop",
    state_sets = list(low = c(1, 2)),
    times = 1:3,
    absorb = max(case$y_levels)
  )
  comparison <- inferences(
    comparison,
    method = "delta",
    vcov = case$covariance
  )

  expect_identical(attr(comparison, "conf_type"), "wald")
  expect_equal(nrow(get_jacobian(comparison)), nrow(comparison))
  expect_equal(
    comparison$conf.low,
    comparison$estimate - stats::qnorm(0.975) * comparison$std.error
  )
  expect_equal(
    comparison$conf.high,
    comparison$estimate + stats::qnorm(0.975) * comparison$std.error
  )
})

test_that("unconditional targets reject explicit coefficient covariance", {
  object <- data.frame(
    estimand = "sop",
    time = 1,
    state_set = "1",
    comparison = "difference",
    estimate = 0
  )
  class(object) <- c("markov_avg_comparisons", "data.frame")
  attr(object, "comparison_args") <- list(
    estimand = "sop",
    comparison = "difference"
  )
  attr(object, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = NULL
  )

  condition <- tryCatch(
    delta_validate_comparison_scope(
      object,
      target = "unconditional",
      vcov = diag(1),
      conf_type = "wald"
    ),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "cannot be supplied", fixed = TRUE)
})

test_that("public fixed VGLM Jacobians match coefficient replay", {
  case <- delta_public_vglm_case()
  point <- sops(
    case$model,
    newdata = case$baseline[1:4, , drop = FALSE],
    times = 1:3,
    absorb = max(case$y_levels)
  )

  set.seed(8111)
  rng_state <- .Random.seed
  inferred <- inferences(
    point,
    method = "delta",
    vcov = case$covariance,
    seed = 9191
  )
  numeric_jacobian <- central_public_sops_jacobian(point, case$model)
  analytical_jacobian <- get_jacobian(inferred)

  expect_identical(.Random.seed, rng_state)
  expect_lt(
    max(abs(unname(analytical_jacobian) - unname(numeric_jacobian))),
    1e-7
  )
  expected_vcov <- analytical_jacobian %*%
    case$covariance[
      colnames(analytical_jacobian),
      colnames(analytical_jacobian)
    ] %*%
    t(analytical_jacobian)
  expect_equal(
    unname(stats::vcov(inferred, rows = 1:5)),
    unname(expected_vcov[1:5, 1:5]),
    tolerance = 1e-12
  )
  expect_identical(attr(inferred, "target"), "fixed")
  expect_null(attr(inferred, "draws"))
  expect_null(attr(inferred, "n_draws"))
})

test_that("public fixed ORM Jacobians match coefficient replay", {
  skip_on_cran()
  case <- delta_public_orm_case()
  point <- sops(
    case$model,
    newdata = case$baseline[1:3, , drop = FALSE],
    times = 1:3,
    absorb = max(case$y_levels)
  )
  inferred <- inferences(
    point,
    method = "delta",
    vcov = case$covariance
  )

  numeric_jacobian <- central_public_sops_jacobian(point, case$model)
  expect_lt(
    max(abs(unname(get_jacobian(inferred)) - unname(numeric_jacobian))),
    1e-7
  )
})

test_that("public empirical ORM averages use complete named covariance", {
  case <- delta_public_orm_case()
  avg <- avg_sops(
    case$model,
    newdata = case$baseline,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = max(case$y_levels)
  )
  inferred <- inferences(
    avg,
    method = "delta",
    vcov = case$covariance
  )

  jacobian <- get_jacobian(inferred)
  aligned_vcov <- case$covariance[colnames(jacobian), colnames(jacobian)]
  expected_variance <- rowSums((jacobian %*% aligned_vcov) * jacobian)
  expect_equal(
    unname(inferred$std.error^2),
    unname(expected_variance),
    tolerance = 1e-12
  )
  expect_identical(attr(inferred, "target"), "empirical")
  expect_identical(attr(inferred, "covariance_source"), "explicit")
})

test_that("fitted-cohort unconditional averages expose stacked influence", {
  skip_on_cran()
  case <- delta_public_unconditional_case()
  avg <- avg_sops(
    case$model,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    absorb = max(case$y_levels)
  )
  inferred <- inferences(avg, method = "delta", vcov = "unconditional")

  avg_args <- attr(avg, "avg_args")
  newdata <- attr(avg, "newdata_pred")
  model_plan <- case$model$vglm_fit
  plan <- delta_compile_plan(
    object = avg,
    model_plan = model_plan,
    newdata = newdata,
    times = avg_args$times,
    output = "average"
  )
  engine <- delta_validate_plan_result(
    run_sop_delta_plan(plan, model_plan),
    n_profiles = nrow(newdata),
    n_times = length(avg_args$times),
    n_states = length(case$y_levels)
  )
  grid <- do.call(expand.grid, avg_args$variables)
  n_each <- nrow(newdata) / nrow(grid)
  reduced <- delta_average_arrays(
    probabilities = engine$probabilities,
    jacobian = engine$jacobian,
    grid = grid,
    times = avg_args$times,
    y_levels = case$y_levels,
    variables = avg_args$variables,
    n_each = n_each
  )
  index <- delta_match_cells(
    avg,
    reduced$cells,
    c("time", "state", names(avg_args$variables))
  )
  average_jacobian <- reduced$jacobian[index, , drop = FALSE]
  colnames(average_jacobian) <- engine$coefficient_names
  individual_values <- reduced$individual_values[, index, drop = FALSE]
  cell_names <- paste0("cell_", seq_len(ncol(individual_values)))
  colnames(individual_values) <- cell_names
  rownames(average_jacobian) <- cell_names
  baseline <- newdata[seq_len(n_each), , drop = FALSE]
  profile_ids <- delta_profile_ids(baseline, markov_result_settings(avg)$id_var)
  direct <- delta_stacked_influence(
    individual_values = individual_values,
    average_jacobian = average_jacobian,
    profile_ids = profile_ids,
    score_components = get_delta_score_components(
      case$model,
      cluster = NULL
    )
  )
  analytical <- attr(inferred, "analytical")

  expect_equal(
    unname(analytical$influence),
    unname(direct$influence),
    tolerance = 1e-12
  )
  expect_equal(
    unname(get_jacobian(inferred)),
    unname(average_jacobian),
    tolerance = 1e-12
  )
  expect_equal(
    unname(inferred$std.error),
    unname(direct$std.error),
    tolerance = 1e-12
  )
  expect_equal(
    stats::vcov(inferred, rows = 1:5),
    direct$covariance[1:5, 1:5],
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_identical(analytical$profile_ids, profile_ids)
  expect_identical(attr(inferred, "target"), "unconditional")
  metadata <- attr(inferred, "covariance_metadata")
  expect_equal(attr(inferred, "covariance_source"), "stacked_patient_influence")
  expect_equal(metadata$n_patients, length(profile_ids))
  expect_equal(metadata$finite_sample_method, "patient_level_sample_covariance")
  expect_equal(
    metadata$finite_sample_factor,
    length(profile_ids) / (length(profile_ids) - 1L)
  )
  expect_equal(metadata$bread_source, "robcov_vglm_observed")
})

test_that("backend HC settings affect empirical but not unconditional inference", {
  skip_on_cran()
  case <- delta_public_unconditional_case()
  raw <- case$model$vglm_fit
  data <- case$data
  hc0 <- robcov_vglm(
    raw,
    cluster = attr(raw, "markov_data")$id,
    type = "HC0",
    cadjust = FALSE
  )
  hc1 <- robcov_vglm(
    raw,
    cluster = attr(raw, "markov_data")$id,
    type = "HC1",
    cadjust = TRUE
  )

  point_hc0 <- avg_sops(
    hc0,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    absorb = max(case$y_levels)
  )
  point_hc1 <- avg_sops(
    hc1,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    absorb = max(case$y_levels)
  )
  empirical_hc0 <- inferences(point_hc0, method = "delta", vcov = "conditional")
  empirical_hc1 <- inferences(point_hc1, method = "delta", vcov = "conditional")
  super_hc0 <- inferences(
    point_hc0,
    method = "delta",
    vcov = "unconditional"
  )
  super_hc1 <- inferences(
    point_hc1,
    method = "delta",
    vcov = "unconditional"
  )

  expect_gt(max(abs(empirical_hc0$std.error - empirical_hc1$std.error)), 1e-10)
  expect_equal(super_hc0$estimate, super_hc1$estimate, tolerance = 1e-12)
  expect_equal(super_hc0$std.error, super_hc1$std.error, tolerance = 1e-12)
  expect_equal(
    attr(super_hc0, "covariance_metadata")$backend_hc_type_ignored,
    "HC0"
  )
  expect_equal(
    attr(super_hc1, "covariance_metadata")$backend_hc_type_ignored,
    "HC1"
  )
  expect_identical(
    attr(super_hc0, "covariance_metadata")$backend_cadjust_ignored,
    FALSE
  )
  expect_identical(
    attr(super_hc1, "covariance_metadata")$backend_cadjust_ignored,
    TRUE
  )
})

test_that("public delta scope enforces fixed targets and patient clustering", {
  skip_on_cran()
  local_reproducible_output(width = 80)
  case <- delta_public_vglm_case()
  fixed <- sops(
    case$model,
    newdata = case$baseline[1:3, , drop = FALSE],
    times = 1:2,
    absorb = max(case$y_levels)
  )

  expect_snapshot(
    inferences(
      fixed,
      method = "delta",
      vcov = "unconditional"
    ),
    error = TRUE
  )
  # Conditional variance uses the patient-clustered covariance that
  # `mostm()` stores on the fit.
  conditional <- inferences(fixed, method = "delta", vcov = "conditional")
  explicit <- inferences(
    fixed,
    method = "delta",
    vcov = robcov_vglm(
      case$model$vglm_fit,
      cluster = attr(case$model, "markov_data")$id
    )$var
  )
  expect_true(all(is.finite(conditional$std.error)))
  expect_equal(conditional$std.error, explicit$std.error, tolerance = 1e-10)

  unconditional_case <- delta_public_unconditional_case()
  avg <- avg_sops(
    unconditional_case$model,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    absorb = max(unconditional_case$y_levels)
  )
  expect_snapshot(
    inferences(
      avg,
      method = "delta",
      vcov = c("conditional", "unconditional")
    ),
    error = TRUE
  )
})

test_that("delta vcov choices preserve calculations and select class defaults", {
  skip_on_cran()
  case <- delta_public_unconditional_case()
  args <- list(
    model = case$model,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    absorb = max(case$y_levels)
  )
  objects <- list(
    do.call(sops, args),
    do.call(avg_sops, args),
    do.call(avg_comparisons, c(args, list(estimand = "sop", state_sets = "1")))
  )
  for (object in objects) {
    individual <- inherits(object, "markov_sops")
    conditional <- inferences(object, method = "delta", vcov = "conditional")
    matrix_result <- inferences(object, method = "delta", vcov = case$model$var)
    expect_equal(
      conditional$std.error,
      matrix_result$std.error,
      tolerance = 1e-12
    )
    expect_equal(
      conditional$conf.low,
      matrix_result$conf.low,
      tolerance = 1e-12
    )
    expect_identical(
      attr(conditional, "target"),
      if (individual) "fixed" else "empirical"
    )
    default <- inferences(object, method = "delta")
    explicit <- if (individual) {
      conditional
    } else {
      inferences(object, method = "delta", vcov = "unconditional")
    }
    expect_equal(default, explicit)
    expect_equal(inferences(object, method = "delta", vcov = NULL), explicit)
    expect_identical(
      attr(default, "target"),
      if (individual) "fixed" else "unconditional"
    )
  }
})

test_that("delta vcov validates strings and does not silently change targets", {
  skip_on_cran()
  object <- make_delta_average_fixture()
  local_reproducible_output(width = 80)
  expect_snapshot(inferences(object, method = "delta"), error = TRUE)
  expect_snapshot(
    inferences(object, method = "delta", vcov = NA_character_),
    error = TRUE
  )
  expect_snapshot(
    inferences(object, method = "delta", vcov = character()),
    error = TRUE
  )
  expect_snapshot(
    inferences(object, method = "delta", vcov = "cond"),
    error = TRUE
  )
  expect_snapshot(
    inferences(object, method = "delta", target = "empirical"),
    error = TRUE
  )
  expect_error(
    mostr:::inferences_score_bootstrap(object, vcov = "conditional"),
    "Character `vcov` choices are only available with `method = \"delta\"`.",
    fixed = TRUE
  )
})

test_that("variance and covariance-matrix checks respect covariance units", {
  # The same checks must hold whatever the units of the coefficients, so
  # every tolerance has to scale with the covariance.
  for (scale in c(1e-20, 1, 1e20)) {
    label <- paste("covariance scale", scale)

    # Variances from a Jacobian and a covariance matrix.
    expect_equal(
      delta_jacobian_variance(matrix(1), matrix(scale)) / scale,
      1,
      info = label
    )
    expect_match(
      tryCatch(
        delta_jacobian_variance(matrix(1), matrix(-scale)),
        error = conditionMessage
      ),
      "negative variance",
      info = label
    )
    cancellation <- diag(c(1, -1 - .Machine$double.eps)) * scale
    expect_equal(
      delta_jacobian_variance(matrix(c(1, 1), nrow = 1), cancellation),
      0,
      info = label
    )

    # Validation of a named covariance matrix.
    covariance <- diag(c(1, 2)) * scale
    dimnames(covariance) <- list(c("a", "b"), c("a", "b"))
    expect_equal(
      delta_validate_named_matrix(
        covariance,
        c("a", "b"),
        "vcov",
        positive_definite = TRUE
      ) /
        scale,
      covariance / scale,
      info = label
    )
    invalid <- covariance
    invalid[1, 1] <- -scale
    expect_match(
      tryCatch(
        delta_validate_named_matrix(
          invalid,
          c("a", "b"),
          "vcov",
          positive_semidefinite = TRUE
        ),
        error = conditionMessage
      ),
      "not numerically positive semidefinite",
      info = label
    )
    invalid <- covariance
    invalid[1, 2] <- scale / 10
    expect_match(
      tryCatch(
        delta_validate_named_matrix(invalid, c("a", "b"), "vcov"),
        error = conditionMessage
      ),
      "not numerically symmetric",
      info = label
    )
  }

  # A zero matrix is positive semidefinite but not positive definite.
  zero <- covariance * 0
  expect_equal(
    delta_validate_named_matrix(
      zero,
      c("a", "b"),
      "vcov",
      positive_semidefinite = TRUE
    ),
    zero
  )
  expect_match(
    tryCatch(
      delta_validate_named_matrix(
        zero,
        c("a", "b"),
        "bread",
        positive_definite = TRUE
      ),
      error = conditionMessage
    ),
    "not numerically positive definite"
  )
})

test_that("SOP replay allows roundoff while detecting changed estimates", {
  point <- data.frame(time = 1:2, estimate = c(0, 0.8))
  replay <- point
  replay$estimate <- replay$estimate + c(5e-13, 5e-11)
  expect_identical(delta_match_cells(point, replay, "time"), 1:2)
  replay$estimate[2] <- replay$estimate[2] + 1e-8
  expect_match(
    tryCatch(
      delta_match_cells(point, replay, "time"),
      error = conditionMessage
    ),
    "did not reproduce"
  )
})
