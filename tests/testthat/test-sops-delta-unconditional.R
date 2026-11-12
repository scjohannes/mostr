local_delta_vglm_fit <- local({
  fit <- NULL
  function() {
    skip_if_not_installed("VGAM")
    if (is.null(fit)) {
      data <- make_test_data(
        n_patients = 35,
        follow_up_time = 6,
        seed = 2401
      )
      fit <<- mostm(
        ordered(y) ~ time_lin + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        engine = "vglm"
      )
    }
    fit
  }
})

local_delta_orm_fit <- local({
  fit <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(fit)) {
      data <- make_test_data(
        n_patients = 35,
        follow_up_time = 6,
        seed = 2402
      )
      fit <<- mostm(
        ordered(y) ~ time + tx + yprev,
        data = data,
        id_var = "id"
      )
    }
    fit
  }
})

test_that("explicit zero ORM penalties preserve unpenalized inference", {
  skip_on_cran()
  skip_if_not_installed("rms")
  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 2402)
  reference <- local_delta_orm_fit()
  infer <- function(model, covariance) {
    avg <- avg_sops(
      model,
      variables = list(tx = c(0, 1)),
      times = 1:3,
      absorb = 6
    )
    inferences(avg, method = "delta", vcov = covariance)
  }
  for (variance in c("simple", "sandwich")) {
    model <- mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      penalty = 0,
      var.penalty = variance
    )
    expect_type(model$penalty, "list")
    expect_equal(orm_model_bread(model)$bread, orm_model_bread(reference)$bread)
    for (covariance in c("conditional", "unconditional")) {
      actual <- infer(model, covariance)
      expected <- infer(reference, covariance)
      expect_equal(actual$estimate, expected$estimate)
      expect_equal(stats::vcov(actual), stats::vcov(expected))
    }
  }
})

delta_weighted_sop_estimate <- function(
  model,
  coefficient,
  prediction_data,
  profile_weights,
  variables,
  times,
  y_levels,
  absorb
) {
  weighted_model <- set_coef(model, coefficient)
  probabilities <- soprob_markov(
    model = weighted_model,
    newdata = prediction_data,
    times = times,
    y_levels = y_levels,
    absorb = absorb
  )
  n <- length(profile_weights)
  unlist(
    lapply(seq_len(nrow(do.call(expand.grid, variables))), function(i) {
      rows <- (i - 1L) * n + seq_len(n)
      as.vector(apply(
        probabilities[rows, , , drop = FALSE],
        c(2, 3),
        function(x) {
          stats::weighted.mean(x, profile_weights)
        }
      ))
    }),
    use.names = FALSE
  )
}

delta_weighted_refit_oracle <- function(
  model,
  fit_data,
  profile_data,
  refit,
  patient_ids,
  variables = list(tx = c(0, 1)),
  times = 1:2,
  y_levels = 1:6,
  absorb = 6,
  step = 1e-2
) {
  point <- avg_sops(
    model,
    variables = variables,
    times = times,
    y_levels = y_levels,
    absorb = absorb
  )
  inferred <- inferences(
    point,
    method = "delta",
    vcov = "unconditional"
  )
  analytical <- attr(inferred, "analytical")
  profile_ids <- as.character(profile_data$id)
  n <- length(profile_ids)
  prediction_data <- attr(point, "newdata_pred")
  base_weights <- rep(1, n)
  components <- get_delta_score_components(model)

  numerical_coefficient <- matrix(
    NA_real_,
    nrow = length(patient_ids),
    ncol = length(stats::coef(model)),
    dimnames = list(patient_ids, names(stats::coef(model)))
  )
  numerical_sop <- matrix(
    NA_real_,
    nrow = length(patient_ids),
    ncol = nrow(point),
    dimnames = list(patient_ids, NULL)
  )
  profile_only_sop <- numerical_sop

  for (patient_id in patient_ids) {
    fit_patient <- as.character(fit_data$id) == patient_id
    profile_patient <- profile_ids == patient_id
    upper_fit_weights <- lower_fit_weights <- rep(1, nrow(fit_data))
    upper_fit_weights[fit_patient] <- 1 + step
    lower_fit_weights[fit_patient] <- 1 - step
    upper_profile_weights <- lower_profile_weights <- base_weights
    upper_profile_weights[profile_patient] <- 1 + step
    lower_profile_weights[profile_patient] <- 1 - step

    upper_coefficient <- refit(upper_fit_weights)
    lower_coefficient <- refit(lower_fit_weights)
    numerical_coefficient[patient_id, ] <- n *
      (upper_coefficient - lower_coefficient) /
      (2 * step)

    upper_sop <- delta_weighted_sop_estimate(
      model = model,
      coefficient = upper_coefficient,
      prediction_data = prediction_data,
      profile_weights = upper_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    lower_sop <- delta_weighted_sop_estimate(
      model = model,
      coefficient = lower_coefficient,
      prediction_data = prediction_data,
      profile_weights = lower_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    numerical_sop[patient_id, ] <- n *
      (upper_sop - lower_sop) /
      (2 * step)

    upper_profile_only <- delta_weighted_sop_estimate(
      model = model,
      coefficient = stats::coef(model),
      prediction_data = prediction_data,
      profile_weights = upper_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    lower_profile_only <- delta_weighted_sop_estimate(
      model = model,
      coefficient = stats::coef(model),
      prediction_data = prediction_data,
      profile_weights = lower_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    profile_only_sop[patient_id, ] <- n *
      (upper_profile_only - lower_profile_only) /
      (2 * step)
  }

  score_rows <- match(patient_ids, components$ids)
  expected_coefficient <- components$scores[score_rows, , drop = FALSE] %*%
    t(components$Ainv)
  rownames(expected_coefficient) <- patient_ids
  influence_rows <- match(patient_ids, analytical$profile_ids)

  list(
    inferred = inferred,
    numerical_coefficient = numerical_coefficient,
    expected_coefficient = expected_coefficient,
    numerical_sop = numerical_sop,
    profile_only_sop = profile_only_sop,
    expected_sop = analytical$influence[influence_rows, , drop = FALSE]
  )
}

# Checks shared by the vglm and orm patient-weight refit tests: the refit
# derivatives agree with the analytical influence, and the covariance built
# from them agrees with the unconditional covariance of the averaged SOPs.
expect_patient_weight_refits_match <- function(oracle, patient_ids) {
  expect_equal(
    oracle$numerical_coefficient,
    oracle$expected_coefficient,
    tolerance = 2e-3
  )
  expect_equal(
    oracle$numerical_sop,
    oracle$expected_sop,
    tolerance = 2e-3
  )
  expect_gt(max(abs(oracle$numerical_sop - oracle$profile_only_sop)), 1e-3)
  # Each derivative changes fitting and averaging weights together. Using every
  # patient checks the joint contribution, including cross-covariances, without
  # constructing scores, sensitivity, or coefficient derivatives in the oracle.
  n <- length(patient_ids)
  centered <- sweep(oracle$numerical_sop, 2L, colMeans(oracle$numerical_sop))
  covariance <- crossprod(centered) / (n * (n - 1))
  expect_lt(
    max(abs(stats::vcov(oracle$inferred) - covariance)) / max(abs(covariance)),
    2e-3
  )
  expect_equal(
    oracle$inferred$std.error,
    sqrt(diag(covariance)),
    tolerance = 2e-3
  )
  without_cross_term <- (stats::cov(oracle$profile_only_sop) +
    stats::cov(oracle$numerical_sop - oracle$profile_only_sop)) /
    n
  expect_gt(
    max(abs(covariance - without_cross_term)) / max(abs(covariance)),
    2e-3
  )
}

test_that("patient clusters resolve explicitly or from stored fitting metadata", {
  skip_on_cran()
  fit <- local_delta_vglm_fit()
  fitted_data <- attr(fit, "markov_data")
  explicit <- paste0("patient-", fitted_data$id)
  cases <- list(
    list(
      label = "stored patient-ID column",
      cluster = NULL,
      expected = fitted_data$id,
      source = "stored_id_var"
    ),
    list(
      label = "explicit cluster vector",
      cluster = explicit,
      expected = explicit,
      source = "explicit_vector"
    ),
    list(
      label = "explicit cluster formula",
      cluster = ~id,
      expected = fitted_data$id,
      source = "explicit_formula"
    )
  )
  for (case in cases) {
    resolved <- resolve_delta_cluster(fit, cluster = case$cluster)
    expect_equal(resolved$cluster, case$expected, info = case$label)
    expect_equal(resolved$metadata$source, case$source, info = case$label)
  }
  stored <- resolve_delta_cluster(fit)
  expect_equal(stored$ids, unique(as.character(fitted_data$id)))
  expect_equal(stored$metadata$cluster_level, "patient")

  # Without the stored fitting data and patient-ID column, the clusters never
  # fall back to the fitting rows.
  attr(fit, "markov_data") <- NULL
  attr(fit, "markov_id_var") <- NULL
  attr(fit$vglm_fit, "markov_data") <- NULL
  attr(fit$vglm_fit, "markov_id_var") <- NULL
  expect_error(
    resolve_delta_cluster(fit),
    "the model does not record its fitting data and patient-ID column",
    fixed = TRUE
  )
})

test_that("explicit coefficient covariance has precedence and strict validation", {
  skip_on_cran()
  fit <- local_delta_vglm_fit()
  attr(fit, "markov_data") <- NULL
  attr(fit, "markov_id_var") <- NULL
  attr(fit$vglm_fit, "markov_data") <- NULL
  attr(fit$vglm_fit, "markov_id_var") <- NULL
  coefficient_names <- names(fit$coefficients)
  reversed_names <- rev(coefficient_names)
  covariance <- diag(seq_along(reversed_names))
  dimnames(covariance) <- list(reversed_names, reversed_names)

  resolved <- get_delta_cluster_vcov(fit, vcov = covariance)
  expect_identical(rownames(resolved$vcov), coefficient_names)
  expect_identical(colnames(resolved$vcov), coefficient_names)
  expect_equal(resolved$metadata$source, "explicit")
  expect_null(resolved$cluster_info)

  unnamed <- covariance
  dimnames(unnamed) <- NULL
  expect_error(
    get_delta_cluster_vcov(fit, vcov = unnamed),
    "must be uniquely named",
    fixed = TRUE
  )

  nonfinite <- covariance
  nonfinite[1L, 1L] <- Inf
  expect_error(
    get_delta_cluster_vcov(fit, vcov = nonfinite),
    "contains non-finite values",
    fixed = TRUE
  )

  nonsymmetric <- covariance
  nonsymmetric[1L, 2L] <- 1
  expect_error(
    get_delta_cluster_vcov(fit, vcov = nonsymmetric),
    "not numerically symmetric",
    fixed = TRUE
  )

  non_psd <- covariance
  non_psd[1L, 1L] <- -1
  expect_error(
    get_delta_cluster_vcov(fit, vcov = non_psd),
    "not numerically positive semidefinite",
    fixed = TRUE
  )
})

test_that("vglm covariance and score components preserve backend conventions", {
  skip_on_cran()
  fit <- local_delta_vglm_fit()
  fitted_data <- attr(fit, "markov_data")
  cluster_ids <- unique(as.character(fitted_data$id))
  expected_scores <- rowsum(
    fit$scores,
    factor(as.character(fitted_data$id), levels = cluster_ids),
    reorder = FALSE
  )
  rownames(expected_scores) <- cluster_ids

  # The robust fit that `mostm()` stores, and the raw `vglm` fit inside it with
  # the patient clusters given explicitly. Both give the stored patient
  # covariance and the same patient scores.
  cases <- list(
    list(
      label = "robcov_vglm fit from mostm()",
      model = fit,
      cluster = NULL,
      vcov_tolerance = testthat_tolerance(),
      bread_source = "robcov_vglm_observed"
    ),
    list(
      label = "raw vglm fit with explicit clusters",
      model = fit$vglm_fit,
      cluster = fitted_data$id,
      vcov_tolerance = 1e-10,
      bread_source = "vglm_observed_information"
    )
  )
  for (case in cases) {
    covariance <- get_delta_cluster_vcov(case$model, cluster = case$cluster)
    components <- get_delta_score_components(
      case$model,
      cluster = case$cluster
    )
    expect_equal(
      covariance$vcov,
      fit$var,
      tolerance = case$vcov_tolerance,
      info = case$label
    )
    expect_equal(components$ids, cluster_ids, info = case$label)
    expect_equal(components$scores, expected_scores, info = case$label)
    expect_equal(components$n, length(cluster_ids), info = case$label)
    expect_equal(
      components$Ainv,
      length(cluster_ids) * fit$bread,
      info = case$label
    )
    expect_false(
      isTRUE(all.equal(components$Ainv, fit$var)),
      info = case$label
    )
    expect_equal(
      components$metadata$bread_source,
      case$bread_source,
      info = case$label
    )
    expect_equal(
      components$metadata$n_patients,
      length(cluster_ids),
      info = case$label
    )
  }

  # The stored fit keeps its backend settings in the metadata; the raw fit's
  # covariance is computed.
  stored <- get_delta_cluster_vcov(fit)
  expect_equal(stored$metadata$type, fit$type)
  expect_identical(stored$metadata$cadjust, fit$cadjust)
  expect_equal(stored$metadata$adjustment_factor, fit$adjustment_factor)
  stored_components <- get_delta_score_components(fit)
  expect_equal(stored_components$metadata$backend_hc_type_ignored, fit$type)
  expect_identical(
    stored_components$metadata$backend_cadjust_ignored,
    fit$cadjust
  )
  computed <- get_delta_cluster_vcov(fit$vglm_fit, cluster = fitted_data$id)
  expect_equal(computed$metadata$source, "computed_robcov_vglm")
})

test_that("orm covariance and score components use inverse total information", {
  fit <- local_delta_orm_fit()
  fitted_data <- attr(fit, "markov_data")
  covariance <- get_delta_cluster_vcov(fit)

  expect_equal(covariance$vcov, fit$var)
  expect_equal(covariance$metadata$type, "HC0")
  expect_true(covariance$metadata$cadjust)
  expect_equal(covariance$metadata$source, "stored_markov_robcov_orm")

  recomputed <- get_delta_cluster_vcov(fit, cluster = fitted_data$id)
  expect_equal(recomputed$vcov, fit$var, tolerance = 1e-10)
  expect_equal(recomputed$metadata$source, "computed_markov_robcov_orm")

  components <- get_delta_score_components(fit)
  row_scores <- compute_scores_orm(fit)
  cluster_ids <- unique(as.character(fitted_data$id))
  expected_scores <- rowsum(
    row_scores,
    factor(as.character(fitted_data$id), levels = cluster_ids),
    reorder = FALSE
  )
  rownames(expected_scores) <- cluster_ids
  expected_bread <- orm_model_bread(fit)$bread

  expect_equal(components$ids, cluster_ids)
  expect_equal(components$scores, expected_scores)
  expect_equal(components$Ainv, length(cluster_ids) * expected_bread)
  expect_equal(
    components$metadata$bread_source,
    "rms_inverse_total_information"
  )
})

test_that("orm delta covariance recomputes after stored covariance mutation", {
  fit <- local_delta_orm_fit()
  fitted_data <- attr(fit, "markov_data")
  stored <- get_delta_cluster_vcov(fit)

  mutated <- fit
  mutated$var[1L, 1L] <- mutated$var[1L, 1L] + 1
  recomputed <- get_delta_cluster_vcov(mutated)
  explicit <- get_delta_cluster_vcov(fit, cluster = fitted_data$id)

  expect_equal(recomputed$metadata$source, "computed_markov_robcov_orm")
  expect_equal(recomputed$vcov, stored$vcov, tolerance = 1e-12)
  expect_equal(explicit$metadata$source, "computed_markov_robcov_orm")
  expect_equal(explicit$vcov, stored$vcov, tolerance = 1e-12)
})

test_that("vglm patient-weight refits independently validate unconditional covariance", {
  skip_on_cran()
  fit <- local_delta_vglm_fit()
  fit_data <- attr(fit, "markov_data")
  profile_data <- attr(fit, "markov_starting_profile_data")
  patient_ids <- as.character(profile_data$id)

  oracle <- delta_weighted_refit_oracle(
    model = fit,
    fit_data = fit_data,
    profile_data = profile_data,
    refit = function(weights) {
      weighted <- VGAM::vglm(
        ordered(y) ~ time_lin + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = fit_data,
        weights = weights,
        coefstart = stats::coef(fit),
        control = VGAM::vglm.control(epsilon = 1e-10, maxit = 100)
      )
      stats::coef(weighted)
    },
    patient_ids = patient_ids
  )

  expect_patient_weight_refits_match(oracle, patient_ids)
})

test_that("orm patient-weight refits independently validate unconditional covariance", {
  skip_on_cran()
  fit <- local_delta_orm_fit()
  fit_data <- attr(fit, "markov_data")
  profile_data <- attr(fit, "markov_starting_profile_data")
  patient_ids <- as.character(profile_data$id)

  oracle <- delta_weighted_refit_oracle(
    model = fit,
    fit_data = fit_data,
    profile_data = profile_data,
    refit = function(weights) {
      weighted <- suppressWarnings(rms::orm(
        ordered(y) ~ time + tx + yprev,
        data = fit_data,
        weights = weights,
        x = TRUE,
        y = TRUE
      ))
      stats::coef(weighted)
    },
    patient_ids = patient_ids
  )

  expect_patient_weight_refits_match(oracle, patient_ids)
})

test_that("unconditional score machinery rejects partial proportional odds", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  data <- make_test_data(
    n_patients = 35,
    follow_up_time = 6,
    seed = 2403
  )
  fit <- mostm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(
      reverse = TRUE,
      parallel = FALSE ~ tx
    ),
    data = data,
    id_var = "id",
    engine = "vglm"
  )

  expect_error(
    resolve_delta_cluster(fit, cluster = data$id),
    "requires a full proportional-odds vglm model",
    fixed = TRUE
  )
})

test_that("stacked influence matches fitted patients and retains the cross term", {
  values <- rbind(
    p1 = c(cell_a = 0.2, cell_b = 0.8),
    p2 = c(cell_a = 0.5, cell_b = 0.5)
  )
  jacobian <- rbind(
    cell_a = c(beta_1 = 1, beta_2 = 0.5),
    cell_b = c(beta_1 = -1, beta_2 = -0.5)
  )
  scores <- rbind(
    p2 = c(beta_1 = -0.3, beta_2 = 0.4),
    p1 = c(beta_1 = 0.2, beta_2 = -0.1)
  )
  Ainv <- diag(c(2, 1))
  dimnames(Ainv) <- list(c("beta_1", "beta_2"), c("beta_1", "beta_2"))
  components <- list(
    ids = rownames(scores),
    scores = scores,
    Ainv = Ainv,
    n = nrow(scores)
  )

  result <- delta_stacked_influence(
    individual_values = values,
    average_jacobian = jacobian,
    profile_ids = rownames(values),
    score_components = components
  )

  aligned_scores <- rbind(
    p1 = scores["p1", ],
    p2 = scores["p2", ]
  )
  coefficient_influence <- aligned_scores %*% t(Ainv) %*% t(jacobian)
  profile_component <- sweep(values, 2L, colMeans(values), "-")
  expected_influence <- profile_component + coefficient_influence
  expected_covariance <- stats::cov(expected_influence) / nrow(values)
  centered <- sweep(
    expected_influence,
    2L,
    colMeans(expected_influence),
    "-"
  )
  centered_hc0 <- crossprod(centered) / nrow(values)^2
  raw_uncentered <- crossprod(expected_influence) / nrow(values)^2
  influence_mean <- colMeans(expected_influence)

  expect_equal(result$scores, aligned_scores)
  expect_equal(result$Ainv, Ainv)
  expect_equal(result$coefficient_influence, coefficient_influence)
  expect_equal(result$influence, expected_influence)
  expect_equal(result$covariance, expected_covariance)
  expect_equal(result$std.error, sqrt(diag(expected_covariance)))
  expect_equal(
    result$covariance,
    crossprod(centered) / (nrow(values) * (nrow(values) - 1L))
  )
  expect_equal(
    result$covariance,
    centered_hc0 * nrow(values) / (nrow(values) - 1L)
  )
  expect_equal(
    result$covariance,
    nrow(values) /
      (nrow(values) - 1L) *
      raw_uncentered -
      tcrossprod(influence_mean) / (nrow(values) - 1L)
  )
  expect_equal(
    result$metadata$finite_sample_method,
    "patient_level_sample_covariance"
  )
  expect_equal(result$metadata$finite_sample_factor, 2)
  expect_equal(result$metadata$n_patients, 2)
  expect_equal(result$metadata$n_score_clusters, 2)
  expect_null(result$metadata$missing_score_profiles)
  expect_null(result$metadata$sensitivity_scale)

  additive_covariance <- stats::cov(profile_component) /
    nrow(values) +
    stats::cov(coefficient_influence) / nrow(values)
  expect_gt(max(abs(result$covariance - additive_covariance)), 1e-8)
})

test_that("stacked influence rejects patients missing from profiles or scores", {
  jacobian <- matrix(1, nrow = 1L, dimnames = list("cell", "beta"))
  Ainv <- matrix(1, dimnames = list("beta", "beta"))
  cases <- list(
    list(
      label = "profile patient without likelihood scores",
      profile_ids = c("p1", "p2", "p3"),
      score_ids = c("p1", "p2"),
      message = "Target-profile patients have no usable likelihood transition: p3"
    ),
    list(
      label = "scored patient absent from the profiles",
      profile_ids = c("p1", "p2"),
      score_ids = c("p1", "p3"),
      message = "Likelihood-score patients are absent from the target profiles: p3"
    )
  )
  for (case in cases) {
    values <- matrix(
      c(0.2, 0.4, 0.6)[seq_along(case$profile_ids)],
      ncol = 1L,
      dimnames = list(case$profile_ids, "cell")
    )
    scores <- matrix(
      c(0.1, -0.1),
      ncol = 1L,
      dimnames = list(case$score_ids, "beta")
    )
    expect_error(
      delta_stacked_influence(
        values,
        jacobian,
        rownames(values),
        list(ids = rownames(scores), scores = scores, Ainv = Ainv, n = 2L)
      ),
      case$message,
      fixed = TRUE,
      info = case$label
    )
  }
})

test_that("empirical covariance propagates named raw coefficients", {
  jacobian <- rbind(
    first = c(beta_2 = 2, beta_1 = 1),
    second = c(beta_2 = -1, beta_1 = 0.5)
  )
  covariance <- matrix(
    c(4, 1, 1, 2),
    nrow = 2L,
    dimnames = list(c("beta_1", "beta_2"), c("beta_1", "beta_2"))
  )
  expected <- jacobian %*%
    covariance[colnames(jacobian), colnames(jacobian)] %*%
    t(jacobian)

  expect_equal(delta_empirical_covariance(jacobian, covariance), expected)
})
