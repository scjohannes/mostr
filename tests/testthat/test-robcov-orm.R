local_robcov_orm_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(value)) {
      data <- suppressWarnings(make_test_data(
        n_patients = 28,
        times = 1:5,
        seed = 2811
      ))
      fit <- suppressWarnings(rms::orm(
        ordered(y) ~ time + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      ))
      value <<- list(data = data, fit = fit)
    }
    value
  }
})

orm_test_row_loglik <- function(model, row, coefficient) {
  gamma <- get_effective_coefs(model, beta = coefficient)
  design <- cbind("(Intercept)" = 1, as.matrix(model$x))
  design <- design[, colnames(gamma), drop = FALSE]
  eta <- design[row, , drop = FALSE] %*% t(gamma)
  probabilities <- lp_to_probs(eta, nrow(gamma))
  response_index <- match(
    as_state_labels(model$y[row]),
    as_state_labels(model$yunique)
  )
  weight <- orm_case_weights(model, nrow(model$x))[row]
  weight * log(probabilities[1L, response_index])
}

test_that("orm HC0 covariance equals a directly assembled cluster sandwich", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  cluster <- case$data$id
  robust <- suppressWarnings(robcov_orm(
    fit,
    cluster = cluster,
    type = "HC0",
    cadjust = FALSE
  ))

  scores <- compute_scores_orm(fit)
  cluster_ids <- unique(as.character(cluster))
  clustered_scores <- rowsum(
    scores,
    factor(as.character(cluster), levels = cluster_ids),
    reorder = FALSE
  )
  bread <- orm_model_bread(fit)$bread
  expected <- crossprod(clustered_scores %*% bread)
  dimnames(expected) <- dimnames(robust$var)

  expect_equal(robust$var, expected, tolerance = 1e-12)
  expect_equal(
    attr(robust, "markov_robust_covariance")$adjustment_factor,
    1
  )
})

test_that("orm robust covariance is preserved when time units change", {
  case <- local_robcov_orm_case()
  days <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    id_var = "id",
    absorb = 6
  )
  data <- case$data
  data$minutes <- data$time * 1440
  minutes <- mostm(
    ordered(y) ~ minutes + tx + yprev,
    data = data,
    id_var = "id",
    absorb = 6
  )
  units <- ifelse(names(stats::coef(minutes)) == "minutes", 1440, 1)

  expect_equal(
    unname(minutes$orig.var * outer(units, units)),
    unname(days$orig.var),
    tolerance = 1e-7
  )
  expect_equal(
    unname(minutes$var * outer(units, units)),
    unname(days$var),
    tolerance = 1e-7
  )
})

test_that("orm bread validation rejects nonpositive and singular matrices", {
  skip_on_cran()
  fit <- local_robcov_orm_case()$fit
  coefficient_names <- names(stats::coef(fit))
  fit$orig.var <- diag(length(coefficient_names))
  dimnames(fit$orig.var) <- list(coefficient_names, coefficient_names)

  fit$orig.var[1L, 1L] <- 0
  expect_snapshot(error = TRUE, orm_model_bread(fit))
  fit$orig.var[1L, 1L] <- -1
  expect_snapshot(error = TRUE, orm_model_bread(fit))
  fit$orig.var[1L, 1L] <- 1
  fit$orig.var[1L, 2L] <- fit$orig.var[2L, 1L] <- 1
  expect_snapshot(error = TRUE, orm_model_bread(fit))
  fit$orig.var[1L, 2L] <- fit$orig.var[2L, 1L] <- 2
  expect_snapshot(error = TRUE, orm_model_bread(fit))
})

test_that("orm small-sample and cluster corrections count only represented rows and clusters", {
  skip_on_cran()
  case <- local_robcov_orm_case()

  # Weighted case: two whole clusters and one further row get weight zero.
  weights <- rep(1, nrow(case$data))
  zero_cluster_ids <- unique(case$data$id)[1:2]
  weights[case$data$id %in% zero_cluster_ids] <- 0
  weights[which(!case$data$id %in% zero_cluster_ids)[1L]] <- 0
  weighted_data <- case$data
  weighted_data$case_weight <- weights
  weighted <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = weighted_data,
    weights = case_weight,
    x = TRUE,
    y = TRUE
  ))

  cases <- list(
    list(label = "unweighted", fit = case$fit, excluded_clusters = integer()),
    list(
      label = "zero-weight rows and clusters",
      fit = weighted,
      excluded_clusters = zero_cluster_ids
    )
  )

  for (correction_case in cases) {
    label <- correction_case$label
    fit <- correction_case$fit
    cluster <- case$data$id
    represented <- orm_case_weights(fit, nrow(fit$x)) > 0
    represented_clusters <- unique(cluster[represented])
    n_represented <- sum(represented)
    n_clusters <- length(represented_clusters)
    p <- length(stats::coef(fit))

    hc0 <- robcov_orm(fit, cluster, type = "HC0", cadjust = FALSE)
    hc1 <- robcov_orm(fit, cluster, type = "HC1", cadjust = FALSE)
    cluster_adjusted <- robcov_orm(fit, cluster, type = "HC0", cadjust = TRUE)
    combined <- robcov_orm(fit, cluster, type = "HC1", cadjust = TRUE)
    hc0_metadata <- attr(hc0, "markov_robust_covariance")
    combined_metadata <- attr(combined, "markov_robust_covariance")

    hc1_factor <- (n_represented - 1) / (n_represented - p)
    cluster_factor <- n_clusters / (n_clusters - 1)
    expect_equal(hc1$var, hc0$var * hc1_factor, tolerance = 1e-12, info = label)
    expect_equal(
      cluster_adjusted$var,
      hc0$var * cluster_factor,
      tolerance = 1e-12,
      info = label
    )
    expect_equal(
      combined$var,
      hc0$var * hc1_factor * cluster_factor,
      tolerance = 1e-12,
      info = label
    )
    expect_equal(
      combined_metadata$adjustment_factor,
      hc1_factor * cluster_factor,
      info = label
    )
    expect_equal(hc0_metadata$n_rows, nrow(fit$x), info = label)
    expect_equal(hc0_metadata$n_represented_rows, n_represented, info = label)
    expect_equal(hc0_metadata$n_clusters, n_clusters, info = label)
    expect_identical(hc0_metadata$cluster, cluster, info = label)
    expect_false(
      any(correction_case$excluded_clusters %in% represented_clusters),
      info = label
    )
  }
})

test_that("orm robust covariance ignores empty cluster factor levels", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  with_empty_level <- factor(
    case$data$id,
    levels = c(unique(case$data$id), 9999)
  )
  robust <- robcov_orm(fit, with_empty_level, type = "HC1", cadjust = TRUE)
  without_empty_level <- robcov_orm(
    fit,
    case$data$id,
    type = "HC1",
    cadjust = TRUE
  )

  expect_equal(robust$var, without_empty_level$var, tolerance = 1e-12)
  expect_equal(
    attr(robust, "markov_robust_covariance")$n_clusters,
    length(unique(case$data$id))
  )
})

test_that("the orm engine stores package robust covariance metadata by default", {
  case <- local_robcov_orm_case()
  fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    id_var = "id",
    absorb = 6
  )
  metadata <- attr(fit, "markov_robust_covariance")
  n_clusters <- length(unique(case$data$id))

  expect_s3_class(fit, "orm")
  expect_equal(metadata$implementation, "mostr::robcov_orm")
  expect_equal(metadata$backend, "orm")
  expect_equal(metadata$type, "HC0")
  expect_true(metadata$cadjust)
  expect_equal(metadata$n_clusters, n_clusters)
  expect_equal(metadata$adjustment_factor, n_clusters / (n_clusters - 1))
  expect_identical(metadata$covariance_identity, fit$var)
})

test_that("weighted orm scores equal numerical weighted row likelihood derivatives", {
  case <- local_robcov_orm_case()
  weights <- rep(c(0.5, 1, 2), length.out = nrow(case$data))
  weighted <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    weights = weights,
    x = TRUE,
    y = TRUE
  ))
  scores <- compute_scores_orm(weighted)
  coefficient <- stats::coef(weighted)
  rows <- unique(c(1L, nrow(weighted$x) %/% 2L, nrow(weighted$x)))
  numerical <- central_difference_jacobian(
    function(coefficient) {
      vapply(
        rows,
        function(row) orm_test_row_loglik(weighted, row, coefficient),
        numeric(1)
      )
    },
    coefficient
  )
  rownames(numerical) <- rownames(scores)[rows]

  expect_equal(scores[rows, ], numerical, tolerance = 1e-7)

  robust <- robcov_orm(
    weighted,
    cluster = case$data$id,
    type = "HC0",
    cadjust = TRUE
  )
  expect_true(attr(robust, "markov_robust_covariance")$weighted)
  expect_true(all(is.finite(robust$var)))
  expect_true(all(diag(robust$var) >= 0))

  weighted_data <- case$data
  weighted_data$case_weight <- weights
  wrapped <- suppressWarnings(mostm(
    ordered(y) ~ time + tx + yprev,
    data = weighted_data,
    weights = case_weight,
    id_var = "id",
    absorb = 6
  ))
  average <- avg_sops(
    wrapped,
    variables = list(tx = c(0, 1)),
    times = 1:2
  )
  inferred <- inferences(average, method = "delta", vcov = "conditional")
  expect_true(all(is.finite(inferred$std.error)))
  expect_true(all(inferred$std.error >= 0))
})

test_that("orm scores stay exact when the observed state has probability below 1e-16", {
  # Issue #48: for an observed state with fitted probability below
  # .Machine$double.eps, 1 - plogis(eta) cancels to 0 and the old score
  # formula returned 0 instead of about -1 (lowest state) or +1 (highest
  # state). Gross outliers in a strong covariate produce such rows, and they
  # are the rows a robust variance should weigh most.
  skip_if_not_installed("rms")
  n <- 1000L
  x <- seq(-2, 2, length.out = n)
  # Deterministic logistic noise, so that the fit does not depend on a seed.
  noise <- stats::qlogis(stats::ppoints(n))[(seq_len(n) * 7919L) %% n + 1L]
  y <- cut(6 * x + noise, c(-Inf, -1, 1, Inf), labels = FALSE)
  # One outlier in each state, each in the state its covariate makes least
  # likely: lowest state at large x, highest state at very negative x, middle
  # state at large x (where both cut-point probabilities round to 1).
  outliers <- n + 1:3
  x <- c(x, 15, -15, 15)
  y <- c(y, 1L, 3L, 2L)
  fit <- rms::orm(y ~ x, x = TRUE, y = TRUE)
  coefficient <- stats::coef(fit)

  # Row log-likelihoods written without cancellation, as a reference that does
  # not share the package's probability code.
  row_loglik <- function(beta, row) {
    eta <- beta[1:2] + beta[[3]] * x[row]
    switch(
      y[row],
      stats::plogis(eta[1], lower.tail = FALSE, log.p = TRUE),
      stats::plogis(eta[1], log.p = TRUE) +
        stats::plogis(eta[2], lower.tail = FALSE, log.p = TRUE) +
        log(-expm1(eta[2] - eta[1])),
      stats::plogis(eta[2], log.p = TRUE)
    )
  }
  reference <- central_difference_jacobian(
    function(beta) {
      vapply(seq_along(y), function(row) row_loglik(beta, row), numeric(1))
    },
    coefficient
  )

  # The premise: each outlier's observed state has probability below 1e-16.
  observed_log_probability <- vapply(
    outliers,
    function(row) row_loglik(coefficient, row),
    numeric(1)
  )
  expect_true(all(observed_log_probability < log(1e-16)))

  scores <- compute_scores_orm(fit)
  expect_equal(
    unname(scores[outliers, ]),
    unname(reference[outliers, ]),
    tolerance = 1e-6
  )
  expect_equal(unname(scores), unname(reference), tolerance = 1e-6)

  # The robust covariance uses the exact scores.
  robust <- robcov_orm(fit, type = "HC0")
  bread <- orm_model_bread(fit)$bread
  expected <- bread %*% crossprod(reference) %*% bread
  expect_equal(unname(robust$var), unname(expected), tolerance = 1e-6)
})

test_that("orm scores refuse cut-points that cross for an observed middle state", {
  # A fitted orm model has ordered intercepts; reversed ones give a middle
  # state a probability of zero or less, where its score is undefined. The
  # old code returned a huge finite score instead.
  case <- local_robcov_orm_case()
  fit <- case$fit
  intercepts <- seq_len(length(fit$yunique) - 1L)
  fit$coefficients[intercepts] <- rev(fit$coefficients[intercepts])
  expect_error(compute_scores_orm(fit), "cut-points cross", class = "mostr_error")
})

test_that("orm robust covariance rejects invalid fits and inputs", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  cluster <- case$data$id

  missing_cluster <- cluster
  missing_cluster[1L] <- NA
  expect_error(robcov_orm(fit, integer()), "does not match")
  expect_error(robcov_orm(fit, missing_cluster), "contains missing values")
  expect_error(
    robcov_orm(fit, rep(1, length(cluster))),
    "at least two clusters"
  )

  invalid_weights <- fit
  invalid_weights$weights <- rep(1, nrow(fit$x))
  invalid_weights$weights[1L] <- -1
  expect_error(robcov_orm(invalid_weights, cluster), "finite and nonnegative")

  missing_x <- fit
  missing_x$x <- NULL
  expect_error(
    robcov_orm(missing_x, cluster),
    "x = TRUE, y = TRUE",
    fixed = TRUE
  )
  missing_y <- fit
  missing_y$y <- NULL
  expect_error(
    robcov_orm(missing_y, cluster),
    "x = TRUE, y = TRUE",
    fixed = TRUE
  )

  non_logit <- fit
  non_logit$family <- "probit"
  expect_error(robcov_orm(non_logit, cluster), "only ordinary ordinal logistic")

  too_small <- fit
  p <- length(stats::coef(too_small))
  too_small$x <- too_small$x[seq_len(p), , drop = FALSE]
  too_small$y <- too_small$y[seq_len(p)]
  too_small$weights <- rep(1, p)
  expect_error(
    robcov_orm(too_small, cluster = seq_len(p), type = "HC1"),
    "requires more observations than parameters"
  )

  malformed_bread <- fit
  malformed_bread$orig.var <- diag(2)
  expect_error(robcov_orm(malformed_bread, cluster), "orm model bread")
})

test_that("orm robust covariance explains failed fits", {
  skip_if_not_installed("rms")

  # When its fitting algorithm fails, rms 8.2 and later return only this field.
  failed <- structure(list(fail = TRUE), class = "orm")
  expect_error(robcov_orm(failed, 1:2), "could not fit the model", fixed = TRUE)
  expect_error(
    get_vcov_robust(failed, cluster = 1:2),
    "could not fit the model",
    fixed = TRUE
  )
})

test_that("orm robust covariance rejects fits that hold random effects", {
  case <- local_robcov_orm_case()

  for (field in c("sigma", "sigma1", "sigma2", "gamma")) {
    with_random_effects <- case$fit
    with_random_effects[[field]] <- 1
    expect_error(
      robcov_orm(with_random_effects, case$data$id),
      "random effects",
      fixed = TRUE
    )
  }

  skip_if_not_installed("rms", "8.2-0")
  cluster <- function(x) x
  random_intercept <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev + cluster(id),
    data = case$data,
    x = TRUE,
    y = TRUE
  ))
  expect_error(
    robcov_orm(random_intercept, case$data$id),
    "random effects",
    fixed = TRUE
  )
  expect_error(
    get_vcov_robust(random_intercept, cluster = case$data$id),
    "random effects",
    fixed = TRUE
  )
})

test_that("stored orm covariance integrity detects mutation", {
  case <- local_robcov_orm_case()
  robust <- robcov_orm(case$fit, case$data$id)
  expect_true(orm_stored_covariance_valid(robust, case$data$id))

  mutated <- robust
  mutated$var[1L, 1L] <- mutated$var[1L, 1L] + 1
  expect_false(orm_stored_covariance_valid(mutated, case$data$id))
  expect_false(orm_stored_covariance_valid(robust, rev(case$data$id)))
})

test_that("penalized sandwich orm uses its dedicated inverse-sensitivity bread", {
  case <- local_robcov_orm_case()
  penalized <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    penalty = 0.1,
    var.penalty = "sandwich",
    x = TRUE,
    y = TRUE
  ))
  expect_false(is.null(penalized$var.from.info.matrix))

  robust <- robcov_orm(
    penalized,
    cluster = case$data$id,
    type = "HC0",
    cadjust = FALSE
  )
  expected_bread <- delta_normalize_orm_backend_matrix(
    penalized$var.from.info.matrix,
    penalized,
    names(stats::coef(penalized))
  )
  metadata <- attr(robust, "markov_robust_covariance")

  expect_equal(robust$orig.var, expected_bread, tolerance = 1e-12)
  expect_equal(metadata$bread_convention, "var.from.info.matrix")

  missing_bread <- penalized
  missing_bread$var.from.info.matrix <- NULL
  expect_error(
    robcov_orm(missing_bread, cluster = case$data$id),
    "does not retain the required inverse-sensitivity matrix",
    fixed = TRUE
  )

  malformed_bread <- penalized
  malformed_bread$var.from.info.matrix <- diag(2)
  expect_error(
    robcov_orm(malformed_bread, cluster = case$data$id),
    "inverse-sensitivity matrix stored in `var.from.info.matrix` is invalid",
    fixed = TRUE
  )
})

test_that("penalized orm sandwich subtracts each patient's share of the penalty gradient", {
  case <- local_robcov_orm_case()
  penalized <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    penalty = 2,
    x = TRUE,
    y = TRUE
  ))
  cluster <- case$data$id
  beta <- stats::coef(penalized)
  penalty_matrix <- penalized$penalty.matrix
  # rms penalizes the non-intercept coefficients, which come last.
  penalized_index <- length(beta) - nrow(penalty_matrix) + seq_len(nrow(penalty_matrix))
  penalty_gradient <- stats::setNames(numeric(length(beta)), names(beta))
  penalty_gradient[penalized_index] <- penalty_matrix %*% beta[penalized_index]

  # A penalized fit maximizes loglik - 0.5 * beta' P beta, so at its estimate
  # the log-likelihood scores sum to the penalty gradient P beta, not to zero.
  scores <- compute_scores_orm(penalized)
  expect_gt(max(abs(penalty_gradient)), 1)
  expect_equal(colSums(scores), penalty_gradient, tolerance = 1e-6)

  # Each patient's contribution to the penalized estimating equations is the
  # patient's score sum minus an equal share of P beta; the contributions sum
  # to zero. Leaving P beta in the scores adds (P beta)(P beta)' / G to the
  # meat, which inflated coefficient SEs by up to 2 times at penalty 100 (#49).
  patient_scores <- rowsum(scores, cluster, reorder = FALSE)
  contributions <- sweep(
    patient_scores,
    2,
    penalty_gradient / nrow(patient_scores)
  )
  bread <- orm_model_bread(penalized)$bread
  expected <- bread %*% crossprod(contributions) %*% bread
  expected <- (expected + t(expected)) / 2
  robust <- robcov_orm(
    penalized,
    cluster = cluster,
    type = "HC0",
    cadjust = FALSE
  )
  expect_equal(unname(robust$var), unname(expected), tolerance = 1e-10)

  # Without a cluster each row is its own unit and carries 1 / n of P beta.
  row_contributions <- sweep(scores, 2, penalty_gradient / nrow(scores))
  expected_rows <- bread %*% crossprod(row_contributions) %*% bread
  expected_rows <- (expected_rows + t(expected_rows)) / 2
  robust_rows <- robcov_orm(penalized, type = "HC0", cadjust = FALSE)
  expect_equal(unname(robust_rows$var), unname(expected_rows), tolerance = 1e-10)
})

test_that("get_vcov_robust inherits and explicitly overrides stored orm corrections", {
  case <- local_robcov_orm_case()
  stored <- robcov_orm(
    case$fit,
    cluster = case$data$id,
    type = "HC1",
    cadjust = FALSE
  )

  expect_identical(get_vcov_robust(stored), stored$var)
  expect_equal(
    get_vcov_robust(stored, cluster = case$data$id),
    stored$var,
    tolerance = 1e-12
  )

  expected_type_override <- robcov_orm(
    case$fit,
    cluster = case$data$id,
    type = "HC0",
    cadjust = FALSE
  )$var
  expect_equal(
    get_vcov_robust(
      stored,
      cluster = case$data$id,
      type = "HC0"
    ),
    expected_type_override,
    tolerance = 1e-12
  )

  expected_cadjust_override <- robcov_orm(
    case$fit,
    cluster = case$data$id,
    type = "HC1",
    cadjust = TRUE
  )$var
  expect_equal(
    get_vcov_robust(
      stored,
      cluster = case$data$id,
      cadjust = TRUE
    ),
    expected_cadjust_override,
    tolerance = 1e-12
  )
})
