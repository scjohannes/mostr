markov_sandwich_adjustment <- function(
  n,
  p,
  n_clusters,
  type = c("HC0", "HC1"),
  cadjust = FALSE
) {
  type <- match.arg(type)
  if (!is.logical(cadjust) || length(cadjust) != 1L || is.na(cadjust)) {
    mostr_abort("`cadjust` must be TRUE or FALSE.")
  }

  adjustment_factor <- 1
  if (type == "HC1") {
    if (n <= p) {
      mostr_abort(
        "`type = \"HC1\"` requires more observations than parameters."
      )
    }
    adjustment_factor <- adjustment_factor * ((n - 1) / (n - p))
  }
  if (cadjust) {
    adjustment_clusters <- n_clusters %||% n
    if (adjustment_clusters < 2L) {
      mostr_abort(
        "`cadjust = TRUE` requires at least two independent clusters."
      )
    }
    adjustment_factor <- adjustment_factor *
      (adjustment_clusters / (adjustment_clusters - 1))
  }

  adjustment_factor
}

markov_assemble_sandwich <- function(
  sandwich_scores,
  bread,
  n,
  p,
  n_clusters,
  type = c("HC0", "HC1"),
  cadjust = FALSE
) {
  type <- match.arg(type)
  adjustment_factor <- markov_sandwich_adjustment(
    n = n,
    p = p,
    n_clusters = n_clusters,
    type = type,
    cadjust = cadjust
  )
  meat <- crossprod(sandwich_scores) * adjustment_factor
  transformed_scores <- sandwich_scores %*% bread
  covariance <- crossprod(transformed_scores) * adjustment_factor
  covariance <- (covariance + t(covariance)) / 2

  list(
    meat = meat,
    covariance = covariance,
    adjustment_factor = adjustment_factor
  )
}

# TRUE when the symmetric matrix `x` (a covariance, bread or information) is
# positive definite with room to spare for inversion. The orm and vglm robust
# covariances both use this rule, so that the two engines accept and refuse the
# same matrices.
markov_correlation_positive_definite <- function(x) {
  if (any(!is.finite(x))) {
    return(FALSE)
  }
  x <- x / 2 + t(x) / 2
  variances <- diag(x)
  # A positive definite matrix has positive diagonal entries. A zero or
  # negative one cannot be scaled to a correlation matrix (cov2cor() would
  # divide by its square root), and it already shows that the matrix is not
  # positive definite.
  if (!all(variances > 0)) {
    return(FALSE)
  }
  eigenvalues <- eigen(
    stats::cov2cor(x),
    symmetric = TRUE,
    only.values = TRUE
  )$values
  # The correlation matrix is checked so that covariate units do not matter: a
  # raw check refuses a well-conditioned covariance whose variances are all
  # about 1e-16. On the correlation scale the largest eigenvalue lies between 1
  # and p. Refuse when the smallest is at most p * 1e-10 times the largest:
  # inverting it would lose about 10 of 16 digits. On simulated trial fits the
  # ratio is 2e-4 to 3e-3 (measured Nov 2026), about 1e5 times above the
  # cut-off; only nearly duplicated covariates reach it.
  tolerance <- max(1, max(abs(eigenvalues))) * length(variances) * 1e-10
  min(eigenvalues) > tolerance
}

orm_case_weights <- function(model, n) {
  weights <- model$weights
  if (is.null(weights) || length(weights) == 0L) {
    return(rep(1, n))
  }
  if (!is.numeric(weights) || length(weights) != n) {
    mostr_abort(
      "Fitted orm case weights must be a numeric vector with one value per fitted row."
    )
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    mostr_abort("Fitted orm case weights must be finite and nonnegative.")
  }
  as.numeric(weights)
}

# Since rms 8.2-0, the "LM" algorithm relaxes its damping slowly (by 1.5 instead
# of 10 after each accepted step), so it often needs more than orm()'s default
# 30 iterations: suggest it only together with a larger `maxit`. Measured Nov
# 2026 on simulated ACTT-2 data: with `maxit = 100` LM converged in every case
# under rms 8.1-1 and 8.2-1, at most 98 iterations. Small data sets where a
# state never follows some previous state have no finite estimates, which is
# why merging rare states is suggested too.
stop_failed_orm_fit <- function() {
  mostr_abort(
    "`rms::orm()` could not fit the model: its fitting algorithm stopped ",
    "without finding the coefficient estimates. You can try another ",
    "algorithm with more iterations, `opt_method = \"LM\", maxit = 100`, ",
    "scale the covariates with `scale = TRUE`, or simplify the model, for ",
    "example by removing covariates that are nearly collinear or by merging ",
    "states that are rarely observed."
  )
}

validate_robcov_orm_fit <- function(fit) {
  if (!inherits(fit, "orm")) {
    mostr_abort("'fit' must be an orm object.")
  }
  if (isTRUE(fit$fail)) {
    stop_failed_orm_fit()
  }
  if (orm_has_random_effects(fit)) {
    stop_unsupported_orm_random_effects()
  }
  if (!identical(fit$family, "logistic")) {
    mostr_abort(
      "`robcov_orm()` supports only ordinary ordinal logistic orm fits."
    )
  }
  if (is.null(fit$x) || is.null(fit$y)) {
    mostr_abort(
      "orm robust covariance requires `x = TRUE, y = TRUE` in the fitted model."
    )
  }

  coefficients <- stats::coef(fit)
  if (
    !is.numeric(coefficients) ||
      !length(coefficients) ||
      is.null(names(coefficients)) ||
      anyNA(names(coefficients)) ||
      any(names(coefficients) == "") ||
      anyDuplicated(names(coefficients)) ||
      any(!is.finite(coefficients))
  ) {
    mostr_abort(
      "orm coefficients must be finite and have unique, non-missing names."
    )
  }
  coefficients
}

orm_model_bread <- function(fit) {
  coefficients <- stats::coef(fit)
  penalty_matrix <- fit$penalty.matrix
  penalty <- unlist(fit$penalty, use.names = FALSE)
  penalized <-
    (!is.null(penalty_matrix) &&
      length(penalty_matrix) > 0L &&
      any(is.finite(penalty_matrix) & penalty_matrix != 0)) ||
    (!is.null(penalty) &&
      length(penalty) > 0L &&
      any(is.finite(penalty) & penalty != 0))
  requested_penalty_variance <- fit$call$var.penalty
  requested_sandwich <- !is.null(requested_penalty_variance) &&
    identical(as.character(requested_penalty_variance)[1L], "sandwich")
  sandwich_penalty <- penalized &&
    (!is.null(fit$var.from.info.matrix) || requested_sandwich)

  if (sandwich_penalty && is.null(fit$var.from.info.matrix)) {
    mostr_abort(
      "The penalized orm fit used `var.penalty = \"sandwich\"` but does not ",
      "retain the required inverse-sensitivity matrix in ",
      "`var.from.info.matrix`."
    )
  }

  source <- "stats::vcov(fit, intercepts = \"all\")"
  if (sandwich_penalty) {
    bread <- fit$var.from.info.matrix
    source <- "var.from.info.matrix"
  } else {
    bread <- fit$orig.var
    if (!is.null(bread)) {
      source <- "orig.var"
    } else {
      bread <- tryCatch(
        stats::vcov(fit, intercepts = "all"),
        error = function(e) NULL
      )
    }
  }
  if (is.null(bread)) {
    mostr_abort(
      "Could not obtain the full model-based covariance for the orm fit. ",
      "Refit with `x = TRUE, y = TRUE`."
    )
  }

  bread <- delta_normalize_orm_backend_matrix(bread, fit, names(coefficients))
  bread <- tryCatch(
    validate_coef_vcov(coefficients, bread, arg = "orm model bread"),
    error = function(error) {
      if (sandwich_penalty) {
        mostr_abort(
          "The penalized orm inverse-sensitivity matrix stored in ",
          "`var.from.info.matrix` is invalid: ",
          conditionMessage(error)
        )
      }
      mostr_abort(error)
    }
  )
  symmetry_error <- max(abs(bread - t(bread))) / max(1, max(abs(bread)))
  if (!is.finite(symmetry_error) || symmetry_error > 1e-8) {
    mostr_abort(
      if (sandwich_penalty) {
        "The penalized orm `var.from.info.matrix` is not numerically symmetric."
      } else {
        "`orm model bread` is not numerically symmetric."
      }
    )
  }
  if (!markov_correlation_positive_definite(bread)) {
    mostr_abort(
      if (sandwich_penalty) {
        "The penalized orm `var.from.info.matrix` must be positive definite."
      } else {
        "`orm model bread` must be positive definite."
      }
    )
  }

  list(bread = bread, source = source)
}

orm_robust_metadata <- function(model) {
  attr(model, "markov_robust_covariance", exact = TRUE)
}

orm_stored_covariance_valid <- function(model, cluster = NULL) {
  metadata <- orm_robust_metadata(model)
  if (
    !inherits(model, "orm") ||
      is.null(metadata) ||
      !identical(metadata$backend, "orm") ||
      !identical(metadata$implementation, "mostr::robcov_orm") ||
      is.null(model$var) ||
      !identical(model$var, metadata$covariance_identity)
  ) {
    return(FALSE)
  }
  if (!is.null(cluster) && !identical(cluster, metadata$cluster)) {
    return(FALSE)
  }
  TRUE
}

robcov_orm <- function(
  fit,
  cluster = NULL,
  type = c("HC0", "HC1"),
  cadjust = NULL
) {
  coefficients <- validate_robcov_orm_fit(fit)
  type <- match.arg(type)
  n <- nrow(fit$x)
  p <- length(coefficients)
  if (!identical(length(fit$y), n)) {
    mostr_abort(
      "Stored orm `x` and `y` components have inconsistent row counts."
    )
  }

  scores <- compute_scores_orm(fit)
  if (
    !is.matrix(scores) ||
      !identical(dim(scores), c(n, p)) ||
      !identical(colnames(scores), names(coefficients)) ||
      any(!is.finite(scores))
  ) {
    mostr_abort("Observation-level orm scores are incomplete or misaligned.")
  }
  weights <- orm_case_weights(fit, n)
  weighted <- any(weights != 1)
  represented <- weights > 0
  n_represented <- sum(represented)
  if (n_represented == 0L) {
    mostr_abort(
      "orm robust covariance requires at least one positive-weight fitted row."
    )
  }

  if (!is.null(cluster)) {
    cluster <- align_cluster_orm(fit, cluster, n)
    if (anyNA(cluster)) {
      mostr_abort(
        "'cluster' contains missing values after alignment with the fitted data."
      )
    }
    represented_cluster <- cluster[represented]
    cluster_group <- as.factor(represented_cluster)
    sandwich_scores <- rowsum(
      scores[represented, , drop = FALSE],
      cluster_group,
      reorder = FALSE
    )
    n_clusters <- nrow(sandwich_scores)
    if (n_clusters < 2L) {
      mostr_abort(
        "Cluster-robust covariance requires at least two clusters represented ",
        "by positive-weight rows."
      )
    }
  } else {
    sandwich_scores <- scores[represented, , drop = FALSE]
    n_clusters <- NULL
  }
  # By default a patient-clustered variance is multiplied by G / (G - 1), with
  # G the number of patients, as sandwich::vcovCL() does (`cadjust = TRUE`,
  # HC0 for models other than lm); rms::robcov() applies no factor. With HC0
  # this also makes the variance equal cov(influence) / G, the convention of
  # the unconditional variance. Without a cluster no factor is applied.
  if (is.null(cadjust)) {
    cadjust <- !is.null(cluster)
  }
  if (!is.logical(cadjust) || length(cadjust) != 1L || is.na(cadjust)) {
    mostr_abort("`cadjust` must be TRUE or FALSE.")
  }

  bread_result <- orm_model_bread(fit)
  assembled <- markov_assemble_sandwich(
    sandwich_scores = sandwich_scores,
    bread = bread_result$bread,
    n = n_represented,
    p = p,
    n_clusters = n_clusters,
    type = type,
    cadjust = cadjust
  )
  dimnames(assembled$meat) <- list(names(coefficients), names(coefficients))
  dimnames(assembled$covariance) <- list(
    names(coefficients),
    names(coefficients)
  )

  result <- fit
  result$orig.var <- bread_result$bread
  result$var <- assembled$covariance
  result$clusterInfo <- list(
    name = if (is.null(cluster)) NULL else "cluster",
    n = n_clusters %||% n_represented
  )
  attr(result, "markov_robust_covariance") <- list(
    implementation = "mostr::robcov_orm",
    backend = "orm",
    cluster = cluster,
    type = type,
    cadjust = cadjust,
    adjustment_factor = assembled$adjustment_factor,
    n_rows = n,
    n_represented_rows = n_represented,
    n_parameters = p,
    n_clusters = n_clusters,
    bread_convention = bread_result$source,
    weighted = weighted,
    covariance_identity = assembled$covariance
  )
  result
}

compute_scores_orm <- function(model) {
  if (!inherits(model, "orm")) {
    mostr_abort("'model' must be an orm object.")
  }
  if (is.null(model$x) || is.null(model$y)) {
    mostr_abort(
      "orm score calculations require the model to be fitted with ",
      "`x = TRUE, y = TRUE`."
    )
  }

  beta <- stats::coef(model)
  Gamma <- get_effective_coefs(model, beta = beta)
  X <- cbind("(Intercept)" = 1, as.matrix(model$x))
  X <- X[, colnames(Gamma), drop = FALSE]

  eta <- X %*% t(Gamma)

  y_levels <- as_state_labels(model$yunique)
  y_idx <- match(as_state_labels(model$y), y_levels)
  if (anyNA(y_idx)) {
    mostr_abort("Could not match orm response values to fitted outcome levels.")
  }

  n <- nrow(X)
  M <- nrow(Gamma)
  d_eta <- matrix(0, nrow = n, ncol = M)

  # Derivatives of each row's log-likelihood with respect to eta, where
  # F = plogis and F(eta_k) = P(Y > state k). The textbook form
  # F(eta)(1 - F(eta)) / P(Y = observed state) is not used: when the observed
  # state's probability is below .Machine$double.eps, `1 - plogis(eta)` or
  # `plogis(eta_{k-1}) - plogis(eta_k)` cancels to 0, and the old floor on the
  # probability turned the score of exactly these outlying rows into 0
  # (issue #48). The forms below cancel the common factor algebraically, use
  # upper tails (`lower.tail = FALSE`) instead of `1 - plogis()`, and take
  # ratios of tail probabilities on the log scale so that they stay exact even
  # when both tails underflow.
  lowest <- y_idx == 1L
  d_eta[lowest, 1L] <- -stats::plogis(eta[lowest, 1L])
  highest <- y_idx == M + 1L
  d_eta[highest, M] <- stats::plogis(eta[highest, M], lower.tail = FALSE)

  middle <- which(y_idx > 1L & y_idx <= M)
  if (length(middle) > 0L) {
    k <- y_idx[middle]
    upper <- cbind(middle, k - 1L)
    lower <- cbind(middle, k)
    a <- eta[upper]
    b <- eta[lower]
    # eta_{k-1} - eta_k from the coefficient differences rather than a - b,
    # which would lose precision when |eta| is large. With
    # P(Y = k) = F(a) F(-b) (1 - exp(b - a)), the scores are
    # F(-a) / (F(-b) (1 - exp(b - a))) and -F(b) / (F(a) (1 - exp(b - a))).
    gap <- rowSums(
      X[middle, , drop = FALSE] *
        (Gamma[k - 1L, , drop = FALSE] - Gamma[k, , drop = FALSE])
    )
    # Crossed cut-points give the observed state a probability of zero or
    # less; its log-likelihood and score are then undefined.
    crossed <- !is.na(gap) & gap <= 0
    if (any(crossed)) {
      mostr_abort(
        "Cannot compute orm scores: the fitted cut-points cross for ",
        sum(crossed),
        " row(s), so their observed states have a fitted probability of ",
        "zero or less."
      )
    }
    one_minus_ratio <- -expm1(-gap)
    d_eta[upper] <- exp(
      stats::plogis(a, lower.tail = FALSE, log.p = TRUE) -
        stats::plogis(b, lower.tail = FALSE, log.p = TRUE)
    ) /
      one_minus_ratio
    d_eta[lower] <- -exp(
      stats::plogis(b, log.p = TRUE) - stats::plogis(a, log.p = TRUE)
    ) /
      one_minus_ratio
  }

  intercept_scores <- d_eta
  slope_scores <- as.matrix(model$x) * rowSums(d_eta)
  scores <- cbind(intercept_scores, slope_scores)
  weights <- orm_case_weights(model, n)
  scores <- scores * weights
  colnames(scores) <- names(beta)
  scores
}

align_cluster_orm <- function(model, cluster, n_scores) {
  if (inherits(cluster, "formula")) {
    mostr_abort("`cluster` for an orm fit must be a row-level vector.")
  }

  if (length(cluster) == n_scores) {
    return(cluster)
  }

  omitted <- na_action_omitted_rows(model$na.action)
  if (
    length(omitted) > 0L &&
      length(cluster) > n_scores &&
      max(omitted) <= length(cluster)
  ) {
    keep <- seq_along(cluster)[-omitted]
    if (length(keep) == n_scores) {
      return(cluster[keep])
    }
  }

  mostr_abort(
    "Length of `cluster` (",
    length(cluster),
    ") does not match the number of orm score rows (",
    n_scores,
    "). Supply the row-level patient ID vector used to fit the model."
  )
}
