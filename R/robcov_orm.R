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

stop_failed_orm_fit <- function() {
  mostr_abort(
    "`rms::orm()` could not fit the model: its fitting algorithm stopped ",
    "without finding the coefficient estimates. You can try another ",
    "algorithm with `opt_method = \"LM\"`, scale the covariates with ",
    "`scale = TRUE`, allow more iterations with a larger `maxit`, or ",
    "simplify the model, for example by removing covariates that are ",
    "nearly collinear."
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
  eigenvalues <- if (all(diag(bread) > 0)) {
    eigen(
      stats::cov2cor(bread / 2 + t(bread) / 2),
      symmetric = TRUE,
      only.values = TRUE
    )$values
  } else {
    0
  }
  tolerance <- max(1, max(abs(eigenvalues))) * length(coefficients) * 1e-10
  if (min(eigenvalues) <= tolerance) {
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
