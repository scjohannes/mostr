# Patient-cluster covariance and unconditional delta-method helpers.

delta_validate_full_po_model <- function(model) {
  fit <- sop_delta_backend_model(model)

  if (inherits(fit, "blrm")) {
    mostr_abort(
      "Analytical unconditional inference supports only frequentist full ",
      "proportional-odds orm and vglm models; blrm models are not supported."
    )
  }
  if (!inherits(fit, c("orm", "vglm"))) {
    mostr_abort(
      "Analytical unconditional inference supports only full proportional-odds ",
      "orm and vglm models."
    )
  }

  if (inherits(fit, "vglm")) {
    constraints <- VGAM::constraints(fit)
    constraint_names <- names(constraints)
    slope_index <- if (is.null(constraint_names)) {
      seq_along(constraints)[-1L]
    } else {
      which(constraint_names != "(Intercept)")
    }
    is_common_slope <- vapply(
      constraints[slope_index],
      function(constraint) {
        constraint <- as.matrix(constraint)
        ncol(constraint) == 1L &&
          all(is.finite(constraint)) &&
          max(abs(constraint[, 1L] - constraint[1L, 1L])) <= 1e-12
      },
      logical(1)
    )
    if (length(is_common_slope) && !all(is_common_slope)) {
      mostr_abort(
        "Analytical inference currently requires a full proportional-odds ",
        "vglm model (`parallel = TRUE`); partial or nonproportional-odds ",
        "constraints are not supported."
      )
    }
  }

  invisible(fit)
}

delta_validate_named_matrix <- function(
  value,
  names,
  label,
  positive_semidefinite = FALSE,
  positive_definite = FALSE
) {
  if (methods::is(value, "Matrix")) {
    value <- as.matrix(value)
  }
  q <- length(names)
  if (
    !is.matrix(value) || !is.numeric(value) || !identical(dim(value), c(q, q))
  ) {
    mostr_abort(
      "`",
      label,
      "` must be a numeric [",
      q,
      " x ",
      q,
      "] matrix with one row and column per raw model coefficient."
    )
  }
  row_names <- rownames(value)
  column_names <- colnames(value)
  valid_names <- function(x) {
    !is.null(x) &&
      length(x) == q &&
      !anyNA(x) &&
      !any(x == "") &&
      !anyDuplicated(x)
  }
  if (
    !valid_names(row_names) ||
      !valid_names(column_names) ||
      !setequal(row_names, names) ||
      !setequal(column_names, names)
  ) {
    mostr_abort(
      "Rows and columns of `",
      label,
      "` must be uniquely named and match the complete raw coefficient vector."
    )
  }
  value <- value[names, names, drop = FALSE]
  if (any(!is.finite(value))) {
    mostr_abort("`", label, "` contains non-finite values.")
  }

  scale <- max(abs(value))
  symmetry_error <- max(abs(value - t(value)))
  if (!is.finite(symmetry_error) || symmetry_error > 1e-8 * scale) {
    mostr_abort("`", label, "` is not numerically symmetric.")
  }
  value <- value / 2 + t(value) / 2

  if (positive_semidefinite || positive_definite) {
    eigenvalues <- eigen(value, symmetric = TRUE, only.values = TRUE)$values
    eigen_scale <- max(abs(eigenvalues))
    semidefinite_tolerance <- sqrt(.Machine$double.eps) * eigen_scale
    definite_tolerance <- 100 * .Machine$double.eps * q * eigen_scale
    if (
      positive_semidefinite &&
        min(eigenvalues) < -semidefinite_tolerance
    ) {
      mostr_abort(
        "`",
        label,
        "` is not numerically positive semidefinite."
      )
    }
    if (positive_definite && min(eigenvalues) <= definite_tolerance) {
      mostr_abort(
        "`",
        label,
        "` is not numerically positive definite."
      )
    }
  }

  value
}

delta_fitting_row_count <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    if (!is.matrix(model$scores)) {
      mostr_abort(
        "The `robcov_vglm` object does not contain row-aligned score ",
        "contributions. Refit it with `robcov_vglm()`."
      )
    }
    return(nrow(model$scores))
  }

  fit <- sop_delta_backend_model(model)
  if (inherits(fit, "vglm")) {
    return(as.integer(stats::nobs(fit, type = "lm")))
  }
  if (inherits(fit, "orm")) {
    if (!is.null(fit$x)) {
      return(nrow(fit$x))
    }
    if (!is.null(fit$y)) {
      return(length(fit$y))
    }
    n <- tryCatch(stats::nobs(fit), error = function(e) NULL)
    if (!is.null(n) && length(n) == 1L && is.finite(n)) {
      return(as.integer(n))
    }
    mostr_abort(
      "Could not determine the fitted orm row count. Refit with ",
      "`x = TRUE, y = TRUE`."
    )
  }

  mostr_abort(
    "Analytical unconditional inference supports only full proportional-odds ",
    "orm and vglm models."
  )
}

delta_align_explicit_cluster <- function(model, cluster, n_rows) {
  fit <- sop_delta_backend_model(model)
  if (inherits(fit, "vglm")) {
    return(align_cluster_vglm(fit, cluster, n_rows))
  }
  align_cluster_orm(fit, cluster, n_rows)
}

delta_validate_cluster <- function(cluster, n_rows) {
  if (length(cluster) != n_rows) {
    mostr_abort(
      "The patient cluster vector has length ",
      length(cluster),
      ", but the fitted model has ",
      n_rows,
      " score rows."
    )
  }
  if (is.list(cluster) || is.data.frame(cluster) || is.matrix(cluster)) {
    mostr_abort(
      "`cluster` must be a one-dimensional patient ID vector."
    )
  }
  if (anyNA(cluster)) {
    mostr_abort(
      "The patient cluster vector contains missing values after alignment ",
      "with the fitted rows."
    )
  }
  cluster_character <- as.character(cluster)
  if (anyNA(cluster_character) || any(cluster_character == "")) {
    mostr_abort(
      "Patient cluster IDs must be non-missing and non-empty."
    )
  }
  ids <- unique(cluster_character)
  if (length(ids) < 2L) {
    mostr_abort(
      "Patient-cluster inference requires at least two independent patients."
    )
  }
  list(cluster = cluster, ids = ids)
}

# Resolve the patient of each fitting row for analytical inference with
# conditional or unconditional variance. Explicit row-level IDs take
# precedence; fitting rows are never implicit clusters. The data stored by
# `mostm()` are the fitted rows in the order of the score rows, so an ID column
# is read from them as it is. The data frame the user fitted to is not looked
# up again: it may have been reordered, filtered or removed since the fit.
resolve_delta_cluster <- function(model, cluster = NULL) {
  delta_validate_full_po_model(model)
  n_rows <- delta_fitting_row_count(model)
  source <- NULL
  id_var <- NULL

  if (!is.null(cluster)) {
    if (inherits(cluster, "formula")) {
      cluster_vars <- all.vars(cluster)
      if (length(cluster_vars) != 1L) {
        mostr_abort(
          "A `cluster` formula must select exactly one patient ID column, ",
          "for example `~ id`."
        )
      }
      data <- markov_model_data(model)
      id_var <- cluster_vars[[1L]]
      if (is.null(data) || !is.data.frame(data) || !id_var %in% names(data)) {
        mostr_abort(
          "Explicit cluster formula `~ ",
          id_var,
          "` requires that column in the model's stored fitting data."
        )
      }
      cluster <- data[[id_var]]
      source <- "explicit_formula"
    } else {
      cluster <- delta_align_explicit_cluster(model, cluster, n_rows)
      source <- "explicit_vector"
    }
  } else {
    id_var <- markov_model_id_var(model)
    data <- markov_model_data(model)
    if (is.null(id_var) || !is.data.frame(data)) {
      mostr_abort(
        "Analytical standard errors need the patient of each fitting row, ",
        "but the model does not record its fitting data and patient-ID ",
        "column. Refit it with `mostm(..., id_var = ...)`, or give the ",
        "patient of each fitting row in `cluster`. Rows are never treated as ",
        "separate patients."
      )
    }
    if (!id_var %in% names(data)) {
      mostr_abort(
        "Stored patient ID column `",
        id_var,
        "` is absent from the stored fitting data."
      )
    }
    cluster <- data[[id_var]]
    source <- "stored_id_var"
  }

  checked <- delta_validate_cluster(cluster, n_rows)
  list(
    cluster = checked$cluster,
    ids = checked$ids,
    metadata = list(
      source = source,
      id_var = id_var,
      cluster_level = "patient",
      n_rows = n_rows,
      n_clusters = length(checked$ids)
    )
  )
}

delta_vglm_covariance_metadata <- function(robust, source) {
  list(
    source = source,
    covariance = "patient_cluster_robust",
    backend = "vglm",
    type = robust$type,
    cadjust = robust$cadjust,
    adjustment_factor = robust$adjustment_factor,
    bread_type = robust$bread_type,
    n_clusters = robust$n_clusters
  )
}

delta_orm_backend_aliases <- function(model, coefficient_names) {
  threshold_count <- as.integer(model$non.slopes %||% 0L)
  slope_count <- length(coefficient_names) - threshold_count
  if (
    threshold_count < 1L ||
      slope_count < 0L ||
      threshold_count + slope_count != length(coefficient_names)
  ) {
    return(list(coefficient_names))
  }

  threshold_names <- coefficient_names[seq_len(threshold_count)]
  design <- model$Design
  slope_aliases <- list(
    design$colnames,
    design$mmcolnames,
    colnames(model$x)
  )
  aliases <- c(
    list(coefficient_names),
    lapply(slope_aliases, function(value) {
      if (is.null(value) || length(value) != slope_count) {
        return(NULL)
      }
      c(threshold_names, as.character(value))
    })
  )
  Filter(
    function(value) {
      length(value) == length(coefficient_names) &&
        !anyNA(value) &&
        !any(value == "") &&
        !anyDuplicated(value)
    },
    aliases
  )
}

# ORM covariance matrices are constructed in coefficient order, but penalized
# fits can lose dimnames and inline rms transforms can use model-matrix labels
# instead of the display labels in coef(). Normalize only matrices produced
# internally from this exact orm fit.
delta_normalize_orm_backend_matrix <- function(
  value,
  model,
  coefficient_names
) {
  if (methods::is(value, "Matrix")) {
    value <- as.matrix(value)
  }
  q <- length(coefficient_names)
  if (!is.matrix(value) || !identical(dim(value), c(q, q))) {
    return(value)
  }

  aliases <- delta_orm_backend_aliases(model, coefficient_names)
  normalize_axis <- function(axis_names) {
    if (is.null(axis_names)) {
      return(coefficient_names)
    }
    for (candidate in aliases) {
      index <- match(axis_names, candidate)
      if (!anyNA(index) && !anyDuplicated(index)) {
        return(coefficient_names[index])
      }
    }
    axis_names
  }

  rownames(value) <- normalize_axis(rownames(value))
  colnames(value) <- normalize_axis(colnames(value))
  value
}

# The case weights the model was fitted with (`mostm(weights = )`), or NULL.
# vglm fits store them as a one-column matrix, orm fits as a vector.
delta_fit_case_weights <- function(model) {
  fit <- sop_delta_backend_model(model)
  weights <- if (inherits(fit, "vglm")) {
    tryCatch(methods::slot(fit, "prior.weights"), error = function(e) NULL)
  } else {
    tryCatch(fit$weights, error = function(e) NULL)
  }
  if (is.null(weights) || length(weights) == 0L) {
    return(NULL)
  }
  as.numeric(weights)
}

# Unconditional variance has not been checked with case weights for any
# engine, so it is refused whenever a fitted row has a weight other than 1.
# This also covers the default `vcov = NULL` for averages.
delta_reject_weighted_unconditional <- function(model) {
  weights <- delta_fit_case_weights(model)
  if (is.null(weights) || all(is.finite(weights) & weights == 1)) {
    return(invisible(NULL))
  }

  mostr_abort(
    "Models fitted with case weights (`weights` in `mostm()`) do not support ",
    "unconditional variance, because its calculation has not been checked ",
    "with weights. Unconditional variance is the default for averages; pass ",
    "`vcov = \"conditional\"` to `inferences()` for conditional variance, ",
    "which accounts for coefficient estimation and treats the patients' ",
    "starting states and covariates as given."
  )
}

delta_reject_penalized_orm_unconditional <- function(model) {
  penalty_matrix <- model$penalty.matrix
  penalty <- unlist(model$penalty, use.names = FALSE)
  penalized <-
    (!is.null(penalty_matrix) &&
      length(penalty_matrix) > 0L &&
      any(is.finite(penalty_matrix) & penalty_matrix != 0)) ||
    (!is.null(penalty) &&
      length(penalty) > 0L &&
      any(is.finite(penalty) & penalty != 0))
  if (!penalized) {
    return(invisible(NULL))
  }

  mostr_abort(
    "Unconditional delta inference does not currently support penalized ",
    "orm likelihoods because the patient scores and sensitivity have not ",
    "been shown to incorporate the penalty consistently. Ordinary spline ",
    "transformations remain supported with `vcov = \"conditional\"`."
  )
}

# Resolve a complete raw-coefficient covariance, preferring an explicit matrix
# and otherwise using the backend's patient-cluster robust convention.
get_delta_cluster_vcov <- function(model, cluster = NULL, vcov = NULL) {
  delta_validate_full_po_model(model)
  beta <- sop_delta_raw_coef(model)

  if (!is.null(vcov)) {
    covariance <- delta_validate_named_matrix(
      vcov,
      names(beta),
      "vcov",
      positive_semidefinite = TRUE
    )
    return(list(
      vcov = covariance,
      metadata = list(
        source = "explicit",
        covariance = "user_supplied",
        backend = class(sop_delta_backend_model(model))[[1L]],
        type = NA_character_,
        cadjust = NA,
        adjustment_factor = NA_real_,
        bread_type = NA_character_,
        n_clusters = NA_integer_
      ),
      cluster_info = NULL
    ))
  }

  cluster_info <- resolve_delta_cluster(model, cluster = cluster)
  fit <- sop_delta_backend_model(model)

  if (inherits(fit, "vglm")) {
    stored_cluster_matches <- inherits(model, "robcov_vglm") &&
      !is.null(model$cluster) &&
      identical(
        as.character(model$cluster),
        as.character(cluster_info$cluster)
      )
    use_stored <- inherits(model, "robcov_vglm") &&
      is.null(cluster) &&
      !is.null(model$var) &&
      stored_cluster_matches
    robust <- if (use_stored) {
      model
    } else {
      robcov_vglm(
        fit,
        cluster = cluster_info$cluster,
        bread = if (inherits(model, "robcov_vglm")) {
          model$bread_type %||% "observed"
        } else {
          "observed"
        },
        type = if (inherits(model, "robcov_vglm")) {
          model$type %||% "HC0"
        } else {
          "HC0"
        },
        cadjust = if (inherits(model, "robcov_vglm")) {
          model$cadjust %||% TRUE
        } else {
          TRUE
        }
      )
    }
    covariance <- delta_validate_named_matrix(
      robust$var,
      names(beta),
      "patient-cluster robust covariance",
      positive_semidefinite = TRUE
    )
    metadata <- delta_vglm_covariance_metadata(
      robust,
      if (use_stored) "stored_robcov_vglm" else "computed_robcov_vglm"
    )
  } else {
    stored_metadata <- orm_robust_metadata(model)
    use_stored <- is.null(cluster) &&
      orm_stored_covariance_valid(model, cluster_info$cluster)
    robust <- if (use_stored) {
      model
    } else {
      robcov_orm(
        fit,
        cluster = cluster_info$cluster,
        type = stored_metadata$type %||% "HC0",
        cadjust = stored_metadata$cadjust %||% TRUE
      )
    }
    backend_covariance <- delta_normalize_orm_backend_matrix(
      robust$var,
      fit,
      names(beta)
    )
    covariance <- delta_validate_named_matrix(
      backend_covariance,
      names(beta),
      "patient-cluster robust covariance",
      positive_semidefinite = TRUE
    )
    metadata <- list(
      source = if (use_stored) {
        "stored_markov_robcov_orm"
      } else {
        "computed_markov_robcov_orm"
      },
      covariance = "patient_cluster_robust",
      backend = "orm",
      type = orm_robust_metadata(robust)$type,
      cadjust = orm_robust_metadata(robust)$cadjust,
      adjustment_factor = orm_robust_metadata(robust)$adjustment_factor,
      bread_type = orm_robust_metadata(robust)$bread_convention,
      n_clusters = orm_robust_metadata(robust)$n_clusters
    )
  }

  metadata$cluster_source <- cluster_info$metadata$source
  metadata$cluster_level <- "patient"
  list(
    vcov = covariance,
    metadata = metadata,
    cluster_info = cluster_info
  )
}

delta_validate_scores <- function(scores, coefficient_names, n_rows) {
  scores <- as.matrix(scores)
  if (
    !is.numeric(scores) ||
      !identical(dim(scores), c(n_rows, length(coefficient_names)))
  ) {
    mostr_abort(
      "Row score contributions must be a numeric [",
      n_rows,
      " x ",
      length(coefficient_names),
      "] matrix."
    )
  }
  score_names <- colnames(scores)
  if (
    is.null(score_names) ||
      anyNA(score_names) ||
      any(score_names == "") ||
      anyDuplicated(score_names) ||
      !setequal(score_names, coefficient_names)
  ) {
    mostr_abort(
      "Score columns must uniquely match the complete raw coefficient vector."
    )
  }
  scores <- scores[, coefficient_names, drop = FALSE]
  if (any(!is.finite(scores))) {
    mostr_abort("Row score contributions contain non-finite values.")
  }
  scores
}

# Aggregate raw likelihood scores by patient and construct inverse
# per-participant sensitivity for the stacked influence function.
get_delta_score_components <- function(model, cluster = NULL) {
  delta_validate_full_po_model(model)
  beta <- sop_delta_raw_coef(model)
  fit <- sop_delta_backend_model(model)
  delta_reject_weighted_unconditional(model)
  if (inherits(fit, "orm")) {
    delta_reject_penalized_orm_unconditional(fit)
  }
  cluster_info <- resolve_delta_cluster(model, cluster = cluster)
  n_rows <- cluster_info$metadata$n_rows

  if (inherits(model, "robcov_vglm")) {
    row_scores <- model$scores
    bread <- model$bread
    component_metadata <- list(
      bread_source = paste0("robcov_vglm_", model$bread_type %||% "observed"),
      backend_hc_type_ignored = model$type %||% NA_character_,
      backend_cadjust_ignored = model$cadjust %||% NA
    )
  } else if (inherits(fit, "vglm")) {
    row_scores <- compute_scores_vglm(fit)
    bread <- compute_observed_bread_vglm(fit)$bread
    component_metadata <- list(
      bread_source = "vglm_observed_information",
      backend_hc_type_ignored = NA_character_,
      backend_cadjust_ignored = NA
    )
  } else {
    row_scores <- compute_scores_orm(fit)
    bread <- delta_normalize_orm_backend_matrix(
      orm_model_bread(fit)$bread,
      fit,
      names(beta)
    )
    component_metadata <- list(
      bread_source = "rms_inverse_total_information",
      backend_hc_type_ignored = NA_character_,
      backend_cadjust_ignored = NA
    )
  }

  row_scores <- delta_validate_scores(row_scores, names(beta), n_rows)
  bread <- delta_validate_named_matrix(
    bread,
    names(beta),
    "model-based bread",
    positive_definite = TRUE
  )

  cluster_character <- as.character(cluster_info$cluster)
  cluster_factor <- factor(cluster_character, levels = cluster_info$ids)
  scores <- rowsum(row_scores, cluster_factor, reorder = FALSE)
  scores <- as.matrix(scores)
  rownames(scores) <- cluster_info$ids
  colnames(scores) <- names(beta)

  n <- nrow(scores)
  Ainv <- n * bread
  dimnames(Ainv) <- list(names(beta), names(beta))

  component_metadata$cluster_source <- cluster_info$metadata$source
  component_metadata$n_patients <- n
  list(
    ids = cluster_info$ids,
    scores = scores,
    Ainv = Ainv,
    n = n,
    metadata = component_metadata
  )
}

delta_validate_profile_matrix <- function(value, label) {
  value <- as.matrix(value)
  if (!is.numeric(value) || length(value) == 0L || any(!is.finite(value))) {
    mostr_abort(
      "`",
      label,
      "` must be a non-empty finite numeric matrix."
    )
  }
  value
}

# Combine centered profile functionals and fitted-coefficient influence for the
# fitted-cohort unconditional target.
delta_stacked_influence <- function(
  individual_values,
  average_jacobian,
  profile_ids,
  score_components
) {
  values <- delta_validate_profile_matrix(
    individual_values,
    "individual_values"
  )
  G <- delta_validate_profile_matrix(average_jacobian, "average_jacobian")
  n <- nrow(values)
  m <- ncol(values)

  # The variance below is a sample covariance over patients, which needs at
  # least two of them.
  if (n < 2L) {
    mostr_abort(
      "Stacked unconditional inference requires at least two profiles."
    )
  }
  if (length(profile_ids) != n || anyNA(profile_ids)) {
    mostr_abort(
      "`profile_ids` must contain one non-missing patient ID per profile row."
    )
  }
  profile_ids <- as.character(profile_ids)
  if (any(profile_ids == "") || anyDuplicated(profile_ids)) {
    mostr_abort("`profile_ids` must be non-empty and unique.")
  }
  if (!is.list(score_components)) {
    mostr_abort(
      "`score_components` must come from `get_delta_score_components()`."
    )
  }
  required <- c("ids", "scores", "Ainv", "n")
  missing_required <- required[vapply(
    required,
    function(name) is.null(score_components[[name]]),
    logical(1)
  )]
  if (length(missing_required)) {
    mostr_abort(
      "`score_components` is missing: ",
      paste(missing_required, collapse = ", "),
      "."
    )
  }

  score_ids <- as.character(score_components$ids)
  scores <- as.matrix(score_components$scores)
  if (
    anyNA(score_ids) ||
      any(score_ids == "") ||
      anyDuplicated(score_ids) ||
      nrow(scores) != length(score_ids)
  ) {
    mostr_abort(
      "Score-component IDs must be non-missing, unique, and row-aligned."
    )
  }
  if (
    length(score_components$n) != 1L ||
      !is.numeric(score_components$n) ||
      is.na(score_components$n) ||
      score_components$n != length(score_ids)
  ) {
    mostr_abort("`score_components$n` does not match its patient score rows.")
  }

  unmatched_scores <- setdiff(score_ids, profile_ids)
  unmatched_profiles <- setdiff(profile_ids, score_ids)
  if (length(unmatched_scores)) {
    mostr_abort(
      "Likelihood-score patients are absent from the target profiles: ",
      paste(utils::head(unmatched_scores, 5L), collapse = ", "),
      if (length(unmatched_scores) > 5L) " ..." else "",
      ". Fitted-cohort unconditional inference requires every score patient to ",
      "appear exactly once in `profile_ids`."
    )
  }
  if (length(unmatched_profiles)) {
    mostr_abort(
      "Target-profile patients have no usable likelihood transition: ",
      paste(utils::head(unmatched_profiles, 5L), collapse = ", "),
      if (length(unmatched_profiles) > 5L) " ..." else "",
      ". Unconditional inference includes only fitted patients and does not ",
      "insert zero patient-score rows."
    )
  }
  if (score_components$n != n) {
    mostr_abort(
      "The starting-profile and patient-score cohorts have different sizes."
    )
  }

  cell_names <- colnames(values)
  if (
    is.null(cell_names) ||
      anyNA(cell_names) ||
      any(cell_names == "") ||
      anyDuplicated(cell_names)
  ) {
    mostr_abort(
      "`individual_values` columns must have unique, non-missing cell names."
    )
  }
  jacobian_rows <- rownames(G)
  if (
    is.null(jacobian_rows) ||
      anyNA(jacobian_rows) ||
      any(jacobian_rows == "") ||
      anyDuplicated(jacobian_rows) ||
      !setequal(jacobian_rows, cell_names)
  ) {
    mostr_abort(
      "Jacobian row names must uniquely match all `individual_values` cells."
    )
  }
  if (nrow(G) != m) {
    mostr_abort("`average_jacobian` must have one row per estimand cell.")
  }
  G <- G[cell_names, , drop = FALSE]

  coefficient_names <- colnames(scores)
  if (
    is.null(coefficient_names) ||
      anyNA(coefficient_names) ||
      any(coefficient_names == "") ||
      anyDuplicated(coefficient_names)
  ) {
    mostr_abort(
      "Score columns must have unique, non-missing coefficient names."
    )
  }
  jacobian_columns <- colnames(G)
  if (
    is.null(jacobian_columns) ||
      anyNA(jacobian_columns) ||
      any(jacobian_columns == "") ||
      anyDuplicated(jacobian_columns) ||
      !setequal(jacobian_columns, coefficient_names)
  ) {
    mostr_abort("Score and Jacobian coefficient names do not match.")
  }
  G <- G[, coefficient_names, drop = FALSE]
  score_Ainv <- delta_validate_named_matrix(
    score_components$Ainv,
    coefficient_names,
    "score_components$Ainv",
    positive_definite = TRUE
  )
  if (any(!is.finite(scores))) {
    mostr_abort("Patient score contributions contain non-finite values.")
  }

  score_order <- match(profile_ids, score_ids)
  aligned_scores <- scores[
    score_order,
    coefficient_names,
    drop = FALSE
  ]
  rownames(aligned_scores) <- profile_ids
  Ainv <- score_Ainv

  estimate <- colMeans(values)
  profile_component <- sweep(values, 2L, estimate, "-")
  coefficient_influence <- aligned_scores %*% t(Ainv) %*% t(G)
  influence <- profile_component + coefficient_influence
  dimnames(influence) <- list(profile_ids, cell_names)

  # Each estimate is an average over the n patients, so its variance is the
  # variance of the patients' influence values divided by n. The sample variance
  # divides by n - 1, as cov() does, which gives the same n / (n - 1) factor as
  # the default conditional variance (HC0 with `cadjust = TRUE`);
  # test-sops-delta-unconditional-identity.R relies on this. With
  # `cadjust = FALSE` only the conditional variance loses the factor. Centring
  # at the column means is harmless: summed scores are 0 at the fit (vglm stops
  # early, but that moves standard errors by less than 1e-11).
  #
  # Only the column variances are computed here. The full covariance has one
  # row and column per cell, so vcov() builds it on request from `influence`,
  # under the `mostr.delta_max_bytes` limit.
  centered_influence <- sweep(influence, 2L, colMeans(influence), "-")
  variance <- colSums(centered_influence^2) / (n - 1) / n
  names(variance) <- cell_names

  list(
    influence = influence,
    std.error = sqrt(pmax(variance, 0)),
    metadata = c(
      list(
        source = "stacked_patient_influence",
        covariance_source = "stacked_patient_influence",
        target = "unconditional",
        cohort = "fitted_patients_with_starting_profiles",
        covariance_convention = "sample_covariance_of_influence_divided_by_n",
        finite_sample_method = "patient_level_sample_covariance",
        finite_sample_factor = n / (n - 1),
        n_patients = n,
        n_score_clusters = score_components$n
      ),
      score_components$metadata %||% list()
    )
  )
}
