# Analytical delta-method inference for SOP objects.

# Averages whose scenarios set the previous-state variable, or for
# second-order models the state two visits back (`p2_var`), start every
# patient in the set history. Their baseline anchor is then a fixed point
# mass, and only conditional variance (coefficient uncertainty) is defined for
# them.
sets_starting_state <- function(object) {
  if (
    !inherits(
      object,
      c("markov_avg_sops", "markov_avg_state_time", "markov_avg_comparisons")
    )
  ) {
    return(FALSE)
  }
  settings <- markov_result_settings(object)
  state_vars <- c(settings$p_var, settings$p2_var)
  any(state_vars %in% names(attr(object, "avg_args")$variables %||% list()))
}

# Delta inference refuses second-order models, so for them only MVN draws
# are offered.
stop_starting_state_conditional_only <- function(object) {
  advice <- if (is.null(markov_result_settings(object)$p2_var)) {
    'use `method = "delta", vcov = "conditional"` or `method = "mvn"`.'
  } else {
    'use `method = "mvn"`.'
  }
  mostr_abort(
    "Setting the starting state with `variables` supports only conditional ",
    "variance: ",
    advice
  )
}

# Results with one row per prediction profile support only conditional
# variance: there is no average over sampled patients. Returns the subject
# of the error that says so, naming the functions that make such results.
delta_individual_subject <- function(object) {
  if (inherits(object, "markov_comparisons")) {
    return(
      "Delta inference for `comparisons()` and `avg_comparisons(by = FALSE)`"
    )
  }
  if (inherits(object, "markov_state_time")) {
    return(
      "Delta inference for `state_time()` and `avg_state_time(by = FALSE)`"
    )
  }
  if (inherits(object, "markov_sops")) {
    return("`sops()` delta inference")
  }
  NULL
}

stop_individual_conditional_only <- function(object) {
  mostr_abort(
    delta_individual_subject(object),
    ' supports only `vcov = "conditional"` or a coefficient covariance matrix.'
  )
}

delta_resolve_vcov <- function(object, vcov) {
  individual <- !is.null(delta_individual_subject(object))
  choice <- if (is.null(vcov)) {
    if (individual) "conditional" else "unconditional"
  } else if (is.character(vcov)) {
    if (
      length(vcov) != 1L ||
        is.na(vcov) ||
        !vcov %in% c("conditional", "unconditional")
    ) {
      mostr_abort(
        '`vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.'
      )
    }
    vcov
  } else {
    "conditional"
  }
  if (individual && identical(choice, "unconditional")) {
    stop_individual_conditional_only(object)
  }
  if (identical(choice, "unconditional") && sets_starting_state(object)) {
    stop_starting_state_conditional_only(object)
  }
  if (identical(choice, "unconditional") && !is.null(attr(object, "model"))) {
    delta_reject_weighted_unconditional(attr(object, "model"))
  }
  target <- if (identical(choice, "unconditional")) {
    "unconditional"
  } else if (individual) {
    "fixed"
  } else {
    "empirical"
  }
  list(
    target = delta_validate_target(object, target),
    vcov = if (is.character(vcov)) NULL else vcov
  )
}

delta_validate_target <- function(object, target) {
  if (!is.null(delta_individual_subject(object))) {
    if (!identical(target, "fixed")) {
      stop_individual_conditional_only(object)
    }
    return(target)
  }

  if (!target %in% c("empirical", "unconditional")) {
    mostr_abort(
      "Averaged delta inference supports `vcov = \"conditional\"`, ",
      "`vcov = \"unconditional\"`, or a coefficient covariance matrix."
    )
  }
  if (
    identical(target, "unconditional") &&
      isTRUE(attr(object, "newdata_supplied"))
  ) {
    mostr_abort(
      "`vcov = \"unconditional\"` requires the stored fitted-patient ",
      "cohort. For user-supplied `newdata`, use `vcov = \"conditional\"`."
    )
  }
  target
}

delta_validate_sop_scope <- function(object, model, by) {
  if (!is.null(by)) {
    mostr_abort(
      "Analytical delta inference with `by` is not yet supported."
    )
  }
  if (!is.null(markov_result_settings(object)$p2_var)) {
    mostr_abort(
      "Analytical delta inference for second-order Markov models is not yet ",
      "supported."
    )
  }
  model_plan <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  if (!inherits(model_plan, c("orm", "vglm")) || inherits(model_plan, "blrm")) {
    mostr_abort(
      "Analytical delta inference currently supports frequentist `orm`, ",
      "`vglm`, and `robcov_vglm` models."
    )
  }
  invisible(model_plan)
}

delta_compile_plan <- function(object, model_plan, newdata, times) {
  plan <- compile_sop_execution_plan(
    model = model_plan,
    newdata = newdata,
    times = times,
    # Only first-order models reach analytical inference
    # (`delta_validate_sop_scope()`).
    settings = markov_result_settings(object)
  )

  Gamma <- get_effective_coefs(model_plan)
  po <- markov_po_structure(
    Gamma,
    plan$components$col_names,
    plan$components$X_init
  )
  if (is.null(po)) {
    mostr_abort(
      "Analytical delta inference currently requires a full ",
      "proportional-odds model; partial proportional odds are not yet ",
      "supported."
    )
  }
  plan
}

delta_validate_plan_result <- function(result, n_profiles, n_times, n_states) {
  if (
    !is.list(result) ||
      is.null(result$probabilities) ||
      is.null(result$jacobian)
  ) {
    mostr_abort("The analytical SOP engine returned an invalid result.")
  }
  expected_prob <- c(n_profiles, n_times, n_states)
  if (!identical(dim(result$probabilities), expected_prob)) {
    mostr_abort("The analytical SOP probabilities have unexpected dimensions.")
  }
  jac_dim <- dim(result$jacobian)
  if (
    length(jac_dim) != 4L ||
      !identical(jac_dim[seq_len(3L)], expected_prob)
  ) {
    mostr_abort("The analytical SOP Jacobian has unexpected dimensions.")
  }

  coefficient_names <- result$coef_names %||% names(result$coef)
  if (is.null(coefficient_names) || length(coefficient_names) != jac_dim[4L]) {
    mostr_abort("The analytical SOP Jacobian is missing raw coefficient names.")
  }
  if (anyDuplicated(coefficient_names)) {
    mostr_abort("The analytical SOP Jacobian coefficient names must be unique.")
  }
  list(
    probabilities = result$probabilities,
    jacobian = result$jacobian,
    coefficient_names = coefficient_names,
    individual_probabilities = result$individual_probabilities
  )
}

delta_flatten_jacobian <- function(jacobian, coefficient_names) {
  out <- matrix(
    as.numeric(jacobian),
    nrow = prod(dim(jacobian)[seq_len(3L)]),
    ncol = dim(jacobian)[4L]
  )
  colnames(out) <- coefficient_names
  out
}

delta_match_cells <- function(object, calculated, key_cols) {
  object_key <- sop_draw_cell_key(object, key_cols)
  calculated_key <- sop_draw_cell_key(calculated, key_cols)
  index <- match(object_key, calculated_key)
  if (anyNA(index)) {
    mostr_abort(
      "Analytical SOP cells do not align with the original point estimate."
    )
  }
  object_estimate <- object$estimate
  calculated_estimate <- calculated$estimate[index]
  # inferences() recomputes the SOPs with the code that also computes their
  # derivatives. If these differ from the stored estimates (an edited
  # `estimate` column, or a bug in either path), the standard errors would
  # belong to other numbers.
  if (
    !is.numeric(object_estimate) ||
      !is.numeric(calculated_estimate) ||
      length(object_estimate) != length(calculated_estimate) ||
      any(!is.finite(object_estimate)) ||
      any(!is.finite(calculated_estimate)) ||
      estimates_differ(object_estimate, calculated_estimate)
  ) {
    mostr_abort(
      "Analytical SOP replay did not reproduce the stored point estimates."
    )
  }
  index
}

# Whether estimates calculated again differ from the stored ones by more than
# rounding. Two routes to the same estimate agree to about 1e-14 on the test
# grid (measured Nov 2026), so 1e-10 relative leaves ample headroom and still
# catches an edited estimate; the 1e-12 absolute part covers estimates near 0.
estimates_differ <- function(calculated, stored) {
  any(
    abs(calculated - stored) >
      1e-12 + 1e-10 * pmax(abs(calculated), abs(stored))
  )
}

delta_align_coefficient_vcov <- function(result, coefficient_names) {
  V <- result$vcov
  if (!is.matrix(V)) {
    mostr_abort("The analytical coefficient covariance must be a matrix.")
  }
  if (
    is.null(rownames(V)) ||
      is.null(colnames(V)) ||
      anyDuplicated(rownames(V)) ||
      anyDuplicated(colnames(V)) ||
      !setequal(coefficient_names, rownames(V)) ||
      !setequal(coefficient_names, colnames(V))
  ) {
    mostr_abort(
      "The analytical coefficient covariance names must match every raw ",
      "model coefficient."
    )
  }
  V <- V[coefficient_names, coefficient_names, drop = FALSE]
  if (any(!is.finite(V))) {
    mostr_abort(
      "The analytical coefficient covariance must contain finite values."
    )
  }
  V / 2 + t(V) / 2
}

delta_coefficient_state <- function(model, cluster, vcov, jacobian) {
  covariance <- get_delta_cluster_vcov(
    model = model,
    cluster = cluster,
    vcov = vcov
  )
  V <- delta_align_coefficient_vcov(covariance, colnames(jacobian))
  variance <- delta_jacobian_variance(jacobian, V)
  list(
    std.error = sqrt(variance),
    coefficient_vcov = V,
    metadata = covariance$metadata %||% list(source = "coefficient_vcov")
  )
}

delta_jacobian_variance <- function(jacobian, coefficient_vcov) {
  working_bytes <- as.double(nrow(jacobian)) * ncol(jacobian) * 8
  delta_assert_bytes(working_bytes, "Analytical standard-error calculation")
  JV <- jacobian %*% coefficient_vcov
  variance <- rowSums(JV * jacobian)
  tolerance <- sqrt(.Machine$double.eps) * rowSums(abs(JV * jacobian))
  variance[variance < 0 & abs(variance) <= tolerance] <- 0
  if (any(variance < 0)) {
    mostr_abort("Analytical propagation produced a negative variance.")
  }
  variance
}

delta_profile_ids <- function(baseline_data, id_var) {
  if (is.null(id_var) || !id_var %in% names(baseline_data)) {
    mostr_abort(
      "Unconditional delta inference requires the stored starting-profile ",
      "ID column."
    )
  }
  ids <- as.character(baseline_data[[id_var]])
  if (anyNA(ids) || anyDuplicated(ids)) {
    mostr_abort(
      "Unconditional delta inference requires exactly one non-missing ",
      "starting profile per fitted patient."
    )
  }
  ids
}

delta_unconditional_state <- function(
  individual_values,
  average_jacobian,
  profile_ids,
  score_components
) {
  cell_names <- paste0("cell_", seq_len(ncol(individual_values)))
  colnames(individual_values) <- cell_names
  rownames(average_jacobian) <- cell_names
  stacked <- delta_stacked_influence(
    individual_values = individual_values,
    average_jacobian = average_jacobian,
    profile_ids = profile_ids,
    score_components = score_components
  )
  influence <- as.matrix(stacked$influence)
  if (
    nrow(influence) != length(profile_ids) ||
      ncol(influence) != nrow(average_jacobian)
  ) {
    mostr_abort("The stacked SOP influence values have unexpected dimensions.")
  }
  colnames(influence) <- NULL
  standard_error <- stacked$std.error
  metadata <- stacked$metadata %||%
    score_components$metadata %||%
    list()
  metadata$source <- metadata$source %||% "stacked_influence"
  list(
    std.error = as.numeric(standard_error),
    influence = influence,
    average_jacobian = average_jacobian,
    metadata = metadata
  )
}

delta_interval_bounds <- function(
  estimate,
  standard_error,
  conf_level,
  conf_type
) {
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (identical(conf_type, "wald")) {
    return(list(
      conf.low = estimate - critical * standard_error,
      conf.high = estimate + critical * standard_error
    ))
  }
  if (!identical(conf_type, "logit")) {
    mostr_abort(
      "Delta SOP inference supports `conf_type = \"logit\"` or `\"wald\"`."
    )
  }

  lower <- upper <- rep(NA_real_, length(estimate))
  interior <- is.finite(estimate) & estimate > 0 & estimate < 1
  if (any(interior)) {
    logit_se <- standard_error[interior] /
      (estimate[interior] * (1 - estimate[interior]))
    center <- stats::qlogis(estimate[interior])
    lower[interior] <- stats::plogis(center - critical * logit_se)
    upper[interior] <- stats::plogis(center + critical * logit_se)
  }

  # The recursion starts with fitted ordinal probabilities and supplies no
  # structural-boundary metadata for the resulting occupancy probabilities.
  # Saturation can make both its probability and derivative numerically exact.
  undefined <- !interior
  if (any(undefined)) {
    warning(
      "Logit-delta limits are undefined for nonstructural boundary SOP ",
      "estimates; their confidence limits were set to NA.",
      call. = FALSE
    )
  }
  list(conf.low = lower, conf.high = upper)
}

delta_copy_result_attrs <- function(result, object) {
  for (name in names(attributes(object))) {
    if (!name %in% c("names", "row.names", "class")) {
      attr(result, name) <- attr(object, name)
    }
  }
  class(result) <- class(object)
  result
}

delta_finalize_result <- function(
  object,
  standard_error,
  conf_level,
  conf_type,
  target,
  analytical
) {
  # Check the retained analytical state before copying the result data frame or
  # attaching any of its potentially large attributes.
  analytical <- delta_attach_byte_accounting(analytical)
  bounds <- delta_interval_bounds(
    estimate = object$estimate,
    standard_error = standard_error,
    conf_level = conf_level,
    conf_type = conf_type
  )
  result <- as.data.frame(object)
  result$conf.low <- bounds$conf.low
  result$conf.high <- bounds$conf.high
  result$std.error <- standard_error
  result <- delta_copy_result_attrs(result, object)

  keys <- delta_row_key_frame(result)
  analytical$version <- 1L
  analytical$target <- target
  analytical$interval <- conf_type
  analytical$row_keys <- keys
  analytical$row_key <- delta_row_key(keys)
  if (anyDuplicated(analytical$row_key)) {
    mostr_abort("Analytical result-row keys are not unique.")
  }

  attr(result, "analytical") <- analytical
  attr(result, "method") <- "delta"
  attr(result, "target") <- target
  attr(result, "conf_level") <- conf_level
  attr(result, "conf_type") <- conf_type
  attr(
    result,
    "covariance_source"
  ) <- analytical$covariance_metadata$source %||%
    analytical$representation
  attr(result, "covariance_metadata") <- analytical$covariance_metadata
  result <- order_estimate_columns(result)
  class(result) <- class(object)
  result
}

delta_average_arrays <- function(
  probabilities,
  jacobian,
  grid,
  times,
  y_levels,
  variables,
  n_each,
  individual_probabilities = NULL
) {
  n_cf <- nrow(grid)
  q <- dim(jacobian)[4L]
  if (
    !identical(
      dim(probabilities),
      c(n_cf, length(times), length(y_levels))
    ) ||
      !identical(
        dim(jacobian)[seq_len(3L)],
        c(n_cf, length(times), length(y_levels))
      )
  ) {
    mostr_abort("Grouped analytical SOP arrays have unexpected dimensions.")
  }
  cells <- marginalize_sops_array(
    sops_array = probabilities,
    grid = grid,
    times = times,
    y_levels = y_levels,
    variables = variables,
    n_cf = n_cf,
    n_each = 1L
  )
  average_jacobian <- matrix(NA_real_, nrow = nrow(cells), ncol = q)
  individual_values <- if (is.null(individual_probabilities)) {
    NULL
  } else {
    matrix(NA_real_, nrow = n_each, ncol = nrow(cells))
  }
  cursor <- 1L
  for (scenario in seq_len(n_cf)) {
    scenario_jac <- jacobian[scenario, , , , drop = FALSE]
    n_cells <- length(times) * length(y_levels)
    output_rows <- cursor:(cursor + n_cells - 1L)
    average_jacobian[output_rows, ] <- matrix(
      scenario_jac,
      nrow = n_cells,
      ncol = q
    )
    if (!is.null(individual_values)) {
      rows <- ((scenario - 1L) * n_each + 1L):(scenario * n_each)
      scenario_prob <- individual_probabilities[rows, , , drop = FALSE]
      individual_values[, output_rows] <- matrix(
        scenario_prob,
        nrow = n_each,
        ncol = n_cells
      )
    }
    cursor <- cursor + n_cells
  }
  list(
    cells = cells,
    jacobian = average_jacobian,
    individual_values = individual_values
  )
}

inferences_delta_sops <- function(
  object,
  target,
  vcov,
  cluster,
  conf_level,
  conf_type
) {
  is_average <- inherits(object, "markov_avg_sops")
  model <- attr(object, "model")
  if (is.null(model)) {
    mostr_abort("Model not stored in object. Cannot perform delta inference.")
  }
  avg_args <- attr(object, "avg_args")
  call_args <- attr(object, "call_args")
  by <- if (is_average) avg_args$by else NULL
  model_plan <- delta_validate_sop_scope(object, model, by)
  score_components <- if (identical(target, "unconditional")) {
    get_delta_score_components(
      model = model,
      cluster = cluster
    )
  } else {
    NULL
  }
  newdata <- attr(object, "newdata_pred")
  if (is.null(newdata) || !is.data.frame(newdata)) {
    mostr_abort("Stored SOP prediction data are required for delta inference.")
  }
  times <- if (is_average) avg_args$times else call_args$times
  y_levels <- markov_result_settings(object)$y_levels
  if (is_average) {
    variables <- avg_args$variables
    grid <- create_counterfactual_grid(variables)
    n_cf <- nrow(grid)
    if (n_cf < 1L || nrow(newdata) < 1L || nrow(newdata) %% n_cf != 0L) {
      mostr_abort(
        "Stored prediction data are not aligned with the counterfactual grid."
      )
    }
    n_each <- as.integer(nrow(newdata) / n_cf)
    # Native averaging needs contiguous scenarios with identical profile order.
    for (column in union(names(newdata), names(grid))) {
      expected <- if (column %in% names(grid)) {
        rep(grid[[column]], each = n_each)
      } else {
        newdata[
          rep(seq_len(n_each), times = n_cf),
          column,
          drop = FALSE
        ][[1L]]
      }
      if (!identical(unname(newdata[[column]]), unname(expected))) {
        mostr_abort(
          "Stored prediction data must contain counterfactual grid blocks ",
          "in grid order, with the same starting-profile order in every block."
        )
      }
    }
  } else {
    variables <- grid <- NULL
    n_cf <- n_each <- NULL
  }
  plan <- delta_compile_plan(
    object = object,
    model_plan = model_plan,
    newdata = newdata,
    times = times
  )
  engine_result <- delta_validate_plan_result(
    run_sop_delta_plan(
      plan,
      model_plan,
      average_group_size = if (is_average) n_each else NULL,
      retain_individual_probabilities = identical(target, "unconditional")
    ),
    n_profiles = if (is_average) n_cf else nrow(newdata),
    n_times = length(times),
    n_states = length(y_levels)
  )

  if (!is_average) {
    calculated <- array_to_df_individual(
      engine_result$probabilities,
      times,
      y_levels,
      newdata
    )
    index <- delta_match_cells(object, calculated, c("rowid", "time", "state"))
    jacobian <- delta_flatten_jacobian(
      engine_result$jacobian,
      engine_result$coefficient_names
    )[index, , drop = FALSE]
    coefficient <- delta_coefficient_state(model, cluster, vcov, jacobian)
    analytical <- list(
      representation = "coefficient",
      jacobian = jacobian,
      coefficient_vcov = coefficient$coefficient_vcov,
      coefficient_names = colnames(jacobian),
      covariance_metadata = coefficient$metadata
    )
    return(delta_finalize_result(
      object,
      standard_error = coefficient$std.error,
      conf_level = conf_level,
      conf_type = conf_type,
      target = target,
      analytical = analytical
    ))
  }

  reduced <- delta_average_arrays(
    probabilities = engine_result$probabilities,
    jacobian = engine_result$jacobian,
    grid = grid,
    times = times,
    y_levels = y_levels,
    variables = variables,
    n_each = n_each,
    individual_probabilities = engine_result$individual_probabilities
  )
  key_cols <- c("time", "state", names(variables))
  index <- delta_match_cells(object, reduced$cells, key_cols)
  average_jacobian <- reduced$jacobian[index, , drop = FALSE]
  colnames(average_jacobian) <- engine_result$coefficient_names

  if (identical(target, "unconditional")) {
    if (is.null(reduced$individual_values)) {
      mostr_abort(
        "Unconditional analytical SOP inference lost individual ",
        "probabilities."
      )
    }
    individual_values <- reduced$individual_values[, index, drop = FALSE]
    baseline_data <- newdata[seq_len(n_each), , drop = FALSE]
    id_var <- markov_result_settings(object)$id_var
    profile_ids <- delta_profile_ids(baseline_data, id_var)
    unconditional <- delta_unconditional_state(
      individual_values = individual_values,
      average_jacobian = average_jacobian,
      profile_ids = profile_ids,
      score_components = score_components
    )
    analytical <- list(
      representation = "influence",
      average_jacobian = unconditional$average_jacobian,
      influence = unconditional$influence,
      profile_ids = profile_ids,
      coefficient_names = colnames(average_jacobian),
      covariance_metadata = unconditional$metadata
    )
    standard_error <- unconditional$std.error
  } else {
    coefficient <- delta_coefficient_state(
      model,
      cluster,
      vcov,
      average_jacobian
    )
    analytical <- list(
      representation = "coefficient",
      jacobian = average_jacobian,
      coefficient_vcov = coefficient$coefficient_vcov,
      coefficient_names = colnames(average_jacobian),
      covariance_metadata = coefficient$metadata
    )
    standard_error <- coefficient$std.error
  }

  delta_finalize_result(
    object,
    standard_error = standard_error,
    conf_level = conf_level,
    conf_type = conf_type,
    target = target,
    analytical = analytical
  )
}
