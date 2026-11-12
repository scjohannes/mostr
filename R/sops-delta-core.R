# Analytic first-order derivatives for compiled SOP execution plans.

delta_max_bytes <- function() {
  value <- getOption("mostr.delta_max_bytes", 256 * 1024^2)
  if (
    length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0
  ) {
    mostr_abort(
      "Option `mostr.delta_max_bytes` must be a positive finite number."
    )
  }
  as.double(value)
}

stop_sop_delta_too_large <- function(required_bytes, max_bytes, purpose) {
  message <- paste0(
    purpose,
    " requires ",
    format(required_bytes, scientific = FALSE, trim = TRUE),
    " bytes, above the analytical delta-method limit of ",
    format(max_bytes, scientific = FALSE, trim = TRUE),
    " bytes. Increase option `mostr.delta_max_bytes` only if the ",
    "required allocation is acceptable."
  )
  mostr_abort(
    message,
    class = "mostr_delta_too_large",
    data = list(
      required_bytes = as.double(required_bytes),
      max_bytes = as.double(max_bytes),
      purpose = purpose
    )
  )
}

delta_assert_bytes <- function(
  bytes,
  purpose = "Analytical delta calculation"
) {
  bytes <- as.double(bytes)
  max_bytes <- delta_max_bytes()
  if (length(bytes) != 1L || is.na(bytes) || !is.finite(bytes) || bytes < 0) {
    mostr_abort(
      "Analytical allocation size must be a finite non-negative number."
    )
  }
  if (bytes > max_bytes) {
    stop_sop_delta_too_large(bytes, max_bytes, purpose)
  }
  invisible(bytes)
}

sop_delta_backend_model <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    if (is.null(model$vglm_fit) || !inherits(model$vglm_fit, "vglm")) {
      mostr_abort(
        "A `robcov_vglm` model must retain its underlying `vglm_fit`."
      )
    }
    return(model$vglm_fit)
  }
  model
}

sop_delta_raw_coef <- function(model) {
  coef <- get_coef(model)
  if (
    !is.numeric(coef) ||
      length(coef) == 0L ||
      anyNA(coef) ||
      any(!is.finite(coef))
  ) {
    mostr_abort(
      "Analytical SOP inference requires finite numeric coefficients."
    )
  }
  if (
    is.null(names(coef)) ||
      anyNA(names(coef)) ||
      any(names(coef) == "") ||
      anyDuplicated(names(coef))
  ) {
    mostr_abort(
      "Analytical SOP inference requires uniquely named raw coefficients."
    )
  }
  coef
}

# Return A such that as.vector(Gamma) = A %*% coef(model). The vectorization is
# column-major, matching base R matrices and the execution-plan coefficient
# layout.
get_effective_coef_map <- function(model) {
  coef <- sop_delta_raw_coef(model)
  backend <- sop_delta_backend_model(model)

  if (inherits(backend, "vglm")) {
    constraints <- VGAM::constraints(backend)
    if (
      !is.list(constraints) ||
        length(constraints) == 0L ||
        is.null(names(constraints)) ||
        any(names(constraints) == "")
    ) {
      mostr_abort(
        "Could not determine the `vglm` coefficient constraints."
      )
    }
    threshold_count <- unique(vapply(constraints, nrow, integer(1)))
    if (length(threshold_count) != 1L || threshold_count < 1L) {
      mostr_abort(
        "The `vglm` constraint matrices have incompatible dimensions."
      )
    }
    M <- unname(threshold_count)
    P <- length(constraints)
    map <- matrix(
      0,
      nrow = M * P,
      ncol = length(coef),
      dimnames = list(NULL, names(coef))
    )
    coefficient_index <- 1L
    for (term_index in seq_along(constraints)) {
      constraint <- constraints[[term_index]]
      if (!is.numeric(constraint) || nrow(constraint) != M) {
        mostr_abort(
          "The `vglm` constraint matrices must be numeric matrices."
        )
      }
      coefficient_count <- ncol(constraint)
      if (coefficient_count < 1L) {
        mostr_abort(
          "A `vglm` constraint matrix has no coefficient columns."
        )
      }
      coefficient_rows <- coefficient_index + seq_len(coefficient_count) - 1L
      if (max(coefficient_rows) > length(coef)) {
        mostr_abort(
          "The `vglm` constraints require more coefficients than the model stores."
        )
      }
      gamma_rows <- (term_index - 1L) * M + seq_len(M)
      map[gamma_rows, coefficient_rows] <- constraint
      coefficient_index <- coefficient_index + coefficient_count
    }
    if (coefficient_index - 1L != length(coef)) {
      mostr_abort(
        "The `vglm` constraints do not consume every raw model coefficient."
      )
    }
    gamma_dimnames <- list(
      paste0("eta", seq_len(M)),
      names(constraints)
    )
  } else if (inherits(backend, "orm")) {
    M <- backend$non.slopes
    if (
      length(M) != 1L ||
        !is.numeric(M) ||
        is.na(M) ||
        M < 1L ||
        M != as.integer(M)
    ) {
      mostr_abort(
        "Cannot determine the number of `orm` threshold coefficients."
      )
    }
    M <- as.integer(M)
    if (length(coef) < M) {
      mostr_abort(
        "The `orm` coefficient vector is shorter than its thresholds."
      )
    }
    slope_count <- length(coef) - M
    P <- slope_count + 1L
    map <- matrix(
      0,
      nrow = M * P,
      ncol = length(coef),
      dimnames = list(NULL, names(coef))
    )
    map[seq_len(M), seq_len(M)] <- diag(M)
    if (slope_count > 0L) {
      for (slope_index in seq_len(slope_count)) {
        gamma_rows <- slope_index * M + seq_len(M)
        map[gamma_rows, M + slope_index] <- 1
      }
    }
    gamma_dimnames <- list(
      paste0("eta", seq_len(M)),
      c("(Intercept)", names(coef)[M + seq_len(slope_count)])
    )
  } else {
    mostr_abort(
      "Analytical SOP inference supports only `vglm`, `robcov_vglm`, and `orm` models."
    )
  }

  gamma <- if (inherits(backend, "vglm")) {
    t(stats::coef(backend, matrix = TRUE))
  } else {
    get_effective_coefs(backend, beta = coef)
  }
  mapped_gamma <- matrix(
    drop(map %*% coef),
    nrow = M,
    ncol = P,
    dimnames = gamma_dimnames
  )
  if (
    !identical(dim(gamma), dim(mapped_gamma)) ||
      !isTRUE(all.equal(
        unname(mapped_gamma),
        unname(gamma),
        tolerance = 1e-12,
        check.attributes = FALSE
      ))
  ) {
    mostr_abort(
      "The raw-to-effective coefficient map does not reproduce the fitted model."
    )
  }

  rownames(map) <- as.vector(outer(
    gamma_dimnames[[1L]],
    gamma_dimnames[[2L]],
    function(eta, term) paste0(term, "::", eta)
  ))
  attr(map, "gamma_dim") <- c(M, P)
  attr(map, "gamma_dimnames") <- gamma_dimnames
  map
}

sop_delta_condition <- function(message, class, ...) {
  mostr_abort(message, class = class, data = list(...))
}

sop_delta_validate_plan <- function(plan, model, map, coef) {
  if (!inherits(plan, "markov_sop_exec_plan") || !identical(plan$version, 1L)) {
    mostr_abort("Invalid or unsupported SOP execution plan.")
  }
  if (!identical(plan$recursion_order, 1L)) {
    sop_delta_condition(
      "Analytical SOP inference currently supports first-order Markov plans only.",
      "mostr_delta_unsupported_order",
      recursion_order = plan$recursion_order
    )
  }
  if (!identical(plan$link, "logit")) {
    sop_delta_condition(
      "Analytical SOP inference currently supports the reverse cumulative logit link only.",
      "mostr_delta_unsupported_link",
      link = plan$link
    )
  }

  backend <- sop_delta_backend_model(model)
  if (!inherits(backend, c("vglm", "orm")) || inherits(backend, "blrm")) {
    mostr_abort(
      "Analytical SOP inference supports only frequentist `vglm` and `orm` models."
    )
  }
  validate_markov_model(backend)
  if (!identical(class(backend)[1L], plan$model_class)) {
    mostr_abort(
      "The execution plan and fitted model use different backends."
    )
  }

  components <- plan$components
  required_component_names <- c(
    "X_init",
    "X_transition",
    "n_pat",
    "n_states",
    "M",
    "y_levels",
    "col_names",
    "transition_layout",
    "transition_origins"
  )
  if (
    !is.list(components) ||
      any(!required_component_names %in% names(components))
  ) {
    mostr_abort(
      "The SOP execution plan is missing first-order design components."
    )
  }
  if (!identical(components$transition_layout, "stacked")) {
    mostr_abort(
      "Analytical SOP inference requires stacked transition designs."
    )
  }
  if (
    !is.matrix(components$X_init) ||
      nrow(components$X_init) != components$n_pat ||
      ncol(components$X_init) != length(components$col_names)
  ) {
    mostr_abort(
      "The SOP initial design matrix is inconsistent with the plan."
    )
  }
  if (length(plan$times) < 1L || components$n_states != length(plan$y_levels)) {
    mostr_abort(
      "The SOP execution plan has inconsistent time or state metadata."
    )
  }
  if (components$M != components$n_states - 1L) {
    mostr_abort(
      "The SOP execution plan has an inconsistent threshold count."
    )
  }

  gamma_dim <- attr(map, "gamma_dim", exact = TRUE)
  gamma_dimnames <- attr(map, "gamma_dimnames", exact = TRUE)
  if (
    length(gamma_dim) != 2L ||
      gamma_dim[1L] != components$M ||
      ncol(map) != length(coef)
  ) {
    mostr_abort(
      "The fitted coefficient map is incompatible with the SOP plan."
    )
  }
  gamma_columns <- gamma_dimnames[[2L]]
  column_index <- match(components$col_names, gamma_columns)
  if (anyNA(column_index)) {
    mostr_abort(
      "The fitted coefficient map is missing execution-plan design columns: ",
      paste(components$col_names[is.na(column_index)], collapse = ", ")
    )
  }
  gamma <- matrix(
    drop(map %*% coef),
    nrow = gamma_dim[1L],
    ncol = gamma_dim[2L],
    dimnames = gamma_dimnames
  )
  gamma <- gamma[, column_index, drop = FALSE]
  map_rows <- unlist(lapply(column_index, function(index) {
    (index - 1L) * components$M + seq_len(components$M)
  }))
  plan_map <- map[map_rows, , drop = FALSE]

  intercept <- grep("Intercept", components$col_names, fixed = TRUE)
  if (
    length(intercept) != 1L ||
      any(components$X_init[, intercept] != 1)
  ) {
    sop_delta_condition(
      "Analytical SOP inference requires one constant intercept design column.",
      "mostr_delta_not_full_po"
    )
  }
  slope_columns <- setdiff(seq_along(components$col_names), intercept)
  for (column in slope_columns) {
    rows <- (column - 1L) * components$M + seq_len(components$M)
    block <- plan_map[rows, , drop = FALSE]
    reference <- matrix(
      block[1L, ],
      nrow = components$M,
      ncol = ncol(block),
      byrow = TRUE
    )
    if (any(abs(block - reference) > 1e-12)) {
      sop_delta_condition(
        paste0(
          "Analytical SOP inference requires full proportional odds; design ",
          "column `",
          components$col_names[column],
          "` has threshold-specific slopes."
        ),
        "mostr_delta_not_full_po",
        design_column = components$col_names[column]
      )
    }
  }

  list(
    components = components,
    Gamma = gamma,
    map = plan_map,
    intercept = intercept
  )
}

sop_delta_preflight <- function(
  plan,
  coefficient_count,
  average_group_size = NULL,
  retain_individual_probabilities = FALSE
) {
  components <- plan$components
  n <- as.double(components$n_pat)
  time_count <- as.double(length(plan$times))
  state_count <- as.double(components$n_states)
  threshold_count <- as.double(components$M)
  q <- as.double(coefficient_count)

  grouped <- !is.null(average_group_size)
  if (grouped) {
    if (
      length(average_group_size) != 1L ||
        !is.numeric(average_group_size) ||
        is.na(average_group_size) ||
        average_group_size < 1L ||
        average_group_size != as.integer(average_group_size) ||
        components$n_pat %% as.integer(average_group_size) != 0L
    ) {
      mostr_abort("Invalid analytical SOP averaging group size.")
    }
    output_profiles <- n / as.double(average_group_size)
  } else {
    if (isTRUE(retain_individual_probabilities)) {
      mostr_abort(
        "Individual analytical probabilities are already retained for fixed profiles."
      )
    }
    output_profiles <- n
  }

  output_cells <- output_profiles * time_count * state_count * (q + 1)
  individual_probability_cells <- if (
    grouped && isTRUE(retain_individual_probabilities)
  ) {
    n * time_count * state_count
  } else {
    0
  }
  recursion_cells <- 2 * n * state_count * (q + 1)
  category_workspace_cells <-
    threshold_count * (q + 1) + state_count * (q + 1) + q
  required_bytes <- (output_cells +
    individual_probability_cells +
    recursion_cells +
    category_workspace_cells) *
    8
  delta_assert_bytes(required_bytes, "The analytical SOP Jacobian calculation")
  required_bytes
}

delta_native_failure <- function(result, plan) {
  status <- as.integer(result$status)
  if (length(status) != 3L || anyNA(status)) {
    mostr_abort(
      "The native analytical SOP engine returned invalid status metadata."
    )
  }
  code <- status[1L]
  if (code == 0L) {
    return(invisible(NULL))
  }
  visit <- status[2L]
  origin <- status[3L]
  value <- as.numeric(result$value)[1L]
  context <- if (visit == 1L) {
    "the first SOP time point"
  } else if (
    visit >= 2L &&
      origin >= 1L &&
      origin <= length(plan$y_levels)
  ) {
    paste0(
      "SOP time point ",
      visit,
      ", origin state ",
      plan$y_levels[origin]
    )
  } else {
    paste0("SOP time point ", visit)
  }

  if (code == 1L) {
    sop_delta_condition(
      paste0(
        "Non-finite raw category probabilities occurred at ",
        context,
        "."
      ),
      "mostr_delta_invalid_probabilities",
      context = context
    )
  }
  if (code == 2L) {
    sop_delta_condition(
      paste0(
        "Raw ordinal category probabilities are crossed at ",
        context,
        " (minimum ",
        format(value, digits = 6),
        ")."
      ),
      "mostr_delta_crossed_probabilities",
      context = context,
      minimum_probability = value
    )
  }
  if (code == 3L) {
    sop_delta_condition(
      paste0("Raw category probabilities do not sum to one at ", context, "."),
      "mostr_delta_invalid_probabilities",
      context = context
    )
  }
  if (code == 4L) {
    mostr_abort(
      "Internal analytical category derivatives do not sum to zero."
    )
  }
  if (code == 5L) {
    mostr_abort(
      "Internal analytical SOP probabilities do not sum to one."
    )
  }
  if (code == 6L) {
    mostr_abort(
      "Internal analytical SOP derivatives do not sum to zero."
    )
  }
  mostr_abort(
    "The native analytical SOP engine returned an unknown failure status."
  )
}

run_sop_delta_plan <- function(
  plan,
  model,
  average_group_size = NULL,
  retain_individual_probabilities = FALSE
) {
  coef <- sop_delta_raw_coef(model)
  map <- get_effective_coef_map(model)
  validated <- sop_delta_validate_plan(plan, model, map, coef)
  components <- validated$components
  Gamma <- validated$Gamma
  map <- validated$map
  sop_delta_preflight(
    plan,
    length(coef),
    average_group_size = average_group_size,
    retain_individual_probabilities = retain_individual_probabilities
  )

  n <- components$n_pat
  K <- components$n_states
  absorb <- which(
    as.character(plan$y_levels) %in% as.character(plan$absorb)
  )
  non_absorb <- setdiff(seq_len(K), absorb)
  origin_positions <- match(non_absorb, components$transition_origins)
  if (anyNA(origin_positions)) {
    mostr_abort("Transition design is missing a required origin state.")
  }

  native <- cpp_run_sop_delta(
    initial_design = components$X_init,
    transition_designs = components$X_transition,
    gamma = Gamma,
    coefficient_map = map,
    origin_positions = origin_positions,
    non_absorb = non_absorb,
    absorb = absorb,
    intercept = validated$intercept,
    group_size = average_group_size %||% 0L,
    retain_individual_probabilities = isTRUE(retain_individual_probabilities)
  )
  delta_native_failure(native, plan)

  probabilities <- native$probabilities
  jacobian <- native$jacobian
  dimnames(probabilities) <- list(NULL, NULL, plan$y_levels)
  dimnames(jacobian) <- list(NULL, NULL, plan$y_levels, names(coef))
  out <- list(
    probabilities = probabilities,
    jacobian = jacobian,
    coef = coef,
    coef_names = names(coef)
  )
  if (isTRUE(retain_individual_probabilities)) {
    individual <- native$individual_probabilities
    if (!identical(dim(individual), c(n, length(plan$times), K))) {
      mostr_abort(
        "The native analytical SOP individual probabilities have unexpected dimensions."
      )
    }
    dimnames(individual) <- list(NULL, NULL, plan$y_levels)
    out$individual_probabilities <- individual
  }
  out
}
