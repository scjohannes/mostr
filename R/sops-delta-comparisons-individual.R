# Analytical propagation for results with one row per prediction profile:
# comparisons from `comparisons()` and times from `state_time()`.

inferences_delta_profile_results <- function(
  object,
  target,
  vcov,
  cluster,
  conf_level,
  conf_type
) {
  args <- attr(object, "comparison_args")
  if (is.null(args)) {
    stop("Stored result arguments are required for delta inference.")
  }
  if (!identical(conf_type, "wald")) {
    stop(
      "Analytical intervals for one row per profile require ",
      "`conf_type = \"wald\"`.",
      call. = FALSE
    )
  }
  is_time <- identical(profile_result_kind(object), "time")
  if (!is_time && !identical(args$comparison, "difference")) {
    stop(
      "Analytical comparisons currently support differences only; ratios ",
      "are not yet supported.",
      call. = FALSE
    )
  }
  if (!is_time && !args$estimand %in% c("sop", "time_in_state")) {
    stop(
      "Analytical comparisons currently support `estimand = \"sop\"` or ",
      "`\"time_in_state\"`; time benefit is not yet supported.",
      call. = FALSE
    )
  }
  model <- attr(object, "model")
  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform delta inference.")
  }
  model_plan <- delta_validate_sop_scope(object, model, by = NULL)
  setup <- setup_from_avg_comparison_object(object)
  index <- if (is_time) {
    delta_time_row_index(object, setup, args)
  } else {
    delta_comparison_row_index(object, setup, args)
  }
  time_weights <- delta_comparison_time_weights(setup, args)

  block_values <- function(block) {
    delta_comparison_block(
      object = object,
      model_plan = model_plan,
      setup = setup,
      block = block,
      state_sets = args$state_sets,
      time_weights = time_weights
    )
  }
  jacobian <- NULL
  estimate <- rep(NA_real_, nrow(object))
  start_jacobian <- function(coefficient_names) {
    delta_assert_bytes(
      as.double(nrow(object)) * length(coefficient_names) * 8,
      "Analytical derivatives for one row per profile"
    )
    matrix(
      NA_real_,
      nrow = nrow(object),
      ncol = length(coefficient_names),
      dimnames = list(NULL, coefficient_names)
    )
  }

  # Scenarios are differentiated one block of profiles at a time, so the
  # derivative arrays of at most two scenarios are held at once.
  if (is_time) {
    # A profile's time under one scenario uses that scenario's derivatives.
    for (scenario in unique(index$scenario)) {
      block <- block_values(scenario)
      jacobian <- jacobian %||% start_jacobian(block$coefficient_names)
      for (set in unique(index$set[index$scenario == scenario])) {
        rows <- which(index$scenario == scenario & index$set == set)
        cells <- index$profile[rows]
        jacobian[rows, ] <- block$jacobian[[set]][cells, , drop = FALSE]
        estimate[rows] <- block$value[[set]][cells]
      }
    }
  } else {
    # A comparison's derivative is the comparison scenario's minus the
    # reference scenario's, for the same profile.
    reference <- block_values(1L)
    jacobian <- start_jacobian(reference$coefficient_names)
    for (level in unique(index$level)) {
      comparison <- block_values(level)
      for (set in unique(index$set[index$level == level])) {
        rows <- which(index$level == level & index$set == set)
        cells <- index$profile[rows] + setup$n_each * (index$time[rows] - 1L)
        jacobian[rows, ] <- comparison$jacobian[[set]][cells, , drop = FALSE] -
          reference$jacobian[[set]][cells, , drop = FALSE]
        estimate[rows] <- comparison$value[[set]][cells] -
          reference$value[[set]][cells]
      }
    }
  }

  stored <- object$estimate
  if (
    any(!is.finite(estimate)) ||
      any(
        abs(estimate - stored) >
          1e-12 + 1e-10 * pmax(abs(estimate), abs(stored))
      )
  ) {
    stop(
      "The analytical calculation did not reproduce the stored point ",
      "estimates.",
      call. = FALSE
    )
  }

  coefficient <- delta_coefficient_state(model, cluster, vcov, jacobian)
  delta_finalize_result(
    object = object,
    standard_error = coefficient$std.error,
    conf_level = conf_level,
    conf_type = "wald",
    target = target,
    analytical = list(
      representation = "coefficient",
      jacobian = jacobian,
      coefficient_vcov = coefficient$coefficient_vcov,
      coefficient_names = colnames(jacobian),
      covariance_metadata = coefficient$metadata
    )
  )
}

# Positions of each result row's profile, time, state set, and scenario.
delta_comparison_row_index <- function(object, setup, args) {
  values <- setup$variables[[1]]
  index <- list(
    profile = as.integer(object$rowid),
    time = if (identical(args$estimand, "sop")) {
      match(as.character(object$time), as.character(setup$times))
    } else {
      rep(1L, nrow(object))
    },
    set = match(as.character(object$state_set), names(args$state_sets)),
    level = match(as.character(object$comparison_level), as.character(values))
  )
  if (
    anyNA(unlist(index, use.names = FALSE)) ||
      any(index$profile < 1L | index$profile > setup$n_each) ||
      any(index$level < 2L)
  ) {
    stop("Comparison rows do not align with the stored prediction profiles.")
  }
  index
}

# Positions of each `state_time()` row's profile, scenario, and state set.
# `rowid` numbers the scenario-profile pairs as the stacked prediction data do.
delta_time_row_index <- function(object, setup, args) {
  index <- list(
    profile = as.integer(profile_index(object)),
    scenario = as.integer((object$rowid - 1L) %/% setup$n_each + 1L),
    set = match(as.character(object$state_set), names(args$state_sets))
  )
  if (
    anyNA(unlist(index, use.names = FALSE)) ||
      any(index$profile < 1L | index$profile > setup$n_each) ||
      any(index$scenario > setup$n_cf)
  ) {
    stop("Time rows do not align with the stored prediction profiles.")
  }
  index
}

# Weights that turn visit SOPs into the reported metric. Real-time totals
# are linear in the visit SOPs plus the starting state at `baseline_time`,
# which has no coefficient derivative.
delta_comparison_time_weights <- function(setup, args) {
  if (identical(args$estimand, "sop")) {
    return(NULL)
  }
  if (is.null(args$time_map)) {
    return(list(visit = rep(1, length(setup$times)), baseline = 0))
  }
  weights <- delta_real_time_weights(data.frame(time = setup$times), args)
  list(
    visit = unname(weights$visit[as.character(setup$times)]),
    baseline = weights$baseline
  )
}

# Values and coefficient derivatives of one scenario's state-set metric for
# every profile. SOP cells are ordered profile fastest, then time.
delta_comparison_block <- function(
  object,
  model_plan,
  setup,
  block,
  state_sets,
  time_weights
) {
  n_each <- setup$n_each
  n_times <- length(setup$times)
  rows <- ((block - 1L) * n_each + 1L):(block * n_each)
  newdata <- setup$newdata_pred[rows, , drop = FALSE]
  plan <- delta_compile_plan(
    object = object,
    model_plan = model_plan,
    newdata = newdata,
    times = setup$times,
    output = "individual"
  )
  result <- delta_validate_plan_result(
    run_sop_delta_plan(plan, model_plan),
    n_profiles = n_each,
    n_times = n_times,
    n_states = length(setup$y_levels)
  )
  n_coef <- length(result$coefficient_names)
  state_labels <- as_state_labels(setup$y_levels)
  start <- as.character(newdata[[setup$p_var]])

  value <- jacobian <- vector("list", length(state_sets))
  for (set in seq_along(state_sets)) {
    states <- as_state_labels(state_sets[[set]])
    probability <- matrix(0, nrow = n_each, ncol = n_times)
    derivative <- array(0, dim = c(n_each, n_times, n_coef))
    for (state in match(states, state_labels)) {
      probability <- probability +
        matrix(result$probabilities[, , state], nrow = n_each, ncol = n_times)
      derivative <- derivative +
        array(
          result$jacobian[, , state, , drop = FALSE],
          dim = c(n_each, n_times, n_coef)
        )
    }
    if (is.null(time_weights)) {
      value[[set]] <- as.vector(probability)
      jacobian[[set]] <- matrix(
        derivative,
        nrow = n_each * n_times,
        ncol = n_coef
      )
    } else {
      value[[set]] <- as.vector(probability %*% time_weights$visit) +
        time_weights$baseline * as.numeric(start %in% states)
      jacobian[[set]] <- matrix(
        matrix(
          aperm(derivative, c(1L, 3L, 2L)),
          nrow = n_each * n_coef,
          ncol = n_times
        ) %*%
          time_weights$visit,
        nrow = n_each,
        ncol = n_coef
      )
    }
  }
  list(
    value = value,
    jacobian = jacobian,
    coefficient_names = result$coefficient_names
  )
}
