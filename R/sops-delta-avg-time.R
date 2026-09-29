# Analytical propagation for average time in states.

inferences_delta_avg_time <- function(
  object,
  target,
  vcov,
  cluster,
  conf_level,
  conf_type
) {
  args <- attr(object, "time_args")
  avg_args <- attr(object, "avg_args")
  if (is.null(args) || is.null(avg_args)) {
    stop("Stored average-time arguments are required for delta inference.")
  }
  if (!identical(conf_type, "wald")) {
    stop(
      'Analytical average time requires `conf_type = "wald"`.',
      call. = FALSE
    )
  }
  target <- delta_validate_target(object, target)
  avg <- replay_avg_time_sops(object, conf_level, return_draws = FALSE)
  avg <- inferences_delta_sops(
    object = avg,
    target = target,
    vcov = vcov,
    cluster = cluster,
    conf_level = conf_level,
    conf_type = "wald"
  )

  state_sets <- args$state_sets
  state_index <- match(as.character(object$state_set), names(state_sets))
  if (anyNA(state_index) || anyDuplicated(names(state_sets))) {
    stop("Average-time state sets do not align with the replayed average SOPs.")
  }
  source <- as.data.frame(avg)
  variables <- names(avg_args$variables)
  if (
    !all(c("time", "state", "estimate", variables) %in% names(source)) ||
      !all(variables %in% names(object))
  ) {
    stop("The replayed average SOP cells are missing average-time identifiers.")
  }
  source_time <- as.character(source$time)
  integration <- if (is.null(args$time_map)) {
    list(
      visit = stats::setNames(
        rep(1, length(unique(source_time))),
        unique(source_time)
      ),
      baseline = 0
    )
  } else {
    delta_real_time_weights(avg, args)
  }
  time_weight <- unname(integration$visit[match(
    source_time,
    names(integration$visit)
  )])
  if (anyNA(time_weight)) {
    stop("Average-time visits do not align with the replayed average SOPs.")
  }
  p_var <- attr(avg, "p_var") %||% "yprev"
  # Scenarios that set the starting state anchor at that state (a constant).
  set_start <- p_var %in% variables
  baseline <- if (integration$baseline != 0 && !set_start) {
    baseline_rows_for_anchor(avg, avg_args$id_var)
  } else {
    NULL
  }
  starting_state <- if (!is.null(baseline)) {
    if (!p_var %in% names(baseline)) {
      stop("Previous-state variable `", p_var, "` not found in baseline data.")
    }
    as.character(baseline[[p_var]])
  } else {
    NULL
  }

  n_result <- nrow(object)
  stored_bytes <- as.double(n_result) * 24
  delta_assert_bytes(
    stored_bytes + as.double(nrow(source)) * 64,
    "The analytical average-time operator"
  )
  operator <- vector("list", n_result)
  baseline_mean <- numeric(n_result)
  # ponytail: scan source cells per result; index groups if setup time dominates.
  for (i in seq_len(n_result)) {
    states <- as.character(state_sets[[state_index[i]]])
    selected <- as.character(source$state) %in% states
    if (length(variables)) {
      selected <- selected &
        rows_match_values(source, object[i, variables, drop = FALSE])
    }
    weight <- time_weight * selected
    index <- which(weight != 0)
    stored_bytes <- stored_bytes + length(index) * 12
    delta_assert_bytes(
      stored_bytes + as.double(nrow(source)) * 64,
      "The analytical average-time operator"
    )
    operator[[i]] <- list(index = index, weight = weight[index])
    if (!is.null(baseline)) {
      baseline_mean[i] <- mean(starting_state %in% states)
    } else if (set_start && integration$baseline != 0) {
      baseline_mean[i] <- as.numeric(
        as.character(object[[p_var]][i]) %in% states
      )
    }
  }
  calculated <- drop(delta_apply_comparison_operator(
    operator,
    source$estimate
  )) +
    integration$baseline * baseline_mean
  if (
    any(!is.finite(calculated)) ||
      any(
        abs(calculated - object$estimate) >
          1e-12 + 1e-10 * pmax(abs(calculated), abs(object$estimate))
      )
  ) {
    stop(
      "The analytical average-time operator did not reproduce the stored point estimates.",
      call. = FALSE
    )
  }

  propagated <- delta_propagate_comparison_state(avg, operator)
  if (
    identical(propagated$analytical$representation, "influence") &&
      !is.null(baseline)
  ) {
    ids <- delta_profile_ids(baseline, avg_args$id_var %||% attr(avg, "id_var"))
    profile_ids <- propagated$analytical$profile_ids
    index <- match(as.character(profile_ids), as.character(ids))
    if (anyNA(index) || length(ids) != length(profile_ids)) {
      stop(
        "Baseline-anchor patients do not align with the analytical patient influences."
      )
    }
    starting_state <- starting_state[index]
    n <- length(profile_ids)
    delta_assert_bytes(
      (as.double(n) * n_result + n * 3) * 8,
      "Analytical baseline-anchor propagation"
    )
    for (i in seq_len(n_result)) {
      indicator <- as.numeric(starting_state %in% state_sets[[state_index[i]]])
      influence <- propagated$analytical$influence[, i] +
        integration$baseline * (indicator - baseline_mean[i])
      propagated$analytical$influence[, i] <- influence
      propagated$std.error[i] <- sqrt(stats::var(influence) / n)
    }
  }
  delta_finalize_result(
    object = object,
    standard_error = propagated$std.error,
    conf_level = conf_level,
    conf_type = "wald",
    target = target,
    analytical = propagated$analytical
  )
}
