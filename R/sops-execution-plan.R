# Internal serializable execution plans for frequentist SOPs.

validate_sop_threshold_count <- function(Gamma, y_levels) {
  expected <- length(y_levels) - 1L
  actual <- nrow(Gamma)
  if (!identical(actual, expected)) {
    mostr_abort(
      "`y_levels` defines ",
      length(y_levels),
      " states, but the fitted model defines ",
      actual + 1L,
      " states through ",
      actual,
      " threshold coefficients."
    )
  }
  invisible(NULL)
}

# Caps the retained designs of all visits: (profiles + (visits - 1) x profiles x
# non-absorbing states) x columns x 8 bytes. Above it, prediction uses the R
# engine and delta inference stops. 256 MiB was a choice, not a measurement. Nov
# 2026, 25k profiles x 29 visits: 10 MVN draws took 69 s in the R engine against
# 19 s compiled (0.29 vs 0.72 GB peak), so it is on the conservative side.
execution_plan_max_bytes <- function() {
  value <- getOption("mostr.execution_plan_max_bytes", 256 * 1024^2)
  if (
    length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      value <= 0
  ) {
    mostr_abort(
      "Option `mostr.execution_plan_max_bytes` must be a positive number."
    )
  }
  as.double(value)
}

stop_execution_plan_too_large <- function(required_bytes, max_bytes) {
  message <- paste0(
    "Compiled SOP designs require ",
    format(required_bytes, scientific = FALSE),
    " bytes, above the execution-plan limit of ",
    format(max_bytes, scientific = FALSE),
    " bytes."
  )
  mostr_abort(
    message,
    class = "mostr_execution_plan_too_large",
    data = list(
      required_bytes = required_bytes,
      max_bytes = max_bytes
    )
  )
}

preflight_execution_plan_designs <- function(
  n_pat,
  n_times,
  rows_per_transition,
  n_cols
) {
  initial_cells <- as.double(n_pat) * as.double(n_cols)
  transition_cells <-
    as.double(max(n_times - 1L, 0L)) *
    as.double(rows_per_transition) *
    as.double(n_cols)
  design_bytes <- (initial_cells + transition_cells) * 8
  workspace_bytes <- execution_plan_max_bytes()
  if (!is.finite(design_bytes) || design_bytes > workspace_bytes) {
    stop_execution_plan_too_large(design_bytes, workspace_bytes)
  }
  list(
    design_bytes = design_bytes,
    workspace_bytes = workspace_bytes
  )
}

# A plan for the SOPs of `newdata` at `times`: the prediction designs of every
# visit, built once (`markov_msm_build()`), with the columns and states named
# in `settings` (`resolve_markov_settings()`). `run_sop_execution_plan()` runs
# it for any coefficients.
compile_sop_execution_plan <- function(model, newdata, times, settings) {
  validate_markov_model(model)
  newdata <- prepare_prediction_rows(model, newdata, settings)
  components <- markov_msm_build(
    model = model,
    newdata = newdata,
    times = times,
    settings = settings
  )
  second_order <- !is.null(settings$p2_var)
  if (second_order) {
    # What the second-order recursion reads besides the designs: each
    # patient's starting state and, for each block of a later visit's design,
    # the state two visits back (`older`) and one visit back (`current`).
    non_absorb <- components$transition_origins
    components <- c(
      components,
      list(
        previous = match_state_indices(
          newdata[[settings$p_var]],
          components$y_levels,
          settings$p_var
        ),
        older = components$previous_states[[settings$p2_var]],
        current = components$previous_states[[settings$p_var]],
        absorb = setdiff(seq_len(components$n_states), non_absorb),
        non_absorb = non_absorb
      )
    )
  }
  structure(
    list(
      version = 1L,
      model_class = class(model)[1L],
      link = "logit",
      recursion_order = if (second_order) 2L else 1L,
      times = components$times,
      y_levels = components$y_levels,
      absorb = as_state_labels(settings$absorb),
      workspace_bytes = components$workspace_bytes,
      design_bytes = components$design_bytes,
      components = components
    ),
    class = "markov_sop_exec_plan"
  )
}

run_sop_execution_plan <- function(plan, Gamma) {
  if (!inherits(plan, "markov_sop_exec_plan") || plan$version != 1L) {
    mostr_abort("Invalid or unsupported SOP execution plan.")
  }
  validate_sop_threshold_count(Gamma, plan$y_levels)
  if (identical(plan$recursion_order, 2L)) {
    return(run_second_order_execution_plan(plan, Gamma))
  }
  markov_msm_run(
    plan$components,
    Gamma,
    times = plan$times,
    absorb = plan$absorb
  )
}

run_second_order_execution_plan <- function(plan, Gamma) {
  components <- plan$components
  Gamma <- Gamma[, components$col_names, drop = FALSE]
  po <- markov_po_structure(Gamma, components$col_names, components$X_init)
  probabilities <- function(X) {
    if (is.null(po)) {
      return(lp_to_probs(X %*% t(Gamma), nrow(Gamma)))
    }
    lp_to_probs(
      outer(drop(X %*% po$beta), po$cutpoints, "+"),
      length(po$cutpoints)
    )
  }
  initial <- probabilities(components$X_init)
  check_transition_probabilities(initial, "the first SOP time point")
  check_crossed_thresholds(initial, components$n_pat, plan$times[1L])
  output <- array(
    0,
    dim = c(components$n_pat, length(plan$times), components$n_states)
  )
  output[, 1L, ] <- initial
  joint <- array(
    0,
    dim = c(components$n_pat, components$n_states, components$n_states)
  )
  # SOP calls refuse starting rows in an absorbing state before predicting
  # (`stop_absorbing_start()`), so every patient's first visit is predicted
  # from a non-absorbing history.
  for (patient in seq_len(components$n_pat)) {
    joint[patient, components$previous[patient], ] <- initial[patient, ]
  }
  if (length(plan$times) >= 2L) {
    for (visit in 2:length(plan$times)) {
      if (is.null(po)) {
        transition <- probabilities(components$X_transition[[visit]])
        check_crossed_thresholds(
          transition,
          components$n_pat,
          plan$times[visit]
        )
        joint <- markov_update_second_order_native(
          joint,
          transition,
          components$older,
          components$current,
          components$absorb
        )
      } else {
        scalar <- drop(components$X_transition[[visit]] %*% po$beta)
        joint <- markov_update_second_order_po_native(
          joint,
          scalar,
          po$cutpoints,
          components$older,
          components$current,
          components$absorb,
          time = plan$times[visit]
        )
      }
      output[, visit, ] <- apply(joint, c(1L, 3L), sum)
    }
  }
  dimnames(output) <- list(NULL, plan$times, plan$y_levels)
  output
}
