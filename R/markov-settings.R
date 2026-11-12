# One object for the column names and model settings that the SOP,
# comparison, inference and diagnostic code would otherwise pass around one by
# one.
#
# `mostm()` builds it once and stores it on the fit (`markov_set_settings()`);
# `markov_model_settings()` reads it back. A public SOP call merges its own
# arguments into the fit's settings with `resolve_markov_settings()`, which
# refuses a value that contradicts the fit, and passes the result to the
# internal functions. Results keep that call's settings in their `settings`
# attribute (`markov_result_settings()`), so inference replays the same
# calculation.

#' Settings of a fitted model or of one SOP call
#'
#' @param engine The engine that fitted the model (`"orm"`, `"vglm"`, `"blrm"`)
#'   or `NULL` when unknown.
#' @param id_var,time_var,p_var,p2_var,gap_var Column names: patient ID, visit
#'   time, previous state, state two visits back (`NULL` for first order), and
#'   elapsed time since the previous visit.
#' @param y_levels State labels in fitted order.
#' @param absorb Absorbing states.
#' @param time_covariates Lookup table of time-varying covariates by visit.
#' @param first_followup_time The visit of the stored starting profiles.
#' @param from_fit Names of the fields that the fit fixes. A public call may
#'   repeat these but not change them.
#' @return An object of class `markov_settings`; `order` is the Markov order.
#' @noRd
markov_settings <- function(
  engine = NULL,
  id_var = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  y_levels = NULL,
  absorb = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  first_followup_time = NULL,
  from_fit = character()
) {
  for (arg in c("id_var", "time_var", "p_var", "p2_var", "gap_var")) {
    value <- get(arg, inherits = FALSE)
    if (
      !is.null(value) &&
        (!is.character(value) ||
          length(value) != 1L ||
          is.na(value) ||
          !nzchar(value))
    ) {
      mostr_abort("`", arg, "` must be a single column name.")
    }
  }
  if (is.null(time_var) || is.null(p_var)) {
    mostr_abort(
      "`",
      if (is.null(time_var)) "time_var" else "p_var",
      "` must be a single column name."
    )
  }
  if (identical(p2_var, p_var)) {
    mostr_abort(
      "`p2_var` must name the column with the state two visits back, which ",
      "differs from the previous-state column `p_var = \"",
      p_var,
      "\"`."
    )
  }
  if (!is.null(engine) && !engine %in% markov_engines()) {
    mostr_abort("Internal Markov model engine provenance is invalid.")
  }
  from_fit <- intersect(from_fit, markov_settings_fit_fields())
  structure(
    list(
      engine = engine,
      id_var = id_var,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      order = if (is.null(p2_var)) 1L else 2L,
      y_levels = y_levels,
      absorb = absorb,
      gap_var = gap_var,
      time_covariates = time_covariates,
      first_followup_time = first_followup_time,
      from_fit = from_fit
    ),
    class = "markov_settings"
  )
}

# The fields that `mostm()` fixes for every later call on the fit.
markov_settings_fit_fields <- function() {
  c("engine", "id_var", "time_var", "p_var", "p2_var", "y_levels")
}

# A copy of `settings` with some fields replaced, validated again. Bootstrap
# refits use it for the states present in a resample.
markov_settings_update <- function(settings, ...) {
  changes <- list(...)
  fields <- unclass(settings)
  fields$order <- NULL
  fields[names(changes)] <- changes
  do.call(markov_settings, fields)
}

markov_set_settings <- function(model, settings) {
  if (!inherits(settings, "markov_settings")) {
    mostr_abort("Internal Markov model settings are invalid.")
  }
  attr(model, "markov_settings") <- settings
  model
}

# The settings recorded on a fit. Fits made by `mostm()` store them. Other
# objects, such as the hand-made fits of the package tests, get settings
# built from what they record (engine and patient ID), with nothing fixed, so
# calls on them behave as they did before settings were stored.
markov_model_settings <- function(model) {
  settings <- markov_model_metadata_attr(model, "markov_settings")
  if (inherits(settings, "markov_settings")) {
    return(settings)
  }
  markov_settings(
    engine = markov_model_engine(model),
    id_var = markov_model_id_var(model),
    first_followup_time =
      markov_model_starting_profile_metadata(model)$first_followup_time
  )
}

# The settings of a SOP, comparison or time result.
markov_result_settings <- function(object) {
  settings <- attr(object, "settings", exact = TRUE)
  if (!inherits(settings, "markov_settings")) {
    mostr_abort(
      "The result does not record the settings of the call that made it. ",
      "Recompute it with the current version of mostr."
    )
  }
  settings
}

# The settings of one public call: the call's arguments, with every argument
# left `NULL` taken from the fit. An argument that the fit fixes must agree
# with it. `y_levels` is compared by state label and the supplied vector is
# kept, because the result's `state` column takes its type.
resolve_markov_settings <- function(
  model,
  id_var = NULL,
  time_var = NULL,
  p_var = NULL,
  p2_var = NULL,
  y_levels = NULL,
  absorb = NULL,
  gap_var = NULL,
  time_covariates = NULL
) {
  fit <- markov_model_settings(model)
  supplied <- list(
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels
  )
  for (arg in names(supplied)) {
    value <- supplied[[arg]]
    if (is.null(value) || !arg %in% fit$from_fit) {
      next
    }
    if (identical(arg, "y_levels")) {
      # A wrong number of states is reported by the threshold-count check,
      # `validate_sop_threshold_count()`, as before.
      if (
        length(value) == length(fit$y_levels) &&
          !identical(as_state_labels(value), as_state_labels(fit$y_levels))
      ) {
        stop_y_levels_contradict_fit(value, fit$y_levels)
      }
    } else if (!identical(value, fit[[arg]])) {
      stop_setting_contradicts_fit(arg, value, fit[[arg]])
    }
  }
  # `p2_var = NULL` means "as fitted"; on a fit without a fixed order it is
  # first order, as before.
  p2 <- if ("p2_var" %in% fit$from_fit) fit$p2_var else p2_var
  markov_settings(
    engine = fit$engine,
    id_var = id_var %||% fit$id_var,
    time_var = time_var %||% fit$time_var,
    p_var = p_var %||% fit$p_var,
    p2_var = p2,
    y_levels = y_levels %||%
      fit$y_levels %||%
      markov_model_ylevels_or_null(model),
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
    first_followup_time = fit$first_followup_time,
    from_fit = fit$from_fit
  )
}

markov_model_ylevels_or_null <- function(model) {
  tryCatch(markov_model_ylevels(model), error = function(e) NULL)
}

format_setting_value <- function(value) {
  if (is.null(value)) {
    return("NULL")
  }
  paste0("\"", value, "\"")
}

stop_setting_contradicts_fit <- function(arg, value, fitted) {
  advice <- if (identical(arg, "p2_var") && is.null(fitted)) {
    paste0(
      "The model was fitted as a first-order model. For a second-order model, ",
      "refit it with `mostm(..., p2_var = ",
      format_setting_value(value),
      ")`."
    )
  } else {
    paste0(
      "Leave `",
      arg,
      "` out to use the fit's value, or refit the model with `mostm(..., ",
      arg,
      " = ",
      format_setting_value(value),
      ")`."
    )
  }
  mostr_abort(
    "`",
    arg,
    " = ",
    format_setting_value(value),
    "` does not match the fitted model: `mostm()` fitted it with `",
    arg,
    " = ",
    format_setting_value(fitted),
    "`. ",
    advice
  )
}

stop_y_levels_contradict_fit <- function(value, fitted) {
  mostr_abort(
    "`y_levels` must be the fitted model's states in fitted order (",
    paste(as_state_labels(fitted), collapse = ", "),
    "), but it is ",
    paste(as_state_labels(value), collapse = ", "),
    ". Leave `y_levels` out to use the fitted states."
  )
}

# The previous-state columns that `mostm()` records: names of columns of
# `data`, distinct from each other.
validate_mostm_previous_state_vars <- function(p_var, p2_var, data) {
  for (arg in c("p_var", "p2_var")) {
    value <- get(arg, inherits = FALSE)
    if (is.null(value) && identical(arg, "p2_var")) {
      next
    }
    if (
      !is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(value)
    ) {
      mostr_abort("`", arg, "` must be a single column name.")
    }
  }
  if (identical(p_var, p2_var)) {
    mostr_abort(
      "`p2_var` must name the column with the state two visits back, which ",
      "differs from the previous-state column `p_var = \"",
      p_var,
      "\"`."
    )
  }
  if (!is.null(p2_var) && !p2_var %in% names(data)) {
    mostr_abort(
      "The column `",
      p2_var,
      "` given as `p2_var` is not in `data`. Set `p2_var` to the column ",
      "that holds each patient's state two visits back."
    )
  }
  invisible(NULL)
}
