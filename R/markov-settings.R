# One object for the column names and model settings that the SOP,
# comparison, inference and diagnostic code would otherwise pass around one by
# one.
#
# `mostm()` builds it once and stores it on the fit (`markov_set_settings()`);
# `markov_model_settings()` reads it back. A public SOP call reads the fit's
# settings with `resolve_markov_settings()` and passes them to the internal
# functions. The column names, the order, the states, the absorbing states
# and the gap column come from the fit only; the diagnostic plots, which also
# take observed data, may repeat the column names and states but not change
# them. Results keep that call's settings in their `settings` attribute
# (`markov_result_settings()`), so inference replays the same calculation.

#' Settings of a fitted model or of one SOP call
#'
#' @param id_var,time_var,p_var,p2_var,gap_var Column names: patient ID, visit
#'   time, previous state, state two visits back (`NULL` for first order), and
#'   elapsed time since the previous visit.
#' @param y_levels State labels in fitted order.
#' @param absorb Absorbing states, as given to `mostm()`.
#' @param from_fit Names of the fields that the fit fixes. A diagnostic plot
#'   may repeat these but not change them.
#' @return An object of class `markov_settings`; `order` is the Markov order.
#' @noRd
markov_settings <- function(
  id_var = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  y_levels = NULL,
  absorb = NULL,
  gap_var = NULL,
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
  from_fit <- intersect(from_fit, markov_settings_fit_fields())
  structure(
    list(
      id_var = id_var,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      order = if (is.null(p2_var)) 1L else 2L,
      y_levels = y_levels,
      absorb = absorb,
      gap_var = gap_var,
      from_fit = from_fit
    ),
    class = "markov_settings"
  )
}

# The fields that `mostm()` fixes for every later call on the fit.
markov_settings_fit_fields <- function() {
  c(
    "id_var",
    "time_var",
    "p_var",
    "p2_var",
    "y_levels",
    "absorb",
    "gap_var"
  )
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

# The settings recorded on a fit. Fits made by `mostm()` store them, and so do
# the hand-made fits of the package tests (`with_fit_settings()` in
# `tests/testthat/helper-mocks.R`). Other objects get settings built from the
# patient ID they record and the default column names, with nothing fixed.
markov_model_settings <- function(model) {
  settings <- markov_model_metadata_attr(model, "markov_settings")
  if (inherits(settings, "markov_settings")) {
    return(settings)
  }
  markov_settings(id_var = markov_model_id_var(model))
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

# The settings of one public call: the fit's settings. Only the diagnostic
# plots pass `id_var`, `time_var`, `p_var`, `p2_var` and `y_levels`, which they
# also need for observed data; for a fit, every one left `NULL` is taken from
# the fit, and one that the fit fixes must agree with it. `y_levels` is
# compared by state label and the supplied vector is kept, because the plot's
# state labels take its type. The absorbing states and the gap column always
# come from the fit.
resolve_markov_settings <- function(
  model,
  id_var = NULL,
  time_var = NULL,
  p_var = NULL,
  p2_var = NULL,
  y_levels = NULL
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
      # A wrong number of states is refused here too: the plots that pass
      # `y_levels` would otherwise fail later with an unrelated R error.
      if (!identical(as_state_labels(value), as_state_labels(fit$y_levels))) {
        stop_y_levels_contradict_fit(value, fit$y_levels)
      }
    } else if (!identical(value, fit[[arg]])) {
      stop_setting_contradicts_fit(arg, value, fit[[arg]])
    }
  }
  # `p2_var = NULL` means "as fitted"; on an object without stored settings it
  # is first order.
  p2 <- if ("p2_var" %in% fit$from_fit) fit$p2_var else p2_var
  markov_settings(
    id_var = id_var %||% fit$id_var,
    time_var = time_var %||% fit$time_var,
    p_var = p_var %||% fit$p_var,
    p2_var = p2,
    y_levels = y_levels %||%
      fit$y_levels %||%
      markov_model_ylevels_or_null(model),
    absorb = fit$absorb,
    gap_var = fit$gap_var,
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

# Every fit made by `mostm()` records its states; an object without them can't
# label the state occupancy probabilities.
stop_no_fitted_states <- function() {
  mostr_abort(
    "The fitted model does not record its states in fitted order. Refit it ",
    "with `mostm()`."
  )
}

# The arguments that the prediction functions (`sops()`, `avg_sops()`,
# `state_time()`, `avg_state_time()`, `comparisons()`, `avg_comparisons()`,
# `soprob_markov()`) no longer have, because the fit records them, and
# `time_covariates`, which was removed.
markov_removed_fit_arguments <- function() {
  c(
    "id_var",
    "time_var",
    "p_var",
    "p2_var",
    "y_levels",
    "absorb",
    "gap_var",
    "time_covariates"
  )
}

# Refuses a removed argument passed to a prediction function. Their `...`
# would otherwise take it without effect: a call that passes `p2_var` for a
# model fitted without it would silently predict with the first-order
# recursion, and one that passes `absorb` would let patients leave a state
# that the fit does not record as absorbing.
stop_removed_fit_arguments <- function(...) {
  supplied <- intersect(...names(), markov_removed_fit_arguments())
  if (!length(supplied)) {
    return(invisible(NULL))
  }
  arg <- supplied[[1L]]
  if (identical(arg, "time_covariates")) {
    mostr_abort(
      "`time_covariates` was removed. Write time terms inline in the model ",
      "formula given to `mostm()`, for example `rcs(time, 4)`, `log(time)` ",
      "or `pol(time, 2)`; the prediction functions evaluate them at each ",
      "requested time."
    )
  }
  what <- switch(
    arg,
    y_levels = "the fitted model records its states in fitted order",
    absorb = "the fitted model records the absorbing states given to `mostm()`",
    paste0("the fitted model records the `", arg, "` given to `mostm()`")
  )
  advice <- switch(
    arg,
    y_levels = "Leave it out.",
    absorb = paste0(
      "Leave it out; to change them, give them to `mostm(absorb = )` and ",
      "refit."
    ),
    "Leave it out; to use another column, give it to `mostm()` and refit."
  )
  mostr_abort(
    "`",
    arg,
    "` is not an argument of this function: ",
    what,
    ". ",
    advice
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

# The absorbing states given to `mostm()`, checked against the rows it fits
# (`rows`, the rows selected by `subset`). Returns `absorb`, or `NULL` when it
# names no state. Absorbing states are never worked out from the data: a
# state that never occurs as a previous state is often absorbing, but it can
# also be a rare state seen only at a last visit, so such a state only
# gets a warning.
validate_mostm_absorb <- function(
  absorb,
  formula,
  data,
  rows,
  p_var,
  p2_var
) {
  if (!is.null(absorb) && !length(absorb)) {
    absorb <- NULL
  }
  if (!is.null(absorb) && (!is.atomic(absorb) || anyNA(absorb))) {
    mostr_abort(
      "`absorb` must be a vector of state labels without missing values, ",
      "for example `absorb = \"8\"` or `absorb = c(\"home\", \"dead\")`."
    )
  }
  # A response that can't be evaluated here also fails the fit, which then
  # gives the engine's own error.
  response <- tryCatch(
    eval(formula[[2L]], data, environment(formula)),
    error = function(e) NULL
  )
  if (length(formula) != 3L || NROW(response) != nrow(data)) {
    return(absorb)
  }
  response <- response[rows]
  outcome_states <- unique(as_state_labels(response[!is.na(response)]))
  absorb_labels <- unique(as_state_labels(absorb))

  unknown <- setdiff(absorb_labels, outcome_states)
  if (length(unknown)) {
    mostr_abort(
      "`absorb` must name states of the response, but ",
      format_sop_visit_list(paste0("\"", unknown, "\"")),
      if (length(unknown) == 1L) " is not one" else " are not",
      ". The states in the fitted rows of `data` are ",
      paste(state_column_labels_sorted(outcome_states), collapse = ", "),
      "."
    )
  }

  for (column in c(p_var, p2_var)) {
    values <- as_state_labels(data[[column]][rows])
    hit <- !is.na(values) & values %in% absorb_labels
    if (!any(hit)) {
      next
    }
    states <- unique(values[hit])
    what <- if (identical(column, p_var)) {
      "the previous state"
    } else {
      "the state two visits back"
    }
    mostr_abort(
      "`absorb` names ",
      if (length(states) == 1L) "the state " else "the states ",
      format_sop_visit_list(paste0("\"", states, "\"")),
      ", but ",
      sum(hit),
      " of the fitted rows of `data` have ",
      if (length(states) == 1L) "it" else "one of them",
      " as ",
      what,
      " (`",
      column,
      "`). A patient in an absorbing state can't leave it, so these rows ",
      "are not transitions the model can learn from. Remove them before ",
      "fitting, for example with `data[!data$",
      column,
      " %in% ",
      deparse1(absorb_labels),
      ", ]`; `prepare_markov_data(absorbing_state = )` does this for the ",
      "previous state."
    )
  }

  previous <- data[[p_var]][rows]
  previous_states <- unique(as_state_labels(previous[!is.na(previous)]))
  outcome_only <- setdiff(outcome_states, c(previous_states, absorb_labels))
  if (length(outcome_only)) {
    warn_outcome_only_states(outcome_only, p_var, absorb_labels)
  }
  absorb
}

# A state that occurs as an outcome but never as a previous state is often
# absorbing, such as death. If it is not named in `absorb`, a numeric previous
# state (for example `rcs(yprev, 4)`) extrapolates transitions out of it: in a
# simulated trial measured on 2026-11-12, dead patients moved back to state 7,
# and the day-14 probability of death was 0.012 instead of 0.057.
warn_outcome_only_states <- function(states, p_var, absorb_labels) {
  several <- length(states) > 1L
  quoted <- paste0("\"", states, "\"")
  suggestion <- deparse1(c(absorb_labels, states))
  warning(warningCondition(
    paste0(
      if (several) "The states " else "The state ",
      format_sop_visit_list(quoted),
      if (several) " occur" else " occurs",
      " as an outcome but never as the previous state (`",
      p_var,
      "`). ",
      if (several) "They" else "It",
      " may be absorbing, such as death; if so, name ",
      if (several) "them" else "it",
      " in `absorb`, for example `mostm(..., absorb = ",
      suggestion,
      ")`. Otherwise the model may predict transitions out of ",
      if (several) "them" else "it",
      ", so that patients who reach ",
      if (several) "one" else "it",
      " leave again."
    ),
    class = "mostr_warning_outcome_only_states",
    call = NULL
  ))
}

# State labels sorted as numbers when they all are numbers, otherwise as text.
state_column_labels_sorted <- function(labels) {
  numbers <- suppressWarnings(as.numeric(labels))
  if (!anyNA(numbers)) {
    return(labels[order(numbers)])
  }
  sort(labels)
}

# The gap column given to `mostm()`. The prediction functions set it to the
# difference between consecutive requested times (the first requested time
# for the first transition), so the fitted gaps must be the time since each
# patient's previous row, in the units of the time column. A gap in another
# unit, such as days with time in weeks, would be fitted without error and
# then predicted on the wrong scale.
validate_mostm_gap_var <- function(gap_var, data, id_var, time_var) {
  if (is.null(gap_var)) {
    return(invisible(NULL))
  }
  if (
    !is.character(gap_var) ||
      length(gap_var) != 1L ||
      is.na(gap_var) ||
      !nzchar(gap_var)
  ) {
    mostr_abort("`gap_var` must be a single column name.")
  }
  if (!gap_var %in% names(data)) {
    mostr_abort(
      "The column `",
      gap_var,
      "` given as `gap_var` is not in `data`. Set `gap_var` to the column ",
      "that holds the time since each patient's previous visit."
    )
  }
  gap <- data[[gap_var]]
  if (!is.numeric(gap)) {
    mostr_abort(
      "The gap column `",
      gap_var,
      "` must be numeric: it holds the time since each patient's previous ",
      "visit."
    )
  }
  time <- data[[time_var]]
  if (!is.numeric(time)) {
    mostr_abort(
      "`gap_var` needs numeric time, but `",
      time_var,
      "` is a ",
      if (is.factor(time)) "factor" else class(time)[[1L]],
      ". The prediction functions set the gap to the difference between ",
      "consecutive requested times, which visit labels don't give, and ",
      "there is no other way to give each visit's gap. Use a numeric time ",
      "column, or leave `gap_var` out."
    )
  }

  ord <- order(data[[id_var]], time)
  id <- data[[id_var]][ord]
  time <- time[ord]
  gap <- gap[ord]
  later <- c(FALSE, id[-1L] == id[-length(id)])
  elapsed <- c(NA_real_, diff(time))
  # Gaps computed from times in another way, for example as differences of
  # dates converted to days, can differ from `diff(time)` by rounding only.
  tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(elapsed))
  bad <- which(
    later &
      !is.na(gap) &
      !is.na(elapsed) &
      abs(gap - elapsed) > tolerance
  )
  if (!length(bad)) {
    return(invisible(NULL))
  }
  first <- bad[[1L]]
  mostr_abort(
    "The gap column `",
    gap_var,
    "` must hold the time since the patient's previous row, in the units of ",
    "`",
    time_var,
    "`, because the prediction functions set it to the difference between ",
    "consecutive requested times. It differs from that difference in ",
    length(bad),
    " of the rows after a patient's first: for example, patient ",
    format(id[[first]]),
    " at ",
    time_var,
    " ",
    format(time[[first]]),
    " has ",
    gap_var,
    " ",
    format(gap[[first]]),
    ", but the previous row is ",
    format(elapsed[[first]]),
    " earlier. Check that both columns use the same unit."
  )
}
