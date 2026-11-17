# Internal fitted-model plotting helpers.

markov_supported_model <- function(object) {
  inherits(object, c("orm", "blrm", "vglm", "robcov_vglm"))
}

# The settings of one diagnostic plot call. For a fitted model, the call's
# arguments are merged into the fit's settings (`resolve_markov_settings()`):
# an argument left `NULL` takes the fit's value, and one that contradicts the
# fit is an error. Trajectory data have no fit, so `NULL` means the default
# column name; `data_id_var` is the default patient-ID column for data.
plot_markov_settings <- function(
  x,
  id_var = NULL,
  time_var = NULL,
  p_var = NULL,
  p2_var = NULL,
  y_levels = NULL,
  data_id_var = NULL
) {
  if (markov_supported_model(x)) {
    return(resolve_markov_settings(
      x,
      id_var = id_var,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      y_levels = y_levels
    ))
  }
  markov_settings(
    id_var = id_var %||% data_id_var,
    time_var = time_var %||% "time",
    p_var = p_var %||% "yprev",
    p2_var = p2_var,
    y_levels = y_levels
  )
}

markov_model_ylevels <- function(model, y_levels = NULL) {
  if (!is.null(y_levels)) {
    return(y_levels)
  }

  if (inherits(model, "vglm")) {
    out <- tryCatch(model@extra$colnames.y, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  if (inherits(model, "robcov_vglm")) {
    out <- model$extra$colnames.y %||%
      tryCatch(model$vglm_fit@extra$colnames.y, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  if (inherits(model, "blrm")) {
    out <- model$ylevels
    if (!is.null(out)) {
      return(out)
    }
  }
  if (inherits(model, "orm")) {
    out <- model$yunique
    if (!is.null(out)) {
      return(out)
    }
  }

  mostr_abort("`y_levels` cannot be NULL.")
}

validate_plot_counterfactual_variables <- function(
  variables,
  require_two = FALSE
) {
  if (is.null(variables)) {
    return(NULL)
  }
  if (!is.list(variables) || is.null(names(variables))) {
    mostr_abort(
      "`variables` must be a named list, for example `list(tx = c(0, 1))`."
    )
  }
  if (length(variables) != 1L || !nzchar(names(variables)[1])) {
    mostr_abort("`variables` must contain exactly one named variable.")
  }
  if (length(variables[[1]]) < 1L) {
    mostr_abort("`variables[[1]]` must contain at least one value.")
  }
  if (require_two && length(variables[[1]]) < 2L) {
    mostr_abort("`variables[[1]]` must contain at least two values.")
  }
  if (anyDuplicated(as.character(variables[[1]]))) {
    mostr_abort("`variables[[1]]` must not contain duplicate values.")
  }
  variables
}

# The visits the diagnostic recursion runs, one transition per visit, so
# that sparse `times` such as `c(1, 7, 14)` on a daily fit still give the
# probabilities at visits 7 and 14. Visits that `times` leaves out are filled
# in from the fit's own visits, as `sops()` expects them (`R/sops-times.R`);
# filling in every whole number instead ran transitions the model was not
# fitted for (issue #51).
complete_plot_recursion_times <- function(
  times,
  time_info,
  model,
  time_var,
  gap_var = NULL,
  newdata_supplied = FALSE
) {
  if (isTRUE(time_info$is_factor)) {
    return(complete_plot_factor_times(times, time_info, model, newdata_supplied))
  }

  if (!is.numeric(times)) {
    # Character labels carry no order, so the requested order is the visit
    # order, as in `sops()`. Sorting them would run "six" before "three".
    return(unique(times))
  }

  times <- sort(unique(times))
  # With `gap_var`, each transition carries its own elapsed time, so the
  # requested times are the visits and nothing is filled in, as in `sops()`.
  if (!is.null(gap_var)) {
    return(times)
  }

  # Called with no times, this returns the fitted visits, or `NULL` when the
  # fit's data are not stored or its time column is not numeric. Without the
  # fitted visits there is nothing to fill in from.
  visits <- numeric_sop_visit_schedule(model, numeric(), time_var)$visits
  if (is.null(visits)) {
    return(times)
  }
  complete_plot_numeric_times(times, visits, model, newdata_supplied)
}

complete_plot_factor_times <- function(
  times,
  time_info,
  model,
  newdata_supplied
) {
  levels <- time_info$levels
  position <- match(as.character(times), levels)
  if (anyNA(position)) {
    mostr_abort("Could not align plot times with fitted visit levels.")
  }

  # Supplied rows hold the state before the first requested visit; stored
  # profiles hold the state before the first follow-up visit.
  start <- if (newdata_supplied) {
    min(position)
  } else {
    stored_sop_start_position(
      model,
      list(visits = levels, is_factor = TRUE),
      position
    )
  }
  if (!newdata_supplied && any(position < start)) {
    stop_plot_times_before_first_followup(
      levels[position[position < start]],
      levels[start]
    )
  }

  factor(
    levels[seq.int(start, max(position))],
    levels = levels,
    ordered = isTRUE(time_info$ordered)
  )
}

# `times` is sorted and unique; `visits` are the fitted visits, sorted.
complete_plot_numeric_times <- function(
  times,
  visits,
  model,
  newdata_supplied
) {
  position <- match_sop_visit_times(times, visits)
  below <- is.na(position) & times < visits[1L]
  above <- is.na(position) & times > visits[length(visits)]

  # A time between two fitted visits would add a transition the model was not
  # fitted for, and every later probability would be one transition off.
  between <- is.na(position) & !below & !above
  if (any(between)) {
    labels <- format_sop_time_labels(times[between])
    mostr_abort(
      "`times` includes ",
      format_sop_visit_list(labels),
      if (length(labels) == 1L) ", which is not a fitted visit" else
        ", which are not fitted visits",
      ". The fitted visits are ",
      format_sop_visit_list(format_sop_time_labels(visits)),
      ". The plot runs one transition per fitted visit, so a time between ",
      "fitted visits would add a transition the model was not fitted for. ",
      "If the visits are unevenly spaced, name the time since the previous ",
      "visit in `mostm(gap_var = )` and refit."
    )
  }

  requested <- position[!is.na(position)]
  if (newdata_supplied) {
    # Supplied rows hold the state before the first requested time, which may
    # come before the first fitted visit (kept as given) or be a later visit.
    first <- if (any(below)) 1L else min(c(requested, length(visits) + 1L))
  } else {
    first <- stored_sop_start_position(
      model,
      list(visits = visits, is_factor = FALSE),
      requested
    )
    early <- below | (!is.na(position) & position < first)
    if (any(early)) {
      stop_plot_times_before_first_followup(
        format_sop_time_labels(times[early]),
        format_sop_time_labels(visits[first])
      )
    }
  }
  last <- if (any(above)) length(visits) else max(c(requested, 0L))

  # Requested times replace the fitted visits they match, so that the plotted
  # times are found in the schedule by label even when they differ from a
  # fitted visit by rounding error.
  filled <- visits
  filled[requested] <- times[!is.na(position)]
  filled <- if (first <= last) filled[seq.int(first, last)] else numeric()

  sort(unique(c(times[below], filled, times[above])))
}

stop_plot_times_before_first_followup <- function(labels, first_followup) {
  mostr_abort(
    "`times` includes ",
    format_sop_visit_list(labels),
    ", before the first follow-up visit (",
    first_followup,
    ") of the stored starting profiles. These profiles hold each patient's ",
    "state before that visit, so the plot starts there. Request times from ",
    first_followup,
    " on, or supply `newdata` with each patient's state before the first ",
    "requested time."
  )
}
