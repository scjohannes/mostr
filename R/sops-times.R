# Checks that requested SOP times follow the fitted visit schedule.

# The SOP recursion takes one transition per value of `times`. A request that
# skips fitted visits, for example `times = c(1, 4, 8)` with daily visits, runs
# three transitions and labels the third result 8, although the probability at
# visit 8 needs eight transitions. The user-facing SOP functions warn about such
# requests. With `gap_var`, each transition carries its own elapsed time, so
# irregular `times` are legitimate and are not checked.

warn_skipped_sop_visits <- function(
  model,
  data,
  times,
  time_var,
  gap_var,
  newdata_supplied
) {
  issue <- skipped_sop_visit_issue(
    model,
    data,
    times,
    time_var,
    gap_var,
    newdata_supplied
  )
  if (!is.null(issue)) {
    warning(warningCondition(
      skipped_sop_visit_message(issue),
      class = "mostr_warning_skipped_visits",
      call = NULL
    ))
  }
  invisible(NULL)
}

# Internal SOP calls whose request was already checked use this so the
# user sees the warning once.
muffle_skipped_sop_visits <- function(expr) {
  withCallingHandlers(
    expr,
    mostr_warning_skipped_visits = function(w) {
      invokeRestart("muffleWarning")
    }
  )
}

# Returns `NULL` when the request is fine or the fitted visit schedule is
# unknown, otherwise a list describing what is wrong with `times`.
skipped_sop_visit_issue <- function(
  model,
  data,
  times,
  time_var,
  gap_var,
  newdata_supplied
) {
  if (
    !is.null(gap_var) ||
      is.null(time_var) ||
      !length(times) ||
      anyNA(times)
  ) {
    return(NULL)
  }

  time_info <- get_time_info(model, data, time_var)
  schedule <- if (isTRUE(time_info$is_factor)) {
    factor_sop_visit_schedule(time_info$levels, times)
  } else {
    numeric_sop_visit_schedule(model, times, time_var)
  }
  if (is.null(schedule)) {
    return(NULL)
  }

  visits <- schedule$visits
  position <- schedule$position
  requested <- position[!is.na(position)]

  # Supplied rows may hold the state before a later visit, so their start is
  # the earliest requested visit. Stored profiles start at the first
  # follow-up visit.
  start <- if (newdata_supplied) {
    if (length(requested)) min(requested) else NA_integer_
  } else {
    stored_sop_start_position(model, schedule, requested)
  }

  unordered <- is.unsorted(
    if (schedule$is_factor) position else as.numeric(times),
    strictly = TRUE
  )
  early <- !newdata_supplied &&
    (any(schedule$below) || (!is.na(start) && any(requested < start)))

  end <- if (any(schedule$above)) {
    length(visits)
  } else if (length(requested)) {
    max(requested)
  } else {
    NA_integer_
  }
  skipped <- if (!is.na(start) && !is.na(end) && end >= start) {
    setdiff(seq.int(start, end), requested)
  } else {
    integer()
  }

  if (!unordered && !early && !length(skipped)) {
    return(NULL)
  }

  list(
    times = times,
    schedule = schedule,
    start = start,
    end = end,
    unordered = unordered,
    early = early,
    skipped = skipped
  )
}

factor_sop_visit_schedule <- function(levels, times) {
  visits <- as.character(levels)
  position <- match(as.character(times), visits)
  # Unknown levels are reported by `resolve_sop_times()`.
  if (anyNA(position)) {
    return(NULL)
  }
  outside <- rep(FALSE, length(times))
  list(
    visits = visits,
    is_factor = TRUE,
    position = position,
    below = outside,
    above = outside
  )
}

numeric_sop_visit_schedule <- function(model, times, time_var) {
  if (!is.numeric(times) || any(!is.finite(times))) {
    return(NULL)
  }
  fit_data <- markov_model_refit_data(model)
  if (
    !is.data.frame(fit_data) ||
      !time_var %in% names(fit_data) ||
      !is.numeric(fit_data[[time_var]])
  ) {
    return(NULL)
  }
  values <- fit_data[[time_var]]
  visits <- sort(unique(values[is.finite(values)]))
  if (!length(visits)) {
    return(NULL)
  }

  times <- as.numeric(times)
  position <- match_sop_visit_times(times, visits)
  below <- is.na(position) & times < visits[1L]
  above <- is.na(position) & times > visits[length(visits)]
  # A time inside the schedule's range that is not a fitted visit means the
  # request does not follow the fitted schedule, so skipped visits cannot be
  # identified.
  if (any(is.na(position) & !below & !above)) {
    return(NULL)
  }

  list(
    visits = visits,
    is_factor = FALSE,
    position = position,
    below = below,
    above = above
  )
}

# Position of each time among the fitted visits, allowing for rounding error;
# `NA` when a time is not a fitted visit.
match_sop_visit_times <- function(times, visits) {
  vapply(
    as.numeric(times),
    function(x) {
      tolerance <- sqrt(.Machine$double.eps) * max(1, abs(x))
      hit <- which(abs(visits - x) <= tolerance)
      if (length(hit)) hit[1L] else NA_integer_
    },
    integer(1),
    USE.NAMES = FALSE
  )
}

stored_sop_start_position <- function(model, schedule, requested) {
  first <- markov_model_starting_profile_metadata(model)$first_followup_time
  known <- length(first) == 1L && !is.na(first)

  if (schedule$is_factor) {
    position <- if (known) {
      match(as.character(first), schedule$visits)
    } else {
      NA_integer_
    }
    if (is.na(position)) min(requested) else position
  } else {
    position <- if (known && is.numeric(first)) {
      match_sop_visit_times(first, schedule$visits)
    } else {
      NA_integer_
    }
    if (is.na(position)) 1L else position
  }
}

skipped_sop_visit_message <- function(issue) {
  schedule <- issue$schedule
  visits <- format_sop_time_labels(schedule$visits)
  times <- issue$times
  n_times <- length(times)
  last_label <- format_sop_time_labels(times[n_times])

  problems <- c(
    if (issue$unordered) "`times` is not in increasing order.",
    if (issue$early) {
      first_requested <- if (schedule$is_factor) {
        visits[min(schedule$position)]
      } else {
        format_sop_time_labels(min(times))
      }
      paste0(
        "`times` starts at ",
        first_requested,
        ", before the first follow-up visit (",
        visits[issue$start],
        ") of the stored starting profiles."
      )
    },
    if (length(issue$skipped)) {
      paste0(
        "`times` skips fitted ",
        if (length(issue$skipped) == 1L) "visit " else "visits ",
        format_sop_visit_list(visits[issue$skipped]),
        "."
      )
    }
  )

  steps <- if (issue$early || length(issue$skipped)) {
    paste0(
      "SOPs are calculated one transition per value of `times`, so the ",
      "result labelled ",
      last_label,
      " is the probability after ",
      n_times,
      if (n_times == 1L) " transition" else " transitions",
      " from the starting state, not at visit ",
      last_label,
      "."
    )
  } else {
    "SOPs are calculated one transition per value of `times`, in the order given."
  }

  paste(
    c(
      problems,
      steps,
      paste0(
        "To get SOPs at selected visits, ",
        skipped_sop_visit_remedy(issue),
        " and keep the rows you need."
      ),
      paste0(
        "If the model includes the time elapsed since the previous visit, ",
        "name it in `gap_var`."
      )
    ),
    collapse = " "
  )
}

# Describes the consecutive request that runs from the starting visit to the
# last requested time.
skipped_sop_visit_remedy <- function(issue) {
  schedule <- issue$schedule
  visits <- schedule$visits

  if (schedule$is_factor) {
    first <- visits[issue$start]
    last <- visits[max(issue$start, issue$end)]
    if (identical(first, last)) {
      return(paste0("request visit \"", first, "\""))
    }
    return(paste0(
      "request every fitted visit from \"",
      first,
      "\" to \"",
      last,
      "\""
    ))
  }

  times <- as.numeric(issue$times)
  first <- if (is.na(issue$start)) min(times) else visits[issue$start]
  last <- max(first, times)
  span <- visits[visits >= first]
  consecutive <- all(is_whole_sop_time(span)) &&
    all(abs(diff(span) - 1) <= sqrt(.Machine$double.eps))

  if (consecutive && is_whole_sop_time(first) && is_whole_sop_time(last)) {
    request <- if (first == last) {
      format_sop_time_labels(first)
    } else {
      paste0(
        format_sop_time_labels(first),
        ":",
        format_sop_time_labels(last)
      )
    }
    return(paste0("request every visit (`times = ", request, "`)"))
  }

  if (first == last) {
    return(paste0("request visit ", format_sop_time_labels(first)))
  }
  last_visit <- visits[length(visits)]
  if (first > last_visit) {
    return("request the times in increasing order")
  }
  paste0(
    "request every fitted visit from ",
    format_sop_time_labels(first),
    " to ",
    format_sop_time_labels(min(last, last_visit)),
    if (any(issue$schedule$above)) ", then the later times"
  )
}

is_whole_sop_time <- function(x) {
  abs(x - round(x)) <= sqrt(.Machine$double.eps) * pmax(1, abs(x))
}

format_sop_time_labels <- function(x) {
  if (!is.numeric(x)) {
    return(as.character(x))
  }
  vapply(
    x,
    function(value) {
      format(value, trim = TRUE, digits = 10, scientific = FALSE)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# "2", "2 and 3", "2, 3, and 5", or the first `limit` labels and a count.
format_sop_visit_list <- function(labels, limit = 5L) {
  n <- length(labels)
  if (n <= 1L) {
    return(paste(labels, collapse = ""))
  }
  if (n == 2L) {
    return(paste(labels, collapse = " and "))
  }
  if (n > limit) {
    return(paste0(
      paste(labels[seq_len(limit)], collapse = ", "),
      ", and ",
      n - limit,
      " more"
    ))
  }
  paste0(paste(labels[-n], collapse = ", "), ", and ", labels[n])
}
