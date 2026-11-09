# Crossed cut-points in partial proportional-odds models.
#
# In a partial proportional-odds model some predictors have a separate effect
# at each cut-point (threshold), so for some covariate values the cut-points
# can cross. The cumulative probabilities P(Y >= k) then stop decreasing in k
# and the model gives a negative probability for a state. State occupancy
# probabilities (SOPs) are not defined there, so every route that turns linear
# predictors into state probabilities stops with a
# `mostr_crossed_thresholds` error instead of clipping and rescaling.
#
# Rounding alone makes a probability negative by a few multiples of 1e-16.
# Values down to -`crossed_threshold_tolerance()` are treated as rounding: they
# are set to zero and the row is rescaled, as before. The same tolerance is used
# by the analytical (delta-method) recursion in src/sops.cpp, and by the
# compiled kernels (`kCrossedTolerance`), which must stay equal to it.
crossed_threshold_tolerance <- function() {
  1e-10
}

# Stop when `probs` holds a probability below the tolerance.
#
# `probs` is a prediction matrix [row, state] or a posterior array
# [draw, row, state]. Prediction rows stack `n_profiles` profiles once per
# previous state (or pair of previous states), so row r belongs to profile
# (r - 1) %% n_profiles + 1. `rows` gives the stacked row number of each row of
# `probs` when only some rows were predicted. `time` is the visit's time label.
check_crossed_thresholds <- function(probs, n_profiles, time, rows = NULL) {
  tolerance <- crossed_threshold_tolerance()
  crossed_cells <- !is.na(probs) & probs < -tolerance
  if (!any(crossed_cells)) {
    return(invisible(NULL))
  }
  dims <- dim(probs)
  draws <- NULL
  n_draws <- NULL
  if (length(dims) == 3L) {
    crossed <- apply(crossed_cells, c(1L, 2L), any)
    draws <- sum(rowSums(crossed) > 0)
    n_draws <- dims[1L]
    crossed_rows <- which(colSums(crossed) > 0)
    rows <- rows %||% seq_len(dims[2L])
  } else {
    crossed_rows <- which(rowSums(crossed_cells) > 0)
    rows <- rows %||% seq_len(nrow(probs))
  }
  profiles <- unique((rows[crossed_rows] - 1L) %% n_profiles + 1L)
  stop_crossed_thresholds(
    profiles = length(profiles),
    n_profiles = n_profiles,
    time = time,
    minimum = min(probs[crossed_cells]),
    draws = draws,
    n_draws = n_draws
  )
}

# Stop when a compiled update kernel reported crossed cut-points. The kernels
# mark their output with the profiles (1-based) and the lowest probability.
check_native_crossed_thresholds <- function(out, n_profiles, time) {
  profiles <- attr(out, "crossed_profiles", exact = TRUE)
  if (is.null(profiles)) {
    return(out)
  }
  stop_crossed_thresholds(
    profiles = length(profiles),
    n_profiles = n_profiles,
    time = time %||% "unknown",
    minimum = attr(out, "minimum_probability", exact = TRUE)
  )
}

stop_crossed_thresholds <- function(
  profiles,
  n_profiles,
  time,
  minimum,
  draws = NULL,
  n_draws = NULL
) {
  where <- paste0(
    profiles,
    " of the ",
    n_profiles,
    if (n_profiles == 1) " profile" else " profiles",
    " at time ",
    as.character(time)
  )
  if (!is.null(draws)) {
    where <- paste0(where, " in ", draws, " of the ", n_draws, " posterior draws")
  }
  message <- paste0(
    "The fitted cut-points (thresholds) of the partial proportional-odds ",
    "model cross for ",
    where,
    ", so the model gives a negative probability for a state (as low as ",
    format(minimum, digits = 3),
    "). State occupancy probabilities aren't defined there. ",
    crossed_threshold_advice()
  )
  stop(structure(
    list(
      message = message,
      call = NULL,
      profiles = profiles,
      n_profiles = n_profiles,
      time = time,
      minimum_probability = minimum,
      draws = draws,
      n_draws = n_draws
    ),
    class = c("mostr_crossed_thresholds", "error", "condition")
  ))
}

crossed_threshold_advice <- function() {
  paste0(
    "Check the covariate values of those profiles, which are often outside ",
    "the range of the data, or simplify the model's non-proportional ",
    "(partial proportional-odds) terms."
  )
}

# Simulation and bootstrap loops turn a failed draw into a warning and drop it.
# Crossed cut-points must not be dropped silently: call this first in such an
# error handler so the error reaches the loop driver, which counts the draws.
rethrow_crossed_thresholds <- function(condition) {
  if (inherits(condition, "mostr_crossed_thresholds")) {
    stop(condition)
  }
  invisible(NULL)
}

# Run one draw of a simulation or bootstrap loop, recording crossed cut-points
# instead of stopping, so the driver can report how many draws are affected.
run_draw_recording_crossings <- function(fn, ...) {
  tryCatch(
    fn(...),
    mostr_crossed_thresholds = function(e) {
      structure(list(condition = e), class = "mostr_crossed_threshold_draw")
    }
  )
}

# Stop when any draw of a simulation or bootstrap loop had crossed cut-points.
# `unit` names the draws in the plural, for example "simulation draws".
stop_crossed_threshold_draws <- function(results, unit, unit_one) {
  crossed <- which(vapply(
    results,
    inherits,
    logical(1),
    what = "mostr_crossed_threshold_draw"
  ))
  if (length(crossed) == 0L) {
    return(invisible(results))
  }
  first <- results[[crossed[1L]]]$condition
  message <- paste0(
    "The fitted cut-points (thresholds) of the partial proportional-odds ",
    "model cross in ",
    length(crossed),
    " of the ",
    length(results),
    " ",
    unit,
    ", so those ",
    unit,
    " give a negative probability for a state. State occupancy ",
    "probabilities aren't defined there, so no confidence intervals were ",
    "computed. In the first such ",
    unit_one,
    " (number ",
    crossed[1L],
    "), the cut-points cross for ",
    first$profiles,
    " of the ",
    first$n_profiles,
    if (first$n_profiles == 1) " profile" else " profiles",
    " at time ",
    as.character(first$time),
    ". ",
    crossed_threshold_advice()
  )
  stop(structure(
    list(
      message = message,
      call = NULL,
      draws = length(crossed),
      n_draws = length(results),
      first = first
    ),
    class = c("mostr_crossed_thresholds", "error", "condition")
  ))
}
