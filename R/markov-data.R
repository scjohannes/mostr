#' Prepare Markov Data for Proportional Odds Modeling
#'
#' Prepares longitudinal Markov trajectory data for proportional odds regression
#' by removing absorbing states and converting variables to appropriate types.
#' The result can be fitted with [mostm()] and any of its engines.
#'
#' @param data A data frame containing Markov trajectory data with columns
#'   `y` (current state) and `yprev` (previous state).
#' @param absorbing_state Integer or integer vector. The absorbing state(s) to
#'   filter out (default: 6 for death). Observations where `yprev` equals any of
#'   these values are removed. Can be a single value (e.g., 6) or a vector of
#'   multiple absorbing states (e.g., c(1, 6) for both home discharge and death).
#' @param ordered_response Logical. Should `y` be converted to ordered factor?
#'   (default: TRUE). Required for vglm, not required for orm.
#' @param factor_previous Logical. Should `yprev` be converted to factor?
#'   Defaults to TRUE for backward compatibility and categorical previous-state
#'   effects. Set to FALSE when modeling the previous state linearly or with
#'   nonlinear numeric terms such as `rms::rcs(yprev, 6)`.
#'
#' @return A data frame with modified `y` and `yprev` columns and absorbing
#'   states removed.
#'
#' @details
#' This is a convenience function that performs standard data preparation for
#' proportional odds models on Markov data:
#' - Removes rows where patients were in an absorbing state at the previous time
#' - Converts current state to ordered factor (for cumulative models)
#' - Converts previous state to factor by default (for categorical transition
#'   effects)
#' - Can preserve previous state as numeric for linear or spline transition
#'   effects
#'
#' **Engines**: The `"vglm"` engine of [mostm()] needs an ordered response,
#' which `ordered_response = TRUE` provides (or write `ordered(y)` in the
#' formula). The `"orm"` and `"blrm"` engines also accept an integer response.
#'
#' **Multiple absorbing states**: When modeling scenarios with multiple absorbing
#' states (e.g., discharge home and death), specify them as a vector. For example,
#' `absorbing_state = c(1, 6)` will remove all transitions from either state 1
#' (home) or state 6 (death).
#'
#' @examples
#' \dontrun{
#' # Prepare data with single absorbing state (death only)
#' model_data <- prepare_markov_data(trajectories, absorbing_state = 6)
#'
#' # Prepare data with multiple absorbing states (home and death)
#' model_data <- prepare_markov_data(trajectories, absorbing_state = c(1, 6))
#'
#' # Fit a proportional-odds transition model with the default orm engine
#' fit <- mostm(
#'   y ~ tx + rms::rcs(time, 4) + yprev,
#'   data = model_data,
#'   id_var = "id"
#' )
#'
#' # Keep yprev numeric to use one linear or nonlinear previous-state effect
#' model_data <- prepare_markov_data(
#'   trajectories,
#'   absorbing_state = 6,
#'   factor_previous = FALSE
#' )
#' fit <- mostm(
#'   ordered(y) ~ tx + rms::rcs(time, 4) + rms::rcs(yprev, 6),
#'   data = model_data,
#'   id_var = "id",
#'   engine = "vglm"
#' )
#' }
#'
#' @export
prepare_markov_data <- function(
  data,
  absorbing_state = 6,
  ordered_response = TRUE,
  factor_previous = TRUE
) {
  if (!is.null(absorbing_state)) {
    data <- data[!data$yprev %in% absorbing_state, , drop = FALSE]
  }

  if (factor_previous) {
    data$yprev <- factor(data$yprev)
  }

  if (ordered_response) {
    data$y <- ordered(data$y)
  }

  data
}


#' Drop States Absent From a Bootstrap Resample
#'
#' A bootstrap resample can lack a state. Its refit then has no outcome
#' category, and no previous-state level, for that state. This function finds
#' the states absent from `data` and removes them from the factor levels of the
#' state columns. The states keep their original labels; they are never
#' renumbered.
#'
#' @param data A data frame containing the data to refit.
#' @param factor_cols Character vector of the column names holding states: the
#'   response and the previous-state columns. These may be factor, character
#'   or numeric state columns.
#' @param original_data Optional data frame containing the original data before
#'   subsetting. Used to determine which states exist in the full data.
#' @param y_levels Optional vector of original state labels in state order
#'   (e.g., 1:6 or `c("home", "ward", "dead")`). If provided, the labels of the
#'   states present are returned.
#' @param absorb Optional label(s) of the absorbing state(s). If provided, the
#'   labels of those present are returned.
#'
#' @return A list with components:
#'   - data: The data frame, with absent states removed from factor levels
#'   - ylevels: The labels of the states present, in state order (if
#'     `y_levels` was provided), or NULL
#'   - absorb: The labels of the absorbing states present (if `absorb` was
#'     provided), or NULL
#'   - missing_levels: Character vector of the absent state labels
#'
#' @details
#' States are matched by their labels (`as.character()` of the values), never
#' by converting labels to numbers, so letter labels such as `"a"`, `"b"`,
#' `"c"` work as numbers do. The state order is that of `y_levels`, followed
#' by any further labels of the original data in their own order (factor
#' level order; sorted values for numeric and character columns). A state is
#' present when any state column of `data` holds it.
#'
#' When a state is absent, each factor state column keeps only the levels it
#' holds in `data`, in their original order; an ordered factor stays ordered.
#' The refit needs no empty levels, and `rms::orm()` cannot fit with one.
#' Numeric and character columns are left unchanged, so a numeric previous
#' state keeps its numbers and a term such as `rcs(yprev, 3)` sees the same
#' scale as in the original fit. Because the labels are kept, the refit's
#' states, previous-state levels and coefficient names match the original
#' fit's, and supplied prediction profiles need no translation.
#'
#' When no state is absent, `data` is returned unchanged.
#'
#' @noRd
drop_absent_states <- function(
  data,
  factor_cols = c("y", "yprev"),
  original_data = NULL,
  y_levels = NULL,
  absorb = NULL
) {
  # All state labels in state order: `y_levels` first, then any other labels
  # of the original data (or, without either, of `data` itself).
  label_source <- original_data %||% (if (is.null(y_levels)) data)
  original_levels <- union(
    as_state_labels(y_levels),
    state_column_labels(label_source, factor_cols)
  )

  # Labels present in current data, kept in state order. A label unknown to
  # the original data is kept after the known ones.
  present <- unique(unlist(lapply(
    intersect(factor_cols, names(data)),
    function(col) {
      values <- data[[col]]
      as_state_labels(unique(values[!is.na(values)]))
    }
  )))
  states_present <- c(
    original_levels[original_levels %in% present],
    setdiff(present, original_levels)
  )

  missing_levels <- setdiff(original_levels, states_present)

  # If no levels are missing, return data as-is
  if (length(missing_levels) == 0) {
    return(list(
      data = data,
      ylevels = y_levels,
      absorb = absorb,
      missing_levels = character(0)
    ))
  }

  # Factor state columns lose their empty levels. `droplevels()` keeps the
  # level order and whether the factor is ordered. Other columns keep their
  # values: a numeric previous state must keep its scale.
  for (col in intersect(factor_cols, names(data))) {
    if (is.factor(data[[col]])) {
      data[[col]] <- droplevels(data[[col]])
    }
  }

  new_ylevels <- if (!is.null(y_levels)) states_present else NULL

  # The absorbing states present keep their labels.
  new_absorb <- if (!is.null(absorb)) {
    absorb_labels <- as_state_labels(absorb)
    absorb_present <- absorb_labels[absorb_labels %in% states_present]
    absorb_missing <- setdiff(absorb_labels, states_present)
    if (length(absorb_present) > 0) {
      if (length(absorb_missing) > 0) {
        warning(
          "Absorbing state(s) ",
          paste(absorb_missing, collapse = ", "),
          " not present in data and were dropped."
        )
      }
      absorb_present
    } else {
      warning(
        "Absorbing state(s) ",
        paste(absorb, collapse = ", "),
        " not present in data. Set to NULL."
      )
      NULL
    }
  } else {
    NULL
  }

  return(list(
    data = data,
    ylevels = new_ylevels,
    absorb = new_absorb,
    missing_levels = as.character(missing_levels)
  ))
}

# A factor of state labels whose levels keep the given state order; plain
# `factor()` would sort text labels such as "home", "ward", "dead".
state_factor <- function(states) {
  labels <- as_state_labels(states)
  factor(labels, levels = unique(labels))
}

# The state labels of the columns `cols` of `data`, in state order: factor
# levels in level order, and the sorted values of numeric or character
# columns. Numeric columns are sorted together, so that `y` and `yprev`
# holding different numbers still give the numbers in order.
state_column_labels <- function(data, cols) {
  cols <- intersect(cols, names(data))
  if (is.null(data) || length(cols) == 0L) {
    return(character(0))
  }
  columns <- lapply(cols, function(col) data[[col]])
  observed <- function(x) unique(x[!is.na(x)])
  if (all(vapply(columns, is.numeric, logical(1)))) {
    return(as_state_labels(sort(unique(unlist(lapply(columns, observed))))))
  }
  unique(unlist(lapply(columns, function(x) {
    if (is.factor(x)) levels(x) else as_state_labels(sort(observed(x)))
  })))
}
