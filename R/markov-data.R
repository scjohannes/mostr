#' Prepare Markov Data for Proportional Odds Modeling
#'
#' Prepares longitudinal Markov trajectory data for proportional odds regression
#' by removing the rows after an absorbing state and converting the state
#' columns to appropriate types. The result can be fitted with [mostm()] and
#' any of its engines.
#'
#' @param data A data frame containing Markov trajectory data, with a column
#'   for the current state (`y_var`) and one for the previous state (`p_var`).
#' @param absorb Required. The absorbing state or states, such as death, that
#'   a patient never leaves. Rows whose previous state is one of them are
#'   removed. Give one state (`absorb = 6`), several (`absorb = c(1, 6)` for
#'   both discharge home and death), or `NULL` to keep every row when the data
#'   have no absorbing state. There is no default because the absorbing states
#'   depend on the data: a default of `6` would silently remove rows from data
#'   in which state 6 is not death. Give the same states to
#'   `mostm(absorb = )`.
#' @param ordered_response Logical. Should the current state be converted to an
#'   ordered factor? (default: TRUE). Required for vglm, not required for orm.
#' @param factor_previous Logical. Should the previous state be converted to a
#'   factor? Defaults to TRUE for categorical previous-state effects. Set to
#'   FALSE when modeling the previous state linearly or with nonlinear numeric
#'   terms such as `rms::rcs(yprev, 6)`.
#' @param y_var Name of the column with the current state (default `"y"`).
#' @param p_var Name of the column with the previous state (default
#'   `"yprev"`), as in [mostm()].
#'
#' @return `data` with the rows after an absorbing state removed and the
#'   `y_var` and `p_var` columns converted as requested.
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
#' `absorb = c(1, 6)` will remove all transitions from either state 1
#' (home) or state 6 (death).
#'
#' @examples
#' # 200 simulated patients followed for 14 days. States 1 and 2 are outside
#' # hospital, states 3 to 7 in hospital, and state 8 is death. sim_markov()
#' # keeps the rows after death, with state 8 on every later day.
#' trajectories <- sim_markov("actt2", n_patients = 200, times = 1:14, seed = 1)
#'
#' # Death is absorbing: remove the rows that start in state 8
#' model_data <- prepare_markov_data(trajectories, absorb = 8)
#' c(before = nrow(trajectories), after = nrow(model_data))
#' str(model_data)
#'
#' # Fit a transition model to the prepared data
#' fit <- mostm(
#'   y ~ yprev + rms::rcs(time, 3) + tx,
#'   data = model_data,
#'   id_var = "id",
#'   absorb = "8"
#' )
#'
#' # Keep the previous state numeric, for one linear previous-state effect
#' # instead of one coefficient per state
#' numeric_data <- prepare_markov_data(
#'   trajectories,
#'   absorb = 8,
#'   factor_previous = FALSE
#' )
#' fit_linear <- mostm(
#'   y ~ yprev + rms::rcs(time, 3) + tx,
#'   data = numeric_data,
#'   id_var = "id",
#'   absorb = "8"
#' )
#'
#' # Several absorbing states: also treat leaving hospital (states 1 and 2)
#' # as final, as in an analysis of time to discharge
#' discharge_data <- prepare_markov_data(trajectories, absorb = c(1, 2, 8))
#' nrow(discharge_data)
#'
#' @export
prepare_markov_data <- function(
  data,
  absorb,
  ordered_response = TRUE,
  factor_previous = TRUE,
  y_var = "y",
  p_var = "yprev"
) {
  if (missing(absorb)) {
    mostr_abort(
      "`absorb` is required: give the absorbing state(s), such as ",
      "`absorb = 6` for death, or `absorb = NULL` to keep every row."
    )
  }
  # A missing previous-state column would make every row look non-absorbing
  # (`NULL %in% absorb` is empty) and drop them all, so check the names first.
  missing_cols <- setdiff(c(y_var, p_var), names(data))
  if (length(missing_cols) > 0) {
    mostr_abort(
      "Columns not in `data`: ",
      paste0("`", missing_cols, "`", collapse = ", "),
      ". Name the state columns with `y_var` and `p_var`."
    )
  }

  if (!is.null(absorb)) {
    data <- data[!data[[p_var]] %in% absorb, , drop = FALSE]
  }

  if (factor_previous) {
    data[[p_var]] <- factor(data[[p_var]])
  }

  if (ordered_response) {
    data[[y_var]] <- ordered(data[[y_var]])
  }

  data
}


#' Drop States Absent From a Bootstrap Resample
#'
#' A bootstrap resample can lack a state as an outcome, as a previous state,
#' or both. Its refit then has no outcome category, or no previous-state
#' level, for that state. This function finds the outcome states absent from
#' the response of `data` and removes, from each factor state column, the
#' levels that column does not hold. The states keep their original labels;
#' they are never renumbered.
#'
#' @param data A data frame containing the data to refit.
#' @param factor_cols Character vector of the column names holding states: the
#'   response and the previous-state columns. These may be factor, character
#'   or numeric state columns.
#' @param response_col The name of the response column, whose states are the
#'   refit's outcome states. It is treated as a state column even when it is
#'   not in `factor_cols`.
#' @param original_data Optional data frame containing the original data before
#'   subsetting. Used to determine which outcome states exist in the full data.
#' @param y_levels Optional vector of the original outcome state labels in
#'   state order (e.g., 1:6 or `c("home", "ward", "dead")`). If provided, the
#'   labels of the outcome states present are returned.
#' @param absorb Optional label(s) of the absorbing state(s). If provided, the
#'   labels of those present as outcomes are returned.
#'
#' @return A list with components:
#'   - data: The data frame, each factor state column without the levels it
#'     does not hold
#'   - ylevels: The labels of the outcome states present, in state order (if
#'     `y_levels` was provided), or NULL
#'   - absorb: The labels of the absorbing states present (if `absorb` was
#'     provided), or NULL
#'   - missing_levels: Character vector of the absent outcome state labels
#'   - levels_dropped: TRUE when some factor state column lost a level that
#'     it holds in `original_data` (or any level, without `original_data`),
#'     so that the refit has fewer coefficients than the original fit
#'
#' @details
#' States are matched by their labels (`as.character()` of the values), never
#' by converting labels to numbers, so letter labels such as `"a"`, `"b"`,
#' `"c"` work as numbers do. The outcome states of the original fit are
#' `y_levels`, followed by any further labels that the response of the
#' original data holds, in their own order (factor level order; sorted values
#' for numeric and character columns). An outcome state is present when the
#' response of `data` holds it. A state that only the previous-state columns
#' hold is a starting state, not an outcome state of the refit, just as in a
#' fit to the full data. The refit still predicts from it as a previous state
#' and gives it probability 0 as an outcome, its maximum-likelihood value in
#' that resample.
#'
#' Each factor state column keeps only the levels it holds in `data`, in their
#' original order; an ordered factor stays ordered. The response loses the
#' absent outcome states, and each previous-state column loses its own empty
#' levels. The refit needs no empty levels, and `rms::orm()` cannot fit with
#' one. Numeric and character columns are left unchanged, so a numeric
#' previous state keeps its numbers and a term such as `rcs(yprev, 3)` sees
#' the same scale as in the original fit. Because the labels are kept, the
#' refit's states, previous-state levels and coefficient names match the
#' original fit's, and supplied prediction profiles need no translation.
#'
#' When no factor state column has an empty level, `data` is returned
#' unchanged.
#'
#' @noRd
drop_absent_states <- function(
  data,
  factor_cols = c("y", "yprev"),
  response_col = "y",
  original_data = NULL,
  y_levels = NULL,
  absorb = NULL
) {
  if (!response_col %in% names(data)) {
    mostr_abort(
      "The bootstrap data do not contain the response column `",
      response_col,
      "`."
    )
  }

  # The outcome states of the original fit in state order: `y_levels` first,
  # then any other labels that the response of the original data (or, without
  # either, of `data` itself) holds. Labels held only by a previous-state
  # column are starting states, not outcome states, as in a fit to the full
  # data.
  label_source <- original_data %||% (if (is.null(y_levels)) data)
  original_levels <- union(
    as_state_labels(y_levels),
    observed_state_labels(label_source, response_col)
  )

  # Outcome labels present in the resample, kept in state order. A label
  # unknown to the original data is kept after the known ones.
  present <- observed_state_labels(data, response_col)
  states_present <- c(
    original_levels[original_levels %in% present],
    setdiff(present, original_levels)
  )
  missing_levels <- setdiff(original_levels, states_present)

  # Each factor state column loses the levels it does not hold: the response
  # its absent outcome states, a previous-state column the states that no
  # patient of the resample is observed after. `droplevels()` keeps the level
  # order and whether the factor is ordered. Other columns keep their values:
  # a numeric previous state must keep its scale.
  # `levels_dropped` records whether a column lost a level that the original
  # data hold, so that the refit has fewer coefficients than the original fit
  # (and the original coefficients can't start it). A level that was empty
  # in the original data as well changes nothing.
  levels_dropped <- FALSE
  for (col in intersect(union(response_col, factor_cols), names(data))) {
    values <- data[[col]]
    if (!is.factor(values)) {
      next
    }
    lost <- setdiff(levels(values), as_state_labels(unique(values)))
    if (length(lost) > 0L) {
      data[[col]] <- droplevels(values)
      if (is.null(original_data) ||
        any(lost %in% observed_state_labels(original_data, col))) {
        levels_dropped <- TRUE
      }
    }
  }

  if (length(missing_levels) == 0) {
    return(list(
      data = data,
      ylevels = y_levels,
      absorb = absorb,
      missing_levels = character(0),
      levels_dropped = levels_dropped
    ))
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
    missing_levels = as.character(missing_levels),
    levels_dropped = levels_dropped
  ))
}

# The state labels that column `col` of `data` holds, in state order (factor
# level order; sorted values for numeric and character columns). Unlike
# `state_column_labels()`, it leaves out factor levels that no row holds.
observed_state_labels <- function(data, col) {
  if (is.null(data) || !col %in% names(data)) {
    return(character(0))
  }
  values <- data[[col]]
  held <- as_state_labels(unique(values[!is.na(values)]))
  ordered_labels <- state_column_labels(data, col)
  ordered_labels[ordered_labels %in% held]
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
