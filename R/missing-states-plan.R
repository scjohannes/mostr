# Fitting when some states are not observed: reading the response and laying
# out the visits whose likelihood sums over paths.
#
# A visit's state is either observed exactly, known to lie in a range of
# states (`Ocens(lower, upper)`), or not observed at all (`NA`). A patient's
# likelihood is the product of ordinary transition probabilities, for visits
# whose state and previous state are both observed, and, for each run of
# visits where that is not so (a "stretch"), the sum over every sequence of
# states the patient could have passed through:
#
#   f = sum over (z_1, ..., z_T), z_t in range_t, of
#       P(z_1 | known previous state) * P(z_2 | z_1) * ... * P(z_T | z_{T-1})
#
# The sum is a product of transition matrices restricted to the ranges, so it
# is computed with one pass over the visits (`cpp_path_loglik()`), never by
# listing paths. The functions here only work out which rows form stretches
# and which design rows each step needs.

# The response of `formula` as a range of states for each row of `data`.
# Returns the state labels in model order and, for each row, the lowest and
# highest possible state (as positions in the labels). `lo == hi` is an
# observed state; both `NA` is a visit without any record of the state.
#
# mostr reads `Ocens()` itself instead of passing it to rms: `rms::orm()`
# 8.2-1 fits an infinite bound as if the state were the finite bound exactly
# (`Ocens(-Inf, 2)` as state 1, `Ocens(3, Inf)` as state 3; measured Nov 2026
# against the likelihood written out), and a range cannot be passed on anyway
# once the next visit's previous state is unknown.
markov_state_ranges <- function(formula, data, env) {
  lhs <- formula[[2L]]
  n <- nrow(data)
  if (markov_is_ocens_call(lhs)) {
    if (length(lhs) != 3L) {
      mostr_abort(
        "A range of states is written `Ocens(lower, upper)`, with the lowest ",
        "and the highest state the patient can have been in."
      )
    }
    lower <- eval(lhs[[2L]], data, env)
    upper <- eval(lhs[[3L]], data, env)
    if (length(lower) != n || length(upper) != n) {
      mostr_abort(
        "Both arguments of `Ocens()` must have one value per row of `data`."
      )
    }
    labels <- markov_state_labels(lower, upper)
    lo <- markov_match_state(lower, labels, infinite = 1L)
    hi <- markov_match_state(upper, labels, infinite = length(labels))
    one_sided <- xor(is.na(lo), is.na(hi))
    if (any(one_sided)) {
      mostr_abort(
        "`Ocens()` has a missing bound on ",
        sum(one_sided),
        " row",
        if (sum(one_sided) > 1L) "s",
        " where the other bound is given. Use `-Inf` or `Inf` for an open ",
        "end, and `NA` in both bounds for a visit whose state is unknown."
      )
    }
    if (any(lo > hi, na.rm = TRUE)) {
      mostr_abort(
        "`Ocens(lower, upper)` needs `lower <= upper`; ",
        sum(lo > hi, na.rm = TRUE),
        " rows have the bounds the other way round."
      )
    }
    return(list(labels = labels, lo = lo, hi = hi, ocens = TRUE))
  }

  y <- eval(lhs, data, env)
  if (NCOL(y) != 1L || length(y) != n) {
    mostr_abort(
      "The response must be one state per row of `data`, for example ",
      "`ordered(y)`, or a range of states, `Ocens(lower, upper)`."
    )
  }
  labels <- markov_state_labels(y)
  position <- markov_match_state(y, labels)
  list(labels = labels, lo = position, hi = position, ocens = FALSE)
}

markov_is_ocens_call <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  fun <- expr[[1L]]
  if (is.call(fun) && length(fun) == 3L && is.name(fun[[3L]])) {
    fun <- fun[[3L]]
  }
  is.name(fun) && identical(as.character(fun), "Ocens")
}

# State labels in the order the engines use: the levels of a factor, else the
# sorted distinct values. Infinite bounds are not states.
markov_state_labels <- function(...) {
  values <- list(...)
  if (all(vapply(values, is.factor, logical(1)))) {
    level_sets <- lapply(values, levels)
    if (!all(vapply(level_sets, identical, logical(1), level_sets[[1L]]))) {
      mostr_abort(
        "The two bounds of `Ocens()` must be factors with the same levels."
      )
    }
    present <- unique(unlist(lapply(values, function(v) as.character(v))))
    return(level_sets[[1L]][level_sets[[1L]] %in% present])
  }
  if (any(vapply(values, is.factor, logical(1)))) {
    mostr_abort(
      "The two bounds of `Ocens()` must both be factors or both be numbers."
    )
  }
  all_values <- unlist(values, use.names = FALSE)
  if (is.numeric(all_values)) {
    all_values <- all_values[is.finite(all_values)]
  }
  as.character(sort(unique(all_values[!is.na(all_values)])))
}

markov_match_state <- function(values, labels, infinite = NA_integer_) {
  position <- match(as.character(values), labels)
  if (is.numeric(values) && !is.na(infinite)) {
    position[is.infinite(values)] <- infinite
  }
  position
}

# TRUE when some row's state or previous state is not observed, so that the
# fit has to sum over paths. `rows` are the rows `subset` keeps.
markov_has_unobserved_states <- function(ranges, previous, rows) {
  lo <- ranges$lo[rows]
  hi <- ranges$hi[rows]
  anyNA(lo) || any(lo != hi, na.rm = TRUE) || anyNA(previous[rows])
}

# Lay out the stretches of `data` (the rows `subset` keeps).
#
# `ranges` comes from `markov_state_ranges()` for the same rows, `absorb` holds
# positions of the absorbing states, `weights` the row weights or NULL.
# `design` turns rows of `data` into design-matrix rows of the model.
#
# The result holds what `cpp_path_loglik()` needs. Design rows and whole
# stretches that occur more than once are kept once: two patients with the
# same covariates, the same known state before the gap and the same record
# after it have the same likelihood term, which is computed once and counted
# with their summed weight.
markov_path_plan <- function(
  data,
  ranges,
  absorb,
  id_var,
  time_var,
  p_var,
  weights,
  design
) {
  n <- nrow(data)
  labels <- ranges$labels
  n_states <- length(labels)
  lo <- ranges$lo
  hi <- ranges$hi

  # An unrecorded visit is known not to be in an absorbing state: these are
  # recorded when they happen (a death has a date), and the data hold no rows
  # after them. So `NA` means "one of the other states".
  unrecorded <- is.na(lo)
  if (any(unrecorded)) {
    possible <- setdiff(seq_len(n_states), absorb)
    if (!identical(possible, seq.int(possible[[1L]], possible[[length(possible)]]))) {
      mostr_abort(
        "A visit with a missing state can be in any state that is not ",
        "absorbing, and here those states are not next to each other in the ",
        "order of the states (",
        paste(labels[possible], collapse = ", "),
        "). Give the possible states of each such visit as a range with ",
        "`Ocens(lower, upper)`."
      )
    }
    lo[unrecorded] <- possible[[1L]]
    hi[unrecorded] <- possible[[length(possible)]]
  }
  exact <- lo == hi

  previous <- match(as.character(data[[p_var]]), labels)
  unknown_previous_label <- is.na(previous) & !is.na(data[[p_var]])
  previous_known <- !is.na(data[[p_var]])

  # Visits in time order within patient. Other time types keep the order of
  # the rows.
  time <- data[[time_var]]
  id <- as.character(data[[id_var]])
  ord <- if (is.numeric(time) || is.factor(time)) {
    order(match(id, unique(id)), time)
  } else {
    order(match(id, unique(id)))
  }
  s_id <- id[ord]
  s_exact <- exact[ord]
  s_prev_known <- previous_known[ord]
  same_patient <- c(FALSE, s_id[-1L] == s_id[-n])
  # A row continues a stretch when the visit before it, of the same patient,
  # has no exactly observed state.
  continues <- same_patient & c(FALSE, !s_exact[-n])

  # The previous state of a row must be missing exactly when the visit before
  # it is not observed exactly; anything else means the rows are not
  # consecutive visits or `p_var` was filled in some other way.
  filled <- continues & s_prev_known
  if (any(filled)) {
    mostr_abort(
      "`",
      p_var,
      "` is recorded on ",
      sum(filled),
      " row",
      if (sum(filled) > 1L) "s",
      " whose previous visit has no exactly observed state (first: `",
      id_var,
      "` = ",
      s_id[which(filled)[[1L]]],
      "). Set `",
      p_var,
      "` to `NA` on the row after each visit with a missing state or a ",
      "range of states: the model sums over the states the patient can have ",
      "been in, so a filled-in value, such as the last observed state, would ",
      "be wrong."
    )
  }
  dangling <- !continues & !s_prev_known & same_patient
  if (any(dangling)) {
    mostr_abort(
      "`",
      p_var,
      "` is missing on ",
      sum(dangling),
      " row",
      if (sum(dangling) > 1L) "s",
      " although the patient's previous row has an observed state (first: `",
      id_var,
      "` = ",
      s_id[which(dangling)[[1L]]],
      "). Each patient needs one row per visit, in time order, with `",
      p_var,
      "` equal to the state on the row before; add a row with a missing ",
      "state for every visit that was not observed, or use ",
      "`missing_states = \"drop\"` to leave out all rows with a missing value."
    )
  }

  in_stretch <- continues | !s_exact
  starts <- in_stretch & !continues
  run <- cumsum(starts)
  run[!in_stretch] <- 0L
  # A patient's first rows can't be used until a previous state is known:
  # the model has no distribution for a state without a previous state.
  orphan_runs <- run[starts & !s_prev_known]
  in_stretch[run %in% orphan_runs] <- FALSE

  # Steps at the end of a stretch that allow every state say nothing: the
  # probabilities of all states sum to 1.
  uninformative <- in_stretch & lo[ord] == 1L & hi[ord] == n_states
  if (any(uninformative)) {
    last_informative <- tapply(
      seq_len(n)[in_stretch & !uninformative],
      run[in_stretch & !uninformative],
      max
    )
    cut_after <- last_informative[as.character(run)]
    cut_after[is.na(cut_after)] <- 0L
    in_stretch <- in_stretch & seq_len(n) <= cut_after
  }

  step_position <- which(in_stretch)
  step_data_row <- ord[step_position]
  n_steps <- length(step_position)
  step_run <- match(run[step_position], unique(run[step_position]))
  n_instances <- if (n_steps) max(step_run) else 0L
  is_first <- !duplicated(step_run)

  if (any(unknown_previous_label[step_data_row[is_first]])) {
    mostr_abort(
      "`",
      p_var,
      "` holds values that are not states of the response."
    )
  }

  step_lo <- lo[step_data_row]
  step_hi <- hi[step_data_row]
  step_weight <- if (is.null(weights)) rep(1, n_steps) else weights[step_data_row]
  if (n_steps) {
    weight_range <- tapply(step_weight, step_run, function(w) diff(range(w)))
    if (any(weight_range != 0)) {
      mostr_abort(
        "Row weights must be the same on all rows of a run of visits with ",
        "unobserved states, because the model has one likelihood term for ",
        "the whole run. Use one weight per patient."
      )
    }
  }

  # Design rows: a first step uses its own row; a later step needs one row
  # for every state the patient can have been in at the step before, except
  # absorbing states, which the patient can't leave.
  later <- which(!is_first)
  candidate_step <- integer()
  candidate_state <- integer()
  if (length(later)) {
    before_lo <- step_lo[later - 1L]
    before_hi <- step_hi[later - 1L]
    width <- before_hi - before_lo + 1L
    candidate_step <- rep(later, width)
    candidate_state <- sequence(width, from = before_lo)
    keep <- !candidate_state %in% absorb
    candidate_step <- candidate_step[keep]
    candidate_state <- candidate_state[keep]
  }
  first_steps <- which(is_first)
  expanded <- data[
    c(step_data_row[first_steps], step_data_row[candidate_step]),
    ,
    drop = FALSE
  ]
  if (length(candidate_step)) {
    rows <- length(first_steps) + seq_along(candidate_step)
    expanded[[p_var]][rows] <- markov_previous_state_values(
      labels[candidate_state],
      data[[p_var]]
    )
  }
  x_all <- if (nrow(expanded)) design(expanded) else NULL
  if (!is.null(x_all) && anyNA(x_all)) {
    bad <- unique(c(step_data_row[first_steps], step_data_row[candidate_step])[
      rowSums(is.na(x_all)) > 0
    ])
    mostr_abort(
      "A visit with an unobserved state needs every predictor of the model, ",
      "because the model predicts the transitions into and out of it. ",
      length(bad),
      " such row",
      if (length(bad) > 1L) "s have" else " has",
      " a missing predictor (first: `",
      id_var,
      "` = ",
      id[bad[[1L]]],
      "). Fill in the predictors, or remove the patient's rows from that ",
      "visit on."
    )
  }

  if (!n_steps) {
    return(markov_empty_path_plan(labels, absorb))
  }

  # Keep each distinct design row once.
  rank <- data.table::frankv(
    data.table::as.data.table(x_all),
    ties.method = "dense"
  )
  design_row <- match(rank, unique(rank))
  x_gap <- x_all[!duplicated(design_row), , drop = FALSE]
  rownames(x_gap) <- NULL

  first_row <- design_row[seq_along(first_steps)]
  step_rows <- matrix(0L, n_steps, n_states)
  if (length(candidate_step)) {
    step_rows[cbind(candidate_step, candidate_state)] <-
      design_row[length(first_steps) + seq_along(candidate_step)]
  }

  # Keep each distinct stretch once: same design rows, same ranges.
  step_key <- do.call(
    paste,
    c(list(step_lo, step_hi), asplit(step_rows, 2L), sep = ",")
  )
  instance_key <- paste(
    first_row,
    vapply(split(step_key, step_run), paste, character(1), collapse = "|")
  )
  stretch <- match(instance_key, unique(instance_key))
  representative <- which(!duplicated(stretch))
  instance_steps <- split(seq_len(n_steps), step_run)
  kept_steps <- unlist(instance_steps[representative], use.names = FALSE)
  n_per <- lengths(instance_steps[representative])
  instance_weight <- step_weight[first_steps]

  list(
    labels = labels,
    absorb = as.integer(absorb),
    x_gap = x_gap,
    stretch_first_row = as.integer(first_row[representative]),
    stretch_first_step = as.integer(cumsum(c(1L, n_per[-length(n_per)]))),
    stretch_n_steps = as.integer(n_per),
    stretch_weight = as.numeric(rowsum(instance_weight, stretch)[, 1L]),
    step_lo = as.integer(step_lo[kept_steps]),
    step_hi = as.integer(step_hi[kept_steps]),
    step_rows = step_rows[kept_steps, , drop = FALSE],
    instances = data.frame(
      stretch = stretch,
      id = id[step_data_row[first_steps]],
      weight = instance_weight,
      stringsAsFactors = FALSE
    ),
    step_data_rows = step_data_row,
    n_steps_total = n_steps
  )
}

markov_empty_path_plan <- function(labels, absorb) {
  list(
    labels = labels,
    absorb = as.integer(absorb),
    x_gap = NULL,
    stretch_first_row = integer(),
    stretch_first_step = integer(),
    stretch_n_steps = integer(),
    stretch_weight = numeric(),
    step_lo = integer(),
    step_hi = integer(),
    step_rows = matrix(0L, 0L, length(labels)),
    instances = data.frame(
      stretch = integer(),
      id = character(),
      weight = numeric(),
      stringsAsFactors = FALSE
    ),
    step_data_rows = integer(),
    n_steps_total = 0L
  )
}

# State labels as values of the previous-state column, which the model reads
# as a factor, a number or text.
markov_previous_state_values <- function(labels, prototype) {
  if (is.factor(prototype)) {
    missing_levels <- setdiff(labels, levels(prototype))
    if (length(missing_levels)) {
      mostr_abort(
        "The previous-state column has no level for the state",
        if (length(missing_levels) > 1L) "s " else " ",
        paste(missing_levels, collapse = ", "),
        ", which a patient can have been in at an unobserved visit."
      )
    }
    return(factor(
      labels,
      levels = levels(prototype),
      ordered = is.ordered(prototype)
    ))
  }
  if (is.numeric(prototype)) {
    return(as.numeric(labels))
  }
  labels
}
