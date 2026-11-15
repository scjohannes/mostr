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

# Lay out the stretches of `data` and build their design rows.
#
# `ranges` comes from `markov_state_ranges()` for the same rows, `absorb` holds
# positions of the absorbing states, `weights` the row weights or NULL.
# `design` turns rows of `data` into design-matrix rows of the model; without
# it only the layout is returned (`markov_path_layout()`), which tells whether
# any stretch is left.
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
  design = NULL,
  additive = TRUE,
  cluster = NULL
) {
  layout <- markov_path_layout(
    data = data,
    ranges = ranges,
    absorb = absorb,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    weights = weights
  )
  if (is.null(design)) {
    return(layout)
  }
  markov_path_plan_design(layout, design, additive = additive, cluster = cluster)
}

# Which rows form stretches, without any design rows. The fit needs this
# before the engine's fit (to know whether there is a stretch at all) and the
# design rows only after it (they come from the fitted model), so the two
# halves are separate and the layout is computed once per fit.
markov_path_layout <- function(
  data,
  ranges,
  absorb,
  id_var,
  time_var,
  p_var,
  weights
) {
  n <- nrow(data)
  labels <- ranges$labels
  n_states <- length(labels)
  lo <- ranges$lo
  hi <- ranges$hi

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
  same_patient <- c(FALSE, s_id[-1L] == s_id[-n])

  # A visit without any record of the state (`NA`, or `NA` in both bounds of
  # `Ocens()`) means one of two things. Between recorded visits (a later row
  # of the patient has a state or a range) it is one of the states that are
  # not absorbing: the patient is seen alive later, or the absorbing state is
  # recorded later on its own visit, so they were not absorbed during the
  # gap. After the patient's last recorded visit it says nothing, every state
  # is possible: the patient may have been lost to follow-up with the rows
  # left in the data, and taking them as "alive" would bias the fit towards
  # survival. Such rows are trimmed below as steps that allow every state, so
  # the fit is the same as without them. A user who knows that the patient
  # was alive at those visits gives the range with `Ocens()`.
  unrecorded <- is.na(lo)
  if (any(unrecorded)) {
    s_recorded <- !unrecorded[ord]
    patient_start <- !same_patient
    recorded_so_far <- cumsum(s_recorded)
    patient_total <- recorded_so_far[c(which(patient_start)[-1L] - 1L, n)]
    recorded_later <- patient_total[cumsum(patient_start)] - recorded_so_far > 0L
    between <- logical(n)
    between[ord] <- recorded_later
    between <- between & unrecorded
    lo[unrecorded] <- 1L
    hi[unrecorded] <- n_states
    if (any(between)) {
      possible <- setdiff(seq_len(n_states), absorb)
      if (!identical(possible, seq.int(possible[[1L]], possible[[length(possible)]]))) {
        mostr_abort(
          "A visit with a missing state between recorded visits can be in any ",
          "state that is not absorbing, and here those states are not next to ",
          "each other in the order of the states (",
          paste(labels[possible], collapse = ", "),
          "). Give the possible states of each such visit as a range with ",
          "`Ocens(lower, upper)`."
        )
      }
      lo[between] <- possible[[1L]]
      hi[between] <- possible[[length(possible)]]
    }
  }
  # An unrecorded visit is never "observed exactly", even when only one state
  # is left for it (two states, one absorbing): its row has no response for
  # the engine's fit and the next row's previous state is missing, so it
  # must be a step of a stretch.
  exact <- !unrecorded & lo == hi

  s_exact <- exact[ord]
  s_prev_known <- previous_known[ord]
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
  is_first <- !duplicated(step_run)

  if (any(!is.na(data[[p_var]][step_data_row[is_first]]) &
    is.na(match(as.character(data[[p_var]][step_data_row[is_first]]), labels)))) {
    mostr_abort(
      "`",
      p_var,
      "` holds values that are not states of the response."
    )
  }

  step_weight <- if (is.null(weights)) rep(1, n_steps) else weights[step_data_row]
  # Compared with the weight of the run's first step: one vectorised
  # comparison instead of a function call per run, which took 0.13 s for
  # 7,500 runs.
  if (n_steps && any(step_weight != step_weight[is_first][step_run])) {
    mostr_abort(
      "Row weights must be the same on all rows of a run of visits with ",
      "unobserved states, because the model has one likelihood term for ",
      "the whole run. Use one weight per patient."
    )
  }

  list(
    data = data,
    labels = labels,
    absorb = as.integer(absorb),
    id_var = id_var,
    p_var = p_var,
    id = id,
    exact = exact,
    step_lo = lo[step_data_row],
    step_hi = hi[step_data_row],
    step_weight = step_weight,
    step_run = step_run,
    is_first = is_first,
    step_data_rows = step_data_row,
    n_steps_total = n_steps
  )
}

# The design rows of the stretches laid out by `markov_path_layout()`, and
# the plan for `cpp_path_loglik()`.
#
# The design rows come in one of two layouts (see
# `inst/include/mostr/path_likelihood.hpp`). When the previous state enters
# the design only through columns that depend on nothing else
# (`markov_additive_design()`), as in `y ~ yprev + time + age`, `x_gap` holds
# one row per step with those columns set to 0 and `x_prev` the K rows of
# those columns, so that the design row for previous state j is
# `x_gap[step_rows[u, j], ] + x_prev[j, ]`: with continuous covariates,
# where no two rows are equal, this stores a fifth as many rows and lets the
# likelihood code compute one linear predictor per step. Otherwise (`yprev`
# interacting with another predictor, or `additive = FALSE`) `x_gap` holds
# every distinct design row and `x_prev` is NULL. Both layouts give the same
# likelihood up to rounding.
#
# `cluster`, for a model with random intercepts, holds the cluster number of
# each row of `layout$data`. A stretch's likelihood term then depends on its
# cluster's intercept, so two stretches are stored once only when they also
# belong to the same cluster, and the plan gives each stored stretch its
# cluster (`stretch_cluster`). Design rows (`x_gap`) are still shared between
# clusters: the intercept is added per stretch, not per design row.
markov_path_plan_design <- function(layout, design, additive = TRUE, cluster = NULL) {
  data <- layout$data
  labels <- layout$labels
  absorb <- layout$absorb
  n_states <- length(labels)
  p_var <- layout$p_var
  id <- layout$id
  step_lo <- layout$step_lo
  step_hi <- layout$step_hi
  step_run <- layout$step_run
  is_first <- layout$is_first
  step_data_row <- layout$step_data_rows
  n_steps <- layout$n_steps_total

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
  # Rows are taken column by column: `[.data.frame` makes the row names of
  # repeated rows unique, which took 0.17 s for 150,000 rows.
  take <- c(step_data_row[first_steps], step_data_row[candidate_step])
  expanded <- lapply(data, function(column) {
    if (is.null(dim(column))) column[take] else column[take, , drop = FALSE]
  })
  attr(expanded, "row.names") <- .set_row_names(length(take))
  class(expanded) <- "data.frame"
  if (length(candidate_step)) {
    rows <- length(first_steps) + seq_along(candidate_step)
    expanded[[p_var]][rows] <- markov_previous_state_values(
      labels[candidate_state],
      data[[p_var]]
    )
  }
  x_all <- if (length(take)) design(expanded) else NULL
  if (!is.null(x_all) && anyNA(x_all)) {
    bad <- unique(take[rowSums(is.na(x_all)) > 0])
    mostr_abort(
      "A visit with an unobserved state needs every predictor of the model, ",
      "because the model predicts the transitions into and out of it. ",
      length(bad),
      " such row",
      if (length(bad) > 1L) "s have" else " has",
      " a missing predictor (first: `",
      layout$id_var,
      "` = ",
      id[bad[[1L]]],
      "). Fill in the predictors, or remove the patient's rows from that ",
      "visit on."
    )
  }

  if (!n_steps) {
    plan <- markov_empty_path_plan(labels, absorb)
    if (!is.null(cluster)) {
      plan$stretch_cluster <- integer()
    }
    return(plan)
  }

  # The previous state of every row of `x_all`: the known one for a first
  # step, the candidate for a later step.
  first_prev <- match(
    as.character(data[[p_var]][step_data_row[first_steps]]),
    labels
  )
  row_step <- c(first_steps, candidate_step)
  row_prev <- c(first_prev, candidate_state)
  split_design <- if (additive) {
    markov_additive_design(x_all, row_step, row_prev, n_states)
  }

  first_row <- integer(length(first_steps))
  step_rows <- matrix(0L, n_steps, n_states)
  if (is.null(split_design)) {
    # General layout: keep each distinct design row once.
    x_prev <- NULL
    design_row <- markov_distinct_rows(x_all)
    x_gap <- x_all[!duplicated(design_row), , drop = FALSE]
    first_row <- design_row[seq_along(first_steps)]
    if (length(candidate_step)) {
      step_rows[cbind(candidate_step, candidate_state)] <-
        design_row[length(first_steps) + seq_along(candidate_step)]
    }
    start_key <- first_row
  } else {
    # Additive layout: one row per step without the previous state's
    # columns, kept once if equal, and those columns by previous state.
    x_prev <- split_design$prev
    base_row <- markov_distinct_rows(split_design$base)
    x_gap <- split_design$base[!duplicated(base_row), , drop = FALSE]
    step_base <- integer(n_steps)
    step_base[split_design$base_step] <- base_row
    first_row <- step_base[first_steps]
    if (length(candidate_step)) {
      step_rows[cbind(candidate_step, candidate_state)] <-
        step_base[candidate_step]
    }
    # The first step's row no longer says which state came before it. The
    # key is an integer: `as.character()` of a double took 26 ms for 7,440
    # stretches, against under 1 ms for an integer.
    start_key <- (first_row - 1L) * n_states + first_prev
  }
  rownames(x_gap) <- NULL
  first_cluster <- NULL
  if (!is.null(cluster)) {
    # Stretches of different clusters have different likelihood terms, so
    # the cluster becomes part of the key that decides which stretches are
    # equal. Without `cluster` the key, and with it the plan, is unchanged.
    first_cluster <- as.integer(cluster[step_data_row[first_steps]])
    start_key <- data.table::frankv(
      list(start_key, first_cluster),
      ties.method = "dense"
    )
  }

  # Keep each distinct stretch once: same design rows, same ranges. Two
  # stretches can only be equal if they share the design row of their first
  # step, so the text keys that compare whole stretches are built for those
  # alone. With a continuous covariate hardly any stretch is shared, and
  # building a key for every step would cost more than the rest of the plan.
  instance_key <- as.character(start_key)
  shared <- start_key %in% start_key[duplicated(start_key)]
  if (any(shared)) {
    in_shared <- shared[step_run]
    step_key <- do.call(
      paste,
      c(
        list(step_lo[in_shared], step_hi[in_shared]),
        asplit(step_rows[in_shared, , drop = FALSE], 2L),
        sep = ","
      )
    )
    instance_key[shared] <- paste(
      start_key[shared],
      vapply(
        split(step_key, step_run[in_shared]),
        paste,
        character(1),
        collapse = "|"
      )
    )
  }
  stretch <- match(instance_key, unique(instance_key))
  representative <- which(!duplicated(stretch))
  instance_steps <- split(seq_len(n_steps), step_run)
  kept_steps <- unlist(instance_steps[representative], use.names = FALSE)
  n_per <- lengths(instance_steps[representative])
  instance_weight <- layout$step_weight[first_steps]

  plan <- list(
    labels = labels,
    absorb = as.integer(absorb),
    layout = if (is.null(x_prev)) "general" else "additive",
    x_gap = x_gap,
    x_prev = x_prev,
    stretch_first_row = as.integer(first_row[representative]),
    stretch_first_prev = as.integer(first_prev[representative]),
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
  if (!is.null(first_cluster)) {
    plan$stretch_cluster <- first_cluster[representative]
  }
  plan
}

# Whether the design rows `x` split exactly into a part that depends on the
# step only and a part that depends on the previous state only. `row_step`
# and `row_prev` give the step and the previous state of each row.
#
# The previous state's columns are those that differ between the rows of one
# step. The split holds when each of these columns has the same value in
# every row with the same previous state, at every step; then
#   x[i, ] = base[row_step[i], ] + prev[row_prev[i], ]
# exactly, because in every column one of the two terms is 0. An
# interaction of the previous state with time or a covariate puts a column
# into the first set that fails the second test, and NULL is returned.
# Comparisons are exact: a column computed from the previous state alone
# (dummy variables, a spline of a numeric previous state) has identical
# values for identical states.
#
# Returns NULL, or `base` (one row per step that has rows, previous-state
# columns 0), `base_step` (the step of each row of `base`) and `prev` (K
# rows, other columns 0; rows of states without rows are 0).
markov_additive_design <- function(x, row_step, row_prev, n_states) {
  first_of_step <- match(row_step, row_step)
  varies <- colSums(x != x[first_of_step, , drop = FALSE]) > 0
  # Nothing depends on the previous state: all rows of a step are equal and
  # the general layout stores them once.
  if (!any(varies)) {
    return(NULL)
  }
  first_of_prev <- match(row_prev, row_prev)
  if (any(x[, varies, drop = FALSE] != x[first_of_prev, varies, drop = FALSE])) {
    return(NULL)
  }
  base_rows <- which(!duplicated(row_step))
  base <- x[base_rows, , drop = FALSE]
  base[, varies] <- 0
  prev <- matrix(0, n_states, ncol(x))
  prev_rows <- which(!duplicated(row_prev))
  prev[row_prev[prev_rows], varies] <- x[prev_rows, varies, drop = FALSE]
  list(base = unname(base), base_step = row_step[base_rows], prev = prev)
}

markov_empty_path_plan <- function(labels, absorb) {
  list(
    labels = labels,
    absorb = as.integer(absorb),
    layout = "general",
    x_gap = NULL,
    x_prev = NULL,
    stretch_first_row = integer(),
    stretch_first_prev = integer(),
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

# For each row of the numeric matrix `x`, the number of the distinct row it
# equals, numbered in order of first appearance.
#
# Equal rows are found from one number per row, a weighted sum of its
# columns with fixed irregular weights: equal rows give exactly the same
# sum, and matching one number per row is several times faster than ordering
# the rows on all columns (0.3 s of a 2 s fit with 158,000 rows and 12
# columns, Nov 2026). Two different rows could give the same sum, so the
# result is checked against `x` and the rows are ordered on all columns if
# it is ever wrong. The sum is formed column by column in R, not with
# `x %*% weights`: an optimised BLAS may add the columns of different rows in
# different orders (SIMD blocks and remainders), and equal rows would then
# get keys that differ in the last bit and silently stop being merged.
markov_distinct_rows <- function(x) {
  weights <- sqrt(seq_len(ncol(x)) + 0.5)
  key <- numeric(nrow(x))
  for (j in seq_len(ncol(x))) {
    key <- key + x[, j] * weights[[j]]
  }
  group <- match(key, unique(key))
  first <- which(!duplicated(group))
  if (!identical(unname(x[first[group], , drop = FALSE]), unname(x))) {
    rank <- data.table::frankv(
      data.table::as.data.table(x),
      ties.method = "dense"
    )
    group <- match(rank, unique(rank))
  }
  group
}
