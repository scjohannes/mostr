# Data set-up for maximum likelihood with missing and range-coded states.
#
# The first-order transition model is treated as a hidden Markov model whose
# emissions are either exact (state observed) or uninformative within a set of
# allowed states (state missing, or known only up to a range with
# rms::Ocens()). Each row carries a 0/1 mask over the K states. Rows are
# grouped into chains: the consecutive rows of one patient. A chain starts at
# a patient's first row and after any row that is not used (a missing
# covariate, or a missing response when `missing_states` is NULL); its first
# row's previous state must be observed. Within a chain, the previous state of
# a row is the state of the row before it, so only the first row's `p_var` is
# used.

# ---------------------------------------------------------------------------
# Response parsing and state scale
# ---------------------------------------------------------------------------

# TRUE when the response of `formula` is an Ocens() call.
markov_formula_has_ocens <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    return(FALSE)
  }
  lhs <- formula[[2L]]
  if (!is.call(lhs)) {
    return(FALSE)
  }
  fn <- lhs[[1L]]
  identical(fn, as.name("Ocens")) ||
    (is.call(fn) &&
      length(fn) == 3L &&
      (identical(fn[[1L]], as.name("::")) ||
        identical(fn[[1L]], as.name(":::"))) &&
      identical(fn[[3L]], as.name("Ocens")))
}

# Evaluate the response of `formula` on `data`. Returns the lower and upper
# endpoints on a common scale plus the information needed to build the state
# scale: `kind` is "numeric" (numeric values, possibly infinite endpoints) or
# "labels" (ordered category labels, compared by position in `levels`).
markov_ml_parse_response <- function(formula, data) {
  lhs <- formula[[2L]]
  env <- new.env(parent = environment(formula) %||% baseenv())
  env$Ocens <- rms::Ocens
  response <- eval(lhs, data, env)
  n <- nrow(data)

  if (inherits(response, "Ocens") && is.matrix(response)) {
    lower <- unclass(response)[, 1L]
    upper <- unclass(response)[, 2L]
    levels <- attr(response, "levels")
    if (length(lower) != n) {
      stop("The `Ocens()` response must have one row per row of `data`.")
    }
    if (!is.null(levels)) {
      return(list(
        kind = "labels",
        lower = as.numeric(lower),
        upper = as.numeric(upper),
        levels = as.character(levels),
        response_type = "factor",
        interval = TRUE
      ))
    }
    return(list(
      kind = "numeric",
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      levels = NULL,
      response_type = "numeric",
      interval = TRUE
    ))
  }

  if (is.matrix(response) || is.data.frame(response)) {
    stop(
      "The response of `formula` must be a single state column or ",
      "`Ocens(lower, upper)`.",
      call. = FALSE
    )
  }
  if (length(response) != n) {
    stop("The response must have one value per row of `data`.")
  }

  if (is.factor(response)) {
    levels <- levels(response)
    position <- as.numeric(match(as.character(response), levels))
    return(list(
      kind = "labels",
      lower = position,
      upper = position,
      levels = levels,
      response_type = "factor",
      interval = FALSE
    ))
  }
  if (is.character(response)) {
    levels <- sort(unique(response[!is.na(response)]))
    position <- as.numeric(match(response, levels))
    return(list(
      kind = "labels",
      lower = position,
      upper = position,
      levels = levels,
      response_type = "factor",
      interval = FALSE
    ))
  }
  if (is.numeric(response) || is.logical(response)) {
    value <- as.numeric(response)
    if (any(is.infinite(value))) {
      stop("A plain state response cannot contain infinite values.")
    }
    return(list(
      kind = "numeric",
      lower = value,
      upper = value,
      levels = NULL,
      response_type = "numeric",
      interval = FALSE
    ))
  }

  stop(
    "The response must be numeric, factor, character, or `Ocens()`.",
    call. = FALSE
  )
}

# Map parsed endpoints to lower and upper state indices on the 1..K scale.
# The scale consists of the exactly observed states. Endpoint rules follow
# rms: two finite endpoints give a closed interval, a lower endpoint of -Inf
# means states below the upper endpoint, and an upper endpoint of Inf means
# states above the lower endpoint. A range keeps the states of the scale it
# contains, so finite endpoints beyond the lowest or highest state are
# clipped to it; a range without any state of the scale is an error.
markov_ml_state_scale <- function(parsed, row_ids) {
  lower <- parsed$lower
  upper <- parsed$upper
  missing <- is.na(lower) & is.na(upper)
  if (any(xor(is.na(lower), is.na(upper)))) {
    stop(
      "`Ocens()` responses must be missing at both endpoints together.",
      call. = FALSE
    )
  }
  if (any(is.infinite(lower) & is.infinite(upper))) {
    bad <- which(is.infinite(lower) & is.infinite(upper))
    stop(
      "`Ocens()` responses may not have two infinite endpoints (patient(s) ",
      markov_format_patient_ids(row_ids[bad]),
      ").",
      call. = FALSE
    )
  }
  exact <- !missing & lower == upper & is.finite(lower)

  if (identical(parsed$kind, "labels")) {
    used <- sort(unique(lower[exact]))
    values <- used
    labels <- parsed$levels[used]
  } else {
    values <- sort(unique(lower[exact]))
    labels <- as.character(values)
  }
  if (length(values) < 2L) {
    stop(
      "At least two different states must be observed exactly.",
      call. = FALSE
    )
  }
  K <- length(values)

  lo <- rep(NA_integer_, length(lower))
  hi <- rep(NA_integer_, length(lower))
  present <- which(!missing)
  a <- lower[present]
  b <- upper[present]
  # First and last state of the scale inside the interval. findInterval()
  # counts the states below (`left.open = TRUE`) or at most at a value, so an
  # endpoint below the lowest state gives `first = 1` and one above the
  # highest state gives `last = K`.
  first <- ifelse(
    is.infinite(a),
    1L,
    ifelse(
      is.infinite(b),
      findInterval(a, values) + 1L,
      findInterval(a, values, left.open = TRUE) + 1L
    )
  )
  last <- ifelse(
    is.infinite(a),
    findInterval(b, values, left.open = TRUE),
    ifelse(is.infinite(b), K, findInterval(b, values))
  )
  empty <- first > last | b < a
  lo[present] <- ifelse(empty, NA_integer_, as.integer(first))
  hi[present] <- ifelse(empty, -1L, as.integer(last))
  invalid <- which(!is.na(hi) & hi < 0L)
  if (length(invalid)) {
    stop(
      "`Ocens()` range(s) contain none of the states observed exactly (",
      paste(labels, collapse = ", "),
      ") for patient(s) ",
      markov_format_patient_ids(row_ids[invalid]),
      ". Each range must contain at least one state that is observed exactly ",
      "somewhere in the data.",
      call. = FALSE
    )
  }

  list(
    labels = labels,
    values = values,
    K = K,
    lo = lo,
    hi = hi,
    missing = missing,
    exact = !missing & lo == hi
  )
}

# Indices on the state scale of the states a missing response allows.
markov_ml_missing_index <- function(missing_states, labels) {
  index <- match(as_state_labels(missing_states), labels)
  if (anyNA(index)) {
    stop(
      "`missing_states` value(s) ",
      paste(missing_states[is.na(index)], collapse = ", "),
      " are not among the states observed exactly: ",
      paste(labels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  sort(unique(index))
}

# Compatibility masks: one 0/1 row of length K per element of `lo`. Missing
# responses (NA `lo`) allow the states `missing_index`.
markov_ml_masks <- function(lo, hi, K, missing_index) {
  n <- length(lo)
  fill <- numeric(K)
  fill[missing_index] <- 1
  C <- matrix(rep(fill, each = n), nrow = n, ncol = K)
  present <- which(!is.na(lo))
  if (length(present)) {
    C[present, ] <- 0
    width <- hi[present] - lo[present] + 1L
    C[cbind(
      rep(present, width),
      sequence(width, from = lo[present])
    )] <- 1
  }
  C
}

# Internal response and weight columns of the start and final fits.
markov_ml_response_column <- function() {
  ".mostr_state"
}

markov_ml_weight_column <- function() {
  ".mostr_em_weight"
}

# Formula with the response replaced by the internal exact-state column.
markov_ml_working_formula <- function(formula) {
  out <- formula
  out[[2L]] <- as.name(markov_ml_response_column())
  out
}

# Exact-state values for the internal response: numeric for numeric responses
# (rms::orm() then keeps numeric state values), otherwise an ordered factor
# with the scale's labels.
markov_ml_state_values <- function(index, scale, response_type) {
  if (identical(response_type, "numeric")) {
    return(scale$values[index])
  }
  factor(scale$labels[index], levels = scale$labels, ordered = TRUE)
}

# ---------------------------------------------------------------------------
# Set-up
# ---------------------------------------------------------------------------

# Rows, chains, masks, and the split into fully observed transitions and runs
# of rows that need the forward-backward recursion. `missing_states` is NULL
# (rows with a missing response are not used) or the states a missing
# response may be. `extra_vars` are covariates of the model outside `formula`
# (the `ppo` terms of the blrm engine). The design matrices of the possible
# previous states are added by the engine (markov_ml_add_design()).
markov_ml_prepare <- function(
  formula,
  data,
  id_var,
  time_var,
  p_var,
  missing_states,
  extra_vars = character()
) {
  if (!nrow(data)) {
    stop("`data` has no rows.", call. = FALSE)
  }
  if (anyNA(data[[id_var]])) {
    stop("Patient IDs in `", id_var, "` must not be missing.", call. = FALSE)
  }
  if (!time_var %in% names(data)) {
    stop("Time column `", time_var, "` not found in `data`.", call. = FALSE)
  }
  # Time orders each patient's rows into a chain. Character time would sort
  # as text ("day 10" before "day 2"); a factor sorts by its levels.
  time <- data[[time_var]]
  if (!is.numeric(time) && !is.factor(time)) {
    stop(
      "With `missing_states` or an `Ocens()` response, `time_var` must name ",
      "a numeric column or a factor whose levels are in time order, because ",
      "it orders each patient's visits; `",
      time_var,
      "` is of class ",
      class(time)[[1L]],
      ".",
      call. = FALSE
    )
  }
  if (anyNA(time)) {
    stop("The time column `", time_var, "` must not be missing.", call. = FALSE)
  }

  ord <- order(data[[id_var]], time)
  sd <- data[ord, , drop = FALSE]
  n_all <- nrow(sd)
  ids_all <- as.character(sd[[id_var]])

  rhs_vars <- all.vars(stats::delete.response(stats::terms(formula)))
  covariate_vars <- setdiff(
    intersect(unique(c(rhs_vars, extra_vars)), names(sd)),
    p_var
  )
  missing_covariate <- if (length(covariate_vars)) {
    !stats::complete.cases(sd[covariate_vars])
  } else {
    logical(n_all)
  }
  parsed <- markov_ml_parse_response(formula, sd)
  missing_response <- is.na(parsed$lower) & is.na(parsed$upper)

  # Rows not used, and the chain starts they create. A chain's first row
  # needs an observed previous state; without one it is not used either, and
  # the next row starts the chain.
  drop_covariate <- missing_covariate
  drop_response <- !missing_covariate &
    missing_response &
    is.null(missing_states)
  drop <- drop_covariate | drop_response
  same_patient <- c(FALSE, ids_all[-1L] == ids_all[-n_all])
  prev_missing <- is.na(sd[[p_var]])
  repeat {
    kept <- !drop
    starts <- kept & !(same_patient & c(FALSE, kept[-n_all]))
    bad_start <- starts & prev_missing
    if (!any(bad_start)) {
      break
    }
    drop <- drop | bad_start
  }
  drop_previous <- drop & !drop_covariate & !drop_response
  if (!any(kept)) {
    stop(
      "No rows can be used: every row has a missing covariate, response, or ",
      "previous state.",
      call. = FALSE
    )
  }

  sd <- sd[kept, , drop = FALSE]
  source_rows <- ord[kept]
  n <- nrow(sd)
  ids <- ids_all[kept]
  starts <- starts[kept]
  parsed$lower <- parsed$lower[kept]
  parsed$upper <- parsed$upper[kept]
  patient_ids <- unique(ids)
  row_patient <- match(ids, patient_ids)
  row_chain <- cumsum(starts)

  scale <- markov_ml_state_scale(parsed, ids)
  K <- scale$K
  labels <- scale$labels
  missing_index <- if (is.null(missing_states)) {
    integer()
  } else {
    markov_ml_missing_index(missing_states, labels)
  }

  # Rows whose observation rules out at least one state.
  n_allowed <- scale$hi - scale$lo + 1L
  n_allowed[scale$missing] <- length(missing_index)
  informative <- n_allowed < K

  # Previous states: observed at chain starts, from the chain afterwards.
  prototype <- sd[[p_var]]
  prev_label <- ifelse(
    is.na(prototype),
    NA_character_,
    as_state_labels(prototype)
  )
  prev_index <- match(prev_label, labels)
  start_rows <- which(starts)
  baseline <- prev_index[start_rows]
  if (anyNA(baseline)) {
    bad <- start_rows[is.na(baseline)]
    stop(
      "The previous state (`",
      p_var,
      "`) at the start of a patient's rows is not one of the states observed ",
      "exactly (",
      paste(labels, collapse = ", "),
      ") for patient(s) ",
      markov_format_patient_ids(ids[bad]),
      ".",
      call. = FALSE
    )
  }

  exact_row <- scale$exact
  exact_state <- ifelse(exact_row, scale$lo, NA_integer_)
  prev_exact <- logical(n)
  prev_state <- rep(NA_integer_, n)
  prev_exact[start_rows] <- TRUE
  prev_state[start_rows] <- baseline
  later <- which(!starts)
  prev_exact[later] <- exact_row[later - 1L]
  prev_state[later] <- exact_state[later - 1L]

  check <- later[prev_exact[later] & !is.na(prev_label[later])]
  inconsistent <- check[
    is.na(prev_index[check]) | prev_index[check] != prev_state[check]
  ]
  if (length(inconsistent)) {
    stop(
      "`",
      p_var,
      "` disagrees with the state observed on the previous row for ",
      "patient(s) ",
      markov_format_patient_ids(ids[inconsistent]),
      ".",
      call. = FALSE
    )
  }

  complete <- exact_row & prev_exact
  if (!any(complete)) {
    stop("No transition has both states observed exactly.", call. = FALSE)
  }
  # The model describes transitions from the states that occur as the
  # previous state of a fully observed transition.
  origins <- sort(unique(prev_state[complete]))
  unseen_outcomes <- setdiff(seq_len(K), exact_state[complete])
  if (length(unseen_outcomes)) {
    stop(
      "State(s) ",
      paste(labels[unseen_outcomes], collapse = ", "),
      " never occur as the outcome of a transition with both states ",
      "observed exactly, so their thresholds cannot be estimated.",
      call. = FALSE
    )
  }

  # The chain splits at every exactly observed state (Markov property). A
  # row whose own and previous states are both observed is an ordinary
  # transition; only segments that contain missing or range rows, each
  # starting from the last observed state, need the recursion.
  chain_rows <- which(!complete)
  segment <- cumsum(prev_exact)[chain_rows]
  segment <- match(segment, unique(segment))
  first_of_segment <- !duplicated(segment)
  segment_start <- chain_rows[first_of_segment]
  segment_baseline <- prev_state[segment_start]
  chain_pos <- stats::ave(seq_along(chain_rows), segment, FUN = seq_along)
  segment_length <- tabulate(segment, length(segment_start))
  C <- markov_ml_masks(
    scale$lo[chain_rows],
    scale$hi[chain_rows],
    K,
    missing_index
  )

  # The recursion leaves a hidden state through the transition model, so a
  # state that a segment may pass through before its last row must be one
  # the model can leave.
  not_origin <- setdiff(seq_len(K), origins)
  bad_baseline <- !segment_baseline %in% origins
  if (any(bad_baseline)) {
    stop(
      "State(s) ",
      paste(unique(labels[segment_baseline[bad_baseline]]), collapse = ", "),
      " are followed by a missing or range-coded state for patient(s) ",
      markov_format_patient_ids(ids[segment_start[bad_baseline]]),
      ", but never occur as the previous state of a transition with both ",
      "states observed, so the model cannot describe leaving them.",
      call. = FALSE
    )
  }
  inner <- chain_pos < segment_length[segment]
  if (length(not_origin) && any(inner)) {
    hidden_bad <- inner & rowSums(C[, not_origin, drop = FALSE]) > 0
    if (any(hidden_bad)) {
      bad_states <- not_origin[colSums(C[hidden_bad, not_origin, drop = FALSE]) > 0]
      stop(
        "A missing or range-coded state that is followed by further rows of ",
        "the same patient may be state(s) ",
        paste(labels[bad_states], collapse = ", "),
        " (patient(s) ",
        markov_format_patient_ids(ids[chain_rows[hidden_bad]]),
        "). These states never occur as the previous state of a transition ",
        "with both states observed, so the model cannot describe leaving ",
        "them; usually they are absorbing states such as death (`absorb` in ",
        "`sops()`). A patient who may have entered such a state is lost to ",
        "follow-up from then on: end the patient's rows at the last visit ",
        "before it, or exclude these states from `missing_states` or the ",
        "`Ocens()` range.",
        call. = FALSE
      )
    }
  }

  # Design rows for every possible previous state are keyed by the unique
  # covariate rows and the origin. Only keys that the likelihood uses are
  # kept: each transition's own previous state, and every origin for rows
  # that need the recursion.
  if (length(covariate_vars)) {
    key <- data.table::as.data.table(sd[covariate_vars])
    cov_id <- data.table::frankv(key, ties.method = "dense")
  } else {
    cov_id <- rep(1L, n)
  }
  n_origin <- length(origins)
  comp_rows <- which(complete)
  comp_key <- (cov_id[comp_rows] - 1L) * n_origin +
    match(prev_state[comp_rows], origins)
  chain_key <- outer(
    (cov_id[chain_rows] - 1L) * n_origin,
    seq_len(n_origin),
    "+"
  )
  used <- sort(unique(c(comp_key, as.vector(chain_key))))
  used_cov <- (used - 1L) %/% n_origin + 1L
  used_origin <- (used - 1L) %% n_origin + 1L
  design_data <- sd[match(used_cov, cov_id), , drop = FALSE]
  design_data[[p_var]] <- coerce_previous_state_values(
    labels[origins[used_origin]],
    prototype,
    p_var
  )
  rownames(design_data) <- NULL
  DI <- matrix(match(chain_key, used), nrow = length(chain_rows))
  storage.mode(DI) <- "integer"
  comp_design <- match(comp_key, used)
  comp_state <- exact_state[comp_rows]

  # Fully observed transitions, aggregated by (design row, state).
  counts <- rowsum(
    rep(1, length(comp_rows)),
    (comp_design - 1L) * K + comp_state
  )
  keys <- as.integer(rownames(counts)) - 1L
  complete_part <- list(
    rows = comp_rows,
    origin = match(prev_state[comp_rows], origins),
    design = comp_design,
    state = comp_state,
    patient = row_patient[comp_rows],
    agg_design = keys %/% K + 1L,
    agg_state = keys %% K + 1L,
    agg_count = counts[, 1L]
  )

  rows_by_pos <- split(seq_along(chain_rows), chain_pos)
  chain <- list(
    rows = chain_rows,
    n = length(chain_rows),
    N = length(segment_start),
    baseline = segment_baseline,
    segment = segment,
    patient = row_patient[chain_rows],
    ids = ids[chain_rows],
    C = C,
    DI = DI,
    rows_by_pos = rows_by_pos,
    pats_by_pos = lapply(rows_by_pos, function(r) segment[r])
  )

  start_data <- sd[comp_rows, , drop = FALSE]
  start_data[[p_var]] <- coerce_previous_state_values(
    labels[prev_state[comp_rows]],
    prototype,
    p_var
  )

  list(
    data = sd,
    source_rows = source_rows,
    ids = ids,
    patient_ids = patient_ids,
    row_patient = row_patient,
    row_chain = row_chain,
    complete_part = complete_part,
    chain = chain,
    N = length(patient_ids),
    n = n,
    K = K,
    labels = labels,
    scale = scale,
    response_type = parsed$response_type,
    interval_response = parsed$interval,
    missing_index = missing_index,
    origins = origins,
    n_origin = n_origin,
    complete = complete,
    informative = informative,
    p_var = p_var,
    prototype = prototype,
    design_data = design_data,
    start_data = start_data,
    working_formula = markov_ml_working_formula(formula),
    dropped = list(
      covariate = sum(drop_covariate),
      response = sum(drop_response),
      previous = sum(drop_previous)
    )
  )
}

# Add the internal response column of exact states to rows of the problem.
markov_ml_with_state <- function(problem, data, state) {
  data[[markov_ml_response_column()]] <- markov_ml_state_values(
    state,
    problem$scale,
    problem$response_type
  )
  data
}

markov_ml_missing_counts <- function(problem) {
  exact <- problem$scale$exact
  missing <- problem$scale$missing
  n_range <- sum(!exact & !missing)
  # Runs of rows without an exact state; a run is a gap when an exact row
  # follows it within the same chain and trailing otherwise.
  next_exact <- c(exact[-1L], FALSE)
  next_same <- c(problem$row_chain[-1L], NA) == problem$row_chain
  next_same[is.na(next_same)] <- FALSE
  run_end <- !exact & (!next_same | next_exact)
  closes_gap <- run_end & next_same & next_exact
  list(
    n_missing = sum(missing),
    n_range = n_range,
    n_gaps = sum(closes_gap),
    n_trailing = sum(run_end & !next_same),
    n_patients = problem$N,
    n_rows = problem$n
  )
}
