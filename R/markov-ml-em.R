# Observed-data likelihood and EM fitting for `orm_markov_ml()`.
#
# The first-order proportional-odds transition model is treated as a hidden
# Markov model whose emissions are either exact (state observed) or
# uninformative within a set of allowed states (state missing or known only up
# to an interval). Each row of the long data carries a 0/1 compatibility mask
# over the K states. A scaled forward pass gives the log-likelihood; the
# matching backward pass gives the posterior transition weights used by the
# EM M-step and by the exact (Fisher-identity) score.

# ---------------------------------------------------------------------------
# Response parsing and state scale
# ---------------------------------------------------------------------------

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
# states above the lower endpoint.
markov_ml_state_scale <- function(parsed, row_ids) {
  lower <- parsed$lower
  upper <- parsed$upper
  missing <- is.na(lower) & is.na(upper)
  if (any(xor(is.na(lower), is.na(upper)))) {
    stop("Interval responses must be missing at both endpoints together.")
  }
  if (any(is.infinite(lower) & is.infinite(upper))) {
    bad <- which(is.infinite(lower) & is.infinite(upper))
    stop(
      "Interval responses may not have two infinite endpoints (rows ",
      paste(utils::head(bad, 5L), collapse = ", "),
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
  for (r in present) {
    a <- lower[[r]]
    b <- upper[[r]]
    allowed <- if (is.infinite(a)) {
      values < b
    } else if (is.infinite(b)) {
      values > a
    } else {
      values >= a & values <= b
    }
    if (
      (is.finite(a) && (a < values[[1L]] || a > values[[K]])) ||
        (is.finite(b) && (b < values[[1L]] || b > values[[K]]))
    ) {
      allowed <- FALSE
    }
    if (!any(allowed) || b < a) {
      lo[[r]] <- NA_integer_
      hi[[r]] <- -1L
      next
    }
    idx <- which(allowed)
    lo[[r]] <- min(idx)
    hi[[r]] <- max(idx)
  }
  invalid <- which(!is.na(hi) & hi < 0L)
  if (length(invalid)) {
    stop(
      "Interval response(s) are empty or outside the observed states ",
      paste(labels, collapse = ", "),
      " for patient(s) ",
      markov_format_patient_ids(row_ids[invalid]),
      ". Each interval must contain at least one state on the model's scale.",
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

# Compatibility masks: one 0/1 row of length K per data row.
markov_ml_masks <- function(lo, hi, K, absorbing, missing_states) {
  n <- length(lo)
  C <- matrix(0, nrow = n, ncol = K)
  present <- !is.na(lo)
  states <- matrix(seq_len(K), nrow = n, ncol = K, byrow = TRUE)
  C[present, ] <- (states[present, , drop = FALSE] >= lo[present]) &
    (states[present, , drop = FALSE] <= hi[present])
  if (any(!present)) {
    fill <- rep(1, K)
    if (identical(missing_states, "non_absorbing")) {
      fill[absorbing] <- 0
    }
    C[!present, ] <- matrix(fill, nrow = sum(!present), ncol = K, byrow = TRUE)
  }
  storage.mode(C) <- "double"
  C
}

# ---------------------------------------------------------------------------
# Set-up
# ---------------------------------------------------------------------------

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

# Exact-state values formatted like the user's response: numeric for numeric
# responses, otherwise an ordered factor with the scale's labels.
markov_ml_state_values <- function(index, scale, response_type) {
  if (identical(response_type, "numeric")) {
    return(scale$values[index])
  }
  factor(scale$labels[index], levels = scale$labels, ordered = TRUE)
}

markov_ml_setup <- function(
  formula,
  data,
  id_var,
  time_var,
  p_var,
  absorbing_state,
  missing_states,
  dots_fit
) {
  n <- nrow(data)
  if (n < 1L) {
    stop("`data` has no rows.")
  }
  ids_raw <- data[[id_var]]
  if (anyNA(ids_raw)) {
    stop("Patient IDs in `", id_var, "` must not be missing.", call. = FALSE)
  }
  if (anyNA(data[[time_var]])) {
    stop("The time column `", time_var, "` must not be missing.", call. = FALSE)
  }

  ord <- order(ids_raw, data[[time_var]])
  sd <- data[ord, , drop = FALSE]
  ids <- as.character(sd[[id_var]])
  patient_ids <- unique(ids)
  row_patient <- match(ids, patient_ids)
  N <- length(patient_ids)
  row_pos <- stats::ave(seq_len(n), row_patient, FUN = seq_along)
  first_row <- row_pos == 1L

  # Covariates other than the previous state must be observed.
  rhs_vars <- all.vars(stats::delete.response(stats::terms(formula)))
  covariate_vars <- setdiff(intersect(rhs_vars, names(sd)), p_var)
  for (v in covariate_vars) {
    if (anyNA(sd[[v]])) {
      stop(
        "Covariate `",
        v,
        "` is missing for patient(s) ",
        markov_format_patient_ids(ids[is.na(sd[[v]])]),
        ". Only the outcome state may be missing.",
        call. = FALSE
      )
    }
  }

  parsed <- markov_ml_parse_response(formula, sd)
  scale <- markov_ml_state_scale(parsed, ids)
  K <- scale$K
  labels <- scale$labels

  absorbing <- match(as_state_labels(absorbing_state), labels)
  if (anyNA(absorbing)) {
    stop(
      "`absorbing_state` value(s) ",
      paste(absorbing_state[is.na(absorbing)], collapse = ", "),
      " are not among the observed states ",
      paste(labels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  absorbing <- sort(unique(absorbing))
  origins <- setdiff(seq_len(K), absorbing)
  if (!length(origins)) {
    stop("At least one state must be non-absorbing.", call. = FALSE)
  }

  C <- markov_ml_masks(scale$lo, scale$hi, K, absorbing, missing_states)

  # Previous-state column: baseline and consistency checks.
  prototype <- sd[[p_var]]
  prev_label <- ifelse(is.na(prototype), NA_character_, as_state_labels(prototype))
  prev_index <- match(prev_label, labels)
  baseline_rows <- which(first_row)
  baseline <- prev_index[baseline_rows]
  bad_baseline <- is.na(baseline)
  if (any(bad_baseline)) {
    stop(
      "The baseline state (`",
      p_var,
      "` on each patient's first row) is missing or not an observed state ",
      "for patient(s) ",
      markov_format_patient_ids(ids[baseline_rows[bad_baseline]]),
      ".",
      call. = FALSE
    )
  }
  if (any(baseline %in% absorbing)) {
    stop(
      "Patient(s) ",
      markov_format_patient_ids(ids[baseline_rows[baseline %in% absorbing]]),
      " start in an absorbing state; remove rows after absorption with ",
      "`prepare_markov_data()`.",
      call. = FALSE
    )
  }

  exact_row <- scale$exact
  exact_state <- ifelse(exact_row, scale$lo, NA_integer_)
  prev_exact <- logical(n)
  prev_state <- rep(NA_integer_, n)
  prev_exact[first_row] <- TRUE
  prev_state[first_row] <- baseline
  later <- which(!first_row)
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

  complete <- exact_row & prev_exact & !(prev_state %in% absorbing)
  if (!any(complete)) {
    stop("No transition has both states observed exactly.", call. = FALSE)
  }
  start_data <- sd[complete, , drop = FALSE]
  start_data[[p_var]] <- coerce_previous_state_values(
    labels[prev_state[complete]],
    prototype,
    p_var
  )
  start_data[[markov_ml_response_column()]] <- markov_ml_state_values(
    exact_state[complete],
    scale,
    parsed$response_type
  )

  unseen_outcomes <- setdiff(seq_len(K), exact_state[complete])
  if (length(unseen_outcomes)) {
    stop(
      "State(s) ",
      paste(labels[unseen_outcomes], collapse = ", "),
      " never occur as an outcome among transitions with both states ",
      "observed exactly, so their thresholds cannot be estimated.",
      call. = FALSE
    )
  }
  unseen_origins <- setdiff(origins, prev_state[complete])
  if (length(unseen_origins)) {
    stop(
      "Non-absorbing state(s) ",
      paste(labels[unseen_origins], collapse = ", "),
      " never occur as the previous state among transitions with both ",
      "states observed exactly.",
      call. = FALSE
    )
  }

  working_formula <- markov_ml_working_formula(formula)
  start_fit <- markov_ml_fit_orm(
    working_formula,
    start_data,
    dots_fit,
    weights = NULL,
    initial = NULL
  )
  if (
    length(start_fit$yunique) != K ||
      !identical(as_state_labels(start_fit$yunique), labels)
  ) {
    stop("Internal error: the starting fit does not use every state.")
  }
  if (isTRUE(start_fit$fail)) {
    stop("The starting fit on fully observed transitions did not converge.")
  }

  # Design rows for every possible origin, built once from the unique
  # covariate rows with the starting fit's stored design.
  if (length(covariate_vars)) {
    key <- data.table::as.data.table(sd[covariate_vars])
    cov_id <- data.table::frankv(key, ties.method = "dense")
  } else {
    cov_id <- rep(1L, n)
  }
  n_cov <- max(cov_id)
  unique_rows <- match(seq_len(n_cov), cov_id)
  n_origin <- length(origins)
  nd <- sd[rep(unique_rows, each = n_origin), , drop = FALSE]
  nd[[p_var]] <- coerce_previous_state_values(
    labels[rep(origins, times = n_cov)],
    prototype,
    p_var
  )
  rownames(nd) <- NULL
  X <- orm_model_matrix(start_fit, nd, include_intercept = FALSE)
  if (any(!is.finite(X))) {
    stop("The design matrix for possible previous states is not finite.")
  }
  DI <- (cov_id - 1L) * n_origin +
    matrix(seq_len(n_origin), nrow = n, ncol = n_origin, byrow = TRUE)
  storage.mode(DI) <- "integer"

  # The chain splits at every exactly observed state (Markov property). A row
  # whose own and previous states are both observed is an ordinary
  # proportional-odds row; only segments that contain missing or range rows,
  # each starting from the last observed state, need the recursion.
  origin_index <- match(prev_state, origins)
  comp_rows <- which(complete)
  comp_origin <- origin_index[comp_rows]
  comp_design <- DI[cbind(comp_rows, comp_origin)]
  comp_state <- exact_state[comp_rows]
  W_complete <- matrix(0, nrow = nrow(X), ncol = K)
  counts <- rowsum(
    rep(1, length(comp_rows)),
    (comp_design - 1L) * K + comp_state
  )
  keys <- as.integer(rownames(counts)) - 1L
  W_complete[cbind(keys %/% K + 1L, keys %% K + 1L)] <- counts[, 1L]
  agg <- which(W_complete > 0, arr.ind = TRUE)
  complete_part <- list(
    rows = comp_rows,
    origin = comp_origin,
    design = comp_design,
    state = comp_state,
    patient = row_patient[comp_rows],
    W = W_complete,
    agg_design = agg[, 1L],
    agg_state = agg[, 2L],
    agg_count = W_complete[agg]
  )

  chain_rows <- which(!complete)
  segment <- cumsum(prev_exact)[chain_rows]
  segment <- match(segment, unique(segment))
  segment_start <- chain_rows[!duplicated(segment)]
  chain_pos <- stats::ave(seq_along(chain_rows), segment, FUN = seq_along)
  rows_by_pos <- split(seq_along(chain_rows), chain_pos)
  chain <- list(
    rows = chain_rows,
    n = length(chain_rows),
    N = length(segment_start),
    baseline = prev_state[segment_start],
    segment = segment,
    patient = row_patient[chain_rows],
    ids = ids[chain_rows],
    C = C[chain_rows, , drop = FALSE],
    DI = DI[chain_rows, , drop = FALSE],
    rows_by_pos = rows_by_pos,
    pats_by_pos = lapply(rows_by_pos, function(r) segment[r])
  )

  informative <- rowSums(C) < K
  coefficients <- stats::coef(start_fit)

  list(
    data = sd,
    order = ord,
    ids = ids,
    patient_ids = patient_ids,
    row_patient = row_patient,
    row_pos = row_pos,
    complete_part = complete_part,
    chain = chain,
    N = N,
    n = n,
    K = K,
    labels = labels,
    scale = scale,
    response_type = parsed$response_type,
    interval_response = parsed$interval,
    absorbing = absorbing,
    origins = origins,
    n_origin = n_origin,
    baseline = baseline,
    C = C,
    complete = complete,
    informative = informative,
    p_var = p_var,
    prototype = prototype,
    cov_id = cov_id,
    design_data = nd,
    X = X,
    DI = DI,
    n_int = K - 1L,
    coef_names = names(coefficients),
    start_fit = start_fit,
    working_formula = working_formula
  )
}

# Fit rms::orm() with an optional weight column and starting values, passing
# through the user's extra arguments. Data and weights are referenced by name
# so the stored call stays small.
markov_ml_fit_orm <- function(formula, data, dots_fit, weights, initial) {
  fit_call <- as.call(c(
    list(quote(rms::orm)),
    list(formula = formula, data = quote(data), x = TRUE, y = TRUE),
    dots_fit
  ))
  if (!is.null(weights)) {
    fit_call$weights <- as.name(weights)
  }
  if (!is.null(initial)) {
    fit_call$initial <- unname(initial)
  }
  suppress_orm_bootstrap_weight_warning(
    structure(list(), class = "orm"),
    eval(fit_call)
  )
}

# ---------------------------------------------------------------------------
# Transition probabilities
# ---------------------------------------------------------------------------

# Probabilities for every unique design row. Returns NULL when `strict = FALSE`
# and the coefficients give invalid probabilities (used by step halving).
markov_ml_transition_probs <- function(problem, theta, strict = TRUE) {
  K <- problem$K
  n_int <- problem$n_int
  alpha <- theta[seq_len(n_int)]
  beta <- theta[-seq_len(n_int)]
  eta <- drop(problem$X %*% beta)
  cum <- stats::plogis(outer(eta, alpha, "+"))
  if (!is.matrix(cum)) {
    cum <- matrix(cum, ncol = n_int)
  }
  P <- cbind(1, cum) - cbind(cum, 0)
  if (any(P < -1e-10) || any(abs(rowSums(P) - 1) > 1e-12)) {
    if (!strict) {
      return(NULL)
    }
    stop(
      "Transition probabilities are invalid: thresholds are not decreasing.",
      call. = FALSE
    )
  }
  P[P < 0] <- 0
  list(P = P, dens = cum * (1 - cum), theta = theta)
}

# ---------------------------------------------------------------------------
# Forward-backward recursion
# ---------------------------------------------------------------------------

# Scaled forward pass over the segments that need the recursion. `FP` holds
# the forward vector before each segment row (the filtered distribution of the
# previous state) and `s` the scaling constants; sum(log(s)) is the segments'
# log-likelihood.
markov_ml_forward <- function(problem, probs) {
  chain <- problem$chain
  K <- problem$K
  P <- probs$P
  origins <- problem$origins
  absorbing <- problem$absorbing

  cur <- matrix(0, nrow = chain$N, ncol = K)
  cur[cbind(seq_len(chain$N), chain$baseline)] <- 1
  FP <- matrix(0, nrow = chain$n, ncol = K)
  s <- numeric(chain$n)

  for (k in seq_along(chain$rows_by_pos)) {
    r <- chain$rows_by_pos[[k]]
    i <- chain$pats_by_pos[[k]]
    fp <- cur[i, , drop = FALSE]
    FP[r, ] <- fp
    new <- matrix(0, nrow = length(r), ncol = K)
    for (ai in seq_along(origins)) {
      new <- new + fp[, origins[[ai]]] * P[chain$DI[r, ai], , drop = FALSE]
    }
    if (length(absorbing)) {
      new[, absorbing] <- new[, absorbing] + fp[, absorbing]
    }
    new <- new * chain$C[r, , drop = FALSE]
    sk <- rowSums(new)
    s[r] <- sk
    cur[i, ] <- new / sk
  }

  list(FP = FP, s = s, loglik = sum(log(s)), valid = all(s > 0))
}

markov_ml_stop_impossible <- function(ids) {
  stop(
    "Observed states have zero probability under the transition model for ",
    "patient(s) ",
    markov_format_patient_ids(ids),
    ". Check for impossible transitions, such as leaving an absorbing state.",
    call. = FALSE
  )
}

# Scaled backward pass. Returns Q = c * g / s for each segment row, so the
# posterior of a transition is xi(a, b) = FP[a] * P[a, b] * Q[b].
markov_ml_backward <- function(problem, probs, fwd) {
  chain <- problem$chain
  K <- problem$K
  P <- probs$P
  origins <- problem$origins
  absorbing <- problem$absorbing

  Q <- matrix(0, nrow = chain$n, ncol = K)
  gcur <- matrix(1, nrow = chain$N, ncol = K)
  for (k in rev(seq_along(chain$rows_by_pos))) {
    r <- chain$rows_by_pos[[k]]
    i <- chain$pats_by_pos[[k]]
    q <- chain$C[r, , drop = FALSE] * gcur[i, , drop = FALSE] / fwd$s[r]
    Q[r, ] <- q
    gp <- matrix(0, nrow = length(r), ncol = K)
    for (ai in seq_along(origins)) {
      gp[, origins[[ai]]] <- rowSums(P[chain$DI[r, ai], , drop = FALSE] * q)
    }
    if (length(absorbing)) {
      gp[, absorbing] <- q[, absorbing]
    }
    gcur[i, ] <- gp
  }
  Q
}

# Posterior transition weights of the segment rows for origin index `ai`.
markov_ml_chain_xi <- function(problem, state, ai) {
  a <- problem$origins[[ai]]
  state$fwd$FP[, a] *
    state$probs$P[problem$chain$DI[, ai], , drop = FALSE] *
    state$Q
}

# Posterior transition weights of every data row for origin index `ai`
# (n x K); fully observed transitions are point masses.
markov_ml_xi <- function(problem, state, ai) {
  out <- matrix(0, nrow = problem$n, ncol = problem$K)
  if (problem$chain$n) {
    out[problem$chain$rows, ] <- markov_ml_chain_xi(problem, state, ai)
  }
  comp <- problem$complete_part
  sel <- comp$origin == ai
  out[cbind(comp$rows[sel], comp$state[sel])] <- 1
  out
}

# Forward-backward pass at `theta`. With `strict = FALSE` an invalid
# coefficient vector or a zero-probability observation returns NULL instead of
# an error; callers use this only for step control.
markov_ml_estep <- function(problem, theta, strict = TRUE, backward = TRUE) {
  probs <- markov_ml_transition_probs(problem, theta, strict = strict)
  if (is.null(probs)) {
    return(NULL)
  }
  comp <- problem$complete_part
  p_complete <- probs$P[cbind(comp$agg_design, comp$agg_state)]
  if (any(!(p_complete > 0))) {
    if (!strict) {
      return(NULL)
    }
    p_rows <- probs$P[cbind(comp$design, comp$state)]
    markov_ml_stop_impossible(problem$ids[comp$rows[!(p_rows > 0)]])
  }
  fwd <- markov_ml_forward(problem, probs)
  if (!fwd$valid) {
    if (!strict) {
      return(NULL)
    }
    markov_ml_stop_impossible(problem$chain$ids[!(fwd$s > 0)])
  }
  loglik <- fwd$loglik + sum(comp$agg_count * log(p_complete))
  out <- list(theta = theta, probs = probs, fwd = fwd, loglik = loglik)
  if (backward) {
    out$Q <- markov_ml_backward(problem, probs, fwd)
  }
  out
}

markov_ml_loglik <- function(problem, theta, strict = TRUE) {
  state <- markov_ml_estep(problem, theta, strict = strict, backward = FALSE)
  if (is.null(state)) {
    return(NA_real_)
  }
  state$loglik
}

# Expected complete-data counts aggregated by (unique design row, destination
# state), as integer keys.
markov_ml_aggregate_weights <- function(problem, state) {
  W <- problem$complete_part$W
  if (problem$chain$n) {
    for (ai in seq_len(problem$n_origin)) {
      xi <- markov_ml_chain_xi(problem, state, ai)
      agg <- rowsum(xi, problem$chain$DI[, ai], reorder = TRUE)
      rows <- as.integer(rownames(agg))
      W[rows, ] <- W[rows, ] + agg
    }
  }
  idx <- which(W > 0, arr.ind = TRUE)
  list(design = idx[, 1L], state = idx[, 2L], weight = W[idx])
}

# Weighted rms::orm.fit() M-step on the aggregated rows. orm.fit() always takes
# at least one full Newton-Raphson step from `initial` and checks convergence
# on the gradient before that step, so its default tolerances are adequate.
markov_ml_mstep <- function(problem, triples, theta) {
  fit <- suppress_orm_bootstrap_weight_warning(
    structure(list(), class = "orm"),
    rms::orm.fit(
      x = problem$X[triples$design, , drop = FALSE],
      y = triples$state,
      weights = triples$weight,
      initial = unname(theta),
      compstats = FALSE,
      maxit = 50L
    )
  )
  if (isTRUE(fit$fail) || length(stats::coef(fit)) != length(theta)) {
    stop("The weighted orm M-step failed to converge.", call. = FALSE)
  }
  out <- unname(stats::coef(fit))
  names(out) <- problem$coef_names
  out
}

# Patient-level posterior triples (row, origin, destination, weight) for the
# final weighted refit.
markov_ml_patient_triples <- function(problem, state) {
  comp <- problem$complete_part
  pieces <- list(data.frame(
    row = comp$rows,
    origin = problem$origins[comp$origin],
    state = comp$state,
    weight = 1
  ))
  if (problem$chain$n) {
    for (ai in seq_len(problem$n_origin)) {
      xi <- markov_ml_chain_xi(problem, state, ai)
      idx <- which(xi > 0, arr.ind = TRUE)
      pieces[[length(pieces) + 1L]] <- data.frame(
        row = problem$chain$rows[idx[, 1L]],
        origin = problem$origins[[ai]],
        state = idx[, 2L],
        weight = xi[idx]
      )
    }
  }
  out <- do.call(rbind, pieces)
  out[order(out$row, out$origin, out$state), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# EM and Newton finish
# ---------------------------------------------------------------------------

markov_ml_em <- function(problem, theta, control) {
  state <- markov_ml_estep(problem, theta)
  trace <- state$loglik
  converged <- FALSE
  iterations <- 0L
  for (iter in seq_len(control$maxit)) {
    iterations <- iter
    triples <- markov_ml_aggregate_weights(problem, state)
    theta_new <- markov_ml_mstep(problem, triples, theta)
    state_new <- markov_ml_estep(problem, theta_new)
    trace <- c(trace, state_new$loglik)
    dll <- state_new$loglik - state$loglik
    dtheta <- max(abs(theta_new - theta))
    if (isTRUE(control$trace)) {
      message(sprintf(
        "EM iteration %d: loglik = %.10f, max |change| = %.3g",
        iter,
        state_new$loglik,
        dtheta
      ))
    }
    theta <- theta_new
    state <- state_new
    if (
      abs(dll) <= control$reltol * (abs(state$loglik) + control$reltol) ||
        dtheta < control$coef_tol
    ) {
      converged <- TRUE
      break
    }
  }
  list(
    theta = theta,
    state = state,
    trace = trace,
    iterations = iterations,
    converged = converged
  )
}

markov_ml_newton <- function(problem, theta, loglik, information, control) {
  steps <- 0L
  for (step in seq_len(control$newton_steps)) {
    score <- colSums(markov_ml_patient_scores(problem, theta))
    if (max(abs(score)) < 1e-12) {
      break
    }
    direction <- tryCatch(
      drop(solve(information, score)),
      error = function(e) NULL
    )
    if (is.null(direction) || any(!is.finite(direction))) {
      break
    }
    t <- 1
    accepted <- FALSE
    while (t >= 1 / 1024) {
      candidate <- theta + t * direction
      ll <- markov_ml_loglik(problem, candidate, strict = FALSE)
      if (is.finite(ll) && ll >= loglik - 1e-12 * max(1, abs(loglik))) {
        accepted <- TRUE
        break
      }
      t <- t / 2
    }
    if (!accepted) {
      break
    }
    theta <- candidate
    loglik <- ll
    steps <- steps + 1L
  }
  list(theta = theta, loglik = loglik, steps = steps)
}
