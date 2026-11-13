# SOP fast-path matrix prediction.

# =============================================================================
# FAST PATH HELPER FUNCTIONS
# =============================================================================
# These internal functions provide the fast path for orm and vglm models: the
# design matrices of all visits are built once, and the recursion then runs
# for any coefficients without calling the models' predict() functions.

#' Build Fast Markov Simulation Components (Internal)
#'
#' Pre-calculates design matrices for fast Markov simulation. This is the
#' "heavy lifting" step that should be done ONCE before running thousands of
#' simulation draws.
#'
#' The first visit is predicted from each patient's own previous state(s).
#' Every later visit needs the transition probabilities from each previous
#' state a patient may be in: each non-absorbing state for a first-order model,
#' and each pair of non-absorbing states (two visits back, one visit back) for
#' a second-order model (`settings$p2_var`). The design of a later visit stacks
#' one block of patient rows per such combination, in the order of
#' `previous_states`.
#'
#' @param model Fitted vglm or orm model
#' @param newdata Baseline data frame (one row per patient)
#' @param times Time points to predict
#' @param settings The resolved settings of the call
#'   (`resolve_markov_settings()`): the state labels `y_levels`, the absorbing
#'   states `absorb` (their transition matrices are not needed, because
#'   absorbing states keep their prior probability mass), and the columns
#'   `time_var`, `p_var`, `p2_var` and `gap_var`.
#'
#' @return A list containing pre-calculated matrices and metadata.
#'
#' @keywords internal
markov_msm_build <- function(model, newdata, times = NULL, settings) {
  if (inherits(model, "robcov_vglm")) {
    model <- model$vglm_fit
  }
  is_orm <- inherits(model, "orm")
  if (!is_orm && !inherits(model, "vglm")) {
    mostr_abort("markov_msm_build() supports only vglm and orm models.")
  }
  Gamma <- get_effective_coefs(model)
  validate_sop_threshold_count(Gamma, settings$y_levels)
  time_res <- resolve_sop_times(
    model,
    newdata,
    times,
    settings$time_var,
    default = "fast"
  )
  times <- time_res$times

  n_pat <- nrow(newdata)
  ylevel_names <- as_state_labels(settings$y_levels)
  absorb_idx <- which(ylevel_names %in% as_state_labels(settings$absorb))
  non_absorb_idx <- setdiff(seq_along(ylevel_names), absorb_idx)
  # One row per combination of previous states, as state indices; with two
  # columns the state two visits back changes fastest.
  state_vars <- c(settings$p2_var, settings$p_var)
  previous_states <- expand.grid(
    rep(list(non_absorb_idx), length(state_vars)),
    KEEP.OUT.ATTRS = FALSE
  )
  names(previous_states) <- state_vars
  budget <- preflight_execution_plan_designs(
    n_pat = n_pat,
    n_times = length(times),
    rows_per_transition = n_pat * nrow(previous_states),
    n_cols = ncol(Gamma)
  )
  missing_vars <- setdiff(state_vars, names(newdata))
  if (length(missing_vars) > 0L) {
    mostr_abort(
      "Previous-state variable ",
      paste0("`", missing_vars, "`", collapse = ", "),
      " not found in `newdata`."
    )
  }

  # The design matrix columns used by the effective coefficient matrix. For
  # inline transforms such as rms::rcs(), model.matrix() uses the fitted terms
  # object's predvars, including stored knot locations.
  terms <- if (!is_orm) stats::delete.response(stats::terms(model))
  design <- function(data) {
    X <- if (is_orm) {
      orm_model_matrix(model, data, include_intercept = TRUE)
    } else {
      vglm_prediction_design(model, data, terms)
    }
    missing_cols <- setdiff(colnames(Gamma), colnames(X))
    if (length(missing_cols) > 0L) {
      mostr_abort(
        "Design matrix columns missing from the prediction design: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    X[, colnames(Gamma), drop = FALSE]
  }
  visit_rows <- function(index) {
    assign_sop_visit(
      newdata,
      settings = settings,
      times = times,
      index = index,
      time_info = time_res$time_info
    )
  }
  # The patient rows once per combination in `combinations` (rows of
  # `previous_states`), with the previous-state columns set to it.
  with_previous_states <- function(data, combinations) {
    if (length(combinations) > 1L) {
      data <- data[
        rep(seq_len(n_pat), times = length(combinations)),
        ,
        drop = FALSE
      ]
    }
    for (column in state_vars) {
      data[[column]] <- make_previous_state_column(
        states = ylevel_names[previous_states[[column]][combinations]],
        prototype = newdata[[column]],
        n = n_pat,
        p_var = column
      )
    }
    data
  }

  # A visit's design is built in one piece up to this many rows, and one
  # combination of previous states at a time above it, with identical results.
  # Measured Nov 2026 (orm, 11 design columns, Ryzen AI 9 HX 370): building in
  # one piece is faster only below about 25k rows; from 50k it is 1.3 to 1.9
  # times slower with up to 1.8 times the peak memory (at 2e6, 8.6 s and
  # 0.93 GB against 5.2 s and 0.51 GB). The default was 2e6 until then. See
  # dev/decisions.md.
  in_one_piece <- n_pat * nrow(previous_states) <=
    getOption("mostr.max_batched_design_rows", 5e4)
  pieces <- if (in_one_piece) {
    list(seq_len(nrow(previous_states)))
  } else {
    as.list(seq_len(nrow(previous_states)))
  }

  X_init <- design(visit_rows(1L))
  X_transition <- vector("list", length(times))
  for (index in seq_along(times)[-1L]) {
    data_visit <- visit_rows(index)
    X_transition[[index]] <- do.call(
      rbind,
      lapply(pieces, function(combinations) {
        design(with_previous_states(data_visit, combinations))
      })
    )
  }

  list(
    X_init = X_init,
    X_transition = X_transition,
    times = times,
    previous_states = previous_states,
    n_pat = n_pat,
    n_states = length(ylevel_names),
    M = length(ylevel_names) - 1,
    y_levels = ylevel_names,
    col_names = colnames(X_init),
    transition_layout = "stacked",
    transition_origins = non_absorb_idx,
    workspace_bytes = budget$workspace_bytes,
    design_bytes = budget$design_bytes
  )
}


#' Run Fast Markov Simulation (Internal)
#'
#' Runs the Markov loop using pre-calculated components and a specific
#' coefficient vector.
#'
#' @param components List returned by `markov_msm_build`
#' @param Gamma Effective coefficient matrix with one row per linear predictor
#'   and one column per small-design column, including the intercept.
#' @param times Vector of time points
#' @param absorb Absorbing state(s). Can be NULL for no absorbing states.
#'
#' @return Array of state probabilities with dimensions n_pat x n_times x n_states
#'
#' @keywords internal
markov_msm_run <- function(components, Gamma, times, absorb = NULL) {
  # Unpack
  X_init <- components$X_init
  X_transition <- components$X_transition
  n_pat <- components$n_pat
  n_states <- components$n_states
  M <- components$M
  y_levels <- components$y_levels
  col_names <- components$col_names
  n_times <- length(times)

  # Align Gamma to components
  Gamma <- Gamma[, col_names, drop = FALSE]
  Gamma_t <- t(Gamma) # [P x M]

  # Helper: X * Gamma_t -> LP [N x M]
  calc_lp <- function(X) {
    X %*% Gamma_t
  }

  absorb_idx <- if (!is.null(absorb)) {
    which(as.character(y_levels) %in% as.character(absorb))
  } else {
    integer(0)
  }
  non_absorb_idx <- setdiff(1:n_states, absorb_idx)

  po <- markov_po_structure(Gamma, col_names, X_init)
  initial <- if (is.null(po)) {
    lp_to_probs(calc_lp(X_init), M)
  } else {
    scalar <- drop(X_init %*% po$beta)
    lp_to_probs(outer(scalar, po$cutpoints, "+"), M)
  }
  check_transition_probabilities(initial, "the first SOP time point")
  check_crossed_thresholds(initial, n_pat, times[1L])
  P_out <- array(0, dim = c(n_pat, n_times, n_states))
  P_out[, 1L, ] <- initial
  previous <- initial

  if (n_times >= 2L) {
    for (t_idx in 2:n_times) {
      X_visit <- X_transition[[t_idx]]
      if (is.null(po)) {
        logits <- calc_lp(X_visit)
        check_transition_probabilities(logits, paste0("SOP time point ", t_idx))
        previous <- markov_update_logits_native(
          previous,
          logits,
          as.integer(non_absorb_idx),
          as.integer(absorb_idx),
          time = times[t_idx]
        )
      } else {
        scalar <- drop(X_visit %*% po$beta)
        check_transition_probabilities(scalar, paste0("SOP time point ", t_idx))
        previous <- markov_update_po_native(
          previous,
          scalar,
          po$cutpoints,
          as.integer(non_absorb_idx),
          as.integer(absorb_idx),
          time = times[t_idx]
        )
      }
      P_out[, t_idx, ] <- previous
    }
  }
  dimnames(P_out) <- list(NULL, NULL, y_levels)
  P_out
}

markov_po_structure <- function(Gamma, col_names, X_init) {
  intercept <- grep("Intercept", col_names, fixed = TRUE)
  if (length(intercept) != 1L || any(X_init[, intercept] != 1)) {
    return(NULL)
  }
  slopes <- setdiff(seq_along(col_names), intercept)
  if (
    length(slopes) &&
      any(vapply(
        slopes,
        function(j) !all(Gamma[, j] == Gamma[1L, j]),
        logical(1)
      ))
  ) {
    return(NULL)
  }
  beta <- Gamma[1L, ]
  beta[intercept] <- 0
  list(
    beta = beta,
    cutpoints = as.numeric(Gamma[, intercept]),
    intercept = intercept
  )
}


#' Convert Cumulative Log-Odds to Probabilities (Internal)
#'
#' Converts a matrix of cumulative log-odds to category probabilities.
#' Negative values from rounding (down to `-crossed_threshold_tolerance()`)
#' are set to zero and the row is rescaled. Lower values mean that the
#' cut-points of a partial proportional-odds model cross; they are kept so that
#' the caller reports them with `check_crossed_thresholds()`.
#'
#' @param eta Matrix with one row per prediction and one column per ordinal
#'   threshold.
#' @param M Number of thresholds (one less than the number of categories)
#' @return Matrix with one row per prediction and one column per outcome state.
#'
#' @keywords internal
lp_to_probs <- function(eta, M) {
  cum_probs <- stats::plogis(eta) # [N x M]
  n_rows <- nrow(eta)
  probs <- matrix(0, nrow = n_rows, ncol = M + 1)

  probs[, 1] <- 1 - cum_probs[, 1]

  if (M > 1) {
    probs[, 2:M] <- cum_probs[, 1:(M - 1)] - cum_probs[, 2:M]
  }

  probs[, M + 1] <- cum_probs[, M]
  rounding <- !is.na(probs) &
    probs < 0 &
    probs >= -crossed_threshold_tolerance()
  probs[rounding] <- 0

  normalize_probability_rows(probs)
}

normalize_probability_rows <- function(probs) {
  totals <- rowSums(probs, na.rm = TRUE)
  valid <- is.finite(totals) & totals > 0
  probs[valid, ] <- probs[valid, , drop = FALSE] / totals[valid]
  probs
}

