# Internal transition traces for model-based diagnostic plots.

markov_transition_trace <- function(
  object,
  data,
  times = NULL,
  settings,
  return_kernels = FALSE,
  ...
) {
  dots <- list(...)
  unknown_dots <- setdiff(
    names(dots),
    c(".draw_indices", ".gamma_draws", ".n_draws_total")
  )
  if (length(unknown_dots) > 0) {
    mostr_abort("Unused arguments: ", paste(unknown_dots, collapse = ", "))
  }
  draw_indices <- dots$.draw_indices
  gamma_draws_arg <- dots$.gamma_draws
  # A caller that traces posterior draws in chunks passes the call's number of
  # draws, which the crossed cut-point error reports.
  n_draws_total <- dots$.n_draws_total

  cl <- if (inherits(object, "blrm")) {
    "blrm"
  } else if (inherits(object, "robcov_vglm")) {
    "robcov_vglm"
  } else if (inherits(object, "orm")) {
    "orm"
  } else if (inherits(object, "vglm")) {
    "vglm"
  } else {
    class(object)[1]
  }
  ftypes <- c(
    orm = "rms",
    blrm = "rmsb",
    vglm = "vglm",
    robcov_vglm = "robcov"
  )
  ftype <- ftypes[cl]

  if (is.na(ftype)) {
    mostr_abort("Object class not supported")
  }

  p_var <- settings$p_var
  p2_var <- settings$p2_var
  if (!p_var %in% names(data)) {
    mostr_abort("Previous-state variable `", p_var, "` not found in `data`.")
  }
  if (!is.null(p2_var) && !p2_var %in% names(data)) {
    mostr_abort(
      "Second previous-state variable `",
      p2_var,
      "` not found in `data`."
    )
  }

  time_res <- resolve_sop_times(
    object,
    data,
    times,
    settings$time_var,
    default = "unique"
  )
  times <- time_res$times
  time_info <- time_res$time_info

  prd <- switch(
    ftype,
    rms = function(obj, d) predict_orm_response_markov(obj, d),
    vglm = function(obj, d) predict_vglm_response_markov(obj, d),
    rmsb = function(obj, d) {
      predict_blrm_response_markov(
        obj,
        d,
        include_re = FALSE,
        id_var = settings$id_var,
        draw_indices = draw_indices,
        gamma_draws = gamma_draws_arg
      )
    },
    robcov = function(obj, d) {
      if (is.null(obj$vglm_fit)) {
        mostr_abort(
          "robcov_vglm object does not contain the original vglm fit. ",
          "Please re-run robcov_vglm() with the latest version of the package."
        )
      }
      predict_vglm_response_markov(obj$vglm_fit, d)
    }
  )

  n_pat <- nrow(data)
  n_times <- length(times)
  ylevel_names <- as_state_labels(settings$y_levels)
  absorb_names <- as_state_labels(settings$absorb)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  nd <- if (ftype == "rmsb" && length(object$draws)) {
    length(draw_indices)
  } else {
    0L
  }

  # Every array has the posterior draws as its first dimension. A fit without
  # draws is traced as a single draw, and `markov_transition_trace_result()`
  # drops that dimension at the end.
  n_layers <- max(nd, 1L)
  trace_array <- function(n_state_dims) {
    array(
      0,
      dim = c(n_layers, n_pat, n_times, rep(n_states, n_state_dims)),
      dimnames = c(
        list(if (nd > 0) draw_indices, rownames(data), times),
        rep(list(ylevel_names), n_state_dims)
      )
    )
  }
  # [draw, patient, time, state]
  P <- trace_array(1L)
  # [draw, patient, time, previous state, state]
  transitions <- trace_array(2L)
  # As `transitions`; a second-order model has the state before the previous
  # one in front of the previous state.
  kernels <- if (isTRUE(return_kernels)) {
    trace_array(if (is.null(p2_var)) 2L else 3L)
  }

  data <- assign_sop_visit(
    data,
    settings = settings,
    times = times,
    index = 1L,
    time_info = time_info
  )

  prev_idx <- match_state_indices(data[[p_var]], ylevel_names, p_var)
  predict_rows <- !prev_idx %in% absorb_idx

  if (any(predict_rows)) {
    p_t1 <- prd(object, data[predict_rows, , drop = FALSE])
    check_transition_probabilities(p_t1, "the first transition time point")
    check_crossed_thresholds(
      p_t1,
      n_pat,
      times[1L],
      rows = which(predict_rows),
      n_draws = n_draws_total
    )
    p_t1 <- trace_draw_layers(p_t1, nd)
    P[, predict_rows, 1, ] <- p_t1
    rows <- which(predict_rows)
    for (pos in seq_along(rows)) {
      row <- rows[pos]
      transitions[, row, 1L, prev_idx[row], ] <- p_t1[, pos, ]
    }
  }
  # A patient who starts in an absorbing state stays there.
  for (row in which(!predict_rows)) {
    state <- prev_idx[row]
    P[, row, 1L, state] <- 1
    transitions[, row, 1L, state, state] <- 1
  }

  trace_later_visits <- if (is.null(p2_var)) {
    markov_transition_trace_first_order
  } else {
    markov_transition_trace_second_order
  }
  trace_later_visits(
    P = P,
    transitions = transitions,
    kernels = kernels,
    object = object,
    data = data,
    prd = prd,
    nd = nd,
    times = times,
    ylevel_names = ylevel_names,
    absorb_names = absorb_names,
    settings = settings,
    time_info = time_info,
    n_draws_total = n_draws_total
  )
}

# Predicted probabilities as [draw, row, state]: a fit without posterior draws
# (`nd` of 0) returns [row, state], which becomes one draw.
trace_draw_layers <- function(probs, nd) {
  if (nd > 0) {
    return(probs)
  }
  array(probs, dim = c(1L, dim(probs)))
}

# The trace as its callers read it: with the draws as first dimension for a
# fit with posterior draws, and without that dimension otherwise.
markov_transition_trace_result <- function(P, transitions, kernels, nd) {
  drop_layer <- function(x) {
    if (nd > 0 || is.null(x)) {
      return(x)
    }
    array(x, dim = dim(x)[-1L], dimnames = dimnames(x)[-1L])
  }
  list(
    sops = drop_layer(P),
    transitions = drop_layer(transitions),
    kernels = drop_layer(kernels)
  )
}

markov_transition_trace_first_order <- function(
  P,
  transitions,
  kernels,
  object,
  data,
  prd,
  nd,
  times,
  ylevel_names,
  absorb_names,
  settings,
  time_info,
  n_draws_total = NULL
) {
  p_var <- settings$p_var
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)
  yna <- ylevel_names[non_absorb_idx]

  if (n_times < 2) {
    return(markov_transition_trace_result(P, transitions, kernels, nd))
  }

  edata_base <- data[rep(seq_len(n_pat), times = length(yna)), , drop = FALSE]
  edata_base[[p_var]] <- make_previous_state_column(
    states = yna,
    prototype = data[[p_var]],
    n = n_pat,
    p_var = p_var
  )

  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      settings = settings,
      times = times,
      index = it,
      time_info = time_info
    )

    trans_probs <- prd(object, edata_base)
    check_transition_probabilities(
      trans_probs,
      paste0("transition time point ", it)
    )
    check_crossed_thresholds(
      trans_probs,
      n_pat,
      times[it],
      n_draws = n_draws_total
    )
    trans_probs <- trace_draw_layers(trans_probs, nd)

    # [draw, patient, previous state, state]
    transition_current <- array(
      0,
      dim = c(dim(P)[1L], n_pat, n_states, n_states)
    )
    for (i in seq_along(yna)) {
      prev_state_name <- yna[i]
      prev_state_idx <- match(prev_state_name, ylevel_names)
      prob_prev <- P[,, it - 1, prev_state_name]
      row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
      probs_transition <- trans_probs[, row_indices, , drop = FALSE]
      if (!is.null(kernels)) {
        kernels[,, it, prev_state_idx, ] <- probs_transition
      }
      for (k in seq_len(n_states)) {
        transition_current[,, prev_state_idx, k] <-
          probs_transition[,, k] * prob_prev
      }
    }
    for (a_state in absorb_idx) {
      transition_current[,, a_state, a_state] <- P[,, it - 1, a_state]
      if (!is.null(kernels)) {
        kernels[,, it, a_state, a_state] <- 1
      }
    }
    transitions[,, it, , ] <- transition_current
    P[,, it, ] <- apply(transition_current, c(1, 2, 4), sum)
  }

  markov_transition_trace_result(P, transitions, kernels, nd)
}

markov_transition_trace_second_order <- function(
  P,
  transitions,
  kernels,
  object,
  data,
  prd,
  nd,
  times,
  ylevel_names,
  absorb_names,
  settings,
  time_info,
  n_draws_total = NULL
) {
  p_var <- settings$p_var
  p2_var <- settings$p2_var
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)

  if (n_times < 2) {
    return(markov_transition_trace_result(P, transitions, kernels, nd))
  }

  # [draw, patient, state before the previous one, previous state]
  joint_dim <- c(dim(P)[1L], n_pat, n_states, n_states)
  joint_prev <- array(transitions[,, 1, , ], dim = joint_dim)

  pair_grid <- expand.grid(
    h = seq_len(n_states),
    j = seq_len(n_states),
    KEEP.OUT.ATTRS = FALSE
  )
  n_pairs <- nrow(pair_grid)
  predictable_pair <- pair_grid$h %in%
    non_absorb_idx &
    pair_grid$j %in% non_absorb_idx
  block_rows <- lapply(seq_len(n_pairs), function(i) {
    ((i - 1) * n_pat + 1):(i * n_pat)
  })

  edata_base <- data[rep(seq_len(n_pat), times = n_pairs), , drop = FALSE]
  edata_base[[p2_var]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$h],
    prototype = data[[p2_var]],
    n = n_pat,
    p_var = p2_var
  )
  edata_base[[p_var]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$j],
    prototype = data[[p_var]],
    n = n_pat,
    p_var = p_var
  )

  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      settings = settings,
      times = times,
      index = it,
      time_info = time_info
    )

    predict_rows <- unlist(block_rows[predictable_pair], use.names = FALSE)
    trans_probs <- prd(object, edata_base[predict_rows, , drop = FALSE])
    check_transition_probabilities(
      trans_probs,
      paste0("second-order transition time point ", it)
    )
    check_crossed_thresholds(
      trans_probs,
      n_pat,
      times[it],
      rows = predict_rows,
      n_draws = n_draws_total
    )
    trans_probs <- trace_draw_layers(trans_probs, nd)
    cursor <- 1L

    joint_current <- array(0, dim = joint_dim)
    for (pair_i in seq_len(n_pairs)) {
      h <- pair_grid$h[pair_i]
      j <- pair_grid$j[pair_i]
      if (j %in% absorb_idx) {
        joint_current[,, j, j] <- joint_current[,, j, j] +
          joint_prev[,, h, j]
        if (!is.null(kernels)) {
          kernels[,, it, h, j, j] <- 1
        }
      } else if (predictable_pair[pair_i]) {
        rows <- cursor:(cursor + n_pat - 1L)
        prob_prev <- joint_prev[,, h, j]
        transition <- trans_probs[, rows, , drop = FALSE]
        if (!is.null(kernels)) {
          kernels[,, it, h, j, ] <- transition
        }
        if (any(prob_prev > 0)) {
          for (l in seq_len(n_states)) {
            joint_current[,, j, l] <- joint_current[,, j, l] +
              transition[,, l] * prob_prev
          }
        }
        cursor <- cursor + n_pat
      }
    }
    transitions[,, it, , ] <- joint_current
    P[,, it, ] <- apply(joint_current, c(1, 2, 4), sum)
    joint_prev <- joint_current
  }

  markov_transition_trace_result(P, transitions, kernels, nd)
}
