# Test-only pure-R references for the native kernels in src/sops.cpp.
#
# Each reference is written from the mathematics of the calculation and never
# calls package code that reaches C++. They must never be used as production
# fallbacks.

# Reverse cumulative logit category probabilities.
#
# `eta` has one row per prediction and one column per threshold j = 1..M, with
# P(Y >= j + 1) = plogis(eta[, j]). Category probabilities are differences of
# the survival function P(Y >= l), with P(Y >= 1) = 1 and P(Y >= K + 1) = 0.
reference_category_probabilities <- function(eta) {
  eta <- as.matrix(eta)
  survival <- cbind(1, stats::plogis(eta), 0)
  survival[, -ncol(survival), drop = FALSE] - survival[, -1L, drop = FALSE]
}

# Point-estimate and posterior kernels clip crossed probabilities at zero and
# rescale every row to unit mass.
reference_clipped_probabilities <- function(eta) {
  probabilities <- reference_category_probabilities(eta)
  probabilities[probabilities < 0] <- 0
  probabilities / rowSums(probabilities)
}

# One first-order Markov step, written as S_t = S_{t-1} %*% T_t per patient.
#
# `previous` is n x K. `transitions[[k]]` is the n x K matrix of
# P(Y_t = l | Y_{t-1} = k) for each non-absorbing origin k. Absorbing states
# have identity rows in T_t, so their mass is carried over exactly once.
reference_markov_step <- function(previous, transitions, absorb) {
  n <- nrow(previous)
  K <- ncol(previous)
  out <- matrix(0, nrow = n, ncol = K)
  for (i in seq_len(n)) {
    transition_matrix <- diag(K)
    for (origin in setdiff(seq_len(K), absorb)) {
      transition_matrix[origin, ] <- transitions[[origin]][i, ]
    }
    out[i, ] <- previous[i, ] %*% transition_matrix
  }
  out
}

# Split stacked transition rows into per-origin blocks. Row block b contains
# the n rows for origin `origins[b]`.
reference_origin_blocks <- function(stacked, origins, n) {
  blocks <- vector("list", max(c(origins, 0L)))
  for (position in seq_along(origins)) {
    rows <- (position - 1L) * n + seq_len(n)
    blocks[[origins[position]]] <- stacked[rows, , drop = FALSE]
  }
  blocks
}

# One second-order step on the joint probability J(h, j) of the states two
# visits ago (h) and one visit ago (j).
#
# J_t(j, l) = sum_h J_{t-1}(h, j) T_t(l | h, j) for non-absorbing j, and
# J_t(j, j) = sum_h J_{t-1}(h, j) for absorbing j. `transition_for(h, j)`
# returns the n x K transition probabilities for history (h, j).
reference_second_order_step <- function(joint, transition_for, absorb) {
  n <- dim(joint)[1L]
  K <- dim(joint)[2L]
  out <- array(0, dim = dim(joint))
  for (j in seq_len(K)) {
    if (j %in% absorb) {
      out[, j, j] <- rowSums(matrix(joint[, , j], nrow = n))
      next
    }
    for (h in seq_len(K)) {
      mass <- joint[, h, j]
      if (all(mass == 0)) {
        next
      }
      out[, j, ] <- out[, j, ] + mass * transition_for(h, j)
    }
  }
  out
}

# Normalize every (draw, observation) probability vector over states. Missing
# entries are ignored; vectors with infinite entries or non-positive mass are
# returned unchanged.
reference_normalize_probability_array <- function(values) {
  out <- values
  dims <- dim(values)
  for (draw in seq_len(dims[1L])) {
    for (observation in seq_len(dims[2L])) {
      vector <- values[draw, observation, ]
      present <- !is.na(vector)
      total <- sum(vector[present])
      if (all(is.finite(vector[present])) && total > 0) {
        out[draw, observation, present] <- vector[present] / total
      }
    }
  }
  out
}

# Average draw-level SOP arrays [draw, observation, time, state] within groups.
reference_group_means <- function(values, groups, group_count) {
  dims <- dim(values)
  out <- array(0, dim = c(dims[1L], group_count, dims[3L], dims[4L]))
  for (group in seq_len(group_count)) {
    members <- which(groups == group)
    if (length(members) > 0L) {
      out[, group, , ] <- apply(
        values[, members, , , drop = FALSE],
        c(1L, 3L, 4L),
        mean
      )
    }
  }
  out
}

# Inverse-CDF categorical sampling: with U ~ Uniform(0, 1), select the
# smallest category k whose cumulative weight reaches U times the row total.
# The uniforms come from R's generator in row order. Partial sums are
# accumulated left to right in double precision.
reference_sample_categorical_rows <- function(probabilities, uniforms) {
  vapply(
    seq_len(nrow(probabilities)),
    function(i) {
      cumulative <- Reduce(`+`, probabilities[i, ], accumulate = TRUE)
      target <- uniforms[i] * cumulative[length(cumulative)]
      selected <- which(target <= cumulative)
      if (length(selected) == 0L) {
        length(cumulative)
      } else {
        selected[1L]
      }
    },
    integer(1)
  )
}

# Random probability rows with strictly positive entries.
reference_random_stochastic <- function(rows, K) {
  weights <- matrix(stats::rexp(rows * K), nrow = rows, ncol = K)
  weights / rowSums(weights)
}

# Effective coefficients Gamma (thresholds x design columns) read directly from
# the backend, without the package's coefficient helpers.
reference_backend_gamma <- function(model, col_names) {
  if (inherits(model, "robcov_vglm")) {
    model <- model$vglm_fit
  }
  if (inherits(model, "vglm")) {
    gamma <- t(VGAM::coef(model, matrix = TRUE))
  } else if (inherits(model, "orm")) {
    coefficients <- stats::coef(model)
    M <- model$non.slopes
    slopes <- coefficients[-seq_len(M)]
    gamma <- cbind(
      "(Intercept)" = coefficients[seq_len(M)],
      matrix(slopes, nrow = M, ncol = length(slopes), byrow = TRUE)
    )
    colnames(gamma) <- c("(Intercept)", names(slopes))
  } else {
    stop("Unsupported reference backend.")
  }
  gamma[, col_names, drop = FALSE]
}

# First-order SOPs from compiled design matrices and effective coefficients.
# Pass `category_probabilities = reference_clipped_probabilities` for partial
# proportional-odds fits, whose point estimates clip crossed probabilities.
reference_first_order_plan_sops <- function(
  plan,
  gamma,
  category_probabilities = reference_category_probabilities
) {
  components <- plan$components
  n <- components$n_pat
  K <- components$n_states
  times <- length(plan$times)
  absorb <- which(plan$y_levels %in% plan$absorb)
  probabilities <- array(0, dim = c(n, times, K))
  current <- category_probabilities(
    components$X_init %*% t(gamma)
  )
  probabilities[, 1L, ] <- current
  if (times >= 2L) {
    for (visit in 2:times) {
      transitions <- reference_origin_blocks(
        category_probabilities(
          components$X_transition[[visit]] %*% t(gamma)
        ),
        components$transition_origins,
        n
      )
      current <- reference_markov_step(current, transitions, absorb)
      probabilities[, visit, ] <- current
    }
  }
  probabilities
}

# Second-order SOPs from a compiled second-order plan. Patients whose most
# recent observed state is absorbing start with all mass in that state.
reference_second_order_plan_sops <- function(plan, gamma) {
  components <- plan$components
  n <- components$n_pat
  K <- components$n_states
  times <- length(plan$times)
  absorb <- components$absorb
  initial <- reference_category_probabilities(components$X_init %*% t(gamma))
  probabilities <- array(0, dim = c(n, times, K))
  joint <- array(0, dim = c(n, K, K))
  for (i in seq_len(n)) {
    previous <- components$previous[i]
    if (previous %in% absorb) {
      joint[i, previous, previous] <- 1
    } else {
      joint[i, previous, ] <- initial[i, ]
    }
  }
  probabilities[, 1L, ] <- apply(joint, c(1L, 3L), sum)
  if (times >= 2L) {
    for (visit in 2:times) {
      transition <- reference_category_probabilities(
        components$X_transition[[visit]] %*% t(gamma)
      )
      pair_rows <- function(h, j) {
        pair <- which(components$older == h & components$current == j)
        if (length(pair) != 1L) {
          stop("Missing second-order history in the reference plan.")
        }
        transition[(pair - 1L) * n + seq_len(n), , drop = FALSE]
      }
      joint <- reference_second_order_step(joint, pair_rows, absorb)
      probabilities[, visit, ] <- apply(joint, c(1L, 3L), sum)
    }
  }
  probabilities
}
