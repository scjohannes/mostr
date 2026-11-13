run_sop_delta_r_oracle <- function(plan, model) {
  coef <- sop_delta_raw_coef(model)
  map <- get_effective_coef_map(model)
  validated <- sop_delta_validate_plan(plan, model, map, coef)
  components <- validated$components
  Gamma <- validated$Gamma
  map <- validated$map

  n <- components$n_pat
  n_times <- length(plan$times)
  n_states <- components$n_states
  n_coef <- length(coef)
  absorb <- which(
    as.character(plan$y_levels) %in% as.character(plan$absorb)
  )
  non_absorb <- setdiff(seq_len(n_states), absorb)

  initial <- oracle_raw_probabilities(
    components$X_init,
    Gamma,
    map,
    "the first SOP time point"
  )
  probabilities <- array(0, dim = c(n, n_times, n_states))
  jacobian <- array(0, dim = c(n, n_times, n_states, n_coef))
  probabilities[, 1L, ] <- initial$probabilities
  jacobian[, 1L, , ] <- initial$derivative
  current_probability <- initial$probabilities
  current_jacobian <- initial$derivative

  if (n_times >= 2L) {
    for (visit in 2L:n_times) {
      X_visit <- oracle_visit_design(components, visit, non_absorb)
      next_probability <- matrix(0, nrow = n, ncol = n_states)
      next_jacobian <- array(0, dim = c(n, n_states, n_coef))

      for (origin in absorb) {
        next_probability[, origin] <- current_probability[, origin]
        next_jacobian[, origin, ] <- current_jacobian[, origin, ]
      }
      for (origin_position in seq_along(non_absorb)) {
        origin <- non_absorb[origin_position]
        rows <- (origin_position - 1L) * n + seq_len(n)
        transition <- oracle_raw_probabilities(
          X_visit[rows, , drop = FALSE],
          Gamma,
          map,
          paste0(
            "SOP time point ",
            visit,
            ", origin state ",
            plan$y_levels[origin]
          )
        )
        origin_probability <- current_probability[, origin]
        for (destination in seq_len(n_states)) {
          destination_probability <- transition$probabilities[, destination]
          next_probability[, destination] <-
            next_probability[, destination] +
            origin_probability * destination_probability
          for (coefficient in seq_len(n_coef)) {
            next_jacobian[, destination, coefficient] <-
              next_jacobian[, destination, coefficient] +
              current_jacobian[, origin, coefficient] *
                destination_probability +
              origin_probability *
                transition$derivative[, destination, coefficient]
          }
        }
      }
      probabilities[, visit, ] <- next_probability
      jacobian[, visit, , ] <- next_jacobian
      current_probability <- next_probability
      current_jacobian <- next_jacobian
    }
  }

  dimnames(probabilities) <- list(NULL, NULL, plan$y_levels)
  dimnames(jacobian) <- list(NULL, NULL, plan$y_levels, names(coef))
  list(probabilities = probabilities, jacobian = jacobian)
}

# Independent matrix-product reference for the analytical SOP recursion.
#
# Unlike run_sop_delta_r_oracle(), this reference uses no package helpers: it
# takes the raw inputs of the native routine. `X_init` is n x P,
# `X_transition[[t]]` (t >= 2) stacks one n-row block per origin state in
# `origins`, `gamma` is the M x P matrix of effective coefficients, and `map` is
# the (M * P) x q matrix with vec(gamma) = map %*% theta (column-major).
#
# With eta_j = x' gamma[j, ] and F_j = plogis(eta_j) = P(Y >= j + 1), category
# probabilities are p_l = F_{l-1} - F_l with F_0 = 1 and F_K = 0. Their
# derivatives are dp_l / dtheta_r = dF_{l-1} - dF_l, where
# dF_j = F_j (1 - F_j) x' dgamma_r[j, ] and dgamma_r is column r of `map`.
#
# For each patient the SOP row vector follows S_1 = p(x_init) and
# S_t = S_{t-1} T_t, where row k of T_t holds the transition probabilities from
# non-absorbing state k and absorbing rows are unit vectors. The product rule
# gives J_t[, r] = J_{t-1}[, r] T_t + S_{t-1} dT_t[, , r], with zero derivative
# rows for absorbing states. Grouped output averages contiguous blocks of
# `group_size` patients.
reference_sop_delta <- function(
  X_init,
  X_transition,
  gamma,
  map,
  origins,
  absorb,
  group_size = 0L
) {
  n <- nrow(X_init)
  M <- nrow(gamma)
  K <- M + 1L
  P <- ncol(gamma)
  q <- ncol(map)
  times <- length(X_transition)
  non_absorb <- setdiff(seq_len(K), absorb)

  category <- function(X) {
    eta <- X %*% t(gamma)
    survival <- stats::plogis(eta)
    probabilities <- cbind(1, survival, 0)
    probabilities <- probabilities[, -(K + 1L), drop = FALSE] -
      probabilities[, -1L, drop = FALSE]
    derivative <- array(0, dim = c(nrow(X), K, q))
    for (r in seq_len(q)) {
      gamma_r <- matrix(map[, r], nrow = M, ncol = P)
      d_survival <- cbind(0, survival * (1 - survival) * (X %*% t(gamma_r)), 0)
      derivative[, , r] <- d_survival[, -(K + 1L), drop = FALSE] -
        d_survival[, -1L, drop = FALSE]
    }
    list(probabilities = probabilities, derivative = derivative)
  }

  probabilities <- array(0, dim = c(n, times, K))
  jacobian <- array(0, dim = c(n, times, K, q))
  initial <- category(X_init)
  S <- initial$probabilities
  J <- initial$derivative
  probabilities[, 1L, ] <- S
  jacobian[, 1L, , ] <- J

  if (times >= 2L) {
    for (visit in 2:times) {
      X <- X_transition[[visit]]
      by_origin <- vector("list", K)
      for (origin in non_absorb) {
        position <- match(origin, origins)
        if (is.na(position)) {
          stop("The reference transition design lacks an origin block.")
        }
        rows <- (position - 1L) * n + seq_len(n)
        by_origin[[origin]] <- category(X[rows, , drop = FALSE])
      }
      S_next <- matrix(0, nrow = n, ncol = K)
      J_next <- array(0, dim = c(n, K, q))
      for (i in seq_len(n)) {
        transition <- diag(K)
        d_transition <- array(0, dim = c(K, K, q))
        for (origin in non_absorb) {
          transition[origin, ] <- by_origin[[origin]]$probabilities[i, ]
          d_transition[origin, , ] <- by_origin[[origin]]$derivative[i, , ]
        }
        S_next[i, ] <- S[i, ] %*% transition
        for (r in seq_len(q)) {
          J_next[i, , r] <- J[i, , r] %*% transition +
            S[i, ] %*% d_transition[, , r]
        }
      }
      S <- S_next
      J <- J_next
      probabilities[, visit, ] <- S
      jacobian[, visit, , ] <- J
    }
  }

  average <- function(x) {
    if (group_size == 0L) {
      return(x)
    }
    groups <- rep(seq_len(n / group_size), each = group_size)
    dims <- dim(x)
    means <- rowsum(matrix(x, nrow = dims[1L]), groups) / group_size
    array(means, dim = c(nrow(means), dims[-1L]))
  }
  list(
    probabilities = average(probabilities),
    jacobian = average(jacobian),
    individual_probabilities = probabilities
  )
}

# Effective coefficients and their raw-coefficient map, read from the fitted
# backend without the package's map construction. VGLM maps are recovered by
# evaluating VGAM's coefficient matrix at each raw-coefficient basis vector;
# ORM maps follow from eta_j = alpha_j + x' beta.
reference_delta_model_inputs <- function(model, col_names) {
  if (inherits(model, "robcov_vglm")) {
    model <- model$vglm_fit
  }
  if (inherits(model, "vglm")) {
    theta <- model@coefficients
    full_gamma <- t(VGAM::coef(model, matrix = TRUE))
    map <- vapply(
      seq_along(theta),
      function(r) {
        basis <- model
        basis@coefficients[] <- 0
        basis@coefficients[r] <- 1
        as.vector(t(VGAM::coef(basis, matrix = TRUE)))
      },
      numeric(length(full_gamma))
    )
  } else if (inherits(model, "orm")) {
    theta <- stats::coef(model)
    M <- model$non.slopes
    slope_names <- names(theta)[-seq_len(M)]
    full_gamma <- matrix(
      0,
      nrow = M,
      ncol = 1L + length(slope_names),
      dimnames = list(NULL, c("(Intercept)", slope_names))
    )
    map <- matrix(0, nrow = length(full_gamma), ncol = length(theta))
    map[seq_len(M), seq_len(M)] <- diag(M)
    for (slope in seq_along(slope_names)) {
      map[slope * M + seq_len(M), M + slope] <- 1
    }
  } else {
    stop("Unsupported reference backend.")
  }
  M <- nrow(full_gamma)
  columns <- match(col_names, colnames(full_gamma))
  if (anyNA(columns)) {
    stop("The reference coefficient map lacks plan design columns.")
  }
  rows <- unlist(lapply(columns, function(column) {
    (column - 1L) * M + seq_len(M)
  }))
  map <- map[rows, , drop = FALSE]
  list(
    theta = theta,
    map = map,
    gamma = matrix(drop(map %*% theta), nrow = M, ncol = length(columns))
  )
}

# Category probabilities and their derivatives with respect to the model
# coefficients for one design matrix, and the transition rows of one visit.
# Only run_sop_delta_r_oracle() uses them; the package computes both in the
# compiled analytical recursion.
oracle_raw_probabilities <- function(X, Gamma, map, context) {
  n <- nrow(X)
  M <- nrow(Gamma)
  K <- M + 1L
  q <- ncol(map)
  P <- ncol(Gamma)

  eta <- X %*% t(Gamma)
  cumulative <- stats::plogis(eta)
  probabilities <- matrix(0, nrow = n, ncol = K)
  probabilities[, 1L] <- 1 - cumulative[, 1L]
  if (M > 1L) {
    probabilities[, 2L:M] <-
      cumulative[, seq_len(M - 1L), drop = FALSE] -
      cumulative[, 2L:M, drop = FALSE]
  }
  probabilities[, K] <- cumulative[, M]

  if (any(!is.finite(probabilities))) {
    sop_delta_condition(
      paste0(
        "Non-finite raw category probabilities occurred at ",
        context,
        "."
      ),
      "mostr_delta_invalid_probabilities",
      context = context
    )
  }
  minimum <- min(probabilities)
  if (minimum < -1e-10) {
    sop_delta_condition(
      paste0(
        "Raw ordinal category probabilities are crossed at ",
        context,
        " (minimum ",
        format(minimum, digits = 6),
        ")."
      ),
      "mostr_delta_crossed_probabilities",
      context = context,
      minimum_probability = minimum
    )
  }
  if (max(abs(rowSums(probabilities) - 1)) > 1e-12) {
    sop_delta_condition(
      paste0("Raw category probabilities do not sum to one at ", context, "."),
      "mostr_delta_invalid_probabilities",
      context = context
    )
  }

  derivative <- array(0, dim = c(n, K, q))
  logistic_derivative <- cumulative * (1 - cumulative)
  for (coefficient in seq_len(q)) {
    gamma_derivative <- matrix(
      map[, coefficient],
      nrow = M,
      ncol = P
    )
    eta_derivative <- X %*% t(gamma_derivative)
    cumulative_derivative <- logistic_derivative * eta_derivative
    derivative[, 1L, coefficient] <- -cumulative_derivative[, 1L]
    if (M > 1L) {
      derivative[, 2L:M, coefficient] <-
        cumulative_derivative[, seq_len(M - 1L), drop = FALSE] -
        cumulative_derivative[, 2L:M, drop = FALSE]
    }
    derivative[, K, coefficient] <- cumulative_derivative[, M]
  }
  derivative_mass <- apply(derivative, c(1L, 3L), sum)
  derivative_scale <- apply(abs(derivative), c(1L, 3L), sum)
  if (any(abs(derivative_mass) > 1e-12 + 1e-12 * derivative_scale)) {
    stop(
      "Internal analytical category derivatives do not sum to zero.",
      call. = FALSE
    )
  }

  list(probabilities = probabilities, derivative = derivative)
}

oracle_visit_design <- function(components, visit, non_absorb) {
  X <- components$X_transition[[visit]]
  if (!is.matrix(X)) {
    stop("Missing transition design for SOP visit ", visit, ".", call. = FALSE)
  }
  origins <- components$transition_origins
  origin_positions <- match(non_absorb, origins)
  if (anyNA(origin_positions)) {
    stop("Transition design is missing a required origin state.", call. = FALSE)
  }
  rows <- unlist(lapply(origin_positions, function(position) {
    start <- (position - 1L) * components$n_pat + 1L
    start:(start + components$n_pat - 1L)
  }))
  if (nrow(X) < max(rows) || ncol(X) != length(components$col_names)) {
    stop(
      "The transition design has dimensions inconsistent with the SOP plan.",
      call. = FALSE
    )
  }
  X[rows, , drop = FALSE]
}

# R reference for the compiled averaging in analytical inference: averages
# per-profile SOPs and their Jacobians over the profiles of each scenario.
# `probabilities` and `jacobian` hold one row per scenario-profile pair, in
# grid blocks of `n_each` profiles.
delta_average_profile_arrays <- function(
  probabilities,
  jacobian,
  grid,
  times,
  y_levels,
  variables,
  n_each
) {
  n_cf <- nrow(grid)
  q <- dim(jacobian)[4L]
  cells <- marginalize_sops_array(
    sops_array = probabilities,
    grid = grid,
    times = times,
    y_levels = y_levels,
    variables = variables,
    n_cf = n_cf,
    n_each = n_each
  )
  average_jacobian <- matrix(NA_real_, nrow = nrow(cells), ncol = q)
  individual_values <- matrix(NA_real_, nrow = n_each, ncol = nrow(cells))
  n_cells <- length(times) * length(y_levels)
  cursor <- 1L
  for (scenario in seq_len(n_cf)) {
    rows <- ((scenario - 1L) * n_each + 1L):(scenario * n_each)
    output_rows <- cursor:(cursor + n_cells - 1L)
    average_jacobian[output_rows, ] <- matrix(
      apply(jacobian[rows, , , , drop = FALSE], c(2L, 3L, 4L), mean),
      nrow = n_cells,
      ncol = q
    )
    individual_values[, output_rows] <- matrix(
      probabilities[rows, , , drop = FALSE],
      nrow = n_each,
      ncol = n_cells
    )
    cursor <- cursor + n_cells
  }
  list(
    cells = cells,
    jacobian = average_jacobian,
    individual_values = individual_values
  )
}
