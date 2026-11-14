# Data and a plain-R reference for fits with unobserved states.
#
# The reference, `ms_reference()`, takes the data frame itself (not the
# plan that `markov_path_plan()` builds), writes the design rows out by hand
# and sums over every sequence of states each patient can have passed through,
# by listing the sequences. It calls no package code, so it checks the plan
# builder, the kernel and the maximisation together.
#
# Model, for K ordered states 1..K and a transition into the state at visit t
# from the state k at the visit before:
#   x_t = (1(k == 2), ..., 1(k == K - 1), time_t, tx)
#   P(Y_t >= j + 1) = plogis(alpha_j + x_t'beta), j = 1..K-1,
#   P(Y_t = y) = P(Y_t >= y) - P(Y_t >= y + 1).
# State K is absorbing when `absorb = K`: from it the next state is K with
# probability 1. A patient who was never in state K before has no previous
# state K, so the previous-state dummies run over 2..K-1 only.
#
# A patient's states are s_0 (the previous state of the first row, `yprev`),
# s_1, ..., s_n (the states of the rows, in time order). A row's state is known
# to lie in lo_t..hi_t: lo == hi is an observed state; `NA` in the data means
# "any state that is not absorbing". The patient's likelihood conditions on the
# first exactly known state s_j (s_0 when `yprev` is given) and is
#   L = sum over (s_{j+1}, ..., s_n) in their ranges of
#       prod_{t > j} P(s_t | s_{t-1}, x_t).
# With K absorbing and s_{t-1} = K, the factor is 1(s_t == K).

# Simulated visits of `n` patients over `n_visits` days with K = 4 states, state
# 4 (death) absorbing, and no rows after death. Columns: id, time, tx, yprev
# (factor with levels 1..3), y (integer).
ms_simulate <- function(n = 60, n_visits = 6, seed = 20261114) {
  withr::with_seed(seed, {
    alpha <- c(1.5, -0.8, -3.2)
    beta <- c(yprev2 = 1.2, yprev3 = 2.6, time = -0.15, tx = -0.5)
    rows <- vector("list", n)
    for (i in seq_len(n)) {
      tx <- i %% 2
      state <- sample(1:3, 1, prob = c(0.4, 0.4, 0.2))
      out <- NULL
      for (t in seq_len(n_visits)) {
        x <- c(state == 2, state == 3, t, tx)
        cumulative <- c(1, stats::plogis(alpha + sum(x * beta)), 0)
        prob <- cumulative[1:4] - cumulative[2:5]
        new_state <- sample(1:4, 1, prob = prob)
        out <- rbind(out, data.frame(id = i, time = t, tx = tx, yprev = state, y = new_state))
        state <- new_state
        if (state == 4) break
      }
      rows[[i]] <- out
    }
  })
  data <- do.call(rbind, rows)
  data$yprev <- factor(data$yprev, levels = 1:3)
  rownames(data) <- NULL
  data
}

# Simulated visits with K states (K absorbing), three continuous covariates
# (age, bmi, sofa) and a fraction `miss` of the non-absorbing states set to
# `NA`, with the previous state of the next row set to `NA` as well. With
# many states and much missing, the starting values (the fit to the fully
# observed transitions) can lie where the log-likelihood with unobserved
# states is not concave.
ms_simulate_many_states <- function(n, n_visits, K, miss, seed) {
  withr::with_seed(seed, {
    alpha <- seq(3, -5, length.out = K - 1)
    b_prev <- seq(0, 4.5, length.out = K - 1)
    rows <- vector("list", n)
    for (i in seq_len(n)) {
      tx <- i %% 2
      age <- stats::rnorm(1)
      bmi <- stats::rnorm(1)
      sofa <- stats::rnorm(1)
      state <- sample(2:min(5, K - 1), 1)
      out <- vector("list", n_visits)
      for (t in seq_len(n_visits)) {
        eta <- b_prev[state] - 0.03 * t - 0.4 * tx + 0.3 * age - 0.2 * bmi +
          0.25 * sofa
        cumulative <- c(1, stats::plogis(alpha + eta), 0)
        new <- sample.int(K, 1, prob = cumulative[1:K] - cumulative[2:(K + 1)])
        out[[t]] <- data.frame(
          id = i, time = t, tx = tx, age = age, bmi = bmi, sofa = sofa,
          yprev = state, y = new
        )
        state <- new
        if (state == K) break
      }
      rows[[i]] <- do.call(rbind, out)
    }
    data <- do.call(rbind, rows)
    hide <- data$y != K & stats::runif(nrow(data)) < miss
  })
  data$y[hide] <- NA
  after <- c(FALSE, hide[-nrow(data)]) &
    c(FALSE, data$id[-1] == data$id[-nrow(data)])
  data$yprev[after] <- NA
  data$yprev <- factor(data$yprev, levels = 1:(K - 1))
  rownames(data) <- NULL
  data
}

# Set the state of the visits (ids[k], times[k]) to `NA` and the previous
# state of the visit after each to `NA`.
ms_hide <- function(data, ids, times) {
  for (k in seq_along(ids)) {
    row <- which(data$id == ids[[k]] & data$time == times[[k]])
    stopifnot(length(row) == 1L)
    data$y[row] <- NA
    nxt <- which(data$id == ids[[k]] & data$time == times[[k]] + 1)
    if (length(nxt)) data$yprev[nxt] <- NA
  }
  data
}

# The standard data set with unobserved visits: every kind of stretch the
# tests need. Patient numbers refer to `ms_simulate()` with its defaults.
ms_gap_data <- function(data = ms_simulate()) {
  last_time <- tapply(data$time, data$id, max)
  died <- tapply(data$y == 4, data$id, any)
  alive_long <- as.integer(names(last_time)[last_time == 6 & !died])
  dead_late <- as.integer(names(last_time)[died & last_time >= 4])
  # Interior gaps of length 1, 2 and 3 followed by an observed state.
  hide <- data.frame(id = alive_long[1], time = 3)
  hide <- rbind(hide, data.frame(id = alive_long[2], time = 2:3))
  hide <- rbind(hide, data.frame(id = alive_long[3], time = 2:4))
  # Two gaps in one patient.
  hide <- rbind(hide, data.frame(id = alive_long[4], time = c(2, 4)))
  # A gap that runs to the patient's last row.
  hide <- rbind(hide, data.frame(id = alive_long[5], time = 5:6))
  # Gaps followed by death.
  for (k in seq_len(min(2, length(dead_late)))) {
    j <- dead_late[k]
    hide <- rbind(hide, data.frame(id = j, time = last_time[[as.character(j)]] - 1))
  }
  # A gap of length 2 before death.
  if (length(dead_late) >= 3) {
    j <- dead_late[3]
    lt <- last_time[[as.character(j)]]
    hide <- rbind(hide, data.frame(id = j, time = c(lt - 2, lt - 1)))
  }
  out <- ms_hide(data, hide$id, hide$time)
  # A patient whose first two rows have unknown previous states: the first
  # previous state and the first state are not recorded.
  first <- alive_long[6]
  rows <- which(out$id == first)
  out$yprev[rows[1:2]] <- NA
  out$y[rows[1]] <- NA
  attr(out, "hidden") <- hide
  attr(out, "special") <- list(
    alive_long = alive_long,
    dead_late = dead_late,
    unknown_start = first,
    trailing = alive_long[5]
  )
  out
}

ms_fit <- function(data, ..., formula = ordered(y) ~ yprev + time + tx) {
  mostm(formula, data = data, id_var = "id", absorb = 4, ...)
}

# Shared fixtures, computed once per test run: a fit takes about 0.3 s, and
# most tests read the same data and fit. Tests must not modify them.
ms_fixtures <- new.env(parent = emptyenv())
ms_fixture <- function(name) {
  if (is.null(ms_fixtures[[name]])) {
    ms_fixtures[[name]] <- switch(
      name,
      complete = ms_simulate(),
      gaps = ms_gap_data(ms_fixture("complete")),
      gap_fit = ms_fit(ms_fixture("gaps")),
      # Without the patient whose first previous state is unknown, who has
      # no starting row for predictions.
      downstream = {
        data <- ms_fixture("gaps")
        data[data$id != attr(data, "special")$unknown_start, ]
      },
      downstream_fit = ms_fit(ms_fixture("downstream")),
      stop("unknown fixture ", name)
    )
  }
  ms_fixtures[[name]]
}

# Ranges lo..hi per row from `y` or from two columns `lower`, `upper` (with
# -Inf, Inf and NA as in `Ocens()`), and whether the row's state is recorded
# at all. A row without a record (`NA`, or `NA` in both bounds) is any state
# that is not absorbing when a later row of the patient is recorded, and any
# state at all after the patient's last recorded row.
ms_row_ranges <- function(data, K = 4, absorb = 4, lower = NULL, upper = NULL) {
  top <- if (length(absorb)) K - 1 else K
  if (is.null(lower)) {
    y <- as.integer(as.character(data$y))
    recorded <- !is.na(y)
    lo <- ifelse(recorded, y, NA)
    hi <- lo
  } else {
    recorded <- !(is.na(lower) & is.na(upper))
    lo <- ifelse(lower == -Inf, 1, lower)
    hi <- ifelse(upper == Inf, K, upper)
  }
  recorded_later <- vapply(seq_len(nrow(data)), function(i) {
    any(recorded & data$id == data$id[i] & data$time > data$time[i])
  }, logical(1))
  lo[!recorded] <- 1
  hi[!recorded] <- ifelse(recorded_later[!recorded], top, K)
  list(lo = as.integer(lo), hi = as.integer(hi), recorded = recorded)
}

# The reference log-likelihood for a data frame, as a function of `theta`
# (named like `coef(fit)`: intercepts "y>=2", ..., "y>=K", then "yprev=2", ...,
# "yprev=K-1", "time", "tx"; without the "yprev" terms when
# `previous_effect = FALSE`, for the formula `ordered(y) ~ time + tx`). The
# paths are listed once, here; the returned function
# `f(theta, per_patient = FALSE)` evaluates every transition of every path,
# multiplies along each path, adds the paths of each patient and takes the
# log. `per_patient = TRUE` returns one value per patient (named by id).
# `weights`: one weight per row, the same on all rows of a patient.
ms_reference <- function(
  data,
  K = 4,
  absorb = 4,
  ranges = ms_row_ranges(data, K, absorb),
  weights = NULL,
  previous_effect = TRUE
) {
  ids <- unique(data$id)
  steps <- list()
  path_count <- 0L
  patient_weight <- rep(1, length(ids))
  recorded <- ranges$recorded %||% rep(TRUE, nrow(data))
  for (i in seq_along(ids)) {
    rows <- which(data$id == ids[[i]])
    rows <- rows[order(data$time[rows])]
    if (!is.null(weights)) patient_weight[i] <- weights[rows[1]]
    n <- length(rows)
    s0 <- as.integer(as.character(data$yprev[rows[1]]))
    lo <- c(s0, ranges$lo[rows])
    hi <- c(s0, ranges$hi[rows])
    # A state is known exactly only on a recorded row (or as `yprev`); an
    # unrecorded row with a single possible state is still summed over.
    known <- which(!is.na(lo) & lo == hi & c(TRUE, recorded[rows]))
    # Nothing to condition on, or nothing after the first known state.
    if (!length(known) || known[1] > n) next
    positions <- known[1]:(n + 1) # 1 is s_0, t + 1 is row t
    paths <- as.matrix(expand.grid(lapply(positions, function(p) lo[p]:hi[p])))
    for (m in seq_along(positions)[-1]) {
      row <- rows[positions[m] - 1]
      steps[[length(steps) + 1L]] <- data.frame(
        patient = i,
        path = path_count + seq_len(nrow(paths)),
        previous = paths[, m - 1],
        current = paths[, m],
        time = data$time[row],
        tx = data$tx[row]
      )
    }
    path_count <- path_count + nrow(paths)
  }
  steps <- do.call(rbind, steps)
  path_patient <- integer(path_count)
  path_patient[steps$path] <- steps$patient
  # Design rows: previous-state dummies for states 2..K-1, time, tx.
  dummies <- if (previous_effect && K > 2) {
    outer(steps$previous, 2:(K - 1), "==") * 1
  }
  x <- cbind(dummies, steps$time, steps$tx)
  from_absorbing <- steps$previous %in% absorb

  function(theta, per_patient = FALSE) {
    alpha <- unname(theta[seq_len(K - 1)])
    beta <- unname(theta[-seq_len(K - 1)])
    eta <- as.numeric(x %*% beta)
    cumulative <- cbind(1, stats::plogis(outer(eta, alpha, "+")), 0)
    prob <- cumulative[cbind(seq_along(eta), steps$current)] -
      cumulative[cbind(seq_along(eta), steps$current + 1L)]
    prob[from_absorbing] <- as.numeric(
      steps$current[from_absorbing] == steps$previous[from_absorbing]
    )
    path_prob <- exp(rowsum(log(prob), steps$path, reorder = TRUE)[, 1])
    patient_lik <- rowsum(path_prob, path_patient, reorder = TRUE)[, 1]
    out <- numeric(length(ids))
    names(out) <- as.character(ids)
    used <- as.integer(names(patient_lik))
    out[used] <- patient_weight[used] * log(patient_lik)
    if (per_patient) out else sum(out)
  }
}

ms_reference_loglik <- function(theta, data, ...) {
  ms_reference(data, ...)(theta)
}

# The same with death allowed at unobserved visits (a wrong model, used to show
# that the reference distinguishes the two).
ms_reference_loglik_death_allowed <- function(theta, data) {
  y <- as.integer(as.character(data$y))
  ranges <- list(
    lo = ifelse(is.na(y), 1L, y),
    hi = ifelse(is.na(y), 4L, y),
    recorded = !is.na(y)
  )
  ms_reference_loglik(theta, data, ranges = ranges)
}

# Second derivatives by central differences with step h (the Hessian of a
# function of a few coefficients).
ms_numeric_hessian <- function(f, x, h = 1e-4) {
  p <- length(x)
  out <- matrix(0, p, p)
  f0 <- f(x)
  for (a in seq_len(p)) {
    for (b in a:p) {
      ea <- replace(numeric(p), a, h)
      eb <- replace(numeric(p), b, h)
      if (a == b) {
        out[a, a] <- (f(x + ea) - 2 * f0 + f(x - ea)) / h^2
      } else {
        out[a, b] <- out[b, a] <-
          (f(x + ea + eb) - f(x + ea - eb) - f(x - ea + eb) + f(x - ea - eb)) / (4 * h^2)
      }
    }
  }
  out
}

# The log-likelihood an orm fit reports (path fits and plain fits alike).
ms_loglik <- function(fit) {
  -fit$deviance[[2L]] / 2
}

# The fit maximises `reference`: a Newton step built from the numerical
# gradient of the reference and the fit's inverse information moves no
# coefficient by more than 1e-4 of its standard error. (The gradient itself
# is not near 0 on an absolute scale: the fit stops when the remaining gain in
# log-likelihood is below 1e-10, which allows gradients of about 1e-5 for
# well-determined coefficients.)
ms_expect_maximum <- function(fit, reference, info = NULL) {
  gradient <- central_difference_jacobian(reference, stats::coef(fit))[1, ]
  inverse_information <- fit$orig.var %||% stats::vcov(fit, intercepts = "all")
  step <- as.numeric(inverse_information %*% gradient)
  se <- sqrt(diag(inverse_information))
  testthat::expect_lt(max(abs(step / se)), 1e-4, label = info)
}

# The numeric pieces of an orm fit, without the call or the stored data.
ms_numbers <- function(fit) {
  list(
    coef = stats::coef(fit),
    vcov = stats::vcov(fit, intercepts = "all"),
    var = fit$var,
    orig_var = fit$orig.var
  )
}
