# The readable reference for state occupancy probabilities (SOPs).
#
# Everything here is plain R loops: one patient, one visit, one previous state
# at a time, written directly from the Markov recursion. It calls no mostr
# function. `predict_transition(row)` is the only input from the model. It takes
# one data row and returns that row's probabilities of each current state, and
# in practice it wraps the fitting package's own predict() (see
# `slow_predictor()`), so these loops share no prediction code with mostr either.
#
# It is deliberately slow and repetitive so that it can be checked by eye against
# the formulas below. Please keep it that way: do not vectorize it, and do not
# call package code from it. Otherwise it stops being an independent check.
#
# Profiles must start outside the absorbing states. The first- and second-order
# paths treat a start in an absorbing state differently, and which behaviour is
# correct is still undecided, so these loops do not take a position on it.

# The row used for a transition at one visit: the patient's own covariates, with
# the time variable, and any visit-specific columns such as a gap, set for that
# visit.
slow_row_at_visit <- function(patient, visit, times, visit_columns = NULL) {
  row <- patient
  row$time <- if (is.factor(patient$time)) {
    factor(times[visit], levels = levels(patient$time))
  } else {
    times[visit]
  }
  for (column in names(visit_columns)) {
    row[[column]] <- visit_columns[visit, column]
  }
  row
}

slow_set_state <- function(row, column, state) {
  row[[column]] <- if (is.factor(row[[column]])) {
    factor(state, levels = levels(row[[column]]))
  } else {
    as.numeric(state)
  }
  row
}

slow_check_start <- function(patient, absorb) {
  if (as.character(patient$yprev) %in% as.character(absorb)) {
    stop("The slow loop needs a profile that starts outside `absorb`.")
  }
}

# First-order recursion. At the first visit the patient's previous state is
# known, so S_1(l) is the model's transition row. At later visits:
#
#   S_t(l) = sum over non-absorbing k of S_(t-1)(k) * P(Y_t = l | Y_(t-1) = k)
#            + S_(t-1)(l) if l is absorbing
#
# Returns a matrix with one row per visit and one column per state.
slow_first_order_sops <- function(
  predict_transition,
  patient,
  times,
  states,
  absorb = NULL,
  visit_columns = NULL
) {
  slow_check_start(patient, absorb)
  occupancy <- matrix(0, nrow = length(times), ncol = length(states))

  occupancy[1, ] <- predict_transition(
    slow_row_at_visit(patient, 1, times, visit_columns)
  )

  for (visit in seq_along(times)[-1]) {
    for (from in seq_along(states)) {
      was_here <- occupancy[visit - 1, from]

      if (states[from] %in% absorb) {
        # Once in an absorbing state, the patient stays there.
        occupancy[visit, from] <- occupancy[visit, from] + was_here
      } else {
        row <- slow_row_at_visit(patient, visit, times, visit_columns)
        row <- slow_set_state(row, "yprev", states[from])
        transition <- predict_transition(row)

        for (to in seq_along(states)) {
          occupancy[visit, to] <- occupancy[visit, to] + was_here * transition[to]
        }
      }
    }
  }

  occupancy
}

# Second-order recursion. The current state now depends on the previous state j
# and the state before it, h, so the loop tracks J_t(j, l), the probability of
# being in j at visit t - 1 and in l at visit t:
#
#   J_t(j, l) = sum over h of J_(t-1)(h, j) * P(Y_t = l | Y_(t-1) = j, Y_(t-2) = h)
#   S_t(l)    = sum over j of J_t(j, l)
#
# A patient in an absorbing state j stays there, whichever h came before.
slow_second_order_sops <- function(
  predict_transition,
  patient,
  times,
  states,
  absorb = NULL
) {
  slow_check_start(patient, absorb)
  n_states <- length(states)
  occupancy <- matrix(0, nrow = length(times), ncol = n_states)

  first <- predict_transition(slow_row_at_visit(patient, 1, times))
  start <- match(as.character(patient$yprev), as.character(states))
  joint <- matrix(0, nrow = n_states, ncol = n_states)
  joint[start, ] <- first
  occupancy[1, ] <- first

  for (visit in seq_along(times)[-1]) {
    updated <- matrix(0, nrow = n_states, ncol = n_states)

    for (previous in seq_len(n_states)) {
      for (older in seq_len(n_states)) {
        was_here <- joint[older, previous]
        if (was_here == 0) {
          next
        }

        if (states[previous] %in% absorb) {
          updated[previous, previous] <- updated[previous, previous] + was_here
        } else {
          row <- slow_row_at_visit(patient, visit, times)
          row <- slow_set_state(row, "ypprev", states[older])
          row <- slow_set_state(row, "yprev", states[previous])
          transition <- predict_transition(row)

          for (to in seq_len(n_states)) {
            updated[previous, to] <- updated[previous, to] +
              was_here * transition[to]
          }
        }
      }
    }

    joint <- updated
    occupancy[visit, ] <- colSums(joint)
  }

  occupancy
}

# A row-at-a-time predictor from the fitting package's own predict().
slow_predictor <- function(fit) {
  if (inherits(fit, "robcov_vglm")) {
    fit <- fit$vglm_fit
  }
  if (inherits(fit, "vglm")) {
    return(function(row) {
      as.vector(VGAM::predict(fit, row, type = "response"))
    })
  }
  function(row) {
    as.vector(stats::predict(fit, newdata = row, type = "fitted.ind"))
  }
}
