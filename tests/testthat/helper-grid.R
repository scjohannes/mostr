# A shared grid of model set-ups, so that one set of checks can run on every
# combination rather than on one fit. Fits are cached, because several test files
# use the same cell and most of the cost is fitting.

grid_cache <- new.env(parent = emptyenv())

grid_n_states <- 4L

# One simulated data set: 120 patients, 6 visits, 4 states, state 4 absorbing.
# Unequal `marker` groups make a wrong average easy to see.
grid_data <- function() {
  if (is.null(grid_cache$data)) {
    grid_cache$data <- withr::with_seed(2024, {
      data <- make_markov_trajectories(
        n_patients = 120,
        follow_up_time = 6,
        n_states = grid_n_states,
        thresholds = c(-1.5, 0, 1.5),
        allowed_start_state = 1:3,
        absorbing_state = grid_n_states,
        treatment_effect = 0.5,
        seed = 2024
      ) |>
        prepare_markov_data(absorbing_state = grid_n_states)
      data <- data[order(data$id, data$time), , drop = FALSE]
      rownames(data) <- NULL

      # Before the first visit only one state is known, so the state two visits
      # back repeats the previous state there.
      data$ypprev <- stats::ave(
        as.integer(as.character(data$yprev)),
        data$id,
        FUN = function(x) c(x[1L], utils::head(x, -1L))
      )
      data$ypprev <- factor(data$ypprev, levels = levels(data$yprev))

      ids <- unique(data$id)
      data$age <- round(stats::rnorm(length(ids), 60, 10))[match(data$id, ids)]
      data$w <- sample(1:3, length(ids), TRUE)[match(data$id, ids)]
      data$marker <- as.integer(data$id %% 3 == 0)
      data$gap <- 1 + (data$time %% 3)
      data$gap[data$time == 1] <- 1

      basis <- rms::rcs(data$time, 3)
      data$time_lin <- as.vector(basis[, 1])
      data$time_nlin <- as.vector(basis[, 2])
      data
    })
  }
  grid_cache$data
}

# Set-ups as a list of cells. `previous`: how the previous state is coded.
# `time`: linear, spline, factor, or an explicit spline basis given to prediction
# through `time_covariates`. `visit`: "gap" adds a gap variable and uneven visits.
grid_cell <- function(
  engine,
  order = 1L,
  time = "linear",
  weighted = FALSE,
  odds = "full",
  previous = "factor",
  visit = "regular"
) {
  list(
    engine = engine,
    order = order,
    time = time,
    weighted = weighted,
    odds = odds,
    previous = previous,
    visit = visit
  )
}

grid_label <- function(cell) {
  paste(
    cell$engine,
    paste0("order ", cell$order),
    paste0(cell$time, " time"),
    if (cell$weighted) "weighted" else "unweighted",
    paste0(cell$odds, " odds"),
    paste0(cell$previous, " previous state"),
    paste0(cell$visit, " visits"),
    sep = ", "
  )
}

# Engine x order x time, plus partial odds for vglm: the three axes that decide
# which code path runs, so they stay fully crossed.
grid_core_cells <- function(times = c("linear", "rcs", "factor")) {
  cells <- list()
  for (engine in c("orm", "vglm")) {
    for (order in 1:2) {
      for (time in times) {
        cells[[length(cells) + 1L]] <- grid_cell(engine, order, time)
      }
    }
  }
  for (order in 1:2) {
    for (time in intersect(times, c("linear", "factor"))) {
      cells[[length(cells) + 1L]] <- grid_cell(
        "vglm",
        order,
        time,
        odds = "partial"
      )
    }
  }
  cells
}

grid_cell_data <- function(cell) {
  data <- grid_data()
  if (cell$time == "factor") {
    data$time <- factor(data$time, levels = sort(unique(data$time)))
  }
  if (cell$previous == "numeric") {
    data$yprev <- as.numeric(as.character(data$yprev))
  }
  data
}

grid_fit <- function(cell) {
  key <- grid_label(cell)
  if (!is.null(grid_cache$fits[[key]])) {
    return(grid_cache$fits[[key]])
  }

  data <- grid_cell_data(cell)
  time_term <- switch(
    cell$time,
    rcs = "rms::rcs(time, 3)",
    basis = "time_lin + time_nlin",
    "time"
  )
  rhs <- paste(
    c(
      time_term,
      "tx",
      "age",
      "yprev",
      if (cell$order == 2L) "ypprev",
      if (cell$visit == "gap") "gap"
    ),
    collapse = " + "
  )
  first_followup <- if (cell$time == "factor") "1" else 1
  family <- VGAM::cumulative(
    reverse = TRUE,
    parallel = if (cell$odds == "partial") FALSE ~ tx else TRUE
  )

  fit <- suppressWarnings(
    if (cell$engine == "orm") {
      formula <- stats::as.formula(paste("y ~", rhs))
      if (cell$weighted) {
        mostm(
          formula,
          data = data,
          id_var = "id",
          first_followup_time = first_followup,
          weights = w
        )
      } else {
        mostm(
          formula,
          data = data,
          id_var = "id",
          first_followup_time = first_followup
        )
      }
    } else {
      formula <- stats::as.formula(paste("ordered(y) ~", rhs))
      if (cell$weighted) {
        mostm(
          formula,
          data = data,
          id_var = "id",
          first_followup_time = first_followup,
          engine = "vglm",
          family = family,
          weights = w
        )
      } else {
        mostm(
          formula,
          data = data,
          id_var = "id",
          first_followup_time = first_followup,
          engine = "vglm",
          family = family
        )
      }
    }
  )

  if (is.null(grid_cache$fits)) {
    grid_cache$fits <- list()
  }
  grid_cache$fits[[key]] <- fit
  fit
}

grid_times <- function(cell) {
  if (cell$visit == "gap") {
    return(c(1, 2, 4, 5))
  }
  if (cell$time == "factor") {
    return(as.character(1:4))
  }
  1:4
}

# Columns that change by visit, as the slow loop needs them. A gap is the time
# since the previous requested visit, and the first visit's gap is its own time.
grid_visit_columns <- function(cell, times = grid_times(cell)) {
  if (cell$visit == "gap") {
    return(data.frame(gap = c(times[1], diff(times))))
  }
  if (cell$time == "basis") {
    data <- grid_data()
    rows <- data[match(times, data$time), c("time_lin", "time_nlin")]
    rownames(rows) <- NULL
    return(rows)
  }
  NULL
}

# Arguments that select the same set-up in the public functions.
grid_public_args <- function(cell, absorb = grid_n_states, times = grid_times(cell)) {
  args <- list(times = times, absorb = absorb)
  if (cell$order == 2L) {
    args$p2_var <- "ypprev"
  }
  if (cell$visit == "gap") {
    args$gap_var <- "gap"
  }
  if (cell$time == "basis") {
    table <- grid_visit_columns(cell, times)
    table$time <- times
    args$time_covariates <- table
  }
  args
}

# One patient per starting state that is outside `absorb`, with the state two
# visits back shifted so that it differs from the previous state.
grid_profiles <- function(cell, absorb = grid_n_states) {
  data <- grid_cell_data(cell)
  start <- as.character(data$time) == "1"
  baseline <- data[start, , drop = FALSE]
  states <- setdiff(1:(grid_n_states - 1L), absorb)
  rows <- match(states, as.integer(as.character(baseline$yprev)))
  profiles <- baseline[rows, , drop = FALSE]
  profiles$ypprev <- factor(
    rev(states),
    levels = levels(data$ypprev)
  )
  rownames(profiles) <- NULL
  profiles
}

# SOPs for each profile from the readable loop, as a list of visit-by-state
# matrices.
grid_slow_sops <- function(cell, fit, profiles, absorb, times = grid_times(cell)) {
  predict_transition <- slow_predictor(fit)
  states <- seq_len(grid_n_states)
  lapply(seq_len(nrow(profiles)), function(i) {
    patient <- profiles[i, , drop = FALSE]
    if (cell$order == 2L) {
      slow_second_order_sops(predict_transition, patient, times, states, absorb)
    } else {
      slow_first_order_sops(
        predict_transition,
        patient,
        times,
        states,
        absorb,
        grid_visit_columns(cell, times)
      )
    }
  })
}

# The largest absolute difference between a tidy result and visit-by-state
# matrices. `which_matrix` gives the matrix for each row of the result.
grid_largest_gap <- function(result, matrices, which_matrix, times) {
  gaps <- vapply(
    seq_len(nrow(result)),
    function(row) {
      reference <- matrices[[which_matrix[row]]]
      visit <- match(as.character(result$time[row]), as.character(times))
      state <- as.integer(as.character(result$state[row]))
      abs(result$estimate[row] - reference[visit, state])
    },
    numeric(1)
  )
  if (length(gaps) == 0L) {
    return(Inf)
  }
  max(gaps)
}

# Tests that are too slow for the default run. Set MOSTR_RUN_SLOW=true to run them.
skip_unless_slow_tier <- function() {
  skip_if_not(
    identical(Sys.getenv("MOSTR_RUN_SLOW"), "true"),
    "Set MOSTR_RUN_SLOW=true to run slow checks."
  )
}
