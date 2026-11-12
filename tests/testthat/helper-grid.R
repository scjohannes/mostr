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
# `id`: how patient IDs are stored, as numbers, character strings ("p1", ...) or
# a factor of those strings, whose sorted levels differ from the order in which
# patients appear.
grid_cell <- function(
  engine,
  order = 1L,
  time = "linear",
  weighted = FALSE,
  odds = "full",
  previous = "factor",
  visit = "regular",
  id = "numeric"
) {
  list(
    engine = engine,
    order = order,
    time = time,
    weighted = weighted,
    odds = odds,
    previous = previous,
    visit = visit,
    id = id
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
    paste0(cell$id, " IDs"),
    sep = ", "
  )
}

# TRUE when tests run on CRAN, by the rule `skip_on_cran()` uses: not
# interactive, and NOT_CRAN not set to "true" (devtools and CI set it).
grid_on_cran <- function() {
  !interactive() && !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
}

# Every set-up for the consistency checks: the core crossing plus one cell per
# engine for each extra axis. On CRAN only one cell per engine and Markov order
# runs, so that the checks there stay cheap.
grid_consistency_cells <- function() {
  if (grid_on_cran()) {
    return(list(
      grid_cell("orm"),
      grid_cell("orm", order = 2L, time = "factor"),
      grid_cell("vglm", time = "rcs"),
      grid_cell("vglm", order = 2L, odds = "partial")
    ))
  }
  extras <- list()
  for (engine in c("orm", "vglm")) {
    extras <- c(
      extras,
      list(
        grid_cell(engine, weighted = TRUE),
        grid_cell(engine, visit = "gap"),
        grid_cell(engine, time = "basis"),
        grid_cell(engine, previous = "numeric"),
        grid_cell(engine, id = "character"),
        grid_cell(engine, time = "factor", id = "factor")
      )
    )
  }
  c(grid_core_cells(), extras)
}

# Analytical (delta) inference supports first-order, full proportional-odds
# fits. Unconditional variance also needs a fit without case weights.
grid_has_delta <- function(cell) {
  cell$order == 1L && cell$odds == "full"
}

grid_has_unconditional <- function(cell) {
  grid_has_delta(cell) && !cell$weighted
}

# Elapsed times of the requested visits, unevenly spaced, for real-time
# summaries. The span from baseline time 0 is 14.
grid_time_map <- function(cell, times = grid_times(cell)) {
  stats::setNames(c(2, 5, 9, 14)[seq_along(times)], as.character(times))
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
  if (cell$id != "numeric") {
    data$id <- paste0("p", data$id)
  }
  if (cell$id == "factor") {
    data$id <- factor(data$id)
  }
  data
}

grid_fit <- function(cell) {
  grid_refit(cell, "unchanged data", identity)
}

# Fits `cell` to `change(grid_cell_data(cell))` and caches the fit under
# `change_name`, so that a changed data set is refitted once however many tests
# use it.
grid_refit <- function(cell, change_name, change) {
  key <- paste(grid_label(cell), change_name, sep = " | ")
  if (is.null(grid_cache$fits[[key]])) {
    if (is.null(grid_cache$fits)) {
      grid_cache$fits <- list()
    }
    grid_cache$fits[[key]] <- grid_fit_to(cell, change(grid_cell_data(cell)))
  }
  grid_cache$fits[[key]]
}

# A fixture that several tests build with the same arguments, made once by
# `make()` and kept under `key` for the rest of the test run. Test files share
# it, so the result must not depend on anything a test changes.
grid_fixture <- function(key, make) {
  if (is.null(grid_cache$fixtures)) {
    grid_cache$fixtures <- list()
  }
  if (is.null(grid_cache$fixtures[[key]])) {
    grid_cache$fixtures[[key]] <- make()
  }
  grid_cache$fixtures[[key]]
}

# `p2_var` names the column with the state two visits back that the fit
# records; by default the column of second-order cells.
grid_fit_to <- function(
  cell,
  data,
  p2_var = if (cell$order == 2L) "ypprev"
) {
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

  suppressWarnings(
    if (cell$engine == "orm") {
      formula <- stats::as.formula(paste("y ~", rhs))
      if (cell$weighted) {
        mostm(
          formula,
          data = data,
          id_var = "id",
          first_followup_time = first_followup,
          p2_var = p2_var,
          weights = w
        )
      } else {
        mostm(
          formula,
          data = data,
          id_var = "id",
          first_followup_time = first_followup,
          p2_var = p2_var
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
          p2_var = p2_var,
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
          p2_var = p2_var,
          engine = "vglm",
          family = family
        )
      }
    }
  )
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

# One string per row from the `keys` columns, for matching rows of two results.
grid_key <- function(x, keys) {
  x <- as.data.frame(x)
  if (length(keys) == 0L) {
    return(rep("", nrow(x)))
  }
  do.call(
    paste,
    c(lapply(unname(as.list(x[keys])), as.character), sep = "|")
  )
}

# The largest absolute difference between column `x_col` of `x` and `y_col`
# of `y`, after matching rows on `keys`. Inf unless the rows match one to one,
# so that a lost or duplicated row fails the check. With `all_rows = FALSE`,
# `y` may have rows that `x` lacks.
grid_gap <- function(x, y, keys, x_col = "estimate", y_col = x_col, all_rows = TRUE) {
  key_x <- grid_key(x, keys)
  key_y <- grid_key(y, keys)
  index <- match(key_x, key_y)
  if (
    nrow(as.data.frame(x)) == 0L ||
      (all_rows && length(key_x) != length(key_y)) ||
      anyDuplicated(key_x) > 0L ||
      anyDuplicated(key_y) > 0L ||
      anyNA(index)
  ) {
    return(Inf)
  }
  gaps <- abs(x[[x_col]] - y[[y_col]][index])
  if (anyNA(gaps)) {
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

# ---- Shared by the consistency checks ----

grid_treatments <- list(tx = c(0, 1))

# The public functions, called on a cell's fit with that cell's arguments.
# Results are cached, because several checks start from the same result.
grid_call <- function(fun, cell, ..., args = grid_public_args(cell)) {
  # The order of named arguments does not matter for the key.
  dots <- list(...)
  if (length(dots) > 0L) {
    dots <- dots[order(names(dots))]
  }
  key <- paste(grid_label(cell), rlang::hash(list(fun, dots, args)))
  if (is.null(grid_cache$results[[key]])) {
    if (is.null(grid_cache$results)) {
      grid_cache$results <- list()
    }
    grid_cache$results[[key]] <- do.call(fun, c(list(grid_fit(cell), ...), args))
  }
  grid_cache$results[[key]]
}

# Analytical inference on a result, cached in the same way.
grid_delta <- function(result, vcov) {
  key <- paste("delta", vcov, rlang::hash(result))
  if (is.null(grid_cache$results[[key]])) {
    grid_cache$results[[key]] <- inferences(result, method = "delta", vcov = vcov)
  }
  grid_cache$results[[key]]
}

grid_vcov_types <- function(cell) {
  c("conditional", if (grid_has_unconditional(cell)) "unconditional")
}

# The setting whose turn it is in the `index`-th set-up.
grid_turn <- function(choices, index) {
  choices[[(index - 1L) %% length(choices) + 1L]]
}

# Real-time time benefit is slow to compute, so every other set-up checks it,
# and the others check real-time time in state instead.
grid_real_time_benefit_turn <- function(index) {
  index %% 2L == 1L
}

# Real-time totals for each profile are slow to compute, so each set-up
# computes them with one choice of baseline time, shared by the checks.
grid_profile_real_time <- function(cell, index) {
  list(
    time_map = grid_time_map(cell),
    baseline_time = grid_turn(list(0, NULL), index)
  )
}

# The proportion of stored profiles starting in each state: the averaged
# occupancy at baseline time, the same in every treatment scenario.
grid_start_distribution <- function(cell) {
  first <- grid_call(sops, cell)
  first <- first[as.character(first$time) == as.character(grid_times(cell)[1]) &
    as.character(first$state) == "1", , drop = FALSE]
  start <- as.character(first$yprev)
  vapply(
    as.character(seq_len(grid_n_states)),
    function(state) mean(start == state),
    numeric(1)
  )
}

grid_trapezoid <- function(x, y) {
  sum(diff(x) * (utils::head(y, -1L) + utils::tail(y, -1L)) / 2)
}

# Largest fall between consecutive times of the occupancy of `state`, within
# each group of `groups`. `times` gives the order of the time values.
grid_largest_fall <- function(x, state, groups, times) {
  x <- as.data.frame(x)
  x <- x[as.character(x$state) == as.character(state), , drop = FALSE]
  position <- match(as.character(x$time), as.character(times))
  x <- x[order(grid_key(x, groups), position), , drop = FALSE]
  falls <- unlist(lapply(
    split(x$estimate, grid_key(x, groups)),
    function(values) -diff(values)
  ))
  max(c(0, falls))
}
