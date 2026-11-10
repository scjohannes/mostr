# Checks that changes to the fitting data which should not matter leave the
# results alone: the order of the rows, how patient IDs are written, how
# treatment is coded, the units of a covariate, columns the model does not use,
# and case weights against duplicated patients. Each check refits one base
# set-up of the grid (first order, full proportional odds) and compares the
# refit with the base fit. test-grid-data-changes.R runs them on half of the
# twelve base set-ups, and test-slow-grid-data-changes.R on the other half.

change_treatments <- list(tx = c(0, 1))

# Engine x time x case weights. `half = "default"` keeps engine x time crossed
# and alternates the weights, so each engine is checked with and without them;
# `half = "slow"` is the other six.
change_base_cells <- function(half = c("default", "slow")) {
  half <- match.arg(half)
  cells <- list()
  weighted <- half == "slow"
  for (engine in c("orm", "vglm")) {
    for (time in c("linear", "rcs", "factor")) {
      cells[[length(cells) + 1L]] <- grid_cell(engine, time = time, weighted = weighted)
      weighted <- !weighted
    }
  }
  cells
}

# Two fits of the same likelihood converge to the same maximum only up to the
# fitting algorithms' convergence criteria. On these fits the largest
# difference measured in SOPs, time totals and standard errors was 3e-12; the
# tolerance leaves a wide margin above that and still catches any real change.
refit_tolerance <- 1e-8

# The same rows in the same order: the fits are the same computation.
same_rows_tolerance <- 1e-10

# Results of a fit that the changes should leave alone, with the keys that
# identify each row. `id_back` maps the refit's patient IDs back to the base
# data's IDs. Weighted fits refuse unconditional variance, so their averages
# get conditional standard errors only.
change_results <- function(cell, fit, id_back = identity) {
  args <- grid_public_args(cell)
  time_map <- grid_time_map(cell)
  call <- function(fun, ...) do.call(fun, c(list(fit, ...), args))
  each <- as.data.frame(call(sops))
  each$id <- id_back(each$id)
  average <- call(avg_sops, variables = change_treatments)
  vcov <- if (cell$weighted) "conditional" else "unconditional"
  list(
    each = list(each, c("id", "time", "state")),
    average = list(average, c("tx", "time", "state")),
    real_time = list(
      call(avg_state_time, variables = change_treatments, time_map = time_map),
      c("tx", "state_set")
    ),
    time_benefit = list(
      call(
        avg_comparisons,
        variables = change_treatments,
        estimand = "time_benefit",
        time_map = time_map
      ),
      "contrast"
    ),
    conditional_se = list(
      inferences(average, method = "delta", vcov = "conditional"),
      c("tx", "time", "state"),
      "std.error"
    ),
    average_se = list(
      inferences(average, method = "delta", vcov = vcov),
      c("tx", "time", "state"),
      "std.error"
    )
  )
}

# The largest difference between two sets of results, for each result.
change_gaps <- function(base, changed) {
  vapply(names(base), function(name) {
    column <- if (length(base[[name]]) > 2L) base[[name]][[3]] else "estimate"
    grid_gap(changed[[name]][[1]], base[[name]][[1]], base[[name]][[2]], column)
  }, numeric(1))
}

expect_same_results <- function(base, changed, tolerance, label) {
  gaps <- change_gaps(base, changed)
  for (name in names(gaps)) {
    expect_lt(gaps[[name]], tolerance, label = paste(label, name))
  }
}

shuffle_rows <- function(data) {
  withr::with_seed(5, data[sample(nrow(data)), , drop = FALSE])
}

# Patient IDs written another way, and how to map them back. Each keeps the
# order in which patients first appear.
id_changes <- list(
  "relabelled numeric IDs" = list(
    change = function(data) {
      data$id <- 10000L - data$id
      data
    },
    back = function(id) 10000L - id
  ),
  "character IDs" = list(
    change = function(data) {
      data$id <- paste0("p", data$id)
      data
    },
    back = function(id) as.integer(sub("^p", "", as.character(id)))
  ),
  "factor IDs" = list(
    change = function(data) {
      data$id <- factor(paste0("p", data$id))
      data
    },
    back = function(id) as.integer(sub("^p", "", as.character(id)))
  )
)

# tx and 1 - tx with the compared values swapped, and tx as a factor with
# labels: two values saturate the treatment term, so every coding is the same
# model. `back` maps the recoded values back to 0 and 1.
treatment_codings <- list(
  "reversed treatment" = list(
    change = function(data) {
      data$tx <- 1 - data$tx
      data
    },
    variables = list(tx = c(1, 0)),
    back = function(tx) 1 - tx
  ),
  "labelled treatment" = list(
    change = function(data) {
      data$tx <- factor(ifelse(data$tx == 1, "treated", "control"))
      data
    },
    variables = list(tx = c("control", "treated")),
    back = function(tx) as.numeric(as.character(tx) == "treated")
  )
)

age_in_months <- function(data) {
  data$age <- data$age * 12
  data
}

add_unused_column <- function(data) {
  data$noise <- withr::with_seed(7, stats::rnorm(nrow(data)))
  data$noise[data$id %% 10 == 0] <- NA
  data
}

# Patient i, with weight w_i, becomes w_i patients with new IDs.
duplicate_patients <- function(data) {
  rows <- rep(seq_len(nrow(data)), data$w)
  copy <- stats::ave(seq_along(rows), rows, FUN = seq_along)
  out <- data[rows, , drop = FALSE]
  out$id <- out$id * 10L + copy
  out <- out[order(out$id, out$time), , drop = FALSE]
  rownames(out) <- NULL
  out
}

constant_weight <- function(data) {
  data$w <- 2
  data
}

check_row_order <- function(cell) {
  base <- change_results(cell, grid_fit(cell))
  changed <- change_results(cell, grid_refit(cell, "shuffled rows", shuffle_rows))
  expect_same_results(base, changed, refit_tolerance, paste(grid_label(cell), "shuffled rows"))
}

check_id_writing <- function(cell) {
  base <- change_results(cell, grid_fit(cell))
  for (name in names(id_changes)) {
    fit <- grid_refit(cell, name, id_changes[[name]]$change)
    changed <- change_results(cell, fit, id_changes[[name]]$back)
    expect_same_results(base, changed, same_rows_tolerance, paste(grid_label(cell), name))
  }
}

check_treatment_coding <- function(cell) {
  args <- grid_public_args(cell)
  keys <- c("tx", "time", "state")
  summaries <- function(fit, variables, back) {
    average <- do.call(avg_sops, c(list(fit, variables = variables), args))
    with_se <- as.data.frame(inferences(average, method = "delta", vcov = "conditional"))
    with_se$tx <- back(with_se$tx)
    benefit <- do.call(
      avg_comparisons,
      c(list(fit, variables = variables, estimand = "time_benefit"), args)
    )
    list(average = with_se, benefit = benefit$estimate)
  }
  base <- summaries(grid_fit(cell), change_treatments, identity)
  for (name in names(treatment_codings)) {
    coding <- treatment_codings[[name]]
    changed <- summaries(grid_refit(cell, name, coding$change), coding$variables, coding$back)
    label <- paste(grid_label(cell), name)
    expect_lt(grid_gap(changed$average, base$average, keys), refit_tolerance, label = label)
    expect_lt(
      grid_gap(changed$average, base$average, keys, "std.error"),
      refit_tolerance,
      label = paste(label, "standard error")
    )
    expect_lt(
      abs(changed$benefit - base$benefit),
      refit_tolerance,
      label = paste(label, "time benefit")
    )
  }
}

check_covariate_units <- function(cell) {
  base <- change_results(cell, grid_fit(cell))
  changed <- change_results(cell, grid_refit(cell, "age in months", age_in_months))
  expect_same_results(base, changed, refit_tolerance, paste(grid_label(cell), "age in months"))
}

check_unused_column <- function(cell) {
  label <- paste(grid_label(cell), "unused column")
  base <- change_results(cell, grid_fit(cell))
  fit <- grid_refit(cell, "unused column with missing values", add_unused_column)
  expect_same_results(base, change_results(cell, fit), same_rows_tolerance, label)

  # Per-profile real-time totals keep every patient, so their mean is the
  # average total.
  args <- grid_public_args(cell)
  time_map <- grid_time_map(cell)
  each <- as.data.frame(do.call(
    state_time,
    c(list(fit, variables = change_treatments, time_map = time_map), args)
  ))
  expect_false(anyNA(each$estimate), label = label)
  means <- stats::aggregate(each["estimate"], each[c("tx", "state_set")], mean)
  expect_lt(
    grid_gap(means, base$real_time[[1]], c("tx", "state_set")),
    same_rows_tolerance,
    label = paste(label, "mean of state_time()")
  )
  totals <- time_in_state(
    do.call(sops, c(list(fit), args)),
    target_states = "1",
    time_map = time_map
  )
  expect_identical(nrow(totals), 120L, label = label)

  # The same column, with a missing value, in supplied profiles.
  profiles <- grid_profiles(cell)
  noisy <- profiles
  noisy$noise <- c(NA, 1, 2)
  average <- function(newdata) {
    do.call(
      avg_sops,
      c(list(grid_fit(cell), newdata = newdata, variables = change_treatments), args)
    )
  }
  expect_lt(
    grid_gap(average(noisy), average(profiles), c("tx", "time", "state")),
    same_rows_tolerance,
    label = paste(label, "in newdata")
  )
}

# Weight w on every row of a patient gives the same likelihood as w copies of
# the patient with new IDs. So the coefficients, and SOPs for fixed profiles,
# agree. Averages over the stored patients do not: the weighted fit counts each
# patient once, the duplicated fit w times. The w-weighted mean of the weighted
# fit's per-patient SOPs equals the duplicated fit's average instead. Robust
# standard errors are not compared: they treat one patient with weight w and w
# independent patients differently by design.
check_weights_against_duplicates <- function(cell) {
  label <- grid_label(cell)
  unweighted <- cell
  unweighted$weighted <- FALSE
  weighted_fit <- grid_fit(cell)
  duplicated_fit <- grid_refit(unweighted, "patients duplicated by weight", duplicate_patients)
  args <- grid_public_args(cell)

  profiles <- grid_profiles(cell)
  expect_lt(
    grid_gap(
      do.call(sops, c(list(weighted_fit, newdata = profiles), args)),
      do.call(sops, c(list(duplicated_fit, newdata = profiles), args)),
      c("rowid", "time", "state")
    ),
    refit_tolerance,
    label = paste(label, "fixed profiles")
  )

  each <- as.data.frame(do.call(
    sops,
    c(list(weighted_fit, variables = change_treatments), args)
  ))
  data <- grid_cell_data(cell)
  weight <- data$w[match(each$id, data$id)]
  weighted_mean <- stats::aggregate(
    data.frame(estimate = each$estimate * weight, weight = weight),
    each[c("tx", "time", "state")],
    sum
  )
  weighted_mean$estimate <- weighted_mean$estimate / weighted_mean$weight
  expect_lt(
    grid_gap(
      weighted_mean,
      do.call(avg_sops, c(list(duplicated_fit, variables = change_treatments), args)),
      c("tx", "time", "state")
    ),
    refit_tolerance,
    label = paste(label, "weighted mean of stored patients")
  )
}

# A weight c on every row multiplies every patient's score by c and the
# information by c, so the robust (sandwich) covariance is unchanged, and so
# are the conditional standard errors.
check_constant_weight <- function(cell) {
  weighted <- cell
  weighted$weighted <- TRUE
  unweighted <- cell
  unweighted$weighted <- FALSE
  base <- change_results(unweighted, grid_fit(unweighted))
  changed <- change_results(weighted, grid_refit(weighted, "constant weight 2", constant_weight))
  base$average_se <- NULL
  changed$average_se <- NULL
  expect_same_results(base, changed, refit_tolerance, paste(grid_label(cell), "constant weight"))
}
