# Unconditional variance splits exactly into three parts (catalogue R15):
#
#   unconditional = conditional
#                   + cov(F) / n
#                   + (cov(F, C) + cov(C, F)) / n
#
# F has one row per patient and one column per estimate: the patient's own
# prediction (SOP, contrast or time in state), whose average over patients is
# the estimate. C is the same shape: each patient's influence on the
# coefficients, n * B s_i, carried to the estimates by their Jacobian G, so
# C = n * S B G'. Here S holds the patients' summed likelihood scores, one row
# per patient, and B the model-based bread (inverse information). cov() is the
# sample covariance, with divisor n - 1.
#
# Why exact: the code takes cov(F + C) / n. The C part alone is
# G B S'S B G' * n / (n - 1) less a term in the mean score, which is zero at the
# maximum. That equals G V G', the conditional variance, for the default
# patient-cluster sandwich V with type = "HC0" and cadjust = TRUE, whose
# adjustment factor is n / (n - 1). So the identity holds to rounding, and a
# part that is off by 1% shows. The three parts are built here from per-patient
# predictions from the public functions, row scores and bread from their named
# helpers, and the Jacobian of the conditional result, not from the code that
# adds them up.

identity_cells <- function() {
  list(grid_cell("orm"), grid_cell("vglm"))
}

identity_times <- 1:4

identity_state_sets <- list(low = c("1", "2"), dead = "4")

# Summed likelihood scores per patient (rows in starting-profile order) and the
# bread B, so that the coefficient influence of patient i is n * B s_i.
identity_coefficient_parts <- function(fit) {
  if (inherits(fit, "robcov_vglm")) {
    row_scores <- compute_scores_vglm(fit$vglm_fit)
    bread <- compute_observed_bread_vglm(fit$vglm_fit)$bread
  } else {
    row_scores <- compute_scores_orm(fit)
    bread <- orm_model_bread(fit)$bread
  }
  fitted_rows <- attr(fit, "markov_data")
  stopifnot(nrow(fitted_rows) == nrow(row_scores))
  ids <- as.character(attr(fit, "markov_starting_profile_data")$id)
  scores <- rowsum(row_scores, as.character(fitted_rows$id), reorder = FALSE)
  scores <- scores[ids, , drop = FALSE]
  list(scores = scores, bread = bread[colnames(scores), colnames(scores)])
}

# A patient-by-estimate matrix from a per-patient result `x`: column j holds,
# for every patient in `ids`, the row of `x` whose `keys` match row j of
# `estimates`.
identity_patient_matrix <- function(x, estimates, keys, ids) {
  x <- as.data.frame(x)
  x_key <- paste(as.character(x$id), grid_key(x, keys))
  sapply(seq_len(nrow(estimates)), function(j) {
    wanted <- paste(ids, grid_key(estimates[j, , drop = FALSE], keys))
    x$estimate[match(wanted, x_key)]
  })
}

# The right-hand side of the identity, from its parts.
identity_sum <- function(conditional, values, coefficient) {
  n <- nrow(values)
  jacobian <- get_jacobian(conditional)
  influence <- n *
    coefficient$scores %*%
    coefficient$bread %*%
    t(jacobian[, colnames(coefficient$scores), drop = FALSE])
  cross <- stats::cov(values, influence)
  unname(stats::vcov(conditional)) +
    stats::cov(values) / n +
    (cross + t(cross)) / n
}

# Checks the identity for `point`, an average from the default starting
# profiles, given the per-patient predictions `values`.
expect_unconditional_identity <- function(point, values, coefficient, label) {
  expect_false(anyNA(values), label = label)
  expect_equal(colMeans(values), point$estimate, tolerance = 1e-10, info = label)

  conditional <- inferences(point, method = "delta", vcov = "conditional")
  unconditional <- inferences(point, method = "delta", vcov = "unconditional")
  actual <- unname(stats::vcov(unconditional))
  expected <- identity_sum(conditional, values, coefficient)
  scale <- max(abs(actual))
  expect_lt(max(abs(actual - expected)) / scale, 1e-8, label = label)

  # The cross term is far above rounding here (0.3% to 2.4% of the largest
  # entry), so leaving it out would fail.
  n <- nrow(values)
  without_cross <- unname(stats::vcov(conditional)) + stats::cov(values) / n
  expect_gt(max(abs(actual - without_cross)) / scale, 1e-3, label = label)
}

test_that("unconditional variance is conditional variance plus the patient terms", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  for (cell in identity_cells()) {
    fit <- grid_fit(cell)
    label <- grid_label(cell)
    profiles <- attr(fit, "markov_starting_profile_data")
    ids <- as.character(profiles$id)
    coefficient <- identity_coefficient_parts(fit)

    # Each patient's SOPs under each treatment, from the public `sops()`.
    patient_sops <- lapply(c(0, 1), function(value) {
      profiles$tx <- value
      sops(fit, newdata = profiles, times = identity_times, absorb = 4)
    })
    sop_matrix <- function(estimates, keys) {
      identity_patient_matrix(
        rbind(patient_sops[[1]], patient_sops[[2]]),
        estimates,
        keys,
        ids
      )
    }

    # Average SOPs under each treatment.
    point <- avg_sops(
      fit,
      variables = list(tx = c(0, 1)),
      times = identity_times,
      absorb = 4
    )
    expect_unconditional_identity(
      point,
      sop_matrix(as.data.frame(point), c("time", "state", "tx")),
      coefficient,
      paste(label, "avg_sops")
    )

    # The treatment contrast of average SOPs.
    contrast <- avg_comparisons(
      fit,
      variables = list(tx = c(0, 1)),
      times = identity_times,
      absorb = 4
    )
    cells <- as.data.frame(contrast)
    cells$state <- cells$state_set
    treated <- cells
    treated$tx <- 1
    untreated <- cells
    untreated$tx <- 0
    keys <- c("time", "state", "tx")
    expect_unconditional_identity(
      contrast,
      sop_matrix(treated, keys) - sop_matrix(untreated, keys),
      coefficient,
      paste(label, "avg_comparisons")
    )

    # Average time in state on elapsed time from baseline time 0, so that each
    # patient's starting state also enters the totals.
    time_args <- list(
      fit,
      variables = list(tx = c(0, 1)),
      times = identity_times,
      absorb = 4,
      state_sets = identity_state_sets,
      time_map = grid_time_map(cell)
    )
    average_time <- do.call(avg_state_time, time_args)
    patient_time <- do.call(state_time, time_args)
    expect_unconditional_identity(
      average_time,
      identity_patient_matrix(
        patient_time,
        as.data.frame(average_time),
        c("state_set", "tx"),
        ids
      ),
      coefficient,
      paste(label, "avg_state_time")
    )
  }
})
