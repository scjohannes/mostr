# Generate test trajectories with a proportional-odds Markov model.
make_markov_trajectories <- function(
  n_patients = 50,
  follow_up_time = 20,
  n_states = 6,
  thresholds = seq(-3, 3, length.out = n_states - 1),
  allowed_start_state = seq_len(n_states - 1),
  absorbing_state = n_states,
  treatment_prob = 0.5,
  treatment_effect = 0,
  seed = 123
) {
  set.seed(seed)
  baseline <- data.frame(
    id = seq_len(n_patients),
    yprev = sample(allowed_start_state, n_patients, replace = TRUE),
    tx = rbinom(n_patients, 1, treatment_prob)
  )
  sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = follow_up_time,
    intercepts = thresholds,
    absorbing_states = absorbing_state,
    covariate_names = "tx",
    parameter = treatment_effect,
    extra_params = c(center = (n_states + 1) / 2),
    lp_function = function(yprev, t, tx, parameter, extra_params) {
      0.4 *
        (as.numeric(as.character(yprev)) - extra_params["center"]) -
        0.03 * t +
        parameter * tx
    }
  )
}

# Helper functions for simulating data and fitting models in tests

#' Generate standard test data using Markov transition simulation
#'
#' @param n_patients Number of patients
#' @param follow_up_time Follow-up time
#' @param seed Random seed
#' @param treatment_effect Treatment effect (treatment_effect)
#'
#' @return A prepared markov data frame
make_test_data <- function(
  n_patients = 50,
  follow_up_time = 20,
  seed = 123,
  treatment_effect = 0
) {
  raw_data <- make_markov_trajectories(
    n_patients = n_patients,
    follow_up_time = follow_up_time,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = seed,
    treatment_effect = treatment_effect
  )

  data <- prepare_markov_data(raw_data)

  # Generate manual spline basis
  # We use rcs() to get the basis, then extract columns
  time_spl_m <- rms::rcs(data$time, 3)
  knots <- attr(time_spl_m, "parms")
  data$time_lin <- as.vector(time_spl_m[, 1])
  data$time_nlin_1 <- as.vector(time_spl_m[, 2])

  return(data)
}

make_time_covariates <- function(data, time_col = "time", ...) {
  cols <- c(time_col, ...)
  out <- unique(data[cols])
  out <- out[order(out[[time_col]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Fit a standard VGLM model for testing
#'
#' @param data Data frame from make_test_data
#' @param robust Logical, whether to return a robust covariance model
#' @param cluster Cluster variable for robust covariance
#'
#' @return A fitted vglm or robcov_vglm object
make_test_model <- function(data, robust = FALSE, cluster = NULL) {
  if (robust && !is.null(cluster) && !identical(cluster, data$id)) {
    stop("`cluster` must match `data$id` for wrapper-based test models.")
  }

  suppressWarnings(vglm_markov(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = if (robust) "id" else NULL
  ))
}

make_score_bootstrap_case <- function(
  seed,
  n_patients = 60,
  follow_up_time = 10
) {
  data <- make_test_data(
    n_patients = n_patients,
    seed = seed,
    follow_up_time = follow_up_time
  )

  list(
    data = data,
    baseline = data[data$time == 1, , drop = FALSE],
    model = vglm_markov(
      ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id"
    ),
    times = seq_len(follow_up_time),
    y_levels = 1:6,
    absorb = 6
  )
}
