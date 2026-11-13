# Shared coefficient-draw helpers for simulation inference.

generate_sop_coefficient_draws <- function(model, n_draws, vcov) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    mostr_abort(
      "Package 'mvtnorm' is required for MVN simulation inference.\n",
      "Install with: install.packages('mvtnorm')"
    )
  }

  beta_hat <- get_coef(model)
  Sigma <- if (!is.null(vcov)) {
    validate_coef_vcov(beta_hat, vcov, arg = "vcov")
  } else {
    validate_coef_vcov(beta_hat, get_vcov_robust(model), arg = "model vcov")
  }
  mvtnorm::rmvnorm(n_draws, mean = beta_hat, sigma = Sigma)
}

apply_sop_simulation_draws <- function(
  n_draws,
  analysis_fn,
  workers,
  packages = c("rms", "VGAM", "stats", "mostr")
) {
  use_parallel <- !is.null(workers) && workers > 1
  draw_ids <- seq_len(n_draws)
  # Draws whose cut-points cross are counted and reported together.
  draw_fn <- function(i) run_draw_recording_crossings(analysis_fn, i)

  if (!use_parallel) {
    results <- lapply(draw_ids, draw_fn)
    stop_crossed_threshold_draws(results, "simulation draws", "draw")
    return(results)
  }
  if (!requireNamespace("furrr", quietly = TRUE)) {
    mostr_abort("Package 'furrr' is required for parallel processing")
  }

  # `plan(new)` returns the user's whole plan, nested levels included, so that
  # is what gets restored; `plan()` alone returns only the outer level.
  old_plan <- future::plan(future.callr::callr, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)

  results <- furrr::future_map(
    draw_ids,
    draw_fn,
    # `draw_fn` reaches everything it uses through its own environment, which
    # is sent with it, so no other object is exported to the workers.
    .options = furrr::furrr_options(
      seed = FALSE,
      globals = FALSE,
      packages = packages
    ),
    .progress = FALSE
  )
  stop_crossed_threshold_draws(results, "simulation draws", "draw")
  results
}
