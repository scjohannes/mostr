# Shared fixture for `comparisons()` tests: a first-order full
# proportional-odds `vglm` fit with three states, state 3 absorbing.
make_comparisons_case <- function() {
  skip_if_not_installed("VGAM")
  withr::local_seed(4417)
  data <- make_markov_trajectories(
    n_patients = 60,
    follow_up_time = 5,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = 1:2,
    absorbing_state = 3,
    treatment_effect = 0.4,
    seed = 4417
  ) |>
    prepare_markov_data(absorbing_state = 3)
  data$marker <- data$id %% 2
  model <- vglm_markov(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id"
  )
  list(
    model = model,
    data = data,
    profiles = data[data$time == 1, , drop = FALSE],
    time_map = c("1" = 1, "2" = 2, "3" = 4, "4" = 7)
  )
}

# SOPs for every profile with `var` set to `value`, in `sops()` row order:
# profile fastest, then visit, then state.
sops_with_value <- function(case, var, value, times = 1:4, ...) {
  profiles <- case$profiles
  profiles[[var]] <- value
  sops(case$model, newdata = profiles, times = times, absorb = 3, ...)
}
