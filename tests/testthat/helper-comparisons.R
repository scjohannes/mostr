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
  # A group of two patients, often absent from an ordinary bootstrap resample.
  data$rare <- data$id <= 2
  model <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
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

# The same trajectories with visits stored as a factor, so the model has one
# coefficient per visit and `first_followup_time` must be given explicitly.
make_factor_comparisons_case <- function() {
  case <- make_comparisons_case()
  data <- case$data
  data$time <- factor(data$time, levels = sort(unique(data$time)))
  model <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = "1",
    id_var = "id",
    engine = "vglm"
  )
  list(
    model = model,
    data = data,
    profiles = data[data$time == "1", , drop = FALSE],
    time_map = case$time_map
  )
}

# The same trajectories with the state two visits back, `ypprev`, for a
# second-order model. Before the first visit only one state is known, so
# `ypprev` repeats `yprev` there.
make_second_order_comparisons_case <- function() {
  case <- make_comparisons_case()
  data <- case$data[order(case$data$id, case$data$time), , drop = FALSE]
  data$ypprev <- stats::ave(
    as.integer(as.character(data$yprev)),
    data$id,
    FUN = function(x) c(x[1L], utils::head(x, -1L))
  )
  data$ypprev <- factor(data$ypprev, levels = levels(data$yprev))
  model <- mostm(
    ordered(y) ~ time + tx + yprev + ypprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
  )
  list(
    model = model,
    data = data,
    profiles = data[data$time == 1, , drop = FALSE],
    time_map = case$time_map
  )
}
