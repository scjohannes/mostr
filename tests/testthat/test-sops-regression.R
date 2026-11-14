describe("SOP regression baselines", {
  add_patient_age <- function(data) {
    data$age <- NULL

    age_lookup <- data.frame(
      id = sort(unique(data$id)),
      age = 50 + ((sort(unique(data$id)) * 7) %% 23)
    )

    merge(data, age_lookup, by = "id", all.x = TRUE, sort = FALSE)
  }

  fit_inline_spline_model <- function(data, absorb) {
    suppressWarnings(
      mostm(
        ordered(y) ~ rms::rcs(time, 4) * tx + yprev + age,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        absorb = absorb,
        engine = "vglm"
      )
    )
  }

  sops_signature <- function(
    sops_array,
    patient_idx = 1:3,
    time_idx = c(1, 2, 4, 6, 10),
    state_idx = c(1, 3, dim(sops_array)[3])
  ) {
    patient_idx <- patient_idx[patient_idx <= dim(sops_array)[1]]
    time_idx <- time_idx[time_idx <= dim(sops_array)[2]]
    state_idx <- state_idx[state_idx <= dim(sops_array)[3]]

    out <- expand.grid(
      patient = patient_idx,
      time = time_idx,
      state = state_idx,
      KEEP.OUT.ATTRS = FALSE
    )
    out$prob <- round(
      mapply(
        function(i, t, s) sops_array[i, t, s],
        out$patient,
        out$time,
        out$state
      ),
      8
    )

    out[order(out$patient, out$time, out$state), , drop = FALSE]
  }

  avg_signature <- function(
    result,
    time_idx = c(1, 2, 4, 6, 10),
    state_idx = c(1, 3, max(as.integer(as.character(result$state))))
  ) {
    out <- subset(result, time %in% time_idx & state %in% state_idx)
    out <- out[
      order(out$tx, out$time, out$state),
      c("tx", "time", "state", "estimate")
    ]
    out$state <- as.integer(as.character(out$state))
    out$estimate <- round(out$estimate, 8)
    rownames(out) <- NULL
    out
  }

  draw_signature <- function(
    draws,
    draw_id = 1L,
    time_idx = c(1, 2, 4, 6, 8),
    state_idx = c(1, 3, max(as.integer(as.character(draws$state))))
  ) {
    out <- subset(
      draws,
      draw_id == !!draw_id & time %in% time_idx & state %in% state_idx
    )
    out <- out[
      order(out$tx, out$time, out$state),
      c("draw_id", "tx", "time", "state", "draw")
    ]
    out$state <- as.integer(as.character(out$state))
    out$draw <- round(out$draw, 8)
    rownames(out) <- NULL
    out
  }

  build_generic_markov_case <- function() {
    raw_data <- make_markov_trajectories(
      n_patients = 80,
      follow_up_time = 10,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 1001,
      treatment_effect = 0.2
    )

    data <- mostr::prepare_markov_data(raw_data)
    data <- add_patient_age(data)

    list(
      data = data,
      baseline = subset(data, time == 1),
      y_levels = 1:6,
      absorb = "6"
    )
  }

  build_markov_case <- function() {
    raw_data <- mostr::sim_markov(
      "violet",
      n_patients = 80,
      follow_up_time = 10,
      treatment_effect = log(0.85),
      seed = 1004
    )

    # Death, state 4, is absorbing: its rows as the previous state are removed.
    data <- mostr::prepare_markov_data(raw_data, absorbing_state = 4)
    data <- add_patient_age(data)

    list(
      data = data,
      baseline = subset(data, time == 1),
      y_levels = 1:4,
      absorb = "4"
    )
  }

  test_that("generic Markov SOP signatures remain stable", {
    skip_on_cran() # The pass depends on the seed.
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    case <- build_generic_markov_case()
    model <- suppressWarnings(fit_inline_spline_model(case$data, case$absorb))

    sops_result <- mostr::soprob_markov(
      model = model,
      newdata = case$baseline,
      times = 1:10
    )

    avg_result <- mostr::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:10
    )

    withr::local_seed(5001)
    sim_result <- mostr::inferences(
      avg_result,
      method = "mvn",
      n_draws = 2,
      return_draws = TRUE
    )

    expect_snapshot_value(
      list(
        sops = sops_signature(sops_result),
        avg = avg_signature(avg_result),
        sim_draw = draw_signature(mostr::get_draws(sim_result))
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("markov SOP signatures remain stable", {
    skip_on_cran() # The pass depends on the seed.
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_markov_case()
    # The simulated VIOLET states move by at most one level per day, so this
    # proportional-odds fit drifts towards infinite coefficients and has no
    # usable patient-clustered covariance. The point-estimate signatures need
    # only the coefficients, so the model is fitted as bootstrap refits are.
    model <- mostr:::markov_with_coefficients_only_refit(
      suppressWarnings(fit_inline_spline_model(case$data, case$absorb))
    )

    sops_result <- mostr::soprob_markov(
      model = model,
      newdata = case$baseline,
      times = 1:10
    )

    avg_result <- mostr::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:10
    )

    expect_snapshot_value(
      list(
        sops = sops_signature(sops_result),
        avg = avg_signature(avg_result)
      ),
      style = "json2",
      cran = TRUE
    )
  })
})
