# Without `newdata`, every bootstrap and FWB draw must start each patient
# from their row at `first_followup_time`, as the point estimate does. A
# draw that picked the earliest row by sorting `time` started from the wrong
# visit when time is character (alphabetical order puts "Week 10" before
# "Week 2") or a factor with a level before `first_followup_time` (issue
# #47). Each fit below is the same model on the same transitions as the
# numeric-time fit, so the same resamples must give the same draws.

build_start_row_cases <- function() {
  data <- make_markov_trajectories(
    n_patients = 60,
    times = 1:3,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = 1:2,
    absorb = 3,
    treatment_effect = 0.4,
    seed = 4417
  ) |>
    prepare_markov_data(absorb = 3)
  data <- data[order(data$id, data$time), , drop = FALSE]
  rownames(data) <- NULL

  fit <- function(data, first_followup_time) {
    mostm(
      y ~ tx + yprev,
      data = data,
      id_var = "id",
      absorb = 3,
      first_followup_time = first_followup_time
    )
  }

  # Character time whose alphabetical order is not the visit order.
  week_labels <- c("Week 2", "Week 6", "Week 10")
  character_data <- data
  character_data$time <- week_labels[character_data$time]

  # A factor with a screening visit "0" before the first follow-up visit.
  # Its rows have no outcome, so the fit drops them, but their predictors are
  # complete and their previous state differs from the one at visit "1": a
  # draw that started there would run without error and give other SOPs.
  screening <- data[data$time == 1, , drop = FALSE]
  screening$time <- 0
  screening$y <- NA
  screening$yprev <- ifelse(screening$yprev == 1, 2, 1)
  screening_data <- rbind(screening, data)
  screening_data <- screening_data[
    order(screening_data$id, screening_data$time),
    ,
    drop = FALSE
  ]
  screening_data$time <- factor(screening_data$time, levels = 0:3)
  rownames(screening_data) <- NULL

  list(
    numeric = list(model = fit(data, 1), times = 1:3),
    character = list(
      model = fit(character_data, "Week 2"),
      times = week_labels
    ),
    screening = list(
      model = fit(screening_data, "1"),
      times = c("1", "2", "3")
    )
  )
}

start_row_cases <- function() {
  grid_fixture("bootstrap start row cases", build_start_row_cases)
}

start_row_intervals <- function(case, method, estimand) {
  result <- if (identical(estimand, "sops")) {
    avg_sops(case$model, variables = list(tx = c(0, 1)), times = case$times)
  } else {
    avg_comparisons(
      case$model,
      variables = list(tx = c(0, 1)),
      estimand = "time_benefit",
      times = case$times
    )
  }
  result <- inferences(
    result,
    method = method,
    n_draws = 5,
    workers = 1,
    seed = 2147
  )
  as.data.frame(result)[c("estimate", "conf.low", "conf.high")]
}

test_that("bootstrap draws start from first_followup_time for character time", {
  skip_on_cran() # Slow: bootstrap refits. CI still checks #47.
  cases <- start_row_cases()
  expect_equal(
    start_row_intervals(cases$character, "bootstrap", "sops"),
    start_row_intervals(cases$numeric, "bootstrap", "sops"),
    tolerance = 1e-8
  )
})

test_that("bootstrap and FWB draws start from first_followup_time", {
  skip_on_cran() # Repeats the CRAN check above for the other routes.
  cases <- start_row_cases()
  for (case_name in c("character", "screening")) {
    for (method in c("bootstrap", "fwb")) {
      for (estimand in c("sops", "time_benefit")) {
        expect_equal(
          start_row_intervals(cases[[case_name]], method, estimand),
          start_row_intervals(cases$numeric, method, estimand),
          tolerance = 1e-8,
          label = paste(case_name, method, estimand)
        )
      }
    }
  }
})
