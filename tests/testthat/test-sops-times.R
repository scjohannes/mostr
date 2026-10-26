# `times` lists consecutive visits: the SOP recursion takes one transition per
# value. Requests that skip fitted visits silently change what each result
# means, so the user-facing SOP functions warn about them.

make_daily_visit_case <- function(time_as_factor = FALSE) {
  skip_if_not_installed("rms")
  data <- make_test_data(n_patients = 80, follow_up_time = 8, seed = 2612)
  if (time_as_factor) {
    data$time <- factor(data$time, levels = as.character(1:8))
  }
  list(
    data = data,
    model = mostm(
      y ~ time + tx + yprev,
      data = data,
      id_var = "id",
      first_followup_time = if (time_as_factor) "1"
    )
  )
}

# Visits on days 1, 2, 4, and 8, with `gap` holding the days since the
# previous visit (day 0 for the first visit).
make_irregular_visit_case <- function(visit_days = c(1, 2, 4, 8)) {
  skip_if_not_installed("rms")
  raw <- make_markov_trajectories(
    n_patients = 120,
    follow_up_time = max(visit_days),
    seed = 2613
  )
  visits <- raw[raw$time %in% visit_days, , drop = FALSE]
  visits <- visits[order(visits$id, visits$time), , drop = FALSE]
  first <- !duplicated(visits$id)
  previous <- c(NA, seq_len(nrow(visits) - 1L))
  visits$yprev[!first] <- visits$y[previous[!first]]
  visits$gap <- ifelse(
    first,
    visits$time,
    visits$time - visits$time[previous]
  )
  prepare_markov_data(visits)
}

count_skipped_visit_warnings <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(
    expr,
    mostr_warning_skipped_visits = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, n = n)
}

death_at <- function(x, time, tx) {
  x$estimate[x$state == "6" & x$time == time & x$tx == tx]
}

test_that("avg_sops() warns when `times` skips fitted visits", {
  case <- make_daily_visit_case()
  variables <- list(tx = c(0, 1))

  expect_no_warning(
    full <- avg_sops(case$model, variables = variables, times = 1:8, absorb = "6")
  )
  expect_warning(
    sparse <- avg_sops(
      case$model,
      variables = variables,
      times = c(1, 4, 8),
      absorb = "6"
    ),
    class = "mostr_warning_skipped_visits"
  )

  # The result labelled 8 is the probability after three transitions, not at
  # visit 8; subsetting the consecutive request gives the visit-8 estimate.
  expect_equal(death_at(sparse, 1, 0), death_at(full, 1, 0))
  expect_gt(death_at(full, 8, 0) - death_at(sparse, 8, 0), 0.1)
})

test_that("the skipped-visit warning explains `times` and the remedy", {
  case <- make_daily_visit_case()

  warning <- tryCatch(
    avg_sops(case$model, times = c(1, 4, 8), absorb = "6"),
    mostr_warning_skipped_visits = function(w) w
  )

  message <- conditionMessage(warning)
  expect_match(message, "skips fitted visits 2, 3, 5, 6, and 7", fixed = TRUE)
  expect_match(message, "labelled 8 is the probability after 3 transitions")
  expect_match(message, "`times = 1:8`", fixed = TRUE)
  expect_match(message, "`gap_var`", fixed = TRUE)
})

test_that("sops() warns when stored starting profiles do not start at the first visit", {
  case <- make_daily_visit_case()

  expect_warning(
    sops(case$model, times = 8, absorb = "6"),
    "skips fitted visits 1, 2, 3, 4, 5, and 2 more"
  )
  expect_warning(
    sops(case$model, times = 4:8, absorb = "6"),
    class = "mostr_warning_skipped_visits"
  )
  expect_warning(
    sops(case$model, times = 0:8, absorb = "6"),
    "starts at 0, before the first follow-up visit (1)",
    fixed = TRUE
  )
  expect_warning(
    sops(case$model, times = c(1, 3, 2, 4:8), absorb = "6"),
    "not in increasing order"
  )
  expect_no_warning(sops(case$model, times = 1:8, absorb = "6"))
})

test_that("supplied newdata may start later but must not skip visits", {
  case <- make_daily_visit_case()
  profiles <- case$data[case$data$time == 1, , drop = FALSE]

  # Each row's previous state is the state before the first requested visit.
  expect_no_warning(
    sops(case$model, newdata = profiles, times = 4:8, absorb = "6")
  )
  expect_warning(
    sops(case$model, newdata = profiles, times = c(4, 8), absorb = "6"),
    "skips fitted visits 5, 6, and 7"
  )
  expect_warning(
    avg_sops(case$model, newdata = profiles, times = c(1, 4, 8), absorb = "6"),
    class = "mostr_warning_skipped_visits"
  )
})

test_that("times beyond the fitted visits are not treated as skipped visits", {
  case <- make_daily_visit_case()

  expect_no_warning(avg_sops(case$model, times = 1:12, absorb = "6"))
  expect_warning(
    avg_sops(case$model, times = c(1, 4, 8, 12), absorb = "6"),
    class = "mostr_warning_skipped_visits"
  )
})

test_that("explicit time covariates do not change the skipped-visit check", {
  case <- make_daily_visit_case()
  fit <- mostm(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = case$data,
    id_var = "id",
    engine = "vglm"
  )
  time_covariates <- make_time_covariates(
    case$data,
    "time",
    "time_lin",
    "time_nlin_1"
  )

  expect_no_warning(avg_sops(
    fit,
    times = 1:8,
    absorb = "6",
    time_covariates = time_covariates
  ))
  expect_warning(
    avg_sops(
      fit,
      times = c(1, 4, 8),
      absorb = "6",
      time_covariates = time_covariates[c(1, 4, 8), , drop = FALSE]
    ),
    class = "mostr_warning_skipped_visits"
  )
})

test_that("irregular visit schedules are checked against the fitted visits", {
  data <- make_irregular_visit_case()
  fit <- mostm(y ~ time + tx + yprev, data = data, id_var = "id")

  expect_no_warning(avg_sops(fit, times = c(1, 2, 4, 8), absorb = "6"))
  expect_warning(
    avg_sops(fit, times = c(1, 4, 8), absorb = "6"),
    "skips fitted visit 2\\."
  )
  expect_warning(
    avg_sops(fit, times = c(1, 4, 8, 12), absorb = "6"),
    "request every fitted visit from 1 to 8, then the later times",
    fixed = TRUE
  )
})

test_that("`gap_var` allows times that skip fitted visits", {
  data <- make_irregular_visit_case()
  fit <- mostm(y ~ time + gap + tx + yprev, data = data, id_var = "id")

  expect_no_warning(
    avg_sops(fit, times = c(1, 4, 8), absorb = "6", gap_var = "gap")
  )
  expect_no_warning(
    sops(fit, times = c(1, 4, 8), absorb = "6", gap_var = "gap")
  )
})

test_that("factor visit times warn when levels are skipped", {
  case <- make_daily_visit_case(time_as_factor = TRUE)

  expect_no_warning(avg_sops(case$model, times = as.character(1:8), absorb = "6"))
  expect_warning(
    avg_sops(case$model, times = c("1", "4", "8"), absorb = "6"),
    "skips fitted visits 2, 3, 5, 6, and 7"
  )
})

test_that("summaries built on avg_sops() warn once and not again during inference", {
  case <- make_daily_visit_case()
  variables <- list(tx = c(0, 1))

  time_total <- count_skipped_visit_warnings(
    avg_state_time(
      case$model,
      variables = variables,
      times = c(1, 4, 8),
      absorb = "6"
    )
  )
  comparison <- count_skipped_visit_warnings(avg_comparisons(
    case$model,
    variables = variables,
    times = c(1, 4, 8),
    absorb = "6"
  ))
  benefit <- count_skipped_visit_warnings(avg_comparisons(
    case$model,
    variables = variables,
    estimand = "time_benefit",
    times = c(1, 4, 8),
    absorb = "6"
  ))

  expect_equal(c(time_total$n, comparison$n, benefit$n), c(1L, 1L, 1L))

  expect_no_warning(inferences(time_total$value, method = "delta"))
  expect_no_warning(inferences(comparison$value, method = "delta"))
  expect_no_warning(
    inferences(benefit$value, method = "mvn", n_draws = 5)
  )
  expect_no_warning(
    avg_state_time(case$model, variables = variables, times = 1:8, absorb = "6")
  )
})

test_that("comparisons() and state_time() check `times` against the stored starting profiles", {
  case <- make_daily_visit_case()

  # Without `newdata`, each patient starts from the state before visit 1, so
  # `times` must begin at visit 1 here too.
  expect_warning(
    comparisons(case$model, variables = "tx", times = 4:8, absorb = "6"),
    "skips fitted visits 1, 2, and 3"
  )
  expect_warning(
    comparisons(
      case$model,
      variables = "tx",
      estimand = "time_benefit",
      times = 4:8,
      absorb = "6"
    ),
    "skips fitted visits 1, 2, and 3"
  )
  expect_warning(
    state_time(case$model, times = 4:8, absorb = "6"),
    "skips fitted visits 1, 2, and 3"
  )
  expect_warning(
    state_time(case$model, variables = "tx", times = 0:8, absorb = "6"),
    "starts at 0, before the first follow-up visit (1)",
    fixed = TRUE
  )

  # Supplied rows may hold the state before a later visit.
  profiles <- case$data[case$data$time == 1, , drop = FALSE]
  expect_no_warning(comparisons(
    case$model,
    newdata = profiles,
    variables = "tx",
    times = 4:8,
    absorb = "6"
  ))
  expect_no_warning(
    state_time(case$model, newdata = profiles, times = 4:8, absorb = "6")
  )
  expect_warning(
    state_time(case$model, newdata = profiles, times = c(4, 8), absorb = "6"),
    "skips fitted visits 5, 6, and 7"
  )
})

test_that("comparisons() and state_time() warn once and not again during inference", {
  case <- make_daily_visit_case()

  late <- count_skipped_visit_warnings(
    comparisons(case$model, variables = "tx", times = 4:8, absorb = "6")
  )
  sparse <- count_skipped_visit_warnings(
    comparisons(case$model, variables = "tx", times = c(1, 4, 8), absorb = "6")
  )
  time_total <- count_skipped_visit_warnings(
    state_time(case$model, times = 4:8, absorb = "6")
  )

  expect_equal(c(late$n, sparse$n, time_total$n), c(1L, 1L, 1L))

  expect_no_warning(inferences(late$value, method = "delta"))
  expect_no_warning(inferences(time_total$value, method = "mvn", n_draws = 5))
  expect_no_warning(
    comparisons(case$model, variables = "tx", times = 1:8, absorb = "6")
  )
  expect_no_warning(state_time(case$model, times = 1:8, absorb = "6"))
})
