# A patient cannot start in an absorbing state. SOP calls refuse a profile whose
# previous state, or for second-order models the state two visits back, is one
# of the states in `absorb`. Without the refusal, a first-order fit predicts a
# transition out of the absorbing state from a previous state it never saw (a
# factor `yprev` fails with an unrelated error, a numeric `yprev` extends the
# fitted slope), while a second-order fit silently keeps the patient there.

# The error names the absorbing state, the column and how many profiles start
# there.
expect_absorbing_start_error <- function(expr, column, count, info) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "error")
  if (!inherits(err, "error")) {
    return(invisible())
  }
  message <- conditionMessage(err)
  expect_match(message, "absorbing", info = info)
  expect_match(message, paste0("`", column, "`"), fixed = TRUE, info = info)
  expect_match(message, paste0("\\b", count, "\\b"), info = info)
}

# Profiles of the comparisons case with the first `count` patients started in
# the absorbing state 3.
# `prepare_markov_data()` drops the transitions out of the absorbing state, so
# a factor previous state may lack the level "3"; it is added first, or the
# assignment would give missing values instead.
absorbing_start_profiles <- function(profiles, count, column = "yprev") {
  value <- if (is.factor(profiles[[column]])) {
    levels(profiles[[column]]) <- union(levels(profiles[[column]]), "3")
    factor("3", levels = levels(profiles[[column]]))
  } else {
    3
  }
  profiles[[column]][seq_len(count)] <- value
  profiles
}

# Every public SOP call type on the same profiles. `...` takes further
# arguments for every call, such as `p2_var` for a second-order model.
absorbing_start_calls <- function(model, profiles, ...) {
  base <- list(model, newdata = profiles, times = 1:4, absorb = 3, ...)
  compare <- c(base, list(variables = list(tx = c(0, 1))))
  list(
    sops = function() do.call(sops, base),
    avg_sops = function() do.call(avg_sops, base),
    state_time = function() {
      do.call(state_time, c(base, list(state_sets = "1")))
    },
    avg_state_time = function() {
      do.call(avg_state_time, c(base, list(state_sets = "1")))
    },
    comparisons = function() {
      do.call(comparisons, c(compare, list(state_sets = "1")))
    },
    avg_comparisons_sop = function() {
      do.call(avg_comparisons, c(compare, list(state_sets = "1")))
    },
    avg_comparisons_time_benefit = function() {
      do.call(avg_comparisons, c(compare, list(estimand = "time_benefit")))
    }
  )
}

test_that("every SOP call refuses profiles starting in an absorbing state", {
  skip_on_cran()
  case <- make_comparisons_case()
  second <- make_second_order_comparisons_case()

  # A numeric previous state is compared with `absorb` as a label.
  numeric_data <- case$data
  numeric_data$yprev <- as.numeric(as.character(numeric_data$yprev))
  numeric_model <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = numeric_data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
  )

  # The most recent state is absorbing, as is the state two visits back.
  both <- absorbing_start_profiles(second$profiles, 2)
  both <- absorbing_start_profiles(both, 2, column = "ypprev")

  # `column` and `count` are the column and number of profiles the error names.
  cases <- list(
    list(
      label = "first order, factor yprev",
      model = case$model,
      profiles = absorbing_start_profiles(case$profiles, 2),
      column = "yprev",
      count = 2
    ),
    list(
      label = "first order, numeric yprev",
      model = numeric_model,
      profiles = absorbing_start_profiles(
        numeric_data[numeric_data$time == 1, , drop = FALSE],
        3
      ),
      column = "yprev",
      count = 3
    ),
    list(
      label = "second order, both previous states",
      model = second$model,
      profiles = both,
      p2_var = "ypprev",
      column = "yprev",
      count = 2
    ),
    # Only the state two visits back is absorbing: an impossible history, since
    # the patient would have had to leave the absorbing state.
    list(
      label = "second order, ypprev only",
      model = second$model,
      profiles = absorbing_start_profiles(second$profiles, 1, column = "ypprev"),
      p2_var = "ypprev",
      column = "ypprev",
      count = 1
    )
  )

  for (row in cases) {
    extra <- if (!is.null(row$p2_var)) list(p2_var = row$p2_var)
    calls <- do.call(absorbing_start_calls, c(list(row$model, row$profiles), extra))
    for (name in names(calls)) {
      expect_absorbing_start_error(
        calls[[name]](),
        row$column,
        row$count,
        paste0(row$label, ": ", name)
      )
    }
  }
})

test_that("SOPs from the profiles stored by mostm() refuse an absorbing start", {
  case <- make_comparisons_case()
  data <- case$data
  # Patients 1 and 2 are recorded in the absorbing state before their first
  # follow-up visit.
  levels(data$yprev) <- union(levels(data$yprev), "3")
  data$yprev[data$id %in% 1:2 & data$time == 1] <- "3"
  model <- suppressWarnings(mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
  ))
  err <- tryCatch(
    avg_sops(model, times = 1:4, absorb = 3),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(
    conditionMessage(err),
    "`yprev` is the absorbing state 3 for 2 of the 60 starting profiles stored by `mostm()`",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "(`id` 1, 2)", fixed = TRUE)
})

test_that("soprob_markov() refuses prediction rows with an absorbing start", {
  case <- make_comparisons_case()
  rows <- absorbing_start_profiles(case$profiles, 2)
  expect_error(
    soprob_markov(
      case$model,
      newdata = rows,
      times = 1:4,
      y_levels = 1:3,
      absorb = 3
    ),
    "`yprev` is the absorbing state 3 for 2 of the 60 prediction rows",
    fixed = TRUE
  )
  # `variables` that set the previous state are checked on the stacked rows.
  expect_error(
    avg_sops(
      case$model,
      newdata = case$profiles,
      variables = list(yprev = factor(c("2", "3"), levels = 1:3)),
      times = 1:4,
      absorb = 3
    ),
    "if `variables` sets `yprev`, use non-absorbing states only",
    fixed = TRUE
  )
})
