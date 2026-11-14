# A `by` column with missing values among the averaged profiles is refused.
# Dropping those patients silently changes who is averaged, so the user must
# drop or impute them first. Missing values in rows that are not averaged,
# such as later visits of the fitted data, do not matter.

# The comparisons profiles with a grouping column that is missing for the
# first ten patients.
profiles_with_missing_group <- function(case) {
  profiles <- case$profiles
  profiles$grp <- ifelse(
    profiles$id <= 10,
    NA_character_,
    ifelse(profiles$id %% 2 == 0, "a", "b")
  )
  profiles
}

# The error names the column and the number of affected profiles.
expect_missing_by_error <- function(expr, info) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "error")
  if (!inherits(err, "error")) {
    return(invisible())
  }
  message <- conditionMessage(err)
  expect_match(message, "`grp`", fixed = TRUE, info = info)
  expect_match(message, "\\b10\\b", info = info)
  expect_match(message, "missing", info = info)
}

test_that("averages refuse a `by` column with missing values in newdata", {
  case <- make_comparisons_case()
  profiles <- profiles_with_missing_group(case)
  base <- list(case$model, newdata = profiles, times = 1:4, by = "grp")
  compare <- c(base, list(variables = list(tx = c(0, 1))))
  calls <- list(
    avg_sops = function() do.call(avg_sops, base),
    sops = function() do.call(sops, base),
    avg_state_time = function() {
      do.call(avg_state_time, c(base, list(state_sets = "1")))
    },
    state_time = function() do.call(state_time, c(base, list(state_sets = "1"))),
    avg_comparisons_sop = function() {
      do.call(avg_comparisons, c(compare, list(state_sets = "1")))
    },
    avg_comparisons_time_in_state = function() {
      do.call(
        avg_comparisons,
        c(compare, list(estimand = "time_in_state", state_sets = "1"))
      )
    },
    avg_comparisons_time_benefit = function() {
      do.call(avg_comparisons, c(compare, list(estimand = "time_benefit")))
    },
    comparisons = function() {
      do.call(comparisons, c(compare, list(state_sets = "1")))
    }
  )
  for (name in names(calls)) {
    expect_missing_by_error(calls[[name]](), name)
  }
})

test_that("averages over the stored patients refuse missing `by` values in starting rows", {
  case <- make_comparisons_case()
  data <- case$data
  data$grp <- ifelse(data$id <= 10, NA_character_, ifelse(data$id %% 2 == 0, "a", "b"))
  model <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    absorb = 3,
    engine = "vglm"
  )
  expect_missing_by_error(
    avg_sops(model, times = 1:4, by = "grp"),
    "stored avg_sops"
  )
  expect_missing_by_error(
    avg_comparisons(
      model,
      variables = list(tx = c(0, 1)),
      times = 1:4,
      by = "grp",
      state_sets = "1"
    ),
    "stored avg_comparisons"
  )
})

test_that("missing `by` values outside the starting rows do not matter", {
  case <- make_comparisons_case()
  data <- case$data
  data$grp <- ifelse(data$id %% 2 == 0, "a", "b")
  # Missing only at later visits, which are not prediction profiles.
  data$grp[data$time > 1 & data$id <= 10] <- NA_character_
  model <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    absorb = 3,
    engine = "vglm"
  )
  out <- avg_sops(model, times = 1:4, by = "grp")
  expect_setequal(unique(out$grp), c("a", "b"))
  starts <- data[data$time == 1, , drop = FALSE]
  for (value in c("a", "b")) {
    expected <- avg_sops(
      model,
      newdata = starts[starts$grp == value, , drop = FALSE],
      times = 1:4
    )
    got <- out[out$grp == value, , drop = FALSE]
    key <- function(x) paste(x$time, x$state)
    expect_equal(got$estimate[match(key(expected), key(got))], expected$estimate)
  }
})

test_that("complete `by` columns, `by = TRUE`, and missing values in other columns still work", {
  case <- make_comparisons_case()
  profiles <- profiles_with_missing_group(case)
  profiles$grp2 <- ifelse(profiles$id %% 2 == 0, "a", "b")
  expect_no_error(
    avg_sops(case$model, newdata = profiles, times = 1:4, by = "grp2")
  )
  expect_no_error(
    avg_sops(case$model, newdata = profiles, times = 1:4, by = TRUE)
  )
  expect_no_error(
    avg_comparisons(
      case$model,
      newdata = profiles,
      variables = list(tx = c(0, 1)),
      times = 1:4,
      by = "tx",
      state_sets = "1"
    )
  )
})
