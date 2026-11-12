# Summaries of SOP tables (`time_in_state()`, `interpolate_sops()`,
# `state_time()`) keep one group per profile, treatment value, or `by` group.
# Columns that only describe a row, such as a profile column the model does
# not use, or the test columns added by `inferences(null = )`, must not split
# or drop those groups.

# The comparisons profiles with a column the model does not use, missing for
# the first five profiles.
profiles_with_unused_note <- function(case) {
  profiles <- case$profiles
  profiles$note <- ifelse(seq_len(nrow(profiles)) <= 5L, NA_character_, "ok")
  profiles
}

test_that("a missing value in an unused profile column keeps the profile", {
  skip_on_cran()
  case <- make_comparisons_case()
  plain <- case$profiles
  noted <- profiles_with_unused_note(case)
  n <- nrow(plain)
  s_plain <- sops(case$model, newdata = plain, times = 1:4, absorb = 3)
  s_noted <- sops(case$model, newdata = noted, times = 1:4, absorb = 3)

  for (time_map in list(NULL, case$time_map)) {
    expected <- time_in_state(s_plain, target_states = "1", time_map = time_map)
    out <- time_in_state(s_noted, target_states = "1", time_map = time_map)
    expect_identical(nrow(out), n)
    out <- out[order(out$rowid), , drop = FALSE]
    expected <- expected[order(expected$rowid), , drop = FALSE]
    expect_equal(out$total_time, expected$total_time)
    expect_identical(out$note, noted$note)
  }

  expected <- interpolate_sops(s_plain, time_map = case$time_map, target_times = c(1, 3, 5, 7))
  out <- interpolate_sops(s_noted, time_map = case$time_map, target_times = c(1, 3, 5, 7))
  expect_identical(nrow(out), nrow(expected))
  key <- function(x) paste(x$rowid, x$time, x$state)
  expect_equal(out$estimate[match(key(expected), key(out))], expected$estimate)
})

test_that("inferences(null = ) does not change interpolate_sops() or time_in_state() groups", {
  skip_on_cran()
  case <- make_comparisons_case()
  averaged <- avg_sops(
    case$model,
    times = 1:4,
    absorb = 3,
    variables = list(tx = c(0, 1))
  )
  target_times <- c(1, 3, 5, 7)
  # Delta-method results have no draws to rebuild intervals from, and the
  # summaries say so in a warning that is beside the point here.
  without_draws_warning <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("no stored draws", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  for (method in c("delta", "mvn")) {
    extra <- if (method == "mvn") list(n_draws = 30, seed = 1) else list()
    without <- do.call(inferences, c(list(averaged, method = method), extra))
    with <- do.call(inferences, c(list(averaged, method = method, null = 0.5), extra))

    expected <- without_draws_warning(
      interpolate_sops(without, time_map = case$time_map, target_times = target_times)
    )
    out <- without_draws_warning(
      interpolate_sops(with, time_map = case$time_map, target_times = target_times)
    )
    # 2 treatment values x 3 states x 4 target times.
    expect_identical(nrow(out), nrow(expected))
    expect_identical(nrow(out), 2L * 3L * length(target_times))
    key <- function(x) paste(x$tx, x$state, x$time)
    index <- match(key(expected), key(out))
    expect_false(anyNA(index))
    for (column in intersect(c("estimate", "std.error", "conf.low", "conf.high"), names(expected))) {
      expect_equal(out[[column]][index], expected[[column]], info = paste(method, column))
    }
    # The tests refer to the interpolated estimates, against the same null.
    tested <- !is.na(out$std.error) & out$std.error > 0
    expect_equal(
      out$statistic[tested],
      (out$estimate[tested] - 0.5) / out$std.error[tested],
      info = method
    )
    expect_equal(out$p.value[tested], 2 * stats::pnorm(-abs(out$statistic[tested])), info = method)

    expected <- without_draws_warning(
      time_in_state(without, target_states = "1", time_map = case$time_map)
    )
    out <- without_draws_warning(
      time_in_state(with, target_states = "1", time_map = case$time_map)
    )
    expect_identical(nrow(out), 2L)
    out <- out[order(out$tx), , drop = FALSE]
    expected <- expected[order(expected$tx), , drop = FALSE]
    for (column in intersect(c("total_time", "std.error", "conf.low", "conf.high"), names(expected))) {
      expect_equal(out[[column]], expected[[column]], info = paste(method, column))
    }
    # A null value for probabilities says nothing about a total time.
    expect_false(any(c("statistic", "p.value", "s.value") %in% names(out)), info = method)
  }
})
