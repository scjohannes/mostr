# Each profile's SOPs for a state set, summed over visits, from `sops()` on
# the profiles with `var` set to `value`. Returned in profile order.
profile_visit_totals <- function(case, var, value, states, ...) {
  ind <- sops_with_value(case, var, value, ...)
  keep <- as.character(ind$state) %in% states
  totals <- stats::aggregate(estimate ~ rowid, ind[keep, ], sum)
  totals$estimate[order(totals$rowid)]
}

test_that("state_time() totals each profile's SOPs under each value", {
  case <- make_comparisons_case()
  n <- nrow(case$profiles)
  out <- state_time(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    state_sets = list(alive = c("1", "2")),
    times = 1:4,
    absorb = 3
  )
  expect_s3_class(out, "markov_state_time")
  expect_identical(nrow(out), 2L * n)
  expect_identical(out$rowid, seq_len(2L * n))
  expect_identical(out$rowidcf, rep(seq_len(n), 2L))
  expect_true(all(out$state_set == "alive"))
  for (value in c(0, 1)) {
    rows <- out[out$tx == value, , drop = FALSE]
    expect_equal(
      rows$estimate,
      profile_visit_totals(case, "tx", value, c("1", "2"))
    )
    expect_equal(rows$marker, case$profiles$marker)
  }
  expect_false(any(startsWith(names(out), ".markov_")))
  expect_false("time" %in% names(out))
  expect_identical(names(out)[ncol(out)], "estimate")

  plain <- state_time(case$model, times = 1:4, absorb = 3)
  expect_false("rowidcf" %in% names(plain))
  expect_setequal(plain$state_set, as.character(1:3))
})

test_that("posterior state_time() summarizes each profile's total per draw", {
  coefficients <- rbind(
    c(0.5, -0.5, 0.4, 1.2),
    c(2, 1, -0.8, -1.1),
    c(-0.4, -2, 1.5, 0.2)
  )
  colnames(coefficients) <- c("y>=2", "y>=3", "tx", "time")
  model <- fake_blrm(coefficients)
  profiles <- data.frame(
    id = c("a", "b"),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1
  )
  local_mocked_bindings(blrm_design_matrix = fake_blrm_design(c("tx", "time")))
  out <- state_time(
    model,
    newdata = profiles,
    state_sets = "1",
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "median",
    return_draws = TRUE
  )
  ind <- sops(
    model,
    newdata = profiles,
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    return_draws = TRUE
  )
  draws <- attr(ind, "draws")
  totals <- stats::aggregate(
    estimate ~ rowid + draw_id,
    draws[as.character(draws$state) == "1", ],
    sum
  )
  for (i in seq_len(nrow(out))) {
    values <- totals$estimate[totals$rowid == out$rowid[i]]
    expect_equal(out$estimate[i], stats::median(values))
    expect_equal(out$std.error[i], stats::sd(values))
  }
  expect_identical(nrow(get_draws(out)), 3L * nrow(out))
})
