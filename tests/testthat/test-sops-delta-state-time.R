test_that("delta state_time() uses each profile's time derivatives", {
  case <- make_comparisons_case()
  for (time_map in list(NULL, case$time_map)) {
    args <- list(
      case$model,
      variables = list(tx = c(0, 1)),
      state_sets = list(alive = c("1", "2")),
      times = 1:4,
      absorb = 3,
      time_map = time_map
    )
    unit <- inferences(do.call(state_time, args), method = "delta")
    expect_identical(attr(unit, "method"), "delta")
    expect_identical(attr(unit, "conf_type"), "wald")
    expect_null(attr(unit, "null"))
    average <- inferences(
      do.call(avg_state_time, args),
      method = "delta",
      vcov = "conditional"
    )
    unit_jacobian <- get_jacobian(unit)
    average_jacobian <- get_jacobian(average)
    for (value in c(0, 1)) {
      expect_equal(
        colMeans(unit_jacobian[unit$tx == value, , drop = FALSE]),
        average_jacobian[average$tx == value, ],
        tolerance = 1e-10
      )
    }
    V <- attr(unit, "analytical")$coefficient_vcov
    expect_equal(
      unit$std.error,
      unname(sqrt(rowSums((unit_jacobian %*% V) * unit_jacobian))),
      tolerance = 1e-10
    )
    block <- vcov(unit, rows = 1:2)
    expect_equal(unname(sqrt(diag(block))), unit$std.error[1:2])
  }

  plain <- inferences(
    state_time(case$model, state_sets = "1", times = 1:4, absorb = 3),
    method = "delta"
  )
  expect_equal(
    colMeans(get_jacobian(plain)),
    get_jacobian(inferences(
      avg_state_time(case$model, state_sets = "1", times = 1:4, absorb = 3),
      method = "delta",
      vcov = "conditional"
    ))[1, ],
    tolerance = 1e-10
  )
  expect_error(
    inferences(
      state_time(case$model, times = 1:4, absorb = 3),
      method = "delta",
      vcov = "unconditional"
    ),
    "`state_time()` and `avg_state_time(by = FALSE)`",
    fixed = TRUE
  )
})
