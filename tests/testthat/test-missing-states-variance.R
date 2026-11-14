# Standard errors of fits with unobserved states: one score per patient,
# summed over the patient's fully observed transitions and stretches, and the
# patient-clustered sandwich built from them.

test_that("patient scores are the gradients of each patient's log-likelihood", {
  skip_on_cran() # Slow: numerical gradients of the reference.
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  info <- attr(fit, "markov_path_fit")
  scores <- info$scores
  expect_identical(rownames(scores), info$ids)
  expect_identical(colnames(scores), names(coef(fit)))
  # One row per patient of the fit, including the patient whose first rows
  # have no previous state.
  expect_setequal(info$ids, as.character(unique(data$id)))
  # At the maximum the scores add up to the total score, which is ~0.
  expect_equal(colSums(scores), fit$u, tolerance = 1e-8)
  expect_lt(max(abs(colSums(scores))), 1e-6)

  special <- attr(data, "special")
  patients <- as.character(c(
    special$alive_long[3], # gap of length 3
    special$alive_long[4], # two gaps
    special$dead_late[1], # gap before death
    special$trailing, # gap at the end of follow-up
    special$alive_long[10] # no unobserved visit
  ))
  reference <- ms_reference(data)
  numeric <- central_difference_jacobian(
    function(theta) reference(theta, per_patient = TRUE)[patients],
    coef(fit)
  )
  expect_equal(
    unname(scores[patients, ]),
    unname(numeric),
    tolerance = 1e-6
  )
})

test_that("the robust covariance is the patient-clustered sandwich", {
  data <- ms_fixture("gaps")
  fit <- ms_fixture("gap_fit")
  info <- attr(fit, "markov_path_fit")
  scores <- info$scores
  G <- nrow(scores)
  bread <- fit$orig.var
  expect_equal(unname(bread), unname(solve(info$information)), tolerance = 1e-10)
  sandwich <- G / (G - 1) * bread %*% crossprod(scores) %*% bread
  expect_equal(unname(fit$var), unname(sandwich), tolerance = 1e-10)
  expect_equal(
    unname(vcov(fit, intercepts = "all")),
    unname(fit$var),
    tolerance = 1e-12
  )
  # Without the small-sample factor.
  fit_raw <- ms_fit(data, cadjust = FALSE)
  expect_equal(
    unname(fit_raw$var),
    unname(bread %*% crossprod(scores) %*% bread),
    tolerance = 1e-10
  )
})

test_that("delta-method SEs work with both variances", {
  skip_on_cran() # Slow: predictions with both variances.
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  out <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:3)
  conditional <- inferences(out, method = "delta", vcov = "conditional")
  unconditional <- inferences(out, method = "delta", vcov = "unconditional")
  for (result in list(conditional, unconditional)) {
    expect_true(all(is.finite(result$std.error)))
    expect_true(all(result$std.error > 0))
    expect_equal(result$estimate, out$estimate)
  }
  # The unconditional variance adds the variation of the starting states, so
  # it is not the conditional one.
  expect_false(isTRUE(all.equal(conditional$std.error, unconditional$std.error)))
})

test_that("another `cluster` can't be supplied to a fit with unobserved states", {
  skip_on_cran() # Slow: predictions.
  data <- ms_fixture("downstream")
  fit <- ms_fixture("downstream_fit")
  out <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:2)
  expect_error(
    inferences(out, method = "delta", vcov = "unconditional", cluster = data$id),
    "another `cluster` can't be supplied",
    class = "mostr_error"
  )
  fitted_rows <- attr(fit, "markov_path_fit")$n_observed_transitions
  expect_error(
    inferences(
      out,
      method = "delta",
      vcov = "conditional",
      cluster = rep(1:2, length.out = fitted_rows)
    ),
    "another `cluster` can't be supplied",
    class = "mostr_error"
  )
})
