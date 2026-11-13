test_that("avg_sops averages observed covariates overall and within groups", {
  skip_on_cran()
  case <- make_observed_average_case()
  individual <- sops(case$model, times = case$times)
  for (by in list(TRUE, "tx")) {
    average <- avg_sops(
      case$model,
      by = by,
      times = case$times
    )
    keys <- c("time", "state", if (is.character(by)) by)
    manual <- stats::aggregate(individual$estimate, individual[keys], mean)
    names(manual)[names(manual) == "x"] <- "manual"
    compared <- merge(average, manual, by = keys)
    expect_equal(compared$estimate, compared$manual, tolerance = 1e-12)
    expect_named(average, c(keys, "estimate"))
    expect_null(attr(average, "avg_args")$variables)
  }
  empty_by <- avg_sops(
    case$model,
    by = character(),
    times = case$times
  )
  explicit <- avg_sops(
    case$model,
    newdata = attr(individual, "newdata_pred"),
    variables = NULL,
    times = case$times
  )
  expect_equal(empty_by$estimate, explicit$estimate)
  expect_null(attr(empty_by, "avg_args")$by)
})

test_that("observed averages retain conditional and unconditional covariance", {
  skip_on_cran()
  case <- make_observed_average_case()
  average <- avg_sops(case$model, times = case$times)
  individual <- sops(case$model, times = case$times) |>
    inferences(method = "delta")
  conditional <- inferences(average, method = "delta", vcov = "conditional")
  cell <- match(
    paste(individual$time, individual$state),
    paste(average$time, average$state)
  )
  summed <- rowsum(get_jacobian(individual), cell)
  jacobian <- summed[as.character(seq_len(nrow(average))), , drop = FALSE] /
    length(unique(individual$rowid))
  expected <- jacobian %*% stats::vcov(case$model) %*% t(jacobian)
  expect_equal(unname(stats::vcov(conditional)), unname(expected))
  for (vcov in c("conditional", "unconditional")) {
    result <- inferences(average, method = "delta", vcov = vcov)
    expect_equal(result$estimate, average$estimate)
    expect_equal(unname(diag(stats::vcov(result))), result$std.error^2)
    expect_equal(all(is.finite(result$std.error)), TRUE)
    sums <- rowsum(stats::vcov(result), result$time)
    expect_equal(
      unname(sums),
      matrix(0, nrow(sums), ncol(sums)),
      tolerance = 1e-10
    )
  }
})

test_that("observed averages support both refit bootstrap methods", {
  skip_on_cran()
  case <- make_observed_average_case()
  average <- avg_sops(case$model, times = case$times)
  for (method in c("bootstrap", "fwb")) {
    result <- inferences(average, method = method, n_draws = 3, seed = 5211)
    draws <- get_draws(result)
    expect_named(
      result,
      c(
        "time",
        "state",
        "estimate",
        "std.error",
        "conf.low",
        "conf.high"
      )
    )
    expect_equal(all(is.finite(result$std.error)), TRUE)
    expect_equal(nrow(draws), nrow(average) * 3L)
    totals <- stats::aggregate(draw ~ time + draw_id, draws, sum)
    expect_equal(totals$draw, rep(1, nrow(totals)), tolerance = 1e-12)
    anchors <- attr(result, "baseline_anchor_draws")
    expect_equal(sort(unique(anchors$draw_id)), 1:3)
  }
})
