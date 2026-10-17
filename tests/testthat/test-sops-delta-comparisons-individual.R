test_that("delta intervals difference each profile's SOP derivatives", {
  case <- make_comparisons_case()
  levels <- c(0, 1, 2)
  state_sets <- list(alive = c("1", "2"), dead = "3")
  cmp <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = levels),
    state_sets = state_sets,
    times = 1:4,
    absorb = 3
  )
  inferred <- inferences(cmp, method = "delta")
  expect_s3_class(inferred, "markov_comparisons")
  expect_identical(attr(inferred, "method"), "delta")
  expect_identical(attr(inferred, "conf_type"), "wald")
  expect_identical(attr(inferred, "null"), 0)
  expect_null(attr(inferred, "draws"))
  expect_equal(inferred$estimate, cmp$estimate)

  # Reference derivatives: analytical `sops()` on the profiles stacked once
  # per value of `tx`.
  stacked <- do.call(rbind, lapply(levels, function(value) {
    profiles <- case$profiles
    profiles$tx <- value
    profiles
  }))
  ind <- inferences(
    sops(case$model, newdata = stacked, times = 1:4, absorb = 3),
    method = "delta"
  )
  ind_jacobian <- get_jacobian(ind)
  ind_key <- paste(ind$rowid, ind$time, ind$state)
  n <- nrow(case$profiles)
  set_jacobian <- function(rowid, time, states) {
    rows <- match(paste(rowid, time, states), ind_key)
    colSums(ind_jacobian[rows, , drop = FALSE])
  }
  expected <- t(vapply(
    seq_len(nrow(cmp)),
    function(i) {
      states <- state_sets[[cmp$state_set[i]]]
      offset <- (match(cmp$comparison_level[i], levels) - 1L) * n
      set_jacobian(cmp$rowid[i] + offset, cmp$time[i], states) -
        set_jacobian(cmp$rowid[i], cmp$time[i], states)
    },
    numeric(ncol(ind_jacobian))
  ))
  V <- attr(ind, "analytical")$coefficient_vcov
  expect_equal(unname(get_jacobian(inferred)), unname(expected), tolerance = 1e-10)
  expect_equal(
    inferred$std.error,
    sqrt(rowSums((expected %*% V) * expected)),
    tolerance = 1e-10
  )
  expect_equal(
    inferred$conf.low,
    inferred$estimate - stats::qnorm(0.975) * inferred$std.error
  )

  block <- vcov(inferred, rows = 1:3)
  expect_equal(
    unname(block),
    unname(expected[1:3, ] %*% V %*% t(expected[1:3, ])),
    tolerance = 1e-10
  )
  expect_equal(unname(sqrt(diag(block))), inferred$std.error[1:3])

  scaled <- inferences(cmp, method = "delta", vcov = 4 * V)
  expect_equal(scaled$std.error, 2 * inferred$std.error, tolerance = 1e-10)

  after_draws <- inferences(
    inferences(cmp, n_draws = 3, seed = 5),
    method = "delta"
  )
  expect_null(attr(after_draws, "draws"))
  expect_equal(after_draws$std.error, inferred$std.error)
})

test_that("averaged delta time-in-state derivatives match avg_comparisons()", {
  case <- make_comparisons_case()
  starts <- factor(c("1", "2"), levels = levels(case$profiles$yprev))
  settings <- list(
    list(variables = list(tx = c(0, 1)), time_map = NULL),
    list(variables = list(tx = c(0, 1)), time_map = case$time_map),
    list(variables = list(yprev = starts), time_map = case$time_map)
  )
  for (setting in settings) {
    args <- list(
      case$model,
      variables = setting$variables,
      estimand = "time_in_state",
      state_sets = list(alive = c("1", "2")),
      times = 1:4,
      absorb = 3,
      time_map = setting$time_map
    )
    unit <- inferences(do.call(comparisons, args), method = "delta")
    average <- inferences(
      do.call(avg_comparisons, args),
      method = "delta",
      vcov = "conditional"
    )
    expect_equal(
      colMeans(get_jacobian(unit)),
      get_jacobian(average)[1, ],
      tolerance = 1e-10
    )
  }
})

test_that("delta comparison inference keeps its restrictions", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    times = 1:4,
    absorb = 3
  )
  cmp <- do.call(comparisons, args)
  expect_error(
    inferences(cmp, method = "delta", vcov = "unconditional"),
    "`comparisons()` delta inference supports only",
    fixed = TRUE
  )
  expect_error(
    inferences(cmp, method = "delta", conf_type = "logit"),
    "conf_type = \"wald\""
  )
  ratio <- do.call(comparisons, c(args, list(comparison = "ratio")))
  expect_error(inferences(ratio, method = "delta"), "differences only")
  benefit <- do.call(comparisons, c(args, list(estimand = "time_benefit")))
  expect_error(inferences(benefit, method = "delta"), "time benefit")

  # Comparisons for each profile are not averaged over patients, so a set
  # starting state keeps the patient-weighted draw methods.
  starts <- factor(c("1", "2"), levels = levels(case$profiles$yprev))
  set_start <- comparisons(
    case$model,
    variables = list(yprev = starts),
    estimand = "time_in_state",
    state_sets = "1",
    times = 1:4,
    absorb = 3,
    time_map = case$time_map
  )
  scored <- inferences(
    set_start,
    method = "score_bootstrap",
    n_draws = 3,
    seed = 21
  )
  expect_true(all(is.finite(scored$std.error)))
})

test_that("delta comparisons handle factor visit times", {
  case <- make_factor_comparisons_case()
  cmp <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:4,
    absorb = 3
  )
  inferred <- inferences(cmp, method = "delta")
  lo <- case$profiles
  lo$tx <- 0
  hi <- case$profiles
  hi$tx <- 1
  ind <- inferences(
    sops(case$model, newdata = rbind(lo, hi), times = 1:4, absorb = 3),
    method = "delta"
  )
  ind_jacobian <- get_jacobian(ind)
  ind_key <- paste(ind$rowid, ind$time, ind$state)
  n <- nrow(case$profiles)
  expected <- ind_jacobian[
    match(paste(cmp$rowid + n, cmp$time, "1"), ind_key),
    ,
    drop = FALSE
  ] -
    ind_jacobian[match(paste(cmp$rowid, cmp$time, "1"), ind_key), , drop = FALSE]
  expect_equal(
    unname(get_jacobian(inferred)),
    unname(expected),
    tolerance = 1e-10
  )

  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    state_sets = "1",
    times = 1:4,
    absorb = 3,
    time_map = case$time_map
  )
  unit <- inferences(do.call(comparisons, args), method = "delta")
  average <- inferences(
    do.call(avg_comparisons, args),
    method = "delta",
    vcov = "conditional"
  )
  expect_equal(
    colMeans(get_jacobian(unit)),
    get_jacobian(average)[1, ],
    tolerance = 1e-10
  )
})
