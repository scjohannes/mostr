test_that("delta intervals difference each profile's SOP derivatives", {
  # Each case: a fit, the compared values of `tx`, and the state sets, given
  # as passed to `comparisons()` (`state_sets`) and as the states of each named
  # set in the result (`members`).
  cases <- list(
    list(
      label = "numeric visits, three treatment levels",
      case = make_comparisons_case(),
      levels = c(0, 1, 2),
      state_sets = list(alive = c("1", "2"), dead = "3"),
      members = list(alive = c("1", "2"), dead = "3")
    ),
    list(
      label = "factor visits",
      case = make_factor_comparisons_case(),
      levels = c(0, 1),
      state_sets = "1",
      members = list("1" = "1")
    )
  )
  for (setting in cases) {
    label <- setting$label
    case <- setting$case
    levels <- setting$levels
    cmp <- comparisons(
      case$model,
      newdata = case$profiles,
      variables = list(tx = levels),
      state_sets = setting$state_sets,
      times = 1:4,
      absorb = 3
    )
    inferred <- inferences(cmp, method = "delta")
    expect_s3_class(inferred, "markov_comparisons")
    expect_identical(attr(inferred, "method"), "delta", info = label)
    expect_identical(attr(inferred, "conf_type"), "wald", info = label)
    expect_identical(attr(inferred, "null"), 0, info = label)
    expect_null(attr(inferred, "draws"), info = label)
    expect_equal(inferred$estimate, cmp$estimate, info = label)

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
        states <- setting$members[[as.character(cmp$state_set[i])]]
        offset <- (match(cmp$comparison_level[i], levels) - 1L) * n
        set_jacobian(cmp$rowid[i] + offset, cmp$time[i], states) -
          set_jacobian(cmp$rowid[i], cmp$time[i], states)
      },
      numeric(ncol(ind_jacobian))
    ))
    V <- attr(ind, "analytical")$coefficient_vcov
    expect_equal(
      unname(get_jacobian(inferred)),
      unname(expected),
      tolerance = 1e-10,
      info = label
    )
    expect_equal(
      inferred$std.error,
      sqrt(rowSums((expected %*% V) * expected)),
      tolerance = 1e-10,
      info = label
    )
    expect_equal(
      inferred$conf.low,
      inferred$estimate - stats::qnorm(0.975) * inferred$std.error,
      info = label
    )

    block <- vcov(inferred, rows = 1:3)
    expect_equal(
      unname(block),
      unname(expected[1:3, ] %*% V %*% t(expected[1:3, ])),
      tolerance = 1e-10,
      info = label
    )
    expect_equal(
      unname(sqrt(diag(block))),
      inferred$std.error[1:3],
      info = label
    )

    scaled <- inferences(cmp, method = "delta", vcov = 4 * V)
    expect_equal(
      scaled$std.error,
      2 * inferred$std.error,
      tolerance = 1e-10,
      info = label
    )

    after_draws <- inferences(
      inferences(cmp, n_draws = 3, seed = 5),
      method = "delta"
    )
    expect_null(attr(after_draws, "draws"), info = label)
    expect_equal(after_draws$std.error, inferred$std.error, info = label)
  }
})

test_that("averaged delta derivatives of each profile match the averages", {
  case <- make_comparisons_case()
  factor_case <- make_factor_comparisons_case()
  starts <- factor(c("1", "2"), levels = levels(case$profiles$yprev))
  alive <- list(alive = c("1", "2"))
  # Each case: the function for each profile and its average-level
  # counterpart, their arguments, the column whose values are averaged
  # separately (`group`, NULL when the average has a single row) and the
  # `null` attribute of the result for each profile.
  comparison_case <- function(label, model, variables, state_sets, time_map) {
    list(
      label = label,
      unit = comparisons,
      average = avg_comparisons,
      args = list(
        model,
        variables = variables,
        estimand = "time_in_state",
        state_sets = state_sets,
        times = 1:4,
        absorb = 3,
        time_map = time_map
      ),
      group = NULL,
      null = 0
    )
  }
  state_time_case <- function(label, args, group) {
    list(
      label = label,
      unit = state_time,
      average = avg_state_time,
      args = c(list(case$model), args, list(times = 1:4, absorb = 3)),
      group = group,
      null = NULL
    )
  }
  cases <- list(
    comparison_case(
      "comparisons of tx, visit times",
      case$model,
      list(tx = c(0, 1)),
      alive,
      NULL
    ),
    comparison_case(
      "comparisons of tx, real time",
      case$model,
      list(tx = c(0, 1)),
      alive,
      case$time_map
    ),
    comparison_case(
      "comparisons of a set starting state, real time",
      case$model,
      list(yprev = starts),
      alive,
      case$time_map
    ),
    comparison_case(
      "comparisons of tx, factor visits, real time",
      factor_case$model,
      list(tx = c(0, 1)),
      "1",
      factor_case$time_map
    ),
    state_time_case(
      "state_time() by tx, visit times",
      list(variables = list(tx = c(0, 1)), state_sets = alive, time_map = NULL),
      "tx"
    ),
    state_time_case(
      "state_time() by tx, real time",
      list(
        variables = list(tx = c(0, 1)),
        state_sets = alive,
        time_map = case$time_map
      ),
      "tx"
    ),
    state_time_case(
      "state_time() without variables",
      list(state_sets = "1"),
      NULL
    )
  )
  for (setting in cases) {
    label <- setting$label
    unit <- inferences(do.call(setting$unit, setting$args), method = "delta")
    average <- inferences(
      do.call(setting$average, setting$args),
      method = "delta",
      vcov = "conditional"
    )
    expect_identical(attr(unit, "method"), "delta", info = label)
    expect_identical(attr(unit, "conf_type"), "wald", info = label)
    expect_identical(attr(unit, "null"), setting$null, info = label)
    unit_jacobian <- get_jacobian(unit)
    average_jacobian <- get_jacobian(average)
    if (is.null(setting$group)) {
      expect_equal(
        colMeans(unit_jacobian),
        average_jacobian[1, ],
        tolerance = 1e-10,
        info = label
      )
    } else {
      group <- setting$group
      for (value in unique(unit[[group]])) {
        expect_equal(
          colMeans(unit_jacobian[unit[[group]] == value, , drop = FALSE]),
          average_jacobian[average[[group]] == value, ],
          tolerance = 1e-10,
          info = paste(label, value)
        )
      }
    }
    V <- attr(unit, "analytical")$coefficient_vcov
    expect_equal(
      unit$std.error,
      unname(sqrt(rowSums((unit_jacobian %*% V) * unit_jacobian))),
      tolerance = 1e-10,
      info = label
    )
    block <- vcov(unit, rows = 1:2)
    expect_equal(
      unname(sqrt(diag(block))),
      unit$std.error[1:2],
      info = label
    )
  }

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
    paste(
      "Delta inference for `comparisons()` and `avg_comparisons(by = FALSE)`",
      "supports only `vcov = \"conditional\"`"
    ),
    fixed = TRUE
  )

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
  scored <- mostr:::inferences_score_bootstrap(
    set_start,
    n_draws = 3,
    seed = 21
  )
  expect_true(all(is.finite(scored$std.error)))
})
