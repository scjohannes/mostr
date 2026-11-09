# `by` that names the compared variable groups patients by the value they were
# observed with, as the `comparisons()` help promises ("results can be
# summarized, for example, among treated patients"). Each group's average must
# equal the same call with `newdata` restricted to that group's profiles.

# Key columns that identify one estimate, besides the `by` group.
by_estimate_key <- function(x, cols) {
  do.call(paste, c(unname(as.list(x[cols])), sep = "|"))
}

# `avg_comparisons()` within each observed treatment group, from `newdata`
# restricted to that group, stacked with the group's value in `tx`.
comparisons_within_observed_tx <- function(case, args) {
  out <- lapply(c(0, 1), function(value) {
    profiles <- case$profiles[case$profiles$tx == value, , drop = FALSE]
    res <- do.call(avg_comparisons, c(list(case$model, newdata = profiles), args))
    res$tx <- value
    res
  })
  do.call(rbind, lapply(out, as.data.frame))
}

expect_matches_observed_groups <- function(out, expected, columns, info) {
  expect_true("tx" %in% names(out), info = info)
  if (!"tx" %in% names(out)) {
    return(invisible())
  }
  expect_setequal(unique(out$tx), c(0, 1))
  expect_identical(nrow(out), nrow(expected), info = info)
  cols <- intersect(
    c("time", "state_set", "state", "contrast"),
    intersect(names(out), names(expected))
  )
  key_out <- paste(out$tx, by_estimate_key(out, cols))
  key_expected <- paste(expected$tx, by_estimate_key(expected, cols))
  index <- match(key_expected, key_out)
  expect_false(anyNA(index), info = info)
  for (column in columns) {
    expect_equal(
      out[[column]][index],
      expected[[column]],
      tolerance = 1e-10,
      info = paste(info, column)
    )
  }
}

comparison_estimand_args <- function(estimand) {
  args <- list(
    variables = list(tx = c(0, 1)),
    estimand = estimand,
    times = 1:4,
    absorb = 3
  )
  if (estimand != "time_benefit") {
    args$state_sets <- "1"
  }
  args
}

test_that("avg_comparisons(by = compared variable) averages within observed groups", {
  case <- make_comparisons_case()
  for (estimand in c("sop", "time_in_state", "time_benefit")) {
    args <- comparison_estimand_args(estimand)
    out <- do.call(
      avg_comparisons,
      c(list(case$model, newdata = case$profiles, by = "tx"), args)
    )
    expected <- comparisons_within_observed_tx(case, args)
    expect_matches_observed_groups(out, expected, "estimate", estimand)

    # The same as the mean of the per-profile comparisons in each group.
    each <- do.call(comparisons, c(list(case$model, newdata = case$profiles), args))
    keys <- intersect(c("time", "state_set"), names(each))
    means <- stats::aggregate(
      each["estimate"],
      c(each[keys], list(tx = each$tx)),
      mean
    )
    expect_matches_observed_groups(out, means, "estimate", paste(estimand, "mean"))
  }
})

test_that("time in state with time_map by the compared variable uses each group's own start", {
  case <- make_comparisons_case()
  args <- c(comparison_estimand_args("time_in_state"), list(time_map = case$time_map))
  out <- do.call(
    avg_comparisons,
    c(list(case$model, newdata = case$profiles, by = "tx"), args)
  )
  expected <- comparisons_within_observed_tx(case, args)
  expect_matches_observed_groups(out, expected, "estimate", "time_map")
})

test_that("draw-based intervals by the compared variable match each group's own call", {
  case <- make_comparisons_case()
  for (estimand in c("sop", "time_in_state", "time_benefit")) {
    args <- comparison_estimand_args(estimand)
    out <- inferences(
      do.call(
        avg_comparisons,
        c(list(case$model, newdata = case$profiles, by = "tx"), args)
      ),
      method = "mvn",
      n_draws = 50,
      seed = 11
    )
    expected <- lapply(c(0, 1), function(value) {
      profiles <- case$profiles[case$profiles$tx == value, , drop = FALSE]
      res <- inferences(
        do.call(avg_comparisons, c(list(case$model, newdata = profiles), args)),
        method = "mvn",
        n_draws = 50,
        seed = 11
      )
      res <- as.data.frame(res)
      res$tx <- value
      res
    })
    expected <- do.call(rbind, expected)
    expect_matches_observed_groups(
      out,
      expected,
      c("estimate", "std.error", "conf.low", "conf.high"),
      estimand
    )
  }
})

test_that("the delta method still refuses `by`, including the compared variable", {
  case <- make_comparisons_case()
  averaged <- do.call(
    avg_comparisons,
    c(list(case$model, newdata = case$profiles, by = "tx"), comparison_estimand_args("sop"))
  )
  expect_error(
    inferences(averaged, method = "delta", vcov = "conditional"),
    "Analytical delta inference with `by` is not yet supported.",
    fixed = TRUE
  )
})

test_that("avg_sops() refuses `by` on the variable it sets, and a copied column works", {
  case <- make_comparisons_case()
  args <- list(times = 1:4, absorb = 3, variables = list(tx = c(0, 1)))
  expect_error(
    do.call(avg_sops, c(list(case$model, newdata = case$profiles, by = "tx"), args)),
    "copy",
    info = "the message explains the copied-column workaround"
  )

  # The workaround: group by a copy of the observed treatment.
  profiles <- case$profiles
  profiles$observed_tx <- profiles$tx
  out <- do.call(
    avg_sops,
    c(list(case$model, newdata = profiles, by = "observed_tx"), args)
  )
  for (value in c(0, 1)) {
    expected <- do.call(
      avg_sops,
      c(list(case$model, newdata = case$profiles[case$profiles$tx == value, ]), args)
    )
    got <- out[out$observed_tx == value, , drop = FALSE]
    key <- function(x) paste(x$tx, x$time, x$state)
    expect_equal(got$estimate[match(key(expected), key(got))], expected$estimate)
  }
})

test_that("avg_sops() by treatment without `variables` still groups by observed value", {
  case <- make_comparisons_case()
  out <- avg_sops(case$model, newdata = case$profiles, times = 1:4, absorb = 3, by = "tx")
  for (value in c(0, 1)) {
    expected <- avg_sops(
      case$model,
      newdata = case$profiles[case$profiles$tx == value, ],
      times = 1:4,
      absorb = 3
    )
    got <- out[out$tx == value, , drop = FALSE]
    key <- function(x) paste(x$time, x$state)
    expect_equal(got$estimate[match(key(expected), key(got))], expected$estimate)
  }
})

test_that("`by` with the compared variable and another column, and ratios, use observed groups", {
  case <- make_comparisons_case()
  for (estimand in c("sop", "time_in_state")) {
    args <- c(comparison_estimand_args(estimand), list(comparison = "ratio"))
    out <- do.call(
      avg_comparisons,
      c(list(case$model, newdata = case$profiles, by = c("tx", "marker")), args)
    )
    expect_true(all(c("tx", "marker") %in% names(out)), info = estimand)
    for (tx in c(0, 1)) {
      for (marker in c(0, 1)) {
        keep <- case$profiles$tx == tx & case$profiles$marker == marker
        expected <- do.call(
          avg_comparisons,
          c(list(case$model, newdata = case$profiles[keep, , drop = FALSE]), args)
        )
        got <- out[out$tx == tx & out$marker == marker, , drop = FALSE]
        cols <- intersect(c("time", "state_set"), names(expected))
        index <- match(by_estimate_key(expected, cols), by_estimate_key(got, cols))
        expect_false(anyNA(index), info = estimand)
        # A ratio of the group's averaged summaries, not an average of ratios.
        expect_equal(
          got$estimate[index],
          expected$estimate,
          tolerance = 1e-10,
          info = paste(estimand, tx, marker)
        )
      }
    }
  }
})

test_that("resampling by the compared variable matches grouping by a copied column", {
  case <- make_comparisons_case()
  # The reference: the same model fitted to data with a copy of the observed
  # treatment, grouped by that copy as an ordinary column.
  data <- case$data
  data$observed_tx <- data$tx
  copied <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
  )
  args <- c(
    comparison_estimand_args("time_in_state"),
    list(time_map = case$time_map)
  )
  for (method in c("fwb", "bootstrap")) {
    out <- inferences(
      do.call(avg_comparisons, c(list(case$model, by = "tx"), args)),
      method = method,
      n_draws = 8,
      seed = 5
    )
    expected <- as.data.frame(inferences(
      do.call(avg_comparisons, c(list(copied, by = "observed_tx"), args)),
      method = method,
      n_draws = 8,
      seed = 5
    ))
    names(expected)[names(expected) == "observed_tx"] <- "tx"
    expect_matches_observed_groups(
      out,
      expected,
      c("estimate", "std.error", "conf.low", "conf.high"),
      method
    )
  }
})

test_that("avg_state_time() refuses `by` on the variable it sets", {
  case <- make_comparisons_case()
  expect_error(
    avg_state_time(
      case$model,
      newdata = case$profiles,
      variables = list(tx = c(0, 1)),
      by = "tx",
      times = 1:4,
      absorb = 3
    ),
    "copy"
  )
  # Without `variables`, `by = "tx"` groups by the observed treatment.
  out <- avg_state_time(
    case$model,
    newdata = case$profiles,
    by = "tx",
    times = 1:4,
    absorb = 3
  )
  expect_setequal(unique(out$tx), c(0, 1))
})
