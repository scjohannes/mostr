# Checks of inference results that need no known answer: limits that hold for
# any data, and two routes to one number. test-grid-inference.R runs them on
# the full crossing of engine, time coding and method, and
# test-grid-inference-array.R on a pairwise covering array of the other
# settings. The IDs (L9, R13, N4, ...) refer to the catalogue of claims behind
# the test overhaul.
#
# Draws are few (MVN 20, bootstrap and fractional weighted bootstrap 3), so
# the checks are about how results are put together, not about whether the
# intervals cover. None of them depends on the seed. Deliberately absent, as
# they are false: percentile intervals need not contain the estimate,
# unconditional standard errors need not exceed conditional ones, and
# weighted means of draws reproduce averaged draws for differences only.

inference_seed <- 11L

inference_n_draws <- function(method) {
  if (identical(method, "mvn")) 20L else 3L
}

inference_is_draw_method <- function(method) {
  method %in% c("mvn", "bootstrap", "fwb")
}

# One inference set-up: what is estimated, by which method, on which fit.
# Settings that do not apply are "-" (`by` for results with one row per
# profile, `scale` for SOPs, `comparison` for anything but comparisons).
inference_setting <- function(
  result = "avg_sops",
  method = "mvn",
  engine = "orm",
  time = "linear",
  order = "1",
  odds = "full",
  weights = "none",
  profiles = "stored",
  by = "none",
  scale = "-",
  comparison = "-",
  vcov = NULL
) {
  setting <- list(
    result = result,
    method = method,
    engine = engine,
    time = time,
    order = order,
    odds = odds,
    weights = weights,
    profiles = profiles,
    by = by,
    scale = scale,
    comparison = comparison
  )
  # Unconditional variance needs the stored patients of an unweighted fit and
  # an average; otherwise analytical inference uses conditional variance.
  if (method == "delta" && is.null(vcov)) {
    vcov <- if (
      startsWith(result, "avg_") && profiles == "stored" && weights == "none"
    ) {
      "unconditional"
    } else {
      "conditional"
    }
  }
  setting$vcov <- vcov
  setting
}

inference_label <- function(setting) {
  shown <- setting[setting != "-"]
  paste(names(shown), unlist(shown), sep = " = ", collapse = ", ")
}

inference_cell <- function(setting) {
  grid_cell(
    setting$engine,
    order = as.integer(setting$order),
    time = setting$time,
    weighted = setting$weights == "integer",
    odds = setting$odds
  )
}

# Twelve stored starting profiles in `newdata`, from both `marker` groups.
inference_profiles <- function(cell) {
  data <- grid_cell_data(cell)
  start <- data[as.character(data$time) == "1", , drop = FALSE]
  start <- start[seq_len(12), , drop = FALSE]
  rownames(start) <- NULL
  start
}

inference_is_sop <- function(setting) {
  setting$result %in% c("sops", "avg_sops")
}

inference_is_time <- function(setting) {
  setting$result %in% c("state_time", "avg_state_time")
}

inference_is_comparison <- function(setting) {
  setting$result %in% c("comparisons", "avg_comparisons")
}

inference_estimand <- function(setting) {
  if (identical(setting$comparison, "time_benefit")) {
    return("time_benefit")
  }
  if (identical(setting$scale, "time_map")) "time_in_state" else "sop"
}

# Arguments of the public call for a setting, apart from the fit and the
# grid's visits and absorbing state.
inference_call_args <- function(setting, variables = grid_treatments) {
  cell <- inference_cell(setting)
  args <- list(variables = variables)
  if (setting$profiles == "newdata") {
    args$newdata <- inference_profiles(cell)
  }
  if (setting$by == "marker") {
    args$by <- "marker"
  }
  if (setting$scale == "time_map") {
    args$time_map <- grid_time_map(cell)
  }
  if (inference_is_comparison(setting)) {
    args$estimand <- inference_estimand(setting)
    if (setting$comparison == "ratio") {
      args$comparison <- "ratio"
    }
  }
  args
}

inference_point <- function(setting, variables = grid_treatments) {
  do.call(
    grid_call,
    c(
      list(get(setting$result), inference_cell(setting)),
      inference_call_args(setting, variables)
    )
  )
}

# The averaged comparison of treatment arms that matches an averaged SOP or
# time result: same profiles, groups and time scale.
inference_matching_contrast <- function(setting) {
  contrast <- setting
  contrast$result <- "avg_comparisons"
  contrast$comparison <- "difference"
  if (inference_is_time(setting)) {
    args <- inference_call_args(contrast)
    args$estimand <- "time_in_state"
    return(do.call(
      grid_call,
      c(list(avg_comparisons, inference_cell(setting)), args)
    ))
  }
  contrast$scale <- "visit"
  inference_point(contrast)
}

# Inference on `x` by the setting's method. Draw methods use a fixed seed.
# `...` takes further arguments of inferences(). The warning that supplied
# `newdata` profiles get no bootstrap draw weights is expected here.
inference_run <- function(x, setting, ...) {
  if (setting$method == "delta") {
    return(inferences(x, method = "delta", vcov = setting$vcov, ...))
  }
  withCallingHandlers(
    inferences(
      x,
      method = setting$method,
      n_draws = inference_n_draws(setting$method),
      seed = inference_seed,
      ...
    ),
    warning = function(w) {
      if (grepl("treated as fixed", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Columns that identify a row of a result or of its draws.
inference_keys <- function(x) {
  intersect(
    c("rowid", "tx", "marker", "estimand", "time", "state", "state_set"),
    names(x)
  )
}

inference_draws <- function(x) {
  as.data.frame(get_draws(x))
}

# The value a Wald test is run against: "no effect" for comparisons, and a
# value inside the range of the estimates otherwise.
inference_null <- function(setting) {
  if (inference_is_comparison(setting)) {
    return("auto")
  }
  if (inference_is_sop(setting)) 0.25 else 1
}

# Rows whose summaries across all states are fixed: SOPs add to one at each
# visit, times in all states add to the follow-up, and differences of those
# add to zero. Returns the grouping of rows into such partitions, or NULL.
inference_partition_groups <- function(x, setting) {
  if (inference_is_comparison(setting)) {
    if (setting$comparison != "difference") {
      return(NULL)
    }
  }
  member <- intersect(c("state", "state_set"), names(x))
  if (length(member) != 1L) {
    return(NULL)
  }
  states <- sort(unique(as.character(x[[member]])))
  if (!identical(states, as.character(seq_len(grid_n_states)))) {
    return(NULL)
  }
  grid_key(x, setdiff(inference_keys(x), member))
}

# For results with one row per profile, the analytical covariance of every
# pair of rows is large, so the covariance checks use the rows of the first
# three profiles.
inference_vcov_rows <- function(x) {
  if (!"rowid" %in% names(x)) {
    return(seq_len(nrow(x)))
  }
  which(x$rowid %in% sort(unique(x$rowid))[1:3])
}

# Runs every check that applies to `setting`. Each run uses the same seed:
# A at the 95% level with draws kept, B at the 80% level, and C with Wald
# intervals and a null test, without draws.
check_inference_setting <- function(setting) {
  label <- inference_label(setting)
  draw_method <- inference_is_draw_method(setting$method)
  point <- inference_point(setting)
  keys <- inference_keys(point)
  null <- inference_null(setting)

  # N4: the caller's random number stream is left as it was.
  random_state <- function() get(".Random.seed", envir = globalenv())
  set.seed(99)
  before <- random_state()
  a <- inference_run(point, setting, null = null)
  expect_identical(random_state(), before, label = paste(label, "N4 random state"))

  b <- inference_run(point, setting, conf_level = 0.8, null = null)
  c <- inference_run(
    point,
    setting,
    conf_type = "wald",
    return_draws = FALSE,
    null = null
  )
  a_df <- as.data.frame(a)
  b_df <- as.data.frame(b)
  c_df <- as.data.frame(c)

  # L9 and L13.
  for (run in list(a = a_df, b = b_df, c = c_df)) {
    expect_false(anyNA(run$std.error), label = paste(label, "L9 complete"))
    expect_gte(min(run$std.error), 0, label = paste(label, "L9"))
    expect_true(
      all(run$conf.low <= run$conf.high + 1e-12),
      label = paste(label, "L13")
    )
  }

  # N21: neither the level, the interval type nor keeping draws changes the
  # estimates or standard errors.
  for (run in list(b = b_df, c = c_df)) {
    expect_lt(grid_gap(run, a_df, keys), 1e-14, label = paste(label, "N21 estimate"))
    expect_lt(
      grid_gap(run, a_df, keys, "std.error"),
      1e-14,
      label = paste(label, "N21 standard error")
    )
  }
  expect_null(attr(c, "draws"), label = paste(label, "N21 no draws"))

  # L10: Wald intervals are centred on the estimate, z standard errors wide.
  z <- stats::qnorm(0.975)
  expect_lt(
    max(abs((c_df$conf.low + c_df$conf.high) / 2 - c_df$estimate)),
    1e-12,
    label = paste(label, "L10 centre")
  )
  expect_lt(
    max(abs(c_df$conf.high - c_df$estimate - z * c_df$std.error)),
    1e-12,
    label = paste(label, "L10 width")
  )

  # L24: the Wald interval excludes the null exactly when p < 0.05.
  if ("p.value" %in% names(c_df)) {
    null_value <- if (identical(null, "auto")) {
      if (setting$comparison == "ratio") 1 else 0
    } else {
      null
    }
    tested <- c_df$std.error > 0 & !is.na(c_df$p.value)
    excluded <- c_df$conf.low > null_value | c_df$conf.high < null_value
    expect_identical(
      excluded[tested],
      c_df$p.value[tested] < 0.05,
      label = paste(label, "L24")
    )
    expect_true(any(tested), label = paste(label, "L24 has tests"))
  } else {
    expect_true(FALSE, label = paste(label, "L24 null test missing"))
  }

  # L23: the 80% interval lies within the 95% interval. For percentile
  # intervals this needs the same draws, which the same seed gives.
  expect_true(
    all(b_df$conf.low >= a_df$conf.low[match(grid_key(b_df, keys), grid_key(a_df, keys))] - 1e-12),
    label = paste(label, "L23 lower")
  )
  expect_true(
    all(b_df$conf.high <= a_df$conf.high[match(grid_key(b_df, keys), grid_key(a_df, keys))] + 1e-12),
    label = paste(label, "L23 upper")
  )

  # L11: logit limits of analytical SOPs stay in [0, 1].
  if (setting$method == "delta" && inference_is_sop(setting)) {
    expect_gte(min(a_df$conf.low), 0, label = paste(label, "L11 lower"))
    expect_lte(max(a_df$conf.high), 1, label = paste(label, "L11 upper"))
  }

  # L18: percentile limits of an absorbing state never fall over visits,
  # because each draw's occupancy of it never falls.
  if (draw_method && inference_is_sop(setting)) {
    groups <- setdiff(keys, c("time", "state"))
    times <- grid_times(inference_cell(setting))
    for (limit in c("conf.low", "conf.high")) {
      limits <- a_df
      limits$estimate <- limits[[limit]]
      expect_lt(
        grid_largest_fall(limits, grid_n_states, groups, times),
        1e-12,
        label = paste(label, "L18", limit)
      )
    }
  }

  partition <- inference_partition_groups(a_df, setting)
  if (draw_method) {
    draws <- inference_draws(a)

    # L29: bootstrap bookkeeping.
    n_ok <- attr(a, "n_successful")
    expect_true(
      !is.null(n_ok) && n_ok <= inference_n_draws(setting$method) && n_ok >= 2L,
      label = paste(label, "L29 successful draws")
    )
    expect_identical(length(unique(draws$draw_id)), as.integer(n_ok), label = paste(label, "L29 draw IDs"))
    expect_identical(nrow(draws), as.integer(n_ok) * nrow(a_df), label = paste(label, "L29 draw rows"))

    # N4: the same seed gives the same draws, whatever the level.
    # The stored draw tables are compared directly; get_draws() adds the
    # summary columns, which differ with the level.
    stored <- function(x) {
      table <- attr(x, "draws")
      if ("draw" %in% names(table)) table$draw else table$estimate
    }
    expect_identical(stored(b), stored(a), label = paste(label, "N4 draws"))
    expect_identical(length(stored(a)), nrow(draws), label = paste(label, "N4 has draws"))

    # R14 for draws: each standard error is the SD of that row's draws.
    spread <- tapply(draws$draw, grid_key(draws, keys), stats::sd)
    expect_identical(length(spread), nrow(a_df), label = paste(label, "R14 rows"))
    expect_lt(
      max(abs(spread[grid_key(a_df, keys)] - a_df$std.error)),
      1e-12,
      label = paste(label, "R14 draws")
    )

    # L21 for draws: every draw's partition adds to the same total.
    if (!is.null(partition)) {
      cell <- grid_key(draws, setdiff(keys, c("state", "state_set")))
      totals <- tapply(draws$draw, paste(cell, draws$draw_id), sum)
      cell_of_total <- cell[match(names(totals), paste(cell, draws$draw_id))]
      spread <- tapply(totals, cell_of_total, stats::sd)
      expect_lt(max(spread), 1e-10, label = paste(label, "L21 draws"))
    }

    # R13 for draws: the treatment contrast's standard error is the SD of
    # the difference of the two arms' draws, under the same seed.
    if (setting$result %in% c("avg_sops", "avg_state_time")) {
      contrast <- inference_run(inference_matching_contrast(setting), setting)
      contrast_df <- as.data.frame(contrast)
      member <- intersect(c("state", "state_set"), names(draws))
      arm_keys <- setdiff(keys, c("tx", member))
      treated <- draws[draws$tx == 1, , drop = FALSE]
      control <- draws[draws$tx == 0, , drop = FALSE]
      match_keys <- c(arm_keys, member, "draw_id")
      treated$draw <- treated$draw -
        control$draw[match(grid_key(treated, match_keys), grid_key(control, match_keys))]
      treated$state_set <- as.character(treated[[member]])
      spread <- stats::aggregate(
        treated["draw"],
        treated[c(arm_keys, "state_set")],
        stats::sd
      )
      names(spread)[names(spread) == "draw"] <- "std.error"
      expect_lt(
        grid_gap(contrast_df, spread, c(arm_keys, "state_set"), "std.error"),
        1e-12,
        label = paste(label, "R13 draws")
      )
    }
  } else {
    # R14 and L21 from the analytical covariance.
    rows <- inference_vcov_rows(a_df)
    covariance <- stats::vcov(a, rows = rows)
    expect_lt(
      max(abs(sqrt(diag(covariance)) - a_df$std.error[rows])),
      1e-12,
      label = paste(label, "R14")
    )
    if (!is.null(partition)) {
      groups <- partition[rows]
      sums <- vapply(unique(groups), function(group) {
        w <- as.numeric(groups == group)
        drop(t(w) %*% covariance %*% w)
      }, numeric(1))
      expect_lt(
        max(abs(sums)),
        1e-12 * max(1, max(diag(covariance))),
        label = paste(label, "L21")
      )
    }
  }

  # N1: swapping the compared values flips the contrast, its standard error
  # and its limits; a ratio becomes its reciprocal.
  if (inference_is_comparison(setting)) {
    swapped <- inference_run(inference_point(setting, list(tx = c(1, 0))), setting, null = null)
    swapped <- as.data.frame(swapped)
    swap_keys <- setdiff(keys, "tx")
    if (setting$comparison == "ratio") {
      swapped$estimate <- 1 / swapped$estimate
      expect_lt(grid_gap(a_df, swapped, swap_keys), 1e-10, label = paste(label, "N1 ratio"))
    } else {
      flipped <- swapped
      flipped$estimate <- -swapped$estimate
      flipped$conf.low <- -swapped$conf.high
      flipped$conf.high <- -swapped$conf.low
      for (column in c("estimate", "std.error", "conf.low", "conf.high")) {
        expect_lt(
          grid_gap(a_df, flipped, swap_keys, column),
          1e-12,
          label = paste(label, "N1", column)
        )
      }
    }
  }
  invisible(NULL)
}
