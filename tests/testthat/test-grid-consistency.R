# Consistency checks on every model set-up of the shared grid: limits, two
# routes to one number, changes that should not matter, and orm against vglm.
# Each check computes one number in two ways, or checks a limit that holds
# whatever the model, so it needs no known answer. The IDs in the test names
# (L1, R1, N1, ...) refer to the catalogue of claims behind the test overhaul.
# Time totals, comparisons and time benefit are in
# test-grid-consistency-time.R; the grid and shared helpers in helper-grid.R.
#
# On CRAN only one set-up per engine and Markov order runs
# (`grid_consistency_cells()`); everything here is deterministic.

test_that("SOPs lie in [0, 1] and add to one at every visit (L1, L2)", {
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    average <- grid_call(avg_sops, cell, variables = grid_treatments)
    results <- list(
      stored = list(grid_call(sops, cell), c("rowid", "time")),
      average = list(average, c("tx", "time")),
      inner_absorbing = list(
        grid_call(
          sops,
          grid_inner_absorbing(cell),
          newdata = grid_profiles(grid_inner_absorbing(cell))
        ),
        c("rowid", "time")
      ),
      real_time = list(
        interpolate_sops(average, time_map = grid_time_map(cell)),
        c("tx", "time")
      )
    )
    for (name in names(results)) {
      x <- results[[name]][[1]]
      expect_gte(min(x$estimate), -1e-12, label = paste(label, name))
      expect_lte(max(x$estimate), 1 + 1e-12, label = paste(label, name))
      totals <- tapply(x$estimate, grid_key(x, results[[name]][[2]]), sum)
      expect_lt(max(abs(totals - 1)), 1e-10, label = paste(label, name))
    }
  }
})

test_that("occupancy of an absorbing state never falls (L4)", {
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    times <- grid_times(cell)
    average <- grid_call(avg_sops, cell, variables = grid_treatments)
    expect_lt(
      grid_largest_fall(grid_call(sops, cell), grid_n_states, "rowid", times),
      1e-12,
      label = paste(label, "stored")
    )
    expect_lt(
      grid_largest_fall(average, grid_n_states, "tx", times),
      1e-12,
      label = paste(label, "average")
    )
    real <- interpolate_sops(average, time_map = grid_time_map(cell))
    expect_lt(
      grid_largest_fall(real, grid_n_states, "tx", sort(unique(real$time))),
      1e-12,
      label = paste(label, "real time from baseline")
    )
    inner <- grid_inner_absorbing(cell)
    several <- grid_call(sops, inner, newdata = grid_profiles(inner))
    for (state in inner$absorb) {
      expect_lt(
        grid_largest_fall(several, state, "rowid", times),
        1e-12,
        label = paste(label, "absorbing state", state)
      )
    }
  }
})

test_that("avg_sops() is the mean of sops(), overall and within groups (R1)", {
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    each <- as.data.frame(grid_call(sops, cell, variables = grid_treatments))
    for (by in list(TRUE, "marker")) {
      groups <- c("time", "state", "tx", if (is.character(by)) by)
      average <- grid_call(avg_sops, cell, variables = grid_treatments, by = by)
      means <- stats::aggregate(each["estimate"], each[groups], mean)
      expect_lt(
        grid_gap(average, means, groups),
        1e-10,
        label = paste(label, "by", by)
      )
    }
  }
})

test_that("contrasts are differences of the averaged summaries (R3)", {
  skip_on_cran()
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    average <- as.data.frame(grid_call(avg_sops, cell, variables = grid_treatments))
    treated <- average[average$tx == 1, , drop = FALSE]
    control <- average[average$tx == 0, , drop = FALSE]
    keys <- c("time", "state")
    expected <- data.frame(
      time = treated$time,
      state_set = as.character(treated$state),
      estimate = treated$estimate -
        control$estimate[match(grid_key(treated, keys), grid_key(control, keys))]
    )
    contrast <- grid_call(
      avg_comparisons,
      cell,
      variables = grid_treatments,
      estimand = "sop"
    )
    expect_lt(
      grid_gap(contrast, expected, c("time", "state_set")),
      1e-12,
      label = paste(label, "sop")
    )

    for (extra in list(list(), list(time_map = grid_time_map(cell)))) {
      totals <- as.data.frame(do.call(
        grid_call,
        c(list(avg_state_time, cell, variables = grid_treatments), extra)
      ))
      treated <- totals[totals$tx == 1, , drop = FALSE]
      control <- totals[totals$tx == 0, , drop = FALSE]
      expected <- data.frame(
        state_set = treated$state_set,
        estimate = treated$estimate -
          control$estimate[match(treated$state_set, control$state_set)]
      )
      contrast <- do.call(
        grid_call,
        c(
          list(avg_comparisons, cell, variables = grid_treatments, estimand = "time_in_state"),
          extra
        )
      )
      expect_lt(
        grid_gap(contrast, expected, "state_set"),
        1e-12,
        label = paste(label, "time in state", names(extra))
      )
    }
  }
})

test_that("a state set is the sum of its members, with matching variance (R6)", {
  skip_on_cran()
  sets <- list(low = c("1", "2"), one = "1", two = "2")
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    for (estimand in c("sop", "time_in_state")) {
      contrast <- as.data.frame(grid_call(
        avg_comparisons,
        cell,
        variables = grid_treatments,
        estimand = estimand,
        state_sets = sets
      ))
      keys <- intersect("time", names(contrast))
      low <- contrast[contrast$state_set == "low", , drop = FALSE]
      one <- contrast[contrast$state_set == "one", , drop = FALSE]
      two <- contrast[contrast$state_set == "two", , drop = FALSE]
      one$estimate <- one$estimate +
        two$estimate[match(grid_key(one, keys), grid_key(two, keys))]
      expect_lt(grid_gap(low, one, keys), 1e-12, label = paste(label, estimand))
    }

    time_map <- grid_time_map(cell)
    members <- grid_call(
      avg_state_time,
      cell,
      variables = grid_treatments,
      time_map = time_map
    )
    lumped <- grid_call(
      avg_state_time,
      cell,
      variables = grid_treatments,
      time_map = time_map,
      state_sets = list(low = c("1", "2"))
    )
    in_set <- as.data.frame(members)[members$state_set %in% c("1", "2"), , drop = FALSE]
    expected <- stats::aggregate(in_set["estimate"], in_set["tx"], sum)
    expect_lt(grid_gap(lumped, expected, "tx"), 1e-12, label = paste(label, "time"))

    if (!grid_has_delta(cell)) next
    for (vcov in grid_vcov_types(cell)) {
      members_inf <- grid_delta(members, vcov)
      lumped_inf <- grid_delta(lumped, vcov)
      covariance <- stats::vcov(members_inf)
      for (tx in c(0, 1)) {
        w <- as.numeric(members_inf$state_set %in% c("1", "2") & members_inf$tx == tx)
        expect_equal(
          lumped_inf$std.error[lumped_inf$tx == tx]^2,
          drop(t(w) %*% covariance %*% w),
          tolerance = 1e-8,
          label = paste(label, vcov, "variance of the set, tx", tx)
        )
      }
    }
  }
})

test_that("interpolated SOPs pass through the visits and the baseline distribution (R9)", {
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    time_map <- grid_time_map(cell)
    result <- grid_call(avg_sops, cell, variables = grid_treatments)
    real <- interpolate_sops(result, time_map = time_map)
    average <- as.data.frame(result)
    average$time <- unname(time_map[as.character(average$time)])

    at_visits <- real[real$time > 0, , drop = FALSE]
    expect_lt(
      grid_gap(at_visits, average, c("time", "state", "tx")),
      1e-12,
      label = paste(label, "visits")
    )
    start <- grid_start_distribution(cell)
    at_baseline <- real[real$time == 0, , drop = FALSE]
    expect_identical(nrow(at_baseline), 2L * grid_n_states, label = label)
    expect_lt(
      max(abs(at_baseline$estimate - start[as.character(at_baseline$state)])),
      1e-12,
      label = paste(label, "baseline")
    )

    # Halfway between the first two visits: the mean of their SOPs.
    between <- interpolate_sops(
      result,
      time_map = time_map,
      target_times = mean(time_map[1:2])
    )
    first <- average[average$time == time_map[[1]], , drop = FALSE]
    second <- average[average$time == time_map[[2]], , drop = FALSE]
    keys <- c("state", "tx")
    first$estimate <- (first$estimate +
      second$estimate[match(grid_key(first, keys), grid_key(second, keys))]) / 2
    expect_lt(grid_gap(between, first, keys), 1e-12, label = paste(label, "halfway"))
  }
})

test_that("analytical covariances give the standard errors and contrast variances (R13, R14)", {
  skip_on_cran()
  for (cell in Filter(grid_has_delta, grid_consistency_cells())) {
    label <- grid_label(cell)
    time_map <- grid_time_map(cell)
    for (vcov in grid_vcov_types(cell)) {
      what <- paste(label, vcov)
      sop <- grid_delta(grid_call(avg_sops, cell, variables = grid_treatments), vcov)
      totals <- grid_delta(
        grid_call(avg_state_time, cell, variables = grid_treatments, time_map = time_map),
        vcov
      )
      sop_contrast <- grid_delta(
        grid_call(avg_comparisons, cell, variables = grid_treatments, estimand = "sop"),
        vcov
      )
      time_contrast <- grid_delta(
        grid_call(
          avg_comparisons,
          cell,
          variables = grid_treatments,
          estimand = "time_in_state",
          time_map = time_map
        ),
        vcov
      )
      for (x in list(sop, totals, sop_contrast, time_contrast)) {
        expect_lt(
          max(abs(sqrt(diag(stats::vcov(x))) - x$std.error)),
          1e-12,
          label = paste(what, "R14", class(x)[1])
        )
      }

      # Var(a - b) = Var(a) + Var(b) - 2 Cov(a, b).
      contrast_se <- function(summary, contrast, keys) {
        covariance <- stats::vcov(summary)
        key_summary <- grid_key(summary, keys)
        vapply(seq_len(nrow(contrast)), function(i) {
          key <- grid_key(contrast[i, , drop = FALSE], keys)
          hi <- which(key_summary == key & summary$tx == 1)
          lo <- which(key_summary == key & summary$tx == 0)
          sqrt(covariance[hi, hi] + covariance[lo, lo] - 2 * covariance[hi, lo])
        }, numeric(1))
      }
      sop$state_set <- as.character(sop$state)
      expect_lt(
        max(abs(
          contrast_se(sop, sop_contrast, c("time", "state_set")) -
            sop_contrast$std.error
        )),
        1e-10,
        label = paste(what, "R13 sop")
      )
      expect_lt(
        max(abs(
          contrast_se(totals, time_contrast, "state_set") - time_contrast$std.error
        )),
        1e-10,
        label = paste(what, "R13 time in state")
      )
    }
  }
})

test_that("restarting the first-order recursion part-way gives the same SOPs (R27)", {
  # S_4(l) = sum over non-absorbing k of S_2(k) * P(state l at visit 4 | state
  # k at visit 2), plus S_2(l) if l is absorbing. A derived gap starts again
  # at the first requested visit, so the gap set-up is left out.
  cells <- Filter(
    function(cell) cell$order == 1L && cell$visit != "gap",
    grid_consistency_cells()
  )
  for (cell in cells) {
    label <- grid_label(cell)
    times <- grid_times(cell)
    profiles <- grid_profiles(cell)
    full <- grid_call(sops, cell, newdata = profiles)
    later <- grid_public_args(cell, times = times[3:4])
    for (i in seq_len(nrow(profiles))) {
      mine <- full[full$rowid == i, , drop = FALSE]
      at <- function(visit) {
        rows <- mine[as.character(mine$time) == as.character(times[visit]), ]
        rows$estimate[match(seq_len(grid_n_states), as.integer(as.character(rows$state)))]
      }
      second <- at(2)
      restarted <- c(rep(0, grid_n_states - 1L), second[grid_n_states])
      for (k in seq_len(grid_n_states - 1L)) {
        start <- profiles[i, , drop = FALSE]
        start$yprev <- if (is.factor(start$yprev)) {
          factor(k, levels = levels(start$yprev))
        } else {
          k
        }
        step <- grid_call(sops, cell, newdata = start, args = later)
        step <- step[as.character(step$time) == as.character(times[4]), ]
        restarted <- restarted + second[k] *
          step$estimate[match(seq_len(grid_n_states), as.integer(as.character(step$state)))]
      }
      expect_lt(max(abs(restarted - at(4))), 1e-10, label = paste(label, "profile", i))
    }
  }
})

test_that("swapping the compared values flips the contrast (N1)", {
  skip_on_cran()
  swapped <- list(tx = c(1, 0))
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    settings <- list(
      list(estimand = "sop"),
      list(estimand = "time_in_state"),
      list(estimand = "time_benefit")
    )
    for (setting in settings) {
      what <- paste(label, setting$estimand)
      forward <- do.call(
        grid_call,
        c(list(avg_comparisons, cell, variables = grid_treatments), setting)
      )
      backward <- do.call(
        grid_call,
        c(list(avg_comparisons, cell, variables = swapped), setting)
      )
      keys <- intersect(c("time", "state_set"), names(forward))
      flipped <- backward
      flipped$estimate <- -flipped$estimate
      expect_lt(grid_gap(forward, flipped, keys), 1e-12, label = what)

      if (grid_has_delta(cell) && setting$estimand != "time_benefit") {
        expect_lt(
          grid_gap(
            grid_delta(forward, "conditional"),
            grid_delta(backward, "conditional"),
            keys,
            "std.error"
          ),
          1e-12,
          label = paste(what, "standard error")
        )
      }
    }

    ratio <- function(variables) {
      grid_call(
        avg_comparisons,
        cell,
        variables = variables,
        estimand = "sop",
        comparison = "ratio"
      )
    }
    forward <- ratio(grid_treatments)
    backward <- ratio(swapped)
    backward$estimate <- 1 / backward$estimate
    expect_lt(
      grid_gap(forward, backward, c("time", "state_set")),
      1e-10,
      label = paste(label, "ratio")
    )
  }
})

test_that("a variable that is not in the model changes nothing (N24)", {
  skip_on_cran()
  markers <- list(marker = c(0, 1))
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    average <- grid_call(avg_sops, cell, variables = markers)
    expect_lt(
      grid_gap(
        average[average$marker == 0, , drop = FALSE],
        average[average$marker == 1, , drop = FALSE],
        c("time", "state")
      ),
      1e-14,
      label = label
    )
    for (estimand in c("sop", "time_in_state", "time_benefit")) {
      contrast <- grid_call(
        avg_comparisons,
        cell,
        variables = markers,
        estimand = estimand
      )
      expect_lt(max(abs(contrast$estimate)), 1e-14, label = paste(label, estimand))
      # Such a comparison is 0 by construction, so its uncertainty is not
      # reported, as in marginaleffects.
      if (grid_has_delta(cell) && estimand != "time_benefit") {
        for (vcov in grid_vcov_types(cell)) {
          inferred <- grid_delta(contrast, vcov)
          for (column in c("std.error", "statistic", "p.value", "conf.low", "conf.high")) {
            expect_true(
              all(is.na(inferred[[column]])),
              label = paste(label, estimand, vcov, column)
            )
          }
        }
      }
    }
  }
})

test_that("results at early visits do not depend on later visits (L28)", {
  for (cell in grid_consistency_cells()) {
    label <- grid_label(cell)
    times <- grid_times(cell)
    early <- grid_public_args(cell, times = times[1:3])
    expect_lt(
      grid_gap(
        grid_call(sops, cell, args = early),
        grid_call(sops, cell),
        c("rowid", "time", "state"),
        all_rows = FALSE
      ),
      1e-12,
      label = paste(label, "sops")
    )

    # Adding a visit adds that visit's SOP to the total.
    short <- as.data.frame(grid_call(
      avg_state_time,
      cell,
      variables = grid_treatments,
      args = early
    ))
    long <- grid_call(avg_state_time, cell, variables = grid_treatments)
    last <- grid_call(avg_sops, cell, variables = grid_treatments)
    last <- last[as.character(last$time) == as.character(times[4]), , drop = FALSE]
    keys <- c("tx", "state_set")
    short$estimate <- short$estimate +
      last$estimate[match(grid_key(short, keys), paste(last$tx, last$state, sep = "|"))]
    expect_lt(grid_gap(short, long, keys), 1e-12, label = paste(label, "totals"))
  }
})

test_that("orm and vglm fits of the same proportional-odds model give the same results (R8)", {
  skip_on_cran()
  # Both engines maximize the same likelihood, so results differ only by
  # convergence error. At the default convergence criteria the coefficients
  # differ by up to 4e-5 here (factor time), which moves SOPs by up to 2e-6
  # and standard errors by about 1e-7. The tolerance is 1e-5 on the
  # probability scale, and the same relative to the follow-up for time totals
  # (times 4 visits, or 14 in real time), so that it scales with the totals.
  orm_cells <- Filter(
    function(cell) cell$engine == "orm" && cell$odds == "full",
    if (grid_on_cran()) list(grid_cell("orm")) else grid_consistency_cells()
  )
  for (orm_cell in orm_cells) {
    vglm_cell <- orm_cell
    vglm_cell$engine <- "vglm"
    label <- grid_label(orm_cell)
    time_map <- grid_time_map(orm_cell)
    both <- function(fun, ...) {
      list(grid_call(fun, orm_cell, ...), grid_call(fun, vglm_cell, ...))
    }
    checks <- list(
      sops = list(
        both(sops, newdata = grid_profiles(orm_cell)),
        c("rowid", "time", "state"),
        1e-5
      ),
      avg_sops = list(both(avg_sops, variables = grid_treatments), c("tx", "time", "state"), 1e-5),
      state_time = list(
        both(avg_state_time, variables = grid_treatments, time_map = time_map),
        c("tx", "state_set"),
        14e-5
      ),
      sop_contrast = list(
        both(avg_comparisons, variables = grid_treatments, estimand = "sop"),
        c("time", "state_set"),
        1e-5
      ),
      time_contrast = list(
        both(avg_comparisons, variables = grid_treatments, estimand = "time_in_state"),
        "state_set",
        length(grid_times(orm_cell)) * 1e-5
      ),
      time_benefit = list(
        both(
          avg_comparisons,
          variables = grid_treatments,
          estimand = "time_benefit",
          time_map = time_map
        ),
        character(),
        14e-5
      )
    )
    for (name in names(checks)) {
      pair <- checks[[name]][[1]]
      keys <- checks[[name]][[2]]
      tolerance <- checks[[name]][[3]]
      expect_lt(grid_gap(pair[[1]], pair[[2]], keys), tolerance, label = paste(label, name))
      if (!grid_has_delta(orm_cell) || name == "time_benefit") next
      types <- if (name == "sops") "conditional" else grid_vcov_types(orm_cell)
      for (vcov in types) {
        expect_lt(
          grid_gap(
            grid_delta(pair[[1]], vcov),
            grid_delta(pair[[2]], vcov),
            keys,
            "std.error"
          ),
          tolerance,
          label = paste(label, name, vcov, "standard error")
        )
      }
    }
  }
})

# Free the cached results; other files do not use them.
grid_cache$results <- NULL
