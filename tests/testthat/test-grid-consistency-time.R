# Consistency checks on every model set-up of the shared grid, for time in
# state, comparisons and time benefit; see test-grid-consistency.R for the
# rest. The IDs in the test names (L6, R2, ...) refer to the catalogue of
# claims behind the test overhaul.
#
# On CRAN only one set-up per engine and Markov order runs
# (`grid_consistency_cells()`); everything here is deterministic. Settings
# that only vary the request, such as the `by` groups or the real-time grid,
# take turns across set-ups (`grid_turn()`) rather than being crossed with
# every set-up, which keeps the run time near a minute.

test_that("time totals lie between 0 and the follow-up and add up to it (L6, L7, L8)", {
  skip_on_cran()
  lumps <- list(low = c("1", "2"), high = c("3", "4"))
  cells <- grid_consistency_cells()
  for (index in seq_along(cells)) {
    cell <- cells[[index]]
    time_map <- grid_time_map(cell)
    # Follow-up: one unit per visit, 14 from baseline time 0, 12 from the
    # first mapped visit (time 2).
    scale <- grid_turn(
      list(
        list(name = "visits", extra = list(), span = length(grid_times(cell))),
        list(name = "real time", extra = list(time_map = time_map), span = 14),
        list(
          name = "real time from the first visit",
          extra = list(time_map = time_map, baseline_time = NULL),
          span = 12
        )
      ),
      index
    )
    label <- paste(grid_label(cell), scale$name)
    for (sets in list(NULL, lumps)) {
      what <- paste(label, if (is.null(sets)) "states" else "lumps")
      average <- do.call(
        grid_call,
        c(list(avg_state_time, cell, variables = grid_treatments, state_sets = sets), scale$extra)
      )
      expect_gte(min(average$estimate), -1e-10, label = what)
      expect_lte(max(average$estimate), scale$span + 1e-10, label = what)
      expect_lt(
        max(abs(tapply(average$estimate, average$tx, sum) - scale$span)),
        1e-9,
        label = what
      )
    }
    real <- grid_profile_real_time(cell, index)
    profile_scales <- list(
      list(extra = list(), span = length(grid_times(cell))),
      list(extra = real, span = if (is.null(real$baseline_time)) 12 else 14)
    )
    for (profile_scale in profile_scales) {
      each <- do.call(
        grid_call,
        c(list(state_time, cell, variables = grid_treatments), profile_scale$extra)
      )
      what <- paste(grid_label(cell), "each profile, span", profile_scale$span)
      expect_gte(min(each$estimate), -1e-10, label = what)
      expect_lte(max(each$estimate), profile_scale$span + 1e-10, label = what)
      expect_lt(
        max(abs(
          tapply(each$estimate, grid_key(each, c("rowid", "tx")), sum) -
            profile_scale$span
        )),
        1e-9,
        label = what
      )
    }

    # A state set and its complement: their differences cancel.
    lumped <- grid_call(
      avg_comparisons,
      cell,
      variables = grid_treatments,
      estimand = "sop",
      state_sets = lumps
    )
    low <- lumped[lumped$state_set == "low", , drop = FALSE]
    high <- lumped[lumped$state_set == "high", , drop = FALSE]
    expect_lte(max(abs(lumped$estimate)), 1, label = grid_label(cell))
    expect_lt(
      grid_gap(low, transform(high, estimate = -estimate), "time"),
      1e-12,
      label = paste(grid_label(cell), "lump and complement")
    )
  }
})

test_that("time benefit lies within plus or minus the follow-up (L19)", {
  skip_on_cran()
  cells <- grid_consistency_cells()
  for (index in seq_along(cells)) {
    cell <- cells[[index]]
    label <- grid_label(cell)
    for (real in c(FALSE, grid_real_time_benefit_turn(index))) {
      span <- if (real) 14 else length(grid_times(cell))
      extra <- if (real) list(time_map = grid_time_map(cell)) else list()
      for (fun in list(comparisons, avg_comparisons)) {
        benefit <- do.call(
          grid_call,
          c(list(fun, cell, variables = grid_treatments, estimand = "time_benefit"), extra)
        )
        expect_lte(
          max(abs(benefit$estimate)),
          span,
          label = paste(label, if (real) "real time" else "visits")
        )
      }
    }
  }
})

test_that("avg_comparisons() is the mean of comparisons(), overall and within groups (R2)", {
  skip_on_cran()
  # by = "tx" names the compared variable, so it groups by the treatment each
  # patient was observed with. Every set-up averages overall and takes its
  # turn with one kind of group and one real-time estimand.
  cells <- grid_consistency_cells()
  for (index in seq_along(cells)) {
    cell <- cells[[index]]
    label <- grid_label(cell)
    settings <- list(
      list(estimand = "sop"),
      list(estimand = "time_in_state"),
      list(estimand = "time_benefit"),
      list(
        estimand = if (grid_real_time_benefit_turn(index)) "time_benefit" else "time_in_state",
        time_map = grid_time_map(cell)
      )
    )
    for (setting in settings) {
      extra <- c(list(variables = grid_treatments), setting)
      each <- as.data.frame(do.call(grid_call, c(list(comparisons, cell), extra)))
      # Pairs of set-ups share a group, so that each group meets both
      # real-time estimands.
      for (by in list(TRUE, grid_turn(list("tx", "marker"), (index + 1L) %/% 2L))) {
        average <- do.call(
          grid_call,
          c(list(avg_comparisons, cell), extra, if (is.character(by)) list(by = by))
        )
        keys <- intersect(c("time", "state_set", "contrast"), names(average))
        groups <- c(keys, if (is.character(by)) by)
        means <- stats::aggregate(each["estimate"], each[groups], mean)
        expect_lt(
          grid_gap(average, means, groups),
          1e-10,
          label = paste(
            label,
            setting$estimand,
            if (!is.null(setting$time_map)) "real time",
            "by",
            by
          )
        )
      }
    }
  }
})

test_that("time in state is the sum or area of the SOPs (R5)", {
  skip_on_cran()
  cells <- grid_consistency_cells()
  for (index in seq_along(cells)) {
    cell <- cells[[index]]
    label <- grid_label(cell)
    times <- grid_times(cell)
    time_map <- grid_time_map(cell)

    # Visits: each requested visit counts one unit.
    each <- as.data.frame(grid_call(sops, cell, variables = grid_treatments))
    sums <- stats::aggregate(each["estimate"], each[c("rowid", "tx", "state")], sum)
    sums$state_set <- as.character(sums$state)
    per_profile <- grid_call(state_time, cell, variables = grid_treatments)
    expect_lt(
      grid_gap(per_profile, sums, c("rowid", "tx", "state_set")),
      1e-10,
      label = paste(label, "each profile")
    )
    average <- as.data.frame(grid_call(avg_sops, cell, variables = grid_treatments))
    sums <- stats::aggregate(average["estimate"], average[c("tx", "state")], sum)
    sums$state_set <- as.character(sums$state)
    totals <- grid_call(avg_state_time, cell, variables = grid_treatments)
    expect_lt(
      grid_gap(totals, sums, c("tx", "state_set")),
      1e-10,
      label = paste(label, "average")
    )

    # Real time: straight lines between the baseline distribution and the
    # mapped visits, and the trapezoid area under them.
    start <- grid_start_distribution(cell)
    for (baseline in list(0, NULL)) {
      expected <- do.call(rbind, lapply(c(0, 1), function(tx) {
        do.call(rbind, lapply(seq_len(grid_n_states), function(state) {
          rows <- average[average$tx == tx & as.character(average$state) == state, ]
          rows <- rows[match(as.character(times), as.character(rows$time)), ]
          x <- c(baseline, unname(time_map))
          y <- c(if (!is.null(baseline)) start[[state]], rows$estimate)
          data.frame(tx = tx, state_set = as.character(state), estimate = grid_trapezoid(x, y))
        }))
      }))
      real <- grid_call(
        avg_state_time,
        cell,
        variables = grid_treatments,
        time_map = time_map,
        baseline_time = baseline
      )
      what <- paste(label, "real time from", if (is.null(baseline)) "first visit" else "baseline")
      expect_lt(grid_gap(real, expected, c("tx", "state_set")), 1e-10, label = what)
    }

    # The average of the per-profile areas is the area of the average.
    real_time <- grid_profile_real_time(cell, index)
    baseline <- real_time$baseline_time
    real <- do.call(
      grid_call,
      c(list(avg_state_time, cell, variables = grid_treatments), real_time)
    )
    per_profile <- as.data.frame(do.call(
      grid_call,
      c(list(state_time, cell, variables = grid_treatments), real_time)
    ))
    means <- stats::aggregate(
      per_profile["estimate"],
      per_profile[c("tx", "state_set")],
      mean
    )
    what <- paste(label, "mean of real-time totals, baseline", format(baseline))
    expect_equal(nrow(per_profile), nrow(each) / length(times), label = what)
    expect_lt(grid_gap(real, means, c("tx", "state_set")), 1e-10, label = what)
  }
})

test_that("time_in_state() on SOPs matches avg_state_time() and state_time() (R11)", {
  skip_on_cran()
  cells <- grid_consistency_cells()
  for (index in seq_along(cells)) {
    cell <- cells[[index]]
    label <- grid_label(cell)
    time_map <- grid_time_map(cell)
    average <- grid_call(avg_sops, cell, variables = grid_treatments)
    grid <- grid_turn(
      list(
        list(baseline_time = 0, target_times = NULL),
        list(baseline_time = 0, target_times = c(1, 3, 6, 10, 14)),
        list(baseline_time = NULL, target_times = NULL)
      ),
      index
    )
    what <- paste(
      label,
      "baseline",
      format(grid$baseline_time),
      "targets",
      paste(grid$target_times, collapse = " ")
    )
    for (set in list("4", c("1", "2"))) {
      expected <- grid_call(
        avg_state_time,
        cell,
        variables = grid_treatments,
        state_sets = set,
        time_map = time_map,
        baseline_time = grid$baseline_time,
        target_times = grid$target_times
      )
      totals <- time_in_state(
        average,
        target_states = set,
        time_map = time_map,
        baseline_time = grid$baseline_time,
        target_times = grid$target_times
      )
      expect_lt(
        grid_gap(totals, expected, "tx", "total_time", "estimate"),
        1e-10,
        label = paste(what, "states", paste(set, collapse = " "))
      )
    }

    # One total per profile from sops() on the stored profiles; every third
    # set-up, as it is slow.
    if (index %% 3L == 1L) {
      per_profile <- time_in_state(
        grid_call(sops, cell),
        target_states = "1",
        time_map = time_map
      )
      expected <- grid_call(state_time, cell, state_sets = "1", time_map = time_map)
      expect_lt(
        grid_gap(per_profile, expected, "rowid", "total_time", "estimate"),
        1e-10,
        label = paste(label, "each profile")
      )
    }

    # Null-hypothesis columns from inferences() are not grouping columns.
    if (!grid_has_delta(cell)) next
    tested <- inferences(average, method = "delta", vcov = "conditional", null = 0.5)
    expect_lt(
      grid_gap(
        suppressWarnings(time_in_state(tested, target_states = "1", time_map = time_map)),
        time_in_state(average, target_states = "1", time_map = time_map),
        "tx",
        "total_time"
      ),
      1e-12,
      label = paste(label, "after a null test")
    )
    expect_identical(
      nrow(suppressWarnings(interpolate_sops(tested, time_map = time_map))),
      nrow(interpolate_sops(average, time_map = time_map)),
      label = paste(label, "interpolated after a null test")
    )
  }
})

test_that("time benefit is the formula applied to each profile's paired SOPs (R18)", {
  skip_on_cran()
  # At each visit, with reference distribution P0 and comparison P1 (states
  # ordered from best to worst): sum over s of P0(s) * (P1(below s) - P1(above s)).
  cells <- grid_consistency_cells()
  for (index in seq_along(cells)) {
    cell <- cells[[index]]
    label <- grid_label(cell)
    times <- grid_times(cell)
    each <- as.data.frame(grid_call(sops, cell, variables = grid_treatments))
    patients <- unique(as.character(each$id))
    probabilities <- lapply(c(0, 1), function(tx) {
      rows <- each[each$tx == tx, , drop = FALSE]
      out <- array(0, c(length(patients), length(times), grid_n_states))
      out[cbind(
        match(as.character(rows$id), patients),
        match(as.character(rows$time), as.character(times)),
        as.integer(as.character(rows$state))
      )] <- rows$estimate
      out
    })
    reference <- probabilities[[1]]
    comparison <- probabilities[[2]]
    benefit <- matrix(0, length(patients), length(times))
    for (state in seq_len(grid_n_states)) {
      below <- apply(comparison[, , seq_len(state - 1L), drop = FALSE], c(1, 2), sum)
      above <- apply(
        comparison[, , setdiff(seq_len(grid_n_states), seq_len(state)), drop = FALSE],
        c(1, 2),
        sum
      )
      benefit <- benefit + reference[, , state] * (below - above)
    }

    visits <- grid_call(
      comparisons,
      cell,
      variables = grid_treatments,
      estimand = "time_benefit"
    )
    patient <- match(as.character(visits$id), patients)
    expect_lt(
      max(abs(visits$estimate - rowSums(benefit)[patient])),
      1e-10,
      label = paste(label, "visits")
    )

    # In real time, both scenarios start in the profile's own state, where the
    # benefit is 0.
    if (!grid_real_time_benefit_turn(index)) next
    time_map <- grid_time_map(cell)
    real <- grid_call(
      comparisons,
      cell,
      variables = grid_treatments,
      estimand = "time_benefit",
      time_map = time_map
    )
    x <- c(0, unname(time_map))
    areas <- apply(benefit, 1, function(row) grid_trapezoid(x, c(0, row)))
    expect_lt(
      max(abs(real$estimate - areas[match(as.character(real$id), patients)])),
      1e-10,
      label = paste(label, "real time")
    )
  }
})

# Free the cached results; other files do not use them.
grid_cache$results <- NULL
