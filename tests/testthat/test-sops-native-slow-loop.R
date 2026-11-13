# Are the fast paths right? These tests compare public results with the readable
# loops in helper-slow-loop.R, which use only the fitting package's predict(), so
# a bug anywhere between the fit and the output has nothing to hide behind. That
# includes how prediction rows are built, which comparisons on prebuilt designs
# cannot see.
#
# The file name keeps these tests in the native-safety CI jobs.

test_that("sops() matches the readable loop on every model set-up", {
  # On CRAN, check the four set-ups whose fits the consistency checks already
  # cache, so that one comparison with the readable loop runs there cheaply.
  cells <- if (grid_on_cran()) grid_consistency_cells() else grid_core_cells()
  for (cell in cells) {
    fit <- grid_fit(cell)
    profiles <- grid_profiles(cell)
    result <- do.call(
      sops,
      c(list(fit, newdata = profiles), grid_public_args(cell))
    )
    reference <- grid_slow_sops(cell, fit, profiles, grid_n_states)

    expect_identical(nrow(result), nrow(profiles) * 4L * grid_n_states)
    expect_lt(
      grid_largest_gap(result, reference, result$rowid, grid_times(cell)),
      1e-9,
      label = grid_label(cell)
    )
  }
})

test_that("sops() matches the loop with several absorbing states, inner ones included", {
  skip_on_cran()
  for (cell in lapply(grid_core_cells(), grid_inner_absorbing)) {
    fit <- grid_fit(cell)
    profiles <- grid_profiles(cell)
    result <- do.call(
      sops,
      c(list(fit, newdata = profiles), grid_public_args(cell))
    )
    reference <- grid_slow_sops(cell, fit, profiles, cell$absorb)

    expect_lt(
      grid_largest_gap(result, reference, result$rowid, grid_times(cell)),
      1e-9,
      label = grid_label(cell)
    )
  }
})

test_that("sops() matches the loop for gaps, log time, numeric states and weights", {
  skip_on_cran()
  cells <- list(
    grid_cell("orm", visit = "gap"),
    grid_cell("vglm", visit = "gap"),
    grid_cell("orm", time = "log"),
    grid_cell("vglm", time = "log"),
    grid_cell("orm", previous = "numeric"),
    grid_cell("vglm", previous = "numeric"),
    grid_cell("orm", weighted = TRUE),
    grid_cell("vglm", weighted = TRUE)
  )
  for (cell in cells) {
    fit <- grid_fit(cell)
    profiles <- grid_profiles(cell)
    result <- suppressWarnings(do.call(
      sops,
      c(list(fit, newdata = profiles), grid_public_args(cell))
    ))
    reference <- grid_slow_sops(cell, fit, profiles, grid_n_states)

    expect_lt(
      grid_largest_gap(result, reference, result$rowid, grid_times(cell)),
      1e-9,
      label = grid_label(cell)
    )
  }
})

test_that("avg_sops() matches averaging the loop, within groups and under each treatment", {
  skip_on_cran()
  for (cell in grid_core_cells(times = c("linear", "factor"))) {
    fit <- grid_fit(cell)
    times <- grid_times(cell)
    profiles <- grid_profiles(cell)
    # Unequal groups: a wrong weighting would not cancel.
    profiles$marker <- c(1L, 1L, 0L)

    for (tx in c(0, 1)) {
      treated <- profiles
      treated$tx <- tx
      reference <- grid_slow_sops(cell, fit, treated, grid_n_states)
      average <- list(Reduce(`+`, reference) / length(reference))
      groups <- list(
        Reduce(`+`, reference[profiles$marker == 0L]) /
          sum(profiles$marker == 0L),
        Reduce(`+`, reference[profiles$marker == 1L]) /
          sum(profiles$marker == 1L)
      )

      overall <- do.call(
        avg_sops,
        c(
          list(fit, newdata = profiles, variables = list(tx = tx)),
          grid_public_args(cell)
        )
      )
      expect_lt(
        grid_largest_gap(overall, average, rep(1L, nrow(overall)), times),
        1e-9,
        label = grid_label(cell)
      )

      grouped <- do.call(
        avg_sops,
        c(
          list(fit, newdata = treated, by = "marker"),
          grid_public_args(cell)
        )
      )
      expect_lt(
        grid_largest_gap(grouped, groups, grouped$marker + 1L, times),
        1e-9,
        label = paste(grid_label(cell), "by group")
      )
    }
  }
})

test_that("sops() on the stored patients matches the loop on their baseline rows", {
  skip_on_cran()
  for (cell in grid_core_cells(times = c("linear", "factor"))) {
    fit <- grid_fit(cell)
    times <- grid_times(cell)
    result <- do.call(sops, c(list(fit), grid_public_args(cell)))

    data <- grid_cell_data(cell)
    baseline <- data[as.character(data$time) == "1", , drop = FALSE]
    checked <- baseline[1:3, , drop = FALSE]
    reference <- grid_slow_sops(cell, fit, checked, grid_n_states)
    keep <- result[result$id %in% checked$id, , drop = FALSE]

    expect_identical(nrow(keep), 3L * 4L * grid_n_states)
    expect_lt(
      grid_largest_gap(keep, reference, match(keep$id, checked$id), times),
      1e-9,
      label = grid_label(cell)
    )
  }
})

test_that("second-order recursion on a first-order fit reproduces first-order SOPs", {
  skip_on_cran()
  # A first-order model ignores the older state, so the same model fitted
  # with `p2_var`, which makes the predictions use second-order recursion,
  # must give the same SOPs.
  cells <- Filter(function(cell) cell$order == 1L, grid_core_cells())
  for (cell in cells) {
    fit <- grid_fit(cell)
    second_fit <- grid_fit_to(cell, grid_cell_data(cell), p2_var = "ypprev")
    profiles <- grid_profiles(cell)
    args <- grid_public_args(cell)
    first <- do.call(sops, c(list(fit, newdata = profiles), args))
    second <- do.call(sops, c(list(second_fit, newdata = profiles), args))
    expect_identical(markov_model_settings(second_fit)$order, 2L)

    expect_lt(
      max(abs(first$estimate - second$estimate)),
      1e-9,
      label = grid_label(cell)
    )
  }
})

test_that("analytical derivatives match numerical differences of public results", {
  skip_on_cran()
  # Analytical inference supports first-order, full proportional-odds fits.
  cells <- Filter(
    function(cell) cell$order == 1L && cell$odds == "full",
    grid_core_cells()
  )
  # The age column is large, which inflates the error of a central difference,
  # so this tolerance is looser than the loop's. Planted derivative bugs move
  # these numbers by 0.1 or more.
  tolerance <- 1e-5

  for (cell in cells) {
    fit <- grid_fit(cell)
    profiles <- grid_profiles(cell)
    args <- grid_public_args(cell)

    per_profile <- function(model) {
      do.call(sops, c(list(model, newdata = profiles), args))$estimate
    }
    point <- do.call(sops, c(list(fit, newdata = profiles), args))
    inferred <- inferences(point, method = "delta", vcov = "conditional")
    analytical <- get_jacobian(inferred)
    numerical <- numeric_jacobian(fit, per_profile)
    expect_lt(
      max(abs(analytical - numerical[, colnames(analytical), drop = FALSE])),
      tolerance,
      label = paste(grid_label(cell), "sops")
    )

    averaged <- function(model) {
      do.call(
        avg_sops,
        c(list(model, newdata = profiles, variables = list(tx = c(0, 1))), args)
      )$estimate
    }
    point <- do.call(
      avg_sops,
      c(list(fit, newdata = profiles, variables = list(tx = c(0, 1))), args)
    )
    inferred <- inferences(point, method = "delta", vcov = "conditional")
    analytical <- get_jacobian(inferred)
    numerical <- numeric_jacobian(fit, averaged)
    expect_lt(
      max(abs(analytical - numerical[, colnames(analytical), drop = FALSE])),
      tolerance,
      label = paste(grid_label(cell), "avg_sops")
    )
  }
})
