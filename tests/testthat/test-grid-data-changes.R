# Changes to the fitting data that should not change any result, refitted on
# six base set-ups: engine x time, with and without case weights. The checks
# are in helper-grid-changes.R; test-slow-grid-data-changes.R runs them on the
# other six. The IDs in the test names refer to the catalogue of claims behind
# the test overhaul.
#
# These refit many models, so they do not run on CRAN.

test_that("the order of the data rows does not matter (N13)", {
  skip_on_cran()
  for (cell in change_base_cells("default")) check_row_order(cell)
})

test_that("how patient IDs are written does not matter (N14)", {
  skip_on_cran()
  for (cell in change_base_cells("default")) check_id_writing(cell)
})

test_that("recoding the treatment does not matter (N15)", {
  skip_on_cran()
  for (cell in change_base_cells("default")) check_treatment_coding(cell)
})

test_that("the units of a covariate do not matter (N3)", {
  skip_on_cran()
  for (cell in change_base_cells("default")) check_covariate_units(cell)
})

test_that("a column the model does not use does not matter, even with missing values (N20)", {
  skip_on_cran()
  for (cell in change_base_cells("default")) check_unused_column(cell)
})

test_that("integer case weights equal duplicated patients where they should (N9)", {
  skip_on_cran()
  for (cell in Filter(function(cell) cell$weighted, change_base_cells("default"))) {
    check_weights_against_duplicates(cell)
  }
})

test_that("a constant case weight changes neither estimates nor robust standard errors (N9)", {
  skip_on_cran()
  for (cell in Filter(function(cell) cell$weighted, change_base_cells("default"))) {
    check_constant_weight(cell)
  }
})

test_that("bootstrap draws do not depend on how patient IDs are written (N14)", {
  skip_on_cran()
  # Patients are drawn by their position in the order in which they first
  # appear, which these relabellings keep, so the same seed gives the same
  # draws. Shuffling the rows changes that order, so it is not compared here.
  cells <- list(grid_cell("orm"), grid_cell("vglm", time = "factor"))
  for (cell in cells) {
    args <- grid_public_args(cell)
    with_draws <- function(fit, method) {
      average <- do.call(avg_sops, c(list(fit, variables = change_treatments), args))
      as.data.frame(inferences(average, method = method, n_draws = 3, seed = 42))
    }
    for (method in c("bootstrap", "fwb")) {
      base <- with_draws(grid_fit(cell), method)
      for (name in names(id_changes)) {
        changed <- with_draws(grid_refit(cell, name, id_changes[[name]]$change), method)
        for (column in c("std.error", "conf.low", "conf.high")) {
          expect_lt(
            grid_gap(changed, base, c("tx", "time", "state"), column),
            same_rows_tolerance,
            label = paste(grid_label(cell), name, method, column)
          )
        }
      }
    }
  }
})
