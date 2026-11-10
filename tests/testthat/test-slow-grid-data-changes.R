# The checks of test-grid-data-changes.R on the other six base set-ups, which
# swap which engine-and-time cells have case weights. Together the two files
# cover engine x time x case weights in full. Set MOSTR_RUN_SLOW=true to run.

test_that("data changes that should not matter, on the remaining base set-ups (N3, N9, N13, N14, N15, N20)", {
  skip_on_cran()
  skip_unless_slow_tier()
  for (cell in change_base_cells("slow")) {
    check_row_order(cell)
    check_id_writing(cell)
    check_treatment_coding(cell)
    check_covariate_units(cell)
    check_unused_column(cell)
    if (cell$weighted) {
      check_weights_against_duplicates(cell)
      check_constant_weight(cell)
    }
  }
})
