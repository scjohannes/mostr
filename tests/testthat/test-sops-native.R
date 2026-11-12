test_that("PO structure detection requires threshold-invariant slopes", {
  X <- cbind("(Intercept)" = 1, x = c(-1, 0, 1))
  po <- mostr:::markov_po_structure(
    rbind(c(1, 0.5), c(-1, 0.5)),
    colnames(X),
    X
  )
  non_po <- mostr:::markov_po_structure(
    rbind(c(1, 0.5), c(-1, 0.7)),
    colnames(X),
    X
  )

  expect_equal(po$cutpoints, c(1, -1))
  expect_equal(unname(po$beta), c(0, 0.5))
  expect_null(non_po)
})

test_that("vgam fits are no longer supported", {
  model <- with_fit_settings(structure(list(), class = "vgam"), y_levels = 1:2)
  expect_snapshot(
    error = TRUE,
    soprob_markov(
      model,
      newdata = data.frame(time = 1, yprev = 1),
      times = 1
    )
  )
})
