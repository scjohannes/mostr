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
