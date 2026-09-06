test_that("prepare_markov_data can keep previous state numeric", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(2, 3, 4, 5),
    yprev = c(1, 2, 4, 6)
  )

  prepared <- prepare_markov_data(
    data,
    absorbing_state = 6,
    factor_previous = FALSE
  )

  expect_type(prepared$yprev, "double")
  expect_s3_class(prepared$y, "ordered")
  expect_equal(prepared$yprev, c(1, 2, 4))
})
