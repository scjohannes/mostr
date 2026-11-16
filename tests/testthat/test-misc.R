test_that("prepare_markov_data can keep previous state numeric", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(2, 3, 4, 5),
    yprev = c(1, 2, 4, 6)
  )

  prepared <- prepare_markov_data(
    data,
    absorb = 6,
    factor_previous = FALSE
  )

  expect_type(prepared$yprev, "double")
  expect_s3_class(prepared$y, "ordered")
  expect_equal(prepared$yprev, c(1, 2, 4))
})

test_that("prepare_markov_data requires absorb and drops nothing with NULL", {
  data <- data.frame(y = c(2, 3, 6, 6), yprev = c(1, 2, 3, 6))

  # A default of 6 would drop the last row even where 6 is not death.
  expect_error(
    prepare_markov_data(data),
    "`absorb` is required",
    fixed = TRUE,
    class = "mostr_error"
  )
  expect_equal(nrow(prepare_markov_data(data, absorb = NULL)), 4L)
  expect_equal(nrow(prepare_markov_data(data, absorb = 6)), 3L)
})

test_that("prepare_markov_data uses the state columns named by y_var and p_var", {
  data <- data.frame(
    status = c(2, 3, 6, 6),
    previous = c(1, 2, 3, 6),
    yprev = c(6, 6, 6, 6)
  )

  prepared <- prepare_markov_data(
    data,
    absorb = 6,
    y_var = "status",
    p_var = "previous"
  )

  expect_equal(as.character(prepared$previous), c("1", "2", "3"))
  expect_s3_class(prepared$status, "ordered")
  # A column that happens to be called `yprev` is left alone.
  expect_identical(prepared$yprev, c(6, 6, 6))
  expect_error(
    prepare_markov_data(data, absorb = 6),
    "Columns not in `data`: `y`",
    fixed = TRUE,
    class = "mostr_error"
  )
})
