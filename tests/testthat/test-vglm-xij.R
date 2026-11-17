# With `xij`, VGAM lets a predictor take a different value in each linear
# predictor; `form2` builds those columns. The SOP predictions build only the
# ordinary design matrix, so they used the same column in every linear
# predictor and returned wrong probabilities without a word (issue #45). With
# xij = list(z ~ z1 + z2), the probability of state 3 was 0.00664 where
# VGAM::predict() gave 0.00786. mostm() therefore refuses both arguments.

xij_data <- function() {
  trajectories <- sim_markov("actt2", n_patients = 300, times = 1:6, seed = 1)
  collapse <- function(s) ifelse(s <= 2, 1L, ifelse(s <= 7, 2L, 3L))
  trajectories$y <- collapse(trajectories$y)
  trajectories$yprev <- collapse(trajectories$yprev)
  data <- prepare_markov_data(trajectories, absorb = 3)
  withr::with_seed(2, {
    data$z1 <- stats::rnorm(nrow(data))
    data$z2 <- data$z1 + stats::rnorm(nrow(data))
  })
  data$z <- data$z1
  data
}

test_that("mostm(engine = 'vglm') refuses xij, given directly or in control", {
  skip_if_not_installed("VGAM")
  data <- xij_data()

  expect_error(
    mostm(
      ordered(y) ~ yprev + time + tx + z,
      data = data,
      id_var = "id",
      absorb = "3",
      engine = "vglm",
      xij = list(z ~ z1 + z2),
      form2 = ~ yprev + time + tx + z + z1 + z2
    ),
    "`xij`",
    class = "mostr_error"
  )
  expect_error(
    mostm(
      ordered(y) ~ yprev + time + tx + z,
      data = data,
      id_var = "id",
      absorb = "3",
      engine = "vglm",
      control = VGAM::vglm.control(xij = list(z ~ z1 + z2))
    ),
    "`xij`",
    class = "mostr_error"
  )
})

test_that("mostm(engine = 'vglm') refuses form2 without xij", {
  skip_if_not_installed("VGAM")
  data <- xij_data()

  expect_error(
    mostm(
      ordered(y) ~ yprev + time + tx + z,
      data = data,
      id_var = "id",
      absorb = "3",
      engine = "vglm",
      form2 = ~ yprev + time + tx + z + z1 + z2
    ),
    "`form2`",
    class = "mostr_error"
  )
})
