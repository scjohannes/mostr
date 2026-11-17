# The vglm engine builds its prediction design with stats::model.frame(),
# which leaves a numeric column numeric even when the fit used it as a factor;
# model.matrix() then stops with "contrasts apply only to factors". The orm
# engine turns such columns into factors with the fitted levels, so the same
# call must give the same answer with both engines (issue #44).

# Three states (1-2 -> 1, 3-7 -> 2, 8 -> 3), with state 3 absorbing, so
# `yprev` is a factor with levels "1" and "2".
vglm_factor_data <- function(factor_time = FALSE) {
  trajectories <- sim_markov("actt2", n_patients = 300, times = 1:6, seed = 1)
  collapse <- function(s) ifelse(s <= 2, 1L, ifelse(s <= 7, 2L, 3L))
  trajectories$y <- collapse(trajectories$y)
  trajectories$yprev <- collapse(trajectories$yprev)
  data <- prepare_markov_data(trajectories, absorb = 3)
  data$txf <- factor(data$tx)
  if (factor_time) {
    data$time <- factor(data$time)
  }
  data
}

# The fits are shared by the tests below, so each is made only once.
vglm_factor_cache <- new.env(parent = emptyenv())

vglm_factor_fits <- function(factor_time = FALSE) {
  key <- if (factor_time) "factor_time" else "numeric_time"
  if (!is.null(vglm_factor_cache[[key]])) {
    return(vglm_factor_cache[[key]])
  }
  data <- vglm_factor_data(factor_time)
  fit <- function(engine) {
    args <- list(
      ordered(y) ~ yprev + time + txf,
      data = data,
      id_var = "id",
      absorb = "3",
      engine = engine
    )
    if (factor_time) {
      args$first_followup_time <- "1"
    }
    do.call(mostm, args)
  }
  vglm_factor_cache[[key]] <- list(
    data = data,
    orm = fit("orm"),
    vglm = fit("vglm")
  )
  vglm_factor_cache[[key]]
}

test_that("vglm predictions accept numeric values of factor predictors", {
  skip_if_not_installed("VGAM")
  fits <- vglm_factor_fits()
  newdata <- data.frame(
    id = 1:4,
    yprev = factor(c(1, 2, 2, 1)),
    time = 1,
    txf = factor(c(0, 0, 1, 1))
  )

  as_factor <- avg_sops(
    fits$vglm,
    newdata = newdata,
    variables = list(txf = factor(c(0, 1))),
    times = 1:3
  )
  as_number <- avg_sops(
    fits$vglm,
    newdata = newdata,
    variables = list(txf = c(0, 1)),
    times = 1:3
  )
  from_orm <- avg_sops(
    fits$orm,
    newdata = newdata,
    variables = list(txf = c(0, 1)),
    times = 1:3
  )
  expect_equal(as_number$estimate, as_factor$estimate)
  expect_equal(as_number$estimate, from_orm$estimate, tolerance = 1e-6)

  numeric_rows <- transform(newdata, yprev = c(1, 2, 2, 1), txf = c(0, 0, 1, 1))
  expect_equal(
    sops(fits$vglm, newdata = numeric_rows, times = 1:3)$estimate,
    sops(fits$vglm, newdata = newdata, times = 1:3)$estimate
  )
})

test_that("plot_transitions() takes numeric previous states for a vglm fit", {
  skip_if_not_installed("VGAM")
  fits <- vglm_factor_fits()
  tiles <- function(model, values) {
    plot_transitions(
      model,
      times = 2:3,
      variables = list(yprev = values),
      comparison = "difference"
    )$data$estimate
  }

  as_number <- tiles(fits$vglm, c(1, 2))
  expect_equal(as_number, tiles(fits$vglm, factor(c(1, 2))))
  expect_equal(as_number, tiles(fits$orm, c(1, 2)), tolerance = 1e-6)
})

test_that("vglm predictions accept numeric or character values of factor time", {
  skip_if_not_installed("VGAM")
  fits <- vglm_factor_fits(factor_time = TRUE)
  rows <- function(time) {
    data.frame(id = 1:2, yprev = factor(c(1, 2)), txf = factor(0), time = time)
  }
  time_levels <- levels(fits$data$time)
  reference <- sops(
    fits$vglm,
    newdata = rows(factor("1", levels = time_levels)),
    times = 1:3
  )$estimate

  expect_equal(
    sops(fits$vglm, newdata = rows(1), times = 1:3)$estimate,
    reference
  )
  expect_equal(
    sops(fits$vglm, newdata = rows("1"), times = 1:3)$estimate,
    reference
  )
  expect_equal(
    sops(fits$orm, newdata = rows(1), times = 1:3)$estimate,
    reference,
    tolerance = 1e-6
  )
})

test_that("vglm predictions refuse values outside the fitted levels", {
  skip_if_not_installed("VGAM")
  fits <- vglm_factor_fits()
  newdata <- data.frame(id = 1:2, yprev = factor(c(1, 2)), time = 1, txf = 0)

  expect_error(
    avg_sops(
      fits$vglm,
      newdata = newdata,
      variables = list(txf = c(0, 5)),
      times = 1:2
    ),
    "`txf`",
    class = "mostr_error"
  )
  expect_error(
    sops(fits$vglm, newdata = transform(newdata, txf = "a"), times = 1:2),
    "`txf`",
    class = "mostr_error"
  )
})
