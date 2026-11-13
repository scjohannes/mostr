# MVN draws against conditional analytical (delta) variance, given the same
# coefficient covariance.
#
# MVN pushes coefficient draws through the nonlinear SOP recursion; delta
# linearizes the recursion. With the model's covariance V they differ slightly
# however many draws are taken, so they are not compared there. With
# vcov = eps * V and small eps the draws stay so close to the estimate that
# the recursion is linear for them, and the MVN standard deviation divided by
# the delta standard error tends to 1. Its Monte Carlo error with B draws is
# about 1 / sqrt(2 B) relative, 0.5% for B = 20000, so 2% allows for the
# largest of the cells. With eps = 1e-3 the nonlinearity is of order eps,
# far below that.
#
# Delta inference supports first-order fits only, so there is no second-order
# set-up. Both methods use the same 40 supplied starting profiles, which delta
# treats conditionally.
#
# About 4 minutes in all (most of it the 20,000 draws, about 1.5 minutes per
# set-up), so this runs only with MOSTR_RUN_SLOW=true.

mvn_delta_eps <- 1e-3
mvn_delta_draws <- 20000L

# The SD of the MVN draws per result row, computed here from the draws
# (`get_draws()` keeps each draw's value in `draw`).
mvn_draw_sd <- function(result, keys) {
  draws <- as.data.frame(get_draws(result))
  sds <- tapply(draws$draw, grid_key(draws, keys), stats::sd)
  unname(sds[grid_key(result, keys)])
}

expect_mvn_matches_delta <- function(point, covariance, keys, label) {
  mvn <- inferences(
    point,
    method = "mvn",
    vcov = covariance,
    n_draws = mvn_delta_draws,
    seed = 20261111,
    return_draws = TRUE
  )
  delta <- inferences(point, method = "delta", vcov = covariance)
  expect_identical(attr(delta, "target"), "empirical", label = label)

  ratio <- mvn_draw_sd(mvn, keys) / delta$std.error
  expect_false(anyNA(ratio), label = label)
  expect_lt(max(abs(ratio - 1)), 0.02, label = label)
}

test_that("MVN draws match conditional delta variance with a scaled covariance", {
  skip_unless_slow_tier()
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("mvtnorm")

  set_ups <- list(
    list(cell = grid_cell("orm"), estimand = "avg_sops"),
    list(cell = grid_cell("vglm"), estimand = "avg_sops"),
    list(cell = grid_cell("vglm"), estimand = "avg_comparisons")
  )
  for (set_up in set_ups) {
    fit <- grid_fit(set_up$cell)
    label <- paste(grid_label(set_up$cell), set_up$estimand)
    profiles <- utils::head(attr(fit, "markov_starting_profile_data"), 40L)
    covariance <- mvn_delta_eps * fit$var

    if (set_up$estimand == "avg_sops") {
      point <- avg_sops(
        fit,
        newdata = profiles,
        variables = list(tx = c(0, 1)),
        times = 1:4
      )
      keys <- c("time", "state", "tx")
    } else {
      point <- avg_comparisons(
        fit,
        newdata = profiles,
        variables = list(tx = c(0, 1)),
        times = 1:4
      )
      keys <- c("time", "state_set")
    }
    expect_mvn_matches_delta(point, covariance, keys, label)
  }
})
