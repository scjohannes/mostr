# Unconditional variance is refused for every fit with non-unit case weights:
# whether the formula holds with weights has not been checked, so no engine
# reports it. Conditional variance stays available.

# Made once per test run.
make_weighted_case <- function() {
  case <- make_comparisons_case()
  grid_fixture("weighted comparisons case", function() build_weighted_case(case))
}

build_weighted_case <- function(case) {
  data <- case$data
  ids <- unique(data$id)
  patient_weight <- (ids %% 3) + 1
  data$w <- patient_weight[match(data$id, ids)]
  vglm_fit <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm",
    weights = w
  )
  orm_fit <- suppressWarnings(mostm(
    y ~ time + tx + marker + yprev,
    data = data,
    first_followup_time = 1,
    id_var = "id",
    weights = w
  ))
  list(vglm = vglm_fit, orm = orm_fit)
}

expect_weighted_unconditional_error <- function(expr, info) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "error")
  if (!inherits(err, "error")) {
    return(invisible())
  }
  message <- conditionMessage(err)
  expect_match(message, "weight", info = info)
  expect_match(message, "unconditional", info = info)
  expect_match(message, "conditional\"", fixed = TRUE, info = info)
}

test_that("unconditional variance is refused for weighted fits of every engine", {
  fits <- make_weighted_case()
  for (engine in names(fits)) {
    estimates <- avg_sops(fits[[engine]], times = 1:2, absorb = 3)
    expect_weighted_unconditional_error(
      inferences(estimates, method = "delta", vcov = "unconditional"),
      paste(engine, "explicit")
    )
    # Unconditional variance is the default for averages.
    expect_weighted_unconditional_error(
      inferences(estimates, method = "delta"),
      paste(engine, "default")
    )
    contrast <- avg_comparisons(
      fits[[engine]],
      variables = list(tx = c(0, 1)),
      times = 1:2,
      absorb = 3,
      state_sets = "1"
    )
    expect_weighted_unconditional_error(
      inferences(contrast, method = "delta", vcov = "unconditional"),
      paste(engine, "avg_comparisons")
    )
  }
})

test_that("conditional variance still works for weighted fits", {
  fits <- make_weighted_case()
  for (engine in names(fits)) {
    out <- inferences(
      avg_sops(fits[[engine]], times = 1:2, absorb = 3),
      method = "delta",
      vcov = "conditional"
    )
    expect_true(all(is.finite(out$std.error)), info = engine)
  }
})

test_that("inferences() no longer offers the score bootstrap", {
  case <- make_comparisons_case()
  estimates <- avg_sops(case$model, times = 1:2, absorb = 3)
  expect_error(
    inferences(estimates, method = "score_bootstrap", n_draws = 10, seed = 1),
    "score_bootstrap"
  )
  expect_false("score_bootstrap" %in% eval(formals(inferences)$method))
})
