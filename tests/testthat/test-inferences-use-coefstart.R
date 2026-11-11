# `use_coefstart = TRUE` starts each refit from the original coefficients. Only
# vglm refits take starting coefficients, so on other fits the argument would
# do nothing; it is an error there instead.

test_that("use_coefstart = TRUE is refused for orm fits", {
  case <- make_comparisons_case()
  orm_fit <- mostm(
    y ~ time + tx + marker + yprev,
    data = case$data,
    first_followup_time = 1,
    id_var = "id"
  )
  estimates <- avg_sops(orm_fit, times = 1:2, absorb = 3)
  contrast <- avg_comparisons(
    orm_fit,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    absorb = 3,
    state_sets = "1"
  )
  for (method in c("bootstrap", "fwb")) {
    for (x in list(estimates, contrast)) {
      err <- tryCatch(
        inferences(x, method = method, n_draws = 2, seed = 1, use_coefstart = TRUE),
        error = function(e) e
      )
      expect_s3_class(err, "error")
      if (inherits(err, "error")) {
        expect_match(conditionMessage(err), "`use_coefstart`", fixed = TRUE)
        expect_match(conditionMessage(err), "vglm")
      }
    }
  }
})

test_that("use_coefstart = TRUE still works for vglm fits", {
  case <- make_comparisons_case()
  estimates <- avg_sops(case$model, times = 1:2, absorb = 3)
  expect_no_error(
    inferences(estimates, method = "bootstrap", n_draws = 2, seed = 1, use_coefstart = TRUE)
  )
})

test_that("use_coefstart = TRUE is refused by bootstrap_model_coefs() for orm fits", {
  case <- make_comparisons_case()
  orm_fit <- mostm(
    y ~ time + tx + marker + yprev,
    data = case$data,
    first_followup_time = 1,
    id_var = "id"
  )
  expect_error(
    bootstrap_model_coefs(orm_fit, n_boot = 2, use_coefstart = TRUE),
    "`use_coefstart` can be `TRUE` only for models fitted with `engine = \"vglm\"`",
    fixed = TRUE
  )
})

# The simulation and analytical methods never refit, so starting values would
# have no effect there either.
test_that("use_coefstart = TRUE is refused for methods that don't refit", {
  case <- make_comparisons_case()
  estimates <- avg_sops(case$model, times = 1:2, absorb = 3)
  for (method in c("mvn", "delta")) {
    expect_error(
      inferences(estimates, method = method, n_draws = 2, use_coefstart = TRUE),
      "does not refit it",
      fixed = TRUE,
      info = method
    )
  }
  expect_error(
    inferences(estimates, method = "bootstrap", n_draws = 2, use_coefstart = NA),
    "`use_coefstart` must be `TRUE` or `FALSE`.",
    fixed = TRUE
  )
})
