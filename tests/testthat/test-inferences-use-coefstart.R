# `use_coefstart = TRUE` starts each refit from the original coefficients. Only
# vglm refits take starting coefficients, so on other fits the argument would
# do nothing; it is an error there instead.

test_that("use_coefstart = TRUE still works for vglm fits", {
  case <- make_comparisons_case()
  estimates <- avg_sops(case$model, times = 1:2)
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
    id_var = "id",
    absorb = 3
  )
  expect_error(
    bootstrap_model_coefs(orm_fit, n_boot = 2, use_coefstart = TRUE),
    "`use_coefstart` can be `TRUE` only for models fitted with `engine = \"vglm\"`",
    fixed = TRUE
  )
})
