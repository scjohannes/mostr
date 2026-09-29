fit_reverse_vglm <- function(data, family) {
  vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = family,
    data = data,
    id_var = "id"
  )
}

avg_sops_reverse_case <- function(fit) {
  avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    y_levels = 1:6,
    absorb = 6
  )
}

test_that("the reverse check accepts `reverse` set from a local variable", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5401)
  fit_in_closure <- function(data, rev_flag) {
    vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = rev_flag, parallel = TRUE),
      data = data,
      id_var = "id"
    )
  }
  fit <- fit_in_closure(data, rev_flag = TRUE)

  expect_no_error(mostr:::validate_markov_model(fit))
  expect_s3_class(avg_sops_reverse_case(fit), "data.frame")
})

test_that("the reverse check accepts a family object passed by name", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5402)
  fam <- VGAM::cumulative(reverse = TRUE, parallel = TRUE)
  fit <- fit_reverse_vglm(data, fam)

  expect_no_error(mostr:::validate_markov_model(fit))
  expect_s3_class(avg_sops_reverse_case(fit), "data.frame")
})

test_that("the reverse check uses the fitted family, not later variable values", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5403)
  rev_flag <- TRUE
  fit <- vglm_markov(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = rev_flag, parallel = TRUE),
    data = data,
    id_var = "id"
  )
  rev_flag <- FALSE

  expect_no_error(mostr:::validate_markov_model(fit))
})

test_that("the reverse check rejects reverse = FALSE family objects", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5404)
  fam <- VGAM::cumulative(reverse = FALSE, parallel = TRUE)
  fit <- suppressWarnings(fit_reverse_vglm(data, fam))

  expect_error(
    mostr:::validate_markov_model(fit),
    "must use reverse = TRUE",
    fixed = TRUE
  )
})

test_that("bootstrap refits pass the reverse check for family objects", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5405)
  fit_in_closure <- function(data, rev_flag) {
    fam <- VGAM::cumulative(reverse = rev_flag, parallel = TRUE)
    fit_reverse_vglm(data, fam)
  }
  fit <- fit_in_closure(data, rev_flag = TRUE)
  out <- avg_sops_reverse_case(fit)

  withr::local_seed(5406)
  inferred <- inferences(out, method = "bootstrap", n_draws = 2)
  expect_equal(attr(inferred, "n_successful"), 2L)
})
