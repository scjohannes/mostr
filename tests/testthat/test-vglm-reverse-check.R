fit_reverse_vglm <- function(data, family) {
  mostm(
    ordered(y) ~ time + tx + yprev,
    family = family,
    data = data,
    id_var = "id",
    engine = "vglm"
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

test_that("the reverse check reads `reverse` from the fitted family however it was passed", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5401)

  # Each `fit()` returns a vglm fit whose family was spelled differently.
  cases <- list(
    list(
      label = "`reverse` set from a local variable",
      fit = function() {
        fit_in_closure <- function(data, rev_flag) {
          mostm(
            ordered(y) ~ time + tx + yprev,
            family = VGAM::cumulative(reverse = rev_flag, parallel = TRUE),
            data = data,
            id_var = "id",
            engine = "vglm"
          )
        }
        fit_in_closure(data, rev_flag = TRUE)
      },
      error = NULL
    ),
    list(
      label = "a family object passed by name",
      fit = function() {
        fam <- VGAM::cumulative(reverse = TRUE, parallel = TRUE)
        fit_reverse_vglm(data, fam)
      },
      error = NULL
    ),
    list(
      label = "the variable changes after the fit",
      fit = function() {
        rev_flag <- TRUE
        fit <- mostm(
          ordered(y) ~ time + tx + yprev,
          family = VGAM::cumulative(reverse = rev_flag, parallel = TRUE),
          data = data,
          id_var = "id",
          engine = "vglm"
        )
        rev_flag <- FALSE
        fit
      },
      error = NULL
    ),
    list(
      label = "a reverse = FALSE family object",
      fit = function() {
        fam <- VGAM::cumulative(reverse = FALSE, parallel = TRUE)
        suppressWarnings(fit_reverse_vglm(data, fam))
      },
      error = "must use reverse = TRUE"
    )
  )

  for (case in cases) {
    fit <- case$fit()
    if (is.null(case$error)) {
      expect_error(mostr:::validate_markov_model(fit), NA, info = case$label)
    } else {
      expect_error(
        mostr:::validate_markov_model(fit),
        case$error,
        fixed = TRUE,
        info = case$label
      )
    }
  }
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
