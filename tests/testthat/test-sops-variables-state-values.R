# State labels are names. When `variables` sets the previous state (`yprev`,
# or `ypprev` in a second-order model), a value is the state whose label it
# matches, whatever its type: 2, 2L and "2" all mean state "2", whether the
# fitted model stores the previous state as a factor or as a number. A value
# that matches no state is an error that names the column and the states.

numeric_yprev_model <- function(case) {
  data <- case$data
  data$yprev <- as.integer(as.character(data$yprev))
  mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
  )
}

test_that("numbers set a factor previous state by its labels", {
  case <- make_comparisons_case()
  levels <- levels(case$profiles$yprev)
  by_factor <- avg_sops(
    case$model,
    variables = list(yprev = factor(c("1", "2"), levels = levels)),
    times = 1:3,
    absorb = 3
  )
  for (values in list(1:2, c(1, 2), c("1", "2"))) {
    out <- avg_sops(case$model, variables = list(yprev = values), times = 1:3, absorb = 3)
    expect_equal(out$estimate, by_factor$estimate, tolerance = 1e-12, info = class(values))
  }
  contrast <- avg_comparisons(case$model, variables = list(yprev = 1:2), times = 1:3, absorb = 3)
  expected <- avg_comparisons(
    case$model,
    variables = list(yprev = factor(c("1", "2"), levels = levels)),
    times = 1:3,
    absorb = 3
  )
  expect_equal(contrast$estimate, expected$estimate, tolerance = 1e-12)
})

test_that("character labels set a numeric previous state", {
  skip_on_cran()
  case <- make_comparisons_case()
  model <- numeric_yprev_model(case)
  by_number <- avg_sops(model, variables = list(yprev = 1:2), times = 1:3, absorb = 3)
  by_label <- avg_sops(model, variables = list(yprev = c("1", "2")), times = 1:3, absorb = 3)
  expect_equal(by_label$estimate, by_number$estimate, tolerance = 1e-12)
  contrast <- avg_comparisons(model, variables = list(yprev = c("1", "2")), times = 1:3, absorb = 3)
  expected <- avg_comparisons(model, variables = list(yprev = 1:2), times = 1:3, absorb = 3)
  expect_equal(contrast$estimate, expected$estimate, tolerance = 1e-12)
})

test_that("the state two visits back is matched by label too", {
  case <- make_second_order_comparisons_case()
  levels <- levels(case$data$ypprev)
  by_factor <- avg_sops(
    case$model,
    variables = list(ypprev = factor(c("1", "2"), levels = levels)),
    times = 1:3,
    absorb = 3
  )
  by_number <- avg_sops(
    case$model,
    variables = list(ypprev = 1:2),
    times = 1:3,
    absorb = 3
  )
  expect_equal(by_number$estimate, by_factor$estimate, tolerance = 1e-12)
})

test_that("a previous-state value that is no state is an error naming the states", {
  case <- make_comparisons_case()
  for (values in list(c("1", "5"), c(1, 5), c(1, 1.5))) {
    expect_error(
      avg_sops(case$model, variables = list(yprev = values), times = 1:3, absorb = 3),
      "`yprev`.*not a state.*1, 2, 3",
      info = paste(values, collapse = " ")
    )
  }
  model <- numeric_yprev_model(case)
  expect_error(
    avg_sops(model, variables = list(yprev = c("1", "x")), times = 1:3, absorb = 3),
    "`yprev`.*not a state.*1, 2, 3"
  )
})
