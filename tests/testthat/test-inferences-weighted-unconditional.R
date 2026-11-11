# Unconditional variance is refused for every fit with non-unit case weights:
# whether the formula holds with weights has not been checked, so no engine
# reports it (rows I36 to I38 of test-refusals.R). Conditional variance stays
# available.

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
