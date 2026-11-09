# Only rmsb::blrm() fits can have random effects, so `include_re = TRUE` on any
# other fit is an error rather than being ignored.

expect_include_re_error <- function(expr, info) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "error")
  if (!inherits(err, "error")) {
    return(invisible())
  }
  message <- conditionMessage(err)
  expect_match(message, "`include_re`", fixed = TRUE, info = info)
  expect_match(message, "blrm", info = info)
}

test_that("include_re = TRUE is refused for vglm and orm fits", {
  case <- make_comparisons_case()
  orm_fit <- mostm(
    y ~ time + tx + marker + yprev,
    data = case$data,
    first_followup_time = 1,
    id_var = "id"
  )
  fits <- list(vglm = case$model, orm = orm_fit)
  for (engine in names(fits)) {
    model <- fits[[engine]]
    args <- list(model, times = 1:2, absorb = 3, include_re = TRUE)
    calls <- list(
      sops = function() do.call(sops, args),
      avg_sops = function() do.call(avg_sops, args),
      state_time = function() do.call(state_time, c(args, list(state_sets = "1"))),
      avg_state_time = function() {
        do.call(avg_state_time, c(args, list(state_sets = "1")))
      },
      comparisons = function() {
        do.call(
          comparisons,
          c(args, list(variables = list(tx = c(0, 1)), state_sets = "1"))
        )
      },
      avg_comparisons = function() {
        do.call(
          avg_comparisons,
          c(args, list(variables = list(tx = c(0, 1)), state_sets = "1"))
        )
      }
    )
    for (name in names(calls)) {
      expect_include_re_error(calls[[name]](), paste(engine, name))
    }
    # The default, FALSE, is fine.
    expect_no_error(avg_sops(model, times = 1:2, absorb = 3))
  }
})
