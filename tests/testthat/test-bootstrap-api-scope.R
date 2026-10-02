# Refit bootstraps always resample the data stored on a wrapper fit. The public
# API therefore does not accept separate refit data, and the low-level
# bootstrap helpers are internal.

test_that("SOP functions do not accept separate refit data", {
  for (fun in list(sops, avg_sops, avg_comparisons, avg_time)) {
    expect_false("refit_data" %in% names(formals(fun)))
  }
})

test_that("bootstrap_model_coefs resamples the stored refit data only", {
  expect_false("data" %in% names(formals(bootstrap_model_coefs)))
})

test_that("low-level bootstrap helpers are not exported", {
  # `load_all()` exports every function, so read the declared exports.
  namespace_lines <- readLines(system.file("NAMESPACE", package = "mostr"))
  exports <- sub(
    "^export\\((.*)\\)$",
    "\\1",
    grep("^export\\(", namespace_lines, value = TRUE)
  )
  helpers <- c(
    "bootstrap_analysis_wrapper",
    "apply_to_bootstrap",
    "fast_group_bootstrap",
    "materialize_bootstrap_sample_indexed",
    "relevel_factors_consecutive"
  )
  expect_length(intersect(helpers, exports), 0L)
  expect_true(all(c("inferences", "bootstrap_model_coefs") %in% exports))
})
