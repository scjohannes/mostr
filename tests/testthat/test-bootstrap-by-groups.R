test_that("ordinary bootstrap stops when a resample lacks a `by` group", {
  case <- make_comparisons_case()
  # With 60 patients, a resample misses both patients of the `rare` group
  # with probability (58/60)^60, about 0.13, so 40 resamples almost surely
  # include such a resample. The check runs before any refit.
  expect_error(
    inferences(
      avg_sops(case$model, by = "rare", times = 1:3, absorb = 3),
      method = "bootstrap",
      n_draws = 40,
      seed = 11
    ),
    "rare = TRUE.*method = \"fwb\""
  )
  expect_error(
    inferences(
      avg_state_time(case$model, by = "rare", times = 1:3, absorb = 3),
      method = "bootstrap",
      n_draws = 40,
      seed = 11
    ),
    "rare = TRUE"
  )
  for (estimand in c("sop", "time_benefit")) {
    expect_error(
      inferences(
        avg_comparisons(
          case$model,
          variables = list(tx = c(0, 1)),
          estimand = estimand,
          by = "rare",
          times = 1:3,
          absorb = 3
        ),
        method = "bootstrap",
        n_draws = 40,
        seed = 11
      ),
      "rare = TRUE"
    )
  }
})

test_that("ordinary bootstrap keeps groups that every resample contains", {
  case <- make_comparisons_case()
  boot <- inferences(
    avg_sops(case$model, by = "marker", times = 1:2, absorb = 3),
    method = "bootstrap",
    n_draws = 3,
    seed = 12
  )
  expect_true(all(is.finite(boot$std.error)))

  # Fractional weighted bootstrap keeps every patient, so small groups stay.
  fwb <- inferences(
    avg_sops(case$model, by = "rare", times = 1:2, absorb = 3),
    method = "fwb",
    n_draws = 2,
    seed = 13
  )
  expect_true(all(is.finite(fwb$std.error)))
})

test_that("group check reports each missing group and its resample count", {
  baseline <- data.frame(
    id = c("a", "b", "c", "d"),
    sex = c("f", "f", "m", NA),
    site = c(1, 2, 1, 1)
  )
  resamples <- list(
    data.frame(original_id = c("a", "a", "b", "d")),
    data.frame(original_id = c("a", "c", "c", "d")),
    data.frame(original_id = c("b", "b", "b", "b"))
  )
  expect_silent(check_bootstrap_by_groups(resamples[2], baseline, "id", "sex"))
  expect_error(
    check_bootstrap_by_groups(resamples, baseline, "id", "sex"),
    "`sex = m` in 2 of 3 resamples"
  )
  expect_error(
    check_bootstrap_by_groups(resamples, baseline, "id", c("sex", "site")),
    "`sex = f, site = 1`"
  )
  expect_null(check_bootstrap_by_groups(resamples, baseline, "id", NULL))
})
