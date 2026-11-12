# A score-bootstrap fit and its baseline rows, made once for the draw-generator
# tests below.
score_draw_case <- function() {
  grid_fixture("score bootstrap draw case", function() {
    make_score_bootstrap_case(seed = 2006, n_patients = 35, follow_up_time = 8)
  })
}

test_that("score-bootstrap simulation adds CIs, metadata and draws for avg_sops", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2001)
  withr::local_seed(2001)
  result <- avg_sops(
    model = case$model,
    variables = list(tx = c(0, 1)),
    times = case$times,
    absorb = case$absorb
  ) |>
    mostr:::inferences_score_bootstrap(
      n_draws = 40,
      return_draws = TRUE
    )

  expect_equal(attr(result, "method"), "score_bootstrap")
  expect_equal(attr(result, "n_draws"), 40)
  expect_true(attr(result, "n_successful") <= 40)
  expect_true(attr(result, "n_successful") > 0)
  expect_inference_intervals(result, require_positive_std_error = TRUE)

  draws <- get_draws(result)

  expect_s3_class(draws, "data.frame")
  expect_contains(names(draws), c("draw_id", "estimate", "time", "state", "tx"))
  expect_equal(length(unique(draws$draw_id)), attr(result, "n_successful"))
  expect_all_true(is.finite(draws$draw_id))
  expect_all_true(is.finite(draws$estimate))
})

test_that("score-bootstrap treats supplied newdata as fixed profiles", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2010)
  profiles <- case$baseline[seq_len(5), , drop = FALSE]
  profiles <- profiles[, setdiff(names(profiles), "id"), drop = FALSE]

  avg <- avg_sops(
    model = case$model,
    newdata = profiles,
    variables = list(tx = c(0, 1)),
    times = case$times[1:2],
    absorb = case$absorb
  )

  withr::local_seed(2010)
  expect_warning(
    result <- mostr:::inferences_score_bootstrap(
      avg,
      n_draws = 3,
      return_draws = TRUE
    ),
    "baseline weights have been set to NULL"
  )

  draws <- get_draws(result)

  expect_equal(attr(result, "n_successful"), 3L)
  expect_equal(nrow(attr(avg, "newdata_pred")), nrow(profiles) * 2L)
  expect_contains(names(draws), c("draw_id", "estimate", "time", "state", "tx"))
  expect_false("id" %in% names(draws))
})

test_that("generate_score_bootstrap_draws returns valid dimensions and normalized weights", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- score_draw_case()

  draws <- generate_score_bootstrap_draws(
    model = case$model,
    baseline_data = case$baseline,
    id_var = "id",
    n_draws = 25
  )

  expect_true(is.matrix(draws$beta_draws))
  expect_true(is.matrix(draws$baseline_weights))
  expect_equal(nrow(draws$beta_draws), 25)
  expect_equal(ncol(draws$beta_draws), length(case$model$coefficients))
  expect_equal(nrow(draws$baseline_weights), 25)
  expect_equal(ncol(draws$baseline_weights), nrow(case$baseline))
  expect_equal(
    as.numeric(rowSums(draws$baseline_weights)),
    rep(1, 25),
    tolerance = 1e-10
  )
  expect_true(all(is.finite(draws$beta_draws)))
  expect_true(all(is.finite(draws$baseline_weights)))
  expect_true(all(draws$baseline_weights >= 0))
})

test_that("generate_score_bootstrap_draws is reproducible with a fixed seed", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- score_draw_case()

  withr::local_seed(9999)
  draws_a <- generate_score_bootstrap_draws(
    model = case$model,
    baseline_data = case$baseline,
    id_var = "id",
    n_draws = 15
  )

  withr::local_seed(9999)
  draws_b <- generate_score_bootstrap_draws(
    model = case$model,
    baseline_data = case$baseline,
    id_var = "id",
    n_draws = 15
  )

  expect_equal(draws_a$beta_draws, draws_b$beta_draws)
  expect_equal(draws_a$baseline_weights, draws_b$baseline_weights)
})

test_that("score bootstrap draws are invariant to working-memory chunks", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- score_draw_case()

  withr::local_options(mostr.score_bootstrap_working_bytes = 1)
  withr::local_seed(777)
  small_chunks <- generate_score_bootstrap_draws(
    model = case$model,
    baseline_data = case$baseline,
    id_var = "id",
    n_draws = 17
  )

  withr::local_options(mostr.score_bootstrap_working_bytes = Inf)
  withr::local_seed(777)
  one_chunk <- generate_score_bootstrap_draws(
    model = case$model,
    baseline_data = case$baseline,
    id_var = "id",
    n_draws = 17
  )

  expect_equal(small_chunks$beta_draws, one_chunk$beta_draws, tolerance = 1e-13)
  expect_identical(small_chunks$baseline_weights, one_chunk$baseline_weights)
})

test_that("generate_score_bootstrap_draws rejects bare fits, unknown baseline ids and missing robust components", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- score_draw_case()
  unknown_ids <- case$baseline
  unknown_ids$id <- as.character(unknown_ids$id)
  unknown_ids$id[1] <- "missing_cluster_id"
  no_scores <- case$model
  no_scores$scores <- NULL

  cases <- list(
    # `mostm(engine = "vglm")` always returns a robcov_vglm object, so only
    # the bare vglm fit inside it lacks the patient scores.
    list(
      label = "bare vglm fit",
      model = case$model$vglm_fit,
      baseline = case$baseline,
      error = "requires a 'robcov_vglm' model"
    ),
    list(
      label = "baseline ids not among the clusters",
      model = case$model,
      baseline = unknown_ids,
      error = "were not found in the cluster variable"
    ),
    list(
      label = "missing scores",
      model = no_scores,
      baseline = case$baseline,
      error = "missing required components.*scores"
    )
  )

  for (row in cases) {
    withr::local_seed(2008)
    expect_error(
      generate_score_bootstrap_draws(
        model = row$model,
        baseline_data = row$baseline,
        id_var = "id",
        n_draws = 10
      ),
      row$error,
      info = row$label
    )
  }
})
