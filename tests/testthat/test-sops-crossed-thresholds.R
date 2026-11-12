# In a partial proportional-odds model the cut-points can cross for some
# patients, so the model implies a negative probability for a state. SOPs are
# then not defined, and every route stops with an error instead of clipping the
# negative probability to zero and rescaling the rest.

# A partial proportional-odds vglm fit with a non-proportional effect of a
# continuous `score`. At score = -20 the fitted cut-points cross (state 2 gets
# a predicted probability of about -0.06); at score = 0 they don't. Made once
# per test run.
make_crossed_threshold_case <- function() {
  case <- make_comparisons_case()
  grid_fixture("crossed threshold case", function() {
    build_crossed_threshold_case(case)
  })
}

build_crossed_threshold_case <- function(case) {
  data <- case$data
  withr::local_seed(1)
  data$score <- stats::rnorm(nrow(data))
  model <- mostm(
    ordered(y) ~ time + tx + score + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ score),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    engine = "vglm"
  )
  profiles <- data[data$time == 1, , drop = FALSE][1:3, , drop = FALSE]
  profiles$score <- c(0, 1, -20)
  list(model = model, data = data, profiles = profiles)
}

expect_crossed_threshold_error <- function(expr, info) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "error")
  if (!inherits(err, "error")) {
    return(invisible())
  }
  expect_match(conditionMessage(err), "cross", info = info)
  expect_match(conditionMessage(err), "negative", info = info)
}

test_that("the fixture's thresholds really cross at score = -20", {
  skip_on_cran() # The pass depends on the seed.
  case <- make_crossed_threshold_case()
  fit <- case$model$vglm_fit
  probabilities <- VGAM::predict(fit, newdata = case$profiles, type = "response")
  expect_true(all(probabilities[1:2, ] >= 0))
  expect_lt(min(probabilities[3, ]), 0)
})

test_that("SOP calls stop when the fitted thresholds cross", {
  skip_on_cran() # The pass depends on the seed.
  case <- make_crossed_threshold_case()
  base <- list(case$model, newdata = case$profiles, times = 1:3, absorb = 3)
  compare <- c(base, list(variables = list(tx = c(0, 1)), state_sets = "1"))
  calls <- list(
    sops = function() do.call(sops, base),
    avg_sops = function() do.call(avg_sops, base),
    state_time = function() do.call(state_time, c(base, list(state_sets = "1"))),
    comparisons = function() do.call(comparisons, compare),
    avg_comparisons = function() do.call(avg_comparisons, compare)
  )
  for (name in names(calls)) {
    expect_crossed_threshold_error(calls[[name]](), name)
  }
})

test_that("the error says how many profiles cross, and when", {
  skip_on_cran() # The pass depends on the seed.
  case <- make_crossed_threshold_case()
  err <- tryCatch(
    sops(case$model, newdata = case$profiles, times = 1:3, absorb = 3),
    error = function(e) e
  )
  expect_s3_class(err, "mostr_crossed_thresholds")
  expect_match(conditionMessage(err), "for 1 of the 3 profiles at time 1")
  expect_match(conditionMessage(err), "outside the range of the data")

  # Comparisons predict each profile once per compared value.
  expect_error(
    comparisons(
      case$model,
      newdata = case$profiles,
      times = 1:3,
      absorb = 3,
      variables = list(tx = c(0, 1)),
      state_sets = "1"
    ),
    "for 2 of the 6 profiles at time 1",
    class = "mostr_crossed_thresholds"
  )
})

test_that("sops() agrees with the readable loop on both sides of the crossing", {
  skip_on_cran() # The pass depends on the seed.
  case <- make_crossed_threshold_case()
  predict_transition <- slow_predictor(case$model)
  loop <- function(i) {
    slow_first_order_sops(
      predict_transition,
      case$profiles[i, , drop = FALSE],
      times = 1:3,
      states = 1:3,
      absorb = 3
    )
  }

  out <- sops(
    case$model,
    newdata = case$profiles[1:2, , drop = FALSE],
    times = 1:3,
    absorb = 3
  )
  for (i in 1:2) {
    rows <- out[out$rowid == i, , drop = FALSE]
    expected <- loop(i)[cbind(rows$time, as.integer(as.character(rows$state)))]
    expect_equal(rows$estimate, expected, tolerance = 1e-8)
  }
  expect_error(loop(3L), "cut-points cross")
})

test_that("the R engine, which predicts with VGAM, also stops", {
  skip_on_cran() # The pass depends on the seed.
  case <- make_crossed_threshold_case()
  expect_error(
    soprob_markov_r_engine(
      case$model,
      newdata = case$profiles,
      times = 1:3,
      y_levels = 1:3,
      absorb = 3
    ),
    "cross for 1 of the 3 profiles at time 1",
    class = "mostr_crossed_thresholds"
  )
})

test_that("rounding negatives are set to zero instead of stopping", {
  # Nearly equal cut-points: the middle probability is about -2e-13.
  eta <- matrix(c(0.3, 0.3 + 1e-12), nrow = 1L)
  raw <- c(1 - stats::plogis(eta[1L]), -diff(stats::plogis(eta[1L, ])))
  expect_lt(raw[2L], 0)
  probabilities <- lp_to_probs(eta, 2L)
  expect_true(all(probabilities >= 0))
  expect_equal(sum(probabilities), 1)
  expect_no_error(check_crossed_thresholds(probabilities, 1L, 1))
  expect_no_error(
    check_crossed_thresholds(matrix(c(0.5, -1e-11, 0.5), nrow = 1L), 1L, 1)
  )
  expect_error(
    check_crossed_thresholds(matrix(c(0.5, -1e-9, 0.5), nrow = 1L), 1L, 1),
    class = "mostr_crossed_thresholds"
  )
})

# At score = -6 the estimated cut-points don't cross, but they do for some
# coefficient draws around the estimate and for some bootstrap refits.
test_that("inference stops when cut-points cross in some draws or refits", {
  skip_on_cran() # The pass depends on the seed.
  case <- make_crossed_threshold_case()
  profiles <- case$profiles[1:2, , drop = FALSE]
  profiles$score[2L] <- -6
  cases <- list(
    list(
      label = "sops(), simulation draws",
      result = sops(case$model, newdata = profiles, times = 1:3, absorb = 3),
      method = "mvn",
      n_draws = 20,
      seed = 1,
      phrase = "cross in [0-9]+ of the 20 simulation draws"
    ),
    list(
      label = "avg_comparisons(), simulation draws",
      result = avg_comparisons(
        case$model,
        newdata = profiles,
        times = 1:3,
        absorb = 3,
        variables = list(tx = c(0, 1)),
        state_sets = "1"
      ),
      method = "mvn",
      n_draws = 20,
      seed = 2,
      phrase = "simulation draws"
    ),
    list(
      label = "avg_sops(), bootstrap refits",
      result = avg_sops(case$model, newdata = profiles, times = 1:3, absorb = 3),
      method = "bootstrap",
      n_draws = 10,
      seed = 3,
      phrase = "cross in [0-9]+ of the 10 bootstrap samples"
    )
  )
  for (row in cases) {
    withr::with_seed(row$seed, {
      expect_error(
        inferences(row$result, method = row$method, n_draws = row$n_draws),
        row$phrase,
        class = "mostr_crossed_thresholds",
        info = row$label
      )
    })
  }
})

# A hand-made blrm fit with a partial proportional-odds effect of time,
# P(Y >= k) = plogis(alpha_k + x * beta + time * tau * (k - 1)), as in
# test-sops-native-blrm-slow-loop.R. With tau = 0.5 the first two cut-points
# cross once time * tau > 1.1, that is from time 3 on.
crossed_fake_blrm <- function() {
  draws <- rbind(
    c(1.2, 0.1, -1.3, 0.40, 0.08, -0.60, 0.05),
    c(1.2, 0.1, -1.3, 0.10, 0.03, -0.45, 0.50),
    c(1.2, 0.1, -1.3, 0.70, 0.12, -0.75, 0.50),
    c(1.2, 0.1, -1.3, 0.25, 0.05, -0.55, 0.02)
  )
  colnames(draws) <- c("y>=2", "y>=3", "y>=4", "tx", "time", "yprev", "tau_time")
  fake_blrm(draws, tau = "tau_time", cppo = function(y) as.numeric(y) - 1)
}

crossed_fake_design <- fake_blrm_design(c("tx", "time", "yprev"), tau_var = "time")

test_that("blrm partial proportional-odds fits stop where cut-points cross", {
  local_mocked_bindings(blrm_design_matrix = crossed_fake_design)
  model <- crossed_fake_blrm()
  profiles <- data.frame(
    id = c("a", "b"),
    time = 1,
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:4)
  )

  before <- sops(model, newdata = profiles, times = 1:2, absorb = 4L, n_draws = NULL)
  expect_true(all(before$estimate >= 0))

  expect_error(
    sops(model, newdata = profiles, times = 1:3, absorb = 4L, n_draws = NULL),
    "cross for 2 of the 2 profiles at time 3 in 2 of the 4 posterior draws",
    class = "mostr_crossed_thresholds"
  )
})
