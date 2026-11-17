# Scores of a constrained partial proportional-odds blrm fit. rmsb::blrm()
# computes them once on all state levels, `cppo(ylevels)`, and its Stan model
# uses element j for the cut-point "y >= state j" (pposcore[2..k]). Calling
# `cppo` on one state at a time gives other scores for a `cppo` that uses the
# other states' values. On a real 10-state fit with `cppo = rank`, those
# scores missed Stan's own log-likelihood by up to 2.7, while
# `cppo(ylevels)[-1]` matched it to 5e-14 (issue #46).

# Rescales the states to (0, 1] by the highest one, so the score of each state
# depends on the other states: called on one state, it always returns 1.
cppo_by_highest <- function(y) as.numeric(y) / max(as.numeric(y))

partial_odds_fake_blrm <- function(cppo) {
  draws <- rbind(
    c(1.2, 0.1, -1.3, 0.40, 0.08, -0.60, 0.50),
    c(0.9, -0.2, -1.1, 0.10, 0.03, -0.45, -0.80)
  )
  colnames(draws) <- c("y>=2", "y>=3", "y>=4", "tx", "time", "yprev", "tau_time")
  fake_blrm(draws, tau = "tau_time", cppo = cppo, cluster = FALSE)
}

partial_odds_rows <- data.frame(
  id = c("a", "b"),
  time = c(1, 3),
  tx = c(0, 1),
  yprev = factor(c(1, 2), levels = 1:4)
)

test_that("blrm cppo scores are computed on all states at once, as rmsb fits them", {
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time", "yprev"), tau_var = "time")
  )
  fit <- partial_odds_fake_blrm(cppo_by_highest)
  probs <- predict_blrm_response_markov(fit, partial_odds_rows)

  # P(Y >= k) = plogis(alpha_k + x * beta + time * tau * score_k), with the
  # scores of states 2 to 4 out of 1 to 4: 0.5, 0.75 and 1.
  scores <- cppo_by_highest(1:4)[-1]
  for (draw in 1:2) {
    coefficient <- fit$draws[draw, ]
    for (i in 1:2) {
      row <- partial_odds_rows[i, ]
      linear <- row$tx * coefficient[["tx"]] +
        row$time * coefficient[["time"]] +
        as.numeric(as.character(row$yprev)) * coefficient[["yprev"]]
      at_least <- stats::plogis(
        coefficient[1:3] + linear + row$time * coefficient[["tau_time"]] * scores
      )
      expected <- -diff(c(1, at_least, 0))
      expect_equal(unname(probs[draw, i, ]), unname(expected), tolerance = 1e-12)
    }
  }
})

test_that("a cppo that returns no number for each state is refused", {
  local_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(c("tx", "time", "yprev"), tau_var = "time")
  )
  one_score <- partial_odds_fake_blrm(function(y) 1)
  expect_error(
    predict_blrm_response_markov(one_score, partial_odds_rows),
    class = "mostr_error"
  )
  text_scores <- partial_odds_fake_blrm(function(y) as.character(y))
  expect_error(
    predict_blrm_response_markov(text_scores, partial_odds_rows),
    class = "mostr_error"
  )
})
