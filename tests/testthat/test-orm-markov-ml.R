# orm_markov_ml(): observed-data maximum likelihood with missing states.

ml_fit <- function(data, formula = y ~ rms::rcs(time, 3) + tx + yprev, ...) {
  orm_markov_ml(
    formula,
    data = data,
    id_var = "id",
    absorbing_state = 6,
    ...
  )
}

# --- Likelihood and posteriors against brute-force enumeration --------------

test_that("the likelihood and posteriors equal enumeration over completions", {
  data <- ml_oracle_data()
  for (rule in c("any", "non_absorbing")) {
    problem <- ml_oracle_problem(data, missing_states = rule)
    theta0 <- stats::coef(problem$start_fit)
    for (seed in 1:3) {
      theta <- ml_random_theta(theta0, problem$n_int, seed)
      state <- markov_ml_estep(problem, theta)
      oracle <- ml_oracle(
        data,
        problem$start_fit,
        theta,
        missing_states = rule
      )
      expect_equal(state$loglik, oracle$loglik, tolerance = 1e-12)

      xi <- ml_problem_xi(problem, state)
      keys <- names(xi)
      expect_setequal(keys, names(oracle$xi))
      diff <- max(vapply(
        keys,
        function(k) max(abs(xi[[k]] - oracle$xi[[k]])),
        numeric(1)
      ))
      expect_lt(diff, 1e-12)
    }
  }
})

test_that("without an absorbing state the likelihood equals enumeration", {
  data <- ml_oracle_data_no_absorbing()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  # Without absorbing states, both rules let a missing state be any state.
  for (rule in c("any", "non_absorbing")) {
    fit <- orm_markov_ml(formula, data = data, id_var = "id",
      missing_states = rule)
    expect_length(fit$ml$absorbing_state, 0L)
    expect_gt(fit$ml$n_missing, 0)
    expect_true(fit$ml$converged)
    problem <- markov_ml_setup(formula, data, "id", "time", "yprev", NULL,
      rule, list())
    oracle <- ml_oracle(data, problem$start_fit, coef(fit), K = 3L,
      absorbing = integer(), missing_states = rule)
    expect_equal(fit$ml$loglik, oracle$loglik, tolerance = 1e-10)
  }
})

test_that("posteriors sum to one and match the forward-backward identity", {
  data <- ml_oracle_data()
  problem <- ml_oracle_problem(data)
  state <- markov_ml_estep(problem, stats::coef(problem$start_fit))
  xi <- ml_problem_xi(problem, state)
  expect_equal(unname(vapply(xi, sum, numeric(1))), rep(1, problem$n),
    tolerance = 1e-12
  )
})

# --- Equivalence with orm_markov() when nothing is missing -------------------

test_that("with nothing missing the fit reproduces orm_markov()", {
  data <- make_test_data(n_patients = 80, follow_up_time = 10)
  formula <- ordered(y) ~ rms::rcs(time, 3) + tx + yprev
  reference <- orm_markov(formula, data = data, id_var = "id")
  fit <- ml_fit(data, formula)

  expect_s3_class(fit, "orm_markov_ml")
  expect_s3_class(fit, "orm")
  expect_lte(fit$ml$iterations, 1L)
  expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
  expect_equal(fit$ml$loglik, -reference$deviance[[2L]] / 2, tolerance = 1e-10)
  expect_equal(
    fit$ml$var_model,
    reference$orig.var,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(fit$ml$var_robust, reference$var, tolerance = 1e-6)
  expect_equal(vcov(fit, intercepts = "all"), fit$ml$var_robust)
  expect_lt(fit$ml$max_abs_score, 1e-6)

  hc1 <- orm_markov(formula, data = data, id_var = "id", type = "HC1")
  fit_hc1 <- ml_fit(data, formula, type = "HC1")
  expect_equal(fit_hc1$var, hc1$var, tolerance = 1e-6)

  a_ref <- avg_sops(reference, variables = list(tx = c(0, 1)), times = 1:10,
    absorb = "6")
  a_ml <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:10,
    absorb = "6")
  expect_equal(a_ml$estimate, a_ref$estimate, tolerance = 1e-8)
  for (target in c("conditional", "unconditional")) {
    expect_equal(
      inferences(a_ml, method = "delta", vcov = target)$std.error,
      inferences(a_ref, method = "delta", vcov = target)$std.error,
      tolerance = 1e-6
    )
  }
  s_ref <- sops(reference, times = 1:4, absorb = "6")
  s_ml <- sops(fit, times = 1:4, absorb = "6")
  expect_equal(s_ml$estimate, s_ref$estimate, tolerance = 1e-8)
})

test_that("with nothing missing the fit costs no forward-backward passes", {
  data <- make_test_data(n_patients = 80, follow_up_time = 10)
  formula <- ordered(y) ~ rms::rcs(time, 3) + tx + yprev
  reference <- orm_markov(formula, data = data, id_var = "id")
  original <- markov_ml_estep
  passes <- 0L
  local_mocked_bindings(markov_ml_estep = function(...) {
    passes <<- passes + 1L
    original(...)
  })
  fit <- ml_fit(data, formula)

  # Nothing needs imputing: the orm fit on the observed transitions is the
  # maximum likelihood fit, and its exact information is the observed
  # information. Forward-backward passes and numerical derivatives would only
  # add time, one pass per coefficient and direction.
  expect_identical(passes, 0L)
  expect_equal(coef(fit), coef(reference), tolerance = 1e-10)
  expect_equal(
    fit$ml$var_model,
    reference$orig.var,
    tolerance = 1e-9,
    ignore_attr = TRUE
  )
  expect_identical(fit$ml$missing_info_fraction, 0)
})

test_that("without an absorbing state the fit reproduces orm_markov()", {
  data <- make_test_data(n_patients = 120, follow_up_time = 10)
  # Patients who never die: states 1 to 5, none of them absorbing.
  data <- droplevels(data[!data$id %in% data$id[data$y == 6], , drop = FALSE])
  formula <- ordered(y) ~ rms::rcs(time, 3) + tx + yprev
  reference <- orm_markov(formula, data = data, id_var = "id")
  fit <- orm_markov_ml(formula, data = data, id_var = "id")

  expect_length(fit$ml$absorbing_state, 0L)
  expect_equal(coef(fit), coef(reference), tolerance = 1e-10)
  expect_equal(fit$var, reference$var, tolerance = 1e-9, ignore_attr = TRUE)
  expect_equal(
    avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:5)$estimate,
    avg_sops(reference, variables = list(tx = c(0, 1)), times = 1:5)$estimate,
    tolerance = 1e-10
  )
})

# --- Ocens() responses ------------------------------------------------------

test_that("ranges only on last rows reproduce rms::orm() with Ocens()", {
  raw <- make_markov_trajectories(
    n_patients = 150,
    follow_up_time = 6,
    n_states = 4,
    seed = 3
  )
  data <- prepare_markov_data(raw, absorbing_state = 4, ordered_response = FALSE)
  data$y <- as.numeric(as.character(data$y))
  last <- !duplicated(data$id, fromLast = TRUE)
  set.seed(2)
  # Right-censoring only. With left- or interval-censoring present, rms
  # recodes the response through icenReg and moves finite endpoints one state
  # inward (Ocens(1, Inf) becomes "above 2"), so it no longer fits the same
  # model; left-censored masks are checked against hand-coded tables below.
  above <- last & data$y %in% c(2, 3) & stats::runif(nrow(data)) < 0.6
  data$lo <- data$y
  data$hi <- data$y
  data$lo[above] <- 1
  data$hi[above] <- Inf
  formula <- rms::Ocens(lo, hi) ~ rms::rcs(time, c(1, 3, 5)) + tx + yprev

  reference <- rms::orm(formula, data = data)
  fit <- orm_markov_ml(formula, data = data, id_var = "id", absorbing_state = 4)
  expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance = 1e-6)
  expect_equal(fit$ml$loglik, -reference$deviance[[2L]] / 2, tolerance = 1e-8)
})

test_that("Ocens() endpoints map to the documented masks", {
  data <- data.frame(
    id = 1,
    time = 1:6,
    lo = c(1, 2, -Inf, 2, NA, 3),
    hi = c(1, 3, 3, Inf, NA, 3),
    yprev = 1
  )
  parsed <- markov_ml_parse_response(
    rms::Ocens(lo, hi) ~ yprev,
    data
  )
  # The scale comes from exact states 1 and 3 here, so add the others.
  parsed$lower <- c(parsed$lower, 2, 4)
  parsed$upper <- c(parsed$upper, 2, 4)
  scale <- markov_ml_state_scale(parsed, rep("1", 8))
  masks <- markov_ml_masks(scale$lo, scale$hi, scale$K, 4L, "non_absorbing")
  expected <- rbind(
    c(1, 0, 0, 0),
    c(0, 1, 1, 0),
    c(1, 1, 0, 0),
    c(0, 0, 1, 1),
    c(1, 1, 1, 0),
    c(0, 0, 1, 0),
    c(0, 1, 0, 0),
    c(0, 0, 0, 1)
  )
  expect_equal(masks, expected)

  # Factor endpoints map through their levels.
  states <- factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  fdata <- data.frame(
    id = 1,
    lo = states[c(1, 1, 2, 3)],
    hi = states[c(1, 2, 2, 3)],
    yprev = 1
  )
  fparsed <- markov_ml_parse_response(rms::Ocens(lo, hi) ~ yprev, fdata)
  fscale <- markov_ml_state_scale(fparsed, rep("1", 4))
  expect_equal(fscale$labels, c("a", "b", "c"))
  expect_equal(fscale$lo, c(1L, 1L, 2L, 3L))
  expect_equal(fscale$hi, c(1L, 2L, 2L, 3L))

  # An all-exact Ocens() response is returned as a plain vector.
  plain <- markov_ml_parse_response(
    rms::Ocens(lo, lo) ~ yprev,
    data.frame(lo = c(1, 2, 3), yprev = 1)
  )
  expect_false(plain$interval)
  expect_equal(plain$lower, c(1, 2, 3))
})

test_that("full ranges give the same fit as missing states", {
  data <- ml_oracle_data()
  plain <- orm_markov_ml(
    y ~ time + tx + yprev,
    data = data[!data$id %in% c(107, 108), ],
    id_var = "id",
    absorbing_state = 4
  )
  any_range <- data[!data$id %in% c(107, 108), ]
  missing <- is.na(any_range$y)
  any_range$y_lo[missing] <- 1
  any_range$y_hi[missing] <- 4
  ranged <- orm_markov_ml(
    rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = any_range,
    id_var = "id",
    absorbing_state = 4
  )
  expect_equal(coef(ranged), coef(plain), tolerance = 1e-10)

  alive <- data[!data$id %in% c(105, 107, 108), ]
  plain_alive <- orm_markov_ml(
    y ~ time + tx + yprev,
    data = alive,
    id_var = "id",
    absorbing_state = 4,
    missing_states = "non_absorbing"
  )
  missing <- is.na(alive$y)
  alive$y_lo[missing] <- 1
  alive$y_hi[missing] <- 3
  ranged_alive <- orm_markov_ml(
    rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = alive,
    id_var = "id",
    absorbing_state = 4
  )
  expect_equal(coef(ranged_alive), coef(plain_alive), tolerance = 1e-10)
})

# --- Scores and information --------------------------------------------------

ml_numeric_gradient <- function(problem, theta, h = 1e-6) {
  vapply(
    seq_along(theta),
    function(k) {
      e <- replace(numeric(length(theta)), k, h)
      (markov_ml_loglik(problem, theta + e) -
        markov_ml_loglik(problem, theta - e)) /
        (2 * h)
    },
    numeric(1)
  )
}

test_that("patient scores sum to the gradient of the log-likelihood", {
  data <- ml_oracle_data()
  problem <- ml_oracle_problem(data)
  theta <- ml_random_theta(stats::coef(problem$start_fit), problem$n_int, 9)
  score <- markov_ml_total_score(problem, theta)
  expect_equal(
    unname(score),
    ml_numeric_gradient(problem, theta),
    tolerance = 1e-6
  )

  fit <- orm_markov_ml(
    rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorbing_state = 4
  )
  expect_lt(fit$ml$max_abs_score, 1e-6)
  expect_equal(nrow(fit$ml$patient_scores), length(unique(data$id)))
  expect_equal(
    unname(fit$ml$patient_scores["106", ]),
    rep(0, length(coef(fit)))
  )
})

test_that("po_row_scores() matches compute_scores_orm()", {
  data <- make_test_data(n_patients = 40, follow_up_time = 8)
  fit <- orm_markov(ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    data = data, id_var = "id")
  expected <- compute_scores_orm(fit)
  actual <- po_row_scores(
    fit$x,
    as.integer(fit$y),
    coef(fit),
    fit$non.slopes
  )
  expect_equal(unname(actual), unname(expected), tolerance = 1e-10)
})

test_that("exact proportional-odds information and grouped scores are exact", {
  data <- make_test_data(n_patients = 40, follow_up_time = 8)
  fit <- rms::orm(ordered(y) ~ rms::rcs(time, 3) + tx + yprev, data = data,
    x = TRUE, y = TRUE)
  y <- as.integer(fit$y)
  k <- fit$non.slopes
  J <- po_information(fit$x, y, coef(fit), k)
  expect_equal(solve(J), vcov(fit, intercepts = "all"), tolerance = 1e-10,
    ignore_attr = TRUE)

  # Weighted sums within groups equal sums of the row scores.
  w <- seq_along(y) %% 3 + 1
  group <- rep(1:5, length.out = length(y))
  expect_equal(
    po_group_scores(fit$x, y, coef(fit), k, group, 5L, weights = w),
    rowsum(po_row_scores(fit$x, y, coef(fit), k) * w, group),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )

  # Away from the estimate, the weighted information is the negative
  # derivative of the weighted score.
  theta <- coef(fit) + 0.1
  score <- function(t) drop(po_group_scores(fit$x, y, t, k, weights = w))
  H <- vapply(seq_along(theta), function(j) {
    e <- replace(numeric(length(theta)), j, 1e-5)
    -(score(theta + e) - score(theta - e)) / 2e-5
  }, numeric(length(theta)))
  expect_equal(po_information(fit$x, y, theta, k, weights = w), H,
    tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("the observed information matches a numerical Hessian", {
  data <- ml_oracle_data()
  fit <- orm_markov_ml(
    rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorbing_state = 4
  )
  problem <- ml_oracle_problem(data)
  theta <- coef(fit)
  h <- 1e-3
  p <- length(theta)
  H <- matrix(0, p, p)
  for (j in seq_len(p)) {
    for (k in seq_len(p)) {
      ej <- replace(numeric(p), j, h)
      ek <- replace(numeric(p), k, h)
      H[j, k] <- (markov_ml_loglik(problem, theta + ej + ek) -
        markov_ml_loglik(problem, theta + ej - ek) -
        markov_ml_loglik(problem, theta - ej + ek) +
        markov_ml_loglik(problem, theta - ej - ek)) /
        (4 * h^2)
    }
  }
  expect_equal(unname(fit$ml$information), -H, tolerance = 1e-4)
})

test_that("the observed information differentiates only rows that need the recursion", {
  data <- make_test_data(n_patients = 60, follow_up_time = 8)
  data <- data[order(data$id, data$time), ]
  # Two untreated patients miss their day-4 state, so only their day-4 and
  # day-5 rows need the forward-backward recursion.
  day4 <- data$time == 4 & data$tx == 0 & data$y != 6
  ids <- utils::head(intersect(data$id[day4], data$id[data$time == 5]), 2L)
  data$y[data$id %in% ids & data$time == 4] <- NA
  data$yprev[data$id %in% ids & data$time == 5] <- NA
  problem <- markov_ml_setup(y ~ rms::rcs(time, 3) + tx + yprev, data, "id",
    "time", "yprev", 6, "any", list())
  expect_identical(problem$chain$n, 4L)
  theta <- stats::coef(problem$start_fit)

  # Central differences of the score summed over every row.
  p <- length(theta)
  reference <- matrix(0, p, p)
  for (k in seq_len(p)) {
    h <- 1e-4 * max(1, abs(theta[[k]]))
    e <- replace(numeric(p), k, h)
    reference[, k] <- -(markov_ml_total_score(problem, theta + e) -
      markov_ml_total_score(problem, theta - e)) / (2 * h)
  }
  reference <- (reference + t(reference)) / 2

  full_design <- nrow(problem$X)
  n_complete <- length(problem$complete_part$agg_count)
  original_probs <- markov_ml_transition_probs
  original_scores <- po_row_scores
  probs_all_designs <- 0L
  scores_all_complete <- 0L
  local_mocked_bindings(
    markov_ml_transition_probs = function(problem, theta, strict = TRUE) {
      if (nrow(problem$X) == full_design) {
        probs_all_designs <<- probs_all_designs + 1L
      }
      original_probs(problem, theta, strict = strict)
    },
    po_row_scores = function(x, ...) {
      if (NROW(x) == n_complete) {
        scores_all_complete <<- scores_all_complete + 1L
      }
      original_scores(x, ...)
    }
  )
  J <- markov_ml_information(problem, theta)

  # Fully observed transitions add their exact proportional-odds information,
  # so the cost of the numerical derivatives grows with the rows that need
  # the recursion, not with all rows.
  expect_equal(unname(J), unname(reference), tolerance = 1e-6)
  expect_lte(probs_all_designs, 1L)
  expect_lte(scores_all_complete, 1L)
})

# --- Fitting -----------------------------------------------------------------

test_that("the estimate maximises the likelihood and EM is monotone", {
  data <- ml_sim_missing(n_patients = 60, follow_up_time = 8, missing = 0.2)
  formula <- y ~ rms::rcs(time, 3) + tx + yprev
  fit <- ml_fit(data, formula)
  expect_true(fit$ml$converged)
  expect_true(all(diff(fit$ml$trace) > -1e-10))
  expect_lt(fit$ml$newton_decrement, 1e-10)

  problem <- markov_ml_setup(
    formula, data, "id", "time", "yprev", 6, "any", list()
  )
  start <- fit$ml$start_coefficients
  optimum <- stats::optim(
    start,
    fn = function(theta) {
      ll <- markov_ml_loglik(problem, theta, strict = FALSE)
      if (is.na(ll)) Inf else -ll
    },
    gr = function(theta) -markov_ml_total_score(problem, theta),
    method = "BFGS",
    control = list(reltol = 1e-15, maxit = 5000)
  )
  expect_equal(unname(coef(fit)), unname(optimum$par), tolerance = 1e-5)
  expect_equal(fit$ml$loglik, -optimum$value, tolerance = 1e-8)
  expect_gt(fit$ml$n_missing, 0)
  expect_gt(fit$ml$missing_info_fraction, 0)
  expect_lt(fit$ml$missing_info_fraction, 1)
})

test_that("the fit does not depend on row order, ID labels, or empty rows", {
  data <- ml_sim_missing(n_patients = 50, follow_up_time = 8, missing = 0.2)
  fit <- ml_fit(data)

  set.seed(4)
  shuffled <- data[sample(nrow(data)), ]
  expect_equal(coef(ml_fit(shuffled)), coef(fit), tolerance = 1e-10)

  relabelled <- data
  relabelled$id <- paste0("p", 1000 - relabelled$id)
  expect_equal(coef(ml_fit(relabelled)), coef(fit), tolerance = 1e-10)

  last <- data[!duplicated(data$id, fromLast = TRUE), ]
  last <- last[last$y != "6" | is.na(last$y), ]
  extra <- last
  extra$time <- extra$time + 1
  extra$y[] <- NA
  extra$yprev[] <- NA
  appended <- rbind(data, extra)
  expect_equal(coef(ml_fit(appended)), coef(fit), tolerance = 1e-10)
})

test_that("both rules for missing states agree without deaths in gaps or trailing runs", {
  data <- ml_oracle_data()
  keep <- !data$id %in% c(104, 105, 106, 107, 108)
  any <- orm_markov_ml(y ~ time + tx + yprev, data = data[keep, ],
    id_var = "id", absorbing_state = 4)
  alive <- orm_markov_ml(y ~ time + tx + yprev, data = data[keep, ],
    id_var = "id", absorbing_state = 4, missing_states = "non_absorbing")
  expect_equal(coef(alive), coef(any), tolerance = 1e-10)

  trailing_any <- orm_markov_ml(y ~ time + tx + yprev, data = data,
    id_var = "id", absorbing_state = 4)
  trailing_alive <- orm_markov_ml(y ~ time + tx + yprev, data = data,
    id_var = "id", absorbing_state = 4, missing_states = "non_absorbing")
  expect_false(isTRUE(all.equal(coef(trailing_any), coef(trailing_alive))))
})

# --- Errors ------------------------------------------------------------------

test_that("invalid input gives informative errors", {
  data <- ml_oracle_data()
  fit_ml <- function(...) {
    orm_markov_ml(
      y ~ time + tx + yprev,
      data = data,
      id_var = "id",
      absorbing_state = 4,
      ...
    )
  }
  expect_error(
    orm_markov_ml(y ~ time + tx, data = data, id_var = "id",
      absorbing_state = 4),
    "must appear on the right-hand side"
  )
  expect_error(fit_ml(p2_var = "ypprev"), "first-order")
  expect_error(fit_ml(weights = rep(1, nrow(data))), "weights")
  expect_error(fit_ml(subset = data$id > 0), "subset")
  expect_error(fit_ml(offset = rep(0, nrow(data))), "offset")
  expect_error(
    orm_markov_ml(y ~ time + tx + yprev, data = data, absorbing_state = 4),
    "`id_var` is required"
  )
  no_baseline <- data
  no_baseline$yprev[no_baseline$id == 101 & no_baseline$time == 1] <- NA
  expect_error(
    orm_markov_ml(y ~ time + tx + yprev, data = no_baseline, id_var = "id",
      absorbing_state = 4),
    "baseline state.*101"
  )

  inconsistent <- data
  inconsistent$yprev[inconsistent$id == 101 & inconsistent$time == 5] <- "3"
  expect_error(
    orm_markov_ml(y ~ time + tx + yprev, data = inconsistent, id_var = "id",
      absorbing_state = 4),
    "disagrees.*101"
  )

  impossible <- rbind(
    data,
    data.frame(id = 999, time = 1:3, y_lo = c(2, 4, 2), y_hi = c(2, 4, 2),
      yprev = factor(c(1, 2, NA), levels = 1:3), tx = 0, y = c(2, 4, 2))
  )
  expect_error(
    orm_markov_ml(y ~ time + tx + yprev, data = impossible, id_var = "id",
      absorbing_state = 4),
    "zero probability.*999"
  )

  out_of_range <- data
  out_of_range$y_lo[1] <- 5
  out_of_range$y_hi[1] <- 6
  expect_error(
    orm_markov_ml(rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
      data = out_of_range, id_var = "id", absorbing_state = 4),
    "empty or outside"
  )
  empty <- data
  empty$y_lo[1] <- 3
  empty$y_hi[1] <- 2
  expect_error(
    orm_markov_ml(rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
      data = empty, id_var = "id", absorbing_state = 4),
    "empty or outside"
  )

  fit <- fit_ml()
  expect_error(
    sops(fit, times = 1:3, absorb = "4", p2_var = "yprev"),
    "first-order"
  )
  expect_error(
    avg_sops(fit, variables = list(tx = 0:1), times = 1:3, absorb = "4",
      p2_var = "yprev"),
    "first-order"
  )
  expect_error(robcov_orm(fit, cluster = data$id), "complete-data sandwich")
  expect_error(compute_scores_orm(fit), "complete-data scores")
})

test_that("given absorbing states are still checked", {
  data <- ml_oracle_data()
  fit_ml <- function(absorbing_state) {
    orm_markov_ml(y ~ time + tx + yprev, data = data, id_var = "id",
      absorbing_state = absorbing_state)
  }
  expect_error(fit_ml(NA), "without missing values")
  expect_error(fit_ml(c(4, NA)), "without missing values")
  expect_error(fit_ml(list(4)), "without missing values")
  expect_error(fit_ml(7), "7 are not among the observed states")
  # State 4 is never left in these data; omitting it is caught.
  expect_error(fit_ml(NULL), "give it in `absorbing_state`")
  expect_error(fit_ml(integer()), "give it in `absorbing_state`")
})

# --- Downstream inference ----------------------------------------------------

test_that("prediction and inference run on an ML fit without robcov_orm()", {
  data <- ml_sim_missing(n_patients = 80, follow_up_time = 8, missing = 0.15)
  fit <- ml_fit(data)
  local_mocked_bindings(
    robcov_orm = function(...) stop("robcov_orm() must not be called"),
    compute_scores_orm = function(...) stop("compute_scores_orm() must not be called")
  )

  a <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:8, absorb = "6")
  for (target in c("conditional", "unconditional")) {
    out <- inferences(a, method = "delta", vcov = target)
    expect_true(all(is.finite(out$std.error)))
  }
  set.seed(1)
  expect_true(all(is.finite(inferences(a, method = "mvn", n_draws = 20)$std.error)))
  set.seed(1)
  expect_true(all(is.finite(
    inferences(a, method = "score_bootstrap", n_draws = 20)$std.error
  )))

  cmp <- avg_comparisons(fit, variables = list(tx = c(0, 1)), times = 1:8,
    absorb = "6")
  expect_true(all(is.finite(
    inferences(cmp, method = "delta", vcov = "unconditional")$std.error
  )))
  tm <- avg_time(fit, variables = list(tx = c(0, 1)), times = 1:8, absorb = "6")
  expect_true(all(is.finite(inferences(tm, method = "delta")$std.error)))
  s <- sops(fit, times = 1:3, absorb = "6")
  expect_true(all(is.finite(inferences(s, method = "delta")$std.error)))

  expect_error(inferences(a, method = "fwb", n_draws = 5), "not available")
})

test_that("without an absorbing state inference and refits keep none", {
  data <- ml_oracle_data_no_absorbing()
  fit <- orm_markov_ml(rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = data, id_var = "id")
  expect_null(fit$call$absorbing_state)

  a <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:4)
  for (target in c("conditional", "unconditional")) {
    out <- inferences(a, method = "delta", vcov = target)
    expect_true(all(is.finite(out$std.error)))
  }

  boot <- rbind(
    transform(data, new_id = paste0(id, "_1")),
    transform(data[data$id %in% unique(data$id)[1:3], ],
      new_id = paste0(id, "_2"))
  )
  refit <- update_bootstrap_model(fit, boot)
  expect_length(refit$ml$absorbing_state, 0L)
  expect_equal(refit$ml$n_patients, length(unique(boot$new_id)))
  set.seed(3)
  out <- suppressWarnings(inferences(a, method = "bootstrap", n_draws = 4))
  expect_true(all(is.finite(out$std.error)))
})

# --- Printing ----------------------------------------------------------------

test_that("print() reports the ML fit, not the weighted refit", {
  data <- ml_sim_missing(n_patients = 60, follow_up_time = 6, missing = 0.2)
  fit <- orm_markov_ml(y ~ time + tx + yprev, data = data, id_var = "id",
    absorbing_state = 6)
  out <- utils::capture.output(returned <- print(fit))
  expect_identical(returned, fit)
  text <- paste(out, collapse = "\n")

  # Statistics of the refit on expected transition counts describe neither the
  # rows nor the patients, so they are not shown.
  for (refit_only in c("Frequencies of Responses", "LR chi2", "Discrimination",
    "Dxy", "Distinct Y")) {
    expect_false(grepl(refit_only, text, fixed = TRUE), info = refit_only)
  }
  expect_match(text, format(fit$ml$loglik, nsmall = 4), fixed = TRUE)
  expect_match(text, paste0("Patients: ", fit$ml$n_patients), fixed = TRUE)
  expect_match(text, "Wald Z", fixed = TRUE)

  # One row per coefficient with the estimate and the robust standard error.
  se <- sqrt(diag(vcov(fit, intercepts = "all")))
  for (name in names(coef(fit))) {
    line <- out[startsWith(out, paste0(name, " "))]
    expect_length(line, 1L)
    printed <- as.numeric(
      strsplit(trimws(substring(line, nchar(name) + 1L)), " +")[[1L]][1:2]
    )
    expect_lt(max(abs(printed - c(coef(fit)[[name]], se[[name]]))), 1e-3)
  }

  no_variance <- orm_markov_ml(y ~ time + tx + yprev, data = data,
    id_var = "id", absorbing_state = 6,
    control = orm_markov_ml_control(variance = FALSE))
  out <- utils::capture.output(print(no_variance))
  expect_true(any(grepl("not computed", out, fixed = TRUE)))
  expect_true(any(startsWith(out, "tx ")))
})

test_that("refit bootstrap refits orm_markov_ml() on whole patients", {
  data <- ml_oracle_data()
  fit <- orm_markov_ml(
    rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorbing_state = 4
  )
  boot <- rbind(
    transform(data, new_id = paste0(id, "_1")),
    transform(data[data$id %in% c(101, 106, 107), ], new_id = paste0(id, "_2"))
  )
  refit <- update_bootstrap_model(fit, boot)
  expect_s3_class(refit, "orm_markov_ml")
  expect_identical(markov_model_fit_wrapper(refit), "orm_markov_ml")
  expect_equal(nrow(markov_model_data(refit)), nrow(boot))
  expect_equal(refit$ml$n_patients, length(unique(boot$new_id)))
  expect_gt(refit$ml$n_missing, sum(is.na(data$y_lo)))
  expect_null(refit$ml$var_robust)

  a <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:4, absorb = "4")
  set.seed(3)
  out <- suppressWarnings(inferences(a, method = "bootstrap", n_draws = 4))
  expect_true(all(is.finite(out$std.error)))
})

test_that("influence components come from the stored ML scores", {
  data <- ml_oracle_data()
  fit <- orm_markov_ml(
    rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorbing_state = 4
  )
  components <- get_delta_score_components(fit)
  expect_equal(components$scores, fit$ml$patient_scores)
  expect_equal(components$Ainv, nrow(fit$ml$patient_scores) * fit$ml$bread)

  profile_ids <- as.character(
    markov_validate_starting_profiles(fit)$id
  )
  expect_setequal(rownames(fit$ml$patient_scores), profile_ids)
  expect_identical(
    rownames(fit$ml$patient_scores),
    attr(fit, "markov_starting_profile_metadata")$fitted_ids
  )

  U <- fit$ml$patient_scores[fit$ml$contributing, , drop = FALSE]
  n <- nrow(U)
  hand <- fit$ml$bread %*% crossprod(U) %*% fit$ml$bread * n / (n - 1)
  expect_equal(fit$var, hand, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(get_vcov_robust(fit), fit$var)
  expect_equal(
    get_vcov_robust(fit, type = "HC1"),
    hand * (fit$ml$n_informative_rows - 1) /
      (fit$ml$n_informative_rows - length(coef(fit))),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )

  sb <- score_bootstrap_components(fit)
  expect_equal(sb$scores, fit$ml$patient_scores)
  expect_equal(sb$bread, fit$ml$bread)
})

test_that("near-separated data fall back to Levenberg-Marquardt", {
  full <- prepare_markov_data(
    sim_actt2_markov(n_patients = 400, follow_up_time = 28,
      treatment_effect = log(0.8), seed = 2026),
    absorbing_state = 8
  )
  full <- full[order(full$id, full$time), ]
  set.seed(1)
  drop <- full$y != "8" & stats::runif(nrow(full)) < 0.25
  full$y[drop] <- NA
  same <- c(FALSE, full$id[-1] == full$id[-nrow(full)])
  full$yprev[c(FALSE, drop[-nrow(full)]) & same] <- NA
  # A bootstrap resample whose starting fit is nearly separated.
  ids <- unique(full$id)
  set.seed(35)
  picked <- sample(ids, replace = TRUE)
  boot <- do.call(rbind, lapply(seq_along(picked), function(k) {
    transform(full[full$id == picked[[k]], ], new_id = paste0(picked[[k]], "_", k))
  }))
  formula <- y ~ rms::rcs(time, 4) + tx + yprev
  problem <- markov_ml_setup(formula, boot, "new_id", "time", "yprev", 8,
    "non_absorbing", list())
  theta <- stats::coef(problem$start_fit)
  state <- markov_ml_estep(problem, theta)
  triples <- markov_ml_aggregate_weights(problem, state)

  nr <- suppressMessages(rms::orm.fit(
    x = problem$X[triples$design, , drop = FALSE], y = triples$state,
    weights = triples$weight, initial = unname(theta), compstats = FALSE
  ))
  expect_true(isTRUE(nr$fail))
  updated <- markov_ml_mstep(problem, triples, theta)
  expect_gt(markov_ml_loglik(problem, updated), state$loglik)

  fit <- expect_no_warning(orm_markov_ml(formula, data = boot,
    id_var = "new_id", absorbing_state = 8, missing_states = "non_absorbing"))
  expect_true(fit$ml$converged)
  expect_lt(fit$ml$newton_decrement, 1e-6)
})

test_that("the M-step uses a user-supplied opt_method as given", {
  data <- ml_sim_missing(n_patients = 50, follow_up_time = 8, missing = 0.2)
  formula <- y ~ rms::rcs(time, 3) + tx + yprev
  original <- rms::orm.fit
  used <- character()
  local_mocked_bindings(
    orm.fit = function(..., opt_method = "NR") {
      used <<- c(used, opt_method)
      original(..., opt_method = opt_method)
    },
    .package = "rms"
  )

  lm_problem <- markov_ml_setup(formula, data, "id", "time", "yprev", 6,
    "any", list(opt_method = "LM"))
  default_problem <- markov_ml_setup(formula, data, "id", "time", "yprev", 6,
    "any", list())
  theta <- stats::coef(default_problem$start_fit)
  state <- markov_ml_estep(default_problem, theta)
  triples <- markov_ml_aggregate_weights(default_problem, state)

  used <- character()
  markov_ml_mstep(lm_problem, triples, theta)
  expect_identical(used, "LM")

  used <- character()
  markov_ml_mstep(default_problem, triples, theta)
  expect_identical(used, "NR")

  fit_lm <- ml_fit(data, formula, opt_method = "LM")
  expect_equal(coef(fit_lm), coef(ml_fit(data, formula)), tolerance = 1e-6)
})
