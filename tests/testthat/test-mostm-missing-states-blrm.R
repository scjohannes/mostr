# mostm(engine = "blrm") with missing or range-coded states: the posterior of
# the rmsb::blrm() transition model with the observed-data likelihood
# (R/markov-ml-bayes.R and inst/stan/markov_ml.stan).
#
# Tests that compile the Stan program or sample from it run only when the
# environment variable MOSTR_TEST_STAN is "true", because compiling takes
# minutes. The others check the data passed to Stan and run by default.

skip_without_stan <- function() {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("rmsb")
  skip_if_not(
    identical(Sys.getenv("MOSTR_TEST_STAN"), "true"),
    "Set MOSTR_TEST_STAN=true to run tests that compile the Stan program."
  )
}

bayes_cppo <- function(y) as.numeric(y)

# Problem, fit object, and Stan data of a blrm fit with missing states, set up
# as mostm() does. A `cluster(id)` term in `formula` adds a random intercept
# per patient.
bayes_setup <- function(
  data,
  formula,
  missing_states = 1:3,
  ppo = NULL,
  conc = 0.7,
  normcppo = FALSE,
  pcontrast = NULL,
  npcontrast = NULL,
  psigma = 1,
  rsdmean = if (psigma == 1) 0 else 1,
  rsdsd = 1
) {
  setup <- list(
    formula = formula,
    data = data,
    dots = NULL,
    subset_rows = NULL,
    prior_weight_column = NULL,
    eval_env = environment(),
    data_name = "data"
  )
  problem <- mostm_ml_problem(
    setup,
    engine = "blrm",
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    missing_states = missing_states,
    engine_args = list(ppo = ppo)
  )
  if (is.null(problem)) {
    stop("No row of `data` needs the recursion.")
  }
  problem <- markov_bayes_problem(problem, setup)
  fit <- markov_bayes_fit_object(
    problem,
    formula,
    ppo,
    if (!is.null(ppo)) bayes_cppo,
    normcppo
  )
  input <- markov_bayes_stan_data(
    problem,
    fit,
    conc,
    pcontrast,
    npcontrast,
    random_effects = problem$random_effects,
    psigma = psigma,
    rsdmean = rsdmean,
    rsdsd = rsdsd
  )
  list(problem = problem, fit = fit, input = input)
}

# Stan-scale parameters as a one-draw array with Stan's parameter names.
bayes_stan_array <- function(alpha, beta, tau = numeric()) {
  names <- c(
    paste0("alpha[", seq_along(alpha), "]"),
    paste0("beta[", seq_along(beta), "]"),
    if (length(tau)) paste0("tau[", seq_along(tau), "]")
  )
  array(
    c(alpha, beta, tau),
    dim = c(1L, 1L, length(names)),
    dimnames = list(NULL, NULL, names)
  )
}

# Observed-data log-likelihood written as a plain loop over the runs of rows
# that need the recursion, independent of the index lists passed to Stan,
# with transition probabilities built from original-scale coefficients
# (including partial proportional-odds terms) and each patient's random
# intercept `gamma`. Every run is followed to its last row.
bayes_r_loglik <- function(setup, coefficients, gamma = NULL) {
  problem <- setup$problem
  K <- problem$K
  n_int <- K - 1L
  p <- ncol(problem$X)
  eta <- drop(problem$X %*% coefficients[n_int + seq_len(p)])
  tau <- coefficients[-seq_len(n_int + p)]
  zeta <- if (length(tau)) {
    Z <- blrm_design_matrix(setup$fit, problem$design_data, second = TRUE)
    drop(Z %*% tau)
  } else {
    0 * eta
  }
  score <- setup$input$pposcore
  gamma <- gamma %||% rep(0, problem$N)
  probs <- function(d, shift) {
    cum <- vapply(
      2:K,
      function(j) {
        stats::plogis(coefficients[[j - 1L]] + eta[d] + shift +
          score[[j]] * zeta[d])
      },
      numeric(length(d))
    )
    cum <- matrix(cum, ncol = n_int)
    cbind(1, cum) - cbind(cum, 0)
  }

  comp <- problem$complete_part
  P_obs <- probs(comp$design, gamma[comp$patient])
  loglik <- sum(log(P_obs[cbind(seq_along(comp$design), comp$state)]))

  chain <- problem$chain
  for (g in unique(chain$segment)) {
    rows <- which(chain$segment == g)
    shift <- gamma[chain$patient[rows[[1L]]]]
    f <- numeric(K)
    f[chain$baseline[[g]]] <- 1
    for (r in rows) {
      next_f <- numeric(K)
      for (ai in seq_along(problem$origins)) {
        next_f <- next_f +
          f[problem$origins[[ai]]] * probs(chain$DI[r, ai], shift)[1L, ]
      }
      next_f <- next_f * chain$C[r, ]
      loglik <- loglik + log(sum(next_f))
      f <- next_f / sum(next_f)
    }
  }
  loglik
}

# Log density of the Stan program at unconstrained parameters `u` without the
# Jacobian of the transformations, its original-scale coefficients, and its
# constrained parameters.
bayes_stan_density <- function(stanfit, setup, u) {
  pars <- rstan::constrain_pars(stanfit, u)
  list(
    log_density = rstan::log_prob(stanfit, u, adjust_transform = FALSE),
    coefficients = markov_bayes_draws_from_array(
      bayes_stan_array(pars$alpha, pars$beta, pars$tau),
      setup$input
    )[1L, ],
    pars = pars
  )
}

# Log density of the Dirichlet prior of the intercepts.
bayes_dirichlet_prior <- function(pi, conc) {
  K <- length(pi)
  lgamma(K * conc) - K * lgamma(conc) + (conc - 1) * sum(log(pi))
}

# ml_oracle_data() without the patients who reach state 4: states 1 to 3, each
# of which occurs as the previous state of a fully observed transition, so a
# row before further rows of the same patient may allow every state.
bayes_data_all_left <- function() {
  data <- ml_oracle_data()
  reaches_4 <- (!is.na(data$y_lo) & data$y_lo >= 4) |
    (!is.na(data$y_hi) & is.finite(data$y_hi) & data$y_hi >= 4)
  data[!data$id %in% data$id[reaches_4], , drop = FALSE]
}

# `data` with one more visit for a patient alive at the end of follow-up,
# whose state is missing. With every state in `missing_states`, that visit
# allows every state, so it adds nothing to the likelihood: the fit uses
# missing states, but the Stan data hold only fully observed transitions.
bayes_add_missing_visit <- function(data) {
  alive <- data$time == max(data$time) & as.character(data$y) != "6"
  extra <- data[which(alive)[[1L]], , drop = FALSE]
  extra$yprev[] <- as.character(extra$y)
  extra$time <- extra$time + 1L
  for (column in intersect(c("y", "y_factor", "y_character"), names(extra))) {
    extra[[column]][] <- NA
  }
  rbind(data, extra)
}

# --- Data set-up without Stan -------------------------------------------------

test_that("the predictor transformation matches rmsb's selectedQr()", {
  skip_if_not_installed("rmsb")
  data <- ml_sim_missing(n_patients = 40, follow_up_time = 6, seed = 7)
  setup <- bayes_setup(
    data,
    y ~ rms::rcs(time, 3) + tx + yprev,
    missing_states = 1:5,
    ppo = ~time
  )
  comp <- setup$problem$complete_part
  rows <- rep(comp$agg_design, comp$agg_count)
  qr_fun <- utils::getFromNamespace("selectedQr", "rmsb")

  reference <- qr_fun(setup$problem$X[rows, , drop = FALSE], center = TRUE)
  expect_equal(setup$input$xbar, reference$xbar, ignore_attr = TRUE)
  expect_equal(unname(setup$input$Rinv), unname(reference$Rinv))

  Z <- blrm_design_matrix(setup$fit, setup$problem$design_data, second = TRUE)
  reference_z <- qr_fun(Z[rows, , drop = FALSE], center = TRUE)
  expect_equal(setup$input$zbar, reference_z$xbar, ignore_attr = TRUE)
  expect_equal(unname(setup$input$Rinv_z), unname(reference_z$Rinv))
})

test_that("back-transformed draws give the same linear predictors", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data()
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev
  for (ppo in list(NULL, ~time)) {
    setup <- bayes_setup(data, formula, ppo = ppo)
    input <- setup$input
    stan <- input$data
    set.seed(3)
    alpha <- sort(stats::rnorm(stan$k - 1L), decreasing = TRUE)
    beta <- stats::rnorm(stan$p)
    tau <- stats::rnorm(stan$q)
    coefficients <- markov_bayes_draws_from_array(
      bayes_stan_array(alpha, beta, tau),
      input
    )[1L, ]
    expect_named(coefficients, c(input$coef_names, input$tau_names))

    # Linear predictor of P(Y >= j) for the design rows Stan sees.
    used <- input$design_rows
    expect_identical(nrow(stan$X), length(used))
    n_int <- stan$k - 1L
    original_x <- setup$problem$X[used, , drop = FALSE] %*%
      coefficients[n_int + seq_len(stan$p)]
    stan_x <- stan$X %*% beta
    if (stan$q > 0) {
      Z <- blrm_design_matrix(setup$fit, setup$problem$design_data,
        second = TRUE)[used, , drop = FALSE]
      original_z <- Z %*% coefficients[-seq_len(n_int + stan$p)]
      stan_z <- stan$Z %*% tau
    } else {
      original_z <- stan_z <- 0 * stan_x
    }
    for (j in 2:stan$k) {
      expect_equal(
        drop(coefficients[[j - 1L]] + original_x +
          input$pposcore[[j]] * original_z),
        drop(alpha[[j - 1L]] + stan_x + input$pposcore[[j]] * stan_z),
        tolerance = 1e-10,
        ignore_attr = TRUE
      )
    }
  }
})

test_that("random intercepts aggregate observed transitions within patients", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data()
  setup <- bayes_setup(
    data,
    Ocens(y_lo, y_hi) ~ time + tx + yprev + cluster(id)
  )
  stan <- setup$input$data
  problem <- setup$problem
  comp <- problem$complete_part
  expect_true(problem$random_effects)
  expect_identical(stan$Nc, problem$N)
  expect_equal(sum(stan$obs_count), length(comp$rows))
  expect_length(stan$obs_cluster, stan$n_obs)
  # Each (design row, state, patient) appears once with its count.
  key <- paste(stan$obs_design, stan$obs_state, stan$obs_cluster)
  expect_false(anyDuplicated(key) > 0)
  expect_equal(
    tapply(stan$obs_count, stan$obs_cluster, sum),
    tapply(comp$patient, comp$patient, length),
    ignore_attr = TRUE
  )
  expect_identical(as.integer(stan$chain_cluster), problem$chain$patient)
  expect_length(stan$rsdmean, 1L)
  expect_length(stan$rsdsd, 1L)

  # Without random intercepts the transitions are aggregated over patients.
  plain <- bayes_setup(data, Ocens(y_lo, y_hi) ~ time + tx + yprev)$input$data
  expect_identical(plain$Nc, 0L)
  expect_length(plain$obs_cluster, 0L)
  expect_length(plain$chain_cluster, 0L)
  expect_length(plain$rsdmean, 0L)
  expect_equal(sum(plain$obs_count), length(comp$rows))
  expect_lt(plain$n_obs, stan$n_obs)
})

test_that("the recursion visits the previous states that can have probability", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data(long_runs = TRUE)
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev
  # A missing state may be any state the model leaves (1 to 3), or only 1 or
  # 3, so that rows with a missing state exclude a state between allowed ones.
  for (missing_states in list(1:3, c(1, 3))) {
    setup <- bayes_setup(data, formula, missing_states = missing_states)
    stan <- setup$input$data
    problem <- setup$problem
    chain <- problem$chain
    K <- problem$K
    allowed <- chain$C > 0
    expect_null(stan$n_absorbing)
    expect_identical(as.integer(stan$row_lo),
      unname(apply(allowed, 1L, function(a) min(which(a)))))
    expect_identical(as.integer(stan$row_hi),
      unname(apply(allowed, 1L, function(a) max(which(a)))))
    expect_identical(
      as.logical(stan$row_gappy),
      unname(rowSums(allowed) < as.vector(stan$row_hi - stan$row_lo) + 1L)
    )
    expect_identical(
      as.logical(stan$row_full),
      unname(rowSums(allowed) == K)
    )
    expect_identical(
      any(stan$row_gappy == 1L),
      identical(missing_states, c(1, 3))
    )

    # Each run stops at its last row that excludes a state; patient 108's run
    # has no such row (its only row allows every state) and is left out.
    informative <- which(rowSums(allowed) < K)
    last <- tapply(informative, chain$segment[informative], max)
    expect_identical(as.integer(stan$seg_end), as.integer(last))
    expect_identical(
      as.integer(stan$seg_start),
      match(as.integer(names(last)), chain$segment)
    )
    expect_identical(
      as.integer(stan$seg_baseline),
      chain$baseline[as.integer(names(last))]
    )
    expect_false("108" %in% chain$ids[stan$seg_start])
    in_run <- unlist(Map(seq, stan$seg_start, stan$seg_end))
    expect_true(all(
      stan$cell_end[-in_run] < stan$cell_start[-in_run]
    ))

    # Forward through each run: the cells of a row are the states with
    # probability on the row before (the run's baseline on its first row),
    # each with the design row of the transition from it.
    used <- setup$input$design_rows
    for (g in seq_len(stan$n_seg)) {
      reach <- seq_len(K) == stan$seg_baseline[[g]]
      for (t in stan$seg_start[[g]]:stan$seg_end[[g]]) {
        cells <- seq_len(stan$cell_end[[t]] - stan$cell_start[[t]] + 1L) +
          stan$cell_start[[t]] - 1L
        origins <- which(reach)
        expect_true(all(origins %in% problem$origins))
        expect_identical(as.integer(stan$cell_origin[cells]), origins)
        design <- used[stan$cd_rows[stan$cell_cd[cells]]]
        expect_identical(
          design,
          chain$DI[cbind(rep(t, length(origins)),
            match(origins, problem$origins))]
        )
        expect_true(all(stan$cd_lo[stan$cell_cd[cells]] <= stan$row_lo[[t]]))
        expect_true(all(stan$cd_hi[stan$cell_cd[cells]] >= stan$row_hi[[t]]))
        reach <- allowed[t, ]
      }
    }
  }
})

test_that("the R reference recursion equals enumeration over completions", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data(long_runs = TRUE)
  for (missing_states in list(1:3, c(1, 3))) {
    setup <- bayes_setup(
      data,
      Ocens(y_lo, y_hi) ~ yprev + time + tx,
      missing_states = missing_states
    )
    problem <- setup$problem
    for (seed in 1:2) {
      theta <- ml_random_theta(problem, problem$start_coefficients, seed)
      expect_equal(
        bayes_r_loglik(setup, theta),
        ml_oracle_loglik(data, problem, theta, missing_states = missing_states),
        tolerance = 1e-10
      )
    }
  }
})

test_that("data that no sequence of states fits stop before sampling", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data()
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev

  # Patient 106 has four missing states in a row. If the second could only be
  # state 4, which the model never leaves, the third could not be reached.
  setup <- bayes_setup(data, formula)
  problem <- setup$problem
  rows <- which(problem$chain$ids == "106")
  problem$chain$C[rows[[2L]], ] <- c(0, 0, 0, 1)
  expect_error(
    markov_bayes_stan_data(problem, setup$fit, 0.7, NULL, NULL),
    "zero probability.*patient\\(s\\) 106\\."
  )

  # mostm() checks the data before compiling or sampling: a missing state
  # followed by further rows may not be a state that is never left, and such
  # a state may not be followed by a missing state.
  expect_error(
    mostm(formula, data = data, engine = "blrm", id_var = "id",
      missing_states = 1:4),
    "never occur as the previous state"
  )
  after_death <- data.frame(
    id = 998,
    time = 1:3,
    y_lo = c(2, 4, NA),
    y_hi = c(2, 4, NA),
    yprev = factor(c(1, 2, NA), levels = 1:3),
    tx = 0
  )
  after_death$y <- after_death$y_lo
  expect_error(
    mostm(formula, data = rbind(data, after_death), engine = "blrm",
      id_var = "id", missing_states = 1:3),
    paste0(
      "State\\(s\\) 4 are followed by a missing or range-coded state for ",
      "patient\\(s\\) 998,"
    )
  )
})

test_that("cluster() terms are checked and removed for set-up", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data()
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev + cluster(id)
  setup <- bayes_setup(data, formula)
  expect_true(setup$problem$random_effects)
  expect_identical(
    attr(stats::terms(setup$problem$formula), "term.labels"),
    c("time", "tx", "yprev")
  )
  expect_identical(environment(setup$problem$formula), environment(formula))

  fit_bayes <- function(formula, engine = "blrm") {
    mostm(formula, data = data, engine = engine, id_var = "id",
      missing_states = 1:3)
  }
  expect_error(
    fit_bayes(Ocens(y_lo, y_hi) ~ time + tx + yprev + cluster(tx)),
    "must name the patient ID column: use `cluster(id)`",
    fixed = TRUE
  )
  expect_error(
    fit_bayes(Ocens(y_lo, y_hi) ~ time + yprev + tx * cluster(id)),
    "one `cluster()` term, not inside an interaction",
    fixed = TRUE
  )
  for (engine in c("orm", "vglm")) {
    expect_error(
      fit_bayes(formula, engine = engine),
      "`engine = \"blrm\"` supports a random intercept per patient",
      fixed = TRUE
    )
  }
})

test_that("the Stan data have no absorbing states", {
  skip_if_not_installed("rmsb")
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev
  # State 4 is never left in ml_oracle_data(); every state is left in the
  # data without the patients who reach it. A row before further rows of the
  # same patient then allows every state, and its probabilities add up to 1.
  with_4 <- bayes_setup(ml_oracle_data(), formula)$input$data
  all_left <- bayes_setup(bayes_data_all_left(), formula)
  expect_identical(all_left$problem$origins, 1:3)
  for (stan in list(with_4, all_left$input$data)) {
    expect_null(stan$n_absorbing)
  }
  stan <- all_left$input$data
  inner <- unlist(Map(
    function(start, end) if (end > start) start:(end - 1L),
    stan$seg_start,
    stan$seg_end
  ))
  expect_true(any(stan$row_full[inner] == 1L))
})

test_that("normcppo centres and scales the scores over observed transitions", {
  skip_if_not_installed("rmsb")
  data <- ml_sim_missing(n_patients = 40, follow_up_time = 6, seed = 8)
  setup <- bayes_setup(
    data,
    y ~ time + tx + yprev,
    missing_states = 1:5,
    ppo = ~time,
    normcppo = TRUE
  )
  problem <- setup$problem
  comp <- problem$complete_part
  observed <- rep(comp$agg_state, comp$agg_count)
  scores <- setup$fit$cppo(problem$labels[observed])
  expect_equal(mean(scores), 0, tolerance = 1e-12)
  expect_equal(stats::sd(scores), 1, tolerance = 1e-12)
  expect_equal(setup$input$pposcore, setup$fit$cppo(problem$labels))

  expect_error(
    markov_bayes_fit_object(problem, y ~ time + tx + yprev, ~time,
      function(y) rep(1, length(y)), TRUE),
    "scores that vary"
  )
})

test_that("invalid input gives informative errors before sampling", {
  skip_if_not_installed("rmsb")
  data <- ml_oracle_data()
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev
  fit_bayes <- function(...) {
    mostm(formula, data = data, engine = "blrm", id_var = "id",
      missing_states = 1:3, ...)
  }
  expect_error(
    suppressWarnings(
      mostm(formula, data = data, engine = "blrm", missing_states = 1:3)
    ),
    "`id_var` is required"
  )
  expect_error(
    mostm(Ocens(y_lo, y_hi) ~ time + tx, data = data, engine = "blrm",
      id_var = "id", missing_states = 1:3),
    "must appear on the right-hand side"
  )
  expect_error(fit_bayes(type = "HC1"), "posterior draws")
  expect_error(fit_bayes(weights = rep(1, nrow(data))), "row weights")
  expect_error(fit_bayes(backend = "cmdstan"), "does not accept `backend`")
  expect_error(
    fit_bayes(loo = FALSE, standata = TRUE),
    "does not accept `loo`, `standata`"
  )
  expect_error(fit_bayes(psigma = 3), "`psigma`")
  expect_error(fit_bayes(psigma = 2, rsdmean = 0), "`rsdmean` must be positive")
  expect_error(fit_bayes(rsdsd = 0), "`rsdsd` a positive number")
  expect_error(fit_bayes(ppo = ~time), "constrained partial proportional odds")
  expect_error(fit_bayes(cppo = bayes_cppo), "require `ppo`")
  expect_error(fit_bayes(normcppo = TRUE), "require `ppo`")
  expect_error(
    fit_bayes(ppo = ~time, cppo = bayes_cppo, normcppo = NA),
    "`normcppo` must be TRUE or FALSE"
  )
  expect_error(
    fit_bayes(ppo = ~age, cppo = bayes_cppo),
    "`age` in `ppo` are not columns"
  )
  expect_error(fit_bayes(iter = 1), "`iter`")
  expect_error(fit_bayes(iter = 100, warmup = 100), "`warmup`")
  expect_error(fit_bayes(iter = 100, warmup = 10.5), "`warmup`")
  # Without `warmup`, an odd `iter` warms up for floor(iter / 2) iterations,
  # as rstan::sampling() does, so the next check is reached.
  expect_error(fit_bayes(iter = 301, conc = -1), "`conc`")
  expect_error(fit_bayes(chains = 0), "`chains`")
  expect_error(fit_bayes(refresh = -1), "`refresh`")
  expect_error(fit_bayes(conc = -1), "`conc`")
  expect_error(fit_bayes(sampling.args = 1), "`sampling.args` must be a list")
  expect_error(fit_bayes(seed = 1), "does not use `seed`.*`sampling.args`")
  expect_error(fit_bayes(pcontrast = 1), "`pcontrast` must be a list")
})

test_that("with nothing missing, missing_states gives rmsb's usual fit", {
  skip_on_cran()
  skip_if_not_installed("rmsb")
  data <- make_test_data(n_patients = 30, follow_up_time = 6, seed = 12)
  states <- sort(unique(as.character(data$y)))
  formula <- ordered(y) ~ time + tx + yprev
  setup <- list(
    formula = formula,
    data = data,
    dots = NULL,
    subset_rows = NULL,
    prior_weight_column = NULL,
    eval_env = environment(),
    data_name = "data"
  )
  expect_null(mostm_ml_problem(setup, "blrm", "id", "time", "yprev",
    missing_states = states, engine_args = list()))

  fit_blrm <- function(...) {
    suppressMessages(suppressWarnings(mostm(
      formula,
      data = data,
      engine = "blrm",
      id_var = "id",
      iter = 100,
      chains = 1,
      sampling.args = list(seed = 1),
      ...
    )))
  }
  with_states <- fit_blrm(missing_states = states)
  expect_s3_class(with_states, "blrm")
  expect_false(inherits(with_states, "mostm_ml"))
  expect_null(with_states$ml)
  expect_equal(with_states$draws, fit_blrm()$draws)
})

# --- Stan program ---------------------------------------------------------------

test_that("with nothing missing the log density equals rmsb's", {
  skip_without_stan()
  # rmsb::blrm() evaluates `cluster(id)`. The identity function is exported
  # by rmsb up to 1.1-2 and by rms from 8.2-0; define it locally so the test
  # works with either.
  cluster <- function(x) x
  model <- markov_bayes_stan_model()
  data <- make_test_data(n_patients = 50, follow_up_time = 6, seed = 3)
  data$y_factor <- factor(as.character(data$y))
  data$y_character <- as.character(data$y)
  states <- sort(unique(as.character(data$y)))
  extended <- bayes_add_missing_visit(data)
  # Contrast priors hold unnamed predictors at their datadist settings.
  withr::local_options(datadist = rms::datadist(data))
  rhs <- ~ rms::rcs(time, 3) + tx + yprev
  tx_prior <- list(sd = 1, c1 = list(tx = 1), c2 = list(tx = 0),
    contrast = expression(c1 - c2))
  time_prior <- list(sd = 0.5, mu = 0.1, c1 = list(time = 5),
    c2 = list(time = 2), contrast = expression(c1 - c2))
  cases <- list(
    list(),
    list(ppo = ~time),
    list(ppo = ~time, re = TRUE),
    list(conc = 0.3, pcontrast = tx_prior),
    list(ppo = ~time, normcppo = TRUE, npcontrast = time_prior),
    list(re = TRUE, psigma = 2, rsdmean = 0.7),
    list(re = TRUE, psigma = 1, rsdmean = 0.4, rsdsd = 0.6),
    # Responses other than ordered factors; rmsb::blrm() gets the factor.
    list(response = "y_factor", ppo = ~ rms::rcs(time, 3)),
    list(response = "y_character", rmsb_response = "y_factor", ppo = ~time)
  )
  for (case in cases) {
    response <- case$response %||% "y"
    re <- isTRUE(case$re)
    terms <- c(
      attr(stats::terms(rhs), "term.labels"),
      if (re) "cluster(id)"
    )
    formula <- stats::reformulate(terms, response = response,
      env = environment())
    conc <- case$conc %||% (1 / (0.8 + 0.35 * 6))
    psigma <- case$psigma %||% 1
    rsdmean <- case$rsdmean %||% if (psigma == 1) 0 else 1
    rsdsd <- case$rsdsd %||% 1
    normcppo <- isTRUE(case$normcppo)
    setup <- bayes_setup(extended, formula, missing_states = states,
      ppo = case$ppo, conc = conc, normcppo = normcppo,
      pcontrast = case$pcontrast, npcontrast = case$npcontrast,
      psigma = psigma, rsdmean = rsdmean, rsdsd = rsdsd)
    expect_identical(setup$input$data$n_seg, 0L)
    expect_identical(setup$input$data$Nc, if (re) 50L else 0L)
    ours <- suppressMessages(
      rstan::sampling(model, data = setup$input$data, chains = 0)
    )

    # rmsb::blrm() returns its Stan data before adding contrast priors, so
    # take the data from a short run instead.
    rmsb_formula <- stats::reformulate(terms,
      response = case$rmsb_response %||% response, env = environment())
    rmsb_args <- list(rmsb_formula, data = data, conc = conc,
      pcontrast = case$pcontrast, psigma = psigma, rsdmean = rsdmean,
      rsdsd = rsdsd, iter = 10, warmup = 5, chains = 1, refresh = 0,
      loo = FALSE, sampling.args = list(seed = 1))
    if (!is.null(case$ppo)) {
      rmsb_args <- c(rmsb_args, list(ppo = case$ppo, cppo = bayes_cppo,
        normcppo = normcppo, npcontrast = case$npcontrast))
    }
    theirs <- suppressMessages(suppressWarnings(
      do.call(rmsb::blrm, rmsb_args)
    ))$rstan

    n_par <- rstan::get_num_upars(ours)
    expect_identical(n_par, rstan::get_num_upars(theirs))
    expect_identical(setup$input$data$cn, length(case$pcontrast$sd))
    expect_identical(setup$input$data$cn2, length(case$npcontrast$sd))
    set.seed(1)
    for (i in 1:3) {
      u <- stats::rnorm(n_par, sd = 0.3)
      expect_equal(
        rstan::log_prob(ours, u),
        rstan::log_prob(theirs, u),
        tolerance = 1e-10
      )
    }
  }
})

test_that("draws of an rmsb fit give the same SOPs through either fit", {
  skip_without_stan()
  # mostm() lets rmsb::blrm() find `cluster()` on the search path.
  withr::local_package("rmsb")
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5)
  # Patient IDs whose sorted order differs from their numeric order, so the
  # random-intercept columns are mapped to patients by ID.
  data$id <- paste0("p", data$id)
  formula <- ordered(y) ~ rms::rcs(time, 3) + tx + yprev + cluster(id)
  rmsb_fit <- suppressMessages(suppressWarnings(mostm(
    formula,
    data = data,
    engine = "blrm",
    id_var = "id",
    ppo = ~time,
    cppo = bayes_cppo,
    iter = 100,
    chains = 2,
    refresh = 0,
    loo = FALSE,
    sampling.args = list(seed = 1)
  )))
  expect_false(inherits(rmsb_fit, "mostm_ml"))

  # The Stan-scale draws of rmsb's own program, transformed as the fit with
  # missing states transforms its draws, reproduce rmsb's draws.
  states <- sort(unique(as.character(data$y)))
  extended <- bayes_add_missing_visit(data)
  setup <- bayes_setup(extended, formula, missing_states = states,
    ppo = ~time, conc = rmsb_fit$conc)
  stan_draws <- rstan::extract(rmsb_fit$rstan, pars = c("alpha", "beta", "tau"),
    permuted = FALSE)
  draws <- markov_bayes_draws_from_array(stan_draws, setup$input)
  expect_identical(colnames(draws), colnames(rmsb_fit$draws))
  expect_equal(draws, rmsb_fit$draws, tolerance = 1e-10)

  # A fit with missing states carrying the same draws, including the
  # patients' random intercepts as that fit labels them.
  ml_fit <- suppressMessages(suppressWarnings(mostm(
    formula,
    data = extended,
    engine = "blrm",
    id_var = "id",
    missing_states = states,
    ppo = ~time,
    cppo = bayes_cppo,
    iter = 20,
    chains = 1,
    sampling.args = list(seed = 1)
  )))
  expect_s3_class(ml_fit, "mostm_ml")
  ml_fit$draws <- draws
  ml_fit$param <- rmsb_fit$param
  ml_fit$gamma_draws <- markov_bayes_add_random_effects(
    ml_fit,
    rmsb_fit$rstan,
    setup$problem,
    "id"
  )$gamma_draws

  args <- list(
    variables = list(tx = c(0, 1)),
    times = 1:4,
    absorb = "6",
    n_draws = 30,
    seed = 1,
    return_draws = TRUE
  )
  results <- list()
  for (include_re in c(FALSE, TRUE)) {
    sops_rmsb <- do.call(avg_sops, c(list(rmsb_fit), args,
      include_re = include_re))
    sops_ml <- do.call(avg_sops, c(list(ml_fit), args, include_re = include_re))
    expect_equal(as.data.frame(sops_ml), as.data.frame(sops_rmsb),
      tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(attr(sops_ml, "draws"), attr(sops_rmsb, "draws"),
      tolerance = 1e-8, ignore_attr = TRUE)

    cmp_args <- c(args, list(state_sets = list(death = "6"),
      include_re = include_re))
    cmp_rmsb <- do.call(avg_comparisons, c(list(rmsb_fit), cmp_args))
    cmp_ml <- do.call(avg_comparisons, c(list(ml_fit), cmp_args))
    expect_equal(as.data.frame(cmp_ml), as.data.frame(cmp_rmsb),
      tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(attr(cmp_ml, "draws"), attr(cmp_rmsb, "draws"),
      tolerance = 1e-8, ignore_attr = TRUE)
    results[[length(results) + 1L]] <- sops_ml$estimate
  }
  # The random intercepts change the predictions, so the comparison with
  # include_re = TRUE checks that they reach the same patients.
  expect_false(isTRUE(all.equal(results[[1L]], results[[2L]])))
})

test_that("the likelihood equals the R recursion with missing and range states", {
  skip_without_stan()
  model <- markov_bayes_stan_model()
  data <- ml_oracle_data(long_runs = TRUE)
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev
  conc <- 0.7
  for (missing_states in list(1:3, c(1, 3))) {
    for (ppo in list(NULL, ~time)) {
      setup <- bayes_setup(data, formula, missing_states = missing_states,
        ppo = ppo, conc = conc)
      stan <- setup$input$data
      stanfit <- suppressMessages(
        rstan::sampling(model, data = stan, chains = 0)
      )
      set.seed(2)
      for (i in 1:3) {
        u <- stats::rnorm(stan$p + stan$q + stan$k - 1L, sd = 0.3)
        density <- bayes_stan_density(stanfit, setup, u)
        expect_equal(
          density$log_density - bayes_dirichlet_prior(density$pars$pi, conc),
          bayes_r_loglik(setup, density$coefficients),
          tolerance = 1e-10
        )
      }
    }
  }
})

test_that("when every state is left the density matches the R recursion", {
  skip_without_stan()
  model <- markov_bayes_stan_model()
  data <- bayes_data_all_left()
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev
  conc <- 0.7
  setup <- bayes_setup(data, formula, conc = conc)
  stan <- setup$input$data
  expect_true(any(stan$row_full == 1L))
  stanfit <- suppressMessages(rstan::sampling(model, data = stan, chains = 0))
  set.seed(3)
  for (i in 1:3) {
    u <- stats::rnorm(stan$p + stan$q + stan$k - 1L, sd = 0.3)
    density <- bayes_stan_density(stanfit, setup, u)
    expect_equal(
      density$log_density - bayes_dirichlet_prior(density$pars$pi, conc),
      bayes_r_loglik(setup, density$coefficients),
      tolerance = 1e-10
    )
  }

  fit <- suppressWarnings(mostm(
    y ~ time + tx + yprev,
    data = data,
    engine = "blrm",
    id_var = "id",
    missing_states = 1:3,
    iter = 301,
    chains = 2,
    sampling.args = list(seed = 1)
  ))
  expect_s3_class(fit, "mostm_ml")
  # An odd `iter` warms up for floor(iter / 2) iterations.
  expect_identical(fit$warmup, 150)
  expect_identical(nrow(fit$draws), 2L * 151L)
  expect_true(all(is.finite(fit$draws)))
  expect_identical(fit$ml$missing_states, c("1", "2", "3"))
})

test_that("with random intercepts the density matches the R recursion", {
  skip_without_stan()
  model <- markov_bayes_stan_model()
  data <- ml_oracle_data(long_runs = TRUE)
  formula <- Ocens(y_lo, y_hi) ~ time + tx + yprev + cluster(id)
  conc <- 0.7
  for (missing_states in list(1:3, c(1, 3))) {
    for (ppo in list(NULL, ~time)) {
      setup <- bayes_setup(data, formula, missing_states = missing_states,
        ppo = ppo, conc = conc)
      stan <- setup$input$data
      stanfit <- suppressMessages(
        rstan::sampling(model, data = stan, chains = 0)
      )
      set.seed(4)
      # Priors written with `~` drop constants, so compare up to a constant.
      gaps <- vapply(1:3, function(i) {
        u <- stats::rnorm(stan$p + stan$q + stan$k - 1L + stan$Nc + 1L,
          sd = 0.3)
        density <- bayes_stan_density(stanfit, setup, u)
        pars <- density$pars
        log_prior <- bayes_dirichlet_prior(pars$pi, conc) -
          0.5 * sum(pars$gamma_raw^2) -
          2.5 * log1p(pars$sigmag^2 / 4)
        density$log_density - log_prior -
          bayes_r_loglik(setup, density$coefficients, gamma = pars$gamma)
      }, numeric(1))
      expect_equal(gaps, rep(gaps[[1L]], 3L), tolerance = 1e-10)
    }
  }
})

test_that("a fitted model works with the SOP functions", {
  skip_without_stan()
  data <- ml_sim_missing(n_patients = 60, follow_up_time = 6, missing = 0.2)
  fit <- suppressWarnings(mostm(
    y ~ time + tx + yprev,
    data = data,
    engine = "blrm",
    id_var = "id",
    missing_states = 1:5,
    ppo = ~time,
    cppo = bayes_cppo,
    iter = 300,
    chains = 2,
    sampling.args = list(seed = 1)
  ))
  expect_s3_class(fit, c("mostm_ml", "blrm", "rmsb", "rms"), exact = TRUE)
  expect_identical(markov_model_engine(fit), "blrm")
  expect_equal(nrow(fit$draws), 300L)
  expect_identical(
    colnames(fit$draws),
    c(names(fit$ml$start_coefficients), "time x f(y)")
  )
  expect_identical(rownames(fit$diagnostics), colnames(fit$draws))
  expect_equal(stats::coef(fit), colMeans(fit$draws))
  expect_identical(fit$ml$missing_states, as.character(1:5))
  expect_gt(fit$ml$n_missing, 0L)
  expect_output(print(fit), "Bayesian Markov model with missing states")
  expect_output(print(fit), "A missing state may be: 1, 2, 3, 4, 5")

  avg <- avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:4,
    absorb = "6", n_draws = 50, seed = 1, posterior_summary = "mean")
  totals <- stats::aggregate(estimate ~ tx + time, data = avg, FUN = sum)
  expect_equal(totals$estimate, rep(1, nrow(totals)), tolerance = 1e-8)
  expect_true(all(c("conf.low", "conf.high") %in% names(avg)))
  expect_error(
    sops(fit, times = 1:3, absorb = "6", p2_var = "yprev"),
    "first-order"
  )
})

test_that("random intercepts reach the SOP functions with include_re", {
  skip_without_stan()
  data <- ml_sim_missing(n_patients = 60, follow_up_time = 6, missing = 0.2)
  fit <- suppressWarnings(mostm(
    y ~ time + tx + yprev + cluster(id),
    data = data,
    engine = "blrm",
    id_var = "id",
    missing_states = 1:5,
    iter = 300,
    chains = 2,
    sampling.args = list(seed = 2)
  ))
  ids <- unique(as.character(data$id[order(data$id)]))
  expect_identical(dim(fit$gamma_draws), c(300L, length(ids)))
  expect_identical(colnames(fit$gamma_draws), ids)
  expect_identical(colnames(fit$omega), "sigmag")
  expect_identical(fit$clusterInfo$name, "id")
  expect_true("sigmag" %in% rownames(fit$diagnostics))
  expect_output(print(fit), "Random intercepts for 60 patients")

  args <- list(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:4,
    absorb = "6",
    n_draws = 40,
    seed = 1,
    posterior_summary = "mean"
  )
  typical <- do.call(avg_sops, args)
  own <- do.call(avg_sops, c(args, include_re = TRUE))
  for (result in list(typical, own)) {
    totals <- stats::aggregate(estimate ~ tx + time, data = result, FUN = sum)
    expect_equal(totals$estimate, rep(1, nrow(totals)), tolerance = 1e-8)
  }
  expect_false(isTRUE(all.equal(typical$estimate, own$estimate)))

  comparison <- avg_comparisons(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:4,
    state_sets = list(death = "6"),
    absorb = "6",
    include_re = TRUE,
    n_draws = 40,
    seed = 1
  )
  expect_true(all(is.finite(comparison$estimate)))
  expect_true(all(comparison$conf.low <= comparison$conf.high))
})
