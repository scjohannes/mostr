# blrm_markov_ml(): Bayesian fit with the observed-data likelihood.
#
# Tests that need the compiled Stan program run only when the environment
# variable MOSTR_TEST_STAN is "true", because compiling takes minutes.

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

# Set-up shared by the density tests: problem, fit object, and Stan data.
bayes_setup <- function(
  data,
  formula,
  ppo = NULL,
  missing_states = "any",
  absorbing_state = 6,
  conc = 0.7,
  random_effects = FALSE,
  normcppo = FALSE,
  pcontrast = NULL,
  npcontrast = NULL,
  psigma = 1,
  rsdmean = if (psigma == 1) 0 else 1,
  rsdsd = 1
) {
  problem <- markov_ml_setup(
    formula,
    data,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    absorbing_state = absorbing_state,
    missing_states = missing_states,
    dots_fit = list(),
    extra_vars = if (!is.null(ppo)) all.vars(ppo) else character()
  )
  fit <- markov_bayes_fit_object(
    problem,
    formula,
    ppo,
    if (!is.null(ppo)) bayes_cppo,
    normcppo = normcppo
  )
  input <- markov_bayes_stan_data(problem, fit, conc, pcontrast, npcontrast,
    random_effects = random_effects, psigma = psigma, rsdmean = rsdmean,
    rsdsd = rsdsd)
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

# Observed-data log-likelihood written as a plain loop over patients' runs of
# rows, independent of the vectorised recursion, with transition probabilities
# built from original-scale coefficients (including partial proportional-odds
# terms) and each patient's random intercept `gamma`.
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
      next_f[problem$absorbing] <- next_f[problem$absorbing] +
        f[problem$absorbing]
      next_f <- next_f * chain$C[r, ]
      loglik <- loglik + log(sum(next_f))
      f <- next_f / sum(next_f)
    }
  }
  loglik
}

# Rows with an unknown state allow every state but 2, so their allowed states
# have a gap, as with an absorbing state inside the scale and
# `missing_states = "non_absorbing"`.
bayes_gappy_setup <- function(setup, conc, random_effects = FALSE) {
  problem <- setup$problem
  open <- rowSums(problem$chain$C) == problem$K
  problem$chain$C[open, 2L] <- 0
  setup$problem <- problem
  setup$input <- markov_bayes_stan_data(problem, setup$fit, conc, NULL, NULL,
    random_effects = random_effects)
  stopifnot(any(setup$input$data$row_gappy == 1L))
  setup
}

# --- Data set-up without Stan -------------------------------------------------

test_that("the predictor transformation matches rmsb's selectedQr()", {
  skip_if_not_installed("rmsb")
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 7)
  setup <- bayes_setup(data, y ~ rms::rcs(time, 3) + tx + yprev, ppo = ~time)
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
  data <- ml_oracle_data()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  for (ppo in list(NULL, ~time)) {
    setup <- bayes_setup(data, formula, ppo = ppo, absorbing_state = 4)
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
  data <- ml_oracle_data()
  setup <- bayes_setup(data, rms::Ocens(y_lo, y_hi) ~ time + tx + yprev,
    absorbing_state = 4, random_effects = TRUE)
  stan <- setup$input$data
  problem <- setup$problem
  comp <- problem$complete_part
  expect_identical(stan$Nc, problem$N)
  expect_equal(sum(stan$obs_count), length(comp$rows))
  expect_length(stan$obs_cluster, stan$n_obs)
  # Each (design row, state, patient) appears once with its count.
  key <- paste(stan$obs_design, stan$obs_state, stan$obs_cluster)
  expect_false(anyDuplicated(key) > 0)
  expect_identical(as.integer(stan$chain_cluster), problem$chain$patient)
  expect_length(stan$rsdmean, 1L)
  expect_length(stan$rsdsd, 1L)
})

test_that("the recursion visits the previous states that can have probability", {
  data <- ml_oracle_data()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  for (rule in c("any", "non_absorbing")) {
    setup <- bayes_setup(data, formula, missing_states = rule,
      absorbing_state = 4)
    stan <- setup$input$data
    problem <- setup$problem
    chain <- problem$chain
    K <- problem$K
    allowed <- chain$C > 0
    expect_identical(as.integer(stan$row_lo),
      unname(apply(allowed, 1L, function(a) min(which(a)))))
    expect_identical(as.integer(stan$row_hi),
      unname(apply(allowed, 1L, function(a) max(which(a)))))
    expect_identical(
      as.logical(stan$row_full),
      unname(rowSums(allowed) == K)
    )

    # Forward pass with equal transition probabilities: the cells of a row are
    # the non-absorbing states with probability on the row before.
    used <- setup$input$design_rows
    for (g in seq_len(stan$n_seg)) {
      f <- as.numeric(seq_len(K) == stan$seg_baseline[[g]])
      for (t in stan$seg_start[[g]]:stan$seg_end[[g]]) {
        cells <- seq_len(stan$cell_end[[t]] - stan$cell_start[[t]] + 1L) +
          stan$cell_start[[t]] - 1L
        origins <- which(f > 0 & !seq_len(K) %in% problem$absorbing)
        expect_identical(as.integer(stan$cell_origin[cells]), origins)
        design <- used[stan$cd_rows[stan$cell_cd[cells]]]
        expect_identical(
          design,
          chain$DI[cbind(rep(t, length(origins)),
            match(origins, problem$origins))]
        )
        expect_true(all(stan$cd_lo[stan$cell_cd[cells]] <= stan$row_lo[[t]]))
        expect_true(all(stan$cd_hi[stan$cell_cd[cells]] >= stan$row_hi[[t]]))
        f <- (sum(f[origins]) / K + f * seq_len(K) %in% problem$absorbing) *
          allowed[t, ]
      }
    }
  }
})

test_that("data that no sequence of states fits stop before sampling", {
  data <- ml_oracle_data()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  # Both patients die (state 4) at time 2. Patient 998's state at time 3 is
  # missing, and patient 999 is recorded in state 2 at time 3.
  after_death <- data.frame(
    id = rep(c(998, 999), each = 3L),
    time = rep(1:3, 2L),
    y_lo = c(2, 4, NA, 2, 4, 2),
    y_hi = c(2, 4, NA, 2, 4, 2),
    yprev = factor(c(1, 2, NA, 1, 2, NA), levels = 1:3),
    tx = 0
  )
  after_death$y <- after_death$y_lo
  impossible <- rbind(data, after_death)

  # With missing states restricted to the non-absorbing ones, neither patient
  # can be alive at time 3 after dying at time 2.
  expect_error(
    bayes_setup(impossible, formula, missing_states = "non_absorbing",
      absorbing_state = 4),
    "zero probability.*998, 999\\."
  )
  # A missing state may be death, so only patient 999 has no possible states.
  expect_error(
    bayes_setup(impossible, formula, absorbing_state = 4),
    "zero probability.*patient\\(s\\) 999\\."
  )

  # The fitting function stops with this error before compiling or sampling.
  skip_if_not_installed("rstan")
  expect_error(
    blrm_markov_ml(formula, data = impossible, id_var = "id",
      absorbing_state = 4, missing_states = "non_absorbing"),
    "zero probability.*998, 999\\."
  )
})

test_that("cluster() terms are checked and removed for set-up", {
  data <- ml_oracle_data()
  formula <- y ~ time + tx + yprev + cluster(id)
  checked <- markov_ml_check_inputs(
    formula,
    data,
    id_var = "id",
    absorbing_state = 4,
    time_var = "time",
    p_var = "yprev",
    random_effects = TRUE
  )
  expect_true(checked$random_effects)
  expect_identical(
    attr(stats::terms(checked$formula), "term.labels"),
    c("time", "tx", "yprev")
  )
  expect_identical(environment(checked$formula), environment(formula))
  expect_error(
    markov_ml_check_inputs(y ~ time + tx + yprev + cluster(tx), data, "id",
      4, "time", "yprev", random_effects = TRUE),
    "must name the patient ID column"
  )
  expect_error(
    markov_ml_check_inputs(y ~ time + yprev + tx * cluster(id), data, "id",
      4, "time", "yprev", random_effects = TRUE),
    "one `cluster\\(\\)` term"
  )
})

test_that("normcppo centres and scales the scores over observed transitions", {
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 8)
  formula <- y ~ time + tx + yprev
  problem <- markov_ml_setup(formula, data, "id", "time", "yprev", 6, "any",
    list(), extra_vars = "time")
  fit <- markov_bayes_fit_object(problem, formula, ~time, bayes_cppo, TRUE)
  comp <- problem$complete_part
  observed <- rep(comp$agg_state, comp$agg_count)
  scores <- fit$cppo(problem$labels[observed])
  expect_equal(mean(scores), 0, tolerance = 1e-12)
  expect_equal(stats::sd(scores), 1, tolerance = 1e-12)
})

test_that("invalid input gives informative errors before sampling", {
  skip_if_not_installed("rstan")
  data <- ml_oracle_data()
  fit_bayes <- function(...) {
    blrm_markov_ml(y ~ time + tx + yprev, data = data, id_var = "id",
      absorbing_state = 4, ...)
  }
  expect_error(
    blrm_markov_ml(y ~ time + tx + yprev, data = data, absorbing_state = 4),
    "`id_var` is required"
  )
  expect_error(
    blrm_markov_ml(y ~ time + tx + yprev, data = data, id_var = "id"),
    "`absorbing_state` is required"
  )
  expect_error(
    blrm_markov_ml(y ~ time + tx, data = data, id_var = "id",
      absorbing_state = 4),
    "must appear on the right-hand side"
  )
  expect_error(fit_bayes(psigma = 3), "`psigma`")
  expect_error(fit_bayes(psigma = 2, rsdmean = 0), "`rsdmean` must be positive")
  expect_error(fit_bayes(p2_var = "ypprev"), "first-order")
  expect_error(fit_bayes(ppo = ~time), "constrained partial proportional odds")
  expect_error(fit_bayes(cppo = bayes_cppo), "require `ppo`")
  expect_error(
    fit_bayes(ppo = ~age, cppo = bayes_cppo),
    "`age` in `ppo` are not columns"
  )
  expect_error(fit_bayes(weights = rep(1, nrow(data))), "weights")
  expect_error(fit_bayes(backend = "cmdstan"), "backend")
  expect_error(fit_bayes(iter = 100, warmup = 100), "`warmup`")
  expect_error(fit_bayes(conc = -1), "`conc`")
  expect_error(
    orm_markov_ml(y ~ time + tx + yprev + cluster(id), data = data,
      id_var = "id", absorbing_state = 4),
    "Random effects"
  )
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
    formula <- stats::reformulate(
      attr(stats::terms(rhs), "term.labels"),
      response = response,
      env = environment()
    )
    conc <- case$conc %||% (1 / (0.8 + 0.35 * 6))
    re <- isTRUE(case$re)
    psigma <- case$psigma %||% 1
    rsdmean <- case$rsdmean %||% if (psigma == 1) 0 else 1
    rsdsd <- case$rsdsd %||% 1
    normcppo <- isTRUE(case$normcppo)
    setup <- bayes_setup(data, formula, ppo = case$ppo, conc = conc,
      random_effects = re, normcppo = normcppo, pcontrast = case$pcontrast,
      npcontrast = case$npcontrast, psigma = psigma, rsdmean = rsdmean,
      rsdsd = rsdsd)
    ours <- suppressMessages(
      rstan::sampling(model, data = setup$input$data, chains = 0)
    )

    # rmsb::blrm() returns its Stan data before adding contrast priors, so
    # take the data from a short run instead.
    rmsb_formula <- stats::reformulate(
      c(attr(stats::terms(rhs), "term.labels"), if (re) "cluster(id)"),
      response = case$rmsb_response %||% response,
      env = environment()
    )
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
  # blrm_markov() lets rmsb::blrm() find `cluster()` on the search path.
  withr::local_package("rmsb")
  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 5)
  # Patient IDs whose sorted order differs from their numeric order, so the
  # random-intercept columns are mapped to patients by ID.
  data$id <- paste0("p", data$id)
  formula <- y ~ rms::rcs(time, 3) + tx + yprev + cluster(id)
  rmsb_fit <- suppressMessages(suppressWarnings(blrm_markov(
    formula,
    ppo = ~time,
    cppo = bayes_cppo,
    data = data,
    id_var = "id",
    iter = 100,
    chains = 2,
    refresh = 0,
    loo = FALSE,
    sampling.args = list(seed = 1)
  )))

  # The Stan-scale draws of rmsb's own program, transformed as
  # blrm_markov_ml() transforms its draws, reproduce rmsb's draws.
  setup <- bayes_setup(data, y ~ rms::rcs(time, 3) + tx + yprev, ppo = ~time,
    conc = rmsb_fit$conc, random_effects = TRUE)
  stan_draws <- rstan::extract(rmsb_fit$rstan, pars = c("alpha", "beta", "tau"),
    permuted = FALSE)
  draws <- markov_bayes_draws_from_array(stan_draws, setup$input)
  expect_identical(colnames(draws), colnames(rmsb_fit$draws))
  expect_equal(draws, rmsb_fit$draws, tolerance = 1e-10)

  # A blrm_markov_ml() fit carrying the same draws, including the patients'
  # random intercepts as blrm_markov_ml() would label them.
  ml_fit <- suppressMessages(suppressWarnings(blrm_markov_ml(
    formula,
    ppo = ~time,
    cppo = bayes_cppo,
    data = data,
    id_var = "id",
    absorbing_state = 6,
    iter = 20,
    chains = 1,
    seed = 1
  )))
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
  data <- ml_oracle_data()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  conc <- 0.7
  for (rule in c("any", "non_absorbing", "gappy")) {
    for (ppo in list(NULL, ~time)) {
      setup <- bayes_setup(data, formula, ppo = ppo,
        missing_states = if (rule == "gappy") "any" else rule,
        absorbing_state = 4, conc = conc)
      if (rule == "gappy") setup <- bayes_gappy_setup(setup, conc)
      stan <- setup$input$data
      stanfit <- suppressMessages(
        rstan::sampling(model, data = stan, chains = 0)
      )
      K <- stan$k
      set.seed(2)
      for (i in 1:3) {
        u <- stats::rnorm(stan$p + stan$q + K - 1L, sd = 0.3)
        log_density <- rstan::log_prob(stanfit, u, adjust_transform = FALSE)
        pars <- rstan::constrain_pars(stanfit, u)
        coefficients <- markov_bayes_draws_from_array(
          bayes_stan_array(pars$alpha, pars$beta, pars$tau),
          setup$input
        )[1L, ]
        log_prior <- lgamma(K * conc) - K * lgamma(conc) +
          (conc - 1) * sum(log(pars$pi))
        expect_equal(
          log_density - log_prior,
          bayes_r_loglik(setup, coefficients),
          tolerance = 1e-10
        )
      }
    }
  }
})

test_that("with random intercepts the density matches the R recursion", {
  skip_without_stan()
  model <- markov_bayes_stan_model()
  data <- ml_oracle_data()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  conc <- 0.7
  cases <- expand.grid(rule = c("any", "non_absorbing", "gappy"),
    ppo = c(FALSE, TRUE), stringsAsFactors = FALSE)
  for (i_case in seq_len(nrow(cases))) {
    rule <- cases$rule[[i_case]]
    ppo <- if (cases$ppo[[i_case]]) ~time
    setup <- bayes_setup(data, formula, ppo = ppo,
      missing_states = if (rule == "gappy") "any" else rule,
      absorbing_state = 4, conc = conc, random_effects = TRUE)
    if (rule == "gappy") {
      setup <- bayes_gappy_setup(setup, conc, random_effects = TRUE)
    }
    stan <- setup$input$data
    stanfit <- suppressMessages(
      rstan::sampling(model, data = stan, chains = 0)
    )
    K <- stan$k
    set.seed(4)
    # Priors written with `~` drop constants, so compare up to a constant.
    gaps <- vapply(1:3, function(i) {
      u <- stats::rnorm(stan$p + stan$q + K - 1L + stan$Nc + 1L, sd = 0.3)
      log_density <- rstan::log_prob(stanfit, u, adjust_transform = FALSE)
      pars <- rstan::constrain_pars(stanfit, u)
      coefficients <- markov_bayes_draws_from_array(
        bayes_stan_array(pars$alpha, pars$beta, pars$tau),
        setup$input
      )[1L, ]
      log_prior <- lgamma(K * conc) - K * lgamma(conc) +
        (conc - 1) * sum(log(pars$pi)) -
        0.5 * sum(pars$gamma_raw^2) -
        2.5 * log1p(pars$sigmag^2 / 4)
      log_density - log_prior -
        bayes_r_loglik(setup, coefficients, gamma = pars$gamma)
    }, numeric(1))
    expect_equal(gaps, rep(gaps[[1L]], 3L), tolerance = 1e-10)
  }
})

test_that("a fitted model works with the SOP functions", {
  skip_without_stan()
  data <- ml_sim_missing(n_patients = 60, follow_up_time = 6, missing = 0.2)
  fit <- suppressWarnings(blrm_markov_ml(
    y ~ time + tx + yprev,
    ppo = ~time,
    cppo = bayes_cppo,
    data = data,
    id_var = "id",
    absorbing_state = 6,
    iter = 300,
    chains = 2,
    seed = 1
  ))
  expect_s3_class(fit, c("blrm_markov_ml", "blrm"))
  expect_identical(markov_model_fit_wrapper(fit), "blrm_markov_ml")
  expect_equal(nrow(fit$draws), 300L)
  expect_identical(
    colnames(fit$draws),
    c(names(fit$ml$start_coefficients), "time x f(y)")
  )
  expect_identical(rownames(fit$diagnostics), colnames(fit$draws))
  expect_equal(coef(fit), colMeans(fit$draws))
  expect_output(print(fit), "blrm_markov_ml")

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
  fit <- suppressWarnings(blrm_markov_ml(
    y ~ time + tx + yprev + cluster(id),
    data = data,
    id_var = "id",
    absorbing_state = 6,
    iter = 300,
    chains = 2,
    seed = 2
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
