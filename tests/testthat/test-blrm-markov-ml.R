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
  random_effects = FALSE
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
    normcppo = FALSE
  )
  input <- markov_bayes_stan_data(problem, fit, conc, NULL, NULL,
    random_effects = random_effects)
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
    used <- sort(unique(c(
      setup$problem$complete_part$agg_design,
      as.vector(setup$problem$chain$DI)
    )))
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
  expect_identical(
    as.integer(stan$seg_cluster),
    problem$chain$patient[stan$seg_start]
  )
  expect_length(stan$rsdmean, 1L)
  expect_length(stan$rsdsd, 1L)
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
  formula <- y ~ rms::rcs(time, 3) + tx + yprev
  cases <- list(
    list(ppo = NULL, re = FALSE),
    list(ppo = ~time, re = FALSE),
    list(ppo = ~time, re = TRUE)
  )
  for (case in cases) {
    ppo <- case$ppo
    setup <- bayes_setup(data, formula, ppo = ppo,
      conc = 1 / (0.8 + 0.35 * 6), random_effects = case$re)
    ours <- suppressMessages(
      rstan::sampling(model, data = setup$input$data, chains = 0)
    )
    rmsb_formula <- if (case$re) {
      stats::update(formula, . ~ . + cluster(id))
    } else {
      formula
    }
    standata <- if (is.null(ppo)) {
      rmsb::blrm(rmsb_formula, data = data, standata = TRUE)
    } else {
      rmsb::blrm(rmsb_formula, ppo = ppo, cppo = bayes_cppo, data = data,
        standata = TRUE)
    }
    stanmodels <- utils::getFromNamespace("stanmodels", "rmsb")
    theirs <- suppressMessages(
      rstan::sampling(stanmodels$lrmconppo, data = standata, chains = 0)
    )
    n_par <- setup$input$data$p + setup$input$data$q + setup$problem$K - 1L +
      if (case$re) setup$problem$N + 1L else 0L
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

test_that("the likelihood equals the R recursion with missing and range states", {
  skip_without_stan()
  model <- markov_bayes_stan_model()
  data <- ml_oracle_data()
  formula <- rms::Ocens(y_lo, y_hi) ~ time + tx + yprev
  conc <- 0.7
  for (rule in c("any", "non_absorbing")) {
    for (ppo in list(NULL, ~time)) {
      setup <- bayes_setup(data, formula, ppo = ppo, missing_states = rule,
        absorbing_state = 4, conc = conc)
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
  for (ppo in list(NULL, ~time)) {
    setup <- bayes_setup(data, formula, ppo = ppo, absorbing_state = 4,
      conc = conc, random_effects = TRUE)
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
