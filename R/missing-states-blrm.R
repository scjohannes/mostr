# Fitting when some states are not observed: `engine = "blrm"`.
#
# rmsb's Stan programs take one row per transition, so they can't sum over
# the states of unobserved visits. mostr samples from its own Stan program
# (`inst/stan/po_paths.stan`) with the same priors, and puts the draws into
# a `blrm` object so that everything downstream reads an ordinary rmsb fit.
# The program samples the intercepts and coefficients on a scale taken from
# the maximum likelihood fit of the model (`markov_path_stan_scale()`), not
# in rmsb's parameters, because that needs far fewer gradient evaluations.

mostm_fit_blrm_paths <- function(setup, ranges, path, ...) {
  args <- list(...)
  backend <- args$backend %||% getOption("rmsb.backend", "rstan")
  backend <- match.arg(backend, c("rstan", "cmdstan"))
  labels <- ranges$labels
  layout <- mostm_path_layout(setup, ranges, path)
  # No stretch: rmsb's own fit to the observed transitions is the fit (see
  # `mostm_fit_orm_paths()`), and nothing needs to be compiled.
  if (!layout$n_steps_total) {
    return(mostm_fit_blrm(markov_observed_transition_setup(setup, ranges), ...))
  }
  set_up <- markov_path_setup(setup, ranges, path)
  mostm_check_setup_states(ranges, set_up, layout, path$p_var)
  observed_setup <- set_up$setup
  # Rows of `setup$data` that hold the first step of each stretch.
  stretch_data_rows <- set_up$subset_index[layout$step_data_rows[layout$is_first]]
  cluster <- markov_path_cluster_values(setup)
  if (!is.null(cluster)) {
    markov_check_path_clusters(
      cluster,
      setup$data[[path$id_var]],
      c(which(set_up$use), set_up$subset_index[layout$step_data_rows]),
      path$id_var
    )
  }

  # The model structure comes from `rmsb::blrm()` itself, fitted to the
  # set-up rows (`markov_setup_rows()`): design matrix, state labels and
  # every field the prediction code reads. Its draws are replaced below, so a
  # few iterations of one chain are enough; they take about a second with
  # rmsb's compiled program, which is why the backend is rstan here whatever
  # the user chose.
  structure_args <- utils::modifyList(
    args,
    list(
      iter = 20,
      warmup = 10,
      chains = 1,
      loo = FALSE,
      method = "sampling",
      refresh = 0,
      progress = "",
      backend = "rstan"
    )
  )
  structure_args[c("sampling.args", "file", "ppairs", "inits")] <- NULL
  structure_call <- c(
    list(
      formula = observed_setup$formula,
      data = observed_setup$data,
      x = TRUE,
      y = TRUE
    ),
    structure_args
  )
  structure_call$subset <- observed_setup$subset_rows
  # Ten warm-up iterations always draw sampler warnings (R-hat, effective
  # sample size); they describe draws that are thrown away.
  fit <- suppressWarnings(do.call(rmsb::blrm, structure_call))
  markov_check_observed_states(fit$ylevels, labels)
  fit_frame <- markov_rms_model_frame(
    observed_setup$formula,
    observed_setup$data,
    setup$dots,
    subset = observed_setup$subset_rows,
    eval_env = setup$eval_env
  )
  fit_data <- markov_align_model_data(setup$data, fit_frame)
  keep <- markov_setup_recorded_rows(
    fit_data,
    setup$data,
    set_up$recorded,
    nrow(fit$x)
  )
  fit_data <- fit_data[keep, , drop = FALSE]

  # Random intercepts: one cluster number per row, over the clusters of the
  # fully observed transitions and of the stretches. The structure fit
  # numbered only the clusters of its own rows, and a patient whose
  # transitions all touch unobserved visits has none there.
  clusters <- NULL
  if (!is.null(cluster)) {
    observed_data_rows <- match(rownames(fit_data), rownames(setup$data))
    clusters <- markov_path_cluster_codes(
      cluster$values,
      c(observed_data_rows, stretch_data_rows)
    )
    clusters$obs <- clusters$code[observed_data_rows]
  }

  plan <- markov_path_plan_design(
    layout,
    function(rows) blrm_design_matrix(fit, rows),
    cluster = if (!is.null(clusters)) clusters$code[set_up$subset_index]
  )

  x_setup <- as.matrix(fit$x)
  x <- x_setup[keep, , drop = FALSE]
  n_states <- length(labels)
  observed_state <- match(
    as.character(markov_blrm_observed_states(fit))[keep],
    labels
  )
  assembled <- markov_blrm_path_stan_data(
    x_setup = x_setup,
    keep = keep,
    observed_state = observed_state,
    plan = plan,
    n_states = n_states,
    contrast = fit$Contrast,
    conc = args$conc %||% fit$conc %||% (1 / (0.8 + 0.35 * max(n_states, 3))),
    pcontrast = args$pcontrast,
    clusters = clusters,
    psigma = fit$psigma %||% 1,
    rsdmean = fit$rsdmean %||% 0,
    rsdsd = fit$rsdsd %||% 1
  )
  stan_data <- assembled$data
  transform <- assembled$transform
  rinv <- assembled$rinv

  # The scale the program samples the intercepts and coefficients on
  # (`theta_hat`, `theta_L`), from the maximum likelihood fit of the model
  # without random intercepts. Newton's method starts from the structure
  # fit's draws and, if that fails, from the state frequencies.
  scale_time <- system.time(
    stan_scale <- markov_path_stan_scale(
      stan_data,
      starts = list(
        structure_fit = markov_path_stan_start_from_draws(
          fit$draws,
          n_states,
          rinv,
          transform$xbar
        ),
        state_frequencies = markov_path_stan_marginal_start(stan_data)
      )
    )
  )
  stan_data$theta_hat <- stan_scale$theta_hat
  stan_data$theta_L <- stan_scale$theta_L

  iter <- args$iter %||% 2000
  warmup <- args$warmup %||% (iter / 2)
  chains <- args$chains %||% 4
  refresh <- args$refresh %||% 0
  model <- markov_path_stan_model(backend)
  pars <- c("alpha", "beta", if (!is.null(clusters)) c("gamma", "sigmag"))
  sampler_args <- markov_path_sampler_args(
    backend,
    stan_data = stan_data,
    pars = pars,
    iter = iter,
    warmup = warmup,
    chains = chains,
    refresh = refresh,
    sampling_args = args$sampling.args,
    random_intercepts = !is.null(clusters)
  )
  elapsed <- system.time(
    sampled <- if (backend == "rstan") {
      do.call(rstan::sampling, c(list(model), sampler_args))
    } else {
      do.call(model$sample, sampler_args)
    }
  )
  if (backend == "rstan" && sampled@mode == 2L) {
    mostr_abort("The Stan sampler failed for the model with unobserved states.")
  }

  pattern <- paste0("^(", paste(pars, collapse = "|"), ")\\[")
  if (backend == "rstan") {
    raw <- as.matrix(sampled)
    parameter_names <- colnames(raw)
  } else {
    all_draws <- sampled$draws()
    parameter_names <- dimnames(all_draws)[[3L]]
    parameter_names <- parameter_names[grepl(pattern, parameter_names)]
    raw <- matrix(
      all_draws[, , parameter_names],
      ncol = length(parameter_names),
      dimnames = list(NULL, parameter_names)
    )
  }
  alpha_names <- parameter_names[grepl("^alpha\\[", parameter_names)]
  beta_names <- parameter_names[grepl("^beta\\[", parameter_names)]
  betas <- raw[, beta_names, drop = FALSE] %*% t(rinv)
  # The intercepts were sampled for centred predictors.
  alphas <- raw[, alpha_names, drop = FALSE] -
    as.numeric(betas %*% transform$xbar)
  colnames(alphas) <- colnames(fit$draws)[seq_len(n_states - 1L)]
  colnames(betas) <- colnames(x)
  draws <- cbind(alphas, betas)

  fit$draws <- draws
  # Downstream code reads `x` and `y` as the fully observed transitions, as
  # for the orm fit (`markov_orm_keep_rows()`).
  fit$x <- x
  fit$y <- markov_keep_response_rows(fit$y, keep)
  fit$freq <- stats::setNames(
    as.integer(table(factor(
      markov_blrm_observed_states(fit),
      levels = fit$ylevels
    ))),
    fit$ylevels
  )
  fit$N <- nrow(x)
  fit$param <- rbind(
    mean = colMeans(draws),
    median = apply(draws, 2L, stats::median)
  )
  cluster_parameters <- character()
  if (!is.null(clusters)) {
    # As `rmsb::blrm()` stores them: the posterior draws of the SD of the
    # random intercepts as `omega`, the posterior medians of the intercepts
    # as `gammas` (named "gamma[c]"), and the clusters. `clusterInfo$cluster`
    # holds the cluster of each fully observed transition (the rows of `x`)
    # as a factor whose levels are the clusters in the order of `gamma`, so
    # that a cluster seen only in stretches, which has no row of `x`, keeps
    # its place: `include_re = TRUE` matches clusters to `gamma` by these
    # levels.
    gamma_names <- parameter_names[grepl("^gamma\\[", parameter_names)]
    fit$omega <- cbind(sigmag = raw[, "sigmag[1]"])
    fit$gammas <- apply(raw[, gamma_names, drop = FALSE], 2L, stats::median)
    fit$clusterInfo <- list(
      cluster = factor(clusters$levels[clusters$obs], levels = clusters$levels),
      n = clusters$n,
      name = fit$clusterInfo$name %||% cluster$name
    )
    cluster_parameters <- "sigmag[1]"
  }
  fit$diagnostics <- if (backend == "rstan") {
    tryCatch(
      rstan::summary(
        sampled,
        pars = c(alpha_names, beta_names, cluster_parameters),
        probs = NULL
      )$summary[, c("n_eff", "Rhat")],
      error = function(e) NULL
    )
  } else {
    list(
      diagnostic_summary = sampled$diagnostic_summary(),
      summary = as.data.frame(sampled$summary(
        variables = setdiff(pars, "gamma")
      ))
    )
  }
  fit$alphas <- alpha_names
  fit$betas <- beta_names
  fit$xbar <- transform$xbar
  fit$iter <- iter
  fit$chains <- chains
  fit$backend <- backend
  fit$sampling_time <- unname(elapsed[["elapsed"]])
  fit$stancode <- markov_path_stan_code()
  fit$rstan <- sampled
  # Leave-one-out and posterior-mode results of the structure fit describe
  # other draws; none are computed for this model.
  fit["loo"] <- list(NULL)
  fit["opt"] <- list(NULL)
  fit["datahash"] <- list(NULL)

  patient_ids <- unique(c(
    as.character(fit_data[[path$id_var]]),
    plan$instances$id
  ))
  attr(fit, "markov_path_fit") <- list(
    ids = patient_ids,
    n_observed_transitions = nrow(x),
    n_stretches = nrow(plan$instances),
    n_distinct_stretches = length(plan$stretch_n_steps),
    n_unobserved_visits = plan$n_steps_total,
    # The sampler's parametrisation: coefficients beta_s of the centred,
    # transformed predictors, beta = rinv %*% beta_s, intercepts for
    # predictors centred at `fit$xbar`. It comes from the set-up rows, which
    # are not stored, so it is kept for checking the sampled density.
    rinv = rinv,
    # The scale the intercepts and coefficients were sampled on
    # (`markov_path_stan_scale()`): whether it came from the maximum
    # likelihood fit (`maximum_likelihood`; if not, the sampler ran on a
    # rough scale and was only slower), which start Newton's method
    # converged from, its iterations and time, and `theta_hat` and
    # `theta_L`, which the tests need to map points to the sampled
    # parameters.
    stan_scale = c(
      stan_scale,
      list(seconds = unname(scale_time[["elapsed"]]))
    )
  )
  list(fit = fit, fit_data = fit_data, fitted_ids = patient_ids)
}

# The data of `inst/stan/po_paths.stan`, with the transformation of the
# predictors it is sampled in.
#
# `x_setup` is the design matrix of the set-up rows (`fit$x` of the structure
# fit), `keep` marks its fully observed transitions and `observed_state`
# gives their states as positions in the labels. `plan` comes from
# `markov_path_plan_design()`. `clusters`, for a model with random
# intercepts, holds `n` (the number of clusters) and `obs` (the cluster of
# each fully observed transition); the stretches' clusters are
# `plan$stretch_cluster`. `psigma`, `rsdmean` and `rsdsd` are rmsb's prior on
# the SD of the random intercepts.
#
# Returns `data` (the Stan data), `transform` (`rmsb::selectedQr()` of the
# set-up rows) and `rinv` (its `Rinv`, scaled as below).
markov_blrm_path_stan_data <- function(
  x_setup,
  keep,
  observed_state,
  plan,
  n_states,
  contrast = NULL,
  conc,
  pcontrast = NULL,
  clusters = NULL,
  psigma = 1,
  rsdmean = 0,
  rsdsd = 1
) {
  x <- x_setup[keep, , drop = FALSE]
  n_predictors <- ncol(x)
  # Stan samples coefficients of centred, orthogonalised predictors, as rmsb
  # does (`rmsb::selectedQr()`); the stretches' design rows go through the
  # same transformation, and the draws are transformed back by the caller.
  # The transformation is computed from all set-up rows, which have a
  # full-rank design whenever the set-up fit worked; the fully observed
  # transitions alone need not (a previous state that is always followed by
  # an unobserved visit has a column of zeros there). The Dirichlet prior on
  # the intercepts is therefore placed at the mean predictors of the set-up
  # rows.
  transform <- rmsb::selectedQr(x_setup, center = TRUE)
  # `selectedQr()` returns orthonormal columns (length 1, so SD
  # 1 / sqrt(n - 1)), whose coefficients have posterior SDs of order 1 to 10
  # while the intercepts' unconstrained parameters have SDs of order
  # n^-1/2. Stan starts from a unit metric and needed trajectories of 100 or
  # more leapfrog steps until it had adapted: on 300 patients x 14 visits,
  # 6 states, 4 chains x 1000 iterations, warm-up took about 15,900 leapfrog
  # steps and 56 s per chain with these columns, and about 4,300 steps and
  # 16 s with columns of SD 1 (Nov 2026).
  # So the columns are multiplied by sqrt(n - 1), and `Rinv` by the same
  # factor: with Q the orthonormal columns and R^-1 their `Rinv`, the linear
  # predictor (x - xbar)' beta = Q beta_q with beta = R^-1 beta_q equals
  # (c Q) (beta_q / c), and beta = (c R^-1) (beta_q / c). The sampled
  # coefficients are beta_q / c; the scaled `Rinv` turns them, the
  # stretches' design rows and the contrasts of the contrast prior into the
  # original scale. The model has no prior on these coefficients other than
  # the contrast prior, which is written in the original coefficients
  # (`C beta`) and so is unchanged, and the Dirichlet prior on the
  # intercepts is defined at the mean predictors, which the scaling does not
  # move: the posterior of the original coefficients is the same, only the
  # sampler's parametrisation differs.
  column_scale <- sqrt(max(nrow(x_setup) - 1, 1))
  rinv <- transform$Rinv * column_scale
  # The design row for previous state j is x_gap + x_prev[j] (additive
  # layout) or x_gap alone; the transformation is linear, so the centring
  # goes with x_gap and x_prev only gets `rinv`.
  x_gap <- if (is.null(plan$x_gap)) {
    matrix(0, 0L, n_predictors)
  } else {
    sweep(unname(plan$x_gap), 2L, transform$xbar) %*% rinv
  }
  x_prev <- if (is.null(plan$x_prev)) {
    matrix(0, 0L, n_predictors)
  } else {
    unname(plan$x_prev) %*% rinv
  }
  n_contrasts <- if (is.null(contrast)) 0L else nrow(contrast)
  # Observed transitions with the same predictors and state are counted.
  # With random intercepts, transitions of different clusters have different
  # linear predictors, so the cluster is part of what must be equal.
  observed_row <- markov_distinct_rows(
    if (is.null(clusters)) x else cbind(x, clusters$obs)
  )
  n_observed_rows <- if (length(observed_row)) max(observed_row) else 0L
  cell_count <- table(
    factor(observed_row, levels = seq_len(n_observed_rows)),
    factor(observed_state, levels = seq_len(n_states))
  )
  cells <- which(cell_count > 0)
  cell_row <- (cells - 1L) %% n_observed_rows + 1L
  first_of_row <- !duplicated(observed_row)
  x_observed <- unname(unclass(transform$X))[keep, , drop = FALSE][
    first_of_row,
    ,
    drop = FALSE
  ] * column_scale
  # The Stan program reads matrices as arrays stored by row.
  by_row <- function(m) as.array(as.numeric(t(m)))
  n_clusters <- if (is.null(clusters)) 0L else as.integer(clusters$n)
  data <- list(
    p = n_predictors,
    k = n_states,
    n_obs = length(cells),
    x_obs = by_row(x_observed[cell_row, , drop = FALSE]),
    y_obs = as.array((cells - 1L) %/% n_observed_rows + 1L),
    w_obs = as.array(as.numeric(cell_count[cells])),
    G = nrow(x_gap),
    x_gap = by_row(x_gap),
    has_prev = as.integer(nrow(x_prev) > 0L),
    x_prev = by_row(x_prev),
    S = length(plan$stretch_n_steps),
    T = length(plan$step_lo),
    first_row = as.array(plan$stretch_first_row),
    first_prev = as.array(plan$stretch_first_prev),
    first_step = as.array(plan$stretch_first_step),
    n_steps = as.array(plan$stretch_n_steps),
    weight = as.array(plan$stretch_weight),
    lo = as.array(plan$step_lo),
    hi = as.array(plan$step_hi),
    step_rows = as.array(as.integer(t(plan$step_rows))),
    absorbing = as.array(as.integer(seq_len(n_states) %in% plan$absorb)),
    Nc = n_clusters,
    obs_cluster = if (n_clusters) {
      as.array(as.integer(clusters$obs[first_of_row][cell_row]))
    } else {
      integer()
    },
    stretch_cluster = if (n_clusters) {
      as.array(as.integer(plan$stretch_cluster))
    } else {
      integer()
    },
    psigma = as.integer(psigma),
    rsdmean = if (n_clusters) as.array(rsdmean) else numeric(),
    rsdsd = if (n_clusters && psigma == 1) as.array(rsdsd) else numeric(),
    conc = conc,
    cn = n_contrasts,
    C = if (n_contrasts) {
      contrast %*% rinv
    } else {
      matrix(0, 0L, n_predictors)
    },
    cmus = as.array(rep_len(pcontrast$mu %||% 0, n_contrasts)),
    csds = as.array(rep_len(pcontrast$sd %||% 100, n_contrasts))
  )
  list(data = data, transform = transform, rinv = rinv)
}

# The scale on which `inst/stan/po_paths.stan` samples the intercepts and
# coefficients: theta = theta_hat + theta_L z, z sampled. theta holds the
# log-odds u[j] = log(pi[j] / pi[k]), j < k, of the state probabilities pi at
# the centred predictors, then the coefficients beta of the Stan data's
# transformed predictors. theta_hat is the maximum likelihood estimate of the
# model without random intercepts, computed with the same C++ likelihood
# from the Stan data's own design rows, and theta_L the lower Cholesky
# factor of its covariance (the inverse of the observed information, moved
# from the intercepts to u with `markov_path_log_odds()`).
#
# Why: in rmsb's parameters Stan starts from a unit metric while the
# coefficients' posterior SDs are about 20 times smaller than those of the
# random intercepts and strongly correlated with the intercepts, so early
# warm-up runs with tiny steps, and a diagonal metric never learns the
# correlations. On 1,000 patients x 30 visits, half the states unobserved,
# 2 chains x (1000 + 1000), this scale cut the gradient evaluations from
# 92,000 to 28,000 without `cluster()` and from 175,000 to 81,000 with
# `cluster(id)`, both at adapt_delta 0.8 (dev/decisions.md, 2026-11-15).
#
# The posterior is unchanged, also with weights or with the contrast prior:
# z -> theta is affine, so its Jacobian is a constant, and the program adds
# the Jacobian of u -> pi for the Dirichlet prior. Weights enter the Stan
# data's counts (`w_obs`, `weight`) as before, and the contrast prior is
# written on beta, which the program computes from theta.
#
# Newton's method (`markov_newton_maximise()`) runs from each of `starts`
# in turn (lists of `alpha`, intercepts at the centred predictors, and
# `beta`; NULL entries are skipped) until it converges. If none converges,
# the information at the maximum is not positive definite, or anything is
# not finite, theta_hat is the first start that gives finite log-odds and
# theta_L the identity: the program then samples the same posterior on a
# rough scale, only more slowly, so the fit goes ahead without a warning.
# Returns `theta_hat`, `theta_L`, `maximum_likelihood` (whether the scale
# came from the maximum likelihood fit), `start` (the name of the start
# used) and `iterations`.
markov_path_stan_scale <- function(data, starts, maxit = 50L) {
  n_intercepts <- data$k - 1L
  n_theta <- n_intercepts + data$p
  intercepts <- seq_len(n_intercepts)
  evaluate <- markov_path_stan_loglik_function(data)
  weights <- list(stretch_weight = as.numeric(data$weight))
  score_of <- function(result) markov_path_total_score(result, weights)
  if (is.null(names(starts))) {
    names(starts) <- paste0("start_", seq_along(starts))
  }
  starts <- Filter(Negate(is.null), starts)
  for (name in names(starts)) {
    theta <- c(starts[[name]]$alpha, starts[[name]]$beta)
    scale <- tryCatch(
      markov_path_stan_ml_scale(theta, evaluate, score_of, n_intercepts, maxit),
      error = function(e) NULL
    )
    if (!is.null(scale)) {
      return(c(scale, list(maximum_likelihood = TRUE, start = name)))
    }
  }
  for (name in names(starts)) {
    log_odds <- markov_path_log_odds(starts[[name]]$alpha)
    theta_hat <- c(log_odds$u, starts[[name]]$beta)
    if (length(theta_hat) == n_theta && all(is.finite(theta_hat))) {
      return(list(
        theta_hat = theta_hat,
        theta_L = diag(n_theta),
        iterations = NA_integer_,
        maximum_likelihood = FALSE,
        start = name
      ))
    }
  }
  # No start gave finite log-odds: equal state probabilities and no effects.
  list(
    theta_hat = numeric(n_theta),
    theta_L = diag(n_theta),
    iterations = NA_integer_,
    maximum_likelihood = FALSE,
    start = "none"
  )
}

# The maximum likelihood scale from one start, or NULL when it can't be
# formed (see `markov_path_stan_scale()`).
markov_path_stan_ml_scale <- function(theta, evaluate, score_of, n_intercepts, maxit) {
  if (length(theta) == 0L || !all(is.finite(theta))) {
    return(NULL)
  }
  current <- evaluate(theta, 2L)
  if (!is.finite(current$loglik)) {
    return(NULL)
  }
  maximum <- markov_newton_maximise(
    theta,
    evaluate,
    score_of,
    current = current,
    maxit = maxit
  )
  if (!maximum$converged || !all(is.finite(maximum$theta))) {
    return(NULL)
  }
  information <- -maximum$result$hessian
  if (!all(is.finite(information))) {
    return(NULL)
  }
  # A maximum must have a positive definite information; without one there
  # is no covariance to take the scale from.
  factor <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(factor)) {
    return(NULL)
  }
  covariance <- chol2inv(factor)
  intercepts <- seq_len(n_intercepts)
  log_odds <- markov_path_log_odds(maximum$theta[intercepts])
  jacobian <- diag(length(theta))
  jacobian[intercepts, intercepts] <- log_odds$jacobian
  covariance <- jacobian %*% covariance %*% t(jacobian)
  # Rounding leaves the product slightly asymmetric, and `chol()` reads only
  # the upper triangle.
  covariance <- (covariance + t(covariance)) / 2
  theta_hat <- c(log_odds$u, maximum$theta[-intercepts])
  upper <- tryCatch(chol(covariance), error = function(e) NULL)
  if (is.null(upper) || !all(is.finite(upper)) || !all(is.finite(theta_hat))) {
    return(NULL)
  }
  list(
    theta_hat = theta_hat,
    theta_L = t(upper),
    iterations = maximum$iterations
  )
}

# The log-likelihood of the model without random intercepts as a function
# of theta = (alpha, beta) on the scale of the Stan data `data`, with its
# derivatives (`cpp_path_loglik()`, the program's own C++ likelihood). The
# observed rows of a model with random intercepts are split by cluster in
# the Stan data, which changes nothing here.
markov_path_stan_loglik_function <- function(data) {
  p <- data$p
  k <- data$k
  by_row <- function(values, ncol) {
    matrix(as.numeric(values), ncol = ncol, byrow = TRUE)
  }
  x_obs <- by_row(data$x_obs, p)
  x_gap <- by_row(data$x_gap, p)
  x_prev <- if (data$has_prev == 1L) by_row(data$x_prev, p) else matrix(0, 0L, p)
  step_rows <- matrix(as.integer(data$step_rows), ncol = k, byrow = TRUE)
  n_intercepts <- k - 1L
  function(theta, order) {
    cpp_path_loglik(
      alpha = theta[seq_len(n_intercepts)],
      beta = theta[-seq_len(n_intercepts)],
      x_obs = x_obs,
      y_obs = as.integer(data$y_obs),
      w_obs = as.numeric(data$w_obs),
      group_obs = rep(1L, data$n_obs),
      n_groups = 1L,
      x_gap = x_gap,
      x_prev = x_prev,
      absorb = which(as.integer(data$absorbing) == 1L),
      stretch_first_row = as.integer(data$first_row),
      stretch_first_prev = as.integer(data$first_prev),
      stretch_first_step = as.integer(data$first_step),
      stretch_n_steps = as.integer(data$n_steps),
      stretch_weight = as.numeric(data$weight),
      step_lo = as.integer(data$lo),
      step_hi = as.integer(data$hi),
      step_rows = step_rows,
      cluster_shift = double(),
      obs_cluster = integer(),
      stretch_cluster = integer(),
      order = as.integer(order)
    )
  }
}

# The log-odds u[j] = log(pi[j] / pi[k]), j < k, of the state probabilities
# pi given by intercepts alpha[j] = logit(P(Y >= j + 1)), and the matrix of
# their derivatives, jacobian[i, j] = du[i] / dalpha[j].
#
# pi[i] = plogis(a) - plogis(b), with a = alpha[i - 1] (Inf for i = 1) and
# b = alpha[i] (-Inf for i = k), is computed on the log scale as
# log(plogis(a)) + log(plogis(-b)) + log(1 - exp(b - a)), because at the
# centred predictors a state can have pi of about 1e-7 and the difference of
# two probabilities near 1 would lose most of its digits. With
# c(a) = plogis(a) and c'(a) = c(a) (1 - c(a)), alpha[j] is the lower end of
# pi[j] and the upper end of pi[j + 1]:
#   dlog(pi[j]) / dalpha[j] = -c'(alpha[j]) / pi[j],
#   dlog(pi[j + 1]) / dalpha[j] = c'(alpha[j]) / pi[j + 1],
# and du[i] / dalpha[j] = dlog(pi[i]) / dalpha[j] - dlog(pi[k]) / dalpha[j].
# Intercepts that do not decrease give a pi <= 0 and non-finite values.
markov_path_log_odds <- function(alpha) {
  alpha <- as.numeric(alpha)
  m <- length(alpha)
  k <- m + 1L
  upper <- c(Inf, alpha)
  lower <- c(alpha, -Inf)
  log_pi <- suppressWarnings(
    stats::plogis(upper, log.p = TRUE) +
      stats::plogis(-lower, log.p = TRUE) +
      log(-expm1(lower - upper))
  )
  log_slope <- stats::plogis(alpha, log.p = TRUE) +
    stats::plogis(-alpha, log.p = TRUE)
  d_log_pi <- matrix(0, k, m)
  for (j in seq_len(m)) {
    d_log_pi[j, j] <- -exp(log_slope[j] - log_pi[j])
    d_log_pi[j + 1L, j] <- exp(log_slope[j] - log_pi[j + 1L])
  }
  list(
    u = log_pi[-k] - log_pi[k],
    jacobian = d_log_pi[-k, , drop = FALSE] -
      matrix(d_log_pi[k, ], m, m, byrow = TRUE)
  )
}

# A start for `markov_path_stan_scale()` from the structure fit's draws: their
# medians, moved to the Stan data's scale (intercepts at the predictors
# centred at `xbar`, coefficients of the transformed predictors, beta_s =
# rinv^-1 beta). The structure fit has only a few short draws, but Newton's
# method needs only a point near the maximum. NULL when the draws do not
# have the expected columns.
markov_path_stan_start_from_draws <- function(draws, n_states, rinv, xbar) {
  n_intercepts <- n_states - 1L
  if (is.null(draws) || ncol(draws) != n_intercepts + ncol(rinv) || !nrow(draws)) {
    return(NULL)
  }
  center <- apply(draws, 2L, stats::median)
  alpha <- unname(center[seq_len(n_intercepts)])
  beta <- unname(center[-seq_len(n_intercepts)])
  beta_scaled <- tryCatch(as.numeric(solve(rinv, beta)), error = function(e) NULL)
  if (is.null(beta_scaled)) {
    return(NULL)
  }
  list(alpha = alpha + sum(beta * xbar), beta = beta_scaled)
}

# A start that needs no fit: no effects, and intercepts from the frequencies
# of the states of the fully observed transitions, with half a count added
# to each state so that a state seen only after unobserved visits gets a
# finite intercept. This is the usual start of an ordinal regression.
markov_path_stan_marginal_start <- function(data) {
  k <- data$k
  counts <- vapply(
    seq_len(k),
    function(state) sum(as.numeric(data$w_obs)[as.integer(data$y_obs) == state]),
    numeric(1)
  ) + 0.5
  share <- counts / sum(counts)
  at_least <- rev(cumsum(rev(share)))[-1L]
  list(alpha = stats::qlogis(at_least), beta = numeric(data$p))
}

# The arguments of `rstan::sampling()` (without the model) or of cmdstanr's
# `$sample()`, with the user's `sampling.args` added.
#
# With random intercepts, `adapt_delta` (the target acceptance rate of the
# step size adaptation) defaults to 0.9 instead of Stan's 0.8. With a true
# SD of the random intercepts near 0, 300 patients and Stan's 0.8, three
# seeds gave 0, 5 and 2 divergent transitions after warm-up, and one chain's
# step size collapsed; with 0.9, 0, 1 and 0, at 40-60% more gradient
# evaluations than with 0.8, still fewer than rmsb's parametrisation needed;
# at 1,000 patients x 30 visits with SD 0.5, 136,000 against 82,000 at 0.8
# and 175,000 for rmsb's parametrisation at 0.8 (dev/decisions.md,
# 2026-11-15). Without random intercepts Stan's default
# stays. A value the user gives (`sampling.args = list(control =
# list(adapt_delta = ))` for rstan, `list(adapt_delta = )` for cmdstan) is
# kept, and the rest of a user's `control` list too.
#
# The number of chains run in parallel is left to the backend, as
# `rmsb::blrm()` leaves it: `rstan::sampling()` uses
# `getOption("mc.cores", 1)` cores and cmdstanr's `$sample()` as many
# parallel chains, and `sampling.args` (`cores`, `parallel_chains`)
# overrides either.
markov_path_sampler_args <- function(
  backend,
  stan_data,
  pars,
  iter,
  warmup,
  chains,
  refresh,
  sampling_args = NULL,
  random_intercepts = FALSE
) {
  sampling_args <- as.list(sampling_args)
  adapt_delta <- 0.9
  if (backend == "rstan") {
    if (random_intercepts) {
      control <- as.list(sampling_args$control)
      if (is.null(control$adapt_delta)) {
        control$adapt_delta <- adapt_delta
      }
      sampling_args$control <- control
    }
    return(c(
      list(
        data = stan_data,
        pars = pars,
        iter = iter,
        warmup = warmup,
        chains = chains,
        refresh = refresh
      ),
      sampling_args
    ))
  }
  if (random_intercepts && is.null(sampling_args$adapt_delta)) {
    sampling_args$adapt_delta <- adapt_delta
  }
  c(
    list(
      data = stan_data,
      iter_sampling = iter - warmup,
      iter_warmup = warmup,
      chains = chains,
      refresh = refresh
    ),
    sampling_args
  )
}

# The cluster variable of a `cluster()` term in the formula, evaluated on
# every row of `setup$data`, or NULL without one. `mostm()` passes only a
# bare `cluster()` to rmsb (`mostm_check_random_effects()`).
markov_path_cluster_values <- function(setup) {
  calls <- formula_rhs_calls(setup$formula, "cluster")
  if (!length(calls)) {
    return(NULL)
  }
  expr <- calls[[1L]][[2L]]
  values <- eval(expr, setup$data, setup$eval_env)
  if (length(values) != nrow(setup$data)) {
    mostr_abort(
      "The variable in `cluster()` must have one value per row of `data`."
    )
  }
  list(name = deparse1(expr), values = values)
}

# A run of unobserved visits is one likelihood term with one random
# intercept, and the stretches and the fully observed transitions of a
# patient are laid out per patient, so the cluster must be known and the same
# on every row of a patient that the fit uses (`rows`). A cluster coarser
# than the patient, such as a site, satisfies this; a cluster that changes
# within a patient's follow-up does not.
markov_check_path_clusters <- function(cluster, id, rows, id_var) {
  values <- cluster$values[rows]
  if (anyNA(values)) {
    mostr_abort(
      "The variable in `cluster(",
      cluster$name,
      ")` is missing on ",
      sum(is.na(values)),
      " row",
      if (sum(is.na(values)) > 1L) "s",
      " used by the model with unobserved states. Every visit needs its ",
      "cluster, because the random intercept enters the probabilities of ",
      "the unobserved visits too."
    )
  }
  patient <- as.character(id[rows])
  pairs <- unique(data.frame(patient = patient, value = as.character(values)))
  varying <- unique(pairs$patient[duplicated(pairs$patient)])
  if (length(varying)) {
    mostr_abort(
      "The variable in `cluster(",
      cluster$name,
      ")` changes within a patient (first: `",
      id_var,
      "` = ",
      varying[[1L]],
      "). With unobserved states, a patient's visits must all belong to one ",
      "cluster, because a run of unobserved visits has one likelihood term ",
      "and so one random intercept."
    )
  }
  invisible(NULL)
}

# Cluster numbers as `rmsb::blrm()` assigns them,
# `as.integer(as.factor(cluster))` over the fitted rows, here the rows `used`
# (the fully observed transitions and the stretches): clusters are numbered
# in the order of `as.factor()` (sorted values, or the factor's levels), and
# there are as many as the largest number. Returns `code` (the number of
# every element of `values`, `NA` for clusters not used), `levels` and `n`.
markov_path_cluster_codes <- function(values, used) {
  used_factor <- as.factor(values[used])
  n <- max(as.integer(used_factor))
  cluster_levels <- levels(used_factor)[seq_len(n)]
  list(
    code = match(as.character(values), cluster_levels),
    levels = cluster_levels,
    n = n
  )
}

# The observed state of each fitted row of a blrm fit. rmsb stores the
# response as a two-column matrix of positions in the state labels (the ends
# of a range), which are equal here.
markov_blrm_observed_states <- function(fit) {
  y <- fit$y
  if (is.matrix(y)) {
    y <- y[, 1L]
  }
  if (is.numeric(y)) {
    return(fit$ylevels[y])
  }
  as.character(y)
}

markov_path_stan_code <- function() {
  file <- system.file("stan", "po_paths.stan", package = "mostr")
  paste(readLines(file, warn = FALSE), collapse = "
")
}

# The C++ headers compiled with the Stan program: the likelihood shared with
# the orm fit, its Stan interface, and the definition of the program's
# `path_loglik()`, which must be included inside the model's namespace.
markov_path_stan_headers <- function() {
  files <- system.file(
    "include",
    "mostr",
    c("path_likelihood.hpp", "stan_path_loglik.hpp", "stan_path_loglik_functions.hpp"),
    package = "mostr"
  )
  if (length(files) != 3L || !all(nzchar(files))) {
    mostr_abort(
      "The C++ files of mostr's Stan program are missing from the ",
      "installation; reinstall mostr."
    )
  }
  stats::setNames(
    normalizePath(files, winslash = "/"),
    c("likelihood", "interface", "functions")
  )
}

# Everything the compiled program is made from, for the cache key: a change
# to the likelihood's C++ code must compile the program again just as a
# change to the Stan code does.
markov_path_stan_source <- function(headers = markov_path_stan_headers()) {
  c(
    markov_path_stan_code(),
    vapply(
      headers,
      function(file) paste(readLines(file, warn = FALSE), collapse = "
"),
      character(1)
    )
  )
}

# Compiled Stan programs of this session, by backend.
markov_stan_models <- new.env(parent = emptyenv())

# The compiled Stan program. mostr does not ship it compiled, which would
# make rstan and a C++ toolchain build requirements of the whole package for
# one optional engine. It is compiled on first use, which takes a minute or
# more, and kept in the user's cache directory so that later sessions load it
# in well under a second. `options(mostr.stan_cache = )` names another
# directory, and `FALSE` keeps the compiled program for this session only.
#
# The program declares `path_loglik()` without a body (`allow_undefined`);
# the C++ definition comes from mostr's headers (`markov_path_stan_headers()`),
# which are included inside the model's namespace: by rstan, which inserts
# `includes` there, and for CmdStan through a user header that opens the
# namespace itself.
markov_path_stan_model <- function(backend) {
  cached <- markov_stan_models[[backend]]
  if (!is.null(cached)) {
    return(cached)
  }
  cache <- markov_stan_cache_dir(getOption("mostr.stan_cache", NULL))
  announce <- function() message(markov_stan_compile_message(cache))
  headers <- markov_path_stan_headers()
  source <- markov_path_stan_source(headers)

  if (backend == "cmdstan") {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      mostr_abort(
        "`backend = \"cmdstan\"` needs the cmdstanr package."
      )
    }
    model <- markov_cmdstan_model(cache, source, headers, announce)
  } else {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      mostr_abort("`engine = \"blrm\"` needs the rstan package.")
    }
    model <- markov_cached_stan_model(
      cache,
      key = markov_stan_cache_key(source),
      compile = function() {
        announce()
        rstan::stan_model(
          model_code = source[[1L]],
          model_name = "mostr_po_paths",
          allow_undefined = TRUE,
          includes = paste0(
            "
#include \"",
            headers[["interface"]],
            "\"
#include \"",
            headers[["functions"]],
            "\"
"
          ),
          auto_write = FALSE
        )
      }
    )
  }
  markov_stan_models[[backend]] <- model
  model
}

# The program compiled by CmdStan. cmdstanr compiles a `.stan` file into an
# executable next to it and recompiles only when the file is newer than the
# executable, so a change to the C++ headers alone would not be picked up:
# each version of the program (Stan code, headers, CmdStan version) gets its
# own directory, named by its key, and the directories of other versions are
# deleted after a successful compile, as for rstan
# (`markov_cached_stan_model()`). The model's C++ namespace is named after
# the file, `po_paths_model_namespace`. Without a usable cache directory the
# program is compiled in a temporary directory for this session.
markov_cmdstan_model <- function(cache, source, headers, announce) {
  key <- markov_stan_cache_key(
    source,
    versions = c(cmdstan = as.character(cmdstanr::cmdstan_version()))
  )
  root <- if (cache$use) cache$dir else tempdir()
  dir <- file.path(root, paste0("po_paths_cmdstan_", key))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  stan_file <- file.path(dir, "po_paths.stan")
  header_file <- file.path(dir, "po_paths_functions.hpp")
  if (!file.exists(stan_file)) {
    writeLines(source[[1L]], stan_file)
  }
  # The header names mostr's installed headers by absolute path, but the key
  # covers their contents, not their location. Written only when absent, a
  # header left by a failed or interrupted first compile would keep pointing
  # to a library that may be gone (a `devtools::load_all()` tree, an
  # `R CMD check` library) and every later compile would fail until the
  # cache was deleted. It is rewritten whenever its text differs, and left
  # untouched otherwise so that its time stamp does not trigger a recompile.
  header_text <- c(
    "#include <stan/model/model_header.hpp>",
    paste0("#include \"", headers[["interface"]], "\""),
    "namespace po_paths_model_namespace {",
    paste0("#include \"", headers[["functions"]], "\""),
    "}"
  )
  current <- if (file.exists(header_file)) {
    tryCatch(readLines(header_file, warn = FALSE), error = function(e) NULL)
  }
  if (!identical(current, header_text)) {
    writeLines(header_text, header_file)
  }
  executable <- file.path(
    dir,
    paste0("po_paths", if (.Platform$OS.type == "windows") ".exe")
  )
  if (!file.exists(executable)) {
    announce()
  }
  # cmdstanr names the model after the file (stanc's `--name=po_paths_model`),
  # which gives the namespace the header opens.
  model <- suppressMessages(cmdstanr::cmdstan_model(
    stan_file,
    user_header = header_file
  ))
  if (cache$use) {
    outdated <- setdiff(
      list.files(cache$dir, "^po_paths_cmdstan_"),
      basename(dir)
    )
    unlink(file.path(cache$dir, outdated), recursive = TRUE)
  }
  model
}

# The cache directory from `options(mostr.stan_cache = )`: `NULL` is the
# user's cache directory, a path is that directory, anything else (`FALSE`)
# turns the cache off. `use` says whether compiled programs are written there:
# not when it is off or can't be written to.
markov_stan_cache_dir <- function(option) {
  dir <- if (is.null(option)) tools::R_user_dir("mostr", "cache") else option
  if (!is.character(dir) || length(dir) != 1L || is.na(dir)) {
    return(list(dir = NULL, use = FALSE, status = "off"))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dir) || !markov_dir_writable(dir)) {
    return(list(dir = dir, use = FALSE, status = "unwritable"))
  }
  list(dir = dir, use = TRUE, status = "ok")
}

# Whether a file can be created in `dir`. `file.access(dir, 2)` is not used:
# on Windows it reports only whether the directory exists, and in a
# read-only directory `saveRDS()` then fails after the compile.
markov_dir_writable <- function(dir) {
  probe <- tempfile("write_test_", tmpdir = dir)
  on.exit(unlink(probe), add = TRUE)
  isTRUE(suppressWarnings(file.create(probe, showWarnings = FALSE)))
}

markov_stan_compile_message <- function(cache) {
  paste0(
    "Compiling the Stan program for models with unobserved states. This ",
    "takes a minute or more and is done once",
    switch(
      cache$status,
      ok = paste0("; the result is kept in ", cache$dir),
      off = " per session",
      unwritable = paste0(
        " per session, because the cache directory ",
        cache$dir,
        " can't be written to (set `options(mostr.stan_cache = )` to another ",
        "directory)"
      )
    ),
    "."
  )
}

# The name of the cached program. It changes with the program, with R, and
# with the versions of rstan and of the packages whose code the compiled
# program links (StanHeaders, Rcpp, RcppEigen, RcppParallel; rstan's own
# version does not change when they are updated), so an outdated compiled
# program is never loaded.
markov_stan_cache_key <- function(code, versions = markov_stan_toolchain_versions()) {
  rlang::hash(list(
    code,
    versions,
    R.version$version.string,
    R.version$platform
  ))
}

markov_stan_toolchain_versions <- function() {
  packages <- c("rstan", "StanHeaders", "Rcpp", "RcppEigen", "RcppParallel")
  installed <- packages[nzchar(vapply(
    packages,
    function(package) system.file(package = package),
    character(1)
  ))]
  vapply(
    installed,
    function(package) as.character(utils::packageVersion(package)),
    character(1)
  )
}

# Load the program compiled under `key` from the cache, or compile it with
# `compile()` and save it there. The file is written under a temporary name
# in the same directory and then renamed, so another R session never reads a
# half-written file. After a successful save the files of other keys
# (earlier programs or toolchains) are deleted: CRAN allows a package to keep
# files in the user's cache directory only if it removes outdated ones, and
# each is about 1.6 MB.
markov_cached_stan_model <- function(cache, key, compile) {
  if (!cache$use) {
    return(compile())
  }
  file <- file.path(cache$dir, paste0("po_paths_", key, ".rds"))
  if (file.exists(file)) {
    model <- tryCatch(readRDS(file), error = function(e) NULL)
    if (inherits(model, "stanmodel")) {
      return(model)
    }
  }
  model <- compile()
  temporary <- tempfile("po_paths_", tmpdir = cache$dir, fileext = ".tmp")
  saved <- tryCatch(
    {
      saveRDS(model, temporary)
      file.rename(temporary, file)
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  unlink(temporary)
  if (isTRUE(saved)) {
    outdated <- setdiff(
      list.files(cache$dir, "^po_paths_.*[.]rds$"),
      basename(file)
    )
    unlink(file.path(cache$dir, outdated))
  }
  model
}
