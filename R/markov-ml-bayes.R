# Bayesian fits of the blrm engine with missing states.
#
# With `missing_states` or an Ocens() response, mostm(engine = "blrm") samples
# the posterior of the rmsb::blrm() transition model with the observed-data
# likelihood. The data set-up is shared with the frequentist engines
# (markov_ml_prepare()); the Stan program in `inst/stan/markov_ml.stan` keeps
# the parameterisation and priors of rmsb's `lrmconppo.stan` and replaces its
# likelihood by the forward recursion over the rows that need it.

# Fit the blrm engine with missing states. `engine_args` are the rmsb::blrm()
# arguments given to mostm().
mostm_fit_ml_blrm <- function(setup, problem, engine_args) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop(
      "Package 'rstan' is required for `engine = \"blrm\"` with ",
      "`missing_states` or an `Ocens()` response.",
      call. = FALSE
    )
  }
  defaults <- list(
    ppo = NULL,
    cppo = NULL,
    normcppo = FALSE,
    conc = NULL,
    pcontrast = NULL,
    npcontrast = NULL,
    psigma = 1,
    rsdsd = 1,
    iter = 2000,
    chains = 4,
    refresh = 0,
    sampling.args = NULL
  )
  # rmsb::blrm() hands unknown arguments to the Stan optimizer, which this
  # sampler does not run.
  unknown <- setdiff(names(engine_args), c(names(defaults), "warmup", "rsdmean"))
  if (length(unknown)) {
    stop(
      "With `missing_states` or an `Ocens()` response, `engine = \"blrm\"` ",
      "does not use ",
      paste0("`", unknown, "`", collapse = ", "),
      ". Settings of the Stan sampler, such as `seed` or `cores`, go in ",
      "`sampling.args`.",
      call. = FALSE
    )
  }
  args <- utils::modifyList(defaults, engine_args, keep.null = TRUE)
  # Without `warmup`, half the iterations, rounded down as rstan::sampling()
  # does, so that an odd `iter` is accepted.
  args$warmup <- engine_args$warmup %||%
    (if (is.numeric(args$iter)) floor(args$iter / 2) else NA)
  args$rsdmean <- engine_args$rsdmean %||% (if (args$psigma == 1) 0 else 1)
  markov_bayes_check_sampling(args$iter, args$warmup, args$chains, args$refresh)
  if (
    !is.null(args$conc) &&
      !(is.numeric(args$conc) &&
        length(args$conc) == 1L &&
        is.finite(args$conc) &&
        args$conc > 0)
  ) {
    stop("`conc` must be a positive number.", call. = FALSE)
  }
  if (!is.null(args$sampling.args) && !is.list(args$sampling.args)) {
    stop("`sampling.args` must be a list.", call. = FALSE)
  }
  cppo <- markov_bayes_check_ppo(
    args$ppo,
    args$cppo,
    args$normcppo,
    args$npcontrast,
    problem$selected
  )
  markov_bayes_check_sigma_prior(args$psigma, args$rsdmean, args$rsdsd)

  problem <- markov_bayes_problem(problem, setup)
  random_effects <- problem$random_effects
  conc <- args$conc %||% (1 / (0.8 + 0.35 * max(problem$K, 3)))
  fit <- markov_bayes_fit_object(
    problem,
    setup$formula,
    args$ppo,
    cppo,
    args$normcppo
  )
  stan_input <- markov_bayes_stan_data(
    problem,
    fit,
    conc = conc,
    pcontrast = args$pcontrast,
    npcontrast = args$npcontrast,
    random_effects = random_effects,
    psigma = args$psigma,
    rsdmean = args$rsdmean,
    rsdsd = args$rsdsd
  )

  model <- markov_bayes_stan_model()
  sampling_call <- c(
    list(
      model,
      data = stan_input$data,
      pars = c(
        "alpha",
        "beta",
        if (stan_input$data$q > 0) "tau",
        if (random_effects) c("gamma", "sigmag")
      ),
      iter = args$iter,
      warmup = args$warmup,
      chains = args$chains,
      refresh = args$refresh
    ),
    args$sampling.args
  )
  sampling_time <- system.time(
    stanfit <- do.call(rstan::sampling, sampling_call)
  )[["elapsed"]]
  if (stanfit@mode != 0L) {
    stop("Stan sampling failed.", call. = FALSE)
  }
  draws_array <- markov_bayes_draws(stanfit, stan_input)
  draws <- matrix(
    draws_array,
    ncol = dim(draws_array)[[3L]],
    dimnames = list(NULL, dimnames(draws_array)[[3L]])
  )

  fit$draws <- draws
  fit$param <- rbind(
    mean = colMeans(draws),
    median = apply(draws, 2L, stats::median)
  )
  fit$diagnostics <- markov_bayes_diagnostics(draws_array)
  if (random_effects) {
    fit <- markov_bayes_add_random_effects(fit, stanfit, problem, problem$id_var)
  }
  fit$xbar <- stan_input$xbar
  fit$zbar <- stan_input$zbar
  fit$conc <- conc
  fit$pcontrast <- args$pcontrast
  fit$npcontrast <- args$npcontrast
  fit$Contrast <- stan_input$contrast
  fit$nContrast <- stan_input$ncontrast
  fit$psigma <- args$psigma
  fit$rsdmean <- args$rsdmean
  fit$rsdsd <- args$rsdsd
  fit$iter <- args$iter
  fit$warmup <- args$warmup
  fit$chains <- args$chains
  fit$sampling_time <- sampling_time
  fit$rstan <- stanfit
  fit$ml <- c(
    list(
      states = problem$labels,
      missing_states = problem$labels[problem$missing_index],
      p_var = problem$p_var,
      id_var = problem$id_var,
      start_coefficients = problem$start_coefficients
    ),
    markov_ml_missing_counts(problem)
  )
  list(fit = fit)
}

# Add the proportional-odds design of the rms fit to the transitions with
# both states observed: the design rows `X` (without intercept) of the
# possible previous states, the coefficient names, and the start fit, whose
# design metadata the returned blrm object uses for prediction.
markov_bayes_problem <- function(problem, setup) {
  start <- markov_ml_start(problem, "blrm", setup, list())
  problem$start_fit <- start$fit
  problem$lp <- start$lp
  problem$X <- start$lp$X[, -1L, drop = FALSE]
  problem$coef_names <- names(start$theta)
  problem$start_coefficients <- start$theta
  problem
}

markov_bayes_check_sampling <- function(iter, warmup, chains, refresh) {
  count <- function(x, min) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= min && x == round(x)
  }
  if (!count(iter, 2)) {
    stop("`iter` must be an integer of at least 2.", call. = FALSE)
  }
  if (!count(warmup, 1) || warmup >= iter) {
    stop("`warmup` must be a positive integer below `iter`.", call. = FALSE)
  }
  if (!count(chains, 1)) {
    stop("`chains` must be a positive integer.", call. = FALSE)
  }
  if (!count(refresh, 0)) {
    stop("`refresh` must be a non-negative integer.", call. = FALSE)
  }
  invisible(NULL)
}

markov_bayes_check_sigma_prior <- function(psigma, rsdmean, rsdsd) {
  if (!identical(as.numeric(psigma), 1) && !identical(as.numeric(psigma), 2)) {
    stop("`psigma` must be 1 (t prior) or 2 (exponential prior).", call. = FALSE)
  }
  non_negative <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  }
  if (!non_negative(rsdmean) || !non_negative(rsdsd) || rsdsd == 0) {
    stop(
      "`rsdmean` must be a non-negative number and `rsdsd` a positive number.",
      call. = FALSE
    )
  }
  if (psigma == 2 && rsdmean == 0) {
    stop("`rsdmean` must be positive when `psigma = 2`.", call. = FALSE)
  }
  invisible(NULL)
}

# Posterior draws of the random intercepts, stored as `blrm()` fits with
# `cluster()` store them, plus a draws x patients matrix named by patient ID
# that the SOP functions use with `include_re = TRUE`.
markov_bayes_add_random_effects <- function(fit, stanfit, problem, id_var) {
  gamma <- as.matrix(stanfit, pars = "gamma")
  sigmag <- as.matrix(stanfit, pars = "sigmag")[, 1L]
  colnames(gamma) <- problem$patient_ids
  fit$gamma_draws <- gamma
  fit$gammas <- apply(gamma, 2L, stats::median)
  fit$omega <- cbind(sigmag = sigmag)
  fit$clusterInfo <- list(
    cluster = problem$patient_ids[problem$row_patient],
    n = problem$N,
    name = id_var
  )
  sims <- rstan::extract(stanfit, pars = "sigmag", permuted = FALSE)
  sims <- matrix(sims, ncol = dim(sims)[[2L]])
  fit$diagnostics <- rbind(
    fit$diagnostics,
    sigmag = c(
      Rhat = rstan::Rhat(sims),
      ess_bulk = rstan::ess_bulk(sims),
      ess_tail = rstan::ess_tail(sims)
    )
  )
  fit
}

# Validate the partial proportional-odds arguments and return `cppo` as a
# function (or NULL).
markov_bayes_check_ppo <- function(ppo, cppo, normcppo, npcontrast, data) {
  if (!is.logical(normcppo) || length(normcppo) != 1L || is.na(normcppo)) {
    stop("`normcppo` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(ppo)) {
    if (!is.null(cppo) || !is.null(npcontrast) || isTRUE(normcppo)) {
      stop(
        "`cppo`, `normcppo`, and `npcontrast` require `ppo`.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (!inherits(ppo, "formula") || length(ppo) != 2L) {
    stop("`ppo` must be a one-sided formula such as `~ time`.", call. = FALSE)
  }
  if (is.null(cppo)) {
    stop(
      "Only constrained partial proportional odds is supported: give `cppo`, ",
      "a function scoring the states, together with `ppo`.",
      call. = FALSE
    )
  }
  cppo <- resolve_blrm_cppo(cppo, envir = environment(ppo) %||% globalenv())
  absent <- setdiff(all.vars(ppo), names(data))
  if (length(absent)) {
    stop(
      "Variable(s) ",
      paste0("`", absent, "`", collapse = ", "),
      " in `ppo` are not columns of `data`.",
      call. = FALSE
    )
  }
  cppo
}

# A `blrm`-like fit object without draws. The design metadata of `formula`
# come from the starting rms fit on transitions with both states observed; the
# `ppo` design is built from the same transitions with rms::Design(), as
# rmsb::blrm() does, so prediction rebuilds the design rows used in fitting.
# The response column is the internal `.mostr_state`.
markov_bayes_fit_object <- function(problem, formula, ppo, cppo, normcppo) {
  start_fit <- problem$start_fit
  K <- problem$K
  comp <- problem$complete_part
  observed_state <- rep(comp$agg_state, comp$agg_count)
  freq <- tabulate(observed_state, nbins = K)
  names(freq) <- problem$labels
  p <- ncol(problem$X)

  fit <- list(
    call = NULL,
    fitter = "markov_ml",
    link = "logistic",
    yname = paste(deparse(formula[[2L]]), collapse = ""),
    ylevels = problem$labels,
    freq = freq,
    pppo = 0L,
    cppo = NULL,
    Design = start_fit$Design,
    zDesign = NULL,
    scale.pred = c("log odds", "Odds Ratio"),
    terms = start_fit$terms,
    assign = rms::DesignAssign(start_fit$Design, K - 1L, start_fit$terms),
    zassign = NULL,
    na.action = NULL,
    fail = FALSE,
    non.slopes = K - 1L,
    interceptRef = start_fit$interceptRef %||% 1L,
    sformula = start_fit$sformula,
    zsformula = NULL,
    N = problem$n,
    Ncens = c(left = 0, right = 0, interval = 0),
    p = p,
    tauInfo = NULL,
    iprior = 0,
    backend = "rstan"
  )

  if (!is.null(ppo)) {
    rows <- problem$design_data[rep(comp$agg_design, comp$agg_count), ,
      drop = FALSE]
    rows[[markov_ml_response_column()]] <- markov_ml_state_values(
      observed_state,
      problem$scale,
      problem$response_type
    )
    model_data <- rms::modelData(
      rows,
      problem$working_formula,
      ppo,
      callenv = environment(ppo) %||% globalenv()
    )
    zframe <- rms::Design(attr(model_data, "data2"), formula = ppo)
    zattr <- attributes(zframe)
    zterms <- zattr$terms
    attr(zterms, "formula") <- ppo
    if (length(zattr$Design$colnames) < 1L) {
      stop("`ppo` must contain at least one term.", call. = FALSE)
    }
    if (isTRUE(normcppo)) {
      scores <- cppo(problem$labels[observed_state])
      cppo <- markov_bayes_normalized_cppo(
        cppo,
        center = mean(scores),
        scale = stats::sd(scores)
      )
    }
    fit$pppo <- length(zattr$Design$colnames)
    fit$cppo <- cppo
    fit$zDesign <- zattr$Design
    fit$zsformula <- zattr$sformula
    fit$zassign <- rms::DesignAssign(zattr$Design, 0L, zterms)
    fit$tauInfo <- data.frame(
      name = paste(zattr$Design$colnames, "x f(y)"),
      x = zattr$Design$colnames
    )
  }
  # rms prediction treats a fit as Bayesian only when it has draws.
  fit$param <- markov_bayes_param_placeholder(problem$coef_names, fit$tauInfo)
  fit$draws <- fit$param[1L, , drop = FALSE]
  class(fit) <- c("mostm_ml", "blrm", "rmsb", "rms")
  fit
}

markov_bayes_normalized_cppo <- function(fun, center, scale) {
  force(fun)
  force(center)
  if (!is.finite(scale) || scale <= 0) {
    stop(
      "`normcppo = TRUE` needs `cppo` scores that vary across the observed ",
      "states.",
      call. = FALSE
    )
  }
  function(y) (fun(y) - center) / scale
}

markov_bayes_param_placeholder <- function(coef_names, tau_info) {
  names <- c(coef_names, tau_info$name)
  matrix(
    0,
    nrow = 2L,
    ncol = length(names),
    dimnames = list(c("mean", "median"), names)
  )
}

# Centre the rows of `X` at their weighted mean and return the inverse R
# factor of the QR decomposition of the weighted, centered rows. With integer
# weights this equals rmsb's selectedQr() on the rows repeated `weight` times.
markov_bayes_qr <- function(X, weight) {
  p <- ncol(X)
  xbar <- colSums(X * weight) / sum(weight)
  qr_obs <- qr(sweep(X, 2L, xbar) * sqrt(weight))
  if (qr_obs$rank < p || !identical(qr_obs$pivot, seq_len(p))) {
    stop(
      "The design of the transitions with both states observed is not of ",
      "full rank, so some coefficients are not identified.",
      call. = FALSE
    )
  }
  R <- qr.R(qr_obs)
  R <- sweep(R, 1L, sign(diag(R)), "*")
  list(xbar = xbar, Rinv = backsolve(R, diag(p)))
}

# Stan data. Predictors are centered at their averages over the transitions
# with both states observed and multiplied by the inverse R factor of the QR
# decomposition of those centered rows, as rmsb::blrm() does for the rows it
# fits. The Dirichlet prior then refers to a patient with average predictors,
# and with nothing missing the parameters match blrm()'s one to one. Only the
# design rows that the likelihood uses are passed.
markov_bayes_stan_data <- function(
  problem,
  fit,
  conc,
  pcontrast,
  npcontrast,
  random_effects = FALSE,
  psigma = 1,
  rsdmean = 0,
  rsdsd = 1
) {
  comp <- problem$complete_part
  chain <- problem$chain
  K <- problem$K
  X <- problem$X
  p <- ncol(X)

  X_check <- blrm_design_matrix(fit, problem$design_data)
  if (
    !identical(dim(X_check), dim(X)) ||
      max(abs(X_check - X)) > 1e-10 * max(1, max(abs(X)))
  ) {
    stop(
      "Internal error: prediction would build different design rows from ",
      "the ones used for fitting.",
      call. = FALSE
    )
  }

  q <- fit$pppo
  if (q > 0) {
    Z <- blrm_design_matrix(fit, problem$design_data, second = TRUE)
    if (ncol(Z) != q || any(!is.finite(Z))) {
      stop("Internal error: the partial proportional-odds design is invalid.")
    }
  }

  # Runs of rows that need the forward recursion, and for each row the
  # previous states that can have probability.
  if (chain$n) {
    if (any(diff(chain$segment) < 0L)) {
      stop("Internal error: chain rows are not grouped by run.")
    }
    seg_start <- which(!duplicated(chain$segment))
    # Rows that allow every state add nothing at the end of a run (their
    # probabilities add up to 1), so each run stops at its last row that
    # excludes a state, and runs without such a row are left out.
    informative <- which(rowSums(chain$C > 0) < K)
    seg_end <- tapply(informative, chain$segment[informative], max)
    keep <- as.integer(names(seg_end))
    seg_start <- seg_start[keep]
    seg_end <- as.integer(seg_end)
    seg_baseline <- chain$baseline[keep]
  } else {
    seg_start <- seg_end <- seg_baseline <- integer()
  }
  recursion <- markov_bayes_recursion_index(problem, seg_start, seg_end,
    seg_baseline)
  # The recursion sums over states, so the cumulative probabilities of every
  # possible previous state on every recursion row must not cross. Whether
  # they cross depends only on the partial proportional-odds design row, so
  # one row per distinct `ppo` design row is checked (one row without `ppo`).
  check_rows <- unique(as.vector(chain$DI))
  check_rows <- if (q > 0) {
    check_rows[!duplicated(Z[check_rows, , drop = FALSE])]
  } else {
    utils::head(check_rows, 1L)
  }

  weight <- comp$agg_count
  qr_x <- markov_bayes_qr(X[comp$agg_design, , drop = FALSE], weight)
  used <- sort(unique(c(comp$agg_design, recursion$cell_design, check_rows)))
  X_used <- sweep(X[used, , drop = FALSE], 2L, qr_x$xbar) %*% qr_x$Rinv

  if (q > 0) {
    qr_z <- markov_bayes_qr(Z[comp$agg_design, , drop = FALSE], weight)
    Z_used <- sweep(Z[used, , drop = FALSE], 2L, qr_z$xbar) %*% qr_z$Rinv
    pposcore <- as.numeric(fit$cppo(problem$labels))
    if (length(pposcore) != K || any(!is.finite(pposcore))) {
      stop(
        "`cppo` must return one finite score for each of the ",
        K,
        " states.",
        call. = FALSE
      )
    }
  } else {
    qr_z <- list(xbar = NULL, Rinv = matrix(0, 0L, 0L))
    Z_used <- matrix(0, nrow = length(used), ncol = 0L)
    pposcore <- rep(0, K)
  }

  # With random intercepts each patient has its own transition probabilities,
  # so fully observed transitions are aggregated within patients.
  if (random_effects) {
    key <- paste(comp$design, comp$state, comp$patient)
    first <- !duplicated(key)
    obs <- list(
      design = comp$design[first],
      state = comp$state[first],
      count = tabulate(match(key, key[first])),
      cluster = comp$patient[first]
    )
  } else {
    obs <- list(
      design = comp$agg_design,
      state = comp$agg_state,
      count = comp$agg_count,
      cluster = integer()
    )
  }

  # Design rows used by the cells, with the range of states that the cells
  # of each need (for probabilities shared by rows without random effects).
  cd_rows <- sort(unique(recursion$cell_design))
  cell_cd <- match(recursion$cell_design, cd_rows)
  cd_factor <- factor(cell_cd, levels = seq_along(cd_rows))
  cd_lo <- as.integer(tapply(recursion$row_lo[recursion$cell_row], cd_factor,
    min))
  cd_hi <- as.integer(tapply(recursion$row_hi[recursion$cell_row], cd_factor,
    max))

  contrast <- markov_bayes_contrast_prior(
    pcontrast,
    problem$start_fit,
    qr_x$Rinv,
    "pcontrast"
  )
  ncontrast <- if (q > 0) {
    markov_bayes_contrast_prior(
      npcontrast,
      markov_bayes_ppo_contrast_fit(fit, problem$K),
      qr_z$Rinv,
      "npcontrast"
    )
  } else {
    markov_bayes_contrast_prior(NULL, NULL, qr_z$Rinv, "npcontrast")
  }

  data <- list(
    k = K,
    p = p,
    q = q,
    D = length(used),
    X = unname(X_used),
    Z = unname(Z_used),
    pposcore = as.array(pposcore),
    Nc = if (random_effects) problem$N else 0L,
    n_obs = length(obs$design),
    obs_design = as.array(match(obs$design, used)),
    obs_state = as.array(as.integer(obs$state)),
    obs_count = as.array(as.numeric(obs$count)),
    obs_cluster = as.array(as.integer(obs$cluster)),
    n_seg = length(seg_start),
    seg_baseline = as.array(as.integer(seg_baseline)),
    seg_start = as.array(as.integer(seg_start)),
    seg_end = as.array(as.integer(seg_end)),
    n_chain = chain$n,
    chain_mask = unname(chain$C),
    chain_cluster = as.array(
      if (random_effects) as.integer(chain$patient) else integer()
    ),
    row_lo = as.array(recursion$row_lo),
    row_hi = as.array(recursion$row_hi),
    row_gappy = as.array(recursion$row_gappy),
    row_full = as.array(recursion$row_full),
    cell_start = as.array(recursion$cell_start),
    cell_end = as.array(recursion$cell_end),
    n_cell = length(cell_cd),
    cell_origin = as.array(recursion$cell_origin),
    cell_cd = as.array(cell_cd),
    n_cd = length(cd_rows),
    cd_rows = as.array(match(cd_rows, used)),
    cd_lo = as.array(cd_lo),
    cd_hi = as.array(cd_hi),
    n_check = length(check_rows),
    check_rows = as.array(match(check_rows, used)),
    conc = conc,
    cn = nrow(contrast$C),
    C = contrast$C,
    cmus = as.array(contrast$mu),
    csds = as.array(contrast$sd),
    cn2 = nrow(ncontrast$C),
    C2 = ncontrast$C,
    cmus2 = as.array(ncontrast$mu),
    csds2 = as.array(ncontrast$sd),
    psigma = as.integer(psigma),
    rsdmean = if (random_effects) as.array(rsdmean) else numeric(),
    rsdsd = if (random_effects && psigma == 1) as.array(rsdsd) else numeric()
  )
  list(
    data = data,
    design_rows = used,
    xbar = qr_x$xbar,
    Rinv = qr_x$Rinv,
    zbar = qr_z$xbar,
    Rinv_z = qr_z$Rinv,
    pposcore = pposcore,
    coef_names = problem$coef_names,
    tau_names = fit$tauInfo$name,
    contrast = contrast$contrast,
    ncontrast = ncontrast$contrast
  )
}

# Index lists for the forward recursion over the runs `seg_start[g]` to
# `seg_end[g]` of chain rows, starting from state `seg_baseline[g]`, derived
# from the masks before sampling:
# * `row_lo`, `row_hi`: the lowest and highest state allowed on each row; the
#   recursion computes probabilities of these states only. `row_gappy` marks
#   rows that also exclude states in between, so the mask is still applied.
#   `row_full` marks rows that allow every state, whose summed probability is
#   1.
# * cells: for each row, the previous states that can have probability (the
#   run's baseline on its first row; after that, the states the previous row
#   allows, as every state the model can leave can move to every state), with
#   the design row of the transition from each. Cells of
#   row t are `cell_start[t]` to `cell_end[t]`.
# Stops, naming the patients, when a run has no possible sequence of states.
markov_bayes_recursion_index <- function(
  problem,
  seg_start,
  seg_end,
  seg_baseline
) {
  chain <- problem$chain
  K <- problem$K
  allowed <- chain$C > 0
  is_origin <- seq_len(K) %in% problem$origins
  from <- vector("list", chain$n)
  impossible <- integer()
  for (g in seq_along(seg_start)) {
    reach <- seq_len(K) == seg_baseline[[g]]
    for (t in seg_start[[g]]:seg_end[[g]]) {
      from[[t]] <- which(reach & is_origin)
      reach <- allowed[t, ] & length(from[[t]]) > 0L
      # No state can be reached on this row, such as a state other than death
      # after a recorded death, so no coefficients fit this run.
      if (!any(reach)) {
        impossible <- c(impossible, t)
        break
      }
    }
  }
  if (length(impossible) > 0L) {
    markov_ml_stop_impossible(chain$ids[impossible])
  }
  n_cells <- lengths(from)
  cell_row <- rep(seq_len(chain$n), n_cells)
  cell_origin <- as.integer(unlist(from))
  cell_design <- chain$DI[cbind(cell_row, match(cell_origin, problem$origins))]
  row_lo <- max.col(allowed * 1, ties.method = "first")
  row_hi <- max.col(allowed * 1, ties.method = "last")
  row_gappy <- rowSums(allowed) < row_hi - row_lo + 1L
  cell_start <- cumsum(c(1L, n_cells))[seq_len(chain$n)]
  list(
    row_lo = as.integer(row_lo),
    row_hi = as.integer(row_hi),
    row_gappy = as.integer(row_gappy),
    row_full = as.integer(row_lo == 1L & row_hi == K & !row_gappy),
    cell_start = as.integer(cell_start),
    cell_end = as.integer(cell_start + n_cells - 1L),
    cell_row = cell_row,
    cell_origin = cell_origin,
    cell_design = as.integer(cell_design)
  )
}

# An object whose design is the `ppo` design, so rms::Xcontrast() builds
# contrast rows for the partial proportional-odds coefficients (as
# rmsb::blrm() does for `npcontrast`).
markov_bayes_ppo_contrast_fit <- function(fit, K) {
  out <- fit
  out$Design <- fit$zDesign
  out$assign <- fit$zassign
  out$sformula <- fit$zsformula
  out$pppo <- 0L
  out$cppo <- NULL
  out$zDesign <- NULL
  out$param <- matrix(
    0,
    nrow = 2L,
    ncol = K - 1L + fit$pppo,
    dimnames = list(c("mean", "median"), NULL)
  )
  out$draws <- out$param[1L, , drop = FALSE]
  class(out) <- c("blrm", "rmsb", "rms")
  out
}

# Normal priors for contrasts of the coefficients, following rmsb::blrm():
# settings are turned into design-row differences with rms::Xcontrast() and
# mapped to the QR-transformed coefficients.
markov_bayes_contrast_prior <- function(spec, fit, Rinv, arg) {
  p <- nrow(Rinv)
  if (is.null(spec)) {
    return(list(
      C = matrix(0, nrow = 0L, ncol = p),
      mu = numeric(),
      sd = numeric(),
      contrast = NULL
    ))
  }
  if (!is.list(spec)) {
    stop("`", arg, "` must be a list, as in `rmsb::blrm()`.", call. = FALSE)
  }
  weights <- spec$weights %||% "equal"
  expand <- spec$expand %||% TRUE
  mu <- spec$mu %||% 0
  sd <- spec$sd %||% 100
  contrast <- spec$contrast
  if (!length(contrast) || !is.expression(contrast)) {
    stop(
      "`",
      arg,
      "` must contain `contrast`, an `expression()` of the named settings.",
      call. = FALSE
    )
  }
  settings <- spec
  settings[c("weights", "expand", "mu", "sd", "contrast", "ycut")] <- NULL
  if (!length(settings) || is.null(names(settings)) ||
    any(names(settings) == "")) {
    stop(
      "Predictor settings in `",
      arg,
      "` must be named lists.",
      call. = FALSE
    )
  }
  rows <- lapply(settings, function(a) {
    rms::Xcontrast(fit, a = a, weights = weights, expand = expand,
      Zmatrix = FALSE)
  })
  C <- do.call(rbind, lapply(contrast, eval, rows))
  cn <- nrow(C)
  if (length(mu) == 1L) mu <- rep(mu, cn)
  if (length(sd) == 1L) sd <- rep(sd, cn)
  if (length(mu) != cn || length(sd) != cn) {
    stop(
      "`",
      arg,
      "` defines ",
      cn,
      " contrast(s) but gives ",
      length(mu),
      " mean(s) and ",
      length(sd),
      " standard deviation(s).",
      call. = FALSE
    )
  }
  if (any(!is.finite(sd)) || any(sd <= 0)) {
    stop("`", arg, "$sd` must be positive.", call. = FALSE)
  }
  if (ncol(C) != p) {
    stop("Internal error: contrast rows do not match the design.")
  }
  list(
    C = unname(C %*% Rinv),
    mu = as.numeric(mu),
    sd = as.numeric(sd),
    contrast = C
  )
}

# The compiled Stan program, cached for the session and, unless
# `options(mostr.stan_cache_dir = FALSE)`, on disk.
markov_bayes_cache <- new.env(parent = emptyenv())

markov_bayes_stan_model <- function() {
  source_file <- system.file("stan", "markov_ml.stan", package = "mostr")
  if (!nzchar(source_file)) {
    stop("Internal error: the Stan program for missing states is missing.")
  }
  hash <- unname(tools::md5sum(source_file))
  key <- paste(hash, utils::packageVersion("rstan"))
  if (identical(markov_bayes_cache$key, key)) {
    return(markov_bayes_cache$model)
  }

  cache_dir <- getOption("mostr.stan_cache_dir")
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("mostr", "cache")
  }
  model <- NULL
  if (!isFALSE(cache_dir)) {
    if (!is.character(cache_dir) || length(cache_dir) != 1L ||
      is.na(cache_dir)) {
      stop(
        "`options(mostr.stan_cache_dir = )` must be a directory path or FALSE.",
        call. = FALSE
      )
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    target <- file.path(
      cache_dir,
      paste0("markov_ml_", substr(hash, 1L, 12L), ".stan")
    )
    if (!file.exists(target)) {
      file.copy(source_file, target)
    }
    if (file.exists(target)) {
      if (!file.exists(sub("\\.stan$", ".rds", target))) {
        message(
          "Compiling the Stan program for missing states with ",
          "`engine = \"blrm\"`. This is done once and can take a few minutes."
        )
      }
      model <- rstan::stan_model(
        file = target,
        model_name = "markov_ml",
        auto_write = TRUE
      )
    }
  }
  if (is.null(model)) {
    message(
      "Compiling the Stan program for missing states with ",
      "`engine = \"blrm\"`. This can take a few minutes."
    )
    model <- rstan::stan_model(
      model_code = paste(readLines(source_file), collapse = "\n"),
      model_name = "markov_ml"
    )
  }
  markov_bayes_cache$key <- key
  markov_bayes_cache$model <- model
  model
}

# Posterior draws as an iterations x chains x coefficients array on the
# original scale: undo the QR transformations and move the intercepts from the
# average patient back to predictors at zero, as rmsb::blrm() does.
markov_bayes_draws <- function(stanfit, stan_input) {
  pars <- c("alpha", "beta", if (stan_input$data$q > 0) "tau")
  arr <- rstan::extract(stanfit, pars = pars, permuted = FALSE,
    inc_warmup = FALSE)
  out <- markov_bayes_draws_from_array(arr, stan_input)
  array(
    out,
    dim = c(dim(arr)[[1L]], dim(arr)[[2L]], ncol(out)),
    dimnames = list(NULL, NULL, colnames(out))
  )
}

# Transform an iterations x chains x parameters array of Stan draws to a
# matrix with one row per draw on the original scale.
markov_bayes_draws_from_array <- function(arr, stan_input) {
  names <- dimnames(arr)[[3L]]
  dims <- dim(arr)
  flat <- matrix(arr, ncol = dims[[3L]])
  betas <- flat[, grep("^beta\\[", names), drop = FALSE] %*%
    t(stan_input$Rinv)
  alphas <- flat[, grep("^alpha\\[", names), drop = FALSE] -
    drop(betas %*% stan_input$xbar)
  taus <- NULL
  if (stan_input$data$q > 0) {
    taus <- flat[, grep("^tau\\[", names), drop = FALSE] %*%
      t(stan_input$Rinv_z)
    alphas <- alphas -
      outer(drop(taus %*% stan_input$zbar), stan_input$pposcore[-1L])
    colnames(taus) <- stan_input$tau_names
  }
  out <- cbind(alphas, betas, taus)
  colnames(out) <- c(stan_input$coef_names, stan_input$tau_names)
  out
}

# Rhat and bulk and tail effective sample sizes of each coefficient on the
# original scale.
markov_bayes_diagnostics <- function(draws_array) {
  out <- t(vapply(
    seq_len(dim(draws_array)[[3L]]),
    function(j) {
      sims <- matrix(draws_array[, , j], ncol = dim(draws_array)[[2L]])
      c(
        Rhat = rstan::Rhat(sims),
        ess_bulk = rstan::ess_bulk(sims),
        ess_tail = rstan::ess_tail(sims)
      )
    },
    numeric(3)
  ))
  rownames(out) <- dimnames(draws_array)[[3L]]
  out
}
