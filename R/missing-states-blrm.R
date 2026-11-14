# Fitting when some states are not observed: `engine = "blrm"`.
#
# rmsb's Stan programs take one row per transition, so they can't sum over
# the states of unobserved visits. mostr samples from its own Stan program
# (`inst/stan/po_paths.stan`) with the same parameters and priors, and puts
# the draws into a `blrm` object so that everything downstream reads an
# ordinary rmsb fit.

mostm_fit_blrm_paths <- function(setup, ranges, path, ...) {
  args <- list(...)
  backend <- args$backend %||% getOption("rmsb.backend", "rstan")
  backend <- match.arg(backend, c("rstan", "cmdstan"))
  labels <- ranges$labels
  observed_setup <- markov_observed_transition_setup(setup, ranges)

  # The model structure comes from `rmsb::blrm()` itself, fitted to the fully
  # observed transitions: design matrix, state labels and every field the
  # prediction code reads. Its draws are replaced below, so a few iterations
  # of one chain are enough; they take about a second with rmsb's compiled
  # program, which is why the backend is rstan here whatever the user chose.
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
  if (!is.null(setup$subset_rows)) {
    structure_call$subset <- setup$subset_rows
  }
  # Ten warm-up iterations always draw sampler warnings (R-hat, effective
  # sample size); they describe draws that are thrown away.
  fit <- suppressWarnings(do.call(rmsb::blrm, structure_call))
  markov_check_observed_states(fit$ylevels, labels)
  fit_frame <- markov_rms_model_frame(
    observed_setup$formula,
    observed_setup$data,
    setup$dots,
    subset = setup$subset_rows,
    eval_env = setup$eval_env
  )
  fit_data <- markov_align_model_data(setup$data, fit_frame)

  rows <- mostm_subset_row_index(setup$subset_rows, nrow(setup$data))
  data <- setup$data[rows, , drop = FALSE]
  plan <- markov_path_plan(
    data = data,
    ranges = list(
      labels = labels,
      lo = ranges$lo[rows],
      hi = ranges$hi[rows]
    ),
    absorb = match(as_state_labels(path$absorb), labels),
    id_var = path$id_var,
    time_var = path$time_var,
    p_var = path$p_var,
    weights = NULL,
    design = function(rows) blrm_design_matrix(fit, rows)
  )

  # Stan samples coefficients of centred, orthogonalised predictors, as rmsb
  # does (`rmsb::selectedQr()`); the stretches' design rows go through the
  # same transformation, and the draws are transformed back below.
  x <- as.matrix(fit$x)
  n_predictors <- ncol(x)
  transform <- rmsb::selectedQr(x, center = TRUE)
  x_gap <- if (is.null(plan$x_gap)) {
    matrix(0, 0L, n_predictors)
  } else {
    sweep(unname(plan$x_gap), 2L, transform$xbar) %*% transform$Rinv
  }
  n_states <- length(labels)
  contrast <- fit$Contrast
  n_contrasts <- if (is.null(contrast)) 0L else nrow(contrast)
  stan_data <- list(
    N = nrow(x),
    p = n_predictors,
    k = n_states,
    X = unname(unclass(transform$X))[, , drop = FALSE],
    y = match(as.character(markov_blrm_observed_states(fit)), labels),
    G = nrow(x_gap),
    Xg = x_gap,
    S = length(plan$stretch_n_steps),
    T = length(plan$step_lo),
    first_row = as.array(plan$stretch_first_row),
    first_step = as.array(plan$stretch_first_step),
    n_steps = as.array(plan$stretch_n_steps),
    weight = as.array(plan$stretch_weight),
    lo = as.array(plan$step_lo),
    hi = as.array(plan$step_hi),
    step_rows = plan$step_rows,
    absorbing = as.array(as.integer(seq_len(n_states) %in% plan$absorb)),
    conc = args$conc %||% fit$conc %||% (1 / (0.8 + 0.35 * max(n_states, 3))),
    cn = n_contrasts,
    C = if (n_contrasts) {
      contrast %*% transform$Rinv
    } else {
      matrix(0, 0L, n_predictors)
    },
    cmus = as.array(rep_len(args$pcontrast$mu %||% 0, n_contrasts)),
    csds = as.array(rep_len(args$pcontrast$sd %||% 100, n_contrasts))
  )
  attr(stan_data$X, "scaled:center") <- NULL

  iter <- args$iter %||% 2000
  warmup <- args$warmup %||% (iter / 2)
  chains <- args$chains %||% 4
  refresh <- args$refresh %||% 0
  model <- markov_path_stan_model(backend)
  elapsed <- system.time(
    sampled <- if (backend == "rstan") {
      do.call(
        rstan::sampling,
        c(
          list(
            model,
            data = stan_data,
            pars = c("alpha", "beta"),
            iter = iter,
            warmup = warmup,
            chains = chains,
            refresh = refresh
          ),
          args$sampling.args
        )
      )
    } else {
      do.call(
        model$sample,
        c(
          list(
            data = stan_data,
            iter_sampling = iter - warmup,
            iter_warmup = warmup,
            chains = chains,
            refresh = refresh
          ),
          args$sampling.args
        )
      )
    }
  )
  if (backend == "rstan" && sampled@mode == 2L) {
    mostr_abort("The Stan sampler failed for the model with unobserved states.")
  }

  if (backend == "rstan") {
    raw <- as.matrix(sampled)
    parameter_names <- colnames(raw)
  } else {
    all_draws <- sampled$draws()
    parameter_names <- dimnames(all_draws)[[3L]]
    parameter_names <- parameter_names[grepl("^(alpha|beta)\\[", parameter_names)]
    raw <- matrix(
      all_draws[, , parameter_names],
      ncol = length(parameter_names),
      dimnames = list(NULL, parameter_names)
    )
  }
  alpha_names <- parameter_names[grepl("^alpha\\[", parameter_names)]
  beta_names <- parameter_names[grepl("^beta\\[", parameter_names)]
  betas <- raw[, beta_names, drop = FALSE] %*% t(transform$Rinv)
  # The intercepts were sampled for centred predictors.
  alphas <- raw[, alpha_names, drop = FALSE] -
    as.numeric(betas %*% transform$xbar)
  colnames(alphas) <- colnames(fit$draws)[seq_len(n_states - 1L)]
  colnames(betas) <- colnames(x)
  draws <- cbind(alphas, betas)

  fit$draws <- draws
  fit$param <- rbind(
    mean = colMeans(draws),
    median = apply(draws, 2L, stats::median)
  )
  fit$diagnostics <- if (backend == "rstan") {
    tryCatch(
      rstan::summary(
        sampled,
        pars = c(alpha_names, beta_names),
        probs = NULL
      )$summary[, c("n_eff", "Rhat")],
      error = function(e) NULL
    )
  } else {
    list(
      diagnostic_summary = sampled$diagnostic_summary(),
      summary = as.data.frame(sampled$summary(variables = c("alpha", "beta")))
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
    n_unobserved_visits = plan$n_steps_total
  )
  list(fit = fit, fit_data = fit_data, fitted_ids = patient_ids)
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
  paste(readLines(file, warn = FALSE), collapse = "\n")
}

# Compiled Stan programs of this session, by backend.
markov_stan_models <- new.env(parent = emptyenv())

# The compiled Stan program. mostr does not ship it compiled, which would
# make rstan and a C++ toolchain build requirements of the whole package for
# one optional engine. It is compiled on first use, which takes a minute or
# more, and kept in the user's cache directory so that later sessions load it
# in well under a second. `options(mostr.stan_cache = )` names another
# directory, and `FALSE` keeps the compiled program for this session only.
markov_path_stan_model <- function(backend) {
  cached <- markov_stan_models[[backend]]
  if (!is.null(cached)) {
    return(cached)
  }
  cache <- getOption("mostr.stan_cache", NULL)
  if (is.null(cache)) {
    cache <- tools::R_user_dir("mostr", "cache")
  }
  use_cache <- is.character(cache) && length(cache) == 1L && !is.na(cache)
  if (use_cache) {
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    use_cache <- dir.exists(cache)
  }
  code <- markov_path_stan_code()
  announce <- function() {
    message(
      "Compiling the Stan program for models with unobserved states. This ",
      "takes a minute or more and is done once",
      if (use_cache) paste0("; the result is kept in ", cache) else
        " per session",
      "."
    )
  }

  if (backend == "cmdstan") {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      mostr_abort(
        "`backend = \"cmdstan\"` needs the cmdstanr package."
      )
    }
    directory <- if (use_cache) cache else tempdir()
    announce()
    model <- suppressMessages(cmdstanr::cmdstan_model(
      system.file("stan", "po_paths.stan", package = "mostr"),
      dir = directory
    ))
  } else {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      mostr_abort("`engine = \"blrm\"` needs the rstan package.")
    }
    # The key changes with the program, rstan and R, so an outdated compiled
    # program is never loaded.
    key <- rlang::hash(list(
      code,
      as.character(utils::packageVersion("rstan")),
      R.version$version.string,
      R.version$platform
    ))
    file <- if (use_cache) file.path(cache, paste0("po_paths_", key, ".rds"))
    model <- if (use_cache && file.exists(file)) {
      tryCatch(readRDS(file), error = function(e) NULL)
    }
    if (!inherits(model, "stanmodel")) {
      announce()
      model <- rstan::stan_model(
        model_code = code,
        model_name = "mostr_po_paths",
        auto_write = FALSE
      )
      if (use_cache) {
        tryCatch(saveRDS(model, file), error = function(e) NULL)
      }
    }
  }
  markov_stan_models[[backend]] <- model
  model
}
