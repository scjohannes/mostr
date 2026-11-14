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
  layout <- mostm_path_layout(setup, ranges, path)
  # No stretch: rmsb's own fit to the observed transitions is the fit (see
  # `mostm_fit_orm_paths()`), and nothing needs to be compiled.
  if (!layout$n_steps_total) {
    return(mostm_fit_blrm(markov_observed_transition_setup(setup, ranges), ...))
  }
  set_up <- markov_path_setup(setup, ranges, path)
  mostm_check_setup_states(ranges, set_up, layout, path$p_var)
  observed_setup <- set_up$setup

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

  plan <- markov_path_plan_design(
    layout,
    function(rows) blrm_design_matrix(fit, rows)
  )

  # Stan samples coefficients of centred, orthogonalised predictors, as rmsb
  # does (`rmsb::selectedQr()`); the stretches' design rows go through the
  # same transformation, and the draws are transformed back below. The
  # transformation is computed from all set-up rows, which have a full-rank
  # design whenever the set-up fit worked; the fully observed transitions
  # alone need not (a previous state that is always followed by an
  # unobserved visit has a column of zeros there). The Dirichlet prior on
  # the intercepts is therefore placed at the mean predictors of the set-up
  # rows.
  x_setup <- as.matrix(fit$x)
  x <- x_setup[keep, , drop = FALSE]
  n_predictors <- ncol(x)
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
  n_states <- length(labels)
  contrast <- fit$Contrast
  n_contrasts <- if (is.null(contrast)) 0L else nrow(contrast)
  # Observed transitions with the same predictors and state are counted.
  observed_row <- markov_distinct_rows(x)
  n_observed_rows <- if (length(observed_row)) max(observed_row) else 0L
  observed_state <- match(
    as.character(markov_blrm_observed_states(fit))[keep],
    labels
  )
  cell_count <- table(
    factor(observed_row, levels = seq_len(n_observed_rows)),
    factor(observed_state, levels = seq_len(n_states))
  )
  cells <- which(cell_count > 0)
  x_observed <- unname(unclass(transform$X))[keep, , drop = FALSE][
    !duplicated(observed_row),
    ,
    drop = FALSE
  ] * column_scale
  # The Stan program reads matrices as arrays stored by row.
  by_row <- function(m) as.array(as.numeric(t(m)))
  stan_data <- list(
    p = n_predictors,
    k = n_states,
    n_obs = length(cells),
    x_obs = by_row(x_observed[(cells - 1L) %% n_observed_rows + 1L, , drop = FALSE]),
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
    conc = args$conc %||% fit$conc %||% (1 / (0.8 + 0.35 * max(n_states, 3))),
    cn = n_contrasts,
    C = if (n_contrasts) {
      contrast %*% rinv
    } else {
      matrix(0, 0L, n_predictors)
    },
    cmus = as.array(rep_len(args$pcontrast$mu %||% 0, n_contrasts)),
    csds = as.array(rep_len(args$pcontrast$sd %||% 100, n_contrasts))
  )

  iter <- args$iter %||% 2000
  warmup <- args$warmup %||% (iter / 2)
  chains <- args$chains %||% 4
  refresh <- args$refresh %||% 0
  model <- markov_path_stan_model(backend)
  # The number of chains run in parallel is left to the backend, as
  # `rmsb::blrm()` leaves it: `rstan::sampling()` uses
  # `getOption("mc.cores", 1)` cores and cmdstanr's `$sample()` as many
  # parallel chains, and `sampling.args` (`cores`, `parallel_chains`)
  # overrides either.
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
    n_unobserved_visits = plan$n_steps_total,
    # The sampler's parametrisation: coefficients beta_s of the centred,
    # transformed predictors, beta = rinv %*% beta_s, intercepts for
    # predictors centred at `fit$xbar`. It comes from the set-up rows, which
    # are not stored, so it is kept for checking the sampled density.
    rinv = rinv
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
  if (!file.exists(header_file)) {
    writeLines(
      c(
        "#include <stan/model/model_header.hpp>",
        paste0("#include \"", headers[["interface"]], "\""),
        "namespace po_paths_model_namespace {",
        paste0("#include \"", headers[["functions"]], "\""),
        "}"
      ),
      header_file
    )
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
