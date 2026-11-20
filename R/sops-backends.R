# SOP model adapters and time helpers.

# Only `rmsb::blrm()` fits have random effects, so `include_re = TRUE` on any
# other fit would otherwise be ignored without notice. Every public function
# with an `include_re` argument calls this before doing any work.
stop_include_re_without_blrm <- function(model, include_re) {
  if (!isTRUE(include_re)) {
    return(invisible(NULL))
  }
  fit <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  if (inherits(fit, "blrm")) {
    return(invisible(NULL))
  }
  engine <- markov_model_engine(model) %||% class(fit)[1]
  mostr_abort(
    "`include_re = TRUE` adds each patient's random effect to the ",
    "predictions, and only models fitted with `rmsb::blrm()` ",
    "(`mostm(engine = \"blrm\")`) have random effects. This model was fitted ",
    "with engine \"",
    engine,
    "\"; drop `include_re` or set it to FALSE."
  )
}

#' Validate Model for Markov Simulation
#'
#' Internal helper that checks whether a model is compatible with Markov
#' simulation. For vglm models, verifies that the family is cumulative with
#' reverse = TRUE.
#'
#' @param object A fitted model object
#' @return NULL (invisibly) if valid, otherwise throws an error
#' @keywords internal
validate_markov_model <- function(object) {
  # Unwrap robcov_vglm if needed
  model_chk <- if (inherits(object, "robcov_vglm")) {
    object$vglm_fit
  } else {
    object
  }

  expected_engine <- if (inherits(model_chk, "blrm")) {
    "blrm"
  } else if (inherits(model_chk, "vglm")) {
    "vglm"
  } else if (inherits(model_chk, "orm")) {
    "orm"
  } else {
    NULL
  }
  if (
    !is.null(expected_engine) &&
      !identical(markov_model_engine(object), expected_engine)
  ) {
    mostr_abort(
      "Model-based Markov workflows require fits created by `mostm()`. ",
      "Refit this model with `mostm(engine = \"",
      expected_engine,
      "\")`."
    )
  }

  if (inherits(model_chk, "blrm")) {
    return(invisible(NULL))
  }

  # Check if it's a vglm model
  if (inherits(model_chk, "vglm")) {
    # Extract family object
    fam <- model_chk@family

    # Check if it's cumulative family
    if (!grepl("cumulative", fam@vfamily[1], ignore.case = TRUE)) {
      mostr_abort(
        "For vglm models, only the cumulative family is supported.\n",
        "Current family: ",
        fam@vfamily[1],
        "\n",
        "Please refit your model with: family = cumulative(reverse = TRUE, ...)"
      )
    }

    link <- unique(as.character(model_chk@misc$link))
    if (length(link) != 1L || !identical(tolower(link), "logitlink")) {
      stop_unsupported_markov_link(link)
    }

    # Read the reverse coding from the fitted family object. Its `infos`
    # function reports the value used at fit time, however the family was
    # supplied (inline call, local variable, or family object).
    reverse_param <- vglm_family_reverse(fam)

    # Reject fits whose reverse coding cannot be verified.
    if (is.null(reverse_param)) {
      mostr_abort(
        "Cannot determine 'reverse' parameter from the fitted vglm family.\n",
        "For safety, please refit with: family = cumulative(reverse = TRUE, ...)"
      )
    }

    # Check that reverse = TRUE
    if (!isTRUE(reverse_param)) {
      mostr_abort(
        "For vglm models, the cumulative family must use reverse = TRUE.\n",
        "Current setting: reverse = ",
        reverse_param,
        "\n",
        "Please refit your model with: family = cumulative(reverse = TRUE, ...)"
      )
    }
  }

  if (inherits(model_chk, "orm")) {
    if (orm_has_random_effects(model_chk)) {
      stop_unsupported_orm_random_effects()
    }
    family <- model_chk$family
    if (
      is.null(family) ||
        length(family) != 1L ||
        !identical(tolower(as.character(family)), "logistic")
    ) {
      stop_unsupported_markov_link(family)
    }
  }

  invisible(NULL)
}

# Return the `reverse` setting of a fitted VGAM family, or NULL when the
# family does not report it as one non-missing logical value.
vglm_family_reverse <- function(family) {
  infos <- tryCatch(methods::slot(family, "infos"), error = function(e) NULL)
  if (!is.function(infos)) {
    return(NULL)
  }
  reverse <- tryCatch(infos()$reverse, error = function(e) NULL)
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    return(NULL)
  }
  reverse
}

stop_unsupported_markov_link <- function(link = NULL) {
  message <- paste0(
    "Only cumulative-logit models are supported; ",
    "refit the model with a logit link."
  )
  mostr_abort(
    message,
    class = "mostr_unsupported_link",
    data = list(
      link = if (length(link)) as.character(link) else NA_character_
    )
  )
}

# The vglm design matrix of `newdata`, one row per row of `newdata`.
# stats::model.matrix() given a data frame builds its model frame with the
# default `na.action`, stats::na.omit(), which drops rows with a missing value
# without a word; the recursion then matched the remaining rows' probabilities
# to the wrong profiles. The frame is built with stats::na.pass() instead, and
# a missing value in the design is refused. The public routes refuse missing
# predictor values earlier, naming the column (`stop_missing_predictors()`);
# this keeps any other route from returning shifted probabilities.
vglm_prediction_design <- function(
  model,
  newdata,
  terms = stats::delete.response(stats::terms(model))
) {
  newdata <- normalize_vglm_prediction_data(model, newdata)
  frame <- stats::model.frame(
    terms,
    newdata,
    na.action = stats::na.pass,
    xlev = model@xlevels
  )
  X <- stats::model.matrix(
    terms,
    frame,
    contrasts.arg = if (length(model@contrasts)) model@contrasts else NULL
  )
  if (anyNA(X)) {
    columns <- colnames(X)[colSums(is.na(X)) > 0]
    mostr_abort(
      "The vglm design matrix has missing values in ",
      sum(rowSums(is.na(X)) > 0),
      " of the ",
      nrow(X),
      " prediction rows, in columns ",
      paste0("`", columns, "`", collapse = ", "),
      ". Supply prediction rows with a value of every column the model uses."
    )
  }
  X
}

predict_vglm_response_markov <- function(object, newdata) {
  X <- vglm_prediction_design(object, newdata)
  attr(X, "assign") <- object@assign

  lm2vlm <- utils::getFromNamespace("lm2vlm.model.matrix", "VGAM")
  X_vlm <- lm2vlm(
    X,
    Hlist = object@constraints,
    M = object@misc$M,
    xij = object@control$xij,
    Xm2 = NULL
  )

  eta <- matrix(
    X_vlm %*% VGAM::coefvlm(object),
    nrow = nrow(X),
    ncol = object@misc$M,
    byrow = TRUE
  )

  object@family@linkinv(eta = eta, extra = object@extra)
}

predict_orm_response_markov <- function(object, newdata) {
  Gamma <- get_effective_coefs(object)
  X <- orm_model_matrix(object, newdata, include_intercept = TRUE)
  X <- X[, colnames(Gamma), drop = FALSE]
  probs <- lp_to_probs(X %*% t(Gamma), nrow(Gamma))
  colnames(probs) <- as_state_labels(object$yunique)
  probs
}

select_posterior_draws <- function(model, n_draws = 100L, seed = NULL) {
  n_available <- nrow(model$draws)
  if (is.null(n_available) || n_available < 1) {
    mostr_abort("`blrm` model does not contain posterior draws.")
  }

  if (is.null(n_draws)) {
    n <- n_available
  } else {
    if (
      !is.numeric(n_draws) ||
        length(n_draws) != 1 ||
        is.na(n_draws) ||
        n_draws < 1
    ) {
      mostr_abort("`n_draws` must be a positive integer or NULL.")
    }
    n <- min(as.integer(n_draws), n_available)
  }

  if (n == n_available) {
    return(seq_len(n_available))
  }

  if (!is.null(seed)) {
    old_seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    old_seed <- if (old_seed_exists) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit(
      {
        if (old_seed_exists) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }

  sample.int(n_available, n)
}

blrm_design_matrix <- function(model, newdata, second = FALSE) {
  if (!requireNamespace("rms", quietly = TRUE)) {
    mostr_abort("Package 'rms' is required for `blrm` prediction.")
  }
  X <- rms::predictrms(model, newdata = newdata, type = "x", second = second)
  normalize_blrm_design_names(model, X, second = second)
}

normalize_blrm_design_names <- function(model, X, second = FALSE) {
  if (second) {
    return(X)
  }

  design <- model$Design
  rms_names <- design$colnames
  mm_names <- design$mmcolnames
  if (
    !is.null(rms_names) &&
      !is.null(mm_names) &&
      length(rms_names) == ncol(X) &&
      length(mm_names) == ncol(X) &&
      identical(colnames(X), mm_names)
  ) {
    colnames(X) <- rms_names
  }

  X
}

# The posterior draws of the random intercepts, one row per requested draw and
# one column per fitted cluster, with every column named by its cluster value
# (`name_blrm_gamma_draws()`). `model$gamma_draws` stands in for the Stan fit
# in tests.
get_blrm_gamma_draws <- function(model, draw_indices) {
  if (!is.null(model$gamma_draws)) {
    gamma <- name_blrm_gamma_draws(model$gamma_draws, model)
    return(gamma[draw_indices, , drop = FALSE])
  }

  if (!requireNamespace("rmsb", quietly = TRUE)) {
    mostr_abort("Package 'rmsb' is required for `include_re = TRUE`.")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    mostr_abort("Package 'posterior' is required for `include_re = TRUE`.")
  }

  stan_fit <- rmsb::stanGet(model)
  gamma <- tryCatch(
    {
      posterior::as_draws_matrix(stan_fit$draws(variables = "gamma"))
    },
    error = function(e) {
      mat <- as.matrix(stan_fit)
      mat[, grep("^gamma\\[", colnames(mat)), drop = FALSE]
    }
  )

  if (ncol(gamma) == 0) {
    mostr_abort(
      "No random-effect draws named `gamma` were found in the `blrm` fit."
    )
  }

  # A plain matrix, so that the columns can be renamed to cluster values and
  # selected by name without the checks a `draws_matrix` makes on its
  # variable names.
  gamma <- matrix(
    as.vector(gamma),
    nrow = nrow(gamma),
    ncol = ncol(gamma),
    dimnames = list(NULL, colnames(gamma))
  )
  gamma <- name_blrm_gamma_draws(gamma, model)
  gamma[draw_indices, , drop = FALSE]
}

# How rmsb::blrm() numbered the clusters of the fitted data. It passes
# `as.integer(as.factor(cluster))` to Stan with `Nc`, the largest of these
# numbers, as the number of clusters. So cluster number `k` is the `k`-th
# factor level (`levels`), and Stan has draws for the numbers 1 to `n`. `n`
# is smaller than the number of levels when the cluster variable is a factor
# whose last levels do not occur in the fitted data.
blrm_cluster_numbering <- function(model) {
  cluster <- model$clusterInfo$cluster
  if (is.null(cluster) || length(cluster) == 0) {
    mostr_abort(
      "`include_re = TRUE` requires a `blrm` model fitted with `cluster()` ",
      "and stored cluster information."
    )
  }
  cluster <- as.factor(cluster)
  list(levels = levels(cluster), n = max(as.integer(cluster), na.rm = TRUE))
}

# Names every column of the random-intercept draws by the cluster value it
# belongs to, so that prediction rows are matched to columns by that value
# only (`blrm_random_effect_matrix()`).
#
# Stan names the draws `gamma[1]`, ..., `gamma[Nc]`, and `gamma[k]` belongs to
# the `k`-th factor level of the fitted cluster values, not to the cluster
# whose value is `k`. Matching patient IDs against the number in the brackets
# gave patients the draws of other patients without any message (#39):
# character IDs "1" to "10" sort as "1", "10", "2", ..., so `gamma[2]` is
# patient "10", and with numeric IDs 2 to 11 `gamma[2]` is patient 3.
#
# - Columns named as Stan names them are renamed through the number in the
#   brackets, not through their position.
# - Columns without names are taken to be in the order of the factor levels.
# - Columns with any other names are taken to be named by cluster value
#   already and are returned unchanged.
name_blrm_gamma_draws <- function(gamma, model) {
  gamma_names <- colnames(gamma)
  stan_pattern <- "^gamma\\[([0-9]+)\\]$"
  stan_named <- !is.null(gamma_names) && all(grepl(stan_pattern, gamma_names))
  if (!is.null(gamma_names) && !stan_named) {
    return(gamma)
  }

  numbering <- blrm_cluster_numbering(model)
  stan_index <- if (stan_named) {
    as.integer(sub(stan_pattern, "\\1", gamma_names))
  } else {
    seq_len(ncol(gamma))
  }
  # Every cluster number of the fit needs exactly one column. Anything else
  # means the draws and the stored cluster values come from different fits,
  # and naming the columns anyway would attach draws to the wrong clusters.
  if (!identical(sort(stan_index), seq_len(numbering$n))) {
    mostr_abort(
      "The `blrm` fit has random-effect draws for ",
      ncol(gamma),
      if (ncol(gamma) == 1L) " cluster" else " clusters",
      if (stan_named) {
        paste0(
          " (`",
          paste(utils::head(gamma_names, 3), collapse = "`, `"),
          if (length(gamma_names) > 3L) "`, ...)" else "`)"
        )
      },
      ", but the cluster values stored with the fit are numbered 1 to ",
      numbering$n,
      ", so the draws cannot be matched to clusters."
    )
  }
  colnames(gamma) <- numbering$levels[stan_index]
  gamma
}

cache_blrm_gamma_draws <- function(model, draw_indices, include_re) {
  if (!isTRUE(include_re)) {
    return(NULL)
  }

  gamma <- get_blrm_gamma_draws(model, draw_indices)
  rownames(gamma) <- as.character(draw_indices)
  gamma
}

subset_cached_draw_matrix <- function(x, all_draw_indices, draw_indices) {
  if (is.null(x)) {
    return(NULL)
  }

  idx <- match(draw_indices, all_draw_indices)
  if (anyNA(idx)) {
    mostr_abort(
      "Cached posterior draw matrix is missing requested draw indices."
    )
  }

  x[idx, , drop = FALSE]
}

# The column of the prediction rows that holds each row's cluster value: the
# variable inside `cluster()` of the fitted formula, which rmsb records as
# `clusterInfo$name`. It is usually the patient ID, but for a fit with, say,
# `cluster(site)` the random intercepts belong to sites, and looking them up
# by patient ID would give wrong intercepts or none. `id_var` is used only
# when the fit does not record the name.
blrm_cluster_var <- function(model, id_var = NULL) {
  cluster_var <- model$clusterInfo$name
  if (
    is.character(cluster_var) &&
      length(cluster_var) == 1L &&
      !is.na(cluster_var) &&
      nzchar(cluster_var)
  ) {
    return(cluster_var)
  }
  id_var
}

# One column of random-intercept draws for each row of `newdata`: the draws
# of the cluster that the row belongs to. Rows are matched to draw columns by
# cluster value and column name only, never by Stan's cluster number; see
# `name_blrm_gamma_draws()`.
blrm_random_effect_matrix <- function(
  model,
  newdata,
  id_var,
  draw_indices,
  gamma_draws = NULL
) {
  cluster <- model$clusterInfo$cluster
  if (is.null(cluster) || length(cluster) == 0) {
    mostr_abort(
      "`include_re = TRUE` requires a `blrm` model fitted with `cluster()` ",
      "and stored cluster information."
    )
  }
  cluster_var <- blrm_cluster_var(model, id_var)
  if (!cluster_var %in% names(newdata)) {
    mostr_abort(
      "`newdata` must contain `",
      cluster_var,
      "`, the variable in `cluster()` of the fitted model, when ",
      "`include_re = TRUE`: each row gets the random effect of its value of ",
      "that variable."
    )
  }

  cluster_values <- as.character(newdata[[cluster_var]])

  gamma <- gamma_draws %||% get_blrm_gamma_draws(model, draw_indices)
  if (nrow(gamma) != length(draw_indices)) {
    mostr_abort(
      "Random-effect draw matrix does not match requested posterior draws."
    )
  }
  # Draws from `get_blrm_gamma_draws()` are named already; this names a
  # matrix that a caller passed without names or with Stan's names.
  gamma <- name_blrm_gamma_draws(gamma, model)

  gamma_index <- match(cluster_values, colnames(gamma))
  if (anyNA(gamma_index)) {
    missing_values <- unique(cluster_values[is.na(gamma_index)])
    mostr_abort(
      "`newdata` contains IDs not present in the fitted `blrm` clusters: ",
      paste(utils::head(missing_values, 5), collapse = ", "),
      if (length(missing_values) > 5) " ..." else "",
      ". These are values of `",
      cluster_var,
      "`, the variable in `cluster()` of the fitted model; only clusters ",
      "that were in the fitted data have an estimated random effect."
    )
  }

  gamma[, gamma_index, drop = FALSE]
}

resolve_blrm_cppo <- function(cppo, envir = parent.frame()) {
  if (is.function(cppo)) {
    return(cppo)
  }

  if (
    !is.character(cppo) || length(cppo) != 1L || is.na(cppo) || !nzchar(cppo)
  ) {
    mostr_abort(
      "Only constrained partial proportional odds `blrm` models are supported."
    )
  }

  is_syntactic_name <- make.names(cppo) == cppo && !grepl("^\\.[0-9]", cppo)
  if (!is_syntactic_name) {
    mostr_abort(
      "String `cppo` values must name an available function. ",
      "Inline expressions are not supported."
    )
  }

  tryCatch(
    get(cppo, mode = "function", envir = envir, inherits = TRUE),
    error = function(e) {
      mostr_abort(
        "String `cppo` values must name an available function. ",
        "Could not find function `",
        cppo,
        "`."
      )
    }
  )
}

# Scores that scale the partial proportional-odds effects at the cut-points
# "y >= state j", j = 2..k, of a constrained blrm fit. rmsb::blrm() computes
# `pposcore <- cppo(ylevels)` once on all states, and its Stan model uses
# `pposcore[j]` for cut-point j, so the scores are `cppo(ylevels)[-1]`.
# Calling `cppo` on one state at a time gives other scores whenever a state's
# score depends on the other states: on a real 10-state fit with `cppo = rank`
# these reproduced Stan's stored log-likelihood to 5e-14, while per-state calls
# missed it by 2.7 (issue #46). rmsb::predict.blrm() itself uses
# `cppo(ylevels[-1])`, which missed its own fit by 0.35 for such functions, so
# it is not followed here.
blrm_cppo_scores <- function(cppo, ylevels) {
  scores <- cppo(ylevels)
  n_states <- length(ylevels)
  valid <- (is.numeric(scores) || is.logical(scores)) &&
    length(scores) == n_states
  if (valid) {
    scores <- as.double(scores)
    # The first state's score is never used: it belongs to no cut-point.
    valid <- all(is.finite(scores[-1L]))
  }
  if (!valid) {
    mostr_abort(
      "`cppo` must return one finite number for each state when it is given ",
      "the labels of all ",
      n_states,
      " states at once, as `rmsb::blrm()` calls it. For a factor response ",
      "the labels are character strings, so convert them, for example with ",
      "`cppo = function(y) as.numeric(y)`."
    )
  }
  scores[-1L]
}

predict_blrm_response_markov <- function(
  object,
  newdata,
  include_re = FALSE,
  id_var = NULL,
  draw_indices = NULL,
  gamma_draws = NULL,
  prediction_cache = NULL
) {
  if (is.null(draw_indices)) {
    draw_indices <- seq_len(nrow(object$draws))
  }
  cache_index <- NULL
  if (!is.null(prediction_cache)) {
    prediction_cache$cursor <- prediction_cache$cursor + 1L
    cache_index <- prediction_cache$cursor
  }
  X <- if (!is.null(cache_index) && length(prediction_cache$X) >= cache_index) {
    prediction_cache$X[[cache_index]]
  } else {
    value <- blrm_design_matrix(object, newdata = newdata, second = FALSE)
    if (!is.null(cache_index)) {
      prediction_cache$X[[cache_index]] <- value
    }
    value
  }
  draws <- object$draws[draw_indices, , drop = FALSE]
  ndraws <- nrow(draws)
  ns <- object$non.slopes
  cn <- colnames(draws)
  tauinfo <- object$tauInfo
  tau_names <- if (!is.null(tauinfo) && length(tauinfo$name) > 0L) {
    tauinfo$name
  } else {
    character()
  }

  intercept_draws <- draws[, seq_len(ns), drop = FALSE]
  beta_names <- setdiff(cn, c(cn[seq_len(ns)], tau_names))
  if (length(beta_names) == 0) {
    xb <- matrix(0, nrow = ndraws, ncol = nrow(newdata))
  } else {
    beta_draws <- draws[, beta_names, drop = FALSE]
    missing_beta <- setdiff(beta_names, colnames(X))
    if (length(missing_beta) > 0) {
      mostr_abort(
        "Could not match `blrm` posterior coefficients to prediction ",
        "design columns. Missing columns: ",
        paste(utils::head(missing_beta, 10), collapse = ", "),
        if (length(missing_beta) > 10) " ..." else "",
        "."
      )
    }
    X <- X[, beta_names, drop = FALSE]
    xb <- beta_draws %*% t(X)
  }

  pppo <- object$pppo %||% 0L
  has_npo <- pppo > 0
  if (has_npo) {
    cppo <- resolve_blrm_cppo(object$cppo, envir = parent.frame())
    Z <- if (
      !is.null(cache_index) && length(prediction_cache$Z) >= cache_index
    ) {
      prediction_cache$Z[[cache_index]]
    } else {
      value <- blrm_design_matrix(object, newdata = newdata, second = TRUE)
      if (!is.null(cache_index)) {
        prediction_cache$Z[[cache_index]] <- value
      }
      value
    }
    tau_draws <- draws[, tau_names, drop = FALSE]
    zt <- tau_draws %*% t(Z)
    cppos <- blrm_cppo_scores(cppo, object$ylevels)
  }

  u_draws <- matrix(0, nrow = ndraws, ncol = nrow(newdata))
  if (include_re) {
    id_var <- markov_require_id_var(object, id_var)
    u_draws <- blrm_random_effect_matrix(
      object,
      newdata,
      id_var,
      draw_indices,
      gamma_draws = gamma_draws
    )
  }

  ynam <- if (!is.null(object$yname)) {
    paste(object$yname, object$ylevels, sep = "=")
  } else {
    as_state_labels(object$ylevels)
  }
  base_eta <- xb + u_draws
  out <- blrm_probabilities_native(
    base_eta,
    intercept_draws,
    threshold_eta = if (has_npo) zt else numeric(),
    threshold_scale = if (has_npo) cppos else numeric()
  )
  dimnames(out) <- list(draw_indices, rownames(newdata), ynam)
  out
}

orm_model_matrix <- function(model, newdata, include_intercept = FALSE) {
  if (!inherits(model, "orm")) {
    mostr_abort("orm_model_matrix() requires an orm model.")
  }
  if (is.null(model$Design) || is.null(model$sformula)) {
    mostr_abort(
      "Fast orm prediction requires a model fitted by rms::orm() with ",
      "stored Design metadata."
    )
  }

  design <- model$Design
  newdata <- normalize_orm_prediction_data(model, newdata)

  oldopts <- options(
    contrasts = c(factor = "contr.treatment", ordered = "contr.treatment"),
    Design.attr = design
  )
  on.exit(
    {
      options(contrasts = oldopts$contrasts)
      options(Design.attr = oldopts$Design.attr)
    },
    add = TRUE
  )

  formulano <- add_rms_formula_helpers(model$sformula)
  Terms <- stats::delete.response(
    stats::terms(formulano, specials = "strat")
  )
  attr(Terms, "intercept") <- 1L

  strata_terms <- attr(Terms, "specials")$strat
  Terms_ns <- if (length(strata_terms)) {
    Terms[-strata_terms]
  } else {
    Terms
  }

  mf <- stats::model.frame(Terms, newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(Terms_ns, mf)[, -1L, drop = FALSE]

  expected_cols <- design$colnames
  if (length(expected_cols) != ncol(X)) {
    mostr_abort(
      "Could not reconstruct the orm design matrix. Expected ",
      length(expected_cols),
      " columns, got ",
      ncol(X),
      "."
    )
  }
  colnames(X) <- expected_cols

  if (include_intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  X
}

normalize_orm_prediction_data <- function(model, newdata) {
  design <- model$Design
  categorical <- design$name[design$assume.code %in% c(5L, 8L)]
  fitted_levels <- lapply(
    stats::setNames(nm = categorical),
    function(nm) as.character(design$parms[[nm]])
  )
  factors_with_fitted_levels(newdata, fitted_levels, "orm")
}

# The prediction rows with each predictor that the vglm fit used as a factor
# turned into a factor with the fitted levels. stats::model.frame() with
# `xlev` leaves a numeric column numeric, so numbers such as
# `variables = list(txf = c(0, 1))` or the previous states of
# plot_transitions() made model.matrix() stop with "contrasts apply only to
# factors" (issue #44), while the orm engine accepts them. The fitted levels
# come from `model@xlevels`, which names the variables of the model frame:
# a term such as factor(tx) is listed as "factor(tx)", which no prediction
# column is named, so its column `tx` stays as given and is converted by the
# formula.
normalize_vglm_prediction_data <- function(model, newdata) {
  fitted_levels <- tryCatch(model@xlevels, error = function(e) NULL)
  factors_with_fitted_levels(newdata, fitted_levels, "vglm")
}

# Each column of `newdata` named in `fitted_levels` (a list of level vectors)
# as a factor with those levels, whatever its type: numbers and characters
# are matched to the levels by their labels. Values that are not fitted levels
# are refused, naming the column, because they would otherwise become missing
# values or an error from deep inside the engine. Missing values stay missing,
# for the refusals that name them later. An ordered factor stays ordered; the
# fitted contrasts are passed to model.matrix(), so a column that arrives
# plain gets the same design.
factors_with_fitted_levels <- function(newdata, fitted_levels, engine) {
  for (nm in intersect(names(fitted_levels), names(newdata))) {
    values <- newdata[[nm]]
    coerced <- factor(
      as.character(values),
      levels = as.character(fitted_levels[[nm]]),
      ordered = is.ordered(values)
    )
    bad <- is.na(coerced) & !is.na(values)
    if (any(bad)) {
      bad_values <- unique(as.character(values[bad]))
      mostr_abort(
        "Values in `",
        nm,
        "` are not among fitted ",
        engine,
        " levels: ",
        paste(utils::head(bad_values, 5), collapse = ", "),
        if (length(bad_values) > 5) " ..." else ""
      )
    }
    newdata[[nm]] <- coerced
  }

  newdata
}

as_state_labels <- function(x) {
  as.character(x)
}

coerce_previous_state_values <- function(values, prototype, p_var) {
  labels <- as_state_labels(values)

  if (is.factor(prototype)) {
    return(factor(
      labels,
      levels = levels(prototype),
      ordered = is.ordered(prototype)
    ))
  }

  if (is.integer(prototype)) {
    out <- suppressWarnings(as.integer(labels))
    if (anyNA(out) && any(!is.na(labels))) {
      mostr_abort(
        "The fitted states must be integer-compatible when the ",
        "previous-state column `",
        p_var,
        "` holds integers."
      )
    }
    return(out)
  }

  if (is.numeric(prototype)) {
    out <- suppressWarnings(as.numeric(labels))
    if (anyNA(out) && any(!is.na(labels))) {
      mostr_abort(
        "The fitted states must be numeric-compatible when the ",
        "previous-state column `",
        p_var,
        "` holds numbers."
      )
    }
    return(out)
  }

  if (is.character(prototype)) {
    return(labels)
  }

  labels
}

# The function takes the previous-state column a user supplied for prediction
# and converts it to the type and levels the fitted model expects, refusing any value
# that is not a fitted state.
# `values` of the previous-state column `column` (`p_var`, or `p2_var` for the
# state two visits back) in the form the model reads them, whatever their own
# type or factor levels. The recursion builds the previous states of later
# visits with the type and levels of the starting rows' column
# (`make_previous_state_column()`), so a factor holding only the profiles' own
# states turned every other state into NA: orm then stopped with an error
# about missing transition probabilities, vglm with a dimension error from
# src/sops.cpp. A column the model reads as a factor gets its fitted levels in
# fitted order; one it reads as a number gets numbers; with neither recorded
# (a fit without its metadata) a factor gets the fitted states as levels.
# Values that are not fitted states are refused, naming the column, because
# vglm would otherwise stop in stats::model.frame() ("factor yprev has new
# levels 9") and blrm would predict from a state the model never saw.
# Missing values are left for stop_missing_predictors().

fitted_state_values <- function(model, values, column, states) {
  labels <- as_state_labels(values)
  fitted <- markov_model_refit_data(model)[[column]]
  levels <- get_model_factor_levels(model, column) %||%
    (if (is.factor(fitted)) levels(fitted))
  allowed <- levels %||% as_state_labels(states)
  unknown <- !is.na(values) & !labels %in% allowed
  if (any(unknown)) {
    shown <- unique(labels[unknown])
    mostr_abort(
      "`",
      column,
      "` holds ",
      paste(utils::head(shown, 5L), collapse = ", "),
      if (length(shown) > 5L) ", ...",
      " in ",
      sum(unknown),
      " of the ",
      length(values),
      " prediction rows, which ",
      if (length(shown) == 1L) "is not a state" else "are not states",
      " the model can predict from. The states of `",
      column,
      "` in the fitted model are ",
      paste(allowed, collapse = ", "),
      "."
    )
  }
  if (!is.null(levels)) {
    return(factor(
      labels,
      levels = levels,
      ordered = is.ordered(if (is.null(fitted)) values else fitted)
    ))
  }
  if (is.numeric(fitted) && !is.numeric(values)) {
    return(if (is.integer(fitted)) as.integer(labels) else as.numeric(labels))
  }
  if (is.null(fitted) && is.factor(values)) {
    return(factor(labels, levels = allowed, ordered = is.ordered(values)))
  }
  values
}

# The prediction rows with their previous-state columns in the form the model
# reads them (`fitted_state_values()`).
normalize_prediction_states <- function(model, data, settings) {
  for (column in c(settings$p_var, settings$p2_var)) {
    if (column %in% names(data)) {
      data[[column]] <- fitted_state_values(
        model,
        data[[column]],
        column,
        settings$y_levels
      )
    }
  }
  data
}

make_previous_state_column <- function(states, prototype, n, p_var) {
  values <- coerce_previous_state_values(states, prototype, p_var)
  rep(values, each = n)
}

get_model_factor_levels <- function(model, varname) {
  model_chk <- if (inherits(model, "robcov_vglm")) {
    model$vglm_fit
  } else {
    model
  }

  if (inherits(model_chk, "vglm")) {
    xlevels <- tryCatch(
      model_chk@xlevels,
      error = function(e) NULL
    )
    if (!is.null(xlevels[[varname]])) {
      return(as.character(xlevels[[varname]]))
    }
  }

  if (inherits(model_chk, c("orm", "blrm"))) {
    design <- tryCatch(model_chk$Design, error = function(e) NULL)
    # rms keeps other settings in `parms` too, such as the knots of
    # rcs(yprev, 4), which are not levels. Only categorical terms (the codes
    # that normalize_orm_prediction_data() turns into factors) have levels.
    categorical <- design$name[design$assume.code %in% c(5L, 8L)]
    if (varname %in% categorical && !is.null(design$parms[[varname]])) {
      return(as.character(design$parms[[varname]]))
    }
  }

  NULL
}

get_time_info <- function(model, data, time_var) {
  values <- if (!is.null(time_var) && time_var %in% names(data)) {
    data[[time_var]]
  } else {
    NULL
  }
  model_levels <- if (!is.null(model)) {
    get_model_factor_levels(model, time_var)
  } else {
    NULL
  }
  fitted_time <- if (!is.null(model) && !is.null(time_var)) {
    markov_model_refit_data(model)[[time_var]]
  }

  # When the fit's time column is a factor, the fitted levels decide, whatever
  # the type of the prediction rows' time column: a vglm fit with factor time
  # given numeric or character times otherwise treated time as numeric, so
  # requested `times` were neither checked against the fitted levels nor
  # turned into factors (issue #44). Only then: the engines record levels for
  # character time in the formula as well, sorted alphabetically, and taking
  # them would run character visits in that order instead of the requested
  # one (issue #51). So with any other fitted time, or a fit without its
  # stored data, the rows' time column decides: a factor, or no column, takes
  # the fitted levels; numbers and characters stay as they are.
  levels <- NULL
  ordered <- FALSE
  if (is.factor(fitted_time) || is.factor(values) || is.null(values)) {
    levels <- model_levels %||% (if (is.factor(values)) levels(values))
    ordered <- is.ordered(values)
  }

  list(
    is_factor = !is.null(levels),
    levels = levels,
    ordered = ordered
  )
}

default_sop_times <- function(data, time_var, default) {
  if (default == "fast") {
    if (!is.null(time_var) && time_var %in% names(data)) {
      return(sort(unique(data[[time_var]])))
    }
    return(1)
  }

  if (is.null(time_var) || !time_var %in% names(data)) {
    mostr_abort("`times` must be specified if `time_var` is not in data.")
  }

  sort(unique(data[[time_var]]))
}

coerce_visit_times <- function(times, time_info, time_var) {
  labels <- as.character(times)
  bad <- is.na(match(labels, time_info$levels))
  if (any(bad)) {
    bad_values <- unique(labels[bad])
    mostr_abort(
      "Requested visit values are not among fitted factor levels for `",
      time_var,
      "`: ",
      paste(utils::head(bad_values, 5), collapse = ", "),
      if (length(bad_values) > 5) " ..." else ""
    )
  }

  factor(
    labels,
    levels = time_info$levels,
    ordered = isTRUE(time_info$ordered)
  )
}

resolve_sop_times <- function(
  model,
  data,
  times,
  time_var,
  default = c("unique", "fast") #* needed only for testing
) {
  default <- match.arg(default)
  time_info <- get_time_info(model, data, time_var)

  if (is.null(times)) {
    times <- if (time_info$is_factor) {
      time_info$levels
    } else {
      default_sop_times(data, time_var, default)
    }
  }

  if (time_info$is_factor) {
    times <- coerce_visit_times(times, time_info, time_var)
  }

  list(times = times, time_info = time_info)
}

assign_sop_time <- function(data, time_var, value, time_info) {
  if (isTRUE(time_info$is_factor)) {
    data[[time_var]] <- factor(
      as.character(value),
      levels = time_info$levels,
      ordered = isTRUE(time_info$ordered)
    )
  } else {
    data[[time_var]] <- value
  }
  data
}

# The gap since the previous requested visit, for a model fitted with
# `gap_var`. The first transition's gap is its own time, the time since time
# 0; `?mostm` and the "Choosing visits" section of `?sops` document this rule.
# `mostm()` accepts `gap_var` only with numeric time, so the times can be
# subtracted.
assign_sop_gap <- function(data, gap_var, times, index) {
  if (is.null(gap_var)) {
    return(data)
  }
  data[[gap_var]] <- if (index == 1L) {
    times[index]
  } else {
    times[index] - times[index - 1L]
  }
  data
}

# Set the columns of visit `index` of `times`: the time column and the gap
# since the previous visit, as named in `settings`. Time terms of the formula,
# such as `rcs(time, 4)`, are evaluated from the time column.
assign_sop_visit <- function(data, settings, times, index, time_info) {
  data <- assign_sop_time(data, settings$time_var, times[index], time_info)
  assign_sop_gap(data, settings$gap_var, times, index)
}
