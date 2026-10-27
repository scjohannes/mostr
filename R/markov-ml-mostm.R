# mostm() integration of maximum likelihood fits with missing states.

# Refits for the bootstrap need only the coefficients. update_bootstrap_model()
# sets this flag while it refits, so mostm() skips the covariances.
markov_ml_flags <- new.env(parent = emptyenv())

markov_ml_coefficients_only <- function() {
  isTRUE(markov_ml_flags$coefficients_only)
}

with_markov_ml_coefficients_only <- function(expr) {
  old <- markov_ml_flags$coefficients_only
  markov_ml_flags$coefficients_only <- TRUE
  on.exit(markov_ml_flags$coefficients_only <- old, add = TRUE)
  expr
}

# Check the arguments of a fit with missing states and set up its data.
# Returns NULL when neither `missing_states` nor an Ocens() response is given,
# or when `missing_states` is given but no row needs the recursion: the
# engine's ordinary fit then maximises the same likelihood.
mostm_ml_problem <- function(
  setup,
  engine,
  id_var,
  time_var,
  p_var,
  missing_states,
  engine_args
) {
  formula <- setup$formula
  ocens <- markov_formula_has_ocens(formula)
  if (is.null(missing_states) && !ocens) {
    return(NULL)
  }
  if (
    !is.null(missing_states) &&
      (!is.atomic(missing_states) ||
        !length(missing_states) ||
        anyNA(missing_states))
  ) {
    stop(
      "`missing_states` must give the states a missing response may be, ",
      "without missing values, or be NULL.",
      call. = FALSE
    )
  }
  if (is.null(id_var)) {
    stop(
      "`id_var` is required with `missing_states` or an `Ocens()` response: ",
      "each patient's rows form one Markov chain.",
      call. = FALSE
    )
  }
  if (!is.null(setup$prior_weight_column)) {
    stop(
      "`weights` cannot be combined with `missing_states` or an `Ocens()` ",
      "response.",
      call. = FALSE
    )
  }
  if (
    !is.character(p_var) ||
      length(p_var) != 1L ||
      is.na(p_var) ||
      !p_var %in% names(setup$data)
  ) {
    stop("`p_var` must name the previous-state column of `data`.", call. = FALSE)
  }
  mt <- stats::terms(formula, specials = "cluster", data = setup$data)
  if (!p_var %in% all.vars(stats::delete.response(mt))) {
    stop(
      "The previous-state variable `",
      p_var,
      "` must appear on the right-hand side of `formula`.",
      call. = FALSE
    )
  }
  random_effects <- !is.null(attr(mt, "specials")$cluster)
  if (random_effects) {
    if (engine != "blrm") {
      stop(
        "Random effects (`cluster()`) cannot be combined with ",
        "`missing_states` or an `Ocens()` response for `engine = \"",
        engine,
        "\"`; `engine = \"blrm\"` supports a random intercept per patient.",
        call. = FALSE
      )
    }
    formula <- markov_ml_drop_cluster(mt, formula, id_var)
  }

  selected <- markov_subset_source_data(
    setup$data,
    setup$subset_rows,
    setup$eval_env
  )$data
  ppo_vars <- if (inherits(engine_args$ppo, "formula")) {
    all.vars(engine_args$ppo)
  } else {
    character()
  }
  problem <- markov_ml_prepare(
    formula = formula,
    data = selected,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    missing_states = missing_states,
    extra_vars = ppo_vars
  )
  if (!problem$chain$n && !ocens) {
    return(NULL)
  }
  markov_ml_reject_engine_args(engine, names(engine_args))
  if (
    engine == "orm" &&
      !is.null(engine_args$family) &&
      !identical(tolower(as.character(engine_args$family)), "logistic")
  ) {
    stop_unsupported_markov_link(engine_args$family)
  }
  # Spline knots from the transitions with both states observed, shared by
  # every fit of this problem.
  problem$working_formula <- markov_ml_lock_knots(
    problem$working_formula,
    problem$start_data
  )
  problem$formula <- formula
  problem$random_effects <- random_effects
  problem$selected <- selected
  problem$id_var <- id_var
  problem
}

# `formula` without its cluster(<id_var>) term.
markov_ml_drop_cluster <- function(mt, formula, id_var) {
  labels <- attr(mt, "term.labels")
  cluster_term <- grepl("^cluster\\(", labels)
  if (
    sum(cluster_term) != 1L ||
      any(grepl("cluster\\(", labels[!cluster_term]))
  ) {
    stop(
      "`formula` may contain one `cluster()` term, not inside an interaction.",
      call. = FALSE
    )
  }
  cluster_call <- str2lang(labels[cluster_term])
  if (
    length(cluster_call) != 2L ||
      !identical(paste(deparse(cluster_call[[2L]]), collapse = ""), id_var)
  ) {
    stop(
      "`cluster()` must name the patient ID column: use `cluster(",
      id_var,
      ")` for a random intercept per patient.",
      call. = FALSE
    )
  }
  fixed <- stats::formula(stats::drop.terms(
    mt,
    dropx = which(cluster_term),
    keep.response = TRUE
  ))
  environment(fixed) <- environment(formula)
  fixed
}

# Warn that rows with missing values were not used.
markov_warn_dropped_rows <- function(
  n_rows,
  response,
  previous,
  covariate,
  p_var,
  missing_states_given
) {
  total <- response + previous + covariate
  if (total <= 0) {
    return(invisible(NULL))
  }
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)
  parts <- c(
    if (response) paste(fmt(response), "with a missing response"),
    if (previous) {
      paste0(fmt(previous), " with a missing previous state (`", p_var, "`)")
    },
    if (covariate) paste(fmt(covariate), "with missing covariates")
  )
  if (length(parts) > 1L) {
    parts <- paste0(
      paste(parts[-length(parts)], collapse = ", "),
      " and ",
      parts[[length(parts)]]
    )
  }
  message <- paste0(
    "mostm() did not use ",
    fmt(total),
    " of ",
    fmt(n_rows),
    " rows because of missing values: ",
    parts,
    "."
  )
  if (response && !missing_states_given) {
    message <- paste0(
      message,
      " Without the transitions into and out of visits with a missing ",
      "state, the estimates can be biased, for example when deaths are always ",
      "recorded but other states are not. To use these visits, give the ",
      "states a patient may be in at a visit without a recorded state in ",
      "`missing_states`, such as `missing_states = 1:7` when state 8 is ",
      "death; see `?mostm`."
    )
  }
  warning(
    structure(
      class = c("mostr_dropped_rows_warning", "warning", "condition"),
      list(message = message, call = NULL)
    )
  )
}

# Counts of rows of `selected` that a complete-row fit does not use, by
# reason, and the warning about them.
markov_warn_dropped_complete_rows <- function(
  selected,
  fit_data,
  formula,
  p_var,
  missing_states_given
) {
  n_dropped <- nrow(selected) - NROW(fit_data)
  if (n_dropped <= 0L) {
    return(invisible(NULL))
  }
  vars <- intersect(all.vars(formula), names(selected))
  response_vars <- intersect(all.vars(formula[[2L]]), names(selected))
  missing_response <- if (length(response_vars)) {
    !stats::complete.cases(selected[response_vars])
  } else {
    logical(nrow(selected))
  }
  missing_previous <- !missing_response &
    (if (p_var %in% vars) is.na(selected[[p_var]]) else FALSE)
  others <- setdiff(vars, c(response_vars, p_var))
  missing_covariate <- !missing_response &
    !missing_previous &
    (if (length(others)) {
      !stats::complete.cases(selected[others])
    } else {
      FALSE
    })
  counts <- c(
    sum(missing_response),
    sum(missing_previous),
    sum(missing_covariate)
  )
  if (sum(counts) != n_dropped) {
    # Rows removed for another reason (for example by a custom `na.action`):
    # report the total as rows with missing covariates.
    counts <- c(0, 0, n_dropped)
  }
  markov_warn_dropped_rows(
    nrow(selected),
    response = counts[[1L]],
    previous = counts[[2L]],
    covariate = counts[[3L]],
    p_var = p_var,
    missing_states_given = missing_states_given
  )
}

# Add the maximum likelihood components to the engine fit and give it the
# patient-clustered covariance. `fitted_ids` are the patients of the stored
# starting profiles, the order of the stored patient scores.
markov_ml_finish <- function(fit, ml, problem, engine, id_var, fitted_ids) {
  order <- match(fitted_ids, rownames(ml$patient_scores))
  if (anyNA(order)) {
    stop("Internal error: the patient scores do not match the fitted patients.")
  }
  ml$patient_scores <- ml$patient_scores[order, , drop = FALSE]
  ml$contributing <- ml$contributing[order]
  ml$id_var <- id_var
  n_clusters <- sum(ml$contributing)

  if (engine == "orm") {
    fit$ml <- ml
    fit <- markov_ml_fit_statistics(fit, problem)
    if (!is.null(ml$var_robust)) {
      fit$var <- ml$var_robust
      fit$orig.var <- ml$var_model
      fit$clusterInfo <- list(name = id_var, n = n_clusters)
      attr(fit, "markov_robust_covariance") <- list(
        implementation = "mostr::mostm_ml",
        backend = "orm",
        cluster = NULL,
        type = ml$type,
        cadjust = ml$cadjust,
        adjustment_factor = ml$adjustment_factor,
        n_rows = ml$n_informative_rows,
        n_represented_rows = ml$n_informative_rows,
        n_parameters = length(stats::coef(fit)),
        n_clusters = n_clusters,
        bread_convention = "observed_data_information",
        weighted = FALSE,
        covariance_identity = ml$var_robust
      )
    } else {
      fit$var <- NULL
    }
    class(fit) <- unique(c("mostm_ml", class(fit)))
    return(fit)
  }

  # vglm: the patient-clustered covariance is held by a `robcov_vglm` object,
  # as for fits without missing states, with one score row per patient. The
  # model information is copied from the vglm fit as robcov_vglm() does.
  coefficients <- stats::coef(fit)
  se <- if (!is.null(ml$var_robust)) sqrt(diag(ml$var_robust))
  result <- list(
    coefficients = coefficients,
    var = ml$var_robust,
    se = se,
    z = if (!is.null(se)) coefficients / se,
    pvalues = if (!is.null(se)) 2 * stats::pnorm(-abs(coefficients / se)),
    original_vcov = ml$var_model,
    original_se = if (!is.null(ml$var_model)) sqrt(diag(ml$var_model)),
    scores = ml$patient_scores,
    bread = ml$var_model,
    bread_type = "observed_data_information",
    cluster = rownames(ml$patient_scores),
    n = ml$n_informative_rows,
    n_clusters = n_clusters,
    type = ml$type,
    cadjust = ml$cadjust,
    adjustment_factor = ml$adjustment_factor,
    adjust = ml$cadjust,
    vglm_fit = fit,
    ml = ml
  )
  for (name in c(
    "fitted.values",
    "residuals",
    "family",
    "predictors",
    "prior.weights",
    "df.residual",
    "df.total",
    "rank",
    "criterion",
    "constraints",
    "control",
    "extra",
    "misc",
    "model",
    "x",
    "y",
    "terms",
    "assign",
    "xlevels",
    "na.action",
    "call"
  )) {
    result[name] <- list(tryCatch(
      methods::slot(fit, name),
      error = function(e) NULL
    ))
  }
  class(result) <- c("mostm_ml", "robcov_vglm")
  for (name in c("markov_vglm", "markov_split_assign", "markov_basis_terms")) {
    attr(result, name) <- attr(fit, name, exact = TRUE)
  }
  result
}

# rms::orm() computed `stats`, `freq`, and `deviance` of the returned fit on
# the expected transition counts. Replace them with values that describe the
# data: `Obs` is nobs(), the frequencies count exactly observed states, and
# the deviance is -2 times the observed-data log-likelihood (there is no
# intercept-only fit, so the null deviance is NA). Statistics that compare
# with an intercept-only model (likelihood-ratio and score tests, R2) or
# compare predictions with observed states row by row (rank correlations)
# are NA; the R2 variants named after the refit's sample sizes are dropped.
markov_ml_fit_statistics <- function(fit, problem) {
  scale <- problem$scale
  if (length(fit$freq) != scale$K) {
    stop("Internal error: the orm fit does not have one level per state.")
  }
  observed <- scale$lo[scale$exact]
  freq <- tabulate(observed, nbins = scale$K)
  names(freq) <- names(fit$freq)
  fit$freq <- freq

  stats <- fit$stats
  stats <- stats[!startsWith(names(stats), "R2(")]
  stats[] <- NA_real_
  stats[c("Obs", "Distinct Y", "Median Y", "Max Deriv", "d.f.")] <- c(
    fit$ml$n_informative_rows,
    scale$K,
    stats::median(scale$values[observed]),
    fit$ml$max_abs_score,
    length(stats::coef(fit)) - (scale$K - 1L)
  )
  fit$stats <- stats
  fit$deviance <- c(NA_real_, -2 * fit$ml$loglik)
  fit
}

# Maximum likelihood components of an orm or vglm fit with missing states, or
# NULL for other fits.
markov_ml_components <- function(model) {
  if (!inherits(model, "mostm_ml") || inherits(model, "blrm")) {
    return(NULL)
  }
  model$ml
}

markov_is_ml_fit <- function(model) {
  !is.null(markov_ml_components(model))
}

# TRUE when stored `mostm()` arguments ask for missing states: their refits
# must keep the patients' chains.
markov_refit_uses_ml <- function(refit_args) {
  !is.null(refit_args$missing_states) ||
    markov_formula_has_ocens(refit_args$formula)
}

stop_markov_ml_fwb <- function(model) {
  if (
    inherits(model, "mostm_ml") ||
      markov_refit_uses_ml(markov_model_refit_args(model))
  ) {
    stop(
      "The fractional weighted bootstrap is not available for fits with ",
      "`missing_states` or an `Ocens()` response, which do not accept row ",
      "weights. Use `method = \"bootstrap\"`, `\"delta\"`, `\"mvn\"`, or ",
      "`\"score_bootstrap\"`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

# The engine print methods are not used: their statistics (frequencies,
# likelihood-ratio tests, R2, rank correlations) would describe the fit to
# expected transition counts, not the data.
#' @export
print.mostm_ml <- function(x, digits = 4, ...) {
  ml <- x$ml
  if (inherits(x, "blrm")) {
    return(markov_ml_print_bayes(x, digits = digits))
  }
  cat(
    "Markov model fitted by maximum likelihood with missing states ",
    "(engine \"",
    markov_model_engine(x) %||% "?",
    "\")\n\n",
    sep = ""
  )
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat(
    "  Log-likelihood: ",
    format(ml$loglik, nsmall = 4),
    "; largest |score|: ",
    format(signif(ml$max_abs_score, 3)),
    "\n  Patients: ",
    ml$n_patients,
    "; patient clusters for robust standard errors: ",
    sum(ml$contributing),
    "\n  Rows: ",
    ml$n_rows,
    " (",
    ml$n_missing,
    " missing, ",
    ml$n_range,
    " with a range of states)",
    if (length(ml$missing_states)) {
      paste0(
        "\n  A missing state may be: ",
        paste(ml$missing_states, collapse = ", ")
      )
    },
    "\n  Share of information lost to missing states: ",
    format(signif(ml$missing_info_fraction, 3)),
    "\n\n",
    sep = ""
  )
  theta <- get_coef(x)
  if (is.null(ml$var_robust)) {
    cat("Coefficients (standard errors not computed):\n")
    print(cbind(Coef = round(theta, digits)))
    return(invisible(x))
  }
  se <- sqrt(diag(ml$var_robust))
  z <- theta / se
  table <- cbind(
    Coef = theta,
    `Robust S.E.` = se,
    `Wald Z` = z,
    `Pr(>|Z|)` = 2 * stats::pnorm(-abs(z))
  )
  cat(
    "Coefficients (patient-clustered robust standard errors, ",
    ml$type,
    if (isTRUE(ml$cadjust)) " with n/(n-1)",
    "):\n",
    sep = ""
  )
  if (!isTRUE(ml$information_positive_definite)) {
    cat(
      "  (unavailable: the observed information is not positive definite)\n"
    )
  }
  stats::printCoefmat(
    table,
    digits = digits,
    signif.stars = FALSE,
    P.values = TRUE,
    has.Pvalue = TRUE,
    eps.Pvalue = 1e-4,
    na.print = "NA"
  )
  invisible(x)
}

markov_ml_print_bayes <- function(x, digits = 4) {
  ml <- x$ml
  cat(
    "Bayesian Markov model with missing states (engine \"blrm\")\n",
    "  Posterior draws: ",
    nrow(x$draws),
    " (",
    x$chains,
    " chains); sampling time: ",
    format(round(x$sampling_time, 1)),
    "s\n  Rows: ",
    ml$n_rows,
    " (",
    ml$n_missing,
    " missing, ",
    ml$n_range,
    " with a range of states)",
    if (length(ml$missing_states)) {
      paste0(
        "\n  A missing state may be: ",
        paste(ml$missing_states, collapse = ", ")
      )
    },
    "\n  Dirichlet prior concentration for the intercepts: ",
    format(signif(x$conc, 3)),
    if (!is.null(x$Contrast) || !is.null(x$nContrast)) {
      "; normal priors for contrasts given"
    },
    if (x$pppo > 0) "\n  Constrained partial proportional odds",
    if (!is.null(x$clusterInfo)) {
      paste0("\n  Random intercepts for ", x$clusterInfo$n, " patients")
    },
    "\n\n",
    sep = ""
  )
  draws <- x$draws
  table <- cbind(
    Mean = colMeans(draws),
    SD = apply(draws, 2L, stats::sd),
    `Lower 0.95` = apply(draws, 2L, stats::quantile, probs = 0.025),
    `Upper 0.95` = apply(draws, 2L, stats::quantile, probs = 0.975),
    `Pr(>0)` = colMeans(draws > 0)
  )
  diagnostics <- x$diagnostics[colnames(draws), , drop = FALSE]
  table <- cbind(
    round(table, digits),
    Rhat = round(diagnostics[, "Rhat"], 3),
    ess_bulk = round(diagnostics[, "ess_bulk"]),
    ess_tail = round(diagnostics[, "ess_tail"])
  )
  print(table)
  if (!is.null(x$omega)) {
    sigmag <- x$omega[, "sigmag"]
    cat(
      "\nSD of the patient random intercepts (sigmag): mean ",
      format(round(mean(sigmag), digits)),
      ", 95% interval [",
      format(round(stats::quantile(sigmag, 0.025), digits)),
      ", ",
      format(round(stats::quantile(sigmag, 0.975), digits)),
      "], Rhat ",
      format(round(x$diagnostics["sigmag", "Rhat"], 3)),
      ", ess_bulk ",
      round(x$diagnostics["sigmag", "ess_bulk"]),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

# The maximised observed-data log-likelihood, with one degree of freedom per
# coefficient; AIC() and BIC() use it.
#' @export
logLik.mostm_ml <- function(object, ...) {
  if (inherits(object, "blrm")) {
    return(NextMethod())
  }
  structure(
    object$ml$loglik,
    nobs = stats::nobs(object),
    df = length(get_coef(object)),
    class = "logLik"
  )
}

# Rows whose observation rules out at least one state: the transitions a
# complete-row fit counts when nothing is missing, and `m` of the HC1
# correction.
#' @export
nobs.mostm_ml <- function(object, ...) {
  if (inherits(object, "blrm")) {
    return(NextMethod())
  }
  object$ml$n_informative_rows
}

# `missing_states` of a bootstrap refit whose resample lacks some states, so
# that its states were renumbered consecutively by `state_mapping` (old label
# to new number). States absent from the resample are dropped. Ocens() ranges
# are not renumbered, so such a refit stops; the bootstrap then drops that
# resample with a warning.
markov_ml_refit_missing_states <- function(refit_args, state_mapping) {
  if (markov_formula_has_ocens(refit_args$formula)) {
    stop(
      "A bootstrap resample lacks state(s) of the data, and the states of an ",
      "`Ocens()` response cannot be renumbered.",
      call. = FALSE
    )
  }
  states <- refit_args$missing_states
  if (is.null(states) || is.null(state_mapping)) {
    return(states)
  }
  mapped <- unname(state_mapping[as_state_labels(states)])
  mapped <- mapped[!is.na(mapped)]
  if (!length(mapped)) {
    stop(
      "A bootstrap resample contains none of the states in `missing_states`.",
      call. = FALSE
    )
  }
  mapped
}
