# Draw-based inference for SOP comparisons of each prediction profile.

inferences_comparisons <- function(
  object,
  method,
  n_draws,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  update_datadist,
  use_coefstart
) {
  if (identical(method, "bootstrap")) {
    stop(
      "Standard refit bootstrap is not supported for comparisons of each ",
      "profile, from `comparisons()` or `avg_comparisons(by = FALSE)`, ",
      "because a resample can lack outcome states that the fixed profiles ",
      "need. Use `method = \"fwb\"`, or average the comparisons (`by = TRUE` ",
      "or grouping columns) for ordinary bootstrap inference.",
      call. = FALSE
    )
  }
  model <- attr(object, "model")
  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform inference.")
  }
  validate_comparison_draw_cells(n_draws, nrow(object))
  setup <- setup_from_avg_comparison_object(object)
  reduce_draw <- comparison_draw_reducer(object, setup)
  check_comparison_replay(object, setup, reduce_draw, model)

  draws <- if (identical(method, "fwb")) {
    comparison_fwb_draws(
      model = model,
      object = object,
      setup = setup,
      reduce_draw = reduce_draw,
      n_boot = n_draws,
      workers = workers,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart,
      return_draws = return_draws
    )
  } else {
    comparison_simulation_draws(
      model = model,
      setup = setup,
      reduce_draw = reduce_draw,
      engine = method,
      n_draws = n_draws,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      return_draws = return_draws
    )
  }

  summary_stats <- summarize_sop_draw_matrix(
    draw_values = draws$values,
    conf_level = conf_level,
    conf_type = conf_type,
    point_estimates = object$estimate
  )
  result <- clear_comparison_inference(object)
  result$conf.low <- summary_stats$conf.low
  result$conf.high <- summary_stats$conf.high
  result$std.error <- summary_stats$std.error
  for (name in names(draws$metadata)) {
    attr(result, name) <- draws$metadata[[name]]
  }
  if (isTRUE(return_draws)) {
    keys <- as.data.frame(object)[,
      comparison_key_columns(object),
      drop = FALSE
    ]
    draws_df <- posterior_draw_matrix_to_df(
      draws$values,
      keys,
      draws$draw_ids
    )
    # Each row gets its patient's weight in that draw, in the same
    # draw-by-row order as the estimates.
    if (isTRUE(draws$metadata$draw_weights_attached)) {
      draws_df[[draws$metadata$draw_weight_col]] <- as.vector(t(
        draws$weights[, keys$rowid, drop = FALSE]
      ))
    }
    attr(result, "draws") <- draws_df
  }
  result
}

# The draw values are held as a draws x result-rows matrix, and returned
# draws expand it further, so the size is checked before any draw is run.
validate_comparison_draw_cells <- function(n_draws, n_rows) {
  max_cells <- getOption("mostr.max_sops_draw_cells", 5e7)
  n_cells <- as.double(n_draws) * n_rows
  if (n_cells > max_cells) {
    stop(
      "Inference for `comparisons()` would keep ",
      format(n_cells, scientific = FALSE),
      " draw values (draws x result rows), above the current limit of ",
      format(max_cells, scientific = FALSE),
      ". Lower `n_draws`; request fewer `times`, state sets, or profiles; ",
      "use `avg_comparisons()` for averaged comparisons; or increase option ",
      "`mostr.max_sops_draw_cells` if this allocation is intentional.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

comparison_key_columns <- function(object) {
  intersect(
    c(
      "rowid",
      "estimand",
      "time",
      "state_set",
      "term",
      "reference_level",
      "comparison_level",
      "contrast",
      "comparison",
      "time_unit"
    ),
    names(object)
  )
}

clear_comparison_inference <- function(object) {
  result <- object
  for (col in setdiff(inference_columns(), "estimate")) {
    result[[col]] <- NULL
  }
  for (name in c(
    "draws",
    "null",
    "analytical",
    "target",
    "covariance_source",
    "covariance_metadata",
    "n_boot",
    "n_successful",
    "engine",
    "score_weight_dist",
    "fwb_weight_type",
    "fwb_weight_scale",
    "draw_weights_attached",
    "draw_weight_col",
    "draw_weight_omission_reason",
    "baseline_anchor_draws"
  )) {
    attr(result, name) <- NULL
  }
  result
}

# Returns a function that reduces one stacked SOP array
# [scenario x profile, time, state] to the comparison estimates, in the row
# order of `object`.
comparison_draw_reducer <- function(object, setup) {
  args <- attr(object, "comparison_args")
  real_time <- !is.null(args$time_map) && !identical(args$estimand, "sop")
  # Interpolation needs each row's starting state only; the other profile
  # columns would be copied into every interpolated row.
  slim_setup <- setup
  slim_setup$newdata_pred <- setup$newdata_pred[,
    intersect(setup$p_var, names(setup$newdata_pred)),
    drop = FALSE
  ]

  function(sops_array) {
    out <- if (real_time) {
      comparisons_from_individual_sops(
        x = time_benefit_sops_from_array(sops_array, slim_setup),
        setup = setup,
        estimand = args$estimand,
        state_sets = args$state_sets,
        comparison = args$comparison,
        time_map = args$time_map,
        baseline_time = args$baseline_time,
        target_times = args$target_times,
        time_unit = args$time_unit,
        posterior_summary = "mean",
        conf_level = 0.95,
        return_draws = FALSE
      )
    } else {
      individual_comparisons_from_array(
        sops_array = stacked_to_counterfactual_array(
          sops_array,
          n_cf = setup$n_cf,
          n_each = setup$n_each
        ),
        setup = setup,
        estimand = args$estimand,
        state_sets = args$state_sets,
        comparison = args$comparison,
        real_time = FALSE,
        time_unit = args$time_unit
      )
    }
    out$estimate
  }
}

# [scenario x profile, time, state], with scenarios in consecutive blocks of
# `n_each` profiles -> [scenario, profile, time, state].
stacked_to_counterfactual_array <- function(sops_array, n_cf, n_each) {
  dims <- dim(sops_array)
  if (dims[1L] != n_cf * n_each) {
    stop("Stacked SOP rows are not aligned with the comparison scenarios.")
  }
  aperm(
    array(sops_array, dim = c(n_each, n_cf, dims[2L], dims[3L])),
    c(2L, 1L, 3L, 4L)
  )
}

# Draws are reduced without their row labels, so first check that reducing
# the fitted model's predictions reproduces the stored estimates.
check_comparison_replay <- function(object, setup, reduce_draw, model) {
  replayed <- reduce_draw(soprob_markov(
    model = model,
    newdata = setup$newdata_pred,
    times = setup$times,
    y_levels = setup$y_levels,
    absorb = setup$absorb,
    time_var = setup$time_var,
    p_var = setup$p_var,
    p2_var = setup$p2_var,
    gap_var = setup$gap_var,
    time_covariates = setup$time_covariates
  ))
  stored <- object$estimate
  same_missing <- identical(is.na(replayed), is.na(stored))
  ok <- !is.na(stored)
  if (
    length(replayed) != length(stored) ||
      !same_missing ||
      any(
        abs(replayed[ok] - stored[ok]) >
          1e-12 + 1e-10 * pmax(abs(replayed[ok]), abs(stored[ok]))
      )
  ) {
    stop(
      "Replaying `comparisons()` did not reproduce the stored estimates.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Returns a function from a coefficient vector to the stacked SOP array.
comparison_coefficient_predictor <- function(model, setup) {
  model_plan <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  if (inherits(model_plan, c("vglm", "orm")) && !inherits(model_plan, "blrm")) {
    plan <- tryCatch(
      compile_sop_execution_plan(
        model = model_plan,
        newdata = setup$newdata_pred,
        time_covariates = setup$time_covariates,
        times = setup$times,
        y_levels = setup$y_levels,
        absorb = setup$absorb,
        time_var = setup$time_var,
        p_var = setup$p_var,
        p2_var = setup$p2_var,
        gap_var = setup$gap_var,
        builder = "streamed",
        output = "individual"
      ),
      error = function(e) {
        notify_sop_reference_fallback(conditionMessage(e))
        NULL
      }
    )
    if (!is.null(plan)) {
      return(function(beta) {
        run_sop_execution_plan(
          plan,
          get_effective_coefs(model_plan, beta = beta)
        )
      })
    }
  }
  function(beta) {
    soprob_markov(
      model = set_coef(model, beta),
      newdata = setup$newdata_pred,
      times = setup$times,
      y_levels = setup$y_levels,
      absorb = setup$absorb,
      time_var = setup$time_var,
      p_var = setup$p_var,
      p2_var = setup$p2_var,
      gap_var = setup$gap_var,
      time_covariates = setup$time_covariates
    )
  }
}

comparison_simulation_draws <- function(
  model,
  setup,
  reduce_draw,
  engine,
  n_draws,
  vcov,
  cluster,
  workers,
  return_draws
) {
  # Score-bootstrap patient weights are not used here, because nothing is
  # averaged, but they are returned with the draws as `score_weight`, as for
  # `sops()`, so the profiles' comparisons can be averaged afterwards.
  simulation_draws <- generate_sop_coefficient_draws(
    model = model,
    engine = engine,
    score_weight_dist = "exponential",
    n_draws = n_draws,
    vcov = vcov,
    cluster = cluster,
    baseline_data = setup$baseline_data,
    id_var = setup$id_var,
    newdata_supplied = setup$newdata_supplied,
    is_avg = FALSE,
    by = NULL,
    return_draws = return_draws,
    prediction_data = setup$newdata_pred
  )
  beta_draws <- simulation_draws$beta_draws
  predict_array <- comparison_coefficient_predictor(model, setup)

  analysis_fn <- function(i) {
    sops_array <- tryCatch(
      predict_array(beta_draws[i, ]),
      error = function(e) {
        warning(
          "SOP prediction failed in draw ",
          i,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }
    reduce_draw(sops_array)
  }
  results <- apply_sop_simulation_draws(
    n_draws = nrow(beta_draws),
    analysis_fn = analysis_fn,
    workers = workers,
    globals = c("beta_draws", "predict_array", "reduce_draw")
  )
  successful <- which(!vapply(results, is.null, logical(1)))
  if (length(successful) == 0L) {
    stop("All simulation draws failed.")
  }
  weights <- simulation_draws$baseline_weights_draws
  weight_col <- simulation_draws$draw_weight_col
  list(
    values = do.call(rbind, results[successful]),
    draw_ids = successful,
    weights = if (!is.null(weight_col) && !is.null(weights)) {
      weights[successful, , drop = FALSE]
    },
    metadata = c(
      list(n_successful = length(successful)),
      if (identical(engine, "score_bootstrap")) {
        list(
          draw_weights_attached = simulation_draws$draw_weights_attached,
          draw_weight_col = weight_col,
          draw_weight_omission_reason =
            simulation_draws$draw_weight_omission_reason
        )
      }
    )
  )
}

comparison_fwb_draws <- function(
  model,
  object,
  setup,
  reduce_draw,
  n_boot,
  workers,
  update_datadist,
  use_coefstart,
  return_draws
) {
  refit_data <- stored_refit_bootstrap_data(object, model)
  validate_refit_bootstrap_data(refit_data, setup$id_var, setup$time_var)
  factor_cols <- intersect(c("y", setup$p_var, setup$p2_var), names(refit_data))

  # The patient weights enter the comparisons only through the refits, but,
  # as for `sops()`, the stored patients' weights are returned with the draws
  # as `fwb_weight` so the comparisons can be averaged afterwards. Supplied
  # profiles cannot be matched to the bootstrap patients.
  weight_col <- if (setup$newdata_supplied) NULL else "fwb_weight"
  if (is.null(weight_col)) {
    if (isTRUE(return_draws)) {
      warn_fixed_profile_draw_weights("FWB")
    }
  } else {
    validate_prediction_weight_ids(setup$baseline_data, setup$id_var)
    validate_draw_weight_column_available(setup$newdata_pred, weight_col)
  }

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = setup$y_levels,
      absorb = setup$absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart,
      fit_weights = boot_data$fwb_weight
    )
    if (is.null(boot_result$model)) {
      return(NULL)
    }
    if (length(boot_result$missing_states) > 0L) {
      warning(
        "Fractional weighted bootstrap unexpectedly dropped outcome-state ",
        "support; skipping this draw.",
        call. = FALSE
      )
      return(NULL)
    }
    sops_array <- tryCatch(
      soprob_markov(
        model = boot_result$model,
        newdata = setup$newdata_pred,
        times = setup$times,
        y_levels = factor(boot_result$ylevels %||% boot_result$y_levels),
        absorb = boot_result$absorb,
        time_var = setup$time_var,
        p_var = setup$p_var,
        p2_var = setup$p2_var,
        gap_var = setup$gap_var,
        time_covariates = setup$time_covariates
      ),
      error = function(e) {
        warning("soprob_markov failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }
    list(
      estimate = reduce_draw(sops_array),
      weight = if (!is.null(weight_col)) {
        fwb_baseline_weights(fwb_weights, setup$baseline_data, setup$id_var)
      }
    )
  }

  fwb_samples <- generate_fwb_bootstrap_weights(
    data = refit_data,
    id_var = setup$id_var,
    n_boot = n_boot
  )
  results <- apply_to_fwb_bootstrap(
    fwb_samples = fwb_samples,
    analysis_fn = analysis_fn,
    data = refit_data,
    id_var = setup$id_var,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "stats"),
    globals = c(
      "model",
      "setup",
      "refit_data",
      "factor_cols",
      "update_datadist",
      "use_coefstart",
      "reduce_draw",
      "weight_col"
    )
  )
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    stop("All bootstrap iterations failed.")
  }
  list(
    values = do.call(rbind, lapply(results, `[[`, "estimate")),
    draw_ids = seq_along(results),
    weights = if (!is.null(weight_col)) {
      do.call(rbind, lapply(results, `[[`, "weight"))
    },
    metadata = list(
      n_boot = n_boot,
      n_successful = length(results),
      draw_weights_attached = !is.null(weight_col) && isTRUE(return_draws),
      draw_weight_col = weight_col,
      draw_weight_omission_reason = if (is.null(weight_col)) {
        "user_supplied_newdata"
      },
      fwb_weight_type = "exponential",
      fwb_weight_scale = "cluster_mean_1"
    )
  )
}
