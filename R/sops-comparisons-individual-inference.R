# Draw-based inference for results with one row per prediction profile:
# comparisons from `comparisons()` and times from `state_time()`.

# "comparison" or "time", selecting the reduction in `profile_array_reducer()`.
profile_result_kind <- function(object) {
  if (inherits(object, "markov_state_time")) "time" else "comparison"
}

# The functions that return such a result, for error messages.
profile_result_source <- function(object) {
  if (identical(profile_result_kind(object), "time")) {
    "`state_time()` or `avg_state_time(by = FALSE)`"
  } else {
    "`comparisons()` or `avg_comparisons(by = FALSE)`"
  }
}

# The profile (patient) of each row: `rowidcf` when scenarios repeat the
# profiles, otherwise `rowid`.
profile_index <- function(data) {
  if ("rowidcf" %in% names(data)) data$rowidcf else data$rowid
}

inferences_profile_results <- function(
  object,
  method,
  n_draws,
  vcov,
  workers,
  conf_level,
  conf_type,
  return_draws,
  use_coefstart
) {
  if (identical(method, "bootstrap")) {
    mostr_abort(
      "Standard refit bootstrap is not supported for results with one row ",
      "per profile, from ",
      profile_result_source(object),
      ", because a resample can lack outcome states that the fixed profiles ",
      "need. Use `method = \"fwb\"`, or average over the profiles ",
      "(`by = TRUE` or grouping columns) for ordinary bootstrap inference."
    )
  }
  model <- attr(object, "model")
  if (is.null(model)) {
    mostr_abort("Model not stored in object. Cannot perform inference.")
  }
  validate_comparison_draw_cells(n_draws, nrow(object))
  setup <- setup_from_avg_comparison_object(object)
  reduce_draw <- profile_draw_reducer(object, setup)
  check_comparison_replay(object, setup, reduce_draw, model)

  draws <- if (identical(method, "fwb")) {
    comparison_fwb_draws(
      model = model,
      object = object,
      setup = setup,
      reduce_draw = reduce_draw,
      n_boot = n_draws,
      workers = workers,
      use_coefstart = use_coefstart,
      return_draws = return_draws
    )
  } else {
    comparison_simulation_draws(
      model = model,
      setup = setup,
      reduce_draw = reduce_draw,
      n_draws = n_draws,
      vcov = vcov,
      workers = workers
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
        draws$weights[, profile_index(keys), drop = FALSE]
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
    mostr_abort(
      "Inference for one row per profile would keep ",
      format(n_cells, scientific = FALSE),
      " draw values (draws x result rows), above the current limit of ",
      format(max_cells, scientific = FALSE),
      ". Lower `n_draws`; request fewer `times`, state sets, or profiles; ",
      "average over the profiles (`by = TRUE`); or increase option ",
      "`mostr.max_sops_draw_cells` if this allocation is intentional."
    )
  }
  invisible(NULL)
}

comparison_key_columns <- function(object) {
  intersect(
    c(
      "rowid",
      "rowidcf",
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
    "n_successful",
    "engine",
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
# [scenario x profile, time, state] to the estimates, in the row order of
# `object`.
profile_draw_reducer <- function(object, setup) {
  args <- attr(object, "comparison_args")
  reducer <- profile_array_reducer(profile_result_kind(object), args)
  # Interpolation needs each row's starting state only; the other profile
  # columns would be copied into every interpolated row.
  slim_setup <- setup
  slim_setup$newdata_pred <- setup$newdata_pred[,
    intersect(setup$settings$p_var, names(setup$newdata_pred)),
    drop = FALSE
  ]

  function(sops_array) {
    out <- if (reducer$real_time) {
      profile_results_from_individual_sops(
        x = time_benefit_sops_from_array(sops_array, slim_setup),
        setup = setup,
        reducer = reducer,
        time_map = args$time_map,
        baseline_time = args$baseline_time,
        target_times = args$target_times,
        posterior_summary = "mean",
        conf_level = 0.95,
        return_draws = FALSE
      )
    } else {
      reducer$reduce_array(
        stacked_to_counterfactual_array(
          sops_array,
          n_cf = setup$n_cf,
          n_each = setup$n_each
        ),
        setup = setup
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
    mostr_abort(
      "Stacked SOP rows are not aligned with the comparison scenarios."
    )
  }
  aperm(
    array(sops_array, dim = c(n_each, n_cf, dims[2L], dims[3L])),
    c(2L, 1L, 3L, 4L)
  )
}

# Draws are reduced without their row labels, so first check that reducing
# the fitted model's predictions reproduces the stored estimates.
check_comparison_replay <- function(object, setup, reduce_draw, model) {
  replayed <- reduce_draw(predict_markov_sops(
    model = model,
    newdata = setup$newdata_pred,
    times = setup$times,
    settings = setup$settings
  ))
  stored <- object$estimate
  # Each draw goes through the same reducer, which drops row labels, so check
  # that a replay at the fitted coefficients rebuilds the stored rows in order.
  # The values are identical when both came from the compiled code and differ by
  # about 1e-14 when the estimates came from the R engine; a larger gap means
  # edited estimates or rows the draws would not line up with.
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
    mostr_abort(
      "Replaying `comparisons()` did not reproduce the stored estimates."
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
        times = setup$times,
        settings = setup$settings
      ),
      error = function(e) {
        notify_sop_r_engine_fallback(conditionMessage(e))
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
    predict_markov_sops(
      model = set_coef(model, beta),
      newdata = setup$newdata_pred,
      times = setup$times,
      settings = setup$settings
    )
  }
}

comparison_simulation_draws <- function(
  model,
  setup,
  reduce_draw,
  n_draws,
  vcov,
  workers
) {
  beta_draws <- generate_sop_coefficient_draws(
    model = model,
    n_draws = n_draws,
    vcov = vcov
  )
  predict_array <- comparison_coefficient_predictor(model, setup)

  analysis_fn <- function(i) {
    sops_array <- tryCatch(
      predict_array(beta_draws[i, ]),
      error = function(e) {
        rethrow_crossed_thresholds(e)
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
    workers = workers
  )
  successful <- which(!vapply(results, is.null, logical(1)))
  if (length(successful) == 0L) {
    mostr_abort("All simulation draws failed.")
  }
  list(
    values = do.call(rbind, results[successful]),
    draw_ids = successful,
    metadata = list(n_successful = length(successful))
  )
}

comparison_fwb_draws <- function(
  model,
  object,
  setup,
  reduce_draw,
  n_boot,
  workers,
  use_coefstart,
  return_draws
) {
  refit_data <- stored_refit_bootstrap_data(object, model)
  validate_refit_bootstrap_data(
    refit_data,
    setup$settings$id_var,
    setup$settings$time_var
  )
  factor_cols <- bootstrap_state_columns(
    model,
    refit_data,
    setup$settings$p_var,
    setup$settings$p2_var
  )

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
    validate_prediction_weight_ids(setup$baseline_data, setup$settings$id_var)
    validate_draw_weight_column_available(setup$newdata_pred, weight_col)
  }

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = setup$settings$y_levels,
      absorb = setup$settings$absorb,
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
      predict_markov_sops(
        model = boot_result$model,
        newdata = setup$newdata_pred,
        times = setup$times,
        settings = markov_settings_update(
          setup$settings,
          y_levels = state_factor(boot_result$ylevels),
          absorb = boot_result$absorb
        )
      ),
      error = function(e) {
        rethrow_crossed_thresholds(e)
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
        fwb_baseline_weights(
          fwb_weights,
          setup$baseline_data,
          setup$settings$id_var
        )
      }
    )
  }

  fwb_samples <- generate_fwb_bootstrap_weights(
    data = refit_data,
    id_var = setup$settings$id_var,
    n_boot = n_boot
  )
  results <- apply_to_fwb_bootstrap(
    fwb_samples = fwb_samples,
    analysis_fn = analysis_fn,
    data = refit_data,
    id_var = setup$settings$id_var,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "stats")
  )
  results <- keep_successful_draws(results)
  list(
    values = do.call(rbind, lapply(results, `[[`, "estimate")),
    draw_ids = seq_along(results),
    weights = if (!is.null(weight_col)) {
      do.call(rbind, lapply(results, `[[`, "weight"))
    },
    metadata = list(
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
