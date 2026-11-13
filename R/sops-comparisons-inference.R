# Average comparison inference helpers.

# `method` is the `method` of `inferences()`: "mvn", "bootstrap" or "fwb".
inferences_avg_comparisons <- function(
  object,
  method,
  n_draws,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  use_coefstart
) {
  args <- attr(object, "comparison_args")
  estimand <- args$estimand

  if (estimand %in% c("sop", "time_in_state")) {
    return(inferences_avg_comparisons_linear(
      object = object,
      method = method,
      n_draws = n_draws,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      use_coefstart = use_coefstart
    ))
  }

  if (method == "mvn") {
    return(inferences_avg_comparisons_time_benefit_simulation(
      object = object,
      n_draws = n_draws,
      vcov = vcov,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws
    ))
  }

  inferences_avg_comparisons_time_benefit_bootstrap(
    object = object,
    method = method,
    n_boot = n_draws,
    workers = workers,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = return_draws,
    use_coefstart = use_coefstart
  )
}

inferences_avg_comparisons_linear <- function(
  object,
  method,
  n_draws,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  use_coefstart
) {
  args <- attr(object, "comparison_args")
  avg_args <- attr(object, "avg_args")
  newdata <- if (isTRUE(attr(object, "newdata_supplied"))) {
    attr(object, "newdata_orig")
  } else {
    NULL
  }
  avg <- avg_comparison_replay_avg_sops(
    model = attr(object, "model"),
    newdata = newdata,
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    settings = markov_result_settings(object),
    include_re = args$include_re,
    n_draws = args$n_draws,
    seed = args$seed,
    posterior_summary = args$posterior_summary,
    conf_level = conf_level,
    return_draws = FALSE,
    extra_args = args$extra_args,
    warn_skipped_visits = FALSE
  )
  avg <- inferences_internal(
    avg,
    method = method,
    n_draws = n_draws,
    vcov = vcov,
    cluster = cluster,
    workers = workers,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = TRUE,
    use_coefstart = use_coefstart
  )

  state_sets <- normalize_comparison_state_sets(
    args$state_sets,
    markov_result_settings(avg)$y_levels
  )
  result <- avg_comparison_from_avg_sops(
    avg,
    estimand = args$estimand,
    state_sets = state_sets,
    comparison = args$comparison,
    time_map = args$time_map,
    baseline_time = args$baseline_time,
    target_times = args$target_times,
    time_unit = args$time_unit,
    return_draws = return_draws
  )
  result <- restore_avg_comparison_inference_attrs(result, object, avg)
  if (!isTRUE(return_draws)) {
    attr(result, "draws") <- NULL
  }
  result
}

restore_avg_comparison_inference_attrs <- function(result, object, inferred) {
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  for (a in c(
    "n_draws",
    "n_successful",
    "conf_level",
    "conf_type",
    "method",
    "engine",
    "fwb_weight_type",
    "fwb_weight_scale",
    "draw_weights_attached",
    "draw_weight_col",
    "draw_weight_omission_reason"
  )) {
    if (!is.null(attr(inferred, a))) {
      attr(result, a) <- attr(inferred, a)
    }
  }
  class(result) <- class(object)
  result
}

inferences_avg_comparisons_time_benefit_simulation <- function(
  object,
  n_draws,
  vcov,
  workers,
  conf_level,
  conf_type,
  return_draws
) {
  setup <- setup_from_avg_comparison_object(object)
  args <- attr(object, "comparison_args")
  model <- attr(object, "model")

  beta_draws <- generate_sop_coefficient_draws(
    model = model,
    n_draws = n_draws,
    vcov = vcov
  )

  analysis_fn <- function(i) {
    model_i <- set_coef(model, beta_draws[i, ])
    sops_array <- tryCatch(
      predict_markov_sops(
        model = model_i,
        newdata = setup$newdata_pred,
        times = setup$times,
        settings = setup$settings
      ),
      error = function(e) {
        rethrow_crossed_thresholds(e)
        warning("soprob_markov failed in draw ", i, ": ", e$message)
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }
    draw_i <- reduce_time_benefit_array_for_setup(
      sops_array = sops_array,
      setup = setup,
      comparison = args$comparison,
      weights = NULL,
      time_unit = args$time_unit,
      time_map = args$time_map,
      baseline_time = args$baseline_time,
      target_times = args$target_times
    )
    draw_i$draw_id <- i
    draw_i
  }

  draw_results <- apply_sop_simulation_draws(
    n_draws = nrow(beta_draws),
    analysis_fn = analysis_fn,
    workers = workers
  )

  draw_results <- Filter(Negate(is.null), draw_results)
  if (length(draw_results) == 0L) {
    mostr_abort("All simulation draws failed.")
  }
  draws_df <- bind_rows_fill(draw_results)
  finalize_avg_comparison_draw_inference(
    object = object,
    draws_df = draws_df,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = return_draws,
    metadata = list(
      n_draws = n_draws,
      n_successful = length(draw_results),
      conf_level = conf_level,
      conf_type = conf_type,
      method = "mvn",
      engine = "mvn"
    )
  )
}

# `method` is "bootstrap", which resamples patients, or "fwb", which reweights
# them.
inferences_avg_comparisons_time_benefit_bootstrap <- function(
  object,
  method,
  n_boot,
  workers,
  conf_level,
  conf_type,
  return_draws,
  use_coefstart
) {
  method <- match.arg(method, choices = c("bootstrap", "fwb"))
  setup <- setup_from_avg_comparison_object(object)
  args <- attr(object, "comparison_args")
  model <- attr(object, "model")
  refit_data <- stored_refit_bootstrap_data(object, model)
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))

  validate_refit_bootstrap_data(
    refit_data,
    setup$settings$id_var,
    setup$settings$time_var
  )
  if (newdata_supplied && method == "fwb") {
    warn_fixed_profile_bootstrap_weights("FWB")
  }

  factor_cols <- bootstrap_state_columns(
    model,
    refit_data,
    setup$settings$p_var,
    setup$settings$p2_var
  )

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    fit_weights <- if (method == "fwb") boot_data$fwb_weight else NULL
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = setup$settings$y_levels,
      absorb = setup$settings$absorb,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )
    m_boot <- boot_result$model
    if (is.null(m_boot)) {
      return(NULL)
    }

    if (newdata_supplied) {
      newdata_cf <- setup$newdata_pred
      baseline_data <- setup$baseline_data
      baseline_weights <- NULL
    } else if (method == "bootstrap") {
      baseline_data <- resolve_markov_prediction_data(
        boot_result$data,
        id_var = "new_id",
        time_var = setup$settings$time_var,
        data_label = "bootstrap data"
      )
      baseline_data <- ensure_markov_rowid(baseline_data)
      newdata_cf <- create_counterfactual_data(
        baseline_data,
        setup$grid,
        setup$variables
      )
      baseline_weights <- NULL
    } else {
      baseline_data <- resolve_markov_prediction_data(
        boot_result$data,
        id_var = setup$settings$id_var,
        time_var = setup$settings$time_var,
        data_label = "bootstrap data"
      )
      baseline_data <- ensure_markov_rowid(baseline_data)
      newdata_cf <- create_counterfactual_data(
        baseline_data,
        setup$grid,
        setup$variables
      )
      baseline_weights <- fwb_baseline_weights(
        fwb_weights = fwb_weights,
        baseline_data = baseline_data,
        id_var = setup$settings$id_var
      )
    }

    sops_array <- tryCatch(
      predict_markov_sops(
        model = m_boot,
        newdata = newdata_cf,
        times = setup$times,
        settings = markov_settings_update(
          setup$settings,
          y_levels = state_factor(boot_result$ylevels),
          absorb = boot_result$absorb
        )
      ),
      error = function(e) {
        rethrow_crossed_thresholds(e)
        warning("soprob_markov failed: ", e$message)
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }

    # The refit keeps the original state labels, so its states are matched
    # to the original ones by label; a missing state gets probability 0.
    sops_array <- complete_comparison_sops_array_states(
      sops_array = sops_array,
      boot_ylevels = boot_result$ylevels,
      y_levels = setup$settings$y_levels
    )
    boot_setup <- setup
    boot_setup$baseline_data <- baseline_data
    boot_setup$newdata_pred <- newdata_cf
    boot_setup$n_each <- nrow(baseline_data)
    reduce_time_benefit_array_for_setup(
      sops_array = sops_array,
      setup = boot_setup,
      comparison = args$comparison,
      weights = baseline_weights,
      time_unit = args$time_unit,
      time_map = args$time_map,
      baseline_time = args$baseline_time,
      target_times = args$target_times
    )
  }

  if (method == "bootstrap") {
    boot_ids <- fast_group_bootstrap(
      data = refit_data,
      id_var = setup$settings$id_var,
      n_boot = n_boot
    )
    # Supplied profiles keep their starting states in every refit, so each
    # resample must hold those states.
    check_bootstrap_starting_states(
      boot_ids,
      refit_data = refit_data,
      id_var = setup$settings$id_var,
      state_cols = factor_cols,
      starting_states = bootstrap_starting_states(
        profiles = if (newdata_supplied) setup$newdata_pred,
        state_vars = c(setup$settings$p_var, setup$settings$p2_var)
      )
    )
    # Supplied profiles are fixed and always contain every group.
    if (!newdata_supplied) {
      check_bootstrap_by_groups(
        boot_ids,
        baseline_data = setup$baseline_data,
        id_var = setup$settings$id_var,
        by = setup$by
      )
    }
    boot_results <- apply_to_bootstrap(
      boot_samples = boot_ids,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = setup$settings$id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats")
    )
  } else {
    fwb_samples <- generate_fwb_bootstrap_weights(
      data = refit_data,
      id_var = setup$settings$id_var,
      n_boot = n_boot
    )
    boot_results <- apply_to_fwb_bootstrap(
      fwb_samples = fwb_samples,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = setup$settings$id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats")
    )
  }

  boot_results <- keep_successful_draws(boot_results)
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  finalize_avg_comparison_draw_inference(
    object = object,
    draws_df = boot_df,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = return_draws,
    metadata = list(
      n_boot = n_boot,
      n_successful = length(boot_results),
      conf_level = conf_level,
      method = "bootstrap",
      # The stored name of the resampling scheme.
      engine = if (method == "fwb") "fwb" else "standard",
      fwb_weight_type = if (method == "fwb") "exponential" else NULL,
      fwb_weight_scale = if (method == "fwb") "cluster_mean_1" else NULL
    )
  )
}

setup_from_avg_comparison_object <- function(object) {
  avg_args <- attr(object, "avg_args")
  call_args <- attr(object, "call_args")
  list(
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    newdata_orig = attr(object, "newdata_orig"),
    refit_data = attr(object, "refit_data"),
    newdata_supplied = isTRUE(attr(object, "newdata_supplied")),
    baseline_data = attr(object, "comparison_baseline_data"),
    grid = attr(object, "comparison_grid"),
    newdata_pred = attr(object, "newdata_pred"),
    n_cf = nrow(attr(object, "comparison_grid")),
    n_each = attr(object, "comparison_n_each"),
    call_args = call_args,
    settings = markov_result_settings(object)
  )
}

# The averaged time benefit of one stacked SOP array
# [scenario x profile, time, state]: each profile's time benefit, as in
# `comparisons()`, averaged over the profiles of each `setup$by` group.
reduce_time_benefit_array_for_setup <- function(
  sops_array,
  setup,
  comparison,
  weights,
  time_unit,
  time_map = NULL,
  baseline_time = 0,
  target_times = NULL
) {
  args <- list(
    estimand = "time_benefit",
    comparison = comparison,
    time_map = time_map,
    baseline_time = baseline_time,
    target_times = target_times
  )
  out <- average_time_benefit(
    profile_rows_from_stacked_sops(
      sops_array,
      setup = setup,
      reducer = profile_array_reducer("comparison", args),
      args = args
    ),
    setup = setup,
    weights = weights
  )
  if (!is.null(time_unit)) {
    out$time_unit <- time_unit
  }
  out
}

time_benefit_sops_from_array <- function(sops_array, setup) {
  newdata <- setup$newdata_pred
  if (is.null(newdata)) {
    mostr_abort(
      "`setup$newdata_pred` is required for real-time time-benefit inference."
    )
  }
  newdata$rowid <- seq_len(nrow(newdata))

  out <- array_to_df_individual(
    sops_array = sops_array,
    times = setup$times,
    y_levels = setup$settings$y_levels,
    newdata = newdata
  )
  attr(out, "settings") <- setup$settings
  attr(out, "call_args") <- setup$call_args %||% sops_call_args(setup$times)
  attr(out, "newdata_orig") <- setup$newdata_orig %||% setup$baseline_data
  attr(out, "newdata_pred") <- newdata
  attr(out, "newdata_supplied") <- setup$newdata_supplied
  class(out) <- c("markov_sops", class(out))
  out
}

complete_comparison_sops_array_states <- function(
  sops_array,
  boot_ylevels,
  y_levels
) {
  boot_labels <- as_state_labels(boot_ylevels)
  target_labels <- as_state_labels(y_levels)
  if (identical(boot_labels, target_labels)) {
    return(sops_array)
  }

  full <- array(
    0,
    dim = c(dim(sops_array)[1], dim(sops_array)[2], length(target_labels))
  )
  for (i in seq_along(boot_labels)) {
    target <- match(boot_labels[i], target_labels)
    if (!is.na(target)) {
      full[,, target] <- sops_array[,, i]
    }
  }
  full
}

finalize_avg_comparison_draw_inference <- function(
  object,
  draws_df,
  conf_level,
  conf_type,
  return_draws,
  metadata
) {
  group_cols <- setdiff(names(draws_df), c("estimate", "draw_id"))
  ci <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type,
    point_estimates = as.data.frame(object)
  )
  result <- merge_comparison_ci(as.data.frame(object), ci, group_cols)
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  for (nm in names(metadata)) {
    if (!is.null(metadata[[nm]])) {
      attr(result, nm) <- metadata[[nm]]
    }
  }
  if (isTRUE(return_draws)) {
    attr(result, "draws") <- draws_df
  }
  class(result) <- class(object)
  result
}
