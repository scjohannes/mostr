# SOP refit-bootstrap inference.

# =============================================================================
# BOOTSTRAP-BASED INFERENCE
# =============================================================================

#' Bootstrap-Based Inference for SOPs
#'
#' Internal function that implements bootstrap inference for SOPs. Resamples
#' patients with replacement, refits the model, and computes SOPs for each
#' bootstrap sample.
#'
#' @param object A `markov_avg_sops` object, or a `markov_sops` object when
#'   `method = "fwb"`.
#' @param method `"bootstrap"`, which resamples patients, or `"fwb"`, which
#'   reweights them.
#' @param n_boot Number of bootstrap iterations.
#' @param workers Number of parallel workers. If NULL or 1, uses sequential
#'   processing. If > 1, uses parallel processing.
#' @param conf_level Confidence level.
#' @param return_draws Store individual bootstrap draws?
#' @param use_coefstart Use original coefficients as starting values?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_bootstrap <- function(
  object,
  method,
  n_boot,
  workers,
  conf_level,
  conf_type,
  return_draws,
  use_coefstart
) {
  if (!inherits(object, "markov_avg_sops")) {
    if (method != "fwb") {
      mostr_abort(
        "Standard refit bootstrap is not supported for individual-level ",
        "SOPs because ordinary resampling can drop outcome-state support ",
        "needed by fixed prediction rows. Use `method = \"fwb\"`, or call ",
        "`avg_sops()` for marginal bootstrap inference."
      )
    }
    return(inferences_bootstrap_sops_fwb(
      object = object,
      n_boot = n_boot,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      use_coefstart = use_coefstart
    ))
  }

  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  prediction_data <- attr(object, "newdata_pred")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))
  avg_args <- attr(object, "avg_args")

  settings <- markov_result_settings(object)
  y_levels <- settings$y_levels

  variables <- avg_args$variables
  by <- avg_args$by
  times <- avg_args$times
  id_var <- settings$id_var

  refit_data <- stored_refit_bootstrap_data(object, model)
  validate_refit_bootstrap_data(refit_data, id_var, settings$time_var)

  if (newdata_supplied && method == "fwb") {
    warn_fixed_profile_bootstrap_weights("FWB")
  }

  # Dimensions for result expansion
  n_times <- length(times)
  n_states <- length(y_levels)

  # --- 2. Identify State Columns ---
  factor_cols <- bootstrap_state_columns(
    model,
    refit_data,
    settings$p_var,
    settings$p2_var
  )

  # --- 3. Define Analysis Function ---
  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    baseline_anchor <- NULL
    if (method == "bootstrap") {
      # A patient sampled more than once represents distinct bootstrap
      # patients, so the ID column takes the unique bootstrap identifier
      # instead of the repeated source identifier.
      boot_data[[id_var]] <- boot_data$new_id
    }

    fit_weights <- NULL
    if (method == "fwb") {
      fit_weights <- boot_data$fwb_weight
    }

    # A. Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = y_levels,
      absorb = settings$absorb,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }

    # Get states present in original numbering (for mapping back later)
    ylevel_names <- as_state_labels(y_levels)
    states_present <- ylevel_names[
      !ylevel_names %in% as_state_labels(missing_states)
    ]

    grid <- create_counterfactual_grid(variables)
    if (newdata_supplied) {
      newdata_cf <- prediction_data
      n_cf <- nrow(grid)
      if (is.null(newdata_cf) || nrow(newdata_cf) %% n_cf != 0L) {
        mostr_abort(
          "Stored prediction data is not aligned with counterfactual grid."
        )
      }
      baseline_weights <- NULL
    } else {
      # B. Compute standardized SOPs using G-computation on bootstrap data
      # Extract baseline data (one row per patient)
      if (method == "bootstrap") {
        baseline_boot <- resolve_markov_prediction_data(
          boot_data,
          id_var = "new_id",
          time_var = settings$time_var,
          data_label = "bootstrap data"
        )
        baseline_weights <- NULL
      } else {
        baseline_boot <- resolve_markov_prediction_data(
          boot_data,
          id_var = id_var,
          time_var = settings$time_var,
          data_label = "bootstrap data"
        )
        baseline_weights <- fwb_baseline_weights(
          fwb_weights = fwb_weights,
          baseline_data = baseline_boot,
          id_var = id_var
        )
      }

      baseline_anchor <- empirical_baseline_anchor_from_data(
        x = object,
        baseline = baseline_boot,
        baseline_time = 0,
        weights = baseline_weights
      )

      # Create counterfactual datasets for each variable value
      newdata_cf <- create_counterfactual_data(
        baseline_boot,
        grid,
        variables,
        by = by
      )
    }

    # Compute individual SOPs for counterfactual data. The refit keeps the
    # original state labels, so supplied profiles are predicted as they are;
    # `check_bootstrap_starting_states()` has made sure that their starting
    # states are in every resample.
    sops_array <- tryCatch(
      predict_markov_sops(
        model = m_boot,
        newdata = newdata_cf,
        times = times,
        settings = markov_settings_update(
          settings,
          y_levels = state_factor(boot_ylevels),
          absorb = boot_absorb
        )
      ),
      error = function(e) {
        rethrow_crossed_thresholds(e)
        warning("soprob_markov failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(sops_array)) {
      return(NULL)
    }

    # C. Marginalize (average) across patients for each counterfactual
    # sops_array is [n_pat, n_times, n_boot_states]
    n_cf <- nrow(grid)
    n_each <- nrow(newdata_cf) / n_cf
    boot_avg <- marginalize_sops_array(
      sops_array = sops_array,
      grid = grid,
      times = times,
      y_levels = state_factor(boot_ylevels),
      variables = variables,
      n_cf = n_cf,
      n_each = n_each,
      weights = baseline_weights,
      by = by,
      newdata = newdata_cf
    )

    if (!is.null(by)) {
      if (length(missing_states) > 0) {
        boot_avg <- complete_bootstrap_sop_states(
          boot_avg,
          times = times,
          y_levels = y_levels,
          variables = variables,
          by = by
        )
      }
      attr(boot_avg, "baseline_anchor") <- baseline_anchor
      return(boot_avg)
    }

    avg_sops_list <- bootstrap_avg_df_to_matrices(
      boot_avg = boot_avg,
      grid = grid,
      times = times,
      states = state_factor(boot_ylevels),
      variables = variables
    )

    # D. Expand back to original state space with zero-padding: the refit's
    # states are the present states, in state order and with their labels.
    if (length(missing_states) > 0) {
      for (cf_i in seq_len(n_cf)) {
        avg_sops_mat <- avg_sops_list[[cf_i]]
        full_mat <- matrix(0, nrow = n_times, ncol = n_states)

        # Fill in states that were present
        for (i in seq_along(states_present)) {
          original_state <- states_present[i]
          original_idx <- which(ylevel_names == as_state_labels(original_state))
          full_mat[, original_idx] <- avg_sops_mat[, i]
        }

        avg_sops_list[[cf_i]] <- full_mat
      }
    }

    # E. Format as data frame
    result_list <- vector("list", n_cf)
    for (cf_i in seq_len(n_cf)) {
      avg_sops_mat <- avg_sops_list[[cf_i]]

      df <- expand.grid(time = times, state = y_levels)
      df$estimate <- as.vector(avg_sops_mat)

      # Add variable values for this counterfactual
      for (v in names(grid)) {
        df[[v]] <- grid[cf_i, v]
      }

      result_list[[cf_i]] <- df
    }

    result <- bind_rows_fill(result_list)
    attr(result, "baseline_anchor") <- baseline_anchor
    result
  }

  # --- 4. Apply to Bootstrap Samples ---
  if (method == "bootstrap") {
    boot_ids <- fast_group_bootstrap(
      data = refit_data,
      id_var = id_var,
      n_boot = n_boot
    )
    # Supplied profiles keep their starting states in every refit, so each
    # resample must hold those states.
    check_bootstrap_starting_states(
      boot_ids,
      refit_data = refit_data,
      id_var = id_var,
      state_cols = factor_cols,
      starting_states = bootstrap_starting_states(
        profiles = if (newdata_supplied) prediction_data,
        state_vars = c(settings$p_var, settings$p2_var)
      )
    )
    # Supplied profiles are fixed and always contain every group. Otherwise
    # each resample's profiles come from its own rows, as in `analysis_fn`.
    if (!newdata_supplied && !is.null(by)) {
      check_bootstrap_by_groups(
        boot_ids,
        baseline_data = add_observed_by_columns(
          resolve_markov_prediction_data(
            refit_data,
            id_var = id_var,
            time_var = settings$time_var,
            data_label = "refit data"
          ),
          variables,
          by
        ),
        id_var = id_var,
        by = by
      )
    }

    boot_results <- apply_to_bootstrap(
      boot_samples = boot_ids,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats")
    )
  } else {
    fwb_samples <- generate_fwb_bootstrap_weights(
      data = refit_data,
      id_var = id_var,
      n_boot = n_boot
    )

    boot_results <- apply_to_fwb_bootstrap(
      fwb_samples = fwb_samples,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats")
    )
  }

  # --- 5. Combine and Compute Summary Statistics ---
  boot_results <- keep_successful_draws(boot_results)
  baseline_anchor_draws <- combine_baseline_anchor_draws(
    lapply(boot_results, attr, "baseline_anchor"),
    seq_along(boot_results)
  )

  # Add draw_id and combine
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  # Compute confidence intervals
  group_cols <- c("time", "state", names(variables))
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = boot_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type,
    point_estimates = as.data.frame(object)
  )

  # Merge with original object
  object$.sop_order <- seq_len(nrow(object))
  final_result <- merge(
    object,
    summary_stats,
    by = group_cols,
    all.x = TRUE,
    sort = FALSE
  )
  final_result <- final_result[order(final_result$.sop_order), , drop = FALSE]
  final_result$.sop_order <- NULL
  rownames(final_result) <- NULL

  final_result <- restore_sops_attrs(final_result, object)

  # Add bootstrap metadata
  attr(final_result, "n_successful") <- length(boot_results)
  if (method == "fwb") {
    attr(final_result, "fwb_weight_type") <- "exponential"
    attr(final_result, "fwb_weight_scale") <- "cluster_mean_1"
    attr(final_result, "draw_weights_attached") <- FALSE
    attr(final_result, "draw_weight_omission_reason") <- "averaged_sops"
  }

  # Store full bootstrap draws if requested
  if (return_draws) {
    attr(final_result, "draws") <- boot_df
    if (!is.null(baseline_anchor_draws)) {
      attr(final_result, "baseline_anchor_draws") <- baseline_anchor_draws
    }
  }

  final_result
}

inferences_bootstrap_sops_fwb <- function(
  object,
  n_boot,
  workers,
  conf_level,
  conf_type,
  return_draws,
  use_coefstart
) {
  model <- attr(object, "model")
  prediction_data <- attr(object, "newdata_pred")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))
  call_args <- attr(object, "call_args")

  settings <- markov_result_settings(object)
  y_levels <- settings$y_levels
  times <- call_args$times
  id_var <- markov_require_id_var(model, settings$id_var)
  refit_data <- stored_refit_bootstrap_data(object, model)
  validate_refit_bootstrap_data(refit_data, id_var, settings$time_var)

  prediction_data <- ensure_markov_rowid(prediction_data)
  if (!newdata_supplied) {
    # Copies of a profile from `variables` share their patient's weight.
    validate_prediction_weight_ids(
      prediction_weight_profiles(prediction_data),
      id_var
    )
  } else if (isTRUE(return_draws)) {
    warn_fixed_profile_draw_weights("FWB")
  }

  draw_weight_col <- if (newdata_supplied) NULL else "fwb_weight"
  validate_draw_weight_column_available(prediction_data, draw_weight_col)
  draw_weights_attached <- !is.null(draw_weight_col) && isTRUE(return_draws)
  draw_weight_omission_reason <- if (newdata_supplied) {
    "user_supplied_newdata"
  } else {
    NULL
  }

  factor_cols <- bootstrap_state_columns(
    model,
    refit_data,
    settings$p_var,
    settings$p2_var
  )

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    fit_weights <- boot_data$fwb_weight

    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = y_levels,
      absorb = settings$absorb,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )

    m_boot <- boot_result$model
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }
    if (length(missing_states) > 0) {
      warning(
        "Fractional weighted bootstrap unexpectedly dropped outcome-state ",
        "support; skipping this draw.",
        call. = FALSE
      )
      return(NULL)
    }

    sops_array <- tryCatch(
      predict_markov_sops(
        model = m_boot,
        newdata = prediction_data,
        times = times,
        settings = markov_settings_update(
          settings,
          y_levels = state_factor(boot_ylevels),
          absorb = boot_absorb
        )
      ),
      error = function(e) {
        rethrow_crossed_thresholds(e)
        warning("soprob_markov failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(sops_array)) {
      return(NULL)
    }

    prediction_weights <- if (newdata_supplied) {
      NULL
    } else {
      fwb_baseline_weights(
        fwb_weights = fwb_weights,
        baseline_data = prediction_data,
        id_var = id_var
      )
    }

    array_to_df_individual(
      sops_array,
      times,
      state_factor(boot_ylevels),
      prediction_data,
      weights = prediction_weights,
      weight_col = draw_weight_col
    )
  }

  fwb_samples <- generate_fwb_bootstrap_weights(
    data = refit_data,
    id_var = id_var,
    n_boot = n_boot
  )

  boot_results <- apply_to_fwb_bootstrap(
    fwb_samples = fwb_samples,
    analysis_fn = analysis_fn,
    data = refit_data,
    id_var = id_var,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "stats")
  )

  boot_results <- keep_successful_draws(boot_results)

  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  group_cols <- c("rowid", "time", "state")

  summary_stats <- compute_ci_from_draws(
    draws_df = boot_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type,
    point_estimates = as.data.frame(object)
  )

  object$.sop_order <- seq_len(nrow(object))
  final_result <- merge(
    object,
    summary_stats,
    by = group_cols,
    all.x = TRUE,
    sort = FALSE
  )
  final_result <- final_result[order(final_result$.sop_order), , drop = FALSE]
  final_result$.sop_order <- NULL
  rownames(final_result) <- NULL
  final_result <- restore_sops_attrs(final_result, object)

  attr(final_result, "n_successful") <- length(boot_results)
  attr(final_result, "fwb_weight_type") <- "exponential"
  attr(final_result, "fwb_weight_scale") <- "cluster_mean_1"
  attr(final_result, "draw_weights_attached") <- draw_weights_attached
  attr(final_result, "draw_weight_col") <- draw_weight_col
  if (!is.null(draw_weight_omission_reason)) {
    attr(final_result, "draw_weight_omission_reason") <-
      draw_weight_omission_reason
  }

  if (return_draws) {
    attr(final_result, "draws") <- boot_df
  }

  final_result
}


# =============================================================================
# HELPER FUNCTIONS FOR INFERENCE
# =============================================================================

# Refit bootstraps resample the refit data that `mostm()` stores on the
# model: they already exclude rows removed by `subset` and carry the row
# weights. SOP results keep that data in their "refit_data" attribute.
stored_refit_bootstrap_data <- function(object, model) {
  refit_data <- attr(object, "refit_data", exact = TRUE) %||%
    markov_model_refit_data(model)
  if (is.null(refit_data)) {
    mostr_abort(
      "Refit bootstrap inference resamples the data stored on the fitted ",
      "model, but this model has no stored data. Fit it with `mostm()`."
    )
  }
  refit_data
}

# The ordinary bootstrap resamples patients, so a small `by` group can be
# absent from a resample. That resample then has no estimate for the group,
# and the group's interval would rest on fewer resamples without notice. The
# sampled IDs are checked against the stored starting profiles before any
# refit. A missing `by` value among the starting profiles is refused by
# `stop_missing_by_values()` before the point estimate, so the complete-case
# filter below only guards against a profile that belongs to no group.
check_bootstrap_by_groups <- function(boot_ids, baseline_data, id_var, by) {
  if (is.null(by)) {
    return(invisible(NULL))
  }
  complete <- stats::complete.cases(baseline_data[, by, drop = FALSE])
  ids <- as.character(baseline_data[[id_var]])[complete]
  labels <- do.call(
    paste,
    c(
      lapply(by, function(col) {
        paste(
          observed_by_label(col),
          "=",
          as.character(baseline_data[[col]][complete])
        )
      }),
      sep = ", "
    )
  )
  groups <- unique(labels)
  missing_counts <- stats::setNames(integer(length(groups)), groups)
  for (sample in boot_ids) {
    present <- labels[match(as.character(sample$original_id), ids)]
    absent <- setdiff(groups, present)
    missing_counts[absent] <- missing_counts[absent] + 1L
  }
  missing_counts <- missing_counts[missing_counts > 0L]
  if (length(missing_counts) == 0L) {
    return(invisible(NULL))
  }
  shown <- utils::head(missing_counts, 5L)
  mostr_abort(
    "The ordinary bootstrap left out every patient of ",
    paste0(
      "`",
      names(shown),
      "` in ",
      shown,
      " of ",
      length(boot_ids),
      " resamples",
      collapse = "; "
    ),
    if (length(missing_counts) > length(shown)) "; and other groups",
    ", so those groups' intervals would rest on fewer resamples. Use ",
    "`method = \"fwb\"`, which keeps every patient in every refit, or ",
    "combine small groups."
  )
}

# The ordinary bootstrap can also leave out every patient who was ever in a
# rare state. That resample's refit knows no transitions from the state, so it
# can't predict for a prediction profile that starts there; with supplied
# `newdata`, such a profile is in every resample's predictions. The sampled
# IDs are checked against the refit data before any refit, as in
# `check_bootstrap_by_groups()`. A state counts as present when a sampled
# patient holds it in any state column (`state_cols`: the response and the
# previous-state columns), as in `drop_absent_states()`. A state that no
# profile starts from may be missing: it gets probability 0 in that draw.
#
# `starting_states` is a named list from `bootstrap_starting_states()`: for
# each previous-state column, the state labels that prediction starts from.
check_bootstrap_starting_states <- function(
  boot_ids,
  refit_data,
  id_var,
  state_cols,
  starting_states
) {
  starting <- unique(unlist(
    lapply(starting_states, as_state_labels),
    use.names = FALSE
  ))
  starting <- starting[!is.na(starting)]
  if (length(starting) == 0L) {
    return(invisible(NULL))
  }
  ids <- as.character(refit_data[[id_var]])
  state_cols <- intersect(state_cols, names(refit_data))
  missing_counts <- vapply(
    starting,
    function(state) {
      hit <- Reduce(
        `|`,
        lapply(state_cols, function(col) {
          as_state_labels(refit_data[[col]]) %in% state
        }),
        FALSE
      )
      holders <- unique(ids[hit])
      sum(vapply(
        boot_ids,
        function(sample) !any(as.character(sample$original_id) %in% holders),
        logical(1)
      ))
    },
    integer(1)
  )
  names(missing_counts) <- starting
  missing_counts <- missing_counts[missing_counts > 0L]
  if (length(missing_counts) == 0L) {
    return(invisible(NULL))
  }
  states <- names(missing_counts)
  columns <- vapply(
    names(starting_states),
    function(column) {
      hits <- intersect(as_state_labels(starting_states[[column]]), states)
      if (length(hits) == 0L) {
        return(NA_character_)
      }
      paste0("`", column, "` = ", paste(hits, collapse = ", "))
    },
    character(1)
  )
  columns <- columns[!is.na(columns)]
  one <- length(states) == 1L
  mostr_abort(
    "The ordinary bootstrap left out every patient who was ever in ",
    paste0(
      "state ",
      states,
      " in ",
      missing_counts,
      " of the ",
      length(boot_ids),
      " resamples",
      collapse = "; "
    ),
    ", but some prediction profiles start from ",
    if (one) "that state" else "those states",
    " (",
    paste(columns, collapse = "; "),
    "). The refit of a resample without a state has no estimate for leaving ",
    "it, so it can't predict for those profiles. Use `method = \"fwb\"`, ",
    "which keeps every patient in every refit, or remove those profiles ",
    "from `newdata`."
  )
}

# The states that prediction starts from, for each previous-state column in
# `state_vars`: the values of supplied prediction `profiles` (NULL when the
# profiles come from each resample's own patients). `variables` never sets a
# previous-state column here: `inferences()` refuses that for the bootstrap.
bootstrap_starting_states <- function(profiles, state_vars) {
  out <- list()
  for (column in unique(state_vars)) {
    values <- if (is.data.frame(profiles) && column %in% names(profiles)) {
      as_state_labels(profiles[[column]])
    }
    values <- unique(values[!is.na(values)])
    if (length(values) > 0L) {
      out[[column]] <- values
    }
  }
  out
}

validate_refit_bootstrap_data <- function(refit_data, id_var, time_var) {
  validate_markov_id_var(
    id_var,
    refit_data,
    "the refit data stored on the fitted model"
  )

  if (!is.null(time_var) && time_var %in% names(refit_data)) {
    rows_per_patient <- tabulate(match(
      refit_data[[id_var]],
      unique(refit_data[[id_var]])
    ))

    if (all(rows_per_patient == 1)) {
      mostr_abort(
        "Bootstrap inference requires full longitudinal data (all time ",
        "points), but the refit data stored on the fitted model contain one ",
        "row per patient. Fit the model to the full longitudinal data.\n",
        "For simulation: baseline-only prediction profiles are supported."
      )
    }
  }

  invisible(NULL)
}

bootstrap_avg_df_to_matrices <- function(
  boot_avg,
  grid,
  times,
  states,
  variables
) {
  out <- vector("list", nrow(grid))

  for (cf_i in seq_len(nrow(grid))) {
    keep <- rep(TRUE, nrow(boot_avg))
    for (v in names(variables)) {
      keep <- keep & boot_avg[[v]] == grid[[v]][cf_i]
    }

    df <- boot_avg[keep, , drop = FALSE]
    mat <- matrix(
      df$estimate,
      nrow = length(times),
      ncol = length(states),
      byrow = FALSE
    )
    out[[cf_i]] <- mat
  }

  out
}

complete_bootstrap_sop_states <- function(
  boot_avg,
  times,
  y_levels,
  variables,
  by = NULL
) {
  group_cols <- unique(c(names(variables), by))
  group_values <- unique(boot_avg[, group_cols, drop = FALSE])
  base <- expand.grid(
    time = times,
    state = y_levels,
    KEEP.OUT.ATTRS = FALSE
  )

  pieces <- vector("list", nrow(group_values))
  for (i in seq_len(nrow(group_values))) {
    piece <- base
    for (col in group_cols) {
      piece[[col]] <- group_values[[col]][i]
    }
    pieces[[i]] <- piece
  }

  full <- bind_rows_fill(pieces)
  keys <- unique(c("time", "state", group_cols))
  out <- merge(full, boot_avg, by = keys, all.x = TRUE, sort = FALSE)
  out$estimate[is.na(out$estimate)] <- 0
  out
}
