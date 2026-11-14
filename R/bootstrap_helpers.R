#' Fast Group Bootstrap ID Sampler
#'
#' Generates bootstrap samples by resampling group IDs (e.g., patient IDs) with
#' replacement. Returns only the sampled IDs, not the full data, for memory
#' efficiency. This is much faster and more memory-efficient than
#' rsample::group_bootstraps() for datasets with many groups.
#'
#' @param data A data frame containing the patient data
#' @param id_var Character string specifying the name of the ID variable for
#'   group bootstrap (default "id")
#' @param n_boot Number of bootstrap samples to generate
#'
#' @return A list of length \code{n_boot}, where each element is a data frame
#'   (ID lookup table) with columns:
#'   \itemize{
#'     \item original_id: The sampled group ID from the original data
#'     \item new_id: Unique identifier for resampled groups (e.g., "1_1", "1_2")
#'     \item boot_id: Bootstrap iteration number
#'   }
#'
#' @details
#' The ID lookup tables are later joined with the original data on-demand using
#' \code{\link{materialize_bootstrap_sample_indexed}}, so each parallel worker
#' only materializes the data it needs for analysis.
#'
#' The new_id column ensures that if a group is sampled multiple times in
#' the same bootstrap iteration, each instance gets a unique identifier
#' (e.g., "1_1", "1_2"), which is necessary for refitting models that expect
#' unique group identifiers.
#'
#' This approach is dramatically faster than rsample::group_bootstraps() for
#' datasets with many groups or unbalanced group sizes, because rsample uses
#' an inefficient "oversampling then trimming" strategy. See GitHub issue
#' tidymodels/rsample#357 for details.
#'
#' @importFrom stats ave
#' @noRd

fast_group_bootstrap <- function(data, id_var = "id", n_boot) {
  # Input validation
  if (!id_var %in% names(data)) {
    mostr_abort("id_var '", id_var, "' not found in data")
  }

  # Get unique group IDs
  unique_ids <- unique(data[[id_var]])
  n_groups <- length(unique_ids)

  # Generate bootstrap ID samples (IDs only, not full data)
  bootstrap_id_samples <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    # Sample group positions with replacement. Indexing by position keeps the
    # random stream of sample(unique_ids, ...) and also works for a single
    # numeric ID, which sample() would treat as 1:id.
    sampled_pos <- sample.int(n_groups, size = n_groups, replace = TRUE)
    sampled_ids <- unique_ids[sampled_pos]

    # Number repeated draws of the same group 1, 2, ... by counting integer
    # positions. Counting inside the ID vector itself would fail for factor
    # IDs, because ave() returns a factor that cannot hold the counts.
    occurrence <- stats::ave(
      seq_along(sampled_pos),
      sampled_pos,
      FUN = seq_along
    )

    id_table <- data.frame(
      original_id = sampled_ids,
      new_id = paste(as.character(sampled_ids), occurrence, sep = "_"),
      boot_id = i
    )

    bootstrap_id_samples[[i]] <- id_table
  }

  return(bootstrap_id_samples)
}


bootstrap_row_plan <- function(data, id_var) {
  ids <- as.character(data[[id_var]])
  unique_ids <- unique(ids)
  group <- match(ids, unique_ids)
  list(
    ids = unique_ids,
    rows = split(seq_len(nrow(data)), group, drop = TRUE)
  )
}

#' Materialize Bootstrap Sample from ID Lookup
#'
#' Joins bootstrap ID lookup table with the original data to create a full
#' bootstrap dataset. This is used internally by bootstrap functions to
#' materialize data on-demand for memory efficiency.
#'
#' @param boot_ids A data frame with columns: original_id, new_id, boot_id
#'   (output from \code{\link{fast_group_bootstrap}})
#' @param data The original data frame
#' @param id_var Name of the grouping variable in data (e.g., "id")
#' @param row_plan The rows of each group in `data`, from
#'   `bootstrap_row_plan(data, id_var)`. It is computed once and reused for
#'   every bootstrap sample.
#'
#' @return A data frame containing the bootstrap sample with columns from
#'   the original data plus new_id and boot_id
#'
#' @details
#' Each row of `boot_ids` is expanded to all rows of its group in `data`,
#' preserving the sampling order. A group sampled several times appears once
#' per draw, each copy with its own `new_id`.
#'
#' @noRd
materialize_bootstrap_sample_indexed <- function(
  boot_ids,
  data,
  id_var,
  row_plan
) {
  # This function is also serialized to workers independently of the namespace.
  .datatable.aware <- TRUE
  required <- c("original_id", "new_id", "boot_id")
  if (!all(required %in% names(boot_ids))) {
    mostr_abort(
      "`boot_ids` must contain original_id, new_id, and boot_id columns."
    )
  }

  group <- match(as.character(boot_ids$original_id), row_plan$ids)
  row_blocks <- lapply(group, function(index) {
    if (is.na(index)) {
      return(NA_integer_)
    }
    row_plan$rows[[index]]
  })
  counts <- lengths(row_blocks)
  data_rows <- unlist(row_blocks, use.names = FALSE)
  boot_rows <- rep.int(seq_len(nrow(boot_ids)), counts)

  boot_ids_renamed <- boot_ids
  names(boot_ids_renamed)[names(boot_ids_renamed) == "original_id"] <- id_var
  # A matrix column can be a model predictor; as.data.table would split it.
  nested <- function(x) {
    any(vapply(x, function(col) is.list(col) || !is.null(dim(col)), logical(1)))
  }
  if (nested(data) || nested(boot_ids_renamed)) {
    left <- as.data.frame(boot_ids_renamed)[boot_rows, , drop = FALSE]
    right <- as.data.frame(data)[
      data_rows,
      setdiff(names(data), id_var),
      drop = FALSE
    ]
    out <- cbind(left, right)
    rownames(out) <- NULL
    return(out)
  }
  left <- data.table::as.data.table(boot_ids_renamed)[boot_rows]
  right <- data.table::as.data.table(data)[
    data_rows,
    setdiff(names(data), id_var),
    with = FALSE
  ]
  as.data.frame(cbind(left, right))
}

generate_fwb_bootstrap_weights <- function(
  data,
  id_var = "id",
  n_boot,
  weight_dist = "exponential"
) {
  if (!id_var %in% names(data)) {
    mostr_abort("id_var '", id_var, "' not found in data")
  }

  weight_dist <- match.arg(weight_dist, choices = "exponential")
  cluster_ids <- unique(data[[id_var]])
  n_clusters <- length(cluster_ids)
  fwb_samples <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    raw_weights <- stats::rexp(n_clusters, rate = 1)
    weight_mean <- mean(raw_weights)
    if (!is.finite(weight_mean) || weight_mean <= 0) {
      mostr_abort("Invalid fractional bootstrap weights generated.")
    }

    fwb_samples[[i]] <- data.frame(
      original_id = cluster_ids,
      fwb_weight = raw_weights / weight_mean,
      boot_id = i
    )
  }

  fwb_samples
}

materialize_fwb_bootstrap_sample <- function(fwb_weights, data, id_var) {
  if (!all(c("original_id", "fwb_weight", "boot_id") %in% names(fwb_weights))) {
    mostr_abort(
      "`fwb_weights` must contain original_id, fwb_weight, and boot_id columns."
    )
  }

  idx <- match(
    as.character(data[[id_var]]),
    as.character(fwb_weights$original_id)
  )
  if (anyNA(idx)) {
    missing_ids <- unique(as.character(data[[id_var]])[is.na(idx)])
    mostr_abort(
      "Some data IDs were not found in the fractional bootstrap weights: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else "",
      "."
    )
  }

  boot_data <- data
  boot_data$fwb_weight <- fwb_weights$fwb_weight[idx]
  boot_data$boot_id <- fwb_weights$boot_id[idx]
  boot_data
}

fwb_baseline_weights <- function(fwb_weights, baseline_data, id_var) {
  if (!id_var %in% names(baseline_data)) {
    mostr_abort("ID variable '", id_var, "' not found in baseline data.")
  }

  idx <- match(
    as.character(baseline_data[[id_var]]),
    as.character(fwb_weights$original_id)
  )
  if (anyNA(idx)) {
    missing_ids <- unique(as.character(baseline_data[[id_var]])[is.na(idx)])
    mostr_abort(
      "Some baseline IDs were not found in the fractional bootstrap weights: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else "",
      "."
    )
  }

  fwb_weights$fwb_weight[idx]
}


#' Apply Function to Bootstrap Samples with Parallelization
#'
#' Helper function to apply an analysis function to bootstrap samples with
#' optional parallel processing. Handles just-in-time materialization of
#' bootstrap data for memory efficiency.
#'
#' @param boot_samples List of bootstrap ID lookups from
#'   \code{\link{fast_group_bootstrap}}. Each element should be a data frame
#'   with columns: original_id, new_id, boot_id.
#' @param analysis_fn Function to apply to each bootstrap sample. Should accept
#'   a data frame and return a result (can be any type).
#' @param data The original data frame (before bootstrap sampling). This will
#'   be joined with boot_samples on each worker to materialize the data.
#' @param id_var Name of the grouping variable in data (e.g., "id")
#' @param workers Number of workers for parallel processing. If NULL (default)
#'   or 1, sequential processing is used. If > 1, uses parallel processing with
#'   future.callr::callr strategy.
#' @param packages Character vector of packages needed by \code{analysis_fn}.
#'   Default is c("rms", "VGAM", "Hmisc", "stats", "dplyr")
#' @param globals Character vector of global variables/functions needed by
#'   \code{analysis_fn}. These will be exported to each worker.
#'
#' @return A list of length \code{length(boot_samples)} containing the results
#'   from applying \code{analysis_fn} to each bootstrap sample.
#'
#' @details
#' This function implements a memory-efficient just-in-time bootstrap approach:
#' \enumerate{
#'   \item Bootstrap ID lookups are distributed to workers (minimal memory)
#'   \item Each worker receives a copy of the original data
#'   \item Worker materializes bootstrap data by joining IDs with original data
#'   \item Worker runs analysis and returns only the result
#'   \item Full bootstrap datasets are never stored in memory simultaneously
#' }
#'
#' If \code{workers} is NULL or 1, uses sequential processing with \code{lapply()}.
#' If \code{workers > 1}, uses future.callr::callr strategy for parallel
#' processing, which provides isolated R processes for each worker.
#'
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan
#' @importFrom future.callr callr
#' @noRd

apply_to_bootstrap <- function(
  boot_samples,
  analysis_fn,
  data,
  id_var,
  workers = NULL,
  packages = c("rms", "VGAM", "Hmisc", "stats"),
  globals = character(0)
) {

  # Determine whether to parallelize based on workers
  use_parallel <- !is.null(workers) && workers > 1
  row_plan <- bootstrap_row_plan(data, id_var)
  materialize_indexed <- materialize_bootstrap_sample_indexed
  record_crossings <- run_draw_recording_crossings
  record_failure <- run_draw_recording_failure

  # Wrapper function that materializes data then runs analysis
  wrapper_fn <- function(boot_ids) {
    # Materialize bootstrap sample from ID lookup
    boot_data <- materialize_indexed(
      boot_ids,
      data,
      id_var,
      row_plan
    )

    # Run user's analysis function. Samples whose cut-points cross are
    # counted and reported together below, and so are failed samples.
    record_failure(function() record_crossings(analysis_fn, boot_data))
  }

  # Apply function to bootstrap samples
  if (use_parallel) {
    # `plan(new)` returns the user's whole plan; `plan()` alone returns only
    # its outer level, and restoring from that dropped any nested levels.
    old_plan <- plan(callr, workers = workers)
    on.exit(plan(old_plan), add = TRUE)

    results <- furrr::future_map(
      boot_samples,
      wrapper_fn,
      .options = furrr::furrr_options(
        packages = c(packages),
        globals = c(
          globals,
          "data",
          "id_var",
          "analysis_fn",
          "row_plan",
          "materialize_indexed",
          "record_crossings",
          "record_failure"
        )
      )
    )
  } else {
    results <- lapply(boot_samples, wrapper_fn)
  }
  stop_crossed_threshold_draws(results, "bootstrap samples", "sample")

  return(collect_failed_draws(results))
}

apply_to_fwb_bootstrap <- function(
  fwb_samples,
  analysis_fn,
  data,
  id_var,
  workers = NULL,
  packages = c("rms", "VGAM", "Hmisc", "stats"),
  globals = character(0)
) {
  use_parallel <- !is.null(workers) && workers > 1
  materialize_fwb <- materialize_fwb_bootstrap_sample
  record_crossings <- run_draw_recording_crossings
  record_failure <- run_draw_recording_failure

  wrapper_fn <- function(fwb_weights) {
    boot_data <- materialize_fwb(fwb_weights, data, id_var)
    record_failure(function() {
      record_crossings(analysis_fn, boot_data, fwb_weights)
    })
  }

  if (use_parallel) {
    # `plan(new)` returns the user's whole plan; `plan()` alone returns only
    # its outer level, and restoring from that dropped any nested levels.
    old_plan <- plan(callr, workers = workers)
    on.exit(plan(old_plan), add = TRUE)

    results <- furrr::future_map(
      fwb_samples,
      wrapper_fn,
      .options = furrr::furrr_options(
        packages = c(packages),
        globals = c(
          globals,
          "data",
          "id_var",
          "analysis_fn",
          "materialize_fwb",
          "record_crossings",
          "record_failure"
        )
      )
    )
  } else {
    results <- lapply(fwb_samples, wrapper_fn)
  }
  stop_crossed_threshold_draws(results, "bootstrap samples", "sample")

  collect_failed_draws(results)
}

# A bootstrap draw fails when its analysis returns NULL, after a warning that
# says why (for example "Bootstrap model fitting failed: ..."). This runs one
# draw (`run()`, without arguments) and, if it fails, returns a
# `mostr_failed_draw` marker holding the message of the draw's last warning,
# so that `keep_successful_draws()` can report the reason. The warnings
# themselves still reach the user.
run_draw_recording_failure <- function(run) {
  reason <- NULL
  result <- withCallingHandlers(
    run(),
    warning = function(w) {
      reason <<- conditionMessage(w)
    }
  )
  if (is.null(result)) {
    return(structure(list(reason = reason), class = "mostr_failed_draw"))
  }
  result
}

# Turns the failure markers of `run_draw_recording_failure()` back into NULL
# results and keeps their reasons in the "failure_reasons" attribute (NA for
# a successful draw or an unknown reason).
collect_failed_draws <- function(results) {
  failed <- vapply(results, inherits, logical(1), what = "mostr_failed_draw")
  reasons <- rep(NA_character_, length(results))
  reasons[failed] <- vapply(
    results[failed],
    function(marker) marker$reason %||% NA_character_,
    character(1)
  )
  results[failed] <- list(NULL)
  attr(results, "failure_reasons") <- reasons
  results
}

# The successful draws of a bootstrap loop. Each failed draw has already
# warned on its own; this adds one warning that says how many draws failed,
# since the intervals then rest on fewer draws, and the first reason. When
# every draw failed, it stops.
keep_successful_draws <- function(results) {
  n_draws <- length(results)
  failed <- which(vapply(results, is.null, logical(1)))
  reasons <- attr(results, "failure_reasons")
  first_reason <- if (length(failed) > 0L && length(reasons) == n_draws) {
    reasons[failed[1L]]
  } else {
    NA_character_
  }
  first_text <- if (length(failed) > 0L) {
    paste0(
      " The first failure (draw ",
      failed[1L],
      ")",
      if (is.na(first_reason)) "." else paste0(": ", first_reason)
    )
  }
  if (length(failed) == n_draws) {
    mostr_abort("All bootstrap iterations failed.", first_text)
  }
  if (length(failed) > 0L) {
    n_left <- n_draws - length(failed)
    warning(
      length(failed),
      " of the ",
      n_draws,
      " bootstrap draws failed and ",
      if (length(failed) == 1L) "was" else "were",
      " left out, so the confidence intervals rest on the other ",
      n_left,
      if (n_left == 1L) " draw." else " draws.",
      first_text,
      call. = FALSE
    )
  }
  results[setdiff(seq_len(n_draws), failed)]
}


#' Bootstrap Analysis Wrapper
#'
#' Common bootstrap analysis steps used across multiple bootstrap functions.
#' Drops states absent from the sample and refits the model.
#'
#' @param boot_data Bootstrap sample data frame
#' @param model Original fitted model to update with bootstrap data
#' @param factor_cols Character vector of the state columns (e.g.,
#'   c("y", "yprev"))
#' @param original_data Original data frame before bootstrap sampling
#' @param y_levels Original state labels (e.g., 1:6)
#' @param absorb Absorbing state label(s)
#' @param use_coefstart Logical indicating whether to use starting coefficients
#'   from the original model when refitting (only for vglm models). Starting
#'   values (these or a `coefstart` stored with the fit) are used only when
#'   no state is missing from the sample; a failed refit with starting values
#'   is retried once without them.
#' @param fit_weights Optional non-negative row weights for the bootstrap fit.
#'
#' @return A list with components:
#'   \itemize{
#'     \item model: Refitted model on bootstrap data (or NULL if fitting failed)
#'     \item data: Bootstrap data, absent states removed from factor levels
#'     \item ylevels: The original labels of the states present (or NULL if
#'       `y_levels` was not provided)
#'     \item absorb: The labels of the absorbing states present (or NULL)
#'     \item missing_states: Character vector of missing state labels
#'   }
#'
#' @details
#' This function wraps common bootstrap operations:
#' 1. Drop the states absent from the sample from the factor levels of the
#'    state columns, keeping the original labels (`drop_absent_states()`)
#' 2. Refit the model on the bootstrap data, with the rms `datadist` option
#'    unset (see `use_coefstart` for starting values)
#'
#' @seealso `drop_absent_states()`, `fast_group_bootstrap()`,
#'   `apply_to_bootstrap()`
#'
#' @noRd
bootstrap_analysis_wrapper <- function(
  boot_data,
  model,
  factor_cols,
  original_data,
  y_levels = NULL,
  absorb = NULL,
  use_coefstart = FALSE,
  fit_weights = NULL
) {
  fit_model <- bootstrap_refit_model(model)

  prior_weight_column <- markov_model_prior_weight_column(fit_model)
  if (!is.null(prior_weight_column)) {
    # Fail early rather than refitting without the stored row weights.
    bootstrap_prior_weights(boot_data, prior_weight_column)
  }

  if (!is.null(fit_weights)) {
    if (
      !is.numeric(fit_weights) ||
        length(fit_weights) != nrow(boot_data) ||
        any(!is.finite(fit_weights)) ||
        any(fit_weights < 0)
    ) {
      mostr_abort("`fit_weights` must be non-negative finite row weights.")
    }

    original_weights <- bootstrap_original_fit_weights(fit_model, boot_data)
    boot_data$.mostr_fit_weight <- original_weights * fit_weights
  }

  # A sample can lack a state; its refit then knows only the other states,
  # under their original labels.
  present <- drop_absent_states(
    data = boot_data,
    factor_cols = factor_cols,
    original_data = original_data,
    y_levels = y_levels,
    absorb = absorb
  )

  boot_data <- present$data
  boot_ylevels <- present$ylevels
  boot_absorb <- present$absorb
  missing_states <- present$missing_levels

  # rms fitting functions read the user's `datadist` option, and stop when it
  # names an object that does not exist. The coefficients do not depend on it,
  # so refits run with the option unset; it is restored on exit.
  old_options <- options(datadist = NULL)
  on.exit(options(old_options), add = TRUE)

  # Starting coefficients, whether the original fit's (`use_coefstart`) or a
  # `coefstart` stored with the fit, only match a sample with every state:
  # a missing state removes an intercept. Such samples are refitted without
  # any starting values. Other refits with starting values are retried once
  # without them if they fail, e.g. when a covariate level is absent.
  coefstart <- NULL
  if (use_coefstart && inherits(fit_model, "vglm")) {
    coefstart <- stats::coef(fit_model)
  }
  uses_coefstart <- length(missing_states) == 0 &&
    (!is.null(coefstart) ||
      !is.null(bootstrap_model_call(fit_model)$coefstart))
  refit <- function(with_coefstart) {
    update_bootstrap_model(
      fit_model,
      boot_data,
      fit_weights = fit_weights,
      coefstart = if (with_coefstart) coefstart else NULL,
      drop_coefstart = !with_coefstart
    )
  }

  m_boot <- tryCatch(
    if (uses_coefstart) {
      tryCatch(refit(TRUE), error = function(e) refit(FALSE))
    } else {
      refit(FALSE)
    },
    error = function(e) {
      warning("Bootstrap model fitting failed: ", e$message)
      return(NULL)
    }
  )

  return(list(
    model = m_boot,
    data = boot_data,
    ylevels = boot_ylevels,
    absorb = boot_absorb,
    missing_states = missing_states
  ))
}

# `use_coefstart = TRUE` starts each refit from the original coefficients.
# Only `vglm` refits take starting coefficients, and only the refit methods
# refit, so anywhere else the argument would silently do nothing; it is an
# error instead. `method = NULL` skips the method check (for
# `bootstrap_model_coefs()`, which always refits).
check_use_coefstart <- function(use_coefstart, model, method = NULL) {
  if (
    !is.logical(use_coefstart) ||
      length(use_coefstart) != 1L ||
      is.na(use_coefstart)
  ) {
    mostr_abort("`use_coefstart` must be `TRUE` or `FALSE`.")
  }
  if (!use_coefstart) {
    return(invisible(NULL))
  }
  if (!is.null(method) && !method %in% c("bootstrap", "fwb")) {
    mostr_abort(
      "`use_coefstart` applies only to `method = \"bootstrap\"` and ",
      "`method = \"fwb\"`, which refit the model; `method = \"",
      method,
      "\"` does not refit it. Use `use_coefstart = FALSE`, the default."
    )
  }
  if (!is.null(model) && !inherits(model, c("vglm", "robcov_vglm"))) {
    engine <- markov_model_engine(model) %||% class(model)[1L]
    mostr_abort(
      "`use_coefstart` can be `TRUE` only for models fitted with ",
      "`engine = \"vglm\"`, because only vglm refits start from the ",
      "original coefficients; this model was fitted with `engine = \"",
      engine,
      "\"`. Use `use_coefstart = FALSE`, the default."
    )
  }
  invisible(NULL)
}

bootstrap_refit_model <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    if (is.null(model$vglm_fit)) {
      mostr_abort(
        "The supplied robcov_vglm object does not contain the original vglm fit."
      )
    }
    return(model$vglm_fit)
  }

  model
}

# The data column holding the fitted model's response, for example `state` in
# `ordered(state) ~ yprev + time`. A bootstrap refit drops a state that its
# resample lacks from the levels of this column and of the previous-state
# columns.
# Every `mostm()` fit has a readable formula; `"y"` is used only for objects
# without one, such as the stand-in models of the package tests.
bootstrap_response_column <- function(model) {
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  formula <- tryCatch(
    if (isS4(model_chk)) {
      stats::formula(model_chk)
    } else {
      model_chk$sformula %||% stats::formula(model_chk)
    },
    error = function(e) NULL
  )
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    return("y")
  }
  response <- all.vars(formula[[2L]])
  if (length(response) != 1L) {
    return("y")
  }
  response
}

# The state columns of a bootstrap refit: the response and the previous-state
# columns that are present in the refit data.
bootstrap_state_columns <- function(model, refit_data, p_var, p2_var = NULL) {
  intersect(
    c(bootstrap_response_column(model), p_var, p2_var),
    names(refit_data)
  )
}

bootstrap_original_fit_weights <- function(model, data) {
  # `mostm()` fits store their row weights as a column of the refit data, so
  # resampled rows carry their own weights. Unweighted fits have unit weights.
  prior_weight_column <- markov_model_prior_weight_column(model)
  if (is.null(prior_weight_column)) {
    return(rep(1, nrow(data)))
  }
  bootstrap_prior_weights(data, prior_weight_column)
}

bootstrap_prior_weights <- function(data, column) {
  if (!column %in% names(data)) {
    mostr_abort(
      "The bootstrap data do not contain the model's row weights (column `",
      column,
      "`). Resample the refit data stored on the fitted model."
    )
  }
  weights <- data[[column]]
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    mostr_abort(
      "Original model weights must be non-negative finite row weights."
    )
  }
  weights
}

update_bootstrap_model <- function(
  model,
  boot_data,
  fit_weights = NULL,
  coefstart = NULL,
  drop_coefstart = FALSE
) {
  refit_args <- markov_model_refit_args(model)

  # The stored call is the `mostm()` call, so the refit calls `mostm()` with
  # the original engine.
  refit_call <- bootstrap_model_call(model)
  refit_call$data <- quote(boot_data)

  # Stored refit data already exclude the rows removed by `subset`, and
  # starting values are aligned with the original rows, so neither applies to
  # resampled data.
  refit_call$subset <- NULL
  refit_call$etastart <- NULL
  refit_call$mustart <- NULL

  # Row weights come from columns of the resampled data: FWB refits use the
  # original weights times the random weights; ordinary refits use the
  # weights `mostm()` stored with each row.
  prior_weight_column <- markov_model_prior_weight_column(model)
  if (!is.null(fit_weights)) {
    refit_call$weights <- quote(.mostr_fit_weight)
  } else if (!is.null(prior_weight_column)) {
    refit_call$weights <- as.name(prior_weight_column)
  }

  # Starting coefficients supplied for this refit replace a stored
  # `coefstart`; `drop_coefstart` removes the stored one as well.
  if (!is.null(coefstart)) {
    refit_call$coefstart <- quote(coefstart)
  } else if (drop_coefstart) {
    refit_call$coefstart <- NULL
  }

  # `mostm()` stores its argument values from fit time. Inline them so the
  # refit does not re-evaluate expressions that may refer to variables local
  # to the user's fitting frame or changed after fitting.
  inlined <- intersect(names(refit_args), setdiff(names(refit_call), ""))
  inlined <- setdiff(inlined, c("data", "weights"))
  if (!is.null(coefstart)) {
    inlined <- setdiff(inlined, "coefstart")
  }
  for (nm in inlined) {
    refit_call[nm] <- list(refit_args[[nm]])
  }

  # Refits need only the coefficients, so `mostm()` skips the
  # patient-clustered covariance and the starting profiles.
  suppress_expected_refit_warnings(
    model,
    markov_with_coefficients_only_refit(eval(refit_call))
  )
}

bootstrap_model_call <- function(model) {
  if (isS4(model)) {
    return(tryCatch(methods::slot(model, "call"), error = function(e) NULL))
  }
  model$call
}

# State of a bootstrap refit in progress. While `coefficients_only` is TRUE,
# `mostm()` fits only the coefficients: it skips the patient-clustered
# covariance, the starting profiles, and the duplicate-visit warning.
markov_refit_state <- new.env(parent = emptyenv())
markov_refit_state$coefficients_only <- FALSE

markov_refit_coefficients_only <- function() {
  isTRUE(markov_refit_state$coefficients_only)
}

markov_with_coefficients_only_refit <- function(expr) {
  previous <- markov_refit_state$coefficients_only
  markov_refit_state$coefficients_only <- TRUE
  on.exit(markov_refit_state$coefficients_only <- previous, add = TRUE)
  expr
}

# A refit of an orm model runs `rms::orm()` again, so it hides the same
# warning as the original fit in `mostm()`.
suppress_expected_refit_warnings <- function(model, expr) {
  if (!inherits(model, "orm")) {
    return(expr)
  }
  suppress_orm_weights_warning(expr)
}

# Whenever it is given weights, `rms::orm()` warns "currently weights are
# ignored in model validation and bootstrapping orm fits". The warning is
# about rms's own `validate()` and `bootcov()`, which mostr never calls;
# mostr's bootstrap and fractional weighted bootstrap refit with each row's
# weight, so the warning would tell users that their weights are lost when
# they are not. Only this exact message is muffled, so every other warning
# from the fit, such as one about convergence, still reaches the user.
suppress_orm_weights_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (
        grepl(
          "currently weights are ignored in model validation and bootstrapping orm fits",
          conditionMessage(w),
          fixed = TRUE
        )
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
