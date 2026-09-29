# Coefficient bootstrap summaries.

#' Bootstrap confidence intervals for model coefficients
#'
#' Resamples patients with replacement, refits the model to each resample, and
#' returns the refitted coefficients (intercepts and slopes). Percentiles of
#' each column give bootstrap confidence intervals for that coefficient.
#'
#' @param model A model fitted with [orm_markov()] or [vglm_markov()]. These
#'   wrappers store the data used for fitting, and the bootstrap resamples
#'   patients from those stored data.
#' @param n_boot Number of bootstrap samples.
#' @param workers Number of workers for parallel processing. If NULL (default)
#'   or 1, sequential processing is used. If > 1, uses parallel processing.
#' @param parallel (Deprecated) Logical indicating whether to use parallel
#'   processing. Use the \code{workers} parameter instead.
#' @param id_var Name of the patient ID column in the stored data. The default,
#'   `NULL`, uses the `id_var` given when the model was fitted. Supply it only
#'   when the model was fitted without `id_var`.
#' @param use_coefstart Logical. If `TRUE`, `vglm` refits start from the
#'   original model coefficients, which can speed up convergence. Starting
#'   values, including a `coefstart` passed to [vglm_markov()], are used only
#'   when a resample contains every outcome state; a refit with starting values
#'   that fails is retried once without them. Default is `FALSE`.
#'
#' @return A data frame with one row per bootstrap sample:
#'   \itemize{
#'     \item boot_id: Bootstrap iteration number
#'     \item One column per model coefficient (intercepts and slopes). Failed
#'       refits give missing values.
#'   }
#'
#' @details
#' Patients, not rows, are resampled, so every resampled patient keeps all of
#' their visits and the within-patient correlation is preserved. A patient
#' drawn more than once enters the refit as separate patients. The refit uses
#' the same fitting rows as the original model: rows excluded by `subset` stay
#' excluded, and fitting weights stay attached to their rows.
#'
#' For each bootstrap iteration:
#' \enumerate{
#'   \item Draws patients with replacement from the stored data
#'   \item Renumbers outcome and previous-state levels to consecutive integers
#'     when a state is absent from the resample
#'   \item Updates the `rms` datadist for `orm` models (safe with
#'     future.callr workers)
#'   \item Refits the model and extracts its coefficients
#' }
#'
#' When a state is missing from a resample, the refit has fewer intercepts than
#' the original model, so the intercept columns of that row refer to renumbered
#' states.
#'
#' Parallelization is handled via future.callr::callr strategy, which provides
#' isolated R processes for each worker, making it safe to modify the global
#' environment (for datadist) without conflicts.
#'
#' \strong{Previous-state type:} \code{yprev} may be a factor for categorical
#' previous-state effects or numeric for linear/spline effects. Bootstrap
#' releveling preserves the column type used to fit the original model.
#'
#' @keywords bootstrap coefficients confidence intervals
#'
#' @importFrom stats update coef
#'
#' @examples
#' \dontrun{
#' fit <- vglm_markov(
#'   ordered(y) ~ tx + yprev + time,
#'   family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
#'   data = my_data,
#'   id_var = "id"
#' )
#'
#' # Resample the patients stored on the fit.
#' bs_coefs <- bootstrap_model_coefs(fit, n_boot = 1000)
#'
#' # 95% percentile confidence intervals
#' apply(
#'   bs_coefs[, -1],
#'   2,
#'   stats::quantile,
#'   probs = c(0.025, 0.975),
#'   na.rm = TRUE
#' )
#' }
#' @export

bootstrap_model_coefs <- function(
  model,
  n_boot,
  workers = NULL,
  parallel = NULL,
  id_var = NULL,
  use_coefstart = FALSE
) {
  # Handle deprecated parallel parameter
  if (!is.null(parallel)) {
    warning(
      "The 'parallel' argument is deprecated. ",
      "Please use 'workers' instead.\n",
      "  - For sequential processing: workers = NULL or workers = 1\n",
      "  - For parallel processing: workers = N (e.g., workers = 8)"
    )
    if (parallel && is.null(workers)) {
      workers <- parallel::detectCores() - 1
    }
  }

  # Check model class
  if (!inherits(model, c("orm", "vglm", "robcov_vglm"))) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  data <- markov_model_refit_data(model)
  if (is.null(data)) {
    stop(
      "The model has no stored data to resample. Fit it with `orm_markov()` ",
      "or `vglm_markov()`.",
      call. = FALSE
    )
  }

  id_var <- markov_model_id_var(model, id_var)
  if (is.null(id_var)) {
    stop(
      "No patient ID column is known for this model. Fit it with `id_var`, ",
      "or supply `id_var`.",
      call. = FALSE
    )
  }
  if (!id_var %in% names(data)) {
    stop("id_var '", id_var, "' not found in the stored model data")
  }

  # Identify state columns that may need releveling in bootstrap samples.
  # Numeric yprev is preserved for linear/spline previous-state models.
  factor_cols <- unique(c(names(data)[sapply(data, is.factor)], "yprev"))
  factor_cols <- intersect(factor_cols, names(data))

  # Generate bootstrap ID samples using fast helper (memory-efficient JIT approach)
  boot_ids <- fast_group_bootstrap(
    data = data,
    id_var = id_var,
    n_boot = n_boot
  )

  # Define analysis function
  analysis_fn <- function(boot_data) {
    # A patient drawn more than once enters the refit as separate patients.
    if (!is.null(boot_data$new_id)) {
      boot_data[[id_var]] <- boot_data$new_id
    }

    # Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = data,
      y_levels = NULL,
      absorb = NULL,
      update_datadist = inherits(model, "orm"),
      use_coefstart = use_coefstart
    )

    m_boot <- boot_result$model

    # Extract coefficients
    if (!is.null(m_boot)) {
      coefs <- coef(m_boot)
      return(as.list(coefs))
    } else {
      return(NULL)
    }
  }

  # Apply analysis function to bootstrap samples with JIT materialization
  bs_coefs_list <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = data,
    id_var = id_var,
    workers = workers,
    packages = c("rms", "VGAM", "stats"),
    globals = c(
      "model",
      "factor_cols",
      "id_var"
    )
  )

  # Convert to a data frame with boot_id and one column per coefficient.
  result <- named_list_to_wide(bs_coefs_list, id = seq_len(n_boot))

  return(result)
}
