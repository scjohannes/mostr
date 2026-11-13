# Coefficient bootstrap summaries.

#' Bootstrap confidence intervals for model coefficients
#'
#' Resamples patients with replacement, refits the model to each resample, and
#' returns the refitted coefficients (intercepts and slopes). Percentiles of
#' each column give bootstrap confidence intervals for that coefficient.
#'
#' @param model A model fitted with [mostm()] and the `"orm"` or `"vglm"`
#'   engine. [mostm()] stores the data used for fitting, and the bootstrap
#'   resamples patients from those stored data.
#' @param n_boot Number of bootstrap samples.
#' @param workers Number of workers for parallel processing. If NULL (default)
#'   or 1, sequential processing is used. If > 1, uses parallel processing.
#' @param parallel (Deprecated) Logical indicating whether to use parallel
#'   processing. Use the \code{workers} parameter instead.
#' @param use_coefstart Logical. If `TRUE`, each refit starts from the
#'   original model coefficients, which can speed up convergence. Only models
#'   fitted with `engine = "vglm"` accept starting coefficients, so `TRUE` with
#'   an `orm` model is an error. Starting values, including a `coefstart`
#'   passed to [mostm()], are used only when a resample contains every outcome
#'   state; a refit with starting values that fails is retried once without
#'   them. Default is `FALSE`.
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
#'   \item When a state is absent from the resample, removes it from the
#'     factor levels of the response column (the variable on the left of the
#'     model formula) and of \code{yprev}. The other states keep their labels;
#'     they are not renumbered. States are matched by their labels, so numbers
#'     and text labels such as \code{"home"} work alike. Other factor columns
#'     are covariates and are left as they are.
#'   \item Refits the model by calling [mostm()] again with the original
#'     engine and arguments, and extracts its coefficients
#' }
#'
#' When a state is missing from a resample, the refit has fewer intercepts than
#' the original model. Coefficients named after a state, such as the
#' intercept \code{y>=2} of an \code{"orm"} fit or the previous-state
#' coefficient \code{yprev=2}, are then \code{NA} in that row. The intercepts of
#' a \code{"vglm"} fit are numbered (\code{(Intercept):1}, \code{(Intercept):2},
#' ...) in order of the states present, so in that row they belong to other
#' cut-points than in the rows with every state.
#'
#' The refits do not use an rms `datadist`: a `datadist` you have set with
#' `options(datadist = ...)` is ignored during the refits and left as it is.
#'
#' With `workers > 1`, the refits run in separate R processes started by
#' the `future.callr` package.
#'
#' \strong{Previous-state type:} \code{yprev} may be a factor for categorical
#' previous-state effects or numeric for linear/spline effects. Refits keep
#' the column type and the state numbers used to fit the original model.
#'
#' @keywords bootstrap coefficients confidence intervals
#'
#' @importFrom stats update coef
#'
#' @examples
#' \dontrun{
#' fit <- mostm(
#'   ordered(y) ~ tx + yprev + time,
#'   data = my_data,
#'   id_var = "id",
#'   engine = "vglm"
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
    mostr_abort(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }
  check_use_coefstart(use_coefstart, model)

  data <- markov_model_refit_data(model)
  if (is.null(data)) {
    mostr_abort(
      "The model has no stored data to resample. Fit it with `mostm()`."
    )
  }

  id_var <- markov_require_id_var(model)
  if (!id_var %in% names(data)) {
    mostr_abort("id_var '", id_var, "' not found in the stored model data")
  }

  # The state columns whose factor levels lose a state that a resample lacks:
  # the response and the previous-state columns the fit records (`p_var`, and
  # `p2_var` for a second-order fit). Other factor columns are covariates,
  # whose levels are not states. A numeric previous state keeps its values for
  # linear/spline previous-state models.
  settings <- markov_model_settings(model)
  factor_cols <- bootstrap_state_columns(
    model,
    data,
    settings$p_var,
    settings$p2_var
  )

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
