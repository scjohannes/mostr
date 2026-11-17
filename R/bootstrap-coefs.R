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
#' @param n_draws Number of bootstrap refits, as `n_draws` in [inferences()].
#' @param workers Number of workers for parallel processing. If NULL (default)
#'   or 1, sequential processing is used. If > 1, uses parallel processing.
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
#' The coefficients of a factor previous state compare each state with the
#' first one, the reference level: \code{yprev=3} is the effect of having been
#' in state 3 rather than state 1 at the previous visit, and the intercepts
#' belong to a patient whose previous state is state 1. A resample can hold no
#' row whose previous state is the reference level. Its refit would then
#' compare with the next level instead, so \code{yprev=3} would mean state 3
#' rather than state 2, under the same name. In such a row the intercepts and
#' the previous-state coefficients are therefore \code{NA}, the other
#' coefficients (for example those of time and treatment) are kept, and a
#' warning gives the number of rows affected.
#'
#' A model fitted with explicit \code{constraints} (\code{engine = "vglm"})
#' is refused before any refit when some resample lacks an outcome state: its
#' constraint matrices have one row per cut-point of the full model and do
#' not fit a refit with fewer cut-points. Merge rarely observed states, or
#' give partial proportional odds through \code{family}, for example
#' \code{VGAM::cumulative(parallel = FALSE ~ tx)}, which builds the matrices
#' for each refit.
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
#' # 200 simulated patients followed for 14 days; state 8, death, is
#' # absorbing
#' trajectories <- sim_markov("actt2", n_patients = 200, times = 1:14, seed = 1)
#' model_data <- prepare_markov_data(trajectories, absorb = 8)
#' fit <- mostm(
#'   ordered(y) ~ yprev + rms::rcs(time, 3) + tx,
#'   data = model_data,
#'   id_var = "id",
#'   absorb = "8"
#' )
#'
#' # Resample the patients stored on the fit and refit the model 10 times;
#' # use 1000 or more refits in an analysis
#' set.seed(1)
#' boot_coefs <- bootstrap_model_coefs(fit, n_draws = 10)
#' dim(boot_coefs)
#'
#' # Percentile intervals for the treatment effect (log odds ratio)
#' stats::quantile(boot_coefs$tx, c(0.025, 0.975), na.rm = TRUE)
#'
#' @export

bootstrap_model_coefs <- function(
  model,
  n_draws,
  workers = NULL,
  use_coefstart = FALSE
) {

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

  response_col <- bootstrap_response_column(model)
  # The factor previous-state columns and their reference levels: the first
  # level that the stored data hold, which the full fit compares with.
  previous_cols <- intersect(c(settings$p_var, settings$p2_var), names(data))
  previous_cols <- previous_cols[vapply(
    previous_cols,
    function(col) is.factor(data[[col]]),
    logical(1)
  )]
  reference_levels <- vapply(
    previous_cols,
    function(col) observed_state_labels(data, col)[1L],
    character(1)
  )

  # Generate bootstrap ID samples using fast helper (memory-efficient JIT approach)
  boot_ids <- fast_group_bootstrap(
    data = data,
    id_var = id_var,
    n_boot = n_draws
  )
  check_bootstrap_constraint_outcomes(
    boot_ids,
    refit_data = data,
    id_var = id_var,
    response_col = response_col,
    model = model,
    fwb = FALSE
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
      use_coefstart = use_coefstart,
      response_col = response_col
    )

    m_boot <- boot_result$model

    if (is.null(m_boot)) {
      return(NULL)
    }
    coefs <- coef(m_boot)

    # A resample without the reference level of a factor previous-state
    # column refits that column against another level. Its coefficients,
    # and the intercepts (which belong to a patient in the reference level),
    # then estimate other quantities than the same-named coefficients of the
    # full fit, so they are set to NA. Coefficients of other terms keep their
    # meaning and are kept. Only a column present in the refit data can lack
    # a level; a refit result without data has nothing to check.
    lacks_reference <- previous_cols[vapply(
      previous_cols,
      function(col) {
        col %in% names(boot_result$data) &&
          !reference_levels[[col]] %in% as_state_labels(boot_result$data[[col]])
      },
      logical(1)
    )]
    if (length(lacks_reference) > 0L) {
      na_positions <- unique(unlist(lapply(
        lacks_reference,
        function(col) bootstrap_coef_positions(m_boot, col)
      )))
      coefs[na_positions] <- NA_real_
      return(structure(as.list(coefs), lacks_reference = TRUE))
    }
    as.list(coefs)
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

  n_lacking <- sum(vapply(
    bs_coefs_list,
    function(x) isTRUE(attr(x, "lacks_reference")),
    logical(1)
  ))
  if (n_lacking > 0L) {
    warning(
      "In ",
      n_lacking,
      " of the ",
      n_draws,
      " resamples, no patient was observed after being in the reference ",
      "(first) level of the previous state (",
      paste0(
        "`",
        previous_cols,
        "` = ",
        reference_levels,
        collapse = " or "
      ),
      "). The refits of those resamples compare the previous states with ",
      "another level, so their intercepts and previous-state coefficients ",
      "would mean something different from those in the other rows; they ",
      "are NA. The ",
      "other coefficients of those rows are kept.",
      call. = FALSE
    )
  }

  # Convert to a data frame with boot_id and one column per coefficient.
  result <- named_list_to_wide(bs_coefs_list, id = seq_len(n_draws))

  return(result)
}

# The positions in `coef(fit)` of the intercepts and of the coefficients of
# every term that involves the column `column` (its main effect and its
# interactions), read from the fit's own term information rather than from
# coefficient names. `rms` fits list each term's coefficient positions in
# `$assign`, after the `non.slopes` intercepts. `vglm` fits map each term to
# model matrix columns in `@assign`, and each column to as many coefficients
# as its constraint matrix has columns.
bootstrap_coef_positions <- function(fit, column) {
  fit <- bootstrap_refit_model(fit)
  term_factors <- attr(stats::terms(fit), "factors")
  involves <- vapply(
    rownames(term_factors),
    function(variable) column %in% all.vars(str2lang(variable)),
    logical(1)
  )
  terms_with_column <- colnames(term_factors)[
    colSums(term_factors[involves, , drop = FALSE]) > 0
  ]

  if (inherits(fit, "vglm")) {
    constraint_list <- VGAM::constraints(fit)
    widths <- vapply(constraint_list, ncol, integer(1))
    ends <- cumsum(widths)
    positions <- lapply(seq_along(widths), function(j) {
      seq.int(ends[j] - widths[j] + 1L, length.out = widths[j])
    })
    names(positions) <- names(constraint_list)
    model_columns <- unlist(fit@assign[terms_with_column], use.names = FALSE)
    return(c(
      positions[["(Intercept)"]],
      unlist(positions[model_columns], use.names = FALSE)
    ))
  }

  c(
    seq_len(fit$non.slopes),
    unlist(fit$assign[terms_with_column], use.names = FALSE)
  )
}
