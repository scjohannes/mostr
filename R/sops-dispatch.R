#' Calculate State Occupancy Probabilities
#'
#' Computes patient-level state occupancy probabilities from a supported
#' Markov transition model. Frequentist first-order models use the optimized
#' native engine; other supported workflows use an engine written in R.
#'
#' @param model A Markov transition model fitted by [mostm()], including the
#'   `robcov_vglm` object that the vglm engine returns. Models
#'   fitted directly with rms, VGAM, or rmsb are not supported.
#' @param newdata A data frame containing one baseline prediction row per
#'   patient or profile.
#' @param times Required consecutive visits, in increasing order. The
#'   recursion takes one transition per value, so leaving out visits does not
#'   skip time; see [sops()]. Unlike [sops()], this function does not warn
#'   when `times` leaves out fitted visits.
#' @param include_re For `blrm` fits, include fitted random-effect draws for
#'   known IDs, which `newdata` gives in the patient ID column (`id_var` of
#'   [mostm()]). Only `blrm` fits have random effects, so `TRUE` gives an
#'   error for other fits.
#' @param n_draws Number of `blrm` posterior draws, or `NULL` for all draws.
#' @param seed Optional seed used when selecting posterior draws.
#' @param ... Reserved for internal posterior draw batching.
#' @return An array of state probabilities. Frequentist fits return
#'   `[patient, time, state]`; Bayesian fits return
#'   `[draw, patient, time, state]`. The state dimension is named by the
#'   fitted states.
#' @inheritSection sops Settings taken from the fit
#' @export
soprob_markov <- function(
  model,
  newdata,
  times = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  ...
) {
  stop_removed_fit_arguments(...)
  stop_include_re_without_blrm(model, include_re)
  settings <- resolve_markov_settings(model)
  if (is.null(settings$y_levels)) {
    stop_no_fitted_states()
  }
  validate_markov_model(model)
  predict_markov_sops(
    model,
    newdata,
    times = times,
    settings = settings,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    ...
  )
}

# The SOP array for `newdata` under resolved `settings` (from
# `resolve_markov_settings()`). Package code calls this rather than
# `soprob_markov()`, because its settings are already resolved and bootstrap
# refits pass only the states present in their resample. The model is checked
# once by the public function that was called (`validate_markov_model()`), and
# the prediction rows are prepared once, here.
predict_markov_sops <- function(
  model,
  newdata,
  times,
  settings,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  ...
) {
  dots <- list(...)
  # The public SOP functions already refused absorbing starts among the user's
  # profiles; this also covers direct calls and `variables` that set the
  # previous state.
  stop_absorbing_start(newdata, settings, source = "rows")
  newdata <- prepare_prediction_rows(model, newdata, settings)
  run_r_engine <- function() {
    soprob_markov_r_engine(
      model = model,
      newdata = newdata,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      ...
    )
  }
  model_fast <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  vglm_fast <- inherits(model_fast, "vglm") &&
    methods::is(model_fast, "vglm") &&
    isTRUE(tryCatch(
      !is.null(methods::slot(model_fast, "constraints")),
      error = function(e) FALSE
    ))
  orm_fast <- inherits(model_fast, "orm") && !inherits(model_fast, "blrm")
  use_fast <- (vglm_fast || orm_fast) && length(dots) == 0L

  if (!use_fast) {
    if (!inherits(model_fast, "blrm")) {
      reason <- if (length(dots) > 0L) {
        "additional backend arguments need the R engine."
      } else {
        "the fitted model or configuration has no compiled execution plan."
      }
      notify_sop_r_engine_fallback(reason)
    }
    return(run_r_engine())
  }

  validate_fast_previous_state_support(model_fast, settings)
  fallback_condition <- NULL
  plan <- tryCatch(
    compile_sop_execution_plan(
      model = model_fast,
      newdata = newdata,
      times = times,
      settings = settings,
      rows_prepared = TRUE
    ),
    mostr_execution_plan_too_large = function(e) {
      fallback_condition <<- e
      NULL
    }
  )
  if (is.null(plan)) {
    notify_sop_r_engine_fallback(conditionMessage(fallback_condition))
    return(run_r_engine())
  }
  resolved_times <- plan$times
  out <- run_sop_execution_plan(plan, get_effective_coefs(model_fast))
  dimnames(out) <- list(
    rownames(newdata),
    resolved_times,
    as_state_labels(settings$y_levels)
  )
  out
}

# A patient in an absorbing state can't leave it, so state occupancy
# probabilities must start from a non-absorbing state. This refuses starting
# rows whose previous state (`p_var`) is in `absorb`, and, for second-order
# models, rows whose state two visits back (`p2_var`) is: that history is
# impossible, since the patient would still be in the absorbing state. States
# are compared as labels, so factor, character and numeric columns all work.
# The columns and absorbing states come from `settings`.
#
# `source` says where the rows come from, which sets the advice in the error:
# - "newdata": the user's `newdata`, one starting profile per row.
# - "stored": the starting profiles that `mostm()` stored, one per patient;
#   `settings$id_var` names their ID column.
# - "rows": prediction rows given to `soprob_markov()`, which hold one copy of
#   each profile for each combination of `variables`.
stop_absorbing_start <- function(
  data,
  settings,
  source = c("newdata", "stored", "rows")
) {
  source <- match.arg(source)
  absorb <- settings$absorb
  p_var <- settings$p_var
  id_var <- if (identical(source, "stored")) settings$id_var
  if (length(absorb) == 0L || !is.data.frame(data) || nrow(data) == 0L) {
    return(invisible(NULL))
  }
  absorb_labels <- as_state_labels(absorb)
  for (column in c(p_var, settings$p2_var)) {
    if (!column %in% names(data)) {
      next
    }
    values <- as_state_labels(data[[column]])
    hit <- !is.na(values) & values %in% absorb_labels
    if (!any(hit)) {
      next
    }
    states <- unique(values[hit])
    state_text <- if (length(states) == 1L) {
      paste0("the absorbing state ", states)
    } else {
      paste0("one of the absorbing states ", paste(states, collapse = ", "))
    }
    rows_text <- switch(
      source,
      newdata = "starting profiles in `newdata`",
      stored = "starting profiles stored by `mostm()`",
      rows = "prediction rows"
    )
    reason <- if (identical(column, p_var)) {
      paste0(
        "A patient in an absorbing state can't leave it, so state occupancy ",
        "probabilities must start from a non-absorbing state."
      )
    } else {
      paste0(
        "`",
        column,
        "` is the state two visits back, and a patient in an absorbing state ",
        "can't leave it, so such a patient would still be in that state at ",
        "the start. State occupancy probabilities must start from a ",
        "non-absorbing history."
      )
    }
    advice <- switch(
      source,
      newdata = "Remove those patients from `newdata`.",
      stored = paste0(
        "Supply `newdata` without those patients",
        if (!is.null(id_var) && id_var %in% names(data)) {
          paste0(
            " (`",
            id_var,
            "` ",
            markov_format_patient_ids(data[[id_var]][hit]),
            ")"
          )
        },
        "."
      ),
      rows = paste0(
        "Remove those rows from `newdata`, and if `variables` sets `",
        column,
        "`, use non-absorbing states only."
      )
    )
    mostr_abort(
      "`",
      column,
      "` is ",
      state_text,
      " for ",
      sum(hit),
      " of the ",
      nrow(data),
      " ",
      rows_text,
      ". ",
      reason,
      " ",
      advice
    )
  }
  invisible(NULL)
}

# Every column the model uses needs a value in every prediction row. Without
# this check a missing value reached the engines: vglm's
# stats::model.matrix() dropped the row, so with two profiles the other
# profile's probabilities filled both rows and were returned without a word,
# and orm and blrm stopped with errors about invalid probabilities that named
# no column. The columns are those of the fitted formula
# (`markov_model_variables()`); a model whose columns can't be read is not
# checked, and columns the model doesn't use may hold missing values. The visit
# time and the gap since the last visit are set at every step of the recursion
# (`assign_sop_visit()`), and the patient ID only looks up blrm random
# effects, so missing values there never reach a prediction. `set_columns`
# names further columns the caller replaces before predicting, and `data_arg`
# the argument the rows came from, for the advice.
stop_missing_predictors <- function(
  model,
  data,
  settings,
  set_columns = character(),
  data_arg = "newdata"
) {
  used <- markov_model_variables(model)
  if (is.null(used) || !is.data.frame(data) || nrow(data) == 0L) {
    return(invisible(NULL))
  }
  columns <- setdiff(
    intersect(used, names(data)),
    c(settings$time_var, settings$gap_var, settings$id_var, set_columns)
  )
  counts <- vapply(
    columns,
    function(column) {
      missing <- is.na(data[[column]])
      if (!is.null(dim(missing))) {
        missing <- rowSums(missing) > 0
      }
      sum(missing)
    },
    numeric(1)
  )
  counts <- counts[counts > 0]
  if (length(counts) == 0L) {
    return(invisible(NULL))
  }
  mostr_abort(
    "The model's predictors have missing values: ",
    paste0("`", names(counts), "` in ", counts, collapse = ", "),
    " of the ",
    nrow(data),
    " prediction rows. Predictions need a value of every column the model ",
    "uses; fill in those values or remove those rows from `",
    data_arg,
    "`."
  )
}

# The rows the engines predict from, checked and with their previous states
# in the form the model reads them. Both are needed wherever prediction rows
# enter an engine: `predict_markov_sops()`, which every SOP route passes
# through, and `compile_sop_execution_plan()` when the inference paths call
# it directly on the rows of an earlier result.
prepare_prediction_rows <- function(model, data, settings) {
  stop_missing_predictors(model, data, settings)
  normalize_prediction_states(model, data, settings)
}

validate_fast_previous_state_support <- function(model, settings) {
  p_var <- settings$p_var
  model_data <- markov_model_data(model)
  if (
    is.null(model_data) ||
      !p_var %in% names(model_data) ||
      !is.factor(model_data[[p_var]])
  ) {
    return(invisible(NULL))
  }
  observed <- markov_fitted_previous_states(model, p_var)
  unsupported <- setdiff(as_state_labels(settings$y_levels), observed)
  unsupported <- setdiff(unsupported, as_state_labels(settings$absorb))
  if (length(unsupported) > 0L) {
    mostr_abort(
      "Model prediction requires transitions from state level(s) not ",
      "represented in the fitted transition data: ",
      paste(unsupported, collapse = ", "),
      ". No fitted row has such a state as its previous state. If it is ",
      "absorbing, such as death, name it in `mostm(absorb = )` and refit."
    )
  }
  invisible(NULL)
}
