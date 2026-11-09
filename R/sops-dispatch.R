#' Calculate State Occupancy Probabilities
#'
#' Computes patient-level state occupancy probabilities from a supported
#' Markov transition model. Frequentist first-order models use the optimized
#' native engine; other supported workflows use the reference engine.
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
#' @param y_levels Ordered outcome-state levels. The number of levels must equal
#'   the fitted threshold count plus one.
#' @param absorb Optional absorbing state levels: states that a patient can't
#'   leave. A row of `newdata` whose `p_var` or `p2_var` value is absorbing
#'   gives an error: such a patient can't move to another state, and the model
#'   has no transitions out of an absorbing state to predict from.
#' @param time_var Name of the model time variable.
#' @param p_var Name of the previous-state variable.
#' @param p2_var Optional second previous-state variable. A non-`NULL` value
#'   selects second-order recursion.
#' @param gap_var Optional elapsed-time-gap variable.
#' @param time_covariates Optional data frame with one row per requested visit
#'   containing precomputed time-varying covariates.
#' @param include_re For `blrm` fits, include fitted random-effect draws for
#'   known IDs. Only `blrm` fits have random effects, so `TRUE` gives an error
#'   for other fits.
#' @param id_var ID column used for `blrm` random effects.
#' @param n_draws Number of `blrm` posterior draws, or `NULL` for all draws.
#' @param seed Optional seed used when selecting posterior draws.
#' @param ... Reserved for internal posterior draw batching.
#' @return An array of state probabilities. Frequentist fits return
#'   `[patient, time, state]`; Bayesian fits return
#'   `[draw, patient, time, state]`.
#' @export
soprob_markov <- function(
  model,
  newdata,
  times = NULL,
  y_levels,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL,
  ...
) {
  stop_include_re_without_blrm(model, include_re)
  dots <- list(...)
  # The public SOP functions already refused absorbing starts among the user's
  # profiles; this also covers direct calls and `variables` that set the
  # previous state.
  stop_absorbing_start(
    newdata,
    absorb = absorb,
    p_var = p_var,
    p2_var = p2_var,
    source = "rows"
  )
  run_reference <- function() {
    soprob_markov_reference(
      model = model,
      newdata = newdata,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      include_re = include_re,
      id_var = id_var,
      n_draws = n_draws,
      seed = seed,
      ...
    )
  }
  model_fast <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  reference_supported <- inherits(
    model,
    c("vglm", "orm", "blrm", "robcov_vglm")
  )
  vglm_fast <- inherits(model_fast, "vglm") &&
    methods::is(model_fast, "vglm") &&
    isTRUE(tryCatch(
      !is.null(methods::slot(model_fast, "constraints")),
      error = function(e) FALSE
    ))
  orm_fast <- inherits(model_fast, "orm") && !inherits(model_fast, "blrm")
  use_fast <- (vglm_fast || orm_fast) &&
    !inherits(model_fast, "blrm") &&
    length(dots) == 0L

  if (!use_fast) {
    if (reference_supported && !inherits(model_fast, "blrm")) {
      reason <- if (length(dots) > 0L) {
        "additional backend arguments require reference evaluation."
      } else {
        "the fitted model or configuration has no compiled execution plan."
      }
      notify_sop_reference_fallback(reason)
    }
    return(run_reference())
  }

  validate_markov_model(model_fast)
  validate_fast_previous_state_support(
    model_fast,
    p_var = p_var,
    y_levels = y_levels,
    absorb = absorb
  )
  fallback_condition <- NULL
  plan <- tryCatch(
    compile_sop_execution_plan(
      model = model_fast,
      newdata = newdata,
      time_covariates = time_covariates,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      builder = "batched",
      output = "array"
    ),
    mostr_execution_plan_too_large = function(e) {
      fallback_condition <<- e
      NULL
    }
  )
  if (is.null(plan)) {
    notify_sop_reference_fallback(conditionMessage(fallback_condition))
    return(run_reference())
  }
  resolved_times <- plan$times
  out <- run_sop_execution_plan(plan, get_effective_coefs(model_fast))
  dimnames(out) <- list(
    rownames(newdata),
    resolved_times,
    as_state_labels(y_levels)
  )
  out
}

# A patient in an absorbing state can't leave it, so state occupancy
# probabilities must start from a non-absorbing state. This refuses starting
# rows whose previous state (`p_var`) is in `absorb`, and, for second-order
# models, rows whose state two visits back (`p2_var`) is: that history is
# impossible, since the patient would still be in the absorbing state. States
# are compared as labels, so factor, character and numeric columns all work.
#
# `source` says where the rows come from, which sets the advice in the error:
# - "newdata": the user's `newdata`, one starting profile per row.
# - "stored": the starting profiles that `mostm()` stored, one per patient;
#   `id_var` names their ID column.
# - "rows": prediction rows given to `soprob_markov()`, which hold one copy of
#   each profile for each combination of `variables`.
stop_absorbing_start <- function(
  data,
  absorb,
  p_var,
  p2_var = NULL,
  source = c("newdata", "stored", "rows"),
  id_var = NULL
) {
  source <- match.arg(source)
  if (length(absorb) == 0L || !is.data.frame(data) || nrow(data) == 0L) {
    return(invisible(NULL))
  }
  absorb_labels <- as_state_labels(absorb)
  for (column in c(p_var, p2_var)) {
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
    stop(
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
      advice,
      call. = FALSE
    )
  }
  invisible(NULL)
}

validate_fast_previous_state_support <- function(
  model,
  p_var,
  y_levels,
  absorb
) {
  model_data <- markov_model_data(model)
  if (
    is.null(model_data) ||
      !p_var %in% names(model_data) ||
      !is.factor(model_data[[p_var]])
  ) {
    return(invisible(NULL))
  }
  observed <- unique(as.character(model_data[[p_var]]))
  unsupported <- setdiff(as_state_labels(y_levels), observed)
  unsupported <- setdiff(unsupported, as_state_labels(absorb))
  if (length(unsupported) > 0L) {
    stop(
      "Model prediction requires transitions from state level(s) not ",
      "represented in the fitted transition data: ",
      paste(unsupported, collapse = ", "),
      ". If a level is absorbing, pass it via `absorb`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
