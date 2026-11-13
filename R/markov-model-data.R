# Markov model data metadata helpers.

# Rebuild the rows used by an rms fit. `subset` is the already evaluated row
# selection and `weight_column` names the stored row-weight column of `data`.
markov_rms_model_frame <- function(
  formula,
  data,
  dots,
  subset = NULL,
  weight_column = NULL,
  eval_env = parent.frame()
) {
  na_action <- if (!is.null(dots$na.action)) {
    eval(dots$na.action, eval_env)
  } else {
    utils::getFromNamespace("na.delete", "Hmisc")
  }
  args <- list(formula = formula, data = data)
  if (!is.null(subset)) {
    args$subset <- subset
  }
  if (!is.null(weight_column)) {
    args$weights <- as.name(weight_column)
  }
  args$na.action <- na_action
  args$drop.unused.levels <- TRUE
  eval(as.call(c(list(quote(stats::model.frame)), args)))
}

validate_markov_id_var <- function(id_var, data, data_arg = "data") {
  if (!is.character(id_var) || length(id_var) != 1 || is.na(id_var)) {
    mostr_abort("`id_var` must be a single column name.")
  }
  if (!id_var %in% names(data)) {
    mostr_abort("ID variable '", id_var, "' not found in ", data_arg, ".")
  }
  id_var
}

stop_missing_mostm_id_var <- function() {
  mostr_abort(
    "`mostm()` needs `id_var`, the name of the column in `data` that ",
    "identifies patients, for example `id_var = \"id\"`. Standard errors that ",
    "allow for correlation between the transitions of one patient, the ",
    "stored starting row of each patient, and the bootstrap all need to know ",
    "which rows belong to the same patient."
  )
}

warn_duplicate_markov_id_time <- function(
  data,
  id_var,
  time_var = "time",
  wrapper = "mostm()"
) {
  if (!time_var %in% names(data)) {
    return(invisible(NULL))
  }

  key <- data[c(id_var, time_var)]
  if (anyDuplicated(key) == 0L) {
    return(invisible(NULL))
  }

  warning(
    wrapper,
    " received data with duplicate `",
    id_var,
    "` x `",
    time_var,
    "` combinations. This is usually accidental in Markov transition data; ",
    "automatic SOP prediction from stored data will use one prediction row ",
    "per `",
    id_var,
    "`.",
    call. = FALSE
  )

  invisible(NULL)
}

# Capture the values of the arguments a user passed to `mostm()`.
# Bootstrap and FWB refits reuse these values instead of re-evaluating the
# stored call's argument expressions, which may refer to variables that exist
# only in the user's fitting frame or that changed after fitting. Arguments
# evaluated inside `data` (or replaced at refit time) are not captured.
markov_capture_refit_args <- function(call, env, formal_names, ...) {
  skip <- c("data", "subset", "weights", "etastart", "mustart", "offset")
  call_names <- setdiff(names(call)[-1L], "")
  out <- list()

  for (nm in intersect(call_names, setdiff(formal_names, c("...", skip)))) {
    out[nm] <- list(get(nm, envir = env, inherits = FALSE))
  }

  dot_names <- ...names()
  for (i in seq_len(...length())) {
    nm <- dot_names[[i]]
    if (is.null(nm) || is.na(nm) || !nzchar(nm) || nm %in% skip) {
      next
    }
    out[nm] <- list(...elt(i))
  }

  out
}

# Name of the column in which `mostm()` stores row weights. Stored refit
# data carry this column, so bootstrap resampling keeps each row's weight.
markov_prior_weight_column <- function() {
  ".mostr_prior_weight"
}

# Evaluate a row-weight expression once, with the data and enclosing
# environment the engine would use, and validate the result.
markov_eval_prior_weights <- function(expr, data, eval_env) {
  weights <- eval(expr, data, eval_env)
  if (is.null(weights)) {
    return(NULL)
  }
  if (is.data.frame(weights)) {
    weights <- as.matrix(weights)
  }
  if (is.matrix(weights) && ncol(weights) == 1L) {
    weights <- weights[, 1L]
  }
  n_weights <- if (is.matrix(weights)) nrow(weights) else length(weights)
  if (!is.numeric(weights) || n_weights != nrow(data)) {
    mostr_abort(
      "`weights` must be numeric with one value per row of `data`."
    )
  }
  if (any(!is.finite(weights))) {
    mostr_abort("`weights` must be finite.")
  }
  if (any(weights < 0)) {
    mostr_abort("negative weights not allowed")
  }
  weights
}

markov_attach_model_data <- function(
  model,
  data,
  id_var = NULL,
  refit_data = NULL,
  starting_profile_data = NULL,
  starting_profile_metadata = NULL,
  refit_args = NULL,
  prior_weight_column = NULL
) {
  if (!is.null(refit_args)) {
    attr(model, "markov_refit_args") <- refit_args
  }
  if (!is.null(prior_weight_column)) {
    attr(model, "markov_prior_weight_column") <- prior_weight_column
  }
  if (!is.null(data)) {
    attr(model, "markov_data") <- data
  }
  if (!is.null(id_var)) {
    attr(model, "markov_id_var") <- id_var
  }
  if (!is.null(refit_data)) {
    attr(model, "markov_refit_data") <- refit_data
  }
  if (!is.null(starting_profile_data)) {
    attr(model, "markov_starting_profile_data") <- starting_profile_data
  }
  if (!is.null(starting_profile_metadata)) {
    attr(model, "markov_starting_profile_metadata") <-
      starting_profile_metadata
  }
  model
}

# Provenance: `mostm()` records the engine that fitted a model. The SOP and
# inference workflows accept only models with this record.
markov_engines <- function() {
  c("orm", "vglm", "blrm")
}

markov_set_engine <- function(model, engine) {
  if (
    !is.character(engine) ||
      length(engine) != 1L ||
      is.na(engine) ||
      !engine %in% markov_engines()
  ) {
    mostr_abort("Internal Markov model engine provenance is invalid.")
  }
  attr(model, "markov_engine") <- engine
  model
}

markov_model_engine <- function(model) {
  engine <- attr(model, "markov_engine", exact = TRUE)
  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    inner <- attr(model$vglm_fit, "markov_engine", exact = TRUE)
    if (!is.null(engine) && !is.null(inner) && !identical(engine, inner)) {
      return(NULL)
    }
    engine <- engine %||% inner
  }
  if (
    !is.character(engine) ||
      length(engine) != 1L ||
      is.na(engine) ||
      !engine %in% markov_engines()
  ) {
    return(NULL)
  }
  engine
}

# Names of the data columns on the right-hand side of the fitted model's
# formula, including those inside terms such as `rms::rcs(age, 3)`, `poly()`,
# interactions and a blrm `cluster(id)`, and those of a blrm partial proportional
# odds formula (`ppo`). They come from the stored formula, never from the
# coefficients. `NULL` when the variables cannot be read from the object, for
# example a formula with `.`, so that callers can leave such models alone.
markov_model_variables <- function(model) {
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  formulas <- tryCatch(
    if (isS4(model_chk)) {
      list(stats::formula(model_chk))
    } else {
      list(
        model_chk$sformula %||% stats::formula(model_chk),
        model_chk$zsformula
      )
    },
    error = function(e) NULL
  )
  formulas <- Filter(Negate(is.null), formulas)
  if (
    length(formulas) == 0L ||
      !inherits(formulas[[1L]], "formula") ||
      !all(vapply(formulas, inherits, logical(1), "formula"))
  ) {
    return(NULL)
  }
  # An object that only looks like a fit, such as a robcov_vglm without its
  # vglm fit, can give an empty formula, whose right-hand side can't be read.
  variables <- tryCatch(
    unique(unlist(lapply(formulas, function(formula) {
      all.vars(formula[[length(formula)]])
    }))),
    error = function(e) NULL
  )
  if (is.null(variables) || "." %in% variables) {
    return(NULL)
  }
  # rmsb::blrm() keeps `cluster(id)` out of its stored formula. The random
  # effect of each patient is looked up by that column.
  cluster <- if (inherits(model_chk, "blrm")) model_chk$clusterInfo$name
  unique(c(variables, if (is.character(cluster)) cluster))
}

# A model derived from `source`, such as a refit or the inner vglm fit, keeps
# its engine record and, when the derived model has none of its own, its
# settings.
markov_inherit_engine <- function(model, source) {
  engine <- markov_model_engine(source)
  if (is.null(engine)) {
    return(model)
  }
  settings <- markov_model_metadata_attr(source, "markov_settings")
  if (
    inherits(settings, "markov_settings") &&
      is.null(attr(model, "markov_settings", exact = TRUE))
  ) {
    model <- markov_set_settings(model, settings)
  }
  markov_set_engine(model, engine)
}

markov_model_data <- function(model) {
  data <- attr(model, "markov_data", exact = TRUE)
  if (!is.null(data)) {
    return(data)
  }

  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    data <- attr(model$vglm_fit, "markov_data", exact = TRUE)
    if (!is.null(data)) {
      return(data)
    }
  }

  NULL
}

# The states that occur as the previous state in the data the model was fitted
# to, as labels, or NULL when the fit does not carry those data. An absorbing
# state never occurs there: a patient in it has no further transition, so the
# model has no linear predictor for it as a previous state.
markov_fitted_previous_states <- function(model, p_var) {
  model_data <- markov_model_data(model)
  if (is.null(model_data) || !p_var %in% names(model_data)) {
    return(NULL)
  }
  values <- model_data[[p_var]]
  unique(as_state_labels(values[!is.na(values)]))
}

markov_model_id_var <- function(model, id_var = NULL) {
  if (!is.null(id_var)) {
    return(id_var)
  }

  out <- attr(model, "markov_id_var", exact = TRUE)
  if (!is.null(out)) {
    return(out)
  }

  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    out <- attr(model$vglm_fit, "markov_id_var", exact = TRUE)
    if (!is.null(out)) {
      return(out)
    }
  }

  NULL
}

# The patient-ID column: `id_var` when given, otherwise the column named in the
# `mostm()` call. Every `mostm()` fit records one, so a model without it was
# not made by `mostm()` or lost its metadata; no column is guessed.
markov_require_id_var <- function(model, id_var = NULL) {
  id_var <- markov_model_id_var(model, id_var)
  if (is.null(id_var)) {
    mostr_abort(
      "The model does not record which column of its data identifies ",
      "patients. Refit it with `mostm(..., id_var = ...)`."
    )
  }
  id_var
}

markov_model_metadata_attr <- function(model, name) {
  value <- attr(model, name, exact = TRUE)
  if (!is.null(value)) {
    return(value)
  }

  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    return(attr(model$vglm_fit, name, exact = TRUE))
  }

  NULL
}

markov_model_refit_args <- function(model) {
  markov_model_metadata_attr(model, "markov_refit_args")
}

markov_model_prior_weight_column <- function(model) {
  markov_model_metadata_attr(model, "markov_prior_weight_column")
}

markov_model_refit_data <- function(model) {
  markov_model_metadata_attr(model, "markov_refit_data") %||%
    markov_model_data(model)
}

markov_model_starting_profile_metadata <- function(model) {
  markov_model_metadata_attr(model, "markov_starting_profile_metadata")
}

markov_format_patient_ids <- function(ids, limit = 5L) {
  ids <- unique(as.character(ids))
  paste0(
    paste(utils::head(ids, limit), collapse = ", "),
    if (length(ids) > limit) " ..." else ""
  )
}

markov_validate_starting_profiles <- function(
  model,
  required_variables = NULL
) {
  profiles <- markov_model_metadata_attr(model, "markov_starting_profile_data")
  metadata <- markov_model_starting_profile_metadata(model)
  if (is.null(profiles) || is.null(metadata)) {
    mostr_abort(
      "Provide newdata: the model has no stored starting row for each ",
      "patient. Fits made by `mostm()` store these rows."
    )
  }

  id_var <- metadata$id_var
  fitted_ids <- as.character(metadata$fitted_ids)
  if (
    is.null(id_var) ||
      !id_var %in% names(profiles) ||
      anyNA(fitted_ids) ||
      any(fitted_ids == "") ||
      anyDuplicated(fitted_ids)
  ) {
    mostr_abort("Stored starting-profile patient metadata is invalid.")
  }

  profile_ids <- as.character(profiles[[id_var]])
  missing_ids <- setdiff(fitted_ids, profile_ids)
  if (length(missing_ids)) {
    mostr_abort(
      "Fitted patient(s) lack a row at the designated starting visit `",
      as.character(metadata$first_followup_time),
      "`: ",
      markov_format_patient_ids(missing_ids),
      ". A later likelihood row is not substituted for the starting profile."
    )
  }

  duplicated_ids <- unique(profile_ids[
    duplicated(profile_ids) |
      duplicated(profile_ids, fromLast = TRUE)
  ])
  if (length(duplicated_ids)) {
    mostr_abort(
      "Fitted patient(s) have multiple rows at the designated starting visit ",
      "`",
      as.character(metadata$first_followup_time),
      "`: ",
      markov_format_patient_ids(duplicated_ids),
      ". Exactly one starting profile is required per fitted patient."
    )
  }

  required <- unique(c(metadata$required_variables, required_variables))
  missing_variables <- setdiff(required, names(profiles))
  if (length(missing_variables)) {
    mostr_abort(
      "Stored starting profiles are missing required prediction variables: ",
      paste(missing_variables, collapse = ", "),
      "."
    )
  }
  incomplete <- !stats::complete.cases(profiles[, required, drop = FALSE])
  if (any(incomplete)) {
    mostr_abort(
      "Fitted patient(s) have incomplete starting-profile predictors: ",
      markov_format_patient_ids(profile_ids[incomplete]),
      ". The transition response is not required, but ID, modeled predictors, ",
      "time metadata, and the previous state must be complete."
    )
  }

  fitting_data <- markov_model_data(model)
  factor_predictors <- intersect(
    metadata$predictor_variables %||% character(),
    names(profiles)[vapply(profiles, is.factor, logical(1))]
  )
  for (variable in factor_predictors) {
    if (is.null(fitting_data) || !variable %in% names(fitting_data)) {
      next
    }
    supported <- unique(as.character(fitting_data[[variable]]))
    unsupported <- !as.character(profiles[[variable]]) %in% supported
    if (any(unsupported)) {
      mostr_abort(
        "Fitted patient(s) have unsupported factor value(s) for starting ",
        "predictor `",
        variable,
        "`: ",
        markov_format_patient_ids(profile_ids[unsupported]),
        "."
      )
    }
  }

  unexpected_ids <- setdiff(profile_ids, fitted_ids)
  if (length(unexpected_ids)) {
    mostr_abort(
      "Stored starting profiles contain patient(s) with no usable likelihood ",
      "transition: ",
      markov_format_patient_ids(unexpected_ids),
      "."
    )
  }

  profiles <- profiles[match(fitted_ids, profile_ids), , drop = FALSE]
  rownames(profiles) <- NULL
  profiles
}

markov_subset_source_data <- function(data, subset, eval_env) {
  if (is.null(subset)) {
    return(list(data = data, rows = seq_len(nrow(data))))
  }

  selection <- eval(subset, data, eval_env)
  rows <- seq_len(nrow(data))[selection]
  rows <- rows[!is.na(rows)]
  list(data = data[rows, , drop = FALSE], rows = rows)
}

# Earlier fitting functions accepted these starting-profile arguments.
markov_reject_legacy_args <- function(dots) {
  legacy <- intersect(c("start_time", "origin_time", "time_map"), names(dots))
  if (!length(legacy)) {
    return(invisible(NULL))
  }

  mostr_abort(
    "mostm() does not accept the former starting-profile argument(s) ",
    paste0("`", legacy, "`", collapse = ", "),
    ". Use `first_followup_time` to select the stored starting profiles; ",
    "supply baseline and time-mapping arguments to downstream SOP summary ",
    "functions."
  )
}

markov_first_followup_time <- function(time, first_followup_time) {
  if (is.numeric(time)) {
    scheduled <- time[!is.na(time)]
    if (!length(scheduled)) {
      mostr_abort(
        "The numeric time variable has no non-missing scheduled values."
      )
    }
    if (any(!is.finite(scheduled))) {
      mostr_abort("The numeric time variable must contain only finite values.")
    }
    if (any(scheduled < 1)) {
      mostr_abort(
        "The numeric time variable contains a value below 1; follow-up time ",
        "must begin at 1 or later."
      )
    }

    if (is.null(first_followup_time)) {
      if (!1 %in% scheduled) {
        mostr_abort(
          "The default `first_followup_time = 1` is not observed in the ",
          "numeric time variable."
        )
      }
      return(1)
    }
    if (
      !is.numeric(first_followup_time) ||
        length(first_followup_time) != 1L ||
        is.na(first_followup_time) ||
        !is.finite(first_followup_time)
    ) {
      mostr_abort(
        "For numeric time, `first_followup_time` must be one finite numeric ",
        "value."
      )
    }
    if (!first_followup_time %in% scheduled) {
      mostr_abort(
        "`first_followup_time = ",
        as.character(first_followup_time),
        "` is not observed in the numeric time variable."
      )
    }
    earliest <- min(scheduled)
    if (!identical(as.numeric(first_followup_time), as.numeric(earliest))) {
      mostr_abort(
        "For numeric time, `first_followup_time` must equal the earliest ",
        "observed scheduled time (`",
        as.character(earliest),
        "`)."
      )
    }
    return(first_followup_time)
  }

  if (!is.factor(time) && !is.character(time)) {
    mostr_abort(
      "The time variable must be numeric, factor, or character."
    )
  }
  if (is.null(first_followup_time)) {
    mostr_abort(
      "Factor and character time require an explicit `first_followup_time`."
    )
  }
  if (length(first_followup_time) != 1L || is.na(first_followup_time)) {
    mostr_abort(
      "`first_followup_time` must contain one non-missing value."
    )
  }

  resolved <- as.character(first_followup_time)
  observed <- unique(as.character(time[!is.na(time)]))
  if (!resolved %in% observed) {
    mostr_abort(
      "`first_followup_time = ",
      resolved,
      "` is not observed in the factor or character time variable."
    )
  }
  resolved
}

markov_prepare_stored_data <- function(
  data,
  formula,
  subset,
  eval_env,
  fit_data,
  id_var,
  time_var,
  first_followup_time,
  starting_profiles = TRUE
) {
  if (!is.data.frame(data)) {
    return(list(
      refit_data = NULL,
      starting_profile_data = NULL,
      starting_profile_metadata = NULL
    ))
  }
  source <- markov_subset_source_data(data, subset, eval_env)
  refit_data <- source$data
  if (!starting_profiles) {
    return(list(
      refit_data = refit_data,
      starting_profile_data = NULL,
      starting_profile_metadata = NULL
    ))
  }
  scheduled_start <- markov_first_followup_time(
    refit_data[[time_var]],
    first_followup_time
  )
  validate_markov_id_var(id_var, fit_data, "fitted data")

  fitted_ids <- unique(as.character(fit_data[[id_var]]))
  if (anyNA(fitted_ids) || any(fitted_ids == "")) {
    mostr_abort("Fitted patient IDs must be non-missing and non-empty.")
  }
  # The fitted cohort is the patients with at least one row in the fit. A
  # patient whose every row the engine dropped (a model variable is missing on
  # each) is left out of the stored data as well, so that bootstrap resamples,
  # and the cohort each draw standardises over, hold the same patients as the
  # point estimate. Bootstrap refits are unchanged: the engine dropped these
  # rows anyway. Patients excluded by `subset` were removed above and are not
  # counted as dropped.
  in_fit <- as.character(refit_data[[id_var]]) %in% fitted_ids
  dropped_ids <- unique(as.character(refit_data[[id_var]][!in_fit]))
  refit_data <- refit_data[in_fit, , drop = FALSE]
  source_rows <- source$rows[in_fit]
  at_start <- as.character(refit_data[[time_var]]) %in%
    as.character(scheduled_start)
  starting_profile_data <- refit_data[at_start, , drop = FALSE]
  # Stored row weights are fitting metadata, not prediction variables.
  starting_profile_data[[markov_prior_weight_column()]] <- NULL
  # A bookkeeping column, so it starts with `.markov_`; see
  # `markov_internal_column()`.
  starting_profile_data$.markov_source_row <- source_rows[at_start]

  predictor_terms <- stats::delete.response(stats::terms(formula, data = data))
  predictor_variables <- all.vars(predictor_terms)
  required <- unique(c(
    id_var,
    predictor_variables,
    time_var
  ))

  list(
    refit_data = refit_data,
    dropped_ids = dropped_ids,
    starting_profile_data = starting_profile_data,
    starting_profile_metadata = list(
      id_var = id_var,
      time_var = time_var,
      first_followup_time = scheduled_start,
      required_variables = required,
      predictor_variables = predictor_variables,
      fitted_ids = fitted_ids,
      factor_levels = lapply(
        refit_data[vapply(refit_data, is.factor, logical(1))],
        levels
      )
    )
  )
}

# Columns that mostr adds to the stored starting profiles for its own
# bookkeeping start with `.markov_`, such as `.markov_source_row` (the row of
# the source data each profile came from). They are not patient variables, so
# no result shows them, whichever function or inference method made it.
# Results match profiles to patients through the ID column and `rowid`, not
# through these columns.
markov_internal_column <- function(names) {
  startsWith(names, ".markov_")
}

drop_markov_internal_columns <- function(data) {
  for (col in names(data)[markov_internal_column(names(data))]) {
    data[[col]] <- NULL
  }
  data
}

markov_align_model_data <- function(data, model_or_frame) {
  if (!is.data.frame(data)) {
    return(NULL)
  }

  rn <- NULL
  if (is.data.frame(model_or_frame)) {
    rn <- rownames(model_or_frame)
  } else {
    mf <- tryCatch(stats::model.frame(model_or_frame), error = function(e) NULL)
    if (!is.null(mf)) {
      rn <- rownames(mf)
    }
  }

  if (!is.null(rn) && length(rn) > 0 && all(rn %in% rownames(data))) {
    return(data[rn, , drop = FALSE])
  }

  model_n <- tryCatch(stats::nobs(model_or_frame), error = function(e) NULL)
  if (!is.null(model_n) && nrow(data) != model_n) {
    model_vars <- tryCatch(
      all.vars(stats::formula(model_or_frame)),
      error = function(e) character()
    )
    model_vars <- intersect(model_vars, names(data))
    if (length(model_vars) > 0) {
      keep <- stats::complete.cases(data[, model_vars, drop = FALSE])
      if (sum(keep) == model_n) {
        return(data[keep, , drop = FALSE])
      }
    }
  }

  data
}

# The data that SOP calls start from: the supplied `newdata`, or the starting
# profiles stored by `mostm()`. With `settings$absorb`, profiles that start in
# an absorbing state are refused here, before anything is predicted, so the
# error counts the user's profiles rather than the stacked prediction rows.
# Without `settings` (bare variable names look up values) nothing is checked.
resolve_markov_source_data <- function(model, newdata, settings = NULL) {
  time_var <- settings$time_var
  newdata_supplied <- !is.null(newdata)
  refit_data <- markov_model_refit_data(model)
  source_data <- if (newdata_supplied) {
    newdata
  } else {
    metadata <- markov_model_starting_profile_metadata(model)
    if (
      !is.null(time_var) &&
        !is.null(metadata$time_var) &&
        !identical(time_var, metadata$time_var)
    ) {
      mostr_abort(
        "Automatic SOP prediction requested `time_var = \"",
        time_var,
        "\"`, but `mostm()` stored the starting profiles using `",
        metadata$time_var,
        "`. Refit with the intended `time_var` or supply explicit `newdata`."
      )
    }
    markov_validate_starting_profiles(
      model,
      required_variables = c(time_var, settings$p_var)
    )
  }

  if (!is.data.frame(source_data)) {
    mostr_abort("`newdata` must be a data frame.")
  }
  if (!is.null(settings)) {
    stop_absorbing_start(
      source_data,
      settings,
      source = if (newdata_supplied) "newdata" else "stored"
    )
  }

  list(
    source_data = source_data,
    refit_data = refit_data,
    newdata_supplied = newdata_supplied
  )
}

resolve_markov_prediction_data <- function(
  data,
  id_var,
  time_var,
  data_label = "newdata"
) {
  if (is.null(id_var) || !id_var %in% names(data)) {
    return(data)
  }

  id <- data[[id_var]]
  if (!anyDuplicated(id)) {
    return(data)
  }

  if (time_var %in% names(data)) {
    work <- data.table::data.table(
      group = factor(id),
      time = data[[time_var]],
      row = seq_along(id)
    )
    data.table::setorderv(work, c("group", "time"), na.last = TRUE)
    baseline_idx <- work$row[!is.na(work$group) & !duplicated(work$group)]
    return(data[baseline_idx, , drop = FALSE])
  }

  warning(
    data_label,
    " contains repeated IDs but no `",
    time_var,
    "` column; using the first row per ID as the prediction row.",
    call. = FALSE
  )
  data[!duplicated(id), , drop = FALSE]
}

ensure_markov_rowid <- function(data) {
  data$rowid <- seq_len(nrow(data))
  data
}
