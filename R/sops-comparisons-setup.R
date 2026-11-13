# Average comparison setup helpers.


validate_comparison_times <- function(times, fn) {
  if (missing(times) || is.null(times)) {
    mostr_abort("`times` must be supplied to `", fn, "`.")
  }
  if (anyDuplicated(times)) {
    mostr_abort(
      "`times` must not contain duplicate visits; each visit is counted once."
    )
  }
  invisible(NULL)
}

# Data in which bare variable names are looked up: the data stored with the
# model, or `newdata` for a model without stored data, as marginaleffects
# uses the model's data.
bare_variable_data <- function(model, newdata) {
  data_res <- resolve_markov_source_data(model, newdata)
  data_res$refit_data %||% data_res$source_data
}

# Turns a bare variable name into the named list of values to compare, then
# validates the list. Values of the previous-state columns are matched to
# states by label (`match_previous_state_variables()`).
resolve_comparison_variables <- function(
  variables,
  model,
  newdata,
  fn,
  settings
) {
  if (!missing(variables) && is.character(variables)) {
    if (length(variables) != 1L || is.na(variables) || !nzchar(variables)) {
      mostr_abort("`", fn, "` supports exactly one named variable.")
    }
    data <- bare_variable_data(model, newdata)
    if (!variables %in% names(data)) {
      mostr_abort("Variables not in data: ", variables)
    }
    variables <- stats::setNames(
      list(comparison_variable_values(data[[variables]], variables)),
      variables
    )
  }
  variables <- validate_avg_comparison_variables(variables, fn)
  match_previous_state_variables(variables, model, newdata, settings)
}

# Turns bare variable names for predictions (`sops()`, `avg_sops()`,
# `state_time()`, `avg_state_time()`) into a named list of values. Named
# lists are kept, except that values of the previous-state columns are
# matched to states by label (`match_previous_state_variables()`).
resolve_prediction_variables <- function(variables, model, newdata, settings) {
  if (is.null(variables)) {
    return(variables)
  }
  if (is.character(variables)) {
    if (
      anyNA(variables) || any(!nzchar(variables)) || anyDuplicated(variables)
    ) {
      mostr_abort(
        "Bare variable names must be unique, non-missing, and non-empty."
      )
    }
    data <- bare_variable_data(model, newdata)
    missing_vars <- setdiff(variables, names(data))
    if (length(missing_vars) > 0L) {
      mostr_abort(
        "Variables not in data: ",
        paste(missing_vars, collapse = ", ")
      )
    }
    variables <- stats::setNames(
      lapply(variables, function(name) {
        prediction_variable_values(data[[name]], name)
      }),
      variables
    )
  }
  match_previous_state_variables(variables, model, newdata, settings)
}

# State labels are names: when `variables` sets the previous state (`p_var`)
# or the state two visits back (`p2_var`), a value means the state whose
# label it matches, so 2, 2L and "2" are all state "2". A value that matches
# no state is an error. Matched values are converted to the type of that
# column in the data the model was fitted to (or `newdata`, or the stored
# starting profiles), so every later step (the native engine, the reference
# engine, inference replays) sees values the model's design can use. A
# factor keeps its levels; a state that is not among them, such as an
# absorbing state that no transition starts from, is added as a level, so
# the absorbing-start refusal can still name it.
match_previous_state_variables <- function(
  variables,
  model,
  newdata,
  settings
) {
  if (!is.list(variables) || is.null(names(variables))) {
    return(variables)
  }
  p_var <- settings$p_var
  columns <- intersect(c(p_var, settings$p2_var), names(variables))
  if (length(columns) == 0L) {
    return(variables)
  }
  states <- tryCatch(
    as_state_labels(markov_model_ylevels(model, settings$y_levels)),
    error = function(e) NULL
  )
  for (column in columns) {
    values <- variables[[column]]
    labels <- as_state_labels(values)
    if (!is.null(states)) {
      unknown <- unique(labels[is.na(labels) | !labels %in% states])
      if (length(unknown) > 0L) {
        role <- if (identical(column, p_var)) {
          "the previous state"
        } else {
          "the state two visits back"
        }
        mostr_abort(
          "`",
          column,
          "` is ",
          role,
          ", so `variables` must set it to states. ",
          paste(unknown, collapse = ", "),
          if (length(unknown) == 1L) " is" else " are",
          " not a state; the states are ",
          paste(states, collapse = ", "),
          "."
        )
      }
    }
    prototype <- previous_state_prototype(model, newdata, column)
    if (is.null(prototype)) {
      next
    }
    variables[[column]] <- if (is.factor(prototype)) {
      factor(
        labels,
        levels = union(levels(prototype), intersect(states, labels)),
        ordered = is.ordered(prototype)
      )
    } else {
      coerce_previous_state_values(labels, prototype, column)
    }
  }
  variables
}

# The previous-state column whose type `variables` values take: from the
# data the model was fitted to, else `newdata`, else the starting profiles
# stored by `mostm()`.
previous_state_prototype <- function(model, newdata, column) {
  sources <- list(
    markov_model_refit_data(model),
    newdata,
    markov_model_metadata_attr(model, "markov_starting_profile_data")
  )
  for (data in sources) {
    if (is.data.frame(data) && column %in% names(data)) {
      return(data[[column]])
    }
  }
  NULL
}

# Observed values of a variable named without values, chosen as in
# marginaleffects: factor levels in level order, sorted character values,
# FALSE then TRUE, or 0 then 1. `numeric` decides other numeric variables:
# predictions use their five-number summary, as marginaleffects'
# predictions() does; comparisons stop, because marginaleffects' default for
# them (each row's value plus one) is not a comparison of fixed values.
bare_variable_values <- function(x, name, numeric = c("error", "fivenum")) {
  numeric <- match.arg(numeric)
  explicit <- paste0(
    "list the values instead, for example `variables = list(",
    name,
    " = c(...))`."
  )
  observed <- x[!is.na(x)]
  if (is.factor(x)) {
    return(factor(
      levels(x)[levels(x) %in% as.character(observed)],
      levels = levels(x)
    ))
  }
  if (is.logical(x) || is.character(x)) {
    return(sort(unique(observed)))
  }
  if (!is.numeric(x) || is.object(x)) {
    mostr_abort(
      "A bare variable name supports factor, character, logical, or numeric ",
      "variables; ",
      explicit
    )
  }
  unique_values <- sort(unique(observed))
  if (all(unique_values %in% c(0, 1))) {
    return(unique_values)
  }
  if (identical(numeric, "fivenum")) {
    return(unique(stats::fivenum(observed)))
  }
  mostr_abort(
    "`",
    name,
    "` is numeric with values other than 0 and 1, so a bare name cannot ",
    "choose values to compare; ",
    explicit
  )
}

# The first value is the reference of every comparison.
comparison_variable_values <- function(x, name) {
  values <- bare_variable_values(x, name, numeric = "error")
  if (length(values) < 2L) {
    mostr_abort(
      "`",
      name,
      "` needs at least two observed values to compare."
    )
  }
  values
}

prediction_variable_values <- function(x, name) {
  values <- bare_variable_values(x, name, numeric = "fivenum")
  if (length(values) < 1L) {
    mostr_abort("`", name, "` needs at least one observed value.")
  }
  values
}

validate_avg_comparison_variables <- function(
  variables,
  fn = "avg_comparisons()"
) {
  if (missing(variables) || is.null(variables)) {
    mostr_abort(
      "`variables` must be a named list with one counterfactual variable."
    )
  }
  if (!is.list(variables) || is.null(names(variables))) {
    mostr_abort(
      "`variables` must be a named list, for example `list(tx = c(0, 1))`."
    )
  }
  if (length(variables) != 1L || !nzchar(names(variables)[1])) {
    mostr_abort("`", fn, "` supports exactly one named variable.")
  }

  values <- variables[[1]]
  if (length(values) < 2L) {
    mostr_abort("`variables[[1]]` must contain at least two values.")
  }
  if (anyDuplicated(as.character(values))) {
    mostr_abort("`variables[[1]]` must not contain duplicate values.")
  }

  variables
}

validate_avg_comparison_estimand <- function(estimand, state_sets, comparison) {
  if (estimand == "time_benefit") {
    if (!is.null(state_sets)) {
      mostr_abort(
        "`state_sets` is not used with `estimand = \"time_benefit\"`."
      )
    }
    if (comparison != "difference") {
      mostr_abort(
        "`estimand = \"time_benefit\"` only supports `comparison = \"difference\"`."
      )
    }
  }
  invisible(NULL)
}

# The profiles of a comparison or time call, stacked once per value of
# `variables`, with the resolved `settings` of the call.
avg_comparison_setup <- function(
  model,
  newdata,
  variables,
  by,
  times,
  settings,
  include_re
) {
  data_res <- resolve_markov_source_data(model, newdata, settings)
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied

  id_var <- markov_require_id_var(model, settings$id_var)
  if (inherits(model, "blrm") && isTRUE(include_re)) {
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  }

  if (!is.null(refit_data)) {
    validate_markov_id_var(
      id_var,
      refit_data,
      "the refit data stored on the fitted model"
    )
  }
  if (!newdata_supplied) {
    validate_markov_id_var(id_var, newdata_orig, "stored model data")
  }

  missing_vars <- setdiff(names(variables), names(newdata_orig))
  if (length(missing_vars) > 0L) {
    mostr_abort("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  baseline_data <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = settings$time_var
    )
  }
  validate_sops_by(by, baseline_data)
  stop_missing_by_values(baseline_data, by, variables)
  baseline_data <- ensure_markov_rowid(baseline_data)

  grid <- create_counterfactual_grid(variables)
  newdata_pred <- create_counterfactual_data(baseline_data, grid, variables)

  list(
    variables = variables,
    by = by,
    times = times,
    settings = settings,
    newdata_orig = newdata_orig,
    refit_data = refit_data,
    newdata_supplied = newdata_supplied,
    baseline_data = baseline_data,
    grid = grid,
    newdata_pred = newdata_pred,
    n_cf = nrow(grid),
    n_each = nrow(baseline_data)
  )
}

avg_comparison_setup_from_sops <- function(x) {
  avg_args <- attr(x, "avg_args")
  variables <- avg_args$variables
  grid <- create_counterfactual_grid(variables)
  newdata_pred <- attr(x, "newdata_pred")
  n_each <- if (!is.null(newdata_pred) && nrow(grid) > 0L) {
    nrow(newdata_pred) / nrow(grid)
  } else {
    NA_real_
  }

  list(
    variables = variables,
    by = avg_args$by,
    times = avg_args$times,
    newdata_orig = attr(x, "newdata_orig"),
    refit_data = attr(x, "refit_data"),
    newdata_supplied = isTRUE(attr(x, "newdata_supplied")),
    baseline_data = if (!is.null(newdata_pred) && is.finite(n_each)) {
      newdata_pred[seq_len(n_each), , drop = FALSE]
    } else {
      NULL
    },
    grid = grid,
    newdata_pred = newdata_pred,
    n_cf = nrow(grid),
    n_each = as.integer(n_each),
    call_args = attr(x, "call_args"),
    settings = markov_result_settings(x)
  )
}

# Average SOPs again with the settings stored on a result, for inference and
# the reductions of comparisons and times.
avg_comparison_replay_avg_sops <- function(
  model,
  newdata,
  variables,
  by,
  times,
  settings,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  extra_args,
  warn_skipped_visits = TRUE
) {
  args <- c(
    list(
      model = model,
      newdata = newdata,
      variables = variables,
      by = by,
      times = times,
      settings = settings,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws
    ),
    extra_args
  )
  # Inference replays the user's request, which already warned about `times`.
  if (!isTRUE(warn_skipped_visits)) {
    return(muffle_skipped_sop_visits(do.call(avg_sops_impl, args)))
  }
  do.call(avg_sops_impl, args)
}

normalize_comparison_state_sets <- function(states, y_levels) {
  state_labels <- as_state_labels(y_levels)
  if (is.null(states)) {
    out <- as.list(state_labels)
    names(out) <- state_labels
    return(out)
  }

  if (is.list(states)) {
    out <- lapply(states, as_state_labels)
    if (is.null(names(out))) {
      names(out) <- vapply(out, state_set_label, character(1))
    }
  } else {
    out <- list(as_state_labels(states))
    names(out) <- state_set_label(out[[1]])
  }

  missing_states <- setdiff(
    unique(unlist(out, use.names = FALSE)),
    state_labels
  )
  if (length(missing_states) > 0L) {
    mostr_abort(
      "State(s) not found in SOP output: ",
      paste(missing_states, collapse = ", ")
    )
  }

  out
}

state_set_label <- function(states) {
  paste(as_state_labels(states), collapse = "+")
}


set_avg_comparison_attrs <- function(
  result,
  model,
  setup,
  estimand,
  state_sets,
  comparison,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  time_map,
  baseline_time,
  target_times,
  time_unit,
  extra_args,
  class_name = "markov_avg_comparisons"
) {
  result <- as.data.frame(result)
  attrs <- list(
    model = model,
    avg_args = list(
      variables = setup$variables,
      by = setup$by,
      times = setup$times
    ),
    comparison_args = list(
      estimand = estimand,
      state_sets = state_sets,
      comparison = comparison,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      time_map = time_map,
      baseline_time = baseline_time,
      target_times = target_times,
      time_unit = time_unit,
      extra_args = extra_args
    ),
    newdata_orig = setup$newdata_orig,
    newdata_pred = setup$newdata_pred,
    refit_data = setup$refit_data,
    newdata_supplied = setup$newdata_supplied,
    comparison_baseline_data = setup$baseline_data,
    comparison_grid = setup$grid,
    comparison_n_each = setup$n_each,
    call_args = setup$call_args,
    settings = setup$settings
  )

  for (nm in names(attrs)) {
    if (!is.null(attrs[[nm]])) {
      attr(result, nm) <- attrs[[nm]]
    }
  }

  if (inherits(model, "blrm")) {
    attr(result, "method") <- "posterior"
    attr(result, "engine") <- "posterior"
    attr(result, "conf_level") <- conf_level
    attr(result, "posterior_summary") <- posterior_summary
  }

  result <- order_estimate_columns(result)
  class(result) <- c(class_name, "data.frame")
  # Posterior intervals of blrm comparisons are returned here.
  drop_unused_comparison_uncertainty(result)
}
