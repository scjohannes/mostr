# Average comparison setup helpers.

stop_legacy_baseline_args <- function(extra_args) {
  legacy_args <- intersect(names(extra_args), c("origin_time", "origin"))
  if (length(legacy_args) > 0L) {
    stop(
      "`",
      legacy_args[1L],
      "` is no longer supported; use `baseline_time` instead.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

validate_comparison_times <- function(times, fn) {
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied to `", fn, "`.")
  }
  if (anyDuplicated(times)) {
    stop(
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
# validates the list.
resolve_comparison_variables <- function(variables, model, newdata, fn) {
  if (!missing(variables) && is.character(variables)) {
    if (length(variables) != 1L || is.na(variables) || !nzchar(variables)) {
      stop("`", fn, "` supports exactly one named variable.")
    }
    data <- bare_variable_data(model, newdata)
    if (!variables %in% names(data)) {
      stop("Variables not in data: ", variables)
    }
    variables <- stats::setNames(
      list(comparison_variable_values(data[[variables]], variables)),
      variables
    )
  }
  validate_avg_comparison_variables(variables, fn)
}

# Turns bare variable names for predictions (`sops()`, `avg_sops()`,
# `state_time()`, `avg_state_time()`) into a named list of values. `NULL`
# and named lists are returned unchanged.
resolve_prediction_variables <- function(variables, model, newdata) {
  if (is.null(variables) || !is.character(variables)) {
    return(variables)
  }
  if (anyNA(variables) || any(!nzchar(variables)) || anyDuplicated(variables)) {
    stop("Bare variable names must be unique, non-missing, and non-empty.")
  }
  data <- bare_variable_data(model, newdata)
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0L) {
    stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }
  stats::setNames(
    lapply(variables, function(name) {
      prediction_variable_values(data[[name]], name)
    }),
    variables
  )
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
    stop(
      "A bare variable name supports factor, character, logical, or numeric ",
      "variables; ",
      explicit,
      call. = FALSE
    )
  }
  unique_values <- sort(unique(observed))
  if (all(unique_values %in% c(0, 1))) {
    return(unique_values)
  }
  if (identical(numeric, "fivenum")) {
    return(unique(stats::fivenum(observed)))
  }
  stop(
    "`",
    name,
    "` is numeric with values other than 0 and 1, so a bare name cannot ",
    "choose values to compare; ",
    explicit,
    call. = FALSE
  )
}

# The first value is the reference of every comparison.
comparison_variable_values <- function(x, name) {
  values <- bare_variable_values(x, name, numeric = "error")
  if (length(values) < 2L) {
    stop(
      "`",
      name,
      "` needs at least two observed values to compare.",
      call. = FALSE
    )
  }
  values
}

prediction_variable_values <- function(x, name) {
  values <- bare_variable_values(x, name, numeric = "fivenum")
  if (length(values) < 1L) {
    stop("`", name, "` needs at least one observed value.", call. = FALSE)
  }
  values
}

validate_avg_comparison_variables <- function(
  variables,
  fn = "avg_comparisons()"
) {
  if (missing(variables) || is.null(variables)) {
    stop("`variables` must be a named list with one counterfactual variable.")
  }
  if (!is.list(variables) || is.null(names(variables))) {
    stop("`variables` must be a named list, for example `list(tx = c(0, 1))`.")
  }
  if (length(variables) != 1L || !nzchar(names(variables)[1])) {
    stop("`", fn, "` supports exactly one named variable.")
  }

  values <- variables[[1]]
  if (length(values) < 2L) {
    stop("`variables[[1]]` must contain at least two values.")
  }
  if (anyDuplicated(as.character(values))) {
    stop("`variables[[1]]` must not contain duplicate values.")
  }

  variables
}

validate_avg_comparison_estimand <- function(estimand, state_sets, comparison) {
  if (estimand == "time_benefit") {
    if (!is.null(state_sets)) {
      stop("`state_sets` is not used with `estimand = \"time_benefit\"`.")
    }
    if (comparison != "difference") {
      stop(
        "`estimand = \"time_benefit\"` only supports `comparison = \"difference\"`."
      )
    }
  }
  invisible(NULL)
}

avg_comparison_setup <- function(
  model,
  newdata,
  variables,
  by,
  times,
  y_levels,
  absorb,
  time_var,
  p_var,
  id_var,
  p2_var,
  gap_var,
  time_covariates,
  include_re,
  ...
) {
  validate_markov_model(model)
  markov_reject_ml_second_order(model, p2_var)

  data_res <- resolve_markov_source_data(model, newdata)
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied

  if (inherits(model, "blrm") && isTRUE(include_re)) {
    id_var <- resolve_blrm_id_var(model, newdata_orig, id_var)
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  } else {
    id_var <- markov_model_id_var(model, id_var) %||% "id"
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
    stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  baseline_data <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = time_var
    )
  }
  validate_sops_by(by, baseline_data)
  baseline_data <- ensure_markov_rowid(baseline_data)

  grid <- create_counterfactual_grid(variables)
  newdata_pred <- create_counterfactual_data(baseline_data, grid, variables)

  list(
    variables = variables,
    by = by,
    times = times,
    time_var = time_var,
    p_var = p_var,
    id_var = id_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
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
    id_var = avg_args$id_var,
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
    time_var = attr(x, "time_var"),
    p_var = attr(x, "p_var"),
    p2_var = attr(x, "p2_var"),
    y_levels = attr(x, "y_levels"),
    absorb = attr(x, "absorb"),
    gap_var = attr(x, "gap_var"),
    time_covariates = attr(x, "time_covariates")
  )
}

avg_comparison_replay_avg_sops <- function(
  model,
  newdata,
  variables,
  by,
  times,
  y_levels,
  absorb,
  time_var,
  p_var,
  id_var,
  p2_var,
  gap_var,
  time_covariates,
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
      # Stored averages keep `by = NULL` for the overall average, which
      # `avg_sops()` spells `by = TRUE`.
      by = by %||% TRUE,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      id_var = id_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
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
    return(muffle_skipped_sop_visits(do.call(avg_sops, args)))
  }
  do.call(avg_sops, args)
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
    stop(
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
      times = setup$times,
      id_var = setup$id_var
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
    id_var = setup$id_var
  )

  for (nm in names(attrs)) {
    if (!is.null(attrs[[nm]])) {
      attr(result, nm) <- attrs[[nm]]
    }
  }

  for (nm in c(
    "call_args",
    "time_var",
    "p_var",
    "p2_var",
    "y_levels",
    "absorb",
    "gap_var",
    "time_covariates"
  )) {
    if (!is.null(setup[[nm]])) {
      attr(result, nm) <- setup[[nm]]
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
  result
}
