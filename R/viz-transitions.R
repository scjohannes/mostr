# Transition heatmaps.

#' Plot Empirical or Model-Based Transition Proportions
#'
#' `plot_transitions()` creates heatmaps of population transition proportions
#' from observed trajectory data or from transition probabilities implied by a
#' fitted Markov model. Model-based plots are population averaged over the
#' supplied or stored patient profiles. When `comparison = "difference"`,
#' counterfactual transition summaries are computed for the levels in
#' `variables`, and the heatmap shows the difference in transition proportions
#' relative to the first level.
#'
#' @details
#' Minimum inputs depend on the object type. For trajectory data, `object` must
#' contain the current time, current state, and previous-state columns named by
#' `time_var`, `y_var`, and `p_var`. For model-based plots, pass a fitted
#' Markov model and `times`; use `newdata` for explicit prediction profiles. The
#' column names and states default to those of the fit, and the absorbing
#' states are those given to [mostm()]. Difference plots also require
#' `comparison = "difference"` and `variables = list(var = c(reference,
#' comparison))`.
#'
#' @param x A trajectory data frame or a fitted Markov model.
#' @param newdata Optional data frame of prediction profiles for model-based
#'   plots. If `NULL`, models fitted with [mostm()] use their stored data and
#'   extract
#'   one prediction row per ID.
#' @param variables Optional named list with one counterfactual variable. For
#'   `comparison = "difference"`, at least two values are required and the
#'   first value is the reference.
#' @param times Optional current-time values to plot. Required for model-based
#'   plots.
#' @param y_levels The states in order. For a fitted model, `NULL` (the
#'   default) uses the states of the model; a supplied vector must hold the
#'   same states in the same order, otherwise it is an error. For trajectory
#'   data, `NULL` uses the states found in the `y_var` and `p_var` columns,
#'   in numeric order when they are numbers.
#' @param comparison `"none"` for transition proportions or `"difference"` for
#'   counterfactual differences in transition proportions.
#' @param time_var Name of the time column. For a fitted model, `NULL` (the
#'   default) uses the `time_var` given to [mostm()]; a different name is an
#'   error. For trajectory data, `NULL` means `"time"`.
#' @param y_var Character name of the current-state column in trajectory data.
#' @param p_var Name of the previous-state column. For a fitted model, `NULL`
#'   (the default) uses the `p_var` given to [mostm()], `"yprev"` unless set
#'   there; a different name is an error. For trajectory data, `NULL` means
#'   `"yprev"`.
#' @param p2_var Name of the column with the state two visits back, for
#'   model-based plots. `NULL` (the default) uses the `p2_var` given to
#'   [mostm()]: a model fitted with `p2_var` is plotted with a second-order
#'   recursion, and one fitted without it with a first-order recursion. A
#'   different name is an error, so to plot a second-order model, give
#'   `p2_var` to [mostm()]. Not used for trajectory data.
#' @param id_var Name of the patient ID column, used to take one starting row
#'   per patient from the data stored on a fitted model. `NULL` (the default)
#'   uses the `id_var` given to [mostm()]; a different name is an error. Not
#'   used for trajectory data.
#' @param facet_var Optional observed grouping variable. The plot always facets
#'   by time; `facet_var` adds grouping to the time panels.
#' @param seed Optional random seed for `blrm` posterior draw selection.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param show_values Logical. If `TRUE`, print rounded values in each tile.
#' @param digits Number of digits used for tile labels. Defaults to 2 for
#'   proportions and 3 for differences.
#' @param fill_limits Optional fill-scale limits. Difference plots use symmetric
#'   limits around zero when this is `NULL`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_transitions(data, time_var = "week", p_var = "yprev")
#'
#' plot_transitions(
#'   fit,
#'   times = 1:7,
#'   variables = list(tx = c(0, 1)),
#'   comparison = "difference"
#' )
#' }
#' @export
plot_transitions <- function(
  x,
  newdata = NULL,
  variables = NULL,
  times = NULL,
  y_levels = NULL,
  comparison = c("none", "difference"),
  time_var = NULL,
  y_var = "y",
  p_var = NULL,
  p2_var = NULL,
  id_var = NULL,
  facet_var = NULL,
  seed = NULL,
  n_draws = 100L,
  show_values = TRUE,
  digits = NULL,
  fill_limits = NULL
) {
  object <- x
  comparison <- match.arg(comparison)
  if (
    !is.logical(show_values) || length(show_values) != 1L || is.na(show_values)
  ) {
    mostr_abort("`show_values` must be TRUE or FALSE.")
  }
  if (!is.null(digits)) {
    plot_validate_scalar(digits, "digits", lower = 0)
    digits <- floor(digits)
  }
  if (is.null(digits)) {
    digits <- if (identical(comparison, "difference")) 3L else 2L
  }

  if (markov_supported_model(object)) {
    settings <- plot_markov_settings(
      object,
      id_var = id_var,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      y_levels = y_levels
    )
    data <- plot_transitions_model_data(
      model = object,
      newdata = newdata,
      variables = variables,
      times = times,
      comparison = comparison,
      settings = settings,
      facet_var = facet_var,
      seed = seed,
      n_draws = n_draws
    )
  } else {
    if (!is.data.frame(object)) {
      mostr_abort("`object` must be a data frame or a supported Markov model.")
    }
    if (!is.null(newdata)) {
      mostr_abort("`newdata` is only used for model-based plots.")
    }
    if (!is.null(variables)) {
      mostr_abort("`variables` is only used for model-based plots.")
    }
    if (!identical(comparison, "none")) {
      mostr_abort(
        "`comparison = \"difference\"` requires a fitted Markov model."
      )
    }
    # Trajectory data use only the time and previous-state columns and the
    # state levels; the other settings are ignored.
    settings <- plot_markov_settings(
      object,
      time_var = time_var,
      p_var = p_var,
      y_levels = y_levels
    )
    data <- plot_transitions_empirical_data(
      data = object,
      times = times,
      settings = settings,
      y_var = y_var,
      facet_var = facet_var
    )
  }

  plot_transitions_heatmap(
    data = data,
    comparison = comparison,
    show_values = show_values,
    digits = digits,
    fill_limits = fill_limits
  )
}

plot_transitions_heatmap <- function(
  data,
  comparison,
  show_values,
  digits,
  fill_limits
) {
  data$.label <- plot_value_label(data$estimate, digits)
  p <- ggplot2::ggplot(data) +
    ggplot2::aes(
      x = .data$previous_state,
      y = .data$state,
      fill = .data$estimate
    ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::facet_wrap(ggplot2::vars(.data$.panel)) +
    ggplot2::labs(
      x = "Previous State",
      y = "Current State",
      fill = if (identical(comparison, "difference")) {
        "Difference in transition proportion"
      } else {
        "Transition proportion"
      }
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      aspect.ratio = 1
    )

  if (show_values) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$.label),
        color = "black",
        size = 2.5
      )
  }

  if (identical(comparison, "difference")) {
    if (is.null(fill_limits)) {
      max_abs <- max(abs(data$estimate), na.rm = TRUE)
      if (!is.finite(max_abs) || max_abs == 0) {
        max_abs <- 1
      }
      fill_limits <- c(-max_abs, max_abs)
    }
    return(
      p +
        ggplot2::scale_fill_gradient2(
          low = "#d73027",
          mid = "white",
          high = "#2c7bb6",
          midpoint = 0,
          limits = fill_limits,
          na.value = "grey90"
        )
    )
  }

  p +
    ggplot2::scale_fill_viridis_c(
      limits = fill_limits,
      na.value = "grey90"
    )
}
