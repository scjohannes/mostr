# Correlation heatmaps and variograms.

#' Plot Empirical or Model-Implied Correlations Over Time
#'
#' `plot_correlation()` computes correlations between ordinal state scores at
#' pairs of follow-up times. With a fitted Markov model, correlations are
#' computed from model-implied pairwise moments rather than sampled paths.
#'
#' @details
#' Minimum inputs depend on the object type. For trajectory data, `object` must
#' contain the ID, time, and state columns named by `id_var`, `time_var`, and
#' `y_var`. For model-based plots, pass a fitted Markov model and `times`; use
#' `newdata` for explicit prediction profiles. The column names and states
#' default to those of the fit, and the absorbing states are those given to
#' [mostm()]. Model-based correlations are computed from exact first- or
#' second-order Markov recursions. For `blrm` fits, posterior draws are kept
#' separate through the moment calculation and summarized only after
#' correlations are computed for each draw.
#'
#' Spearman correlation is the default: it measures association between patients'
#' ranks at two times, assigning average ranks to tied states. Observed-data
#' calculations use patients observed at both times and rank each pair separately.
#' Model-based calculations use the corresponding ranks from each time's state
#' probabilities, averaged over prediction profiles within each stratum. Ranks
#' are calculated separately for each posterior draw before summarizing.
#' Correlations are undefined when either time has no variation in state.
#'
#' With `method = "pearson"`, correlations describe numerical state scores.
#' Model-based scores are positions in `y_levels` (1 through the number of states).
#' Observed scores use `y_levels` positions when supplied, factor codes for
#' factors, or integer state values. Use the same `y_levels` for comparisons.
#' Pearson is useful when differences between these scores are meaningful.
#'
#' @param method Correlation method: `"spearman"` (default) or `"pearson"`.
#' @param x A trajectory data frame or a fitted Markov model.
#' @param newdata Optional data frame of prediction profiles for model-based
#'   plots. If `NULL`, models fitted with [mostm()] use their stored data and
#'   extract
#'   one prediction row per ID.
#' @param times Optional time values to include. Required for model-based plots.
#' @param y_levels The states in order. For a fitted model, `NULL` (the
#'   default) uses the states of the model; a supplied vector must hold the
#'   same states in the same order, otherwise it is an error. For trajectory
#'   data, supplied states set the order of the state scores; with `NULL`, the
#'   scores are factor codes or the integer state values.
#' @param id_var Name of the patient ID column. For a fitted model, `NULL`
#'   (the default) uses the `id_var` given to [mostm()], which selects one
#'   starting row per patient from the stored data; a different name is an
#'   error. For trajectory data, `NULL` means `"id"`.
#' @param time_var Name of the time column. For a fitted model, `NULL` (the
#'   default) uses the `time_var` given to [mostm()]; a different name is an
#'   error. For trajectory data, `NULL` means `"time"`.
#' @param y_var Character name of the state column in raw trajectory data.
#' @param p_var Name of the previous-state column for model recursions.
#'   `NULL` (the default) uses the `p_var` given to [mostm()], `"yprev"`
#'   unless set there; a different name is an error. Not used for trajectory
#'   data.
#' @param p2_var Name of the column with the state two visits back, for
#'   model-based plots. `NULL` (the default) uses the `p2_var` given to
#'   [mostm()]: a model fitted with `p2_var` uses a second-order recursion,
#'   and one fitted without it a first-order recursion. A different name is
#'   an error, so to plot a second-order model, give `p2_var` to [mostm()].
#'   Supplying it for trajectory data is an error.
#' @param facet_var Optional grouping variable. Correlations are computed
#'   separately within each observed stratum.
#' @param seed Optional random seed for `blrm` posterior draw selection.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param triangle `"upper"` to show one triangle or `"full"` to show both
#'   off-diagonal triangles.
#' @param show_values Logical. If `TRUE`, print rounded correlations in tiles.
#' @param digits Number of digits used for tile labels.
#' @param fill_limits Optional fill-scale limits. By default (`NULL`), use
#'   `c(0, 1)` when all available correlations are nonnegative, and `c(-1, 1)`
#'   when any are negative or all are missing. Supply limits explicitly to use
#'   the same color scale across several plots.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_correlation(data, id_var = "id", time_var = "day")
#' plot_correlation(fit, times = 1:14)
#' }
#' @export
plot_correlation <- function(
  x,
  newdata = NULL,
  times = NULL,
  y_levels = NULL,
  id_var = NULL,
  time_var = NULL,
  y_var = "y",
  p_var = NULL,
  p2_var = NULL,
  facet_var = NULL,
  seed = NULL,
  n_draws = 100L,
  triangle = c("upper", "full"),
  show_values = TRUE,
  digits = 2,
  fill_limits = NULL,
  method = c("spearman", "pearson")
) {
  method <- match.arg(method)
  object <- x
  triangle <- match.arg(triangle)
  if (
    !is.logical(show_values) || length(show_values) != 1L || is.na(show_values)
  ) {
    mostr_abort("`show_values` must be TRUE or FALSE.")
  }
  plot_validate_scalar(digits, "digits", lower = 0)
  settings <- plot_markov_settings(
    object,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    data_id_var = "id"
  )

  corr <- plot_correlation_input_data(
    object = object,
    newdata = newdata,
    times = times,
    settings = settings,
    y_var = y_var,
    facet_var = facet_var,
    seed = seed,
    n_draws = n_draws,
    triangle = triangle,
    method = method
  )

  plot_correlation_heatmap(
    corr,
    show_values = show_values,
    digits = floor(digits),
    fill_limits = fill_limits
  )
}

#' Plot an Empirical or Model-Implied Correlation Variogram
#'
#' `plot_variogram()` summarizes the same time-by-time correlation matrix used
#' by [plot_correlation()] as correlations against absolute time differences.
#' Its minimum input requirements are the same as [plot_correlation()].
#'
#' @inheritParams plot_correlation
#' @param smooth Logical. If `TRUE`, add a loess smooth.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_variogram(data, id_var = "id", time_var = "day")
#' }
#' @export
plot_variogram <- function(
  x,
  newdata = NULL,
  times = NULL,
  y_levels = NULL,
  id_var = NULL,
  time_var = NULL,
  y_var = "y",
  p_var = NULL,
  p2_var = NULL,
  facet_var = NULL,
  seed = NULL,
  n_draws = 100L,
  smooth = TRUE,
  method = c("spearman", "pearson")
) {
  method <- match.arg(method)
  object <- x
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth)) {
    mostr_abort("`smooth` must be TRUE or FALSE.")
  }
  settings <- plot_markov_settings(
    object,
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    data_id_var = "id"
  )

  corr <- plot_correlation_input_data(
    object = object,
    newdata = newdata,
    times = times,
    settings = settings,
    y_var = y_var,
    facet_var = facet_var,
    seed = seed,
    n_draws = n_draws,
    triangle = "upper",
    method = method
  )
  variogram <- plot_variogram_data(corr)

  p <- ggplot2::ggplot(variogram) +
    ggplot2::aes(x = .data$delta, y = .data$correlation) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(
      x = "Absolute Time Difference",
      y = "Correlation"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw()

  if (smooth) {
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )
  }
  if (".panel" %in% names(variogram)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.panel))
  }

  p
}

plot_correlation_heatmap <- function(
  data,
  show_values,
  digits,
  fill_limits
) {
  if (is.null(fill_limits)) {
    correlations <- data$correlation[is.finite(data$correlation)]
    fill_limits <- if (length(correlations) > 0L && all(correlations >= 0)) {
      c(0, 1)
    } else {
      c(-1, 1)
    }
  }
  data$.label <- plot_value_label(data$correlation, digits)
  p <- ggplot2::ggplot(data) +
    ggplot2::aes(x = .data$time_1, y = .data$time_2, fill = .data$correlation) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_viridis_c(limits = fill_limits, na.value = "grey90") +
    ggplot2::labs(x = "Time", y = "Time", fill = "") +
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
  if (".panel" %in% names(data)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.panel))
  }

  p
}
