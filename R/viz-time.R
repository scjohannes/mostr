# Average time-in-state plotting.

#' Plot Average Time in States
#'
#' Plots output from [avg_time()], either as stacked horizontal bars that
#' divide follow-up into the average time spent in each state, or as points
#' with confidence intervals.
#'
#' @param x A `markov_avg_time` object returned by [avg_time()], optionally
#'   after [inferences()].
#' @param type Character string. The kind of plot:
#'   \itemize{
#'     \item `"bar"` (the default): one horizontal bar per scenario, split into
#'       the average time in each state set ("Grotta bars"). Use it to show
#'       how follow-up is shared between states. Confidence intervals cannot
#'       be shown on stacked segments, so this type ignores them.
#'     \item `"pointrange"`: states along the x axis, with one point per
#'       scenario and a vertical confidence interval when `x` has `conf.low`
#'       and `conf.high`. Scenarios are drawn in black and successively
#'       lighter shades of gray, so the plot prints well. Use it to show the
#'       uncertainty of each average.
#'   }
#' @param group_var Character vector or `NULL`. Columns that identify the
#'   scenarios being compared, for example the treatment variable set by
#'   `variables`. They label the bars (`type = "bar"`) or set the point colors
#'   (`type = "pointrange"`). If `NULL`, every column other than `state_set`,
#'   the estimates, and those in `facet_var` is used. Several columns are
#'   combined into one label.
#' @param facet_var Character vector or `NULL`. Optional faceting variable,
#'   for example a `by` variable. Use a length-two vector for
#'   `facet_grid(row ~ column)`.
#' @param scale Character string. `"time"` (the default) plots average time in
#'   the units of `x`. `"proportion"` divides each average by the scenario's
#'   total follow-up, so a bar spans 0 to 1. `"proportion"` requires every
#'   state to belong to exactly one state set.
#' @param show_labels Logical. For `type = "bar"`, print each segment's value
#'   inside it. Segments shorter than 6% of the longest bar are left unlabeled
#'   because the label would not fit.
#' @param digits Number of decimal places for segment labels.
#' @param show_uncertainty Logical. For `type = "pointrange"`, draw confidence
#'   intervals when `conf.low` and `conf.high` exist.
#' @param bar_width Numeric bar thickness for `type = "bar"`.
#' @param point_size Numeric point size for `type = "pointrange"`.
#' @param line_width Numeric line width for interval lines.
#'
#' @return A ggplot object.
#'
#' @details
#' A stacked bar only makes sense when no state is counted twice. With
#' `type = "bar"`, the state sets in `x` must therefore not overlap, which
#' holds for the default `state_sets = NULL` of [avg_time()] (every state
#' separately). When the sets also include every state, each bar's length is
#' the total follow-up: the number of visits, or the elapsed time covered when
#' `avg_time()` used `time_map`.
#'
#' To compare scenarios with a confidence interval for the difference, use
#' `avg_comparisons(estimand = "time_in_state")` and [plot_comparisons()].
#'
#' @examples
#' \dontrun{
#' tis <- avg_time(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:28,
#'   absorb = "6"
#' ) |>
#'   inferences(method = "delta")
#'
#' # How the 28 days are divided between states in each arm
#' plot_time(tis)
#'
#' # Average days in each state with 95% confidence intervals
#' plot_time(tis, type = "pointrange")
#' }
#'
#' @export
plot_time <- function(
  x,
  type = c("bar", "pointrange"),
  group_var = NULL,
  facet_var = NULL,
  scale = c("time", "proportion"),
  show_labels = TRUE,
  digits = 1,
  show_uncertainty = TRUE,
  bar_width = 0.6,
  point_size = 2,
  line_width = 0.7
) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  plot_validate_flag(show_labels, "show_labels")
  plot_validate_flag(show_uncertainty, "show_uncertainty")
  plot_validate_scalar(digits, "digits", lower = 0)
  plot_validate_scalar(bar_width, "bar_width", lower = 0, upper = 1)
  plot_validate_scalar(point_size, "point_size", lower = 0)
  plot_validate_scalar(line_width, "line_width", lower = 0)
  if (!inherits(x, "markov_avg_time")) {
    stop("`x` must be a `markov_avg_time` object returned by `avg_time()`.")
  }

  data <- as.data.frame(x)
  plot_validate_facets(data, facet_var)
  if (is.null(group_var)) {
    group_var <- setdiff(
      names(data),
      c("state_set", "time_unit", inference_columns(), facet_var)
    )
  } else {
    if (!is.character(group_var)) {
      stop("`group_var` must be NULL or a character vector.")
    }
    plot_validate_columns(data, group_var, "`group_var`")
  }

  state_sets <- attr(x, "time_args")$state_sets
  if (identical(type, "bar") || identical(scale, "proportion")) {
    plot_time_validate_partition(state_sets, attr(x, "y_levels"), scale)
  }

  data$state_set <- factor(
    as.character(data$state_set),
    levels = plot_comparisons_state_set_levels(x)
  )
  data$.group <- plot_time_group_labels(data, group_var)
  # One bar per scenario within each facet.
  data$.bar <- interaction(data[, c(".group", facet_var)], drop = TRUE)
  if (identical(scale, "proportion")) {
    data <- plot_time_proportions(data)
  }

  p <- if (identical(type, "bar")) {
    plot_time_bars(data, show_labels, digits, bar_width)
  } else {
    plot_time_points(data, show_uncertainty, point_size, line_width)
  }

  group_title <- plot_time_group_title(group_var)
  x_label <- plot_time_axis_label(scale, unique(data$time_unit))
  p <- p +
    if (identical(type, "bar")) {
      ggplot2::labs(x = x_label, y = group_title, fill = "State")
    } else {
      ggplot2::labs(x = "State", y = x_label, colour = group_title)
    }
  if (identical(type, "bar")) {
    p <- plot_add_default_scales(p)
  }
  plot_add_facets(p, facet_var)
}

plot_validate_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.")
  }
  invisible(NULL)
}

# Bars stack state-set totals, so a state counted in two sets would lengthen
# the bar. Proportions additionally need every state, so that the bar total is
# the follow-up length.
plot_time_validate_partition <- function(state_sets, y_levels, scale) {
  if (is.null(state_sets)) {
    stop(
      "`x` has no stored state sets; recreate it with `avg_time()`.",
      call. = FALSE
    )
  }
  states <- unlist(state_sets, use.names = FALSE)
  if (anyDuplicated(states)) {
    stop(
      "Stacked bars need state sets that do not overlap; state(s) ",
      paste(unique(states[duplicated(states)]), collapse = ", "),
      " appear in more than one set. Use `type = \"pointrange\"` or ",
      "non-overlapping `state_sets` in `avg_time()`.",
      call. = FALSE
    )
  }
  if (identical(scale, "proportion")) {
    missing_states <- setdiff(as_state_labels(y_levels), states)
    if (length(missing_states)) {
      stop(
        "`scale = \"proportion\"` needs every state in a state set; ",
        "missing state(s): ",
        paste(missing_states, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

plot_time_group_labels <- function(data, group_var) {
  if (!length(group_var)) {
    return(factor("All patients"))
  }
  if (length(group_var) == 1) {
    values <- data[[group_var]]
    return(if (is.factor(values)) values else factor(values))
  }
  interaction(data[, group_var, drop = FALSE], sep = " / ", lex.order = TRUE)
}

plot_time_group_title <- function(group_var) {
  if (!length(group_var)) {
    return(NULL)
  }
  paste(group_var, collapse = " / ")
}

plot_time_proportions <- function(data) {
  total <- stats::ave(data$estimate, data$.bar, FUN = sum)
  for (col in intersect(c("estimate", "conf.low", "conf.high"), names(data))) {
    data[[col]] <- data[[col]] / total
  }
  data
}

plot_time_axis_label <- function(scale, time_unit) {
  if (identical(scale, "proportion")) {
    return("Proportion of follow-up")
  }
  time_unit <- time_unit[!is.na(time_unit)]
  if (length(time_unit) == 1) {
    return(paste0("Average time in state (", time_unit, ")"))
  }
  "Average time in state"
}

plot_time_bars <- function(data, show_labels, digits, bar_width) {
  # Reverse the y axis so the first scenario is on top, and the stacking so the
  # first state starts at zero.
  data$.group <- factor(data$.group, levels = rev(levels(data$.group)))
  stack <- ggplot2::position_stack(reverse = TRUE)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[["estimate"]],
      y = .data[[".group"]],
      fill = .data[["state_set"]]
    )
  ) +
    ggplot2::geom_col(width = bar_width, colour = "white", position = stack)
  if (show_labels) {
    # Segments narrower than 6% of the longest bar cannot fit a label.
    width_limit <- 0.06 * max(stats::ave(data$estimate, data$.bar, FUN = sum))
    data$.label <- ifelse(
      data$estimate >= width_limit,
      formatC(data$estimate, format = "f", digits = digits),
      ""
    )
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[[".label"]]),
        data = data,
        position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE),
        size = 3
      )
  }
  p
}

plot_time_points <- function(data, show_uncertainty, point_size, line_width) {
  dodge <- ggplot2::position_dodge(width = 0.5)
  # Scenarios run from black to mid gray. Stopping at 65% gray keeps the last
  # scenario visible on the light grey panel.
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[["state_set"]],
      y = .data[["estimate"]],
      colour = .data[[".group"]],
      group = .data[[".group"]]
    )
  )
  if (show_uncertainty && all(c("conf.low", "conf.high") %in% names(data))) {
    p <- p +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data[["conf.low"]], ymax = .data[["conf.high"]]),
        linewidth = line_width,
        position = dodge
      )
  }
  p <- p +
    ggplot2::geom_point(size = point_size, position = dodge) +
    ggplot2::scale_colour_grey(start = 0, end = 0.65)
  if (nlevels(data$.group) == 1) {
    p <- p + ggplot2::guides(colour = "none")
  }
  p
}
