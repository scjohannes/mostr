# The rows one layer of `plot` draws, with the data values its scales and
# panels show. `keys` maps an aesthetic ("x", "colour", "fill", "linetype") or
# a facet variable to the name of a new character column holding the decoded
# value, so the rows can be matched to the data they should show.
plot_layer_values <- function(plot, layer, keys) {
  # ggplot2 4.0 reports legend titles set for aesthetics that no layer maps.
  built <- suppressMessages(ggplot2::ggplot_build(plot))
  values <- built$data[[layer]]
  panels <- built$layout$layout
  for (aesthetic in names(keys)) {
    if (aesthetic %in% names(panels)) {
      decoded <- panels[[aesthetic]][match(values$PANEL, panels$PANEL)]
    } else if (aesthetic == "x") {
      x_scale <- built$layout$panel_params[[1]]$x
      decoded <- if (x_scale$is_discrete()) {
        x_scale$get_limits()[round(values$x)]
      } else {
        values$x
      }
    } else {
      scale <- built$plot$scales$get_scales(aesthetic)
      limits <- scale$get_limits()
      decoded <- limits[match(values[[aesthetic]], scale$map(limits))]
    }
    values[[keys[[aesthetic]]]] <- as.character(decoded)
  }
  if (all(c("ymin", "ymax") %in% names(values))) {
    values$height <- values$ymax - values$ymin
  }
  values
}

# Every plotted row matches one row of `expected` on the key columns, and the
# plotted columns (names of `columns`) equal the expected ones (its values).
expect_plotted_values <- function(plotted, expected, keys, columns, info) {
  keys <- unname(keys)
  expected[keys] <- lapply(expected[keys], as.character)
  joined <- merge(
    plotted[c(keys, names(columns))],
    expected[c(keys, unname(columns))],
    by = keys,
    suffixes = c(".plotted", ".expected")
  )
  expect_equal(nrow(plotted), nrow(expected), info = info)
  expect_equal(nrow(joined), nrow(expected), info = info)
  for (column in names(columns)) {
    plotted_column <- column
    expected_column <- columns[[column]]
    if (identical(plotted_column, expected_column)) {
      plotted_column <- paste0(column, ".plotted")
      expected_column <- paste0(column, ".expected")
    }
    expect_equal(
      joined[[plotted_column]],
      joined[[expected_column]],
      info = info
    )
  }
}

test_that("plot_sops plots empirical state proportions within each group", {
  data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, 4),
    y = c(1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 3, 3),
    tx = rep(c(0, 1), each = 6),
    site = rep(c("a", "b"), times = 6)
  )
  # Share of the rows in each group and time that are in each observed state,
  # counted directly from `data`.
  proportions <- function(by) {
    cells <- unique(data[c(by, "time", "y")])
    cells$prob <- vapply(
      seq_len(nrow(cells)),
      function(i) {
        same <- data$time == cells$time[i]
        for (column in by) {
          same <- same & data[[column]] == cells[[column]][i]
        }
        mean(data$y[same] == cells$y[i])
      },
      numeric(1)
    )
    names(cells)[names(cells) == "y"] <- "state"
    cells
  }

  cases <- list(
    list(
      label = "bars faceted by tx (default)",
      args = list(geom = "bar"),
      geom = "GeomBar",
      by = "tx",
      keys = c(x = "time", fill = "state", tx = "tx"),
      value = "height",
      n_panels = 2L,
      linetype = NULL
    ),
    list(
      label = "lines without facets (default geom)",
      args = list(facet_var = NULL),
      geom = "GeomLine",
      by = character(),
      keys = c(x = "time", colour = "state"),
      value = "y",
      n_panels = 1L,
      linetype = NULL
    ),
    list(
      label = "lines with tx as linetype",
      args = list(facet_var = NULL, linetype_var = "tx", geom = "line"),
      geom = "GeomLine",
      by = "tx",
      keys = c(x = "time", colour = "state", linetype = "tx"),
      value = "y",
      n_panels = 1L,
      linetype = "tx"
    ),
    list(
      label = "lines faceted by tx with site as linetype",
      args = list(facet_var = "tx", linetype_var = "site", geom = "line"),
      geom = "GeomLine",
      by = c("tx", "site"),
      keys = c(x = "time", colour = "state", linetype = "site", tx = "tx"),
      value = "y",
      n_panels = 2L,
      linetype = "site"
    )
  )

  for (case in cases) {
    label <- case$label
    plot <- do.call(plot_sops, c(list(data), case$args))
    plotted <- plot_layer_values(plot, 1L, case$keys)

    expect_length(plot$layers, 1L)
    expect_true(inherits(plot$layers[[1]]$geom, case$geom), info = label)
    expect_length(unique(plotted$PANEL), case$n_panels)
    expect_identical(
      plot$labels$title,
      "Empirical State Occupancy Probabilities Over Time",
      info = label
    )
    expect_identical(plot$labels$linetype, case$linetype, info = label)
    expect_true(plot$scales$has_scale("colour"), info = label)
    expect_true(plot$scales$has_scale("fill"), info = label)
    expect_plotted_values(
      plotted,
      proportions(case$by),
      case$keys,
      stats::setNames("prob", case$value),
      info = label
    )
  }
})

make_avg_sops_plot_data <- function(with_draws = TRUE) {
  data <- expand.grid(
    time = 1:2,
    state = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
    tx = c(0, 1),
    KEEP.OUT.ATTRS = FALSE
  )
  data$estimate <- c("1" = 0.2, "2" = 0.3, "3" = 0.5)[
    as.character(data$state)
  ]
  data$conf.low <- pmax(0, data$estimate - 0.05)
  data$conf.high <- pmin(1, data$estimate + 0.05)
  class(data) <- c("markov_avg_sops", class(data))

  if (with_draws) {
    draws <- expand.grid(
      draw_id = 1:7,
      time = 1:2,
      state = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
      tx = c(0, 1),
      KEEP.OUT.ATTRS = FALSE
    )
    draws$estimate <- c("1" = 0.18, "2" = 0.31, "3" = 0.51)[
      as.character(draws$state)
    ]
    attr(data, "draws") <- draws
  }

  data
}

test_that("plot_sops draws model-derived estimates and confidence limits", {
  # Distinct estimates, so a value drawn at the wrong time, state or arm shows.
  distinct <- function(data) {
    data$estimate <- seq_len(nrow(data)) / 100
    data$conf.low <- data$estimate - 0.005
    data$conf.high <- data$estimate + 0.01
    data
  }
  with_draws <- distinct(make_avg_sops_plot_data())
  without_draws <- distinct(make_avg_sops_plot_data(with_draws = FALSE))
  plain <- without_draws
  class(plain) <- "data.frame"
  plain$source <- rep(c("a", "b"), length.out = nrow(plain))

  line_keys <- c(x = "time", colour = "state", tx = "tx")
  fill_keys <- c(x = "time", fill = "state", tx = "tx")
  cases <- list(
    list(
      label = "avg_sops lines with ribbons",
      data = with_draws,
      args = list(geom = "line", facet_var = "tx"),
      warning = NULL,
      n_layers = 2L,
      facet = "FacetWrap",
      checks = list(
        list(layer = 1L, keys = fill_keys, columns = c(ymin = "conf.low", ymax = "conf.high")),
        list(layer = 2L, keys = line_keys, columns = c(y = "estimate"))
      )
    ),
    list(
      label = "avg_sops bars without stored draws",
      data = without_draws,
      args = list(geom = "bar", facet_var = "tx"),
      warning = "No stored draws found",
      n_layers = 1L,
      facet = "FacetWrap",
      checks = list(
        list(layer = 1L, keys = fill_keys, columns = c(height = "estimate"))
      )
    ),
    list(
      label = "plain data frame lines in a facet grid",
      data = plain,
      args = list(geom = "line", facet_var = c("tx", "source")),
      warning = NULL,
      n_layers = 2L,
      facet = "FacetGrid",
      checks = list(
        list(
          layer = 2L,
          keys = c(line_keys, source = "source"),
          columns = c(y = "estimate")
        )
      )
    )
  )

  for (case in cases) {
    label <- case$label
    if (is.null(case$warning)) {
      plot <- do.call(plot_sops, c(list(case$data), case$args))
    } else {
      expect_warning(
        plot <- do.call(plot_sops, c(list(case$data), case$args)),
        case$warning,
        fixed = TRUE
      )
    }

    expect_length(plot$layers, case$n_layers)
    expect_identical(class(plot$facet)[1], case$facet, info = label)
    expect_identical(
      plot$labels$title,
      "State Occupancy Probabilities Over Time",
      info = label
    )
    for (check in case$checks) {
      expect_plotted_values(
        plot_layer_values(plot, check$layer, check$keys),
        as.data.frame(case$data),
        check$keys,
        check$columns,
        info = paste(label, "layer", check$layer)
      )
    }
  }
})

test_that("plot_sops respects stored y_levels for model-derived state order", {
  state_levels <- as.character(1:12)
  data <- expand.grid(
    time = 1,
    state = sort(state_levels),
    KEEP.OUT.ATTRS = FALSE
  )
  data$estimate <- 1 / length(state_levels)
  class(data) <- c("markov_avg_sops", class(data))
  attr(data, "y_levels") <- state_levels

  plot <- plot_sops(data, geom = "bar", facet_var = NULL, n_draws = 0)
  built <- ggplot2::ggplot_build(plot)
  fill_scale <- built$plot$scales$get_scales("fill")

  expect_equal(fill_scale$get_limits(), state_levels)
})

test_that("plot_sops overlays a deterministic subset of avg_sops draws", {
  data <- make_avg_sops_plot_data()

  plot <- plot_sops(data, geom = "bar", facet_var = "tx", n_draws = 3)
  built <- ggplot2::ggplot_build(plot)
  draw_layer_ids <- unname(vapply(
    plot$layers[-1],
    function(layer) {
      unique(layer$data$draw_id)
    },
    numeric(1)
  ))
  layer_max_y <- vapply(
    built$data,
    function(layer_data) {
      max(layer_data$y)
    },
    numeric(1)
  )
  draw_layer_nrows <- unname(vapply(
    plot$layers[-1],
    function(layer) nrow(layer$data),
    integer(1)
  ))

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 4)
  expect_equal(draw_layer_ids, c(1, 4, 7))
  expect_equal(draw_layer_nrows, rep(12L, 3))
  expect_equal(layer_max_y, rep(1, 4))
})

test_that("plot_sops rejects ambiguous summary groups", {
  data <- make_avg_sops_plot_data(with_draws = FALSE)
  class(data) <- "data.frame"
  data <- rbind(
    transform(data, scenario = "a"),
    transform(data, scenario = "b")
  )

  line_msg <- tryCatch(
    plot_sops(data, geom = "line", facet_var = "tx"),
    error = conditionMessage
  )
  bar_msg <- tryCatch(
    plot_sops(data, geom = "bar", facet_var = "tx"),
    error = conditionMessage
  )

  expect_match(line_msg, "`facet_var` or `linetype_var`", fixed = TRUE)
  expect_match(bar_msg, "`facet_var`, aggregate", fixed = TRUE)
  expect_no_match(bar_msg, "`linetype_var`", fixed = TRUE)

  plot <- plot_sops(data, geom = "line", facet_var = c("tx", "scenario"))
  expect_s3_class(plot, "ggplot")
})

make_comparison_plot_data <- function(
  comparison = "difference",
  with_time = TRUE
) {
  if (isTRUE(with_time)) {
    data <- expand.grid(
      time = 1:2,
      state_set = c("1", "2"),
      contrast = c("1 - 0", "2 - 0"),
      KEEP.OUT.ATTRS = FALSE
    )
    data$estimand <- "sop"
  } else {
    data <- data.frame(
      state_set = c("1", "2"),
      contrast = c("1 / 0", "1 / 0"),
      estimand = "time_in_state"
    )
  }
  data$variable <- "tx"
  data$reference_level <- 0
  data$comparison_level <- ifelse(data$contrast %in% c("2 - 0", "2 / 0"), 2, 1)
  data$comparison <- comparison
  data$estimate <- seq_len(nrow(data)) / 10
  if (identical(comparison, "ratio")) {
    data$estimate <- data$estimate + 0.8
  }
  data$conf.low <- data$estimate - 0.05
  data$conf.high <- data$estimate + 0.05
  class(data) <- c("markov_avg_comparisons", class(data))
  data
}

test_that("plot_comparisons draws estimates on the difference and ratio scales", {
  custom_colour <- make_comparison_plot_data(
    comparison = "ratio",
    with_time = FALSE
  )
  custom_colour$initial_y <- c("1", "2")

  cases <- list(
    list(
      label = "differences over time",
      data = make_comparison_plot_data(),
      args = list(),
      geom = "GeomLine",
      reference = 0,
      keys = c(x = "time", colour = "state_set", linetype = "contrast"),
      labels = list(
        y = "Difference",
        colour = "State set",
        linetype = "Contrast"
      ),
      colour_var = "state_set"
    ),
    list(
      label = "ratios by state set",
      data = make_comparison_plot_data(comparison = "ratio", with_time = FALSE),
      args = list(),
      geom = "GeomPoint",
      reference = 1,
      keys = c(x = "state_set"),
      labels = list(x = "State set", y = "Ratio"),
      colour_var = NULL
    ),
    list(
      label = "ratios coloured by a custom variable",
      data = custom_colour,
      args = list(color_var = "initial_y"),
      geom = "GeomPoint",
      reference = 1,
      keys = c(x = "state_set", colour = "initial_y"),
      labels = list(colour = "initial_y", fill = "initial_y"),
      colour_var = "initial_y"
    )
  )

  for (case in cases) {
    label <- case$label
    plot <- do.call(plot_comparisons, c(list(case$data), case$args))
    estimates <- plot_layer_values(plot, 3L, case$keys)

    expect_length(plot$layers, 3L)
    expect_true(inherits(plot$layers[[3]]$geom, case$geom), info = label)
    expect_equal(
      unique(plot_layer_values(plot, 1L, character())$yintercept),
      case$reference,
      info = label
    )
    for (name in names(case$labels)) {
      expect_identical(
        plot$labels[[name]],
        case$labels[[name]],
        info = paste(label, name)
      )
    }
    if (!is.null(case$colour_var)) {
      expect_length(
        unique(estimates$colour),
        length(unique(case$data[[case$colour_var]]))
      )
    }
    expect_plotted_values(
      estimates,
      as.data.frame(case$data),
      case$keys,
      c(y = "estimate"),
      info = label
    )
  }
})

test_that("plot_comparisons draws one profile and rejects several", {
  data <- make_comparison_plot_data()
  data$rowid <- 1L
  class(data) <- c("markov_comparisons", "data.frame")

  plot <- plot_comparisons(data)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$layers[[3]]$geom, "GeomLine")
  expect_equal(plot$labels$title, "SOP Comparisons for One Profile")

  several <- rbind(data, transform(data, rowid = 2L))
  expect_error(plot_comparisons(several), "2 profiles")
  expect_error(
    plot_comparisons(as.data.frame(several)),
    "subset(x, rowid == 1)",
    fixed = TRUE
  )
  chosen <- plot_comparisons(several[several$rowid == 2, ])
  expect_equal(nrow(chosen$data), nrow(data))
})

test_that("plot_comparisons orders numeric state sets naturally", {
  data <- data.frame(
    state_set = c("1", "10", "11", "2", "3"),
    contrast = "1 / 0",
    estimand = "time_in_state",
    variable = "tx",
    reference_level = 0,
    comparison_level = 1,
    comparison = "ratio",
    estimate = c(1, 1.1, 1.2, 0.9, 0.8)
  )
  data$conf.low <- data$estimate - 0.05
  data$conf.high <- data$estimate + 0.05
  class(data) <- c("markov_avg_comparisons", class(data))
  attr(data, "y_levels") <- as.character(1:11)

  plot <- plot_comparisons(data)
  built <- ggplot2::ggplot_build(plot)

  expect_equal(
    built$layout$panel_params[[1]]$x$get_labels(),
    c("1", "2", "3", "10", "11")
  )
})

test_that("plot_sops and plot_comparisons reject invalid inputs", {
  sops <- make_avg_sops_plot_data()
  mixed_scales <- make_comparison_plot_data()
  mixed_scales$comparison[1] <- "ratio"

  cases <- list(
    list(
      label = "infinite n_draws",
      call = function() plot_sops(sops, geom = "bar", n_draws = Inf),
      error = "`n_draws`"
    ),
    list(
      label = "missing draw_alpha",
      call = function() plot_sops(sops, geom = "bar", draw_alpha = NA),
      error = "`draw_alpha`"
    ),
    list(
      label = "missing show_uncertainty",
      call = function() plot_sops(sops, geom = "bar", show_uncertainty = NA),
      error = "`show_uncertainty`"
    ),
    list(
      label = "comparison data without required columns",
      call = function() plot_comparisons(data.frame(estimate = 1)),
      error = "missing required comparison column"
    ),
    list(
      label = "comparison data on two scales",
      call = function() plot_comparisons(mixed_scales),
      error = "requires one comparison scale"
    )
  )

  for (case in cases) {
    expect_error(case$call(), case$error, fixed = TRUE, info = case$label)
  }
})
