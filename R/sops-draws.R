# SOP draw summaries and extraction.

# The summary of one estimate's draws `x`: its standard error (their standard
# deviation) and its interval. A percentile interval is given as it is; a Wald
# interval as the two distances from its centre, which the caller adds.
draw_summary_function <- function(conf_level, conf_type) {
  alpha <- 1 - conf_level
  function(x) {
    # Missing and undefined (NaN) draws are left out; infinite draws, such as
    # a ratio with a reference probability of 0, are kept, so a percentile
    # limit can be Inf. The standard deviation of draws that include an
    # infinite value is not finite, so the standard error (and with it a Wald
    # interval) is NA rather than the NaN that `sd()` returns.
    se <- if (any(is.infinite(x))) NA_real_ else stats::sd(x, na.rm = TRUE)
    if (conf_type == "perc") {
      return(c(
        conf.low = unname(stats::quantile(x, alpha / 2, na.rm = TRUE)),
        conf.high = unname(stats::quantile(x, 1 - alpha / 2, na.rm = TRUE)),
        std.error = se
      ))
    }

    critical <- abs(stats::qnorm(alpha / 2))
    c(
      conf.low = -critical * se,
      conf.high = critical * se,
      std.error = se
    )
  }
}

#' Compute Confidence Intervals from Draws
#'
#' Computes confidence intervals from simulation or bootstrap draws.
#'
#' @param draws_df Data frame with draws (must have `estimate` and `draw_id` columns).
#' @param group_cols Character vector of grouping columns.
#' @param conf_level Confidence level (default 0.95).
#' @param conf_type Type of CI: "perc" (percentile) or "wald".
#' @param keep_na_groups Logical. If `FALSE` (the default), draws with a
#'   missing value in `group_cols` are omitted. If `TRUE`, a missing value
#'   forms its own group, as needed for profile columns the model does not use.
#'
#' @return Data frame with summary statistics.
#'
#' @keywords internal
compute_ci_from_draws <- function(
  draws_df,
  group_cols,
  conf_level = 0.95,
  conf_type = "perc",
  point_estimates = NULL,
  keep_na_groups = FALSE
) {
  .datatable.aware <- TRUE
  draws_df <- as.data.frame(draws_df)
  conf_level <- validate_conf_level(conf_level)
  if (
    !is.character(conf_type) ||
      length(conf_type) != 1L ||
      is.na(conf_type) ||
      !conf_type %in% c("perc", "wald")
  ) {
    mostr_abort("conf_type must be 'perc' or 'wald'")
  }
  summarize_values <- draw_summary_function(conf_level, conf_type)
  warn_undefined_draws(draws_df, group_cols, keep_na_groups)

  if (length(group_cols) == 0L) {
    stats <- summarize_values(draws_df$estimate)
    out <- data.frame(
      conf.low = stats[["conf.low"]],
      conf.high = stats[["conf.high"]],
      std.error = stats[["std.error"]]
    )
    if (conf_type == "wald") {
      center <- if (is.null(point_estimates)) {
        mean(draws_df$estimate)
      } else {
        point_estimates$estimate[1]
      }
      out$conf.low <- center + out$conf.low
      out$conf.high <- center + out$conf.high
    }
    return(out)
  }

  work <- data.table::as.data.table(
    draws_df[, c(group_cols, "estimate"), drop = FALSE]
  )
  work <- if (isTRUE(keep_na_groups)) {
    work[!is.na(work$estimate)]
  } else {
    work[stats::complete.cases(work)]
  }
  if (nrow(work) == 0L) {
    mostr_abort("no rows to aggregate")
  }
  summary_stats <- as.data.frame(work[,
    as.list(summarize_values(estimate)),
    by = group_cols
  ])
  rows <- do.call(order, lapply(summary_stats[rev(group_cols)], as.factor))
  summary_stats <- summary_stats[rows, , drop = FALSE]
  rownames(summary_stats) <- NULL
  if (conf_type == "wald") {
    if (is.null(point_estimates)) {
      centers <- aggregate_value(
        draws_df,
        "estimate",
        group_cols,
        mean,
        keep_na_groups = keep_na_groups
      )
    } else {
      centers <- point_estimates[, c(group_cols, "estimate"), drop = FALSE]
    }
    names(centers)[names(centers) == "estimate"] <- ".center"
    summary_stats <- merge(
      summary_stats,
      centers,
      by = group_cols,
      all.x = TRUE,
      sort = FALSE
    )
    summary_stats$conf.low <- summary_stats$.center + summary_stats$conf.low
    summary_stats$conf.high <- summary_stats$.center + summary_stats$conf.high
    summary_stats$.center <- NULL
  }

  summary_stats
}

# Undefined draws (NaN, as 0 / 0 gives for a ratio of two probabilities of 0)
# can't be ordered, so the interval leaves them out. That changes which
# quantiles the limits are, so the user is told how many were left out. With
# several estimates, the message gives the largest count, from the estimate
# that lost the most draws. Missing (NA) draws are left out without a warning,
# as before.
warn_undefined_draws <- function(draws_df, group_cols, keep_na_groups) {
  undefined <- is.nan(draws_df$estimate)
  if (!any(undefined)) {
    return(invisible(NULL))
  }
  key <- split_key(draws_df, group_cols, keep_na = keep_na_groups)
  n_undefined <- tapply(undefined, key, sum)
  n_total <- tapply(undefined, key, length)
  affected <- which(!is.na(n_undefined) & n_undefined > 0)
  if (length(affected) == 0L) {
    return(invisible(NULL))
  }
  worst <- affected[which.max(n_undefined[affected])]
  counts <- paste(n_undefined[[worst]], "of", n_total[[worst]], "draws")
  if (length(affected) == 1L) {
    warning(
      counts,
      " of an estimate were undefined (0 / 0) and were left out of its ",
      "interval.",
      call. = FALSE
    )
  } else {
    warning(
      "Draws of ",
      length(affected),
      " estimates were undefined (0 / 0) and were left out of their ",
      "intervals: up to ",
      counts,
      " of an estimate.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Extract Individual Draws from Inference Objects
#'
#' Extracts the individual draws (bootstrap samples, simulated values from MVN, etc.)
#' from an object returned by `inferences()` with `return_draws = TRUE`.
#' This function joins the draws back to the original point estimate object,
#' preserving all covariates, grouping variables, and summary statistics.
#'
#' @param x A SOP, average-time, or average-comparison object returned with
#'   stored draws.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item draw_id: Simulation or bootstrap iteration number
#'     \item time, state: Time point and state for SOP results
#'     \item state_set: State set for average-time or comparison results
#'     \item draw: Draw-specific estimate of the requested quantity
#'     \item estimate: The original point estimate from the model
#'     \item conf.low, conf.high, std.error: Summary statistics from the point estimate object
#'     \item fwb_weight: Optional draw-specific patient weight
#'       for eligible stored-data `sops()` or [comparisons()] draws
#'     \item Additional columns from the original object (covariates, etc.)
#'   }
#'   Each row represents one draw for a specific result row.
#'
#' @details
#' This function retrieves the draws from objects created by `inferences()`.
#' It performs a join between the draws and the point estimate object to ensure
#' that all metadata is preserved. To avoid conflict, the `estimate` column from
#' the draws is renamed to `draw`.
#'
#' Draw-specific patient weights are included only when they are meaningful for
#' manual downstream averaging: `sops()` and [comparisons()] output
#' from FWB on the stored empirical prediction cohort. For
#' [comparisons()] draws, the weighted mean of a draw's comparisons,
#' \eqn{\sum_i w_i x_i / \sum_i w_i} over patients \eqn{i} with weights
#' \eqn{w_i} and draw values \eqn{x_i}, reproduces the corresponding
#' [avg_comparisons()] draw. They are not included for
#' `avg_sops()` because its draws are already averaged, and
#' they are not included when the original SOP object used user-supplied
#' prediction profiles. The column name `fwb_weight` is reserved for these draw
#' weights when they are attached.
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
#' # Average SOPs under control and under treatment, with 20 draws of the
#' # coefficients; use 1000 or more in an analysis
#' sops_by_arm <- inferences(
#'   avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:14),
#'   method = "mvn",
#'   n_draws = 20,
#'   seed = 1
#' )
#'
#' # One row per draw and result row; the column draw holds the draw's value
#' draws <- get_draws(sops_by_arm)
#' head(draws)
#'
#' # The probability of death on day 14 in each draw, one column per arm
#' dead <- draws[draws$time == 14 & draws$state == "8", ]
#' dead <- tapply(dead$draw, list(dead$draw_id, dead$tx), sum)
#'
#' # A 95% interval for the risk ratio of death on day 14, from the same draws
#' stats::quantile(dead[, "1"] / dead[, "0"], c(0.025, 0.5, 0.975))
#'
#' @seealso [inferences()], [avg_sops()]
#'
#' @export
get_draws <- function(x) {
  if (
    !inherits(
      x,
      c(
        "markov_avg_sops",
        "markov_sops",
        "markov_avg_state_time",
        "markov_state_time",
        "markov_avg_comparisons",
        "markov_comparisons"
      )
    )
  ) {
    mostr_abort(
      "get_draws() requires an object from inferences(). ",
      "Got: ",
      paste(class(x), collapse = ", ")
    )
  }

  # 1. Extract draws attribute
  draws <- attr(x, "draws")

  if (is.null(draws)) {
    method <- attr(x, "method")
    msg <- "No draws found. Run inferences() with return_draws = TRUE."
    if (!is.null(method)) {
      msg <- paste0(msg, " (Method used: '", method, "')")
    }
    mostr_abort(msg)
  }

  # 2. Prepare metadata from the original object
  meta <- as.data.frame(x)

  # Rename 'estimate' in draws to 'draw' to avoid conflict with point estimates
  if ("estimate" %in% names(draws)) {
    names(draws)[names(draws) == "estimate"] <- "draw"
  }

  # 3. Identify join keys
  # Keys are columns present in both draws and meta (excluding estimate/draw)
  common_cols <- intersect(names(draws), names(meta))
  keys <- setdiff(common_cols, c("draw", "estimate"))

  if (length(keys) == 0) {
    return(draws)
  }

  # 4. Join metadata back to draws
  draws <- merge(draws, meta, by = keys, all.x = TRUE, sort = FALSE)

  return(draws)
}
