make_comparisons_case <- function() {
  skip_if_not_installed("VGAM")
  withr::local_seed(4417)
  data <- make_markov_trajectories(
    n_patients = 60,
    follow_up_time = 5,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = 1:2,
    absorbing_state = 3,
    treatment_effect = 0.4,
    seed = 4417
  ) |>
    prepare_markov_data(absorbing_state = 3)
  data$marker <- data$id %% 2
  model <- vglm_markov(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id"
  )
  list(
    model = model,
    data = data,
    profiles = data[data$time == 1, , drop = FALSE],
    time_map = c("1" = 1, "2" = 2, "3" = 4, "4" = 7)
  )
}

# SOPs for every profile with `var` set to `value`, in `sops()` row order:
# profile fastest, then visit, then state.
sops_with_value <- function(case, var, value, times = 1:4, ...) {
  profiles <- case$profiles
  profiles[[var]] <- value
  sops(case$model, newdata = profiles, times = times, absorb = 3, ...)
}

test_that("comparisons() returns paired SOP differences for each profile", {
  case <- make_comparisons_case()
  out <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1)),
    times = 1:4,
    absorb = 3
  )
  lo <- sops_with_value(case, "tx", 0)
  hi <- sops_with_value(case, "tx", 1)
  n <- nrow(case$profiles)

  expect_s3_class(out, "markov_comparisons")
  expect_false(inherits(out, "markov_sops"))
  expect_false(inherits(out, "markov_avg_comparisons"))
  expect_identical(nrow(out), n * 4L * 3L)
  for (state in as.character(1:3)) {
    rows <- out[out$state_set == state, , drop = FALSE]
    keep <- as.character(hi$state) == state
    expect_equal(rows$rowid, hi$rowid[keep])
    expect_equal(rows$time, hi$time[keep])
    expect_equal(rows$estimate, hi$estimate[keep] - lo$estimate[keep])
  }
  expect_true(all(out$estimand == "sop"))
  expect_true(all(out$term == "tx"))
  expect_true(all(out$contrast == "1 - 0"))
  expect_true(all(out$comparison == "difference"))
  # Profile columns hold each patient's observed values.
  expect_equal(out$tx, case$profiles$tx[out$rowid])
  expect_equal(out$id, case$profiles$id[out$rowid])
  expect_equal(out$marker, case$profiles$marker[out$rowid])
  expect_identical(names(out)[ncol(out)], "estimate")
})

test_that("comparisons() lumps state sets and forms ratios for every level", {
  case <- make_comparisons_case()
  out <- comparisons(
    case$model,
    newdata = case$profiles,
    variables = list(tx = c(0, 1, 2)),
    state_sets = list(alive = c("1", "2")),
    comparison = "ratio",
    times = 1:4,
    absorb = 3
  )
  alive <- function(x) {
    keep <- as.character(x$state) %in% c("1", "2")
    stats::aggregate(estimate ~ rowid + time, x[keep, ], sum)
  }
  lo <- alive(sops_with_value(case, "tx", 0))
  expect_identical(unique(out$state_set), "alive")
  expect_identical(unique(out$contrast), c("1 / 0", "2 / 0"))
  for (level in c(1, 2)) {
    hi <- alive(sops_with_value(case, "tx", level))
    rows <- out[out$comparison_level == level, , drop = FALSE]
    expected <- merge(hi, lo, by = c("rowid", "time"))
    expected <- expected[match(
      paste(rows$rowid, rows$time),
      paste(expected$rowid, expected$time)
    ), ]
    expect_equal(rows$estimate, expected$estimate.x / expected$estimate.y)
  }
})

test_that("averaging comparisons() reproduces avg_comparisons()", {
  case <- make_comparisons_case()
  variables <- list(tx = c(0, 1, 2))
  settings <- list(
    list(estimand = "sop", state_sets = list(alive = c("1", "2"), dead = "3")),
    list(estimand = "time_in_state", state_sets = NULL),
    list(
      estimand = "time_in_state",
      state_sets = "1",
      time_map = case$time_map
    ),
    list(estimand = "time_benefit", state_sets = NULL),
    list(estimand = "time_benefit", state_sets = NULL, time_map = case$time_map)
  )
  for (setting in settings) {
    args <- c(
      list(
        case$model,
        variables = variables,
        times = 1:4,
        absorb = 3
      ),
      setting
    )
    unit <- do.call(comparisons, args)
    average <- do.call(avg_comparisons, args)
    keys <- intersect(c("time", "state_set", "contrast"), names(average))
    expected <- stats::aggregate(
      unit$estimate,
      unit[, keys, drop = FALSE],
      mean
    )
    index <- match(
      do.call(paste, average[, keys, drop = FALSE]),
      do.call(paste, expected[, keys, drop = FALSE])
    )
    expect_false(anyNA(index))
    expect_equal(expected$x[index], average$estimate, tolerance = 1e-10)
    expect_identical(
      nrow(unit),
      nrow(average) * length(unique(unit$rowid))
    )
  }
})

test_that("comparisons() starts real-time totals from a set starting state", {
  case <- make_comparisons_case()
  levels <- levels(case$profiles$yprev)
  starts <- factor(c("1", "2"), levels = levels)
  args <- list(
    case$model,
    variables = list(yprev = starts),
    estimand = "time_in_state",
    state_sets = "1",
    times = 1:4,
    absorb = 3,
    time_map = case$time_map,
    target_times = seq(0, 7, by = 0.5)
  )
  unit <- do.call(comparisons, args)
  average <- do.call(avg_comparisons, args)
  expect_equal(mean(unit$estimate), average$estimate, tolerance = 1e-10)
  stored <- attr(unit, "comparison_baseline_data")
  expect_equal(unit$yprev, stored$yprev[unit$rowid])
})

test_that("comparisons(by = ) hands off to avg_comparisons()", {
  case <- make_comparisons_case()
  args <- list(
    case$model,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    comparison = "ratio",
    times = 1:4,
    absorb = 3
  )
  expect_equal(
    do.call(comparisons, c(args, list(by = TRUE))),
    do.call(avg_comparisons, args)
  )
  expect_equal(
    do.call(comparisons, c(args, list(by = "marker"))),
    do.call(avg_comparisons, c(args, list(by = "marker")))
  )
  expect_s3_class(
    do.call(comparisons, c(args, list(by = NULL))),
    "markov_comparisons"
  )
  expect_error(
    do.call(comparisons, c(args, list(by = 1))),
    "`by` must be"
  )
  expect_error(
    do.call(comparisons, c(args, list(by = NA_character_))),
    "`by` must be"
  )
})

test_that("comparisons() validates arguments before predicting", {
  case <- make_comparisons_case()
  local_mocked_bindings(
    sops = function(...) {
      stop("Predictions should not run before argument validation.")
    }
  )
  clash <- case$profiles
  clash$contrast <- "a"
  clash$estimate <- 1
  expect_error(
    comparisons(
      case$model,
      newdata = clash,
      variables = list(tx = c(0, 1)),
      times = 1:4
    ),
    "`contrast`, `estimate`"
  )
  expect_error(
    comparisons(case$model, variables = list(tx = c(0, 1))),
    "`times` must be supplied to `comparisons\\(\\)`"
  )
  expect_error(
    comparisons(case$model, variables = list(tx = 1), times = 1:4),
    "at least two values"
  )
  expect_error(
    comparisons(
      case$model,
      variables = list(tx = c(0, 1), marker = c(0, 1)),
      times = 1:4
    ),
    "exactly one named variable"
  )
  expect_error(
    comparisons(
      case$model,
      variables = list(tx = c(0, 1)),
      times = c(1, 1, 2)
    ),
    "duplicate visits"
  )
  expect_error(
    comparisons(
      case$model,
      variables = list(tx = c(0, 1)),
      estimand = "time_benefit",
      comparison = "ratio",
      times = 1:4
    ),
    "only supports `comparison = \"difference\"`"
  )
})

test_that("posterior comparisons summarize paired draws for each profile", {
  coefficients <- rbind(
    c(0.5, -0.5, 0.4, 1.2),
    c(2, 1, -0.8, -1.1),
    c(-0.4, -2, 1.5, 0.2)
  )
  colnames(coefficients) <- c("y>=2", "y>=3", "tx", "time")
  model <- structure(
    list(
      draws = coefficients,
      non.slopes = 2L,
      pppo = 0L,
      tauInfo = data.frame(name = character()),
      ylevels = 1:3,
      yname = "y",
      clusterInfo = list(name = "id", cluster = c("a", "b"))
    ),
    class = c("blrm", "orm"),
    markov_fit_wrapper = "blrm_markov"
  )
  profiles <- data.frame(
    id = c("a", "b"),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1
  )
  local_mocked_bindings(
    blrm_design_matrix = function(model, newdata, second = FALSE) {
      cbind(tx = newdata$tx, time = newdata$time)
    }
  )
  posterior_sops <- function(value) {
    newdata <- profiles
    newdata$tx <- value
    attr(
      sops(
        model,
        newdata = newdata,
        times = 1:3,
        n_draws = 3,
        seed = 9325,
        return_draws = TRUE
      ),
      "draws"
    )
  }
  out <- comparisons(
    model,
    newdata = profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:3,
    n_draws = 3,
    seed = 9325,
    posterior_summary = "median",
    return_draws = TRUE
  )
  lo <- posterior_sops(0)
  hi <- posterior_sops(1)
  keep <- as.character(lo$state) == "1"
  paired <- data.frame(
    rowid = lo$rowid[keep],
    time = lo$time[keep],
    draw_id = lo$draw_id[keep],
    estimate = hi$estimate[keep] - lo$estimate[keep]
  )
  for (i in seq_len(nrow(out))) {
    values <- paired$estimate[
      paired$rowid == out$rowid[i] & paired$time == out$time[i]
    ]
    expect_equal(out$estimate[i], stats::median(values))
    expect_equal(out$std.error[i], stats::sd(values))
    expect_equal(out$conf.low[i], unname(stats::quantile(values, 0.025)))
    expect_equal(out$conf.high[i], unname(stats::quantile(values, 0.975)))
  }
  expect_equal(out$tx, profiles$tx[out$rowid])
  expect_identical(attr(out, "method"), "posterior")
  draws <- get_draws(out)
  expect_identical(nrow(draws), 3L * nrow(out))
  expect_true(all(c("draw", "draw_id", "contrast", "tx") %in% names(draws)))
  expect_identical(inferences(out, method = "mvn", n_draws = 2), out)

  compact <- comparisons(
    model,
    newdata = profiles,
    variables = list(tx = c(0, 1)),
    state_sets = "1",
    times = 1:3,
    n_draws = 3,
    seed = 9325
  )
  expect_null(attr(compact, "draws"))
  means <- stats::aggregate(estimate ~ rowid + time, paired, mean)
  index <- match(
    paste(compact$rowid, compact$time),
    paste(means$rowid, means$time)
  )
  expect_equal(compact$estimate, means$estimate[index])
})
