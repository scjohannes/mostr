make_mock_avg_sops <- function(by = NULL, tx_levels = c(0, 1, 2)) {
  if (is.null(by)) {
    out <- expand.grid(
      time = 1:2,
      state = c("1", "2"),
      tx = tx_levels,
      KEEP.OUT.ATTRS = FALSE
    )
    values <- c(
      "0" = 0.2,
      "1" = 0.5,
      "2" = 0.1
    )
    out$estimate <- values[as.character(out$tx)]
    out$estimate[out$time == 2 & out$tx == 0] <- 0.4
    out$estimate[out$time == 2 & out$tx == 1] <- 0.7
    out$estimate[out$time == 2 & out$tx == 2] <- 0.2
    out$estimate[out$state == "2"] <- 1 - out$estimate[out$state == "2"]
  } else {
    out <- expand.grid(
      time = 1,
      state = c("1", "2"),
      tx = c(0, 1),
      sex = c("f", "m"),
      KEEP.OUT.ATTRS = FALSE
    )
    out$estimate <- c(0.2, 0.8, 0.6, 0.4, 0.3, 0.7, 0.9, 0.1)
  }
  class(out) <- c("markov_avg_sops", class(out))
  attr(out, "avg_args") <- list(
    variables = list(tx = tx_levels),
    by = by,
    times = sort(unique(out$time)),
    id_var = "id"
  )
  attr(out, "y_levels") <- c("1", "2")
  attr(out, "call_args") <- list(times = sort(unique(out$time)))
  attr(out, "time_var") <- "time"
  attr(out, "p_var") <- "yprev"
  attr(out, "newdata_orig") <- data.frame(
    id = 1:2,
    yprev = c("1", "2")
  )
  attr(out, "newdata_supplied") <- TRUE
  out
}

test_that("avg_comparisons() forms differences, ratios and by strata from averaged SOPs", {
  model <- structure(list(), class = "mock_model", markov_id_var = "id")
  # `order` sorts the result before it is compared with `expected`, column by
  # column.
  cases <- list(
    list(
      label = "SOP differences for two comparison levels",
      avg = make_mock_avg_sops(),
      args = list(
        variables = list(tx = c(0, 1, 2)),
        estimand = "sop",
        state_sets = "1",
        times = 1:2,
        comparison = "difference"
      ),
      order = c("comparison_level", "time"),
      expected = list(
        estimate = c(0.3, 0.3, -0.1, -0.2),
        reference_level = c(0, 0, 0, 0),
        comparison_level = c(1, 1, 2, 2)
      )
    ),
    list(
      label = "time-in-state ratio",
      avg = make_mock_avg_sops(tx_levels = c(0, 1)),
      args = list(
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        comparison = "ratio"
      ),
      order = character(),
      expected = list(
        estimate = (0.5 + 0.7) / (0.2 + 0.4),
        state_set = "1",
        contrast = "1 / 0"
      )
    ),
    list(
      label = "SOP differences within by strata",
      avg = make_mock_avg_sops(by = "sex", tx_levels = c(0, 1)),
      args = list(
        variables = list(tx = c(0, 1)),
        estimand = "sop",
        state_sets = "1",
        times = 1,
        by = "sex"
      ),
      order = "sex",
      expected = list(sex = c("f", "m"), estimate = c(0.4, 0.6))
    )
  )

  for (case in cases) {
    out <- with_mocked_bindings(
      avg_sops = function(...) case$avg,
      do.call(avg_comparisons, c(list(model), case$args))
    )
    expect_true(inherits(out, "markov_avg_comparisons"), info = case$label)
    if (length(case$order) > 0L) {
      out <- out[do.call(order, unname(as.list(out[case$order]))), ]
    }
    for (column in names(case$expected)) {
      values <- out[[column]]
      if (is.factor(values)) {
        values <- as.character(values)
      }
      expect_equal(
        values,
        case$expected[[column]],
        info = paste(case$label, column, sep = ": ")
      )
    }
  }
})

test_that("avg_comparisons() uses real-time AUC for time-in-state", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model", markov_id_var = "id"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        time_map = c("1" = 1, "2" = 3),
        baseline_time = NULL
      )
    }
  )

  expect_equal(out$estimate, 0.6)
})

test_that("real-time comparisons use baseline support without integrating it", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      explicit <- avg_comparisons(
        structure(list(), class = "mock_model", markov_id_var = "id"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        time_map = c("1" = 7, "2" = 14),
        baseline_time = 0,
        target_times = 1:14
      )
      disabled <- avg_comparisons(
        structure(list(), class = "mock_model", markov_id_var = "id"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        time_map = c("1" = 7, "2" = 14),
        baseline_time = NULL
      )
    }
  )

  difference <- stats::approx(
    x = c(0, 7, 14),
    y = c(0, 0.3, 0.3),
    xout = 1:14
  )$y
  expected <- sum(
    utils::head(difference, -1L) + utils::tail(difference, -1L)
  ) /
    2
  expect_equal(explicit$estimate, expected)
  expect_equal(disabled$estimate, 2.1)
})

test_that("avg_comparisons() forwards model-structure arguments to the SOP call", {
  model <- structure(list(), class = "mock_model", markov_id_var = "id")
  structure_args <- list(
    time_var = "visit",
    p_var = "prev",
    p2_var = "preprev",
    gap_var = "gap"
  )

  # The time-benefit estimand pairs the profiles itself, so it calls `sops()`;
  # the other estimands call `avg_sops()`. `forwarded` lists the arguments the
  # mocked call must receive under exactly these names.
  cases <- list(
    list(
      label = "time in state, through avg_sops()",
      time_covariates = data.frame(visit = 1:2, spline = c(0, 1)),
      run = function(captured, time_covariates) {
        avg <- make_mock_avg_sops(tx_levels = c(0, 1))
        with_mocked_bindings(
          avg_sops = function(...) {
            captured$args <- list(...)
            avg
          },
          do.call(
            avg_comparisons,
            c(
              list(
                model,
                variables = list(tx = c(0, 1)),
                estimand = "time_in_state",
                state_sets = "1",
                times = 1:2,
                y_levels = c("1", "2"),
                absorb = "2",
                time_covariates = time_covariates
              ),
              structure_args
            )
          )
        )
      },
      forwarded = c(structure_args, list(y_levels = c("1", "2"), absorb = "2")),
      expected = list()
    ),
    list(
      label = "time benefit from paired profiles, through sops()",
      time_covariates = data.frame(visit = 1, spline = 0),
      run = function(captured, time_covariates) {
        newdata <- data.frame(
          id = 1:2,
          tx = c(0, 1),
          yprev = c(2, 2),
          prev = c(2, 2),
          preprev = c(1, 1)
        )
        with_mocked_bindings(
          validate_markov_model = function(object) NULL,
          sops = function(model, newdata, times, ...) {
            captured$args <- list(...)
            out <- expand.grid(
              rowid = seq_len(nrow(newdata)),
              time = 1,
              state = c("1", "2"),
              KEEP.OUT.ATTRS = FALSE
            )
            out$tx <- rep(newdata$tx, times = 2)
            out$yprev <- rep(newdata$yprev, times = 2)
            out$prev <- rep(newdata$prev, times = 2)
            out$estimate <- c(0, 0, 1, 0, 1, 1, 0, 1)
            class(out) <- c("markov_sops", class(out))
            attr(out, "call_args") <- list(times = 1)
            attr(out, "y_levels") <- c("1", "2")
            attr(out, "time_var") <- captured$args[["time_var"]]
            attr(out, "p_var") <- captured$args[["p_var"]]
            attr(out, "p2_var") <- captured$args[["p2_var"]]
            out
          },
          do.call(
            avg_comparisons,
            c(
              list(
                model,
                newdata = newdata,
                variables = list(tx = c(0, 1)),
                estimand = "time_benefit",
                times = 1,
                y_levels = c("1", "2"),
                time_covariates = time_covariates
              ),
              structure_args
            )
          )
        )
      },
      forwarded = structure_args,
      expected = list(estimate = 0.5, estimand = "time_benefit")
    )
  )

  for (case in cases) {
    captured <- new.env(parent = emptyenv())
    out <- case$run(captured, case$time_covariates)
    expect_true(inherits(out, "markov_avg_comparisons"), info = case$label)
    for (name in names(case$expected)) {
      expect_equal(out[[name]], case$expected[[name]], info = case$label)
    }
    forwarded <- c(case$forwarded, list(time_covariates = case$time_covariates))
    for (name in names(forwarded)) {
      # `[[` matches names exactly, so a misnamed argument cannot pass by
      # partial matching.
      expect_equal(
        captured$args[[name]],
        forwarded[[name]],
        info = paste(case$label, name, sep = ": ")
      )
    }
  }
})

test_that("time-benefit array reducer keeps scenario blocks on visit and real-time scales", {
  mapped_baseline <- data.frame(id = 1, tx = 0, yprev = 2, rowid = 1)
  mapped_grid <- data.frame(tx = c(0, 1))
  # `sops_array` holds SOPs by prediction row, visit and state; the rows come
  # in one block of `n_each` profiles per scenario.
  cases <- list(
    list(
      label = "visit scale, two profiles per scenario",
      setup = list(
        n_cf = 2,
        n_each = 2,
        variables = list(tx = c(0, 1)),
        y_levels = c("1", "2"),
        times = 1,
        by = NULL,
        baseline_data = data.frame(id = 1:2)
      ),
      sops_array = local({
        a <- array(NA_real_, dim = c(4, 1, 2))
        a[1, 1, ] <- c(0, 1)
        a[2, 1, ] <- c(0, 1)
        a[3, 1, ] <- c(1, 0)
        a[4, 1, ] <- c(0, 1)
        a
      }),
      time_args = list(),
      expected = 0.5
    ),
    list(
      label = "real-time scale through time_map",
      setup = list(
        n_cf = 2,
        n_each = 1,
        variables = list(tx = c(0, 1)),
        y_levels = c("1", "2"),
        times = 1:2,
        by = NULL,
        baseline_data = mapped_baseline,
        newdata_orig = mapped_baseline,
        newdata_pred = mostr:::create_counterfactual_data(
          mapped_baseline,
          mapped_grid,
          list(tx = c(0, 1))
        ),
        grid = mapped_grid,
        id_var = "id",
        time_var = "time",
        p_var = "yprev",
        newdata_supplied = TRUE
      ),
      sops_array = local({
        a <- array(NA_real_, dim = c(2, 2, 2))
        a[1, 1, ] <- c(0, 1)
        a[1, 2, ] <- c(0, 1)
        a[2, 1, ] <- c(0, 1)
        a[2, 2, ] <- c(1, 0)
        a
      }),
      time_args = list(time_map = c("1" = 0, "2" = 4), baseline_time = NULL),
      expected = 2
    )
  )

  for (case in cases) {
    out <- do.call(
      mostr:::reduce_time_benefit_array_for_setup,
      c(
        list(
          sops_array = case$sops_array,
          setup = case$setup,
          comparison = "difference",
          weights = NULL,
          time_unit = NULL
        ),
        case$time_args
      )
    )
    expect_equal(out$estimate, case$expected, info = case$label)
  }
})

test_that("linear ordinal-benefit reducer matches the dense score definition", {
  reference <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.7, 0.2, 0.1)
  )
  comparison <- rbind(
    c(0.4, 0.4, 0.2),
    c(0.1, 0.3, 0.6)
  )
  states <- seq_len(ncol(reference))
  score <- outer(states, states, function(ref, cmp) sign(ref - cmp))
  expected <- rowSums((reference %*% score) * comparison)

  expect_equal(
    mostr:::ordinal_pairwise_benefit(reference, comparison),
    expected,
    tolerance = 1e-15
  )
})

test_that("time-benefit simulation inference uses real-time comparison scale", {
  baseline <- data.frame(id = 1, tx = 0, yprev = 2, rowid = 1)
  grid <- data.frame(tx = c(0, 1))
  newdata_pred <- mostr:::create_counterfactual_data(
    baseline,
    grid,
    list(tx = c(0, 1))
  )
  object <- data.frame(
    estimand = "time_benefit",
    term = "tx",
    reference_level = 0,
    comparison_level = 1,
    contrast = "1 - 0",
    comparison = "difference",
    estimate = 2
  )
  class(object) <- c("markov_avg_comparisons", class(object))
  attr(object, "model") <- structure(
    list(coefficients = c(a = 0, b = 0)),
    class = "mock_model"
  )
  attr(object, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = NULL,
    times = 1:2,
    id_var = "id"
  )
  attr(object, "comparison_args") <- list(
    estimand = "time_benefit",
    comparison = "difference",
    time_map = c("1" = 0, "2" = 4),
    baseline_time = NULL,
    target_times = NULL,
    time_unit = NULL
  )
  attr(object, "newdata_orig") <- baseline
  attr(object, "newdata_pred") <- newdata_pred
  attr(object, "newdata_supplied") <- TRUE
  attr(object, "comparison_baseline_data") <- baseline
  attr(object, "comparison_grid") <- grid
  attr(object, "comparison_n_each") <- 1L
  attr(object, "call_args") <- list(times = 1:2, y_levels = c("1", "2"))
  attr(object, "time_var") <- "time"
  attr(object, "p_var") <- "yprev"
  attr(object, "y_levels") <- c("1", "2")

  sops_array <- array(NA_real_, dim = c(2, 2, 2))
  sops_array[1, 1, ] <- c(0, 1)
  sops_array[1, 2, ] <- c(0, 1)
  sops_array[2, 1, ] <- c(0, 1)
  sops_array[2, 2, ] <- c(1, 0)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2) * 1e-8,
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) sops_array,
    {
      withr::local_seed(1)
      out <- mostr:::inferences_avg_comparisons_time_benefit_simulation(
        object = object,
        engine = "mvn",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = NULL,
        cluster = NULL,
        workers = NULL,
        conf_level = 0.95,
        conf_type = "perc",
        return_draws = TRUE
      )
    }
  )

  expect_equal(out$conf.low, 2)
  expect_equal(out$conf.high, 2)
  expect_equal(attr(out, "draws")$estimate, c(2, 2))
})

test_that("avg_comparisons() validates v1 interface boundaries", {
  avg <- make_mock_avg_sops()

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model", markov_id_var = "id"),
          variables = list(tx = c(0, 1))
        ),
        "`times` must be supplied to `avg_comparisons()`.",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model", markov_id_var = "id"),
          variables = list(tx = c(0, 1)),
          times = 1,
          origin = "none"
        ),
        "`origin` is no longer supported; use `baseline_time` instead.",
        fixed = TRUE
      )
    }
  )
})

test_that("comparison inference reduces paired draw-level SOPs", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))
  draws <- avg[rep(seq_len(nrow(avg)), times = 2), ]
  draws$draw_id <- rep(1:2, each = nrow(avg))
  draws$estimate <- avg$estimate
  draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "1"] <-
    draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "1"] +
    0.1
  draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "2"] <-
    draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "2"] -
    0.1
  draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "1"] <-
    draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "1"] -
    0.1
  draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "2"] <-
    draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "2"] +
    0.1

  inferred_avg <- avg
  attr(inferred_avg, "draws") <- draws
  attr(inferred_avg, "conf_level") <- 0.95
  attr(inferred_avg, "conf_type") <- "perc"
  attr(inferred_avg, "method") <- "mvn"

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      cmp <- avg_comparisons(
        structure(list(), class = "mock_model", markov_id_var = "id"),
        variables = list(tx = c(0, 1)),
        estimand = "sop",
        times = 1:2,
        state_sets = "1"
      )
    }
  )

  with_mocked_bindings(
    avg_sops = function(...) avg,
    inferences_internal = function(...) inferred_avg,
    {
      out <- mostr:::inferences_avg_comparisons_linear(
        object = cmp,
        method = "mvn",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = NULL,
        cluster = NULL,
        workers = NULL,
        conf_level = 0.95,
        conf_type = "perc",
        return_draws = TRUE,
        use_coefstart = FALSE
      )
    }
  )

  reduced <- attr(out, "draws")
  reduced <- reduced[order(reduced$draw_id, reduced$time), ]
  rownames(reduced) <- NULL

  expect_contains(names(out), c("conf.low", "conf.high", "std.error"))
  expect_equal(reduced$estimate, c(0.4, 0.4, 0.2, 0.2))
  expect_equal(reduced$draw_id, c(1, 1, 2, 2))
})

test_that("comparison inference centers Wald intervals on point comparisons", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))
  draws <- avg[rep(seq_len(nrow(avg)), times = 2), ]
  draws$draw_id <- rep(1:2, each = nrow(avg))
  # Unequal shifts put the mean of the comparison draws above the point
  # comparison, so a draw-mean center would move the interval.
  shift <- c(0.1, 0.3)[draws$draw_id] * (draws$tx == 1)
  draws$estimate <- draws$estimate + ifelse(draws$state == "1", shift, -shift)

  inferred_avg <- avg
  attr(inferred_avg, "draws") <- draws
  attr(inferred_avg, "conf_level") <- 0.95
  attr(inferred_avg, "conf_type") <- "wald"
  attr(inferred_avg, "method") <- "mvn"

  for (estimand in c("sop", "time_in_state")) {
    with_mocked_bindings(
      avg_sops = function(...) avg,
      {
        cmp <- avg_comparisons(
          structure(list(), class = "mock_model", markov_id_var = "id"),
          variables = list(tx = c(0, 1)),
          estimand = estimand,
          times = 1:2,
          state_sets = "1"
        )
      }
    )

    with_mocked_bindings(
      avg_sops = function(...) avg,
      inferences_internal = function(...) inferred_avg,
      {
        out <- mostr:::inferences_avg_comparisons_linear(
          object = cmp,
          method = "simulation",
          engine = "mvn",
          score_weight_dist = "exponential",
          n_draws = 2,
          vcov = NULL,
          cluster = NULL,
          workers = NULL,
          conf_level = 0.95,
          conf_type = "wald",
          return_draws = FALSE,
          use_coefstart = FALSE
        )
      }
    )

    expect_equal(out$estimate, cmp$estimate, info = estimand)
    expect_equal(
      (out$conf.low + out$conf.high) / 2,
      cmp$estimate,
      info = estimand
    )
    expect_equal(
      out$conf.high - out$estimate,
      stats::qnorm(0.975) * out$std.error,
      info = estimand
    )
    expect_gt(min(out$std.error), 0)
  }
})
