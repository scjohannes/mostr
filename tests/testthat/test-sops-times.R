# `times` lists consecutive visits: the SOP recursion takes one transition per
# value. Requests that skip fitted visits silently change what each result
# means, so the user-facing SOP functions warn about them.

# Made once per test run for each `time_as_factor`.
make_daily_visit_case <- function(time_as_factor = FALSE) {
  skip_if_not_installed("rms")
  grid_fixture(
    paste("daily visit case, factor time", time_as_factor),
    function() build_daily_visit_case(time_as_factor)
  )
}

build_daily_visit_case <- function(time_as_factor) {
  data <- make_test_data(n_patients = 80, follow_up_time = 8, seed = 2612)
  if (time_as_factor) {
    data$time <- factor(data$time, levels = as.character(1:8))
  }
  list(
    data = data,
    model = mostm(
      y ~ time + tx + yprev,
      data = data,
      id_var = "id",
      first_followup_time = if (time_as_factor) "1"
    )
  )
}

# The daily visits fitted by vglm with an explicit spline basis for time, which
# prediction takes from `time_covariates`. Made once per test run.
make_time_covariate_visit_case <- function() {
  case <- make_daily_visit_case()
  grid_fixture("daily visit case, vglm with time covariates", function() {
    list(
      model = mostm(
        ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = case$data,
        id_var = "id",
        engine = "vglm"
      ),
      time_covariates = make_time_covariates(
        case$data,
        "time",
        "time_lin",
        "time_nlin_1"
      )
    )
  })
}

# Visits on days 1, 2, 4, and 8, with `gap` holding the days since the
# previous visit (day 0 for the first visit).
make_irregular_visit_case <- function(visit_days = c(1, 2, 4, 8)) {
  skip_if_not_installed("rms")
  raw <- make_markov_trajectories(
    n_patients = 120,
    follow_up_time = max(visit_days),
    seed = 2613
  )
  visits <- raw[raw$time %in% visit_days, , drop = FALSE]
  visits <- visits[order(visits$id, visits$time), , drop = FALSE]
  first <- !duplicated(visits$id)
  previous <- c(NA, seq_len(nrow(visits) - 1L))
  visits$yprev[!first] <- visits$y[previous[!first]]
  visits$gap <- ifelse(
    first,
    visits$time,
    visits$time - visits$time[previous]
  )
  prepare_markov_data(visits)
}

# The irregular visits fitted without and with the gap, made once per test run.
make_irregular_visit_fits <- function() {
  grid_fixture("irregular visit fits", function() {
    data <- make_irregular_visit_case()
    list(
      plain = mostm(y ~ time + tx + yprev, data = data, id_var = "id"),
      gap = mostm(y ~ time + gap + tx + yprev, data = data, id_var = "id")
    )
  })
}

count_skipped_visit_warnings <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(
    expr,
    mostr_warning_skipped_visits = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, n = n)
}

death_at <- function(x, time, tx) {
  x$estimate[x$state == "6" & x$time == time & x$tx == tx]
}

test_that("avg_sops() warns when `times` skips fitted visits", {
  case <- make_daily_visit_case()
  variables <- list(tx = c(0, 1))

  expect_no_warning(
    full <- avg_sops(case$model, variables = variables, times = 1:8, absorb = "6")
  )
  expect_warning(
    sparse <- avg_sops(
      case$model,
      variables = variables,
      times = c(1, 4, 8),
      absorb = "6"
    ),
    class = "mostr_warning_skipped_visits"
  )

  # The result labelled 8 is the probability after three transitions, not at
  # visit 8; subsetting the consecutive request gives the visit-8 estimate.
  expect_equal(death_at(sparse, 1, 0), death_at(full, 1, 0))
  expect_gt(death_at(full, 8, 0) - death_at(sparse, 8, 0), 0.1)
})

test_that("SOP functions warn when `times` skips fitted visits, and only then", {
  skip_on_cran()
  daily <- make_daily_visit_case()
  daily_factor <- make_daily_visit_case(time_as_factor = TRUE)
  covariates <- make_time_covariate_visit_case()
  irregular <- make_irregular_visit_fits()
  # Each supplied row's previous state is the state before the first requested
  # visit, so supplied profiles may start at a later visit.
  profiles <- daily$data[daily$data$time == 1, , drop = FALSE]
  skipped <- "mostr_warning_skipped_visits"

  # `phrases` must all appear in the warning; `class` is the warning's class.
  # A case with neither expects no warning at all.
  cases <- list(
    list(
      label = "avg_sops(): the warning explains `times` and the remedy",
      fun = avg_sops,
      model = daily$model,
      args = list(times = c(1, 4, 8)),
      class = skipped,
      phrases = c(
        "skips fitted visits 2, 3, 5, 6, and 7",
        "labelled 8 is the probability after 3 transitions",
        "`times = 1:8`",
        "`gap_var`"
      )
    ),
    # Without `newdata`, each patient starts from the state before visit 1, so
    # `times` must begin at visit 1.
    list(
      label = "sops(): stored profiles, one late visit",
      fun = sops,
      model = daily$model,
      args = list(times = 8),
      phrases = "skips fitted visits 1, 2, 3, 4, 5, and 2 more"
    ),
    list(
      label = "sops(): stored profiles, late start",
      fun = sops,
      model = daily$model,
      args = list(times = 4:8),
      class = skipped
    ),
    list(
      label = "sops(): stored profiles, start before the first visit",
      fun = sops,
      model = daily$model,
      args = list(times = 0:8),
      phrases = "starts at 0, before the first follow-up visit (1)"
    ),
    list(
      label = "sops(): unordered times",
      fun = sops,
      model = daily$model,
      args = list(times = c(1, 3, 2, 4:8)),
      phrases = "not in increasing order"
    ),
    list(
      label = "sops(): every fitted visit",
      fun = sops,
      model = daily$model,
      args = list(times = 1:8)
    ),
    list(
      label = "sops(): supplied profiles starting later",
      fun = sops,
      model = daily$model,
      args = list(newdata = profiles, times = 4:8)
    ),
    list(
      label = "sops(): supplied profiles skipping visits",
      fun = sops,
      model = daily$model,
      args = list(newdata = profiles, times = c(4, 8)),
      phrases = "skips fitted visits 5, 6, and 7"
    ),
    list(
      label = "avg_sops(): supplied profiles skipping visits",
      fun = avg_sops,
      model = daily$model,
      args = list(newdata = profiles, times = c(1, 4, 8)),
      class = skipped
    ),
    # Times beyond the fitted visits are not skipped visits.
    list(
      label = "avg_sops(): times beyond the fitted visits",
      fun = avg_sops,
      model = daily$model,
      args = list(times = 1:12)
    ),
    list(
      label = "avg_sops(): skipped visits before times beyond the fitted visits",
      fun = avg_sops,
      model = daily$model,
      args = list(times = c(1, 4, 8, 12)),
      class = skipped
    ),
    # Explicit time covariates do not change the check.
    list(
      label = "avg_sops(): time covariates for every visit",
      fun = avg_sops,
      model = covariates$model,
      args = list(times = 1:8, time_covariates = covariates$time_covariates)
    ),
    list(
      label = "avg_sops(): time covariates for skipped visits",
      fun = avg_sops,
      model = covariates$model,
      args = list(
        times = c(1, 4, 8),
        time_covariates = covariates$time_covariates[c(1, 4, 8), , drop = FALSE]
      ),
      class = skipped
    ),
    # Irregular visit schedules are checked against the fitted visits.
    list(
      label = "avg_sops(): every irregular visit",
      fun = avg_sops,
      model = irregular$plain,
      args = list(times = c(1, 2, 4, 8))
    ),
    list(
      label = "avg_sops(): skipped irregular visit",
      fun = avg_sops,
      model = irregular$plain,
      args = list(times = c(1, 4, 8)),
      phrases = "skips fitted visit 2."
    ),
    list(
      label = "avg_sops(): skipped irregular visit before later times",
      fun = avg_sops,
      model = irregular$plain,
      args = list(times = c(1, 4, 8, 12)),
      phrases = "request every fitted visit from 1 to 8, then the later times"
    ),
    # `gap_var` allows times that skip fitted visits.
    list(
      label = "avg_sops(): gap_var",
      fun = avg_sops,
      model = irregular$gap,
      args = list(times = c(1, 4, 8), gap_var = "gap")
    ),
    list(
      label = "sops(): gap_var",
      fun = sops,
      model = irregular$gap,
      args = list(times = c(1, 4, 8), gap_var = "gap")
    ),
    # Factor visit times warn when levels are skipped.
    list(
      label = "avg_sops(): every factor level",
      fun = avg_sops,
      model = daily_factor$model,
      args = list(times = as.character(1:8))
    ),
    list(
      label = "avg_sops(): skipped factor levels",
      fun = avg_sops,
      model = daily_factor$model,
      args = list(times = c("1", "4", "8")),
      phrases = "skips fitted visits 2, 3, 5, 6, and 7"
    ),
    # comparisons() and state_time() check `times` against the stored starting
    # profiles too.
    list(
      label = "comparisons(): stored profiles, late start",
      fun = comparisons,
      model = daily$model,
      args = list(variables = "tx", times = 4:8),
      phrases = "skips fitted visits 1, 2, and 3"
    ),
    list(
      label = "comparisons(time_benefit): stored profiles, late start",
      fun = comparisons,
      model = daily$model,
      args = list(variables = "tx", estimand = "time_benefit", times = 4:8),
      phrases = "skips fitted visits 1, 2, and 3"
    ),
    list(
      label = "state_time(): stored profiles, late start",
      fun = state_time,
      model = daily$model,
      args = list(times = 4:8),
      phrases = "skips fitted visits 1, 2, and 3"
    ),
    list(
      label = "state_time(): start before the first visit",
      fun = state_time,
      model = daily$model,
      args = list(variables = "tx", times = 0:8),
      phrases = "starts at 0, before the first follow-up visit (1)"
    ),
    list(
      label = "comparisons(): supplied profiles starting later",
      fun = comparisons,
      model = daily$model,
      args = list(newdata = profiles, variables = "tx", times = 4:8)
    ),
    list(
      label = "state_time(): supplied profiles starting later",
      fun = state_time,
      model = daily$model,
      args = list(newdata = profiles, times = 4:8)
    ),
    list(
      label = "state_time(): supplied profiles skipping visits",
      fun = state_time,
      model = daily$model,
      args = list(newdata = profiles, times = c(4, 8)),
      phrases = "skips fitted visits 5, 6, and 7"
    )
  )

  for (case in cases) {
    run <- function() {
      do.call(case$fun, c(list(case$model), case$args, list(absorb = "6")))
    }
    if (is.null(case$class) && length(case$phrases) == 0L) {
      expect_warning(run(), NA, info = case$label)
      next
    }
    warning <- if (is.null(case$class)) {
      expect_warning(run(), case$phrases[1], fixed = TRUE, info = case$label)
    } else {
      expect_warning(run(), class = case$class, info = case$label)
    }
    for (phrase in case$phrases) {
      expect_match(
        conditionMessage(warning),
        phrase,
        fixed = TRUE,
        info = case$label
      )
    }
  }
})

test_that("summaries warn about skipped visits once and not again during inference", {
  skip_on_cran()
  case <- make_daily_visit_case()
  model <- case$model
  variables <- list(tx = c(0, 1))
  delta <- function(x) inferences(x, method = "delta")
  mvn <- function(x) inferences(x, method = "mvn", n_draws = 5)

  # `warnings` is the number of skipped-visit warnings from the call; `infer`,
  # if given, must then run on the result without any warning. A case with no
  # skipped-visit warning must give no warning at all.
  cases <- list(
    list(
      label = "avg_state_time()",
      run = function() {
        avg_state_time(model, variables = variables, times = c(1, 4, 8), absorb = "6")
      },
      warnings = 1L,
      infer = delta
    ),
    list(
      label = "avg_comparisons()",
      run = function() {
        avg_comparisons(model, variables = variables, times = c(1, 4, 8), absorb = "6")
      },
      warnings = 1L,
      infer = delta
    ),
    list(
      label = "avg_comparisons(time_benefit)",
      run = function() {
        avg_comparisons(
          model,
          variables = variables,
          estimand = "time_benefit",
          times = c(1, 4, 8),
          absorb = "6"
        )
      },
      warnings = 1L,
      infer = mvn
    ),
    list(
      label = "comparisons(), late start",
      run = function() {
        comparisons(model, variables = "tx", times = 4:8, absorb = "6")
      },
      warnings = 1L,
      infer = delta
    ),
    list(
      label = "comparisons(), skipped visits",
      run = function() {
        comparisons(model, variables = "tx", times = c(1, 4, 8), absorb = "6")
      },
      warnings = 1L
    ),
    list(
      label = "state_time(), late start",
      run = function() state_time(model, times = 4:8, absorb = "6"),
      warnings = 1L,
      infer = mvn
    ),
    list(
      label = "avg_state_time(), every visit",
      run = function() {
        avg_state_time(model, variables = variables, times = 1:8, absorb = "6")
      },
      warnings = 0L
    ),
    list(
      label = "comparisons(), every visit",
      run = function() {
        comparisons(model, variables = "tx", times = 1:8, absorb = "6")
      },
      warnings = 0L
    ),
    list(
      label = "state_time(), every visit",
      run = function() state_time(model, times = 1:8, absorb = "6"),
      warnings = 0L
    )
  )

  for (case in cases) {
    if (case$warnings == 0L) {
      expect_warning(case$run(), NA, info = case$label)
      next
    }
    counted <- count_skipped_visit_warnings(case$run())
    expect_identical(counted$n, case$warnings, info = case$label)
    if (!is.null(case$infer)) {
      expect_warning(case$infer(counted$value), NA, info = case$label)
    }
  }
})
