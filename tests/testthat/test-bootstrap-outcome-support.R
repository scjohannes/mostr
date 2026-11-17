# A bootstrap resample can hold patients who start in a state (their `yprev`
# at the first visit) without holding any patient who moves into it. The
# state is then a starting state of the resample but not one of its outcome
# states, and the refit has one outcome state fewer. With four states (4
# absorbing), the resamples below hold every patient who never has `y == 2`,
# so state 2 occurs only as a previous state. The refit gives the state
# probability 0 as an outcome, its maximum-likelihood value in that resample;
# a draw that failed instead would leave out exactly the resamples in which
# that probability is 0, and move the lower confidence bound up.
#
# Every test mocks the resampler: draw 1 is the special resample, draw 2
# holds every patient.

outcome_support_case <- function() {
  data <- withr::with_seed(99, make_markov_trajectories(
    n_patients = 300,
    times = 1:5,
    n_states = 4,
    thresholds = c(-1.5, 0, 1.5),
    allowed_start_state = 1:3,
    absorb = 4,
    treatment_effect = 0.5,
    seed = 99
  ) |>
    prepare_markov_data(absorb = 4))
  data <- data[order(data$id, data$time), , drop = FALSE]
  rownames(data) <- NULL
  data$ypprev <- stats::ave(
    as.integer(as.character(data$yprev)),
    data$id,
    FUN = function(x) c(x[1L], utils::head(x, -1L))
  )
  data$ypprev <- factor(data$ypprev, levels = levels(data$yprev))
  ids_all <- unique(data$id)
  list(
    data = data,
    ids_all = ids_all,
    # Patients who never move into state 2; some of them start in it.
    ids_no_entry_2 = setdiff(ids_all, data$id[data$y == 2]),
    # Patients never in state 2 at all.
    ids_no2 = setdiff(ids_all, data$id[data$y == 2 | data$yprev == 2]),
    # Patients never in state 1 as a previous state; a few enter it at their
    # last visit, so state 1 is an outcome but not a previous state.
    ids_no_prev_1 = setdiff(ids_all, data$id[data$yprev == 1])
  )
}

mock_resamples <- function(case, ids) {
  function(data, id_var, n_boot) {
    list(
      data.frame(original_id = ids, new_id = paste0(ids, "_1"), boot_id = 1),
      data.frame(
        original_id = case$ids_all,
        new_id = paste0(case$ids_all, "_2"),
        boot_id = 2
      )
    )
  }
}

fit_support_case <- function(data, engine, order = 1L, ...) {
  rhs <- paste(
    c("time", "tx", "yprev", if (order == 2L) "ypprev"),
    collapse = " + "
  )
  p2_var <- if (order == 2L) "ypprev"
  suppressWarnings(
    if (engine == "orm") {
      mostm(stats::as.formula(paste("y ~", rhs)), data = data, id_var = "id",
        first_followup_time = 1, p2_var = p2_var, absorb = 4)
    } else {
      mostm(stats::as.formula(paste("ordered(y) ~", rhs)), data = data,
        id_var = "id", first_followup_time = 1, p2_var = p2_var, absorb = 4,
        engine = "vglm", ...)
    }
  )
}

vglm_po <- function() VGAM::cumulative(reverse = TRUE, parallel = TRUE)

first_boot_draw <- function(result) {
  draws <- attr(result, "draws")
  draws[draws$draw_id == 1, , drop = FALSE]
}

test_that("the case has starting states in state 2 but no entries into it", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  sample <- case$data[case$data$id %in% case$ids_no_entry_2, ]
  expect_false(any(sample$y == 2))
  expect_gt(sum(sample$time == 1 & sample$yprev == 2), 0)
  no_prev_1 <- case$data[case$data$id %in% case$ids_no_prev_1, ]
  expect_false(any(no_prev_1$yprev == 1))
  expect_gt(sum(no_prev_1$y == 1), 0)
})

test_that("a resample without entries into a starting state keeps its draw, with probability 0", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  local_mocked_bindings(
    fast_group_bootstrap = mock_resamples(case, case$ids_no_entry_2)
  )
  sample <- case$data[case$data$id %in% case$ids_no_entry_2, , drop = FALSE]
  sample$y <- droplevels(sample$y)
  baseline <- sample[sample$time == 1, , drop = FALSE]
  for (engine in c("orm", "vglm")) {
    model <- fit_support_case(case$data, engine, family = vglm_po())
    est <- avg_sops(model, times = 1:3)
    expect_no_warning(
      boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
        return_draws = TRUE)
    )
    expect_identical(attr(boot, "n_successful"), 2L, info = engine)
    draw <- first_boot_draw(boot)
    expect_equal(
      draw$estimate[as.character(draw$state) == "2"], c(0, 0, 0),
      info = engine
    )
    # The other states match a refit made by hand, whose states are 1, 3, 4.
    ref <- as.data.frame(avg_sops(
      fit_support_case(sample, engine, family = vglm_po()),
      newdata = baseline,
      times = 1:3
    ))
    for (s in c("1", "3", "4")) {
      got <- draw[as.character(draw$state) == s, , drop = FALSE]
      want <- ref[as.character(ref$state) == s, , drop = FALSE]
      expect_equal(
        got$estimate[order(got$time)], want$estimate[order(want$time)],
        tolerance = 1e-6, info = paste(engine, "state", s)
      )
    }
    sums <- as.vector(tapply(draw$estimate, draw$time, sum))
    expect_equal(sums, c(1, 1, 1), tolerance = 1e-10, info = engine)
  }
})

test_that("bootstrap comparisons keep a resample without entries into a starting state", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  local_mocked_bindings(
    fast_group_bootstrap = mock_resamples(case, case$ids_no_entry_2)
  )
  model <- fit_support_case(case$data, "vglm", family = vglm_po())
  est <- avg_comparisons(model, variables = list(tx = 0:1), times = 1:2)
  expect_no_warning(
    boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
      return_draws = TRUE)
  )
  expect_identical(attr(boot, "n_successful"), 2L)
  draw <- first_boot_draw(boot)
  expect_equal(draw$estimate[as.character(draw$state_set) == "2"], c(0, 0))
})

test_that("a second-order bootstrap with a starting state it never enters is refused before refitting", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  local_mocked_bindings(
    fast_group_bootstrap = mock_resamples(case, case$ids_no_entry_2)
  )
  for (engine in c("orm", "vglm")) {
    model <- fit_support_case(case$data, engine, order = 2L, family = vglm_po())
    est <- avg_sops(model, times = 1:2)
    refits <- 0L
    local_mocked_bindings(
      update_bootstrap_model = function(...) {
        refits <<- refits + 1L
        stop("no refit expected")
      }
    )
    err <- expect_error(
      inferences(est, method = "bootstrap", n_draws = 2, seed = 1),
      class = "mostr_error"
    )
    expect_match(conditionMessage(err), "state 2", fixed = TRUE, info = engine)
    expect_match(conditionMessage(err), "1 of the 2 resamples", fixed = TRUE,
      info = engine)
    expect_match(conditionMessage(err), "method = \"fwb\"", fixed = TRUE,
      info = engine)
    expect_identical(refits, 0L, info = engine)
  }
})

test_that("explicit vglm constraints with a resample that lacks an outcome state are refused", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  constraints <- list(
    "(Intercept)" = diag(3),
    time = matrix(1, 3, 1),
    tx = cbind(1, c(0, 0, 1)),
    yprev = matrix(1, 3, 1)
  )
  model <- fit_support_case(
    case$data, "vglm",
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    constraints = constraints
  )
  est <- avg_sops(model, times = 1:2)
  # State 2 missing from the outcomes only, and missing altogether.
  for (ids in list(case$ids_no_entry_2, case$ids_no2)) {
    local_mocked_bindings(fast_group_bootstrap = mock_resamples(case, ids))
    refits <- 0L
    local_mocked_bindings(
      update_bootstrap_model = function(...) {
        refits <<- refits + 1L
        stop("no refit expected")
      }
    )
    err <- expect_error(
      inferences(est, method = "bootstrap", n_draws = 2, seed = 1),
      class = "mostr_error"
    )
    expect_match(conditionMessage(err), "constraints", fixed = TRUE)
    expect_match(conditionMessage(err), "method = \"fwb\"", fixed = TRUE)
    expect_error(bootstrap_model_coefs(model, n_draws = 2), "constraints",
      class = "mostr_error")
    expect_identical(refits, 0L)
  }
})

test_that("bootstrap_model_coefs() gives NA for contrasts whose reference level a resample lacks", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  local_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      ids <- list(case$ids_all, case$ids_no_prev_1)
      lapply(1:2, function(b) {
        data.frame(
          original_id = ids[[b]],
          new_id = paste0(ids[[b]], "_", b),
          boot_id = b
        )
      })
    }
  )
  sample <- case$data[case$data$id %in% case$ids_no_prev_1, , drop = FALSE]
  sample$yprev <- droplevels(sample$yprev)
  for (engine in c("orm", "vglm")) {
    model <- fit_support_case(case$data, engine, family = vglm_po())
    warnings <- character()
    coefs <- withCallingHandlers(
      bootstrap_model_coefs(model, n_draws = 2),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_length(warnings, 1L)
    expect_match(warnings, "1 of the 2", fixed = TRUE, info = engine)
    full <- stats::coef(model)
    # Row 1 holds every patient and keeps every coefficient.
    expect_equal(unlist(coefs[1, names(full)]), full, tolerance = 1e-6,
      info = engine)
    # Row 2 lacks the reference level 1 of `yprev`: its intercepts and
    # previous-state coefficients compare with another level, so they are NA.
    prev_cols <- grep("^yprev", names(full), value = TRUE)
    intercept_cols <- grep("Intercept|^y>=", names(full), value = TRUE)
    expect_true(all(is.na(unlist(coefs[2, c(prev_cols, intercept_cols)]))),
      info = engine)
    # Coefficients of other terms are kept and match a refit made by hand.
    ref <- stats::coef(fit_support_case(sample, engine, family = vglm_po()))
    expect_equal(unlist(coefs[2, c("time", "tx")]), ref[c("time", "tx")],
      tolerance = 1e-6, info = engine)
  }
})

test_that("a supplied profile in a state a resample only enters is refused", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  local_mocked_bindings(
    fast_group_bootstrap = mock_resamples(case, case$ids_no_prev_1)
  )
  profiles <- case$data[case$data$time == 1 & case$data$yprev == 1, ][1:20, ]
  model <- fit_support_case(case$data, "vglm", family = vglm_po())
  est <- avg_sops(model, newdata = profiles, times = 1:2)
  expect_error(
    inferences(est, method = "bootstrap", n_draws = 2, seed = 1),
    "state 1 in 1 of the 2 resamples"
  )
})

test_that("the R engine refuses `y_levels` that do not match the fitted states", {
  skip_on_cran() # The pass depends on the seed.
  case <- outcome_support_case()
  # A fit with two outcome states, 3 and 4, against four states in `y_levels`:
  # at one visit the probabilities used to be recycled and summed to 2.
  sample <- case$data[case$data$y %in% c(3, 4), , drop = FALSE]
  sample$y <- droplevels(sample$y)
  fit <- fit_support_case(sample, "vglm", family = vglm_po())
  settings <- markov_settings_update(
    markov_model_settings(fit),
    y_levels = state_factor(1:4),
    absorb = 4
  )
  for (times in list(1, 1:3)) {
    expect_error(
      soprob_markov_r_engine(fit, newdata = sample[sample$time == 1, ][1:3, ],
        times = times, settings = settings),
      "`y_levels` defines 4 states, but the fitted model defines 2 states",
      fixed = TRUE
    )
  }
})
