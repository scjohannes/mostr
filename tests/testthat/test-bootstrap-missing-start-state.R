# An ordinary bootstrap resample can lack a state that is not the absorbing
# one. With four states (4 absorbing) a resample can lack state 2; its refit
# then knows only states 1, 3 and 4. The refit keeps those labels, so supplied
# profiles (`newdata`) are predicted from their own states, and state 2 gets
# probability 0. A refit without state 2 can't predict for a profile that
# starts there, so `inferences()` refuses such a bootstrap before refitting.
#
# Every test mocks the resampler: draw 1 holds only the patients never in
# state 2, draw 2 holds every patient.

missing_state_case <- function() {
  data <- withr::with_seed(99, make_markov_trajectories(
    n_patients = 300,
    follow_up_time = 5,
    n_states = 4,
    thresholds = c(-1.5, 0, 1.5),
    allowed_start_state = 1:3,
    absorbing_state = 4,
    treatment_effect = 0.5,
    seed = 99
  ) |>
    prepare_markov_data(absorbing_state = 4))
  data <- data[order(data$id, data$time), , drop = FALSE]
  rownames(data) <- NULL
  data$ypprev <- stats::ave(
    as.integer(as.character(data$yprev)),
    data$id,
    FUN = function(x) c(x[1L], utils::head(x, -1L))
  )
  data$ypprev <- factor(data$ypprev, levels = levels(data$yprev))
  ids_all <- unique(data$id)
  ids_no2 <- setdiff(ids_all, data$id[data$y == 2 | data$yprev == 2])
  list(data = data, ids_all = ids_all, ids_no2 = ids_no2)
}

resample_without_state_2 <- function(case) {
  function(data, id_var, n_boot) {
    list(
      data.frame(
        original_id = case$ids_no2,
        new_id = paste0(case$ids_no2, "_1"),
        boot_id = 1
      ),
      data.frame(
        original_id = case$ids_all,
        new_id = paste0(case$ids_all, "_2"),
        boot_id = 2
      )
    )
  }
}

fit_case <- function(data, engine, order = 1L) {
  rhs <- paste(c("time", "tx", "yprev", if (order == 2L) "ypprev"), collapse = " + ")
  suppressWarnings(
    if (engine == "orm") {
      mostm(stats::as.formula(paste("y ~", rhs)), data = data, id_var = "id",
        first_followup_time = 1)
    } else {
      mostm(stats::as.formula(paste("ordered(y) ~", rhs)), data = data,
        id_var = "id", first_followup_time = 1, engine = "vglm",
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE))
    }
  )
}

# The refit of draw 1, made by hand with the states renumbered
# 1 -> 1, 3 -> 2, 4 -> 3. For factor previous states this is the same model as
# a refit on the original labels, so the two must give the same SOPs.
renumber <- c("1" = 1L, "3" = 2L, "4" = 3L)
renumber_states <- function(data, cols = c("yprev", "ypprev")) {
  for (col in intersect(cols, names(data))) {
    data[[col]] <- factor(renumber[as.character(data[[col]])], levels = 1:2)
  }
  data
}
hand_refit <- function(case, engine, order = 1L) {
  s <- case$data[case$data$id %in% case$ids_no2, , drop = FALSE]
  s$y <- renumber[as.character(s$y)]
  fit_case(renumber_states(s), engine, order)
}
# Hand-refit SOPs on the original state labels, state 2 having probability 0.
hand_sops <- function(ref_model, profiles, times, p2_var = NULL, by = NULL) {
  out <- as.data.frame(avg_sops(ref_model, newdata = renumber_states(profiles),
    times = times, absorb = 3, p2_var = p2_var, by = by %||% TRUE))
  groups <- if (is.null(by)) list(NULL) else sort(unique(profiles[[by]]))
  pieces <- lapply(groups, function(group) {
    full <- expand.grid(time = times, state = 1:4)
    full$estimate <- 0
    rows <- if (is.null(by)) out else out[out[[by]] == group, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      state <- c(1, 3, 4)[as.integer(as.character(rows$state[i]))]
      full$estimate[full$time == rows$time[i] & full$state == state] <- rows$estimate[i]
    }
    if (!is.null(by)) full[[by]] <- group
    full
  })
  do.call(rbind, pieces)
}

first_draw <- function(result) {
  draws <- attr(result, "draws")
  draws[draws$draw_id == 1, , drop = FALSE]
}

start_profiles <- function(case, without_state) {
  profiles <- case$data[case$data$time == 1, , drop = FALSE]
  profiles <- profiles[profiles$yprev != without_state, , drop = FALSE][1:40, ]
  rownames(profiles) <- NULL
  profiles
}

test_that("supplied profiles in states the resample has are predicted from those states", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  profiles <- start_profiles(case, without_state = 2)
  for (order in 1:2) {
    p2 <- if (order == 2L) "ypprev" else NULL
    model <- fit_case(case$data, "vglm", order)
    est <- avg_sops(model, newdata = profiles, times = 1:2, absorb = 4, p2_var = p2)
    expect_no_warning(
      boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
        return_draws = TRUE)
    )
    expect_identical(attr(boot, "n_successful"), 2L, info = paste("order", order))
    expected <- hand_sops(hand_refit(case, "vglm", order), profiles, 1:2, p2)
    draw <- first_draw(boot)
    draw <- draw[order(draw$state, draw$time), ]
    expect_equal(draw$estimate, expected$estimate, tolerance = 1e-6, info = paste("order", order))
  }
})

test_that("a supplied profile in a state some resample lacks is refused before refitting", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  profiles <- start_profiles(case, without_state = 3)
  for (engine in c("vglm", "orm")) {
    model <- fit_case(case$data, engine)
    est <- avg_sops(model, newdata = profiles, times = 1:2, absorb = 4)
    refits <- 0L
    local_mocked_bindings(
      update_bootstrap_model = function(...) {
        refits <<- refits + 1L
        stop("no refit expected")
      }
    )
    err <- expect_error(
      inferences(est, method = "bootstrap", n_draws = 2, seed = 1),
      "state 2 in 1 of the 2 resamples"
    )
    expect_match(conditionMessage(err), "method = \"fwb\"", fixed = TRUE, info = engine)
    expect_identical(refits, 0L, info = engine)
  }
})

test_that("an orm refit of a resample without a middle state succeeds", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  model <- fit_case(case$data, "orm")
  est <- avg_sops(model, times = 1:2, absorb = 4)
  expect_no_warning(
    boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
      return_draws = TRUE)
  )
  expect_identical(attr(boot, "n_successful"), 2L)
  baseline <- case$data[case$data$id %in% case$ids_no2 & case$data$time == 1, ]
  expected <- hand_sops(hand_refit(case, "orm"), baseline, 1:2)
  draw <- first_draw(boot)
  draw <- draw[order(draw$state, draw$time), ]
  expect_equal(draw$estimate, expected$estimate, tolerance = 1e-5)
})

test_that("a refit without a middle state keeps a numeric previous state's numbers", {
  # Five states, 5 absorbing. Without state 2 the previous states are 1, 3
  # and 4; renumbering them 1, 2, 3 would change the fitted linear yprev
  # effect, so the refit must see the original numbers.
  data <- withr::with_seed(99, make_markov_trajectories(
    n_patients = 300,
    follow_up_time = 5,
    n_states = 5,
    thresholds = c(-2, -0.7, 0.7, 2),
    allowed_start_state = 1:4,
    absorbing_state = 5,
    treatment_effect = 0.5,
    seed = 99
  ) |>
    prepare_markov_data(absorbing_state = 5))
  data$yprev <- as.integer(as.character(data$yprev))
  case <- list(
    ids_all = unique(data$id),
    ids_no2 = setdiff(unique(data$id), data$id[data$y == 2 | data$yprev == 2])
  )
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  model <- fit_case(data, "orm")
  est <- avg_sops(model, times = 1:2, absorb = 5)
  expect_no_warning(
    boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
      return_draws = TRUE)
  )
  sample <- data[data$id %in% case$ids_no2, , drop = FALSE]
  expect_setequal(unique(sample$yprev), c(1L, 3L, 4L))
  sample$y <- droplevels(sample$y)
  ref <- fit_case(sample, "orm")
  baseline <- sample[sample$time == 1, , drop = FALSE]
  expected <- as.data.frame(avg_sops(ref, newdata = baseline, times = 1:2, absorb = 5))
  draw <- first_draw(boot)
  draw <- draw[as.character(draw$state) != "2", , drop = FALSE]
  draw <- draw[order(as.character(draw$state), draw$time), ]
  expected <- expected[order(as.character(expected$state), expected$time), ]
  expect_equal(draw$estimate, expected$estimate, tolerance = 1e-5)
})

test_that("averages within `by` groups keep the labels of a resample without a middle state", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  model <- fit_case(case$data, "vglm")
  est <- avg_sops(model, times = 1:2, absorb = 4, by = "tx")
  expect_no_warning(
    boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
      return_draws = TRUE)
  )
  baseline <- case$data[case$data$id %in% case$ids_no2 & case$data$time == 1, ]
  expected <- hand_sops(hand_refit(case, "vglm"), baseline, 1:2, by = "tx")
  expected <- expected[order(expected$tx, expected$state, expected$time), ]
  draw <- first_draw(boot)
  draw <- draw[order(draw$tx, as.character(draw$state), draw$time), ]
  expect_equal(draw$estimate, expected$estimate, tolerance = 1e-6)
})

test_that("the bootstrap refuses `variables` that set the state two visits back", {
  case <- missing_state_case()
  model <- fit_case(case$data, "vglm", order = 2L)
  levels <- levels(case$data$ypprev)
  for (values in list(c(1, 3), c(1, 2))) {
    est <- avg_sops(model, variables = list(ypprev = factor(values, levels = levels)),
      times = 1, absorb = 4, p2_var = "ypprev")
    for (method in c("bootstrap", "fwb")) {
      expect_error(
        inferences(est, method = method, n_draws = 2, seed = 1),
        "Setting the starting state with `variables` supports only conditional",
        info = paste(method, paste(values, collapse = ","))
      )
    }
  }
})

test_that("the starting-state check also covers states that `variables` sets", {
  case <- missing_state_case()
  resamples <- resample_without_state_2(case)(case$data, "id", 2)
  check <- function(...) {
    check_bootstrap_starting_states(
      resamples,
      refit_data = case$data,
      id_var = "id",
      state_cols = c("y", "yprev"),
      ...
    )
  }
  expect_silent(check(starting_states = list(yprev = c(1, 3))))
  expect_error(
    check(starting_states = list(yprev = factor(c(1, 2), levels = 1:3))),
    "state 2 in 1 of the 2 resamples"
  )
  expect_silent(check(starting_states = list()))
})

test_that("bootstrap comparisons on supplied profiles use the resample's own states", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  profiles <- start_profiles(case, without_state = 2)
  model <- fit_case(case$data, "vglm")
  est <- avg_comparisons(model, newdata = profiles, variables = list(tx = 0:1),
    times = 1:2, absorb = 4)
  expect_no_warning(
    boot <- inferences(est, method = "bootstrap", n_draws = 2, seed = 1,
      return_draws = TRUE)
  )
  expect_identical(attr(boot, "n_successful"), 2L)
  ref <- as.data.frame(avg_comparisons(hand_refit(case, "vglm"),
    newdata = renumber_states(profiles), variables = list(tx = 0:1),
    times = 1:2, absorb = 3))
  draw <- first_draw(boot)
  for (s in c(1, 3, 4)) {
    got <- draw$estimate[as.character(draw$state_set) == as.character(s)]
    want <- ref$estimate[as.character(ref$state_set) == as.character(renumber[[as.character(s)]])]
    expect_equal(got, want, tolerance = 1e-6, info = paste("state", s))
  }
  expect_equal(draw$estimate[as.character(draw$state_set) == "2"], c(0, 0))
})

test_that("bootstrap time benefit checks and keeps the states of supplied profiles", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  model <- fit_case(case$data, "vglm")
  benefit <- function(profiles) {
    avg_comparisons(model, newdata = profiles, variables = list(tx = 0:1),
      times = 1:2, absorb = 4, estimand = "time_benefit")
  }
  expect_error(
    inferences(benefit(start_profiles(case, without_state = 3)),
      method = "bootstrap", n_draws = 2, seed = 1),
    "state 2 in 1 of the 2 resamples"
  )
  profiles <- start_profiles(case, without_state = 2)
  expect_no_warning(
    boot <- inferences(benefit(profiles), method = "bootstrap", n_draws = 2,
      seed = 1, return_draws = TRUE)
  )
  expect_identical(attr(boot, "n_successful"), 2L)
  ref <- as.data.frame(avg_comparisons(hand_refit(case, "vglm"),
    newdata = renumber_states(profiles), variables = list(tx = 0:1),
    times = 1:2, absorb = 3, estimand = "time_benefit"))
  draw <- first_draw(boot)
  expect_equal(draw$estimate, ref$estimate, tolerance = 1e-6)
})

test_that("dropped bootstrap draws are summed up in one warning", {
  case <- missing_state_case()
  local_mocked_bindings(fast_group_bootstrap = resample_without_state_2(case))
  model <- fit_case(case$data, "vglm")
  est <- avg_sops(model, times = 1:2, absorb = 4)
  real_soprob <- soprob_markov
  calls <- 0L
  local_mocked_bindings(
    soprob_markov = function(...) {
      calls <<- calls + 1L
      if (calls == 1L) stop("planted failure")
      real_soprob(...)
    }
  )
  warnings <- character()
  boot <- withCallingHandlers(
    inferences(est, method = "bootstrap", n_draws = 2, seed = 1),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(attr(boot, "n_successful"), 1L)
  expect_length(warnings, 2L)
  expect_match(warnings[1], "soprob_markov failed: planted failure", fixed = TRUE)
  expect_match(warnings[2], "1 of the 2 bootstrap draws failed", fixed = TRUE)
  expect_match(warnings[2], "draw 1", fixed = TRUE)
  expect_match(warnings[2], "planted failure", fixed = TRUE)
})
