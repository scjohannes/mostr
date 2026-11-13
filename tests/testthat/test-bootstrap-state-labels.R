# State labels are names, so relabelling the states (keeping their order) must
# not change any result, including bootstrap refits. Bootstrap resamples can
# lack a state; the refit then drops it and keeps the labels of the rest, which
# must work for letter labels on factor columns as it does for numbers.

letter_case <- function() {
  case <- make_comparisons_case()
  data <- case$data
  labels <- c("1" = "a", "2" = "b", "3" = "c")
  data$y <- factor(labels[as.character(data$y)], levels = c("a", "b", "c"))
  data$yprev <- factor(labels[as.character(data$yprev)], levels = c("a", "b", "c"))
  model <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    absorb = "c",
    engine = "vglm"
  )
  list(numbers = case, letters = list(model = model, data = data))
}

test_that("bootstrap results don't depend on whether states are numbers or letters", {
  skip_on_cran()
  cases <- letter_case()
  numbers <- avg_sops(cases$numbers$model, times = 1:3)
  letters <- avg_sops(cases$letters$model, times = 1:3)
  expect_equal(letters$estimate, numbers$estimate, tolerance = 1e-8)

  skip_on_cran()
  for (method in c("bootstrap", "fwb")) {
    with_numbers <- inferences(numbers, method = method, n_draws = 5, seed = 7)
    expect_no_warning(
      with_letters <- inferences(letters, method = method, n_draws = 5, seed = 7)
    )
    expect_equal(
      with_letters$std.error,
      with_numbers$std.error,
      tolerance = 1e-6,
      info = method
    )
  }
})

# Every resample below holds only the patients who never reach state 3, so
# each refit lacks that state.
resample_never_absorbed <- function(data) {
  ids <- setdiff(unique(data$id), data$id[data$y == 3])
  function(data, id_var, n_boot) {
    replicate(
      n_boot,
      data.frame(original_id = ids, new_id = paste0(ids, "_1"), boot_id = 1),
      simplify = FALSE
    )
  }
}

# Such resamples also lose the absorbing state, which is reported.
muffle_absent_absorbing_state <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Absorbing state(s) 3 not present", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("a resample that lacks a state is relabelled whatever the response is called", {
  skip_on_cran()
  case <- make_comparisons_case()
  renamed <- case$data
  names(renamed)[names(renamed) == "y"] <- "state"
  renamed_model <- mostm(
    ordered(state) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = renamed,
    first_followup_time = 1,
    id_var = "id",
    absorb = 3,
    engine = "vglm"
  )
  with_y <- avg_sops(case$model, times = 1:2)
  with_state <- avg_sops(renamed_model, times = 1:2)
  expect_equal(with_state$estimate, with_y$estimate, tolerance = 1e-8)

  local_mocked_bindings(fast_group_bootstrap = resample_never_absorbed(case$data))
  boot_y <- muffle_absent_absorbing_state(
    inferences(with_y, method = "bootstrap", n_draws = 2, seed = 1)
  )
  boot_state <- muffle_absent_absorbing_state(
    inferences(with_state, method = "bootstrap", n_draws = 2, seed = 1)
  )
  expect_equal(boot_state$conf.low, boot_y$conf.low, tolerance = 1e-8)
  expect_equal(boot_state$conf.high, boot_y$conf.high, tolerance = 1e-8)
  # The absent state has probability 0 in every refit.
  expect_true(all(boot_y$conf.high[boot_y$state == 3] == 0))
})

# Covariate factors are not states: their levels must neither hide a missing
# state nor be renumbered.
test_that("bootstrap_model_coefs() relabels only the response and yprev", {
  skip_on_cran()
  case <- make_comparisons_case()
  data <- case$data
  data$tx <- factor(data$tx)
  fit <- mostm(
    ordered(y) ~ time + tx + marker + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    first_followup_time = 1,
    id_var = "id",
    absorb = 3,
    engine = "vglm"
  )
  local_mocked_bindings(fast_group_bootstrap = resample_never_absorbed(data))
  as_factor <- bootstrap_model_coefs(fit, n_boot = 1)
  as_number <- bootstrap_model_coefs(case$model, n_boot = 1)
  expect_equal(as_factor$tx1, as_number$tx, tolerance = 1e-8)
  # One threshold is left of the original two.
  expect_true(is.finite(as_factor[["(Intercept)"]]))
  expect_false("(Intercept):2" %in% names(as_factor))
})

# The previous-state column is the one named in `mostm(p_var = )`, whatever it
# is called. A resample that lacks a state must drop that state's level from
# it, as from a column called `yprev`; `rms::orm()` cannot fit with an empty
# level.
test_that("bootstrap_model_coefs() drops absent states from a renamed previous-state column", {
  skip_on_cran()
  skip_if_not_installed("rms")
  data <- make_test_data(n_patients = 120, follow_up_time = 6, seed = 9461)
  data$yprev <- factor(data$yprev)
  renamed <- data
  names(renamed)[names(renamed) == "yprev"] <- "prev_state"
  fit_yprev <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    absorb = 6
  )
  fit_renamed <- mostm(
    ordered(y) ~ time + tx + prev_state,
    data = renamed,
    id_var = "id",
    p_var = "prev_state",
    absorb = 6
  )

  # Every resample holds only the patients who were never in state 1.
  ids <- setdiff(unique(data$id), data$id[data$y == 1 | data$yprev == 1])
  local_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      replicate(
        n_boot,
        data.frame(original_id = ids, new_id = paste0(ids, "_1"), boot_id = 1),
        simplify = FALSE
      )
    }
  )
  with_yprev <- bootstrap_model_coefs(fit_yprev, n_boot = 1)
  expect_no_warning(
    with_renamed <- bootstrap_model_coefs(fit_renamed, n_boot = 1)
  )
  expect_equal(
    unname(unlist(with_renamed)),
    unname(unlist(with_yprev)),
    tolerance = 1e-8
  )
})
