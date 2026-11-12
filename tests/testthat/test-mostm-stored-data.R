test_that("each engine stores its fitting data and the patient ID column", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1001)
  # `last_time` is the last visit kept by `subset`, or NULL for all visits.
  # `check` holds the checks of one engine.
  engines <- list(
    list(
      label = "vglm",
      args = list(
        ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        engine = "vglm"
      ),
      last_time = NULL,
      check = function(fit) {
        expect_s3_class(fit, "robcov_vglm")
        expect_s4_class(fit$vglm_fit, "vglm")
        expect_identical(fit$type, "HC0")
        expect_identical(fit$cadjust, TRUE)
        expect_equal(attr(fit$vglm_fit, "markov_id_var"), "id")
        metadata <- attr(fit, "markov_starting_profile_metadata")
        expect_equal(metadata$first_followup_time, 1)
        expect_equal(
          attr(fit$vglm_fit, "markov_starting_profile_metadata"),
          metadata
        )
      }
    ),
    list(
      label = "orm",
      args = list(ordered(y) ~ time + tx + yprev),
      last_time = NULL,
      check = function(fit) {
        expect_s3_class(fit, "orm")
        expect_false(is.null(fit$orig.var))
      }
    ),
    list(
      label = "orm with `subset`",
      args = list(ordered(y) ~ time + tx + yprev),
      last_time = 4,
      check = function(fit) expect_s3_class(fit, "orm")
    )
  )
  if (requireNamespace("rmsb", quietly = TRUE)) {
    # `standata = TRUE` prepares the data for Stan without sampling.
    engines <- c(engines, list(list(
      label = "blrm with `subset`",
      args = list(ordered(y) ~ time + tx + yprev, standata = TRUE, engine = "blrm"),
      last_time = 3,
      check = function(fit) expect_type(fit, "list")
    )))
  }

  for (case in engines) {
    args <- c(case$args, list(data = data, id_var = "id"))
    if (!is.null(case$last_time)) {
      args$subset <- call("<=", quote(time), case$last_time)
    }
    fit <- expect_no_error(do.call(mostm, args))

    case$check(fit)
    expect_equal(attr(fit, "markov_id_var"), "id", info = case$label)
    stored <- attr(fit, "markov_data")
    if (is.null(case$last_time)) {
      expect_equal(nrow(stored), nrow(data), info = case$label)
      expect_equal(nrow(attr(fit, "markov_refit_data")), nrow(data), info = case$label)
      expect_equal(
        nrow(attr(fit, "markov_starting_profile_data")),
        length(unique(data$id)),
        info = case$label
      )
    } else {
      expect_equal(nrow(stored), sum(data$time <= case$last_time), info = case$label)
      expect_true(all(stored$time <= case$last_time), info = case$label)
    }
  }
})

test_that("the vglm engine forwards robust covariance corrections", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1019)
  fit <- mostm(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    type = "HC1",
    cadjust = FALSE,
    engine = "vglm"
  )

  expect_identical(fit$type, "HC1")
  expect_identical(fit$cadjust, FALSE)
  direct <- robcov_vglm(
    fit$vglm_fit,
    cluster = data$id,
    type = "HC1",
    cadjust = FALSE
  )
  expect_equal(fit$var, direct$var)
})

test_that("mostm() validates numeric first follow-up schedules", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  local_reproducible_output(width = 80)

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1016)
  no_one <- data
  no_one$time <- no_one$time + 1
  expect_snapshot(
    mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = no_one,
      id_var = "id",
      engine = "vglm"
    ),
    error = TRUE
  )

  below_one <- data
  below_one$time[[1L]] <- 0
  expect_snapshot(
    mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = below_one,
      id_var = "id",
      engine = "vglm"
    ),
    error = TRUE
  )

  alternative <- data
  alternative$time <- alternative$time + 1
  fit <- mostm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = alternative,
    id_var = "id",
    first_followup_time = 2,
    engine = "vglm"
  )
  expect_equal(
    markov_model_starting_profile_metadata(fit)$first_followup_time,
    2
  )

  expect_snapshot(
    mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      first_followup_time = 2,
      engine = "vglm"
    ),
    error = TRUE
  )
})

test_that("mostm() rejects former starting-profile arguments", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  local_reproducible_output(width = 80)

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1017)
  expect_snapshot(
    mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      start_time = 1,
      origin_time = 0,
      time_map = 1:6,
      engine = "vglm"
    ),
    error = TRUE
  )
})

test_that("the vglm engine aligns stored fitting data after subset and NA drops", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1007)
  data$time_lin[2] <- NA_real_

  fit <- mostm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    subset = time <= 5,
    engine = "vglm"
  )

  stored_data <- attr(fit, "markov_data")
  expected_data <- data[data$time <= 5 & !is.na(data$time_lin), , drop = FALSE]

  expect_s3_class(fit, "robcov_vglm")
  expect_equal(
    stored_data[c("id", "time", "time_lin")],
    expected_data[c("id", "time", "time_lin")]
  )
  expect_equal(
    attr(fit$vglm_fit, "markov_data")[c("id", "time")],
    expected_data[c("id", "time")]
  )
  expect_equal(nrow(fit$vglm_fit@x), nrow(stored_data))
  expect_equal(
    nrow(attr(fit, "markov_refit_data")),
    sum(data$time <= 5)
  )
})

test_that("the orm engine stores fitting data and applies robust covariance", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1002)
  fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id"
  )

  expect_s3_class(fit, "orm")
  expect_equal(attr(fit, "markov_id_var"), "id")
  expect_equal(nrow(attr(fit, "markov_data")), nrow(data))
  expect_equal(nrow(attr(fit, "markov_refit_data")), nrow(data))
  expect_equal(
    nrow(attr(fit, "markov_starting_profile_data")),
    length(unique(data$id))
  )
  expect_false(is.null(fit$orig.var))

  subset_fit <- expect_no_error(
    mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      subset = time <= 4
    )
  )
  expect_equal(nrow(attr(subset_fit, "markov_data")), sum(data$time <= 4))
  expect_true(all(attr(subset_fit, "markov_data")$time <= 4))
})

test_that("the blrm engine stores fitting data without sampling when requested", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("rmsb")

  data <- make_test_data(n_patients = 12, follow_up_time = 6, seed = 1005)
  fit <- mostm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    subset = time <= 3,
    standata = TRUE,
    engine = "blrm"
  )

  expect_type(fit, "list")
  expect_equal(attr(fit, "markov_id_var"), "id")
  expect_equal(nrow(attr(fit, "markov_data")), sum(data$time <= 3))
  expect_true(all(attr(fit, "markov_data")$time <= 3))

  expect_error(
    mostm(
      ordered(y) ~ offset(time) + tx + yprev,
      data = data,
      id_var = "id",
      standata = TRUE,
      engine = "blrm"
    ),
    "offset"
  )
})

test_that("automatic SOP prediction never substitutes a later profile row", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1003)
  data <- data[!(data$id == 1 & data$time == min(data$time)), , drop = FALSE]
  fit <- mostm(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    engine = "vglm"
  )

  local_reproducible_output(width = 80)
  expect_snapshot(
    sops(
      fit,
      times = 1:3,
      absorb = 6
    ),
    error = TRUE
  )
})

test_that("missing values leave complete starting profiles of the fitted patients", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  base <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1010)
  first_time <- min(base$time)
  vglm_args <- list(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    engine = "vglm"
  )
  # Each case changes the data, fits it, and checks the starting profiles.
  cases <- list(
    list(
      label = "missing first outcome",
      change = function(data) {
        data$y[data$id == 1 & data$time == first_time] <- NA
        data
      },
      args = vglm_args,
      check = function(fit, data, profiles) {
        # Patient 1 keeps a complete starting profile, with a missing `y`.
        fitted_data <- attr(fit, "markov_data")
        expect_equal(
          any(fitted_data$id == 1 & fitted_data$time == first_time),
          FALSE
        )
        expect_equal(1 %in% fitted_data$id, TRUE)
        expect_equal(profiles$time[profiles$id == 1], first_time)
        expect_equal(is.na(profiles$y[profiles$id == 1]), TRUE)
      }
    ),
    list(
      label = "missing time at a later visit",
      change = function(data) {
        data$time[which(data$time == 3)[1L]] <- NA_real_
        data
      },
      args = list(y ~ time + tx + yprev),
      check = function(fit, data, profiles) {
        starting_rows <- which(data$time == 1)
        expect_equal(profiles$id, data$id[starting_rows])
        expect_equal(profiles$.markov_source_row, starting_rows)
        expect_equal(profiles$time, rep(1, length(starting_rows)))
      }
    ),
    list(
      label = "a patient with only missing outcomes",
      change = function(data) {
        data$y[data$id == 1] <- NA
        data
      },
      args = vglm_args,
      check = function(fit, data, profiles) {
        # Such a patient is not part of the fitted cohort.
        fitted_ids <- unique(attr(fit, "markov_data")$id)
        expect_equal(1 %in% profiles$id, FALSE)
        expect_setequal(profiles$id, fitted_ids)
        expect_equal(nrow(profiles), length(unique(data$id)) - 1L)
      }
    )
  )

  sop_args <- list(times = 1:2, absorb = 6)
  for (case in cases) {
    data <- case$change(base)
    fit <- do.call(mostm, c(case$args, list(data = data, id_var = "id")))
    profiles <- markov_validate_starting_profiles(fit)
    case$check(fit, data, profiles)

    # The SOP functions predict from exactly these profiles.
    individual <- do.call(sops, c(list(fit), sop_args))
    average <- do.call(
      avg_sops,
      c(list(fit, variables = list(tx = c(0, 1))), sop_args)
    )
    comparison <- do.call(
      avg_comparisons,
      c(
        list(
          fit,
          variables = list(tx = c(0, 1)),
          estimand = "sop",
          state_sets = list(low = 1:2)
        ),
        sop_args
      )
    )
    expect_equal(attr(individual, "newdata_orig")$id, profiles$id, info = case$label)
    expect_equal(attr(average, "newdata_orig")$id, profiles$id, info = case$label)
    expect_equal(attr(comparison, "newdata_orig")$id, profiles$id, info = case$label)
  }
})

test_that("incomplete and duplicated starting profiles fail clearly", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  local_reproducible_output(width = 80)

  incomplete <- make_test_data(
    n_patients = 35,
    follow_up_time = 6,
    seed = 1012
  )
  incomplete$tx[incomplete$id == 1 & incomplete$time == min(incomplete$time)] <-
    NA_real_
  incomplete_fit <- mostm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = incomplete,
    id_var = "id",
    engine = "vglm"
  )
  expect_snapshot(
    avg_sops(
      incomplete_fit,
      variables = list(tx = c(0, 1)),
      times = 1:2,
      absorb = 6
    ),
    error = TRUE
  )

  duplicated <- make_test_data(
    n_patients = 35,
    follow_up_time = 6,
    seed = 1013
  )
  duplicate_row <- which(
    duplicated$id == 1 & duplicated$time == min(duplicated$time)
  )[[1L]]
  duplicated <- rbind(duplicated, duplicated[duplicate_row, , drop = FALSE])
  # mostm() warns about the duplicated row, and avg_sops() refuses it.
  expect_warning(
    duplicated_fit <- mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = duplicated,
      id_var = "id",
      engine = "vglm"
    ),
    "duplicate `id` x `time` combinations"
  )
  expect_snapshot(
    avg_sops(
      duplicated_fit,
      variables = list(tx = c(0, 1)),
      times = 1:2,
      absorb = 6
    ),
    error = TRUE
  )
})

test_that("factor and character time require explicit first follow-up values", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  local_reproducible_output(width = 80)

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1015)
  complete_ids <- as.integer(names(which(table(data$id) == 6L)))
  data <- data[data$id %in% complete_ids, , drop = FALSE]
  visits <- paste0("v", seq_len(6L))
  data$time <- factor(paste0("v", data$time), levels = visits)
  expect_snapshot(
    mostm(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      engine = "vglm"
    ),
    error = TRUE
  )

  fit_factor <- mostm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    first_followup_time = "v1",
    engine = "vglm"
  )
  expect_equal(
    unique(as.character(markov_validate_starting_profiles(fit_factor)$time)),
    "v1"
  )

  data$visit <- as.character(data$time)
  expect_snapshot(
    mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      time_var = "visit",
      engine = "vglm"
    ),
    error = TRUE
  )
  fit_character <- mostm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    time_var = "visit",
    first_followup_time = "v1",
    engine = "vglm"
  )
  expect_equal(
    unique(markov_validate_starting_profiles(fit_character)$visit),
    "v1"
  )

  expect_snapshot(
    mostm(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      first_followup_time = "v7",
      engine = "vglm"
    ),
    error = TRUE
  )
})

test_that("SOP functions keep the stored refit data and number supplied newdata rows", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1014)
  fit <- mostm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    engine = "vglm"
  )
  # Rows of `newdata` without a patient column: one patient's starting row
  # moved to other times, with a `rowid` that the functions replace.
  profiles_without_id <- function(n_rows, rowid) {
    profiles <- data[data$time == 1, , drop = FALSE][seq_len(n_rows), ]
    profiles$time <- seq_len(n_rows)
    profiles$rowid <- rowid
    profiles[, setdiff(names(profiles), "id"), drop = FALSE]
  }
  sop_args <- list(times = 1:2, absorb = 6)

  # `rowid` gives the expected row numbers of the prediction data, and
  # `warning` the warning of a fractional bootstrap, which uses the stored
  # refit data rather than the prediction rows.
  cases <- list(
    list(
      label = "avg_sops() on the stored starting profiles",
      run = function() {
        do.call(avg_sops, c(list(fit, variables = list(tx = c(0, 1))), sop_args))
      },
      rowid = NULL,
      warning = NULL
    ),
    list(
      label = "sops() on newdata, fixed prediction profiles",
      run = function() {
        do.call(sops, c(list(fit, newdata = profiles_without_id(3, 99)), sop_args))
      },
      rowid = seq_len(3),
      warning = "fixed prediction profiles",
      seed = 10081,
      each_row = TRUE
    ),
    list(
      label = "avg_sops() on newdata, standardization profiles",
      run = function() {
        do.call(
          avg_sops,
          c(
            list(
              fit,
              newdata = profiles_without_id(4, 100 + seq_len(4)),
              variables = list(tx = c(0, 1))
            ),
            sop_args
          )
        )
      },
      # Each profile once per treatment value.
      rowid = rep(seq_len(4), times = 2L),
      warning = "baseline weights have been set to NULL",
      seed = 10091
    )
  )

  for (case in cases) {
    out <- case$run()

    # Refit bootstraps resample the stored refit data, not the prediction rows.
    expect_equal(
      attr(out, "refit_data"),
      mostr:::markov_model_refit_data(fit),
      info = case$label
    )
    expect_equal(nrow(attr(out, "refit_data")), nrow(data), info = case$label)
    expect_equal(
      mostr:::markov_result_settings(out)$id_var,
      "id",
      info = case$label
    )
    if (is.null(case$rowid)) {
      next
    }

    prediction_data <- attr(out, "newdata_pred")
    expect_equal(nrow(prediction_data), length(case$rowid), info = case$label)
    expect_equal(prediction_data$rowid, case$rowid, info = case$label)
    each_row <- isTRUE(case$each_row)
    if (each_row) {
      # SOPs for each profile, and their draws, keep the new row numbers.
      expect_equal(sort(unique(out$rowid)), sort(unique(case$rowid)), info = case$label)
    }

    withr::with_seed(case$seed, {
      expect_warning(
        inferred <- inferences(
          out,
          method = "fwb",
          n_draws = 1,
          return_draws = TRUE,
          use_coefstart = FALSE
        ),
        case$warning,
        fixed = TRUE
      )
    })
    expect_equal(attr(inferred, "n_successful"), 1L, info = case$label)
    if (each_row) {
      draws <- get_draws(inferred)
      expect_equal(sort(unique(draws$rowid)), sort(unique(case$rowid)), info = case$label)
    }
  }
})

test_that("FWB inference on SOPs for each profile works with every frequentist engine", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 1006)
  ids <- unique(data$id)
  sex_by_id <- stats::setNames(
    rep(c("female", "male"), length.out = length(ids)),
    ids
  )
  data$sex <- factor(sex_by_id[as.character(data$id)])

  fits <- list(
    vglm = mostm(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      engine = "vglm"
    ),
    orm = mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id")
  )
  # The ordinary bootstrap of SOPs for each profile is row I42 of
  # test-refusals.R.
  sop_args <- list(
    vglm = list(times = 1:2, absorb = 6),
    orm = list(times = 1:2, absorb = "6")
  )
  fwb_args <- list(
    vglm = list(use_coefstart = FALSE),
    orm = list()
  )

  for (engine in names(fits)) {
    fit <- fits[[engine]]
    individual <- do.call(sops, c(list(fit), sop_args[[engine]]))
    withr::with_seed(10041, {
      individual_ci <- do.call(
        inferences,
        c(
          list(individual, method = "fwb", n_draws = 2, return_draws = TRUE),
          fwb_args[[engine]]
        )
      )
    })
    draws <- get_draws(individual_ci)

    expect_s3_class(individual_ci, "markov_sops")
    expect_equal(attr(individual_ci, "method"), "fwb", info = engine)
    expect_equal(attr(individual_ci, "n_successful"), 2L, info = engine)
    expect_true(
      all(c("conf.low", "conf.high", "std.error") %in% names(individual_ci)),
      info = engine
    )
    expect_false(is.null(attr(individual_ci, "draws")), info = engine)
    expect_equal(length(unique(draws$rowid)), length(unique(data$id)), info = engine)
    expect_contains(names(draws), c("draw_id", "rowid", "sex", "draw"))
    expect_equal(anyNA(draws$draw), FALSE, info = engine)
  }

  # SOPs averaged within groups keep the group column in their draws.
  grouped <- sops(
    fits$orm,
    by = "sex",
    times = 1:2,
    absorb = "6"
  )
  withr::local_seed(10061)
  grouped_ci <- inferences(
    grouped,
    method = "fwb",
    n_draws = 2,
    return_draws = TRUE
  )

  expect_equal(attr(grouped_ci, "n_successful"), 2L)
  expect_setequal(grouped_ci$sex, levels(data$sex))
  expect_equal(anyNA(grouped_ci$conf.low), FALSE)
  expect_contains(
    names(attr(grouped_ci, "draws")),
    c("sex", "draw_id")
  )
})

test_that("sops() on the stored patients returns no internal bookkeeping columns", {
  cell <- grid_cell("orm")
  result <- do.call(sops, c(list(grid_fit(cell)), grid_public_args(cell)))
  expect_false(any(startsWith(names(result), ".markov_")))
  inferred <- inferences(result, method = "delta", vcov = "conditional")
  expect_false(any(startsWith(names(inferred), ".markov_")))
})
