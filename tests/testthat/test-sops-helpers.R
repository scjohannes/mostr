test_that("baseline selection keeps factor id order, earliest visits and ties", {
  data <- data.frame(
    id = factor(
      c("a", "z", "a", "z", "a", "m", "m"),
      levels = c("z", "a", "m")
    ),
    visit = c(NA, 2, 1, 2, 1, NA, NA),
    value = 1:7
  )
  before <- serialize(data, NULL)
  expect_identical(
    resolve_markov_prediction_data(data, "id", "visit"),
    data[c(2, 3, 6), ]
  )
  expect_identical(serialize(data, NULL), before)
})

test_that("sops() validates inputs and supports stratified aggregation", {
  model <- structure(list(), class = "mock_model", markov_id_var = "id")
  newdata <- data.frame(
    id = 1:3,
    time = c(1, 2, 1),
    yprev = c(1, 2, 1),
    tx = c(0, 1, 1)
  )
  sops_array <- array(
    c(
      0.1,
      0.2,
      0.3,
      0.4,
      0.5,
      0.6,
      0.9,
      0.8,
      0.7,
      0.6,
      0.5,
      0.4
    ),
    dim = c(3, 2, 2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(
      object,
      data,
      times,
      y_levels,
      absorb,
      time_var,
      p_var,
      gap,
      time_covariates,
      ...
    ) {
      expect_equal(times, c(1, 2))
      expect_equal(y_levels, c(1, 2))
      sops_array
    },
    {
      expect_error(
        sops(model, newdata = newdata, times = 1:2, y_levels = 1:2, by = 1),
        "`by` must be `FALSE`, `TRUE`, or a character vector of column names.",
        fixed = TRUE
      )
      expect_error(
        sops(
          model,
          newdata = newdata,
          times = 1:2,
          y_levels = 1:2,
          by = "missing"
        ),
        "Grouping variables missing: missing",
        fixed = TRUE
      )

      # Averages within groups are `avg_sops()` results.
      result <- sops(
        model,
        newdata = newdata,
        times = 1:2,
        y_levels = 1:2,
        by = "tx"
      )
    }
  )

  expected <- aggregate(
    estimate ~ time + state + tx,
    data = data.frame(
      tx = rep(newdata$tx, times = 4),
      time = rep(rep(1:2, each = 3), times = 2),
      state = rep(1:2, each = 6),
      estimate = as.vector(sops_array)
    ),
    FUN = mean
  )
  result <- result[order(result$time, result$state, result$tx), ]
  expected <- expected[order(expected$time, expected$state, expected$tx), ]
  rownames(result) <- NULL
  rownames(expected) <- NULL

  expect_equal(
    as.data.frame(result[, c("time", "state", "tx", "estimate")]),
    expected[, c("time", "state", "tx", "estimate")]
  )
  expect_s3_class(result, "markov_avg_sops")
  expect_equal(attr(result, "avg_args")$by, "tx")
})

test_that("sops() infers y_levels from supported model containers", {
  newdata <- data.frame(id = 1, time = 1, yprev = 1)
  sops_array <- array(c(0.25, 0.75), dim = c(1, 1, 2))

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(
      object,
      data,
      times,
      y_levels,
      absorb,
      time_var,
      p_var,
      gap,
      time_covariates,
      ...
    ) {
      expect_equal(y_levels, c("1", "2"))
      sops_array
    },
    {
      vglm_model <- methods::new("vglm")
      vglm_model@extra <- list(colnames.y = c("1", "2"))
      vglm_model@family <- methods::new("vglmff")
      vglm_model@call <- quote(VGAM::vglm(
        y ~ x,
        family = VGAM::cumulative(reverse = TRUE)
      ))

      robust_model <- list(extra = list(colnames.y = c("1", "2")))
      class(robust_model) <- "robcov_vglm"

      orm_model <- list(yunique = c("1", "2"))
      class(orm_model) <- "orm"
      # `mostm()` fits record their patient column.
      attr(vglm_model, "markov_id_var") <- "id"
      attr(robust_model, "markov_id_var") <- "id"
      attr(orm_model, "markov_id_var") <- "id"

      expect_equal(nrow(sops(vglm_model, newdata = newdata, times = 1)), 2)
      expect_equal(nrow(sops(robust_model, newdata = newdata, times = 1)), 2)
      expect_equal(nrow(sops(orm_model, newdata = newdata, times = 1)), 2)

      vglm_model@extra <- list(colnames.y = NULL)
      expect_error(
        sops(vglm_model, newdata = newdata, times = 1),
        "`y_levels` cannot be NULL",
        fixed = TRUE
      )

      robust_model$extra <- list(colnames.y = NULL)
      expect_error(
        sops(robust_model, newdata = newdata, times = 1),
        "`y_levels` cannot be NULL",
        fixed = TRUE
      )

      orm_model$yunique <- NULL
      expect_error(
        sops(orm_model, newdata = newdata, times = 1),
        "`y_levels` cannot be NULL",
        fixed = TRUE
      )

      expect_error(
        sops(
          structure(list(), class = "other", markov_id_var = "id"),
          newdata = newdata,
          times = 1
        ),
        "`y_levels` cannot be NULL",
        fixed = TRUE
      )
    }
  )
})

test_that("avg_sops() validates inputs and preserves grouping metadata", {
  model <- structure(list(), class = "mock_model", markov_id_var = "id")
  newdata <- data.frame(
    id = c(1, 2),
    tx = c(0, 1),
    subgroup = c("a", "b"),
    yprev = c(1, 2)
  )
  time_covariates <- data.frame(visit = 1:2, spline = c(0, 1))
  captured <- new.env(parent = emptyenv())

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    sops = function(model, newdata, times, ...) {
      captured$args <- list(...)
      out <- expand.grid(
        time = times,
        state = c(1, 2),
        tx = c(0, 1),
        subgroup = c("a", "b"),
        KEEP.OUT.ATTRS = FALSE
      )
      out$estimate <- seq_len(nrow(out)) / 100
      attr(out, "model") <- model
      attr(out, "call_args") <- list(times = times)
      attr(out, "time_var") <- "time"
      attr(out, "p_var") <- "yprev"
      attr(out, "y_levels") <- c(1, 2)
      attr(out, "absorb") <- NULL
      attr(out, "gap") <- NULL
      attr(out, "time_covariates") <- NULL
      out
    },
    {
      observed <- avg_sops(model, newdata = newdata, times = 1)
      expect_named(observed, c("time", "state", "estimate"))
      expect_null(attr(observed, "avg_args")$variables)
      expect_equal(attr(observed, "newdata_pred")$tx, newdata$tx)
      expect_error(
        avg_sops(model, newdata = NULL, variables = "tx", times = 1),
        "Provide newdata",
        fixed = TRUE
      )
      expect_error(
        avg_sops(model, newdata = newdata, variables = "tx"),
        "`times` must be supplied to `avg_sops()`.",
        fixed = TRUE
      )
      expect_no_error(
        avg_sops(
          model,
          newdata = data.frame(tx = 0),
          variables = "tx",
          times = 1
        )
      )
      expect_error(
        avg_sops(
          model,
          newdata = newdata,
          variables = list(missing = 1),
          times = 1
        ),
        "Variables not in data",
        fixed = TRUE
      )

      result <- avg_sops(
        model,
        newdata = newdata,
        variables = "tx",
        by = "subgroup",
        times = 1:2,
        y_levels = c(1, 2),
        absorb = 2,
        time_var = "visit",
        p_var = "prev",
        p2_var = "preprev",
        gap_var = "gap",
        time_covariates = time_covariates
      )
    }
  )

  expect_true("subgroup" %in% names(result))
  expect_equal(attr(result, "avg_args")$by, "subgroup")
  expect_equal(captured$args$y_levels, c(1, 2))
  expect_equal(captured$args$absorb, 2)
  expect_equal(captured$args$time_var, "visit")
  expect_equal(captured$args$p_var, "prev")
  expect_equal(captured$args$p2_var, "preprev")
  expect_equal(captured$args$gap, "gap")
  expect_equal(captured$args$time_covariates, time_covariates)
})

test_that("avg_sops() reports missing grouping variables from individual SOPs", {
  model <- structure(list(), class = "mock_model", markov_id_var = "id")
  newdata <- data.frame(id = 1, tx = 0, yprev = 1)

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    sops = function(model, newdata, times, ...) {
      out <- data.frame(time = 1, state = 1, tx = 0, estimate = 0.5)
      attr(out, "model") <- model
      attr(out, "call_args") <- list(times = times)
      out
    },
    {
      expect_error(
        avg_sops(
          model,
          newdata = newdata,
          variables = "tx",
          by = "subgroup",
          times = 1
        ),
        "Grouping variables missing",
        fixed = TRUE
      )
    }
  )
})

test_that("get_draws() rejects plain data frames and joins draws on the keys they share", {
  expect_error(
    get_draws(data.frame(x = 1)),
    "get_draws() requires an object from inferences",
    fixed = TRUE
  )

  sop_result <- function(estimates, draws) {
    object <- data.frame(
      time = seq_along(estimates),
      state = 1,
      estimate = estimates
    )
    class(object) <- c("markov_avg_sops", class(object))
    attr(object, "draws") <- draws
    object
  }
  cases <- list(
    list(
      label = "draws without join keys keep only their own columns",
      object = sop_result(
        0.5,
        data.frame(draw_id = 1:2, estimate = c(0.4, 0.6))
      ),
      names = c("draw_id", "draw"),
      draw = c(0.4, 0.6),
      estimate = NULL
    ),
    list(
      label = "draws with time and state get the estimates joined",
      object = sop_result(
        c(0.5, 0.6),
        data.frame(
          draw_id = c(1, 1),
          time = 1:2,
          state = 1,
          estimate = c(0.4, 0.7)
        )
      ),
      names = NULL,
      draw = c(0.4, 0.7),
      estimate = c(0.5, 0.6)
    )
  )

  for (case in cases) {
    draws <- get_draws(case$object)
    if (!is.null(draws$time)) {
      draws <- draws[order(draws$time), ]
      rownames(draws) <- NULL
    }
    if (!is.null(case$names)) {
      expect_equal(names(draws), case$names, info = case$label)
    }
    expect_equal(draws$draw, case$draw, info = case$label)
    if (!is.null(case$estimate)) {
      expect_equal(draws$estimate, case$estimate, info = case$label)
    }
  }
})

test_that("get_draws() preserves covariates for individual SOP draws", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("dplyr")

  withr::local_seed(123)
  data <- make_test_data(n_patients = 20, seed = 123, follow_up_time = 30)
  data$age <- stats::rnorm(nrow(data), 60, 10)

  expect_warning(
    m_vglm_rob <- mostm(
      ordered(y) ~ rms::rcs(time, 3) + tx + yprev + age,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      engine = "vglm"
    ),
    "fewer than 30 clusters"
  )

  sops_result <- sops(
    model = m_vglm_rob,
    newdata = dplyr::filter(data, time == 1),
    times = 1:5,
    y_levels = 1:6,
    absorb = 6,
    time_var = "time",
    p_var = "yprev"
  )

  result_ci <- inferences(
    sops_result,
    method = "mvn",
    n_draws = 10,
    return_draws = TRUE
  )

  draws <- get_draws(result_ci)
  sample_row <- draws[1, ]
  orig_row <- data[sample_row$rowid, ]

  expect_s3_class(draws, "data.frame")
  expect_contains(
    names(draws),
    c("draw_id", "rowid", "time", "state", "estimate")
  )
  expect_contains(names(draws), c("tx", "age"))
  expect_equal(sample_row$tx, orig_row$tx)
  expect_equal(sample_row$age, orig_row$age)
})

test_that("markov_msm_run() propagates probabilities without absorbing states", {
  components <- list(
    X_init = matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x")),
    X_transition = list(
      list(
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x")),
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x"))
      ),
      list(
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x")),
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x"))
      )
    ),
    n_pat = 2L,
    n_states = 2L,
    M = 1L,
    y_levels = 1:2,
    col_names = "x"
  )
  gamma <- matrix(0, nrow = 1, ncol = 1)
  colnames(gamma) <- "x"

  out <- mostr:::markov_msm_run(
    components = components,
    Gamma = gamma,
    times = 1:2,
    absorb = NULL
  )

  expect_equal(dim(out), c(2L, 2L, 2L))
  expect_equal(out[, 1, 1], c(0.5, 0.5))
  expect_equal(out[, 2, 1], c(0.5, 0.5))
})

test_that("time_in_state() sums occupancy over time for bootstrap data frames and arrays", {
  boot <- data.frame(
    boot_id = rep(1:2, each = 4),
    time = rep(1:2, times = 4),
    tx = rep(c(1, 0), each = 2, times = 2),
    state_1 = c(0.8, 0.7, 0.5, 0.4, 0.6, 0.5, 0.3, 0.2),
    state_2 = c(0.2, 0.3, 0.5, 0.6, 0.4, 0.5, 0.7, 0.8)
  )
  arr3 <- array(
    c(0.8, 0.6, 0.2, 0.4, 0.7, 0.5, 0.3, 0.5),
    dim = c(2, 2, 2),
    dimnames = list(NULL, NULL, c("1", "2"))
  )
  arr4 <- array(
    seq_len(2 * 2 * 2 * 2) / 100,
    dim = c(2, 2, 2, 2),
    dimnames = list(NULL, NULL, NULL, c("1", "2"))
  )
  unnamed <- array(rep(c(0.2, 0.8), 4), dim = c(2, 2, 2))

  # `value` picks the part of the result that `expected` gives.
  cases <- list(
    list(
      label = "data frame, treated arm in one state",
      input = boot,
      target_states = 1,
      value = function(out) out$SOP_tx,
      expected = c(1.5, 1.1)
    ),
    list(
      label = "data frame, control arm in one state",
      input = boot,
      target_states = 1,
      value = function(out) out$SOP_ctrl,
      expected = c(0.9, 0.5)
    ),
    list(
      label = "data frame, all states give no difference",
      input = boot,
      target_states = c(1, 2),
      value = function(out) out$delta,
      expected = c(0, 0)
    ),
    list(
      label = "data frame without tx",
      input = data.frame(boot_id = 1, time = 1),
      error = "must contain 'boot_id', 'time', and 'tx'"
    ),
    list(
      label = "data frame without state columns",
      input = data.frame(boot_id = 1, time = 1, tx = 0),
      error = "does not contain any columns starting with 'state_'"
    ),
    list(
      label = "data frame without the target state",
      input = boot,
      target_states = 3,
      error = "target state columns were not found"
    ),
    list(
      label = "3D array with named states",
      input = arr3,
      target_states = "1",
      value = identity,
      expected = c(1.0, 1.0)
    ),
    list(
      label = "4D array of draws",
      input = arr4,
      target_states = c("1", "2"),
      value = dim,
      expected = c(2L, 2L)
    ),
    list(
      label = "3D array with unnamed states",
      input = unnamed,
      target_states = 1,
      value = identity,
      expected = c(0.4, 1.6)
    ),
    list(
      label = "a vector",
      input = 1:3,
      error = "must be an array or data frame"
    ),
    list(
      label = "a matrix",
      input = matrix(1, 2, 2),
      error = "must be 3D or 4D"
    ),
    list(
      label = "3D array without the target state",
      input = arr3,
      target_states = "3",
      error = "Target states not found"
    )
  )

  for (case in cases) {
    args <- list(case$input)
    if (!is.null(case$target_states)) {
      args$target_states <- case$target_states
    }
    if (is.null(case$error)) {
      out <- do.call(time_in_state, args)
      expect_equal(case$value(out), case$expected, info = case$label)
    } else {
      expect_error(
        do.call(time_in_state, args),
        case$error,
        fixed = TRUE,
        info = case$label
      )
    }
  }
})

test_that("low-level SOP helpers validate model families and state metadata", {
  skip_if_not_installed("VGAM")
  loadNamespace("VGAM")

  expect_error(
    mostr:::validate_markov_model(lm(mpg ~ wt, data = mtcars)),
    NA
  )

  model <- methods::new("vglm")
  attr(model, "markov_engine") <- "vglm"
  model@family <- methods::new("vglmff", vfamily = "binomialff")
  expect_error(
    mostr:::validate_markov_model(model),
    "only the cumulative family",
    fixed = TRUE
  )

  # The reverse coding is read from the fitted family's `infos()` report.
  model@family <- methods::new(
    "vglmff",
    vfamily = "cumulative",
    infos = function(...) list()
  )
  model@misc$link <- "logitlink"
  expect_error(
    mostr:::validate_markov_model(model),
    "Cannot determine 'reverse'",
    fixed = TRUE
  )

  model@family <- methods::new(
    "vglmff",
    vfamily = "cumulative",
    infos = function(...) list(reverse = FALSE)
  )
  expect_error(
    mostr:::validate_markov_model(model),
    "must use reverse = TRUE",
    fixed = TRUE
  )

  offset_data <- data.frame(
    y = ordered(rep(1:3, length.out = 30)),
    x = rep(c(0, 1), length.out = 30),
    off = seq(-0.2, 0.2, length.out = 30)
  )
  offset_model <- suppressWarnings(
    VGAM::vglm(
      y ~ x + offset(off),
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = offset_data
    )
  )
  attr(offset_model, "markov_engine") <- "vglm"
  expect_error(
    mostr:::validate_markov_model(offset_model),
    "Model offsets are not supported",
    fixed = TRUE
  )

  null_offset_model <- suppressWarnings(
    VGAM::vglm(
      y ~ x,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = offset_data,
      offset = NULL
    )
  )
  attr(null_offset_model, "markov_engine") <- "vglm"
  expect_error(
    mostr:::validate_markov_model(null_offset_model),
    NA
  )

  fit_with_optional_offset <- function(offset = NULL) {
    fit <- suppressWarnings(
      VGAM::vglm(
        y ~ x,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = offset_data,
        offset = offset
      )
    )
    attr(fit, "markov_engine") <- "vglm"
    fit
  }
  expect_error(
    mostr:::validate_markov_model(fit_with_optional_offset()),
    NA
  )

  expect_error(
    mostr:::score_bootstrap_components(list()),
    "supports robcov_vglm and orm",
    fixed = TRUE
  )

  expect_error(
    mostr:::marginalize_sops_array(
      array(1, dim = c(1, 1, 1)),
      grid = data.frame(tx = 0),
      times = 1,
      y_levels = 1,
      variables = list(tx = 0),
      n_cf = 1,
      n_each = 1,
      weights = 0
    ),
    "positive finite sum",
    fixed = TRUE
  )

  no_rowid <- mostr:::array_to_df_individual(
    array(1, dim = c(1, 1, 1)),
    times = 1,
    y_levels = 1,
    newdata = data.frame(id = 10)
  )
  expect_equal(no_rowid$rowid, 1)
})

test_that("low-level SOP helpers reject orm offsets", {
  skip_if_not_installed("rms")

  offset_data <- data.frame(
    y = ordered(rep(1:3, length.out = 45)),
    x = rep(c(0, 1, 2), length.out = 45),
    off = seq(-0.2, 0.2, length.out = 45)
  )
  offset_model <- rms::orm(
    y ~ x + offset(off),
    data = offset_data,
    x = TRUE,
    y = TRUE
  )
  attr(offset_model, "markov_engine") <- "orm"

  expect_error(
    mostr:::validate_markov_model(offset_model),
    "Model offsets are not supported",
    fixed = TRUE
  )
})
