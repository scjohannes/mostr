test_that("fast_group_bootstrap samples group ids and creates unique ids", {
  data <- data.frame(
    id = c("a", "a", "b", "c", "c"),
    value = 1:5
  )

  withr::local_seed(42)
  boot <- fast_group_bootstrap(data, id_var = "id", n_boot = 3)

  expect_length(boot, 3)
  expect_equal(vapply(boot, nrow, integer(1)), rep(3L, 3))
  expect_named(boot[[1]], c("original_id", "new_id", "boot_id"))
  expect_all_true(boot[[1]]$original_id %in% unique(data$id))
  expect_equal(boot[[2]]$boot_id, rep(2, 3))
  expect_equal(length(unique(boot[[1]]$new_id)), nrow(boot[[1]]))
})

test_that("fast_group_bootstrap numbers repeated factor ids like character ids", {
  ids <- sprintf("p%02d", 1:10)
  character_data <- data.frame(id = rep(ids, each = 2), value = 1:20)
  factor_data <- character_data
  factor_data$id <- factor(factor_data$id)

  withr::local_seed(7)
  expected <- fast_group_bootstrap(character_data, id_var = "id", n_boot = 4)
  withr::local_seed(7)
  expect_no_warning(
    out <- fast_group_bootstrap(factor_data, id_var = "id", n_boot = 4)
  )
  for (i in seq_along(expected)) {
    expect_identical(as.character(out[[i]]$original_id), expected[[i]]$original_id)
    # A patient drawn twice is two bootstrap patients: p03_1 and p03_2.
    expect_identical(as.character(out[[i]]$new_id), expected[[i]]$new_id)
    expect_identical(anyDuplicated(out[[i]]$new_id), 0L)
  }
})

test_that("fast_group_bootstrap validates the grouping column", {
  data <- data.frame(id = 1:3)

  expect_error(
    fast_group_bootstrap(data, id_var = "patient", n_boot = 1),
    "id_var 'patient' not found in data",
    fixed = TRUE
  )
})

test_that("materialize_bootstrap_sample_indexed keeps sampled copies, row data, factor levels and unmatched ids", {
  cases <- list(
    list(
      label = "repeated ids copy every row in order",
      data = data.frame(
        id = c(1, 1, 2, 2),
        time = c(1, 2, 1, 2),
        y = c(0, 1, 1, 0)
      ),
      ids = data.frame(
        original_id = c(2, 1, 2),
        new_id = c("2_1", "1_1", "2_2"),
        boot_id = 1
      ),
      expected = list(
        new_id = rep(c("2_1", "1_1", "2_2"), each = 2),
        time = c(1, 2, 1, 2, 1, 2),
        y = c(1, 0, 0, 1, 1, 0)
      )
    ),
    list(
      label = "factor levels and an id not in the data",
      data = data.frame(
        id = c(1, 2, 1),
        y = ordered(c("low", "high", "high"), levels = c("low", "high"))
      ),
      ids = data.frame(
        original_id = c("1", "missing", "1"),
        new_id = c("a", "b", "c"),
        boot_id = 1L
      ),
      expected = list(
        id = c("1", "1", "missing", "1", "1"),
        y = ordered(c("low", "high", NA, "low", "high"), levels = c("low", "high"))
      )
    )
  )

  for (case in cases) {
    before <- serialize(list(case$data, case$ids), NULL)
    out <- materialize_bootstrap_sample_indexed(
      case$ids,
      case$data,
      "id",
      bootstrap_row_plan(case$data, "id")
    )

    expect_identical(class(out), "data.frame", info = case$label)
    expect_equal(nrow(out), length(case$expected[[1]]), info = case$label)
    for (column in names(case$expected)) {
      expect_identical(out[[column]], case$expected[[column]], info = case$label)
    }
    # The inputs are not modified.
    expect_identical(
      serialize(list(case$data, case$ids), NULL),
      before,
      info = case$label
    )
  }
})

test_that("fractional bootstrap weights are mean-one cluster weights", {
  data <- data.frame(
    id = c("a", "a", "b", "c", "c"),
    value = 1:5
  )

  withr::local_seed(123)
  boot_a <- mostr:::generate_fwb_bootstrap_weights(
    data,
    id_var = "id",
    n_boot = 2
  )
  withr::local_seed(123)
  boot_b <- mostr:::generate_fwb_bootstrap_weights(
    data,
    id_var = "id",
    n_boot = 2
  )

  expect_equal(boot_a, boot_b)
  expect_length(boot_a, 2)
  expect_named(boot_a[[1]], c("original_id", "fwb_weight", "boot_id"))
  expect_equal(boot_a[[1]]$original_id, c("a", "b", "c"))
  expect_true(all(boot_a[[1]]$fwb_weight > 0))
  expect_equal(mean(boot_a[[1]]$fwb_weight), 1, tolerance = 1e-12)
  expect_equal(boot_a[[2]]$boot_id, rep(2, 3))
})

test_that("fractional bootstrap materialization expands cluster weights to rows", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(0, 1, 1, 0)
  )
  fwb_weights <- data.frame(
    original_id = c(1, 2),
    fwb_weight = c(0.25, 1.75),
    boot_id = 1
  )

  result <- mostr:::materialize_fwb_bootstrap_sample(
    fwb_weights,
    data,
    "id"
  )

  expect_equal(nrow(result), nrow(data))
  expect_equal(result$fwb_weight, c(0.25, 0.25, 1.75, 1.75))
  expect_equal(result$boot_id, rep(1, 4))
  expect_equal(
    mostr:::fwb_baseline_weights(
      fwb_weights,
      data[!duplicated(data$id), ],
      "id"
    ),
    c(0.25, 1.75)
  )
})

test_that("apply_to_bootstrap and apply_to_fwb_bootstrap materialize samples before sequential analysis", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = c(1, 2, 10, 20)
  )
  cases <- list(
    list(
      label = "ordinary bootstrap",
      run = function() {
        apply_to_bootstrap(
          boot_samples = list(
            data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1),
            data.frame(original_id = c(2, 2), new_id = c("2_1", "2_2"), boot_id = 2)
          ),
          analysis_fn = function(boot_data) {
            c(rows = nrow(boot_data), total = sum(boot_data$y))
          },
          data = data,
          id_var = "id",
          workers = 1
        )
      },
      expected = list(c(rows = 4, total = 33), c(rows = 4, total = 60))
    ),
    list(
      label = "fractional weighted bootstrap",
      run = function() {
        mostr:::apply_to_fwb_bootstrap(
          fwb_samples = list(
            data.frame(original_id = c(1, 2), fwb_weight = c(0.5, 1.5), boot_id = 1),
            data.frame(original_id = c(1, 2), fwb_weight = c(1.2, 0.8), boot_id = 2)
          ),
          analysis_fn = function(boot_data, fwb_weights) {
            c(
              rows = nrow(boot_data),
              total = sum(boot_data$y * boot_data$fwb_weight),
              n_weights = nrow(fwb_weights)
            )
          },
          data = data,
          id_var = "id",
          workers = 1
        )
      },
      expected = list(
        c(rows = 4, total = 46.5, n_weights = 2),
        c(rows = 4, total = 27.6, n_weights = 2)
      )
    )
  )

  for (case in cases) {
    result <- case$run()
    expect_length(result, length(case$expected))
    for (i in seq_along(case$expected)) {
      expect_equal(result[[i]], case$expected[[i]], info = case$label)
    }
  }
})

test_that("parallel bootstrap helpers return results and restore the future plan after success and errors", {
  skip_on_cran()
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))

  one_sample <- list(
    data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1)
  )
  cases <- list(
    list(
      label = "fractional weighted bootstrap",
      run = function() {
        mostr:::apply_to_fwb_bootstrap(
          fwb_samples = list(
            data.frame(original_id = c(1, 2), fwb_weight = c(0.5, 1.5), boot_id = 1)
          ),
          analysis_fn = function(boot_data, fwb_weights) {
            c(
              total = sum(boot_data$y * boot_data$fwb_weight),
              n_weights = nrow(fwb_weights)
            )
          },
          data = data.frame(id = c(1, 1, 2, 2), y = c(1, 2, 10, 20)),
          id_var = "id",
          workers = 2,
          packages = character(0)
        )
      },
      expected = c(total = 46.5, n_weights = 2)
    ),
    list(
      label = "ordinary bootstrap",
      run = function() {
        apply_to_bootstrap(
          boot_samples = one_sample,
          analysis_fn = function(boot_data) sum(boot_data$y),
          data = data.frame(id = c(1, 2), y = c(3, 4)),
          id_var = "id",
          workers = 2,
          packages = character(0)
        )
      },
      expected = 7
    ),
    list(
      label = "ordinary bootstrap with an analysis error",
      run = function() {
        apply_to_bootstrap(
          boot_samples = one_sample,
          analysis_fn = function(boot_data) stop("boom"),
          data = data.frame(id = c(1, 2), y = c(3, 4)),
          id_var = "id",
          workers = 2,
          packages = character(0)
        )
      },
      error = "boom"
    )
  )

  for (case in cases) {
    if (is.null(case$error)) {
      result <- suppressWarnings(case$run())
      expect_equal(result[[1]], case$expected, info = case$label)
    } else {
      expect_error(suppressWarnings(case$run()), case$error, info = case$label)
    }
    expect_equal(future::nbrOfWorkers(), 1L, info = case$label)
  }
})

test_that("parallel helpers restore every level of a nested future plan after an error", {
  skip_on_cran()
  # `future::plan()` without arguments returns only the outer level, so
  # restoring from it dropped the inner levels of the user's plan. Two
  # sequential levels show this without starting extra workers.
  user_plan <- list(future::sequential, future::sequential)
  withr::defer(future::plan(future::sequential))

  failing <- list(
    bootstrap = function() {
      apply_to_bootstrap(
        boot_samples = list(
          data.frame(original_id = 1, new_id = "1_1", boot_id = 1)
        ),
        analysis_fn = function(boot_data) stop("boom"),
        data = data.frame(id = 1, y = 3),
        id_var = "id",
        workers = 2,
        packages = character(0)
      )
    },
    fwb = function() {
      mostr:::apply_to_fwb_bootstrap(
        fwb_samples = list(
          data.frame(original_id = 1, fwb_weight = 1, boot_id = 1)
        ),
        analysis_fn = function(boot_data, fwb_weights) stop("boom"),
        data = data.frame(id = 1, y = 3),
        id_var = "id",
        workers = 2,
        packages = character(0)
      )
    },
    simulation = function() {
      mostr:::apply_sop_simulation_draws(
        n_draws = 1,
        analysis_fn = function(i) stop("boom"),
        workers = 2,
        globals = TRUE,
        packages = character(0)
      )
    }
  )

  for (label in names(failing)) {
    future::plan(user_plan)
    expect_error(suppressWarnings(failing[[label]]()), info = label)
    expect_length(future::plan("list"), 2L)
  }
})

test_that("apply_to_bootstrap maps deprecated parallel argument to workers", {
  data <- data.frame(id = c(1, 2), y = c(3, 4))
  boot_samples <- list(
    data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1)
  )

  expect_warning(
    result <- apply_to_bootstrap(
      boot_samples = boot_samples,
      analysis_fn = function(boot_data) sum(boot_data$y),
      data = data,
      id_var = "id",
      parallel = FALSE
    ),
    "parallel.*deprecated"
  )
  expect_equal(result[[1]], 7)
})

test_that("bootstrap_analysis_wrapper drops absent states and refits models", {
  original_data <- data.frame(
    y = factor(c(1, 2, 3)),
    x = c(0, 1, 2)
  )
  boot_data <- data.frame(
    y = factor(c(1, 3, 3), levels = c(1, 2, 3)),
    x = c(0, 2, 3)
  )
  model <- lm(as.numeric(y) ~ x, data = original_data)

  result <- bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = model,
    factor_cols = "y",
    original_data = original_data,
    y_levels = 1:3,
    absorb = 3
  )

  expect_s3_class(result$model, "lm")
  # The states keep their labels; the absent state 2 is no longer a level.
  expect_equal(levels(result$data$y), c("1", "3"))
  expect_equal(result$ylevels, c("1", "3"))
  expect_equal(result$absorb, "3")
  expect_equal(result$missing_states, "2")
})

test_that("bootstrap_analysis_wrapper reports model fitting failures", {
  boot_data <- data.frame(y = factor(c(1, 2)), x = c(1, 2))
  model <- structure(list(), class = "not_updateable")

  expect_warning(
    result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = "y",
      original_data = boot_data
    ),
    "Bootstrap model fitting failed"
  )

  expect_null(result$model)
  expect_equal(result$missing_states, character(0))
})

test_that("bootstrap_analysis_wrapper refits with the datadist option unset", {
  skip_if_not_installed("rms")

  # The user's `datadist` option as the refit sees it.
  datadist_seen <- "not refitted"
  update.orm <- function(object, data, ...) {
    datadist_seen <<- getOption("datadist")
    structure(list(data = data), class = "orm")
  }
  assign("update.orm", update.orm, envir = globalenv())
  withr::defer(rm("update.orm", envir = globalenv()))

  old_dd <- structure(list(source = "existing"), class = "datadist")
  assign("dd", old_dd, envir = globalenv())
  options(datadist = "existing_dd")
  withr::defer(options(datadist = NULL))
  withr::defer(rm("dd", envir = globalenv()))

  boot_data <- data.frame(
    y = factor(c(1, 2, 3)),
    x = c(0, 1, 2),
    boot_id = 1,
    fwb_weight = c(0.5, 1, 1.5),
    .mostr_fit_weight = c(0.5, 1, 1.5)
  )
  result <- expect_no_warning(
    bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = structure(list(), class = "orm"),
      factor_cols = "y",
      original_data = boot_data
    )
  )

  expect_s3_class(result$model, "orm")
  expect_null(datadist_seen)
  expect_equal(getOption("datadist"), "existing_dd")
  expect_identical(get("dd", envir = globalenv()), old_dd)
})

test_that("update_bootstrap_model suppresses only known orm weight warning", {
  update.orm <- function(object, data, weights = NULL, ...) {
    warning(
      "currently weights are ignored in model validation and bootstrapping orm fits",
      call. = FALSE
    )
    structure(list(data = data), class = "orm")
  }
  assign("update.orm", update.orm, envir = globalenv())
  withr::defer(rm("update.orm", envir = globalenv()))

  boot_data <- data.frame(y = 1:3, .mostr_fit_weight = c(0.5, 1, 1.5))

  result <- expect_no_warning(
    mostr:::update_bootstrap_model(
      structure(list(), class = "orm"),
      boot_data,
      fit_weights = boot_data$.mostr_fit_weight
    )
  )
  expect_s3_class(result, "orm")

  update.orm <- function(object, data, weights = NULL, ...) {
    warning("important warning", call. = FALSE)
    structure(list(data = data), class = "orm")
  }
  assign("update.orm", update.orm, envir = globalenv())

  expect_warning(
    mostr:::update_bootstrap_model(
      structure(list(), class = "orm"),
      boot_data,
      fit_weights = boot_data$.mostr_fit_weight
    ),
    "important warning"
  )
})

test_that("bootstrap_analysis_wrapper falls back when coefstart update fails", {
  coef.vglm <- function(object, ...) object$coefficients
  update.vglm <- function(object, data, coefstart = NULL, ...) {
    if (!is.null(coefstart)) {
      stop("coefstart failed")
    }
    structure(
      list(data = data, coefficients = object$coefficients),
      class = "vglm"
    )
  }
  assign("coef.vglm", coef.vglm, envir = globalenv())
  assign("update.vglm", update.vglm, envir = globalenv())
  withr::defer(rm("coef.vglm", envir = globalenv()))
  withr::defer(rm("update.vglm", envir = globalenv()))

  boot_data <- data.frame(y = factor(c(1, 2)), x = c(0, 1))
  model <- structure(list(coefficients = c(a = 1)), class = "vglm")

  result <- bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = model,
    factor_cols = "y",
    original_data = boot_data,
    use_coefstart = TRUE
  )

  expect_s3_class(result$model, "vglm")
  expect_equal(result$model$data, boot_data)
})

test_that("bootstrap_analysis_wrapper multiplies fractional weights with call weights and stored fit weights", {
  update.mock_weighted <- function(object, data, weights = NULL, ...) {
    weight_expr <- substitute(weights)
    structure(
      list(
        data = data,
        weights = eval(weight_expr, data, parent.frame())
      ),
      class = "mock_weighted"
    )
  }
  assign("update.mock_weighted", update.mock_weighted, envir = globalenv())
  withr::defer(rm("update.mock_weighted", envir = globalenv()))

  call_data <- data.frame(
    y = factor(c(1, 2, 3)),
    x = c(0, 1, 2),
    w = c(2, 3, 4)
  )
  make_weighted_lm <- function() {
    original_data <- data.frame(
      y = c(1, 2, 3, 4),
      x = c(0, 1, 2, 3)
    )
    original_weights <- c(2, 3, 4, 5)

    list(
      model = stats::lm(
        y ~ x,
        data = original_data,
        weights = original_weights
      ),
      data = original_data
    )
  }
  weighted_lm <- make_weighted_lm()

  cases <- list(
    list(
      label = "weights named in the model call",
      model = structure(
        list(call = quote(mock_fit(y ~ x, data = original_data, weights = w))),
        class = "mock_weighted"
      ),
      data = call_data,
      factor_cols = "y",
      fit_weights = c(0.5, 1, 1.5),
      class = "mock_weighted",
      expected = c(1, 3, 6),
      # The refit also sees the product as a data column.
      stored = function(model) model$data$.mostr_fit_weight
    ),
    list(
      label = "weights stored on an lm fit",
      model = weighted_lm$model,
      data = weighted_lm$data,
      factor_cols = character(0),
      fit_weights = c(0.5, 1, 1.5, 2),
      class = "lm",
      expected = c(1, 3, 6, 10),
      stored = function(model) model$model[["(weights)"]]
    )
  )

  for (case in cases) {
    result <- expect_no_warning(
      bootstrap_analysis_wrapper(
        boot_data = case$data,
        model = case$model,
        factor_cols = case$factor_cols,
        original_data = case$data,
        fit_weights = case$fit_weights
      )
    )

    expect_s3_class(result$model, case$class)
    expect_equal(result$model$weights, case$expected, info = case$label)
    expect_equal(case$stored(result$model), case$expected, info = case$label)
  }
})

test_that("bootstrap_analysis_wrapper decides when refits use starting values", {
  # A stored fit-time `coefstart` appears in the model call.
  model <- structure(
    list(coefficients = c(a = 1), call = quote(fit(coefstart = start))),
    class = "vglm"
  )
  original_data <- data.frame(y = factor(1:3))
  drop_requests <- logical()
  refit_with <- function(boot_data, fail_with_coefstart) {
    drop_requests <<- logical()
    with_mocked_bindings(
      update_bootstrap_model = function(
        model,
        boot_data,
        fit_weights = NULL,
        coefstart = NULL,
        drop_coefstart = FALSE
      ) {
        drop_requests <<- c(drop_requests, drop_coefstart)
        if (fail_with_coefstart && !drop_coefstart) {
          stop("non-conformable arguments")
        }
        model
      },
      bootstrap_analysis_wrapper(
        boot_data = boot_data,
        model = model,
        factor_cols = "y",
        original_data = original_data
      )
    )
  }

  # Every state present: keep the stored starting values.
  result <- refit_with(original_data, fail_with_coefstart = FALSE)
  expect_false(is.null(result$model))
  expect_equal(drop_requests, FALSE)

  # A failed refit with starting values is retried once without them.
  result <- refit_with(original_data, fail_with_coefstart = TRUE)
  expect_false(is.null(result$model))
  expect_equal(drop_requests, c(FALSE, TRUE))

  # A missing state: refit without starting values from the start.
  result <- refit_with(original_data[1:2, , drop = FALSE], FALSE)
  expect_length(result$missing_states, 1L)
  expect_equal(drop_requests, TRUE)
})
