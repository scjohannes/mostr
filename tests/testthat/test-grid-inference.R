# Inference on every engine, time coding and method: orm and vglm, linear,
# spline and factor time, and MVN draws, the ordinary and the fractional
# weighted bootstrap, and analytical (delta) conditional and unconditional
# variance, all on avg_sops(). That is 2 x 3 x 5 = 30 set-ups. The checks are
# in helper-grid-inference.R; test-grid-inference-array.R covers the other
# settings pairwise. Analytical inference runs on these first-order,
# full proportional-odds fits, which is its documented scope.
#
# On CRAN only MVN draws and analytical variance on one set-up per engine
# run; everything with bootstrap refits is skipped there.

inference_crossing_settings <- function() {
  methods <- list(
    list(method = "mvn"),
    list(method = "bootstrap"),
    list(method = "fwb"),
    list(method = "delta", vcov = "conditional"),
    list(method = "delta", vcov = "unconditional")
  )
  settings <- list()
  for (engine in c("orm", "vglm")) {
    for (time in c("linear", "rcs", "factor")) {
      for (method in methods) {
        settings[[length(settings) + 1L]] <- inference_setting(
          engine = engine,
          time = time,
          method = method$method,
          vcov = method$vcov
        )
      }
    }
  }
  settings
}

test_that("MVN and analytical inference hold their limits on every engine and time coding", {
  skip_on_cran()
  # L9, L10, L11, L13, L18, L21, L23, L24, R13, R14, N4, N21.
  settings <- Filter(
    function(setting) !setting$method %in% c("bootstrap", "fwb"),
    inference_crossing_settings()
  )
  if (grid_on_cran()) {
    settings <- Filter(
      function(setting) {
        (setting$engine == "orm" && setting$time == "linear") ||
          (setting$engine == "vglm" && setting$time == "rcs")
      },
      settings
    )
  }
  for (setting in settings) check_inference_setting(setting)
})

test_that("bootstrap inference holds its limits on every engine and time coding", {
  # As above, plus L29 bootstrap bookkeeping.
  skip_on_cran()
  settings <- Filter(
    function(setting) setting$method %in% c("bootstrap", "fwb"),
    inference_crossing_settings()
  )
  for (setting in settings) check_inference_setting(setting)
})

test_that("weighted profile draws of the fractional bootstrap reproduce the averaged draws (R12)", {
  # With the stored patients, each draw of avg_sops() and avg_comparisons()
  # is the mean of that draw's profile results, weighted by `fwb_weight`.
  # This holds for differences; a weighted mean of ratios is not the ratio
  # of weighted means, so ratios are not checked.
  skip_on_cran()
  checks <- list(
    list(cell = grid_cell("orm"), results = c("sop", "time_contrast")),
    list(cell = grid_cell("vglm", time = "factor"), results = c("contrast", "time"))
  )
  for (check in checks) {
    cell <- check$cell
    label <- grid_label(cell)
    run <- function(fun, ...) {
      x <- grid_call(fun, cell, variables = grid_treatments, ...)
      as.data.frame(get_draws(
        inferences(x, method = "fwb", n_draws = 3, seed = inference_seed)
      ))
    }
    pairs <- list(
      sop = function() list(run(sops), run(avg_sops), c("tx", "time", "state")),
      time = function() {
        list(run(state_time), run(avg_state_time), c("tx", "state_set"))
      },
      contrast = function() {
        list(
          run(comparisons, estimand = "sop"),
          run(avg_comparisons, estimand = "sop"),
          c("time", "state_set")
        )
      },
      time_contrast = function() {
        list(
          run(comparisons, estimand = "time_in_state"),
          run(avg_comparisons, estimand = "time_in_state"),
          "state_set"
        )
      }
    )
    for (name in check$results) {
      pair <- pairs[[name]]()
      each <- pair[[1]]
      average <- pair[[2]]
      keys <- c(pair[[3]], "draw_id")
      expect_true("fwb_weight" %in% names(each), label = paste(label, name))
      each$weighted <- each$draw * each$fwb_weight
      sums <- stats::aggregate(each[c("weighted", "fwb_weight")], each[keys], sum)
      sums$draw <- sums$weighted / sums$fwb_weight
      expect_lt(
        grid_gap(average, sums, keys, "draw"),
        1e-12,
        label = paste(label, name)
      )
    }
  }
})

test_that("`type` and `cadjust` change conditional but not unconditional variance (N6)", {
  skip_on_cran()
  # `type` and `cadjust` set the small-sample correction of the patient-
  # clustered coefficient covariance, which conditional variance uses.
  # Unconditional variance computes its own covariance from each patient's
  # contribution, so the fit-time correction must not reach it.
  for (engine in c("orm", "vglm")) {
    cell <- grid_cell(engine)
    data <- grid_cell_data(cell)
    fit_with <- function(type, cadjust) {
      args <- list(
        if (engine == "orm") y ~ time + tx + age + yprev else ordered(y) ~ time + tx + age + yprev,
        data = data,
        id_var = "id",
        first_followup_time = 1,
        absorb = grid_n_states,
        type = type,
        cadjust = cadjust
      )
      if (engine == "vglm") {
        args$engine <- "vglm"
        args$family <- VGAM::cumulative(reverse = TRUE, parallel = TRUE)
      }
      suppressWarnings(do.call(mostm, args))
    }
    variances <- lapply(
      list(hc0 = fit_with("HC0", FALSE), hc1 = fit_with("HC1", TRUE)),
      function(fit) {
        average <- avg_sops(
          fit,
          variables = grid_treatments,
          times = 1:3
        )
        list(
          conditional = inferences(average, method = "delta", vcov = "conditional"),
          unconditional = inferences(average, method = "delta", vcov = "unconditional")
        )
      }
    )
    expect_gt(
      max(abs(
        variances$hc0$conditional$std.error - variances$hc1$conditional$std.error
      )),
      1e-8,
      label = paste(engine, "conditional variance follows the correction")
    )
    expect_equal(
      variances$hc0$unconditional$std.error,
      variances$hc1$unconditional$std.error,
      tolerance = 1e-12,
      info = paste(engine, "unconditional")
    )
  }
})

# Free the cached results; other files do not use them.
grid_cache$results <- NULL
