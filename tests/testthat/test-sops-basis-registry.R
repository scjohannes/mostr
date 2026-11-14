test_that("formula helper injection includes every registered RMS basis and leaves non-formulas alone", {
  skip_if_not_installed("rms")

  stripped <- stats::as.formula(
    "y ~ rcs(time, 4) + lsp(age, c(30, 60))",
    env = baseenv()
  )
  no_environment <- y ~ x
  environment(no_environment) <- NULL

  cases <- list(
    list(
      label = "a formula whose environment lacks rms",
      input = stripped,
      helpers = c("rcs", "lsp", "%ia%")
    ),
    list(
      label = "a formula without an environment",
      input = no_environment,
      helpers = c("rcs", "lsp")
    ),
    list(label = "a number", input = 1, helpers = NULL)
  )

  for (case in cases) {
    out <- mostr:::add_rms_formula_helpers(case$input)
    if (is.null(case$helpers)) {
      expect_equal(out, case$input, info = case$label)
      next
    }
    expect_s3_class(out, "formula")
    for (helper in case$helpers) {
      expect_true(
        exists(helper, envir = environment(out), inherits = TRUE),
        info = paste(case$label, helper)
      )
    }
  }
})

test_that("lsp time bases give the same SOPs through the compiled plan and the R engine", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2471, follow_up_time = 12)
  baseline <- data[data$time == 1, , drop = FALSE][1:8, , drop = FALSE]
  vglm_fit <- suppressWarnings(
    mostr::mostm(
      stats::as.formula(
        "ordered(y) ~ lsp(time, c(3, 7)) + tx + yprev",
        env = baseenv()
      ),
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      absorb = 6,
      engine = "vglm"
    )
  )

  # The vglm engine splits the lsp basis into one column per piece.
  expect_equal(attr(vglm_fit, "markov_split_assign", exact = TRUE), TRUE)
  expect_length(
    grep("lsp\\(", names(vglm_fit$vglm_fit@constraints), value = TRUE),
    3L
  )

  orm_fit <- suppressWarnings(mostm(
    ordered(y) ~ rms::lsp(time, c(3, 7)) + tx + yprev,
    data = data,
    id_var = "id",
    absorb = 6
  ))
  # The plan compiler takes the fitted model inside the vglm engine's wrapper.
  cases <- list(
    list(label = "vglm", fit = vglm_fit, plan_model = vglm_fit$vglm_fit),
    list(label = "orm", fit = orm_fit, plan_model = orm_fit)
  )
  args <- list(newdata = baseline, times = 1:8)

  for (case in cases) {
    internal_args <- list(
      newdata = baseline,
      times = 1:8,
      settings = mostr:::resolve_markov_settings(case$fit)
    )
    reference <- do.call(
      mostr:::soprob_markov_r_engine,
      c(list(case$fit), internal_args)
    )
    expect_equal(
      do.call(mostr::soprob_markov, c(list(case$fit), args)),
      reference,
      tolerance = 1e-11,
      info = case$label
    )

    plan <- do.call(
      mostr:::compile_sop_execution_plan,
      c(list(case$plan_model), internal_args)
    )
    expect_s3_class(plan, "markov_sop_exec_plan")
    planned <- mostr:::run_sop_execution_plan(
      plan,
      mostr:::get_effective_coefs(case$plan_model)
    )
    expect_equal(
      unname(planned),
      unname(reference),
      tolerance = 1e-11,
      info = case$label
    )
  }
})
