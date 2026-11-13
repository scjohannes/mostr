test_that("RMS basis registry compiles rcs and lsp bases", {
  skip_if_not_installed("rms")
  skip_if_not_installed("Hmisc")

  handlers <- mostr:::rms_basis_registry()
  expect_equal(names(handlers), c("rcs", "lsp"))

  values <- c(-1, 0, 1, 3, 6, NA_real_)

  lsp <- mostr:::new_compiled_rms_basis(
    handler = "lsp",
    variable = "time",
    parameters = c(1, 4),
    columns = c("time", "time'", "time''")
  )
  actual_lsp <- mostr:::evaluate_compiled_rms_basis(lsp, values)
  expected_lsp <- cbind(values, pmax(values - 1, 0), pmax(values - 4, 0))
  expect_equal(as.vector(actual_lsp), as.vector(expected_lsp))
  expect_equal(attr(actual_lsp, "nonlinear"), c(FALSE, TRUE, TRUE))

  knots <- c(0, 1, 4, 5)
  rcs <- mostr:::new_compiled_rms_basis(
    handler = "rcs",
    variable = "time",
    parameters = knots,
    columns = c("time", "time'", "time''")
  )
  actual_rcs <- mostr:::evaluate_compiled_rms_basis(rcs, values)
  expected_rcs <- Hmisc::rcspline.eval(values, knots = knots, inclx = TRUE)
  expect_equal(
    as.vector(actual_rcs),
    as.vector(expected_rcs),
    tolerance = 1e-13
  )
})

test_that("a third RMS basis handler extends evaluation without dispatch changes", {
  hinge2 <- list(
    id = "hinge2",
    helper = "hinge2",
    assume = "hinge2",
    assume_code = 99L,
    pattern = "(^|:)\\s*hinge2\\(",
    evaluate = function(values, knots) {
      cbind(values, pmax(values - knots[1L], 0)^2)
    }
  )
  withr::local_options(
    mostr.rms_basis_handlers = list(hinge2 = hinge2)
  )
  basis <- mostr:::new_compiled_rms_basis(
    "hinge2",
    variable = "time",
    parameters = 2,
    columns = c("time", "time_hinge2")
  )
  actual <- mostr:::evaluate_compiled_rms_basis(basis, c(1, 2, 4))

  expect_equal(as.vector(actual), as.vector(cbind(c(1, 2, 4), c(0, 0, 4))))
  expect_equal(attr(actual, "nonlinear"), c(FALSE, TRUE))
  expect_equal(
    mostr:::rms_basis_handler_for_design("hinge2", 99L)$id,
    "hinge2"
  )
  expect_equal(
    mostr:::rms_basis_handlers_for_term("hinge2(time)"),
    "hinge2"
  )
})

test_that("RMS basis registry recognizes fitted Design metadata", {
  expect_equal(
    mostr:::rms_basis_handler_for_design("rcspline", 4L)$id,
    "rcs"
  )
  expect_equal(
    mostr:::rms_basis_handler_for_design("lspline", 3L)$id,
    "lsp"
  )
  expect_null(mostr:::rms_basis_handler_for_design("polynomial", 2L))
})

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

  # The vglm engine splits the lsp basis into one column per piece and records
  # them, so prediction can rebuild the basis.
  basis_terms <- attr(vglm_fit, "markov_basis_terms", exact = TRUE)
  expect_length(basis_terms, 1L)
  expect_equal(basis_terms[[1L]]$handlers, "lsp")
  expect_length(basis_terms[[1L]]$columns, 3L)
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
      c(list(case$plan_model), internal_args, list(builder = "batched"))
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
