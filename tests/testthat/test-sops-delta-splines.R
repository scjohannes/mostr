local_spline_delta_orm_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(value)) {
      state_levels <- as.character(seq_len(8L))
      raw_data <- make_markov_trajectories(
        n_patients = 48,
        follow_up_time = 6,
        n_states = length(state_levels),
        thresholds = stats::qnorm(seq(0.05, 0.95, length.out = 7L)),
        allowed_start_state = as.integer(2:7),
        absorbing_state = NULL,
        seed = 2701
      )
      data <- prepare_markov_data(
        raw_data,
        absorbing_state = NULL,
        factor_previous = FALSE
      )
      data$y <- ordered(data$y, levels = state_levels)

      model <- mostm(
        y ~ rms::rcs(time, 3) + tx + rms::rcs(yprev, 4),
        data = data,
        id_var = "id",
        opt_method = "LM",
        scale = TRUE,
        penalty = 0.1
      )
      value <<- list(data = data, model = model)
    }
    value
  }
})

replay_spline_avg_sops <- function(object, model) {
  avg_args <- attr(object, "avg_args")
  settings <- markov_result_settings(object)
  avg_sops(
    model = model,
    variables = avg_args$variables,
    # Stored averages use `by = NULL` for the overall average.
    by = avg_args$by %||% TRUE,
    times = avg_args$times,
    y_levels = settings$y_levels,
    absorb = settings$absorb,
    id_var = settings$id_var,
    time_var = settings$time_var,
    p_var = settings$p_var,
    p2_var = settings$p2_var,
    gap_var = settings$gap_var,
    time_covariates = settings$time_covariates
  )
}

central_spline_avg_jacobian <- function(object, model) {
  avg_args <- attr(object, "avg_args")
  key_columns <- c("time", "state", names(avg_args$variables), avg_args$by)
  object_key <- sop_draw_cell_key(object, key_columns)
  numeric_jacobian(model, function(model) {
    result <- replay_spline_avg_sops(object, model)
    index <- match(object_key, sop_draw_cell_key(result, key_columns))
    if (anyNA(index)) {
      stop("Finite-difference spline SOP cells did not align.")
    }
    result$estimate[index]
  })
}

test_that("backend orm spline covariance labels preserve coefficient order", {
  fit <- local_spline_delta_orm_case()$model
  coefficient_names <- names(stats::coef(fit))
  q <- length(coefficient_names)
  expected <- diag(seq_len(q))

  unnamed_fit <- fit
  unnamed_fit$var <- expected
  unnamed_fit$orig.var <- expected
  unnamed_metadata <- attr(unnamed_fit, "markov_robust_covariance")
  unnamed_metadata$covariance_identity <- expected
  attr(unnamed_fit, "markov_robust_covariance") <- unnamed_metadata
  unnamed <- get_delta_cluster_vcov(unnamed_fit)
  expect_identical(rownames(unnamed$vcov), coefficient_names)
  expect_identical(colnames(unnamed$vcov), coefficient_names)
  expect_identical(unname(diag(unnamed$vcov)), as.numeric(seq_len(q)))

  rms_names <- c(
    coefficient_names[seq_len(fit$non.slopes)],
    fit$Design$mmcolnames
  )
  aliased_fit <- fit
  aliased_fit$var <- expected
  aliased_fit$orig.var <- expected
  dimnames(aliased_fit$var) <- list(rms_names, rms_names)
  dimnames(aliased_fit$orig.var) <- list(rms_names, rms_names)
  aliased_metadata <- attr(aliased_fit, "markov_robust_covariance")
  aliased_metadata$covariance_identity <- aliased_fit$var
  attr(aliased_fit, "markov_robust_covariance") <- aliased_metadata
  aliased <- get_delta_cluster_vcov(aliased_fit)
  expect_identical(rownames(aliased$vcov), coefficient_names)
  expect_identical(colnames(aliased$vcov), coefficient_names)
  expect_identical(unname(diag(aliased$vcov)), as.numeric(seq_len(q)))

  condition <- tryCatch(
    get_delta_cluster_vcov(fit, vcov = aliased_fit$var),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "must be uniquely named")
})

test_that("unnamed orm bread preserves transformed and interaction coefficient order", {
  skip_on_cran()
  data <- local_spline_delta_orm_case()$data
  data$yprev_factor <- factor(data$yprev)
  formulas <- list(
    y ~ rms::rcs(time, 3) * tx + rms::rcs(yprev, 4),
    y ~ log(time + 1) * tx + yprev_factor
  )

  for (formula in formulas) {
    for (variance in c("simple", "sandwich")) {
      fit <- rms::orm(
        formula,
        data = data,
        x = TRUE,
        y = TRUE,
        penalty = 0.1,
        scale = TRUE,
        var.penalty = variance
      )
      # Re-evaluate information without scaling so the backend retains labels.
      reference <- rms::orm.fit(
        x = fit$x,
        y = fit$y,
        initial = stats::coef(fit),
        penalty.matrix = fit$penalty.matrix,
        scale = FALSE,
        maxit = 1,
        compstats = FALSE
      )
      expect_equal(stats::coef(reference), stats::coef(fit), tolerance = 1e-12)
      expected <- stats::vcov(reference, intercepts = "all")
      coefficient_names <- names(stats::coef(fit))
      expect_identical(rownames(expected), coefficient_names)
      expect_identical(colnames(expected), coefficient_names)
      expect_identical(
        colnames(fit$x),
        coefficient_names[-seq_len(fit$non.slopes)]
      )
      expect_equal(orm_model_bread(fit)$bread, expected, tolerance = 1e-8)
    }
  }
})

test_that("penalized orm previous-state splines support delta inference", {
  skip_on_cran()
  case <- local_spline_delta_orm_case()
  fit <- case$model
  avg <- avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    y_levels = fit$yunique,
    p_var = "yprev"
  )

  plan <- delta_compile_plan(
    object = avg,
    model_plan = fit,
    newdata = attr(avg, "newdata_pred"),
    times = attr(avg, "avg_args")$times,
    output = "average"
  )
  native <- run_sop_delta_plan(plan, fit)
  oracle <- run_sop_delta_r_oracle(plan, fit)
  expect_equal(native$probabilities, oracle$probabilities, tolerance = 1e-12)
  expect_equal(native$jacobian, oracle$jacobian, tolerance = 1e-12)

  inferred <- inferences(avg, method = "delta", vcov = "conditional")
  numerical <- central_spline_avg_jacobian(avg, fit)
  analytical <- get_jacobian(inferred)

  expect_identical(colnames(analytical), names(stats::coef(fit)))
  expect_lt(max(abs(analytical - numerical)), 1e-7)
  expect_all_true(is.finite(inferred$std.error))
  expect_all_true(inferred$std.error >= 0)

  condition <- tryCatch(
    inferences(avg, method = "delta", vcov = "unconditional"),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(
    conditionMessage(condition),
    "does not currently support penalized orm likelihoods",
    fixed = TRUE
  )
})
