delta_vglm_case <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 100,
        follow_up_time = 8,
        seed = 7301
      )
      model <- make_test_model(data)$vglm_fit
      y_levels <- seq_len(nrow(get_effective_coefs(model)) + 1L)
      baseline <- data[!duplicated(data$id), , drop = FALSE][
        1:4,
        ,
        drop = FALSE
      ]
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = baseline,
        times = 1:5,
        y_levels = y_levels,
        absorb = max(y_levels)
      )
      value <<- list(
        data = data,
        model = model,
        baseline = baseline,
        y_levels = y_levels,
        plan = plan
      )
    }
    value
  }
})

delta_orm_case <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 100,
        follow_up_time = 8,
        seed = 7311
      )
      model <- suppressWarnings(mostm(
        y ~ time + tx + yprev,
        data = data,
        id_var = "id"
      ))
      baseline <- data[!duplicated(data$id), , drop = FALSE][
        1:3,
        ,
        drop = FALSE
      ]
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = baseline,
        times = 1:4,
        y_levels = model$yunique,
        absorb = max(model$yunique)
      )
      value <<- list(
        data = data,
        model = model,
        baseline = baseline,
        y_levels = model$yunique,
        plan = plan
      )
    }
    value
  }
})

central_sop_jacobian <- function(plan, model) {
  point <- run_sop_execution_plan(plan, get_effective_coefs(model))
  out <- numeric_jacobian(
    model,
    function(model) run_sop_execution_plan(plan, get_effective_coefs(model))
  )
  states <- dimnames(point)[[3L]]
  dimnames(out) <- list(NULL, NULL, states, names(get_coef(model)))
  out
}

test_that("VGLM analytical SOP Jacobians match raw-coefficient differences", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  case <- delta_vglm_case()

  analytical <- run_sop_delta_plan(case$plan, case$model)
  point <- run_sop_execution_plan(
    case$plan,
    get_effective_coefs(case$model)
  )
  numeric <- central_sop_jacobian(case$plan, case$model)

  expect_named(
    analytical,
    c("probabilities", "jacobian", "coef", "coef_names")
  )
  expect_identical(
    dim(analytical$jacobian),
    c(dim(analytical$probabilities), length(get_coef(case$model)))
  )
  expect_lt(
    max(abs(unname(analytical$probabilities) - unname(point))),
    1e-12
  )
  expect_lt(max(abs(unname(analytical$jacobian) - unname(numeric))), 1e-7)
  expect_lt(
    max(abs(apply(analytical$probabilities, c(1L, 2L), sum) - 1)),
    1e-10
  )
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )
  expect_identical(analytical$coef_names, names(get_coef(case$model)))
  expect_identical(names(analytical$coef), analytical$coef_names)
})

test_that("ORM analytical SOP Jacobians match raw-coefficient differences", {
  skip_if_not_installed("rms")
  case <- delta_orm_case()

  analytical <- run_sop_delta_plan(case$plan, case$model)
  point <- run_sop_execution_plan(
    case$plan,
    get_effective_coefs(case$model)
  )
  numeric <- central_sop_jacobian(case$plan, case$model)

  expect_lt(
    max(abs(unname(analytical$probabilities) - unname(point))),
    1e-12
  )
  expect_lt(max(abs(unname(analytical$jacobian) - unname(numeric))), 1e-7)
  expect_lt(
    max(abs(apply(analytical$probabilities, c(1L, 2L), sum) - 1)),
    1e-10
  )
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )
})

test_that("effective coefficient maps exactly preserve named raw coefficients", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(7307)
  data <- data.frame(
    x = stats::rnorm(240),
    group = factor(rep(letters[1:3], 80))
  )
  data$y <- ordered(cut(
    0.3 * data$x + stats::rlogis(240),
    c(-Inf, -1, 0, 1, Inf)
  ))
  spline_model <- VGAM::vglm(
    y ~ rms::rcs(x, 3) * group,
    VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  equidistant_model <- VGAM::vglm(
    y ~ group + I(x^2) + x,
    VGAM::cumulative(reverse = TRUE, parallel = TRUE, Thresh = "equid"),
    data = data
  )
  # A non-unit common-slope basis is still proportional odds. Reparameterize
  # the fitted model without changing its effective coefficients.
  scaled_model <- equidistant_model
  scaled_model@constraints[["x"]] <- -2 * scaled_model@constraints[["x"]]
  scaled_model@coefficients["x"] <- scaled_model@coefficients["x"] / -2

  cases <- list(
    list(label = "vglm Markov fit", model = delta_vglm_case()$model),
    list(label = "orm Markov fit", model = delta_orm_case()$model),
    list(label = "vglm spline-by-group interaction", model = spline_model),
    list(label = "vglm equidistant thresholds", model = equidistant_model),
    list(label = "vglm scaled common-slope constraint", model = scaled_model)
  )
  for (case in cases) {
    map <- get_effective_coef_map(case$model)
    coef <- get_coef(case$model)
    gamma <- get_effective_coefs(case$model)

    expect_identical(colnames(map), names(coef), info = case$label)
    expect_identical(anyDuplicated(colnames(map)), 0L, info = case$label)
    expect_identical(attr(map, "gamma_dim"), dim(gamma), info = case$label)
    expect_identical(
      attr(map, "gamma_dimnames"),
      dimnames(gamma),
      info = case$label
    )
    expect_equal(
      unname(drop(map %*% coef)),
      unname(as.vector(gamma)),
      tolerance = 1e-12,
      info = case$label
    )
    if (!inherits(case$model, "vglm")) {
      next
    }
    # Each column of the map is the native expansion of one unit coefficient.
    beta <- stats::coef(case$model)
    expect_identical(colnames(map), names(beta), info = case$label)
    expect_identical(
      attr(map, "gamma_dimnames")[[2L]],
      rownames(stats::coef(case$model, matrix = TRUE)),
      info = case$label
    )
    for (j in seq_along(beta)) {
      basis <- case$model
      basis@coefficients[] <- 0
      basis@coefficients[j] <- 1
      expect_equal(
        unname(map[, j]),
        as.vector(t(stats::coef(basis, matrix = TRUE))),
        tolerance = 1e-12,
        info = paste(case$label, "coefficient", j)
      )
    }
  }

  zero_columns <- scaled_model
  zero_columns@constraints[["x"]] <- matrix(0, nrow = 3, ncol = 0)
  expect_error(get_effective_coef_map(zero_columns), "no coefficient columns")
})

test_that("absorbing states propagate probability and derivative mass", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  case <- delta_vglm_case()
  analytical <- run_sop_delta_plan(case$plan, case$model)
  absorbing_state <- length(case$y_levels)

  absorbing_probability <- analytical$probabilities[,, absorbing_state]
  expect_true(all(apply(absorbing_probability, 1L, diff) >= -1e-12))
  expect_gt(max(abs(analytical$jacobian[, -1L, absorbing_state, ])), 0)
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )

  wrapper <- structure(
    list(
      coefficients = get_coef(case$model),
      vglm_fit = case$model
    ),
    class = "robcov_vglm"
  )
  robust_result <- run_sop_delta_plan(case$plan, wrapper)
  expect_equal(
    unname(robust_result$jacobian),
    unname(analytical$jacobian),
    tolerance = 0
  )
})

test_that("changing covariate units preserves analytical probabilities", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  for (case in list(delta_vglm_case(), delta_orm_case())) {
    original <- run_sop_delta_plan(case$plan, case$model)
    plan <- case$plan
    plan$components$X_init[, "tx"] <- plan$components$X_init[, "tx"] * 1e12
    plan$components$X_transition <- lapply(
      plan$components$X_transition,
      function(design) {
        if (!is.null(design)) {
          design[, "tx"] <- design[, "tx"] * 1e12
        }
        design
      }
    )
    coef <- get_coef(case$model)
    coef["tx"] <- coef["tx"] / 1e12
    model <- set_coef(case$model, coef)
    native <- run_sop_delta_plan(plan, model)
    oracle <- run_sop_delta_r_oracle(plan, model)

    expect_equal(
      native$probabilities,
      original$probabilities,
      tolerance = 1e-12
    )
    expect_equal(native$jacobian, oracle$jacobian, tolerance = 1e-12)
    native$jacobian[,,, "tx"] <- native$jacobian[,,, "tx"] / 1e12
    expect_equal(native$jacobian, original$jacobian, tolerance = 1e-12)
  }
})

test_that("grouped recursion retains only the required analytical state", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  case <- delta_vglm_case()
  plan <- case$plan
  plan$components$X_init <- plan$components$X_init[
    rep(seq_len(4L), each = 2L),
    ,
    drop = FALSE
  ]
  plan$components$n_pat <- 8L
  plan$components$X_transition <- lapply(
    plan$components$X_transition,
    function(design) {
      if (is.null(design)) {
        return(NULL)
      }
      origin_count <- nrow(design) / 4L
      rows <- unlist(lapply(seq_len(origin_count), function(origin) {
        source <- (origin - 1L) * 4L + seq_len(4L)
        rep(source, each = 2L)
      }))
      design[rows, , drop = FALSE]
    }
  )

  full <- run_sop_delta_r_oracle(plan, case$model)
  grouped <- run_sop_delta_plan(
    plan,
    case$model,
    average_group_size = 4L
  )
  unconditional <- run_sop_delta_plan(
    plan,
    case$model,
    average_group_size = 4L,
    retain_individual_probabilities = TRUE
  )

  expected_probability <- apply(full$probabilities[1:4, , ], c(2L, 3L), mean)
  expected_jacobian <- apply(full$jacobian[1:4, , , ], c(2L, 3L, 4L), mean)
  expect_identical(
    dim(grouped$probabilities),
    c(2L, 5L, length(case$y_levels))
  )
  expect_identical(
    dim(grouped$jacobian),
    c(2L, 5L, length(case$y_levels), length(get_coef(case$model)))
  )
  expect_null(grouped$individual_probabilities)
  expect_equal(
    unname(grouped$probabilities[1L, , ]),
    unname(expected_probability),
    tolerance = 1e-15
  )
  expect_equal(
    unname(grouped$jacobian[1L, , , ]),
    unname(expected_jacobian),
    tolerance = 1e-15
  )
  expect_equal(
    unname(grouped$probabilities[2L, , ]),
    unname(apply(full$probabilities[5:8, , ], c(2L, 3L), mean)),
    tolerance = 1e-15
  )
  expect_equal(
    unname(grouped$jacobian[2L, , , ]),
    unname(apply(full$jacobian[5:8, , , ], c(2L, 3L, 4L), mean)),
    tolerance = 1e-15
  )
  expect_identical(
    dim(unconditional$individual_probabilities),
    dim(full$probabilities)
  )
  expect_equal(
    unname(unconditional$individual_probabilities),
    unname(full$probabilities),
    tolerance = 1e-12
  )

  individual_bytes <- sop_delta_preflight(plan, length(get_coef(case$model)))
  grouped_bytes <- sop_delta_preflight(
    plan,
    length(get_coef(case$model)),
    average_group_size = 4L
  )
  unconditional_bytes <- sop_delta_preflight(
    plan,
    length(get_coef(case$model)),
    average_group_size = 4L,
    retain_individual_probabilities = TRUE
  )
  expect_lt(grouped_bytes, unconditional_bytes)
  expect_lt(unconditional_bytes, individual_bytes)
})

test_that("analytical SOP validation rejects unsupported model structures", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  local_reproducible_output(width = 80)
  case <- delta_vglm_case()

  second_order <- case$plan
  second_order$recursion_order <- 2L
  expect_snapshot(
    run_sop_delta_plan(second_order, case$model),
    error = TRUE
  )

  unsupported_link <- case$plan
  unsupported_link$link <- "probit"
  expect_snapshot(
    run_sop_delta_plan(unsupported_link, case$model),
    error = TRUE
  )

  unsupported_model <- stats::lm(mpg ~ wt, data = mtcars)
  expect_snapshot(
    run_sop_delta_plan(case$plan, unsupported_model),
    error = TRUE
  )

  # This small non-parallel fit has a singular model covariance, so it has no
  # patient-clustered covariance; the check needs only its coefficients.
  partial_model <- markov_with_coefficients_only_refit(suppressWarnings(mostm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
    data = case$data,
    id_var = "id",
    engine = "vglm"
  )))
  partial_plan <- compile_sop_execution_plan(
    model = partial_model,
    newdata = case$baseline,
    times = 1:3,
    y_levels = seq_len(nrow(get_effective_coefs(partial_model)) + 1L),
    absorb = max(case$y_levels)
  )
  expect_snapshot(
    run_sop_delta_plan(partial_plan, partial_model),
    error = TRUE
  )
})

test_that("crossed raw ordinal probabilities are rejected without clipping", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  local_reproducible_output(width = 80)
  case <- delta_vglm_case()
  constraints <- VGAM::constraints(case$model)
  intercept <- grep("Intercept", names(constraints), fixed = TRUE)
  preceding <- if (intercept == 1L) {
    0L
  } else {
    sum(vapply(constraints[seq_len(intercept - 1L)], ncol, integer(1)))
  }
  intercept_coef <- preceding + seq_len(ncol(constraints[[intercept]]))
  coef <- get_coef(case$model)
  coef[intercept_coef] <- seq(-4, 4, length.out = length(intercept_coef))
  crossed_model <- set_coef(case$model, coef)

  condition <- rlang::catch_cnd(
    run_sop_delta_plan(case$plan, crossed_model)
  )
  expect_s3_class(condition, "mostr_delta_crossed_probabilities")
  expect_lt(condition$minimum_probability, 0)
  expect_snapshot(
    run_sop_delta_plan(case$plan, crossed_model),
    error = TRUE
  )
})

# Methods for test-only classes that fail when the analytical state is subset
# or the result is copied, so the allocation guards must run before either.
`[.delta_subset_trap` <- function(x, ...) {
  stop(
    "Analytical state was subset before its allocation guard.",
    call. = FALSE
  )
}

as.data.frame.delta_finalize_trap <- function(x, ...) {
  stop("The result was copied before its retained-state guard.", call. = FALSE)
}

test_that("analytical allocation guard is typed and configurable", {
  skip_on_cran()
  local_reproducible_output(width = 80)
  object <- data.frame(time = 1, state = "1", estimate = 0.5)
  class(object) <- c("markov_avg_sops", "data.frame")
  jacobian <- matrix(1, nrow = 1, dimnames = list(NULL, "beta"))
  covariance <- matrix(1, dimnames = list("beta", "beta"))
  finalize <- function(object) {
    delta_finalize_result(
      object,
      standard_error = 1,
      conf_level = 0.95,
      conf_type = "logit",
      target = "empirical",
      analytical = list(
        representation = "coefficient",
        jacobian = jacobian,
        coefficient_vcov = covariance,
        coefficient_names = "beta",
        covariance_metadata = list(source = "test")
      )
    )
  }
  result <- finalize(object)

  coefficient_result <- result
  analytical <- attr(coefficient_result, "analytical")
  class(analytical$jacobian) <- c(
    "delta_subset_trap",
    class(analytical$jacobian)
  )
  attr(coefficient_result, "analytical") <- analytical

  influence_result <- result
  influence_state <- attr(influence_result, "analytical")
  influence_state$representation <- "influence"
  influence_state$jacobian <- NULL
  influence_state$coefficient_vcov <- NULL
  influence_state$influence <- matrix(1, nrow = 2, ncol = 1)
  class(influence_state$influence) <- c(
    "delta_subset_trap",
    class(influence_state$influence)
  )
  attr(influence_result, "analytical") <- influence_state

  guarded_object <- object
  class(guarded_object) <- c(
    "delta_finalize_trap",
    class(guarded_object)
  )

  withr::local_options(mostr.delta_max_bytes = 1)
  cases <- list(
    list(
      label = "direct byte check",
      call = function() delta_assert_bytes(16, "Test allocation"),
      required_bytes = 16
    ),
    list(
      label = "Jacobian accessor",
      call = function() get_jacobian(coefficient_result),
      required_bytes = 8
    ),
    list(
      label = "covariance rows from coefficients",
      call = function() stats::vcov(coefficient_result, rows = 1L),
      required_bytes = 16
    ),
    list(
      label = "covariance rows from patient influence",
      call = function() stats::vcov(influence_result, rows = 1L),
      required_bytes = 24
    ),
    list(
      label = "retained result state",
      call = function() finalize(guarded_object),
      required_bytes = 16
    )
  )
  for (case in cases) {
    condition <- tryCatch(case$call(), error = identity)
    expect_true(
      inherits(condition, "mostr_delta_too_large"),
      info = case$label
    )
    expect_identical(
      condition$required_bytes,
      case$required_bytes,
      info = case$label
    )
    expect_identical(condition$max_bytes, 1, info = case$label)
  }
  expect_snapshot(
    delta_assert_bytes(16, "Test allocation"),
    error = TRUE
  )
})
