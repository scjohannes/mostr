# Native analytical recursion (cpp_run_sop_delta) versus the independent
# matrix-product reference in helper-sops-delta-oracle.R.

delta_synthetic_case <- function(n = 6L, K = 5L, times = 4L, seed = 20260928L) {
  set.seed(seed)
  M <- K - 1L
  col_names <- c("x1", "group_b", "(Intercept)", "group_c", "x2")
  P <- length(col_names)
  theta <- c(a = 1.1, d = 0.8, s1 = 0.45, s2 = -0.6, s3 = 0.35, s4 = 0.3)
  map <- matrix(0, nrow = M * P, ncol = length(theta))
  colnames(map) <- names(theta)
  rows <- function(column) (column - 1L) * M + seq_len(M)
  map[rows(1L), "s1"] <- 1
  map[rows(2L), "s2"] <- 1
  # Equidistant thresholds: alpha_j = a - (j - 1) d.
  map[rows(3L), "a"] <- 1
  map[rows(3L), "d"] <- -(seq_len(M) - 1)
  # One design column depending on two raw coefficients.
  map[rows(4L), "s2"] <- 0.5
  map[rows(4L), "s3"] <- -1
  # A scaled common slope.
  map[rows(5L), "s4"] <- -2

  group <- sample(c("a", "b", "c"), n, replace = TRUE)
  design <- function(x1, x2) {
    X <- cbind(
      x1,
      as.numeric(group == "b"),
      1,
      as.numeric(group == "c"),
      x2
    )
    colnames(X) <- col_names
    X
  }
  x1 <- stats::rnorm(n)
  X_init <- design(x1, stats::rnorm(n))
  # Transition blocks are stored for every state, in reverse state order.
  origins <- rev(seq_len(K))
  X_transition <- vector("list", times)
  if (times >= 2L) {
    for (visit in 2:times) {
      X_transition[[visit]] <- do.call(
        rbind,
        lapply(origins, function(origin) {
          design(x1 + 0.2 * visit, 0.4 * origin + stats::rnorm(n, sd = 0.3))
        })
      )
    }
  }
  list(
    X_init = X_init,
    X_transition = X_transition,
    theta = theta,
    map = map,
    gamma = matrix(drop(map %*% theta), nrow = M, ncol = P),
    origins = origins,
    intercept = 3L,
    K = K
  )
}

run_delta_native <- function(case, absorb, group_size = 0L, retain = FALSE) {
  non_absorb <- setdiff(seq_len(case$K), absorb)
  cpp_run_sop_delta(
    initial_design = case$X_init,
    transition_designs = case$X_transition,
    gamma = case$gamma,
    coefficient_map = case$map,
    origin_positions = match(non_absorb, case$origins),
    non_absorb = as.integer(non_absorb),
    absorb = as.integer(absorb),
    intercept = case$intercept,
    group_size = as.integer(group_size),
    retain_individual_probabilities = retain
  )
}

run_delta_reference <- function(case, absorb, group_size = 0L) {
  reference_sop_delta(
    case$X_init,
    case$X_transition,
    case$gamma,
    case$map,
    case$origins,
    absorb,
    group_size = group_size
  )
}

delta_absorbing_sets <- list(
  none = integer(),
  last = 5L,
  first_and_last = c(1L, 5L),
  middle = 3L,
  two_inner = c(2L, 4L)
)

test_that("native delta recursion matches the reference on synthetic designs", {
  case <- delta_synthetic_case()

  for (absorb in delta_absorbing_sets) {
    native <- run_delta_native(case, absorb)
    reference <- run_delta_reference(case, absorb)

    expect_identical(native$status, c(0L, 0L, 0L))
    expect_identical(dim(native$probabilities), c(6L, 4L, 5L))
    expect_identical(dim(native$jacobian), c(6L, 4L, 5L, 6L))
    expect_equal(
      native$probabilities,
      reference$probabilities,
      tolerance = 1e-12
    )
    expect_equal(native$jacobian, reference$jacobian, tolerance = 1e-12)
    expect_length(native$individual_probabilities, 0L)
  }
})

test_that("native delta recursion matches the reference for grouped averages", {
  case <- delta_synthetic_case()

  grouped_sets <- delta_absorbing_sets[c("none", "first_and_last", "two_inner")]
  for (absorb in grouped_sets) {
    for (group_size in c(1L, 2L, 3L, 6L)) {
      for (retain in c(FALSE, TRUE)) {
        native <- run_delta_native(case, absorb, group_size, retain)
        reference <- run_delta_reference(case, absorb, group_size)

        expect_identical(native$status, c(0L, 0L, 0L))
        expect_identical(
          dim(native$probabilities),
          as.integer(c(6L / group_size, 4L, 5L))
        )
        expect_equal(
          native$probabilities,
          reference$probabilities,
          tolerance = 1e-12
        )
        expect_equal(native$jacobian, reference$jacobian, tolerance = 1e-12)
        if (retain) {
          expect_equal(
            native$individual_probabilities,
            reference$individual_probabilities,
            tolerance = 1e-12
          )
        } else {
          expect_length(native$individual_probabilities, 0L)
        }
      }
    }
  }
})

test_that("native delta recursion matches the reference for one visit", {
  case <- delta_synthetic_case(times = 1L)
  for (absorb in delta_absorbing_sets[c("none", "first_and_last")]) {
    native <- run_delta_native(case, absorb)
    reference <- run_delta_reference(case, absorb)

    expect_identical(native$status, c(0L, 0L, 0L))
    expect_equal(
      native$probabilities,
      reference$probabilities,
      tolerance = 1e-12
    )
    expect_equal(native$jacobian, reference$jacobian, tolerance = 1e-12)
  }
})

test_that("native delta recursion reports crossed probabilities unclipped", {
  case <- delta_synthetic_case()
  case$theta["d"] <- -0.8
  case$gamma[] <- drop(case$map %*% case$theta)

  native <- run_delta_native(case, absorb = 5L)
  reference <- reference_category_probabilities(case$X_init %*% t(case$gamma))

  expect_lt(min(reference), 0)
  expect_identical(native$status, c(2L, 1L, 0L))
  expect_equal(native$value, min(reference), tolerance = 1e-12)
})

delta_reference_fits <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      data <- make_test_data(n_patients = 80, follow_up_time = 6, seed = 7331)
      orm_fit <- suppressWarnings(orm_markov(
        y ~ time + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      ))
      vglm_fit <- make_test_model(data)
      equidistant_fit <- suppressWarnings(vglm_markov(
        ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
        family = VGAM::cumulative(
          reverse = TRUE,
          parallel = TRUE,
          Thresh = "equid"
        ),
        data = data
      ))
      baseline <- data[!duplicated(data$id), , drop = FALSE]
      baseline <- baseline[1:4, , drop = FALSE]
      value <<- list(
        baseline = baseline,
        models = list(
          orm = orm_fit,
          vglm = vglm_fit,
          vglm_equidistant = equidistant_fit
        )
      )
    }
    value
  }
})

delta_factor_time_data <- function(n_patients = 80, n_visits = 4, seed = 7341) {
  set.seed(seed)
  visits <- as.character(seq_len(n_visits))
  data <- expand.grid(
    id = seq_len(n_patients),
    time = visits,
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$time <- factor(data$time, levels = visits)
  data$tx <- stats::rbinom(n_patients, 1, 0.5)[data$id]
  data$yprev <- factor(
    sample(seq_len(4), nrow(data), replace = TRUE),
    levels = as.character(seq_len(4))
  )
  eta <- -0.35 * data$tx + 0.25 * as.integer(data$time) +
    0.2 * as.integer(as.character(data$yprev))
  y <- 1L + stats::rbinom(nrow(data), 3L, stats::plogis(eta - mean(eta)))
  data$y <- ordered(y, levels = as.character(seq_len(4)))
  data
}

expect_delta_plan_matches_reference <- function(plan, model) {
  components <- plan$components
  inputs <- reference_delta_model_inputs(model, components$col_names)
  expect_equal(
    inputs$gamma,
    unname(reference_backend_gamma(model, components$col_names)),
    tolerance = 1e-12
  )
  absorb <- which(plan$y_levels %in% plan$absorb)

  native <- run_sop_delta_plan(plan, model)
  reference <- reference_sop_delta(
    components$X_init,
    components$X_transition,
    inputs$gamma,
    inputs$map,
    components$transition_origins,
    absorb
  )
  expect_identical(native$coef_names, names(inputs$theta))
  expect_equal(
    unname(native$probabilities),
    reference$probabilities,
    tolerance = 1e-12
  )
  expect_equal(unname(native$jacobian), reference$jacobian, tolerance = 1e-12)

  grouped <- run_sop_delta_plan(
    plan,
    model,
    average_group_size = 2L,
    retain_individual_probabilities = TRUE
  )
  grouped_reference <- reference_sop_delta(
    components$X_init,
    components$X_transition,
    inputs$gamma,
    inputs$map,
    components$transition_origins,
    absorb,
    group_size = 2L
  )
  expect_equal(
    unname(grouped$probabilities),
    grouped_reference$probabilities,
    tolerance = 1e-12
  )
  expect_equal(
    unname(grouped$jacobian),
    grouped_reference$jacobian,
    tolerance = 1e-12
  )
  expect_equal(
    unname(grouped$individual_probabilities),
    grouped_reference$individual_probabilities,
    tolerance = 1e-12
  )
}

test_that("fitted orm and vglm delta plans match the independent reference", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  fits <- delta_reference_fits()

  for (model in fits$models) {
    for (absorb in list(6, c(1, 6), c(3, 6))) {
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = fits$baseline,
        times = 1:4,
        y_levels = 1:6,
        absorb = absorb
      )
      expect_delta_plan_matches_reference(plan, model)
    }
  }
})

test_that("factor-time delta plans match the independent reference", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  data <- delta_factor_time_data()
  baseline <- data[!duplicated(data$id), , drop = FALSE][1:4, , drop = FALSE]
  models <- list(
    orm = suppressWarnings(orm_markov(
      y ~ time + tx + yprev,
      data = data,
      x = TRUE,
      y = TRUE
    )),
    vglm = suppressWarnings(vglm_markov(
      y ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    ))
  )

  for (model in models) {
    for (absorb in list(NULL, 2)) {
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = baseline,
        times = levels(data$time),
        y_levels = seq_len(4),
        absorb = absorb
      )
      expect_delta_plan_matches_reference(plan, model)
    }
  }
})
