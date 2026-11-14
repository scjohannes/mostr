# Native SOP, posterior, and sampling kernels versus the pure-R references in
# helper-sops-native-reference.R.

native_absorbing_sets <- list(
  none = integer(),
  last = 5L,
  first_and_last = c(1L, 5L),
  two_inner = c(2L, 4L)
)

# Starting probabilities with some zero-mass origins.
native_reference_previous <- function(n, K) {
  previous <- reference_random_stochastic(n, K)
  previous[1L, ] <- 0
  previous[1L, 2L] <- 1
  previous[2L, c(1L, 3L)] <- 0
  previous / rowSums(previous)
}

test_that("native logit updates match the reference and report crossed cut-points", {
  skip_on_cran()
  set.seed(4102)
  n <- 4L
  K <- 5L
  previous <- native_reference_previous(n, K)

  for (absorb in native_absorbing_sets) {
    non_absorb <- setdiff(seq_len(K), absorb)
    for (origin_order in list(non_absorb, rev(non_absorb))) {
      rows <- n * length(origin_order)
      logits <- t(apply(
        matrix(stats::rnorm(rows * (K - 1L)), nrow = rows),
        1L,
        function(values) sort(values, decreasing = TRUE)
      ))
      # Nearly equal cut-points: rounding makes one probability about -2e-13,
      # which is within the tolerance and set to zero.
      logits[1L, ] <- c(0.3, 0.3 + 1e-12, -0.2, -1.4)
      expect_lt(min(reference_category_probabilities(logits)), 0)
      actual <- markov_update_logits_native(
        previous,
        logits,
        as.integer(origin_order),
        as.integer(absorb)
      )
      blocks <- reference_origin_blocks(
        reference_checked_probabilities(logits),
        origin_order,
        n
      )
      expected <- reference_markov_step(previous, blocks, absorb)
      expect_equal(actual, expected, tolerance = 1e-12)

      # Increasing logits mean crossed cut-points: for patient 1 at the first
      # origin, which may hold no probability, and patient n at the last.
      logits[1L, ] <- seq(-1, 1, length.out = K - 1L)
      logits[rows, ] <- seq(-0.5, 0.5, length.out = K - 1L)
      expect_error(reference_checked_probabilities(logits), "cross")
      marked <- cpp_markov_update_logits(
        previous,
        logits,
        as.integer(origin_order),
        as.integer(absorb)
      )
      expect_identical(attr(marked, "crossed_profiles"), c(1L, n))
      expect_equal(
        attr(marked, "minimum_probability"),
        min(reference_category_probabilities(logits)),
        tolerance = 1e-12
      )
      expect_error(
        markov_update_logits_native(
          previous,
          logits,
          as.integer(origin_order),
          as.integer(absorb),
          time = 3
        ),
        "cross for 2 of the 4 profiles at time 3",
        class = "mostr_crossed_thresholds"
      )
    }
  }
})

test_that("native PO updates match the reference and report crossed cut-points", {
  set.seed(4103)
  n <- 4L
  K <- 5L
  previous <- native_reference_previous(n, K)

  for (absorb in native_absorbing_sets) {
    non_absorb <- setdiff(seq_len(K), absorb)
    scalar <- stats::rnorm(n * length(non_absorb), sd = 1.5)
    cutpoints <- c(2, 0.5, -0.4, -2.2)
    actual <- markov_update_po_native(
      previous,
      scalar,
      cutpoints,
      as.integer(non_absorb),
      as.integer(absorb)
    )
    blocks <- reference_origin_blocks(
      reference_checked_probabilities(outer(scalar, cutpoints, "+")),
      non_absorb,
      n
    )
    expected <- reference_markov_step(previous, blocks, absorb)
    expect_equal(actual, expected, tolerance = 1e-12)

    # Cut-points out of order cross for every patient.
    crossed <- c(1, -1, 0.5, -2)
    expect_error(
      reference_checked_probabilities(outer(scalar, crossed, "+")),
      "cross"
    )
    marked <- cpp_markov_update_po(
      previous,
      scalar,
      crossed,
      as.integer(non_absorb),
      as.integer(absorb)
    )
    expect_identical(attr(marked, "crossed_profiles"), seq_len(n))
    expect_error(
      markov_update_po_native(
        previous,
        scalar,
        crossed,
        as.integer(non_absorb),
        as.integer(absorb),
        time = 2
      ),
      "cross for 4 of the 4 profiles at time 2",
      class = "mostr_crossed_thresholds"
    )
  }
})

test_that("native posterior-draw updates match the reference per draw", {
  set.seed(4104)
  draws <- 3L
  n <- 4L
  K <- 5L
  previous <- array(0, dim = c(draws, n, K))
  for (draw in seq_len(draws)) {
    previous[draw, , ] <- native_reference_previous(n, K)
  }

  for (absorb in native_absorbing_sets) {
    non_absorb <- setdiff(seq_len(K), absorb)
    transition <- array(0, dim = c(draws, n * length(non_absorb), K))
    for (draw in seq_len(draws)) {
      transition[draw, , ] <- reference_random_stochastic(
        n * length(non_absorb),
        K
      )
    }
    actual <- markov_update_draws_native(
      previous,
      transition,
      non_absorb,
      absorb
    )

    expected <- array(0, dim = c(draws, n, K))
    for (draw in seq_len(draws)) {
      blocks <- reference_origin_blocks(transition[draw, , ], non_absorb, n)
      expected[draw, , ] <- reference_markov_step(
        previous[draw, , ],
        blocks,
        absorb
      )
    }
    expect_equal(actual, expected, tolerance = 1e-14)
  }
})

# Joint probabilities over (state two visits ago, previous state). Mass sits
# only on histories that can occur: an absorbing earlier state implies the same
# absorbing previous state.
native_reference_joint <- function(n, K, absorb) {
  joint <- array(stats::rexp(n * K * K), dim = c(n, K, K))
  for (h in absorb) {
    joint[, h, setdiff(seq_len(K), h)] <- 0
  }
  joint[1L, , ] <- 0
  joint[1L, 3L, 3L] <- 1
  joint / apply(joint, 1L, sum)
}

native_reference_pairs <- function(K, absorb) {
  non_absorb <- setdiff(seq_len(K), absorb)
  pairs <- expand.grid(older = non_absorb, current = non_absorb)
  pairs[sample.int(nrow(pairs)), , drop = FALSE]
}

test_that("native second-order updates match the joint-history reference", {
  set.seed(4105)
  n <- 3L
  K <- 5L

  for (absorb in native_absorbing_sets) {
    joint <- native_reference_joint(n, K, absorb)
    pairs <- native_reference_pairs(K, absorb)
    transition <- reference_random_stochastic(n * nrow(pairs), K)
    transition_for <- function(h, j) {
      pair <- which(pairs$older == h & pairs$current == j)
      transition[(pair - 1L) * n + seq_len(n), , drop = FALSE]
    }
    expected <- reference_second_order_step(joint, transition_for, absorb)

    actual <- markov_update_second_order_native(
      joint,
      transition,
      pairs$older,
      pairs$current,
      absorb
    )
    expect_equal(actual, expected, tolerance = 1e-14)

    # The R engine accumulates pair chunks and absorbing mass
    # separately.
    chunks <- split(seq_len(nrow(pairs)), rep(1:2, length.out = nrow(pairs)))
    chunked <- Reduce(`+`, lapply(chunks, function(chunk) {
      rows <- unlist(lapply(chunk, function(pair) (pair - 1L) * n + seq_len(n)))
      markov_update_second_order_native(
        joint,
        transition[rows, , drop = FALSE],
        pairs$older[chunk],
        pairs$current[chunk],
        integer()
      )
    }))
    if (length(absorb) > 0L) {
      chunked <- chunked + markov_update_second_order_native(
        joint,
        matrix(numeric(), nrow = 0L, ncol = K),
        integer(),
        integer(),
        absorb
      )
    }
    expect_equal(chunked, expected, tolerance = 1e-14)
  }
})

test_that("native second-order PO updates match the joint-history reference", {
  set.seed(4106)
  n <- 3L
  K <- 5L

  for (absorb in native_absorbing_sets) {
    joint <- native_reference_joint(n, K, absorb)
    pairs <- native_reference_pairs(K, absorb)
    scalar <- matrix(stats::rnorm(n * nrow(pairs), sd = 1.5), nrow = n)
    cutpoints <- c(2, 0.5, -0.4, -2.2)
    transition_for <- function(h, j) {
      pair <- which(pairs$older == h & pairs$current == j)
      reference_checked_probabilities(outer(scalar[, pair], cutpoints, "+"))
    }
    expected <- reference_second_order_step(joint, transition_for, absorb)

    actual <- markov_update_second_order_po_native(
      joint,
      scalar,
      cutpoints,
      pairs$older,
      pairs$current,
      absorb
    )
    expect_equal(actual, expected, tolerance = 1e-12)

    # Cut-points out of order cross for every patient.
    crossed <- c(1, -1, 0.5, -2)
    expect_error(
      reference_checked_probabilities(outer(scalar[, 1L], crossed, "+")),
      "cross"
    )
    marked <- cpp_markov_update_second_order_po(
      joint,
      scalar,
      crossed,
      n,
      K,
      as.integer(pairs$older),
      as.integer(pairs$current),
      as.integer(absorb)
    )
    expect_identical(attr(marked, "crossed_profiles"), seq_len(n))
    expect_error(
      markov_update_second_order_po_native(
        joint,
        scalar,
        crossed,
        pairs$older,
        pairs$current,
        absorb,
        time = 4
      ),
      "cross for 3 of the 3 profiles at time 4",
      class = "mostr_crossed_thresholds"
    )
  }
})

test_that("native BLRM probabilities match the reference and keep crossed rows raw", {
  set.seed(4107)
  draws <- 3L
  n <- 4L
  thresholds <- 3L
  base_eta <- matrix(stats::rnorm(draws * n), nrow = draws)
  intercepts <- t(replicate(draws, sort(stats::rnorm(thresholds), TRUE)))
  threshold_eta <- matrix(stats::rnorm(draws * n, sd = 2), nrow = draws)
  threshold_scale <- c(0, 1.5, -1.5)

  reference <- function(partial) {
    out <- array(0, dim = c(draws, n, thresholds + 1L))
    for (draw in seq_len(draws)) {
      eta <- outer(base_eta[draw, ], intercepts[draw, ], "+")
      if (partial) {
        eta <- eta + outer(threshold_eta[draw, ], threshold_scale)
      }
      out[draw, , ] <- reference_posterior_probabilities(eta)
    }
    out
  }

  expect_equal(
    blrm_probabilities_native(base_eta, intercepts),
    reference(FALSE),
    tolerance = 1e-12
  )
  partial <- blrm_probabilities_native(
    base_eta,
    intercepts,
    threshold_eta,
    threshold_scale
  )
  expected <- reference(TRUE)
  # Some cut-points cross here, and those rows keep their negative values for
  # the caller to report.
  expect_lt(min(expected), -1e-10)
  expect_equal(partial, expected, tolerance = 1e-12)
  expect_error(
    check_crossed_thresholds(partial, n_profiles = n, time = 1),
    "posterior draws",
    class = "mostr_crossed_thresholds"
  )
})

test_that("native posterior reduction matches unequal unordered groups", {
  set.seed(4109)
  values <- array(stats::runif(2L * 7L * 3L * 4L), dim = c(2L, 7L, 3L, 4L))
  groups <- c(3L, 1L, 2L, 1L, 3L, 3L, 2L)

  expect_equal(
    reduce_sops_draw_array(values, groups, 3L),
    reference_group_means(values, groups, 3L),
    tolerance = 1e-14
  )
})

test_that("native categorical sampling is inverse-CDF sampling on R uniforms", {
  set.seed(4110)
  random_rows <- reference_random_stochastic(400L, 5L)
  weights <- rbind(
    c(0, 2, 0, 5, 1),
    c(3, 0, 0, 0, 0),
    c(0, 0, 0, 0, 7),
    c(1, 1, 1, 1, 1),
    c(0.25, 0, 0.5, 0, 0.25)
  )
  probabilities <- rbind(random_rows, weights[rep(1:5, 40L), ])

  set.seed(99)
  actual <- sample_categorical_rows(probabilities)
  next_actual <- stats::runif(3L)
  set.seed(99)
  uniforms <- stats::runif(nrow(probabilities))
  next_expected <- stats::runif(3L)
  expected <- reference_sample_categorical_rows(probabilities, uniforms)

  expect_identical(actual, expected)
  expect_identical(next_actual, next_expected)
  expect_false(any(probabilities[cbind(seq_along(actual), actual)] == 0))

  single <- matrix(c(0.3, 2, 1), ncol = 1L)
  set.seed(5)
  expect_identical(sample_categorical_rows(single), rep(1L, 3L))
})

native_plan_fits <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      data <- make_test_data(n_patients = 80, follow_up_time = 6, seed = 7351)
      value <<- list(
        baseline = data[!duplicated(data$id), , drop = FALSE][1:4, ],
        models = list(
          orm = suppressWarnings(mostm(
            y ~ time + tx + yprev,
            data = data,
            id_var = "id",
            absorb = 6
          )),
          vglm_po = make_test_model(data)$vglm_fit,
          vglm_partial_po = suppressWarnings(mostm(
            ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
            family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
            data = data,
            id_var = "id",
            absorb = 6,
            engine = "vglm"
          ))$vglm_fit
        )
      )
    }
    value
  }
})

test_that("compiled first-order SOPs match the reference recursion", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  fits <- native_plan_fits()

  for (model in fits$models) {
    for (absorb in list(6, c(1, 6), c(3, 6))) {
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = fits$baseline,
        times = 1:5,
        settings = markov_settings(
          y_levels = 1:6,
          absorb = absorb
        )
      )
      actual <- run_sop_execution_plan(plan, get_effective_coefs(model))
      expected <- reference_first_order_plan_sops(
        plan,
        reference_backend_gamma(model, plan$components$col_names)
      )
      expect_equal(unname(actual), expected, tolerance = 1e-12)
    }
  }
})

native_second_order_data <- function(n_patients = 40L, seed = 4111) {
  set.seed(seed)
  data <- expand.grid(
    id = seq_len(n_patients),
    time = seq_len(6L),
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$tx <- rep(rep(0:1, length.out = n_patients), each = 6L)
  data$y <- factor(sample(1:3, nrow(data), replace = TRUE), levels = 1:3)
  data$yprev <- stats::ave(as.integer(data$y), data$id, FUN = function(value) {
    c(value[1L], value[-length(value)])
  })
  data$ypprev <- stats::ave(data$yprev, data$id, FUN = function(value) {
    c(value[1L], value[-length(value)])
  })
  data$yprev <- factor(data$yprev, levels = 1:3)
  data$ypprev <- factor(data$ypprev, levels = 1:3)
  data
}

test_that("compiled second-order SOPs match the reference recursion", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  data <- native_second_order_data()
  fit_data <- data[data$time > 2L, , drop = FALSE]
  models <- list(
    orm = suppressWarnings(mostm(
      y ~ tx + time + yprev + ypprev,
      data = fit_data,
      id_var = "id",
      first_followup_time = 3
    )),
    vglm_partial_po = suppressWarnings(mostm(
      ordered(y) ~ tx + time + yprev + ypprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
      data = fit_data,
      id_var = "id",
      first_followup_time = 3,
      engine = "vglm"
    ))$vglm_fit
  )
  # Every pair of starting states. With `absorb`, only the pairs outside the
  # absorbing states are valid starts; SOP calls refuse the others.
  baseline <- data[data$time == 1L, , drop = FALSE][1:9, , drop = FALSE]
  baseline$yprev <- factor(rep(1:3, times = 3L), levels = 1:3)
  baseline$ypprev <- factor(rep(1:3, each = 3L), levels = 1:3)

  for (model in models) {
    for (absorb in list(NULL, 3, 1)) {
      valid <- !baseline$yprev %in% absorb & !baseline$ypprev %in% absorb
      plan <- compile_sop_execution_plan(
        model,
        baseline[valid, , drop = FALSE],
        times = 1:5,
        settings = markov_settings(
          y_levels = 1:3,
          absorb = absorb,
          p2_var = "ypprev"
        )
      )
      actual <- run_sop_execution_plan(plan, get_effective_coefs(model))
      expected <- reference_second_order_plan_sops(
        plan,
        reference_backend_gamma(model, plan$components$col_names)
      )
      expect_identical(plan$recursion_order, 2L)
      expect_equal(unname(actual), expected, tolerance = 1e-12)
    }
  }
})
