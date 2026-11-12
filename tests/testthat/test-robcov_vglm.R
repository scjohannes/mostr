# Unit Tests: Robust Covariance Estimation for vglm
#
# These tests verify that:
# 1. robcov_vglm produces valid sandwich covariance estimates
# 2. Results match rms::robcov for equivalent models
# 3. Clustering works correctly
# 4. Score computations have correct properties
# 5. Model information is preserved
# 6. Z-statistics and p-values are computed correctly
#
# The inputs robcov_vglm() refuses are rows V01 to V04 of test-refusals.R.

vglm_candidate_at <- function(fit, beta) {
  X_vlm <- stats::model.matrix(fit, type = "vlm")
  n <- stats::nobs(fit, type = "lm")
  M <- VGAM::npred(fit)
  fitted_eta <- matrix(
    drop(X_vlm %*% stats::coef(fit)),
    nrow = n,
    ncol = M,
    byrow = TRUE
  )
  eta <- matrix(
    drop(X_vlm %*% beta),
    nrow = n,
    ncol = M,
    byrow = TRUE
  ) +
    fit@predictors -
    fitted_eta
  candidate <- fit
  candidate@coefficients <- beta
  candidate@predictors <- eta
  candidate@fitted.values <- candidate@family@linkinv(
    eta = eta,
    extra = candidate@extra
  )
  candidate
}

finite_difference_loglik_scores <- function(fit) {
  central_difference_jacobian(
    function(beta) {
      candidate <- vglm_candidate_at(fit, beta)
      as.numeric(stats::logLik(candidate, summation = FALSE))
    },
    stats::coef(fit),
    step = 1e-5
  )
}

finite_difference_score_jacobian <- function(fit) {
  evaluator <- mostr:::vglm_score_evaluator(fit)
  central_difference_jacobian(evaluator$evaluate, stats::coef(fit), step = 1e-5)
}

# One proportional-odds fit with three predictors and four outcome levels
# (three intercepts and three slopes) for the checks of the returned object.
# `x2` has a strong effect, so its robust p-value is far below 1e-5. Made once
# and shared by the tests that use it.
robcov_contract_cache <- new.env(parent = emptyenv())

robcov_contract_fixture <- function() {
  if (is.null(robcov_contract_cache$fit)) {
    test_data <- withr::with_seed(707, {
      n <- 400
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      x3 <- rnorm(n)
      lp <- 0.3 * x1 + 1.5 * x2 - 0.2 * x3
      probs <- cbind(
        plogis(-1 - lp),
        plogis(0 - lp) - plogis(-1 - lp),
        plogis(1 - lp) - plogis(0 - lp),
        1 - plogis(1 - lp)
      )
      y <- apply(probs, 1, function(p) sample(1:4, 1, prob = pmax(p, 0.001)))
      data.frame(y = ordered(y), x1 = x1, x2 = x2, x3 = x3)
    })
    robcov_contract_cache$n <- nrow(test_data)
    robcov_contract_cache$fit <- VGAM::vglm(
      y ~ x1 + x2 + x3,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = test_data
    )
  }
  list(fit = robcov_contract_cache$fit, n = robcov_contract_cache$n)
}

describe("robcov_vglm()", {
  it("produces valid covariance estimates", {
    skip_on_cran()
    skip_if_not_installed("VGAM")

    fixture <- robcov_contract_fixture()
    m <- fixture$fit
    n <- fixture$n

    # Compute robust covariance
    result <- robcov_vglm(m)

    # Check output structure
    expect_s3_class(result, "robcov_vglm")

    # Check variance and SE components
    expect_contains(
      names(result),
      c(
        "var",
        "se",
        "scores",
        "influence_beta_obs",
        "bread",
        "bread_type",
        "bread_diagnostics",
        "meat"
      )
    )

    # Check new components
    expect_contains(
      names(result),
      c("coefficients", "z", "pvalues", "original_se")
    )

    # Check model information components
    expect_contains(
      names(result),
      c("family", "constraints", "control", "call", "original_call")
    )

    # Check dimensions: 3 intercepts and 3 slopes
    p <- length(coef(m))
    expect_equal(p, 6)
    expect_equal(dim(result$var), c(p, p))
    expect_length(result$se, p)
    expect_length(result$coefficients, p)
    expect_length(result$z, p)
    expect_length(result$pvalues, p)
    expect_equal(dim(result$scores), c(n, p))
    expect_equal(dim(result$influence_beta_obs), c(n, p))
    expect_equal(dim(result$bread), c(p, p))
    expect_equal(dim(result$meat), c(p, p))
    expect_equal(result$bread_type, "observed")
    expect_equal(result$type, "HC0")
    expect_identical(result$cadjust, FALSE)

    result_vglm_bread <- robcov_vglm(m, bread = "vglm")
    expect_equal(result_vglm_bread$bread, vcov(m))
    expect_equal(result_vglm_bread$bread_type, "vglm")
    expected_vglm_var <- crossprod(
      result_vglm_bread$scores %*% vcov(m)
    )
    expected_vglm_var <- (expected_vglm_var + t(expected_vglm_var)) / 2
    dimnames(expected_vglm_var) <- dimnames(vcov(m))
    expect_equal(result_vglm_bread$var, expected_vglm_var)

    expected_var <- crossprod(result$scores %*% result$bread)
    expected_var <- (expected_var + t(expected_var)) / 2
    dimnames(expected_var) <- dimnames(result$var)
    expect_equal(result$var, expected_var, tolerance = 1e-15)

    # The stored influence contributions use the selected bread.
    expect_equal(
      unname(result$influence_beta_obs),
      unname(result$scores %*% result$bread),
      tolerance = 1e-10
    )

    # Check variance matrix is symmetric
    expect_equal(result$var, t(result$var), tolerance = 1e-15)

    # Check variance matrix is positive semi-definite
    eigenvalues <- eigen(result$var, only.values = TRUE)$values
    expect_all_true(eigenvalues >= -1e-10)

    # Check standard errors are positive
    expect_all_true(result$se > 0)

    # Check n matches
    expect_equal(result$n, n)

    result_row_cluster <- robcov_vglm(m, cluster = seq_len(n), adjust = FALSE)
    expect_equal(result_row_cluster$var, result$var, tolerance = 1e-15)

    # coef() and vcov() return the stored estimates.
    expect_equal(coef(result), result$coefficients)
    expect_length(coef(result), p)
    expect_equal(vcov(result), result$var)
    expect_equal(dim(vcov(result)), c(p, p))
    expect_equal(vcov(result), t(vcov(result)))
  })

  it("computes z-statistics and p-values and keeps the fit's information", {
    skip_on_cran()
    skip_if_not_installed("VGAM")

    m <- robcov_contract_fixture()$fit
    result <- robcov_vglm(m)

    # Check z = coef / se
    expected_z <- result$coefficients / result$se
    expect_equal(result$z, expected_z, tolerance = 1e-10)

    # Check p-values = 2 * pnorm(-abs(z))
    expected_p <- 2 * pnorm(-abs(result$z))
    expect_equal(result$pvalues, expected_p, tolerance = 1e-10)

    # Check that strong effect has small p-value
    expect_lt(result$pvalues["x2"], 0.00001)

    # Check all p-values are in [0, 1]
    expect_all_true(result$pvalues >= 0 & result$pvalues <= 1)

    # Check parameter names
    expect_contains(names(result$se), c("x1", "x2", "x3"))

    # Names should be consistent across all vectors
    expect_equal(names(result$coefficients), names(result$se))
    expect_equal(names(result$coefficients), names(result$z))
    expect_equal(names(result$coefficients), names(result$pvalues))

    # Check model information is preserved
    expect_equal(result$coefficients, coef(m))
    expect_equal(coef(result), coef(m))
    expect_equal(result$fitted.values, m@fitted.values)
    expect_equal(result$residuals, m@residuals)
    expect_equal(result$df.residual, m@df.residual)
    expect_equal(result$df.total, m@df.total)

    # Check family is preserved
    expect_equal(class(result$family), class(m@family))

    # Check constraints are preserved
    expect_equal(result$constraints, m@constraints)

    # Check control is preserved
    expect_equal(result$control, m@control)

    # Check original vcov matches
    expect_equal(result$original_vcov, vcov(m))
    expect_equal(result$original_se, sqrt(diag(vcov(m))))
  })

  it("matches rms::robcov for equivalent binary and ordinal likelihoods", {
    skip_on_cran()
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    # Each data set is fitted once with rms::orm() and once with VGAM::vglm().
    fit_both <- function(formula, data, family) {
      list(
        orm = rms::orm(formula, data = data, x = TRUE, y = TRUE),
        vglm = VGAM::vglm(formula, family = family, data = data)
      )
    }

    binary <- withr::with_seed(456, {
      n <- 500
      x <- rnorm(n)
      lp <- -0.5 + 0.8 * x
      y <- rbinom(n, 1, plogis(lp))
      fit_both(y ~ x, data.frame(y = y, x = x), VGAM::binomialff)
    })

    n_clusters <- 50
    binary_clusters <- rep(1:n_clusters, each = 10)
    binary_clustered <- withr::with_seed(789, {
      n <- length(binary_clusters)
      x <- rnorm(n)
      cluster_effect <- rep(rnorm(n_clusters, sd = 0.3), each = 10)
      lp <- -0.5 + 0.8 * x + cluster_effect
      y <- rbinom(n, 1, plogis(lp))
      fit_both(y ~ x, data.frame(y = y, x = x), VGAM::binomialff)
    })

    ordinal <- withr::with_seed(101, {
      n <- 500
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      lp <- 0.5 * x1 + 0.3 * x2
      cuts <- c(-1, 0, 1)
      probs <- matrix(NA, nrow = n, ncol = 4)
      for (j in 1:4) {
        if (j == 1) {
          probs[, j] <- plogis(cuts[1] - lp)
        } else if (j == 4) {
          probs[, j] <- 1 - plogis(cuts[3] - lp)
        } else {
          probs[, j] <- plogis(cuts[j] - lp) - plogis(cuts[j - 1] - lp)
        }
      }
      y <- apply(probs, 1, function(p) sample(1:4, 1, prob = pmax(p, 0.001)))
      fit_both(
        y ~ x1 + x2,
        data.frame(y = ordered(y), x1 = x1, x2 = x2),
        VGAM::cumulative(reverse = TRUE, parallel = TRUE)
      )
    })

    # rms::robcov() makes no small-sample correction for clusters, so the
    # clustered rows turn it off, once with the legacy `adjust` spelling.
    cases <- list(
      list(label = "binary", fits = binary, cluster = NULL, args = list(), tolerance = 1e-6),
      list(
        label = "binary, clustered",
        fits = binary_clustered,
        cluster = binary_clusters,
        args = list(adjust = FALSE),
        tolerance = 1e-6,
        n_clusters = n_clusters
      ),
      list(label = "ordinal", fits = ordinal, cluster = NULL, args = list(), tolerance = 1e-5),
      list(
        label = "ordinal, clustered",
        fits = ordinal,
        cluster = rep(seq_len(100), each = 5),
        args = list(cadjust = FALSE),
        tolerance = 1e-5
      )
    )

    for (case in cases) {
      orm_robust <- if (is.null(case$cluster)) {
        rms::robcov(case$fits$orm)
      } else {
        rms::robcov(case$fits$orm, cluster = case$cluster)
      }
      vglm_robust <- do.call(
        robcov_vglm,
        c(list(case$fits$vglm, cluster = case$cluster), case$args)
      )
      se_orm <- sqrt(diag(orm_robust$var))

      # Use unname() because names differ between packages
      expect_equal(
        unname(vglm_robust$var),
        unname(orm_robust$var),
        tolerance = case$tolerance,
        info = case$label
      )
      expect_equal(
        unname(vglm_robust$se),
        unname(se_orm),
        tolerance = case$tolerance,
        info = case$label
      )
      expect_lt(max(abs(vglm_robust$se / se_orm - 1)), case$tolerance)
      if (!is.null(case$n_clusters)) {
        expect_equal(vglm_robust$n_clusters, case$n_clusters, info = case$label)
      }
    }
  })

  it("rejects non-converged fits and malformed covariance components", {
    skip_on_cran()
    skip_if_not_installed("VGAM")

    withr::local_seed(405)
    n <- 120
    data <- data.frame(
      y = stats::rbinom(n, 1, 0.5),
      x = stats::rnorm(n)
    )
    fit <- VGAM::vglm(
      y ~ x,
      family = VGAM::binomialff,
      data = data
    )

    nonconverged <- fit
    nonconverged@iter <- nonconverged@control$maxit
    expect_snapshot(error = TRUE, {
      robcov_vglm(nonconverged)
    })

    testthat::local_mocked_bindings(
      compute_scores_vglm = function(fit) {
        scores <- matrix(
          0,
          nrow = stats::nobs(fit, type = "lm"),
          ncol = length(stats::coef(fit))
        )
        colnames(scores) <- rev(names(stats::coef(fit)))
        scores
      },
      .package = "mostr"
    )
    expect_snapshot(error = TRUE, {
      robcov_vglm(fit)
    })
  })
})

describe("Score properties", {
  it("sum to approximately zero", {
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(202)
    n <- 300
    x <- rnorm(n)
    lp <- 0.5 * x
    probs <- cbind(
      plogis(-1 - lp),
      plogis(0 - lp) - plogis(-1 - lp),
      1 - plogis(0 - lp)
    )
    y <- apply(probs, 1, function(p) sample(1:3, 1, prob = pmax(p, 0.001)))

    test_data <- data.frame(y = ordered(y), x = x)

    # Fit vglm
    m <- VGAM::vglm(
      y ~ x,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = test_data
    )

    # Compute robust covariance
    result <- robcov_vglm(m)

    # Score columns should sum to approximately zero (property of MLE)
    score_sums <- colSums(result$scores)
    expect_all_true(abs(score_sums) < 1e-5)
  })

  it("matches the loop-equivalent score calculation", {
    skip_if_not_installed("VGAM")

    withr::local_seed(203)
    n <- 160
    test_data <- data.frame(
      y = ordered(sample(1:4, n, replace = TRUE)),
      x1 = rnorm(n),
      x2 = rnorm(n)
    )
    m <- VGAM::vglm(
      y ~ x1 + x2,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = test_data
    )

    scores <- compute_scores_vglm(m)
    deriv_eta <- VGAM::weights(
      m,
      type = "working",
      deriv.arg = TRUE
    )$deriv
    X_vlm <- stats::model.matrix(m, type = "vlm")
    M <- VGAM::npred(m)

    scores_loop <- matrix(0, nrow = n, ncol = length(coef(m)))
    for (i in seq_len(n)) {
      vlm_rows <- ((i - 1) * M + 1):(i * M)
      scores_loop[i, ] <- as.vector(crossprod(
        X_vlm[vlm_rows, , drop = FALSE],
        deriv_eta[i, ]
      ))
    }
    colnames(scores_loop) <- names(coef(m))
    rownames(scores_loop) <- rownames(scores)

    expect_equal(scores, scores_loop)
  })
})

test_that("compare_se_orm_vglm compares model-based and robust SEs", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(8128)
  n <- 80
  data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x = rnorm(n),
    cluster = rep(seq_len(20), each = 4)
  )

  orm_fit <- rms::orm(y ~ x, data = data, x = TRUE, y = TRUE)
  vglm_fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)

  unclustered <- compare_se_orm_vglm(orm_fit, vglm_fit)
  expect_warning(
    clustered <- compare_se_orm_vglm(orm_fit, vglm_fit, cluster = data$cluster),
    "Cluster-robust p-values",
    fixed = TRUE
  )

  expect_named(
    unclustered,
    c(
      "parameter",
      "se_orm_model",
      "se_vglm_model",
      "se_orm_robust",
      "se_vglm_robust",
      "ratio_model",
      "ratio_robust"
    )
  )
  expect_equal(nrow(clustered), nrow(unclustered))

  # Both packages fit the same logistic likelihood, so their robust standard
  # errors agree with rms::robcov() and robcov_vglm() (clustered without a
  # small-sample correction) and their ratios are 1.
  expected <- list(
    unclustered = list(
      comparison = unclustered,
      orm = rms::robcov(orm_fit),
      vglm = robcov_vglm(vglm_fit)
    ),
    clustered = list(
      comparison = clustered,
      orm = rms::robcov(orm_fit, cluster = data$cluster),
      vglm = suppressWarnings(
        robcov_vglm(vglm_fit, cluster = data$cluster, cadjust = FALSE)
      )
    )
  )
  for (label in names(expected)) {
    case <- expected[[label]]
    expect_equal(
      unname(case$comparison$se_orm_robust),
      unname(sqrt(diag(case$orm$var))),
      tolerance = 1e-6,
      info = label
    )
    expect_equal(
      unname(case$comparison$se_vglm_robust),
      unname(case$vglm$se),
      info = label
    )
    expect_equal(case$comparison$ratio_robust, c(1, 1), tolerance = 1e-5, info = label)
    expect_true(all(is.finite(case$comparison$ratio_model)), info = label)
  }

  # An ordinal orm fit has more intercepts than a logistic vglm fit.
  ordinal_data <- data.frame(
    y_ord = ordered(rep(1:3, length.out = 90)),
    y_bin = rep(0:1, length.out = 90),
    x = withr::with_seed(8129, rnorm(90))
  )
  ordinal_orm <- rms::orm(y_ord ~ x, data = ordinal_data, x = TRUE, y = TRUE)
  binary_vglm <- VGAM::vglm(y_bin ~ x, family = VGAM::binomialff, data = ordinal_data)
  expect_warning(
    comparison <- compare_se_orm_vglm(ordinal_orm, binary_vglm),
    "Number of intercepts differs",
    fixed = TRUE
  )
  expect_s3_class(comparison, "data.frame")
})

describe("Small-sample adjustment", {
  it("separates HC1 and cluster corrections", {
    skip_on_cran()
    skip_if_not_installed("VGAM")

    # Generate clustered data
    withr::local_seed(303)
    n_clusters <- 20
    n_per_cluster <- 15
    n <- n_clusters * n_per_cluster
    cluster <- rep(1:n_clusters, each = n_per_cluster)

    x <- rnorm(n)
    y <- rbinom(n, 1, plogis(0.5 * x))

    test_data <- data.frame(y = y, x = x, cluster = cluster)

    # Fit vglm
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

    # Compute with and without adjustment
    expect_warning(
      result_noadj <- robcov_vglm(m, cluster = cluster, cadjust = FALSE),
      "fewer than 30 clusters"
    )
    expect_warning(
      result_adj <- robcov_vglm(m, cluster = cluster, cadjust = TRUE),
      "fewer than 30 clusters"
    )
    expect_warning(
      result_hc1 <- robcov_vglm(
        m,
        cluster = cluster,
        type = "HC1",
        cadjust = FALSE
      ),
      "fewer than 30 clusters"
    )
    expect_warning(
      result_both <- robcov_vglm(
        m,
        cluster = cluster,
        type = "HC1",
        cadjust = TRUE
      ),
      "fewer than 30 clusters"
    )

    # Adjustment factor is G/(G-1) = 20/19 ≈ 1.0526
    expected_ratio <- sqrt(n_clusters / (n_clusters - 1))

    # SEs with adjustment should be larger by this factor
    se_ratio <- result_adj$se / result_noadj$se
    expect_equal(
      unname(se_ratio),
      rep(expected_ratio, length(se_ratio)),
      tolerance = 1e-10
    )
    expect_equal(
      result_adj$var,
      suppressWarnings(robcov_vglm(m, cluster = cluster))$var
    )
    hc1_factor <- (n - 1) / (n - length(coef(m)))
    expect_equal(result_hc1$var, result_noadj$var * hc1_factor)
    expect_equal(
      result_both$var,
      result_noadj$var * hc1_factor * n_clusters / (n_clusters - 1)
    )
    expect_equal(result_both$adjustment_factor, hc1_factor * 20 / 19)

    expect_warning(
      result_legacy <- robcov_vglm(m, cluster = cluster, adjust = TRUE),
      "fewer than 30 clusters"
    )
    expect_equal(result_legacy$var, result_adj$var)

    expect_snapshot(error = TRUE, {
      robcov_vglm(
        m,
        cluster = cluster,
        adjust = TRUE,
        cadjust = TRUE
      )
    })
  })
})

test_that("analytic scores and observed bread match independent differences", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  withr::local_seed(204)
  n <- 240
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  eta <- 0.4 * x - 0.2 * z
  cuts <- c(-1, 0, 1)
  cumulative <- vapply(
    cuts,
    function(cut) {
      stats::plogis(cut - eta)
    },
    numeric(n)
  )
  probabilities <- cbind(
    cumulative[, 1],
    cumulative[, 2] - cumulative[, 1],
    cumulative[, 3] - cumulative[, 2],
    1 - cumulative[, 3]
  )
  y <- apply(probabilities, 1, function(probability) {
    sample(seq_len(4), 1, prob = probability)
  })
  data <- data.frame(
    y = ordered(y),
    x = x,
    z = z,
    offset = 0.1 * z,
    weight = sample(1:3, n, replace = TRUE)
  )
  fit <- VGAM::vglm(
    y ~ x + offset(offset),
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    weights = weight,
    x = TRUE,
    y = TRUE
  )

  scores <- compute_scores_vglm(fit)
  numerical_scores <- finite_difference_loglik_scores(fit)
  expect_equal(scores, numerical_scores, tolerance = 1e-6, ignore_attr = TRUE)

  numerical_jacobian <- finite_difference_score_jacobian(fit)
  numerical_information <- -(numerical_jacobian + t(numerical_jacobian)) / 2
  numerical_bread <- solve(numerical_information)
  observed <- robcov_vglm(fit)
  expect_equal(
    unname(observed$bread),
    unname(numerical_bread),
    tolerance = 1e-6
  )
})

test_that("logistic scores and observed bread match independent differences with weights or an offset", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  # Non-integer case weights, and an offset given as the `offset` argument
  # rather than in the formula.
  weighted_data <- withr::with_seed(9002, {
    n <- 100
    x <- rnorm(n)
    y <- rbinom(n, 1, 0.5)
    weights <- runif(n, 0.5, 2)
    data.frame(y = y, x = x, w = weights)
  })
  offset_data <- withr::with_seed(9003, {
    n <- 100
    x <- rnorm(n)
    off <- rep(0.5, n) # fixed offset
    y <- rbinom(n, 1, plogis(0.5 * x + off))
    data.frame(y = y, x = x, off = off)
  })
  fits <- list(
    weights = VGAM::vglm(
      y ~ x,
      family = VGAM::binomialff,
      data = weighted_data,
      weights = w
    ),
    offset = VGAM::vglm(
      y ~ x,
      family = VGAM::binomialff,
      data = offset_data,
      offset = off
    )
  )

  for (label in names(fits)) {
    fit <- fits[[label]]
    res <- robcov_vglm(fit)

    expect_s3_class(res, "robcov_vglm")
    expect_all_true(res$se > 0)
    expect_equal(length(res$se), 2, info = label)
    # The offset acts like a fixed intercept shift, not a parameter.
    expect_equal(names(coef(res)), c("(Intercept)", "x"), info = label)

    expect_equal(
      compute_scores_vglm(fit),
      finite_difference_loglik_scores(fit),
      tolerance = 1e-6,
      ignore_attr = TRUE,
      info = label
    )
    numerical_jacobian <- finite_difference_score_jacobian(fit)
    numerical_bread <- solve(-(numerical_jacobian + t(numerical_jacobian)) / 2)
    expect_equal(
      unname(res$bread),
      unname(numerical_bread),
      tolerance = 1e-6,
      info = label
    )
  }
})

test_that("observed bread supports vector-valued inverse links", {
  skip_if_not_installed("VGAM")

  withr::local_seed(206)
  n <- 180
  data <- data.frame(x = stats::rnorm(n))
  data$y <- 0.2 + 0.4 * data$x + stats::rnorm(n)
  fit <- VGAM::vglm(
    y ~ x,
    family = VGAM::uninormal,
    data = data
  )

  inverse_link <- fit@family@linkinv(
    eta = fit@predictors,
    extra = fit@extra
  )
  expect_null(dim(inverse_link))

  evaluator <- mostr:::vglm_score_evaluator(fit)
  evaluated_scores <- evaluator$evaluate(stats::coef(fit))
  expected_scores <- colSums(compute_scores_vglm(fit))
  expect_equal(evaluated_scores, expected_scores, tolerance = 1e-12)

  robust <- robcov_vglm(fit)
  expect_s3_class(robust, "robcov_vglm")
  expect_equal(robust$bread_type, "observed")
  expect_equal(all(is.finite(robust$bread)), TRUE)
  expect_equal(all(is.finite(robust$var)), TRUE)
})

test_that("observed bread supports proportional-odds deviation geometries", {
  skip_on_cran()
  skip_if_not_installed("VGAM")

  withr::local_seed(205)
  n <- 500
  x <- stats::runif(n, -1, 1)
  z <- stats::rnorm(n)
  group <- factor(sample(c("a", "b"), n, replace = TRUE))
  eta <- 0.5 * x - 0.25 * z + 0.1 * (group == "b")
  cuts <- c(-1.2, -0.3, 0.5, 1.3)
  cumulative <- vapply(
    cuts,
    function(cut) {
      stats::plogis(cut - eta)
    },
    numeric(n)
  )
  probabilities <- cbind(
    cumulative[, 1],
    cumulative[, 2] - cumulative[, 1],
    cumulative[, 3] - cumulative[, 2],
    cumulative[, 4] - cumulative[, 3],
    1 - cumulative[, 4]
  )
  y <- apply(probabilities, 1, function(probability) {
    sample(seq_len(5), 1, prob = probability)
  })
  data <- data.frame(y = ordered(y), x = x, z = z, group = group)

  full_nonparallel <- VGAM::vglm(
    y ~ x + z,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
    data = data,
    x = TRUE,
    y = TRUE
  )
  custom_constraints <- VGAM::constraints(full_nonparallel)
  n_thresholds <- VGAM::npred(full_nonparallel)
  custom_constraints[["x"]] <- cbind(
    PO = rep(1, n_thresholds),
    linear = seq_len(n_thresholds) - mean(seq_len(n_thresholds)),
    terminal = c(rep(0, n_thresholds - 1), 1)
  )
  custom_constraints[["z"]] <- cbind(PO = rep(1, n_thresholds))

  # A fully non-parallel fit without `reverse = TRUE`, on its own data: its
  # coefficients are (Intercept):1, (Intercept):2, x:1 and x:2.
  nonreversed_data <- withr::with_seed(9004, {
    x <- stats::rnorm(300)
    data.frame(y = ordered(sample(1:3, 300, replace = TRUE)), x = x)
  })

  fits <- list(
    proportional = VGAM::vglm(
      y ~ x + z + group,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      x = TRUE,
      y = TRUE
    ),
    selected_nonparallel = VGAM::vglm(
      y ~ x + z + group,
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ x),
      data = data,
      x = TRUE,
      y = TRUE
    ),
    fully_nonparallel = full_nonparallel,
    custom_deviations = VGAM::vglm(
      y ~ x + z,
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
      data = data,
      constraints = custom_constraints,
      x = TRUE,
      y = TRUE
    ),
    interactions = VGAM::vglm(
      y ~ x * z + group,
      family = VGAM::cumulative(
        reverse = TRUE,
        parallel = FALSE ~ x + x:z
      ),
      data = data,
      x = TRUE,
      y = TRUE
    ),
    nonreversed_nonparallel = VGAM::vglm(
      y ~ x,
      family = VGAM::cumulative(parallel = FALSE),
      data = nonreversed_data,
      x = TRUE,
      y = TRUE
    )
  )

  expect_length(stats::coef(fits$nonreversed_nonparallel), 4)

  for (label in names(fits)) {
    fit <- fits[[label]]
    n_fit <- stats::nobs(fit, type = "lm")
    robust <- robcov_vglm(
      fit,
      cluster = seq_len(n_fit),
      cadjust = FALSE
    )
    expect_equal(robust$bread_type, "observed", info = label)
    expect_equal(colnames(robust$scores), names(stats::coef(fit)), info = label)
    expect_equal(dim(robust$var), rep(length(stats::coef(fit)), 2), info = label)
    expect_equal(robust$var, t(robust$var), tolerance = 1e-12, info = label)
    expect_gt(min(eigen(robust$bread, symmetric = TRUE)$values), 0)
    expect_gte(min(eigen(robust$var, symmetric = TRUE)$values), -1e-10)
    expect_equal(all(is.finite(robust$se)), TRUE, info = label)
    expect_equal(all(robust$se > 0), TRUE, info = label)

    draws <- generate_score_bootstrap_draws(
      robust,
      baseline_data = data.frame(id = seq_len(n_fit)),
      id_var = "id",
      n_draws = 1
    )
    expect_equal(
      dim(draws$beta_draws),
      c(1L, length(stats::coef(fit))),
      info = label
    )
  }
})

describe("Score bootstrap compatibility", {
  it("uses the selected bread directly for one-step score updates", {
    skip_if_not_installed("VGAM")

    withr::local_seed(304)
    n_clusters <- 30
    n_per_cluster <- 4
    n <- n_clusters * n_per_cluster
    cluster <- rep(sprintf("id%02d", seq_len(n_clusters)), each = n_per_cluster)
    x <- rnorm(n)
    y <- rbinom(n, 1, plogis(-0.2 + 0.4 * x))
    test_data <- data.frame(y = y, x = x, id = cluster)

    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)
    robust <- robcov_vglm(m, cluster = test_data$id, adjust = FALSE)
    baseline_data <- data.frame(id = unique(test_data$id))

    cluster_ids <- unique(as.character(robust$cluster))
    scores_clustered <- rowsum(
      robust$scores,
      factor(as.character(robust$cluster), levels = cluster_ids),
      reorder = FALSE
    )

    withr::local_seed(305)
    w_cluster <- stats::rexp(length(cluster_ids), rate = 1)
    U_star <- as.vector(crossprod(w_cluster - 1, scores_clustered))
    expected_beta <- robust$coefficients + as.vector(robust$bread %*% U_star)

    withr::local_seed(305)
    draws <- generate_score_bootstrap_draws(
      robust,
      baseline_data = baseline_data,
      id_var = "id",
      n_draws = 1
    )

    expect_equal(unname(draws$beta_draws[1, ]), unname(expected_beta))
  })
})

describe("Utility methods", {
  it("print and summary methods work", {
    skip_on_cran()
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(606)
    n <- 100
    test_data <- data.frame(y = sample(0:1, n, replace = TRUE), x = rnorm(n))

    # Fit vglm
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

    # Compute robust covariance
    result <- robcov_vglm(m)

    # Test print method
    expect_snapshot(print(result))

    # Test summary method
    expect_snapshot(summary(result))

    # Test with clustering and no adjustment
    cluster <- rep(1:10, each = 10)
    expect_warning(
      result_cl <- robcov_vglm(m, cluster = cluster, adjust = FALSE),
      "fewer than 30 clusters"
    )
    expect_snapshot(summary(result_cl))

    # Test with clustering and explicit adjustment
    expect_warning(
      result_adj <- robcov_vglm(m, cluster = cluster, adjust = TRUE),
      "fewer than 30 clusters"
    )
    expect_snapshot(summary(result_adj))

    # Test print method returns invisibly
    expect_invisible(print(result))

    # Test summary returns coefficient table
    summary_result <- summary(result)
    expect_contains(names(summary_result), "coefficients")
    expect_s3_class(summary_result$coefficients, "data.frame")
    expect_contains(
      names(summary_result$coefficients),
      c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    )
  })
})

# Stress Tests & Edge Cases for robcov_vglm
#
# These tests focus on "breaking" the function with:
# 1. Missing data (NAs) handling
# 2. Few clusters, fewer clusters than coefficients, and empty factor levels
# 3. Ill-conditioned clustered fits

describe("robcov_vglm() stress tests", {
  it("handles missing data correctly (alignment with clusters)", {
    skip_if_not_installed("VGAM")

    withr::local_seed(9001)
    n <- 200
    x <- rnorm(n)
    y <- rbinom(n, 1, 0.5)
    x[1:10] <- NA
    y[11:20] <- NA

    test_data <- data.frame(y = y, x = x, id = 1:n)
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)
    omitted <- as.integer(unlist(m@na.action, use.names = FALSE))
    aligned_cluster <- test_data$id[-omitted]

    res_original <- robcov_vglm(m, cluster = test_data$id)
    res_aligned <- robcov_vglm(m, cluster = aligned_cluster)

    expect_s3_class(res_original, "robcov_vglm")
    expect_equal(res_original$n, 180) # 200 - 20 dropped
    expect_equal(res_original$cluster, aligned_cluster)
    expect_equal(res_original$var, res_aligned$var)
  })

  it("warns for few clusters and counts only the clusters present", {
    skip_on_cran()
    skip_if_not_installed("VGAM")

    # Five clusters; three clusters for four coefficients, so the clustered
    # middle part of the sandwich has rank below the number of coefficients;
    # and a cluster factor with an unused level "C", which is not a cluster.
    five_clusters <- withr::with_seed(90012, {
      n <- 100
      data.frame(
        y = rbinom(n, 1, 0.5),
        x = rnorm(n),
        cluster = rep(seq_len(5), each = 20)
      )
    })
    three_clusters <- withr::with_seed(9005, {
      n <- 60
      cluster <- rep(1:3, each = 20)
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      x3 <- rnorm(n)
      y <- rbinom(n, 1, 0.5)
      data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, cluster = cluster)
    })
    empty_level <- withr::with_seed(9006, {
      n <- 50
      x <- rnorm(n)
      y <- rbinom(n, 1, 0.5)
      cluster_char <- sample(c("A", "B"), n, replace = TRUE)
      data.frame(
        y = y,
        x = x,
        cluster = factor(cluster_char, levels = c("A", "B", "C"))
      )
    })
    cases <- list(
      list(label = "five clusters", formula = y ~ x, data = five_clusters, n_clusters = 5),
      list(
        label = "more coefficients than clusters",
        formula = y ~ x1 + x2 + x3,
        data = three_clusters,
        n_clusters = 3
      ),
      list(
        label = "a cluster factor with an empty level",
        formula = y ~ x,
        data = empty_level,
        n_clusters = 2
      )
    )

    for (case in cases) {
      m <- VGAM::vglm(case$formula, family = VGAM::binomialff, data = case$data)
      expect_warning(
        res <- robcov_vglm(m, cluster = case$data$cluster),
        "fewer than 30 clusters"
      )

      expect_s3_class(res, "robcov_vglm")
      expect_equal(res$n_clusters, case$n_clusters, info = case$label)
      # Clustered p-values stay normal z-tests.
      expect_equal(res$pvalues, 2 * pnorm(-abs(res$z)), info = case$label)
      expect_true(all(is.finite(res$se)), info = case$label)
      # The covariance has rank at most the number of clusters.
      expect_lte(qr(res$var)$rank, case$n_clusters)
    }
  })

  it("uses a stable PSD construction for ill-conditioned clustered fits", {
    skip_on_cran()
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    trial <- sim_markov(
      "actt2",
      n_patients = 40,
      follow_up_time = 6,
      seed = 1
    )
    data <- prepare_markov_data(trial, absorbing_state = 8)
    fit <- suppressWarnings(mostm(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id",
      engine = "vglm"
    ))
    robust <- robcov_vglm(fit$vglm_fit, cluster = data$id)

    cluster_scores <- rowsum(
      robust$scores,
      as.factor(robust$cluster),
      reorder = FALSE
    )
    expected <- crossprod(cluster_scores %*% robust$bread) *
      robust$adjustment_factor
    dimnames(expected) <- dimnames(robust$var)
    expect_equal(robust$var, expected, tolerance = 1e-12)
    expect_gte(
      min(eigen(robust$var, symmetric = TRUE, only.values = TRUE)$values),
      -1e-10
    )
  })
})
