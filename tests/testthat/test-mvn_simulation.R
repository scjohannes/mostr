# Unit Tests: MVN Simulation-Based Inference for SOPs
#
# These tests verify that:
# 1. set_coef() works correctly for vglm and orm models
# 2. get_vcov_robust() returns correct variance matrices
# 3. inferences(method = "mvn") uses a supplied covariance matrix
# 4. inferences() rejects invalid inputs

# A vglm fit for the cluster-robust covariance checks, made once per test run.
mvn_cluster_case <- function() {
  grid_fixture("mvn cluster case", function() {
    data <- make_test_data(n_patients = 30, seed = 444, follow_up_time = 15)
    list(data = data, model = make_test_model(data))
  })
}

describe("MVN Simulation-Based Inference for SOPs", {
  describe("set_coef()", {
    it("replaces coefficients of list, rms-style, robcov_vglm, vglm and orm models", {
      skip_on_cran()
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      coefs_of <- function(model) {
        if (isS4(model)) stats::coef(model) else model$coefficients
      }
      data <- mvn_cluster_case()$data
      m_vglm <- mvn_cluster_case()$model$vglm_fit
      m_orm <- rms::orm(
        y ~ time_lin + time_nlin_1 + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      )

      cases <- list(
        list(
          label = "default method",
          model = list(coefficients = c(a = 1, b = 2)),
          new = c(10, 20),
          expected = c(a = 10, b = 20)
        ),
        list(
          label = "lrm",
          model = structure(list(coefficients = c(a = 1)), class = "lrm"),
          new = 3,
          expected = c(a = 3)
        ),
        list(
          label = "rms",
          model = structure(list(coefficients = c(a = 1)), class = "rms"),
          new = 4,
          expected = c(a = 4)
        ),
        list(
          label = "robcov_vglm",
          model = structure(
            list(coefficients = c(a = 1, b = 2), vglm_fit = NULL),
            class = "robcov_vglm"
          ),
          new = c(5, 6),
          expected = c(a = 5, b = 6)
        ),
        list(
          label = "vglm",
          model = m_vglm,
          new = stats::coef(m_vglm) * 1.1,
          expected = stats::coef(m_vglm) * 1.1
        ),
        list(
          label = "orm",
          model = m_orm,
          new = stats::coef(m_orm) * 0.9,
          expected = stats::coef(m_orm) * 0.9
        )
      )

      for (case in cases) {
        original <- coefs_of(case$model)
        modified <- set_coef(case$model, case$new)
        expect_equal(coefs_of(modified), case$expected, info = case$label)
        # The model passed in keeps its coefficients.
        expect_equal(coefs_of(case$model), original, info = case$label)
      }

      # Changed vglm coefficients change the predictions.
      pred_orig <- predict(m_vglm, newdata = data[1:5, ], type = "response")
      pred_mod <- predict(
        set_coef(m_vglm, stats::coef(m_vglm) * 1.1),
        newdata = data[1:5, ],
        type = "response"
      )
      expect_false(all(pred_orig == pred_mod))

      expect_error(
        set_coef(list(), 1),
        "Model does not have a 'coefficients' element",
        fixed = TRUE
      )
      expect_error(
        set_coef(list(coefficients = c(a = 1, b = 2)), 1),
        "Length of new_coefs",
        fixed = TRUE
      )
    })
  })

  describe("get_vcov_robust()", {
    it("handles precomputed, fallback, formula, and unsupported cases", {
      robust <- structure(list(var = matrix(2, 1, 1)), class = "robcov_vglm")
      expect_equal(get_vcov_robust(robust), matrix(2, 1, 1))

      precomputed <- structure(
        list(var = matrix(3, 1, 1), orig.var = matrix(1, 1, 1)),
        class = "orm",
        markov_robust_covariance = list(
          backend = "orm",
          implementation = "mostr::robcov_orm",
          covariance_identity = matrix(3, 1, 1)
        )
      )
      expect_equal(get_vcov_robust(precomputed), matrix(3, 1, 1))

      fit <- lm(mpg ~ wt, data = mtcars)
      expect_equal(get_vcov_robust(fit), vcov(fit))
      expect_error(
        get_vcov_robust(fit, cluster = ~cyl),
        "Unsupported model class for robust vcov",
        fixed = TRUE
      )
      expect_error(
        get_vcov_robust(
          structure(list(call = list()), class = "orm"),
          cluster = ~ a + b
        ),
        "cluster formula must specify exactly one variable",
        fixed = TRUE
      )
      expect_error(
        get_vcov_robust(
          structure(list(call = list()), class = "orm"),
          cluster = ~id
        ),
        "Cluster variable 'id' not found",
        fixed = TRUE
      )
    })

    it("returns robust (sandwich) vcov when cluster is NULL for vglm models", {
      skip_on_cran()
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 100, seed = 333, follow_up_time = 15)
      # The vglm fit inside the patient-clustered robcov_vglm object
      m_vglm <- make_test_model(data)$vglm_fit

      # Get vcov without clustering
      V <- get_vcov_robust(m_vglm, cluster = NULL)

      # Should NOT match standard vcov exactly (sandwich vs model-based)
      # They will differ if the model is not perfectly specified
      expect_false(isTRUE(all.equal(V, vcov(m_vglm), tolerance = 1e-8)))

      # Should match robcov_vglm output
      V_rob <- robcov_vglm(m_vglm, cluster = NULL)$var
      expect_equal(V, V_rob)

      V_working <- get_vcov_robust(
        m_vglm,
        cluster = NULL,
        bread = "vglm",
        type = "HC1",
        cadjust = TRUE
      )
      expected_working <- robcov_vglm(
        m_vglm,
        cluster = NULL,
        bread = "vglm",
        type = "HC1",
        cadjust = TRUE
      )$var
      expect_equal(V_working, expected_working)
    })

    it("returns robust (sandwich) vcov when cluster is NULL for orm models", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 100, seed = 333, follow_up_time = 15)
      m_orm <- rms::orm(
        y ~ time_lin + time_nlin_1 + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      )

      # Get vcov without clustering
      V <- get_vcov_robust(m_orm, cluster = NULL)

      m_robust <- rms::robcov(m_orm)

      # Should NOT match standard vcov exactly (sandwich vs model-based)
      # They will differ if the model is not perfectly specified
      expect_false(isTRUE(all.equal(V, m_robust$orig.var, tolerance = 1e-8)))

      # Should match robcov_vglm output
      V_rob <- m_robust$var
      expect_equal(V, V_rob)
    })

    it("computes cluster-robust vcov from a cluster formula on the rows vglm used", {
      skip_on_cran()
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- mvn_cluster_case()$data
      data_with_extra <- rbind(data, data[seq_len(5), , drop = FALSE])
      data_with_extra$y[(nrow(data) + 1):nrow(data_with_extra)] <- NA

      cases <- list(
        list(
          label = "complete data",
          model = mvn_cluster_case()$model$vglm_fit,
          data = data
        ),
        # The rows with a missing outcome are dropped by the fit, so the
        # cluster ids must be aligned to the complete cases.
        list(
          label = "rows with a missing outcome",
          model = make_test_model(data_with_extra)$vglm_fit,
          data = data_with_extra
        )
      )

      for (case in cases) {
        V_formula <- get_vcov_robust(case$model, cluster = ~id, data = case$data)
        expected <- robcov_vglm(case$model, cluster = data$id)$var
        expect_equal(V_formula, expected, info = case$label)
        expect_false(
          isTRUE(all.equal(
            V_formula,
            robcov_vglm(case$model, cluster = NULL)$var
          )),
          info = case$label
        )
      }
    })
  })

  describe("inferences() with simulation method", {
    it("uses a supplied covariance matrix for the MVN draws", {
      skip_on_cran()
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(n_patients = 40, seed = 999, follow_up_time = 10)
      baseline_data <- data[data$time == 1, , drop = FALSE]
      m_vglm <- make_test_model(data)

      result <- avg_sops(
        model = m_vglm,
        newdata = baseline_data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        absorb = 6
      )
      run_mvn <- function(...) {
        inferences(result, method = "mvn", n_draws = 50, seed = 4440, ...)
      }

      default <- run_mvn()
      supplied <- run_mvn(vcov = vcov(m_vglm))
      inflated <- run_mvn(vcov = 2 * vcov(m_vglm))

      # Passing the model's own covariance reproduces the default exactly.
      expect_equal(supplied$std.error, default$std.error)
      expect_equal(supplied$conf.low, default$conf.low)
      expect_equal(supplied$conf.high, default$conf.high)
      # Doubling the covariance widens the draws in every row.
      expect_true(all(inflated$std.error > default$std.error))
    })
  })
})

test_that("seeded inference is invariant to outer worker count", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("future.callr")
  if (
    requireNamespace("pkgload", quietly = TRUE) &&
      pkgload::is_dev_package("mostr")
  ) {
    skip("Callr worker invariance requires an installed package.")
  }

  data <- make_test_data(n_patients = 40, seed = 889, follow_up_time = 6)
  model <- make_test_model(data)
  baseline <- data[data$time == 1, , drop = FALSE]
  point <- avg_sops(
    model,
    newdata = baseline,
    variables = list(tx = 0:1),
    times = 1:4,
    absorb = 6
  )

  sequential <- inferences(
    point,
    method = "mvn",
    n_draws = 7,
    workers = 1,
    seed = 123,
    return_draws = TRUE
  )
  parallel <- inferences(
    point,
    method = "mvn",
    n_draws = 7,
    workers = 2,
    seed = 123,
    return_draws = TRUE
  )

  expect_equal(
    as.data.frame(sequential),
    as.data.frame(parallel),
    tolerance = 1e-11
  )
  expect_equal(
    attr(sequential, "draws"),
    attr(parallel, "draws"),
    tolerance = 1e-11
  )
})
