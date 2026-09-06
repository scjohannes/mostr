test_that("bootstrap_model_coefs validates inputs", {
  expect_error(
    bootstrap_model_coefs(
      lm(mpg ~ wt, data = mtcars),
      data = mtcars,
      n_boot = 1
    ),
    "model must be an orm object",
    fixed = TRUE
  )

  model <- structure(list(), class = "vglm")
  expect_error(
    bootstrap_model_coefs(model, data = NULL, n_boot = 1),
    "No data provided",
    fixed = TRUE
  )

  expect_error(
    bootstrap_model_coefs(model, data = data.frame(patient = 1:2), n_boot = 1),
    "id_var 'id' not found in data",
    fixed = TRUE
  )

  expect_warning(
    expect_error(
      bootstrap_model_coefs(
        lm(mpg ~ wt, data = mtcars),
        data = mtcars,
        n_boot = 1,
        parallel = TRUE
      ),
      "model must be an orm object",
      fixed = TRUE
    ),
    "parallel.*deprecated"
  )
})

test_that("bootstrap_model_coefs orchestrates bootstrap coefficient extraction", {
  model <- structure(
    list(coefficients = c("(Intercept)" = 1, tx = -0.2)),
    class = "vglm"
  )
  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = c(0, 1, 0, 1),
    tx = c(0, 0, 1, 1),
    yprev = factor(c(0, 0, 0, 0))
  )

  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      expect_equal(id_var, "id")
      replicate(
        n_boot,
        data.frame(
          original_id = c(1, 2),
          new_id = c("1_1", "2_1"),
          boot_id = 1
        ),
        simplify = FALSE
      )
    },
    apply_to_bootstrap = function(
      boot_samples,
      analysis_fn,
      data,
      id_var,
      workers,
      packages,
      globals
    ) {
      expect_length(boot_samples, 2)
      expect_equal(workers, 1)
      list(
        list("(Intercept)" = 1, tx = -0.1),
        list("(Intercept)" = 2, tx = -0.3)
      )
    },
    {
      expect_warning(
        result <- bootstrap_model_coefs(
          model,
          data = data,
          n_boot = 2,
          workers = 1,
          parallel = FALSE
        ),
        "parallel.*deprecated"
      )
    }
  )

  expect_equal(result$boot_id, 1:2)
  expect_equal(result$tx, c(-0.1, -0.3))
})

test_that("bootstrap_model_coefs preserves numeric yprev before bootstrapping", {
  model <- structure(
    list(coefficients = c("(Intercept)" = 1, tx = -0.2)),
    class = "vglm"
  )
  data <- data.frame(
    id = c(1, 2),
    y = c(0, 1),
    tx = c(0, 1),
    yprev = c(0, 1)
  )

  observed_is_factor <- NULL
  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      observed_is_factor <<- is.factor(data$yprev)
      list(data.frame(
        original_id = c(1, 2),
        new_id = c("1_1", "2_1"),
        boot_id = 1
      ))
    },
    apply_to_bootstrap = function(...) list(list(tx = 1)),
    {
      result <- bootstrap_model_coefs(model, data = data, n_boot = 1)
    }
  )

  expect_false(observed_is_factor)
  expect_equal(result$tx, 1)
})

test_that("bootstrap_model_coefs records coefficients and failed refits", {
  model <- structure(
    list(coefficients = c("(Intercept)" = 1, tx = 2)),
    class = "orm"
  )
  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = factor(c(1, 2, 1, 2)),
    yprev = factor(c(1, 1, 1, 1)),
    tx = c(0, 0, 1, 1)
  )

  call_id <- 0
  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      replicate(
        n_boot,
        data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1")),
        simplify = FALSE
      )
    },
    apply_to_bootstrap = function(
      boot_samples,
      analysis_fn,
      data,
      id_var,
      workers,
      packages,
      globals
    ) {
      lapply(boot_samples, function(sample) analysis_fn(data))
    },
    bootstrap_analysis_wrapper = function(...) {
      call_id <<- call_id + 1
      if (call_id == 2) {
        list(model = NULL)
      } else {
        list(model = model)
      }
    },
    {
      result <- bootstrap_model_coefs(
        model,
        data = data,
        n_boot = 2,
        workers = 1,
        use_coefstart = TRUE
      )
    }
  )

  expect_equal(result$tx, c(2, NA))
})
