# Bootstrap resamples can lack a state. drop_absent_states() then reports it,
# removes it from the factor levels of the response and the previous state,
# and returns the labels of the states and absorbing states present. Labels
# are never renumbered, and numeric columns keep their values. `warning` is
# the expected warning, or NA for none.
test_that("drop_absent_states drops missing states and keeps the state labels", {
  letters_abc <- c("a", "b", "c")
  cases <- list(
    list(
      label = "integer columns, an absorbing state is dropped",
      args = list(
        data = data.frame(y = c(1L, 2L), yprev = c(1L, 2L)),
        original_data = data.frame(y = factor(1:3), yprev = factor(1:3)),
        y_levels = 1:3,
        absorb = c(2, 3)
      ),
      warning = "were dropped",
      check = function(out, label) {
        expect_identical(out$absorb, "2", info = label)
        expect_identical(out$ylevels, c("1", "2"), info = label)
      }
    ),
    list(
      label = "the only absorbing state is dropped",
      args = list(
        data = data.frame(y = c(1, 2), yprev = c(1, 2)),
        original_data = data.frame(y = 1:3, yprev = 1:3),
        y_levels = 1:3,
        absorb = 3
      ),
      warning = "Set to NULL",
      check = function(out, label) {
        expect_null(out$absorb, info = label)
      }
    ),
    list(
      label = "character columns without y_levels",
      args = list(
        data = data.frame(y = c("a", "b"), yprev = c("a", "b")),
        original_data = data.frame(y = letters_abc, yprev = letters_abc),
        y_levels = NULL,
        absorb = NULL
      ),
      warning = NULL,
      check = function(out, label) {
        expect_equal(out$ylevels, NULL, info = label)
        expect_identical(typeof(out$data$y), "character", info = label)
      }
    ),
    list(
      label = "states inferred from the data alone",
      args = list(
        data = data.frame(y = c(1, 2), yprev = c(1, 2)),
        original_data = NULL,
        y_levels = NULL,
        absorb = NULL
      ),
      warning = NULL,
      check = function(out, label) {
        expect_equal(out$missing_levels, character(0), info = label)
      }
    ),
    list(
      label = "states from y_levels alone",
      args = list(
        data = data.frame(y = c(1L, 3L), yprev = c(1L, 3L)),
        factor_cols = c("y", "yprev"),
        original_data = NULL,
        y_levels = 1:3,
        absorb = 3
      ),
      warning = NULL,
      check = function(out, label) {
        expect_equal(out$missing_levels, "2", info = label)
        expect_identical(out$data$y, c(1L, 3L), info = label)
        expect_identical(out$absorb, "3", info = label)
        expect_identical(out$ylevels, c("1", "3"), info = label)
      }
    ),
    list(
      label = "ordered response with a numeric previous state",
      args = list(
        data = data.frame(
          y = ordered(c(1, 3, 4), levels = 1:4),
          yprev = c(1, 3, 4)
        ),
        factor_cols = c("y", "yprev"),
        original_data = data.frame(
          y = ordered(c(1, 2, 3, 4)),
          yprev = c(1, 2, 3, 4)
        ),
        y_levels = 1:4,
        absorb = 4
      ),
      warning = NULL,
      check = function(out, label) {
        expect_true(is.ordered(out$data$y), info = label)
        expect_identical(levels(out$data$y), c("1", "3", "4"), info = label)
        # A numeric previous state keeps its numbers, so rcs(yprev) or a
        # linear yprev effect sees the original scale.
        expect_identical(out$data$yprev, c(1, 3, 4), info = label)
        expect_identical(out$ylevels, c("1", "3", "4"), info = label)
        expect_identical(out$absorb, "4", info = label)
        expect_equal(out$missing_levels, "2", info = label)
      }
    ),
    list(
      label = "a factor previous state gets no empty level",
      args = list(
        data = data.frame(
          y = factor(c(1, 3, 4), levels = 1:4),
          yprev = factor(c(1, 3, 3), levels = 1:3)
        ),
        factor_cols = c("y", "yprev"),
        original_data = data.frame(
          y = factor(1:4),
          yprev = factor(c(1, 2, 3, 3), levels = 1:3)
        ),
        y_levels = 1:4,
        absorb = 4
      ),
      warning = NA,
      check = function(out, label) {
        # rms::orm() cannot fit a factor with an empty level.
        expect_identical(levels(out$data$yprev), c("1", "3"), info = label)
        expect_identical(as.character(out$data$yprev), c("1", "3", "3"), info = label)
        expect_identical(levels(out$data$y), c("1", "3", "4"), info = label)
      }
    ),
    list(
      label = "letter-labelled factor columns",
      args = list(
        data = data.frame(
          y = factor(c("a", "c", "a"), levels = letters_abc),
          yprev = factor(c("a", "a", "a"), levels = letters_abc)
        ),
        factor_cols = c("y", "yprev"),
        original_data = data.frame(
          y = factor(c("a", "b", "c", "a"), levels = letters_abc),
          yprev = factor(c("a", "a", "b", "b"), levels = letters_abc)
        ),
        y_levels = letters_abc,
        absorb = "c"
      ),
      warning = NA,
      check = function(out, label) {
        expect_identical(as.character(out$missing_levels), "b", info = label)
        expect_identical(levels(out$data$y), c("a", "c"), info = label)
        expect_identical(levels(out$data$yprev), "a", info = label)
        expect_identical(out$ylevels, c("a", "c"), info = label)
        # The absorbing state keeps its label and its place as the last state.
        expect_identical(out$absorb, "c", info = label)
      }
    )
  )

  for (case in cases) {
    relevel <- function() do.call(drop_absent_states, case$args)
    if (is.null(case$warning)) {
      out <- relevel()
    } else if (is.na(case$warning)) {
      expect_no_warning(out <- relevel())
    } else {
      expect_warning(out <- relevel(), case$warning, fixed = TRUE)
    }
    case$check(out, case$label)
  }
})

test_that("validate_coef_vcov reports structural mismatches and accepts Matrix covariance objects", {
  skip_if_not_installed("Matrix")

  beta <- c(a = 1, b = 2)
  named <- function(x, rows = names(beta), cols = names(beta)) {
    dimnames(x) <- list(rows, cols)
    x
  }
  cases <- list(
    list(label = "a vector", vcov = c(1, 2), error = "must be a matrix"),
    list(label = "too many rows", vcov = diag(3), error = "Dimension mismatch"),
    list(
      label = "row names in another order",
      vcov = named(diag(2), rows = c("b", "a")),
      error = "Coefficient names do not match row names"
    ),
    list(
      label = "column names in another order",
      vcov = named(diag(2), cols = c("b", "a")),
      error = "Coefficient names do not match column names"
    ),
    list(
      label = "a dense Matrix object, returned as a base matrix",
      vcov = named(Matrix::Matrix(diag(2), sparse = FALSE)),
      expected = named(diag(2))
    )
  )

  for (case in cases) {
    if (is.null(case$error)) {
      expect_equal(
        mostr:::validate_coef_vcov(beta, case$vcov),
        case$expected,
        info = case$label
      )
    } else {
      expect_error(
        mostr:::validate_coef_vcov(beta, case$vcov),
        case$error,
        fixed = TRUE,
        info = case$label
      )
    }
  }
})

test_that("set_coef() updates nested robcov_vglm fits", {
  skip_if_not_installed("VGAM")

  data <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = c(0, 0, 1, 1, 2, 2))
  fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)

  expect_error(
    set_coef(fit, 1),
    "Length of new_coefs",
    fixed = TRUE
  )

  new_coefs <- stats::coef(fit) + 0.1
  robust <- structure(
    list(coefficients = stats::coef(fit), vglm_fit = fit),
    class = "robcov_vglm"
  )
  robust2 <- set_coef(robust, new_coefs)

  expect_equal(robust2$coefficients, new_coefs)
  expect_equal(stats::coef(robust2$vglm_fit), new_coefs)
})

test_that("get_effective_coefs() reports unsupported and mismatched coefficients", {
  skip_if_not_installed("VGAM")

  expect_error(
    get_effective_coefs(list()),
    "only supports vglm and orm",
    fixed = TRUE
  )

  data <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = c(0, 0, 1, 1, 2, 2))
  fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)
  expect_warning(
    gamma <- mostr:::get_effective_coefs_vglm(
      fit,
      beta = c(stats::coef(fit), extra = 1)
    ),
    "Mismatch in coefficient consumption",
    fixed = TRUE
  )
  expect_equal(ncol(gamma), length(VGAM::constraints(fit)))

  orm <- structure(
    list(coefficients = c("y>=2" = 0, x = 1), non.slopes = NULL),
    class = "orm"
  )
  expect_error(
    mostr:::get_effective_coefs_orm(orm),
    "Cannot determine",
    fixed = TRUE
  )
  orm$non.slopes <- 2
  expect_error(
    mostr:::get_effective_coefs_orm(orm, beta = 1),
    "shorter than",
    fixed = TRUE
  )
  expect_error(
    mostr:::get_effective_coefs_orm(orm, beta = c(1, 2, 3)),
    "Length of beta",
    fixed = TRUE
  )
})

test_that("previous-state coercion respects scalar storage modes", {
  expect_equal(
    mostr:::coerce_previous_state_values(c("1", "2"), 1L, "yprev"),
    c(1L, 2L)
  )
  expect_error(
    mostr:::coerce_previous_state_values("home", 1L, "yprev"),
    "integer-compatible",
    fixed = TRUE
  )

  expect_equal(
    mostr:::coerce_previous_state_values(c("1.5", "2"), 1, "yprev"),
    c(1.5, 2)
  )
  expect_error(
    mostr:::coerce_previous_state_values("home", 1, "yprev"),
    "numeric-compatible",
    fixed = TRUE
  )

  expect_equal(
    mostr:::coerce_previous_state_values(c(1, 2), "1", "yprev"),
    c("1", "2")
  )
  expect_equal(
    mostr:::coerce_previous_state_values(c(1, 2), list(), "yprev"),
    c("1", "2")
  )
})

test_that("orm prediction helpers validate metadata and categorical levels", {
  skip_if_not_installed("rms")

  expect_error(
    mostr:::orm_model_matrix(list(), data.frame(x = 1)),
    "requires an orm model",
    fixed = TRUE
  )
  expect_error(
    mostr:::orm_model_matrix(
      structure(list(Design = NULL, sformula = NULL), class = "orm"),
      data.frame(x = 1)
    ),
    "stored Design metadata",
    fixed = TRUE
  )

  data <- data.frame(
    y = ordered(c(1, 2, 3, 1, 2, 3)),
    z = factor(c("a", "b", "a", "b", "a", "b"))
  )
  fit <- rms::orm(y ~ z, data = data, x = TRUE, y = TRUE)
  X <- mostr:::orm_model_matrix(
    fit,
    data.frame(y = ordered(c(1, 1), levels = 1:3), z = c("a", "b")),
    include_intercept = TRUE
  )

  expect_equal(nrow(X), 2)
  expect_true("(Intercept)" %in% colnames(X))
  expect_error(
    mostr:::orm_model_matrix(
      fit,
      data.frame(y = ordered(1, levels = 1:3), z = "c")
    ),
    "not among fitted orm levels",
    fixed = TRUE
  )

  bad_fit <- fit
  bad_fit$Design$colnames <- c(bad_fit$Design$colnames, "extra")
  expect_error(
    mostr:::orm_model_matrix(
      bad_fit,
      data.frame(y = ordered(1, levels = 1:3), z = "a")
    ),
    "Could not reconstruct",
    fixed = TRUE
  )

  expect_equal(
    mostr:::normalize_orm_prediction_data(
      fit,
      data.frame(y = ordered(1, levels = 1:3))
    ),
    data.frame(y = ordered(1, levels = 1:3))
  )

  strat_env <- new.env(parent = baseenv())
  strat_env$strat <- function(x) x
  strat_fit <- structure(
    list(
      Design = list(
        name = character(0),
        assume.code = integer(0),
        parms = list(),
        colnames = "x"
      ),
      sformula = stats::as.formula("y ~ x + strat(z)", env = strat_env)
    ),
    class = "orm"
  )
  X_strat <- mostr:::orm_model_matrix(
    strat_fit,
    data.frame(x = 1:2, z = factor(c("a", "b")))
  )
  expect_equal(colnames(X_strat), "x")
  expect_equal(unname(X_strat[, "x"]), c(1, 2))
})

test_that("orm score-bootstrap helpers cover validation and alignment branches", {
  skip_if_not_installed("rms")

  expect_error(
    mostr:::compute_scores_orm(list()),
    "must be an orm object",
    fixed = TRUE
  )
  expect_error(
    mostr:::compute_scores_orm(structure(
      list(x = NULL, y = NULL),
      class = "orm"
    )),
    "x = TRUE, y = TRUE",
    fixed = TRUE
  )

  data <- data.frame(
    y = ordered(c(1, 2, 3, 1, 2, 3, 1, 2, 3)),
    x = c(-1, 0, 1, -0.5, 0.5, 1.5, -1.5, 0.2, 1.2)
  )
  fit <- rms::orm(y ~ x, data = data, x = TRUE, y = TRUE)
  scores <- mostr:::compute_scores_orm(fit)
  expect_equal(nrow(scores), nrow(data))
  expect_equal(colnames(scores), names(stats::coef(fit)))

  bad_fit <- fit
  bad_fit$y <- as.character(bad_fit$y)
  bad_fit$y[1] <- 99
  expect_error(
    mostr:::compute_scores_orm(bad_fit),
    "Could not match orm response",
    fixed = TRUE
  )

  expect_error(
    mostr:::align_cluster_orm(fit, ~id, nrow(scores)),
    "row-level vector",
    fixed = TRUE
  )
  expect_equal(
    mostr:::align_cluster_orm(fit, seq_len(nrow(scores)), nrow(scores)),
    seq_len(nrow(scores))
  )

  na_fit <- fit
  na_fit$na.action <- c(2L, 4L)
  expect_equal(
    mostr:::align_cluster_orm(
      na_fit,
      letters[1:(nrow(scores) + 2)],
      nrow(scores)
    ),
    letters[-c(2, 4)][1:nrow(scores)]
  )
  expect_error(
    mostr:::align_cluster_orm(fit, 1:2, nrow(scores)),
    "does not match",
    fixed = TRUE
  )

  expect_true(is.matrix(mostr:::orm_model_bread(fit)$bread))
})

test_that("score bootstrap draw helper reports malformed robust components", {
  baseline <- data.frame(id = c("a", "b"))
  robust <- structure(
    list(
      coefficients = c(a = 0, b = 0),
      bread = diag(2),
      scores = matrix(0, nrow = 2, ncol = 2),
      cluster = c("a", "b")
    ),
    class = "robcov_vglm"
  )

  expect_error(
    mostr:::generate_score_bootstrap_draws(
      robust,
      baseline,
      "missing",
      1
    ),
    "ID variable",
    fixed = TRUE
  )

  bad_scores <- robust
  bad_scores$scores <- data.frame(a = 1)
  expect_error(
    mostr:::generate_score_bootstrap_draws(bad_scores, baseline, "id", 1),
    "must be a matrix",
    fixed = TRUE
  )

  bad_bread <- robust
  bad_bread$bread <- c(1, 2)
  expect_error(
    mostr:::generate_score_bootstrap_draws(bad_bread, baseline, "id", 1),
    "`model$bread` must be a matrix",
    fixed = TRUE
  )

  bad_cluster <- robust
  bad_cluster$cluster <- "a"
  expect_error(
    mostr:::generate_score_bootstrap_draws(
      bad_cluster,
      baseline,
      "id",
      1
    ),
    "does not match",
    fixed = TRUE
  )

  bad_coef <- robust
  bad_coef$coefficients <- c(a = 0)
  expect_error(
    mostr:::generate_score_bootstrap_draws(bad_coef, baseline, "id", 1),
    "Coefficient length",
    fixed = TRUE
  )

  expect_error(
    mostr:::generate_score_bootstrap_draws(
      robust,
      baseline[0, , drop = FALSE],
      "id",
      1
    ),
    "Invalid baseline weights",
    fixed = TRUE
  )
})

test_that("robcov_vglm and score extraction cover malformed internals", {
  skip_if_not_installed("VGAM")

  data <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = c(0, 0, 1, 1, 2, 2))
  fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)
  n <- stats::nobs(fit, type = "lm")
  p <- length(stats::coef(fit))

  expect_error(
    compute_scores_vglm(lm(y ~ x, data = data)),
    "must be a vglm object",
    fixed = TRUE
  )

  scores <- compute_scores_vglm(fit)
  expect_equal(dim(scores), c(n, p))

  with_mocked_bindings(
    compute_scores_vglm = function(fit) matrix(0, nrow = n + 1, ncol = p),
    {
      expect_error(
        robcov_vglm(fit),
        "score rows",
        fixed = TRUE
      )
    }
  )
  with_mocked_bindings(
    compute_scores_vglm = function(fit) matrix(0, nrow = n, ncol = p + 1),
    {
      expect_error(
        robcov_vglm(fit),
        "score columns",
        fixed = TRUE
      )
    }
  )
})

test_that("print.robcov_vglm reports clustered fits", {
  obj <- structure(
    list(
      call = NULL,
      n = 4,
      n_clusters = 2,
      coefficients = c(a = 1),
      se = c(a = 0.5)
    ),
    class = "robcov_vglm"
  )

  expect_output(print(obj), "Number of clusters")
})

test_that("get_vcov_robust aligns a cluster formula with the rows an orm fit used", {
  skip_if_not_installed("rms")

  # The states overlap across `x`, so rms can estimate the coefficients. The
  # missing `x` in row 2 makes orm drop that row, so the cluster IDs from `data`
  # must drop it too; keeping the first 14 IDs instead would shift patients.
  data <- data.frame(
    y = ordered(c(1, 2, 1, 3, 2, 1, 3, 2, 3, 1, 3, 2, 1, 2, 3)),
    x = c(1, NA, 3:15),
    id = rep(1:5, each = 3)
  )
  orm_fit <- rms::orm(y ~ x, data = data, x = TRUE, y = TRUE)
  aligned <- mostr:::robcov_orm(orm_fit, cluster = data$id[-2])$var
  shifted <- mostr:::robcov_orm(orm_fit, cluster = data$id[1:14])$var

  from_formula <- mostr:::get_vcov_robust(orm_fit, cluster = ~id, data = data)
  expect_equal(from_formula, aligned)
  expect_false(isTRUE(all.equal(from_formula, shifted)))
  expect_equal(
    mostr:::get_vcov_robust(orm_fit, cluster = data$id[-2]),
    aligned
  )
})

test_that("markov_msm_build validates model classes and design alignment", {
  skip_if_not_installed("VGAM")

  coefs_named <- function(columns) {
    function(model, beta = NULL) {
      matrix(0, nrow = 1, ncol = length(columns), dimnames = list(NULL, columns))
    }
  }
  # Full columns at the first visit, then a design that lacks `x`.
  shrinking_design <- function() {
    calls <- 0
    function(model, newdata, include_intercept = TRUE) {
      calls <<- calls + 1
      columns <- if (calls == 1) c("(Intercept)", "x") else "(Intercept)"
      matrix(
        1,
        nrow = nrow(newdata),
        ncol = length(columns),
        dimnames = list(NULL, columns)
      )
    }
  }
  vglm_fit <- VGAM::vglm(
    y ~ x,
    family = VGAM::binomialff,
    data = data.frame(
      y = c(0, 1, 1, 0, 0, 1),
      x = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    )
  )

  cases <- list(
    list(
      label = "an unsupported model class",
      mocks = list(get_effective_coefs = coefs_named("(Intercept)")),
      args = list(
        structure(list(), class = "not_supported"),
        data.frame(yprev = 1),
        settings = markov_settings(y_levels = 1:6)
      ),
      error = "supports only vglm and orm"
    ),
    list(
      label = "time_covariates rows that do not match the times",
      mocks = list(),
      args = list(
        structure(list(), class = "orm"),
        data.frame(yprev = 1),
        times = 1,
        settings = markov_settings(
          y_levels = 1:6,
          time_covariates = data.frame(z = 1:2)
        )
      ),
      error = "one row per requested time point"
    ),
    list(
      label = "an orm design that loses a column at a later visit",
      mocks = list(
        get_effective_coefs = coefs_named(c("(Intercept)", "x")),
        orm_model_matrix = shrinking_design()
      ),
      args = list(
        structure(list(), class = "orm"),
        data.frame(id = 1, time = 2, yprev = 1),
        times = c(2, 3),
        settings = markov_settings(y_levels = 1:2)
      ),
      error = "Design matrix columns missing"
    ),
    list(
      label = "a vglm design without the coefficient columns",
      mocks = list(get_effective_coefs = coefs_named("missing_column")),
      args = list(
        model = vglm_fit,
        newdata = data.frame(id = 1:2, x = c(0, 1)),
        settings = markov_settings(y_levels = 1:2)
      ),
      error = "Could not construct a design matrix"
    )
  )

  build <- function(case) {
    if (length(case$mocks) > 0L) {
      do.call(local_mocked_bindings, case$mocks)
    }
    do.call(mostr:::markov_msm_build, case$args)
  }
  for (case in cases) {
    expect_error(build(case), case$error, fixed = TRUE, info = case$label)
  }
})
