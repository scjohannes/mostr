# The vglm engine gives each column of an inline rms basis its own VGAM
# assignment so that `constraints` can name single spline columns. A `parallel`
# formula or a `constraints` list that names the whole basis term must still
# fit the model VGAM::vglm() fits (issue #42): before the fix,
# `parallel = FALSE ~ rms::rcs(time, 4)` matched no column and fitted a
# proportional odds spline without any warning.

# Constraints for `ordered(y) ~ yprev + rms::rcs(time, 4) + tx` that name the
# spline's single columns: a linear trend across thresholds for its linear
# column, and proportional odds for everything else.
vglm_column_constraints <- function(n_thresholds) {
  po <- matrix(1, n_thresholds, 1)
  list(
    "(Intercept)" = diag(n_thresholds),
    "yprev" = po,
    "rms::rcs(time, 4)time" = cbind(1, seq_len(n_thresholds) - 1),
    "rms::rcs(time, 4)time'" = po,
    "rms::rcs(time, 4)time''" = po,
    "tx" = po
  )
}

test_that("vglm fits with parallel formulas or term constraints on rms bases match VGAM::vglm()", {
  skip_on_cran() # Slow: each mostm() vglm fit takes about 2 s here.
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 60, seed = 2471, times = 1:12)
  n_thresholds <- nlevels(data$y) - 1L

  po_constraints <- function(formula) {
    po_fit <- VGAM::vglm(
      formula,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
    VGAM::constraints(po_fit, type = "term")
  }
  rcs_formula <- ordered(y) ~ yprev + rms::rcs(time, 4) + tx
  lsp_formula <- ordered(y) ~ yprev + rms::lsp(time, 6) + tx
  term_constraints <- po_constraints(rcs_formula)
  term_constraints[["rms::rcs(time, 4)"]] <- diag(n_thresholds)

  cases <- list(
    list(
      label = "parallel = FALSE ~ rcs()",
      formula = rcs_formula,
      family = VGAM::cumulative(
        reverse = TRUE,
        parallel = FALSE ~ rms::rcs(time, 4)
      ),
      constraints = NULL,
      # The spline has one coefficient per column and threshold.
      basis_coefs = 3L * n_thresholds
    ),
    list(
      label = "parallel = TRUE ~ -1 + yprev + lsp()",
      formula = lsp_formula,
      family = VGAM::cumulative(
        reverse = TRUE,
        parallel = TRUE ~ -1 + yprev + rms::lsp(time, 6)
      ),
      constraints = NULL,
      # The spline is common to all thresholds; only tx varies.
      basis_coefs = 2L
    ),
    list(
      label = "constraints naming the rcs() term",
      formula = rcs_formula,
      # `parallel = FALSE`, because with `parallel = TRUE` VGAM's cm.VGAM()
      # replaces every supplied constraint matrix with the parallel one, so
      # even VGAM::vglm() would fit a proportional odds spline.
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
      constraints = term_constraints,
      basis_coefs = 3L * n_thresholds
    )
  )

  for (case in cases) {
    reference <- VGAM::vglm(
      case$formula,
      family = case$family,
      data = data,
      constraints = case$constraints
    )
    fit <- suppressWarnings(mostr::mostm(
      case$formula,
      family = case$family,
      data = data,
      id_var = "id",
      absorb = 6,
      constraints = case$constraints,
      engine = "vglm"
    ))$vglm_fit

    expect_length(
      grep("^rms::(rcs|lsp)\\(", names(VGAM::coef(reference))),
      case$basis_coefs
    )
    expect_equal(
      names(VGAM::coef(fit)),
      names(VGAM::coef(reference)),
      info = case$label
    )
    expect_equal(
      VGAM::coef(fit),
      VGAM::coef(reference),
      tolerance = 1e-8,
      info = case$label
    )
    expect_equal(
      as.numeric(stats::logLik(fit)),
      as.numeric(stats::logLik(reference)),
      tolerance = 1e-10,
      info = case$label
    )
  }
})

test_that("an unqualified lsp() term in formula and parallel formula is not parallel", {
  skip_on_cran() # Slow: each mostm() vglm fit takes about 2 s here.
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # The VIOLET specification writes `lsp(day, 2)` without `rms::`; mostm()
  # supplies the helper, so the formula environment need not have rms.
  data <- make_test_data(n_patients = 60, seed = 2471, times = 1:12)
  n_thresholds <- nlevels(data$y) - 1L
  fit <- suppressWarnings(mostr::mostm(
    stats::as.formula("ordered(y) ~ yprev + lsp(time, 6) + tx", env = baseenv()),
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ lsp(time, 6)),
    data = data,
    id_var = "id",
    absorb = 6,
    engine = "vglm"
  ))$vglm_fit
  reference <- VGAM::vglm(
    ordered(y) ~ yprev + rms::lsp(time, 6) + tx,
    family = VGAM::cumulative(
      reverse = TRUE,
      parallel = FALSE ~ rms::lsp(time, 6)
    ),
    data = data
  )

  expect_length(
    grep("^lsp\\(time, 6\\)", names(VGAM::coef(fit))),
    2L * n_thresholds
  )
  expect_equal(
    unname(VGAM::coef(fit)),
    unname(VGAM::coef(reference)),
    tolerance = 1e-8
  )
})

test_that("vglm refuses a parallel formula naming a basis term whose columns constraints name", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # A `constraints` list that names single spline columns makes the fit use one
  # VGAM assignment per column, where a `parallel` formula naming the whole
  # spline term matches nothing. VGAM would then silently replace the user's
  # spline matrices by the parallel one, so mostm() refuses the combination.
  data <- make_test_data(n_patients = 60, seed = 2471, times = 1:12)
  column_constraints <- vglm_column_constraints(nlevels(data$y) - 1L)
  fit_with <- function(parallel) {
    suppressWarnings(mostr::mostm(
      ordered(y) ~ yprev + rms::rcs(time, 4) + tx,
      family = VGAM::cumulative(reverse = TRUE, parallel = parallel),
      data = data,
      id_var = "id",
      absorb = 6,
      constraints = column_constraints,
      engine = "vglm"
    ))
  }

  for (parallel in list(
    FALSE ~ rms::rcs(time, 4),
    TRUE ~ -1 + rms::rcs(time, 4)
  )) {
    expect_error(
      fit_with(parallel),
      regexp = "rms::rcs(time, 4)",
      fixed = TRUE,
      class = "mostr_error"
    )
  }
})

test_that("vglm fits per-column basis constraints with a parallel formula naming other terms", {
  skip_on_cran() # Slow: each mostm() vglm fit takes about 2 s here.
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # A formula that names only other terms, or a logical `parallel`, is not
  # refused: the per-column constraints are fitted as given.
  data <- make_test_data(n_patients = 60, seed = 2471, times = 1:12)
  column_constraints <- vglm_column_constraints(nlevels(data$y) - 1L)
  fit_with <- function(parallel) {
    suppressWarnings(mostr::mostm(
      ordered(y) ~ yprev + rms::rcs(time, 4) + tx,
      family = VGAM::cumulative(reverse = TRUE, parallel = parallel),
      data = data,
      id_var = "id",
      absorb = 6,
      constraints = column_constraints,
      engine = "vglm"
    ))
  }
  expect_no_error(fit_with(FALSE ~ tx))
  fit <- fit_with(FALSE)$vglm_fit
  expect_length(
    grep("^rms::rcs\\(time, 4\\)time:", names(VGAM::coef(fit))),
    2L
  )
})
