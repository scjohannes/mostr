# Regression tests for issue #53: an explicit `cluster` for analytical
# inference must group the fitting rows by patient, `get_vcov_robust()` must
# honour its arguments on vglm fits from `mostm()`, and a cluster vector given
# before `NA` rows were dropped must be aligned when the fit used
# `Hmisc::na.delete`.

cluster_alignment_data <- function() {
  data <- make_test_data(n_patients = 40, times = 1:6, seed = 5301)
  # Five sites, each holding whole patients.
  data$site <- as.integer(factor(data$id)) %% 5L + 1L
  data
}

cluster_alignment_fits <- local({
  fits <- list()
  function(engine) {
    if (engine == "vglm") {
      skip_if_not_installed("VGAM")
    }
    if (is.null(fits[[engine]])) {
      data <- cluster_alignment_data()
      fits[[engine]] <<- mostm(
        ordered(y) ~ time + tx + yprev,
        data = data,
        id_var = "id",
        absorb = 6,
        engine = engine
      )
    }
    fits[[engine]]
  }
})

# Data with missing values in a predictor on follow-up rows only, so every
# patient keeps a complete starting row.
cluster_alignment_na_data <- function() {
  data <- make_test_data(n_patients = 40, times = 1:6, seed = 5302)
  data$z <- withr::with_seed(5303, stats::rnorm(nrow(data)))
  follow_up <- which(data$time > 1)
  data$z[follow_up[seq(3L, length(follow_up), by = 9L)]] <- NA
  data
}

test_that("delta inference refuses a cluster that does not group rows by patient", {
  skip_on_cran()
  for (engine in c("orm", "vglm")) {
    fit <- cluster_alignment_fits(engine)
    fitted <- markov_model_data(fit)
    base <- avg_sops(fit, times = 1:4)

    # The data reordered after the fit: the IDs no longer line up with the
    # fitting rows, so the scores of different patients would be summed.
    scrambled <- withr::with_seed(5304, sample(fitted$id))
    for (vcov in c("conditional", "unconditional")) {
      expect_error(
        inferences(base, method = "delta", vcov = vcov, cluster = scrambled),
        "group the fitting rows",
        class = "mostr_error",
        info = paste(engine, vcov)
      )
    }
    # One cluster per fitting row ignores the correlation within a patient.
    expect_error(
      inferences(
        base,
        method = "delta",
        vcov = "conditional",
        cluster = seq_len(nrow(fitted))
      ),
      "group the fitting rows",
      class = "mostr_error",
      info = engine
    )
    # Sites hold several patients; clusters other than patients are not
    # supported.
    expect_error(
      inferences(base, method = "delta", vcov = "conditional", cluster = ~site),
      "group the fitting rows",
      class = "mostr_error",
      info = engine
    )

    # The patient IDs in fitted order, relabelled, still give the default.
    default <- inferences(base, method = "delta", vcov = "conditional")
    relabelled <- inferences(
      base,
      method = "delta",
      vcov = "conditional",
      cluster = paste0("patient-", fitted$id)
    )
    expect_equal(relabelled$std.error, default$std.error, info = engine)
  }
})

test_that("get_vcov_robust() honours its arguments on vglm fits from mostm()", {
  skip_on_cran()
  fit <- cluster_alignment_fits("vglm")
  raw <- sop_delta_backend_model(fit)
  ids <- markov_model_data(fit)$id

  expect_identical(get_vcov_robust(fit), fit$var)
  # Without `cluster`, a correction override treats the fitting rows as
  # independent, as for orm fits.
  expect_equal(
    get_vcov_robust(fit, type = "HC1"),
    robcov_vglm(raw, cluster = NULL, type = "HC1")$var
  )
  expect_equal(
    get_vcov_robust(fit, cluster = ~id, type = "HC1"),
    robcov_vglm(raw, cluster = ids, type = "HC1")$var
  )
  expect_equal(
    get_vcov_robust(fit, cluster = ~id, cadjust = FALSE),
    robcov_vglm(raw, cluster = ids, cadjust = FALSE)$var
  )
  expect_equal(
    get_vcov_robust(fit, cluster = ~id, bread = "vglm"),
    robcov_vglm(raw, cluster = ids, bread = "vglm")$var
  )
  rows <- seq_along(ids)
  expect_equal(
    get_vcov_robust(fit, cluster = rows),
    robcov_vglm(raw, cluster = rows)$var
  )
  expect_false(isTRUE(all.equal(
    get_vcov_robust(fit, cluster = ~id, type = "HC1"),
    fit$var
  )))
})

test_that("a pre-NA cluster vector aligns with orm fits that used na.delete", {
  skip_on_cran()
  skip_if_not_installed("rms")
  data <- cluster_alignment_na_data()
  fit <- mostm(
    ordered(y) ~ time + tx + yprev + z,
    data = data,
    id_var = "id",
    absorb = 6
  )
  expect_s3_class(fit$na.action, "delete")
  expect_lt(nrow(fit$x), nrow(data))

  expect_equal(robcov_orm(fit, cluster = data$id)$var, fit$var)
  expect_equal(
    get_vcov_robust(fit, cluster = data$id),
    get_vcov_robust(fit, cluster = ~id)
  )
  base <- avg_sops(fit, times = 1:4)
  expect_equal(
    inferences(base, method = "delta", cluster = data$id)$std.error,
    inferences(base, method = "delta")$std.error
  )
})

test_that("a pre-NA cluster vector aligns with vglm fits that used na.delete", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  data <- cluster_alignment_na_data()
  fit <- mostm(
    ordered(y) ~ time + tx + yprev + z,
    data = data,
    id_var = "id",
    absorb = 6,
    engine = "vglm",
    na.action = Hmisc::na.delete
  )
  raw <- sop_delta_backend_model(fit)
  expect_lt(stats::nobs(raw, type = "lm"), nrow(data))

  expect_equal(robcov_vglm(raw, cluster = data$id)$var, fit$var)
  base <- avg_sops(fit, times = 1:4)
  expect_equal(
    inferences(base, method = "delta", cluster = data$id)$std.error,
    inferences(base, method = "delta")$std.error
  )
})

test_that("cluster alignment reads the omitted rows of an na.delete record", {
  # `Hmisc::na.delete()` records the number of missing values per variable
  # next to the omitted rows. Read as row numbers, the counts 0 and 2 would
  # drop row 2 as well.
  delete_record <- structure(
    list(
      nmiss = c(y = 0L, x = 2L),
      omit = c("4" = 4L, "7" = 7L),
      na.detail.response = NULL
    ),
    class = "delete"
  )
  ids <- letters[1:10]
  expected <- ids[-c(4L, 7L)]

  orm_fit <- list(na.action = delete_record)
  expect_equal(mostr:::align_cluster_orm(orm_fit, ids, 8L), expected)

  skip_if_not_installed("VGAM")
  vglm_fit <- methods::new("vglm")
  methods::slot(vglm_fit, "na.action") <- list(delete_record)
  expect_equal(mostr:::align_cluster_vglm(vglm_fit, ids, 8L), expected)
})
