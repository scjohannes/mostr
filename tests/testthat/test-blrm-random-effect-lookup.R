# Which patient gets which random intercept with `include_re = TRUE` (#39).
#
# rmsb::blrm() numbers clusters as `as.integer(as.factor(cluster))`, and Stan
# names the draws `gamma[1]`, ..., `gamma[Nc]`. So `gamma[k]` belongs to
# `levels(as.factor(cluster))[k]`, not to the patient whose ID is `k`. These
# tests use no random numbers: each draw column holds the number of its own
# cluster, so a patient who gets another patient's column is seen at once.

# A fake fit whose random-effect draws are named as Stan names them.
# `cluster`: the cluster value of each row of the fitted data.
# Draw `d` of the cluster with label `x` is `as.numeric(x) + d / 100`.
stan_named_blrm <- function(cluster, cluster_name = "id", n_draws = 3L) {
  draws <- cbind(
    "y>=2" = c(0.50, 1.00, 0.25),
    "y>=3" = c(-0.50, -1.00, -0.75),
    tx = c(0.25, -0.25, 0.10)
  )[seq_len(n_draws), , drop = FALSE]
  cluster_levels <- levels(as.factor(cluster))
  gamma <- outer(
    seq_len(n_draws) / 100,
    as.numeric(cluster_levels),
    `+`
  )
  colnames(gamma) <- paste0("gamma[", seq_along(cluster_levels), "]")
  model <- fake_blrm(draws, gamma_draws = gamma)
  model$clusterInfo <- list(name = cluster_name, cluster = cluster)
  model
}

own_draws <- function(labels, draw_indices) {
  outer(draw_indices / 100, as.numeric(labels), `+`)
}

test_that("character IDs that look like numbers get their own random intercepts", {
  ids <- as.character(1:10)
  # Factor levels sort as "1", "10", "2", ..., so `gamma[2]` is patient "10".
  model <- stan_named_blrm(rep(ids, each = 2))
  newdata <- data.frame(id = c(ids, "10", "2"), stringsAsFactors = FALSE)

  u <- mostr:::blrm_random_effect_matrix(model, newdata, "id", 1:3)

  expect_equal(unname(u), own_draws(newdata$id, 1:3), tolerance = 1e-12)

  # The draws cached once per call by the prediction functions give the same.
  cached <- mostr:::cache_blrm_gamma_draws(model, c(3L, 1L), include_re = TRUE)
  u_cached <- mostr:::blrm_random_effect_matrix(
    model,
    newdata,
    "id",
    c(3L, 1L),
    gamma_draws = cached
  )
  expect_equal(
    unname(u_cached),
    own_draws(newdata$id, c(3, 1)),
    tolerance = 1e-12
  )
})

test_that("numeric IDs that do not start at 1 get their own random intercepts", {
  model <- stan_named_blrm(rep(2:11, each = 2))
  # Every ID in `newdata` is also a Stan index (1 to 10), but `gamma[k]` is
  # the patient with ID `k + 1`.
  newdata <- data.frame(id = 2:10)

  u <- mostr:::blrm_random_effect_matrix(model, newdata, "id", 1:3)

  expect_equal(unname(u), own_draws(newdata$id, 1:3), tolerance = 1e-12)
})

test_that("an ID that was not in the fit is an error, not a Stan index", {
  model <- stan_named_blrm(rep(101:110, each = 2))

  expect_error(
    mostr:::blrm_random_effect_matrix(model, data.frame(id = 7), "id", 1:3),
    "not present",
    fixed = TRUE
  )
  expect_error(
    mostr:::blrm_random_effect_matrix(
      model,
      data.frame(id = c(101, 7)),
      "id",
      1:3
    ),
    "not present",
    fixed = TRUE
  )
})

test_that("draws named by cluster value are still matched by name", {
  ids <- as.character(1:10)
  model <- stan_named_blrm(rep(ids, each = 2))
  # Columns in an order that is neither numeric nor that of the factor levels
  by_name <- own_draws(rev(ids), 1:3)
  colnames(by_name) <- rev(ids)
  model$gamma_draws <- by_name
  newdata <- data.frame(id = ids, stringsAsFactors = FALSE)

  u <- mostr:::blrm_random_effect_matrix(model, newdata, "id", 1:3)

  expect_equal(unname(u), own_draws(ids, 1:3), tolerance = 1e-12)
})

test_that("sops(include_re = TRUE) uses each patient's own random intercept", {
  local_mocked_bindings(blrm_design_matrix = fake_blrm_design("tx"))
  ids <- as.character(1:10)
  stan_named <- stan_named_blrm(rep(ids, each = 2))
  # The same draws with every column named by its patient
  named_by_id <- stan_named
  colnames(named_by_id$gamma_draws) <- levels(as.factor(ids))
  newdata <- data.frame(
    id = ids,
    tx = rep(c(0, 1), 5),
    yprev = factor(rep(c(1, 2), each = 5), levels = 1:3),
    time = 1,
    stringsAsFactors = FALSE
  )
  per_patient <- function(model) {
    sops(
      model,
      newdata = newdata,
      times = 1:2,
      include_re = TRUE,
      n_draws = 3,
      seed = 1,
      posterior_summary = "mean"
    )
  }

  expect_equal(
    per_patient(stan_named)$estimate,
    per_patient(named_by_id)$estimate,
    tolerance = 1e-12
  )
})

test_that("cluster() on another variable looks up that variable, not the patient", {
  # As for `cluster(site)`: three sites, numbered 30, 10 and 20. The patient
  # IDs 1 to 3 are all Stan indices, so matching them would give no error.
  model <- stan_named_blrm(
    c("30", "10", "20", "10", "30", "20"),
    cluster_name = "site"
  )
  newdata <- data.frame(
    id = 1:3,
    site = c("20", "30", "20"),
    stringsAsFactors = FALSE
  )

  u <- mostr:::blrm_random_effect_matrix(model, newdata, "id", 1:3)

  expect_equal(unname(u), own_draws(newdata$site, 1:3), tolerance = 1e-12)

  expect_error(
    mostr:::blrm_random_effect_matrix(model, newdata["id"], "id", 1:3),
    "`site`",
    fixed = TRUE
  )
  unseen <- newdata
  unseen$site[2] <- "2"
  expect_error(
    mostr:::blrm_random_effect_matrix(model, unseen, "id", 1:3),
    "not present",
    fixed = TRUE
  )
})

test_that("sops() and avg_sops() use the intercept of each row's cluster(site) value", {
  local_mocked_bindings(blrm_design_matrix = fake_blrm_design("tx"))
  # As for `mostm(..., id_var = "id")` with `cluster(site)`. As text the sites
  # sort "-1", "0.5", "2". The patient IDs 1 to 3 are all Stan indices, and
  # patient 4 is not a site, so a look-up by patient would be wrong or fail.
  model <- stan_named_blrm(
    c("2", "-1", "0.5", "2", "-1", "0.5"),
    cluster_name = "site"
  )
  expect_identical(mostr:::markov_model_id_var(model), "id")
  newdata <- data.frame(
    id = 1:4,
    site = c("2", "-1", "0.5", "2"),
    tx = c(0, 1, 1, 0),
    yprev = factor(c(1, 2, 1, 2), levels = 1:3),
    time = 1,
    stringsAsFactors = FALSE
  )
  # Hand calculation for the first visit: the mean over the three draws of
  # the state probabilities at eta_k = alpha_k + tx * beta + u_site, where
  # draw `d` of site `x` is `as.numeric(x) + d / 100`.
  expected <- function(site, tx, state) {
    mean(vapply(
      1:3,
      function(d) {
        u <- as.numeric(site) + d / 100
        p <- stats::plogis(
          model$draws[d, c("y>=2", "y>=3")] + tx * model$draws[d, "tx"] + u
        )
        c(1 - p[1], p[1] - p[2], p[2])[state]
      },
      numeric(1)
    ))
  }

  result <- sops(
    model,
    newdata = newdata,
    times = 1,
    include_re = TRUE,
    n_draws = NULL,
    posterior_summary = "mean"
  )

  expect_equal(
    result$estimate,
    mapply(expected, result$site, result$tx, as.integer(result$state)),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )

  # `variables` stacks copies of the profiles; each copy keeps its site.
  averaged <- avg_sops(
    model,
    newdata = newdata,
    variables = list(tx = c(0, 1)),
    times = 1,
    include_re = TRUE,
    n_draws = NULL,
    posterior_summary = "mean"
  )

  expect_equal(
    averaged$estimate,
    mapply(
      function(tx, state) {
        mean(vapply(newdata$site, expected, numeric(1), tx = tx, state = state))
      },
      averaged$tx,
      as.integer(averaged$state)
    ),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )

  # Without the site column the error names it, also through sops().
  expect_error(
    sops(
      model,
      newdata = newdata[names(newdata) != "site"],
      times = 1,
      include_re = TRUE,
      n_draws = NULL
    ),
    "`site`",
    fixed = TRUE
  )
})

test_that("draws that do not map onto the fitted clusters are an error", {
  model <- stan_named_blrm(rep(5:7, each = 2))
  newdata <- data.frame(id = 5)

  # One column fewer than fitted clusters
  fewer <- model
  fewer$gamma_draws <- model$gamma_draws[, 1:2, drop = FALSE]
  expect_error(
    mostr:::blrm_random_effect_matrix(fewer, newdata, "id", 1:3),
    "cannot be matched to clusters",
    fixed = TRUE
  )
  expect_error(
    mostr:::cache_blrm_gamma_draws(fewer, 1:3, include_re = TRUE),
    "cannot be matched to clusters",
    fixed = TRUE
  )

  # The right number of columns, but not Stan's numbers 1 to 3
  shifted <- model
  colnames(shifted$gamma_draws) <- paste0("gamma[", 2:4, "]")
  expect_error(
    mostr:::blrm_random_effect_matrix(shifted, newdata, "id", 1:3),
    "cannot be matched to clusters",
    fixed = TRUE
  )

  # Columns without names: one per cluster is needed to put them in order
  unnamed <- model
  unnamed$gamma_draws <- unname(model$gamma_draws)[, 1:2, drop = FALSE]
  expect_error(
    mostr:::blrm_random_effect_matrix(unnamed, newdata, "id", 1:3),
    "cannot be matched to clusters",
    fixed = TRUE
  )
})

test_that("Stan's draw columns are renamed by their number, not their position", {
  ids <- as.character(1:10)
  model <- stan_named_blrm(rep(ids, each = 2))
  model$gamma_draws <- model$gamma_draws[, 10:1]

  named <- mostr:::get_blrm_gamma_draws(model, 1:3)

  expect_identical(colnames(named), rev(levels(as.factor(ids))))
  expect_equal(
    unname(named),
    own_draws(colnames(named), 1:3),
    tolerance = 1e-12
  )
})

test_that("draw columns without names are in the order of the factor levels", {
  ids <- as.character(1:10)
  model <- stan_named_blrm(rep(ids, each = 2))
  model$gamma_draws <- unname(model$gamma_draws)
  newdata <- data.frame(id = c("10", "2"), stringsAsFactors = FALSE)

  u <- mostr:::blrm_random_effect_matrix(model, newdata, "id", 1:3)

  expect_equal(unname(u), own_draws(newdata$id, 1:3), tolerance = 1e-12)
})

test_that("a fit that does not record the cluster variable is looked up by `id_var`", {
  model <- stan_named_blrm(rep(2:4, each = 2))
  model$clusterInfo$name <- NULL

  u <- mostr:::blrm_random_effect_matrix(model, data.frame(id = 4:3), "id", 1:3)

  expect_equal(unname(u), own_draws(4:3, 1:3), tolerance = 1e-12)
  expect_error(
    mostr:::blrm_random_effect_matrix(model, data.frame(site = 4), "id", 1:3),
    "`id`",
    fixed = TRUE
  )
})

test_that("a missing cluster(site) value matters only with include_re = TRUE", {
  local_mocked_bindings(blrm_design_matrix = fake_blrm_design("tx"))
  model <- stan_named_blrm(c("2", "-1", "0.5", "2"), cluster_name = "site")
  # The formula as rmsb stores it, without the `cluster()` term. mostr reads
  # the model's columns from it and checks them for missing values.
  model$sformula <- y ~ tx
  newdata <- data.frame(
    id = 1:2,
    site = c("2", NA),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1,
    stringsAsFactors = FALSE
  )
  expect_contains(mostr:::markov_model_variables(model), "site")

  # The site only looks up random intercepts, so without them it is not used.
  without_re <- sops(model, newdata = newdata, times = 1, n_draws = NULL)
  expect_false(anyNA(without_re$estimate))

  expect_error(
    sops(
      model,
      newdata = newdata,
      times = 1,
      include_re = TRUE,
      n_draws = NULL
    ),
    "not present",
    fixed = TRUE
  )
})

test_that("a factor cluster variable with unused last levels keeps its numbering", {
  # rmsb numbers clusters by factor level and gives Stan the largest number
  # that occurs, so the unused last level 9 has no draws.
  cluster <- factor(c(6, 5, 6, 5), levels = c(5, 6, 9))
  model <- stan_named_blrm(cluster)
  model$gamma_draws <- model$gamma_draws[, 1:2]

  u <- mostr:::blrm_random_effect_matrix(model, data.frame(id = 6:5), "id", 1:3)

  expect_equal(unname(u), own_draws(6:5, 1:3), tolerance = 1e-12)
  expect_error(
    mostr:::blrm_random_effect_matrix(model, data.frame(id = 9), "id", 1:3),
    "not present",
    fixed = TRUE
  )
})
