# Random intercepts (`cluster()`) for blrm fits with unobserved states: the
# cluster numbers and the data of the Stan program, built without compiling
# or sampling anything. The sampled density is checked in
# test-slow-missing-states-blrm.R.

test_that("clusters are numbered as rmsb numbers them, over all fitted rows", {
  # rmsb: as.integer(as.factor(cluster)) over the fitted rows; here the
  # fully observed transitions and the stretches, so that a patient seen
  # only in stretches (id 33 below) gets a random intercept too.
  values <- c(10, 10, 2, 2, 33, 33, 5)
  codes <- mostr:::markov_path_cluster_codes(values, used = c(1, 3, 4, 6))
  expect_identical(codes$levels, c("2", "10", "33"))
  expect_identical(codes$n, 3L)
  expect_identical(codes$code, c(2L, 2L, 1L, 1L, 3L, 3L, NA))
  # A factor keeps its levels, and there are as many clusters as the
  # largest number used, as in rmsb (here level "c" is never used).
  site <- factor(c("b", "a", "b"), levels = c("a", "b", "c"))
  codes <- mostr:::markov_path_cluster_codes(site, used = 1:3)
  expect_identical(codes$levels, c("a", "b"))
  expect_identical(codes$n, 2L)
  expect_identical(codes$code, c(2L, 1L, 2L))
  expect_identical(
    codes$code,
    as.integer(as.factor(site)),
    label = "rmsb's numbering"
  )
})

test_that("a cluster must be known and the same on all rows of a patient", {
  check <- function(values, rows = seq_along(values)) {
    mostr:::markov_check_path_clusters(
      list(name = "site", values = values),
      id = c(1, 1, 2, 2),
      rows = rows,
      id_var = "id"
    )
  }
  expect_silent(check(c("a", "a", "a", "a")))
  expect_silent(check(c(1, 1, 2, 2)))
  expect_error(
    check(c("a", "b", "a", "a")),
    "changes within a patient (first: `id` = 1)",
    class = "mostr_error",
    fixed = TRUE
  )
  expect_error(
    check(c("a", NA, "b", "b")),
    "`cluster(site)` is missing on 1 row used",
    class = "mostr_error",
    fixed = TRUE
  )
  # Rows the fit does not use are not checked.
  expect_silent(check(c("a", NA, "b", "b"), rows = c(1, 3, 4)))
})

# Two patients with the same predictors at every visit: rows 1-2 are the
# fully observed transitions (one per patient, both into state 2), and the
# plan has one stretch per patient.
cluster_stan_inputs <- function() {
  x_setup <- cbind(
    yprev2 = c(1, 1, 0, 0, 1),
    time = c(2, 2, 3, 4, 1)
  )
  data <- cl_two_patients()
  plan <- cl_plan(data, design = function(rows) {
    previous <- as.integer(as.character(rows$yprev))
    cbind(previous == 2, rows$time) * 1
  })
  list(x_setup = x_setup, plan = plan, data = data)
}

cl_two_patients <- function() {
  data.frame(
    id = rep(c(1L, 2L), each = 4),
    time = rep(1:4, 2),
    tx = 0,
    yprev = factor(c(2, 2, NA, NA, 2, 2, NA, NA), levels = 1:3),
    y = c(2L, NA, NA, 3L, 2L, NA, NA, 3L)
  )
}

cl_plan <- function(data, design, cluster = NULL) {
  ranges <- ms_row_ranges(data)
  ranges$lo[is.na(data$y)] <- NA
  ranges$hi[is.na(data$y)] <- NA
  mostr:::markov_path_plan(
    data = data,
    ranges = list(labels = as.character(1:4), lo = ranges$lo, hi = ranges$hi),
    absorb = 4L,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    weights = NULL,
    design = design,
    cluster = cluster
  )
}

test_that("the Stan data keep transitions and stretches of clusters apart", {
  skip_on_cran() # Slow (about 2 s): loads rmsb for `selectedQr()`.
  skip_if_not_installed("rmsb")
  inputs <- cluster_stan_inputs()
  keep <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  build <- function(clusters, plan = inputs$plan, ...) {
    mostr:::markov_blrm_path_stan_data(
      x_setup = inputs$x_setup,
      keep = keep,
      observed_state = c(2L, 2L),
      plan = plan,
      n_states = 4L,
      conc = 0.5,
      clusters = clusters,
      ...
    )$data
  }
  plain <- build(NULL)
  # Without random intercepts the two equal transitions are one cell of
  # weight 2 and the two equal stretches one stretch of weight 2.
  expect_identical(plain$n_obs, 1L)
  expect_identical(as.numeric(plain$w_obs), 2)
  expect_identical(plain$S, 1L)
  expect_identical(plain$Nc, 0L)
  expect_length(plain$obs_cluster, 0L)
  expect_length(plain$stretch_cluster, 0L)
  expect_length(plain$rsdmean, 0L)
  expect_length(plain$rsdsd, 0L)

  # One cluster per patient: nothing is merged across patients.
  by_patient <- cl_plan(
    inputs$data,
    design = function(rows) cbind(as.integer(as.character(rows$yprev)) == 2, rows$time) * 1,
    cluster = rep(c(1L, 2L), each = 4)
  )
  separate <- build(list(n = 2L, obs = c(1L, 2L)), plan = by_patient)
  expect_identical(separate$n_obs, 2L)
  expect_identical(as.numeric(separate$w_obs), c(1, 1))
  expect_identical(as.integer(separate$obs_cluster), c(1L, 2L))
  expect_identical(separate$S, 2L)
  expect_identical(as.integer(separate$stretch_cluster), c(1L, 2L))
  expect_identical(separate$Nc, 2L)
  # rmsb's prior on the SD of the random intercepts, with its defaults.
  expect_identical(separate$psigma, 1L)
  expect_identical(as.numeric(separate$rsdmean), 0)
  expect_identical(as.numeric(separate$rsdsd), 1)
  exponential <- build(
    list(n = 2L, obs = c(1L, 2L)),
    plan = by_patient,
    psigma = 2,
    rsdmean = 0.5
  )
  expect_identical(exponential$psigma, 2L)
  expect_identical(as.numeric(exponential$rsdmean), 0.5)
  expect_length(exponential$rsdsd, 0L)
  # The design data do not depend on the clusters.
  for (field in c("p", "k", "x_obs", "G", "x_gap", "x_prev", "conc")) {
    expect_identical(separate[[field]][seq_along(plain[[field]])], plain[[field]], info = field)
  }

  # Both patients at one site: merged as without random intercepts.
  same_site <- build(list(n = 1L, obs = c(1L, 1L)), plan = cl_plan(
    inputs$data,
    design = function(rows) cbind(as.integer(as.character(rows$yprev)) == 2, rows$time) * 1,
    cluster = rep(1L, 8)
  ))
  expect_identical(same_site$n_obs, 1L)
  expect_identical(as.numeric(same_site$w_obs), 2)
  expect_identical(same_site$S, 1L)
  expect_identical(as.integer(same_site$obs_cluster), 1L)
  expect_identical(as.integer(same_site$stretch_cluster), 1L)
})
