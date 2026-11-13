describe("avg_sops() and inferences() pipeline", {
  add_patient_age <- function(data) {
    ids <- sort(unique(data$id))
    age_lookup <- data.frame(
      id = ids,
      age = 45 + ((ids * 11) %% 31)
    )

    data$age <- NULL
    merge(data, age_lookup, by = "id", all.x = TRUE, sort = FALSE)
  }

  add_time_basis <- function(data) {
    spl <- rms::rcs(data$time, 4)
    data$time_lin <- as.vector(spl[, 1])
    data$time_nlin_1 <- as.vector(spl[, 2])
    data$time_nlin_2 <- as.vector(spl[, 3])
    data
  }

  # The spline columns at each visit, one row per visit in time order, as the
  # readable loop sets them for a model fitted to the columns themselves.
  get_time_basis_table <- function(data) {
    cols <- c("time", "time_lin", "time_nlin_1", "time_nlin_2")
    out <- unique(data[, cols, drop = FALSE])
    out <- out[order(out$time), , drop = FALSE]
    rownames(out) <- NULL
    out[, c("time_lin", "time_nlin_1", "time_nlin_2"), drop = FALSE]
  }

  build_markov_pipeline_case <- function() {
    raw_data <- make_markov_trajectories(
      n_patients = 70,
      follow_up_time = 8,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 4041,
      treatment_effect = 0.18
    )
    data <- mostr::prepare_markov_data(raw_data)
    data <- add_patient_age(data)

    list(
      data = data,
      baseline = data[data$time == 1, , drop = FALSE],
      y_levels = 1:6
    )
  }

  build_numeric_yprev_case <- function() {
    raw_data <- make_markov_trajectories(
      n_patients = 90,
      follow_up_time = 7,
      n_states = 10,
      thresholds = seq(-3, 3, length.out = 9),
      allowed_start_state = as.integer(2:9),
      absorbing_state = 10,
      seed = 4042
    )
    data <- mostr::prepare_markov_data(
      raw_data,
      absorbing_state = 10,
      factor_previous = FALSE
    )

    list(
      data = data,
      baseline = data[data$time == 1, , drop = FALSE],
      y_levels = 1:10
    )
  }

  # The time spline is written inline. It replaces precomputed columns of
  # `rms::rcs(data$time, 4)` that were passed through the removed
  # `time_covariates`; rms computes the same knots from the same data, so the
  # snapshots below did not change.
  fit_pipeline_model <- function(data) {
    suppressWarnings(
      mostr::mostm(
        ordered(y) ~ rms::rcs(time, 4) * tx + yprev + age,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        absorb = 6,
        engine = "vglm"
      )
    )
  }

  fit_numeric_yprev_model <- function(data) {
    suppressWarnings(
      mostr::mostm(
        ordered(y) ~ rms::rcs(time, 4) * tx + rms::rcs(yprev, 6),
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        absorb = 10,
        engine = "vglm"
      )
    )
  }

  # The pipeline data and its fit, made once and shared by the MVN, bootstrap
  # and FWB tests.
  pipeline_case <- function() {
    grid_fixture("sops-inferences pipeline case", build_markov_pipeline_case)
  }

  pipeline_mostm_fit <- function() {
    grid_fixture("sops-inferences pipeline mostm fit", function() {
      fit_pipeline_model(pipeline_case()$data)
    })
  }

  pipeline_signature <- function(result) {
    keep <- result$time %in% c(1, 3, 5, 8) & result$state %in% c(1, 3, 6)
    out <- result[
      keep,
      c("tx", "time", "state", "estimate", "conf.low", "conf.high", "std.error")
    ]
    out <- out[order(out$tx, out$time, out$state), , drop = FALSE]
    out$state <- as.integer(as.character(out$state))
    numeric_cols <- c("estimate", "conf.low", "conf.high", "std.error")
    out[numeric_cols] <- lapply(out[numeric_cols], round, digits = 8)
    rownames(out) <- NULL
    lapply(out, unname)
  }

  draw_signature <- function(draws) {
    keep <- draws$draw_id %in%
      c(1L, 2L) &
      draws$time %in% c(1, 4, 8) &
      draws$state %in% c(1, 6)
    out <- draws[keep, c("draw_id", "tx", "time", "state", "draw")]
    out <- out[order(out$draw_id, out$tx, out$time, out$state), , drop = FALSE]
    out$state <- as.integer(as.character(out$state))
    out$draw <- round(out$draw, digits = 8)
    rownames(out) <- NULL
    lapply(out, unname)
  }

  expect_draw_specific_baseline_anchor <- function(result, n_draws, times) {
    anchors <- attr(result, "baseline_anchor_draws")
    expect_s3_class(anchors, "data.frame")
    expect_equal(sort(unique(anchors$draw_id)), seq_len(n_draws))
    sums <- stats::aggregate(
      estimate ~ draw_id + tx,
      data = anchors,
      FUN = sum
    )
    expect_equal(sums$estimate, rep(1, nrow(sums)), tolerance = 1e-12)

    interpolated <- interpolate_sops(
      result,
      time_map = stats::setNames(seq_along(times) * 7, times),
      target_times = c(0, 1, 7)
    )
    draws <- attr(interpolated, "draws")
    early <- draws[draws$time < 7, , drop = FALSE]
    expect_false(anyNA(early$estimate))
  }

  test_that("create_counterfactual_data() stacks one baseline copy per scenario", {
    baseline <- data.frame(
      id = 1:2,
      tx = c(0, 1),
      sex = c("F", "M"),
      marker = c(10, 20)
    )
    # 2x2 grid of counterfactual scenarios
    variables <- list(tx = c(0, 1), sex = c("F", "M"))
    grid <- expand.grid(variables, KEEP.OUT.ATTRS = FALSE)

    out <- mostr:::create_counterfactual_data(baseline, grid, variables)

    expect_equal(nrow(out), nrow(baseline) * nrow(grid))
    expect_equal(out$id, rep(baseline$id, times = nrow(grid)))
    expect_equal(out$marker, rep(baseline$marker, times = nrow(grid)))

    # Check that the counterfactual variables are correctly assigned in blocks
    for (i in seq_len(nrow(grid))) {
      rows <- ((i - 1) * nrow(baseline) + 1):(i * nrow(baseline))
      expect_equal(out$tx[rows], rep(grid$tx[i], nrow(baseline)))
      expect_equal(
        as.character(out$sex[rows]),
        rep(as.character(grid$sex[i]), nrow(baseline))
      )
    }
  })

  test_that("marginalize_sops_array() averages blocks, applies patient weights and normalizes them within by strata", {
    block_array <- array(
      seq_len(4 * 2 * 3) / 100,
      dim = c(4, 2, 3),
      dimnames = list(NULL, NULL, as.character(1:3))
    )
    weighted_array <- array(
      c(0.1, 0.7, 0.2, 0.2, 0.6, 0.2, 0.3, 0.4, 0.3, 0.4, 0.3, 0.3),
      dim = c(2, 2, 3)
    )
    patient_weights <- c(0.25, 0.75)
    weighted_matrix <- matrix(
      aperm(weighted_array, c(2, 3, 1)),
      nrow = 6,
      ncol = 2
    )
    weighted_args <- list(
      sops_array = weighted_array,
      grid = data.frame(tx = 0),
      times = c(1, 2),
      y_levels = 1:3,
      variables = list(tx = 0),
      n_cf = 1,
      n_each = 2
    )

    cases <- list(
      list(
        label = "unweighted counterfactual blocks",
        args = list(
          sops_array = block_array,
          grid = data.frame(tx = c(0, 1)),
          times = c(1, 2),
          y_levels = 1:3,
          variables = list(tx = c(0, 1)),
          n_cf = 2,
          n_each = 2
        ),
        expected = list(
          estimate = c(
            as.vector(apply(block_array[1:2, , , drop = FALSE], c(2, 3), mean)),
            as.vector(apply(block_array[3:4, , , drop = FALSE], c(2, 3), mean))
          ),
          tx = rep(c(0, 1), each = 6)
        )
      ),
      list(
        label = "patient weights",
        args = c(weighted_args, list(weights = patient_weights)),
        expected = list(
          estimate = as.vector(weighted_matrix %*% patient_weights)
        )
      ),
      list(
        label = "negative patient weights",
        args = c(weighted_args, list(weights = c(1, -1))),
        error = "non-negative"
      ),
      list(
        label = "weights normalized within by strata",
        args = list(
          sops_array = array(c(0.2, 0.8, 0.1, 0.5), dim = c(4, 1, 1)),
          grid = data.frame(tx = 0),
          times = 1,
          y_levels = 1,
          variables = list(tx = 0),
          n_cf = 1,
          n_each = 4,
          weights = c(1, 3, 2, 6),
          by = "subgroup",
          newdata = data.frame(
            id = 1:4,
            tx = 0,
            subgroup = c("a", "a", "b", "b")
          )
        ),
        order_by = "subgroup",
        expected = list(estimate = c(0.65, 0.4), subgroup = c("a", "b"))
      )
    )

    for (case in cases) {
      if (!is.null(case$error)) {
        expect_error(
          do.call(mostr:::marginalize_sops_array, case$args),
          case$error,
          info = case$label
        )
        next
      }
      out <- do.call(mostr:::marginalize_sops_array, case$args)
      if (!is.null(case$order_by)) {
        out <- out[order(out[[case$order_by]]), , drop = FALSE]
      }
      for (column in names(case$expected)) {
        expect_equal(out[[column]], case$expected[[column]], info = case$label)
      }
    }
  })

  test_that("array_to_df_individual() keeps row identity, copies draw weights and rejects weight column collisions", {
    cases <- list(
      list(
        label = "row identity",
        sops_array = array(seq_len(3 * 2 * 2) / 20, dim = c(3, 2, 2)),
        times = c(1, 2),
        y_levels = 1:2,
        newdata = data.frame(
          id = c(101, 102, 103),
          rowid = c(11, 12, 13),
          tx = c(0, 0, 1)
        )
      ),
      list(
        label = "draw weights",
        sops_array = array(c(0.1, 0.5, 0.9), dim = c(3, 1, 1)),
        times = 1,
        y_levels = 1,
        newdata = data.frame(id = 1:3, rowid = 1:3),
        weights = c(1, 3, 100)
      ),
      list(
        label = "weight column collision",
        sops_array = array(c(0.1, 0.5), dim = c(2, 1, 1)),
        times = 1,
        y_levels = 1,
        newdata = data.frame(id = 1:2, rowid = 1:2, fwb_weight = c(10, 20)),
        weights = c(1, 3),
        error = "reserved for bootstrap draw weights"
      )
    )

    for (case in cases) {
      args <- list(
        sops_array = case$sops_array,
        times = case$times,
        y_levels = case$y_levels,
        newdata = case$newdata
      )
      if (!is.null(case$weights)) {
        args$weights <- case$weights
        args$weight_col <- "fwb_weight"
      }
      if (!is.null(case$error)) {
        expect_error(
          do.call(mostr:::array_to_df_individual, args),
          case$error,
          fixed = TRUE,
          info = case$label
        )
        next
      }
      out <- do.call(mostr:::array_to_df_individual, args)
      first_rows <- seq_len(nrow(case$newdata))
      expect_equal(nrow(out), length(case$sops_array), info = case$label)
      expect_equal(out$rowid[first_rows], case$newdata$rowid, info = case$label)
      expect_equal(out$estimate, as.vector(case$sops_array), info = case$label)
      if (!is.null(case$weights)) {
        expect_equal(out$fwb_weight, case$weights, info = case$label)
      }
    }
  })

  test_that("compute_ci_from_draws() computes percentile and wald intervals", {
    draws <- data.frame(
      draw_id = rep(1:4, times = 2),
      time = rep(1, 8),
      state = rep(c(1, 2), each = 4),
      estimate = c(0.1, 0.2, 0.4, 0.5, 0.5, 0.6, 0.8, 0.9)
    )

    perc <- mostr:::compute_ci_from_draws(
      draws_df = draws,
      group_cols = c("time", "state"),
      conf_level = 0.5,
      conf_type = "perc"
    )
    wald <- mostr:::compute_ci_from_draws(
      draws_df = draws,
      group_cols = c("time", "state"),
      conf_level = 0.5,
      conf_type = "wald"
    )

    expect_equal(perc$conf.low, c(0.175, 0.575))
    expect_equal(perc$conf.high, c(0.425, 0.825))
    expect_equal(
      perc$std.error,
      c(stats::sd(draws$estimate[1:4]), stats::sd(draws$estimate[5:8]))
    )

    critical <- abs(stats::qnorm(0.25))
    expected_wald_low <- c(
      mean(draws$estimate[1:4]) - critical * stats::sd(draws$estimate[1:4]),
      mean(draws$estimate[5:8]) - critical * stats::sd(draws$estimate[5:8])
    )
    expected_wald_high <- c(
      mean(draws$estimate[1:4]) + critical * stats::sd(draws$estimate[1:4]),
      mean(draws$estimate[5:8]) + critical * stats::sd(draws$estimate[5:8])
    )

    expect_equal(wald$conf.low, expected_wald_low)
    expect_equal(wald$conf.high, expected_wald_high)
    expect_error(
      mostr:::compute_ci_from_draws(
        draws,
        c("time", "state"),
        conf_level = 1
      ),
      "between 0 and 1",
      fixed = TRUE
    )
    expect_error(
      mostr:::compute_ci_from_draws(
        draws,
        c("time", "state"),
        conf_type = "other"
      ),
      "conf_type"
    )

    overall <- mostr:::compute_ci_from_draws(
      draws,
      character(),
      conf_level = 0.5
    )
    expect_named(overall, c("conf.low", "conf.high", "std.error"))
    expect_equal(nrow(overall), 1L)
  })

  test_that("lp_to_probs() keeps crossed cut-points for the caller to report", {
    # Row 1 has increasing logits, so its cut-points cross; row 2 is valid.
    eta <- matrix(
      c(
        -1,
        0,
        1,
        0.5,
        -0.5,
        -1
      ),
      nrow = 2,
      byrow = TRUE
    )

    probs <- mostr:::lp_to_probs(eta, M = 3)

    expected <- cbind(
      1 - stats::plogis(eta[, 1]),
      stats::plogis(eta[, 1]) - stats::plogis(eta[, 2]),
      stats::plogis(eta[, 2]) - stats::plogis(eta[, 3]),
      stats::plogis(eta[, 3])
    )

    expect_equal(probs, expected)
    expect_true(all(probs[2, ] >= 0))
    expect_lt(min(probs[1, ]), 0)
    expect_equal(rowSums(probs), rep(1, nrow(probs)))
    expect_error(
      mostr:::check_crossed_thresholds(probs, n_profiles = 2L, time = 1),
      "cross for 1 of the 2 profiles at time 1",
      class = "mostr_crossed_thresholds"
    )
  })

  test_that("seeded Markov avg_sops() and fast simulation inference remain stable", {
    skip_on_cran() # The pass depends on the seed.
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    case <- pipeline_case()
    model <- pipeline_mostm_fit()

    avg <- mostr::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8
    )
    expect_equal(nrow(attr(avg, "newdata_orig")), nrow(case$baseline))

    withr::local_seed(4401)
    inferred <- mostr::inferences(
      avg,
      method = "mvn",
      n_draws = 3,
      return_draws = TRUE
    )

    draws <- mostr::get_draws(inferred)

    expect_equal(attr(inferred, "method"), "mvn")
    expect_equal(attr(inferred, "n_successful"), 3L)
    expect_inference_intervals(inferred)
    expect_equal(sort(unique(draws$draw_id)), 1:3)

    expect_snapshot_value(
      list(
        result = pipeline_signature(inferred),
        draws = draw_signature(draws)
      ),
      style = "json2",
      tolerance = 1e-7,
      cran = TRUE
    )
  })

  test_that("numeric previous-state spline supports MVN and bootstrap inference", {
    skip_on_cran()
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    case <- build_numeric_yprev_case()
    model <- fit_numeric_yprev_model(case$data)
    avg <- mostr::avg_sops(
      model = model,
      variables = "tx",
      times = 1:7
    )

    cases <- list(
      list(method = "mvn", seed = 4413, n_draws = 2, return_draws = TRUE),
      list(method = "bootstrap", seed = 4414, n_draws = 2, return_draws = TRUE)
    )

    for (row in cases) {
      withr::local_seed(row$seed)
      inferred <- inferences(
        avg,
        method = row$method,
        n_draws = row$n_draws,
        return_draws = row$return_draws
      )

      expect_equal(attr(inferred, "method"), row$method, info = row$method)
      expect_equal(
        attr(inferred, "n_successful"),
        as.integer(row$n_draws),
        info = row$method
      )
      expect_inference_intervals(inferred)
      if (row$return_draws) {
        expect_equal(
          sort(unique(mostr::get_draws(inferred)$draw_id)),
          seq_len(row$n_draws),
          info = row$method
        )
      }
    }
  })

  test_that("avg_sops() supports standard bootstrap inference", {
    skip_on_cran() # The pass depends on the seed.
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- pipeline_case()
    model <- pipeline_mostm_fit()

    avg <- mostr::avg_sops(
      model = model,
      variables = "tx",
      times = 1:8
    )
    expect_equal(nrow(attr(avg, "newdata_orig")), nrow(case$baseline))
    expect_equal(nrow(attr(avg, "refit_data")), nrow(case$data))

    withr::local_seed(4402)
    inferred <- mostr::inferences(
      avg,
      method = "bootstrap",
      n_draws = 2,
      return_draws = TRUE
    )

    draws <- mostr::get_draws(inferred)

    expect_equal(attr(inferred, "method"), "bootstrap")
    expect_equal(attr(inferred, "n_boot"), 2)
    expect_equal(attr(inferred, "n_successful"), 2L)
    expect_inference_intervals(inferred)
    expect_equal(sort(unique(draws$draw_id)), 1:2)
    expect_draw_specific_baseline_anchor(inferred, 2, 1:8)

    expect_snapshot_value(
      list(
        result = pipeline_signature(inferred),
        draws = draw_signature(draws)
      ),
      style = "json2",
      tolerance = 1e-7,
      cran = TRUE
    )
  })

  test_that("avg_sops() supports fractional weighted bootstrap inference", {
    skip_on_cran()
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- pipeline_case()
    model <- pipeline_mostm_fit()

    avg <- mostr::avg_sops(
      model = model,
      variables = "tx",
      times = 1:8
    )

    withr::local_seed(4404)
    inferred <- mostr::inferences(
      avg,
      method = "fwb",
      n_draws = 2,
      return_draws = TRUE
    )

    draws <- mostr::get_draws(inferred)

    expect_equal(attr(inferred, "method"), "fwb")
    expect_equal(attr(inferred, "fwb_weight_type"), "exponential")
    expect_equal(attr(inferred, "fwb_weight_scale"), "cluster_mean_1")
    expect_equal(attr(inferred, "n_boot"), 2)
    expect_equal(attr(inferred, "n_successful"), 2L)
    expect_inference_intervals(inferred)
    expect_equal(sort(unique(draws$draw_id)), 1:2)
    expect_draw_specific_baseline_anchor(inferred, 2, 1:8)
  })

  test_that("rerunning inference replaces baseline draws", {
    skip_on_cran()
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    data <- suppressWarnings(make_test_data(
      n_patients = 28,
      follow_up_time = 5,
      seed = 2811
    ))
    model <- mostm(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      absorb = 6
    )
    avg <- avg_sops(model, variables = "tx", times = 1:3)
    previous <- inferences(
      avg,
      method = "fwb",
      n_draws = 4,
      seed = 1
    )
    expect_s3_class(attr(previous, "baseline_anchor_draws"), "data.frame")

    for (n_draws in c(4, 6)) {
      fresh <- inferences(avg, method = "mvn", n_draws = n_draws, seed = 2)
      rerun <- inferences(previous, method = "mvn", n_draws = n_draws, seed = 2)
      expect_null(attr(rerun, "baseline_anchor_draws"))
      interpolated <- interpolate_sops(
        rerun,
        time_map = c("1" = 1, "2" = 2, "3" = 3),
        target_times = 0:3
      )
      expected <- interpolate_sops(
        fresh,
        time_map = c("1" = 1, "2" = 2, "3" = 3),
        target_times = 0:3
      )
      expect_equal(interpolated$std.error, expected$std.error)
      expect_equal(interpolated$conf.low, expected$conf.low)
      expect_equal(interpolated$conf.high, expected$conf.high)
      expect_equal(
        interpolated$std.error[interpolated$time == 0],
        rep(0, sum(interpolated$time == 0))
      )
      expect_equal(
        sort(unique(get_draws(interpolated)$draw_id)),
        seq_len(n_draws)
      )
    }

    replacement <- inferences(
      previous,
      method = "fwb",
      n_draws = 2,
      seed = 3
    )
    expect_equal(
      sort(unique(attr(replacement, "baseline_anchor_draws")$draw_id)),
      1:2
    )
    for (method in c("mvn", "fwb", "delta")) {
      omitted <- inferences(
        previous,
        method = method,
        cluster = data$id,
        n_draws = 2,
        seed = 3,
        return_draws = FALSE
      )
      expect_null(attr(omitted, "draws"))
      expect_null(attr(omitted, "baseline_anchor_draws"))
    }
  })

  test_that("full-PO vglm engine fit with inline rcs matches VGAM::vglm()", {
    skip_on_cran()
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    raw_data <- make_markov_trajectories(
      n_patients = 200,
      follow_up_time = 12,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 2375679,
      treatment_effect = 0
    )
    data <- mostr::prepare_markov_data(raw_data, absorbing_state = 6)
    data <- add_patient_age(data)
    form <- ordered(y) ~ rms::rcs(time, 4) + tx + age + yprev

    ordinary_fit <- suppressWarnings(
      VGAM::vglm(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )
    markov_fit <- suppressWarnings(
      mostr::mostm(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        absorb = 6,
        engine = "vglm"
      )
    )$vglm_fit

    ordinary_fitted <- stats::predict(ordinary_fit, type = "response")
    markov_fitted <- mostr:::predict_vglm_response_markov(
      markov_fit,
      data
    )

    expect_s4_class(markov_fit, "vglm")
    expect_true(isTRUE(attr(markov_fit, "markov_vglm")))
    expect_equal(
      unname(markov_fitted),
      unname(ordinary_fitted),
      tolerance = 1e-10
    )
  })

  test_that("vglm engine fit with constrained PPO matches explicit spline basis in fit, SOPs and covariance", {
    skip_on_cran()
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    raw_data <- make_markov_trajectories(
      n_patients = 200,
      follow_up_time = 12,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 2375679,
      treatment_effect = 0
    )
    data <- mostr::prepare_markov_data(raw_data, absorbing_state = 6)
    data <- add_patient_age(data)
    data <- add_time_basis(data)

    explicit_basis <- unname(as.matrix(data[, c(
      "time_lin",
      "time_nlin_1",
      "time_nlin_2"
    )]))
    inline_basis_raw <- rms::rcs(data$time, 4)
    inline_basis <- matrix(
      as.numeric(inline_basis_raw),
      nrow = nrow(inline_basis_raw),
      ncol = ncol(inline_basis_raw)
    )
    expect_equal(inline_basis, explicit_basis, tolerance = 1e-14)

    baseline <- data[data$time == 1, , drop = FALSE]
    time_basis <- get_time_basis_table(data)
    y_levels <- 1:6
    n_thresholds <- length(y_levels) - 1

    explicit_constraints <- list(
      "(Intercept)" = diag(n_thresholds),
      time_lin = cbind(
        PO_effect = rep(1, n_thresholds),
        linear_deviation = seq_len(n_thresholds) - 1
      ),
      time_nlin_1 = cbind(PO_effect = rep(1, n_thresholds)),
      time_nlin_2 = cbind(PO_effect = rep(1, n_thresholds)),
      tx = cbind(PO_effect = rep(1, n_thresholds)),
      age = cbind(PO_effect = rep(1, n_thresholds)),
      yprev = cbind(PO_effect = rep(1, n_thresholds))
    )

    explicit_fit <- suppressWarnings(
      mostr::mostm(
        ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + age + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
        data = data,
        id_var = "id",
        absorb = 6,
        constraints = explicit_constraints,
        engine = "vglm"
      )
    )$vglm_fit
    explicit_fit@call$constraints <- explicit_constraints

    inline_form <- ordered(y) ~ rms::rcs(time, 4) + tx + age + yprev
    inline_po <- suppressWarnings(
      mostr::mostm(
        inline_form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id",
        absorb = 6,
        engine = "vglm"
      )
    )$vglm_fit

    # Column-level constraints on the inline spline: a partial proportional
    # odds effect for its linear column only.
    inline_constraints <- inline_po@constraints
    inline_time_terms <- grep(
      "^rms::rcs\\(time, 4\\)",
      names(inline_constraints),
      value = TRUE
    )
    expect_length(inline_time_terms, 3)
    inline_linear_time_term <- inline_time_terms[[1]]
    expect_match(inline_linear_time_term, "time$")
    yprev_terms <- grep("^yprev", names(inline_constraints), value = TRUE)
    inline_constraints[yprev_terms] <- NULL
    inline_constraints[["yprev"]] <- cbind(PO_effect = rep(1, n_thresholds))
    inline_constraints[[inline_linear_time_term]] <- cbind(
      PO_effect = rep(1, n_thresholds),
      linear_deviation = seq_len(n_thresholds) - 1
    )

    inline_fit <- suppressWarnings(
      mostr::mostm(
        inline_form,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
        data = data,
        id_var = "id",
        absorb = 6,
        constraints = inline_constraints,
        engine = "vglm"
      )
    )$vglm_fit
    inline_fit@call$constraints <- inline_constraints

    explicit_fitted <- stats::predict(explicit_fit, type = "response")
    inline_fitted <- mostr:::predict_vglm_response_markov(
      inline_fit,
      data
    )
    expect_equal(dim(inline_fitted), dim(explicit_fitted))
    expect_equal(
      unname(inline_fitted),
      unname(explicit_fitted),
      tolerance = 1e-10
    )

    # The SOPs of the inline fit, which mostr evaluates at each visit, match
    # the readable loop on the fit to the precomputed columns, with the
    # columns set to each visit's spline values. A few profiles under each
    # treatment keep the loop quick.
    profiles <- baseline[1:3, , drop = FALSE]
    profiles <- rbind(transform(profiles, tx = 0), transform(profiles, tx = 1))
    inline_sops <- mostr::sops(inline_fit, newdata = profiles, times = 1:12)
    basis_table <- time_basis[1:12, , drop = FALSE]
    reference <- lapply(seq_len(nrow(profiles)), function(i) {
      slow_first_order_sops(
        slow_predictor(explicit_fit),
        profiles[i, , drop = FALSE],
        times = 1:12,
        states = y_levels,
        absorb = 6,
        visit_columns = basis_table
      )
    })
    expect_lt(
      grid_largest_gap(inline_sops, reference, inline_sops$rowid, 1:12),
      1e-10
    )

    # The patient-clustered robust covariance does not depend on how the
    # spline basis was written, so MVN inference on the inline fit uses the
    # same covariance as the fit to the precomputed columns.
    inline_robust <- mostr::robcov_vglm(inline_fit, cluster = data$id)
    explicit_robust <- mostr::robcov_vglm(explicit_fit, cluster = data$id)
    expect_equal(inline_robust$bread_type, "observed")
    expect_true(all(is.finite(inline_robust$se)))
    expect_equal(
      unname(stats::coef(inline_robust)),
      unname(stats::coef(explicit_robust)),
      tolerance = 1e-8
    )
    expect_equal(
      unname(stats::vcov(inline_robust)),
      unname(stats::vcov(explicit_robust)),
      tolerance = 1e-8
    )

    inline_inferred <- mostr::avg_sops(
      inline_robust,
      newdata = baseline,
      variables = "tx",
      times = 1:12
    ) |>
      mostr::inferences(method = "mvn", n_draws = 3, workers = 1, seed = 42)
    expect_s3_class(inline_inferred, "markov_avg_sops")
    for (column in c("estimate", "conf.low", "conf.high", "std.error")) {
      expect_true(all(is.finite(inline_inferred[[column]])), info = column)
    }
  })
})
