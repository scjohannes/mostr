make_avg_sops_object <- function(model, newdata, by = NULL) {
  out <- expand.grid(
    time = 1:2,
    state = 1:2,
    tx = c(0, 1),
    KEEP.OUT.ATTRS = FALSE
  )
  if (!is.null(by)) {
    out[[by]] <- "a"
  }
  out$estimate <- 0.5
  class(out) <- c("markov_avg_sops", class(out))
  attr(out, "model") <- model
  attr(out, "newdata_orig") <- newdata
  attr(out, "refit_data") <- newdata
  attr(out, "call_args") <- list(times = 1:2)
  attr(out, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = by,
    times = 1:2,
    id_var = "id"
  )
  attr(out, "time_var") <- "time"
  attr(out, "p_var") <- "yprev"
  attr(out, "y_levels") <- 1:2
  attr(out, "absorb") <- NULL
  attr(out, "time_covariates") <- NULL
  out
}

make_individual_sops_object <- function(model, newdata) {
  out <- expand.grid(
    rowid = seq_len(nrow(newdata)),
    time = 1:2,
    state = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  out$estimate <- 0.5
  class(out) <- c("markov_sops", class(out))
  attr(out, "model") <- model
  attr(out, "newdata_orig") <- newdata
  attr(out, "call_args") <- list(times = 1:2)
  attr(out, "time_var") <- "time"
  attr(out, "p_var") <- "yprev"
  attr(out, "y_levels") <- 1:2
  attr(out, "absorb") <- NULL
  attr(out, "time_covariates") <- NULL
  out
}

test_that("inferences_simulation() validates stored models, simulation engines and custom vcov objects", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "mock_model")
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  object <- make_avg_sops_object(model, newdata)
  no_model <- object
  attr(no_model, "model") <- NULL

  cases <- list(
    list(
      label = "no stored model",
      object = no_model,
      engine = "mvn",
      vcov = NULL,
      error = "Model not stored"
    ),
    list(
      label = "unknown engine",
      object = object,
      engine = "other",
      vcov = NULL,
      error = "Unknown simulation engine"
    ),
    list(
      label = "vcov given as a data frame",
      object = object,
      engine = "mvn",
      vcov = data.frame(a = c(1, 0), b = c(0, 1)),
      error = "`vcov` must be a matrix."
    )
  )

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    for (case in cases) {
      expect_error(
        mostr:::inferences_simulation(
          case$object,
          engine = case$engine,
          score_weight_dist = "exponential",
          n_draws = 2,
          vcov = case$vcov,
          cluster = NULL,
          conf_level = 0.95,
          conf_type = "perc",
          workers = NULL,
          return_draws = FALSE
        ),
        case$error,
        fixed = TRUE,
        info = case$label
      )
    }
  )
})

test_that("inferences() validates the homogenized method interface", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    yprev = 1,
    tx = 0
  )
  object <- make_avg_sops_object(model, newdata)

  expect_error(inferences(object, method = "other", n_draws = 1), "arg")
  expect_false("engine" %in% names(formals(inferences)))
  expect_false("score_weight_dist" %in% names(formals(inferences)))
  expect_false("..." %in% names(formals(inferences)))
  # Obsolete and misspelled arguments are not silently absorbed.
  expect_error(
    inferences(object, n_boot = 2),
    "unused argument (n_boot = 2)",
    fixed = TRUE
  )

  seen_engine <- NULL
  with_mocked_bindings(
    inferences_bootstrap = function(object, engine, ...) {
      seen_engine <<- engine
      attr(object, "engine") <- engine
      object
    },
    {
      out <- inferences(
        object,
        method = "fwb",
        n_draws = 1
      )
    }
  )

  expect_equal(seen_engine, "fwb")
  expect_equal(attr(out, "method"), "fwb")
})

test_that("inferences() reports repeated R fallbacks only once", {
  object <- structure(
    data.frame(estimate = 0.5),
    class = c("markov_sops", "data.frame")
  )

  with_mocked_bindings(
    inferences_impl = function(...) {
      for (i in seq_len(1000L)) {
        notify_sop_reference_fallback(paste("draw", i))
      }
      object
    },
    {
      messages <- testthat::capture_messages(inferences(object))
    }
  )

  expect_length(messages, 1L)
  expect_match(messages, "Compiled C\\+\\+ SOP calculations were not used")
})

test_that("inferences_simulation() falls back when fast-path setup fails", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "vglm")
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  object <- make_avg_sops_object(model, newdata)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2),
    validate_coef_vcov = function(beta, Sigma, arg = "vcov") Sigma,
    compile_sop_execution_plan = function(...) stop("build failed"),
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      withr::local_seed(1)
      expect_message(
        out <- mostr:::inferences_simulation(
          object,
          engine = "mvn",
          score_weight_dist = "exponential",
          n_draws = 2,
          vcov = NULL,
          cluster = NULL,
          conf_level = 0.8,
          conf_type = "wald",
          workers = NULL,
          return_draws = TRUE
        ),
        "Compiled C\\+\\+ SOP calculations were not used"
      )
    }
  )

  expect_s3_class(out, "markov_avg_sops")
  expect_equal(attr(out, "n_successful"), 2L)
  expect_equal(attr(out, "conf_type"), "wald")
  expect_s3_class(attr(out, "draws"), "data.frame")
})

test_that("inferences_simulation() errors when all fast-path or slow-path draws fail", {
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  shared_mocks <- list(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2),
    validate_coef_vcov = function(beta, Sigma, arg = "vcov") Sigma
  )

  cases <- list(
    # A vglm fit takes the compiled path, where every plan run fails.
    list(
      label = "fast path",
      class = "vglm",
      mocks = list(
        compile_sop_execution_plan = function(...) {
          structure(list(ok = TRUE), class = "markov_sop_exec_plan")
        },
        get_effective_coefs = function(model, beta = NULL) {
          matrix(c(0, 0), nrow = 1)
        },
        run_sop_execution_plan = function(...) stop("run failed")
      ),
      seed = 2,
      vcov = NULL,
      warning = "markov_msm_run failed"
    ),
    # Another model class takes the R path, where every SOP call fails.
    list(
      label = "slow path",
      class = "mock_model",
      mocks = list(
        set_coef = function(model, new_coefs) model,
        soprob_markov = function(...) stop("sop failed")
      ),
      seed = 4,
      vcov = diag(2),
      warning = "soprob_markov failed"
    )
  )

  run_case <- function(case) {
    model <- structure(list(coefficients = c(a = 0, b = 0)), class = case$class)
    object <- make_avg_sops_object(model, newdata)
    local_mocked_bindings(!!!c(shared_mocks, case$mocks))
    withr::local_seed(case$seed)
    expect_warning(
      expect_error(
        mostr:::inferences_simulation(
          object,
          engine = "mvn",
          score_weight_dist = "exponential",
          n_draws = 1,
          vcov = case$vcov,
          cluster = NULL,
          conf_level = 0.95,
          conf_type = "perc",
          workers = NULL,
          return_draws = FALSE
        ),
        "All simulation draws failed",
        fixed = TRUE,
        info = case$label
      ),
      case$warning,
      info = case$label
    )
  }
  for (case in cases) {
    run_case(case)
  }
})

test_that("score-bootstrap and FWB draw weights are normalized within by groups", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  make_by_group_object <- function(newdata, refit_data = NULL) {
    object <- expand.grid(
      time = 1,
      state = 1,
      tx = 0,
      grp = c("a", "b"),
      KEEP.OUT.ATTRS = FALSE
    )
    object$estimate <- 0.5
    class(object) <- c("markov_avg_sops", class(object))
    attr(object, "model") <- model
    attr(object, "newdata_orig") <- newdata
    attr(object, "refit_data") <- refit_data
    attr(object, "call_args") <- list(times = 1)
    attr(object, "avg_args") <- list(
      variables = list(tx = 0),
      by = "grp",
      times = 1,
      id_var = "id"
    )
    attr(object, "time_var") <- "time"
    attr(object, "p_var") <- "yprev"
    attr(object, "y_levels") <- 1
    attr(object, "absorb") <- NULL
    attr(object, "time_covariates") <- NULL
    object
  }
  refit_data <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, times = 4),
    y = factor(rep(c(1, 2), times = 4), levels = 1:2),
    yprev = factor(1, levels = 1:2),
    tx = 0,
    grp = rep(c("a", "a", "b", "b"), each = 2)
  )
  # Patient weights 1, 3 (group a) and 2, 6 (group b) give the group means
  # (0.2 + 3 * 0.8) / 4 = 0.65 and (2 * 0.1 + 6 * 0.5) / 8 = 0.4.
  sops <- function(...) array(c(0.2, 0.8, 0.1, 0.5), dim = c(4, 1, 1))

  cases <- list(
    list(
      label = "score bootstrap",
      object = make_by_group_object(data.frame(
        id = 1:4,
        time = 1,
        yprev = 1,
        tx = 0,
        grp = c("a", "a", "b", "b")
      )),
      mocks = list(
        generate_score_bootstrap_draws = function(...) {
          list(
            beta_draws = matrix(0, nrow = 1, ncol = 1),
            baseline_weights = matrix(c(1, 3, 2, 6), nrow = 1)
          )
        },
        set_coef = function(model, new_coefs) model,
        soprob_markov = sops
      ),
      run = function(object) {
        mostr:::inferences_simulation(
          object,
          engine = "score_bootstrap",
          score_weight_dist = "exponential",
          n_draws = 1,
          vcov = NULL,
          cluster = NULL,
          conf_level = 0.8,
          conf_type = "perc",
          workers = NULL,
          return_draws = TRUE
        )
      },
      weight_col = "score_weight",
      omission_reason = "averaged_sops"
    ),
    list(
      label = "FWB",
      object = make_by_group_object(refit_data, refit_data),
      mocks = list(
        generate_fwb_bootstrap_weights = function(data, id_var, n_boot) {
          list(data.frame(
            original_id = 1:4,
            fwb_weight = c(1, 3, 2, 6),
            boot_id = 1
          ))
        },
        apply_to_fwb_bootstrap = function(
          fwb_samples,
          analysis_fn,
          data,
          id_var,
          workers,
          packages,
          globals
        ) {
          list(analysis_fn(data, fwb_samples[[1]]))
        },
        bootstrap_analysis_wrapper = function(
          boot_data,
          model,
          factor_cols,
          original_data,
          y_levels,
          absorb,
          use_coefstart,
          fit_weights
        ) {
          list(
            model = model,
            data = boot_data,
            y_levels = 1,
            absorb = NULL,
            missing_states = character(0)
          )
        },
        soprob_markov = sops
      ),
      run = function(object) {
        mostr:::inferences_bootstrap(
          object,
          engine = "fwb",
          n_boot = 1,
          workers = NULL,
          conf_level = 0.95,
          return_draws = TRUE,
          use_coefstart = FALSE
        )
      },
      weight_col = "fwb_weight"
    )
  )

  run_case <- function(case) {
    local_mocked_bindings(!!!case$mocks)
    case$run(case$object)
  }
  for (case in cases) {
    out <- run_case(case)
    draws <- attr(out, "draws")
    draws <- draws[order(draws$grp), , drop = FALSE]

    expect_equal(draws$estimate, c(0.65, 0.4), info = case$label)
    expect_false(case$weight_col %in% names(draws), info = case$label)
    if (!is.null(case$omission_reason)) {
      expect_equal(
        attr(out, "draw_weight_omission_reason"),
        case$omission_reason,
        info = case$label
      )
    }
  }
})

test_that("score-bootstrap and FWB individual sops attach each patient's draw weight", {
  individual_sops <- function(...) {
    array(c(0.8, 0.7, 0.2, 0.3, 0.6, 0.5, 0.4, 0.5), dim = c(2, 2, 2))
  }

  score_object <- make_individual_sops_object(
    structure(list(coefficients = c(a = 0, b = 0)), class = "robcov_vglm"),
    data.frame(id = 1:2, yprev = c(1, 2), time = 1)
  )
  attr(score_object, "id_var") <- "id"
  attr(score_object, "newdata_supplied") <- FALSE

  prediction_data <- data.frame(id = 1:2, time = 1, yprev = 1)
  fwb_object <- make_individual_sops_object(
    structure(list(coefficients = c(a = 0)), class = "mock_model"),
    prediction_data
  )
  attr(fwb_object, "newdata_pred") <- prediction_data
  attr(fwb_object, "refit_data") <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:2),
    yprev = factor(c(1, 1, 1, 1), levels = 1:2)
  )
  attr(fwb_object, "id_var") <- "id"
  attr(fwb_object, "newdata_supplied") <- FALSE

  cases <- list(
    # The R path of the score bootstrap, with two draws of patient weights.
    list(
      label = "score bootstrap",
      object = score_object,
      mocks = list(
        generate_score_bootstrap_draws = function(...) {
          list(
            beta_draws = matrix(c(0, 0, 0, 0), nrow = 2, byrow = TRUE),
            baseline_weights = matrix(
              c(0.7, 0.3, 0.4, 0.6),
              nrow = 2,
              byrow = TRUE
            )
          )
        },
        set_coef = function(model, new_coefs) model,
        soprob_markov = individual_sops
      ),
      run = function(object) {
        mostr:::inferences_simulation(
          object,
          engine = "score_bootstrap",
          score_weight_dist = "exponential",
          n_draws = 2,
          vcov = NULL,
          cluster = NULL,
          conf_level = 0.8,
          conf_type = "perc",
          workers = NULL,
          return_draws = TRUE
        )
      },
      method = "score_bootstrap",
      weight_col = "score_weight",
      first_weights = c(0.7, 0.3),
      weights = c(0.3, 0.4, 0.6, 0.7)
    ),
    # The FWB refit on the stored data, with one draw of cluster weights.
    list(
      label = "FWB",
      object = fwb_object,
      mocks = list(
        generate_fwb_bootstrap_weights = function(data, id_var, n_boot) {
          list(data.frame(
            original_id = c(1, 2),
            fwb_weight = c(0.5, 1.5),
            boot_id = 1
          ))
        },
        bootstrap_analysis_wrapper = function(
          boot_data,
          model,
          factor_cols,
          original_data,
          y_levels,
          absorb,
          use_coefstart,
          fit_weights
        ) {
          list(
            model = model,
            data = boot_data,
            y_levels = 1:2,
            absorb = NULL,
            missing_states = character(0)
          )
        },
        soprob_markov = individual_sops
      ),
      run = function(object) {
        mostr:::inferences_bootstrap_sops_fwb(
          object,
          n_boot = 1,
          workers = NULL,
          conf_level = 0.95,
          return_draws = TRUE,
          use_coefstart = FALSE
        )
      },
      method = "fwb",
      weight_col = "fwb_weight",
      first_weights = c(0.5, 1.5),
      weights = c(0.5, 1.5)
    )
  )

  run_case <- function(case) {
    local_mocked_bindings(!!!case$mocks)
    case$run(case$object)
  }
  for (case in cases) {
    out <- run_case(case)
    draws <- attr(out, "draws")
    public_draws <- mostr::get_draws(out)

    expect_s3_class(out, "markov_sops")
    expect_equal(attr(out, "method"), case$method, info = case$label)
    expect_equal(attr(out, "draw_weights_attached"), TRUE, info = case$label)
    expect_equal(attr(out, "draw_weight_col"), case$weight_col, info = case$label)
    expect_equal(
      draws[[case$weight_col]][1:2],
      case$first_weights,
      info = case$label
    )
    expect_equal(
      sort(unique(draws[[case$weight_col]])),
      case$weights,
      info = case$label
    )
    expect_equal(
      sort(unique(public_draws[[case$weight_col]])),
      case$weights,
      info = case$label
    )
    expect_equal(
      public_draws$estimate,
      rep(0.5, nrow(public_draws)),
      info = case$label
    )
  }
})

test_that("inferences_simulation() rejects reserved score weight columns", {
  model <- structure(
    list(coefficients = c(a = 0, b = 0)),
    class = "robcov_vglm"
  )
  newdata <- data.frame(
    id = 1:2,
    yprev = c(1, 2),
    time = 1,
    score_weight = c(10, 20)
  )
  object <- make_individual_sops_object(model, newdata)
  attr(object, "id_var") <- "id"
  attr(object, "newdata_supplied") <- FALSE

  with_mocked_bindings(
    generate_score_bootstrap_draws = function(...) {
      list(
        beta_draws = matrix(c(0, 0), nrow = 1),
        baseline_weights = matrix(c(0.7, 0.3), nrow = 1)
      )
    },
    {
      expect_error(
        mostr:::inferences_simulation(
          object,
          engine = "score_bootstrap",
          score_weight_dist = "exponential",
          n_draws = 1,
          vcov = NULL,
          cluster = NULL,
          conf_level = 0.8,
          conf_type = "perc",
          workers = NULL,
          return_draws = TRUE
        ),
        "reserved for bootstrap draw weights",
        fixed = TRUE
      )
    }
  )
})

test_that("soprob_markov() validates model classes and propagates Bayesian draws", {
  expect_error(
    soprob_markov(
      model = structure(list(), class = "other"),
      newdata = data.frame(yprev = 1),
      times = 1,
      y_levels = 1
    ),
    "Object class not supported",
    fixed = TRUE
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    {
      expect_error(
        soprob_markov(
          model = structure(list(vglm_fit = NULL), class = "robcov_vglm"),
          newdata = data.frame(yprev = 1),
          times = 1,
          y_levels = 1:2
        ),
        "does not contain the original vglm fit",
        fixed = TRUE
      )
      expect_error(
        soprob_markov(
          model = structure(list(), class = "vglm"),
          newdata = data.frame(x = 1),
          times = 1,
          y_levels = 1:2
        ),
        "Previous-state variable",
        fixed = TRUE
      )
    }
  )

  draws <- matrix(stats::qlogis(0.4), nrow = 2, ncol = 1)
  colnames(draws) <- "y>=2"
  blrm <- fake_blrm(draws, cluster = FALSE, class = "blrm")

  with_mocked_bindings(
    blrm_design_matrix = fake_blrm_design(character()),
    {
      out <- soprob_markov(
        model = blrm,
        newdata = data.frame(
          id = 1:2,
          yprev = factor(c(1, 1)),
          gap_var = 0,
          z = 0
        ),
        times = 1:2,
        y_levels = 1:2,
        absorb = 2,
        gap_var = "gap_var",
        time_covariates = data.frame(z = c(1, 2))
      )
    }
  )

  expect_equal(dim(out), c(2L, 2L, 2L, 2L))
  expect_equal(unname(out[,, 2, 2]), array(0.64, dim = c(2, 2)))
})

test_that("inferences_bootstrap() validates inputs and records callback outcomes", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:3),
    yprev = factor(c(1, 1, 1, 1), levels = 1:3),
    tx = c(0, 0, 1, 1),
    grp = "a"
  )
  object <- expand.grid(
    time = 1:2,
    state = 1:3,
    tx = c(0, 1),
    grp = "a",
    KEEP.OUT.ATTRS = FALSE
  )
  object$estimate <- 0.5
  class(object) <- c("markov_avg_sops", class(object))
  attr(object, "model") <- model
  attr(object, "newdata_orig") <- newdata
  attr(object, "refit_data") <- newdata
  attr(object, "avg_args") <- list(
    variables = list(tx = c(0, 1), grp = "a"),
    by = "grp",
    times = 1:2,
    id_var = "id"
  )
  attr(object, "call_args") <- list(times = 1:2)
  attr(object, "time_var") <- "time"
  attr(object, "p_var") <- "yprev"
  attr(object, "y_levels") <- 1:3
  attr(object, "absorb") <- 3
  attr(object, "time_covariates") <- NULL

  no_data <- object
  attr(no_data, "refit_data") <- NULL
  expect_error(
    mostr:::inferences_bootstrap(
      no_data,
      engine = "standard",
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = FALSE,
      use_coefstart = FALSE
    ),
    "this model has no stored data",
    fixed = TRUE
  )

  missing_id <- object
  attr(missing_id, "refit_data") <- data.frame(time = 1:2)
  expect_error(
    mostr:::inferences_bootstrap(
      missing_id,
      engine = "standard",
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = FALSE,
      use_coefstart = FALSE
    ),
    "ID variable 'id' not found",
    fixed = TRUE
  )

  wrapper_call <- 0
  sop_call <- 0
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
      wrapper_call <<- wrapper_call + 1
      if (wrapper_call == 2) {
        return(list(
          model = NULL,
          data = transform(newdata, new_id = id),
          y_levels = 1:3,
          absorb = 3,
          missing_states = character(0)
        ))
      }
      list(
        model = model,
        data = transform(newdata, new_id = id),
        y_levels = 1:2,
        absorb = 2,
        missing_states = 3
      )
    },
    soprob_markov = function(...) {
      sop_call <<- sop_call + 1
      if (sop_call == 2) {
        stop("sop failed")
      }
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      expect_warning(
        out <- mostr:::inferences_bootstrap(
          object,
          engine = "standard",
          n_boot = 3,
          workers = NULL,
          conf_level = 0.8,
          return_draws = TRUE,
          use_coefstart = TRUE
        ),
        "soprob_markov failed"
      )
    }
  )

  expect_s3_class(out, "markov_avg_sops")
  expect_equal(attr(out, "n_successful"), 1L)
  expect_s3_class(attr(out, "draws"), "data.frame")
  expect_equal(as.character(out$grp), rep("a", nrow(out)))

  with_mocked_bindings(
    fast_group_bootstrap = function(...) {
      list(data.frame(original_id = 1, new_id = "1_1"))
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
      list(NULL)
    },
    {
      expect_error(
        mostr:::inferences_bootstrap(
          object,
          engine = "standard",
          n_boot = 1,
          workers = NULL,
          conf_level = 0.95,
          return_draws = FALSE,
          use_coefstart = FALSE
        ),
        "All bootstrap iterations failed",
        fixed = TRUE
      )
    }
  )
})

test_that("inferences_bootstrap() supports fractional weighted refits", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:2),
    yprev = factor(c(1, 1, 1, 1), levels = 1:2),
    tx = c(0, 0, 1, 1)
  )
  object <- make_avg_sops_object(model, newdata)

  fit_weights_seen <- NULL
  with_mocked_bindings(
    generate_fwb_bootstrap_weights = function(data, id_var, n_boot) {
      list(data.frame(
        original_id = c(1, 2),
        fwb_weight = c(0.5, 1.5),
        boot_id = 1
      ))
    },
    bootstrap_analysis_wrapper = function(
      boot_data,
      model,
      factor_cols,
      original_data,
      y_levels,
      absorb,
      use_coefstart,
      fit_weights
    ) {
      fit_weights_seen <<- fit_weights
      list(
        model = model,
        data = boot_data,
        y_levels = 1:2,
        absorb = NULL,
        missing_states = character(0)
      )
    },
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      out <- mostr:::inferences_bootstrap(
        object,
        engine = "fwb",
        n_boot = 1,
        workers = NULL,
        conf_level = 0.95,
        return_draws = TRUE,
        use_coefstart = FALSE
      )
    }
  )

  expect_equal(fit_weights_seen, c(0.5, 0.5, 1.5, 1.5))
  expect_equal(attr(out, "method"), "fwb")
  expect_equal(attr(out, "fwb_weight_type"), "exponential")
  expect_equal(attr(out, "fwb_weight_scale"), "cluster_mean_1")
  expect_s3_class(attr(out, "draws"), "data.frame")
})
