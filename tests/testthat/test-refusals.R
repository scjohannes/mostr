# Every situation that mostr refuses, in one table. Each row names the
# situation, makes one call that should stop, and gives a phrase the error
# message must contain (or a condition class, where the package defines one).
# The rows reuse the shared fits of helper-grid.R, so the table costs few fits.
#
# Rows are grouped by the stage that refuses: fitting with `mostm()`,
# prediction (`sops()`, `avg_sops()`, the time and comparison functions),
# inference (`inferences()` and its accessors), simulation (`sim_markov()`)
# and robust variance (`robcov_vglm()`). Within a stage, rows that need only
# cheap calls run everywhere, CRAN included; rows that need a weighted,
# penalized, second-order or partial proportional-odds fit, bootstrap refits,
# or rmsb loaded, are skipped on CRAN. IDs are F (fitting), P (prediction),
# I (inference), S (simulation) and V (robust variance); the F07 and F10 rows
# are in the blrm group.
#
# Error-text snapshots elsewhere (`_snaps/`) keep the full wording where it
# matters; this table checks that each refusal still happens. Rows marked
# "(internal)" reach code that the public functions no longer offer, such as
# the score bootstrap, through `mostr:::inferences_score_bootstrap()`.

# ---- The table machinery ----

# One row of the table. `call` is a function without arguments that makes the
# refused call. `phrase` is matched literally.
refusal <- function(id, situation, call, phrase = NULL, class = NULL) {
  list(
    id = id,
    situation = situation,
    call = call,
    phrase = phrase,
    class = class
  )
}

# Runs every row, and reports a failure by the row's ID and situation.
# Warnings and printed output on the way to the error are not checked.
check_refusals <- function(rows) {
  ids <- vapply(rows, function(row) row$id, character(1))
  expect_false(anyDuplicated(ids) > 0L, info = "refusal IDs are unique")
  for (row in rows) {
    info <- paste0(row$id, ": ", row$situation)
    condition <- tryCatch(
      {
        utils::capture.output(suppressMessages(suppressWarnings(row$call())))
        NULL
      },
      error = identity
    )
    expect_true(inherits(condition, "error"), info = paste(info, "(no error)"))
    if (!inherits(condition, "error")) {
      next
    }
    if (!is.null(row$phrase)) {
      expect_match(conditionMessage(condition), row$phrase, fixed = TRUE, info = info)
    }
    if (!is.null(row$class)) {
      expect_true(
        inherits(condition, row$class),
        info = paste(info, "(class", row$class, ")")
      )
    }
    # A mostr error names no internal function: either no call or the
    # exported function that was called.
    if (inherits(condition, "mostr_error")) {
      expect_true(
        error_call_name(condition) %in% c(NA, getNamespaceExports("mostr")),
        info = paste(info, "(call", error_call_name(condition), ")")
      )
    }
  }
}

# ---- Shared set-ups ----

refusal_orm <- grid_cell("orm")
refusal_vglm <- grid_cell("vglm")

# A public function on a cell's shared fit, with that cell's arguments
# (visits 1 to 4, state 4 absorbing). `args` replaces any of those arguments.
refusal_call <- function(fun, cell = refusal_orm, ..., args = list()) {
  defaults <- grid_public_args(cell)
  defaults[names(args)] <- args
  do.call(fun, c(list(grid_fit(cell), ...), defaults))
}

refusal_average <- function(cell = refusal_orm, ...) {
  refusal_call(avg_sops, cell, variables = grid_treatments, ...)
}

refusal_contrast <- function(cell = refusal_orm, ...) {
  refusal_call(avg_comparisons, cell, variables = grid_treatments, ...)
}

# Starting profiles in `newdata`: the stored ones of the first ten patients.
refusal_profiles <- function(cell = refusal_orm) {
  data <- grid_cell_data(cell)
  start <- data[as.character(data$time) == "1", , drop = FALSE]
  start <- start[seq_len(10), , drop = FALSE]
  rownames(start) <- NULL
  start
}

refusal_formula_orm <- y ~ time + tx + age + yprev
refusal_formula_vglm <- ordered(y) ~ time + tx + age + yprev

refusal_mostm <- function(data = grid_data(), ..., formula = refusal_formula_orm) {
  mostm(formula, data = data, id_var = "id", first_followup_time = 1, ...)
}

# A fit of changed data, made once and kept with the grid's fits.
refusal_refit <- function(name, make) {
  key <- paste("refusals |", name)
  if (is.null(grid_cache$fits[[key]])) {
    if (is.null(grid_cache$fits)) {
      grid_cache$fits <- list()
    }
    grid_cache$fits[[key]] <- suppressWarnings(make())
  }
  grid_cache$fits[[key]]
}

# A hand-made blrm fit (`fake_blrm()` in helper-mocks.R), as in
# test-sops-native-blrm-slow-loop.R: the draws are real numbers and only the
# design matrix is mocked. `cppo = NULL` with `pppo = 1` is an unconstrained
# partial proportional-odds fit; `cluster = FALSE` is a fit without random
# effects.
refusal_fake_blrm <- function(pppo = 0L, cppo = NULL, cluster = TRUE) {
  draws <- rbind(
    c(1.2, 0.1, -1.3, 0.40, 0.08, -0.60),
    c(0.9, -0.2, -1.1, 0.10, 0.03, -0.45)
  )
  colnames(draws) <- c("y>=2", "y>=3", "y>=4", "tx", "time", "yprev")
  tau <- character()
  if (pppo > 0L) {
    draws <- cbind(draws, tau_time = c(0.05, -0.10))
    tau <- "tau_time"
  }
  fake_blrm(
    draws,
    tau = tau,
    cppo = cppo,
    cluster = cluster,
    gamma_draws = matrix(0, nrow = 2, ncol = 2)
  )
}

refusal_fake_design <- fake_blrm_design(c("tx", "time", "yprev"), tau_var = "time")

refusal_fake_profiles <- function() {
  data.frame(
    id = c("a", "b"),
    time = 1,
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:4)
  )
}

# A partial proportional-odds vglm fit with a non-proportional effect of a
# continuous `score`. At score = -20 the fitted cut-points cross (state 1 gets
# a predicted probability of about -0.22); at score = 0 they don't.
refusal_crossing_fit <- function() {
  refusal_refit("partial odds on a score", function() {
    data <- grid_data()
    data$score <- withr::with_seed(1, stats::rnorm(nrow(data)))
    mostm(
      ordered(y) ~ time + tx + age + score + yprev,
      data = data,
      id_var = "id",
      first_followup_time = 1,
      engine = "vglm",
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ score)
    )
  })
}

# Averages whose `variables` set the starting state. The values are given
# as levels of the factor previous state.
refusal_set_start <- function(fun, cell = refusal_orm) {
  start <- factor(c("1", "2"), levels = levels(grid_data()$yprev))
  refusal_call(fun, cell, variables = list(yprev = start))
}

# ---- Fitting ----

fit_refusals <- list(
  refusal(
    "F01", "blrm with case weights",
    function() refusal_mostm(engine = "blrm", weights = w),
    "`weights` is not available for `engine = \"blrm\"`"
  ),
  refusal(
    "F02", "blrm with a robust-covariance setting",
    function() refusal_mostm(engine = "blrm", type = "HC1"),
    "`type` and `cadjust` set the patient-clustered standard errors"
  ),
  refusal(
    "F03", "`data` that is not a data frame",
    function() refusal_mostm(data = as.list(grid_data())),
    "`data` must be supplied as a data frame."
  ),
  refusal(
    "F04", "an offset term in the formula",
    function() refusal_mostm(formula = y ~ time + tx + yprev + offset(age)),
    "Model offsets are not supported"
  ),
  refusal(
    "F05", "an offset argument",
    function() refusal_mostm(offset = grid_data()$age),
    "Model offsets are not supported"
  ),
  refusal(
    "F06", "orm with factor() written inside the formula",
    function() refusal_mostm(formula = y ~ time + factor(tx) + yprev),
    "cannot fit a variable that is converted to a factor inside the formula"
  ),
  refusal(
    "F08", "orm with a random-effect term",
    function() refusal_mostm(formula = y ~ time + tx + yprev + cluster(id)),
    "does not support random effects"
  ),
  refusal(
    "F09", "vglm with a random-effect term",
    function() {
      refusal_mostm(
        formula = ordered(y) ~ time + tx + yprev + mix_re(id),
        engine = "vglm"
      )
    },
    "does not support random effects"
  ),
  refusal(
    "F11", "an argument of another engine",
    function() refusal_mostm(coefstart = 1),
    "Not an argument of `engine = \"orm\"`"
  ),
  refusal(
    "F12", "asking to keep the design matrix, which mostm() always keeps",
    function() refusal_mostm(x = TRUE),
    "always keeps the design matrix and response"
  ),
  refusal(
    "F13", "orm with a VGAM family",
    function() refusal_mostm(family = VGAM::cumulative(reverse = TRUE)),
    "`family` names the link function"
  ),
  refusal(
    "F14", "a former starting-profile argument",
    function() refusal_mostm(start_time = 1),
    "does not accept the former starting-profile argument(s)"
  ),
  refusal(
    "F15", "numeric time below 1",
    function() {
      data <- grid_data()
      data$time <- data$time - 1
      mostm(refusal_formula_orm, data = data, id_var = "id")
    },
    "contains a value below 1"
  ),
  refusal(
    "F16", "factor time without `first_followup_time`",
    function() mostm(refusal_formula_orm, data = grid_cell_data(grid_cell("orm", time = "factor")), id_var = "id"),
    "Factor and character time require an explicit"
  ),
  refusal(
    "F17", "a negative case weight",
    function() {
      data <- grid_data()
      data$w[1] <- -1
      mostm(refusal_formula_orm, data = data, id_var = "id", weights = w)
    },
    "negative weights not allowed"
  ),
  refusal(
    "F18", "a case weight that is not finite",
    function() {
      data <- grid_data()
      data$w <- as.numeric(data$w)
      data$w[1] <- Inf
      mostm(refusal_formula_orm, data = data, id_var = "id", weights = w)
    },
    "`weights` must be finite."
  ),
  refusal(
    "F19", "case weights of the wrong length",
    function() refusal_mostm(weights = c(1, 2)),
    "`weights` must be numeric with one value per row of `data`."
  ),
  refusal(
    "F20", "an orm fit that rms could not fit",
    function() refusal_mostm(maxit = 0),
    "`rms::orm()` could not fit the model"
  ),
  refusal(
    "F21", "vglm with both `form2` and `subset`",
    function() {
      mostm(
        refusal_formula_vglm,
        data = grid_data(),
        id_var = "id",
        engine = "vglm",
        form2 = ~ time + tx + age + yprev,
        subset = time > 0
      )
    },
    "argument 'subset' cannot be used when argument 'form2' is used"
  ),
  refusal(
    "F22", "no `id_var`",
    function() mostm(refusal_formula_orm, data = grid_data()),
    "`mostm()` needs `id_var`"
  ),
  refusal(
    "F23", "an `id_var` that is not a column of `data`",
    function() mostm(refusal_formula_orm, data = grid_data(), id_var = "patient"),
    "ID variable 'patient' not found in data."
  ),
  refusal(
    "F24", "a `time_var` that is not a column of `data`",
    function() mostm(refusal_formula_orm, data = grid_data(), id_var = "id", time_var = "visit"),
    "The time column `visit` is not in `data`."
  ),
  refusal(
    "F25", "an unnamed argument for the engine",
    function() mostm(refusal_formula_orm, grid_data(), TRUE, id_var = "id"),
    "Arguments passed to the engine through `...` must be named."
  ),
  refusal(
    "F26", "a `first_followup_time` that is not a visit of the data",
    function() {
      mostm(refusal_formula_orm, data = grid_data(), id_var = "id", first_followup_time = 99)
    },
    "`first_followup_time = 99` is not observed in the numeric time variable."
  ),
  refusal(
    "F27", "a numeric `first_followup_time` after the earliest visit",
    function() {
      mostm(refusal_formula_orm, data = grid_data(), id_var = "id", first_followup_time = 2)
    },
    "must equal the earliest observed scheduled time"
  ),
  refusal(
    "F28", "a factor `first_followup_time` that is not a level",
    function() {
      mostm(
        refusal_formula_orm,
        data = grid_cell_data(grid_cell("orm", time = "factor")),
        id_var = "id",
        first_followup_time = "9"
      )
    },
    "is not observed in the factor or character time variable."
  ),
  refusal(
    "F29", "a logical time variable",
    function() {
      data <- grid_data()
      data$time <- data$time > 3
      mostm(y ~ tx + yprev, data = data, id_var = "id")
    },
    "The time variable must be numeric, factor, or character."
  ),
  refusal(
    "F30", "orm with a link other than the logit",
    function() refusal_mostm(family = "probit"),
    "Only cumulative-logit models are supported"
  ),
  refusal(
    "F31", "a `p2_var` column that is not in the data",
    function() refusal_mostm(p2_var = "older"),
    "given as `p2_var` is not in `data`"
  ),
  refusal(
    "F32", "a `p2_var` that names the previous-state column",
    function() refusal_mostm(p2_var = "yprev"),
    "differs from the previous-state column"
  ),
  refusal(
    "F33", "a `p_var` column that is not in the data",
    function() refusal_mostm(p_var = "previous"),
    "is not in `data`. Set `p_var`"
  )
)

test_that("mostm() refuses what it cannot fit", {
  skip_on_cran()
  check_refusals(fit_refusals)
})

# Checking a blrm formula loads rmsb and rstan, which takes seconds.
fit_refusals_blrm <- list(
  refusal(
    "F07", "blrm with factor() written inside the formula",
    function() refusal_mostm(formula = y ~ time + factor(tx) + yprev, engine = "blrm"),
    "cannot fit a variable that is converted to a factor inside the formula"
  ),
  refusal(
    "F10", "blrm with a package-prefixed cluster() term",
    function() {
      refusal_mostm(formula = y ~ time + tx + yprev + rms::cluster(id), engine = "blrm")
    },
    "`rmsb::blrm()` would fit"
  )
)

test_that("mostm() refuses blrm formulas it cannot fit", {
  skip_on_cran()
  skip_if_not_installed("rmsb")
  check_refusals(fit_refusals_blrm)
})

# ---- Prediction ----

prediction_refusals <- list(
  refusal(
    "P01", "a fit that mostm() did not make",
    function() {
      fit <- rms::orm(refusal_formula_orm, data = grid_data(), x = TRUE, y = TRUE)
      sops(fit, times = 1:4, absorb = 4)
    },
    "Model-based Markov workflows require fits created by `mostm()`."
  ),
  refusal(
    "P02", "vglm with a link other than the logit",
    function() {
      fit <- refusal_mostm(
        formula = refusal_formula_vglm,
        engine = "vglm",
        family = VGAM::cumulative(link = "probitlink", reverse = TRUE, parallel = TRUE)
      )
      sops(fit, times = 1:4, absorb = 4)
    },
    "Only cumulative-logit models are supported"
  ),
  refusal(
    "P03", "vglm with `reverse = FALSE`",
    function() {
      fit <- refusal_mostm(
        formula = refusal_formula_vglm,
        engine = "vglm",
        family = VGAM::cumulative(reverse = FALSE, parallel = TRUE)
      )
      sops(fit, times = 1:4, absorb = 4)
    },
    "the cumulative family must use reverse = TRUE"
  ),
  refusal(
    "P04", "vglm with a family other than the cumulative",
    function() {
      fit <- refusal_mostm(
        formula = refusal_formula_vglm,
        engine = "vglm",
        family = VGAM::acat(reverse = TRUE, parallel = TRUE)
      )
      sops(fit, times = 1:4, absorb = 4)
    },
    "only the cumulative family is supported"
  ),
  refusal(
    "P05", "an orm fit that holds random effects (saved before mostm() refused them)",
    function() {
      fit <- grid_fit(refusal_orm)
      fit$sigma <- 1
      sops(fit, times = 1:4, absorb = 4)
    },
    "Orm fits with random effects are not supported."
  ),
  refusal(
    "P06", "no `newdata` and no stored starting profiles",
    function() {
      fit <- grid_fit(refusal_orm)
      attr(fit, "markov_starting_profile_data") <- NULL
      sops(fit, times = 1:4, absorb = 4)
    },
    "Provide newdata: the model has no stored starting row for each patient."
  ),
  refusal(
    "P07", "a fitted patient without a row at the starting visit",
    function() {
      data <- grid_data()
      fit <- refusal_mostm(data = data[!(data$id == data$id[1] & data$time == 1), ])
      sops(fit, times = 1:4, absorb = 4)
    },
    "lack a row at the designated starting visit"
  ),
  refusal(
    "P08", "a visit that is not a level of a factor time",
    function() {
      cell <- grid_cell("orm", time = "factor")
      refusal_call(sops, cell, args = list(times = c("1", "9")))
    },
    "Requested visit values are not among fitted factor levels"
  ),
  refusal(
    "P09", "no `absorb` although a state was never a previous state",
    function() refusal_call(sops, args = list(absorb = NULL)),
    "Model prediction requires transitions from state level(s) not represented"
  ),
  refusal(
    "P10", "a model-based plot given `y_levels` with the wrong number of states",
    function() refusal_call(plot_transitions, y_levels = 1:3),
    "`y_levels` must be the fitted model's states in fitted order"
  ),
  refusal(
    "P11", "`time_covariates` without one row per visit",
    function() refusal_call(sops, time_covariates = data.frame(time = 1:2)),
    "`time_covariates` must have one row per requested time point."
  ),
  refusal(
    "P12", "factor time with a gap variable but no numeric gaps",
    function() {
      fit <- refusal_refit("factor time with a gap", function() {
        data <- grid_cell_data(grid_cell("orm", time = "factor"))
        # The grid's gap follows the visit, which factor time already fits.
        data$gap <- withr::with_seed(3, sample(1:3, nrow(data), TRUE))
        mostm(
          y ~ time + tx + age + yprev + gap,
          data = data,
          id_var = "id",
          first_followup_time = "1"
        )
      })
      sops(fit, times = as.character(1:4), absorb = 4, gap_var = "gap")
    },
    "Factor visit time with `gap_var` requires numeric gap_var values"
  ),
  refusal(
    "P13", "a blrm partial proportional-odds fit without constrained cut-points",
    function() {
      local_mocked_bindings(blrm_design_matrix = refusal_fake_design)
      sops(
        refusal_fake_blrm(pppo = 1L),
        newdata = refusal_fake_profiles(),
        times = 1:2,
        absorb = 4,
        n_draws = NULL
      )
    },
    "Only constrained partial proportional odds `blrm` models are supported."
  ),
  refusal(
    "P14", "random effects of a blrm fit without `cluster()`",
    function() {
      local_mocked_bindings(blrm_design_matrix = refusal_fake_design)
      sops(
        refusal_fake_blrm(cluster = FALSE),
        newdata = refusal_fake_profiles(),
        times = 1:2,
        absorb = 4,
        n_draws = NULL,
        include_re = TRUE
      )
    },
    "`include_re = TRUE` requires a `blrm` model fitted with `cluster()`"
  ),
  refusal(
    "P15", "a comparison of two variables at once",
    function() refusal_call(avg_comparisons, variables = list(tx = 0:1, marker = 0:1)),
    "supports exactly one named variable."
  ),
  refusal(
    "P16", "time benefit with `state_sets`",
    function() refusal_contrast(estimand = "time_benefit", state_sets = "1"),
    "`state_sets` is not used with `estimand = \"time_benefit\"`."
  ),
  refusal(
    "P17", "time benefit as a ratio",
    function() refusal_contrast(estimand = "time_benefit", comparison = "ratio"),
    "only supports `comparison = \"difference\"`."
  ),
  refusal(
    "P18", "a visit given twice to a comparison",
    function() refusal_contrast(args = list(times = c(1, 2, 2))),
    "`times` must not contain duplicate visits"
  ),
  refusal(
    "P19", "a visit given twice to a time total",
    function() refusal_call(avg_state_time, args = list(times = c(1, 2, 2))),
    "`times` must not contain duplicate visits"
  ),
  refusal(
    "P20", "`baseline_time` without `time_map`",
    function() refusal_call(avg_state_time, baseline_time = 0),
    "`baseline_time` requires `time_map`"
  ),
  refusal(
    "P21", "`baseline_time` not before the first mapped visit",
    function() {
      interpolate_sops(
        refusal_average(),
        time_map = grid_time_map(refusal_orm),
        baseline_time = 5
      )
    },
    "`baseline_time` must be earlier than the earliest mapped SOP time."
  ),
  refusal(
    "P22", "a `by` column with a name the time total uses for its output",
    function() {
      profiles <- refusal_profiles()
      profiles$state_set <- profiles$marker
      refusal_call(avg_state_time, newdata = profiles, by = "state_set")
    },
    "Grouping or counterfactual variables use reserved output names"
  ),
  refusal(
    "P23", "a profile column with a name comparisons() uses for its output",
    function() {
      profiles <- refusal_profiles()
      profiles$contrast <- 1
      refusal_call(comparisons, variables = grid_treatments, newdata = profiles)
    },
    "uses the column name(s)"
  ),
  refusal(
    "P24", "`by` with missing values among the averaged profiles",
    function() {
      profiles <- refusal_profiles()
      profiles$marker[2] <- NA
      refusal_call(avg_sops, newdata = profiles, by = "marker")
    },
    "which is missing for 1 of the 10 profiles being averaged"
  ),
  refusal(
    "P25", "`by` naming a variable that `variables` sets",
    function() refusal_average(by = "tx"),
    "which `variables` also sets"
  ),
  refusal(
    "P26", "a profile that starts in an absorbing state",
    function() {
      profiles <- refusal_profiles()
      levels(profiles$yprev) <- union(levels(profiles$yprev), "4")
      profiles$yprev[1] <- "4"
      refusal_call(sops, newdata = profiles)
    },
    "probabilities must start from a non-absorbing state"
  ),
  refusal(
    "P27", "an absorbing start given by `variables`",
    function() {
      start <- factor(c("1", "4"), levels = c(levels(grid_data()$yprev), "4"))
      refusal_call(avg_sops, variables = list(yprev = start))
    },
    "absorbing state 4"
  ),
  refusal(
    "P28", "`include_re = TRUE` on an orm fit, in sops()",
    function() refusal_call(sops, include_re = TRUE),
    "only models fitted with `rmsb::blrm()`"
  ),
  refusal(
    "P29", "`include_re = TRUE` on an orm fit, in avg_sops()",
    function() refusal_call(avg_sops, include_re = TRUE),
    "only models fitted with `rmsb::blrm()`"
  ),
  refusal(
    "P30", "`include_re = TRUE` on a vglm fit, in state_time()",
    function() refusal_call(state_time, refusal_vglm, include_re = TRUE),
    "only models fitted with `rmsb::blrm()`"
  ),
  refusal(
    "P31", "`include_re = TRUE` on a vglm fit, in avg_state_time()",
    function() refusal_call(avg_state_time, refusal_vglm, include_re = TRUE),
    "only models fitted with `rmsb::blrm()`"
  ),
  refusal(
    "P32", "`include_re = TRUE` on an orm fit, in comparisons()",
    function() refusal_call(comparisons, variables = grid_treatments, include_re = TRUE),
    "only models fitted with `rmsb::blrm()`"
  ),
  refusal(
    "P33", "`include_re = TRUE` on a vglm fit, in avg_comparisons()",
    function() refusal_contrast(refusal_vglm, include_re = TRUE),
    "only models fitted with `rmsb::blrm()`"
  ),
  refusal(
    "P34", "no `times`",
    function() sops(grid_fit(refusal_orm), absorb = 4),
    "`times` must be supplied to `sops()`."
  ),
  refusal(
    "P35", "`variables` naming a column that is not in the data",
    function() refusal_call(avg_sops, variables = list(dose = 1:2)),
    "Variables not in data: dose"
  ),
  refusal(
    "P36", "`by` naming a column that is not in the data",
    function() refusal_average(by = "region"),
    "not found in"
  ),
  refusal(
    "P37", "`by` that is neither logical nor column names",
    function() refusal_average(by = 1),
    "`by` must be `FALSE`, `TRUE`, or a character vector of column names."
  ),
  refusal(
    "P38", "a comparison with one compared value",
    function() refusal_call(avg_comparisons, variables = list(tx = 1)),
    "must contain at least two values."
  ),
  refusal(
    "P39", "a comparison with a compared value given twice",
    function() refusal_call(avg_comparisons, variables = list(tx = c(1, 1))),
    "must not contain duplicate values."
  ),
  refusal(
    "P40", "comparison `variables` that are not a named list",
    function() refusal_call(avg_comparisons, variables = list(0:1)),
    "`variables` must be a named list"
  ),
  refusal(
    "P41", "a `time_map` without an entry for a visit",
    function() refusal_call(avg_state_time, time_map = c("1" = 2, "2" = 5)),
    "`time_map` is missing entries for visit value(s)"
  ),
  refusal(
    "P42", "a `time_map` with a visit given twice",
    function() {
      interpolate_sops(refusal_average(), time_map = c("1" = 2, "1" = 3, "2" = 5, "3" = 9, "4" = 14))
    },
    "`time_map` contains duplicate visit entries"
  ),
  refusal(
    "P43", "`target_times` outside the mapped follow-up",
    function() {
      refusal_call(
        avg_state_time,
        time_map = grid_time_map(refusal_orm),
        target_times = c(1, 30)
      )
    },
    "`target_times` must stay within the supported time range"
  ),
  refusal(
    "P44", "`target_times` without `time_map`",
    function() refusal_call(avg_state_time, target_times = c(1, 3)),
    "`target_times` requires `time_map`"
  ),
  refusal(
    "P45", "`newdata` that is not a data frame",
    function() refusal_call(sops, newdata = as.list(refusal_profiles())),
    "`newdata` must be a data frame."
  ),
  refusal(
    "P46", "a state set naming a state the model does not have",
    function() refusal_contrast(state_sets = list(bad = "9")),
    "not found"
  ),
  refusal(
    "P47", "state sets without unique names",
    function() refusal_call(avg_state_time, state_sets = list(a = "1", a = "2")),
    "`state_sets` must contain nonempty state sets with unique, non-missing names."
  ),
  refusal(
    "P48", "the former `origin_time` argument",
    function() refusal_contrast(origin_time = 0),
    "is no longer supported; use `baseline_time` instead."
  ),
  refusal(
    "P49", "a model-based plot given a `time_var` other than the one given to `mostm()`",
    function() refusal_call(plot_transitions, time_var = "gap"),
    "`time_var = \"gap\"` does not match the fitted model"
  ),
  refusal(
    "P53", "a model-based plot given `p2_var` for a model fitted as first order",
    function() refusal_call(plot_transitions, p2_var = "ypprev"),
    "The model was fitted as a first-order model"
  ),
  refusal(
    "P56", "a column name that the fit records, given to a prediction function",
    function() refusal_call(sops, p2_var = "ypprev"),
    "`p2_var` is not an argument of this function"
  ),
  refusal(
    "P54", "a linear-predictor plot faceted by an absorbing previous state",
    function() {
      plot_lp_difference(
        grid_fit(refusal_orm),
        profile = grid_profiles(refusal_orm)[1, , drop = FALSE],
        variables = list(tx = c(0, 1)),
        times = 1:4,
        previous_states = 1:4
      )
    },
    "`previous_states` includes state(s) that never occur in `yprev`"
  ),
  refusal(
    "P55", "a linear-predictor plot of a profile missing factor levels of yprev",
    function() {
      profile <- grid_profiles(refusal_vglm)[1, , drop = FALSE]
      profile$yprev <- factor(as.character(profile$yprev))
      plot_lp_difference(
        grid_fit(refusal_vglm),
        profile = profile,
        variables = list(tx = c(0, 1)),
        times = 1:4
      )
    },
    "The model returned linear predictors for 8 of the 24 plotted rows"
  )
)

test_that("SOP and summary functions refuse what they cannot compute", {
  skip_on_cran()
  check_refusals(prediction_refusals)
})

prediction_refusals_with_extra_fits <- list(
  refusal(
    "P50", "a second-order profile whose state two visits back is absorbing",
    function() {
      cell <- grid_cell("orm", order = 2L)
      profiles <- grid_profiles(cell)
      levels(profiles$ypprev) <- union(levels(profiles$ypprev), "4")
      profiles$ypprev[1] <- "4"
      refusal_call(sops, cell, newdata = profiles)
    },
    "`ypprev` is the state two visits back"
  ),
  refusal(
    "P51", "partial proportional-odds cut-points that cross",
    function() {
      profiles <- grid_profiles(refusal_vglm)
      profiles$score <- c(0, -20, 1)
      sops(refusal_crossing_fit(), newdata = profiles, times = 1:4, absorb = 4)
    },
    "cut-points (thresholds) of the partial proportional-odds model cross",
    class = "mostr_crossed_thresholds"
  ),
  refusal(
    "P52", "partial proportional-odds cut-points that cross, in a comparison",
    function() {
      avg_comparisons(
        refusal_crossing_fit(),
        variables = list(score = c(0, -20)),
        times = 1:4,
        absorb = 4
      )
    },
    class = "mostr_crossed_thresholds"
  )
)

test_that("SOP functions refuse second-order and partial proportional-odds cases", {
  skip_on_cran()
  check_refusals(prediction_refusals_with_extra_fits)
})

# ---- Inference ----

inference_refusals <- list(
  refusal(
    "I01", "analytical inference on a blrm result",
    function() {
      local_mocked_bindings(blrm_design_matrix = refusal_fake_design)
      average <- avg_sops(
        refusal_fake_blrm(),
        newdata = refusal_fake_profiles(),
        times = 1:2,
        absorb = 4,
        n_draws = NULL
      )
      inferences(average, method = "delta")
    },
    "Analytical delta inference is not available for `blrm` models"
  ),
  refusal(
    "I02", "unconditional variance for results with one row per profile",
    function() inferences(refusal_call(sops), method = "delta", vcov = "unconditional"),
    "supports only `vcov = \"conditional\"` or a coefficient covariance matrix."
  ),
  refusal(
    "I03", "unconditional variance for an average over `newdata`",
    function() {
      inferences(
        refusal_average(newdata = refusal_profiles()),
        method = "delta",
        vcov = "unconditional"
      )
    },
    "`vcov = \"unconditional\"` requires the stored fitted-patient cohort."
  ),
  refusal(
    "I04", "unconditional variance when `variables` sets the starting state",
    function() {
      inferences(
        refusal_set_start(avg_sops),
        method = "delta",
        vcov = "unconditional"
      )
    },
    "Setting the starting state with `variables` supports only conditional"
  ),
  refusal(
    "I05", "the ordinary bootstrap with a set starting state",
    function() {
      inferences(
        refusal_set_start(avg_sops),
        method = "bootstrap",
        n_draws = 2
      )
    },
    "Setting the starting state with `variables` supports only conditional"
  ),
  refusal(
    "I06", "the fractional weighted bootstrap with a set starting state",
    function() {
      inferences(
        refusal_set_start(avg_comparisons),
        method = "fwb",
        n_draws = 2
      )
    },
    "Setting the starting state with `variables` supports only conditional"
  ),
  refusal(
    "I07", "the score bootstrap (internal) with a set starting state",
    function() {
      mostr:::inferences_score_bootstrap(
        refusal_set_start(avg_state_time),
        n_draws = 2
      )
    },
    "Setting the starting state with `variables` supports only conditional"
  ),
  refusal(
    "I08", "delta inference with percentile intervals for average comparisons",
    function() inferences(refusal_contrast(), method = "delta", conf_type = "perc"),
    "Analytical average comparisons require `conf_type = \"wald\"`."
  ),
  refusal(
    "I09", "delta inference with percentile intervals for average time",
    function() inferences(refusal_call(avg_state_time), method = "delta", conf_type = "perc"),
    "Analytical average time requires `conf_type = \"wald\"`."
  ),
  refusal(
    "I10", "delta inference with percentile intervals for one row per profile",
    function() {
      inferences(
        refusal_call(comparisons, variables = grid_treatments),
        method = "delta",
        conf_type = "perc"
      )
    },
    "Analytical intervals for one row per profile require `conf_type = \"wald\"`."
  ),
  refusal(
    "I11", "delta inference with percentile intervals for SOPs",
    function() inferences(refusal_average(), method = "delta", conf_type = "perc"),
    "Delta SOP inference supports `conf_type = \"logit\"` or `\"wald\"`."
  ),
  refusal(
    "I12", "a character `vcov` with a draw method",
    function() inferences(refusal_average(), method = "mvn", vcov = "conditional"),
    "Character `vcov` choices are only available with `method = \"delta\"`."
  ),
  refusal(
    "I13", "logit intervals with a draw method",
    function() inferences(refusal_average(), method = "mvn", conf_type = "logit"),
    "`conf_type = \"logit\"` is only available with `method = \"delta\"`."
  ),
  refusal(
    "I14", "a zero null value for a ratio",
    function() inferences(refusal_contrast(comparison = "ratio"), method = "delta", null = 0),
    "Ratio comparisons cannot be tested against `null = 0`."
  ),
  refusal(
    "I15", "the score bootstrap through inferences()",
    function() inferences(refusal_average(), method = "score_bootstrap"),
    "`method = \"score_bootstrap\"` is no longer offered by `inferences()`."
  ),
  refusal(
    "I16", "an object that is not a mostr result",
    function() inferences(data.frame(estimate = 1)),
    "inferences() requires a 'markov_avg_sops'"
  ),
  refusal(
    "I17", "a confidence level outside (0, 1)",
    function() inferences(refusal_average(), method = "delta", conf_level = 95),
    "must be a single finite number between 0 and 1."
  ),
  refusal(
    "I18", "a null value given as text",
    function() inferences(refusal_contrast(), method = "delta", null = "zero"),
    "`null` must be \"auto\", NULL, or a single finite numeric value."
  ),
  refusal(
    "I19", "an unknown variance type",
    function() inferences(refusal_average(), method = "delta", vcov = "robust"),
    "`vcov` must be \"conditional\", \"unconditional\", NULL, or a coefficient covariance matrix."
  ),
  refusal(
    "I20", "a covariance matrix of draw-based results",
    function() stats::vcov(inferences(refusal_average(), n_draws = 5, seed = 1)),
    "No analytical inference state found."
  ),
  refusal(
    "I21", "draws of a result without inference",
    function() get_draws(refusal_average()),
    "No draws found. Run inferences() with return_draws = TRUE."
  ),
  refusal(
    "I22", "a seed that is not a number",
    function() inferences(refusal_average(), n_draws = 5, seed = "a"),
    "`seed` must be a single finite number or NULL."
  ),
  refusal(
    "I23", "a cluster formula with two columns",
    function() {
      inferences(refusal_average(), method = "delta", vcov = "conditional", cluster = ~ id + tx)
    },
    "A `cluster` formula must select exactly one patient ID column"
  ),
  refusal(
    "I24", "a cluster formula naming a column the fit did not store",
    function() {
      inferences(refusal_average(), method = "delta", vcov = "conditional", cluster = ~ region)
    },
    "requires that column in the model's stored fitting data."
  ),
  refusal(
    "I25", "a cluster vector of the wrong length",
    function() {
      inferences(refusal_average(), method = "delta", vcov = "conditional", cluster = 1:3)
    },
    "Length of `cluster` (3) does not match"
  ),
  refusal(
    "I26", "a covariance matrix whose names are not the model's coefficients",
    function() {
      V <- diag(8)
      dimnames(V) <- list(letters[1:8], letters[1:8])
      inferences(refusal_average(), method = "delta", vcov = V)
    },
    "must be uniquely named and match the complete raw coefficient vector."
  ),
  refusal(
    "I27", "delta inference with `by`, for SOPs",
    function() inferences(refusal_average(by = "marker"), method = "delta", vcov = "conditional"),
    "Analytical delta inference with `by` is not yet supported."
  ),
  refusal(
    "I28", "delta inference with `by`, for comparisons",
    function() inferences(refusal_contrast(by = "marker"), method = "delta", vcov = "conditional"),
    "Analytical delta inference with `by` is not yet supported."
  ),
  refusal(
    "I29", "delta inference for an average ratio",
    function() inferences(refusal_contrast(comparison = "ratio"), method = "delta"),
    "currently support differences only"
  ),
  refusal(
    "I30", "delta inference for ratios of each profile",
    function() {
      inferences(
        refusal_call(comparisons, variables = grid_treatments, comparison = "ratio"),
        method = "delta"
      )
    },
    "currently support differences only"
  ),
  refusal(
    "I31", "delta inference for average time benefit",
    function() inferences(refusal_contrast(estimand = "time_benefit"), method = "delta"),
    "time benefit is not yet supported."
  ),
  refusal(
    "I32", "delta inference for the time benefit of each profile",
    function() {
      inferences(
        refusal_call(comparisons, variables = grid_treatments, estimand = "time_benefit"),
        method = "delta"
      )
    },
    "time benefit is not yet supported."
  ),
  refusal(
    "I33", "draws above the memory limit for one row per profile",
    function() {
      withr::local_options(mostr.max_sops_draw_cells = 10)
      inferences(refusal_call(comparisons, variables = grid_treatments), n_draws = 5, seed = 1)
    },
    "Inference for one row per profile would keep"
  ),
  refusal(
    "I34", "the score bootstrap (internal) of an orm fit without clusters",
    function() mostr:::inferences_score_bootstrap(refusal_average(), n_draws = 2),
    "`cluster` is required for orm score-bootstrap inference."
  ),
  refusal(
    "I35", "the score bootstrap (internal) with a covariance matrix",
    function() {
      fit <- grid_fit(refusal_vglm)
      mostr:::inferences_score_bootstrap(
        refusal_average(refusal_vglm),
        n_draws = 2,
        vcov = stats::vcov(fit)
      )
    },
    "`vcov` cannot be supplied when `method = \"score_bootstrap\"`"
  ),
  refusal(
    "I49", "starting refits from the original coefficients on an orm fit",
    function() {
      inferences(refusal_average(), method = "bootstrap", n_draws = 2, use_coefstart = TRUE)
    },
    "`use_coefstart` can be `TRUE` only for models fitted with `engine = \"vglm\"`"
  ),
  refusal(
    "I50", "starting values with a method that doesn't refit",
    function() {
      inferences(refusal_average(refusal_vglm), method = "mvn", n_draws = 2, use_coefstart = TRUE)
    },
    "`use_coefstart` applies only to `method = \"bootstrap\"` and `method = \"fwb\"`"
  ),
  refusal(
    "I51", "`use_coefstart` that is not TRUE or FALSE",
    function() {
      inferences(refusal_average(), method = "bootstrap", n_draws = 2, use_coefstart = NA)
    },
    "`use_coefstart` must be `TRUE` or `FALSE`."
  ),
  refusal(
    "I52", "an ordinary bootstrap resample that lacks the state a supplied profile starts in",
    function() {
      # The one resample holds only the patients never in state 2; the
      # refusal comes before any refit.
      data <- grid_cell_data(refusal_orm)
      ids <- setdiff(unique(data$id), data$id[data$y == 2 | data$yprev == 2])
      local_mocked_bindings(
        fast_group_bootstrap = function(data, id_var, n_boot) {
          list(data.frame(original_id = ids, new_id = paste0(ids, "_1"), boot_id = 1))
        }
      )
      profiles <- refusal_profiles()
      profiles$yprev[1] <- "2"
      average <- refusal_call(avg_sops, newdata = profiles)
      inferences(average, method = "bootstrap", n_draws = 1)
    },
    "The ordinary bootstrap left out every patient who was ever in state 2"
  )
)
test_that("inferences() refuses what it cannot compute", {
  skip_on_cran()
  check_refusals(inference_refusals)
})

inference_refusals_with_extra_fits <- list(
  refusal(
    "I36", "unconditional variance on a weighted orm fit",
    function() {
      inferences(
        refusal_average(grid_cell("orm", weighted = TRUE)),
        method = "delta",
        vcov = "unconditional"
      )
    },
    "Models fitted with case weights (`weights` in `mostm()`) do not support unconditional variance"
  ),
  refusal(
    "I37", "unconditional variance on a weighted vglm fit",
    function() {
      inferences(
        refusal_contrast(grid_cell("vglm", weighted = TRUE)),
        method = "delta",
        vcov = "unconditional"
      )
    },
    "Models fitted with case weights (`weights` in `mostm()`) do not support unconditional variance"
  ),
  refusal(
    "I38", "the default variance of an average on a weighted fit, which is unconditional",
    function() {
      inferences(
        refusal_call(avg_state_time, grid_cell("vglm", weighted = TRUE)),
        method = "delta"
      )
    },
    "Models fitted with case weights (`weights` in `mostm()`) do not support unconditional variance"
  ),
  refusal(
    "I39", "unconditional variance on a penalized orm fit",
    function() {
      fit <- refusal_refit("penalized orm", function() refusal_mostm(penalty = 1))
      average <- do.call(
        avg_sops,
        c(list(fit, variables = grid_treatments), grid_public_args(refusal_orm))
      )
      inferences(average, method = "delta", vcov = "unconditional")
    },
    "does not currently support penalized orm likelihoods"
  ),
  refusal(
    "I40", "delta inference on a second-order fit",
    function() {
      inferences(refusal_average(grid_cell("orm", order = 2L)), method = "delta", vcov = "conditional")
    },
    "for second-order Markov models is not yet supported."
  ),
  refusal(
    "I41", "delta inference on a partial proportional-odds fit",
    function() {
      inferences(
        refusal_average(grid_cell("vglm", odds = "partial")),
        method = "delta",
        vcov = "conditional"
      )
    },
    "currently requires a full proportional-odds model"
  ),
  refusal(
    "I42", "the ordinary bootstrap of SOPs for each profile",
    function() inferences(refusal_call(sops), method = "bootstrap", n_draws = 2),
    "Standard refit bootstrap is not supported for individual-level"
  ),
  refusal(
    "I43", "the ordinary bootstrap of comparisons for each profile",
    function() {
      inferences(
        refusal_call(comparisons, variables = grid_treatments),
        method = "bootstrap",
        n_draws = 2
      )
    },
    "Standard refit bootstrap is not supported for results with one row per profile"
  ),
  refusal(
    "I44", "the ordinary bootstrap of time totals for each profile",
    function() inferences(refusal_call(state_time, refusal_vglm), method = "bootstrap", n_draws = 2),
    "Standard refit bootstrap is not supported for results with one row per profile"
  ),
  refusal(
    "I45", "an ordinary bootstrap resample that leaves out a whole `by` group",
    function() {
      fit <- refusal_refit("one-patient group", function() {
        data <- grid_data()
        data$rare <- as.integer(data$id == data$id[1])
        refusal_mostm(data = data)
      })
      average <- do.call(
        avg_sops,
        c(list(fit, by = "rare"), grid_public_args(refusal_orm))
      )
      inferences(average, method = "bootstrap", n_draws = 20, seed = 1)
    },
    "The ordinary bootstrap left out every patient of"
  ),
  refusal(
    "I46", "refit inference on data with one row per patient",
    function() {
      fit <- refusal_refit("one row per patient", function() {
        data <- grid_data()
        refusal_mostm(data = data[data$time == 1, ], formula = y ~ tx + age + yprev)
      })
      average <- avg_sops(fit, variables = grid_treatments, times = 1:2, absorb = 4)
      inferences(average, method = "bootstrap", n_draws = 2, seed = 1)
    },
    "Bootstrap inference requires full longitudinal data"
  ),
  refusal(
    "I47", "fractional bootstrap draws for profiles that hold a column named `fwb_weight`",
    function() {
      fit <- refusal_refit("fwb_weight column", function() {
        data <- grid_data()
        data$fwb_weight <- 1
        refusal_mostm(data = data)
      })
      each <- do.call(sops, c(list(fit), grid_public_args(refusal_orm)))
      inferences(each, method = "fwb", n_draws = 2, seed = 1)
    },
    "reserved for bootstrap draw weights"
  ),
  refusal(
    "I48", "MVN draws in which partial proportional-odds cut-points cross",
    function() {
      # At score = -10 the fitted cut-points don't cross, but they do for
      # some coefficient draws.
      profiles <- grid_profiles(refusal_vglm)
      profiles$score <- -10
      each <- sops(refusal_crossing_fit(), newdata = profiles, times = 1:2, absorb = 4)
      inferences(each, n_draws = 50, seed = 1)
    },
    "simulation draws, so those simulation draws give a negative probability",
    class = "mostr_crossed_thresholds"
  ),
  refusal(
    "I53", "the ordinary bootstrap when `variables` sets the state two visits back",
    function() {
      start <- factor(c("1", "2"), levels = levels(grid_data()$ypprev))
      average <- refusal_call(
        avg_sops,
        grid_cell("orm", order = 2L),
        variables = list(ypprev = start)
      )
      inferences(average, method = "bootstrap", n_draws = 2)
    },
    "Setting the starting state with `variables` supports only conditional"
  )
)

test_that("inferences() refuses weighted, penalized, second-order and bootstrap cases", {
  skip_on_cran()
  check_refusals(inference_refusals_with_extra_fits)
})

test_that("a covariance matrix selects conditional variance rather than being refused", {
  # The catalogue listed "a coefficient covariance matrix cannot be supplied
  # for unconditional inference" as a refusal. A matrix passed as `vcov`
  # selects conditional variance instead, even for averages, whose default
  # is unconditional, so that message is reachable only from internal code.
  # The public call succeeds and equals conditional variance with the model's
  # own covariance matrix.
  fit <- grid_fit(refusal_orm)
  average <- refusal_average()
  with_matrix <- inferences(
    average,
    method = "delta",
    vcov = stats::vcov(fit, intercepts = "all")
  )
  conditional <- inferences(average, method = "delta", vcov = "conditional")
  expect_equal(with_matrix$std.error, conditional$std.error, tolerance = 1e-10)
})

# ---- Simulation ----

# Patients and a linear predictor for `sim_markov(NULL, ...)`, where the user
# supplies the model.
refusal_sim_patients <- data.frame(id = 1:2, yprev = c(1, 1), tx = c(0, 1))
refusal_sim_lp <- function(yprev, t, tx, treatment_effect, extra_params) 0

# One visit of one patient under a user model whose linear predictor `lp`
# returns the given value.
refusal_sim_threshold_lp <- function(intercepts, value) {
  sim_markov(
    NULL,
    baseline_data = data.frame(id = 1, yprev = 1, tx = 1),
    follow_up_time = 1,
    intercepts = intercepts,
    lp_function = function(yprev, t, tx, treatment_effect, extra_params) value,
    absorbing_states = integer(),
    covariate_names = "tx",
    seed = 1
  )
}

simulation_refusals <- list(
  refusal(
    "S01", "no `preset`",
    function() sim_markov(n_patients = 10),
    "`preset` is required"
  ),
  refusal(
    "S02", "a preset name that doesn't exist",
    function() sim_markov("actt3", n_patients = 10),
    "`preset` must be one of \"actt1\", \"actt2\", \"violet\", or NULL"
  ),
  refusal(
    "S03", "two preset names",
    function() sim_markov(c("actt1", "actt2")),
    "`preset` must be one of"
  ),
  refusal(
    "S04", "a preset with both `baseline_data` and `n_patients`",
    function() {
      sim_markov(
        "actt2",
        n_patients = 2,
        baseline_data = data.frame(id = 1:2, yprev = c(5L, 6L), tx = c(0L, 1L))
      )
    },
    "Supply either `baseline_data` or `n_patients` and `treatment_prob`"
  ),
  refusal(
    "S05", "a preset with both `baseline_data` and `treatment_prob`",
    function() {
      sim_markov(
        "actt2",
        treatment_prob = 1,
        baseline_data = data.frame(id = 1:2, yprev = c(5L, 6L), tx = c(0L, 1L))
      )
    },
    "Supply either `baseline_data` or `n_patients` and `treatment_prob`"
  ),
  refusal(
    "S06", "a user model without `baseline_data`",
    function() sim_markov(NULL, intercepts = 0, lp_function = refusal_sim_lp),
    "With `preset = NULL`, supply `baseline_data`."
  ),
  refusal(
    "S07", "a user model without `intercepts` and `lp_function`",
    function() sim_markov(NULL, baseline_data = refusal_sim_patients),
    "With `preset = NULL`, supply `intercepts`, `lp_function`."
  ),
  refusal(
    "S08", "a user model with both `baseline_data` and `n_patients`",
    function() {
      sim_markov(
        NULL,
        n_patients = 2,
        baseline_data = refusal_sim_patients,
        intercepts = 0,
        lp_function = refusal_sim_lp
      )
    },
    "Supply either `baseline_data` or `n_patients`"
  ),
  refusal(
    "S09", "a user model with `n_patients` instead of `baseline_data`",
    function() {
      sim_markov(NULL, n_patients = 2, intercepts = 0, lp_function = refusal_sim_lp)
    },
    "`n_patients` and `treatment_prob` sample patients from a preset"
  ),
  refusal(
    "S10", "`baseline_data` that is not a data frame",
    function() {
      sim_markov(NULL, baseline_data = list(), intercepts = 0, lp_function = refusal_sim_lp)
    },
    "baseline_data must be a data frame"
  ),
  refusal(
    "S11", "`baseline_data` without a `yprev` column",
    function() {
      sim_markov(
        NULL,
        baseline_data = data.frame(id = 1),
        intercepts = 0,
        lp_function = refusal_sim_lp
      )
    },
    "baseline_data must contain columns: yprev"
  ),
  refusal(
    "S12", "`intercepts` that are not numeric",
    function() {
      sim_markov(
        NULL,
        baseline_data = refusal_sim_patients,
        intercepts = character(),
        lp_function = refusal_sim_lp
      )
    },
    "intercepts must be a numeric vector"
  ),
  refusal(
    "S13", "an `lp_function` that is not a function",
    function() {
      sim_markov(
        NULL,
        baseline_data = refusal_sim_patients,
        intercepts = 0,
        lp_function = "nope"
      )
    },
    "lp_function must be a function"
  ),
  refusal(
    "S14", "a linear predictor with neither one value nor one per threshold",
    function() refusal_sim_threshold_lp(c(-1, 0, 1), c(0, 0)),
    "lp_function must return either length 1 or length(intercepts)"
  ),
  refusal(
    "S15", "a linear predictor that is not finite",
    function() refusal_sim_threshold_lp(c(-1, 0, 1), NA_real_),
    "lp_function must return finite numeric values"
  ),
  refusal(
    "S16", "threshold-specific linear predictors that give negative state probabilities",
    function() refusal_sim_threshold_lp(c(0, 1), c(5, -5)),
    "create decreasing cumulative probabilities"
  )
)

test_that("sim_markov() refuses what it cannot simulate", {
  check_refusals(simulation_refusals)
})

# ---- Robust variance ----

# A small logistic vglm fit with eight clusters of ten rows.
refusal_robcov_data <- function() {
  withr::with_seed(90011, {
    data.frame(
      y = stats::rbinom(80, 1, 0.5),
      x = stats::rnorm(80),
      id = rep(1:8, each = 10)
    )
  })
}

refusal_robcov_fit <- function() {
  refusal_refit("robcov_vglm logistic", function() {
    VGAM::vglm(y ~ x, family = VGAM::binomialff, data = refusal_robcov_data())
  })
}

robcov_refusals <- list(
  refusal(
    "V01", "a fit that is not a vglm fit",
    function() robcov_vglm(stats::lm(y ~ x, data = refusal_robcov_data())),
    "'fit' must be a vglm object"
  ),
  refusal(
    "V02", "a `cluster` of the wrong length",
    function() robcov_vglm(refusal_robcov_fit(), cluster = 1:10),
    "Length of 'cluster' must equal"
  ),
  refusal(
    "V03", "a `cluster` with a missing value",
    function() {
      cluster <- refusal_robcov_data()$id
      cluster[1] <- NA
      robcov_vglm(refusal_robcov_fit(), cluster = cluster)
    },
    "'cluster' contains missing values"
  ),
  refusal(
    "V04", "a single cluster",
    function() robcov_vglm(refusal_robcov_fit(), cluster = rep(1, 80)),
    "Cluster-robust covariance requires at least two clusters."
  )
)

test_that("robcov_vglm() refuses what it cannot compute", {
  skip_if_not_installed("VGAM")
  check_refusals(robcov_refusals)
})
