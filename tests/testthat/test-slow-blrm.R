# Real rmsb::blrm() fits against the readable loops, draw by draw, using only
# predict.blrm(). This is the one place where the blrm design matrix is not mocked.
# Each fit takes about a minute, so these run only with MOSTR_RUN_SLOW=true.

slow_blrm_fit <- function(formula, ...) {
  withr::local_options(mc.cores = 2)
  suppressWarnings(mostm(
    formula,
    data = grid_data(),
    id_var = "id",
    absorb = 4L,
    engine = "blrm",
    iter = 200,
    warmup = 100,
    chains = 2,
    method = "sampling",
    refresh = 0,
    sampling.args = list(seed = 2026),
    ...
  ))
}

# Draw-by-draw predictions from predict.blrm(). One call returns every draw, so
# results are kept by row and the loops then read one draw at a time.
slow_blrm_predictor <- function(fit) {
  saved <- new.env(parent = emptyenv())
  function(draw) {
    function(row) {
      key <- paste(row$time, row$tx, row$age, row$yprev, row$ypprev)
      if (is.null(saved[[key]])) {
        all_draws <- stats::predict(
          fit,
          row,
          type = "fitted.ind",
          posterior.summary = "all"
        )
        size <- dim(all_draws)
        saved[[key]] <- matrix(all_draws[, 1L, ], size[1L], size[3L])
      }
      saved[[key]][draw, ]
    }
  }
}

slow_blrm_check <- function(fit, second_order, label) {
  times <- 1:3
  states <- 1:4
  profiles <- grid_profiles(grid_cell("orm"))
  args <- list(
    fit,
    newdata = profiles,
    times = times,
    n_draws = NULL,
    return_draws = TRUE,
    posterior_summary = "mean"
  )
  result <- do.call(sops, args)
  draws <- attr(result, "draws")
  predictor <- slow_blrm_predictor(fit)

  reference <- list()
  for (key in unique(paste(draws$draw_id, draws$rowid))) {
    parts <- as.integer(strsplit(key, " ")[[1]])
    patient <- profiles[parts[2], , drop = FALSE]
    predict_transition <- predictor(parts[1])
    reference[[key]] <- if (second_order) {
      slow_second_order_sops(predict_transition, patient, times, states, 4L)
    } else {
      slow_first_order_sops(predict_transition, patient, times, states, 4L)
    }
  }

  gaps <- vapply(
    seq_len(nrow(draws)),
    function(row) {
      loop <- reference[[paste(draws$draw_id[row], draws$rowid[row])]]
      abs(
        draws$estimate[row] -
          loop[draws$time[row], as.integer(as.character(draws$state[row]))]
      )
    },
    numeric(1)
  )
  expect_lt(max(gaps), 1e-9, label = label)
}

test_that("real blrm fits match the readable loop draw by draw", {
  skip_unless_slow_tier()
  skip_if_not_installed("rmsb")

  slow_blrm_check(
    slow_blrm_fit(y ~ time + tx + age + yprev),
    second_order = FALSE,
    "first order"
  )
  slow_blrm_check(
    slow_blrm_fit(
      y ~ time + tx + age + yprev,
      ppo = ~time,
      cppo = function(y) as.numeric(y)
    ),
    second_order = FALSE,
    "first order, partial odds"
  )
  slow_blrm_check(
    slow_blrm_fit(y ~ time + tx + age + yprev + ypprev, p2_var = "ypprev"),
    second_order = TRUE,
    "second order"
  )
})
