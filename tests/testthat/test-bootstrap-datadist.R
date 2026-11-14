# mostr needs no rms `datadist`. Bootstrap refits of an orm fit must neither
# depend on a user's `datadist` setting nor change it or the global environment.

# A user's own `dd` object in the global environment, locked so that any write
# to it is an error, and a `datadist` option naming an object that does not
# exist. Both are restored when `env` ends.
local_user_datadist <- function(env = parent.frame()) {
  had_dd <- exists("dd", envir = globalenv(), inherits = FALSE)
  if (had_dd) {
    old_dd <- get("dd", envir = globalenv())
  }
  assign("dd", "user's own object", envir = globalenv())
  lockBinding("dd", globalenv())
  withr::defer(
    {
      unlockBinding("dd", globalenv())
      if (had_dd) {
        assign("dd", old_dd, envir = globalenv())
      } else {
        rm("dd", envir = globalenv())
      }
    },
    envir = env
  )
  withr::local_options(datadist = "no_such_dd", .local_envir = env)
}

global_names <- function() {
  setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed")
}

test_that("inferences() and bootstrap_model_coefs() have no update_datadist argument", {
  expect_false("update_datadist" %in% names(formals(inferences)))
  expect_false("update_datadist" %in% names(formals(bootstrap_model_coefs)))
})

test_that("bootstrap refits of an orm fit ignore and keep a user's datadist", {
  skip_on_cran()
  # A spline of time and a factor previous state: terms whose rms design
  # information a `datadist` could describe.
  cell <- grid_cell("orm", time = "rcs")
  fit <- grid_fit(cell)
  result <- do.call(
    avg_sops,
    c(list(fit, newdata = inference_profiles(cell)), grid_public_args(cell))
  )
  run <- function(method) {
    suppressWarnings(
      inferences(result, method = method, n_draws = 3, seed = 7)
    )
  }
  coefs <- function() {
    withr::with_seed(7, bootstrap_model_coefs(fit, n_boot = 3))
  }

  # Reference results with no `datadist` set.
  withr::local_options(datadist = NULL)
  reference <- list(bootstrap = run("bootstrap"), fwb = run("fwb"))
  reference_coefs <- coefs()

  local_user_datadist()
  names_before <- global_names()

  for (method in c("bootstrap", "fwb")) {
    with_user_datadist <- run(method)
    expect_equal(
      with_user_datadist$std.error,
      reference[[method]]$std.error,
      info = method
    )
    expect_true(all(is.finite(with_user_datadist$std.error)), info = method)
  }
  expect_equal(coefs(), reference_coefs)
  expect_false(anyNA(reference_coefs))

  expect_identical(getOption("datadist"), "no_such_dd")
  expect_identical(get("dd", envir = globalenv()), "user's own object")
  expect_identical(global_names(), names_before)
})
