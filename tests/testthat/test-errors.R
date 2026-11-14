# mostr raises every error through mostr_abort() (R/mostr-abort.R). These tests
# check what users get from it: the class, the literal message, and a call that
# names the function they called rather than an internal one.

test_that("mostr_abort() errors carry class mostr_error, extra classes and data", {
  condition <- tryCatch(
    mostr:::mostr_abort(
      "Part one, ",
      2L,
      " and {three} at 100% of %s.",
      class = "mostr_test_error",
      data = list(size = 3)
    ),
    error = identity
  )

  expect_s3_class(condition, c("mostr_test_error", "mostr_error", "error"))
  expect_identical(
    conditionMessage(condition),
    "Part one, 2 and {three} at 100% of %s."
  )
  expect_identical(condition$size, 3)
  # Called directly, with no exported function on the stack, no call is shown.
  expect_null(conditionCall(condition))
  expect_error(
    mostr:::mostr_abort("x", data = list(parent = 1)),
    "named fields"
  )
})

test_that("mostr_abort() re-raises a condition object and adds mostr_error", {
  original <- simpleError("Original {message}.")
  condition <- tryCatch(mostr:::mostr_abort(original), error = identity)
  expect_s3_class(condition, c("simpleError", "mostr_error", "error"))
  expect_identical(conditionMessage(condition), "Original {message}.")

  raised <- tryCatch(
    mostr:::mostr_abort("Raised.", class = "mostr_test_error"),
    error = identity
  )
  again <- tryCatch(mostr:::mostr_abort(raised), error = identity)
  expect_identical(class(again), class(raised))
  expect_identical(conditionMessage(again), "Raised.")
})

test_that("errors from public calls name the exported function the user called", {
  # This error is raised by an internal helper that validates the arguments.
  exports <- getNamespaceExports("mostr")

  condition <- tryCatch(
    sim_markov(preset = "violet", n_patients = -1),
    error = identity
  )
  expect_s3_class(condition, "mostr_error")
  expect_identical(error_call_name(condition), "sim_markov")
  expect_snapshot(sim_markov(preset = "violet", n_patients = -1), error = TRUE)

  # Written in user code, the inner call is the one reported, even though
  # `inferences()` evaluates it.
  user_env <- new.env(parent = globalenv())
  condition <- tryCatch(
    eval(quote(mostr::inferences(mostr::sops(42))), user_env),
    error = identity
  )
  expect_identical(error_call_name(condition), "sops")
  expect_true(error_call_name(condition) %in% exports)
})

test_that("R code raises errors with mostr_abort(), not stop()", {
  # An error from base stop() names the internal function that raised it, and
  # lacks class mostr_error. See "Errors" in ARCHITECTURE.md.
  namespace <- asNamespace("mostr")
  calls_stop <- vapply(
    ls(namespace, all.names = TRUE),
    function(name) {
      object <- get(name, envir = namespace)
      is.function(object) &&
        !is.primitive(object) &&
        "stop" %in% all.names(body(object))
    },
    logical(1)
  )
  expect_identical(names(calls_stop)[calls_stop], character())
})
