# The bookkeeping of the compiled Stan program's cache, with the compile step
# replaced by a function that returns a stand-in object: no Stan program is
# compiled here.

ms_fake_stan_model <- function(id) {
  structure(list(id = id), class = "stanmodel")
}

test_that("the cache key changes with the program and the toolchain", {
  versions <- c(rstan = "2.32.7", StanHeaders = "2.32.10", Rcpp = "1.1.0")
  key <- markov_stan_cache_key("model {}", versions)
  expect_identical(key, markov_stan_cache_key("model {}", versions))
  expect_false(identical(key, markov_stan_cache_key("model { }", versions)))
  newer <- versions
  newer[["StanHeaders"]] <- "2.32.11"
  expect_false(identical(key, markov_stan_cache_key("model {}", newer)))
  # The installed versions that the compiled program depends on are read.
  installed <- markov_stan_toolchain_versions()
  expect_true(all(names(installed) %in%
    c("rstan", "StanHeaders", "Rcpp", "RcppEigen", "RcppParallel")))
})

test_that("a compiled program is saved once and outdated ones are removed", {
  dir <- withr::local_tempdir()
  # Files of an earlier program, and a file that is not mostr's program.
  file.create(file.path(dir, c("po_paths_old1.rds", "po_paths_old2.rds")))
  file.create(file.path(dir, "markov_ml_other.rds"))
  cache <- markov_stan_cache_dir(dir)
  expect_true(cache$use)
  compiled <- 0L
  compile <- function() {
    compiled <<- compiled + 1L
    ms_fake_stan_model(compiled)
  }
  model <- markov_cached_stan_model(cache, "new", compile)
  expect_identical(model$id, 1L)
  expect_setequal(list.files(dir), c("po_paths_new.rds", "markov_ml_other.rds"))
  # The next session loads it without compiling.
  again <- markov_cached_stan_model(cache, "new", compile)
  expect_identical(compiled, 1L)
  expect_identical(again$id, 1L)
  # A damaged file is compiled again and replaced.
  writeLines("not a model", file.path(dir, "po_paths_new.rds"))
  repaired <- markov_cached_stan_model(cache, "new", compile)
  expect_identical(repaired$id, 2L)
  expect_identical(readRDS(file.path(dir, "po_paths_new.rds"))$id, 2L)
})

test_that("an unwritable or switched-off cache keeps the program for the session", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(markov_dir_writable = function(dir) FALSE)
  cache <- markov_stan_cache_dir(dir)
  expect_false(cache$use)
  expect_match(markov_stan_compile_message(cache), "per session, because the cache")
  expect_match(markov_stan_compile_message(cache), dir, fixed = TRUE)
  model <- markov_cached_stan_model(cache, "key", function() ms_fake_stan_model(1L))
  expect_identical(model$id, 1L)
  expect_identical(list.files(dir), character())

  off <- markov_stan_cache_dir(FALSE)
  expect_false(off$use)
  expect_match(markov_stan_compile_message(off), "done once per session[.]$")
})

test_that("the write check finds a writable directory", {
  dir <- withr::local_tempdir()
  expect_true(markov_dir_writable(dir))
  expect_identical(list.files(dir), character())
  expect_false(markov_dir_writable(file.path(dir, "does-not-exist")))
  expect_match(
    markov_stan_compile_message(markov_stan_cache_dir(dir)),
    "the result is kept in"
  )
})
