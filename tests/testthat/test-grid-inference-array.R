# Inference checks on a pairwise covering array of the remaining settings:
# every pair of values of any two settings appears in at least one row. The
# settings are the result type, the method, the engine, the Markov order,
# full or partial proportional odds, case weights, stored profiles or
# `newdata`, `by`, the visit scale or `time_map`, and the comparison (a
# difference, a ratio or a time benefit). Time is linear throughout; the full
# crossing of engine, time coding and method is in test-grid-inference.R, and
# the checks are in helper-grid-inference.R.
#
# "-" marks a setting that does not apply: `by` to results with one row per
# profile, the time scale to SOPs, and the comparison to anything but
# comparisons. For comparisons, the visit scale means SOP differences or
# ratios, and `time_map` means time in state in real time; time benefit is
# on either scale. Analytical (delta) rows use unconditional variance for
# averages of the stored, unweighted patients and conditional variance
# otherwise. Weighted fits therefore get conditional variance only.
#
# Pairs the package refuses are left out and asserted in test-refusals.R:
# the ordinary bootstrap of results with one row per profile (I42 to I44),
# and analytical inference with second-order fits (I40), partial
# proportional odds (I41), `by` (I27, I28), ratios (I29, I30) and time
# benefit (I31, I32). Partial proportional odds are fitted with vglm only.
#
# The rows were chosen greedily: each adds the row that covers the most pairs
# not yet covered. The last test re-checks that every allowed pair is covered.

inference_array <- utils::read.table(
  header = TRUE,
  colClasses = "character",
  text = "
  result           method     engine  order  odds     weights  profiles  by      scale     comparison
  avg_comparisons  mvn        orm     1      full     none     stored    none    visit     difference
  avg_comparisons  bootstrap  vglm    2      partial  integer  newdata   marker  time_map  ratio
  comparisons      fwb        vglm    2      partial  none     stored    -       visit     time_benefit
  avg_comparisons  fwb        orm     1      full     integer  newdata   none    time_map  time_benefit
  avg_state_time   delta      vglm    1      full     integer  stored    none    visit     -
  avg_state_time   bootstrap  orm     2      full     none     newdata   marker  visit     -
  comparisons      mvn        vglm    1      partial  integer  newdata   -       time_map  difference
  state_time       delta      orm     1      full     none     stored    -       time_map  -
  avg_sops         bootstrap  vglm    1      partial  none     stored    none    -         -
  comparisons      mvn        orm     1      full     none     stored    -       visit     ratio
  avg_sops         mvn        orm     2      full     integer  stored    marker  -         -
  sops             delta      orm     1      full     none     newdata   -       -         -
  state_time       mvn        vglm    2      partial  integer  newdata   -       visit     -
  sops             mvn        vglm    2      partial  integer  stored    -       -         -
  avg_state_time   fwb        vglm    1      partial  none     stored    marker  time_map  -
  avg_comparisons  bootstrap  orm     2      full     none     stored    none    visit     difference
  comparisons      delta      orm     1      full     none     stored    -       visit     difference
  avg_sops         fwb        orm     1      full     none     newdata   none    -         -
  avg_comparisons  fwb        orm     1      full     none     stored    marker  visit     difference
  avg_comparisons  fwb        orm     1      full     none     stored    none    visit     ratio
  avg_comparisons  mvn        orm     1      full     none     stored    marker  visit     time_benefit
  avg_state_time   mvn        orm     1      full     none     stored    none    visit     -
  sops             fwb        orm     1      full     none     stored    -       -         -
  state_time       fwb        orm     1      full     none     stored    -       visit     -
  avg_sops         delta      orm     1      full     none     stored    none    -         -
  avg_comparisons  delta      orm     1      full     none     stored    none    visit     difference
  avg_comparisons  bootstrap  orm     1      full     none     stored    none    visit     time_benefit
"
)

inference_array_settings <- function(rows = seq_len(nrow(inference_array))) {
  lapply(rows, function(i) do.call(inference_setting, as.list(inference_array[i, ])))
}

# Rows for CRAN: MVN or analytical inference on the base orm fit, one each
# for an averaged comparison, time for each profile, SOPs for supplied
# profiles and averaged time.
inference_array_cran_rows <- c(1L, 8L, 12L, 22L)

test_that("MVN and analytical inference hold their limits across the covering array", {
  skip_on_cran()
  # L9, L10, L11, L13, L18, L21, L23, L24, R13, R14, N1, N4, N21.
  rows <- which(inference_array$method %in% c("mvn", "delta"))
  if (grid_on_cran()) {
    rows <- inference_array_cran_rows
  }
  for (setting in inference_array_settings(rows)) check_inference_setting(setting)
})

test_that("bootstrap inference holds its limits across the covering array", {
  # As above, plus L29 bootstrap bookkeeping.
  skip_on_cran()
  rows <- which(inference_array$method %in% c("bootstrap", "fwb"))
  for (setting in inference_array_settings(rows)) check_inference_setting(setting)
})

test_that("the covering array covers every allowed pair of settings", {
  values <- list(
    result = c("sops", "avg_sops", "state_time", "avg_state_time", "comparisons", "avg_comparisons"),
    method = c("mvn", "bootstrap", "fwb", "delta"),
    engine = c("orm", "vglm"),
    order = c("1", "2"),
    odds = c("full", "partial"),
    weights = c("none", "integer"),
    profiles = c("stored", "newdata"),
    by = c("none", "marker"),
    scale = c("visit", "time_map"),
    comparison = c("difference", "ratio", "time_benefit")
  )
  # Every allowed row, with settings that do not apply set to "-".
  all_rows <- expand.grid(values, stringsAsFactors = FALSE)
  averaged <- startsWith(all_rows$result, "avg_")
  compares <- all_rows$result %in% c("comparisons", "avg_comparisons")
  all_rows$by[!averaged] <- "-"
  all_rows$scale[all_rows$result %in% c("sops", "avg_sops")] <- "-"
  all_rows$comparison[!compares] <- "-"
  allowed <- with(
    all_rows,
    !(odds == "partial" & engine == "orm") &
      !(method == "bootstrap" & !averaged) &
      !(method == "delta" & (order == "2" | odds == "partial" | by == "marker" |
        comparison %in% c("ratio", "time_benefit")))
  )
  all_rows <- unique(all_rows[allowed, ])
  pairs <- function(rows) {
    out <- character()
    for (i in seq_along(values)) {
      for (j in seq_along(values)) {
        if (j <= i) next
        a <- rows[[i]]
        b <- rows[[j]]
        keep <- a != "-" & b != "-"
        out <- c(out, paste(names(values)[i], a[keep], names(values)[j], b[keep]))
      }
    }
    unique(out)
  }
  table <- inference_array[names(values)]
  # Each row of the table is allowed, and together they cover every pair.
  expect_true(all(grid_key(table, names(values)) %in% grid_key(all_rows, names(values))))
  expect_identical(setdiff(pairs(all_rows), pairs(table)), character())
})

# Free the cached results; other files do not use them.
grid_cache$results <- NULL
