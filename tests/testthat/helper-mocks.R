# Hand-made stand-ins for fitted models, shared by the test files.

# Records settings on a hand-made fit as `mostm()` records them on a real one
# (`markov_settings()` with every field fixed by the fit). The prediction
# functions take the patient ID, time and previous-state columns, the order
# and the states from the fit, so a hand-made fit that tests other columns, a
# second-order recursion or particular state labels records them here, and
# so does one with absorbing states (`absorb`) or a gap column (`gap_var`).
# `y_levels` defaults to the states the object itself records.
with_fit_settings <- function(
  model,
  id_var = markov_model_id_var(model),
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  y_levels = markov_model_ylevels_or_null(model),
  absorb = NULL,
  gap_var = NULL
) {
  settings <- markov_settings(
    id_var = id_var,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    from_fit = setdiff(
      markov_settings_fit_fields(),
      if (is.null(y_levels)) "y_levels"
    )
  )
  markov_set_settings(model, settings)
}

# A blrm fit as `mostm(engine = "blrm")` makes it. Sampling from Stan is too slow
# for the default run, so the posterior draws are real numbers given here, and
# tests mock only the design matrix (see `fake_blrm_design()`).
# test-slow-blrm.R covers the design matrix with real fits.
#
# `draws`: one row per posterior draw. The intercepts come first, named
#   "y>=2", "y>=3", ..., which also fixes the number of states.
# `tau`: names of the columns of `draws` that hold the partial
#   proportional-odds effects (none for a proportional-odds fit).
# `cppo`: the function that constrains those effects across cut-points, or
#   NULL for an unconstrained partial proportional-odds fit.
# `cluster`: FALSE for a fit without random effects for patients "a" and "b".
# `gamma_draws`: random-effect draws, so that no Stan fit is needed to read
#   them.
# `class`: the class of the fit; c("blrm", "orm") as for rmsb fits.
# `p2_var`: the column with the state two visits back for a second-order fit,
#   recorded as `mostm(p2_var = )` records it.
# `absorb`: the absorbing states, recorded as `mostm(absorb = )` records them.
fake_blrm <- function(
  draws,
  tau = character(),
  cppo = NULL,
  cluster = TRUE,
  gamma_draws = NULL,
  class = c("blrm", "orm"),
  p2_var = NULL,
  absorb = NULL
) {
  non_slopes <- sum(startsWith(colnames(draws), "y>="))
  fit <- list(
    draws = draws,
    non.slopes = non_slopes,
    pppo = length(tau),
    tauInfo = data.frame(name = tau),
    cppo = cppo,
    ylevels = seq_len(non_slopes + 1L),
    yname = "y",
    clusterInfo = if (cluster) list(name = "id", cluster = c("a", "b")),
    gamma_draws = gamma_draws
  )
  fit <- structure(
    Filter(Negate(is.null), fit),
    class = class,
    markov_engine = "blrm",
    markov_id_var = "id"
  )
  with_fit_settings(fit, p2_var = p2_var, absorb = absorb)
}

# A stand-in for `blrm_design_matrix()`: the design matrix is the `columns` of
# `newdata` as numbers (factor levels are read as numbers), and the design of
# the partial proportional-odds effect, asked for with `second = TRUE`, is
# column `tau_var`, named "tau_<tau_var>".
fake_blrm_design <- function(columns, tau_var = NULL) {
  as_number <- function(x) {
    if (is.factor(x) || is.character(x)) as.numeric(as.character(x)) else x
  }
  function(model, newdata, second = FALSE) {
    if (second) {
      return(matrix(
        as_number(newdata[[tau_var]]),
        ncol = 1,
        dimnames = list(NULL, paste0("tau_", tau_var))
      ))
    }
    if (length(columns) == 0L) {
      return(matrix(numeric(0), nrow = nrow(newdata), ncol = 0))
    }
    out <- do.call(cbind, lapply(columns, function(name) as_number(newdata[[name]])))
    colnames(out) <- columns
    out
  }
}
