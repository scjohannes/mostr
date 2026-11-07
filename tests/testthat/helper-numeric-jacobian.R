# How each estimate changes with each model coefficient, by central differences.
# `estimates(model)` returns the numbers to differentiate, in a fixed order. This
# re-runs the whole calculation with nudged coefficients, so it shares no
# derivative code with the analytical recursion.
numeric_jacobian <- function(model, estimates) {
  coefficient <- get_coef(model)
  columns <- lapply(seq_along(coefficient), function(column) {
    step <- .Machine$double.eps^(1 / 3) * max(1, abs(coefficient[column]))
    upper <- coefficient
    lower <- coefficient
    upper[column] <- upper[column] + step
    lower[column] <- lower[column] - step
    (estimates(set_coef(model, upper)) - estimates(set_coef(model, lower))) /
      (2 * step)
  })
  out <- do.call(cbind, columns)
  colnames(out) <- names(coefficient)
  out
}
