# Numerical derivatives by central differences, for checking the analytical
# derivatives. They re-run the whole calculation with nudged coefficients, so
# they share no derivative code with the analytical recursion.

# How each number returned by `f(coefficients)` changes with each coefficient.
# `f` takes a coefficient vector and returns the numbers to differentiate, in a
# fixed order, as a vector or an array. Coefficient j is moved by
# `step * max(1, abs(coefficients[j]))` up and down, and the derivative is
# (f(upper) - f(lower)) / (2 * that move).
#
# Returns a matrix with one row per number and one column per coefficient,
# named after the coefficients. When `f` returns an array, the result is an
# array with the coefficients as its last dimension.
central_difference_jacobian <- function(
  f,
  coefficients,
  step = .Machine$double.eps^(1 / 3)
) {
  columns <- lapply(seq_along(coefficients), function(j) {
    move <- step * max(1, abs(coefficients[[j]]))
    upper <- coefficients
    lower <- coefficients
    upper[j] <- upper[j] + move
    lower[j] <- lower[j] - move
    (f(upper) - f(lower)) / (2 * move)
  })
  value_dim <- dim(columns[[1L]])
  if (is.null(value_dim)) {
    value_dim <- length(columns[[1L]])
  }
  out <- array(
    unlist(lapply(columns, as.vector)),
    dim = c(value_dim, length(coefficients))
  )
  if (length(value_dim) == 1L) {
    colnames(out) <- names(coefficients)
  }
  out
}

# The same for a fitted model: `estimates(model)` returns the numbers to
# differentiate, and the model's coefficients are nudged with `set_coef()`.
numeric_jacobian <- function(model, estimates, step = .Machine$double.eps^(1 / 3)) {
  central_difference_jacobian(
    function(coefficients) estimates(set_coef(model, coefficients)),
    get_coef(model),
    step = step
  )
}
