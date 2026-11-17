# Analytical SOP result accessors.

delta_row_key_frame <- function(x) {
  key_cols <- setdiff(names(x), inference_columns())
  out <- as.data.frame(x)[, key_cols, drop = FALSE]
  rownames(out) <- NULL
  out
}

delta_row_key <- function(keys) {
  if (ncol(keys) == 0L) {
    return(as.character(seq_len(nrow(keys))))
  }
  do.call(
    paste,
    c(
      lapply(keys, function(value) {
        value <- as.character(value)
        value[is.na(value)] <- "<NA>"
        value
      }),
      sep = "\r"
    )
  )
}

delta_storage_bytes <- function(analytical) {
  component_bytes <- c(
    jacobian = if (is.null(analytical$jacobian)) {
      0
    } else {
      length(analytical$jacobian) * 8
    },
    average_jacobian = if (is.null(analytical$average_jacobian)) {
      0
    } else {
      length(analytical$average_jacobian) * 8
    },
    coefficient_vcov = if (is.null(analytical$coefficient_vcov)) {
      0
    } else {
      length(analytical$coefficient_vcov) * 8
    },
    influence = if (is.null(analytical$influence)) {
      0
    } else {
      length(analytical$influence) * 8
    }
  )
  list(
    components = component_bytes,
    stored = unname(sum(component_bytes)),
    limit = delta_max_bytes()
  )
}

delta_attach_byte_accounting <- function(analytical) {
  analytical$bytes <- delta_storage_bytes(analytical)
  delta_assert_bytes(
    analytical$bytes$stored,
    "Stored analytical inference state"
  )
  analytical
}

delta_analytical <- function(x) {
  if (
    !inherits(
      x,
      c(
        "markov_sops",
        "markov_avg_sops",
        "markov_comparisons",
        "markov_avg_comparisons",
        "markov_state_time",
        "markov_avg_state_time"
      )
    )
  ) {
    mostr_abort(
      "Analytical accessors require a SOP, average-time, or comparison object."
    )
  }
  analytical <- attr(x, "analytical", exact = TRUE)
  if (is.null(analytical)) {
    mostr_abort(
      "No analytical inference state found. Run `inferences(method = ",
      "\"delta\")` first."
    )
  }
  analytical
}

# Map each current row of `object` to the stored row that holds its
# derivatives. Base R `[` keeps the "analytical" attribute when a user filters
# or sorts a delta result, so the stored rows can differ from the rows the
# object has now; without this map, vcov() returned covariances of rows the
# object no longer has, in the stored order (issue #43). Keys are rebuilt from
# the stored key column names, so added or reordered columns still match.
delta_current_rows <- function(object, analytical) {
  mismatch <- function() {
    mostr_abort(
      "The rows of this result no longer match the values stored by ",
      "`inferences(method = \"delta\")`. Filtering, sorting, or repeating ",
      "rows is supported, but editing the columns that identify a row is ",
      "not. Rerun `inferences()` on the edited result."
    )
  }
  key_cols <- names(analytical$row_keys)
  stored_key <- analytical$row_key
  if (length(key_cols) == 0L) {
    # Without identifying columns the stored keys are row positions, so a
    # subset cannot be matched. Only the unchanged row count is accepted.
    if (nrow(object) != length(stored_key)) {
      mismatch()
    }
    return(list(key = stored_key, stored = seq_along(stored_key)))
  }
  if (!all(key_cols %in% names(object))) {
    mismatch()
  }
  current_key <- delta_row_key(
    as.data.frame(object)[, key_cols, drop = FALSE]
  )
  stored <- match(current_key, stored_key)
  if (anyNA(stored)) {
    mismatch()
  }
  list(key = current_key, stored = stored)
}

# Resolve `rows` against the current rows of the object and return positions
# among those current rows.
delta_resolve_rows <- function(current, rows = NULL) {
  n <- length(current$key)
  if (is.null(rows)) {
    return(seq_len(n))
  }
  if (is.logical(rows)) {
    if (length(rows) != n || anyNA(rows)) {
      mostr_abort(
        "Logical `rows` must be non-missing and match the number of rows ",
        "of the result."
      )
    }
    return(which(rows))
  }
  if (is.character(rows)) {
    index <- match(rows, current$key)
    if (anyNA(index)) {
      mostr_abort("Some requested analytical row keys were not found.")
    }
    return(index)
  }
  if (is.numeric(rows)) {
    if (
      anyNA(rows) ||
        any(!is.finite(rows)) ||
        any(rows != as.integer(rows)) ||
        any(rows < 1L) ||
        any(rows > n)
    ) {
      mostr_abort("Numeric `rows` must contain valid positive row indices.")
    }
    return(as.integer(rows))
  }
  mostr_abort("`rows` must be NULL, numeric, logical, or analytical row keys.")
}

# Internal inspection helper: the unconditional Jacobian contains only the
# coefficient derivative, not the complete patient influence representation.
get_jacobian <- function(x, rows = NULL) {
  analytical <- delta_analytical(x)
  current <- delta_current_rows(x, analytical)
  selected <- delta_resolve_rows(current, rows)
  index <- current$stored[selected]
  jacobian <- analytical$jacobian %||% analytical$average_jacobian
  if (is.null(jacobian)) {
    mostr_abort("The analytical result does not contain an average Jacobian.")
  }
  required <- as.double(length(index)) * as.double(ncol(jacobian)) * 8
  delta_assert_bytes(required, "The requested Jacobian block")
  out <- jacobian[index, , drop = FALSE]
  rownames(out) <- current$key[selected]
  out
}

delta_vcov <- function(object, rows = NULL) {
  analytical <- delta_analytical(object)
  current <- delta_current_rows(object, analytical)
  selected <- delta_resolve_rows(current, rows)
  index <- current$stored[selected]
  m <- length(index)

  if (identical(analytical$representation, "influence")) {
    n <- nrow(analytical$influence)
    # A sample covariance over patients needs at least two of them.
    if (n < 2L) {
      mostr_abort(
        "Unconditional analytical covariance requires at least two patients."
      )
    }
    required <- (as.double(n) * m + as.double(m) * m) * 8
    delta_assert_bytes(required, "The requested analytical covariance block")
    influence <- analytical$influence[, index, drop = FALSE]
    # The estimator of delta_stacked_influence(): the sample covariance of the
    # patients' influence values divided by n, so the diagonal reproduces the
    # stored standard errors.
    out <- stats::cov(influence) / n
  } else {
    q <- ncol(analytical$jacobian)
    required <- (as.double(m) * q + as.double(m) * m) * 8
    delta_assert_bytes(required, "The requested analytical covariance block")
    jacobian <- analytical$jacobian[index, , drop = FALSE]
    coefficient_vcov <- analytical$coefficient_vcov
    out <- jacobian %*% coefficient_vcov %*% t(jacobian)
    out <- (out + t(out)) / 2
  }

  dimnames(out) <- list(
    current$key[selected],
    current$key[selected]
  )
  out
}

#' Analytical Covariance for SOPs, Average Times, and Comparisons
#'
#' Materializes a selected covariance block from the low-rank analytical state
#' retained by `inferences(method = "delta")`.
#' Returns the conditional or unconditional covariance selected by the `vcov`
#' argument to [inferences()]; it does not change that choice. Conditional
#' covariance accounts for coefficient estimation with prediction covariates
#' treated as given. Unconditional covariance also accounts for sampling the
#' patients used in the average. Bootstrap and MVN results are not supported.
#' With `rows = NULL`, the full dense covariance matrix is requested and its
#' storage grows quadratically with the number of result cells. Use `rows` to
#' extract smaller blocks within the analytical memory limit. Standard errors
#' are computed without materializing the full covariance matrix.
#' You may filter or sort the result before calling `vcov()`: the matrix
#' follows the rows of `object` as you pass it, so the square roots of its
#' diagonal entries equal the `std.error` column.
#'
#' @param object A `markov_sops`, `markov_avg_sops`, `markov_state_time`,
#'   `markov_avg_state_time`, `markov_comparisons`, or `markov_avg_comparisons`
#'   object returned by `inferences(method = "delta")`, possibly with rows
#'   filtered, sorted, or repeated.
#' @param rows Optional selection of rows of `object` as passed: row numbers,
#'   a logical vector with one value per row of `object`, or the row names of
#'   a matrix returned by an earlier `vcov()` call. `NULL` returns every row.
#'   Rows whose identifying columns (such as `time` or `state`) were changed
#'   after `inferences()` raise an error; rerun `inferences()` instead.
#' @param ... Reserved for the `stats::vcov()` generic.
#'
#' @return A symmetric covariance matrix for the selected result cells.
#' @export
vcov.markov_sops <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_avg_sops <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_comparisons <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_avg_comparisons <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_state_time <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_avg_state_time <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}
