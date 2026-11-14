#include <cpp11.hpp>
#include <R_ext/Arith.h>
#include <vector>

#include "mostr/path_likelihood.hpp"

using namespace cpp11;

// The log-likelihood of a proportional-odds Markov model when some visits
// have no exactly observed state, and its exact derivatives. The
// implementation is shared with the Stan program of engine = "blrm" and
// lives in inst/include/mostr/path_likelihood.hpp, which documents the
// model, the layouts of the design rows and the recursion. This file only
// converts R objects: matrices are copied to the row-major layout the shared
// code reads, and its outputs are returned as R matrices.

namespace {

// Row-major copy of an R matrix (stored by column).
std::vector<double> row_major(const doubles_matrix<>& x) {
  const int rows = x.nrow();
  const int cols = x.ncol();
  std::vector<double> out(static_cast<size_t>(rows) * cols);
  if (rows == 0 || cols == 0) return out;
  const double* ptr = REAL(SEXP(x));
  for (int col = 0; col < cols; ++col) {
    const double* column = ptr + static_cast<R_xlen_t>(rows) * col;
    for (int row = 0; row < rows; ++row) {
      out[static_cast<size_t>(row) * cols + col] = column[row];
    }
  }
  return out;
}

std::vector<int> row_major(const integers_matrix<>& x) {
  const int rows = x.nrow();
  const int cols = x.ncol();
  std::vector<int> out(static_cast<size_t>(rows) * cols);
  if (rows == 0 || cols == 0) return out;
  const int* ptr = INTEGER(SEXP(x));
  for (int col = 0; col < cols; ++col) {
    const int* column = ptr + static_cast<R_xlen_t>(rows) * col;
    for (int row = 0; row < rows; ++row) {
      // NA_INTEGER is the most negative int, so the shared validation
      // refuses it as an index.
      out[static_cast<size_t>(row) * cols + col] = column[row];
    }
  }
  return out;
}

const double* real_pointer(SEXP x) {
  return Rf_xlength(x) == 0 ? nullptr : REAL(x);
}

const int* integer_pointer(SEXP x) {
  return Rf_xlength(x) == 0 ? nullptr : INTEGER(x);
}

}  // namespace

// Index conventions (all indices 1-based, as R passes them):
// - y_obs in 1..K; group_obs in 1..n_groups.
// - absorb: absorbing states in 1..K (may be empty).
// - x_gap: design rows of the stretches; x_prev: a K-row matrix added to
//   them by previous state (additive layout), or a matrix with 0 rows
//   (general layout); stretch_first_prev: the known previous state of each
//   stretch's step 1, used only with x_prev.
// - The other stretch and step arguments are those of
//   mostr::paths::Problem; step_rows is an n_steps x K matrix.
[[cpp11::register]] list cpp_path_loglik(
    doubles alpha,
    doubles beta,
    doubles_matrix<> x_obs,
    integers y_obs,
    doubles w_obs,
    integers group_obs,
    int n_groups,
    doubles_matrix<> x_gap,
    doubles_matrix<> x_prev,
    integers absorb,
    integers stretch_first_row,
    integers stretch_first_prev,
    integers stretch_first_step,
    integers stretch_n_steps,
    doubles stretch_weight,
    integers step_lo,
    integers step_hi,
    integers_matrix<> step_rows,
    int order) {
  const int thresholds = alpha.size();
  const int states = thresholds + 1;
  const int p = beta.size();
  const int q = thresholds + p;
  const int n_obs = x_obs.nrow();
  const int n_gap = x_gap.nrow();
  const int m = stretch_first_row.size();
  const int n_steps = step_lo.size();
  const bool has_prev = x_prev.nrow() > 0;

  // Shapes; the shared code checks the indices.
  if (thresholds < 1) cpp11::stop("alpha must have length >= 1.");
  if (order < 0 || order > 2) cpp11::stop("order must be 0, 1 or 2.");
  if (n_obs > 0 && x_obs.ncol() != p) cpp11::stop("x_obs needs length(beta) columns.");
  if (n_gap > 0 && x_gap.ncol() != p) cpp11::stop("x_gap needs length(beta) columns.");
  if (has_prev && (x_prev.nrow() != states || x_prev.ncol() != p)) {
    cpp11::stop("x_prev needs K rows and length(beta) columns, or no rows.");
  }
  if (y_obs.size() != n_obs || w_obs.size() != n_obs || group_obs.size() != n_obs) {
    cpp11::stop("y_obs, w_obs and group_obs must have nrow(x_obs) entries.");
  }
  if (stretch_first_step.size() != m || stretch_n_steps.size() != m ||
      stretch_weight.size() != m ||
      (has_prev && stretch_first_prev.size() != m)) {
    cpp11::stop("Stretch vectors must all have length(stretch_first_row) entries.");
  }
  if (step_hi.size() != n_steps || step_rows.nrow() != n_steps ||
      (n_steps > 0 && step_rows.ncol() != states)) {
    cpp11::stop("step_hi and step_rows must match step_lo; step_rows needs K columns.");
  }
  std::vector<int> is_absorbing(states, 0);
  for (R_xlen_t pos = 0; pos < absorb.size(); ++pos) {
    const int state = absorb[pos];
    if (state == NA_INTEGER || state < 1 || state > states) {
      cpp11::stop("absorb must lie in 1..K.");
    }
    is_absorbing[state - 1] = 1;
  }

  const std::vector<double> x_obs_rows = row_major(x_obs);
  const std::vector<double> x_gap_rows = row_major(x_gap);
  const std::vector<double> x_prev_rows = row_major(x_prev);
  const std::vector<int> step_rows_by_step = row_major(step_rows);

  mostr::paths::Problem d;
  d.thresholds = thresholds;
  d.p = p;
  d.n_obs = n_obs;
  d.x_obs = x_obs_rows.data();
  d.y_obs = integer_pointer(y_obs);
  d.w_obs = real_pointer(w_obs);
  d.group_obs = integer_pointer(group_obs);
  d.n_groups = n_groups;
  d.n_gap = n_gap;
  d.x_gap = x_gap_rows.data();
  d.x_prev = has_prev ? x_prev_rows.data() : nullptr;
  d.absorbing = is_absorbing.data();
  d.m = m;
  d.first_row = integer_pointer(stretch_first_row);
  d.first_prev = has_prev ? integer_pointer(stretch_first_prev) : nullptr;
  d.first_step = integer_pointer(stretch_first_step);
  d.n_steps = integer_pointer(stretch_n_steps);
  d.stretch_weight = real_pointer(stretch_weight);
  d.n_steps_total = n_steps;
  d.lo = integer_pointer(step_lo);
  d.hi = integer_pointer(step_hi);
  d.step_rows = step_rows_by_step.data();
  // y_obs, group_obs, step_lo and step_hi may hold NA_INTEGER, the most
  // negative int, which the range checks refuse like any other bad index.
  try {
    mostr::paths::validate(d, true);
  } catch (const std::exception& e) {
    cpp11::stop("%s", e.what());
  }

  writable::doubles stretch_loglik(m);
  writable::doubles_matrix<> group_scores(order >= 1 ? n_groups : 0, order >= 1 ? q : 0);
  writable::doubles_matrix<> stretch_scores(order >= 1 ? m : 0, order >= 1 ? q : 0);
  writable::doubles_matrix<> hessian(order >= 2 ? q : 0, order >= 2 ? q : 0);

  mostr::paths::Output out;
  out.order = order;
  out.stretch_loglik = m > 0 ? REAL(SEXP(stretch_loglik)) : nullptr;
  if (order >= 1) {
    out.group_scores = n_groups > 0 ? REAL(SEXP(group_scores)) : nullptr;
    out.stretch_scores = m > 0 ? REAL(SEXP(stretch_scores)) : nullptr;
  }
  if (order >= 2) out.hessian = REAL(SEXP(hessian));
  mostr::paths::evaluate(d, real_pointer(alpha), real_pointer(beta), out);

  writable::list result;
  result.push_back(cpp11::named_arg("loglik") = writable::doubles({out.loglik}));
  result.push_back(cpp11::named_arg("stretch_loglik") = stretch_loglik);
  if (order >= 1) {
    result.push_back(cpp11::named_arg("obs_group_scores") = group_scores);
    result.push_back(cpp11::named_arg("stretch_scores") = stretch_scores);
  }
  if (order >= 2) {
    result.push_back(cpp11::named_arg("hessian") = hessian);
  }
  return result;
}
