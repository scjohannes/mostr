#include <cpp11.hpp>
#include <R_ext/Arith.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace cpp11;

// Log-likelihood of a proportional-odds Markov transition model when some
// visits have no exactly observed state. Observed transitions enter as
// ordinary cell probabilities. A run of visits whose states are known only to
// lie in a range of states (a "stretch"; an unobserved visit is the range of
// states it could be in) enters as the probability of all state paths
// consistent with the ranges, summed by a forward recursion over the states
// in each step's range. Derivatives with respect to
// theta = c(alpha, beta) are exact: the recursion carries the gradient and
// Hessian of its forward vector (forward-mode differentiation).
//
// Model: for a design row x and cut-point k = 1..K-1,
//   c_k = plogis(alpha_k + x'beta) = P(Y >= k + 1),
//   P(Y = y) = c_{y-1} - c_y with c_0 = 1 and c_K = 0.
// With g_k = (e_k, x), the derivatives are dc_k = s_k g_k and
// d2c_k = u_k g_k g_k', where s_k = c_k (1 - c_k) and u_k = s_k (1 - 2 c_k).

namespace {

double logistic(const double value) {
  if (value >= 0.0) {
    return 1.0 / (1.0 + std::exp(-value));
  }
  const double exponent = std::exp(value);
  return exponent / (1.0 + exponent);
}

// Cut-point quantities of one design row. Both c_k and 1 - c_k are evaluated
// by their own logistic, because a cell probability c_{y-1} - c_y computed
// from c alone loses all its digits when both are close to 1.
struct CutPoints {
  std::vector<double> upper;  // c_k
  std::vector<double> lower;  // 1 - c_k

  explicit CutPoints(const int thresholds)
      : upper(thresholds), lower(thresholds) {}

  void set(const double* alpha, const int thresholds, const double eta) {
    for (int k = 0; k < thresholds; ++k) {
      upper[k] = logistic(alpha[k] + eta);
      lower[k] = logistic(-(alpha[k] + eta));
    }
  }
};

// P(Y = y) for 1-based y from the cut-point values c (upper) and 1 - c
// (lower) of one row. Whichever form subtracts the smaller numbers is used,
// so the difference keeps its relative precision in both tails.
double cell_probability(
    const double* upper,
    const double* lower,
    const int thresholds,
    const int y) {
  const int lo = y - 1;  // cut-point index of c_{y-1}, 1-based; 0 means c_0
  const int hi = y;      // cut-point index of c_y, 1-based; K means c_K
  const double c_lo = lo >= 1 ? upper[lo - 1] : 1.0;
  const double c_hi = hi <= thresholds ? upper[hi - 1] : 0.0;
  const double q_lo = lo >= 1 ? lower[lo - 1] : 0.0;
  const double q_hi = hi <= thresholds ? lower[hi - 1] : 1.0;
  return c_lo > 0.5 ? q_hi - q_lo : c_lo - c_hi;
}

// s_k = c_k (1 - c_k) and u_k = s_k (1 - 2 c_k) for 1-based k; 0 outside
// 1..K-1, where c_k is the constant 1 or 0.
double slope_at(const double* upper, const double* lower, const int thresholds,
                const int k) {
  if (k < 1 || k > thresholds) {
    return 0.0;
  }
  return upper[k - 1] * lower[k - 1];
}

double curvature_at(const double* upper, const double* lower,
                    const int thresholds, const int k) {
  if (k < 1 || k > thresholds) {
    return 0.0;
  }
  return upper[k - 1] * lower[k - 1] * (lower[k - 1] - upper[k - 1]);
}

const double* real_pointer(SEXP x) {
  return Rf_xlength(x) == 0 ? nullptr : REAL(x);
}

const int* integer_pointer(SEXP x) {
  return Rf_xlength(x) == 0 ? nullptr : INTEGER(x);
}

writable::doubles_matrix<> to_matrix(const std::vector<double>& values,
                                     const int rows, const int cols) {
  writable::doubles_matrix<> out(rows, cols);
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      out(row, col) = values[static_cast<size_t>(row) +
                             static_cast<size_t>(rows) * col];
    }
  }
  return out;
}

}  // namespace

// Index conventions (all indices 1-based, as R passes them):
// - y_obs in 1..K; group_obs in 1..n_groups.
// - absorb: absorbing states in 1..K (may be empty). From an absorbing state
//   the transition is the identity and uses no design row.
// - Stretch s covers the steps stretch_first_step[s] ..
//   stretch_first_step[s] + stretch_n_steps[s] - 1 of the step arrays.
// - stretch_first_row[s]: row of x_gap for step 1, whose previous state is
//   known (and non-absorbing); the step_rows row of a first step is ignored.
// - step_lo[u], step_hi[u]: at step u the state lies in step_lo..step_hi.
// - step_rows[u, j]: row of x_gap for step u when the previous state is j;
//   0 where unused. It must be a valid row for every non-absorbing j in the
//   previous step's range.
// Rows and stretches with weight 0 do not enter loglik or the Hessian, so an
// impossible zero-weight row does not turn the total into NaN.
[[cpp11::register]] list cpp_path_loglik(
    doubles alpha,
    doubles beta,
    doubles_matrix<> x_obs,
    integers y_obs,
    doubles w_obs,
    integers group_obs,
    int n_groups,
    doubles_matrix<> x_gap,
    integers absorb,
    integers stretch_first_row,
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

  // Cheap validation of every index once, so the loops can trust them.
  if (thresholds < 1) cpp11::stop("alpha must have length >= 1.");
  if (order < 0 || order > 2) cpp11::stop("order must be 0, 1 or 2.");
  if (n_groups < 0) cpp11::stop("n_groups must be non-negative.");
  if (n_obs > 0 && x_obs.ncol() != p) cpp11::stop("x_obs needs length(beta) columns.");
  if (n_gap > 0 && x_gap.ncol() != p) cpp11::stop("x_gap needs length(beta) columns.");
  if (y_obs.size() != n_obs || w_obs.size() != n_obs || group_obs.size() != n_obs) {
    cpp11::stop("y_obs, w_obs and group_obs must have nrow(x_obs) entries.");
  }
  if (stretch_first_step.size() != m || stretch_n_steps.size() != m ||
      stretch_weight.size() != m) {
    cpp11::stop("Stretch vectors must all have length(stretch_first_row) entries.");
  }
  if (step_hi.size() != n_steps || step_rows.nrow() != n_steps ||
      (n_steps > 0 && step_rows.ncol() != states)) {
    cpp11::stop("step_hi and step_rows must match step_lo; step_rows needs K columns.");
  }
  const int* y_ptr = integer_pointer(y_obs);
  const int* group_ptr = integer_pointer(group_obs);
  for (int i = 0; i < n_obs; ++i) {
    if (y_ptr[i] == NA_INTEGER || y_ptr[i] < 1 || y_ptr[i] > states) {
      cpp11::stop("y_obs must lie in 1..K.");
    }
    if (group_ptr[i] == NA_INTEGER || group_ptr[i] < 1 || group_ptr[i] > n_groups) {
      cpp11::stop("group_obs must lie in 1..n_groups.");
    }
  }
  std::vector<char> is_absorbing(states, 0);
  for (R_xlen_t pos = 0; pos < absorb.size(); ++pos) {
    const int state = absorb[pos];
    if (state == NA_INTEGER || state < 1 || state > states) {
      cpp11::stop("absorb must lie in 1..K.");
    }
    is_absorbing[state - 1] = 1;
  }
  const int* lo_ptr = integer_pointer(step_lo);
  const int* hi_ptr = integer_pointer(step_hi);
  for (int u = 0; u < n_steps; ++u) {
    if (lo_ptr[u] == NA_INTEGER || hi_ptr[u] == NA_INTEGER || lo_ptr[u] < 1 ||
        hi_ptr[u] < lo_ptr[u] || hi_ptr[u] > states) {
      cpp11::stop("step_lo and step_hi must satisfy 1 <= lo <= hi <= K.");
    }
  }
  const int* step_ptr = integer_pointer(SEXP(step_rows));
  const int* first_row_ptr = integer_pointer(stretch_first_row);
  const int* first_step_ptr = integer_pointer(stretch_first_step);
  const int* n_steps_ptr = integer_pointer(stretch_n_steps);
  for (int s = 0; s < m; ++s) {
    if (first_row_ptr[s] == NA_INTEGER || first_row_ptr[s] < 1 ||
        first_row_ptr[s] > n_gap) {
      cpp11::stop("stretch_first_row must lie in 1..nrow(x_gap).");
    }
    if (n_steps_ptr[s] == NA_INTEGER || n_steps_ptr[s] < 1 ||
        first_step_ptr[s] == NA_INTEGER || first_step_ptr[s] < 1 ||
        first_step_ptr[s] > n_steps - n_steps_ptr[s] + 1) {
      cpp11::stop("stretch_first_step and stretch_n_steps point outside the step arrays.");
    }
    for (int t = 1; t < n_steps_ptr[s]; ++t) {
      const int u = first_step_ptr[s] - 1 + t;
      for (int j = lo_ptr[u - 1]; j <= hi_ptr[u - 1]; ++j) {
        if (is_absorbing[j - 1]) continue;
        const int row = step_ptr[u + static_cast<R_xlen_t>(n_steps) * (j - 1)];
        if (row == NA_INTEGER || row < 1 || row > n_gap) {
          cpp11::stop("step_rows needs a row in 1..nrow(x_gap) for every non-absorbing previous state in range.");
        }
      }
    }
  }

  const double* alpha_ptr = real_pointer(alpha);
  const double* beta_ptr = real_pointer(beta);
  const double* w_ptr = real_pointer(w_obs);
  const double* sw_ptr = real_pointer(stretch_weight);
  const double neg_inf = R_NegInf;

  double loglik = 0.0;
  bool impossible = false;
  std::vector<double> hessian(order >= 2 ? static_cast<size_t>(q) * q : 0, 0.0);
  std::vector<double> group_scores(
      order >= 1 ? static_cast<size_t>(n_groups) * q : 0, 0.0);

  // Observed transitions, closed form. With a_lo = s_{y-1}/P, a_hi = -s_y/P,
  // b_lo = u_{y-1}/P and b_hi = -u_y/P, the score of log P is
  // a_lo g_{y-1} + a_hi g_y and its Hessian is
  // (b_lo - a_lo^2) g_{y-1} g_{y-1}' + (b_hi - a_hi^2) g_y g_y'
  //   - a_lo a_hi (g_{y-1} g_y' + g_y g_{y-1}').
  // Because g_k = (e_k, x), the slope block of the Hessian is a weighted sum
  // of x x', accumulated over rows and expanded once at the end.
  {
    const double* x_ptr = real_pointer(SEXP(x_obs));
    CutPoints cuts(thresholds);
    std::vector<double> xx(order >= 2 ? static_cast<size_t>(p) * p : 0, 0.0);
    std::vector<double> alpha_x(order >= 2 ? static_cast<size_t>(thresholds) * p : 0, 0.0);
    std::vector<double> alpha_alpha(order >= 2 ? static_cast<size_t>(thresholds) * thresholds : 0, 0.0);
    std::vector<double> x_row(p);
    for (int i = 0; i < n_obs; ++i) {
      const double w = w_ptr[i];
      if (w == 0.0) {
        continue;
      }
      double eta = 0.0;
      for (int l = 0; l < p; ++l) {
        x_row[l] = x_ptr[i + static_cast<R_xlen_t>(n_obs) * l];
        eta += x_row[l] * beta_ptr[l];
      }
      cuts.set(alpha_ptr, thresholds, eta);
      const int y = y_ptr[i];
      const double prob = cell_probability(
          cuts.upper.data(), cuts.lower.data(), thresholds, y);
      if (!(prob > 0.0)) {
        impossible = true;
        continue;
      }
      loglik += w * std::log(prob);
      if (order < 1) {
        continue;
      }
      const double* up = cuts.upper.data();
      const double* low = cuts.lower.data();
      const int lo = y - 1;
      const int hi = y;
      const double a_lo = slope_at(up, low, thresholds, lo) / prob;
      const double a_hi = -slope_at(up, low, thresholds, hi) / prob;
      double* score = group_scores.data() + (group_ptr[i] - 1);
      const R_xlen_t stride = n_groups;
      if (lo >= 1) score[stride * (lo - 1)] += w * a_lo;
      if (hi <= thresholds) score[stride * (hi - 1)] += w * a_hi;
      const double slope_coef = w * (a_lo + a_hi);
      for (int l = 0; l < p; ++l) {
        score[stride * (thresholds + l)] += slope_coef * x_row[l];
      }
      if (order < 2) {
        continue;
      }
      const double b_lo = curvature_at(up, low, thresholds, lo) / prob;
      const double b_hi = -curvature_at(up, low, thresholds, hi) / prob;
      const double h_ll = w * (b_lo - a_lo * a_lo);
      const double h_hh = w * (b_hi - a_hi * a_hi);
      const double h_lh = -w * a_lo * a_hi;
      // g_0 and g_K do not exist (c_0, c_K are constants); their terms are 0
      // because a and b vanish there.
      if (lo >= 1) alpha_alpha[(lo - 1) + thresholds * (lo - 1)] += h_ll;
      if (hi <= thresholds) alpha_alpha[(hi - 1) + thresholds * (hi - 1)] += h_hh;
      if (lo >= 1 && hi <= thresholds) {
        alpha_alpha[(lo - 1) + thresholds * (hi - 1)] += h_lh;
        alpha_alpha[(hi - 1) + thresholds * (lo - 1)] += h_lh;
      }
      const double v_lo = h_ll + h_lh;
      const double v_hi = h_hh + h_lh;
      for (int l = 0; l < p; ++l) {
        if (lo >= 1) alpha_x[(lo - 1) + thresholds * l] += v_lo * x_row[l];
        if (hi <= thresholds) alpha_x[(hi - 1) + thresholds * l] += v_hi * x_row[l];
      }
      const double v_xx = h_ll + h_hh + 2.0 * h_lh;
      for (int l2 = 0; l2 < p; ++l2) {
        const double scaled = v_xx * x_row[l2];
        for (int l1 = 0; l1 <= l2; ++l1) {
          xx[l1 + static_cast<size_t>(p) * l2] += scaled * x_row[l1];
        }
      }
    }
    if (order >= 2) {
      for (int k2 = 0; k2 < thresholds; ++k2) {
        for (int k1 = 0; k1 < thresholds; ++k1) {
          hessian[k1 + static_cast<size_t>(q) * k2] += alpha_alpha[k1 + thresholds * k2];
        }
      }
      for (int l = 0; l < p; ++l) {
        for (int k = 0; k < thresholds; ++k) {
          const double value = alpha_x[k + thresholds * l];
          hessian[k + static_cast<size_t>(q) * (thresholds + l)] += value;
          hessian[(thresholds + l) + static_cast<size_t>(q) * k] += value;
        }
      }
      for (int l2 = 0; l2 < p; ++l2) {
        for (int l1 = 0; l1 <= l2; ++l1) {
          const double value = xx[l1 + static_cast<size_t>(p) * l2];
          hessian[(thresholds + l1) + static_cast<size_t>(q) * (thresholds + l2)] += value;
          if (l1 != l2) {
            hessian[(thresholds + l2) + static_cast<size_t>(q) * (thresholds + l1)] += value;
          }
        }
      }
    }
  }

  // Stretches: forward recursion over the states allowed at each step. The
  // forward vector a has one entry per state (index state - 1); only the
  // entries in the current step's range are non-zero, and only those are
  // touched. a_t(k) = sum_j a_{t-1}(j) P(k | j) for k in [lo_t, hi_t].
  writable::doubles stretch_loglik(m);
  std::vector<double> stretch_scores(order >= 1 ? static_cast<size_t>(m) * q : 0, 0.0);
  if (m > 0) {
    // Rows of x_gap are shared by many stretches and steps, so their
    // cut-points are computed once per call; the rows are also stored
    // transposed so that one row is contiguous.
    const double* gap_ptr = real_pointer(SEXP(x_gap));
    std::vector<double> gap_rows(static_cast<size_t>(n_gap) * p);
    std::vector<double> gap_upper(static_cast<size_t>(n_gap) * thresholds);
    std::vector<double> gap_lower(static_cast<size_t>(n_gap) * thresholds);
    {
      CutPoints cuts(thresholds);
      for (int r = 0; r < n_gap; ++r) {
        double eta = 0.0;
        for (int l = 0; l < p; ++l) {
          const double value = gap_ptr[r + static_cast<R_xlen_t>(n_gap) * l];
          gap_rows[static_cast<size_t>(r) * p + l] = value;
          eta += value * beta_ptr[l];
        }
        cuts.set(alpha_ptr, thresholds, eta);
        std::copy(cuts.upper.begin(), cuts.upper.end(),
                  gap_upper.begin() + static_cast<size_t>(r) * thresholds);
        std::copy(cuts.lower.begin(), cuts.lower.end(),
                  gap_lower.begin() + static_cast<size_t>(r) * thresholds);
      }
    }

    const size_t qq = static_cast<size_t>(q) * q;
    std::vector<double> a_in(states), a_out(states);
    std::vector<double> da_in(order >= 1 ? static_cast<size_t>(states) * q : 0);
    std::vector<double> da_out(da_in.size());
    std::vector<double> d2a_in(order >= 2 ? states * qq : 0);
    std::vector<double> d2a_out(d2a_in.size());
    std::vector<double> zero_q(q, 0.0), zero_qq(order >= 2 ? qq : 0, 0.0);
    std::vector<double> dp(q), g_lo(q), g_hi(q);
    std::vector<double> df(q), d2f(qq);
    bool negative_cell = false;

    // Adds the transitions out of one non-absorbing previous state, with
    // forward weight a_j (and its derivatives), through design row `row` to
    // the states out_lo..out_hi.
    auto add_transitions = [&](const int row, const double a_j,
                               const double* da_j, const double* d2a_j,
                               const int out_lo, const int out_hi) {
      const double* up = gap_upper.data() + static_cast<size_t>(row) * thresholds;
      const double* low = gap_lower.data() + static_cast<size_t>(row) * thresholds;
      const double* x = gap_rows.data() + static_cast<size_t>(row) * p;
      for (int y = out_lo; y <= out_hi; ++y) {
        const double prob = cell_probability(up, low, thresholds, y);
        if (prob < 0.0) negative_cell = true;
        a_out[y - 1] += a_j * prob;
        if (order < 1) continue;
        // dP = s_{y-1} g_{y-1} - s_y g_y
        const int lo = y - 1;
        const int hi = y;
        const double s_lo = slope_at(up, low, thresholds, lo);
        const double s_hi = slope_at(up, low, thresholds, hi);
        std::fill(g_lo.begin(), g_lo.end(), 0.0);
        std::fill(g_hi.begin(), g_hi.end(), 0.0);
        if (lo >= 1) g_lo[lo - 1] = 1.0;
        if (hi <= thresholds) g_hi[hi - 1] = 1.0;
        for (int l = 0; l < p; ++l) {
          g_lo[thresholds + l] = lo >= 1 ? x[l] : 0.0;
          g_hi[thresholds + l] = hi <= thresholds ? x[l] : 0.0;
        }
        for (int b = 0; b < q; ++b) dp[b] = s_lo * g_lo[b] - s_hi * g_hi[b];
        // Da'_y += Da_j P + a_j dP
        double* da_y = da_out.data() + static_cast<size_t>(y - 1) * q;
        for (int b = 0; b < q; ++b) da_y[b] += da_j[b] * prob + a_j * dp[b];
        if (order < 2) continue;
        // D2a'_y += D2a_j P + Da_j dP' + dP Da_j'
        //           + a_j (u_{y-1} g_{y-1} g_{y-1}' - u_y g_y g_y')
        const double c_lo = a_j * curvature_at(up, low, thresholds, lo);
        const double c_hi = -a_j * curvature_at(up, low, thresholds, hi);
        double* d2a_y = d2a_out.data() + (y - 1) * qq;
        for (int b2 = 0; b2 < q; ++b2) {
          const double da_b2 = da_j[b2];
          const double dp_b2 = dp[b2];
          const double glo_b2 = c_lo * g_lo[b2];
          const double ghi_b2 = c_hi * g_hi[b2];
          double* col = d2a_y + static_cast<size_t>(q) * b2;
          const double* col_in = d2a_j + static_cast<size_t>(q) * b2;
          // Only the upper triangle (b1 <= b2) is carried; the matrices
          // are symmetric and the stretch total is mirrored at the end.
          for (int b1 = 0; b1 <= b2; ++b1) {
            col[b1] += col_in[b1] * prob + da_j[b1] * dp_b2 + dp[b1] * da_b2 +
                       g_lo[b1] * glo_b2 + g_hi[b1] * ghi_b2;
          }
        }
      }
    };

    for (int s = 0; s < m; ++s) {
      negative_cell = false;
      const int first = first_step_ptr[s] - 1;
      const int steps = n_steps_ptr[s];
      int in_lo = 0;
      int in_hi = 0;
      for (int t = 0; t < steps; ++t) {
        const int u = first + t;
        const int out_lo = lo_ptr[u];
        const int out_hi = hi_ptr[u];
        std::fill(a_out.begin() + (out_lo - 1), a_out.begin() + out_hi, 0.0);
        if (order >= 1) {
          std::fill(da_out.begin() + static_cast<size_t>(out_lo - 1) * q,
                    da_out.begin() + static_cast<size_t>(out_hi) * q, 0.0);
        }
        if (order >= 2) {
          std::fill(d2a_out.begin() + (out_lo - 1) * qq,
                    d2a_out.begin() + out_hi * qq, 0.0);
        }
        if (t == 0) {
          // Step 1 starts from the known previous state: weight 1, no
          // derivatives.
          add_transitions(first_row_ptr[s] - 1, 1.0, zero_q.data(),
                          zero_qq.data(), out_lo, out_hi);
        } else {
          for (int j = in_lo; j <= in_hi; ++j) {
            const double a_j = a_in[j - 1];
            const double* da_j = order >= 1 ? da_in.data() + static_cast<size_t>(j - 1) * q : nullptr;
            const double* d2a_j = order >= 2 ? d2a_in.data() + (j - 1) * qq : nullptr;
            if (is_absorbing[j - 1]) {
              // Identity transition: the state stays j with probability 1.
              if (j < out_lo || j > out_hi) continue;
              a_out[j - 1] += a_j;
              if (order >= 1) {
                double* da_y = da_out.data() + static_cast<size_t>(j - 1) * q;
                for (int b = 0; b < q; ++b) da_y[b] += da_j[b];
              }
              if (order >= 2) {
                double* d2a_y = d2a_out.data() + (j - 1) * qq;
                for (size_t pos = 0; pos < qq; ++pos) d2a_y[pos] += d2a_j[pos];
              }
              continue;
            }
            const int row = step_ptr[u + static_cast<R_xlen_t>(n_steps) * (j - 1)] - 1;
            add_transitions(row, a_j, da_j, d2a_j, out_lo, out_hi);
          }
        }
        std::swap(a_in, a_out);
        std::swap(da_in, da_out);
        std::swap(d2a_in, d2a_out);
        in_lo = out_lo;
        in_hi = out_hi;
      }

      // f = sum of the final forward vector over the last step's range.
      double f = 0.0;
      for (int j = in_lo; j <= in_hi; ++j) f += a_in[j - 1];
      const double weight = sw_ptr[s];
      const bool possible = !negative_cell && f > 0.0;
      const double log_f = possible ? std::log(f) : neg_inf;
      stretch_loglik[s] = log_f;
      if (weight != 0.0) {
        if (possible) {
          loglik += weight * log_f;
        } else {
          impossible = true;
        }
      }
      if (order < 1) continue;
      std::fill(df.begin(), df.end(), 0.0);
      for (int j = in_lo; j <= in_hi; ++j) {
        for (int b = 0; b < q; ++b) df[b] += da_in[static_cast<size_t>(j - 1) * q + b];
      }
      // d log f = Df / f; d2 log f = D2f / f - (Df / f)(Df / f)'.
      for (int b = 0; b < q; ++b) {
        df[b] /= f;
        stretch_scores[s + static_cast<size_t>(m) * b] = df[b];
      }
      if (order < 2 || weight == 0.0) continue;
      std::fill(d2f.begin(), d2f.end(), 0.0);
      for (int j = in_lo; j <= in_hi; ++j) {
        const double* d2a_j = d2a_in.data() + (j - 1) * qq;
        for (size_t pos = 0; pos < qq; ++pos) d2f[pos] += d2a_j[pos];
      }
      for (int b2 = 0; b2 < q; ++b2) {
        for (int b1 = 0; b1 <= b2; ++b1) {
          const size_t pos = b1 + static_cast<size_t>(q) * b2;
          const double value = weight * (d2f[pos] / f - df[b1] * df[b2]);
          hessian[pos] += value;
          if (b1 != b2) hessian[b2 + static_cast<size_t>(q) * b1] += value;
        }
      }
    }
  }

  writable::list out;
  out.push_back(cpp11::named_arg("loglik") =
                    writable::doubles({impossible ? neg_inf : loglik}));
  out.push_back(cpp11::named_arg("stretch_loglik") = stretch_loglik);
  if (order >= 1) {
    out.push_back(cpp11::named_arg("obs_group_scores") =
                      to_matrix(group_scores, n_groups, q));
    out.push_back(cpp11::named_arg("stretch_scores") =
                      to_matrix(stretch_scores, m, q));
  }
  if (order >= 2) {
    out.push_back(cpp11::named_arg("hessian") = to_matrix(hessian, q, q));
  }
  return out;
}
