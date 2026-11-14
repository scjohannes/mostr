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
// theta = c(alpha, beta) are exact: they come from a forward and a backward
// pass over each stretch (see the stretch section below).
//
// Model: for a design row x and cut-point k = 1..K-1,
//   c_k = plogis(alpha_k + x'beta) = P(Y >= k + 1),
//   P(Y = y) = c_{y-1} - c_y with c_0 = 1 and c_K = 0.
// With g_k = (e_k, x), the derivatives are dc_k = s_k g_k and
// d2c_k = u_k g_k g_k', where s_k = c_k (1 - c_k) and u_k = s_k (1 - 2 c_k).

namespace {

// Cut-point quantities of one design row. Both c_k and 1 - c_k are kept
// to full relative precision, because a cell probability c_{y-1} - c_y
// computed from c alone loses all its digits when both are close to 1.
struct CutPoints {
  std::vector<double> upper;  // c_k
  std::vector<double> lower;  // 1 - c_k
  std::vector<double> exp_alpha;      // exp(alpha_k)
  std::vector<double> exp_neg_alpha;  // exp(-alpha_k)
  bool factor = true;

  CutPoints(const double* alpha, const int thresholds)
      : upper(thresholds), lower(thresholds), exp_alpha(thresholds),
        exp_neg_alpha(thresholds) {
    for (int k = 0; k < thresholds; ++k) {
      exp_alpha[k] = std::exp(alpha[k]);
      exp_neg_alpha[k] = std::exp(-alpha[k]);
      if (!(std::fabs(alpha[k]) < kFactorLimit)) factor = false;
    }
  }

  // exp is slow (about 40 ns per call in the Windows runtime) and caching
  // the x_gap rows was dominated by it, so exp(-|alpha_k + eta|) is formed
  // as exp(-+alpha_k) * exp(-+eta): one exp per row instead of K - 1. The
  // product is within a few ulps of the direct exp. Beyond kFactorLimit
  // the factors could overflow to Inf while the product is finite, so
  // large values take the direct exp.
  static constexpr double kFactorLimit = 300.0;

  void set(const double* alpha, const int thresholds, const double eta) {
    // c_k = plogis(z) and 1 - c_k = plogis(-z), each from the form that
    // cannot overflow.
    const bool use_factor = factor && std::fabs(eta) < kFactorLimit;
    const double exp_eta = use_factor ? std::exp(eta) : 0.0;
    const double exp_neg_eta = use_factor ? 1.0 / exp_eta : 0.0;
    for (int k = 0; k < thresholds; ++k) {
      const double z = alpha[k] + eta;
      double e;
      if (use_factor) {
        e = z >= 0.0 ? exp_neg_alpha[k] * exp_neg_eta : exp_alpha[k] * exp_eta;
      } else {
        e = std::exp(-std::fabs(z));
      }
      const double big = 1.0 / (1.0 + e);
      const double small = e / (1.0 + e);
      upper[k] = z >= 0.0 ? big : small;
      lower[k] = z >= 0.0 ? small : big;
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
    CutPoints cuts(alpha_ptr, thresholds);
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

  // Stretches. Write f = a_0' M_1 M_2 ... M_T 1, where M_t is the transition
  // matrix of step t restricted to the ranges (rows: previous states in the
  // previous range, identity rows for absorbing ones; columns: states in the
  // step's range) and a_0 picks the known previous state of step 1. Three
  // passes per stretch:
  // 1. forward: a_t = a_{t-1}' M_t, rescaled by c_t = sum(a_t) at every step
  //    so that long stretches cannot underflow; log f = sum_t log c_t.
  // 2. backward: b_T = 1 on the final range, b_{t-1} = M_t b_t / c_t, scaled
  //    so that sum(a_t * b_t) = 1 at every t.
  // 3. derivatives. With the scaled vectors,
  //      d log f = sum_t a_{t-1}' (dM_t) b_t / c_t,
  //    and, because f is multilinear in the M_t,
  //      D2 f / f = sum_t [a_{t-1}' (D2 M_t) b_t + (Da_{t-1})' (dM_t) b_t
  //                        + its transpose] / c_t,
  //    where Da_t = (Da_{t-1} M_t + a_{t-1}' dM_t) / c_t is the scaled
  //    gradient of the forward vector. Only first derivatives are carried
  //    through the recursion, and D2 log f = D2 f / f - (d log f)(d log f)'.
  // For a non-absorbing previous state j with design row x_j,
  //   sum_y dP(y | j) b(y) = sum_k s_jk (b(k+1) - b(k)) g_jk,
  //   sum_y D2P(y | j) b(y) = sum_k u_jk (b(k+1) - b(k)) g_jk g_jk',
  // with g_jk = (e_k, x_j) and b = 0 outside the step's range, because
  // P(y | j) = c_{y-1} - c_y.
  writable::doubles stretch_loglik(m);
  std::vector<double> stretch_scores(order >= 1 ? static_cast<size_t>(m) * q : 0, 0.0);
  if (m > 0) {
    // Rows of x_gap are shared by many stretches and steps, so their cell
    // probabilities, s_k and u_k are computed once per call; the rows are
    // also stored transposed so that one row is contiguous.
    const double* gap_ptr = real_pointer(SEXP(x_gap));
    std::vector<double> gap_rows(static_cast<size_t>(n_gap) * p);
    std::vector<double> gap_prob(static_cast<size_t>(n_gap) * states);
    std::vector<double> gap_s(static_cast<size_t>(n_gap) * thresholds);
    std::vector<double> gap_u(order >= 2 ? static_cast<size_t>(n_gap) * thresholds : 0);
    {
      CutPoints cuts(alpha_ptr, thresholds);
      for (int r = 0; r < n_gap; ++r) {
        double eta = 0.0;
        for (int l = 0; l < p; ++l) {
          const double value = gap_ptr[r + static_cast<R_xlen_t>(n_gap) * l];
          gap_rows[static_cast<size_t>(r) * p + l] = value;
          eta += value * beta_ptr[l];
        }
        cuts.set(alpha_ptr, thresholds, eta);
        const double* up = cuts.upper.data();
        const double* low = cuts.lower.data();
        for (int y = 1; y <= states; ++y) {
          gap_prob[static_cast<size_t>(r) * states + (y - 1)] =
              cell_probability(up, low, thresholds, y);
        }
        for (int k = 1; k <= thresholds; ++k) {
          gap_s[static_cast<size_t>(r) * thresholds + (k - 1)] =
              slope_at(up, low, thresholds, k);
          if (order >= 2) {
            gap_u[static_cast<size_t>(r) * thresholds + (k - 1)] =
                curvature_at(up, low, thresholds, k);
          }
        }
      }
    }

    // The curvature terms g_jk g_jk' depend on the stretch only through a
    // scalar weight per (x_gap row, k), so they are summed into gap_curv and
    // expanded into the Hessian once at the end.
    std::vector<double> gap_curv(order >= 2 ? static_cast<size_t>(n_gap) * thresholds : 0, 0.0);

    int longest = 1;
    for (int s = 0; s < m; ++s) longest = std::max(longest, n_steps_ptr[s]);
    // fwd[t * K + state - 1] holds the scaled a_t (t = 1..T); bwd the same
    // for b_t; scale[t] = c_t.
    std::vector<double> fwd(static_cast<size_t>(longest + 1) * states, 0.0);
    std::vector<double> bwd(static_cast<size_t>(longest + 1) * states, 0.0);
    std::vector<double> scale(longest + 1, 1.0);
    std::vector<double> diff(thresholds);
    std::vector<double> da_in(order >= 2 ? static_cast<size_t>(states) * q : 0);
    std::vector<double> da_out(da_in.size());
    std::vector<double> cross(order >= 2 ? static_cast<size_t>(q) * q : 0);
    std::vector<double> v(q);

    for (int s = 0; s < m; ++s) {
      const int first = first_step_ptr[s] - 1;
      const int steps = n_steps_ptr[s];
      const int first_gap_row = first_row_ptr[s] - 1;
      const double weight = sw_ptr[s];

      // Pass 1: forward values.
      bool negative_cell = false;
      bool possible = true;
      double log_f = 0.0;
      for (int t = 1; t <= steps && possible; ++t) {
        const int u = first + t - 1;
        const int out_lo = lo_ptr[u];
        const int out_hi = hi_ptr[u];
        double* a_t = fwd.data() + static_cast<size_t>(t) * states;
        std::fill(a_t + (out_lo - 1), a_t + out_hi, 0.0);
        if (t == 1) {
          const double* prob = gap_prob.data() + static_cast<size_t>(first_gap_row) * states;
          for (int y = out_lo; y <= out_hi; ++y) {
            if (prob[y - 1] < 0.0) negative_cell = true;
            a_t[y - 1] = prob[y - 1];
          }
        } else {
          const double* a_prev = fwd.data() + static_cast<size_t>(t - 1) * states;
          for (int j = lo_ptr[u - 1]; j <= hi_ptr[u - 1]; ++j) {
            const double a_j = a_prev[j - 1];
            if (is_absorbing[j - 1]) {
              if (j >= out_lo && j <= out_hi) a_t[j - 1] += a_j;
              continue;
            }
            const int row = step_ptr[u + static_cast<R_xlen_t>(n_steps) * (j - 1)] - 1;
            const double* prob = gap_prob.data() + static_cast<size_t>(row) * states;
            for (int y = out_lo; y <= out_hi; ++y) {
              if (prob[y - 1] < 0.0) negative_cell = true;
              a_t[y - 1] += a_j * prob[y - 1];
            }
          }
        }
        double total = 0.0;
        for (int y = out_lo; y <= out_hi; ++y) total += a_t[y - 1];
        if (!(total > 0.0) || negative_cell) {
          possible = false;
          break;
        }
        scale[t] = total;
        log_f += std::log(total);
        for (int y = out_lo; y <= out_hi; ++y) a_t[y - 1] /= total;
      }
      stretch_loglik[s] = possible ? log_f : neg_inf;
      if (weight != 0.0) {
        if (possible) {
          loglik += weight * log_f;
        } else {
          impossible = true;
        }
      }
      if (order < 1) continue;
      double* score = stretch_scores.data() + s;  // stride m
      const R_xlen_t score_stride = m;
      if (!possible) {
        // The derivatives of log f do not exist; a weighted impossible
        // stretch also makes the Hessian undefined.
        for (int b = 0; b < q; ++b) score[score_stride * b] = R_NaN;
        if (order >= 2 && weight != 0.0) {
          std::fill(hessian.begin(), hessian.end(), R_NaN);
        }
        continue;
      }

      // Pass 2: backward values.
      {
        const int u_last = first + steps - 1;
        double* b_last = bwd.data() + static_cast<size_t>(steps) * states;
        std::fill(b_last, b_last + states, 0.0);
        for (int y = lo_ptr[u_last]; y <= hi_ptr[u_last]; ++y) b_last[y - 1] = 1.0;
        for (int t = steps; t >= 2; --t) {
          const int u = first + t - 1;
          const double* b_t = bwd.data() + static_cast<size_t>(t) * states;
          double* b_prev = bwd.data() + static_cast<size_t>(t - 1) * states;
          std::fill(b_prev, b_prev + states, 0.0);
          const int out_lo = lo_ptr[u];
          const int out_hi = hi_ptr[u];
          for (int j = lo_ptr[u - 1]; j <= hi_ptr[u - 1]; ++j) {
            double value = 0.0;
            if (is_absorbing[j - 1]) {
              if (j >= out_lo && j <= out_hi) value = b_t[j - 1];
            } else {
              const int row = step_ptr[u + static_cast<R_xlen_t>(n_steps) * (j - 1)] - 1;
              const double* prob = gap_prob.data() + static_cast<size_t>(row) * states;
              for (int y = out_lo; y <= out_hi; ++y) value += prob[y - 1] * b_t[y - 1];
            }
            b_prev[j - 1] = value / scale[t];
          }
        }
      }

      // Pass 3: score by the adjoint form and, for order 2, the Hessian.
      for (int b = 0; b < q; ++b) score[score_stride * b] = 0.0;
      const bool hess = order >= 2 && weight != 0.0;
      if (hess) std::fill(cross.begin(), cross.end(), 0.0);
      for (int t = 1; t <= steps; ++t) {
        const int u = first + t - 1;
        const int out_lo = lo_ptr[u];
        const int out_hi = hi_ptr[u];
        const double* b_t = bwd.data() + static_cast<size_t>(t) * states;
        const double inv_scale = 1.0 / scale[t];
        // diff[k - 1] = b(k + 1) - b(k), with b = 0 outside the range.
        for (int k = 1; k <= thresholds; ++k) {
          const double b_k = (k >= out_lo && k <= out_hi) ? b_t[k - 1] : 0.0;
          const double b_k1 = (k + 1 >= out_lo && k + 1 <= out_hi) ? b_t[k] : 0.0;
          diff[k - 1] = b_k1 - b_k;
        }
        // Da_t is needed only as input to a later step.
        const bool propagate = hess && t < steps;
        if (propagate) {
          std::fill(da_out.begin() + static_cast<size_t>(out_lo - 1) * q,
                    da_out.begin() + static_cast<size_t>(out_hi) * q, 0.0);
        }
        const int in_lo = t == 1 ? 1 : lo_ptr[u - 1];
        const int in_hi = t == 1 ? 1 : hi_ptr[u - 1];
        const double* a_prev = fwd.data() + static_cast<size_t>(t - 1) * states;
        for (int j = in_lo; j <= in_hi; ++j) {
          // Step 1 has a single previous state with forward weight 1.
          const double a_j = t == 1 ? 1.0 : a_prev[j - 1];
          const double* da_j = (hess && t > 1) ? da_in.data() + static_cast<size_t>(j - 1) * q : nullptr;
          if (t > 1 && is_absorbing[j - 1]) {
            // Identity transition: no parameters, so only Da is carried on.
            if (propagate && j >= out_lo && j <= out_hi) {
              double* da_y = da_out.data() + static_cast<size_t>(j - 1) * q;
              for (int b = 0; b < q; ++b) da_y[b] += da_j[b] * inv_scale;
            }
            continue;
          }
          const int row = t == 1
              ? first_gap_row
              : step_ptr[u + static_cast<R_xlen_t>(n_steps) * (j - 1)] - 1;
          const double* x = gap_rows.data() + static_cast<size_t>(row) * p;
          const double* s_row = gap_s.data() + static_cast<size_t>(row) * thresholds;
          const double coef = a_j * inv_scale;
          // v = sum_y dP(y | j) b(y) = (m_1..m_{K-1}, (sum_k m_k) x_j) with
          // m_k = s_jk (b(k+1) - b(k)).
          double m_sum = 0.0;
          for (int k = 0; k < thresholds; ++k) {
            v[k] = s_row[k] * diff[k];
            m_sum += v[k];
          }
          for (int l = 0; l < p; ++l) v[thresholds + l] = m_sum * x[l];
          for (int b = 0; b < q; ++b) score[score_stride * b] += coef * v[b];
          if (!hess) continue;
          // Curvature: a_j u_jk (b(k+1) - b(k)) / c_t on g_jk g_jk'.
          const double* u_row = gap_u.data() + static_cast<size_t>(row) * thresholds;
          double* curv = gap_curv.data() + static_cast<size_t>(row) * thresholds;
          const double curv_coef = weight * coef;
          for (int k = 0; k < thresholds; ++k) curv[k] += curv_coef * u_row[k] * diff[k];
          // Cross term: Da_{t-1}(j) v' / c_t (symmetrised per stretch).
          if (t > 1) {
            for (int b2 = 0; b2 < q; ++b2) {
              const double scaled = v[b2] * inv_scale;
              if (scaled == 0.0) continue;
              double* col = cross.data() + static_cast<size_t>(q) * b2;
              for (int b1 = 0; b1 < q; ++b1) col[b1] += da_j[b1] * scaled;
            }
          }
          if (!propagate) continue;
          // Da_t(y) += (Da_{t-1}(j) P(y | j) + a_j dP(y | j)) / c_t, with
          // dP(y | j) = s_{y-1} g_{y-1} - s_y g_y.
          const double* prob = gap_prob.data() + static_cast<size_t>(row) * states;
          for (int y = out_lo; y <= out_hi; ++y) {
            double* da_y = da_out.data() + static_cast<size_t>(y - 1) * q;
            if (da_j != nullptr) {
              const double weight_prob = prob[y - 1] * inv_scale;
              for (int b = 0; b < q; ++b) da_y[b] += da_j[b] * weight_prob;
            }
            const double s_lo = y - 1 >= 1 ? s_row[y - 2] : 0.0;
            const double s_hi = y <= thresholds ? s_row[y - 1] : 0.0;
            if (y - 1 >= 1) da_y[y - 2] += coef * s_lo;
            if (y <= thresholds) da_y[y - 1] -= coef * s_hi;
            const double slope = coef * (s_lo - s_hi);
            for (int l = 0; l < p; ++l) da_y[thresholds + l] += slope * x[l];
          }
        }
        if (propagate) std::swap(da_in, da_out);
      }
      if (!hess) continue;
      // Hessian of log f: cross + cross' - score score' (curvature terms are
      // added at the end); upper triangle, mirrored below.
      for (int b2 = 0; b2 < q; ++b2) {
        const double score_b2 = score[score_stride * b2];
        for (int b1 = 0; b1 <= b2; ++b1) {
          const double value = weight * (cross[b1 + static_cast<size_t>(q) * b2] +
                                         cross[b2 + static_cast<size_t>(q) * b1] -
                                         score[score_stride * b1] * score_b2);
          hessian[b1 + static_cast<size_t>(q) * b2] += value;
          if (b1 != b2) hessian[b2 + static_cast<size_t>(q) * b1] += value;
        }
      }
    }

    if (order >= 2) {
      // Expand the curvature weights W_rk on g_rk g_rk' = (e_k, x_r)(e_k, x_r)':
      // W_rk on the intercept diagonal, W_rk x_r in the cross block and
      // (sum_k W_rk) x_r x_r' in the slope block.
      std::vector<double> xx(static_cast<size_t>(p) * p, 0.0);
      for (int r = 0; r < n_gap; ++r) {
        const double* curv = gap_curv.data() + static_cast<size_t>(r) * thresholds;
        const double* x = gap_rows.data() + static_cast<size_t>(r) * p;
        double total = 0.0;
        for (int k = 0; k < thresholds; ++k) {
          const double w_rk = curv[k];
          if (w_rk == 0.0) continue;
          total += w_rk;
          hessian[k + static_cast<size_t>(q) * k] += w_rk;
          for (int l = 0; l < p; ++l) {
            hessian[k + static_cast<size_t>(q) * (thresholds + l)] += w_rk * x[l];
            hessian[(thresholds + l) + static_cast<size_t>(q) * k] += w_rk * x[l];
          }
        }
        if (total == 0.0) continue;
        for (int l2 = 0; l2 < p; ++l2) {
          const double scaled = total * x[l2];
          for (int l1 = 0; l1 <= l2; ++l1) xx[l1 + static_cast<size_t>(p) * l2] += scaled * x[l1];
        }
      }
      for (int l2 = 0; l2 < p; ++l2) {
        for (int l1 = 0; l1 <= l2; ++l1) {
          const double value = xx[l1 + static_cast<size_t>(p) * l2];
          hessian[(thresholds + l1) + static_cast<size_t>(q) * (thresholds + l2)] += value;
          if (l1 != l2) hessian[(thresholds + l2) + static_cast<size_t>(q) * (thresholds + l1)] += value;
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
