// Log-likelihood of a proportional-odds Markov transition model when some
// visits have no exactly observed state, with its exact first and second
// derivatives.
//
// This header is the one implementation of the likelihood. It is used by the
// package's compiled code (src/paths.cpp, for engine = "orm") and by the Stan
// program of engine = "blrm" (through stan_path_loglik.hpp, compiled on the
// user's machine). It therefore depends on the C++ standard library only:
// no R, no Stan, and errors are thrown as std::invalid_argument. It may be
// included inside another namespace (rstan inserts it into the model's
// namespace) only after the standard headers below have been included.
//
// Observed transitions enter as ordinary cell probabilities. A run of visits
// whose states are known only to lie in a range of states (a "stretch"; an
// unobserved visit is the range of states it could be in) enters as the
// probability of all state paths consistent with the ranges, summed by a
// forward recursion over the states in each step's range. Derivatives with
// respect to theta = c(alpha, beta) come from a forward and a backward pass
// over each stretch (see evaluate_stretches()).
//
// Model: for a design row x and cut-point k = 1..K-1,
//   c_k = plogis(alpha_k + x'beta) = P(Y >= k + 1),
//   P(Y = y) = c_{y-1} - c_y with c_0 = 1 and c_K = 0.
// With g_k = (e_k, x), the derivatives are dc_k = s_k g_k and
// d2c_k = u_k g_k g_k', where s_k = c_k (1 - c_k) and u_k = s_k (1 - 2 c_k).

#ifndef MOSTR_PATH_LIKELIHOOD_HPP
#define MOSTR_PATH_LIKELIHOOD_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mostr {
namespace paths {

// The data of one likelihood. All indices are 1-based, as R passes them.
// Matrices are stored by row (row r of an n x p matrix at x + r * p).
//
// Design rows of the stretches. The design row for step u when the previous
// state is j is
//   x_gap[step_rows[u, j]] + x_prev[j],
// and for step 1 of stretch s, whose previous state first_prev[s] is known,
//   x_gap[first_row[s]] + x_prev[first_prev[s]].
// Two layouts use this:
// - general: x_gap holds every distinct design row and x_prev is null;
// - additive: when the previous state enters the model only through columns
//   that depend on nothing else, x_gap holds one row per step with those
//   columns set to 0, x_prev the K rows of those columns, and every
//   step_rows[u, j] of a step is the same row. The linear predictor is then
//   eta(u) + eta_prev(j), one dot product per step instead of one per
//   (step, previous state), and the score is accumulated per step.
struct Problem {
  int thresholds = 0;  // K - 1
  int p = 0;           // length of beta

  // Observed transitions.
  int n_obs = 0;
  const double* x_obs = nullptr;   // n_obs x p
  const int* y_obs = nullptr;      // in 1..K
  const double* w_obs = nullptr;   // weights; rows of weight 0 are skipped
  const int* group_obs = nullptr;  // in 1..n_groups; may be null when no
                                   // group scores are requested
  int n_groups = 0;

  // Stretches.
  int n_gap = 0;
  const double* x_gap = nullptr;   // n_gap x p
  const double* x_prev = nullptr;  // K x p, or null (general layout)
  const int* absorbing = nullptr;  // K flags: from an absorbing state the
                                   // transition is the identity
  int m = 0;                       // stretches
  const int* first_row = nullptr;  // row of x_gap for step 1
  const int* first_prev = nullptr; // previous state of step 1; needed only
                                   // with x_prev
  const int* first_step = nullptr; // stretch s covers steps first_step[s] ..
  const int* n_steps = nullptr;    //   first_step[s] + n_steps[s] - 1
  const double* stretch_weight = nullptr;
  int n_steps_total = 0;
  const int* lo = nullptr;         // at step u the state lies in lo..hi
  const int* hi = nullptr;
  const int* step_rows = nullptr;  // n_steps_total x K; row of x_gap for
                                   // step u and previous state j, 0 where
                                   // unused. It must be valid for every
                                   // non-absorbing j in the previous step's
                                   // range; the entry of a first step is
                                   // ignored.
};

// What evaluate() writes. Pointers that are null are not computed. Matrices
// are stored by column, as R stores them. Rows and stretches with weight 0
// do not enter loglik, the totals or the Hessian, so an impossible
// zero-weight row does not turn the total into NaN.
struct Output {
  int order = 0;                     // 0: loglik; 1: and scores; 2: and Hessian
  double loglik = 0.0;               // -Inf when a weighted term is impossible
  double* stretch_loglik = nullptr;  // m
  double* total_score = nullptr;     // q = K - 1 + p, needs order >= 1
  double* group_scores = nullptr;    // n_groups x q, needs order >= 1
  double* stretch_scores = nullptr;  // m x q, needs order >= 1
  double* hessian = nullptr;         // q x q, needs order 2
};

namespace detail {

// Cut-point quantities of one design row. Both c_k and 1 - c_k are kept
// to full relative precision, because a cell probability c_{y-1} - c_y
// computed from c alone loses all its digits when both are close to 1.
struct CutPoints {
  std::vector<double> upper;          // c_k
  std::vector<double> lower;          // 1 - c_k
  std::vector<double> exp_alpha;      // exp(alpha_k)
  std::vector<double> exp_neg_alpha;  // exp(-alpha_k)
  bool factor = true;

  // exp is slow (about 40 ns per call in the Windows runtime) and caching
  // the design rows was dominated by it, so exp(-|alpha_k + eta|) is formed
  // as exp(-+alpha_k) * exp(-+eta): one exp per row instead of K - 1. The
  // product is within a few ulps of the direct exp. Beyond kFactorLimit
  // the factors could overflow to Inf while the product is finite, so
  // large values take the direct exp.
  static constexpr double kFactorLimit = 300.0;

  CutPoints(const double* alpha, const int thresholds)
      : upper(thresholds), lower(thresholds), exp_alpha(thresholds),
        exp_neg_alpha(thresholds) {
    for (int k = 0; k < thresholds; ++k) {
      exp_alpha[k] = std::exp(alpha[k]);
      exp_neg_alpha[k] = std::exp(-alpha[k]);
      if (!(std::fabs(alpha[k]) < kFactorLimit)) factor = false;
    }
  }

  void set(const double* alpha, const int thresholds, const double eta) {
    const bool use_factor = factor && std::fabs(eta) < kFactorLimit;
    set(alpha, thresholds, eta, use_factor ? std::exp(eta) : 0.0, use_factor);
  }

  // exp_eta is exp(eta), used only when `factorable`: the caller may form it
  // as a product of cached factors (see evaluate_stretches()).
  void set(const double* alpha, const int thresholds, const double eta,
           const double exp_eta, const bool factorable) {
    // c_k = plogis(z) and 1 - c_k = plogis(-z), each from the form that
    // cannot overflow.
    const bool use_factor =
        factor && factorable && std::fabs(eta) < kFactorLimit;
    const double exp_neg_eta = use_factor ? 1.0 / exp_eta : 0.0;
    for (int k = 0; k < thresholds; ++k) {
      const double z = alpha[k] + eta;
      double e;
      if (use_factor) {
        e = z >= 0.0 ? exp_neg_alpha[k] * exp_neg_eta : exp_alpha[k] * exp_eta;
      } else {
        e = std::exp(-std::fabs(z));
      }
      // One division per cut-point: e / (1 + e) = e * big, equal to within
      // an ulp.
      const double big = 1.0 / (1.0 + e);
      const double small = e * big;
      upper[k] = z >= 0.0 ? big : small;
      lower[k] = z >= 0.0 ? small : big;
    }
  }
};

// P(Y = y) for 1-based y from the cut-point values c (upper) and 1 - c
// (lower) of one row. Whichever form subtracts the smaller numbers is used,
// so the difference keeps its relative precision in both tails.
inline double cell_probability(const double* upper, const double* lower,
                               const int thresholds, const int y) {
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
inline double slope_at(const double* upper, const double* lower,
                       const int thresholds, const int k) {
  if (k < 1 || k > thresholds) return 0.0;
  return upper[k - 1] * lower[k - 1];
}

inline double curvature_at(const double* upper, const double* lower,
                           const int thresholds, const int k) {
  if (k < 1 || k > thresholds) return 0.0;
  return upper[k - 1] * lower[k - 1] * (lower[k - 1] - upper[k - 1]);
}

inline double dot(const double* x, const double* beta, const int p) {
  double eta = 0.0;
  for (int l = 0; l < p; ++l) eta += x[l] * beta[l];
  return eta;
}

// Observed transitions, closed form. With a_lo = s_{y-1}/P, a_hi = -s_y/P,
// b_lo = u_{y-1}/P and b_hi = -u_y/P, the score of log P is
// a_lo g_{y-1} + a_hi g_y and its Hessian is
// (b_lo - a_lo^2) g_{y-1} g_{y-1}' + (b_hi - a_hi^2) g_y g_y'
//   - a_lo a_hi (g_{y-1} g_y' + g_y g_{y-1}').
// Because g_k = (e_k, x), the slope block of the Hessian is a weighted sum
// of x x', accumulated over rows and expanded once at the end.
inline void evaluate_observed(const Problem& d, const double* alpha,
                              const double* beta, Output& out,
                              double& loglik, bool& impossible) {
  const int thresholds = d.thresholds;
  const int p = d.p;
  const int q = thresholds + p;
  const int order = out.order;
  CutPoints cuts(alpha, thresholds);
  std::vector<double> xx(order >= 2 ? static_cast<size_t>(p) * p : 0, 0.0);
  std::vector<double> alpha_x(order >= 2 ? static_cast<size_t>(thresholds) * p : 0, 0.0);
  std::vector<double> alpha_alpha(order >= 2 ? static_cast<size_t>(thresholds) * thresholds : 0, 0.0);
  double* total = out.total_score;
  const size_t stride = static_cast<size_t>(d.n_groups);
  for (int i = 0; i < d.n_obs; ++i) {
    const double w = d.w_obs[i];
    if (w == 0.0) continue;
    const double* x_row = d.x_obs + static_cast<size_t>(i) * p;
    cuts.set(alpha, thresholds, dot(x_row, beta, p));
    const int y = d.y_obs[i];
    const double prob = cell_probability(cuts.upper.data(), cuts.lower.data(),
                                         thresholds, y);
    if (!(prob > 0.0)) {
      impossible = true;
      continue;
    }
    loglik += w * std::log(prob);
    if (order < 1) continue;
    const double* up = cuts.upper.data();
    const double* low = cuts.lower.data();
    const int lo = y - 1;
    const int hi = y;
    const double a_lo = slope_at(up, low, thresholds, lo) / prob;
    const double a_hi = -slope_at(up, low, thresholds, hi) / prob;
    const double slope_coef = w * (a_lo + a_hi);
    if (out.group_scores != nullptr) {
      double* score = out.group_scores + (d.group_obs[i] - 1);
      if (lo >= 1) score[stride * (lo - 1)] += w * a_lo;
      if (hi <= thresholds) score[stride * (hi - 1)] += w * a_hi;
      for (int l = 0; l < p; ++l) {
        score[stride * (thresholds + l)] += slope_coef * x_row[l];
      }
    }
    if (total != nullptr) {
      if (lo >= 1) total[lo - 1] += w * a_lo;
      if (hi <= thresholds) total[hi - 1] += w * a_hi;
      for (int l = 0; l < p; ++l) total[thresholds + l] += slope_coef * x_row[l];
    }
    if (order < 2) continue;
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
  if (order < 2) return;
  double* hessian = out.hessian;
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
//
// With x_j = r + x_prev[j] (r the x_gap row of the step; x_prev null in the
// general layout), the slope part of each of these sums is split into a
// part on r, summed over the previous states that share r before it is
// multiplied out ("runs" below; in the additive layout a whole step is one
// run), and a part on x_prev[j], summed per j over the whole stretch or
// call and multiplied out once.
inline void evaluate_stretches(const Problem& d, const double* alpha,
                               const double* beta, Output& out,
                               double& loglik, bool& impossible) {
  const int thresholds = d.thresholds;
  const int states = thresholds + 1;
  const int p = d.p;
  const int q = thresholds + p;
  const int m = d.m;
  const int order = out.order;
  const double neg_inf = -std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double* x_prev = d.x_prev;
  const bool has_prev = x_prev != nullptr;

  // Linear predictors of the x_gap rows and of x_prev, and their exp, once
  // per call. The design row of a (step, previous state) pair is used three
  // times in one stretch (forward, backward, derivatives) and, without tied
  // covariates, in no other stretch, so its cell probabilities are computed
  // per stretch into small tables that stay in cache: caching them per
  // design row for the whole call wrote and re-read about 17 MB for 100,000
  // rows. exp(eta(u) + eta_prev(j)) is the product of the cached factors,
  // so a pair needs no exp of its own.
  constexpr double limit = CutPoints::kFactorLimit;
  std::vector<double> eta_gap(d.n_gap);
  std::vector<double> exp_gap(d.n_gap);
  for (int r = 0; r < d.n_gap; ++r) {
    const double eta = dot(d.x_gap + static_cast<size_t>(r) * p, beta, p);
    eta_gap[r] = eta;
    exp_gap[r] = std::fabs(eta) < limit ? std::exp(eta) : 0.0;
  }
  std::vector<double> eta_prev(states, 0.0);
  std::vector<double> exp_prev(states, 1.0);
  if (has_prev) {
    for (int j = 0; j < states; ++j) {
      const double eta = dot(x_prev + static_cast<size_t>(j) * p, beta, p);
      eta_prev[j] = eta;
      exp_prev[j] = std::fabs(eta) < limit ? std::exp(eta) : 0.0;
    }
  }

  CutPoints cuts(alpha, thresholds);
  const bool hess_wanted = order >= 2;
  double* hessian = out.hessian;
  // The curvature terms g_jk g_jk' depend on the stretch only through a
  // scalar weight per (x_gap row, k) and per (previous state, k), so they
  // are summed into gap_curv and prev_curv and expanded into the Hessian
  // once at the end; row_prev holds, per x_gap row and previous state, the
  // summed weight of the mixed terms r x_prev[j]'.
  std::vector<double> gap_curv(hess_wanted ? static_cast<size_t>(d.n_gap) * thresholds : 0, 0.0);
  std::vector<double> prev_curv(hess_wanted && has_prev ? static_cast<size_t>(states) * thresholds : 0, 0.0);
  std::vector<double> row_prev(hess_wanted && has_prev ? static_cast<size_t>(d.n_gap) * states : 0, 0.0);
  // cross sums weight * Da_{t-1}(j) v' / c_t over all stretches; the
  // Hessian gets cross + cross' once at the end instead of once per
  // stretch. Its slope columns are built from run_vec (per run, times r)
  // and prev_vec (per previous state, times x_prev[j], once per call).
  std::vector<double> cross(hess_wanted ? static_cast<size_t>(q) * q : 0, 0.0);
  std::vector<double> prev_vec(hess_wanted && has_prev ? static_cast<size_t>(states) * q : 0, 0.0);
  std::vector<double> run_vec(hess_wanted ? q : 0, 0.0);

  int longest = 1;
  for (int s = 0; s < m; ++s) longest = std::max(longest, d.n_steps[s]);
  // Per-stretch tables, indexed by step t = 1..T and previous state j:
  // the x_gap row (0-based), the cell probabilities, s_k and u_k.
  const size_t pairs = static_cast<size_t>(longest + 1) * states;
  std::vector<int> row_tab(pairs, -1);
  std::vector<double> prob_tab(pairs * states);
  std::vector<double> s_tab(pairs * thresholds);
  std::vector<double> u_tab(hess_wanted ? pairs * thresholds : 0);
  // fwd[t * K + state - 1] holds the scaled a_t (t = 1..T); bwd the same
  // for b_t; scale[t] = c_t.
  std::vector<double> fwd(static_cast<size_t>(longest + 1) * states, 0.0);
  std::vector<double> bwd(static_cast<size_t>(longest + 1) * states, 0.0);
  std::vector<double> scale(longest + 1, 1.0);
  std::vector<double> diff(thresholds);
  std::vector<double> da_in(hess_wanted ? static_cast<size_t>(states) * q : 0);
  std::vector<double> da_out(da_in.size());
  std::vector<double> v(thresholds);
  std::vector<double> x_full(p);
  // The stretch's score is summed in a contiguous vector and copied out
  // once; prev_coef holds its part on x_prev[j].
  std::vector<double> score_local(q);
  std::vector<double> prev_coef(states, 0.0);

  for (int s = 0; s < m; ++s) {
    const int first = d.first_step[s] - 1;
    const int steps = d.n_steps[s];
    const double weight = d.stretch_weight[s];
    // The previous state of step 1 matters only through x_prev.
    const int first_prev = has_prev ? d.first_prev[s] : 1;
    auto pair = [&](int t, int j) {
      return static_cast<size_t>(t) * states + (j - 1);
    };
    auto in_range = [&](int t, int& in_lo, int& in_hi) {
      const int u = first + t - 1;
      in_lo = t == 1 ? first_prev : d.lo[u - 1];
      in_hi = t == 1 ? first_prev : d.hi[u - 1];
    };

    // Cell probabilities of every (step, previous state) pair.
    for (int t = 1; t <= steps; ++t) {
      const int u = first + t - 1;
      int in_lo, in_hi;
      in_range(t, in_lo, in_hi);
      for (int j = in_lo; j <= in_hi; ++j) {
        if (t > 1 && d.absorbing[j - 1]) continue;
        const int row = t == 1
            ? d.first_row[s] - 1
            : d.step_rows[static_cast<size_t>(u) * states + (j - 1)] - 1;
        const size_t at = pair(t, j);
        row_tab[at] = row;
        const double eta = eta_gap[row] + eta_prev[j - 1];
        const bool factorable = std::fabs(eta_gap[row]) < limit &&
                                std::fabs(eta_prev[j - 1]) < limit;
        cuts.set(alpha, thresholds, eta, exp_gap[row] * exp_prev[j - 1],
                 factorable);
        const double* up = cuts.upper.data();
        const double* low = cuts.lower.data();
        double* prob = prob_tab.data() + at * states;
        for (int y = 1; y <= states; ++y) {
          prob[y - 1] = cell_probability(up, low, thresholds, y);
        }
        double* s_row = s_tab.data() + at * thresholds;
        for (int k = 1; k <= thresholds; ++k) {
          s_row[k - 1] = slope_at(up, low, thresholds, k);
        }
        if (hess_wanted) {
          double* u_row = u_tab.data() + at * thresholds;
          for (int k = 1; k <= thresholds; ++k) {
            u_row[k - 1] = curvature_at(up, low, thresholds, k);
          }
        }
      }
    }

    // Pass 1: forward values.
    bool negative_cell = false;
    bool possible = true;
    double log_f = 0.0;
    for (int t = 1; t <= steps && possible; ++t) {
      const int u = first + t - 1;
      const int out_lo = d.lo[u];
      const int out_hi = d.hi[u];
      double* a_t = fwd.data() + static_cast<size_t>(t) * states;
      std::fill(a_t + (out_lo - 1), a_t + out_hi, 0.0);
      if (t == 1) {
        const double* prob = prob_tab.data() + pair(1, first_prev) * states;
        for (int y = out_lo; y <= out_hi; ++y) {
          if (prob[y - 1] < 0.0) negative_cell = true;
          a_t[y - 1] = prob[y - 1];
        }
      } else {
        const double* a_prev = fwd.data() + static_cast<size_t>(t - 1) * states;
        for (int j = d.lo[u - 1]; j <= d.hi[u - 1]; ++j) {
          const double a_j = a_prev[j - 1];
          if (d.absorbing[j - 1]) {
            if (j >= out_lo && j <= out_hi) a_t[j - 1] += a_j;
            continue;
          }
          const double* prob = prob_tab.data() + pair(t, j) * states;
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
    if (out.stretch_loglik != nullptr) {
      out.stretch_loglik[s] = possible ? log_f : neg_inf;
    }
    if (weight != 0.0) {
      if (possible) {
        loglik += weight * log_f;
      } else {
        impossible = true;
      }
    }
    if (order < 1) continue;
    if (!possible) {
      // The derivatives of log f do not exist; a weighted impossible
      // stretch also makes the total score and the Hessian undefined.
      if (out.stretch_scores != nullptr) {
        for (int b = 0; b < q; ++b) out.stretch_scores[s + static_cast<size_t>(m) * b] = nan;
      }
      if (weight != 0.0) {
        if (out.total_score != nullptr) {
          std::fill(out.total_score, out.total_score + q, nan);
        }
        if (hess_wanted) std::fill(hessian, hessian + static_cast<size_t>(q) * q, nan);
      }
      continue;
    }

    // Pass 2: backward values.
    {
      const int u_last = first + steps - 1;
      double* b_last = bwd.data() + static_cast<size_t>(steps) * states;
      std::fill(b_last, b_last + states, 0.0);
      for (int y = d.lo[u_last]; y <= d.hi[u_last]; ++y) b_last[y - 1] = 1.0;
      for (int t = steps; t >= 2; --t) {
        const int u = first + t - 1;
        const double* b_t = bwd.data() + static_cast<size_t>(t) * states;
        double* b_prev = bwd.data() + static_cast<size_t>(t - 1) * states;
        std::fill(b_prev, b_prev + states, 0.0);
        const int out_lo = d.lo[u];
        const int out_hi = d.hi[u];
        for (int j = d.lo[u - 1]; j <= d.hi[u - 1]; ++j) {
          double value = 0.0;
          if (d.absorbing[j - 1]) {
            if (j >= out_lo && j <= out_hi) value = b_t[j - 1];
          } else {
            const double* prob = prob_tab.data() + pair(t, j) * states;
            for (int y = out_lo; y <= out_hi; ++y) value += prob[y - 1] * b_t[y - 1];
          }
          b_prev[j - 1] = value / scale[t];
        }
      }
    }

    // Pass 3: score by the adjoint form and, for order 2, the Hessian.
    std::fill(score_local.begin(), score_local.end(), 0.0);
    if (has_prev) std::fill(prev_coef.begin(), prev_coef.end(), 0.0);
    const bool hess = hess_wanted && weight != 0.0;
    // A run: consecutive previous states whose design rows share the x_gap
    // row run_row. Their slope terms are summed in run_coef (score) and
    // run_vec (cross term) and multiplied by the row once.
    int run_row = -1;
    double run_coef = 0.0;
    auto flush = [&]() {
      if (run_row < 0) return;
      const double* r = d.x_gap + static_cast<size_t>(run_row) * p;
      for (int l = 0; l < p; ++l) score_local[thresholds + l] += run_coef * r[l];
      if (hess) {
        for (int l = 0; l < p; ++l) {
          const double r_l = r[l];
          if (r_l == 0.0) continue;
          double* col = cross.data() + static_cast<size_t>(q) * (thresholds + l);
          for (int b = 0; b < q; ++b) col[b] += run_vec[b] * r_l;
        }
        std::fill(run_vec.begin(), run_vec.end(), 0.0);
      }
      run_row = -1;
      run_coef = 0.0;
    };
    for (int t = 1; t <= steps; ++t) {
      const int u = first + t - 1;
      const int out_lo = d.lo[u];
      const int out_hi = d.hi[u];
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
      int in_lo, in_hi;
      in_range(t, in_lo, in_hi);
      const double* a_prev = fwd.data() + static_cast<size_t>(t - 1) * states;
      for (int j = in_lo; j <= in_hi; ++j) {
        // Step 1 has a single previous state with forward weight 1.
        const double a_j = t == 1 ? 1.0 : a_prev[j - 1];
        const double* da_j = (hess && t > 1) ? da_in.data() + static_cast<size_t>(j - 1) * q : nullptr;
        if (t > 1 && d.absorbing[j - 1]) {
          // Identity transition: no parameters, so only Da is carried on.
          if (propagate && j >= out_lo && j <= out_hi) {
            double* da_y = da_out.data() + static_cast<size_t>(j - 1) * q;
            for (int b = 0; b < q; ++b) da_y[b] += da_j[b] * inv_scale;
          }
          continue;
        }
        const size_t at = pair(t, j);
        const int row = row_tab[at];
        const double* s_row = s_tab.data() + at * thresholds;
        const double coef = a_j * inv_scale;
        // sum_y dP(y | j) b(y) = (m_1..m_{K-1}, (sum_k m_k) x_j) with
        // m_k = s_jk (b(k+1) - b(k)).
        double m_sum = 0.0;
        for (int k = 0; k < thresholds; ++k) {
          v[k] = s_row[k] * diff[k];
          m_sum += v[k];
          score_local[k] += coef * v[k];
        }
        if (row != run_row) {
          flush();
          run_row = row;
        }
        const double slope_coef = coef * m_sum;
        run_coef += slope_coef;
        if (has_prev) prev_coef[j - 1] += slope_coef;
        if (!hess) continue;
        // Curvature: a_j u_jk (b(k+1) - b(k)) / c_t on g_jk g_jk'.
        const double* u_row = u_tab.data() + at * thresholds;
        double* curv = gap_curv.data() + static_cast<size_t>(row) * thresholds;
        const double curv_coef = weight * coef;
        double curv_total = 0.0;
        for (int k = 0; k < thresholds; ++k) {
          const double value = curv_coef * u_row[k] * diff[k];
          curv[k] += value;
          if (has_prev) {
            prev_curv[static_cast<size_t>(j - 1) * thresholds + k] += value;
            curv_total += value;
          }
        }
        if (has_prev) row_prev[static_cast<size_t>(row) * states + (j - 1)] += curv_total;
        // Cross term: weight * Da_{t-1}(j) v' / c_t (symmetrised at the
        // end); its slope columns go through run_vec and prev_vec.
        if (t > 1) {
          for (int k = 0; k < thresholds; ++k) {
            const double scaled = weight * v[k] * inv_scale;
            if (scaled == 0.0) continue;
            double* col = cross.data() + static_cast<size_t>(q) * k;
            for (int b = 0; b < q; ++b) col[b] += da_j[b] * scaled;
          }
          const double scaled = weight * m_sum * inv_scale;
          if (scaled != 0.0) {
            for (int b = 0; b < q; ++b) run_vec[b] += da_j[b] * scaled;
            if (has_prev) {
              double* pv = prev_vec.data() + static_cast<size_t>(j - 1) * q;
              for (int b = 0; b < q; ++b) pv[b] += da_j[b] * scaled;
            }
          }
        }
        if (!propagate) continue;
        // Da_t(y) += (Da_{t-1}(j) P(y | j) + a_j dP(y | j)) / c_t, with
        // dP(y | j) = s_{y-1} g_{y-1} - s_y g_y.
        const double* x = d.x_gap + static_cast<size_t>(row) * p;
        if (has_prev) {
          const double* xp = x_prev + static_cast<size_t>(j - 1) * p;
          for (int l = 0; l < p; ++l) x_full[l] = x[l] + xp[l];
          x = x_full.data();
        }
        const double* prob = prob_tab.data() + at * states;
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
      flush();
      if (propagate) std::swap(da_in, da_out);
    }
    if (has_prev) {
      for (int j = 0; j < states; ++j) {
        const double c = prev_coef[j];
        if (c == 0.0) continue;
        const double* xp = x_prev + static_cast<size_t>(j) * p;
        for (int l = 0; l < p; ++l) score_local[thresholds + l] += c * xp[l];
      }
    }
    if (out.stretch_scores != nullptr) {
      for (int b = 0; b < q; ++b) {
        out.stretch_scores[s + static_cast<size_t>(m) * b] = score_local[b];
      }
    }
    if (out.total_score != nullptr && weight != 0.0) {
      for (int b = 0; b < q; ++b) out.total_score[b] += weight * score_local[b];
    }
    if (!hess) continue;
    // Hessian of log f: cross + cross' - score score' (cross and the
    // curvature terms are added at the end); upper triangle, mirrored so
    // that the Hessian stays exactly symmetric.
    for (int b2 = 0; b2 < q; ++b2) {
      const double score_b2 = weight * score_local[b2];
      for (int b1 = 0; b1 <= b2; ++b1) {
        const double value = score_local[b1] * score_b2;
        hessian[b1 + static_cast<size_t>(q) * b2] -= value;
        if (b1 != b2) hessian[b2 + static_cast<size_t>(q) * b1] -= value;
      }
    }
  }

  if (!hess_wanted) return;
  // The slope columns of cross on x_prev[j].
  if (has_prev) {
    for (int j = 0; j < states; ++j) {
      const double* pv = prev_vec.data() + static_cast<size_t>(j) * q;
      const double* xp = x_prev + static_cast<size_t>(j) * p;
      for (int l = 0; l < p; ++l) {
        if (xp[l] == 0.0) continue;
        double* col = cross.data() + static_cast<size_t>(q) * (thresholds + l);
        for (int b = 0; b < q; ++b) col[b] += pv[b] * xp[l];
      }
    }
  }
  for (int b2 = 0; b2 < q; ++b2) {
    for (int b1 = 0; b1 < q; ++b1) {
      hessian[b1 + static_cast<size_t>(q) * b2] +=
          cross[b1 + static_cast<size_t>(q) * b2] +
          cross[b2 + static_cast<size_t>(q) * b1];
    }
  }
  // Expand the curvature weights W on g g' = (e_k, x)(e_k, x)' with
  // x = r + x_prev[j]: W on the intercept diagonal (once per pair, from
  // gap_curv), W x in the cross block and (sum_k W) x x' in the slope block,
  // where x x' = r r' + x_prev[j] x_prev[j]' + r x_prev[j]' + x_prev[j] r'.
  std::vector<double> xx(static_cast<size_t>(p) * p, 0.0);
  // ax[k * p + l]: the intercept-by-slope block, summed in one contiguous
  // row per cut-point and written to both halves of the Hessian once.
  std::vector<double> ax(static_cast<size_t>(thresholds) * p, 0.0);
  std::vector<double> mixed(has_prev ? p : 0);
  for (int r = 0; r < d.n_gap; ++r) {
    const double* curv = gap_curv.data() + static_cast<size_t>(r) * thresholds;
    const double* x = d.x_gap + static_cast<size_t>(r) * p;
    double total = 0.0;
    for (int k = 0; k < thresholds; ++k) {
      const double w_rk = curv[k];
      if (w_rk == 0.0) continue;
      total += w_rk;
      hessian[k + static_cast<size_t>(q) * k] += w_rk;
      double* ax_k = ax.data() + static_cast<size_t>(k) * p;
      for (int l = 0; l < p; ++l) ax_k[l] += w_rk * x[l];
    }
    if (total == 0.0) continue;
    for (int l2 = 0; l2 < p; ++l2) {
      const double scaled = total * x[l2];
      for (int l1 = 0; l1 <= l2; ++l1) xx[l1 + static_cast<size_t>(p) * l2] += scaled * x[l1];
    }
    if (!has_prev) continue;
    // mixed = sum_j W_rj x_prev[j]; adds r mixed' + mixed r'.
    std::fill(mixed.begin(), mixed.end(), 0.0);
    bool any = false;
    for (int j = 0; j < states; ++j) {
      const double w_rj = row_prev[static_cast<size_t>(r) * states + j];
      if (w_rj == 0.0) continue;
      any = true;
      const double* xp = x_prev + static_cast<size_t>(j) * p;
      for (int l = 0; l < p; ++l) mixed[l] += w_rj * xp[l];
    }
    if (!any) continue;
    for (int l2 = 0; l2 < p; ++l2) {
      for (int l1 = 0; l1 <= l2; ++l1) {
        xx[l1 + static_cast<size_t>(p) * l2] += x[l1] * mixed[l2] + mixed[l1] * x[l2];
      }
    }
  }
  if (has_prev) {
    for (int j = 0; j < states; ++j) {
      const double* curv = prev_curv.data() + static_cast<size_t>(j) * thresholds;
      const double* xp = x_prev + static_cast<size_t>(j) * p;
      double total = 0.0;
      for (int k = 0; k < thresholds; ++k) {
        const double w_jk = curv[k];
        if (w_jk == 0.0) continue;
        total += w_jk;
        double* ax_k = ax.data() + static_cast<size_t>(k) * p;
        for (int l = 0; l < p; ++l) ax_k[l] += w_jk * xp[l];
      }
      if (total == 0.0) continue;
      for (int l2 = 0; l2 < p; ++l2) {
        const double scaled = total * xp[l2];
        for (int l1 = 0; l1 <= l2; ++l1) xx[l1 + static_cast<size_t>(p) * l2] += scaled * xp[l1];
      }
    }
  }
  for (int k = 0; k < thresholds; ++k) {
    for (int l = 0; l < p; ++l) {
      const double value = ax[static_cast<size_t>(k) * p + l];
      hessian[k + static_cast<size_t>(q) * (thresholds + l)] += value;
      hessian[(thresholds + l) + static_cast<size_t>(q) * k] += value;
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

}  // namespace detail

// Check every index once, so that evaluate() can trust them. Throws
// std::invalid_argument naming the offending argument.
inline void validate(const Problem& d, const bool need_groups) {
  const int states = d.thresholds + 1;
  if (d.thresholds < 1) throw std::invalid_argument("alpha must have length >= 1.");
  if (d.p < 0 || d.n_obs < 0 || d.n_gap < 0 || d.m < 0 || d.n_steps_total < 0) {
    throw std::invalid_argument("Sizes must be non-negative.");
  }
  if (d.n_groups < 0) throw std::invalid_argument("n_groups must be non-negative.");
  for (int i = 0; i < d.n_obs; ++i) {
    if (d.y_obs[i] < 1 || d.y_obs[i] > states) {
      throw std::invalid_argument("y_obs must lie in 1..K.");
    }
    if (need_groups && (d.group_obs[i] < 1 || d.group_obs[i] > d.n_groups)) {
      throw std::invalid_argument("group_obs must lie in 1..n_groups.");
    }
  }
  for (int u = 0; u < d.n_steps_total; ++u) {
    if (d.lo[u] < 1 || d.hi[u] < d.lo[u] || d.hi[u] > states) {
      throw std::invalid_argument("step_lo and step_hi must satisfy 1 <= lo <= hi <= K.");
    }
  }
  for (int s = 0; s < d.m; ++s) {
    if (d.first_row[s] < 1 || d.first_row[s] > d.n_gap) {
      throw std::invalid_argument("stretch_first_row must lie in 1..nrow(x_gap).");
    }
    if (d.x_prev != nullptr && (d.first_prev[s] < 1 || d.first_prev[s] > states)) {
      throw std::invalid_argument("stretch_first_prev must lie in 1..K.");
    }
    if (d.n_steps[s] < 1 || d.first_step[s] < 1 ||
        d.first_step[s] > d.n_steps_total - d.n_steps[s] + 1) {
      throw std::invalid_argument("stretch_first_step and stretch_n_steps point outside the step arrays.");
    }
    for (int t = 1; t < d.n_steps[s]; ++t) {
      const int u = d.first_step[s] - 1 + t;
      for (int j = d.lo[u - 1]; j <= d.hi[u - 1]; ++j) {
        if (d.absorbing[j - 1]) continue;
        const int row = d.step_rows[static_cast<size_t>(u) * states + (j - 1)];
        if (row < 1 || row > d.n_gap) {
          throw std::invalid_argument("step_rows needs a row in 1..nrow(x_gap) for every non-absorbing previous state in range.");
        }
      }
    }
  }
}

// The log-likelihood and, by out.order, its derivatives at
// theta = c(alpha, beta). Every requested output is overwritten.
inline void evaluate(const Problem& d, const double* alpha, const double* beta,
                     Output& out) {
  if (out.order < 0 || out.order > 2) {
    throw std::invalid_argument("order must be 0, 1 or 2.");
  }
  const int q = d.thresholds + d.p;
  if (out.order < 1) {
    out.total_score = out.group_scores = out.stretch_scores = nullptr;
  }
  if (out.order < 2) out.hessian = nullptr;
  if (out.total_score != nullptr) std::fill(out.total_score, out.total_score + q, 0.0);
  if (out.group_scores != nullptr) {
    std::fill(out.group_scores, out.group_scores + static_cast<size_t>(d.n_groups) * q, 0.0);
  }
  if (out.stretch_scores != nullptr) {
    std::fill(out.stretch_scores, out.stretch_scores + static_cast<size_t>(d.m) * q, 0.0);
  }
  if (out.hessian != nullptr) {
    std::fill(out.hessian, out.hessian + static_cast<size_t>(q) * q, 0.0);
  } else if (out.order >= 2) {
    throw std::invalid_argument("order 2 needs a Hessian to write to.");
  }
  double loglik = 0.0;
  bool impossible = false;
  if (d.n_obs > 0) detail::evaluate_observed(d, alpha, beta, out, loglik, impossible);
  if (d.m > 0) detail::evaluate_stretches(d, alpha, beta, out, loglik, impossible);
  out.loglik = impossible ? -std::numeric_limits<double>::infinity() : loglik;
}

}  // namespace paths
}  // namespace mostr

#endif  // MOSTR_PATH_LIKELIHOOD_HPP
