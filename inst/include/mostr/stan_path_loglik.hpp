// The likelihood of path_likelihood.hpp for the Stan program
// inst/stan/po_paths.stan, with the gradient computed analytically.
//
// Stan differentiated the recursion over the states of each stretch by
// reverse-mode automatic differentiation, which records every product and
// sum of the recursion and then walks them backwards: about 120 ms per
// gradient for 7,440 stretches of 21,362 unobserved visits without tied
// covariates. path_likelihood.hpp computes the same gradient by a forward and
// a backward pass over each stretch, and Stan receives it as precomputed
// partial derivatives of one value.
//
// The functions here take the Stan program's data as flat arrays (matrices
// stored by row) and return the log-likelihood and, when asked, its
// gradient with respect to c(alpha, beta, gamma). The Stan-visible overloads are in
// stan_path_loglik_functions.hpp, which must be included inside the model's
// namespace. This header needs Stan Math (stan/math/rev.hpp) to have been
// included, as it is in every compiled Stan model.

#ifndef MOSTR_STAN_PATH_LOGLIK_HPP
#define MOSTR_STAN_PATH_LOGLIK_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

#include "path_likelihood.hpp"

namespace mostr {
namespace stan_paths {

// The log-likelihood at (alpha, beta, gamma); with `gradient` not null, also
// its gradient with respect to c(alpha, beta, gamma) (a vector of that
// length). gamma holds the random intercepts of the clusters, or is empty
// for a model without them; obs_cluster and stretch_cluster then are empty
// too, else they give the cluster of each observed row and each stretch. The
// sizes of the arrays define the problem: K - 1 = size(alpha),
// p = size(beta), observed rows = size(y_obs), stretches = size(first_row),
// steps = size(lo), clusters = size(gamma). Throws std::invalid_argument for
// inconsistent data.
inline double path_loglik(
    const std::vector<double>& alpha, const std::vector<double>& beta,
    const std::vector<double>& gamma,
    const std::vector<double>& x_obs, const std::vector<int>& y_obs,
    const std::vector<double>& w_obs, const std::vector<double>& x_gap,
    const std::vector<double>& x_prev, const std::vector<int>& absorbing,
    const std::vector<int>& first_row, const std::vector<int>& first_prev,
    const std::vector<int>& first_step, const std::vector<int>& n_steps,
    const std::vector<double>& weight, const std::vector<int>& lo,
    const std::vector<int>& hi, const std::vector<int>& step_rows,
    const std::vector<int>& obs_cluster,
    const std::vector<int>& stretch_cluster,
    std::vector<double>* gradient) {
  const int thresholds = static_cast<int>(alpha.size());
  const int states = thresholds + 1;
  const int p = static_cast<int>(beta.size());
  const int n_clusters = static_cast<int>(gamma.size());
  const size_t n_obs = y_obs.size();
  const size_t m = first_row.size();
  const size_t n_steps_total = lo.size();
  if (p < 1) throw std::invalid_argument("path_loglik: beta must not be empty.");
  if (n_clusters > 0 &&
      (obs_cluster.size() != n_obs || stretch_cluster.size() != m)) {
    throw std::invalid_argument(
        "path_loglik: obs_cluster and stretch_cluster must match y_obs and first_row.");
  }
  if (x_obs.size() != n_obs * p || w_obs.size() != n_obs) {
    throw std::invalid_argument("path_loglik: x_obs and w_obs do not match y_obs.");
  }
  if (x_gap.size() % p != 0) {
    throw std::invalid_argument("path_loglik: x_gap must hold whole rows.");
  }
  if (!x_prev.empty() && x_prev.size() != static_cast<size_t>(states) * p) {
    throw std::invalid_argument("path_loglik: x_prev must be empty or K x p.");
  }
  if (absorbing.size() != static_cast<size_t>(states)) {
    throw std::invalid_argument("path_loglik: absorbing must have K entries.");
  }
  if (first_prev.size() != m || first_step.size() != m || n_steps.size() != m ||
      weight.size() != m) {
    throw std::invalid_argument("path_loglik: the stretch arrays differ in length.");
  }
  if (hi.size() != n_steps_total ||
      step_rows.size() != n_steps_total * static_cast<size_t>(states)) {
    throw std::invalid_argument("path_loglik: the step arrays differ in length.");
  }

  mostr::paths::Problem d;
  d.thresholds = thresholds;
  d.p = p;
  d.n_obs = static_cast<int>(n_obs);
  d.x_obs = x_obs.data();
  d.y_obs = y_obs.data();
  d.w_obs = w_obs.data();
  d.n_gap = static_cast<int>(x_gap.size() / p);
  d.x_gap = x_gap.data();
  d.x_prev = x_prev.empty() ? nullptr : x_prev.data();
  d.absorbing = absorbing.data();
  d.m = static_cast<int>(m);
  d.first_row = first_row.data();
  d.first_prev = first_prev.data();
  d.first_step = first_step.data();
  d.n_steps = n_steps.data();
  d.stretch_weight = weight.data();
  d.n_steps_total = static_cast<int>(n_steps_total);
  d.lo = lo.data();
  d.hi = hi.data();
  d.step_rows = step_rows.data();
  if (n_clusters > 0) {
    d.n_clusters = n_clusters;
    d.cluster_shift = gamma.data();
    d.obs_cluster = obs_cluster.empty() ? nullptr : obs_cluster.data();
    d.stretch_cluster = stretch_cluster.empty() ? nullptr : stretch_cluster.data();
  }
  // About a millisecond for the largest data tried, against 10 ms or more
  // for the evaluation, so the indices are checked at every call rather
  // than trusted.
  mostr::paths::validate(d, false);

  mostr::paths::Output out;
  out.order = gradient != nullptr ? 1 : 0;
  if (gradient != nullptr) {
    // c(alpha, beta) first, then gamma, written by the kernel in place.
    gradient->assign(static_cast<size_t>(thresholds + p + n_clusters), 0.0);
    out.total_score = gradient->data();
    if (n_clusters > 0) out.cluster_scores = gradient->data() + thresholds + p;
  }
  mostr::paths::evaluate(d, alpha.data(), beta.data(), out);
  // An impossible observation makes the log density -Inf and the gradient
  // undefined; Stan rejects the point on the value, and a NaN gradient would
  // only add a second error.
  if (gradient != nullptr && !std::isfinite(out.loglik)) {
    std::fill(gradient->begin(), gradient->end(), 0.0);
  }
  return out.loglik;
}

}  // namespace stan_paths
}  // namespace mostr

#endif  // MOSTR_STAN_PATH_LOGLIK_HPP
