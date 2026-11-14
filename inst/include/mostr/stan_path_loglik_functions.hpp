// The C++ definition of the function `path_loglik()` that
// inst/stan/po_paths.stan declares without a body. Include this file inside
// the namespace of the compiled model, after stan_path_loglik.hpp: Stan
// calls the function unqualified from the model's code, so it must be
// declared in that namespace. It has no include guard on purpose.
//
// Two overloads: `double` arguments, for evaluating the log density alone
// (`rstan::log_prob()`), and `var` arguments, for the gradient during
// sampling. The var overload returns one value whose partial derivatives
// with respect to every element of alpha and beta are those computed by
// path_likelihood.hpp (stan::math::precomputed_gradients()), so Stan's
// automatic differentiation does not see the recursion.
//
// The argument types are those of the Stan declaration
//   real path_loglik(array[] real alpha, array[] real beta,
//                    data array[] real x_obs, data array[] int y_obs, ...)
// called with `to_array_1d(alpha)`, `to_array_1d(beta)` and data arrays. As
// non-template functions they match the call exactly, whatever templates
// stanc may declare for the function.

inline double path_loglik(
    const std::vector<double>& alpha, const std::vector<double>& beta,
    const std::vector<double>& x_obs, const std::vector<int>& y_obs,
    const std::vector<double>& w_obs, const std::vector<double>& x_gap,
    const std::vector<double>& x_prev, const std::vector<int>& absorbing,
    const std::vector<int>& first_row, const std::vector<int>& first_prev,
    const std::vector<int>& first_step, const std::vector<int>& n_steps,
    const std::vector<double>& weight, const std::vector<int>& lo,
    const std::vector<int>& hi, const std::vector<int>& step_rows,
    std::ostream* pstream__) {
  return mostr::stan_paths::path_loglik(
      alpha, beta, x_obs, y_obs, w_obs, x_gap, x_prev, absorbing, first_row,
      first_prev, first_step, n_steps, weight, lo, hi, step_rows, nullptr);
}

inline stan::math::var path_loglik(
    const std::vector<stan::math::var>& alpha,
    const std::vector<stan::math::var>& beta,
    const std::vector<double>& x_obs, const std::vector<int>& y_obs,
    const std::vector<double>& w_obs, const std::vector<double>& x_gap,
    const std::vector<double>& x_prev, const std::vector<int>& absorbing,
    const std::vector<int>& first_row, const std::vector<int>& first_prev,
    const std::vector<int>& first_step, const std::vector<int>& n_steps,
    const std::vector<double>& weight, const std::vector<int>& lo,
    const std::vector<int>& hi, const std::vector<int>& step_rows,
    std::ostream* pstream__) {
  std::vector<double> alpha_value(alpha.size());
  std::vector<double> beta_value(beta.size());
  for (size_t i = 0; i < alpha.size(); ++i) alpha_value[i] = alpha[i].val();
  for (size_t i = 0; i < beta.size(); ++i) beta_value[i] = beta[i].val();
  std::vector<double> gradient;
  const double value = mostr::stan_paths::path_loglik(
      alpha_value, beta_value, x_obs, y_obs, w_obs, x_gap, x_prev, absorbing,
      first_row, first_prev, first_step, n_steps, weight, lo, hi, step_rows,
      &gradient);
  std::vector<stan::math::var> operands;
  operands.reserve(alpha.size() + beta.size());
  operands.insert(operands.end(), alpha.begin(), alpha.end());
  operands.insert(operands.end(), beta.begin(), beta.end());
  return stan::math::precomputed_gradients(value, operands, gradient);
}
