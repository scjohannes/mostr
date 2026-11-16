// Proportional odds Markov model when some states are not observed.
//
// The likelihood has two parts. A visit whose state and previous state are
// both observed contributes the probability of that transition, as in the
// 'lrmconppo' program of the rmsb package. A run of visits whose states are
// not all observed (a "stretch") contributes the sum, over every sequence of
// states the patient can have passed through, of the product of the
// transition probabilities along the sequence. That sum is computed by
// passing a vector of state probabilities from visit to visit.
//
// Both parts are computed in C++ by `path_loglik()`, which mostr compiles
// with this program (inst/include/mostr/stan_path_loglik_functions.hpp, the
// same code as the maximum likelihood fit of engine = "orm"). It returns the
// log-likelihood with its gradient computed analytically, which is much
// faster than differentiating the recursion by automatic differentiation.
//
// The priors are those of rmsb's 'lrmconppo' without partial proportional
// odds: a Dirichlet prior on the state probabilities pi of a patient with
// average predictors, which gives the intercepts, optional normal priors on
// contrasts of the coefficients and, with `cluster()`, a random intercept
// gamma[c] for each cluster c, non-centred as in rmsb (gamma = sigmag *
// gamma_raw, gamma_raw standard normal) with rmsb's prior on sigmag.
// P(Y >= j + 1 | x) = inv_logit(alpha[j] + x * beta + gamma[cluster]).
//
// The parameters are not rmsb's (rmsb declares beta, the simplex pi,
// gamma_raw and sigmag). This program declares z, gamma_raw and sigmag. The
// intercepts and coefficients are sampled as z, with
// theta = theta_hat + theta_L * z, where theta holds the log-odds
// u[j] = log(pi[j] / pi[k]), j < k, and the coefficients, and theta_hat and
// theta_L (the lower Cholesky factor of the covariance) come from the
// maximum likelihood fit of the same model without random intercepts,
// computed by mostr before sampling (`markov_path_stan_scale()`). z then has
// about unit scale and little correlation, which is what Stan's sampler
// assumes before it has adapted its step size and metric; in rmsb's
// parameters the coefficients' posterior SDs differ by a factor of about 20
// from those of the random intercepts and are strongly correlated with the
// intercepts. If that fit fails, theta_hat is a rough start and theta_L the
// identity: the same posterior, only sampled more slowly. Only the
// parametrisation changes, not the posterior.
//
// On 1,000 patients x 30 visits with half the states unobserved, 2 chains of
// 1000 warm-up and 1000 sampling iterations took 28,000 instead of 92,000
// gradient evaluations without `cluster()` (7 instead of 23 leapfrog steps
// per sampling iteration) and 81,000 instead of 175,000 with `cluster(id)`
// (15 instead of 31), both at Stan's adapt_delta 0.8, with the same
// posteriors; at the 0.9 that mostr uses with `cluster()`, 136,000
// (dev/decisions.md, 2026-11-15).

functions {
  // Matrices are passed stored by row. See the data block for the meaning
  // of each argument.
  real path_loglik(array[] real alpha, array[] real beta,
                   array[] real gamma,
                   data array[] real x_obs, data array[] int y_obs,
                   data array[] real w_obs, data array[] real x_gap,
                   data array[] real x_prev, data array[] int absorbing,
                   data array[] int first_row, data array[] int first_prev,
                   data array[] int first_step, data array[] int n_steps,
                   data array[] real weight, data array[] int lo,
                   data array[] int hi, data array[] int step_rows,
                   data array[] int obs_cluster,
                   data array[] int stretch_cluster);
}

data {
  int<lower=1> p;                       // predictors
  int<lower=2> k;                       // states
  // Observed transitions, counted: transitions with the same predictors and
  // the same state have the same probability, so each (design row, state)
  // pair enters once, with the number of transitions it stands for.
  int<lower=0> n_obs;                   // (design row, state) pairs
  array[n_obs * p] real x_obs;          // centered, orthogonalised predictors
  array[n_obs] int<lower=1, upper=k> y_obs;
  array[n_obs] real<lower=0> w_obs;

  // Design rows of the stretches, on the scale of x_obs. The row for a step
  // and previous state j is x_gap[step_rows[step, j]] + x_prev[j]; x_prev is
  // empty when x_gap holds every design row.
  int<lower=0> G;                       // rows of x_gap
  array[G * p] real x_gap;
  int<lower=0, upper=1> has_prev;
  array[has_prev * k * p] real x_prev;
  int<lower=0> S;                       // distinct stretches
  int<lower=0> T;                       // steps of all stretches
  array[S] int<lower=1, upper=G> first_row;   // design row of step 1
  array[S] int<lower=1, upper=k> first_prev;  // known state before step 1
  array[S] int<lower=1, upper=T> first_step;  // position of step 1
  array[S] int<lower=1> n_steps;
  array[S] real<lower=0> weight;        // patients sharing the stretch
  array[T] int<lower=1, upper=k> lo;    // possible states at a step: lo..hi
  array[T] int<lower=1, upper=k> hi;
  // design row of a step given the previous state, by step; 0 where not
  // needed
  array[T * k] int<lower=0, upper=G> step_rows;
  array[k] int<lower=0, upper=1> absorbing;

  // Random intercepts, as rmsb's data of the same names: Nc clusters (0 for
  // none), the cluster of each observed (design row, state) pair and of each
  // stretch, and the prior on their SD: psigma = 1 for
  // student_t(4, rsdmean, rsdsd), 2 for exponential with mean rsdmean.
  int<lower=0> Nc;
  array[Nc == 0 ? 0 : n_obs] int<lower=1, upper=Nc> obs_cluster;
  array[Nc == 0 ? 0 : S] int<lower=1, upper=Nc> stretch_cluster;
  int<lower=1, upper=2> psigma;
  array[Nc == 0 ? 0 : 1] real<lower=0> rsdmean;
  array[Nc == 0 || psigma == 2 ? 0 : 1] real<lower=0> rsdsd;

  real<lower=0> conc;
  int<lower=0> cn;                      // contrasts with a prior
  matrix[cn, p] C;
  vector[cn] cmus;
  vector[cn] csds;

  // The map from z to theta = (u, beta): the maximum likelihood estimate
  // and the lower Cholesky factor of its covariance. The full factor, not
  // only the standard errors: with the standard errors alone (a diagonal
  // theta_L) the program needed 1.6 to 2 times the gradient evaluations of
  // rmsb's parametrisation without `cluster()` (182,000 against 92,000 in
  // the run above), because the intercepts and coefficients are strongly
  // correlated and a diagonal map leaves those correlations to the sampler.
  vector[k - 1 + p] theta_hat;
  matrix[k - 1 + p, k - 1 + p] theta_L;
}

parameters {
  vector[k - 1 + p] z;                           // theta = theta_hat + theta_L * z
  vector[Nc] gamma_raw;                          // unscaled random effects
  array[Nc == 0 ? 0 : 1] real<lower=0> sigmag;   // SD of random effects
}

transformed parameters {
  vector[k - 1 + p] theta = theta_hat + theta_L * z;
  vector[k] ue = append_row(theta[1:(k - 1)], 0);  // log(pi / pi[k])
  vector[k] pi = softmax(ue);                      // sums to 1 by construction
  vector[p] beta = theta[k:(k - 1 + p)];
  vector[k - 1] alpha;
  vector[Nc] gamma;
  // alpha[j - 1] = logit(P(Y >= j)), computed from the log-odds and not as
  // logit(sum(pi[j:k])): at the mean predictors the worst state can have
  // pi[k] of about 1e-7, where the sum would lose most of its digits.
  for (j in 2:k) {
    alpha[j - 1] = log_sum_exp(ue[j:k]) - log_sum_exp(ue[1:(j - 1)]);
  }
  gamma = Nc == 0 ? gamma_raw : sigmag[1] * gamma_raw;
}

model {
  // The Dirichlet prior is a density of q = pi[1:(k - 1)] (pi[k] is
  // 1 - sum(q)), but the sampled scale is u, so the log of the Jacobian
  // determinant of u -> q is added. With q = softmax(u, 0)[1:(k - 1)],
  // dq[i] / du[j] = q[i] * ((i == j) - q[j]), the matrix diag(q) - q q',
  // whose determinant is prod(q) * (1 - sum(q)) (matrix determinant lemma),
  // that is prod(pi): the term is sum(log(pi)), computed as
  // sum(log_softmax(ue)) so that a tiny pi[k] does not underflow. The map
  // from z to theta is affine, so its Jacobian, det(theta_L), is a constant
  // and is left out. Without this term the program would sample another
  // posterior; the slow tests check that its log density minus rmsb's is
  // sum(log(pi)) plus a constant.
  target += sum(log_softmax(ue));
  if (Nc > 0) {
    gamma_raw ~ std_normal();  // implies gamma ~ normal(0, sigmag)
    if (psigma == 1) sigmag ~ student_t(4, rsdmean[1], rsdsd[1]);
    else sigmag ~ exponential(1. / rsdmean[1]);
  }
  target += path_loglik(to_array_1d(alpha), to_array_1d(beta),
                        to_array_1d(gamma), x_obs, y_obs, w_obs, x_gap,
                        x_prev, absorbing, first_row, first_prev, first_step,
                        n_steps, weight, lo, hi, step_rows, obs_cluster,
                        stretch_cluster);
  target += dirichlet_lpdf(pi | rep_vector(conc, k));
  if (cn > 0) target += normal_lpdf(C * beta | cmus, csds);
}
