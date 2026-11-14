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
// Parameters and priors are those of rmsb's 'lrmconppo' without partial
// proportional odds and without random effects: a Dirichlet prior on the
// state probabilities of a patient with average predictors, which gives the
// intercepts, and optional normal priors on contrasts of the coefficients.
// P(Y >= j + 1 | x) = inv_logit(alpha[j] + x * beta).

functions {
  // Matrices are passed stored by row. See the data block for the meaning
  // of each argument.
  real path_loglik(array[] real alpha, array[] real beta,
                   data array[] real x_obs, data array[] int y_obs,
                   data array[] real w_obs, data array[] real x_gap,
                   data array[] real x_prev, data array[] int absorbing,
                   data array[] int first_row, data array[] int first_prev,
                   data array[] int first_step, data array[] int n_steps,
                   data array[] real weight, data array[] int lo,
                   data array[] int hi, data array[] int step_rows);
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

  real<lower=0> conc;
  int<lower=0> cn;                      // contrasts with a prior
  matrix[cn, p] C;
  vector[cn] cmus;
  vector[cn] csds;
}

parameters {
  vector[p] beta;
  simplex[k] pi;
}

transformed parameters {
  vector[k - 1] alpha;
  for (j in 2:k) alpha[j - 1] = logit(sum(pi[j:k]));
}

model {
  target += path_loglik(to_array_1d(alpha), to_array_1d(beta), x_obs, y_obs,
                        w_obs, x_gap, x_prev, absorbing, first_row,
                        first_prev, first_step, n_steps, weight, lo, hi,
                        step_rows);
  target += dirichlet_lpdf(pi | rep_vector(conc, k));
  if (cn > 0) target += normal_lpdf(C * beta | cmus, csds);
}
