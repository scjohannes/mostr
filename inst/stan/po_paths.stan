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
// Parameters and priors are those of rmsb's 'lrmconppo' without partial
// proportional odds and without random effects: a Dirichlet prior on the
// state probabilities of a patient with average predictors, which gives the
// intercepts, and optional normal priors on contrasts of the coefficients.
// P(Y >= j + 1 | x) = inv_logit(alpha[j] + x * beta).

data {
  int<lower=0> N;                       // observed transitions
  int<lower=1> p;                       // predictors
  int<lower=2> k;                       // states
  matrix[N, p] X;                       // centered, orthogonalised predictors
  array[N] int<lower=1, upper=k> y;

  int<lower=0> G;                       // distinct design rows of stretches
  matrix[G, p] Xg;                      // on the scale of X
  int<lower=0> S;                       // distinct stretches
  int<lower=0> T;                       // steps of all stretches
  array[S] int<lower=1, upper=G> first_row;   // design row of step 1
  array[S] int<lower=1, upper=T> first_step;  // position of step 1
  array[S] int<lower=1> n_steps;
  vector<lower=0>[S] weight;            // patients sharing the stretch
  array[T] int<lower=1, upper=k> lo;    // possible states at a step: lo..hi
  array[T] int<lower=1, upper=k> hi;
  // design row of a step given the previous state; 0 where not needed
  array[T, k] int<lower=0, upper=G> step_rows;
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
  // ordered_logistic has P(Y > j) = inv_logit(eta - c[j]), so c = -alpha.
  if (N > 0) target += ordered_logistic_lpmf(y | X * beta, -alpha);

  if (S > 0) {
    // Probability of each state for each distinct design row, once.
    vector[G] eta = Xg * beta;
    matrix[G, k] prob;
    for (g in 1:G) {
      vector[k - 1] c = inv_logit(alpha + eta[g]);
      prob[g, 1] = 1 - c[1];
      for (j in 2:(k - 1)) prob[g, j] = c[j - 1] - c[j];
      prob[g, k] = c[k - 1];
    }
    for (s in 1:S) {
      int u = first_step[s];
      vector[k] a = rep_vector(0, k);
      for (c in lo[u]:hi[u]) a[c] = prob[first_row[s], c];
      for (t in 2:n_steps[s]) {
        vector[k] b = rep_vector(0, k);
        u += 1;
        for (j in lo[u - 1]:hi[u - 1]) {
          if (absorbing[j]) {
            // An absorbing state is never left.
            if (j >= lo[u] && j <= hi[u]) b[j] += a[j];
          } else {
            for (c in lo[u]:hi[u]) b[c] += a[j] * prob[step_rows[u, j], c];
          }
        }
        a = b;
      }
      target += weight[s] * log(sum(a));
    }
  }

  target += dirichlet_lpdf(pi | rep_vector(conc, k));
  if (cn > 0) target += normal_lpdf(C * beta | cmus, csds);
}
