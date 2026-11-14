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
  int<lower=1> p;                       // predictors
  int<lower=2> k;                       // states
  // Observed transitions, counted: transitions with the same predictors and
  // the same state have the same probability, so each distinct design row
  // is evaluated once and each (design row, state) pair enters with the
  // number of transitions it stands for.
  int<lower=0> Go;                      // distinct design rows
  matrix[Go, p] Xo;                     // centered, orthogonalised predictors
  int<lower=0> n_cells;                 // (design row, state) pairs observed
  // position of the pair in the Go x k matrix of probabilities, by column
  array[n_cells] int<lower=1, upper=Go * k> cell;
  vector<lower=0>[n_cells] cell_count;

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

transformed data {
  // Whether the previous step's range of a step includes an absorbing state.
  array[T] int has_absorbing = rep_array(0, T);
  for (s in 1:S) {
    for (t in 2:n_steps[s]) {
      int u = first_step[s] + t - 1;
      for (j in lo[u - 1]:hi[u - 1]) {
        if (absorbing[j]) has_absorbing[u] = 1;
      }
    }
  }
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
  if (n_cells > 0) {
    matrix[Go, k - 1] c = inv_logit(rep_matrix(Xo * beta, k - 1)
                                    + rep_matrix(alpha', Go));
    matrix[Go, k] prob;
    prob[, 1] = 1 - c[, 1];
    for (j in 2:(k - 1)) prob[, j] = c[, j - 1] - c[, j];
    prob[, k] = c[, k - 1];
    target += dot_product(cell_count, log(to_vector(prob)[cell]));
  }

  if (S > 0) {
    // Probability of each state for each distinct design row, once, as
    // whole-matrix operations: design rows are shared between steps and
    // stretches, and few large operations differentiate much faster than
    // many scalar ones.
    matrix[G, k - 1] c = inv_logit(rep_matrix(Xg * beta, k - 1)
                                   + rep_matrix(alpha', G));
    matrix[G, k] prob;
    prob[, 1] = 1 - c[, 1];
    for (j in 2:(k - 1)) prob[, j] = c[, j - 1] - c[, j];
    prob[, k] = c[, k - 1];

    for (s in 1:S) {
      int u = first_step[s];
      // State probabilities at the current step; zero outside its range.
      row_vector[k] a = rep_row_vector(0, k);
      real scale = 0;
      a[lo[u]:hi[u]] = prob[first_row[s], lo[u]:hi[u]];
      for (t in 2:n_steps[s]) {
        row_vector[k] b = rep_row_vector(0, k);
        int from_lo = lo[u];
        int from_hi = hi[u];
        u += 1;
        if (has_absorbing[u]) {
          // A range that includes an absorbing state: that state is never
          // left, so its row of the transition matrix is the identity.
          for (j in from_lo:from_hi) {
            if (absorbing[j]) {
              if (j >= lo[u] && j <= hi[u]) b[j] += a[j];
            } else {
              b[lo[u]:hi[u]] += a[j] * prob[step_rows[u, j], lo[u]:hi[u]];
            }
          }
        } else {
          b[lo[u]:hi[u]] = a[from_lo:from_hi]
                           * prob[step_rows[u, from_lo:from_hi], lo[u]:hi[u]];
        }
        a = b;
        // Rescale now and then, so that a long stretch can't underflow.
        if (t % 16 == 0) {
          real total = sum(a);
          scale += log(total);
          a /= total;
        }
      }
      target += weight[s] * (scale + log(sum(a)));
    }
  }

  target += dirichlet_lpdf(pi | rep_vector(conc, k));
  if (cn > 0) target += normal_lpdf(C * beta | cmus, csds);
}
