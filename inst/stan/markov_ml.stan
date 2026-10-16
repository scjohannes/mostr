// Bayesian first-order Markov transition model for ordinal states with
// missing or range-coded states, used by mostr::blrm_markov_ml().
//
// Parameters and priors follow rmsb's lrmconppo.stan: predictors are centered
// and QR-transformed, the intercepts come from the state probabilities `pi`
// of a patient with average predictors under a Dirichlet prior, contrasts of
// `beta` and `tau` may have normal priors, constrained partial proportional
// odds adds pposcore[j] * Z tau to the linear predictor of P(Y >= j), and
// optional patient random intercepts `gamma` are normal with standard
// deviation `sigmag`. The likelihood is the observed-data likelihood of
// orm_markov_ml(): fully observed transitions contribute the probability of
// the observed state, and every run of rows with a missing or range-coded
// state sums over the states the patient could have been in (scaled forward
// recursion). Given its random intercept, each patient's chain is first order.

functions {
  // Probabilities of the k states for linear predictor `lin` (without the
  // intercepts) and partial PO term `zeta`. The caller ensures that the
  // cumulative probabilities decrease across states.
  row_vector state_probs(vector alpha, real lin, real zeta, vector pposcore) {
    int k = num_elements(alpha) + 1;
    row_vector[k] out;
    vector[k - 1] cum;
    for (j in 2:k) cum[j - 1] = alpha[j - 1] + lin + pposcore[j] * zeta;
    out[1] = exp(log1m_inv_logit(cum[1]));
    for (j in 2:(k - 1)) out[j] = exp(log_inv_logit_diff(cum[j - 1], cum[j]));
    out[k] = inv_logit(cum[k - 1]);
    return out;
  }

  // 1 when P(Y >= j) decreases in j, which does not depend on `lin`.
  int ordered_cum(vector alpha, real zeta, vector pposcore) {
    int k = num_elements(alpha) + 1;
    for (j in 2:(k - 1)) {
      if (alpha[j - 1] + pposcore[j] * zeta <= alpha[j] + pposcore[j + 1] * zeta) {
        return 0;
      }
    }
    return 1;
  }
}

data {
  int<lower = 2> k;                            // number of states
  int<lower = 1> p;                            // number of predictors
  int<lower = 0> q;                            // number of partial PO predictors
  int<lower = 1> D;                            // number of design rows
  matrix[D, p] X;                              // centered, QR-transformed rows
  matrix[D, q] Z;                              // same for partial PO predictors
  vector[k] pposcore;                          // cppo scores of the states
  int<lower = 0> Nc;                           // number of patients with random
                                               // intercepts (0 = none)

  // Fully observed transitions, aggregated by design row and state (and by
  // patient when there are random intercepts).
  int<lower = 0> n_obs;
  array[n_obs] int<lower = 1, upper = D> obs_design;
  array[n_obs] int<lower = 1, upper = k> obs_state;
  vector<lower = 0>[n_obs] obs_count;
  array[Nc == 0 ? 0 : n_obs] int<lower = 1, upper = Nc> obs_cluster;

  // Runs of rows that need the recursion, each starting from an observed
  // state. Rows of a run are consecutive in `chain_design` and `chain_mask`,
  // whose entries index the `n_cd` design rows in `cd_rows`.
  int<lower = 1, upper = k> n_origin;
  array[n_origin] int<lower = 1, upper = k> origins;
  int<lower = 0, upper = k> n_absorbing;
  array[n_absorbing] int<lower = 1, upper = k> absorbing;
  int<lower = 0> n_seg;
  array[n_seg] int<lower = 1, upper = k> seg_baseline;
  array[n_seg] int<lower = 1> seg_start;
  array[n_seg] int<lower = 1> seg_end;
  array[Nc == 0 ? 0 : n_seg] int<lower = 1, upper = Nc> seg_cluster;
  int<lower = 0> n_cd;
  array[n_cd] int<lower = 1, upper = D> cd_rows;
  int<lower = 0> n_chain;
  array[n_chain, n_origin] int<lower = 1, upper = n_cd> chain_design;
  matrix<lower = 0, upper = 1>[n_chain, k] chain_mask;

  // Priors.
  real<lower = 0> conc;
  int<lower = 0> cn;
  matrix[cn, p] C;
  vector[cn] cmus;
  vector<lower = 0>[cn] csds;
  int<lower = 0> cn2;
  matrix[cn2, q] C2;
  vector[cn2] cmus2;
  vector<lower = 0>[cn2] csds2;
  int<lower = 1, upper = 2> psigma;  // 1 = t(4, rsdmean, rsdsd), 2 = exponential
  array[Nc == 0 ? 0 : 1] real<lower = 0> rsdmean;
  array[Nc == 0 || psigma == 2 ? 0 : 1] real<lower = 0> rsdsd;
}

parameters {
  vector[p] beta;          // coefficients on the transformed predictors
  vector[q] tau;           // partial PO coefficients on the transformed predictors
  simplex[k] pi;           // state probabilities for a patient with average predictors
  vector[Nc] gamma_raw;    // unscaled random intercepts
  array[Nc == 0 ? 0 : 1] real<lower = 0> sigmag;   // SD of the random intercepts
}

transformed parameters {
  vector[k - 1] alpha;
  vector[Nc] gamma = Nc == 0 ? gamma_raw : sigmag[1] * gamma_raw;
  for (j in 2:k) alpha[j - 1] = logit(sum(pi[j:k]));
}

model {
  vector[D] eta = X * beta;
  vector[D] zeta = rep_vector(0, D);
  if (q > 0) zeta = Z * tau;

  // Probability of the observed state; as in rmsb, a draw is rejected when
  // the cumulative probabilities around an observed state cross.
  for (i in 1:n_obs) {
    int d = obs_design[i];
    int a = obs_state[i];
    real lin = eta[d];
    real lp;
    if (Nc > 0) lin += gamma[obs_cluster[i]];
    if (a == 1) {
      lp = log1m_inv_logit(alpha[1] + lin + pposcore[2] * zeta[d]);
    } else if (a == k) {
      lp = log_inv_logit(alpha[k - 1] + lin + pposcore[k] * zeta[d]);
    } else {
      real upper = alpha[a - 1] + lin + pposcore[a] * zeta[d];
      real lower = alpha[a] + lin + pposcore[a + 1] * zeta[d];
      lp = upper > lower ? log_inv_logit_diff(upper, lower) : negative_infinity();
    }
    target += obs_count[i] * lp;
  }

  if (n_seg > 0) {
    // Every state's probability enters the sum, so all must be valid.
    int valid = 1;
    for (c in 1:n_cd) {
      if (ordered_cum(alpha, zeta[cd_rows[c]], pposcore) == 0) valid = 0;
    }
    if (valid == 0) {
      target += negative_infinity();
    } else {
      // Without random intercepts the transition rows are shared by patients.
      matrix[Nc == 0 ? n_cd : 0, k] P;
      if (Nc == 0) {
        for (c in 1:n_cd) {
          P[c] = state_probs(alpha, eta[cd_rows[c]], zeta[cd_rows[c]], pposcore);
        }
      }
      for (g in 1:n_seg) {
        real shift = 0;
        row_vector[k] f = rep_row_vector(0, k);
        if (Nc > 0) shift = gamma[seg_cluster[g]];
        f[seg_baseline[g]] = 1;
        for (t in seg_start[g]:seg_end[g]) {
          row_vector[k] next_f = rep_row_vector(0, k);
          real s;
          for (a in 1:n_origin) {
            // Origins without probability contribute nothing.
            if (f[origins[a]] > 0) {
              int c = chain_design[t, a];
              if (Nc == 0) {
                next_f += f[origins[a]] * P[c];
              } else {
                next_f += f[origins[a]] *
                  state_probs(alpha, eta[cd_rows[c]] + shift, zeta[cd_rows[c]], pposcore);
              }
            }
          }
          for (b in 1:n_absorbing) {
            next_f[absorbing[b]] += f[absorbing[b]];
          }
          next_f = next_f .* chain_mask[t];
          s = sum(next_f);
          target += log(s);
          f = next_f / s;
        }
      }
    }
  }

  target += dirichlet_lpdf(pi | rep_vector(conc, k));
  if (cn > 0) target += normal_lpdf(C * beta | cmus, csds);
  if (cn2 > 0) target += normal_lpdf(C2 * tau | cmus2, csds2);
  if (Nc > 0) {
    gamma_raw ~ std_normal();   // implies gamma ~ normal(0, sigmag)
    if (psigma == 1) sigmag ~ student_t(4, rsdmean[1], rsdsd[1]);
    else sigmag ~ exponential(1. / rsdmean[1]);
  }
}
