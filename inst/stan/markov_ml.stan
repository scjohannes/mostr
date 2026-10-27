// Bayesian first-order Markov transition model for ordinal states with
// missing or range-coded states, used by mostr::mostm(engine = "blrm") with
// `missing_states` or an Ocens() response.
//
// Parameters and priors follow rmsb's lrmconppo.stan: predictors are centered
// and QR-transformed, the intercepts come from the state probabilities `pi`
// of a patient with average predictors under a Dirichlet prior, contrasts of
// `beta` and `tau` may have normal priors, constrained partial proportional
// odds adds pposcore[j] * Z tau to the linear predictor of P(Y >= j), and
// optional patient random intercepts `gamma` are normal with standard
// deviation `sigmag`. The likelihood is the observed-data likelihood that
// the frequentist engines maximise: fully observed transitions contribute the
// probability of the observed state, and every run of rows with a missing or
// range-coded state sums over the states the patient could have been in
// (scaled forward recursion). Given its random intercept, each patient's chain is first order.
// The recursion computes only the transition probabilities it uses: from the
// previous states that can have probability (index lists built in R) to the
// states each row allows.

functions {
  // Linear predictor of P(Y >= j), j = 2..k, of a transition: the entry of
  // its design row (built for its row and previous state) in the linear
  // predictors without offset `x0`, plus the offset of its row (0, or the
  // patient's random intercept).
  real cum_logit(int j, row_vector x0, real offset) {
    return x0[j - 1] + offset;
  }

  // Probabilities of states lo..hi of that transition, as differences of
  // the cumulative probabilities P(Y >= j). Each of these is held as its
  // smaller tail, P(Y >= j) or P(Y < j), so that probabilities near 1 are
  // not subtracted from each other. The caller ensures that the cumulative
  // probabilities decrease across states.
  row_vector range_probs(int lo, int hi, row_vector x0, real offset) {
    int k = num_elements(x0) + 1;
    row_vector[hi - lo + 1] out;
    // P(Y >= j) is 1 - tail when `upper` (above 1/2) and tail otherwise.
    real tail = 0;   // P(Y >= 1) = 1
    int upper = 1;
    if (lo > 1) {
      real x = cum_logit(lo, x0, offset);
      upper = x > 0;
      tail = inv_logit(upper ? -x : x);
    }
    for (j in lo:hi) {
      real tail_next = 0;   // P(Y >= k + 1) = 0
      int upper_next = 0;
      if (j < k) {
        real x = cum_logit(j + 1, x0, offset);
        upper_next = x > 0;
        tail_next = inv_logit(upper_next ? -x : x);
      }
      if (upper == upper_next) {
        out[j - lo + 1] = upper ? tail_next - tail : tail - tail_next;
      } else if (upper) {
        out[j - lo + 1] = 1 - tail - tail_next;
      } else {
        out[j - lo + 1] = tail + tail_next - 1;   // crossing
      }
      tail = tail_next;
      upper = upper_next;
    }
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
  // state `seg_baseline` and covering rows seg_start..seg_end.
  int<lower = 0> n_seg;
  array[n_seg] int<lower = 1, upper = k> seg_baseline;
  array[n_seg] int<lower = 1> seg_start;
  array[n_seg] int<lower = 1> seg_end;
  int<lower = 0> n_chain;
  matrix<lower = 0, upper = 1>[n_chain, k] chain_mask;   // allowed states
  array[Nc == 0 ? 0 : n_chain] int<lower = 1, upper = Nc> chain_cluster;
  // Lowest and highest allowed state of each row; rows that exclude states in
  // between are `gappy`, rows that allow every state are `full`.
  array[n_chain] int<lower = 1, upper = k> row_lo;
  array[n_chain] int<lower = 1, upper = k> row_hi;
  array[n_chain] int<lower = 0, upper = 1> row_gappy;
  array[n_chain] int<lower = 0, upper = 1> row_full;
  // Cells: the previous states of a row that can have probability, with the
  // design row (`cd_rows[cell_cd]`) of the transition from each. Cells of row
  // t are cell_start[t]..cell_end[t].
  array[n_chain] int<lower = 1> cell_start;
  array[n_chain] int<lower = 0> cell_end;
  int<lower = 0> n_cell;
  array[n_cell] int<lower = 1, upper = k> cell_origin;
  int<lower = 0> n_cd;
  array[n_cell] int<lower = 1, upper = n_cd> cell_cd;
  array[n_cd] int<lower = 1, upper = D> cd_rows;
  array[n_cd] int<lower = 1, upper = k> cd_lo;  // states needed by the cells
  array[n_cd] int<lower = 1, upper = k> cd_hi;  // of each design row
  // Design rows whose cumulative probabilities must not cross.
  int<lower = 0> n_check;
  array[n_check] int<lower = 1, upper = D> check_rows;

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

transformed data {
  // Without random effects the offset of every recursion row is 0.
  int shared = Nc == 0;
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
    for (c in 1:n_check) {
      if (ordered_cum(alpha, zeta[check_rows[c]], pposcore) == 0) valid = 0;
    }
    if (valid == 0) {
      target += negative_infinity();
    } else {
      // Linear predictors of P(Y >= j), j = 2..k (columns), of the
      // transitions of the design rows the cells use, before the offset of
      // their row. The offset is the patient's random intercept, and the
      // only part of the recursion that random effects change. Without them
      // it is 0 on every row, so each design row's probabilities are computed
      // once and shared.
      matrix[n_cd, k - 1] L = rep_matrix(alpha', n_cd)
        + rep_matrix(eta[cd_rows], k - 1);
      vector[shared ? 0 : n_chain] offset;
      matrix[shared ? n_cd : 0, k] P;
      if (q > 0) L += zeta[cd_rows] * pposcore[2:k]';
      if (shared) {
        for (c in 1:n_cd) {
          P[c, cd_lo[c]:cd_hi[c]] = range_probs(cd_lo[c], cd_hi[c], L[c], 0);
        }
      } else {
        offset = gamma[chain_cluster];
      }
      for (g in 1:n_seg) {
        row_vector[k] f = rep_row_vector(0, k);
        f[seg_baseline[g]] = 1;
        for (t in seg_start[g]:seg_end[g]) {
          int lo = row_lo[t];
          int hi = row_hi[t];
          int first = cell_start[t];
          int last = cell_end[t];
          row_vector[k] next_f = rep_row_vector(0, k);
          if (last >= first) {
            // Transition probabilities from the cells' previous states.
            matrix[last - first + 1, hi - lo + 1] Pt;
            for (c in first:last) {
              if (shared) {
                Pt[c - first + 1] = P[cell_cd[c], lo:hi];
              } else {
                Pt[c - first + 1] = range_probs(lo, hi, L[cell_cd[c]], offset[t]);
              }
            }
            if (t == seg_start[g]) {
              // The run starts with probability 1 at its baseline, its only
              // cell.
              next_f[lo:hi] = Pt[1];
            } else {
              next_f[lo:hi] = f[cell_origin[first:last]] * Pt;
            }
          }
          if (row_gappy[t]) next_f = next_f .* chain_mask[t];
          if (row_full[t]) {
            // The probabilities of all states add up to 1.
            f = next_f;
          } else {
            real s = sum(next_f[lo:hi]);
            target += log(s);
            // The last row's distribution is not used.
            if (t < seg_end[g]) f = next_f / s;
          }
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
