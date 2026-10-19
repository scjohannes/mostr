# sim_markov validates preset inputs

    Code
      sim_markov("actt2", n_patients = 0)
    Condition
      Error in `sim_markov()`:
      ! n_patients must be a positive integer

---

    Code
      sim_markov("actt2", treatment_prob = 2)
    Condition
      Error in `sim_markov()`:
      ! treatment_prob must be a numeric scalar between 0 and 1

---

    Code
      sim_markov("actt2", treatment_effect = Inf)
    Condition
      Error in `sim_markov()`:
      ! treatment_effect must be a finite numeric scalar

---

    Code
      sim_markov("actt2", treatment_effect_decay = -1)
    Condition
      Error in `sim_markov()`:
      ! treatment_effect_decay must be a nonnegative finite numeric scalar

---

    Code
      sim_markov("actt2", follow_up_time = 0)
    Condition
      Error in `sim_markov()`:
      ! follow_up_time must be a positive integer

