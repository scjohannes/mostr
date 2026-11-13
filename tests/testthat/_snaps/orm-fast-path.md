# orm SOPs report a state count that differs from the fit and a missing absorbing state

    Code
      mostr:::predict_markov_sops(fit, baseline, times = 1, settings = fewer_states)
    Condition
      Error:
      ! `y_levels` defines 5 states, but the fitted model defines 6 states through 5 threshold coefficients.

---

    Code
      mostr::soprob_markov(fit, baseline[1L, , drop = FALSE], times = 1)
    Condition
      Error in `mostr::soprob_markov()`:
      ! Model prediction requires transitions from state level(s) not represented in the fitted transition data: 6. No fitted row has such a state as its previous state. If it is absorbing, such as death, name it in `mostm(absorb = )` and refit.

