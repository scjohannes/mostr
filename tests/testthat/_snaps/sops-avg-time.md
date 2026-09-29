# avg_time rejects grouping names that would overwrite output columns

    Code
      avg_time(case$model, newdata = profiles, absorb = 3, times = 1:2, by = "state_set")
    Condition
      Error in `avg_time()`:
      ! Grouping or counterfactual variables use reserved output names: state_set. Rename these variables.

---

    Code
      avg_time(case$model, newdata = profiles, absorb = 3, times = 1:2, variables = list(
        state_set = c("a", "b")))
    Condition
      Error in `avg_time()`:
      ! Grouping or counterfactual variables use reserved output names: state_set. Rename these variables.

---

    Code
      avg_time(case$model, newdata = profiles, absorb = 3, times = 1:2, by = "time_unit",
      time_unit = "days")
    Condition
      Error in `avg_time()`:
      ! Grouping or counterfactual variables use reserved output names: time_unit. Rename these variables.

---

    Code
      avg_time(case$model, newdata = profiles, absorb = 3, times = 1:2, variables = list(
        time_unit = c("a", "b")), time_unit = "days")
    Condition
      Error in `avg_time()`:
      ! Grouping or counterfactual variables use reserved output names: time_unit. Rename these variables.

