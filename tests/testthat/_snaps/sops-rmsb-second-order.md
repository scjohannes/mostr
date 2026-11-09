# grouped blrm summaries reject entirely incomplete grouping rows

    Code
      sops(model, newdata = newdata, times = 1:2, by = "grp", n_draws = 2)
    Condition
      Error:
      ! `by` names `grp`, which is missing for 2 of the 2 profiles being averaged. A profile with a missing `by` value belongs to no group, so it would be left out of every average. Drop or impute those patients first, for example by removing them from `newdata` or from the data given to `mostm()`.

---

    Code
      avg_sops(model, newdata = newdata, variables = list(tx = c(0, 1)), times = 1:2,
      by = "grp", n_draws = 2)
    Condition
      Error:
      ! `by` names `grp`, which is missing for 2 of the 2 profiles being averaged. A profile with a missing `by` value belongs to no group, so it would be left out of every average. Drop or impute those patients first, for example by removing them from `newdata` or from the data given to `mostm()`.

