# public delta scope enforces fixed targets and patient clustering

    Code
      inferences(fixed, method = "delta", vcov = "unconditional")
    Condition
      Error in `inferences()`:
      ! `sops()` delta inference supports only `vcov = "conditional"` or a coefficient covariance matrix.

---

    Code
      inferences(avg, method = "delta", vcov = c("conditional", "unconditional"))
    Condition
      Error in `inferences()`:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

# delta vcov validates strings and does not silently change targets

    Code
      inferences(object, method = "delta")
    Condition
      Error in `inferences()`:
      ! `vcov = "unconditional"` requires the stored fitted-patient cohort. For user-supplied `newdata`, use `vcov = "conditional"`.

---

    Code
      inferences(object, method = "delta", vcov = NA_character_)
    Condition
      Error in `inferences()`:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(object, method = "delta", vcov = character())
    Condition
      Error in `inferences()`:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(object, method = "delta", vcov = "cond")
    Condition
      Error in `inferences()`:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(object, method = "delta", target = "empirical")
    Condition
      Error in `inferences()`:
      ! unused argument (target = "empirical")

