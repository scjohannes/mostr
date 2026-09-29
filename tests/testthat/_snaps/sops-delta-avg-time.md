# average-time analytical inference preserves restrictions and row alignment

    Code
      inferences(supplied, method = "delta")
    Condition
      Error:
      ! `vcov = "unconditional"` requires the stored fitted-patient cohort. For user-supplied `newdata`, use `vcov = "conditional"`.

---

    Code
      inferences(point, method = "delta", conf_type = "logit")
    Condition
      Error:
      ! Analytical average time requires `conf_type = "wald"`.

---

    Code
      inferences(grouped, method = "delta")
    Condition
      Error:
      ! Analytical delta inference with `by` is not yet supported.

