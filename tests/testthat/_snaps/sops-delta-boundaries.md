# logit delta intervals distinguish structural boundaries

    Code
      nonstructural <- delta_interval_bounds(estimate = c(0, 1), standard_error = c(
        0.1, 0.1), conf_level = 0.95, conf_type = "logit")
    Condition
      Warning:
      Logit-delta limits are undefined for nonstructural boundary SOP estimates; their confidence limits were set to NA.

