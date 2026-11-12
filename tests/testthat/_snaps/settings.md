# a setting that contradicts the fit is an error

    Code
      avg_sops(second, times = 1:4, absorb = 4, p2_var = "yprev")
    Condition
      Error:
      ! `p2_var = "yprev"` does not match the fitted model: `mostm()` fitted it with `p2_var = "ypprev"`. Leave `p2_var` out to use the fit's value, or refit the model with `mostm(..., p2_var = "yprev")`.
    Code
      avg_sops(first, times = 1:4, absorb = 4, p2_var = "ypprev")
    Condition
      Error:
      ! `p2_var = "ypprev"` does not match the fitted model: `mostm()` fitted it with `p2_var = NULL`. The model was fitted as a first-order model. For a second-order model, refit it with `mostm(..., p2_var = "ypprev")`.
    Code
      avg_sops(first, times = 1:4, absorb = 4, p_var = "ypprev")
    Condition
      Error:
      ! `p_var = "ypprev"` does not match the fitted model: `mostm()` fitted it with `p_var = "yprev"`. Leave `p_var` out to use the fit's value, or refit the model with `mostm(..., p_var = "ypprev")`.
    Code
      avg_sops(day, times = 1:4, absorb = 4, time_var = "time")
    Condition
      Error:
      ! `time_var = "time"` does not match the fitted model: `mostm()` fitted it with `time_var = "day"`. Leave `time_var` out to use the fit's value, or refit the model with `mostm(..., time_var = "time")`.
    Code
      avg_sops(first, times = 1:4, absorb = 4, id_var = "age")
    Condition
      Error:
      ! `id_var = "age"` does not match the fitted model: `mostm()` fitted it with `id_var = "id"`. Leave `id_var` out to use the fit's value, or refit the model with `mostm(..., id_var = "age")`.
    Code
      avg_sops(first, times = 1:4, absorb = 4, y_levels = c(1, 2, 4, 3))
    Condition
      Error:
      ! `y_levels` must be the fitted model's states in fitted order (1, 2, 3, 4), but it is 1, 2, 4, 3. Leave `y_levels` out to use the fitted states.

# `mostm()` records the settings and checks the previous-state columns

    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", p2_var = "missing")
    Condition
      Error:
      ! The column `missing` given as `p2_var` is not in `data`. Set `p2_var` to the column that holds each patient's state two visits back.
    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", p_var = "yprev", p2_var = "yprev")
    Condition
      Error:
      ! `p2_var` must name the column with the state two visits back, which differs from the previous-state column `p_var = "yprev"`.
    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", p_var = c("a", "b"))
    Condition
      Error:
      ! `p_var` must be a single column name.

