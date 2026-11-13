# the prediction functions refuse the settings that the fit records

    Code
      sops(first, times = 1:4, p2_var = "ypprev")
    Condition
      Error in `sops()`:
      ! `p2_var` is not an argument of this function: the fitted model records the `p2_var` given to `mostm()`. Leave it out; to use another column, give it to `mostm()` and refit.
    Code
      avg_sops(first, times = 1:4, p_var = "yprev")
    Condition
      Error in `avg_sops()`:
      ! `p_var` is not an argument of this function: the fitted model records the `p_var` given to `mostm()`. Leave it out; to use another column, give it to `mostm()` and refit.
    Code
      state_time(first, times = 1:4, time_var = "time")
    Condition
      Error in `state_time()`:
      ! `time_var` is not an argument of this function: the fitted model records the `time_var` given to `mostm()`. Leave it out; to use another column, give it to `mostm()` and refit.
    Code
      comparisons(first, variables = "tx", times = 1:4, id_var = "id")
    Condition
      Error in `comparisons()`:
      ! `id_var` is not an argument of this function: the fitted model records the `id_var` given to `mostm()`. Leave it out; to use another column, give it to `mostm()` and refit.
    Code
      avg_comparisons(first, variables = "tx", times = 1:4, y_levels = 1:4)
    Condition
      Error in `avg_comparisons()`:
      ! `y_levels` is not an argument of this function: the fitted model records its states in fitted order. Leave it out.
    Code
      sops(first, times = 1:4, absorb = 4)
    Condition
      Error in `sops()`:
      ! `absorb` is not an argument of this function: the fitted model records the absorbing states given to `mostm()`. Leave it out; to change them, give them to `mostm(absorb = )` and refit.
    Code
      avg_sops(first, times = 1:4, gap_var = "gap")
    Condition
      Error in `avg_sops()`:
      ! `gap_var` is not an argument of this function: the fitted model records the `gap_var` given to `mostm()`. Leave it out; to use another column, give it to `mostm()` and refit.
    Code
      avg_state_time(first, times = 1:4, time_covariates = data.frame(time = 1:4))
    Condition
      Error in `avg_state_time()`:
      ! `time_covariates` was removed. Write time terms inline in the model formula given to `mostm()`, for example `rcs(time, 4)`, `log(time)` or `pol(time, 2)`; the prediction functions evaluate them at each requested time.

# diagnostic plots refuse a setting that contradicts the fit

    Code
      plot_transitions(second, times = 1:4, p2_var = "yprev")
    Condition
      Error in `plot_transitions()`:
      ! `p2_var = "yprev"` does not match the fitted model: `mostm()` fitted it with `p2_var = "ypprev"`. Leave `p2_var` out to use the fit's value, or refit the model with `mostm(..., p2_var = "yprev")`.
    Code
      plot_transitions(first, times = 1:4, p2_var = "ypprev")
    Condition
      Error in `plot_transitions()`:
      ! `p2_var = "ypprev"` does not match the fitted model: `mostm()` fitted it with `p2_var = NULL`. The model was fitted as a first-order model. For a second-order model, refit it with `mostm(..., p2_var = "ypprev")`.
    Code
      plot_correlation(first, times = 1:4, p_var = "ypprev")
    Condition
      Error in `plot_correlation()`:
      ! `p_var = "ypprev"` does not match the fitted model: `mostm()` fitted it with `p_var = "yprev"`. Leave `p_var` out to use the fit's value, or refit the model with `mostm(..., p_var = "ypprev")`.
    Code
      plot_variogram(day, times = 1:4, time_var = "time")
    Condition
      Error in `plot_variogram()`:
      ! `time_var = "time"` does not match the fitted model: `mostm()` fitted it with `time_var = "day"`. Leave `time_var` out to use the fit's value, or refit the model with `mostm(..., time_var = "time")`.
    Code
      plot_transitions(first, times = 1:4, id_var = "age")
    Condition
      Error in `plot_transitions()`:
      ! `id_var = "age"` does not match the fitted model: `mostm()` fitted it with `id_var = "id"`. Leave `id_var` out to use the fit's value, or refit the model with `mostm(..., id_var = "age")`.
    Code
      plot_correlation(first, times = 1:4, y_levels = c(1, 2, 4, 3))
    Condition
      Error in `plot_correlation()`:
      ! `y_levels` must be the fitted model's states in fitted order (1, 2, 3, 4), but it is 1, 2, 4, 3. Leave `y_levels` out to use the fitted states.

# `mostm()` records the settings and checks the previous-state columns

    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", p2_var = "missing")
    Condition
      Error in `mostm()`:
      ! The column `missing` given as `p2_var` is not in `data`. Set `p2_var` to the column that holds each patient's state two visits back.
    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", p_var = "yprev", p2_var = "yprev")
    Condition
      Error in `mostm()`:
      ! `p2_var` must name the column with the state two visits back, which differs from the previous-state column `p_var = "yprev"`.
    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", p_var = c("a", "b"))
    Condition
      Error in `mostm()`:
      ! `p_var` must be a single column name.

# `mostm()` refuses `absorb` labels that are not states of the response

    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", absorb = 5)
    Condition
      Error in `mostm()`:
      ! `absorb` must name states of the response, but "5" is not one. The states in the fitted rows of `data` are 1, 2, 3, 4.
    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", absorb = c("4", "dead"))
    Condition
      Error in `mostm()`:
      ! `absorb` must name states of the response, but "dead" is not one. The states in the fitted rows of `data` are 1, 2, 3, 4.
    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", absorb = c(4, NA))
    Condition
      Error in `mostm()`:
      ! `absorb` must be a vector of state labels without missing values, for example `absorb = "8"` or `absorb = c("home", "dead")`.

# `mostm()` refuses fitted rows that start in an absorbing state

    Code
      mostm(y ~ time + yprev, data = data, id_var = "id", absorb = c(2, 4))
    Condition
      Error in `mostm()`:
      ! `absorb` names the state "2", but 207 of the fitted rows of `data` have it as the previous state (`yprev`). A patient in an absorbing state can't leave it, so these rows are not transitions the model can learn from. Remove them before fitting, for example with `data[!data$yprev %in% c("2", "4"), ]`; `prepare_markov_data(absorbing_state = )` does this for the previous state.

---

    Code
      mostm(y ~ time + yprev + ypprev, data = no_start_in_2, id_var = "id", p2_var = "ypprev",
      absorb = c(2, 4))
    Condition
      Error in `mostm()`:
      ! `absorb` names the state "2", but 87 of the fitted rows of `data` have it as the state two visits back (`ypprev`). A patient in an absorbing state can't leave it, so these rows are not transitions the model can learn from. Remove them before fitting, for example with `data[!data$ypprev %in% c("2", "4"), ]`; `prepare_markov_data(absorbing_state = )` does this for the previous state.

# `mostm()` warns about a state seen only as an outcome that `absorb` leaves out

    Code
      fit <- mostm(y ~ time + yprev, data = data, id_var = "id")
    Condition
      Warning:
      The state "4" occurs as an outcome but never as the previous state (`yprev`). It may be absorbing, such as death; if so, name it in `absorb`, for example `mostm(..., absorb = "4")`. Otherwise the model may predict transitions out of it, so that patients who reach it leave again.

# `mostm()` checks that `gap_var` is the time since the previous row

    Code
      fit_gap(days)
    Condition
      Error in `mostm()`:
      ! The gap column `gap` must hold the time since the patient's previous row, in the units of `time`, because the prediction functions set it to the difference between consecutive requested times. It differs from that difference in 387 of the rows after a patient's first: for example, patient 1 at time 2 has gap 7, but the previous row is 1 earlier. Check that both columns use the same unit.
    Code
      fit_gap(text)
    Condition
      Error in `mostm()`:
      ! The gap column `gap` must be numeric: it holds the time since each patient's previous visit.
    Code
      mostm(formula, data = data, id_var = "id", absorb = 4, gap_var = "interval")
    Condition
      Error in `mostm()`:
      ! The column `interval` given as `gap_var` is not in `data`. Set `gap_var` to the column that holds the time since each patient's previous visit.
    Code
      mostm(formula, data = data, id_var = "id", absorb = 4, gap_var = c("gap",
        "time"))
    Condition
      Error in `mostm()`:
      ! `gap_var` must be a single column name.
    Code
      mostm(formula, data = factor_time, id_var = "id", absorb = 4, gap_var = "gap",
        first_followup_time = "1")
    Condition
      Error in `mostm()`:
      ! `gap_var` needs numeric time, but `time` is a factor. The prediction functions set the gap to the difference between consecutive requested times, which visit labels don't give, and there is no other way to give each visit's gap. Use a numeric time column, or leave `gap_var` out.

