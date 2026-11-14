# mostm() names the engine an argument belongs to

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", ppo = ~time)
    Condition
      Error in `mostm()`:
      ! Not an argument of `engine = "orm"`: `ppo` (an argument of `engine = "blrm"`). See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", engine = "vglm",
      iter = 100, maxitt = 10)
    Condition
      Error in `mostm()`:
      ! Not an argument of `engine = "vglm"`: `iter` (an argument of `engine = "blrm"`), `maxitt`. See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", engine = "blrm",
      constraints = list())
    Condition
      Error in `mostm()`:
      ! Not an argument of `engine = "blrm"`: `constraints` (an argument of `engine = "vglm"`). See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", backend = "cmdstan")
    Condition
      Error in `mostm()`:
      ! Not an argument of `engine = "orm"`: `backend` (an argument of `engine = "blrm"`). See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", family = VGAM::cumulative(
        reverse = TRUE))
    Condition
      Error in `mostm()`:
      ! For `engine = "orm"`, `family` names the link function of `rms::orm()`, such as "logistic". VGAM family objects apply to `engine = "vglm"`.

# fits without mostm() provenance must be refitted

    Code
      sops(old, times = 1:2)
    Condition
      Error in `sops()`:
      ! Model-based Markov workflows require fits created by `mostm()`. Refit this model with `mostm(engine = "orm")`.

# mostm() refuses bases that rms would evaluate again on new data

    Code
      mostm(y ~ splines::ns(time, 3) + yprev, data = data, id_var = "id", absorb = 4)
    Condition
      Error in `mostm()`:
      ! With `rms::orm()`, `splines::ns(time, 3)` would be evaluated again on the prediction rows, which gives a different basis without any error, so the state occupancy probabilities would be wrong. Use `rms::rcs()` for a spline, for example `rcs(time, 4)`, or `rms::pol()` for a polynomial, for example `pol(time, 2)`: both keep the knots or terms of the fit. `engine = "vglm"` keeps the fitted basis of these terms.
    Code
      mostm(y ~ time + bs(age, 3) + yprev, data = data, id_var = "id", absorb = 4)
    Condition
      Error in `mostm()`:
      ! With `rms::orm()`, `bs(age, 3)` would be evaluated again on the prediction rows, which gives a different basis without any error, so the state occupancy probabilities would be wrong. Use `rms::rcs()` for a spline, for example `rcs(age, 4)`, or `rms::pol()` for a polynomial, for example `pol(age, 2)`: both keep the knots or terms of the fit. `engine = "vglm"` keeps the fitted basis of these terms.
    Code
      mostm(y ~ stats::poly(age, 2) + ns(time, 2) + yprev, data = data, id_var = "id",
      absorb = 4)
    Condition
      Error in `mostm()`:
      ! With `rms::orm()`, `stats::poly(age, 2)` and `ns(time, 2)` would each be evaluated again on the prediction rows, which gives a different basis without any error, so the state occupancy probabilities would be wrong. Use `rms::rcs()` for a spline, for example `rcs(age, 4)`, or `rms::pol()` for a polynomial, for example `pol(age, 2)`: both keep the knots or terms of the fit. `engine = "vglm"` keeps the fitted basis of these terms.

