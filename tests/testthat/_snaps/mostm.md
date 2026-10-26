# mostm() names the engine an argument belongs to

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", ppo = ~time)
    Condition
      Error:
      ! Not an argument of `engine = "orm"`: `ppo` (an argument of `engine = "blrm"`). See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", engine = "vglm",
      iter = 100, maxitt = 10)
    Condition
      Error:
      ! Not an argument of `engine = "vglm"`: `iter` (an argument of `engine = "blrm"`), `maxitt`. See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", engine = "blrm",
      constraints = list())
    Condition
      Error:
      ! Not an argument of `engine = "blrm"`: `constraints` (an argument of `engine = "vglm"`). See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", backend = "cmdstan")
    Condition
      Error:
      ! Not an argument of `engine = "orm"`: `backend` (an argument of `engine = "blrm"`). See `?mostm` for the arguments each engine accepts.

---

    Code
      mostm(ordered(y) ~ time + tx + yprev, data = data, id_var = "id", family = VGAM::cumulative(
        reverse = TRUE))
    Condition
      Error:
      ! For `engine = "orm"`, `family` names the link function of `rms::orm()`, such as "logistic". VGAM family objects apply to `engine = "vglm"`.

# fits without mostm() provenance must be refitted

    Code
      sops(old, times = 1:2, y_levels = 1:6, absorb = 6)
    Condition
      Error:
      ! Model-based Markov workflows require fits created by `mostm()`. Refit this model with `mostm(engine = "orm")`.

