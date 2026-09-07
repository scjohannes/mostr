# mostr

Markov Ordinal State Transition models.

Compute and plot state occupancy probabilities (SOPs) and contrasts from first- and second-order Markov models of ordinal health-state trajectories.

## Installation

Install the development version from GitHub with `pak`:

```r
install.packages("pak")
pak::pak("scjohannes/mostr")
```

## SOP Workflow

This example simulates ordinal outcomes from a viral respiratory disease trials, fits a proportional-odds Markov transition model with `orm_markov()`, estimates marginal treatment-arm SOPs, adds uncertainty, and plots the result.

```r
library(mostr)
library(rms)
library(ggplot2)

set.seed(20260526)

trial <- sim_actt2_markov(
  n_patients = 300,
  follow_up_time = 28,
  treatment_prob = 0.5,
  treatment_effect = log(0.8),
  seed = 20260526
)

markov_data <- prepare_markov_data(trial, absorbing_state = 8)

dd <- datadist(markov_data)
options(datadist = "dd")

fit <- orm_markov(
  y ~ rms::rcs(time, 4) + tx + yprev,
  data = markov_data,
  id_var = "id",
  opt_method = "LM",
  scale = TRUE
)

sop <- avg_sops(
  fit,
  variables = list(tx = c(0, 1)),
  times = 1:28,
  y_levels = fit$yunique,
  absorb = "8"
) |> inferences(
  method = "delta"
)

plot_sops(sop, facet_var = NULL, linetype_var = "tx") +
  labs(
    x = "Day",
    y = "State occupancy probability"
  )
```

Use `orm_markov()`, `vglm_markov()`, or `blrm_markov()` for package model-based
SOP and diagnostic workflows. Raw `rms`, `VGAM`, or `rmsb` fits do not carry the
wrapper provenance and stored-data contracts those workflows require.


## Learn More

After installation, see:

```r
vignette("full-po-sops", package = "mostr")
vignette("analytical-confidence-intervals", package = "mostr")
vignette("partial-po-vglm-sops", package = "mostr")
vignette("second-order-orm-sops", package = "mostr")
vignette("many-levels-previous-state-spline", package = "mostr")
vignette("factor-time-orm-sops", package = "mostr")
vignette("bayesian-rmsb-sops", package = "mostr")
vignette("mvn-vs-bootstrap-orm-sops", package = "mostr")
```

## Data Provenance

The included `violet_baseline` dataset is derived from `Hmisc::simlongord`,
created by Frank Harrell and based on the VIOLET trial.