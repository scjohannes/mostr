# data-raw

This directory contains R scripts used to prepare package datasets.

## Files

- `violet_baseline.R`: Creates `violet_baseline`, the starting states, ages,
  and SOFA scores of the 1,351 VIOLET patients, from Frank Harrell's
  `simlongord` data. It is internal data in `R/sysdata.rda`, used by
  `sim_markov("violet")`, and is not exported.

## Usage

To regenerate the internal data:

```r
# Install required packages
install.packages("Hmisc")
install.packages("usethis")

# Run the data preparation script
source("data-raw/violet_baseline.R")
```

The script downloads `simlongord` with `Hmisc::getHdata()`, stops if the
recovered patients disagree with Harrell's published summaries, and writes
`R/sysdata.rda`.

## Attribution

`violet_baseline` is derived from the `simlongord` dataset created
by Frank Harrell (<https://hbiostat.org/data/repo/simlongord>), available
through the Hmisc package. `simlongord` is based on the VIOLET trial.
