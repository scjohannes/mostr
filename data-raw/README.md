# data-raw

This directory contains R scripts used to prepare package datasets.

## Files

- `violet_baseline.R`: Creates the `violet_baseline` dataset, the starting
  states, ages, and SOFA scores of the 1,351 VIOLET patients, from
  Frank Harrell's `simlongord` data.

## Usage

To regenerate the package data:

```r
# Install required packages
install.packages("Hmisc")
install.packages("usethis")

# Run the data preparation script
source("data-raw/violet_baseline.R")
```

The script downloads `simlongord` with `Hmisc::getHdata()`, stops if the
recovered patients disagree with Harrell's published summaries, and writes
`data/violet_baseline.rda`.

## Attribution

The `violet_baseline` dataset is derived from the `simlongord` dataset created
by Frank Harrell (<https://hbiostat.org/data/repo/simlongord>), available
through the Hmisc package. `simlongord` is based on the VIOLET trial.
