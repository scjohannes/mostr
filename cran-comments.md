## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

Checked locally on Windows 11 with R 4.6.1
(`R CMD check --as-cran --run-donttest` on the built tarball). Update these
results after the win-builder and R-hub checks for the release.

## Examples

All examples run. The `\donttest{}` part of `?mostm` fits a Bayesian model
with 'rmsb', whose Stan sampling takes about 50 seconds; it is guarded by
`requireNamespace("rmsb", quietly = TRUE)`. The other `\donttest{}` parts
refit the model to bootstrap resamples and take a few seconds each.

## Method references

The methods are described in Rohde MD, French B, Stewart TG, Harrell FE Jr
(2024). "Bayesian Transition Models for Ordinal Longitudinal Outcomes."
Statistics in Medicine, 43(18), 3539-3561 <doi:10.1002/sim.10133>, cited in
the Description field.
