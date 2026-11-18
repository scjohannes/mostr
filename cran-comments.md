## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* Local: Windows 11, R 4.6.2 (`R CMD check --as-cran --run-donttest` on the
  built tarball).
* win-builder: R-devel and R-release.
* GitHub Actions: Ubuntu (R-release and R-devel), Windows (R-release) and
  macOS (R-release), with `--as-cran`; the compiled code's tests under the
  address and undefined-behaviour sanitizers and under Valgrind on Ubuntu.

## Examples

All examples run. The `\donttest{}` part of `?mostm` fits a Bayesian model
with 'rmsb', whose Stan sampling takes about 60 seconds; it is guarded by
`requireNamespace("rmsb", quietly = TRUE)`. The other `\donttest{}` parts
refit the model to bootstrap resamples and take a few seconds each.

## Method references

The methods are described in Rohde MD, French B, Stewart TG, Harrell FE Jr
(2024). "Bayesian Transition Models for Ordinal Longitudinal Outcomes."
Statistics in Medicine, 43(18), 3539-3561 <doi:10.1002/sim.10133>, cited in
the Description field.
