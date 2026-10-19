#' Starting states, ages, and SOFA scores of the VIOLET patients
#'
#' The joint distribution of the starting state, age, and SOFA score of the
#' 1,351 patients in the
#' [VIOLET trial](https://doi.org/10.1056/NEJMoa1911124) of early high-dose
#' vitamin D3 versus placebo in critically ill adults. Each row is one
#' combination of the three characteristics, and `n` is the number of VIOLET
#' patients who had it. `sim_markov("violet")` samples simulated patients from
#' these rows (see [violet_preset]).
#'
#' @format A data frame with 858 rows and 4 variables:
#' \describe{
#'   \item{yprev}{State before the first day of follow-up: 2 = in hospital or
#'     another health care facility, 3 = on a ventilator or with acute
#'     respiratory distress syndrome (ARDS). The states used by
#'     `sim_markov("violet")` are 1 = at home, 2 = in hospital or facility,
#'     3 = on a ventilator or with ARDS, and 4 = dead.}
#'   \item{age}{Age in years.}
#'   \item{sofa}{Sequential Organ Failure Assessment (SOFA) score at
#'     enrollment, from 0 to 24; higher scores mean more severe organ
#'     failure.}
#'   \item{n}{Number of VIOLET patients with this starting state, age, and
#'     SOFA score. The counts sum to 1,351.}
#' }
#'
#' @details
#' To draw `size` patients as `sim_markov("violet")` does, sample rows with
#' replacement and with probabilities proportional to `n`:
#' `violet_baseline[sample(nrow(violet_baseline), size, replace = TRUE,
#' prob = violet_baseline$n), ]`. Use `rep(seq_len(nrow(violet_baseline)),
#' violet_baseline$n)` as row indices to get one row per VIOLET patient.
#'
#' The original VIOLET data are not public. Frank Harrell published a simulated
#' trial of 250,000 patients, `simlongord`, whose starting states, ages, and
#' SOFA scores were drawn with replacement from the 1,351 VIOLET patients with
#' a record on the first day (Harrell, 2022). Each VIOLET patient was therefore
#' drawn about 185 times, and dividing how often a combination occurs in
#' `simlongord` by 185 recovers how many VIOLET patients had it. The recovered
#' counts reproduce every summary that Harrell reported for these patients,
#' including the 474 patients on a ventilator or with ARDS, the number of
#' patients with each SOFA score, and the mean and percentiles of age. The
#' script `data-raw/violet_baseline.R` in the package sources performs these
#' steps.
#'
#' @source
#' Derived from `simlongord`, created by Frank Harrell and available through
#' `Hmisc::getHdata(simlongord)`.
#'
#' @references
#' National Heart, Lung, and Blood Institute PETAL Clinical Trials Network
#' (2019). "Early High-Dose Vitamin D3 for Critically Ill, Vitamin
#' D-Deficient Patients." *New England Journal of Medicine*, 381(26),
#' 2529-2540. \doi{10.1056/NEJMoa1911124}
#'
#' Harrell FE Jr (2022). "Simulating a Large Clinical Trial From VIOLET."
#' <https://hbiostat.org/data/repo/simlongord>
#'
#' @examples
#' head(violet_baseline)
#'
#' # Number of patients starting in each state
#' tapply(violet_baseline$n, violet_baseline$yprev, sum)
#'
#' # Mean age of the VIOLET patients
#' weighted.mean(violet_baseline$age, violet_baseline$n)
"violet_baseline"
