## Code to prepare `violet_baseline`, internal data in R/sysdata.rda
##
## `violet_baseline` holds the joint distribution of the starting state, age,
## and SOFA score of the 1,351 patients in the VIOLET trial of early high-dose
## vitamin D3 versus placebo in critically ill adults. Each row is one
## combination of the three characteristics, and `n` is the number of VIOLET
## patients who had it (858 rows; the counts sum to 1,351). Columns:
## - yprev: state before the first day of follow-up, 2 = in hospital or
##   another health care facility, 3 = on a ventilator or with acute
##   respiratory distress syndrome (ARDS). `sim_markov("violet")` uses
##   1 = home, 2 = hospital or facility, 3 = ventilator or ARDS, 4 = dead.
## - age: age in years.
## - sofa: Sequential Organ Failure Assessment (SOFA) score at enrollment,
##   from 0 to 24; higher scores mean more severe organ failure.
## - n: number of VIOLET patients with this combination.
## `violet_markov_parameters()` in R/simulate-markov.R passes the table to
## `sim_markov("violet")`, which draws patients as rows sampled with
## replacement and with probabilities proportional to `n`. Drawing the three
## characteristics together keeps their correlation. The table is internal
## because users have no need for it beyond the preset (see dev/decisions.md).
##
## Provenance. The original VIOLET data are not public. Frank Harrell
## simulated a 250,000-patient trial (`simlongord`) from a Markov model fitted
## to VIOLET. Each simulated patient's starting state, age, and SOFA score were
## drawn with replacement from the 1,351 VIOLET patients with a day-1 record,
## so each VIOLET patient was drawn about 185 times, and dividing how often a
## combination occurs in `simlongord` by 185 recovers how many VIOLET patients
## had it. This script performs that recovery, collapsed to one row per
## distinct combination. The recovered counts reproduce every summary that
## Harrell reported for these patients, including the 474 patients on a
## ventilator or with ARDS, the number of patients with each SOFA score, and
## the mean and percentiles of age; the script stops if any disagrees.
##
## Credit: Frank Harrell created `simlongord`, which is available through
## `Hmisc::getHdata(simlongord)`.
##
## References:
## National Heart, Lung, and Blood Institute PETAL Clinical Trials Network
## (2019). "Early High-Dose Vitamin D3 for Critically Ill, Vitamin D-Deficient
## Patients." New England Journal of Medicine, 381(26), 2529-2540.
## doi:10.1056/NEJMoa1911124
## Harrell FE Jr (2022). "Simulating a Large Clinical Trial From VIOLET."
## https://hbiostat.org/data/repo/simlongord

library(Hmisc)

getHdata(simlongord)

n_violet <- 1351L
baseline <- simlongord[simlongord$time == 1L, c("yprev", "age", "sofa")]
stopifnot(nrow(baseline) == 250000L)

cells <- stats::aggregate(
  list(count = rep(1L, nrow(baseline))),
  by = list(
    state = as.character(baseline$yprev),
    age = baseline$age,
    sofa = baseline$sofa
  ),
  FUN = sum
)

# A VIOLET patient is expected to be drawn 250,000 / 1,351 = 185 times, so a
# combination's count divided by 185 estimates how many VIOLET patients had it.
cells$expected <- cells$count * n_violet / nrow(baseline)
cells$n <- as.integer(round(cells$expected))

# Counts published by Harrell for the 1,351 VIOLET patients.
published_state <- c("Vent/ARDS" = 474L, "In Hospital/Facility" = 877L)
published_sofa <- c(
  87L, 123L, 112L, 136L, 158L, 134L, 135L, 106L, 89L, 76L, 65L, 54L, 32L, 17L,
  11L, 6L, 5L, 3L, 2L
)
names(published_sofa) <- 0:18

margin_excess <- function(n) {
  state <- tapply(n, factor(cells$state, names(published_state)), sum)
  sofa <- tapply(n, factor(cells$sofa, names(published_sofa)), sum)
  list(
    state = state - published_state,
    sofa = sofa - published_sofa
  )
}

loglik <- function(n) {
  sum(stats::dbinom(cells$count, nrow(baseline), n / sum(n), log = TRUE))
}

# Rounding leaves a combination that lies about halfway between two
# multiplicities. Remove a patient from, or add one to, the combination whose
# starting state and SOFA score both disagree with the published counts,
# choosing the change that best explains the simulated counts.
repeat {
  excess <- margin_excess(cells$n)
  if (all(excess$state == 0L) && all(excess$sofa == 0L)) {
    break
  }
  direction <- if (sum(cells$n) > n_violet) -1L else 1L
  eligible <- which(
    excess$state[cells$state] * direction < 0 &
      excess$sofa[as.character(cells$sofa)] * direction < 0 &
      cells$n + direction >= 0L
  )
  stopifnot(length(eligible) > 0L)
  candidate_loglik <- vapply(
    eligible,
    function(i) {
      n <- cells$n
      n[i] <- n[i] + direction
      loglik(n)
    },
    numeric(1)
  )
  best <- eligible[which.max(candidate_loglik)]
  cells$n[best] <- cells$n[best] + direction
}

# Check the recovered records against Harrell's published description of the
# 1,351 VIOLET patients.
ages <- rep(cells$age, cells$n)
sofa <- rep(cells$sofa, cells$n)
stopifnot(
  sum(cells$n) == n_violet,
  all(cells$n > 0L),
  length(unique(ages)) == 78L,
  round(mean(ages), 2) == 56.31,
  round(Hmisc::GiniMd(ages), 2) == 18.16,
  all(
    stats::quantile(ages, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)) ==
      c(27, 33, 46, 58, 67.5, 76, 82)
  ),
  round(mean(sofa), 3) == 5.377,
  round(Hmisc::GiniMd(sofa), 3) == 4.029
)

# mostr numbers states from best to worst: 1 = home, 2 = in hospital or
# facility, 3 = on a ventilator or with ARDS, 4 = dead.
state_code <- c("In Hospital/Facility" = 2L, "Vent/ARDS" = 3L)
violet_baseline <- data.frame(
  yprev = unname(state_code[cells$state]),
  age = as.integer(cells$age),
  sofa = as.integer(cells$sofa),
  n = cells$n
)
violet_baseline <- violet_baseline[
  order(violet_baseline$yprev, violet_baseline$age, violet_baseline$sofa),
]
rownames(violet_baseline) <- NULL

usethis::use_data(
  violet_baseline,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
