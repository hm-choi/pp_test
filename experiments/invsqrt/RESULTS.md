# invSqrt: HE-DAP parameter selection

HEaaN FGb (CPU), log_slots = 15. invSqrt input x = V / V_max (normalized Welch variance) in [V_BAR_MIN, 1] = [0.01, 1],
one domain for every dataset comparison. Chebyshev initial value (coefficients generated from the domain) + Newton
iterations. Search (HE-DAP, Park et al. SAC'26, Algorithms 2/3): degree 2^d - 1 (d = 4..9), input level 4..12, Pre-BTS
c in {0, 1}, i_max = 15, theta = delta = 1, one run per combination. Selection `accuracy` = fastest configuration with
MRE <= m_theta (used by `HEHypothesisTesting` as `select`), `speed` = fastest configuration. Script `hedap_optimizer.py`,
raw per-iteration data `results/hedap_raw.csv`, full selection `results/hedap_optimal.json`.

## Selection per input level

m_theta = 4.00e-08.

| level | accuracy: degree / c / iterations | time (s) | BTS | MRE | MaxRE | speed: degree / c / iterations | time (s) | MRE |
|---|---|---|---|---|---|---|---|---|
| 12 | 63 / 0 / 1 | 3.78 | 1 | 3.96e-08 | 2.85e-06 | 63 / 0 / 1 | 3.78 | 3.96e-08 |
| 11 | 127 / 0 / 1 | 4.31 | 1 | 3.98e-08 | 7.15e-06 | 31 / 0 / 2 | 2.66 | 4.00e-08 |
| 10 | 31 / 0 / 2 | 3.13 | 1 | 3.96e-08 | 4.68e-06 | 15 / 0 / 2 | 2.76 | 4.77e-08 |
| 9 | 31 / 0 / 2 | 2.65 | 1 | 3.93e-08 | 3.84e-06 | 31 / 0 / 2 | 2.65 | 3.93e-08 |
| 8 | 63 / 0 / 1 | 3.06 | 1 | 3.97e-08 | 4.56e-06 | 63 / 0 / 1 | 3.06 | 3.97e-08 |
| 7 | 63 / 0 / 1 | 3.40 | 1 | 3.96e-08 | 4.65e-06 | 63 / 1 / 1 | 3.34 | 7.75e-08 |
| 6 | 15 / 0 / 2 | 4.89 | 2 | 5.90e-08 | 7.26e-06 | 63 / 1 / 1 | 3.28 | 7.75e-08 |
| 5 | 63 / 1 / 1 | 3.12 | 1 | 7.72e-08 | 5.14e-06 | 63 / 1 / 1 | 3.12 | 7.72e-08 |
| 4 | 63 / 1 / 1 | 3.61 | 1 | 7.62e-08 | 5.53e-06 | 63 / 1 / 1 | 3.61 | 7.62e-08 |

## Input level in the Welch pipeline

The input level of V / V_max depends on the group sizes: constants below 1e-4 (1/n for n > 10^4) are applied as two
factors and cost one extra level.

| case | V / V_max | input level | degree / c / iterations |
|---|---|---|---|
| insurance_charges_smoker | 0.1151 | 7 | 63 / 0 / 1 |
| adult_edu_income | 0.1111 | 7 | 63 / 0 / 1 |
| adult_age_gender | 0.1448 | 7 | 63 / 0 / 1 |
| heart_chol_disease | 0.0546 | 7 | 63 / 0 / 1 |
| diabetes_labs_readmitted | 0.0874 | 5 | 63 / 1 / 1 |
| credit_limit_default | 0.0579 | 5 | 63 / 1 / 1 |
| bank_age_subscribed | 0.1167 | 5 | 63 / 1 / 1 |
