# invSqrt: HE-DAP parameter selection

HEaaN FGb (CPU), log_slots = 15. invSqrt input x = V (Welch variance term) in [1e-3, V_max], V_max = R^2/4 (1/(n1-1) + 1/(n2-1));
Chebyshev initial value (coefficients generated from the domain) + Newton iterations. Search (HE-DAP, Park et al. SAC'26,
Algorithms 2/3): degree 2^d - 1 (d = 4..9), input level 4..12, Pre-BTS c in {0, 1}, i_max = 15, theta = delta = 1, one run per
combination. Selection `accuracy` = fastest configuration with MRE <= m_theta (used by `HEHypothesisTesting` as `select`),
`speed` = fastest configuration. Script `hedap_optimizer.py`, raw per-iteration data `results/hedap_raw.csv`, full selection
`results/hedap_optimal.json`.

In the Welch pipeline the invSqrt input arrives at level 10 (all three cases).

| case | V_max | level | accuracy: degree / c / iterations | time (s) | MRE | MaxRE | speed: degree / c / iterations | time (s) | MRE |
|---|---|---|---|---|---|---|---|---|---|
| insurance_charges_smoker | 4.68 | 12 | 63 / 0 / 3 | 3.30 | 1.6e-08 | 1.8e-05 | 63 / 0 / 3 | 3.30 | 1.6e-08 |
| insurance_charges_smoker | 4.68 | 11 | 63 / 0 / 3 | 3.47 | 1.8e-08 | 3.0e-05 | 63 / 0 / 3 | 3.47 | 1.8e-08 |
| insurance_charges_smoker | 4.68 | 10 | 127 / 0 / 2 | 6.05 | 1.5e-08 | 9.6e-06 | 127 / 0 / 2 | 6.05 | 1.5e-08 |
| insurance_charges_smoker | 4.68 | 9 | 127 / 0 / 2 | 6.31 | 1.6e-08 | 1.6e-05 | 127 / 0 / 2 | 6.31 | 1.6e-08 |
| insurance_charges_smoker | 4.68 | 8 | 255 / 0 / 1 | 7.27 | 1.7e-08 | 1.8e-05 | 127 / 1 / 2 | 6.50 | 1.6e-07 |
| insurance_charges_smoker | 4.68 | 7 | 255 / 0 / 1 | 8.46 | 1.6e-08 | 1.3e-05 | 127 / 1 / 2 | 6.61 | 1.6e-07 |
| insurance_charges_smoker | 4.68 | 6 | 127 / 0 / 2 | 11.01 | 1.6e-08 | 1.3e-05 | 127 / 1 / 2 | 6.84 | 1.6e-07 |
| insurance_charges_smoker | 4.68 | 5 | 127 / 1 / 2 | 6.80 | 1.6e-07 | 1.2e-04 | 127 / 1 / 2 | 6.80 | 1.6e-07 |
| insurance_charges_smoker | 4.68 | 4 | 127 / 1 / 2 | 6.08 | 1.6e-07 | 1.6e-04 | 127 / 1 / 2 | 6.08 | 1.6e-07 |
| adult_edu_income | 0.1251 | 12 | 31 / 0 / 1 | 0.67 | 3.7e-07 | 3.7e-05 | 31 / 0 / 1 | 0.67 | 3.7e-07 |
| adult_edu_income | 0.1251 | 11 | 31 / 0 / 1 | 2.68 | 3.7e-07 | 3.7e-05 | 31 / 0 / 1 | 2.68 | 3.7e-07 |
| adult_edu_income | 0.1251 | 10 | 15 / 1 / 2 | 2.45 | 3.8e-07 | 4.4e-05 | 15 / 1 / 2 | 2.45 | 3.8e-07 |
| adult_edu_income | 0.1251 | 9 | 15 / 1 / 2 | 2.44 | 3.7e-07 | 4.7e-05 | 15 / 1 / 2 | 2.44 | 3.7e-07 |
| adult_edu_income | 0.1251 | 8 | 15 / 1 / 2 | 2.55 | 3.7e-07 | 5.6e-05 | 15 / 1 / 2 | 2.55 | 3.7e-07 |
| adult_edu_income | 0.1251 | 7 | 15 / 1 / 2 | 2.48 | 3.8e-07 | 5.3e-05 | 15 / 1 / 2 | 2.48 | 3.8e-07 |
| adult_edu_income | 0.1251 | 6 | 15 / 1 / 2 | 2.61 | 3.7e-07 | 4.0e-05 | 15 / 1 / 2 | 2.61 | 3.7e-07 |
| adult_edu_income | 0.1251 | 5 | 31 / 1 / 1 | 2.66 | 3.9e-07 | 4.1e-05 | 31 / 1 / 1 | 2.66 | 3.9e-07 |
| adult_edu_income | 0.1251 | 4 | 15 / 1 / 2 | 2.45 | 3.7e-07 | 6.0e-05 | 15 / 1 / 2 | 2.45 | 3.7e-07 |
| adult_age_gender | 3.959 | 12 | 63 / 0 / 3 | 3.74 | 1.9e-08 | 1.7e-05 | 63 / 0 / 3 | 3.74 | 1.9e-08 |
| adult_age_gender | 3.959 | 11 | 63 / 0 / 3 | 3.55 | 1.8e-08 | 2.3e-05 | 63 / 0 / 3 | 3.55 | 1.8e-08 |
| adult_age_gender | 3.959 | 10 | 63 / 0 / 3 | 7.55 | 1.9e-08 | 3.0e-05 | 127 / 1 / 2 | 6.27 | 1.6e-07 |
| adult_age_gender | 3.959 | 9 | 127 / 0 / 2 | 6.36 | 2.0e-08 | 3.4e-05 | 127 / 1 / 2 | 6.22 | 1.6e-07 |
| adult_age_gender | 3.959 | 8 | 255 / 0 / 1 | 7.91 | 1.8e-08 | 1.7e-05 | 127 / 1 / 2 | 6.32 | 1.6e-07 |
| adult_age_gender | 3.959 | 7 | 255 / 0 / 1 | 8.48 | 1.9e-08 | 1.6e-05 | 127 / 1 / 2 | 6.37 | 1.6e-07 |
| adult_age_gender | 3.959 | 6 | 31 / 0 / 4 | 11.10 | 2.0e-08 | 4.4e-05 | 15 / 1 / 5 | 6.71 | 1.9e-07 |
| adult_age_gender | 3.959 | 5 | 127 / 1 / 2 | 6.42 | 1.6e-07 | 8.8e-05 | 127 / 1 / 2 | 6.42 | 1.6e-07 |
| adult_age_gender | 3.959 | 4 | 127 / 1 / 2 | 6.56 | 1.6e-07 | 1.2e-04 | 127 / 1 / 2 | 6.56 | 1.6e-07 |
