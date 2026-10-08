# invSqrt: HE-DAP parameter selection

HEaaN FGb (CPU), log_slots = 15. invSqrt input x = V (Welch variance term) in [1e-3, V_max], V_max = R^2/4 (1/(n1-1) + 1/(n2-1));
Chebyshev initial value (coefficients generated from the domain) + Newton iterations. Search (HE-DAP, Park et al. SAC'26,
Algorithms 2/3): degree 2^d - 1 (d = 4..9), input level 4..12, Pre-BTS c in {0, 1}, i_max = 15, theta = delta = 1, one run per
combination. Selection `accuracy` = fastest configuration with MRE <= m_theta (used by `HEHypothesisTesting` as `select`),
`speed` = fastest configuration. Script `hedap_optimizer.py`, raw per-iteration data `results/hedap_raw.csv`, full selection
`results/hedap_optimal.json`.

## Selection at the input level reached in the Welch pipeline

The level depends on how many factors the mapping constants are split into (constants below 1e-4 cost one extra level).

| case | V_max | input level | degree / Pre-BTS / iterations | invSqrt time (s) | MRE | MaxRE |
|---|---|---|---|---|---|---|
| insurance_charges_smoker | 4.68 | 10 | 127 / 0 / 2 | 6.05 | 1.5e-08 | 9.6e-06 |
| adult_edu_income | 0.1251 | 10 | 15 / 1 / 2 | 2.45 | 3.8e-07 | 4.4e-05 |
| adult_age_gender | 3.959 | 10 | 63 / 0 / 3 | 7.55 | 1.9e-08 | 3.0e-05 |
| heart_chol_disease | 1064 | 9 | 255 / 0 / 6 | 16.95 | 3.2e-08 | 3.6e-06 |
| diabetes_labs_readmitted | 0.456 | 8 | 63 / 0 / 1 | 3.27 | 1.2e-07 | 2.9e-05 |
| credit_limit_default | 0.4838 | 9 | 31 / 0 / 2 | 2.81 | 1.2e-07 | 4.8e-05 |
| bank_age_subscribed | 0.4832 | 9 | 31 / 0 / 2 | 3.03 | 1.1e-07 | 3.4e-05 |
| wine_alcohol_color | 4.607 | 9 | 127 / 0 / 2 | 6.37 | 1.5e-08 | 1.5e-05 |

## All levels

| case | level | accuracy: degree / c / iterations | time (s) | MRE | MaxRE | speed: degree / c / iterations | time (s) | MRE |
|---|---|---|---|---|---|---|---|---|
| insurance_charges_smoker | 12 | 63 / 0 / 3 | 3.30 | 1.6e-08 | 1.8e-05 | 63 / 0 / 3 | 3.30 | 1.6e-08 |
| insurance_charges_smoker | 11 | 63 / 0 / 3 | 3.47 | 1.8e-08 | 3.0e-05 | 63 / 0 / 3 | 3.47 | 1.8e-08 |
| insurance_charges_smoker | 10 | 127 / 0 / 2 | 6.05 | 1.5e-08 | 9.6e-06 | 127 / 0 / 2 | 6.05 | 1.5e-08 |
| insurance_charges_smoker | 9 | 127 / 0 / 2 | 6.31 | 1.6e-08 | 1.6e-05 | 127 / 0 / 2 | 6.31 | 1.6e-08 |
| insurance_charges_smoker | 8 | 255 / 0 / 1 | 7.27 | 1.7e-08 | 1.8e-05 | 127 / 1 / 2 | 6.50 | 1.6e-07 |
| insurance_charges_smoker | 7 | 255 / 0 / 1 | 8.46 | 1.6e-08 | 1.3e-05 | 127 / 1 / 2 | 6.61 | 1.6e-07 |
| insurance_charges_smoker | 6 | 127 / 0 / 2 | 11.01 | 1.6e-08 | 1.3e-05 | 127 / 1 / 2 | 6.84 | 1.6e-07 |
| insurance_charges_smoker | 5 | 127 / 1 / 2 | 6.80 | 1.6e-07 | 1.2e-04 | 127 / 1 / 2 | 6.80 | 1.6e-07 |
| insurance_charges_smoker | 4 | 127 / 1 / 2 | 6.08 | 1.6e-07 | 1.6e-04 | 127 / 1 / 2 | 6.08 | 1.6e-07 |
| adult_edu_income | 12 | 31 / 0 / 1 | 0.67 | 3.7e-07 | 3.7e-05 | 31 / 0 / 1 | 0.67 | 3.7e-07 |
| adult_edu_income | 11 | 31 / 0 / 1 | 2.68 | 3.7e-07 | 3.7e-05 | 31 / 0 / 1 | 2.68 | 3.7e-07 |
| adult_edu_income | 10 | 15 / 1 / 2 | 2.45 | 3.8e-07 | 4.4e-05 | 15 / 1 / 2 | 2.45 | 3.8e-07 |
| adult_edu_income | 9 | 15 / 1 / 2 | 2.44 | 3.7e-07 | 4.7e-05 | 15 / 1 / 2 | 2.44 | 3.7e-07 |
| adult_edu_income | 8 | 15 / 1 / 2 | 2.55 | 3.7e-07 | 5.6e-05 | 15 / 1 / 2 | 2.55 | 3.7e-07 |
| adult_edu_income | 7 | 15 / 1 / 2 | 2.48 | 3.8e-07 | 5.3e-05 | 15 / 1 / 2 | 2.48 | 3.8e-07 |
| adult_edu_income | 6 | 15 / 1 / 2 | 2.61 | 3.7e-07 | 4.0e-05 | 15 / 1 / 2 | 2.61 | 3.7e-07 |
| adult_edu_income | 5 | 31 / 1 / 1 | 2.66 | 3.9e-07 | 4.1e-05 | 31 / 1 / 1 | 2.66 | 3.9e-07 |
| adult_edu_income | 4 | 15 / 1 / 2 | 2.45 | 3.7e-07 | 6.0e-05 | 15 / 1 / 2 | 2.45 | 3.7e-07 |
| adult_age_gender | 12 | 63 / 0 / 3 | 3.74 | 1.9e-08 | 1.7e-05 | 63 / 0 / 3 | 3.74 | 1.9e-08 |
| adult_age_gender | 11 | 63 / 0 / 3 | 3.55 | 1.8e-08 | 2.3e-05 | 63 / 0 / 3 | 3.55 | 1.8e-08 |
| adult_age_gender | 10 | 63 / 0 / 3 | 7.55 | 1.9e-08 | 3.0e-05 | 127 / 1 / 2 | 6.27 | 1.6e-07 |
| adult_age_gender | 9 | 127 / 0 / 2 | 6.36 | 2.0e-08 | 3.4e-05 | 127 / 1 / 2 | 6.22 | 1.6e-07 |
| adult_age_gender | 8 | 255 / 0 / 1 | 7.91 | 1.8e-08 | 1.7e-05 | 127 / 1 / 2 | 6.32 | 1.6e-07 |
| adult_age_gender | 7 | 255 / 0 / 1 | 8.48 | 1.9e-08 | 1.6e-05 | 127 / 1 / 2 | 6.37 | 1.6e-07 |
| adult_age_gender | 6 | 31 / 0 / 4 | 11.10 | 2.0e-08 | 4.4e-05 | 15 / 1 / 5 | 6.71 | 1.9e-07 |
| adult_age_gender | 5 | 127 / 1 / 2 | 6.42 | 1.6e-07 | 8.8e-05 | 127 / 1 / 2 | 6.42 | 1.6e-07 |
| adult_age_gender | 4 | 127 / 1 / 2 | 6.56 | 1.6e-07 | 1.2e-04 | 127 / 1 / 2 | 6.56 | 1.6e-07 |
| heart_chol_disease | 12 | 255 / 0 / 6 | 10.65 | 3.1e-08 | 1.3e-06 | 255 / 0 / 6 | 10.65 | 3.1e-08 |
| heart_chol_disease | 11 | 63 / 0 / 9 | 13.06 | 3.2e-08 | 1.4e-05 | 63 / 0 / 9 | 13.06 | 3.2e-08 |
| heart_chol_disease | 10 | 511 / 0 / 4 | 17.11 | 3.3e-08 | 2.8e-05 | 511 / 1 / 3 | 11.94 | 1.1e-06 |
| heart_chol_disease | 9 | 255 / 0 / 6 | 16.95 | 3.2e-08 | 3.6e-06 | 511 / 1 / 3 | 11.85 | 1.1e-06 |
| heart_chol_disease | 8 | 511 / 0 / 4 | 24.60 | 4.6e-08 | 3.9e-05 | 511 / 1 / 3 | 11.74 | 1.1e-06 |
| heart_chol_disease | 7 | 511 / 0 / 4 | 24.35 | 4.5e-08 | 3.0e-05 | 511 / 1 / 3 | 11.71 | 1.2e-06 |
| heart_chol_disease | 6 | 511 / 0 / 4 | 16.39 | 3.3e-08 | 4.2e-05 | 127 / 1 / 6 | 11.79 | 1.3e-06 |
| heart_chol_disease | 5 | 511 / 1 / 3 | 12.23 | 1.2e-06 | 3.3e-02 | 511 / 1 / 3 | 12.23 | 1.2e-06 |
| heart_chol_disease | 4 | 511 / 1 / 3 | 12.16 | 1.1e-06 | 2.9e-02 | 511 / 1 / 3 | 12.16 | 1.1e-06 |
| diabetes_labs_readmitted | 12 | 63 / 0 / 1 | 3.13 | 1.2e-07 | 3.2e-05 | 63 / 0 / 1 | 3.13 | 1.2e-07 |
| diabetes_labs_readmitted | 11 | 31 / 0 / 2 | 3.02 | 1.2e-07 | 3.6e-05 | 31 / 0 / 2 | 3.02 | 1.2e-07 |
| diabetes_labs_readmitted | 10 | 63 / 0 / 1 | 2.81 | 1.2e-07 | 3.3e-05 | 63 / 0 / 1 | 2.81 | 1.2e-07 |
| diabetes_labs_readmitted | 9 | 31 / 0 / 2 | 2.75 | 1.2e-07 | 3.7e-05 | 31 / 0 / 2 | 2.75 | 1.2e-07 |
| diabetes_labs_readmitted | 8 | 63 / 0 / 1 | 3.27 | 1.2e-07 | 2.9e-05 | 63 / 0 / 1 | 3.27 | 1.2e-07 |
| diabetes_labs_readmitted | 7 | 63 / 0 / 1 | 3.61 | 1.2e-07 | 1.9e-05 | 63 / 0 / 1 | 3.61 | 1.2e-07 |
| diabetes_labs_readmitted | 6 | 63 / 1 / 1 | 3.29 | 1.7e-07 | 3.9e-05 | 63 / 1 / 1 | 3.29 | 1.7e-07 |
| diabetes_labs_readmitted | 5 | 63 / 1 / 1 | 3.60 | 1.7e-07 | 3.1e-05 | 63 / 1 / 1 | 3.60 | 1.7e-07 |
| diabetes_labs_readmitted | 4 | 63 / 1 / 1 | 3.21 | 1.7e-07 | 3.5e-05 | 63 / 1 / 1 | 3.21 | 1.7e-07 |
| credit_limit_default | 12 | 63 / 0 / 1 | 3.65 | 1.2e-07 | 2.8e-05 | 63 / 0 / 1 | 3.65 | 1.2e-07 |
| credit_limit_default | 11 | 31 / 0 / 2 | 2.92 | 1.1e-07 | 2.9e-05 | 31 / 0 / 2 | 2.92 | 1.1e-07 |
| credit_limit_default | 10 | 31 / 0 / 2 | 2.54 | 1.1e-07 | 4.1e-05 | 31 / 0 / 2 | 2.54 | 1.1e-07 |
| credit_limit_default | 9 | 31 / 0 / 2 | 2.81 | 1.2e-07 | 4.8e-05 | 31 / 0 / 2 | 2.81 | 1.2e-07 |
| credit_limit_default | 8 | 63 / 0 / 1 | 3.15 | 1.2e-07 | 3.4e-05 | 63 / 0 / 1 | 3.15 | 1.2e-07 |
| credit_limit_default | 7 | 63 / 1 / 1 | 3.36 | 1.7e-07 | 3.7e-05 | 63 / 1 / 1 | 3.36 | 1.7e-07 |
| credit_limit_default | 6 | 63 / 1 / 1 | 3.27 | 1.7e-07 | 3.7e-05 | 63 / 1 / 1 | 3.27 | 1.7e-07 |
| credit_limit_default | 5 | 63 / 1 / 1 | 3.25 | 1.8e-07 | 4.5e-05 | 63 / 1 / 1 | 3.25 | 1.8e-07 |
| credit_limit_default | 4 | 63 / 1 / 1 | 3.73 | 1.8e-07 | 7.1e-05 | 63 / 1 / 1 | 3.73 | 1.8e-07 |
| bank_age_subscribed | 12 | 63 / 0 / 1 | 3.76 | 1.2e-07 | 3.3e-05 | 63 / 0 / 1 | 3.76 | 1.2e-07 |
| bank_age_subscribed | 11 | 31 / 0 / 2 | 3.02 | 1.2e-07 | 3.9e-05 | 31 / 0 / 2 | 3.02 | 1.2e-07 |
| bank_age_subscribed | 10 | 63 / 0 / 1 | 2.81 | 1.2e-07 | 4.2e-05 | 63 / 0 / 1 | 2.81 | 1.2e-07 |
| bank_age_subscribed | 9 | 31 / 0 / 2 | 3.03 | 1.1e-07 | 3.4e-05 | 31 / 0 / 2 | 3.03 | 1.1e-07 |
| bank_age_subscribed | 8 | 63 / 1 / 1 | 3.04 | 1.8e-07 | 3.8e-05 | 63 / 1 / 1 | 3.04 | 1.8e-07 |
| bank_age_subscribed | 7 | 63 / 1 / 1 | 3.35 | 1.7e-07 | 6.0e-05 | 63 / 1 / 1 | 3.35 | 1.7e-07 |
| bank_age_subscribed | 6 | 63 / 1 / 1 | 3.53 | 1.7e-07 | 4.6e-05 | 63 / 1 / 1 | 3.53 | 1.7e-07 |
| bank_age_subscribed | 5 | 63 / 1 / 1 | 3.30 | 1.7e-07 | 4.4e-05 | 63 / 1 / 1 | 3.30 | 1.7e-07 |
| bank_age_subscribed | 4 | 63 / 1 / 1 | 3.50 | 1.8e-07 | 5.1e-05 | 63 / 1 / 1 | 3.50 | 1.8e-07 |
| wine_alcohol_color | 12 | 63 / 0 / 3 | 3.25 | 1.7e-08 | 2.4e-05 | 63 / 0 / 3 | 3.25 | 1.7e-08 |
| wine_alcohol_color | 11 | 63 / 0 / 3 | 3.52 | 1.6e-08 | 1.4e-05 | 63 / 0 / 3 | 3.52 | 1.6e-08 |
| wine_alcohol_color | 10 | 127 / 0 / 2 | 5.78 | 1.6e-08 | 1.4e-05 | 127 / 0 / 2 | 5.78 | 1.6e-08 |
| wine_alcohol_color | 9 | 127 / 0 / 2 | 6.37 | 1.5e-08 | 1.5e-05 | 127 / 0 / 2 | 6.37 | 1.5e-08 |
| wine_alcohol_color | 8 | 255 / 0 / 1 | 8.07 | 1.6e-08 | 1.7e-05 | 127 / 1 / 2 | 6.39 | 1.6e-07 |
| wine_alcohol_color | 7 | 255 / 0 / 1 | 8.41 | 1.6e-08 | 2.9e-05 | 127 / 1 / 2 | 7.11 | 1.6e-07 |
| wine_alcohol_color | 6 | 127 / 0 / 2 | 11.22 | 1.7e-08 | 1.1e-05 | 127 / 1 / 2 | 6.33 | 1.6e-07 |
| wine_alcohol_color | 5 | 127 / 1 / 2 | 7.01 | 1.6e-07 | 1.1e-04 | 127 / 1 / 2 | 7.01 | 1.6e-07 |
| wine_alcohol_color | 4 | 127 / 1 / 2 | 6.24 | 1.6e-07 | 1.1e-04 | 127 / 1 / 2 | 6.24 | 1.6e-07 |
