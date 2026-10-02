# invSqrt: HE-DAP parameter selection

Configuration: HEaaN FGb, CPU (`heaan-stat` container), `log_slots=15`. Search: degrees 2^4-1..2^9-1, input levels 4..12,
Pre-BTS c in {0,1}, i_max=15, theta=delta=1 (HE-DAP, Park et al. SAC'26, Algorithms 2/3). One run per combination.
Domains are defined in `domains.py`: `normalized` (x = V/V_max in [v_min, 1], fixed coefficients) and `raw`
(x = V in [1e-3, V_max] per dataset comparison and test, coefficients generated from the public bounds).

Implementation notes (HEApprox.invSqrt): the Chebyshev output y0 is refreshed with the regular bootstrap (its error for
|y0| > 1 is refined by Newton); y in the Newton steps uses the extended bootstrap. Timing excludes decryption.

Selected configuration `accuracy` = HE-DAP u1 (fastest with MRE <= m_theta). Level 9/10 is the invSqrt input level reached in
`HEHypothesisTesting` (10 for Adult educational-num), level 12 is a freshly bootstrapped input.

## normalized

| domain | lo | hi | level | accuracy deg/c/iter | time (s) | MRE | MaxRE | speed deg/c/iter | time (s) | MRE |
|---|---|---|---|---|---|---|---|---|---|---|
| normalized_vmin0.001 | 0.001 | 1 | 12 | 31/0/3 | 5.44 | 5.79e-08 | 2.30e-05 | 63/0/2 | 3.51 | 6.19e-08 |
| normalized_vmin0.001 | 0.001 | 1 | 10 | 31/0/3 | 7.34 | 5.81e-08 | 2.94e-05 | 63/0/2 | 3.11 | 6.02e-08 |
| normalized_vmin0.001 | 0.001 | 1 | 9 | 255/0/1 | 7.98 | 5.88e-08 | 1.98e-05 | 127/0/1 | 6.37 | 6.31e-08 |
| normalized_vmin0.001 | 0.001 | 1 | 6 | 127/0/1 | 8.42 | 7.88e-08 | 2.85e-05 | 127/1/1 | 6.83 | 1.33e-07 |
| normalized_vmin0.001 | 0.001 | 1 | 4 | 127/1/1 | 6.89 | 1.34e-07 | 4.22e-05 | 127/1/1 | 6.89 | 1.34e-07 |
| normalized_vmin0.0001 | 0.0001 | 1 | 12 | 63/0/4 | 8.10 | 7.91e-08 | 1.86e-04 | 127/0/2 | 4.72 | 8.81e-08 |
| normalized_vmin0.0001 | 0.0001 | 1 | 10 | 63/0/3 | 7.19 | 7.84e-08 | 8.24e-05 | 127/0/2 | 5.92 | 8.79e-08 |
| normalized_vmin0.0001 | 0.0001 | 1 | 9 | 127/0/2 | 6.58 | 8.38e-08 | 1.26e-04 | 127/1/2 | 6.53 | 2.05e-07 |
| normalized_vmin0.0001 | 0.0001 | 1 | 6 | 31/0/5 | 8.76 | 9.61e-08 | 1.00e-04 | 127/1/2 | 5.79 | 1.88e-07 |
| normalized_vmin0.0001 | 0.0001 | 1 | 4 | 127/1/3 | 6.02 | 1.97e-07 | 2.50e-04 | 127/1/3 | 6.02 | 1.97e-07 |
| normalized_vmin1e-05 | 1e-05 | 1 | 12 | 255/0/3 | 6.67 | 1.36e-07 | 9.49e-04 | 255/0/3 | 6.67 | 1.36e-07 |
| normalized_vmin1e-05 | 1e-05 | 1 | 10 | 255/0/3 | 13.19 | 9.32e-08 | 1.98e-04 | 255/1/3 | 8.58 | 3.03e-07 |
| normalized_vmin1e-05 | 1e-05 | 1 | 9 | 255/0/2 | 8.83 | 1.74e-07 | 1.94e-03 | 255/1/3 | 8.29 | 2.93e-07 |
| normalized_vmin1e-05 | 1e-05 | 1 | 6 | 127/0/4 | 8.65 | 1.65e-07 | 1.17e-03 | 255/1/3 | 8.51 | 2.79e-07 |
| normalized_vmin1e-05 | 1e-05 | 1 | 4 | 511/1/2 | 12.33 | 2.30e-07 | 1.01e-03 | 255/1/3 | 8.26 | 4.12e-07 |

## raw

| domain | lo | hi | level | accuracy deg/c/iter | time (s) | MRE | MaxRE | speed deg/c/iter | time (s) | MRE |
|---|---|---|---|---|---|---|---|---|---|---|
| raw_insurance_charges_smoker_welch | 0.001 | 4.68 | 12 | 63/0/3 | 3.30 | 1.63e-08 | 1.84e-05 | 63/0/3 | 3.30 | 1.63e-08 |
| raw_insurance_charges_smoker_welch | 0.001 | 4.68 | 10 | 127/0/2 | 6.05 | 1.52e-08 | 9.65e-06 | 127/0/2 | 6.05 | 1.52e-08 |
| raw_insurance_charges_smoker_welch | 0.001 | 4.68 | 9 | 127/0/2 | 6.31 | 1.63e-08 | 1.65e-05 | 127/0/2 | 6.31 | 1.63e-08 |
| raw_insurance_charges_smoker_welch | 0.001 | 4.68 | 6 | 127/0/2 | 11.01 | 1.61e-08 | 1.34e-05 | 127/1/2 | 6.84 | 1.62e-07 |
| raw_insurance_charges_smoker_welch | 0.001 | 4.68 | 4 | 127/1/2 | 6.08 | 1.60e-07 | 1.59e-04 | 127/1/2 | 6.08 | 1.60e-07 |
| raw_insurance_charges_smoker_f | 0.001 | 1018 | 12 | 255/0/6 | 10.09 | 3.07e-08 | 1.75e-05 | 255/0/6 | 10.09 | 3.07e-08 |
| raw_insurance_charges_smoker_f | 0.001 | 1018 | 10 | 511/0/4 | 16.30 | 3.03e-08 | 4.99e-06 | 127/1/6 | 11.69 | 1.26e-06 |
| raw_insurance_charges_smoker_f | 0.001 | 1018 | 9 | 511/0/4 | 16.83 | 3.03e-08 | 1.30e-06 | 511/1/3 | 11.63 | 1.10e-06 |
| raw_insurance_charges_smoker_f | 0.001 | 1018 | 6 | 511/0/8 | 21.04 | 3.00e-08 | 2.14e-06 | 511/1/3 | 12.01 | 1.11e-06 |
| raw_insurance_charges_smoker_f | 0.001 | 1018 | 4 | 127/1/7 | 15.29 | 9.95e-07 | 2.43e-02 | 255/1/4 | 13.28 | 1.57e-06 |
| raw_adult_edu_income_welch | 0.001 | 0.1251 | 12 | 31/0/1 | 0.67 | 3.71e-07 | 3.72e-05 | 31/0/1 | 0.67 | 3.71e-07 |
| raw_adult_edu_income_welch | 0.001 | 0.1251 | 10 | 15/1/2 | 2.45 | 3.76e-07 | 4.36e-05 | 15/1/2 | 2.45 | 3.76e-07 |
| raw_adult_edu_income_welch | 0.001 | 0.1251 | 9 | 15/1/2 | 2.44 | 3.69e-07 | 4.67e-05 | 15/1/2 | 2.44 | 3.69e-07 |
| raw_adult_edu_income_welch | 0.001 | 0.1251 | 6 | 15/1/2 | 2.61 | 3.74e-07 | 3.97e-05 | 15/1/2 | 2.61 | 3.74e-07 |
| raw_adult_edu_income_welch | 0.001 | 0.1251 | 4 | 15/1/2 | 2.45 | 3.70e-07 | 6.03e-05 | 15/1/2 | 2.45 | 3.70e-07 |
| raw_adult_edu_income_f | 0.001 | 64.06 | 12 | 255/0/3 | 6.39 | 3.60e-09 | 1.24e-05 | 255/0/3 | 6.39 | 3.60e-09 |
| raw_adult_edu_income_f | 0.001 | 64.06 | 10 | 255/0/4 | 11.46 | 3.53e-09 | 5.04e-06 | 255/1/2 | 7.71 | 2.85e-07 |
| raw_adult_edu_income_f | 0.001 | 64.06 | 9 | 511/0/2 | 11.21 | 3.51e-09 | 7.60e-06 | 255/1/2 | 8.21 | 2.86e-07 |
| raw_adult_edu_income_f | 0.001 | 64.06 | 6 | 255/0/3 | 12.72 | 3.97e-09 | 1.54e-05 | 255/1/2 | 8.15 | 2.83e-07 |
| raw_adult_edu_income_f | 0.001 | 64.06 | 4 | 255/1/2 | 8.25 | 2.83e-07 | 2.28e-03 | 255/1/2 | 8.25 | 2.83e-07 |
| raw_adult_age_gender_welch | 0.001 | 3.959 | 12 | 63/0/3 | 3.74 | 1.94e-08 | 1.65e-05 | 63/0/3 | 3.74 | 1.94e-08 |
| raw_adult_age_gender_welch | 0.001 | 3.959 | 10 | 63/0/3 | 7.55 | 1.89e-08 | 3.02e-05 | 127/1/2 | 6.27 | 1.56e-07 |
| raw_adult_age_gender_welch | 0.001 | 3.959 | 9 | 127/0/2 | 6.36 | 1.98e-08 | 3.40e-05 | 127/1/2 | 6.22 | 1.61e-07 |
| raw_adult_age_gender_welch | 0.001 | 3.959 | 6 | 31/0/4 | 11.10 | 1.97e-08 | 4.36e-05 | 15/1/5 | 6.71 | 1.94e-07 |
| raw_adult_age_gender_welch | 0.001 | 3.959 | 4 | 127/1/2 | 6.56 | 1.55e-07 | 1.19e-04 | 127/1/2 | 6.56 | 1.55e-07 |
| raw_adult_age_gender_f | 0.001 | 2027 | 12 | 255/0/9 | 14.63 | 5.98e-08 | 1.12e-06 | 255/0/9 | 14.63 | 5.98e-08 |
| raw_adult_age_gender_f | 0.001 | 2027 | 10 | 63/0/10 | 23.33 | 5.96e-08 | 3.95e-06 | 255/1/5 | 12.78 | 2.71e-06 |
| raw_adult_age_gender_f | 0.001 | 2027 | 9 | 511/0/8 | 26.16 | 5.99e-08 | 9.69e-06 | 255/1/5 | 12.32 | 1.96e-06 |
| raw_adult_age_gender_f | 0.001 | 2027 | 6 | 255/0/7 | 18.20 | 5.96e-08 | 5.07e-06 | 255/1/6 | 13.01 | 1.78e-06 |
| raw_adult_age_gender_f | 0.001 | 2027 | 4 | 255/1/6 | 13.45 | 1.93e-06 | 5.44e-02 | 255/1/6 | 13.45 | 1.93e-06 |

Full per-level results: `results/<domain set>/hedap_optimal.json`; raw per-iteration data: `hedap_raw.csv`.
Standalone repeated measurements (`measure.py`) were not rerun after the y0 bootstrap change.
