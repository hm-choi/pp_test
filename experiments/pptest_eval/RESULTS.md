# Encrypted hypothesis tests: accuracy and performance

HEaaN FGb (CPU container `heaan-stat:1.0.0-cpu`, Intel Xeon Icelake, 16 vCPU), log_slots = 15. Two-sided alpha in
{0.001, 0.01, 0.025, 0.05, 0.1}, 10 repetitions. Statistics are computed once per repetition and every alpha is a separate
decision branch; total time = statistics + one decision branch (encryption and decryption excluded). Script `run_eval.py`,
raw results `results/eval_raw.csv`, summary `results/eval_summary.csv`.

## Setup

| case | groups (variable) | n1 / n2 | R | V | Welch df | ciphertexts |
|---|---|---|---|---|---|---|
| insurance_charges_smoker | Insurance charges/1000: smoker yes vs no | 274 / 1064 | 63.7704 | 0.5199 | 312 | 2 |
| adult_edu_income | Adult educational-num: >50K vs <=50K | 1024 / 1024 | 16 | 0.01221 | 2038 | 2 |
| adult_age_gender | Adult age: Female vs Male | 1024 / 1024 | 90 | 0.3773 | 2035 | 2 |
| heart_chol_disease | Heart Disease (Cleveland) cholesterol: disease vs no disease | 139 / 164 | 564 | 35.04 | 299 | 2 |
| diabetes_labs_readmitted | Diabetes 130 num_lab_procedures - 1: readmitted <30 days vs not readmitted | 11357 / 54864 | 131 | 0.03986 | 16691 | 3 |
| credit_limit_default | Credit card LIMIT_BAL/10^4: default vs no default | 6636 / 23364 | 100 | 0.02748 | 11982 | 2 |
| bank_age_subscribed | Bank Marketing age: subscribed vs not subscribed | 5289 / 39922 | 95 | 0.03704 | 6109 | 3 |
| wine_alcohol_color | Wine Quality alcohol x 10: red vs white | 1599 / 4898 | 149 | 0.1019 | 3100 | 2 |

Adult groups use the first 1,024 records of each group; the other datasets use every record. Data lie in [0, R]; R (public)
is the maximum value. Sources and preprocessing: `experiments/pptest_cases.py`.

Decision scores (positive = reject) and public bounds B (the sign approximation receives s / B):
- Welch: s = (mean1 - mean2)^2 - c^2 V, B = max(R^2, c_max^2 V_max); c^2 from the degree-15 polynomial of u = 1/df,
  c_max^2 at df = min(n1, n2) - 1; V_max = R^2/4 (1/(n1-1) + 1/(n2-1)). invSqrt of V (for 1/df) on [1e-3, V_max].
- F: s = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2), B = max(A, F_L B_v) max(A, F_U B_v), A = n1/(n1-1) R^2/4, B_v = n2/(n2-1) R^2/4.
- Z: s = Z^2 - z^2, B = max(R^2 / V_public, z^2); the sample variances are used as the public variances.

Welch invSqrt configurations: `default` = degree 63, 7 Newton iterations; `hedap` = HE-DAP selection per input level
(experiments/invsqrt/RESULTS.md).

Plaintext reference (SciPy):

| case | Welch T^2 (p) | F (p) | Z (p) |
|---|---|---|---|
| insurance_charges_smoker | 1073 (5.89e-103) | 3.708 (3.64e-52) | 32.75 (2.85e-235) |
| adult_edu_income | 323.4 (3.20e-67) | 0.8807 (4.24e-02) | 17.98 (2.68e-72) |
| adult_age_gender | 11.44 (7.34e-04) | 1.159 (1.83e-02) | -3.382 (7.20e-04) |
| heart_chol_disease | 2.227 (1.37e-01) | 0.857 (3.50e-01) | 1.492 (1.36e-01) |
| diabetes_labs_readmitted | 85.35 (2.80e-20) | 0.9481 (2.85e-04) | 9.238 (2.50e-20) |
| credit_limit_default | 838.2 (3.36e-178) | 0.7683 (5.65e-39) | -28.95 (2.68e-184) |
| bank_age_subscribed | 18.65 (1.60e-05) | 1.761 (1.66e-189) | 4.318 (1.57e-05) |
| wine_alcohol_color | 8.174 (4.28e-03) | 0.7499 (5.95e-12) | -2.859 (4.25e-03) |

Plaintext decisions: reject in 132 and fail to reject in 28 of 160 (case, test, configuration, alpha) combinations.
Fail to reject: adult_age_gender f (alpha 0.001, 0.01); adult_edu_income f (alpha 0.001, 0.01, 0.025); heart_chol_disease f (alpha 0.001, 0.01, 0.025, 0.05, 0.1); heart_chol_disease welch (alpha 0.001, 0.01, 0.025, 0.05, 0.1); heart_chol_disease z (alpha 0.001, 0.01, 0.025, 0.05, 0.1); wine_alcohol_color welch (alpha 0.001); wine_alcohol_color z (alpha 0.001).

## Decision agreement

**160/160 combinations** match the SciPy decision in all 10 repetitions (step output > 0.5 vs. plaintext score > 0; the sign of the decrypted score and the decision from the HE-derived p-value agree as well).

## Intermediate-value accuracy

Maximum over cases, configurations, alphas and repetitions.

| test | statistic rel. err. | 1/df rel. err. | critical value rel. err. | score err. / B | p-value abs. err. | \|s / B\| (min - max) |
|---|---|---|---|---|---|---|
| Welch | 7.5e-07 | 9.0e-03 | 5.8e-04 | 4.9e-07 | 2.8e-08 | 7.0e-06 - 1.37e-01 |
| F | - | - | public | 5.2e-10 | - | 3.8e-06 - 7.74e-03 |
| Z | 3.4e-08 | - | public | 3.7e-10 | 2.5e-09 | 7.1e-06 - 1.37e-01 |

Welch by case (maximum over configurations, alphas and repetitions):

| case | 1/df rel. err. | critical value rel. err. alpha = 0.001 / 0.05 | score err. / B | min \|s / B\| | min margin / error |
|---|---|---|---|---|---|
| insurance_charges_smoker | 8.8e-07 | 4.7e-05 / 2.8e-06 | 1.3e-07 | 1.4e-01 | 1033704 |
| adult_edu_income | 9.0e-03 | 4.8e-04 / 5.5e-06 | 4.9e-07 | 1.5e-02 | 30204 |
| adult_age_gender | 1.4e-05 | 4.7e-04 / 6.0e-06 | 4.8e-07 | 2.7e-05 | 56 |
| heart_chol_disease | 8.7e-07 | 7.3e-05 / 1.6e-04 | 1.8e-07 | 5.5e-05 | 1280 |
| diabetes_labs_readmitted | 5.4e-03 | 5.8e-04 / 6.3e-06 | 3.0e-08 | 1.7e-04 | 5856 |
| credit_limit_default | 8.6e-03 | 5.8e-04 / 6.1e-06 | 3.5e-08 | 2.3e-03 | 65767 |
| bank_age_subscribed | 2.7e-03 | 5.5e-04 / 5.4e-06 | 4.9e-08 | 3.2e-05 | 648 |
| wine_alcohol_color | 1.9e-04 | 5.1e-04 / 2.8e-05 | 5.2e-08 | 7.0e-06 | 238 |

- alpha = 0.001: the critical-value error equals the plaintext error of the degree-15 coefficients near df >= 2000
  (c^2 relative error up to 9.4e-4).
- The reported critical value is decrypted from c^2 / B; for large R (Heart B = 3.2e5) its relative noise is larger
  (1.6e-4), while the score error / B stays <= 1.8e-7.
- 1/df errors of 2e-3 to 9e-3 occur for small V (Adult educational-num, Diabetes, Credit, Bank); df >= 2000 there, so the
  critical value changes by < 1e-6.
- The F test computes no F statistic or p-value (the decision score needs no division).

## Performance

Mean +- standard deviation over 10 repetitions (s) / bootstraps, alpha = 0.05 (times vary by less than 1 s across alphas; all
alphas in `results/eval_summary.csv`).

| case | Welch hedap | Welch default | F | Z | Welch hedap statistics part |
|---|---|---|---|---|---|
| insurance_charges_smoker | 24.8 +- 1.5 / 8 | 31.3 +- 1.2 / 9 | 10.6 +- 0.4 / 4 | 10.5 +- 0.5 / 4 | 12.5 |
| adult_edu_income | 20.0 +- 1.0 / 7 | 31.0 +- 1.5 / 9 | 10.2 +- 0.1 / 4 | 9.8 +- 0.1 / 4 | 8.4 |
| adult_age_gender | 21.8 +- 0.8 / 7 | 29.3 +- 1.0 / 9 | 11.0 +- 0.4 / 4 | 10.2 +- 0.5 / 4 | 9.5 |
| heart_chol_disease | 35.1 +- 1.5 / 10 | 34.3 +- 2.1 / 10 | 10.4 +- 0.6 / 4 | 10.0 +- 0.2 / 4 | 23.0 |
| diabetes_labs_readmitted | 20.0 +- 0.3 / 7 | 45.9 +- 1.8 / 13 | 10.7 +- 0.6 / 4 | 10.1 +- 0.6 / 4 | 8.6 |
| credit_limit_default | 24.4 +- 1.2 / 8 | 45.8 +- 2.6 / 13 | 10.7 +- 0.6 / 4 | 10.0 +- 0.5 / 4 | 12.7 |
| bank_age_subscribed | 24.8 +- 0.9 / 8 | 45.9 +- 1.7 / 13 | 10.7 +- 0.6 / 4 | 10.1 +- 0.5 / 4 | 13.0 |
| wine_alcohol_color | 23.4 +- 1.1 / 8 | 33.6 +- 0.7 / 10 | 10.4 +- 0.4 / 4 | 10.2 +- 0.6 / 4 | 11.8 |

- Welch hedap is faster than default except Heart (V_max = 1064, invSqrt domain ratio 1e6; HE-DAP selects degree 255 with
  6 iterations at input level 9, so both take ~35 s).
- Diabetes, Credit and Bank (large n, small V): the default configuration needs 13 bootstraps (46 s), HE-DAP 7-8 (20-25 s).
  Diabetes and Bank use 3 ciphertexts (n > 32,768).

## Boundary stress test

Script `run_stress.py`, raw results `results/stress_raw.csv`. Design of `experiments/experiment1/test3.py`: two synthetic groups
of 128 samples, x = 10 + 5 linspace(-1, 1, 128) +- gap/2, with the Welch t statistic set to (t_crit + margin) SE, where t_crit is
the two-sided critical value at df = 254 for each alpha in {0.001, 0.01, 0.025, 0.05, 0.1} (t_crit = 3.336, 2.595, 1.969 at
alpha = 0.001, 0.01, 0.05), margin in {+-0.001, +-0.002, +-0.005, +-0.01}. invSqrt `default` (degree 63, 7 iterations),
10 repetitions per (alpha, margin). For every repetition the statistics are computed once; the score is evaluated with the
critical-value polynomial of degree 7, 15, 31, 63 and 127, and the full decision (sign, step) with degree 15.
Public bound B = 234 - 244 (R^2 dominates).

Full decision (degree 15, sign and step): **400/400 match** the plaintext decision (80 per alpha; step = 1.000000 for positive
margins, 0.000000 for negative margins).

Score sign by critical-value degree (80 results per cell; min margin / error = smallest |s| / |HE score error| over the cell):

| alpha | \|s / B\| at margin +-0.001 | degree 7 | degree 15 | degree 31 | degree 63 | degree 127 |
|---|---|---|---|---|---|---|
| 0.001 | 3.6e-06 | 40/80 (0.0x) | 80/80 (2.1x) | 80/80 (90.0x) | 80/80 (87.2x) | 80/80 (87.1x) |
| 0.01 | 2.9e-06 | 40/80 (0.0x) | 80/80 (1191x) | 80/80 (1178x) | 80/80 (1222x) | 80/80 (1415x) |
| 0.025 | 2.5e-06 | 50/80 (0.1x) | 80/80 (2060x) | 80/80 (2407x) | 80/80 (1823x) | 80/80 (1816x) |
| 0.05 | 2.2e-06 | 80/80 (1.7x) | 80/80 (1494x) | 80/80 (1653x) | 80/80 (1641x) | 80/80 (1928x) |
| 0.1 | 1.9e-06 | 80/80 (30.4x) | 80/80 (1592x) | 80/80 (3168x) | 80/80 (1796x) | 80/80 (2468x) |

Maximum critical-value relative error (HE) by degree:

| alpha | degree 7 | 15 | 31 | 63 | 127 |
|---|---|---|---|---|---|
| 0.001 | 1.0e+00 | 1.4e-04 | 3.8e-06 | 3.9e-06 | 3.9e-06 |
| 0.01 | 4.5e-02 | 2.1e-07 | 2.8e-07 | 2.0e-07 | 1.6e-07 |
| 0.025 | 3.0e-03 | 1.5e-07 | 1.2e-07 | 1.0e-07 | 1.6e-07 |
| 0.05 | 3.0e-04 | 1.6e-07 | 1.4e-07 | 1.4e-07 | 9.8e-08 |
| 0.1 | 2.0e-05 | 1.7e-07 | 1.5e-07 | 1.6e-07 | 1.4e-07 |

Statistics 17.8 +- 1.2 s (4 bootstraps), full degree-15 decision 11.3 s (5 bootstraps), total 29.1 +- 1.7 s.
T^2 relative error <= 4.5e-07, 1/df relative error <= 1.1e-05.
Degree 7 fails at every alpha below 0.05 (negative margins). Degree 15 is correct everywhere; its margin over the error is
>1000x except alpha = 0.001 (2.1x, from the degree-15 coefficient error 9.4e-4 for that alpha), where degree 31 gives 87x.
