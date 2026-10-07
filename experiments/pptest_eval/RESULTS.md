# Encrypted hypothesis tests: accuracy and performance

HEaaN FGb (CPU container `heaan-stat:1.0.0-cpu`, Intel Xeon Icelake, 16 vCPU), log_slots = 15. Two-sided alpha in
{0.001, 0.01, 0.025, 0.05, 0.1}, 10 repetitions. Statistics are computed once per repetition and every alpha is a separate
decision branch; total time = statistics + one decision branch (encryption and decryption excluded). Script `run_eval.py`,
raw results `results/eval_raw.csv`, summary `results/eval_summary.csv`.

## Setup

| case | groups | n1 / n2 | R |
|---|---|---|---|
| insurance_charges_smoker | Insurance charges/1000, smoker yes vs no | 274 / 1064 | 63.77 |
| adult_edu_income | Adult educational-num, income >50K vs <=50K | 1024 / 1024 | 16 |
| adult_age_gender | Adult age, Female vs Male | 1024 / 1024 | 90 |

Adult groups use the first 1,024 records of each group. Data lie in [0, R]; R (public) is the maximum value.

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
| insurance_charges_smoker | 1073 (5.89e-103) | see note | 32.75 (2.85e-235) |
| adult_edu_income | 323.4 (3.20e-67) | see note | 17.98 (2.68e-72) |
| adult_age_gender | 11.44 (7.34e-04) | see note | -3.382 (7.20e-04) |

F statistics: Insurance 3.708 (p = 3.6e-52), Adult educational-num 0.8807 (p = 0.0424), Adult age 1.159 (p = 0.0183).
Plaintext decisions: reject in 55 of 60 (case, test, configuration, alpha) combinations; fail to reject for F at
Adult educational-num alpha = 0.001, 0.01, 0.025 and Adult age alpha = 0.001, 0.01.

## Decision agreement

**60/60 combinations** match the SciPy decision in all 10 repetitions (step output > 0.5 vs. plaintext score > 0; the sign of the decrypted score and the decision from the HE-derived p-value agree as well).

## Intermediate-value accuracy

Maximum over cases, configurations, alphas and repetitions.

| test | statistic rel. err. | 1/df rel. err. | critical value rel. err. | score err. / B | p-value abs. err. | \|s / B\| (min - max) |
|---|---|---|---|---|---|---|
| Welch | 7.5e-07 | 9.0e-03 | 4.8e-04 | 4.9e-07 | 1.7e-10 | 2.7e-05 - 1.37e-01 |
| F | - | - | public | 4.4e-10 | - | 9.2e-06 - 7.74e-03 |
| Z | 5.9e-09 | - | public | 3.7e-10 | 2.8e-11 | 2.8e-05 - 1.37e-01 |

Welch by alpha (maximum over cases and configurations):

| alpha | critical value rel. err. | score err. / B | min \|s / B\| |
|---|---|---|---|
| 0.001 | 4.8e-04 | 4.9e-07 | 2.7e-05 |
| 0.01 | 8.4e-06 | 5.7e-09 | 2.2e-04 |
| 0.025 | 6.7e-06 | 3.4e-09 | 3.0e-04 |
| 0.05 | 6.0e-06 | 3.3e-09 | 3.5e-04 |
| 0.1 | 1.0e-05 | 2.8e-09 | 4.1e-04 |

1/df relative error by case and configuration: adult_age_gender default 1.4e-05, adult_age_gender hedap 5.9e-06, adult_edu_income default 6.5e-03, adult_edu_income hedap 9.0e-03, insurance_charges_smoker default 8.8e-07, insurance_charges_smoker hedap 6.7e-07.
The F test computes no F statistic or p-value (the decision score needs no division), so those columns are empty.

## Performance

Mean +- standard deviation over 10 repetitions (s) / bootstraps. Times vary by less than 1 s across alphas; alpha = 0.05 shown,
all alphas in `results/eval_summary.csv`.

| case | Welch hedap | Welch default | F | Z | Welch hedap statistics part |
|---|---|---|---|---|---|
| insurance_charges_smoker | 24.8 +- 1.5 / 8 | 31.3 +- 1.2 / 9 | 10.6 +- 0.4 / 4 | 10.5 +- 0.5 / 4 | 12.5 |
| adult_edu_income | 20.0 +- 1.0 / 7 | 31.0 +- 1.5 / 9 | 10.2 +- 0.1 / 4 | 9.8 +- 0.1 / 4 | 8.4 |
| adult_age_gender | 21.8 +- 0.8 / 7 | 29.3 +- 1.0 / 9 | 11.0 +- 0.4 / 4 | 10.2 +- 0.5 / 4 | 9.5 |

Time ranges over all alphas:
- welch hedap: 19.9 - 24.8 s
- welch default: 29.2 - 31.6 s
- f none: 10.2 - 11.2 s
- z none: 9.7 - 10.6 s

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
