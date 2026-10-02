# Sections 4.3 / 4.6: encrypted hypothesis tests

HEaaN FGb CPU. Cases (`experiments/pptest_cases.py`): Insurance charges/1000 by smoker (274 vs 1064), Adult educational-num by income
and Adult age by gender (first 1024 per group). Two-sided alpha in {0.01, 0.025, 0.05, 0.1}; Welch critical value from the degree-15
InvDF polynomial (targets t and t^2). 10 repetitions per configuration; statistics are computed once per repetition and every alpha
(and target) is a separate decision branch; total time = statistics + one branch.

invSqrt configurations: `normalized_vmin*` (V/V_max in [v_min,1], HE-DAP u1 per input level), `raw_default` (V in [1e-3, V_max],
degree 63, 7 iterations), `raw_hedap` (HE-DAP u1 for the case/test domain). Score bounds (public constants per case/test):
insurance_charges_smoker/welch: 3000, insurance_charges_smoker/f: 20, insurance_charges_smoker/z: 3000, adult_edu_income/welch: 700, adult_edu_income/f: 0.02, adult_edu_income/z: 700, adult_age_gender/welch: 20, adult_age_gender/f: 0.03, adult_age_gender/z: 20.

Decision agreement with SciPy: **192/192 configurations** match in all 10 reps (step, score sign and HE-derived p-value).
Max statistic relative error 7.8e-07, max p-value error 2.2e-07.

## Time and accuracy at alpha = 0.05 (Welch target t^2)

| case | test | invSqrt | total (s) | statistics (s) | BTS | stat rel. err | p-value err | plain decision | step match |
|---|---|---|---|---|---|---|---|---|---|
| insurance_charges_smoker | welch | normalized_vmin0.001 | 24.8 +- 1.4 | 12.9 | 8 | 1.73e-08 | 7.73e-106 | 1 | 100% |
| insurance_charges_smoker | welch | normalized_vmin0.0001 | 23.9 +- 0.6 | 11.7 | 8 | 2.46e-08 | 6.33e-106 | 1 | 100% |
| insurance_charges_smoker | welch | normalized_vmin1e-05 | 24.8 +- 0.9 | 12.9 | 8 | 1.56e-08 | 1.40e-105 | 1 | 100% |
| insurance_charges_smoker | welch | raw_default | 32.9 +- 0.9 | 21.4 | 10 | 1.36e-08 | 2.79e-107 | 1 | 100% |
| insurance_charges_smoker | welch | raw_hedap | 23.2 +- 0.6 | 11.5 | 8 | 1.63e-08 | 2.77e-107 | 1 | 100% |
| insurance_charges_smoker | f | normalized_vmin0.001 | 21.2 +- 0.3 | 12.0 | 7 | 2.06e-08 | 1.43e-57 | 1 | 100% |
| insurance_charges_smoker | f | normalized_vmin0.0001 | 20.5 +- 0.5 | 11.1 | 7 | 2.13e-08 | 1.48e-57 | 1 | 100% |
| insurance_charges_smoker | f | normalized_vmin1e-05 | 21.9 +- 0.6 | 12.4 | 7 | 4.44e-08 | 3.08e-57 | 1 | 100% |
| insurance_charges_smoker | f | raw_default | 31.1 +- 1.0 | 21.4 | 9 | 1.44e-07 | 9.96e-57 | 1 | 100% |
| insurance_charges_smoker | f | raw_hedap | 30.5 +- 0.6 | 20.8 | 8 | 1.63e-07 | 1.13e-56 | 1 | 100% |
| insurance_charges_smoker | z | public_variance | 10.3 +- 0.3 | 0.5 | 4 | 1.19e-09 | 3.64e-241 | 1 | 100% |
| adult_edu_income | welch | normalized_vmin0.001 | 24.5 +- 1.0 | 12.6 | 8 | 2.89e-08 | 3.58e-70 | 1 | 100% |
| adult_edu_income | welch | normalized_vmin0.0001 | 24.5 +- 0.6 | 12.8 | 8 | 2.19e-08 | 3.40e-70 | 1 | 100% |
| adult_edu_income | welch | normalized_vmin1e-05 | 28.7 +- 1.1 | 17.1 | 9 | 3.11e-08 | 3.36e-70 | 1 | 100% |
| adult_edu_income | welch | raw_default | 33.3 +- 1.2 | 21.5 | 10 | 7.85e-08 | 1.05e-68 | 1 | 100% |
| adult_edu_income | welch | raw_hedap | 19.8 +- 0.7 | 8.1 | 7 | 7.83e-07 | 2.00e-68 | 1 | 100% |
| adult_edu_income | f | normalized_vmin0.001 | 21.6 +- 0.8 | 11.9 | 7 | 3.53e-08 | 5.73e-08 | 1 | 100% |
| adult_edu_income | f | normalized_vmin0.0001 | 21.8 +- 0.6 | 12.4 | 7 | 4.37e-08 | 7.10e-08 | 1 | 100% |
| adult_edu_income | f | normalized_vmin1e-05 | 26.4 +- 0.7 | 16.8 | 8 | 1.98e-08 | 3.22e-08 | 1 | 100% |
| adult_edu_income | f | raw_default | 28.0 +- 0.6 | 16.6 | 9 | 5.08e-08 | 8.25e-08 | 1 | 100% |
| adult_edu_income | f | raw_hedap | 26.1 +- 0.5 | 16.6 | 8 | 5.78e-08 | 9.39e-08 | 1 | 100% |
| adult_edu_income | z | public_variance | 10.1 +- 0.4 | 0.5 | 4 | 6.26e-09 | 5.45e-78 | 1 | 100% |
| adult_age_gender | welch | normalized_vmin0.001 | 25.1 +- 0.9 | 13.0 | 8 | 1.84e-08 | 1.42e-09 | 1 | 100% |
| adult_age_gender | welch | normalized_vmin0.0001 | 23.2 +- 0.7 | 11.6 | 8 | 1.07e-08 | 1.41e-09 | 1 | 100% |
| adult_age_gender | welch | normalized_vmin1e-05 | 24.6 +- 0.8 | 13.0 | 8 | 1.48e-08 | 2.12e-09 | 1 | 100% |
| adult_age_gender | welch | raw_default | 33.3 +- 1.1 | 21.4 | 10 | 8.57e-09 | 4.94e-11 | 1 | 100% |
| adult_age_gender | welch | raw_hedap | 23.9 +- 0.9 | 11.7 | 8 | 2.56e-08 | 1.38e-10 | 1 | 100% |
| adult_age_gender | f | normalized_vmin0.001 | 24.4 +- 1.0 | 12.6 | 8 | 9.16e-09 | 7.21e-09 | 1 | 100% |
| adult_age_gender | f | normalized_vmin0.0001 | 22.4 +- 0.6 | 11.1 | 8 | 1.36e-08 | 1.07e-08 | 1 | 100% |
| adult_age_gender | f | normalized_vmin1e-05 | 24.3 +- 1.1 | 12.5 | 8 | 1.76e-08 | 1.39e-08 | 1 | 100% |
| adult_age_gender | f | raw_default | 30.7 +- 0.7 | 21.2 | 9 | 8.56e-08 | 6.74e-08 | 1 | 100% |
| adult_age_gender | f | raw_hedap | 38.8 +- 1.0 | 29.3 | 10 | 2.73e-07 | 2.15e-07 | 1 | 100% |
| adult_age_gender | z | public_variance | 10.4 +- 0.5 | 0.5 | 4 | 3.84e-09 | 3.40e-11 | 1 | 100% |

All alphas and both Welch targets: `results/eval_summary.csv`; per repetition: `results/eval_raw.csv`.

## Section 4.5: boundary stress test (`run_stress.py`)

Design of `experiments/experiment1/test3.py`: two synthetic groups of 128 samples, Welch t = (t_crit + margin) * SE,
margins +-0.001, +-0.002, +-0.005, +-0.01, alpha = 0.05, score bound 0.05. Critical value + score for degrees 7..127 and both targets;
full decision (sign/step) at degree 15. 10 repetitions per margin and invSqrt configuration.

| degree | target | with sign | runs | score sign match | step match | max critical err | max score err | branch time (s) |
|---|---|---|---|---|---|---|---|---|
| 7 | t | no | 320 | 100% | - | 3.05e-06 | 1.33e-05 | 1.98 |
| 7 | t2 | no | 320 | 100% | - | 5.84e-04 | 2.30e-03 | 1.91 |
| 15 | t | no | 320 | 100% | - | 2.83e-07 | 1.63e-06 | 2.07 |
| 15 | t | yes | 320 | 100% | 100% | 2.83e-07 | 1.63e-06 | 13.33 |
| 15 | t2 | no | 320 | 100% | - | 2.85e-07 | 1.63e-06 | 2.05 |
| 15 | t2 | yes | 320 | 100% | 100% | 2.85e-07 | 1.63e-06 | 13.34 |
| 31 | t | no | 320 | 100% | - | 2.84e-07 | 1.63e-06 | 2.17 |
| 31 | t2 | no | 320 | 100% | - | 2.85e-07 | 1.63e-06 | 2.30 |
| 63 | t | no | 320 | 100% | - | 2.82e-07 | 1.62e-06 | 2.37 |
| 63 | t2 | no | 320 | 100% | - | 2.85e-07 | 1.63e-06 | 2.72 |
| 127 | t | no | 320 | 100% | - | 2.84e-07 | 1.63e-06 | 2.54 |
| 127 | t2 | no | 320 | 100% | - | 2.85e-07 | 1.63e-06 | 3.36 |

All 3840 results match the plaintext decision (all invSqrt configurations: normalized v_min 1e-3/1e-4/1e-5 and raw_default).
