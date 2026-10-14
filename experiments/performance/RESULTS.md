# Performance

Script `run_performance.py`; raw results `results/performance_raw.csv`.

HEaaN FGb (GPU, NVIDIA Tesla V100), log_slots = 15. invSqrt domain V in [0.01, 1000] (degree-63 Chebyshev + 7 Newton iterations); ApproxTCrit degree 15. alpha = 0.05, 10 repetitions. Mean (standard deviation) of the time (s); encryption and decryption
are not timed. Statistic = statistics part of Welch's t-test (Algorithms 1-2); ApproxTCrit = the critical-value
polynomial with its input bootstrap on the same encrypted 1/df (part of the decision).

| case | Welch total | Welch statistic | ApproxTCrit | Welch BTS | F total | F BTS | Z total | Z BTS |
|---|---|---|---|---|---|---|---|---|
| insurance_charges_smoker | 5.62 (1.77) | 3.29 (1.13) | 0.40 (0.14) | 10 | 2.11 (0.62) | 4 | 2.08 (0.45) | 4 |
| adult_edu_income | 5.74 (1.11) | 2.98 (0.73) | 0.53 (0.09) | 10 | 2.19 (0.49) | 4 | 2.06 (0.45) | 4 |
| adult_age_gender | 6.37 (1.35) | 3.56 (0.86) | 0.48 (0.13) | 10 | 2.26 (0.54) | 4 | 2.01 (0.43) | 4 |
| diabetes_labs_readmitted | 5.78 (1.21) | 3.07 (0.72) | 0.47 (0.14) | 10 | 2.82 (0.55) | 5 | 1.92 (0.46) | 4 |
| credit_limit_default | 6.00 (0.99) | 3.03 (0.73) | 0.52 (0.11) | 10 | 2.73 (0.60) | 5 | 2.03 (0.48) | 4 |
| bank_age_subscribed | 5.84 (1.10) | 2.95 (0.67) | 0.53 (0.09) | 10 | 2.60 (0.67) | 5 | 2.34 (0.31) | 4 |

Decrypted decisions match the plaintext decisions in 180/180 runs.
