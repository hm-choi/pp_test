# Performance

Script `run_performance.py`; raw results `results/performance_raw.csv`.

HEaaN FGb (CPU), log_slots = 15, alpha = 0.05, 2 repetitions. Mean (standard deviation) of the
time (s) and bootstraps; encryption and decryption are not timed. Statistic = statistics part (Welch: Algorithms 2-3);
ApproxTCrit = the critical-value polynomial alone on the same encrypted 1/df (part of the Welch decision); the rest of
the decision is mainly the step function. The decrypted decisions match the plaintext decisions in 42/42 runs.

| case | Welch total / BTS | Welch statistic | ApproxTCrit | F total / BTS | Z total / BTS |
|---|---|---|---|---|---|
| insurance_charges_smoker | 23.42 (0.72) / 8 | 11.57 (0.17) | 2.19 (0.21) | 10.67 (0.33) / 4 | 11.06 (0.16) / 4 |
| adult_edu_income | 21.46 (0.31) / 8 | 10.16 (0.13) | 2.10 (0.01) | 10.85 (0.21) / 4 | 10.33 (0.25) / 4 |
| adult_age_gender | 24.00 (0.31) / 8 | 11.29 (0.33) | 2.41 (0.10) | 10.48 (0.07) / 4 | 9.90 (0.05) / 4 |
| heart_chol_disease | 21.05 (0.08) / 8 | 9.96 (0.07) | 2.13 (0.02) | 12.93 (0.05) / 5 | 9.86 (0.14) / 4 |
| diabetes_labs_readmitted | 23.00 (0.17) / 9 | 9.98 (0.13) | 2.16 (0.03) | 12.71 (0.49) / 5 | 10.15 (0.71) / 4 |
| credit_limit_default | 23.90 (0.31) / 9 | 10.31 (0.33) | 2.07 (0.07) | 12.39 (0.01) / 5 | 9.86 (0.59) / 4 |
| bank_age_subscribed | 23.32 (0.06) / 9 | 10.37 (0.20) | 2.10 (0.14) | 12.22 (0.09) / 5 | 10.14 (0.68) / 4 |
