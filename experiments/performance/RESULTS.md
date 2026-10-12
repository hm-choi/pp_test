# Performance

Script `run_performance.py`; raw results `results/performance_raw.csv`.

HEaaN FGb (CPU), log_slots = 15, alpha = 0.05, 10 repetitions. Mean (standard deviation) of the
time (s) and bootstraps; encryption and decryption are not timed. Statistic = statistics part (Welch: Algorithms 2-3);
ApproxTCrit = the critical-value polynomial alone on the same encrypted 1/df (part of the Welch decision); the rest of
the decision is mainly the step function. The decrypted decisions match the plaintext decisions in 210/210 runs.

| case | Welch total / BTS | Welch statistic | ApproxTCrit | F total / BTS | Z total / BTS |
|---|---|---|---|---|---|
| insurance_charges_smoker | 23.39 (1.07) / 8 | 11.26 (0.44) | 2.27 (0.22) | 11.47 (0.44) / 4 | 10.88 (0.25) / 4 |
| adult_edu_income | 22.35 (1.15) / 8 | 10.72 (0.45) | 2.10 (0.13) | 11.32 (0.41) / 4 | 10.85 (0.37) / 4 |
| adult_age_gender | 22.01 (1.31) / 8 | 10.56 (0.66) | 2.07 (0.19) | 11.15 (0.37) / 4 | 10.61 (0.41) / 4 |
| heart_chol_disease | 23.12 (1.40) / 8 | 10.93 (0.75) | 2.21 (0.15) | 13.84 (0.79) / 5 | 10.84 (0.56) / 4 |
| diabetes_labs_readmitted | 25.10 (1.95) / 9 | 11.08 (0.98) | 2.21 (0.21) | 13.96 (0.87) / 5 | 10.66 (0.43) / 4 |
| credit_limit_default | 24.86 (1.33) / 9 | 11.22 (0.64) | 2.08 (0.16) | 13.23 (0.87) / 5 | 10.70 (0.59) / 4 |
| bank_age_subscribed | 24.51 (1.29) / 9 | 10.96 (0.53) | 2.18 (0.19) | 13.29 (0.83) / 5 | 10.82 (0.47) / 4 |
