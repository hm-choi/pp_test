# Experiment 4: Performance Evaluation

Section 5.7 and Table 4 of the paper. The decision agreement over 72 configurations reported in Section 5.7 comes from Experiment 1.

## Task

Measure the execution time and the number of bootstrapping operations of the encrypted tests.

## Setting

- Six workloads, Welch's t-test, the F-test and the Z-test, alpha = 0.05, ten repetitions.
- Encryption and decryption are not timed. The statistics part and the decision part are timed separately.
- ApproxTCrit is timed once more on the same encrypted 1/df; this time is part of the Welch total.

## Results (V100, alpha = 0.05, ten repetitions)

Time (s), mean (standard deviation). Statistic is the statistics part of Welch's t-test (Algorithms 1 and 2).

| Workload | Welch total | Welch statistic | ApproxTCrit | Welch #BTS | F total | F #BTS | Z total | Z #BTS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `insurance_charges_smoker` | 5.62 (1.77) | 3.29 (1.13) | 0.40 (0.14) | 10 | 2.11 (0.62) | 4 | 2.08 (0.45) | 4 |
| `adult_edu_income` | 5.74 (1.11) | 2.98 (0.73) | 0.53 (0.09) | 10 | 2.19 (0.49) | 4 | 2.06 (0.45) | 4 |
| `adult_age_gender` | 6.37 (1.35) | 3.56 (0.86) | 0.48 (0.13) | 10 | 2.26 (0.54) | 4 | 2.01 (0.43) | 4 |
| `diabetes_labs_readmitted` | 5.78 (1.21) | 3.07 (0.72) | 0.47 (0.14) | 10 | 2.82 (0.55) | 5 | 1.92 (0.46) | 4 |
| `credit_limit_default` | 6.00 (0.99) | 3.03 (0.73) | 0.52 (0.11) | 10 | 2.73 (0.60) | 5 | 2.03 (0.48) | 4 |
| `bank_age_subscribed` | 5.84 (1.10) | 2.95 (0.67) | 0.53 (0.09) | 10 | 2.60 (0.67) | 5 | 2.34 (0.31) | 4 |

For groups larger than 10^4, the constants 1/n and 1/(n - 1) of the mean and the variance are below 10^-4 and are multiplied in two steps, which takes one more level. For the F-test this adds one bootstrapping operation on the Diabetes, Credit Default and Bank Marketing workloads.

## Run

```
python3 experiments/experiment4/run_performance.py --reps 10
```

## Output

| File | Contents |
| --- | --- |
| `performance_raw.csv` | One row per workload, test and repetition: statistics, decision and ApproxTCrit times, bootstrapping counts of both parts, and the decisions |

The results of the paper are in `results/experiment4/`.
