# Experiment 1: Accuracy of PP-TEST

Section 5.4 and Table 2 of the paper. The decision agreement reported in Section 5.7 also comes from this experiment.

## Task

Compare the decrypted intermediate values of the encrypted tests with their plaintext (SciPy) counterparts, and the encrypted decision with the plaintext decision.

## Setting

- Six workloads (see the main README), Welch's t-test, the Z-test and the F-test, two-sided significance levels 0.01, 0.025, 0.05 and 0.1, ten repetitions.
- For every workload, test and repetition, the encrypted statistics are computed once and every significance level is evaluated as a separate decision.
- The Z-test uses the sample variance of each group as its public population variance.
- The critical value of Welch's t-test is compared as c^2. The critical values of the F-test and the Z-test are public.

## Results (V100, ten repetitions)

Relative error, mean (standard deviation) over every workload, significance level and repetition.

| Test | Statistic RE | 1/df RE | Critical RE | Score RE |
| --- | --- | --- | --- | --- |
| Welch's t | 1.36e-07 (2.93e-07) | 2.13e-06 (2.98e-06) | 1.19e-07 (4.35e-08) | 3.91e-06 (6.20e-06) |
| F | -- | -- | Public | 2.65e-06 (4.97e-06) |
| Z | 6.50e-09 (1.19e-08) | -- | Public | 1.71e-08 (3.22e-08) |

Decision agreement (step output > 0.5 against the plaintext decision): 72 configurations (six workloads, four significance levels, three tests), all 720 runs match.

| Test | Configurations | Runs | Match | Reject / fail to reject |
| --- | --- | --- | --- | --- |
| Welch's t | 24 | 240 | 240 | 24 / 0 |
| F | 24 | 240 | 240 | 21 / 3 |
| Z | 24 | 240 | 240 | 24 / 0 |

## Run

```
python3 experiments/experiment1/run_accuracy.py --reps 10
```

## Output

| File | Contents |
| --- | --- |
| `accuracy_raw.csv` | One row per workload, test, repetition and significance level: plaintext and decrypted statistic, 1/df, critical value, score and p-value, sign and step outputs, and decisions |
| `accuracy_summary.csv` | Per workload, test and significance level: maximum errors over the repetitions and decision match rates |

The results of the paper are in `results/experiment1/`.
