# Experiment 3: Boundary-Sensitive Decision Evaluation

Section 5.6 of the paper.

## Task

Check the encrypted Welch decision when the t statistic lies next to the critical value.

## Setting

- Synthetic pairs of groups whose Welch t statistic equals t_{1-alpha/2, df} + margin, with margins +-0.001, +-0.002, +-0.005 and +-0.01.
- Degrees of freedom 3 to 10, 50, 100, 2,000 and 10,000; two-sided significance levels 0.01, 0.025, 0.05 and 0.1; ten repetitions (80 trials per cell).
- Each group has ceil(df / 2) + 1 samples. The variance ratio of the two groups sets the Welch-Satterthwaite df exactly, and both groups are scaled so that V >= 0.02 lies inside the domain of the inverse square root.
- Degree-15 ApproxTCrit and the step function.

## Results (V100, ten repetitions)

All 3,840 encrypted decisions match the plaintext decisions.

| alpha | df 3 | df 4 | df 5 | df 6 | df 7 | df 8 | df 9 | df 10 | df 50 | df 100 | df 2000 | df 10000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.01 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.025 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.05 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.1 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |

## Run

```
python3 experiments/experiment3/run_boundary_decision.py --reps 10
```

## Output

| File | Contents |
| --- | --- |
| `boundary_decision_raw.csv` | One row per df, significance level, margin and repetition: group sizes, V and whether it lies in the domain of the inverse square root, plaintext and decrypted T^2, 1/df, critical value and score, step output, decisions, times and bootstrapping counts |

The results of the paper are in `results/experiment3/`.
