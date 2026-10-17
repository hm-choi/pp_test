# PP-TEST

Implementation and experiment results of *PP-TEST: Lookup-Free Privacy-Preserving Hypothesis Testing Framework using Homomorphic Encryption* (ACM SAC 2027).

PP-TEST evaluates Welch's t-test, the Z-test and the F-test over CKKS ciphertexts and returns an encrypted reject/fail-to-reject decision. For Welch's t-test, the data-dependent critical value is evaluated as a degree-15 Chebyshev polynomial of the inverse degrees of freedom (ApproxTCrit), without table lookup or conversion to MPC.

## Layout

```
engine/        HEaaN context, keys, bootstrapping and ciphertext arithmetic
operators/     Encrypted hypothesis tests (HEhypo_test.py), inverse square root and
               step function (HEapprox.py)
coeffs/        Chebyshev coefficients of the inverse square root, the critical values
               and the step function, and the notebook that generates them
datasets/      Datasets of the six workloads
experiments/   experiment1 to experiment4, one folder per experiment of the paper,
               invsqrt (parameter search of the inverse square root), common
results/       GPU results reported in the paper and their summary
figure/        Script of Figure 2
aggregate.py   Builds results/SUMMARY.md from the result CSVs
```

## Experiments

| Folder | Section of the paper |
| --- | --- |
| `experiments/experiment1` | 5.4 Experiment 1: Accuracy of PP-TEST (Table 2); decision agreement in 5.7 |
| `experiments/experiment2` | 5.5 Experiment 2: Critical-Boundary Approximation (Table 3) |
| `experiments/experiment3` | 5.6 Experiment 3: Boundary-Sensitive Decision Evaluation |
| `experiments/experiment4` | 5.7 Experiment 4: Performance Evaluation (Table 4) |
| `experiments/invsqrt` | 2.1.1 CryptoInvSqrt: selection of the Chebyshev degree and the number of Newton iterations |
| `figure/plot_tcrit.py` | Figure 2 |

Each experiment folder has its own README with the setting and the results.

## Workloads

The dataset files are in `datasets/`. The workloads are defined in `experiments/common/cases.py`.

| Workload | File | Attribute | Group 1 | Group 2 | n1 / n2 | [L, U] |
| --- | --- | --- | --- | --- | --- | --- |
| `insurance_charges_smoker` | `insurance.csv` | charges / 1000 | Smoker | Non-smoker | 274 / 1,064 | [1.12, 63.77] |
| `adult_edu_income` | `adult_dataset.csv` | educational-num | Income > 50K | Income <= 50K | 1,024 / 1,024 | [1, 16] |
| `adult_age_gender` | `adult_dataset.csv` | age | Female | Male | 1,024 / 1,024 | [17, 90] |
| `diabetes_labs_readmitted` | `diabetes_130_us_hospitals.csv` | num_lab_procedures | Readmitted < 30 days | No readmission | 11,357 / 54,864 | [1, 132] |
| `credit_limit_default` | `default_of_credit_card_clients.csv` | LIMIT_BAL / 10^4 | Default | No default | 6,636 / 23,364 | [1, 100] |
| `bank_age_subscribed` | `bank-full.csv` | age | Subscribed | Not subscribed | 5,289 / 39,922 | [18, 95] |

The Adult workloads use the first 1,024 records of each group. [L, U] are the minimum and maximum of the two groups and are treated as public, together with n1 and n2.

## Environment

| | |
| --- | --- |
| HE library | CryptoLab HEaaN.Stat SDK, HEaaN.Stat-gpu v1.0.0, FGb parameter set (32,768 slots) |
| Image | `cryptolabinc/heaan-stat:1.0.0-gpu` |
| Machine | NVIDIA Tesla V100 (32 GB), two Intel Xeon Silver 4210, 256 GB RAM |

## Run

Inside the HEaaN container, from the repository root:

```
export PYTHONPATH=<repo root>
export PPTEST_DEVICE=gpu                    # cpu or gpu (default: cpu)
export PPTEST_HEAAN_SETTING=<key directory> # default: /heaan_setting/
python3 experiments/experiment1/run_accuracy.py --reps 10
```

The commands of every experiment are in its README. The scripts write to `results/` next to the script, so a new run does not overwrite the results reported in the paper.

## Results

`results/experiment1` to `results/experiment4` hold the GPU results reported in the paper (ten repetitions). `python3 aggregate.py` rebuilds `results/SUMMARY.md`, which contains the numbers of the paper and supplementary breakdowns per significance level and per degrees of freedom.

## License

Apache License 2.0 (see `LICENSE`).
