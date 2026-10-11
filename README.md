# pp_test

PP-TEST: encrypted two-sample hypothesis testing (Welch's t-test, F-test, Z-test) on CKKS with HEaaN. The tests run
entirely over encrypted data and return an encrypted reject / fail-to-reject decision; the critical value of Welch's
t-test is evaluated in the ciphertext by ApproxTCrit, a Chebyshev polynomial of the inverse degrees of freedom.

## Layout

- `engine/`: HEaaN wrapper (`HEengine`, ciphertext / message containers, Chebyshev evaluation).
- `operators/HEapprox.py`: inverse square root (Chebyshev initial value + Newton iterations) and the composite sign
  approximation used as the step function.
- `operators/HEhypo_test.py`: `HEHypothesisTesting` with the Welch, F and Z tests.
- `coeffs/`: coefficient tables (`t_critical_coeffs.py`: squared t critical values on u = 1/df in [0, 1];
  `sign_coeffs.py`) and the notebook that generates them (`approx.ipynb`).
- `experiments/common/`: dataset comparisons (`cases.py`, files in `datasets/`), plaintext reference values and
  helpers shared by the experiments.
- `experiments/invsqrt/`: HE-DAP parameter search for the inverse square root (used by the other experiments).
- `experiments/accuracy/`: intermediate-value accuracy and decision agreement with the plaintext tests.
- `experiments/critical_boundary/`: comparison of critical-value approximation methods (plaintext and HE).
- `experiments/boundary_decision/`: decisions next to the critical value across degrees of freedom.
- `experiments/performance/`: execution time and bootstrapping counts.
- `distribution_table/`: textbook t-distribution table used by the lookup-table baselines.

Each experiment directory has a `RESULTS.md` and a `results/` directory with the raw data.

## Running

The experiments run in the HEaaN CPU container (`cryptolabinc/heaan-stat:1.0.0-cpu`) with the repository mounted as the
working directory. Keys are generated under `/heaan_setting/` on first use. Run one HE experiment at a time.

```
python3 experiments/invsqrt/hedap_optimizer.py                       # HE-DAP selection for the inverse square root
python3 experiments/accuracy/run_accuracy.py                         # accuracy, decision agreement
python3 experiments/critical_boundary/run_critical_boundary.py --part plain
python3 experiments/critical_boundary/run_critical_boundary.py --part he
python3 experiments/boundary_decision/run_boundary_decision.py       # decisions next to the critical value
python3 experiments/performance/run_performance.py                   # time and bootstraps
```

The scripts append to their CSV files and skip combinations that are already complete; `--reps` sets the number
of repetitions.
