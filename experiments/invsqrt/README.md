# Inverse square root: parameter search

Section 2.1.1 of the paper. Optimizer used to select the inverse square root configuration of PP-TEST: a degree-63 Chebyshev initial value followed by 7 Newton iterations on V in [0.01, 1000].

## Method

Port of the HE-DAP parameter search (Algorithms 2 and 3 of Park et al., SAC'26) to `HEApprox.invSqrt` on HEaaN.

- Chebyshev degree 2^d - 1 with d = 4 to 9 (stored coefficients `coeffs/invSqrt_coeffs.py`), every input level from l_BTS + 1 to l_max, with and without a bootstrap of the input before the Chebyshev evaluation, up to 15 Newton iterations.
- For every combination, the MRE and the accumulated time after each iteration are recorded. The iteration count and then the configuration per input level are selected as in HE-DAP.

## Run

```
python3 experiments/invsqrt/hedap_optimizer.py
python3 experiments/invsqrt/measure.py       # repeated measurement of the selected configurations
```

## Output

Written to `results/` next to the scripts.

| File | Contents |
| --- | --- |
| `hedap_raw.csv` | MRE, MaxRE, time and bootstrapping count after every Newton iteration, per input level, degree and input bootstrap |
| `hedap_optimal.json` | Selected configuration per input level |
