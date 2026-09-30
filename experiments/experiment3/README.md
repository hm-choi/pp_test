# Experiment 3 — DF versus InvDF

This experiment compares two ways to obtain the encrypted two-sided
Student-t critical value used by Welch's test.

- `InvDF`: approximate \(c_\alpha(1/\mathrm{df})\) on the universal interval
  \(1/\mathrm{df}\in[0,1]\).
- `DF-global`: approximate \(c_\alpha(\mathrm{df})\) on the fixed interval
  \(\mathrm{df}\in[1,2048]\).
- `DF-group-aware`: approximate \(c_\alpha(\mathrm{df})\) on the tighter interval
  derived from public group sizes. For the insurance comparison this is
  \([273,1336]\).

Run the complete experiment on the HEaaN server:

```bash
python3 experiments/experiment3/test.py
```

Run only the platform-independent approximation sweep:

```bash
python3 experiments/experiment3/test.py --plaintext-only
```

The script writes `approximation_results.csv` and, after an HE run,
`he_results.csv`. The direct-DF HE comparison deliberately uses the tighter
public group-size interval, so it is not disadvantaged by an unnecessarily
wide approximation domain. Its coefficients are consequently specific to
the public group sizes; the InvDF coefficients are reusable.
