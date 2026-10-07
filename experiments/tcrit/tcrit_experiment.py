"""Critical-boundary approximation (ApproxTCrit): comparison of methods.

Target: the squared two-sided t critical value c^2 = t_{1-alpha/2, df}^2,
reference SciPy `t.ppf`, error = relative error (MaxRE, MRE).

Plaintext part (--part plain). Every method is evaluated under two df-range
conditions, each with its own evaluation grid (log-spaced df plus every
integer df in the range):
- global : df in [1, 2000];
- case   : df in [min(n1, n2) - 1, n1 + n2 - 2] for each dataset case (the
           Welch-Satterthwaite df always lies in this public range).

Polynomial methods (Chebyshev interpolation, degrees 2^3-1 .. 2^7-1, fitted on
the same range they are evaluated on):
- invdf : input u = 1/df (the global fit is the stored coefficient table,
          coeffs/t_critical_coeffs.py);
- df    : input df.

Lookup tables (degree-independent):
- textbook_floor : textbook df grid (distribution_table/t_distribution_table.csv:
                   1..30, 40, 50, 60, 80, 100, 120, 1000), largest tabulated
                   df <= df;
- textbook_interp: same grid plus df = infinity (normal quantile), linear
                   interpolation of c^2 in 1/df;
- integer_floor / integer_round : every integer df, floor / round of df;
- integer_interp : every integer df (up to 4096, covering every case range),
                   linear interpolation of c^2 in 1/df.

For the global condition, MaxRE is also reported per df band
(1-2, 2-10, 10-30, 30-2000).

HE part (--part he): HEHypothesisTesting._critical_value_from_inv_df on
encrypted 1/df values (df log-spaced in [1, 2000]), `--reps` repetitions per
(degree, alpha): runtime, level consumption, error against the exact quantile
and against the plaintext polynomial.

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/tcrit/tcrit_experiment.py --part plain
    python3 experiments/tcrit/tcrit_experiment.py --part he
"""

import argparse
import csv
import math
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial import Chebyshev
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "experiments"))

from coeffs.t_critical_coeffs import (
    _T_CRITICAL_SQ_COEFFS,
    _T_CRITICAL_U_MIN,
    _T_CRITICAL_U_MAX,
)


ALPHAS = [0.001, 0.01, 0.025, 0.05, 0.1]
LOG_DEGREES = range(3, 8)  # degree = 2^3-1, ..., 2^7-1
DF_MIN, DF_MAX = 1.0, 2000.0
DF_BANDS = [(1.0, 2.0), (2.0, 10.0), (10.0, 30.0), (30.0, 2000.0)]
REPS = 10

RESULT_DIR = Path(__file__).resolve().parent / "results"


def exact_sq(alpha, df):

    return stats.t.ppf(1.0 - alpha / 2.0, df) ** 2


def df_grid(lo, hi, n=20001):

    # Log-spaced df plus every integer df in [lo, hi].
    return np.unique(np.concatenate([
        np.geomspace(lo, hi, n),
        np.arange(math.ceil(lo), math.floor(hi) + 1, dtype=np.float64),
    ]))


def textbook_df_grid():

    header = (PROJECT_ROOT / "distribution_table" / "t_distribution_table.csv").read_text().splitlines()[0]
    return np.array([float(c.removeprefix("df_")) for c in header.split(",")[1:]])


def lookup_floor(table_df, df):

    idx = np.clip(np.searchsorted(table_df, df, side="right") - 1, 0, len(table_df) - 1)
    return table_df[idx]


def lookup_interp(alpha, table_df, df, with_infinity=False):

    # Linear interpolation of c^2 in u = 1/df between tabulated points
    # (np.interp needs increasing x, so interpolate over decreasing df).
    u_table = 1.0 / table_df[::-1]
    c_table = exact_sq(alpha, table_df[::-1])

    if with_infinity:
        # df = infinity: u = 0, c = normal quantile.
        u_table = np.concatenate([[0.0], u_table])
        c_table = np.concatenate([[stats.norm.ppf(1.0 - alpha / 2.0) ** 2], c_table])

    return np.interp(1.0 / df, u_table, c_table)


def errors(y, y_true):

    re = np.abs(y - y_true) / np.abs(y_true)
    return float(re.max()), float(re.mean())


def ranges():

    from pptest_cases import load_cases

    out = {"global": (DF_MIN, DF_MAX)}

    for name, case in load_cases().items():
        out[name] = (float(min(case["n1"], case["n2"]) - 1), float(case["n1"] + case["n2"] - 2))

    return out


def run_plain(out_dir):

    table_df = textbook_df_grid()
    integer_df = np.arange(1.0, 4097.0)
    rows = []

    for range_name, (lo, hi) in ranges().items():
        df = df_grid(lo, hi)
        u = 1.0 / df

        for alpha in ALPHAS:
            y_true = exact_sq(alpha, df)

            def add(method, degree, y):
                max_re, mre = errors(y, y_true)
                row = {"range": range_name, "df_lo": lo, "df_hi": hi, "method": method,
                       "degree": degree, "alpha": alpha, "max_re": max_re, "mre": mre}

                if range_name == "global":
                    for b_lo, b_hi in DF_BANDS:
                        mask = (df >= b_lo) & (df <= b_hi)
                        row[f"max_re_df{b_lo:g}-{b_hi:g}"] = errors(y[mask], y_true[mask])[0]

                rows.append(row)

            # Lookup tables.
            add("textbook_floor", "", exact_sq(alpha, lookup_floor(table_df, df)))
            add("textbook_interp", "", lookup_interp(alpha, table_df, df, with_infinity=True))
            add("integer_floor", "", exact_sq(alpha, np.floor(df)))
            add("integer_round", "", exact_sq(alpha, np.maximum(1.0, np.round(df))))
            add("integer_interp", "", lookup_interp(alpha, integer_df, df))

            for ld in LOG_DEGREES:
                degree = 2 ** ld - 1

                # 1/df input.
                if range_name == "global":
                    z = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
                    y_inv = np.polynomial.chebyshev.chebval(z, _T_CRITICAL_SQ_COEFFS[degree][alpha])
                else:
                    p = Chebyshev.interpolate(lambda v: exact_sq(alpha, 1.0 / v), degree, domain=[1.0 / hi, 1.0 / lo])
                    y_inv = p(u)

                add("invdf", degree, y_inv)

                # df input.
                p = Chebyshev.interpolate(lambda d: exact_sq(alpha, d), degree, domain=[lo, hi])
                add("df", degree, p(df))

    fields = sorted({k for r in rows for k in r}, key=lambda k: (
        ["range", "df_lo", "df_hi", "method", "degree", "alpha", "max_re", "mre"] + [f"max_re_df{a:g}-{b:g}" for a, b in DF_BANDS]
    ).index(k))
    path = out_dir / "plaintext.csv"

    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved {path} ({len(rows)} rows)")


def run_he(out_dir, reps):

    from engine.HEengine import HEengine
    from engine.HEdata import Message
    from operators.HEhypo_test import HEHypothesisTesting

    engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")
    ht = HEHypothesisTesting(engine)

    df = np.geomspace(DF_MIN, DF_MAX, engine.num_slots)
    u = 1.0 / df

    fields = ["degree", "alpha", "rep", "time", "input_level", "output_level",
              "bootstrap_count", "max_re", "mre", "max_ae_vs_plain"]
    rows = []

    for ld in LOG_DEGREES:
        degree = 2 ** ld - 1

        for alpha in ALPHAS:
            y_true = exact_sq(alpha, df)
            z = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
            y_plain = np.polynomial.chebyshev.chebval(z, _T_CRITICAL_SQ_COEFFS[degree][alpha])

            for rep in range(reps):
                ctxt = engine.enc(Message(u, engine.log_slots), level=12)

                ht.reset_bootstrap_count()
                start = time.perf_counter()
                out = ht._critical_value_from_inv_df(ctxt, alpha, ld)
                elapsed = time.perf_counter() - start

                y = np.array(engine.dec(out)[0], dtype=np.complex128).real[:len(u)]
                max_re, mre = errors(y, y_true)
                rows.append({
                    "degree": degree, "alpha": alpha, "rep": rep,
                    "time": elapsed, "input_level": ctxt.level(), "output_level": out.level(),
                    "bootstrap_count": ht.bootstrap_count(), "max_re": max_re, "mre": mre,
                    "max_ae_vs_plain": float(np.max(np.abs(y - y_plain))),
                })

            g = rows[-reps:]
            t = np.array([r["time"] for r in g])
            print(f"deg={degree:3d} alpha={alpha:<6} time={t.mean():.3f}+-{t.std(ddof=1):.3f}s "
                  f"levels {g[0]['input_level']}->{g[0]['output_level']} "
                  f"MaxRE={max(r['max_re'] for r in g):.2e} vs plain={max(r['max_ae_vs_plain'] for r in g):.2e}",
                  flush=True)

    path = out_dir / "he_raw.csv"

    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved {path}")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--part", choices=["plain", "he"], required=True)
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.part == "plain":
        run_plain(args.out_dir)
    else:
        run_he(args.out_dir, args.reps)


if __name__ == "__main__":
    main()
