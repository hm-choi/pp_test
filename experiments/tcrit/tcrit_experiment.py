"""Section 4.4: critical-boundary approximation (ApproxTCrit).

Plaintext part (--part plain):
    Chebyshev approximation of the two-sided t critical value t_{1-alpha/2, df}
    (target "t") and of its square (target "t2") for degrees 2^3-1 .. 2^7-1 and
    alpha in ALPHAS, compared with the exact quantile (scipy, numerical
    integration of the t distribution) over df in [1, 2000]:
    - invdf          : P(1/df) on 1/df in [1/2000, 1] (PP-TEST, the stored
                       coefficients coeffs/t_critical_coeffs.py);
    - df_global      : P(df) on df in [1, 2000];
    - df_group_aware : P(df) on [min(n1, n2) - 1, n1 + n2 - 2] per dataset case;
    - lookup_table   : table lookup at tabulated df values
                       (distribution_table/t_distribution_table.csv grid), using
                       the largest tabulated df <= df (conservative rule);
                       reported over df in [1, 2000], df in [2, 2000] and each
                       dataset case's df range.
    The invdf rows are also evaluated on df in (2000, 4096] (extrapolation
    beyond the approximation domain; Adult cases have Welch df ~2035-2038).

HE part (--part he):
    HEHypothesisTesting._critical_value_from_inv_df on encrypted 1/df values
    (df log-spaced in [1, 2000]), `--reps` repetitions per (degree, alpha,
    target): runtime, level consumption, error against the exact quantile and
    against the plaintext polynomial.

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
    _T_CRITICAL_COEFFS,
    _T_CRITICAL_SQ_COEFFS,
    _T_CRITICAL_U_MIN,
    _T_CRITICAL_U_MAX,
)


ALPHAS = [0.01, 0.025, 0.05, 0.1]
LOG_DEGREES = range(3, 8)  # degree = 2^3-1, ..., 2^7-1
TARGETS = {"t": 1, "t2": 2}
DF_MIN, DF_MAX = 1.0, 2000.0
REPS = 10

RESULT_DIR = Path(__file__).resolve().parent / "results"


def exact(alpha, df, power):

    return stats.t.ppf(1.0 - alpha / 2.0, df) ** power


def df_grid():

    # Dense log grid plus every integer df.
    return np.unique(np.concatenate([
        np.geomspace(DF_MIN, DF_MAX, 20001),
        np.arange(1, int(DF_MAX) + 1, dtype=np.float64),
    ]))


def df_case_grid(lo, hi):

    return np.unique(np.concatenate([np.geomspace(lo, hi, 5001), np.arange(math.ceil(lo), math.floor(hi) + 1)]))


def cheb_on(lo, hi, func, degree):

    return Chebyshev.interpolate(func, degree, domain=[lo, hi])


def table_df_grid():

    header = (PROJECT_ROOT / "distribution_table" / "t_distribution_table.csv").read_text().splitlines()[0]
    return np.array([float(c.removeprefix("df_")) for c in header.split(",")[1:]])


def errors(y, y_true):

    ae = np.abs(y - y_true)
    re = ae / np.abs(y_true)
    return {"max_ae": float(ae.max()), "max_re": float(re.max()), "mre": float(re.mean())}


def run_plain(out_dir):

    from pptest_cases import load_cases

    df = df_grid()
    u = 1.0 / df
    table_df = table_df_grid()
    cases = load_cases()
    rows = []

    for alpha in ALPHAS:
        for target, power in TARGETS.items():
            y_true = exact(alpha, df, power)

            # Lookup table (degree-independent).
            idx = np.searchsorted(table_df, df, side="right") - 1
            y_table = exact(alpha, table_df[idx], power)
            rows.append({"method": "lookup_table", "domain": "df in [1, 2000]", "target": target,
                         "degree": "", "alpha": alpha, "levels": "", **errors(y_table, y_true)})
            mask = df >= 2.0
            rows.append({"method": "lookup_table", "domain": "df in [2, 2000]", "target": target,
                         "degree": "", "alpha": alpha, "levels": "",
                         **errors(y_table[mask], y_true[mask])})

            for name, case in cases.items():
                lo = float(min(case["n1"], case["n2"]) - 1)
                hi = float(case["n1"] + case["n2"] - 2)
                d_case = df_case_grid(lo, hi)
                y_case = exact(alpha, table_df[np.searchsorted(table_df, d_case, side="right") - 1], power)
                rows.append({"method": "lookup_table", "domain": f"{name}: df in [{lo:g}, {hi:g}]",
                             "target": target, "degree": "", "alpha": alpha, "levels": "",
                             **errors(y_case, exact(alpha, d_case, power))})

            for ld in LOG_DEGREES:
                degree = 2 ** ld - 1

                # PP-TEST InvDF with the stored coefficients.
                table = _T_CRITICAL_COEFFS if target == "t" else _T_CRITICAL_SQ_COEFFS
                z = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
                y_inv = np.polynomial.chebyshev.chebval(z, table[degree][alpha])
                rows.append({"method": "invdf", "domain": "1/df in [1/2000, 1]", "target": target,
                             "degree": degree, "alpha": alpha, "levels": ld + 1,
                             **errors(y_inv, y_true)})

                # Extrapolation beyond df = 2000.
                d_ext = np.geomspace(2000.0, 4096.0, 2001)[1:]
                z_ext = 2.0 * (1.0 / d_ext - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
                rows.append({"method": "invdf_extrapolation", "domain": "df in (2000, 4096]", "target": target,
                             "degree": degree, "alpha": alpha, "levels": ld + 1,
                             **errors(np.polynomial.chebyshev.chebval(z_ext, table[degree][alpha]),
                                      exact(alpha, d_ext, power))})

                # Direct DF on the global interval.
                p = cheb_on(DF_MIN, DF_MAX, lambda d: exact(alpha, d, power), degree)
                rows.append({"method": "df_global", "domain": "df in [1, 2000]", "target": target,
                             "degree": degree, "alpha": alpha, "levels": ld + 1,
                             **errors(p(df), y_true)})

                # Direct DF on a group-size-aware interval.
                for name, case in cases.items():
                    lo = float(min(case["n1"], case["n2"]) - 1)
                    hi = float(case["n1"] + case["n2"] - 2)
                    d_case = df_case_grid(lo, hi)
                    p = cheb_on(lo, hi, lambda d: exact(alpha, d, power), degree)
                    rows.append({"method": "df_group_aware", "domain": f"{name}: df in [{lo:g}, {hi:g}]",
                                 "target": target, "degree": degree, "alpha": alpha, "levels": ld + 1,
                                 **errors(p(d_case), exact(alpha, d_case, power))})

    path = out_dir / "plaintext.csv"

    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    for r in rows:
        if r["method"] in ("invdf", "invdf_extrapolation") or (r["method"] == "lookup_table" and "case" not in r["domain"] and ":" not in r["domain"]) \
                or (r["method"] == "df_global" and r["degree"] == 15):
            print(f"{r['method']:19s} {r['domain'][:16]:16s} {r['target']:2s} deg={str(r['degree']):>3} alpha={r['alpha']:<5} "
                  f"MaxAE={r['max_ae']:.2e} MaxRE={r['max_re']:.2e} MRE={r['mre']:.2e}")

    print(f"saved {path}")


def run_he(out_dir, reps):

    from engine.HEengine import HEengine
    from engine.HEdata import Message
    from operators.HEhypo_test import HEHypothesisTesting

    engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")
    ht = HEHypothesisTesting(engine)

    df = np.geomspace(DF_MIN, DF_MAX, engine.num_slots)
    u = 1.0 / df

    fields = ["target", "degree", "alpha", "rep", "time", "input_level", "output_level",
              "bootstrap_count", "max_ae", "max_re", "mre", "max_ae_vs_plain"]
    rows = []

    for target, power in TARGETS.items():
        table = _T_CRITICAL_COEFFS if target == "t" else _T_CRITICAL_SQ_COEFFS

        for ld in LOG_DEGREES:
            degree = 2 ** ld - 1

            for alpha in ALPHAS:
                y_true = exact(alpha, df, power)
                z = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
                y_plain = np.polynomial.chebyshev.chebval(z, table[degree][alpha])

                for rep in range(reps):
                    ctxt = engine.enc(Message(u, engine.log_slots), level=12)

                    ht.reset_bootstrap_count()
                    start = time.perf_counter()
                    out = ht._critical_value_from_inv_df(ctxt, alpha, ld, target)
                    elapsed = time.perf_counter() - start

                    y = np.array(engine.dec(out)[0], dtype=np.complex128).real[:len(u)]
                    rows.append({
                        "target": target, "degree": degree, "alpha": alpha, "rep": rep,
                        "time": elapsed, "input_level": ctxt.level(), "output_level": out.level(),
                        "bootstrap_count": ht.bootstrap_count(), **errors(y, y_true),
                        "max_ae_vs_plain": float(np.max(np.abs(y - y_plain))),
                    })

                g = rows[-reps:]
                t = np.array([r["time"] for r in g])
                print(f"{target:2s} deg={degree:3d} alpha={alpha:<5} time={t.mean():.3f}+-{t.std(ddof=1):.3f}s "
                      f"levels {g[0]['input_level']}->{g[0]['output_level']} "
                      f"MaxAE={max(r['max_ae'] for r in g):.2e} MaxRE={max(r['max_re'] for r in g):.2e} "
                      f"vs plain={max(r['max_ae_vs_plain'] for r in g):.2e}", flush=True)

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
    args = parser.parse_args()

    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    if args.part == "plain":
        run_plain(RESULT_DIR)
    else:
        run_he(RESULT_DIR, args.reps)


if __name__ == "__main__":
    main()
