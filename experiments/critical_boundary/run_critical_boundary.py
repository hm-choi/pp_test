"""Critical-boundary approximation (ApproxTCrit): comparison of methods.

Target: the squared two-sided t critical value c^2 = t_{1-alpha/2, df}^2,
reference SciPy `t.ppf`, error = relative error (MaxRE, MRE).

Plaintext part (--part plain). Every method is evaluated under two df-range
conditions, each with its own evaluation grid of EVAL_POINTS df (log-spaced
over the range, plus the Welch df of the cases for a case range):
- global : df in [1, 10^4];
- case   : df in [min(n1, n2) - 1, n1 + n2 - 2] for each dataset case (the
           Welch-Satterthwaite df always lies in this public range).

Polynomial methods (Chebyshev interpolation, degrees 2^3-1 .. 2^7-1, fitted on
the same range they are evaluated on):
- invdf : input u = 1/df (the global fit is the stored coefficient table,
          coeffs/t_critical_coeffs.py, fitted on u in [0, 1], i.e. df >= 1);
- df    : input df.

For the global condition, MaxRE is also reported per df band
(1-2, 2-10, 10-30, 30-10^4).

HE part (--part he). The same df points (two ciphertexts of 32768 slots) and
df-range conditions (global and one per distinct case range), encrypted
u = 1/df at level 12, `--reps` repetitions; runtime, bootstraps, levels,
relative error against the exact quantile and absolute error against the
same method in plaintext:
- invdf : HEHypothesisTesting._critical_value_from_inv_df (stored coefficients,
          for every condition) (results/he_raw.csv);
- invdf_range : 1/df polynomial refitted on the case range (case conditions);
- df    : df = 1/u in HE (invSqrt of u * 0.01 df_hi, which lies in the
          invSqrt domain, squared), then the df polynomial fitted on the
          condition range (results/he_raw.csv, stage time of the division in
          nu_time).

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/critical_boundary/run_critical_boundary.py --part plain
    python3 experiments/critical_boundary/run_critical_boundary.py --part he
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

from coeffs.invSqrt_coeffs import _INV_SQRT_MIN
from coeffs.t_critical_coeffs import (
    _T_CRITICAL_COEFFS,
    _T_CRITICAL_U_MIN,
    _T_CRITICAL_U_MAX,
)


ALPHAS = [0.01, 0.025, 0.05, 0.1]
LOG_DEGREES = range(3, 8)  # degree = 2^3-1, ..., 2^7-1
DF_MIN, DF_MAX = 1.0, 10000.0
DF_BANDS = [(1.0, 2.0), (2.0, 10.0), (10.0, 30.0), (30.0, 10000.0)]
REPS = 2
EVAL_POINTS = 65536  # two ciphertexts of 32768 slots

RESULT_DIR = Path(__file__).resolve().parent / "results"


def exact_sq(alpha, df):

    return stats.t.ppf(1.0 - alpha / 2.0, df) ** 2




def errors(y, y_true):

    re = np.abs(y - y_true) / np.abs(y_true)
    return float(re.max()), float(re.mean())


def run_plain(out_dir):

    conditions = he_conditions()
    rows = []

    for range_name, condition in conditions.items():
        lo, hi = condition["range"]
        df = eval_grid(condition)[0]
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

            for ld in LOG_DEGREES:
                degree = 2 ** ld - 1

                # 1/df input.
                if range_name == "global":
                    z = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
                    y_inv = np.polynomial.chebyshev.chebval(z, _T_CRITICAL_COEFFS[alpha][degree])
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


def he_conditions():

    # Global range and every distinct case range, with the Welch df of the
    # cases sharing that range (Adult: both comparisons).
    from common.cases import load_cases

    out = {"global": {"range": (DF_MIN, DF_MAX), "welch": {}}}

    for name, case in load_cases().items():
        lo, hi = float(min(case["n1"], case["n2"]) - 1), float(case["n1"] + case["n2"] - 2)
        a1 = np.var(case["x1"], ddof=1) / case["n1"]
        a2 = np.var(case["x2"], ddof=1) / case["n2"]
        df = (a1 + a2) ** 2 / (a1 ** 2 / (case["n1"] - 1) + a2 ** 2 / (case["n2"] - 1))
        key = next((k for k, c in out.items() if c["range"] == (lo, hi)), name)
        out.setdefault(key, {"range": (lo, hi), "welch": {}})["welch"][name] = float(df)

    return out


def eval_grid(condition):

    # EVAL_POINTS df: log-spaced over the range plus the Welch df of the cases.
    lo, hi = condition["range"]
    welch = np.array(sorted(condition["welch"].values()))
    return np.concatenate([np.geomspace(lo, hi, EVAL_POINTS - len(welch)), welch]), welch


def run_he_poly(engine, ht, conditions, reps, path, alphas):

    from engine.HEdata import Ciphertext, Message

    fields = ["condition", "df_lo", "df_hi", "method", "degree", "alpha", "rep",
              "time", "nu_time", "input_level", "output_level", "bootstrap_count", "nu_bootstrap_count",
              "max_re", "mre", "max_ae_vs_plain", "max_re_welch", "nu_max_re"]
    rows = []

    for name, condition in conditions.items():
        lo, hi = condition["range"]
        df, welch = eval_grid(condition)
        u = 1.0 / df
        welch_mask = np.isin(df, welch)
        z_inv = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
        dec = lambda c: np.concatenate([np.array(m, dtype=np.complex128).real for m in engine.dec(c)])[:len(df)]

        def add(method, degree, alpha, rep, y, y_plain, **kw):
            y_true = exact_sq(alpha, df)
            max_re, mre = errors(y, y_true)
            rows.append({"condition": name, "df_lo": lo, "df_hi": hi, "method": method,
                         "degree": degree, "alpha": alpha, "rep": rep, "max_re": max_re, "mre": mre,
                         "max_ae_vs_plain": float(np.max(np.abs(y - y_plain))),
                         "max_re_welch": errors(y[welch_mask], y_true[welch_mask])[0] if welch_mask.any() else "",
                         **kw})

        for rep in range(reps):
            # 1/df input, stored coefficients.
            for ld in LOG_DEGREES:
                degree = 2 ** ld - 1

                for alpha in alphas:
                    ctxt = engine.enc(Message(u, engine.log_slots), level=12)
                    ht.reset_bootstrap_count()
                    start = time.perf_counter()
                    out = ht._critical_value_from_inv_df(ctxt, alpha, ld)
                    elapsed = time.perf_counter() - start
                    add("invdf", degree, alpha, rep, dec(out),
                        np.polynomial.chebyshev.chebval(z_inv, _T_CRITICAL_COEFFS[alpha][degree]),
                        time=elapsed, input_level=12, output_level=out.level(),
                        bootstrap_count=ht.bootstrap_count())

            # 1/df input, refitted on the case range [1/hi, 1/lo].
            if name != "global":
                for ld in LOG_DEGREES:
                    degree = 2 ** ld - 1

                    for alpha in alphas:
                        p = Chebyshev.interpolate(lambda v: exact_sq(alpha, 1.0 / v), degree, domain=[1.0 / hi, 1.0 / lo])
                        ctxt = engine.enc(Message(u, engine.log_slots), level=12)
                        ht.reset_bootstrap_count()
                        start = time.perf_counter()
                        zr = engine.mult(ctxt, 2.0 / (1.0 / lo - 1.0 / hi))
                        zr = engine.sub(zr, (1.0 / lo + 1.0 / hi) / (1.0 / lo - 1.0 / hi))
                        ht.approx._ensure_cheb_level(zr, degree)
                        out = engine.evaluate_chebyshev(zr, engine.make_cheb_coeffs(p.coef))
                        elapsed = time.perf_counter() - start
                        add("invdf_range", degree, alpha, rep, dec(out), p(u),
                            time=elapsed, input_level=12, output_level=out.level(),
                            bootstrap_count=ht.bootstrap_count())

            # df input: df = 1/u, then the df polynomial on [lo, hi].
            ctxt = engine.enc(Message(u, engine.log_slots), level=12)
            ht.reset_bootstrap_count()
            start = time.perf_counter()
            z = he_df_input(engine, ht, ctxt, lo, hi)
            nu_time = time.perf_counter() - start
            nu_bts = ht.bootstrap_count()
            nu = (dec(z) * (hi - lo) + (hi + lo)) / 2.0
            nu_max_re = errors(nu, df)[0]

            for ld in LOG_DEGREES:
                degree = 2 ** ld - 1

                for alpha in alphas:
                    p = Chebyshev.interpolate(lambda d: exact_sq(alpha, d), degree, domain=[lo, hi])
                    zc = Ciphertext(z)
                    ht.reset_bootstrap_count()
                    start = time.perf_counter()
                    ht.approx._ensure_cheb_level(zc, degree)
                    out = engine.evaluate_chebyshev(zc, engine.make_cheb_coeffs(p.coef))
                    elapsed = time.perf_counter() - start
                    add("df", degree, alpha, rep, dec(out), p(df),
                        time=nu_time + elapsed, nu_time=nu_time, input_level=12, output_level=out.level(),
                        bootstrap_count=nu_bts + ht.bootstrap_count(), nu_bootstrap_count=nu_bts,
                        nu_max_re=nu_max_re)

            g = [r for r in rows if r["condition"] == name and r["rep"] == rep]
            print(f"{name} rep={rep} nu: {nu_time:.2f}s MaxRE={nu_max_re:.1e} | "
                  + " ".join(f"{m}15 a={alphas[-1]} {next(r['max_re'] for r in g if r['method'] == m and r['degree'] == 15 and r['alpha'] == alphas[-1]):.1e}"
                             for m in ["invdf", "df"]), flush=True)

    if set(alphas) != set(ALPHAS) and path.is_file():
        # Keep the stored rows of the other alphas.
        with path.open() as file:
            rows = [r for r in csv.DictReader(file) if float(r["alpha"]) not in alphas] + rows

    write_csv(path, fields, rows)


def he_df_input(engine, ht, ctxt_u, lo, hi):

    # df = 1/u for df in [lo, hi]: w = c u with c = _INV_SQRT_MIN * hi lies in
    # [_INV_SQRT_MIN, _INV_SQRT_MIN * hi / lo], inside the invSqrt domain,
    # y = 1/sqrt(w) (|y| <= 1 / sqrt(_INV_SQRT_MIN)) and df = c y^2. Returns
    # df mapped to [-1, 1]: z = 2 df / (hi - lo) - (hi + lo) / (hi - lo)
    # = (k y)^2 - (hi + lo) / (hi - lo) with k = sqrt(2 c / (hi - lo)).
    c = _INV_SQRT_MIN * hi
    w = engine.mult(ctxt_u, c)

    y = ht._inv_sqrt(w, c / lo)
    ht._ensure_level(y, 5, 1.0 / math.sqrt(_INV_SQRT_MIN))

    yk = engine.mult(y, math.sqrt(2.0 * c / (hi - lo)))

    return engine.sub(engine.mult(yk, yk), (hi + lo) / (hi - lo))


def write_csv(path, fields, rows):

    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved {path}")


def run_he(out_dir, reps, alphas):

    from common.he_utils import make_engine
    from operators.HEhypo_test import HEHypothesisTesting

    engine = make_engine()
    ht = HEHypothesisTesting(engine)
    conditions = he_conditions()

    run_he_poly(engine, ht, conditions, reps, out_dir / "he_raw.csv", alphas)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--part", choices=["plain", "he"], required=True)
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--alphas", default="all")
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.part == "plain":
        run_plain(args.out_dir)
    else:
        run_he(args.out_dir, args.reps,
               ALPHAS if args.alphas == "all" else [float(a) for a in args.alphas.split(",")])


if __name__ == "__main__":
    main()
