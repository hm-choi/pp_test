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
- integer_interp : every integer df up to the largest df of all ranges,
                   linear interpolation of c^2 in 1/df.

For the global condition, MaxRE is also reported per df band
(1-2, 2-10, 10-30, 30-2000).

HE part (--part he). The same df points and df-range conditions (global and
one per distinct case range), encrypted u = 1/df at level 12, `--reps`
repetitions; runtime, bootstraps, levels, relative error against the exact
quantile and absolute error against the same method in plaintext:
- invdf : HEHypothesisTesting._critical_value_from_inv_df (stored coefficients,
          for every condition), every slot one df (results/he_raw.csv);
- invdf_range : 1/df polynomial refitted on the case range (case conditions);
- df    : df = 1/u in HE (invSqrt of u * df_lo, squared), then the df
          polynomial fitted on the condition range (results/he_raw.csv,
          stage time of the division in nu_time);
- lookup tables without division: the comparisons df >= d_j are evaluated as
          u <= 1/d_j, i.e. one encrypted sign of (u - knot_j) / width per table
          entry, packed in the slots (one ciphertext per 32768 entries); floor
          and linear interpolation in 1/df are weighted slot sums of the same
          step (floor) or of (u - knot_j) * step (interpolation), with every
          alpha's weights public. One df per ciphertext, at LOOKUP_POINTS
          interior points of each range plus the Welch df of the cases
          (results/he_lookup_raw.csv).

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/tcrit/tcrit_experiment.py --part plain
    python3 experiments/tcrit/tcrit_experiment.py --part he [--methods poly,lookup]
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

# HE part.
LOOKUP_POINTS = {"global": 20, "case": 3}
# (log_degree, Newton iterations) of invSqrt in the df-input baseline
# (domain [df_lo / df_hi, 1]).
DF_INPUT_INVSQRT = {"global": (6, 4), "case": (5, 1)}

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
    all_ranges = ranges()
    integer_df = np.arange(1.0, max(hi for _, hi in all_ranges.values()) + 1.0)
    rows = []

    for range_name, (lo, hi) in all_ranges.items():
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


def he_conditions():

    # Global range and every distinct case range, with the Welch df of the
    # cases sharing that range (Adult: both comparisons).
    from pptest_cases import load_cases

    out = {"global": {"range": (DF_MIN, DF_MAX), "welch": {}}}

    for name, case in load_cases().items():
        lo, hi = float(min(case["n1"], case["n2"]) - 1), float(case["n1"] + case["n2"] - 2)
        a1 = np.var(case["x1"], ddof=1) / case["n1"]
        a2 = np.var(case["x2"], ddof=1) / case["n2"]
        df = (a1 + a2) ** 2 / (a1 ** 2 / (case["n1"] - 1) + a2 ** 2 / (case["n2"] - 1))
        key = next((k for k, c in out.items() if c["range"] == (lo, hi)), name)
        out.setdefault(key, {"range": (lo, hi), "welch": {}})["welch"][name] = float(df)

    return out


def poly_grid(condition, num_slots):

    lo, hi = condition["range"]
    welch = np.array(sorted(condition["welch"].values()))
    return np.concatenate([np.geomspace(lo, hi, num_slots - len(welch)), welch]), welch


def lookup_points(name, condition):

    lo, hi = condition["range"]
    k = LOOKUP_POINTS["global" if name == "global" else "case"]
    points = [(f"grid{i}", float(d)) for i, d in enumerate(np.geomspace(lo, hi, k + 2)[1:-1])]
    return points + [(f"welch:{case}", df) for case, df in condition["welch"].items()]


def lookup_table(alpha, knot_df, with_infinity):

    # Knots in ascending u = 1/df (u = 0 for df = infinity) and the weights of
    # floor(u) = base + sum_j dF_j step(u - u_j) and
    # interp(u) = base_i + sum_j dS_j max(u - u_j, 0).
    df_desc = np.sort(knot_df)[::-1]
    u = 1.0 / df_desc
    f = exact_sq(alpha, df_desc)

    if with_infinity:
        u = np.concatenate([[0.0], u])
        f_inf = stats.norm.ppf(1.0 - alpha / 2.0) ** 2
        f_all = np.concatenate([[f_inf], f])
    else:
        f_all = f

    # Floor: the infinity knot (if any) has weight 0 (largest df = df_desc[0]).
    floor_w = np.zeros(len(u))
    off = 1 if with_infinity else 0
    floor_w[off:-1] = np.diff(f)
    slopes = np.diff(f_all) / np.diff(u)
    interp_w = np.zeros(len(u))
    interp_w[:-1] = np.diff(np.concatenate([[0.0], slopes]))

    return u, {"floor": (float(f[0]), floor_w), "interp": (float(f_all[0]), interp_w)}


def run_he_poly(engine, ht, conditions, reps, path):

    from engine.HEdata import Ciphertext, Message

    fields = ["condition", "df_lo", "df_hi", "method", "degree", "alpha", "rep",
              "time", "nu_time", "input_level", "output_level", "bootstrap_count", "nu_bootstrap_count",
              "max_re", "mre", "max_ae_vs_plain", "max_re_welch", "nu_max_re"]
    rows = []

    for name, condition in conditions.items():
        lo, hi = condition["range"]
        df, welch = poly_grid(condition, engine.num_slots)
        u = 1.0 / df
        welch_mask = np.isin(df, welch)
        z_inv = 2.0 * (u - _T_CRITICAL_U_MIN) / (_T_CRITICAL_U_MAX - _T_CRITICAL_U_MIN) - 1.0
        dec = lambda c: np.array(engine.dec(c)[0], dtype=np.complex128).real[:len(df)]

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

                for alpha in ALPHAS:
                    ctxt = engine.enc(Message(u, engine.log_slots), level=12)
                    ht.reset_bootstrap_count()
                    start = time.perf_counter()
                    out = ht._critical_value_from_inv_df(ctxt, alpha, ld)
                    elapsed = time.perf_counter() - start
                    add("invdf", degree, alpha, rep, dec(out),
                        np.polynomial.chebyshev.chebval(z_inv, _T_CRITICAL_SQ_COEFFS[degree][alpha]),
                        time=elapsed, input_level=12, output_level=out.level(),
                        bootstrap_count=ht.bootstrap_count())

            # 1/df input, refitted on the case range [1/hi, 1/lo].
            if name != "global":
                for ld in LOG_DEGREES:
                    degree = 2 ** ld - 1

                    for alpha in ALPHAS:
                        p = Chebyshev.interpolate(lambda v: exact_sq(alpha, 1.0 / v), degree, domain=[1.0 / hi, 1.0 / lo])
                        ctxt = engine.enc(Message(u, engine.log_slots), level=12)
                        ht.reset_bootstrap_count()
                        start = time.perf_counter()
                        zr = engine.mult(ctxt, 2.0 / (1.0 / lo - 1.0 / hi))
                        zr = engine.sub(zr, (1.0 / lo + 1.0 / hi) / (1.0 / lo - 1.0 / hi))
                        ht.approx._ensure_cheb_level(zr, degree)
                        out = engine.evaluate_chebyshev(zr, engine._make_cheb_coeffs(p.coef))
                        elapsed = time.perf_counter() - start
                        add("invdf_range", degree, alpha, rep, dec(out), p(u),
                            time=elapsed, input_level=12, output_level=out.level(),
                            bootstrap_count=ht.bootstrap_count())

            # df input: df = 1/u, then the df polynomial on [lo, hi].
            ctxt = engine.enc(Message(u, engine.log_slots), level=12)
            ht.reset_bootstrap_count()
            start = time.perf_counter()
            z = he_df_input(engine, ht.approx, ctxt, lo, hi,
                            *DF_INPUT_INVSQRT["global" if name == "global" else "case"])
            nu_time = time.perf_counter() - start
            nu_bts = ht.bootstrap_count()
            nu = (dec(z) * (hi - lo) + (hi + lo)) / 2.0
            nu_max_re = errors(nu, df)[0]

            for ld in LOG_DEGREES:
                degree = 2 ** ld - 1

                for alpha in ALPHAS:
                    p = Chebyshev.interpolate(lambda d: exact_sq(alpha, d), degree, domain=[lo, hi])
                    zc = Ciphertext(z)
                    ht.reset_bootstrap_count()
                    start = time.perf_counter()
                    ht.approx._ensure_cheb_level(zc, degree)
                    out = engine.evaluate_chebyshev(zc, engine._make_cheb_coeffs(p.coef))
                    elapsed = time.perf_counter() - start
                    add("df", degree, alpha, rep, dec(out), p(df),
                        time=nu_time + elapsed, nu_time=nu_time, input_level=12, output_level=out.level(),
                        bootstrap_count=nu_bts + ht.bootstrap_count(), nu_bootstrap_count=nu_bts,
                        nu_max_re=nu_max_re)

            g = [r for r in rows if r["condition"] == name and r["rep"] == rep]
            print(f"{name} rep={rep} nu: {nu_time:.2f}s MaxRE={nu_max_re:.1e} | "
                  + " ".join(f"{m}15 a=0.05 {next(r['max_re'] for r in g if r['method'] == m and r['degree'] == 15 and r['alpha'] == 0.05):.1e}"
                             for m in ["invdf", "df"]), flush=True)

    write_csv(path, fields, rows)


def he_df_input(engine, approx, ctxt_u, lo, hi, log_degree, iteration):

    # df = 1/u for df in [lo, hi]: w = lo * u lies in [lo / hi, 1], y = 1/sqrt(w)
    # (invSqrt, |y| <= sqrt(hi / lo)) and df = lo * y^2. Returns df mapped to
    # [-1, 1]: z = 2 df / (hi - lo) - (hi + lo) / (hi - lo) = (k y)^2 - c.
    w_lo = lo / hi
    x_half = engine.mult(ctxt_u, lo / 2.0)
    x_cheb = engine.mult(ctxt_u, 2.0 * lo / (1.0 - w_lo))
    x_cheb = engine.sub(x_cheb, (1.0 + w_lo) / (1.0 - w_lo))

    y = approx.invSqrt(x_half, x_cheb, (w_lo, 1.0), log_degree, iteration)
    approx.ensure_level(y, 5, math.sqrt(hi / lo))

    yk = engine.mult(y, math.sqrt(2.0 * lo / (hi - lo)))

    return engine.sub(engine.mult(yk, yk), (hi + lo) / (hi - lo))


def run_he_lookup(engine, ht, conditions, reps, path):

    from engine.HEdata import Message

    fields = ["condition", "df_lo", "df_hi", "point", "df", "table", "mapping", "alpha", "rep",
              "entries", "ciphertexts", "time", "compare_time", "combine_time", "bootstrap_count",
              "he_value", "plain_value", "exact", "re", "ae_vs_plain"]
    done = set()

    if path.is_file():
        with path.open() as file:
            done = {(r["condition"], r["point"], r["table"], r["rep"]) for r in csv.DictReader(file)}

    new_file = not path.is_file()
    file = path.open("a", newline="")
    writer = csv.DictWriter(file, fieldnames=fields)

    if new_file:
        writer.writeheader()

    textbook_df = textbook_df_grid()
    S = engine.num_slots

    for name, condition in conditions.items():
        lo, hi = condition["range"]
        tables = {
            "textbook": (textbook_df, True, 1.0),
            "integer": (np.arange(lo, hi + 1.0), False, 1.0 / lo - 1.0 / hi),
        }

        for point, df in lookup_points(name, condition):
            u = 1.0 / df

            for table, (knot_df, with_infinity, width) in tables.items():
                knots, _ = lookup_table(ALPHAS[0], knot_df, with_infinity)
                m = len(knots)
                parts = -(-m // S)
                pad = lambda a: np.concatenate([a, np.zeros(parts * S - m)])
                plain = {
                    ("floor", a): float(exact_sq(a, lookup_floor(textbook_df, df)) if table == "textbook"
                                        else exact_sq(a, math.floor(df)))
                    for a in ALPHAS
                }
                plain.update({("interp", a): float(lookup_interp(a, knot_df, np.array([df]), with_infinity)[0])
                              for a in ALPHAS})
                weights = {a: lookup_table(a, knot_df, with_infinity)[1] for a in ALPHAS}
                knot_msg = Message(pad(knots), engine.log_slots)
                scale_msg = Message(pad(np.full(m, 1.0 / width)), engine.log_slots)

                for rep in range(reps):
                    if (name, point, table, str(rep)) in done:
                        continue

                    ctxt = engine.enc(Message(np.full(parts * S, u), engine.log_slots), level=12)
                    ht.reset_bootstrap_count()
                    start = time.perf_counter()
                    w = engine.sub(ctxt, knot_msg)
                    sign = ht.approx.sign(engine.mult(w, scale_msg))
                    step = engine.mult(engine.add(sign, 1.0), 0.5)
                    compare_time = time.perf_counter() - start
                    bts = ht.bootstrap_count()

                    start = time.perf_counter()
                    relu = engine.mult(w, step)
                    relu_time = time.perf_counter() - start

                    for alpha in ALPHAS:
                        for mapping, src in (("floor", step), ("interp", relu)):
                            base, wv = weights[alpha][mapping]
                            start = time.perf_counter()
                            out = engine.add(engine.sum(engine.mult(src, Message(pad(wv), engine.log_slots))), base)
                            combine_time = time.perf_counter() - start + (relu_time if mapping == "interp" else 0.0)
                            he = float(np.array(engine.dec(out)[0], dtype=np.complex128).real[0])
                            y_true = float(exact_sq(alpha, df))
                            pv = plain[(mapping, alpha)]
                            writer.writerow({
                                "condition": name, "df_lo": lo, "df_hi": hi, "point": point, "df": df,
                                "table": table, "mapping": mapping, "alpha": alpha, "rep": rep,
                                "entries": m, "ciphertexts": parts,
                                "time": compare_time + combine_time, "compare_time": compare_time,
                                "combine_time": combine_time, "bootstrap_count": bts,
                                "he_value": he, "plain_value": pv, "exact": y_true,
                                "re": abs(he - y_true) / y_true, "ae_vs_plain": abs(he - pv),
                            })

                    file.flush()
                    print(f"{name} {point} df={df:.2f} {table} rep={rep} compare={compare_time:.1f}s", flush=True)

    file.close()


def write_csv(path, fields, rows):

    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved {path}")


def run_he(out_dir, reps, methods):

    from engine.HEengine import HEengine
    from operators.HEhypo_test import HEHypothesisTesting

    engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")
    ht = HEHypothesisTesting(engine)
    conditions = he_conditions()

    if "poly" in methods:
        run_he_poly(engine, ht, conditions, reps, out_dir / "he_raw.csv")

    if "lookup" in methods:
        run_he_lookup(engine, ht, conditions, reps, out_dir / "he_lookup_raw.csv")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--part", choices=["plain", "he"], required=True)
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--methods", default="poly,lookup")
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.part == "plain":
        run_plain(args.out_dir)
    else:
        run_he(args.out_dir, args.reps, args.methods.split(","))


if __name__ == "__main__":
    main()
