"""Sections 4.3 / 4.6: accuracy and performance of encrypted hypothesis tests.

For every dataset case (experiments/pptest_cases.py), test (Welch, F, Z),
invSqrt configuration (Welch and F) and repetition, the encrypted statistics
are computed once and then every significance level (and, for Welch, both
critical-value targets t and t^2) is evaluated as a separate decision branch.
Copies of all intermediate ciphertexts are decrypted afterwards and compared
with the plaintext (SciPy) values:

- Welch : T^2, 1/df, critical value, score, sign, step; p-value from the
          decrypted T^2 and df vs. the SciPy p-value;
- F     : F, score, sign, step; p-value from the decrypted F;
- Z     : Z, score, sign, step; p-value from the decrypted Z.

invSqrt configurations:
- normalized_vmin{1e-3,1e-4,1e-5} : x = V / V_max in [v_min, 1], HE-DAP
    accuracy-oriented choice (u1) per input level
    (experiments/invsqrt/results/normalized/hedap_optimal.json);
- raw_default : x = V in [1e-3, V_max], degree 63, 7 iterations;
- raw_hedap   : x = V in [1e-3, V_max], HE-DAP u1 per input level for the
    case/test domain (experiments/invsqrt/results/raw/hedap_optimal.json).

Score bounds are public constants per case and test chosen as in
experiments/experiment1 (|score / B| <= 0.5 for every alpha), inside the valid
input range of the sign approximation (1.5e-7 < |x| < 1).

Timing: statistics and each decision branch are timed separately (encryption
and decryption excluded); the total for one test at one alpha is
statistics + that branch.

Product-score comparison (--score product, results in results/product/): the
decision scores are multiplied out so that the scores themselves contain no
division by an encrypted variance and have public bounds (Welch still needs
invSqrt for 1/df, so V must stay inside the invSqrt domain):
- Welch : s'  = (mean1 - mean2)^2 - c^2 V,  B = max(R^2, c_max^2 V_max);
- F     : s'' = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2),
          B = max(A, F_L B_v) max(A, F_U B_v)  (no invSqrt in the decision);
- Z     : unchanged score, public bound max(R^2 / V_public, z^2).

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/pptest_eval/run_eval.py [--score product]
"""

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "experiments"))

from pptest_cases import load_cases


ALPHAS = [0.01, 0.025, 0.05, 0.1]
TARGETS = ["t", "t2"]
CRITICAL_LOG_DEGREE = 4  # degree 15
REPS = 10
TESTS = ["welch", "f", "z"]

INVSQRT_DIR = PROJECT_ROOT / "experiments" / "invsqrt" / "results"
RESULT_DIR = Path(__file__).resolve().parent / "results"

RAW_FIELDS = [
    "case", "test", "config", "rep", "alpha", "target",
    "stats_time", "decision_time", "total_time", "stats_bootstraps", "decision_bootstraps",
    "plain_stat", "he_stat", "stat_abs_err", "stat_rel_err",
    "plain_inv_df", "he_inv_df", "inv_df_abs_err",
    "plain_critical", "he_critical", "critical_abs_err",
    "plain_score", "he_score", "score_abs_err", "score_bound",
    "he_sign", "he_step", "plain_decision", "step_decision", "score_decision",
    "plain_p", "he_p", "p_abs_err", "p_decision",
]


# ----------------------------------------------------------------------
# Plaintext reference
# ----------------------------------------------------------------------

def welch_plain(x1, x2, alpha):

    n1, n2 = len(x1), len(x2)
    v1, v2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    a1, a2 = v1 / n1, v2 / n2
    v = a1 + a2
    t2 = (np.mean(x1) - np.mean(x2)) ** 2 / v
    df = v ** 2 / (a1 ** 2 / (n1 - 1) + a2 ** 2 / (n2 - 1))
    crit = stats.t.ppf(1.0 - alpha / 2.0, df)
    p = 2.0 * stats.t.sf(np.sqrt(t2), df)

    return {"stat": t2, "inv_df": 1.0 / df, "critical": crit, "score": t2 - crit ** 2, "p": p,
            "score_product": (np.mean(x1) - np.mean(x2)) ** 2 - crit ** 2 * v}


def f_plain(x1, x2, alpha):

    n1, n2 = len(x1), len(x2)
    f = np.var(x1, ddof=1) / np.var(x2, ddof=1)
    fl = stats.f.ppf(alpha / 2.0, n1 - 1, n2 - 1)
    fu = stats.f.ppf(1.0 - alpha / 2.0, n1 - 1, n2 - 1)
    p = 2.0 * min(stats.f.cdf(f, n1 - 1, n2 - 1), stats.f.sf(f, n1 - 1, n2 - 1))

    v1, v2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    return {"stat": f, "fl": fl, "fu": fu, "score": (f - fl) * (f - fu), "p": p,
            "score_product": (v1 - fl * v2) * (v1 - fu * v2)}


def z_plain(x1, x2, alpha, s1, s2):

    n1, n2 = len(x1), len(x2)
    z = (np.mean(x1) - np.mean(x2)) / np.sqrt(s1 / n1 + s2 / n2)
    crit = stats.norm.ppf(1.0 - alpha / 2.0)
    p = 2.0 * stats.norm.sf(abs(z))

    return {"stat": z, "critical": crit, "score": z ** 2 - crit ** 2, "p": p,
            "score_product": z ** 2 - crit ** 2}


def score_bound(scores):

    # Public normalization bound: 2 * max |score|, rounded up to one
    # significant digit (|score / B| <= 0.5).
    m = 2.0 * max(abs(s) for s in scores)
    e = math.floor(math.log10(m))
    return math.ceil(m / 10 ** e) * 10 ** e


# ----------------------------------------------------------------------
# invSqrt configurations
# ----------------------------------------------------------------------

def select_table(path, domain):

    levels = json.loads(Path(path).read_text())[domain]["levels"]

    return {
        int(level): (entry["accuracy"]["degree"].bit_length(), entry["accuracy"]["iteration"], bool(entry["accuracy"]["pre_bts"]))
        for level, entry in levels.items()
    }


INVSQRT_CONFIG_NAMES = [
    "normalized_vmin0.001", "normalized_vmin0.0001", "normalized_vmin1e-05",
    "raw_default", "raw_hedap",
]


def invsqrt_config(name, case_name, test):

    if name.startswith("normalized_vmin"):
        return {
            "method": "normalized", "v_min": float(name.removeprefix("normalized_vmin")),
            "select": select_table(INVSQRT_DIR / "normalized" / "hedap_optimal.json", name),
        }

    if name == "raw_default":
        return {"method": "raw", "v_min": 1e-3, "log_degree": 6, "iteration": 7, "pre_bts": False}

    if name == "raw_hedap":
        return {
            "method": "raw", "v_min": 1e-3,
            "select": select_table(INVSQRT_DIR / "raw" / "hedap_optimal.json", f"raw_{case_name}_{test}"),
        }

    raise ValueError(f"Unknown invSqrt configuration: {name}")


# ----------------------------------------------------------------------
# HE runs
# ----------------------------------------------------------------------

def first(engine, ctxt):

    return float(np.array(engine.dec(ctxt)[0], dtype=np.complex128).real[0])


def timed(ht, func):

    ht.reset_bootstrap_count()
    start = time.perf_counter()
    out = func()
    return out, time.perf_counter() - start, ht.bootstrap_count()


def run_welch(engine, ht, case, c1, c2, bounds):

    tr = {}

    (t2, inv_df), st, sb = timed(ht, lambda: ht.HE_Welch_statistics(
        c1, c2, case["n1"], case["n2"], case["R"], tr))

    branches = []

    for alpha in ALPHAS:
        for target in TARGETS:
            btr = {}
            step, dt, db = timed(ht, lambda: ht.HE_Welch_decision(
                t2, inv_df, alpha, CRITICAL_LOG_DEGREE, bounds[alpha], target, btr))
            branches.append((alpha, target, step, dt, db, btr))

    he_t2 = first(engine, tr["t_squared"])
    he_inv_df = first(engine, tr["inv_df"])
    rows = []

    for alpha, target, step, dt, db, btr in branches:
        B = bounds[alpha]
        cs = first(engine, btr["critical_scaled"])
        he_crit = cs * math.sqrt(B) if target == "t" else math.sqrt(max(cs * B, 0.0))
        he_df = 1.0 / he_inv_df
        he_p = 2.0 * stats.t.sf(math.sqrt(max(he_t2, 0.0)), he_df)
        rows.append(dict(alpha=alpha, target=target, stats_time=st, decision_time=dt,
                         stats_bootstraps=sb, decision_bootstraps=db,
                         he_stat=he_t2, he_inv_df=he_inv_df, he_critical=he_crit,
                         he_score=first(engine, btr["score_normalized"]) * B, score_bound=B,
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step), he_p=he_p))

    return rows


def run_f(engine, ht, case, c1, c2, bounds, plain):

    tr = {}

    (var1_s, inv_x), st, sb = timed(ht, lambda: ht.HE_F_statistics(
        c1, c2, case["n1"], case["n2"], case["R"], tr))

    branches = []

    for alpha in ALPHAS:
        btr = {}
        p = plain[alpha]
        step, dt, db = timed(ht, lambda: ht.HE_F_decision(
            var1_s, inv_x, p["fl"], p["fu"], bounds[alpha], btr))
        branches.append((alpha, step, dt, db, btr))

    rows = []
    n1, n2 = case["n1"], case["n2"]

    for alpha, step, dt, db, btr in branches:
        he_f = first(engine, btr["f_statistic"])
        he_p = 2.0 * min(stats.f.cdf(he_f, n1 - 1, n2 - 1), stats.f.sf(he_f, n1 - 1, n2 - 1))
        rows.append(dict(alpha=alpha, target="", stats_time=st, decision_time=dt,
                         stats_bootstraps=sb, decision_bootstraps=db,
                         he_stat=he_f, he_score=first(engine, btr["score"]), score_bound=bounds[alpha],
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step), he_p=he_p))

    return rows


def run_z(engine, ht, case, c1, c2, bounds):

    s1, s2 = np.var(case["x1"], ddof=1), np.var(case["x2"], ddof=1)
    tr = {}

    z, st, sb = timed(ht, lambda: ht.HE_Z_statistic(c1, c2, case["n1"], case["n2"], s1, s2, 0.0, tr))

    branches = []

    for alpha in ALPHAS:
        btr = {}
        step, dt, db = timed(ht, lambda: ht.HE_Z_decision(z, alpha, bounds[alpha], btr))
        branches.append((alpha, step, dt, db, btr))

    he_z = first(engine, tr["z"])
    rows = []

    for alpha, step, dt, db, btr in branches:
        rows.append(dict(alpha=alpha, target="", stats_time=st, decision_time=dt,
                         stats_bootstraps=sb, decision_bootstraps=db,
                         he_stat=he_z, he_critical=stats.norm.ppf(1.0 - alpha / 2.0),
                         he_score=first(engine, btr["score"]), score_bound=bounds[alpha],
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step),
                         he_p=2.0 * stats.norm.sf(abs(he_z))))

    return rows


def run_welch_product(engine, ht, case, c1, c2):

    # score="product": s' = (mean1 - mean2)^2 - c^2 V with the public bound.
    tr = {}

    (t2, inv_df, terms), st, sb = timed(ht, lambda: ht.HE_Welch_statistics(
        c1, c2, case["n1"], case["n2"], case["R"], tr, return_terms=True))

    branches = []

    for alpha in ALPHAS:
        for target in TARGETS:
            btr = {}
            step, dt, db = timed(ht, lambda: ht.HE_Welch_decision(
                t2, inv_df, alpha, CRITICAL_LOG_DEGREE, None, target, btr, "product", terms))
            branches.append((alpha, target, step, dt, db, btr))

    he_t2 = first(engine, tr["t_squared"])
    he_inv_df = first(engine, tr["inv_df"])
    rows = []

    for alpha, target, step, dt, db, btr in branches:
        B = btr["score_bound"]
        norm = B * terms["x_scale"]
        cs = first(engine, btr["critical_scaled"])
        he_crit = cs * math.sqrt(norm) if target == "t" else math.sqrt(max(cs * norm, 0.0))
        he_p = 2.0 * stats.t.sf(math.sqrt(max(he_t2, 0.0)), 1.0 / he_inv_df)
        rows.append(dict(alpha=alpha, target=target, stats_time=st, decision_time=dt,
                         stats_bootstraps=sb, decision_bootstraps=db,
                         he_stat=he_t2, he_inv_df=he_inv_df, he_critical=he_crit,
                         he_score=first(engine, btr["score_normalized"]) * B, score_bound=B,
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step), he_p=he_p))

    return rows


def run_f_product(engine, ht, case, c1, c2, plain):

    # score="product": s'' = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2); no invSqrt,
    # so the F statistic itself (and its p-value) is not computed here.
    terms, st, sb = timed(ht, lambda: ht.HE_F_terms(c1, c2, case["n1"], case["n2"], case["R"]))

    branches = []

    for alpha in ALPHAS:
        btr = {}
        p = plain[alpha]
        step, dt, db = timed(ht, lambda: ht.HE_F_product_decision(terms, p["fl"], p["fu"], btr))
        branches.append((alpha, step, dt, db, btr))

    rows = []

    for alpha, step, dt, db, btr in branches:
        B = btr["score_bound"]
        rows.append(dict(alpha=alpha, target="", stats_time=st, decision_time=dt,
                         stats_bootstraps=sb, decision_bootstraps=db,
                         he_stat=None, he_score=first(engine, btr["score_normalized"]) * B, score_bound=B,
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step), he_p=None))

    return rows


def complete(row, plain, alpha):

    p = plain[alpha]
    row["plain_stat"] = p["stat"]
    row["plain_score"] = p["score"]
    row["score_abs_err"] = abs(row["he_score"] - p["score"])
    row["plain_p"] = p["p"]

    if row.get("he_stat") is not None:
        row["stat_abs_err"] = abs(row["he_stat"] - p["stat"])
        row["stat_rel_err"] = row["stat_abs_err"] / abs(p["stat"])

    if row.get("he_p") is not None:
        row["p_abs_err"] = abs(row["he_p"] - p["p"])
        row["p_decision"] = int(row["he_p"] < alpha)

    if "inv_df" in p:
        row["plain_inv_df"] = p["inv_df"]
        row["inv_df_abs_err"] = abs(row["he_inv_df"] - p["inv_df"])

    if "critical" in p and row.get("he_critical") is not None:
        row["plain_critical"] = p["critical"]
        row["critical_abs_err"] = abs(row["he_critical"] - p["critical"])

    row["plain_decision"] = int(p["score"] > 0.0)
    row["step_decision"] = int(row["he_step"] > 0.5)
    row["score_decision"] = int(row["he_score"] > 0.0)
    row["total_time"] = row["stats_time"] + row["decision_time"]

    return row


def summarize(raw_path, summary_path):

    groups = {}

    with raw_path.open() as file:
        for row in csv.DictReader(file):
            key = (row["case"], row["test"], row["config"], row["alpha"], row["target"])
            groups.setdefault(key, []).append(row)

    fields = [
        "case", "test", "config", "alpha", "target", "reps",
        "stats_time_mean", "decision_time_mean", "total_time_mean", "total_time_std",
        "bootstraps", "plain_stat", "stat_rel_err_max", "inv_df_abs_err_max",
        "critical_abs_err_max", "plain_score", "score_abs_err_max", "score_bound",
        "plain_p", "p_abs_err_max", "plain_decision", "step_match_rate",
        "score_decision_match_rate", "p_decision_match_rate", "step_min", "step_max",
    ]

    def col(rs, name):
        vals = [float(r[name]) for r in rs if r[name] not in ("", None)]
        return np.array(vals) if vals else np.array([np.nan])

    rows = []

    for key, rs in groups.items():
        plain_dec = int(rs[0]["plain_decision"])
        total = col(rs, "total_time")
        rows.append({
            **dict(zip(fields[:5], key)),
            "reps": len(rs),
            "stats_time_mean": col(rs, "stats_time").mean(),
            "decision_time_mean": col(rs, "decision_time").mean(),
            "total_time_mean": total.mean(),
            "total_time_std": total.std(ddof=1) if len(total) > 1 else 0.0,
            "bootstraps": int(np.median(col(rs, "stats_bootstraps") + col(rs, "decision_bootstraps"))),
            "plain_stat": rs[0]["plain_stat"],
            "stat_rel_err_max": col(rs, "stat_rel_err").max(),
            "inv_df_abs_err_max": col(rs, "inv_df_abs_err").max(),
            "critical_abs_err_max": col(rs, "critical_abs_err").max(),
            "plain_score": rs[0]["plain_score"],
            "score_abs_err_max": col(rs, "score_abs_err").max(),
            "score_bound": rs[0]["score_bound"],
            "plain_p": rs[0]["plain_p"],
            "p_abs_err_max": col(rs, "p_abs_err").max(),
            "plain_decision": plain_dec,
            "step_match_rate": np.mean([int(r["step_decision"]) == plain_dec for r in rs]),
            "score_decision_match_rate": np.mean([int(r["score_decision"]) == plain_dec for r in rs]),
            "p_decision_match_rate": (np.mean([int(r["p_decision"]) == plain_dec for r in rs if r["p_decision"] != ""])
                                      if any(r["p_decision"] != "" for r in rs) else ""),
            "step_min": col(rs, "he_step").min(),
            "step_max": col(rs, "he_step").max(),
        })

    with summary_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--cases", default="all")
    parser.add_argument("--tests", default="welch,f,z")
    parser.add_argument("--configs", default="all")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--score", choices=["ratio", "product"], default="ratio",
                        help="ratio: T^2 - c^2 (data-chosen bound, results/); "
                             "product: multiplied-out scores with public bounds (results/product/)")
    args = parser.parse_args()

    if args.out_dir is None:
        args.out_dir = RESULT_DIR if args.score == "ratio" else RESULT_DIR / "product"

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / "eval_raw.csv"
    summary_path = args.out_dir / "eval_summary.csv"

    if not args.summary_only:
        sys.path.append(str(PROJECT_ROOT / "experiments" / "invsqrt"))

        from engine.HEengine import HEengine
        from engine.HEdata import Message
        from operators.HEhypo_test import HEHypothesisTesting
        from hedap_optimizer import open_append

        engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")

        cases = load_cases()

        if args.cases != "all":
            cases = {k: v for k, v in cases.items() if k in args.cases.split(",")}

        done = set()

        if raw_path.is_file():
            with raw_path.open() as file:
                counts = {}
                for row in csv.DictReader(file):
                    k = (row["case"], row["test"], row["config"], row["rep"])
                    counts[k] = counts.get(k, 0) + 1
            done = {k for k, n in counts.items() if n == len(ALPHAS) * (len(TARGETS) if k[1] == "welch" else 1)}

        file, writer = open_append(raw_path, RAW_FIELDS)
        meta = {}

        for case_name, case in cases.items():
            x1, x2 = case["x1"], case["x2"]
            s1, s2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
            plain = {
                "welch": {a: welch_plain(x1, x2, a) for a in ALPHAS},
                "f": {a: f_plain(x1, x2, a) for a in ALPHAS},
                "z": {a: z_plain(x1, x2, a, s1, s2) for a in ALPHAS},
            }

            if args.score == "product":
                plain = {t: {a: {**v, "score": v["score_product"]} for a, v in d.items()} for t, d in plain.items()}

            for test in args.tests.split(","):
                if args.score == "ratio":
                    bound = score_bound([plain[test][a]["score"] for a in ALPHAS])
                    bounds = {a: float(bound) for a in ALPHAS}

                    # Valid sign input range: 1.5e-7 < |score / B| < 1.
                    ratios = [abs(plain[test][a]["score"]) / bound for a in ALPHAS]
                    assert 1e-5 < min(ratios) and max(ratios) < 1.0, (case_name, test, ratios)
                    meta[f"{case_name}/{test}"] = {"score_bound": bound, "n1": case["n1"], "n2": case["n2"], "R": case["R"]}
                else:
                    # Public bounds (Z: max(R^2 / V_public, z^2); Welch / F: computed in the decision).
                    bounds = {a: HEHypothesisTesting.z_public_bound(case["R"], s1, s2, case["n1"], case["n2"], a)
                              for a in ALPHAS}
                    meta[f"{case_name}/{test}"] = {"score": "product", "n1": case["n1"], "n2": case["n2"], "R": case["R"]}

                if test == "z":
                    configs = {"public_variance": None}
                elif test == "f" and args.score == "product":
                    configs = {"no_invsqrt": None}
                else:
                    names = INVSQRT_CONFIG_NAMES if args.configs == "all" else args.configs.split(",")
                    configs = {name: invsqrt_config(name, case_name, test) for name in names}

                for config_name, cfg in configs.items():
                    ht = HEHypothesisTesting(engine, cfg)

                    for rep in range(args.reps):
                        if (case_name, test, config_name, str(rep)) in done:
                            continue

                        # Encrypt both groups for this repetition (not timed).
                        c1 = engine.enc(Message(x1, engine.log_slots))
                        c2 = engine.enc(Message(x2, engine.log_slots))

                        if test == "welch":
                            rows = (run_welch(engine, ht, case, c1, c2, bounds) if args.score == "ratio"
                                    else run_welch_product(engine, ht, case, c1, c2))
                        elif test == "f":
                            rows = (run_f(engine, ht, case, c1, c2, bounds, plain["f"]) if args.score == "ratio"
                                    else run_f_product(engine, ht, case, c1, c2, plain["f"]))
                        else:
                            rows = run_z(engine, ht, case, c1, c2, bounds)

                        for row in rows:
                            row = complete(row, plain[test], row["alpha"])
                            writer.writerow({"case": case_name, "test": test, "config": config_name, "rep": rep,
                                             **{k: row.get(k, "") for k in RAW_FIELDS[4:]}})

                        file.flush()

                        r0 = rows[0]
                        print(f"{case_name} {test} {config_name} rep={rep} stats={r0['stats_time']:.1f}s "
                              f"decision={np.mean([r['decision_time'] for r in rows]):.1f}s "
                              f"score_abs_err={max(r['score_abs_err'] for r in rows):.2e} "
                              f"step_match={sum(r['step_decision'] == r['plain_decision'] for r in rows)}/{len(rows)}",
                              flush=True)

        file.close()

        meta_path = args.out_dir / "eval_meta.json"
        old_meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
        meta_path.write_text(json.dumps({**old_meta, **meta}, indent=2))

    rows = summarize(raw_path, summary_path)
    print(f"{len(rows)} summary rows, saved {summary_path}")


if __name__ == "__main__":
    main()
