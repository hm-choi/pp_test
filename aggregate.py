"""Summary tables of the four experiments from the result CSVs in results/.

Writes results/SUMMARY.md: the numbers reported in the paper and the
supplementary breakdowns (per significance level, per degrees of freedom).

Run (cwd = project root):
    python3 aggregate.py [--results results]
"""

import argparse
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "experiments"))

from common.cases import load_cases
from coeffs.invSqrt_coeffs import _INV_SQRT_MIN, _INV_SQRT_MAX

ALPHAS = ["0.01", "0.025", "0.05", "0.1"]
TESTS = ("welch", "f", "z")
f = float


def e(x):
    return f"{x:.2e}"


def t2(x):
    return f"{x:.2f}"


def ms(values, fmt=e):
    # Mean (standard deviation).
    sd = fmt(st.stdev(values)) if len(values) > 1 else "--"
    return f"{fmt(st.mean(values))} ({sd})"


def experiment1(rows, cases):

    rel = lambda r, k: f(r[k + "_abs_err"]) / abs(f(r["plain_" + k]))
    w = [r for r in rows if r["test"] == "welch"]
    fz = {t: [r for r in rows if r["test"] == t] for t in ("f", "z")}
    crit = [abs(f(r["he_critical"]) ** 2 - f(r["plain_critical"]) ** 2) / f(r["plain_critical"]) ** 2 for r in w]

    out = ["## Experiment 1: Accuracy (Section 5.4, Table 2)", "",
           "Relative error against the plaintext (SciPy) values, mean (standard deviation) over every",
           "workload, significance level and repetition. The critical value is compared as c^2.", "",
           "| Test | Statistic RE | 1/df RE | Critical RE | Score RE |", "|---|---|---|---|---|",
           f"| Welch's t | {ms([f(r['stat_rel_err']) for r in w])} "
           f"| {ms([f(r['inv_df_abs_err']) / f(r['plain_inv_df']) for r in w])} "
           f"| {ms(crit)} | {ms([rel(r, 'score') for r in w])} |",
           f"| F | -- | -- | Public | {ms([rel(r, 'score') for r in fz['f']])} |",
           f"| Z | {ms([f(r['stat_rel_err']) for r in fz['z']])} | -- | Public "
           f"| {ms([rel(r, 'score') for r in fz['z']])} |", "",
           "Decision agreement (step output > 0.5 against the plaintext decision), reported in Section 5.7:", "",
           "| Test | Configurations | Runs | Match | Reject / fail to reject |", "|---|---|---|---|---|"]

    for t in TESTS:
        rs = [r for r in rows if r["test"] == t]
        cfg = [r for r in rs if r["rep"] == "0"]
        rej = sum(r["plain_decision"] == "1" for r in cfg)
        match = sum(r["plain_decision"] == r["step_decision"] for r in rs)
        out.append(f"| {t} | {len(cfg)} | {len(rs)} | {match} | {rej} / {len(cfg) - rej} |")

    out += ["", "Workloads:", "", "| Workload | n1 / n2 | R | V |", "|---|---|---|---|"]

    for name, c in cases.items():
        v = np.var(c["x1"], ddof=1) / c["n1"] + np.var(c["x2"], ddof=1) / c["n2"]
        out.append(f"| {name} | {c['n1']} / {c['n2']} | {c['R']:.4g} | {v:.3g} |")

    return out


def experiment2(he, plain):

    out = ["## Experiment 2: Critical-boundary approximation (Section 5.5, Table 3)", "",
           "Squared two-sided t critical value, relative error against SciPy `t.ppf`; 65,536 df points",
           "log-spaced over [3, 10^4] in two ciphertexts, alpha = 0.05, mean (standard deviation) over the",
           "repetitions.", "",
           "| Degree | Direct-DF MaxRE | Direct-DF MRE | Inverse-DF MaxRE | Inverse-DF MRE |", "|---|---|---|---|---|"]

    g = [r for r in he if r["condition"] == "global" and r["alpha"] == "0.05"]

    for d in ["7", "15", "31", "63", "127"]:
        cells = []

        for m in ("df", "invdf"):
            rs = [r for r in g if r["method"] == m and r["degree"] == d]
            cells += [ms([f(r["max_re"]) for r in rs]), ms([f(r["mre"]) for r in rs])]

        out.append(f"| {d} | " + " | ".join(cells) + " |")

    out += ["", "Supplementary: plaintext MaxRE over [3, 10^4] for every significance level.", "",
            "| alpha | Degree | Direct-DF MaxRE | Inverse-DF MaxRE |", "|---|---|---|---|"]

    for a in ALPHAS:
        for d in ["7", "15", "31", "63", "127"]:
            p = {r["method"]: r for r in plain if r["range"] == "global" and r["alpha"] == a and r["degree"] == d}
            out.append(f"| {a} | {d} | {e(f(p['df']['max_re']))} | {e(f(p['invdf']['max_re']))} |")

    return out


def experiment3(rows):

    dfs = sorted({int(r["df"]) for r in rows})
    cnt = defaultdict(lambda: [0, 0])

    for r in rows:
        k = (r["alpha"], int(r["df"]))
        cnt[k][1] += 1
        cnt[k][0] += r["step_decision"] == r["plain_decision"]

    total = sum(v[0] for v in cnt.values())
    out = ["## Experiment 3: Boundary-sensitive decisions (Section 5.6)", "",
           f"Decision agreement {total}/{len(rows)}. Margins +-0.001, +-0.002, +-0.005, +-0.01, ten",
           "repetitions (80 trials per cell).", "",
           "| alpha | " + " | ".join(f"df {d}" for d in dfs) + " |", "|---|" + "---|" * len(dfs)]

    for a in ALPHAS:
        out.append(f"| {a} | " + " | ".join(f"{cnt[(a, d)][0]}/{cnt[(a, d)][1]}" for d in dfs) + " |")

    return out


def experiment4(rows, cases):

    by = defaultdict(list)

    for r in rows:
        by[(r["case"], r["test"])].append(r)

    out = ["## Experiment 4: Performance (Section 5.7, Table 4)", "",
           "alpha = 0.05, time (s) as mean (standard deviation) over the repetitions, encryption and",
           "decryption excluded. ApproxTCrit is included in the Welch total; Statistic is the Welch",
           "statistics part (Algorithms 1-2).", "",
           "| Workload | Welch total | Welch statistic | ApproxTCrit | Welch #BTS | F total | F #BTS | Z total | Z #BTS |",
           "|---|---|---|---|---|---|---|---|---|"]

    for name in cases:
        cells = [name]

        for t in TESTS:
            rs = by[(name, t)]
            bts = sorted({int(f(r["stats_bootstraps"])) + int(f(r["decision_bootstraps"])) for r in rs})
            cells.append(ms([f(r["total_time"]) for r in rs], t2))

            if t == "welch":
                cells += [ms([f(r["stats_time"]) for r in rs], t2), ms([f(r["critical_time"]) for r in rs], t2)]

            cells.append("/".join(map(str, bts)))

        out.append("| " + " | ".join(cells) + " |")

    return out


def read(path):

    with path.open() as file:
        return list(csv.DictReader(file))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=PROJECT_ROOT / "results")
    args = parser.parse_args()

    res = args.results
    cases = load_cases()

    out = ["# Summary of the experiments", "",
           "HEaaN FGb, log_slots = 15. invSqrt domain V in "
           f"[{_INV_SQRT_MIN:g}, {_INV_SQRT_MAX:g}] (degree-63 Chebyshev + 7 Newton iterations);",
           "ApproxTCrit degree 15. Generated by `aggregate.py`.", ""]
    out += experiment1(read(res / "experiment1" / "accuracy_raw.csv"), cases) + [""]
    out += experiment2(read(res / "experiment2" / "he_raw.csv"), read(res / "experiment2" / "plaintext.csv")) + [""]
    out += experiment3(read(res / "experiment3" / "boundary_decision_raw.csv")) + [""]
    out += experiment4(read(res / "experiment4" / "performance_raw.csv"), cases) + [""]

    (res / "SUMMARY.md").write_text("\n".join(out))
    print(f"written {res / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
