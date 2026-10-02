"""HE-DAP parameter search for HEApprox.invSqrt on HEaaN.

Port of Algorithm 2 (GetOptIter) and Algorithm 3 of
HE-DAP (Park et al., SAC'26, https://github.com/hm-choi/he_dap).

For each invSqrt domain, every input ciphertext level l in [l_BTS + 1, l_max],
Chebyshev degree 2^d - 1 (d = 4..9) and Pre-BTS indicator c, one run of
Chebyshev + i_max Newton iterations records the MRE and the accumulated
runtime after each iteration. The optimal iteration per (l, d, c) and the
optimal (d, c, i) per level are then chosen as in HE-DAP.

Domains (see domains.py):
- normalized : x = V / V_max in [v_min, 1] (v_min = 1e-3, 1e-4, 1e-5);
- raw        : x = V in [1e-3, V_max] with V_max from the public bounds of
               each dataset comparison (Welch V_max and F Var2_max).

Differences from the Lattigo reference code:
- degree is 2^d - 1 (pp_test coefficients) instead of 2^d - 2;
- decryption / MRE computation is excluded from the accumulated runtime;
- raw per-iteration results are saved, so the selection can be re-derived;
- the algorithmic error against 1/sqrt(dec(x)) is also recorded, and a
  supplementary selection with the algorithmic MaxRE is reported.

HEaaN-specific choices in HEApprox.invSqrt (they shape the search space):
- y in the Newton steps (in [1/sqrt(hi), 1/sqrt(lo)]) is refreshed with the
  extended bootstrap (input range 2^20, min level 4) instead of Lattigo's
  DoBootstrap(y, 4);
- x_half is bootstrapped only when its level is <= 5, so with c = 0 a low
  input level makes every Newton step end at level 4-5 and bootstrap y
  (this is the trade-off the Pre-BTS indicator c explores);
- Pre-BTS bootstraps x_cheb (in [-1, 1]) and derives x_half from it;
- the Chebyshev output y0 is refreshed with the regular bootstrap (its
  error for |y0| > 1 is refined by the following Newton iterations).

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/invsqrt/hedap_optimizer.py --domains raw
"""

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(Path(__file__).resolve().parent))

from engine.HEengine import HEengine
from engine.HEdata import Message
from operators.HEapprox import HEApprox

from domains import get_domains


D_MIN, D_MAX = 4, 9  # degree = 2^d - 1
I_MAX = 15
THETA = 1.0
DELTA = 1.0

# HEaaN FGb: l_max = 12, l_BTS = 3, l_afterBTS = 12
L_MAX = 12
L_BTS = 3
L_AFTER_BTS = 12

RESULT_DIR = Path(__file__).resolve().parent / "results"

RAW_FIELDS = [
    "domain", "method", "lo", "hi", "level", "degree", "pre_bts", "iteration",
    "mre", "max_re", "mre_alg", "max_re_alg", "time", "bootstrap_count",
]


def gen_data(start, middle, end, n):

    a = np.linspace(start, middle, n // 2)
    b = np.linspace(middle, end, n // 2)
    return np.concatenate((a, b))


def hedap_data(domain, n):

    lo, hi = domain["lo"], domain["hi"]
    return gen_data(lo, (lo + hi) / 2, hi, n)


def decrypt(engine, ctxt, size):

    msg = engine.dec(ctxt)
    return np.concatenate([np.array(m, dtype=np.complex128).real for m in msg])[:size]


def relative_error(engine, ctxt, ans, ans_alg):

    # ans: true 1/sqrt(x); ans_alg: 1/sqrt(dec(x)), excluding the input
    # encryption noise (algorithmic error only).
    y = decrypt(engine, ctxt, len(ans))
    re = np.abs(1.0 - y / ans)
    re_alg = np.abs(1.0 - y / ans_alg)
    return float(np.mean(re)), float(np.max(re)), float(np.mean(re_alg)), float(np.max(re_alg))


def threshold(m_min, hyper):

    # (floor(alpha) + hyper) * 10^ell  with  m_min = alpha * 10^ell
    ell = math.floor(math.log10(m_min))
    alpha = m_min / 10 ** ell
    return (math.floor(alpha) + hyper) * 10 ** ell


def select_iteration(mres, delta=DELTA):

    # Algorithm 2, lines 11-15 (1-based iteration count).
    finite = [m for m in mres if np.isfinite(m)]

    if not finite:
        return None

    mre_delta = threshold(min(finite), delta)

    for i, m in enumerate(mres):
        if np.isfinite(m) and m <= mre_delta:
            return i + 1

    return None


def map_input(engine, ctxt, domain):

    # x -> (x_half, x_cheb) as prepared by HEHypothesisTesting (one level).
    lo, hi = domain["lo"], domain["hi"]

    x_half = engine.mult(ctxt, 0.5)
    x_cheb = engine.mult(ctxt, 2.0 / (hi - lo))
    x_cheb = engine.sub(x_cheb, (hi + lo) / (hi - lo))

    return x_half, x_cheb


def get_opt_iter(approx, ctxt, ans, domain, log_degree, pre_bts, i_max=I_MAX):

    # Algorithm 2 (GetOptIter): one Chebyshev + i_max Newton iterations,
    # recording MRE and accumulated runtime after every iteration.
    engine = approx.engine
    dom = (domain["lo"], domain["hi"])

    x_dec = decrypt(engine, ctxt, len(ans))
    ans_alg = np.maximum(x_dec, np.finfo(float).tiny) ** -0.5

    # One-off plaintext work (coefficients) outside the timer.
    approx.inv_sqrt_coeffs(log_degree, dom, domain["method"])
    approx.reset_bootstrap_count()

    start = time.perf_counter()

    x_half, x_cheb = map_input(engine, ctxt, domain)
    x_half, y = approx.invSqrt_init(x_half, x_cheb, log_degree, dom, domain["method"], pre_bts)

    elapsed = time.perf_counter() - start

    rows = []

    for i in range(i_max):
        start = time.perf_counter()
        y = approx.newton_step(x_half, y)
        elapsed += time.perf_counter() - start

        mre, max_re, mre_alg, max_re_alg = relative_error(engine, y, ans, ans_alg)

        rows.append({
            "iteration": i + 1,
            "mre": mre,
            "max_re": max_re,
            "mre_alg": mre_alg,
            "max_re_alg": max_re_alg,
            "time": elapsed,
            "bootstrap_count": approx.bootstrap_count(),
        })

    return rows


def open_append(path, fields):

    # Append-mode CSV writer; drops a partial last line left by a hard kill.
    if path.is_file():
        data = path.read_bytes()

        if data and not data.endswith(b"\n"):
            path.write_bytes(data[:data.rfind(b"\n") + 1])

    new_file = not path.is_file() or path.stat().st_size == 0
    file = path.open("a", newline="")
    writer = csv.DictWriter(file, fieldnames=fields)

    if new_file:
        writer.writeheader()

    return file, writer


def load_raw(raw_path):

    done = {}

    if not raw_path.is_file():
        return done

    with raw_path.open() as file:
        for row in csv.DictReader(file):
            key = (row["domain"], int(row["level"]), int(row["degree"]), int(row["pre_bts"]))
            done.setdefault(key, []).append({
                "iteration": int(row["iteration"]),
                "mre": float(row["mre"]),
                "max_re": float(row["max_re"]),
                "mre_alg": float(row["mre_alg"]),
                "max_re_alg": float(row["max_re_alg"]),
                "time": float(row["time"]),
                "bootstrap_count": int(row["bootstrap_count"]),
            })

    # Keep only complete runs (last occurrence per iteration, in case a run
    # was interrupted and repeated).
    ret = {}

    for key, rows in done.items():
        by_iter = {r["iteration"]: r for r in rows}

        if sorted(by_iter) == list(range(1, I_MAX + 1)):
            ret[key] = [by_iter[i] for i in range(1, I_MAX + 1)]

    return ret


def cases(level):

    # Algorithm 3, lines 9 and 13.
    ret = []

    if level >= L_BTS + 3:
        ret.append(0)

    if level <= L_AFTER_BTS - 2:
        ret.append(1)

    return ret


def level_tuples(raw, name, level, metric):

    # D[l] of Algorithm 3: GetOptIter result for every (d, c) at level l.
    tuples = []

    for d in range(D_MIN, D_MAX + 1):
        for c in cases(level):
            rows = raw.get((name, level, 2 ** d - 1, c))

            if rows is None:
                continue

            i = select_iteration([r[metric] for r in rows])

            if i is None:
                continue

            r = rows[i - 1]
            tuples.append({
                "degree": 2 ** d - 1, "pre_bts": c, "iteration": i,
                "mre": r["mre"], "max_re": r["max_re"],
                "mre_alg": r["mre_alg"], "max_re_alg": r["max_re_alg"], "time": r["time"],
                "bootstrap_count": r["bootstrap_count"],
            })

    return tuples


def optimize(raw, domains):

    # Algorithm 3, lines 19-27, per domain and level.
    # "accuracy" / "speed" follow HE-DAP (MRE criterion). "accuracy_maxre_alg"
    # is a supplementary variant using the algorithmic MaxRE (against
    # 1/sqrt(dec(x))) instead of MRE in Algorithms 2 and 3, which also checks
    # convergence for small inputs.
    result = {}

    for domain in domains:
        name = domain["name"]
        per_level = {}

        for level in range(L_MAX, L_BTS, -1):
            tuples = level_tuples(raw, name, level, "mre")

            if not tuples:
                continue

            m_theta = threshold(min(t["mre"] for t in tuples), THETA)

            u1 = min((t for t in tuples if t["mre"] <= m_theta), key=lambda t: t["time"])
            u2 = min(tuples, key=lambda t: t["time"])

            tuples_maxre = level_tuples(raw, name, level, "max_re_alg") or tuples
            m_theta_maxre = threshold(min(t["max_re_alg"] for t in tuples_maxre), THETA)
            u1_maxre = min(
                (t for t in tuples_maxre if t["max_re_alg"] <= m_theta_maxre),
                key=lambda t: t["time"],
            )

            per_level[str(level)] = {
                "accuracy": u1,   # u1: fastest under the MRE bound
                "speed": u2,      # u2: fastest regardless of MRE
                "accuracy_maxre_alg": u1_maxre,
                "m_theta": m_theta,
                "m_theta_maxre_alg": m_theta_maxre,
                "candidates": tuples,
                "candidates_maxre_alg": tuples_maxre,
            }

        result[name] = {**domain, "levels": per_level}

    return result


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--domains", default="raw", help="normalized, raw, all or comma-separated names")
    parser.add_argument("--levels", type=int, nargs="*", default=list(range(L_MAX, L_BTS, -1)))
    parser.add_argument("--log-degrees", type=int, nargs="*", default=list(range(D_MIN, D_MAX + 1)))
    parser.add_argument("--select-only", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()

    domains = get_domains(args.domains)
    out_dir = args.out_dir or RESULT_DIR / args.domains.replace(",", "_")
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_path = out_dir / "hedap_raw.csv"
    opt_path = out_dir / "hedap_optimal.json"

    raw = load_raw(raw_path)

    if not args.select_only:
        engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")
        approx = HEApprox(engine)

        file, writer = open_append(raw_path, RAW_FIELDS)

        for domain in domains:
            test = hedap_data(domain, engine.num_slots)
            ans = test ** -0.5

            for level in args.levels:
                for d in args.log_degrees:
                    for c in cases(level):
                        key = (domain["name"], level, 2 ** d - 1, c)

                        if key in raw:
                            continue

                        ctxt = engine.enc(Message(test, engine.log_slots), level=level)
                        rows = get_opt_iter(approx, ctxt, ans, domain, d, bool(c))

                        for r in rows:
                            writer.writerow({
                                "domain": domain["name"], "method": domain["method"],
                                "lo": domain["lo"], "hi": domain["hi"], "level": level,
                                "degree": 2 ** d - 1, "pre_bts": c, **r,
                            })

                        file.flush()
                        raw[key] = rows

                        i = select_iteration([r["mre"] for r in rows])
                        best = rows[i - 1] if i else None
                        print(
                            f"{domain['name']} l={level:2d} deg={2 ** d - 1:3d} c={c} "
                            f"-> iter={i} MRE={best['mre'] if best else float('nan'):.3e} "
                            f"time={best['time'] if best else float('nan'):.2f}s "
                            f"(i_max: {rows[-1]['time']:.2f}s)",
                            flush=True,
                        )

        file.close()

    opt_path.write_text(json.dumps(optimize(load_raw(raw_path), domains), indent=2))
    print(f"saved {opt_path}")


if __name__ == "__main__":
    main()
