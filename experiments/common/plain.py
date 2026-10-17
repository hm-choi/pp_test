"""Plaintext (SciPy) reference values of the Welch, F and Z tests."""

import numpy as np
from scipy import stats


def welch_plain(x1, x2, alpha):

    n1, n2 = len(x1), len(x2)
    a1, a2 = np.var(x1, ddof=1) / n1, np.var(x2, ddof=1) / n2
    v = a1 + a2
    t2 = (np.mean(x1) - np.mean(x2)) ** 2 / v
    df = v ** 2 / (a1 ** 2 / (n1 - 1) + a2 ** 2 / (n2 - 1))
    crit = stats.t.ppf(1.0 - alpha / 2.0, df)
    p = 2.0 * stats.t.sf(np.sqrt(t2), df)

    return {"stat": t2, "inv_df": 1.0 / df, "critical": crit, "p": p,
            "score": t2 - crit ** 2}


def f_plain(x1, x2, alpha):

    n1, n2 = len(x1), len(x2)
    v1, v2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    fl = stats.f.ppf(alpha / 2.0, n1 - 1, n2 - 1)
    fu = stats.f.ppf(1.0 - alpha / 2.0, n1 - 1, n2 - 1)

    return {"fl": fl, "fu": fu, "score": (v1 - fl * v2) * (v1 - fu * v2)}


def z_plain(x1, x2, alpha, s1, s2):

    n1, n2 = len(x1), len(x2)
    z = (np.mean(x1) - np.mean(x2)) / np.sqrt(s1 / n1 + s2 / n2)
    crit = stats.norm.ppf(1.0 - alpha / 2.0)
    p = 2.0 * stats.norm.sf(abs(z))

    return {"stat": z, "critical": crit, "score": z ** 2 - crit ** 2, "p": p}
