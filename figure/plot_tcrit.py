"""Figure 2: two-sided t critical value against nu and against u = 1/nu.

Run (cwd = project root):
    python3 figure/plot_tcrit.py [out_dir]   ->  tcrit_vs_df_invdf.pdf / .png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
ALPHAS = [0.01, 0.025, 0.05, 0.1]
# The curves never cross and the legend follows their vertical order.
COLORS = {0.01: "#2a78d6", 0.025: "#eb6834", 0.05: "#1baf7a", 0.1: "#e34948"}
INK, GRID = "#1a1a19", "#e4e3dd"
LABEL_SIZE = 12

plt.rcParams.update({"font.size": 10, "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "axes.edgecolor": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.8})


def style(ax):
    ax.set_yscale("log")
    ax.set_ylim(1, 1e2)
    ax.minorticks_off()
    ax.set_ylabel(r"$t$-critical value", fontsize=LABEL_SIZE)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_zorder(5)


def legend(ax, loc):
    ax.legend(loc=loc, fontsize=9, frameon=True, fancybox=True, framealpha=1.0, edgecolor="#cccccc",
              labelcolor=INK, handlelength=1.8)


def plot_nu(ax):
    df = np.geomspace(1.0, 1e4, 2000)
    for a in ALPHAS:
        ax.plot(df, stats.t.ppf(1.0 - a / 2.0, df), color=COLORS[a], lw=1.6, solid_capstyle="round",
                label=rf"$\alpha={a:g}$")
    ax.set_xscale("log")
    ax.set_xlim(1, 1e4)
    ax.set_xlabel(r"Degrees of freedom $\nu$", fontsize=LABEL_SIZE)
    style(ax)


def plot_u(ax):
    u = np.linspace(0.0, 1.0, 2001)
    nu = np.divide(1.0, u, out=np.full_like(u, np.inf), where=u > 0)
    for a in ALPHAS:
        ax.plot(u, stats.t.ppf(1.0 - a / 2.0, nu), color=COLORS[a], lw=1.6, solid_capstyle="round",
                label=rf"$\alpha={a:g}$")
    ax.set_xlim(0, 1)
    ax.set_xlabel(r"Inverse degrees of freedom $u=1/\nu$", fontsize=LABEL_SIZE)
    style(ax)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.4), sharey=True)
plot_nu(ax1)
plot_u(ax2)
ax2.set_ylabel("")
ax2.tick_params(labelleft=False)
legend(ax1, "upper right")
legend(ax2, "upper left")
fig.tight_layout(w_pad=1.0)
fig.savefig(OUT / "tcrit_vs_df_invdf.pdf", bbox_inches="tight")
fig.savefig(OUT / "tcrit_vs_df_invdf.png", dpi=600, bbox_inches="tight")
print("saved", OUT)
