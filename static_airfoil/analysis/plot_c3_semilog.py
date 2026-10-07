"""log|c3| against Re (linear Re axis), and against log Re for comparison.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/plot_c3_semilog.py [--step 1]
Reads figs/legendre_vs_re_11re_step<step>/legendre_vs_re.csv (legendre_vs_re.py).
A straight line on the left panel would mean exponential growth in Re; a straight
line on the right one, a power law.  The dashed curve is the Re 200-1000 power law.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lasso_re1000 import ROOT

RANGES = {6.0: "#2a78d6", 8.0: "#eb6834"}

ap = argparse.ArgumentParser()
ap.add_argument("--step", type=float, default=1.0)
step = ap.parse_args().step
d = ROOT / "analysis" / "figs" / f"legendre_vs_re_11re_step{step:g}"
t = pd.read_csv(d / "legendre_vs_re.csv")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8))
g = np.linspace(100, 2000, 400)
for amax, col in RANGES.items():
    s = t[t.amax == amax].sort_values("re")
    re, y = s.re.to_numpy(), np.log(-s.c3.to_numpy())
    m = (re >= 200) & (re <= 1000)
    p, c = np.polyfit(np.log(re[m]), y[m], 1)
    for ax, x, xg in ((a1, re, g), (a2, np.log(re), np.log(g))):
        ax.plot(x, y, "o", color=col, ms=7, label=f"fit range 0-{amax:g}°")
        ax.plot(xg, p * np.log(g) + c, "--", color=col, lw=1.2,
                label=rf"$\ln|c_3| = {p:.3f}\,\ln Re {c:+.3f}$  (Re 200-1000)")
a1.set(xlabel="Re", ylabel=r"$\ln|c_3|$")
a1.set_title(r"$\ln|c_3|$ vs Re", loc="left")
a2.set(xlabel=r"$\ln Re$", ylabel=r"$\ln|c_3|$")
a2.set_title(r"$\ln|c_3|$ vs $\ln Re$ (slope = power-law exponent)", loc="left")
sec = a2.secondary_xaxis("top", functions=(np.exp, np.log))
sec.set_xticks([100, 200, 500, 1000, 2000]); sec.set_xlabel("Re")
for ax in (a1, a2):
    ax.grid(alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle(f"Odd-Legendre cubic coefficient (P1-P9 fit, {step:g}° grid, all CONVERGED)",
             x=0.01, ha="left")
fig.tight_layout()
fig.savefig(d / "ln_c3_vs_re.png", dpi=200)
print(d / "ln_c3_vs_re.png")
