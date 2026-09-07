#!/usr/bin/env python3
"""Steady base-flow lift curve, and what the shedding adds on top of it.

Left:  C_L of the steady (SIMPLE) base flow against the time-accurate mean of
       the saturated limit cycle, with the instantaneous range of the cycle
       shaded.  The shaded band never reaches down to the steady curve.
Right: mean-flow distortion C_L,mean - C_L,steady against the oscillation
       amplitude A.  Stuart-Landau says this is O(A^2); it measures O(A^0.88).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "analysis", "figs", "steady_vs_mean.png")

# AoA: C_L steady (SIMPLE);  * = SIMPLE limit-cycled, no fixed point found
STEADY = {0: 0.000000, 2: 0.112979, 4: 0.218670, 6: 0.309705, 8: 0.379699,
          10: 0.426072, 12: 0.452294, 14: 0.465639, 16: 0.472710,
          18: 0.477540, 20: 0.482076, 22: 0.494523, 24: 0.505358}
NOT_CONVERGED = {22, 24}
# AoA: (C_L mean, C_L p2p, A(f_carrier)) from the saturated limit cycle
UNSTEADY = {14: (0.551601, 0.068222, 0.034050), 16: (0.630309, 0.139020, 0.069449),
            18: (0.711513, 0.217655, 0.108713), 20: (0.786700, 0.281306, 0.140328),
            22: (0.839900, 0.328268, 0.161594), 24: (0.879000, 0.396700, 0.188831)}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.2, 5.4))

a = np.array(sorted(STEADY))
cl = np.array([STEADY[x] for x in a])
ok = np.array([x not in NOT_CONVERGED for x in a])
ax.plot(a[ok], cl[ok], "o-", color="#1f77b4", lw=1.8, ms=6, label="steady base flow (SIMPLE)")
ax.plot(a[~ok], cl[~ok], "o", mfc="none", color="#1f77b4", ms=7,
        label="SIMPLE limit-cycled (no fixed point)")
ax.plot([0, 6], [0, 6 * 0.1097], ls=":", color="0.5", lw=1.2,
        label=r"inviscid $2\pi\alpha$")

au = np.array(sorted(UNSTEADY))
cm = np.array([UNSTEADY[x][0] for x in au])
p2 = np.array([UNSTEADY[x][1] for x in au])
ax.fill_between(au, cm - p2 / 2, cm + p2 / 2, color="#d62728", alpha=0.18,
                label="limit-cycle range")
ax.plot(au, cm, "s-", color="#d62728", lw=1.8, ms=6, label="time-accurate mean")
ax.set_xlabel(r"$\alpha$  (deg)"); ax.set_ylabel(r"$C_L$")
ax.set_title("Re=500 NACA0012: the base flow stalls; the shedding does not\n"
             "the cycle never descends to the steady value")
ax.grid(alpha=0.3); ax.legend(fontsize=8.5, loc="upper left")

A = np.array([UNSTEADY[x][2] for x in au])
sh = cm - np.array([STEADY[x] for x in au])
bx.loglog(A, sh, "o", color="k", ms=7, zorder=5, label="measured")
for x, xa, ys in zip(au, A, sh):
    bx.annotate(f"{x}$^\\circ$", (xa, ys), textcoords="offset points",
                xytext=(7, -4), fontsize=9)
xs = np.linspace(A.min() * 0.8, A.max() * 1.2, 50)
q, c = np.polyfit(np.log(A), np.log(sh), 1)
bx.loglog(xs, np.exp(c) * xs ** q, color="#1f77b4", lw=1.6,
          label=fr"fit: $\propto A^{{{q:.2f}}}$")
c2 = (sh / A ** 2).mean()
bx.loglog(xs, c2 * xs ** 2, ls="--", color="#d62728", lw=1.6,
          label=r"Stuart-Landau: $\propto A^{2}$")
bx.set_xlabel(r"oscillation amplitude  $A(f_{carrier})$")
bx.set_ylabel(r"$C_{L,mean}-C_{L,steady}$")
bx.set_title("Mean-flow distortion is LINEAR in the amplitude,\n"
             "not quadratic as the weakly-nonlinear theory requires")
bx.grid(alpha=0.3, which="both"); bx.legend(fontsize=9)

fig.tight_layout()
os.makedirs(os.path.dirname(FIG), exist_ok=True)
fig.savefig(FIG, dpi=150)
print("wrote", FIG)
