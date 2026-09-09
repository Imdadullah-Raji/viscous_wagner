#!/usr/bin/env python3
"""The bifurcation point, located from directly measured linear growth rates.

sigma is the real part of the leading eigenvalue of the linearisation about the
STEADY base flow.  It is measured by perturbing that base flow slightly and
fitting the amplitude of the f=0.600 mode -- not by extrapolating any property of
the saturated state, which is what failed over 14-24 deg (FINDINGS.md 17-18).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "analysis", "figs", "eigenvalue_onset.png")

# AoA: (sigma, R2, how)
S = {11.0: (-0.07093, 0.99999, "decay, 2.2 decades"),
     11.5: (+0.00674, 1.00000, "growth, 63 t span"),
     12.0: (+0.10970, 0.95,    "growth, early linear window")}

a = np.array(sorted(S)); s = np.array([S[k][0] for k in sorted(S)])
k, b = np.polyfit(a, s, 1)
ac_fit = -b / k
ac_br = 11.0 + 0.5 * (-S[11.0][0]) / (S[11.5][0] - S[11.0][0])

fig, ax = plt.subplots(figsize=(8.6, 5.4))
xs = np.linspace(10.8, 12.2, 100)
ax.axhspan(-0.09, 0, color="#1f77b4", alpha=0.06)
ax.axhspan(0, 0.13, color="#d62728", alpha=0.06)
ax.text(10.85, -0.075, "STABLE\n(no shedding)", fontsize=9, color="#1f77b4", va="center")
ax.text(10.85, 0.105, "UNSTABLE\n(sheds)", fontsize=9, color="#d62728", va="center")
ax.plot(xs, k * xs + b, color="0.4", lw=1.3, ls="--",
        label=fr"linear fit: $\alpha_c$={ac_fit:.2f}$^\circ$, "
              fr"$d\sigma/d\alpha$={k:.3f}/deg")
ax.axhline(0, color="k", lw=0.9)
for x in sorted(S):
    sig, r2, how = S[x]
    ax.plot([x], [sig], "o", ms=9, color=("#1f77b4" if sig < 0 else "#d62728"), zorder=5)
    dx, dy, ha = (10, 8, "left") if x != 11.5 else (12, -10, "left")
    ax.annotate(f"{sig:+.5f}\n{how}", (x, sig), textcoords="offset points",
                xytext=(dx, dy), fontsize=8.5, ha=ha,
                va="bottom" if dy > 0 else "top")
ax.plot([ac_br], [0], marker="v", ms=11, color="k", clip_on=False, zorder=6)
ax.annotate(fr"$\alpha_c$ = {ac_br:.2f}$^\circ$", (ac_br, 0),
            textcoords="offset points", xytext=(-14, 26), ha="right", fontsize=12,
            bbox=dict(fc="w", ec="0.6"))
ax.set_xlim(10.8, 12.2); ax.set_ylim(-0.09, 0.13)
ax.set_xlabel(r"angle of attack  $\alpha$  (deg)")
ax.set_ylabel(r"linear growth rate  $\sigma$   (leading eigenvalue, mode $f=0.600$)")
ax.set_title("Re=500 static NACA0012: shedding onset from measured eigenvalues\n"
             "bracketed -- stable at 11.0$^\\circ$, unstable at 11.5$^\\circ$")
ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
os.makedirs(os.path.dirname(FIG), exist_ok=True)
fig.savefig(FIG, dpi=150)
print("wrote", FIG)
