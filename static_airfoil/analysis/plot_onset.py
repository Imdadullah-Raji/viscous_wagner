#!/usr/bin/env python3
"""Shedding onset from saturated amplitudes -- and why it is NOT pinned.

Near a supercritical Hopf the saturated amplitude should obey the normal form
A^2 = k (alpha - alpha_c), i.e. A ~ (alpha - alpha_c)^(1/2).  Fitting that over
16-24 deg gave alpha_c = 14.83 and predicted NO shedding at 14 deg.  AoA 14
sheds, saturating at A = 0.03413 (FINDINGS.md section 17).

The failure is a wrong functional form, not a stray point: the fitted alpha_c
rises monotonically with the highest angle included, which a correct model would
never do.  This figure shows that directly -- several linear fits that disagree,
against a free-exponent fit that describes all the data with p ~ 0.79 rather than
the 0.5 the normal form requires.

Period-2 cases are excluded throughout: the scaling describes a single Hopf mode
and says nothing once the orbit has doubled.
"""
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = [(14, 60.5), (16, 60.5), (18, 60.5), (20, 60.5), (22, 80.5), (24, 60.5)]
FIG = os.path.join(ROOT, "analysis", "figs", "onset_amplitude.png")


def carrier_amplitude(aoa, t_conv):
    out = subprocess.run(
        [sys.executable, os.path.join(ROOT, "analysis", "growth_rate.py"),
         os.path.join(ROOT, f"re500_aoa{aoa}"), "--tmin", str(t_conv)],
        capture_output=True, text=True).stdout
    d = dict(l.split("=", 1) for l in out.strip().splitlines() if "=" in l)
    return float(d["A_LAST"]), float(d["SIGMA"])


def main():
    a, A = [], []
    for aoa, tc in CASES:
        try:
            amp, sig = carrier_amplitude(aoa, tc)
        except KeyError:
            print(f"AoA {aoa}: no converged record yet, skipping"); continue
        a.append(aoa); A.append(amp)
        print(f"AoA {aoa}: A={amp:.6f}  A^2={amp**2:.7f}  sigma={sig:+.2e}")
    a = np.array(a, float); A = np.array(A); y = A ** 2

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.2, 5.4))

    # ---- left: the linear fits disagree, which is the whole point ----------
    subsets = [("16-24 (the original fit)", a >= 16, "#d62728"),
               ("14-20", a <= 20, "#1f77b4"),
               ("all six", np.ones(len(a), bool), "#2ca02c")]
    xs = np.linspace(11.5, 25, 200)
    for lab, m, col in subsets:
        if m.sum() < 2:
            continue
        k, c = np.polyfit(a[m], y[m], 1)
        ac = -c / k
        ax.plot(xs, k * xs + c, color=col, lw=1.4, alpha=0.9,
                label=fr"{lab}: $\alpha_c$={ac:.2f}$^\circ$")
        ax.plot([ac], [0], marker="v", color=col, ms=9, clip_on=False)
    ax.plot(a, y, "o", color="k", ms=7, zorder=5, label=r"measured $A^2$")
    for x, v in zip(a, y):
        ax.annotate(f"{int(x)}$^\\circ$", (x, v), textcoords="offset points",
                    xytext=(7, -3), fontsize=9)
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlim(11.5, 25); ax.set_ylim(-0.003, 0.040)
    ax.set_xlabel(r"$\alpha$  (deg)"); ax.set_ylabel(r"$A(f_{carrier})^2$")
    ax.set_title("Linear fits disagree with each other\n"
                 r"($\alpha_c$ rises with the highest angle included"
                 " $\\Rightarrow$ wrong form)")
    ax.grid(alpha=0.3); ax.legend(loc="upper left", fontsize=9)

    # ---- right: free exponent vs the normal form's p = 1/2 -----------------
    def prof(ac):
        x = np.log(a - ac)
        p, c = np.polyfit(x, np.log(A), 1)
        pred = np.exp(c) * (a - ac) ** p
        return p, np.sqrt(((A - pred) ** 2).mean())
    acs = np.arange(9.0, 13.9, 0.05)
    ps, rms = np.array([prof(v) for v in acs]).T
    bx.plot(acs, rms, color="#1f77b4", lw=1.6, label="free exponent $p$")
    # p forced to 1/2
    rms_half = []
    for ac in acs:
        x = a - ac
        C = (A * np.sqrt(x)).sum() / x.sum()
        rms_half.append(np.sqrt(((A - C * np.sqrt(x)) ** 2).mean()))
    bx.plot(acs, rms_half, color="#d62728", lw=1.6, ls="--",
            label=r"normal form, $p=1/2$ imposed")
    i = int(np.argmin(rms))
    bx.plot([acs[i]], [rms[i]], "o", color="#1f77b4", ms=8)
    bx.annotate(fr"best: $\alpha_c$={acs[i]:.2f}$^\circ$, $p$={ps[i]:.2f}",
                (acs[i], rms[i]), textcoords="offset points", xytext=(-12, 18),
                fontsize=10, ha="right")
    bx.set_xlabel(r"assumed $\alpha_c$  (deg)")
    bx.set_ylabel("rms residual in $A$")
    bx.set_title(r"$p$ and $\alpha_c$ trade off; $p=1/2$ is disfavoured"
                 "\n(shallow minimum $\\Rightarrow$ onset not pinned)")
    bx.grid(alpha=0.3); bx.legend(fontsize=9)
    ax2 = bx.twinx()
    ax2.plot(acs, ps, color="0.55", lw=1.0, ls=":")
    ax2.set_ylabel("fitted exponent $p$  (dotted)", color="0.45")
    ax2.axhline(0.5, color="0.75", lw=0.8, ls=":")

    fig.suptitle("Re=500 static NACA0012: the shedding onset is BELOW 14 deg "
                 "and is not yet determined", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=150)
    print("wrote", FIG)


if __name__ == "__main__":
    main()
