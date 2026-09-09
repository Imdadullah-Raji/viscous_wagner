#!/usr/bin/env python3
"""Force coefficients against angle of attack, with both bifurcations marked.

  analysis/figs/coeffs_vs_aoa.png     three panels: C_L, C_D, C_L/C_D
  analysis/figs/cl_vs_aoa.png         the C_L panel alone, for use as a figure
  analysis/figs/cd_vs_aoa.png
  analysis/figs/ld_vs_aoa.png

Each panel carries three curves, and the distinction between them is the point
of the figure:

  * the STEADY base flow (SIMPLE).  A fixed point of the steady equations,
    whether or not that fixed point is stable -- above ~11.5 deg it is not, so
    past there this curve is a solution nobody would ever observe.  Past 20 deg
    SIMPLE stops converging at all (the physical instability leaks into the
    iteration), and those points are drawn open.
  * the TIME-ACCURATE MEAN of the saturated limit cycle.
  * the instantaneous MIN and MAX over a whole number of orbits, dotted.

The gap between the dotted minimum and the steady curve is the campaign's
central geometric fact: **the limit cycle never descends to the fixed point**.
The instantaneous C_L stays above the steady value at every angle and every
phase, by 0.44-0.94 of a peak-to-peak, so the cycle is displaced bodily from
the fixed point rather than orbiting it -- which is why every weakly-nonlinear
(Stuart-Landau) scaling fitted to these saturated states came out wrong, and
why the onset had to be found from eigenvalues instead.  See FINDINGS.md 19.

Min/max are taken over an integer number of FUNDAMENTAL periods, not over
whatever the record happened to end on: past 24 deg the orbit is two shedding
cycles with unequal peaks, so a half-orbit of extra record moves the maximum.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign_data import (RC, ROOT, ANGLES, ALPHA_HOPF, ALPHA_PD_LO,   # noqa
                           ALPHA_PD_HI, STEADY_CL, STEADY_CD,
                           STEADY_NOT_CONVERGED, AOA11, all_stats)

FIGS = os.path.join(ROOT, "analysis", "figs")

C_MEAN = "#0072B2"      # time-accurate mean
C_ENV = "#0072B2"       # envelope (same hue, dotted)
C_STEADY = "#D55E00"    # steady base flow
C_HOPF = "#009E73"      # primary bifurcation
C_PD = "#CC79A7"        # period doubling


def mark_bifurcations(ax, label=True):
    """The two measured bifurcations, drawn as what each one actually is: the
    Hopf is a point with an uncertainty band, the period doubling is a bracket
    that was never narrowed below one degree.

    Call this LAST -- it opens headroom at the top of the axes for the labels,
    which would otherwise land on the data."""
    ax.axvspan(ALPHA_HOPF - 0.1, ALPHA_HOPF + 0.1, color=C_HOPF, alpha=0.18,
               lw=0, zorder=0)
    ax.axvline(ALPHA_HOPF, color=C_HOPF, lw=1.3, ls="--", zorder=1)
    ax.axvspan(ALPHA_PD_LO, ALPHA_PD_HI, color=C_PD, alpha=0.20, lw=0, zorder=0)
    if not label:
        return
    y0, y1 = ax.get_ylim()
    y1 = y1 + 0.24 * (y1 - y0)
    ax.set_ylim(y0, y1)
    ax.text(ALPHA_HOPF + 0.5, y1 - 0.02 * (y1 - y0),
            "Hopf\n" rf"$\alpha_c = {ALPHA_HOPF}^\circ\pm0.1$",
            color=C_HOPF, fontsize=9, va="top", ha="left", linespacing=1.35)
    ax.text(0.5 * (ALPHA_PD_LO + ALPHA_PD_HI), y1 - 0.02 * (y1 - y0),
            "period doubling\n" rf"$(24^\circ, 25^\circ]$",
            color=C_PD, fontsize=9, va="top", ha="center", linespacing=1.35)


def panel(ax, D, key, steady, ylabel, title, envelope=True):
    a = np.array(ANGLES, float)
    mean = np.array([D[i][f"{key}_mean"] for i in ANGLES])
    lo = np.array([D[i][f"{key}_min"] for i in ANGLES])
    hi = np.array([D[i][f"{key}_max"] for i in ANGLES])

    # --- steady base flow -------------------------------------------------
    if steady is not None:
        sa = np.array(sorted(steady), float)
        sv = np.array([steady[int(x)] for x in sa])
        ok = np.array([int(x) not in STEADY_NOT_CONVERGED for x in sa])
        ax.plot(sa, sv, "-", color=C_STEADY, lw=1.5, zorder=3)
        ax.plot(sa[ok], sv[ok], "o", color=C_STEADY, ms=4.5, zorder=4)
        ax.plot(sa[~ok], sv[~ok], "o", mfc="white", mec=C_STEADY, mew=1.3,
                ms=5.5, zorder=4)

    # --- limit-cycle envelope --------------------------------------------
    if envelope:
        ax.fill_between(a, lo, hi, color=C_ENV, alpha=0.13, lw=0, zorder=2)
        ax.plot(a, hi, ":", color=C_ENV, lw=1.4, zorder=3)
        ax.plot(a, lo, ":", color=C_ENV, lw=1.4, zorder=3)
    ax.plot(a, mean, "-", color=C_MEAN, lw=1.7, zorder=5)
    # period-2 angles get a filled square, period-1 a filled circle
    for i in ANGLES:
        m = "s" if D[i]["order"] >= 2 else "o"
        ax.plot([i], [D[i][f"{key}_mean"]], m, color=C_MEAN, ms=5.5, zorder=6)

    ax.set_xlim(9.0, 31.0)
    ax.set_xlabel(r"angle of attack  $\alpha$   (deg)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    ax.grid(alpha=0.22, lw=0.5)
    return mean, lo, hi


def add_aoa11(ax, value):
    ax.plot([AOA11["aoa"]], [value], "D", color="#000000", ms=5, zorder=7,
            mfc="white", mec="#000000", mew=1.2)


def legend_handles(with_steady=True):
    h = [Line2D([], [], color=C_MEAN, lw=1.7, marker="o", ms=5.5,
                label="time-accurate mean (period-1)"),
         Line2D([], [], color=C_MEAN, lw=1.7, marker="s", ms=5.5,
                label="time-accurate mean (period-2)"),
         Line2D([], [], color=C_ENV, lw=1.4, ls=":",
                label="instantaneous min / max"),
         Patch(facecolor=C_ENV, alpha=0.13, label="limit-cycle range")]
    if with_steady:
        h += [Line2D([], [], color=C_STEADY, lw=1.5, marker="o", ms=4.5,
                     label="steady base flow (SIMPLE)"),
              Line2D([], [], color=C_STEADY, lw=0, marker="o", ms=5.5,
                     mfc="white", mec=C_STEADY, mew=1.3,
                     label="SIMPLE limit-cycled: no fixed point")]
    h += [Line2D([], [], color="k", lw=0, marker="D", ms=5, mfc="white",
                 mew=1.2, label=r"$\alpha=11^\circ$: settled to a fixed point"),
          Line2D([], [], color=C_HOPF, lw=1.3, ls="--",
                 label=rf"Hopf, $\alpha_c={ALPHA_HOPF}^\circ$ (eigenvalue)"),
          Patch(facecolor=C_PD, alpha=0.20,
                label=r"period doubling, $(24^\circ,25^\circ]$")]
    return h


def main():
    plt.rcParams.update(RC)
    os.makedirs(FIGS, exist_ok=True)
    D = all_stats()

    specs = [
        ("cl", STEADY_CL, r"$C_L$",
         "Lift: the cycle never descends to the steady value", "cl_vs_aoa.png"),
        ("cd", STEADY_CD, r"$C_D$",
         "Drag: monotone in $\\alpha$, steady and unsteady alike",
         "cd_vs_aoa.png"),
        ("ld", None, r"$C_L/C_D$",
         "Efficiency: falls monotonically once shedding starts",
         "ld_vs_aoa.png"),
    ]

    # ---- combined three-panel figure -------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.9))
    for ax, (key, steady, ylab, title, _) in zip(axes, specs):
        panel(ax, D, key, steady, ylab, title)
        if key == "cl":
            add_aoa11(ax, AOA11["cl"])
        elif key == "cd":
            add_aoa11(ax, AOA11["cd"])
        else:
            add_aoa11(ax, AOA11["cl"] / AOA11["cd"])
            # steady L/D, from the two steady tables
            sa = [x for x in sorted(STEADY_CL) if x >= 10]
            ax.plot(sa, [STEADY_CL[x] / STEADY_CD[x] for x in sa], "-",
                    color=C_STEADY, lw=1.5, zorder=3)
            ax.plot([x for x in sa if x not in STEADY_NOT_CONVERGED],
                    [STEADY_CL[x] / STEADY_CD[x] for x in sa
                     if x not in STEADY_NOT_CONVERGED],
                    "o", color=C_STEADY, ms=4.5, zorder=4)
            ax.plot([x for x in sa if x in STEADY_NOT_CONVERGED],
                    [STEADY_CL[x] / STEADY_CD[x] for x in sa
                     if x in STEADY_NOT_CONVERGED],
                    "o", mfc="white", mec=C_STEADY, mew=1.3, ms=5.5, zorder=4)
        mark_bifurcations(ax, label=True)
    # No tight_layout here: with a legend anchored outside the axes it shrinks
    # the axes to a few percent of the canvas.  Reserve the space explicitly.
    fig.subplots_adjust(left=0.055, right=0.99, top=0.86, bottom=0.30,
                        wspace=0.24)
    fig.legend(handles=legend_handles(), loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=3, fontsize=9)
    fig.suptitle("Static NACA0012, Re = 500: force coefficients across two "
                 "bifurcations", y=0.97, fontsize=12)
    fig.savefig(os.path.join(FIGS, "coeffs_vs_aoa.png"))
    plt.close(fig)

    # ---- individual panels ----------------------------------------------
    for key, steady, ylab, title, name in specs:
        fig, ax = plt.subplots(figsize=(6.8, 5.6))
        panel(ax, D, key, steady, ylab, title)
        add_aoa11(ax, AOA11["cl"] if key == "cl" else
                  AOA11["cd"] if key == "cd" else AOA11["cl"] / AOA11["cd"])
        if key == "ld":
            sa = [x for x in sorted(STEADY_CL) if x >= 10]
            ax.plot(sa, [STEADY_CL[x] / STEADY_CD[x] for x in sa], "-o",
                    color=C_STEADY, lw=1.5, ms=4.5, zorder=3)
        mark_bifurcations(ax, label=True)
        fig.subplots_adjust(left=0.12, right=0.97, top=0.91, bottom=0.30)
        fig.legend(handles=legend_handles(with_steady=(key != "ld")),
                   loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=2,
                   fontsize=8.5)
        fig.savefig(os.path.join(FIGS, name))
        plt.close(fig)

    # ---- the numbers behind the figure -----------------------------------
    print(f"{'AoA':>5} {'ord':>3} {'Cl_mean':>9} {'Cl_min':>8} {'Cl_max':>8} "
          f"{'Cl_steady':>9} {'min-steady':>10} {'frac_p2p':>8} | "
          f"{'Cd_mean':>8} {'L/D':>7}")
    print(f"{11.0:5.1f} {'--':>3} {AOA11['cl']:9.4f} {'':>8} {'':>8} "
          f"{'':>9} {'':>10} {'':>8} | {AOA11['cd']:8.4f} "
          f"{AOA11['cl']/AOA11['cd']:7.3f}   (fixed point)")
    for i in ANGLES:
        d = D[i]
        st = STEADY_CL.get(i)
        gap = (d["cl_min"] - st) if st is not None else float("nan")
        print(f"{i:5d} {d['order']:3d} {d['cl_mean']:9.4f} {d['cl_min']:8.4f} "
              f"{d['cl_max']:8.4f} "
              + (f"{st:9.4f} {gap:10.4f} {gap/d['cl_p2p']:8.2f}"
                 if st is not None else f"{'--':>9} {'--':>10} {'--':>8}")
              + f" | {d['cd_mean']:8.4f} {d['ld_of_means']:7.3f}")
    print("\nfigures:", ", ".join(["coeffs_vs_aoa.png"] + [s[-1] for s in specs]))


if __name__ == "__main__":
    main()
