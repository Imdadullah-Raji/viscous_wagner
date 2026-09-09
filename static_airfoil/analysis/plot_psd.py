#!/usr/bin/env python3
"""Amplitude spectra of the force coefficients, one figure per angle.

  analysis/figs/PSD/psd_aoa<NN>.png

Amplitude, not power: the y axis is in coefficient units, so the height of an
isolated peak is the physical amplitude of that component (half its
peak-to-peak) and can be read straight off against C_L itself.  A power
spectrum squares that and hides the factor-of-two structure this campaign is
about.

Three things the figures are built to show:
  * the carrier and its harmonic ladder, marked from `spectral_structure` --
    never from the tallest peak, which is the orbit fundamental at AoA 28;
  * the SUBHARMONIC at f_c/2, absent at 12-24 and present at 26-30.  The
    period-1 angles are plotted with the same f_c/2 guide line so the reader
    sees it land on the noise floor rather than having to take it on trust;
  * the noise floor itself, so "no peak here" is a bounded statement.

Record length is cut to a whole number of orbits (see campaign_data) -- with an
integer number of periods every harmonic sits on a bin centre, and the
subharmonic is not buried under the carrier's leakage.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign_data import (RC, ROOT, amplitude_spectrum, stats,   # noqa: E402
                           ANGLE_COLOUR)

ANGLES = [12, 16, 20, 24, 28]
OUT = os.path.join(ROOT, "analysis", "figs", "PSD")

CL_COLOUR = "#0072B2"
CD_COLOUR = "#D55E00"


def peak_amp(f, Y, target, halfwidth=2):
    k = int(np.argmin(np.abs(f - target)))
    lo = max(1, k - halfwidth)
    return float(Y[lo:k + halfwidth + 1].max())


def is_local_peak(f, Y, target, halfwidth=2, wide=8):
    """Is there a genuine local maximum at `target`, or only the skirt of a
    neighbouring peak?

    An amplitude threshold cannot tell these apart.  At AoA 20 the value at
    f_c/2 is 45x the noise floor purely because it sits on the carrier's
    leakage skirt, and calling that a residual subharmonic would invent a
    physical component out of a window function.  A peak has to be higher than
    both of its shoulders.
    """
    k = int(np.argmin(np.abs(f - target)))
    lo, hi = max(1, k - halfwidth), min(len(Y) - 1, k + halfwidth + 1)
    kk = lo + int(np.argmax(Y[lo:hi]))
    l0 = max(1, kk - wide)
    r1 = min(len(Y), kk + wide + 1)
    return bool(Y[kk] > Y[l0:kk].min() * 1.5 and Y[kk] > Y[kk + 1:r1].min() * 1.5
                and Y[kk] >= Y[max(1, kk - 1)] and Y[kk] >= Y[min(len(Y) - 1, kk + 1)])


def one(aoa):
    s = stats(aoa)
    f_c, f_f, order = s["f_carrier"], s["f_fund"], s["order"]
    fl, Yl, span = amplitude_spectrum(s["t"], s["cl"], f_fund=f_f)
    fd, Yd, _ = amplitude_spectrum(s["t"], s["cd"], f_fund=f_f)

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    fmax = 6.4 * f_c

    # noise floor: median of the upper half of the band, where nothing physical
    # remains -- quoted so an absent peak is an upper bound, not a claim.
    hi = fl > 0.6 * fmax
    floor = float(np.median(Yl[hi & (fl <= fmax)]))

    # harmonic ladder of the carrier, and the subharmonic half-ladder
    for n in range(1, 7):
        ax.axvline(n * f_c, color="0.82", lw=0.7, zorder=0)
    ax.axvline(f_c / 2, color="#CC79A7", lw=1.0, ls=(0, (4, 3)), zorder=1)

    m = (fl > 0) & (fl <= fmax)
    ax.semilogy(fl[m], Yl[m], lw=1.1, color=CL_COLOUR, label=r"$C_L$")
    md = (fd > 0) & (fd <= fmax)
    ax.semilogy(fd[md], Yd[md], lw=0.9, color=CD_COLOUR, alpha=0.75,
                label=r"$C_D$")

    a_car = peak_amp(fl, Yl, f_c)
    a_sub = peak_amp(fl, Yl, f_c / 2)
    ax.plot([f_c], [a_car], "v", ms=6, color=CL_COLOUR, mfc="white", mew=1.2,
            zorder=5)
    ax.annotate(rf"$f_c={f_c:.4f}$" "\n" rf"$A={a_car:.4f}$",
                xy=(f_c, a_car), xytext=(9, 4), textcoords="offset points",
                fontsize=9, color=CL_COLOUR, va="bottom", ha="left")

    if order >= 2:
        ax.plot([f_c / 2], [a_sub], "o", ms=5, color="#CC79A7", mfc="white",
                mew=1.2, zorder=5)
        # to the RIGHT of the marker: f_c/2 sits near the left spine, and a
        # right-aligned label there runs off the canvas.
        ax.annotate(rf"$f_c/2 = {f_c/2:.4f}$" "\n" rf"$A/A_c={a_sub/a_car:.2f}$",
                    xy=(f_c / 2, a_sub), xytext=(-9, 12),
                    textcoords="offset points", fontsize=9, color="#CC79A7",
                    ha="right", va="bottom")
        verdict = (rf"period-2:  orbit closes after two shedding cycles"
                   "\n" rf"$T_{{orbit}}=2/f_c={1/f_f:.3f}$")
    else:
        # At 24 deg there IS a small f_c/2 component (2.6e-4 of the carrier,
        # decaying at sigma = -0.137/t); at 12-20 there is nothing but the
        # carrier's leakage skirt.  Distinguish them by shape, not by height.
        note = ("residual, decaying" if is_local_peak(fl, Yl, f_c / 2)
                else "no peak: carrier leakage only")
        ax.annotate(f"$f_c/2$:  $A/A_c$ = {a_sub / a_car:.1e}\n{note}",
                    xy=(f_c / 2, a_sub), xytext=(9, 14),
                    textcoords="offset points", fontsize=9, color="#CC79A7",
                    ha="left", va="bottom")
        verdict = (rf"period-1:  orbit closes every shedding cycle"
                   "\n" rf"$T_{{orbit}}=1/f_c={1/f_f:.3f}$")

    ax.axhline(floor, color="0.55", lw=0.8, ls=":")
    ax.text(fmax, floor, "  noise floor ", ha="right", va="bottom",
            fontsize=8, color="0.4")

    ax.set_xlim(0, fmax)
    ax.set_ylim(max(floor * 0.25, Yl[m].min() * 0.5),
                max(a_car, a_sub) * 20)
    ax.set_xlabel(r"frequency  $f$   ($U_\infty/c$)")
    ax.set_ylabel("amplitude   (coefficient units)")
    ax.set_title(rf"$\alpha={aoa}^\circ$,  Re$=500$" "\n"
                 rf"{s['n_orbits']} orbits, $t={s['t_conv']:.1f}$–"
                 rf"${s['t_end']:.1f}$", loc="left")
    ax.text(0.985, 0.97, verdict, transform=ax.transAxes, ha="right",
            va="top", fontsize=9, color="0.25", linespacing=1.4)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62))
    ax.grid(alpha=0.25, which="major", axis="y", lw=0.5)

    path = os.path.join(OUT, f"psd_aoa{aoa:02d}.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"alpha={aoa:2d}  order={order}  f_c={f_c:.4f}  A(f_c)={a_car:.5f}  "
          f"A(f_c/2)={a_sub:.3e}  ratio={a_sub/a_car:.4f}  "
          f"floor={floor:.2e}  span={span:.1f}t  -> {os.path.basename(path)}")


def main():
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams.update(RC)
    for a in ANGLES:
        one(a)
    print("figures in", OUT)


if __name__ == "__main__":
    main()
