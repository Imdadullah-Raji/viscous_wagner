#!/usr/bin/env python3
"""One place that says what the campaign measured, so the figures agree.

Every plot that quotes a coefficient reads it from here, computed from the raw
force series rather than transcribed into a literal.  Two rules the numbers
depend on:

  * statistics are taken over the CONVERGED window only, cut to a whole number
    of FUNDAMENTAL periods.  A window holding a fraction of an orbit reports a
    biased mean and a min/max that depend on where the series happened to stop
    -- and past 24 deg the orbit is two shedding cycles, not one.
  * `spectral_structure` decides which frequency is the carrier.  Never the
    tallest peak: at AoA 28 the tallest peak is the orbit fundamental.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from limit_cycle import load_series, spectral_structure   # noqa: E402

# Bifurcation points, both measured.  See FINDINGS.md sections 19 and 11.
ALPHA_HOPF = 11.46          # +-0.1, from eigenvalues about the steady base flow
ALPHA_PD_LO, ALPHA_PD_HI = 24.0, 25.0   # period doubling, bracketed not resolved

# t_conv for the two runs that finished before the .phase2 marker existed.
T_CONV_FALLBACK = {20: 60.5, 30: 120.5}

ANGLES = [12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

# Steady (SIMPLE) base flow, from analysis/steady_summary.py.  22 and 24 did not
# converge -- SIMPLE limit-cycles there, so no fixed point was found.
STEADY_CL = {0: 0.000000, 2: 0.112979, 4: 0.218670, 6: 0.309705, 8: 0.379699,
             10: 0.426072, 12: 0.452294, 14: 0.465639, 16: 0.472710,
             18: 0.477540, 20: 0.482076, 22: 0.494523, 24: 0.505358}
STEADY_CD = {0: 0.174132, 2: 0.175460, 4: 0.179485, 6: 0.186303, 8: 0.195962,
             10: 0.208307, 12: 0.222876, 14: 0.239060, 16: 0.256362,
             18: 0.274533, 20: 0.293524, 22: 0.316797, 24: 0.340897}
STEADY_NOT_CONVERGED = {22, 24}

# The collaborator's AoA 11: a time-accurate run that settled onto a FIXED
# POINT, C_L flat to 2e-5 over the last 10 t.  It is the only unsteady-solver
# point below the Hopf and the one that anchors the left end of the curve.
AOA11 = dict(aoa=11.0, cl=0.441474, cd=0.215378,
             source="collaborator, /home/raji/OpenFOAM/raji-dev/run/Sept/re500_aoa11")


def t_conv_of(case_dir, aoa):
    p = os.path.join(case_dir, ".phase2")
    if os.path.exists(p):
        for line in open(p):
            if line.startswith("T_CONV="):
                return float(line.split("=", 1)[1])
    return T_CONV_FALLBACK[aoa]


def uniform(t, y, dt=None):
    """Resample onto a uniform grid -- dt changes between chunks, so the raw
    series is only piecewise uniform and cannot be handed to an FFT."""
    dt = dt or float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, y)


def amplitude_spectrum(t, y, n_orbits=None, f_fund=None):
    """Single-sided AMPLITUDE spectrum of y(t), in the units of y.

    Hann-windowed and corrected by the window's coherent gain, so the height of
    an isolated peak is the physical amplitude of that sinusoid (half its
    peak-to-peak), not an FFT bin value that depends on the record length.

    The record is trimmed to a whole number of fundamental periods when one is
    given: an integer number of periods puts the carrier and every harmonic on
    a bin centre, which is what keeps the subharmonic from being buried in the
    carrier's leakage skirt.
    """
    tu, yu = uniform(t, y)
    dt = tu[1] - tu[0]
    if f_fund:
        n_per = int((tu[-1] - tu[0]) * f_fund)
        if n_orbits:
            n_per = min(n_per, int(n_orbits))
        if n_per >= 2:
            n = int(round(n_per / f_fund / dt))
            tu, yu = tu[-n:], yu[-n:]
    yu = yu - yu.mean()
    w = np.hanning(len(yu))
    cg = w.sum() / len(yu)
    Y = np.abs(np.fft.rfft(yu * w)) * 2.0 / (len(yu) * cg)
    f = np.fft.rfftfreq(len(yu), dt)
    return f, Y, tu[-1] - tu[0]


def stats(aoa, root=ROOT):
    """Converged-window statistics for one angle."""
    case = os.path.join(root, "NACA0012", f"re500_aoa{aoa}")
    t, cl, cd = load_series(case)
    tc = t_conv_of(case, aoa)
    sel = t >= tc
    t, cl, cd = t[sel], cl[sel], cd[sel]

    f_fund, f_car, order, sub_ratio = spectral_structure(t, cl)
    # cut to a whole number of orbits, taken from the END of the record
    if f_fund > 0:
        n_orb = int((t[-1] - t[0]) * f_fund)
        if n_orb >= 1:
            m = t >= t[-1] - n_orb / f_fund
            t, cl, cd = t[m], cl[m], cd[m]
    else:
        n_orb = 0

    ld = cl / cd
    return dict(aoa=aoa, t_conv=tc, t_end=t[-1], n_orbits=n_orb,
                f_fund=f_fund, f_carrier=f_car, order=order,
                sub_ratio=sub_ratio,
                cl_mean=cl.mean(), cl_min=cl.min(), cl_max=cl.max(),
                cl_p2p=cl.max() - cl.min(),
                cd_mean=cd.mean(), cd_min=cd.min(), cd_max=cd.max(),
                ld_mean=ld.mean(), ld_min=ld.min(), ld_max=ld.max(),
                ld_of_means=cl.mean() / cd.mean(),
                t=t, cl=cl, cd=cd)


def all_stats(angles=ANGLES, root=ROOT):
    return {a: stats(a, root) for a in angles}


# Okabe-Ito, colour-blind safe and legible in greyscale.  Used everywhere so a
# given angle keeps its colour from figure to figure.
OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
             "#56B4E9", "#F0E442", "#000000"]
ANGLE_COLOUR = {12: "#0072B2", 14: "#56B4E9", 16: "#009E73", 18: "#66A61E",
                20: "#E69F00", 22: "#D55E00", 24: "#A6761D", 26: "#CC79A7",
                28: "#7570B3", 30: "#000000"}

RC = {
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


if __name__ == "__main__":
    hdr = (f"{'AoA':>4} {'ord':>3} {'f_car':>8} {'orbits':>6} "
           f"{'Cl_mean':>9} {'Cl_min':>9} {'Cl_max':>9} "
           f"{'Cd_mean':>9} {'Cd_min':>9} {'Cd_max':>9} {'L/D':>7}")
    print(hdr)
    for a in ANGLES:
        s = stats(a)
        print(f"{a:4d} {s['order']:3d} {s['f_carrier']:8.4f} {s['n_orbits']:6d} "
              f"{s['cl_mean']:9.4f} {s['cl_min']:9.4f} {s['cl_max']:9.4f} "
              f"{s['cd_mean']:9.4f} {s['cd_min']:9.4f} {s['cd_max']:9.4f} "
              f"{s['ld_of_means']:7.3f}")
