#!/usr/bin/env python3
"""Decide whether the C_L transient has decayed, and measure the shedding period.

This is the gate between the two run stages.  Stage 1 runs field-free until this
says CONVERGED; stage 2 then writes 32 snapshots per cycle at the period this
reports.  Both answers come from the same per-step C_L series, so the sampling
rate is measured rather than guessed from a Strouhal correlation.

Convergence test: split the trailing part of the series into two consecutive
windows and require the mean, the peak-to-peak amplitude and the dominant
frequency to agree between them to within `--tol`.  Mean alone is not enough --
it settles long before the amplitude does -- and amplitude alone misses a slow
mean drift, which is exactly the "converges slowly toward the limit cycle"
behaviour noted in experiment_plan.md.

Output is KEY=VALUE lines so the shell driver can eval it.
"""

import argparse
import glob
import os
import sys

import numpy as np


def load_series(case):
    """Concatenate every forceCoeffs restart directory into one time series.

    On restart forceCoeffs opens postProcessing/forceCoeffs/<restartTime>/, so a
    staged run leaves the series split across directories with overlapping
    times at the seams.  Sort by time and keep the LAST sample at each time:
    the later run is the one that continued.
    """
    pat = os.path.join(case, "postProcessing", "forceCoeffs", "*", "*.dat")
    files = sorted(glob.glob(pat))
    if not files:
        sys.exit(f"no forceCoeffs output under {case}/postProcessing/forceCoeffs")

    t, cl, cd = [], [], []
    for f in files:
        names = None
        for line in open(f):
            if line.startswith("#"):
                parts = line[1:].split()
                if "Cl" in parts or "Cd" in parts:
                    names = parts
                continue
            if not line.strip():
                continue
            v = line.split()
            if names is None:            # fall back to the documented layout
                # Verified against OpenFOAM-dev build dev-84333c37c7aa output.
                names = ["Time", "Cm", "Cd", "Cl", "Cl(f)", "Cl(r)"]
            try:
                row = [float(x) for x in v]
            except ValueError:
                continue
            idx = {n: i for i, n in enumerate(names)}
            if "Cl" not in idx or "Cd" not in idx or len(row) <= max(idx["Cl"], idx["Cd"]):
                continue
            t.append(row[0]); cl.append(row[idx["Cl"]]); cd.append(row[idx["Cd"]])

    t = np.asarray(t); cl = np.asarray(cl); cd = np.asarray(cd)
    order = np.argsort(t, kind="stable")
    t, cl, cd = t[order], cl[order], cd[order]
    keep = np.ones(len(t), bool)
    keep[:-1] = t[:-1] != t[1:]          # duplicate time -> keep the last
    return t[keep], cl[keep], cd[keep]


def dominant_freq(t, y):
    """Dominant non-zero frequency of y(t), assuming a uniform dt.

    Parabolic interpolation on the three bins around the peak: with only ~30
    cycles in the window the bin width is a few percent of f, which is the same
    order as the convergence tolerance we are testing against.
    """
    if len(t) < 16:
        return 0.0, 0.0
    dt = np.median(np.diff(t))
    y = y - y.mean()
    w = np.hanning(len(y))
    Y = np.abs(np.fft.rfft(y * w))
    f = np.fft.rfftfreq(len(y), dt)
    Y[0] = 0.0
    k = int(np.argmax(Y))
    if k == 0:
        return 0.0, 0.0
    if 0 < k < len(Y) - 1:
        a, b, c = Y[k - 1], Y[k], Y[k + 1]
        denom = a - 2 * b + c
        delta = 0.5 * (a - c) / denom if denom != 0 else 0.0
    else:
        delta = 0.0
    df = f[1] - f[0]
    # Prominence = peak / median spectrum.  A saturated limit cycle has a sharp
    # narrow peak (prominence >> 1); the deceptive quiet plateau described in
    # ../FINDINGS.md Part 5 has a tiny broad-band trace whose "peak" is barely
    # above the noise floor.  Amplitude alone cannot tell those apart.
    floor = float(np.median(Y[1:])) if len(Y) > 2 else 0.0
    prom = float(Y[k] / floor) if floor > 0 else 0.0
    return float(f[k] + delta * df), float(Y[k]), prom


def window_stats(t, y):
    f, amp, prom = dominant_freq(t, y)
    return dict(mean=float(y.mean()),
                p2p=float(np.percentile(y, 99) - np.percentile(y, 1)),
                freq=f, peak=amp, prom=prom)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--tol", type=float, default=0.01,
                    help="relative agreement required between the two windows")
    ap.add_argument("--min-cycles", type=float, default=6.0,
                    help="cycles required in EACH window before judging")
    ap.add_argument("--min-prominence", type=float, default=8.0,
                    help="spectral peak / noise floor required to call it shedding")
    args = ap.parse_args()

    t_all, cl_all, cd_all = load_series(args.case)

    # Anchor the test windows on CYCLES counted back from the end, not on a
    # fraction of the series.  A fraction wastes the whole early transient: once
    # the limit cycle exists, only the last 2*min_cycles periods carry
    # information, and waiting for them to become half the run means integrating
    # roughly twice as long as necessary before a verdict is even possible.
    half = len(t_all) // 2
    f0 = dominant_freq(t_all[half:], cl_all[half:])[0]
    if f0 > 0:
        # 1.15x margin: f0 is a rough first-pass estimate and the half-split is
        # not exact, so sizing the window at exactly 2*min_cycles/f0 leaves
        # window B holding a hair UNDER min_cycles (measured: 2.99 for a
        # requested 3).  The >= test below then fails permanently and the gate
        # never fires, running every case to T1_MAX.
        need = 2.0 * args.min_cycles / f0 * 1.15
        sel = t_all >= max(t_all[0], t_all[-1] - need)
    else:
        sel = t_all >= t_all[0] + 0.5 * (t_all[-1] - t_all[0])
    t, cl, cd = t_all[sel], cl_all[sel], cd_all[sel]

    half = len(t) // 2
    wa, wb = (t[:half], cl[:half]), (t[half:], cl[half:])
    sa, sb = window_stats(*wa), window_stats(*wb)

    f = sb["freq"]
    period = 1.0 / f if f > 0 else 0.0
    span = t[-1] - t[half]
    cycles = span * f

    # A flat series has no limit cycle to converge to. Guard before dividing.
    #
    # The prominence term is the guard against ../FINDINGS.md Part 5: at an angle
    # just past onset the impulsive-start disturbance DECAYS to a quiet, almost
    # constant plateau (alpha=10 sat at Cl~0.353, amplitude ~7e-4, over t=8-12)
    # and only then does the real instability re-emerge and grow to saturation
    # near t~22.  Mean, amplitude and frequency are all stable across that
    # plateau, so the three agreement tests below would happily call it
    # converged.  A run stopped there reports "steady, no shedding" and is wrong
    # by 15%.  The plateau's spectrum is broad-band noise, so requiring a sharp
    # peak rejects it; the driver's T1_MIN is the second line of defence.
    # Amplitude threshold matches analysis/shedding_analysis.py upstream, which
    # classified the documented Re=1000 sweep with amp > 5e-3: alpha=5 and 8 sit
    # at amp 8e-4 and 1.5e-3 (STEADY), alpha=10 at 2.5e-2 (SHEDDING).  Without
    # it, slow drift in a steady run produces a spurious low-frequency "peak"
    # with a high prominence and the driver would proceed on a garbage period.
    shedding = (sb["p2p"] > 2 * 5e-3 and f > 0
                and sb["prom"] >= args.min_prominence)

    def rel(a, b, scale):
        return abs(a - b) / scale if scale > 0 else 0.0

    d_mean = rel(sa["mean"], sb["mean"], max(abs(sb["mean"]), sb["p2p"], 1e-9))
    d_p2p = rel(sa["p2p"], sb["p2p"], max(sb["p2p"], 1e-9))
    d_frq = rel(sa["freq"], sb["freq"], max(sb["freq"], 1e-9))

    ok = (shedding
          and cycles >= args.min_cycles
          and max(d_mean, d_p2p, d_frq) <= args.tol)

    print(f"T_END={t[-1]:.6f}")
    print(f"SHEDDING={'yes' if shedding else 'no'}")
    print(f"FREQ={f:.8f}")
    print(f"PERIOD={period:.8f}")
    print(f"CYCLES_IN_WINDOW={cycles:.2f}")
    print(f"CL_MEAN={sb['mean']:.6f}")
    print(f"CL_P2P={sb['p2p']:.6f}")
    print(f"PROMINENCE={sb['prom']:.2f}")
    print(f"CD_MEAN={cd[half:].mean():.6f}")
    print(f"D_MEAN={d_mean:.5f}")
    print(f"D_P2P={d_p2p:.5f}")
    print(f"D_FREQ={d_frq:.5f}")
    print(f"CONVERGED={'yes' if ok else 'no'}")


if __name__ == "__main__":
    main()
