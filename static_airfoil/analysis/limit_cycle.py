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

PERIOD DOUBLING.  Past some angle between 20 and 30 deg at Re=500 the wake is no
longer period-1: successive shedding cycles alternate, so the orbit closes only
after TWO oscillations of C_L.  Measured at AoA 30: carrier f = 0.3199, a
saturated subharmonic at f/2 = 0.1600 carrying 63% of the carrier's amplitude,
and every resolved spectral peak an exact integer multiple of f/2 (max deviation
0.0000 -- a closed period-2 orbit, not an incommensurate torus).

That breaks two things if the largest spectral peak is taken as "the" period:

  * cycle counting.  20 "cycles" of the carrier is only 10 cycles of the actual
    orbit, so a phase-2 run records half the requested data.
  * the convergence test.  Alternate cycles have genuinely different means and
    amplitudes, so a window holding a non-integer number of FUNDAMENTAL periods
    -- or split into two halves that do -- compares different phases of the
    orbit against each other and reports a difference that is not drift.  AoA 30
    reported CONVERGED=no at the end of a run whose subharmonic amplitude was
    flat to 6 parts in 100 000.

So this script separates two periods and reports both:

  PERIOD         the fundamental -- the orbit's true period, m x the carrier
                 period for a period-m state.  Sets how LONG phase 2 runs and
                 how the convergence window is cut.
  SAMPLE_PERIOD  the carrier -- the fastest oscillation actually in the signal.
                 Sets the field WRITE INTERVAL, so PER_CYCLE samples per carrier
                 period are kept regardless of the subharmonic order.

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


def spectrum(t, y):
    """Hann-windowed amplitude spectrum of y(t), on a uniform time grid.

    The series is only PIECEWISE uniform: run_staged.sh re-tunes dt between
    chunks, so a convergence window that straddles a chunk boundary carries two
    or three different dt.  Handing that straight to rfft with a median dt
    stretches part of the record and smears the peak -- which matters most
    exactly when it is being used to detect a subharmonic.  Resample first.
    """
    dt = float(np.median(np.diff(t)))
    if dt <= 0 or len(t) < 16:
        return None, None, 0.0
    tu = np.arange(t[0], t[-1], dt)
    yu = np.interp(tu, t, y)
    yu = yu - yu.mean()
    Y = np.abs(np.fft.rfft(yu * np.hanning(len(yu))))
    f = np.fft.rfftfreq(len(yu), dt)
    Y[0] = 0.0
    return f, Y, dt


def _refine(f, Y, k):
    """Parabolic interpolation on the three bins around bin k."""
    if not 0 < k < len(Y) - 1:
        return float(f[k])
    a, b, c = Y[k - 1], Y[k], Y[k + 1]
    denom = a - 2 * b + c
    delta = 0.5 * (a - c) / denom if denom != 0 else 0.0
    return float(f[k] + delta * (f[1] - f[0]))


def dominant_freq(t, y):
    """Dominant non-zero frequency of y(t).

    Parabolic interpolation on the three bins around the peak: with only ~30
    cycles in the window the bin width is a few percent of f, which is the same
    order as the convergence tolerance we are testing against.

    This is the CARRIER -- the largest peak.  For a period-doubled state the
    orbit's own frequency is a submultiple of it; see subharmonic_order().
    """
    if len(t) < 16:
        return 0.0, 0.0, 0.0
    f, Y, dt = spectrum(t, y)
    if f is None or len(Y) < 3:
        return 0.0, 0.0, 0.0
    k = int(np.argmax(Y))
    if k == 0:
        return 0.0, 0.0, 0.0
    # Prominence = peak / median spectrum.  A saturated limit cycle has a sharp
    # narrow peak (prominence >> 1); the deceptive quiet plateau described in
    # ../FINDINGS.md Part 5 has a tiny broad-band trace whose "peak" is barely
    # above the noise floor.  Amplitude alone cannot tell those apart.
    floor = float(np.median(Y[1:])) if len(Y) > 2 else 0.0
    prom = float(Y[k] / floor) if floor > 0 else 0.0
    return _refine(f, Y, k), float(Y[k]), prom


def subharmonic_order(t, y, max_order=4, min_ratio=0.05, min_prominence=8.0):
    """Is the orbit period-1, or does it close only after m oscillations?

    Look for a sharp peak at f_carrier/m.  Two conditions, both needed:

      * amplitude at least `min_ratio` of the carrier's.  A period-doubled wake
        is not subtle -- AoA 30 sits at 0.63 -- while a period-1 wake reads
        1e-4, so the threshold is nowhere near either.  It exists to reject
        low-frequency drift, not to make a close call.
      * a genuine local maximum standing `min_prominence` above the median
        spectrum.  Broad-band content in the quiet plateau (../FINDINGS.md
        Part 5) has power everywhere including at f/2; it has no peak there.

    Highest order first, so a period-4 state is not mistaken for period-2.
    Returns (m, amplitude ratio, carrier amplitude).
    """
    f, Y, dt = spectrum(t, y)
    if f is None or len(Y) < 8:
        return 1, 0.0, 0.0
    k_dom = int(np.argmax(Y))
    if k_dom < 2:
        return 1, 0.0, 0.0
    floor = float(np.median(Y[1:]))
    for m in range(int(max_order), 1, -1):
        k = int(round(k_dom / m))
        if k < 2 or k + 2 >= len(Y):
            continue
        lo = k - 1
        kk = lo + int(np.argmax(Y[lo:k + 2]))
        if not (0 < kk < len(Y) - 1 and Y[kk] >= Y[kk - 1] and Y[kk] >= Y[kk + 1]):
            continue                      # not a peak, just spectral skirt
        ratio = float(Y[kk] / Y[k_dom])
        prom = float(Y[kk] / floor) if floor > 0 else 0.0
        if ratio >= min_ratio and prom >= min_prominence:
            return m, ratio, float(Y[k_dom])
    return 1, 0.0, float(Y[k_dom])


def spectral_structure(t, y, max_order=4, min_ratio=0.05, min_prominence=8.0,
                       carrier_frac=0.5):
    """Separate the ORBIT frequency from the SHEDDING frequency.

    The earlier version took the largest spectral peak as the carrier and looked
    BELOW it for a subharmonic.  That assumes the carrier is always the tallest
    peak, and it is not:

        AoA 30   f=0.3202 A=12806 (carrier)   f/2=0.1593 A=8580   sub/car = 0.67
        AoA 28   f=0.1893 A=3006  (FUNDAMENTAL) 2f=0.3781 A=2018  car/sub = 0.67

    Same period-2 structure, opposite dominance.  At AoA 28 the search found
    nothing below the peak and reported period-1 at f=0.189 -- which happens to
    give the right PERIOD but halves SAMPLE_PERIOD, so phase 2 would have written
    32 samples per ORBIT (16 per shedding cycle) instead of 32 per shedding
    cycle, and D_SUB would not have been computed at all.

    A(2f)/A(f) cannot be used to tell the two apart: it reads 0.270 at AoA 30
    (period-2) and 0.195 at AoA 24 (period-1).  So find the two frequencies
    independently instead:

      f_carrier  the HIGHEST peak still carrying `carrier_frac` of the maximum
                 amplitude.  "Resolve everything with at least half the peak's
                 amplitude" -- this is what the sampling rate has to satisfy.
      f_fund     the LOWEST peak above `min_ratio`, which is the orbit frequency.
                 The looser threshold matters near onset, where the subharmonic
                 emerges small and would be missed by carrier_frac.

    Returns (f_fund, f_carrier, order, sub_ratio).
    """
    f, Y, dt = spectrum(t, y)
    if f is None or len(Y) < 8:
        return 0.0, 0.0, 1, 0.0
    floor = float(np.median(Y[1:]))
    amax = float(Y.max())
    k_dom = int(np.argmax(Y))
    peaks = []
    for k in range(2, len(Y) - 1):
        if Y[k] >= Y[k - 1] and Y[k] >= Y[k + 1] and Y[k] >= min_ratio * amax \
           and (floor <= 0 or Y[k] / floor >= min_prominence):
            peaks.append((k, _refine(f, Y, k), float(Y[k])))
    if not peaks:
        return _refine(f, Y, k_dom), _refine(f, Y, k_dom), 1, 0.0

    strong = [p for p in peaks if p[2] >= carrier_frac * amax]
    car = max(strong, key=lambda p: p[1]) if strong else (k_dom, _refine(f, Y, k_dom), amax)
    f_car, a_car = car[1], car[2]
    k_fund, f_low, a_low = min(peaks, key=lambda p: p[1])

    order, sub_ratio = 1, 0.0
    if f_low > 0 and f_car > f_low * 1.5:
        r = f_car / f_low
        m = int(round(r))
        # Reject anything that is not a clean integer ratio: that is a second,
        # incommensurate frequency (a torus), not a subharmonic of this one.
        if 2 <= m <= max_order and abs(r - m) <= 0.15:
            order = m
            # Relative to the CARRIER, not to the largest peak: at AoA 28 the
            # subharmonic IS the largest peak, and a_low/amax would read 1.0 for
            # every such case regardless of how strong the subharmonic actually
            # is.  A(f_fund)/A(f_carrier) is the quantity that means something
            # and is comparable across angles -- it is the bifurcation diagram.
            sub_ratio = a_low / a_car if a_car > 0 else 0.0
    # Derive the fundamental FROM the carrier: the carrier is a strong, sharp
    # peak, so its refined location is the more precise of the two.
    f_fund = f_car / order if f_car > 0 else 0.0
    return f_fund, f_car, order, float(sub_ratio)

def harmonic_deviation(t, y, f_fund, min_ratio=0.02):
    """Are all significant peaks integer multiples of f_fund?

    The test that separates a period-m limit cycle from a torus.  A
    period-doubling bifurcation produces one closed orbit, so every peak sits at
    an exact multiple of the orbit frequency; a Neimark-Sacker bifurcation
    produces a second, incommensurate frequency and the ratios are not integers.
    Returns the largest distance from an integer over all peaks above
    `min_ratio` of the tallest, in units of f_fund.  Measured 0.0000 at AoA 30.
    """
    f, Y, dt = spectrum(t, y)
    if f is None or f_fund <= 0 or len(Y) < 8:
        return 0.0
    thr = min_ratio * Y.max()
    dev = 0.0
    for k in range(2, len(Y) - 1):
        if Y[k] >= thr and Y[k] >= Y[k - 1] and Y[k] >= Y[k + 1]:
            r = _refine(f, Y, k) / f_fund
            dev = max(dev, abs(r - round(r)))
    return float(dev)


def window_stats(t, y, max_order=4, min_ratio=0.05, min_prominence=8.0):
    """Per-window statistics.

    `freq` is the structure-aware CARRIER, never the raw largest peak.  Which
    peak is tallest flips as SUBHARM_RATIO crosses 1, and AoA 26 sits almost
    exactly there: window A read ratio 1.072 (fundamental taller, largest peak
    0.19922) and window B read 0.948 (carrier taller, largest peak 0.39532), so
    comparing "dominant frequencies" between the windows compared two different
    physical quantities and pinned D_FREQ at 0.496 -- a factor of two -- for
    every chunk, while the carriers actually agreed to 0.007.  Left alone the run
    could never converge and would have burned to T1_MAX.
    """
    f_fund, f_car, order, ratio = spectral_structure(
        t, y, max_order=max_order, min_ratio=min_ratio,
        min_prominence=min_prominence)
    _, amp, prom = dominant_freq(t, y)
    return dict(mean=float(y.mean()),
                p2p=float(np.percentile(y, 99) - np.percentile(y, 1)),
                freq=f_car, fund=f_fund, order=order, ratio=ratio,
                peak=amp, prom=prom)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--tol", type=float, default=0.01,
                    help="relative agreement required between the two windows")
    ap.add_argument("--min-cycles", type=float, default=6.0,
                    help="FUNDAMENTAL cycles required in EACH window before judging")
    ap.add_argument("--min-prominence", type=float, default=8.0,
                    help="spectral peak / noise floor required to call it shedding")
    ap.add_argument("--max-order", type=int, default=4,
                    help="highest subharmonic order to test for (1 disables)")
    ap.add_argument("--subharmonic-ratio", type=float, default=0.05,
                    help="A(f/m)/A(f) needed to call the state period-m")
    args = ap.parse_args()

    t_all, cl_all, cd_all = load_series(args.case)

    # First pass on the trailing half: the carrier, and whether the orbit is
    # period-1 or period-m.  Both are needed before the test window can even be
    # sized, because the window is cut in FUNDAMENTAL periods.
    half = len(t_all) // 2
    f_fund0, f_car0, order, sub_ratio = spectral_structure(
        t_all[half:], cl_all[half:], max_order=max(1, args.max_order),
        min_ratio=args.subharmonic_ratio, min_prominence=args.min_prominence)

    # Anchor the test windows on CYCLES counted back from the end, not on a
    # fraction of the series.  A fraction wastes the whole early transient: once
    # the limit cycle exists, only the last 2*min_cycles periods carry
    # information, and waiting for them to become half the run means integrating
    # roughly twice as long as necessary before a verdict is even possible.
    #
    # Cut a WHOLE number of fundamental periods, and split at the exact midpoint
    # in time so each half holds the same whole number.  For a period-2 state
    # alternate cycles differ by construction -- at AoA 30 the per-cycle mean
    # alternates 0.9639 / 1.0623 and the per-cycle peak-to-peak 0.41 / 0.83 --
    # so a half-window ending mid-orbit reports that alternation as a difference
    # between the windows.  That is a windowing artefact, not drift, and it is
    # what made a fully saturated AoA 30 read CONVERGED=no.  (The old code sized
    # the window as 2*min_cycles/f * 1.15; the 1.15 was there because an exact
    # cut left window B a hair UNDER min_cycles and the gate then never fired.
    # Cutting whole periods and counting cycles from the midpoint removes the
    # need for the fudge factor.)
    n_per = max(1, int(np.ceil(args.min_cycles)))
    if f_fund0 > 0:
        need = 2.0 * n_per / f_fund0
        t_cut = t_all[-1] - need
        full_window = t_cut >= t_all[0]          # did we get the window we asked for?
        t_cut = max(t_all[0], t_cut)
    else:
        t_cut = t_all[0] + 0.5 * (t_all[-1] - t_all[0])
        full_window = False
    sel = t_all >= t_cut
    t, cl, cd = t_all[sel], cl_all[sel], cd_all[sel]

    t_mid = 0.5 * (t[0] + t[-1])
    a = t < t_mid
    b = ~a
    if a.sum() < 16 or b.sum() < 16:      # too short to judge at all
        a = np.arange(len(t)) < len(t) // 2
        b = ~a
    wkw = dict(max_order=max(1, args.max_order),
               min_ratio=args.subharmonic_ratio,
               min_prominence=args.min_prominence)
    sa, sb = window_stats(t[a], cl[a], **wkw), window_stats(t[b], cl[b], **wkw)

    # Re-measure on window B, where the signal is closest to converged, and let
    # that be the reported answer.
    f_fund_b, f_car_b, order_b, sub_ratio_b = spectral_structure(
        t[b], cl[b], max_order=max(1, args.max_order),
        min_ratio=args.subharmonic_ratio, min_prominence=args.min_prominence)
    if f_car_b > 0:
        f_fund, f_car, order, sub_ratio = f_fund_b, f_car_b, order_b, sub_ratio_b
    else:
        f_fund, f_car = f_fund0, f_car0

    period = 1.0 / f_fund if f_fund > 0 else 0.0            # orbit period
    sample_period = 1.0 / f_car if f_car > 0 else 0.0       # carrier period
    cycles = (t[-1] - t_mid) * f_fund                       # FUNDAMENTAL cycles

    # Torus or closed orbit?  Only meaningful once there is a subharmonic.
    hdev = harmonic_deviation(t[b], cl[b], f_fund) if order > 1 else 0.0
    quasi = order > 1 and hdev > 0.05

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
    shedding = (sb["p2p"] > 2 * 5e-3 and f_car > 0
                and sb["prom"] >= args.min_prominence)

    def rel(a_, b_, scale):
        return abs(a_ - b_) / scale if scale > 0 else 0.0

    d_mean = rel(sa["mean"], sb["mean"], max(abs(sb["mean"]), sb["p2p"], 1e-9))
    d_p2p = rel(sa["p2p"], sb["p2p"], max(sb["p2p"], 1e-9))
    d_frq = rel(sa["freq"], sb["freq"], max(sb["freq"], 1e-9))

    # A period-m state has one more thing that must have stopped moving: the
    # subharmonic amplitude itself.  It is the slowest quantity in the problem
    # -- at AoA 30 the carrier was within 1% of saturation long before the
    # subharmonic was -- so mean/amplitude/frequency can all agree while the
    # orbit is still splitting.  Compare A(f/m)/A(f) between the two windows.
    d_sub = 0.0
    if order > 1:
        d_sub = rel(sa["ratio"], sub_ratio, max(sub_ratio, 1e-9))

    # Gate on whether the window was actually cut to length, NOT on recounting
    # cycles in it.  The window is 2*n_per fundamental periods BY CONSTRUCTION,
    # measured with the first-pass frequency f_fund0; `cycles` recounts the same
    # span with the window-B frequency f_fund, so it equals n_per * f_fund/f_fund0
    # and lands a hair under n_per whenever the refined estimate comes in slightly
    # LOW.  Measured on AoA 22 at t=60.5: 5.99 against a required 6.0, with
    # d(mean,p2p,freq,sub) = (0.00025, 0.00094, 0.00067, 0.00000) -- every
    # agreement test passed by two orders of magnitude and the run was blocked by
    # one hundredth of a cycle.  Which side of 6.0 it lands on is a coin flip:
    # AoA 20 and AoA 30 got 6.00 and 6.01 and passed.
    #
    # This is the same class of bug as the one recorded in ../FINDINGS.md section 8
    # (the old code sized the window at exactly 2*min_cycles/f0 and window B then
    # held marginally too few cycles, so the gate never fired).  It was fixed there
    # with a 1.15x margin; the margin was dropped during the period-doubling rework
    # on the reasoning that cutting whole periods made it unnecessary.  It did not:
    # the window construction and the cycle count still use different frequency
    # estimates.  Asking the question the check is actually for -- "is the window
    # the full length, or was it clamped by the start of the series?" -- removes
    # the disagreement instead of papering over it with a fudge factor.
    ok = (shedding
          and full_window
          and max(d_mean, d_p2p, d_frq, d_sub) <= args.tol)

    print(f"T_END={t[-1]:.6f}")
    print(f"SHEDDING={'yes' if shedding else 'no'}")
    print(f"ORDER={order}")                      # 1 = period-1, 2 = doubled, ...
    print(f"SUBHARM_RATIO={sub_ratio:.6f}")
    print(f"HARMONIC_DEV={hdev:.4f}")
    print(f"QUASIPERIODIC={'yes' if quasi else 'no'}")
    print(f"FREQ={f_fund:.8f}")                  # fundamental (orbit) frequency
    print(f"PERIOD={period:.8f}")                # orbit period  -> run length
    print(f"CARRIER_FREQ={f_car:.8f}")
    print(f"SAMPLE_PERIOD={sample_period:.8f}")  # carrier period -> write interval
    print(f"CYCLES_IN_WINDOW={cycles:.2f}")
    print(f"FULL_WINDOW={'yes' if full_window else 'no'}")
    print(f"CL_MEAN={sb['mean']:.6f}")
    print(f"CL_P2P={sb['p2p']:.6f}")
    print(f"PROMINENCE={sb['prom']:.2f}")
    print(f"CD_MEAN={cd[b].mean():.6f}")
    print(f"D_MEAN={d_mean:.5f}")
    print(f"D_P2P={d_p2p:.5f}")
    print(f"D_FREQ={d_frq:.5f}")
    print(f"D_SUB={d_sub:.5f}")
    print(f"CONVERGED={'yes' if ok else 'no'}")


if __name__ == "__main__":
    main()
