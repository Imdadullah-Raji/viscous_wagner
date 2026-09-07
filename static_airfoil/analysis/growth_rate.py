#!/usr/bin/env python3
"""Exponential growth/decay rate of the shedding oscillation.

WHY THIS EXISTS.  Bracketing the shedding onset by running angles until one
stops shedding is both expensive and imprecise: near onset the growth rate
passes through zero (critical slowing down), so a case just below onset takes
hundreds of time units to visibly decay and a case just above takes as long to
visibly grow.  The quantity that actually locates the onset is the rate itself.

For a Hopf bifurcation the least-stable eigenvalue's real part sigma varies
LINEARLY with the control parameter near the critical point,

    sigma(alpha) ~ k * (alpha - alpha_c),

so two sub-critical angles with measured sigma < 0 extrapolate to alpha_c
without ever having to run the marginal case.

WHY A DECAYING CASE GIVES A TRUSTWORTHY sigma AND A GROWING ONE DOES NOT.
These runs are warm-started from a saturated shedding state, i.e. from a
FINITE-amplitude perturbation.  If the oscillation decays, its amplitude
shrinks toward zero and the late part of the decay is governed by the
linearisation about the steady state -- that tail is exactly the linear
eigenvalue, and this is the regime the fit targets.  If instead it grows and
saturates, a finite-amplitude seed cannot distinguish "the steady state is
unstable" from "a subcritical shedding branch coexists with a stable steady
state"; the number is then a nonlinear transient rate, not an eigenvalue.
Hence --tmin defaults to skipping the far-field readjustment, and R2 is
reported so a fit contaminated by saturation is visible rather than silent.

Usage:
    python3 analysis/growth_rate.py <caseDir> [--tmin T] [--tmax T] [--fit-from F]
Emits KEY=VALUE lines, eval-able from sh, like analysis/limit_cycle.py.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from limit_cycle import load_series, spectral_structure


def uniform(t, y, dt=None):
    dt = dt or float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, y)


def envelope(t, y, f0, n_per=6.0, step_frac=0.25):
    """Amplitude of the f0 component in a sliding Hann window.

    n_per periods per window: long enough that the f0 bin is resolved, short
    enough that a decay with an e-folding time of a few periods is not averaged
    flat.  The +/-1 bin search absorbs the small frequency drift a decaying
    oscillation shows as it relaxes toward the steady state.
    """
    tu, yu = uniform(t, y)
    dt = tu[1] - tu[0]
    win = int(round(n_per / f0 / dt))
    if win < 16 or win >= len(tu):
        return np.array([]), np.array([])
    step = max(1, int(win * step_frac))
    w = np.hanning(win)
    cg = w.sum() / win                      # coherent gain of the window
    f = np.fft.rfftfreq(win, dt)
    k0 = int(np.argmin(np.abs(f - f0)))
    tc, amp = [], []
    for i in range(0, len(tu) - win, step):
        seg = yu[i:i + win]
        seg = seg - seg.mean()
        Y = np.abs(np.fft.rfft(seg * w)) * 2.0 / (win * cg)
        tc.append(tu[i + win // 2])
        amp.append(float(Y[max(1, k0 - 1):k0 + 2].max()))
    return np.array(tc), np.array(amp)


def fit_exponential(tc, amp, frac):
    """Least squares on log(amp) over the last `frac` of the record."""
    good = amp > 0
    tc, amp = tc[good], amp[good]
    if len(tc) < 4:
        return None
    t0 = tc[0] + (1.0 - frac) * (tc[-1] - tc[0])
    sel = tc >= t0
    if sel.sum() < 4:
        sel = np.ones(len(tc), bool)
    x, ly = tc[sel], np.log(amp[sel])
    sigma, c = np.polyfit(x, ly, 1)
    resid = ly - (sigma * x + c)
    ss_tot = ((ly - ly.mean()) ** 2).sum()
    r2 = 1.0 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else 0.0
    return sigma, float(r2), x[0], x[-1], amp[sel][0], amp[sel][-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--tmin", type=float, default=10.0,
                    help="ignore before this time: a warm start needs a few "
                         "convective times for the far field to swing round to "
                         "the new angle, and that readjustment is not the mode "
                         "being measured")
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--fit-from", type=float, default=0.6,
                    help="fraction of the record (from the end) to fit; the "
                         "tail is where a decay is linear")
    ap.add_argument("--n-per", type=float, default=6.0)
    args = ap.parse_args()

    t, cl, _ = load_series(args.case)
    sel = t >= args.tmin
    if args.tmax is not None:
        sel &= t <= args.tmax
    t, cl = t[sel], cl[sel]
    if len(t) < 100:
        sys.exit(f"too few samples after t>={args.tmin}")

    # Carrier from the EARLY third, where the oscillation is largest.  Measuring
    # it on the tail of a decaying case reads the noise floor instead.
    head = t <= t[0] + 0.34 * (t[-1] - t[0])
    _, f_car, _, _ = spectral_structure(t[head], cl[head])
    if f_car <= 0:
        sys.exit("no carrier frequency in the early window")

    tc, amp = envelope(t, cl, f_car, n_per=args.n_per)
    if len(tc) < 4:
        sys.exit("record too short for a sliding envelope at this frequency")
    fit = fit_exponential(tc, amp, args.fit_from)
    if fit is None:
        sys.exit("envelope fit failed")
    sigma, r2, ta, tb, aa, ab = fit

    print(f"T_MIN={t[0]:.6g}")
    print(f"T_MAX={t[-1]:.6g}")
    print(f"CARRIER_FREQ={f_car:.8g}")
    print(f"A_FIRST={amp[0]:.6g}")
    print(f"A_LAST={amp[-1]:.6g}")
    print(f"A_RATIO={amp[-1] / amp[0] if amp[0] > 0 else 0:.6g}")
    print(f"FIT_T0={ta:.6g}")
    print(f"FIT_T1={tb:.6g}")
    print(f"FIT_A0={aa:.6g}")
    print(f"FIT_A1={ab:.6g}")
    print(f"SIGMA={sigma:.6g}")
    print(f"EFOLD={abs(1.0 / sigma) if sigma != 0 else float('inf'):.6g}")
    print(f"R2={r2:.6g}")
    print(f"TREND={'growing' if sigma > 0 else 'decaying'}")


if __name__ == "__main__":
    main()
