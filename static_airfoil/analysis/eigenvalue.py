#!/usr/bin/env python3
"""Linear growth rate of a perturbation about a STEADY base flow.

The case this reads was started at (or very near) a fixed point, so
    q(t) = C_L(t) - C_L_steady
is a perturbation evolving under the linearised operator L(alpha).  Its late-time
behaviour is e^{sigma t}, and sigma is the real part of the leading eigenvalue.
sigma > 0 unstable, sigma < 0 stable, and the angle where it crosses zero is the
bifurcation point -- with no assumption about what KIND of bifurcation it is.

DELIBERATELY DOES NOT ASSUME THE MODE OSCILLATES.  A Hopf bifurcation has a
complex pair and q(t) rings at omega; a real eigenvalue crossing (a pitchfork or
saddle-node, which is what a non-Hopf onset would look like) has no oscillation
at all and q(t) simply grows or decays.  Both are measured and reported:

    SIGMA_OSC   from the amplitude of the dominant spectral component
    SIGMA_MEAN  from |window mean of q|, the NON-oscillatory part
    OSC_FRAC    rms(oscillatory) / rms(total) -- near 0 means no ringing,
                so the onset is not a Hopf

Report the local d(ln)/dt table too, not just one fitted slope: a constant local
rate over several e-folds is an eigenvalue, a drifting one is not.  CHECK THE
SPAN OF THE FIT WINDOW BEFORE THE R^2 -- an R^2 of 0.999 over 6 time units taught
nothing (FINDINGS.md section 17).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from limit_cycle import load_series


def uniform(t, y):
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, y), dt


def windows(tu, qu, dt, win_t, step_frac=0.1, f_track=None, f_min=0.2):
    # step_frac 0.1, not 0.2: with a 10 t window a t=5-20.5 record gave only 3
    # window positions and the 4-point minimum in fit() returned nan.  Windows
    # overlap heavily, so they are not independent -- that inflates R^2 and is
    # exactly why the SPAN of the fit, not its R^2, is the thing to check.
    win = max(16, int(round(win_t / dt)))
    if win >= len(tu):
        return None
    step = max(1, int(win * step_frac))
    w = np.hanning(win)
    cg = w.sum() / win
    f = np.fft.rfftfreq(win, dt)
    tc, a_osc, a_mean, rms_tot = [], [], [], []
    for i in range(0, len(tu) - win, step):
        seg = qu[i:i + win]
        m = seg.mean()
        d = seg - m
        Y = np.abs(np.fft.rfft(d * w)) * 2.0 / (win * cg)
        Y[0] = 0.0
        tc.append(tu[i + win // 2])
        # Track ONE frequency, not the tallest peak anywhere.  Y.max() reported
        # sigma=0.21 at AoA 11.5 because a slowly-decaying non-oscillatory mode
        # from the seed made a single large excursion that dwarfed the f=0.600
        # eigenmode -- the fit was measuring the transient, not the eigenvalue.
        # Frequencies below f_min are drift filling the lowest FFT bins (an
        # apparent "peak" at 1/window is exactly that) and are excluded.
        if f_track is not None:
            k = int(np.argmin(np.abs(f - f_track)))
            a_osc.append(float(Y[max(1, k - 1):k + 2].max()))
        else:
            band = f >= f_min
            a_osc.append(float(Y[band].max()) if band.any() else 0.0)
        a_mean.append(abs(float(m)))
        rms_tot.append(float(np.sqrt((seg ** 2).mean())))
    return (np.array(tc), np.array(a_osc), np.array(a_mean), np.array(rms_tot),
            f[int(np.argmax(np.abs(np.fft.rfft((qu[:win] - qu[:win].mean()) * w))[1:])) + 1])


def fit(tc, a, frac):
    good = a > 0
    tc, a = tc[good], a[good]
    if len(tc) < 4:
        return None
    t0 = tc[0] + (1 - frac) * (tc[-1] - tc[0])
    sel = tc >= t0
    if sel.sum() < 4:
        sel = np.ones(len(tc), bool)
    x, ly = tc[sel], np.log(a[sel])
    s, c = np.polyfit(x, ly, 1)
    r = ly - (s * x + c)
    ss = ((ly - ly.mean()) ** 2).sum()
    return s, (1 - (r ** 2).sum() / ss if ss > 0 else 0.0), x[0], x[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--cl-steady", type=float, required=True,
                    help="C_L of the base flow this case was started from")
    ap.add_argument("--tmin", type=float, default=5.0,
                    help="skip the first few t: mapFields transfers volume "
                         "fields only, so phi is rebuilt and the first pressure "
                         "solve absorbs a divergence error that is not the mode")
    ap.add_argument("--win", type=float, default=10.0, help="window length, time units")
    ap.add_argument("--fit-from", type=float, default=0.6)
    ap.add_argument("--table", action="store_true", help="print the local-rate table")
    ap.add_argument("--freq", type=float, default=None,
                    help="ALWAYS PASS THIS when you know the mode frequency. "
                         "Tracks the amplitude at this frequency instead of the "
                         "tallest peak. Without it, a transient at another "
                         "frequency hijacks the fit: at AoA 11.5 the tallest "
                         "peak jumped from a decaying excursion at f=0.2 to the "
                         "f=0.6 eigenmode between consecutive windows and the "
                         "fit read sigma=-0.10 where the mode was at +0.0067.")
    ap.add_argument("--f-min", type=float, default=0.2,
                    help="ignore spectral content below this when auto-picking; "
                         "the lowest FFT bins are drift, not a mode")
    args = ap.parse_args()

    t, cl, _ = load_series(args.case)
    sel = t >= args.tmin
    t, cl = t[sel], cl[sel]
    if len(t) < 200:
        sys.exit(f"too few samples after t>={args.tmin}")
    q = cl - args.cl_steady
    tu, qu, dt = uniform(t, q)
    w = windows(tu, qu, dt, args.win, f_track=args.freq, f_min=args.f_min)
    if w is None:
        sys.exit(f"record {tu[-1]-tu[0]:.1f} t is shorter than the {args.win} t window")
    tc, a_osc, a_mean, rms_tot, f0 = w

    print(f"T_MIN={t[0]:.6g}")
    print(f"T_MAX={t[-1]:.6g}")
    print(f"SPAN={t[-1]-t[0]:.6g}")
    print(f"CL_STEADY={args.cl_steady:.6g}")
    print(f"Q_FIRST={q[0]:.6g}")
    print(f"Q_LAST={q[-1]:.6g}")
    print(f"FREQ={f0:.6g}")
    print(f"OSC_FRAC={np.mean(a_osc / np.maximum(rms_tot, 1e-30)):.4g}")
    for name, arr in (("OSC", a_osc), ("MEAN", a_mean), ("RMS", rms_tot)):
        r = fit(tc, arr, args.fit_from)
        if r is None:
            print(f"SIGMA_{name}=nan"); continue
        s, r2, t0, t1 = r
        print(f"SIGMA_{name}={s:.6g}")
        print(f"R2_{name}={r2:.5f}")
        print(f"FITSPAN_{name}={t1-t0:.6g}")
    print(f"A_OSC_FIRST={a_osc[0]:.6g}")
    print(f"A_OSC_LAST={a_osc[-1]:.6g}")

    if args.table:
        print("\n#    t      a_osc      |mean q|     rms      local d(ln rms)/dt")
        for i in range(len(tc)):
            d = ""
            if i >= 2:
                k = np.polyfit(tc[i-2:i+1], np.log(np.maximum(rms_tot[i-2:i+1], 1e-30)), 1)[0]
                d = f"{k:+.5f}"
            print(f"# {tc[i]:7.2f}  {a_osc[i]:.6e}  {a_mean[i]:.6e}  {rms_tot[i]:.6e}  {d}")


if __name__ == "__main__":
    main()
