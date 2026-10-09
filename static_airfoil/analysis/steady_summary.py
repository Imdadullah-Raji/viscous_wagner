#!/usr/bin/env python3
"""Collect the steady (SIMPLE) base-flow solves and judge each one's convergence.

CONVERGENCE IS JUDGED ON THE FORCE TRACE, NOT THE RESIDUALS.  A residual can
stall low while C_L is still creeping, and -- more importantly here -- for a
strongly unstable base flow SIMPLE does not converge at all: the physical
instability leaks into the iteration and C_L oscillates in iteration space with
a constant amplitude.  Measured at AoA 24: the C_L range over the last 50%, 30%,
20% and 10% of iterations is 1.258e-02 in every window.  An identical range
across shrinking windows is the signature -- a converging case's range must fall
as the window shrinks.

Verdicts:
  CONVERGED  range over the last 10% below `tol` (default 1e-4 of |C_L|)
  CREEPING   range still falling as the window shrinks, but not yet below tol
  LIMIT-CYCLE  range flat across shrinking windows -- no fixed point found
"""
import argparse
import glob
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# unsteady saturated means at Re=500 ONLY, for the mean-flow-distortion comparison
UNSTEADY = {14: (0.551601, 0.261941), 16: (0.630309, 0.304484),
            18: (0.711513, 0.356784), 20: (0.786700, 0.415700),
            22: (0.839900, 0.474200), 24: (0.879000, 0.532000)}


def load(case):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", "forceCoeffs",
                                          "*", "*.dat")))
    it, cd, cl = [], [], []
    for f in files:
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            try:
                it.append(float(v[0])); cd.append(float(v[2])); cl.append(float(v[3]))
            except (ValueError, IndexError):
                pass
    if not it:
        return None
    o = np.argsort(it, kind="stable")
    return np.array(it)[o], np.array(cd)[o], np.array(cl)[o]


def verdict(it, cl, tol):
    rng = {}
    for frac in (0.5, 0.2, 0.1):
        m = it >= (1 - frac) * it[-1] + frac * it[0]
        rng[frac] = cl[m].max() - cl[m].min()
    # Absolute floor, not just |C_L|: at AoA 0 the lift is exactly zero by
    # symmetry, and a purely relative tolerance then measures the range against
    # ~1e-6 and calls a solve converged to 2e-10 "creeping".
    scale = max(abs(cl[-1]), 0.01)
    if rng[0.1] <= tol * scale:
        return "CONVERGED", rng
    # a converging trace must shrink its range as the window shrinks
    if rng[0.1] < 0.5 * rng[0.5]:
        return "CREEPING", rng
    return "LIMIT-CYCLE", rng


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=1e-4,
                    help="fractional C_L range over the last 10%% of iterations")
    ap.add_argument("--root", default=os.path.join(ROOT, "NACA0012", "steady_re500"),
                    help="directory holding steady_aoa<A> cases (default NACA0012/steady_re500)")
    args = ap.parse_args()
    # the unsteady comparison exists only for the NACA0012 Re=500 campaign
    is_re500 = os.path.realpath(args.root) == os.path.realpath(os.path.join(ROOT, "NACA0012", "steady_re500"))

    def angle(path):
        try:
            return float(path.split("steady_aoa")[-1])
        except ValueError:
            return None

    print(f"{'AoA':>5} {'iters':>6} {'Cl_steady':>11} {'Cd_steady':>11} "
          f"{'rng(50%)':>10} {'rng(10%)':>10}  verdict")
    rows = []
    cases = [c for c in glob.glob(os.path.join(args.root, "steady_aoa*"))
             if os.path.isdir(c) and angle(c) is not None]
    for case in sorted(cases, key=angle):
        aoa = angle(case)
        aoa = int(aoa) if aoa.is_integer() else aoa
        d = load(case)
        if d is None:
            print(f"{aoa:>5}      -           -           -          -          -  NO DATA")
            continue
        it, cd, cl = d
        v, rng = verdict(it, cl, args.tol)
        print(f"{aoa:>5} {int(it[-1]):>6} {cl[-1]:>11.6f} {cd[-1]:>11.6f} "
              f"{rng[0.5]:>10.3e} {rng[0.1]:>10.3e}  {v}")
        rows.append((aoa, cl[-1], cd[-1], v))

    if not is_re500:
        return
    print(f"\n{'AoA':>4} {'Cl_steady':>11} {'Cl_mean':>10} {'Cl_mean-Cl_steady':>19} "
          f"{'as % of steady':>15}  verdict")
    for aoa, cls, cds, v in rows:
        if aoa not in UNSTEADY:
            continue
        clm = UNSTEADY[aoa][0]
        print(f"{aoa:>4} {cls:>11.6f} {clm:>10.6f} {clm - cls:>19.6f} "
              f"{100 * (clm - cls) / abs(cls):>14.1f}%  {v}")


if __name__ == "__main__":
    main()
