"""Tight vs production C_L for the ZTP, per Re: convergence tail, offset, odd Legendre fits.

  python analysis/tight_compare.py 500 700 900 [--drop]   (research env)
  --drop: leave out unconverged angles individually instead of capping the windows at the first one
"""
import sys; sys.path.insert(0, "/home/raji/Research/Thesis/static_airfoil/analysis")
from pathlib import Path
import numpy as np, pandas as pd
from lasso_re1000 import load
from cubic_stress_test import fit_window
T = Path("/home/raji/Research/Thesis/static_airfoil/zero_thickness_plate/tight")

def cases(re):
    rows = []
    for c in sorted(T.glob(f"re{re}/steady_aoa*"), key=lambda p: float(p.name[10:]) if p.is_dir() else 0):
        if not c.is_dir() or not list(c.glob("log.steady_to*")): continue
        lg = c / "log.steady"
        if not lg.exists() or b"\nEnd" not in lg.read_bytes(): continue
        it, cd, cl = load(str(c))
        it = it.astype(int); i = {k: j for j, k in enumerate(it)}
        starts = sorted(int(d.name) for d in (c / "postProcessing/forceCoeffs").iterdir() if d.name != "0")
        t0, tN = starts[0], it[-1]
        v = np.array([cl[i[k]] for k in range(tN - 2500, tN + 1, 500)])
        d = np.diff(v); q = d[-1] / d[-2] if d[-2] else np.nan
        rem = d[-1] * q / (1 - q) if 0 < q < 0.9 else np.nan
        rows.append(dict(aoa=float(c.name[10:]), cl=cl[-1], prod=cl[i[t0]], t0=t0, tN=tN,
                         d1000=cl[-1] - cl[i[tN - 1000]], d500=d[-1], q=q, rem=rem))
    return pd.DataFrame(rows)

DROP = "--drop" in sys.argv
for re in [a for a in sys.argv[1:] if not a.startswith("--")]:
    df = cases(re)
    df["off"] = df.cl - df["prod"]
    # usable: last-500 change and extrapolated remainder both <= 3e-9 (target ~1e-9)
    df["ok"] = (df.d500.abs() <= 3e-9) & ((df.rem.abs() <= 3e-9) | df.rem.isna() & (df.d500.abs() <= 1e-9))
    print(f"\n######## Re {re}: {len(df)} tight cases, {df.ok.sum()} usable")
    print(" aoa    iters        Cl(final)      tight-prod   d(last1000)  d(last500)  ratio  remain   ok")
    for _, r in df.iterrows():
        print(f" {r.aoa:4g}  {r.t0:5d}->{r.tN:5d}  {r.cl:.12f}  {r.off:+.2e}   {r.d1000:+.1e}    {r.d500:+.1e}  {r.q:+.2f}  {r.rem:+.1e}  {'Y' if r.ok else '-'}")
    U = df[df.ok]; Pd = df.assign(cl=df["prod"])
    amax = U.aoa.max()
    # first unusable angle caps the windows (no holes)
    bad = df[~df.ok].aoa
    cap = amax if DROP else (min(amax, bad.min() - 1e-9) if len(bad) else amax)
    U = U[U.aoa <= cap]; Pd = Pd[df.ok & (Pd.aoa <= cap)]
    step = np.diff(U.aoa).min()
    As = [a for a in np.arange(2.0, cap + 1e-9, 0.5)]
    for name, d in (("TIGHT", U), ("PRODUCTION (same angles)", Pd)):
        print(f"\n== Re {re} {name}: {len(d)} angles {d.aoa.min():g}-{d.aoa.max():g} step {step:g}")
        for degs in ((1, 3, 5, 7), (1, 3, 5, 7, 9), (1, 3, 5, 7, 9, 11)):
            for A in As:
                n = ((d.aoa > 0) & (d.aoa <= A + 1e-9)).sum()
                if n < len(degs) + 2: continue
                a, se, rms, n = fit_window(d, A, degs)
                print(f"  P1-P{max(degs):<2} 0-{A:<4g} n={n:2d} dof={n-len(degs)} rms={rms:.1e} " + " ".join(
                    f"a{k}={a[k]:+.6g}±{se[k]:.1g}" for k in degs) +
                    f" | R35={np.degrees(np.sqrt(abs(a[3]/a[5]))):.1f}° R57={np.degrees(np.sqrt(abs(a[5]/a[7]))):.1f}°")
