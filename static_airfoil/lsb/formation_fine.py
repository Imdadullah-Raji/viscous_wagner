"""Plate LSB formation, Re = 500, fine mesh: the 0.01-deg chop 7.90-8.00 (oct9_workflow.md).

For every steady_aoa<A> with 7.895 <= A <= 8.005 in lsb/plate_fine_re500/, at every written
iteration >= 3000 (the 7.91-7.99 cases write every 1000; 7.9 and 8.0 only at 6000):
min upper-surface C_f on 0.02 < x < 0.9, its location, separation, reattachment, length.
The change between the last two writes is the iteration-convergence error.

Then: alpha_f from a linear fit of C_f,min(alpha) through zero, and the bubble length
against alpha - alpha_f on log-log axes (local exponent: 1/2 for a parabolic C_f minimum,
1 for the linear growth seen from 8.0 deg up).

Writes figs/plate_fine/formation.{txt,csv,png}.   Run with the research env.
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from streamlines_lsb import ROOT, read_patch, upper_cf, crossings  # noqa: E402

FINE = ROOT / "lsb" / "plate_fine_re500"
OUT = ROOT / "lsb" / "figs" / "plate_fine"


def at_time(d):
    Cw = read_patch(d / "C")[:, :2]
    tau = -read_patch(d / "wallShearStress")[:, :2]
    c0 = Cw.mean(axis=0)
    o = np.argsort(np.arctan2(Cw[:, 1] - c0[1], Cw[:, 0] - c0[0]))
    x, cf = upper_cf(Cw[o], tau[o])
    m = (x > 0.02) & (x < 0.9)
    cr = crossings(x[m], cf[m])
    s = [v for v, k in cr if k == "separation"]
    r = [v for v, k in cr if k == "reattachment" and s and v > s[0]]
    xs, xr = (s[0], r[0]) if s and r else (np.nan, np.nan)
    return float(cf[m].min()), float(x[m][np.argmin(cf[m])]), xs, xr


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for p in sorted(FINE.glob("steady_aoa*")):
        m = re.fullmatch(r"steady_aoa([\d.]+)", p.name)
        if not m or not 7.895 <= float(m.group(1)) <= 8.005:
            continue
        times = sorted((d for d in p.iterdir() if d.is_dir() and re.fullmatch(r"\d+", d.name)
                        and int(d.name) >= 3000 and (d / "C").exists()
                        and (d / "wallShearStress").exists()), key=lambda d: int(d.name))
        if not times:
            print(f"skip {p.name}: no post-processed time")
            continue
        hist = [(int(d.name),) + at_time(d) for d in times]
        rows.append((float(m.group(1)), hist))
    rows.sort()

    lines = [f"{'aoa':>5} {'iter':>5} {'Cf_min':>10} {'x(min)':>7} {'x_sep':>7} {'x_reat':>7} {'L':>7}"]
    csv = ["aoa,iter,cf_min,x_cfmin,x_sep,x_reat,length"]
    for a, hist in rows:
        for it, cmin, xm, xs, xr in hist:
            lines.append(f"{a:5.2f} {it:5d} {cmin:10.6f} {xm:7.3f} {xs:7.4f} {xr:7.4f} {xr - xs:7.4f}")
            csv.append(f"{a},{it},{cmin},{xm},{xs},{xr},{xr - xs}")
        if len(hist) > 1:
            d = hist[-1][1] - hist[-2][1]
            dL = (hist[-1][4] - hist[-1][3]) - (hist[-2][4] - hist[-2][3])
            lines.append(f"{'':5} last-write change: dCf_min {d:+.1e}, dL {dL:+.4f}")

    A = np.array([a for a, _ in rows])
    cmin = np.array([h[-1][1] for _, h in rows])
    L = np.array([h[-1][4] - h[-1][3] for _, h in rows])
    p = np.polyfit(A, cmin, 1)
    af = -p[1] / p[0]
    res = cmin - np.polyval(p, A)
    lines += ["", f"C_f,min = {p[0]:.5f} (alpha - {af:.4f}),  fit max |residual| {np.abs(res).max():.1e}",
              f"alpha_f (C_f,min = 0) = {af:.4f} deg"]
    ok = np.isfinite(L) & (A > af)
    if ok.sum() >= 2:
        da, LL = A[ok] - af, L[ok]
        lines.append("local exponent d ln L / d ln(alpha - alpha_f) between successive angles:")
        for i in range(len(da) - 1):
            k = np.log(LL[i + 1] / LL[i]) / np.log(da[i + 1] / da[i])
            lines.append(f"  {A[ok][i]:.2f}-{A[ok][i + 1]:.2f}: {k:.2f}")
    (OUT / "formation.txt").write_text("\n".join(lines) + "\n")
    (OUT / "formation.csv").write_text("\n".join(csv) + "\n")
    print("\n".join(lines))

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(A, cmin, "o-k")
    ax[0].axhline(0, color="0.6", lw=0.8); ax[0].axvline(af, color="r", lw=0.8, ls="--")
    ax[0].set(xlabel=r"$\alpha$ (deg)", ylabel=r"min upper $C_f$", title=f"α_f = {af:.3f}°")
    xs = np.array([h[-1][3] for _, h in rows]); xr = np.array([h[-1][4] for _, h in rows])
    ax[1].plot(A, xs, "v-", color="#c2410c", label="separation")
    ax[1].plot(A, xr, "^-", color="#2a78d6", label="reattachment")
    ax[1].axvline(af, color="r", lw=0.8, ls="--")
    ax[1].set(xlabel=r"$\alpha$ (deg)", ylabel="x/c", title="bubble ends"); ax[1].legend()
    if ok.sum() >= 2:
        ax[2].loglog(da, LL, "o-k", label="fine mesh")
        g = np.geomspace(da.min(), da.max(), 50)
        ax[2].loglog(g, LL[0] * (g / da[0]) ** 0.5, "--", color="#1baf7a", label=r"$\propto(\alpha-\alpha_f)^{1/2}$")
        ax[2].loglog(g, LL[-1] * (g / da[-1]), ":", color="#c2410c", label=r"$\propto(\alpha-\alpha_f)^{1}$")
        ax[2].set(xlabel=r"$\alpha-\alpha_f$ (deg)", ylabel="bubble length / c", title="growth law")
        ax[2].legend(fontsize=8)
    fig.suptitle("Plate, Re = 500, fine mesh: LSB formation (steady)")
    fig.tight_layout(); fig.savefig(OUT / "formation.png", dpi=150)


if __name__ == "__main__":
    main()
