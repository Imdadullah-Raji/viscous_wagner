"""Zoom on the plate's leading-edge bubble, wall-normal axis stretched.

The bubble near its formation angle is a few thousandths of a chord thick, so at equal
aspect (streamlines_lsb.py) it is invisible.  Here y is measured from the flat upper
surface and the vertical axis is stretched.  Streamlines are drawn as contours of the
stream function psi(x, y) = int_wall^y u_x dy', exact over the flat part of the plate
(the wall is the psi = 0 line).  A closed bubble is the region psi < 0; its dividing
streamline is psi = 0.

    python plate_bubble_zoom.py --aoa 7.9 8 [--root plate_fine_re500] [--tag plate_fine]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from streamlines_lsb import ROOT, load, grid_field, upper_cf, crossings  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aoa", nargs="+", required=True)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent / "plate_fine_re500")
    ap.add_argument("--tag", default="plate_fine")
    ap.add_argument("--xlim", type=float, nargs=2, default=(0.02, 0.55))
    ap.add_argument("--h", type=float, default=0.02, help="height above the wall shown")
    a = ap.parse_args()

    out = ROOT / "lsb" / "figs" / a.tag
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(a.aoa), 1, figsize=(11, 2.6 * len(a.aoa) + 0.6),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, aoa in zip(axes, a.aoa):
        C, U, Cw, tau = load(a.root / f"steady_aoa{aoa}")
        yt = Cw[(Cw[:, 0] > 0.1) & (Cw[:, 0] < 0.9), 1].max()    # flat upper surface
        xs = np.linspace(*a.xlim, 900)
        ys = yt + np.linspace(0, a.h, 500)
        X, Y, ux, uy = grid_field(C, U, Cw, xs, ys)
        u = ux.filled(0.0)
        u[0] = 0.0                                                # no slip on the wall row
        dy = ys[1] - ys[0]
        psi = np.concatenate([np.zeros((1, len(xs))),
                              np.cumsum(0.5 * (u[1:] + u[:-1]) * dy, axis=0)])
        Yw = Y - yt
        pc = ax.pcolormesh(X, Yw, u, cmap="RdBu_r", shading="auto", rasterized=True,
                           norm=TwoSlopeNorm(vmin=-0.01, vcenter=0.0, vmax=0.15))
        pos = np.geomspace(2e-6, psi.max(), 22)
        ax.contour(X, Yw, psi, levels=pos, colors="0.2", linewidths=0.5)
        neg = psi.min()
        if neg < 0:
            ax.contour(X, Yw, psi, levels=np.linspace(neg, 0, 7)[:-1] + 1e-12, colors="k",
                       linewidths=0.7, linestyles="solid")
        ax.contour(X, Yw, psi, levels=[0.0], colors="yellow", linewidths=1.6)
        ax.contour(X, Yw, u, levels=[0.0], colors="k", linewidths=1.0, linestyles="dashed")
        x, cf = upper_cf(Cw, tau)
        cr = [c for c in crossings(x, cf) if c[0] < a.xlim[1]]
        for v, k in cr:
            ax.plot(v, 0, "v" if k == "separation" else "^", color="yellow", mec="k",
                    ms=9, zorder=5, clip_on=False)
        aft = (x > 0.02) & (x < 0.9)
        msg = (", ".join(f"{k[:3]}. {v:.3f}" for v, k in cr) or "attached") + \
              f";  min C_f = {cf[aft].min():+.5f} at x = {x[aft][np.argmin(cf[aft])]:.3f}"
        if neg < 0:
            msg += f";  bubble height (psi=0) {np.max(np.where(psi < 0, Yw, 0)):.4f}c"
        ax.set_title(f"α = {aoa}°:  {msg}", loc="left", fontsize=9)
        yc = np.unique(np.round(C[(np.abs(C[:, 0] - 0.25) < 2e-3) & (C[:, 1] > yt), 1] - yt, 6))
        for h in yc[yc < a.h]:
            ax.axhline(h, xmin=0.985, xmax=1.0, color="m", lw=1.2)
        ax.set_ylabel("height above wall / c")
        ax.set_ylim(0, a.h)
    axes[-1].set_xlabel("x/c")
    cb = fig.colorbar(pc, ax=axes, shrink=0.85, pad=0.01)
    cb.set_label(r"$u_x$  (dashed: $u_x=0$; yellow: $\psi=0$)")
    fig.suptitle("Plate, Re = 500, fine mesh: leading-edge bubble (vertical axis stretched, not to scale)", fontsize=10)
    name = "bubble_zoom_" + "_".join(a.aoa) + ".png"
    fig.savefig(out / name, dpi=200, bbox_inches="tight")
    print(out / name)


if __name__ == "__main__":
    main()
