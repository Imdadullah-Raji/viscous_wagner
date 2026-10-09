"""Streamlines and reversed-flow regions at 7 deg, Re = 500: plate, NACA0006, NACA0012.

Run with: /home/raji/miniconda3/envs/research/bin/python lsb/streamlines_lsb.py --aoa 7 [8 9 ...]

Steady (SIMPLE) solutions at the latest iteration.  Needs the cell centres:
    foamPostProcess -func writeCellCentres -latestTime      (in each case)

Geometry stays at zero incidence (AoA is imposed by rotating the freestream), so the
mesh x-axis IS the chord direction and u_x < 0 is flow reversed relative to the chord.

  * velocity from the cell centres plus the wall-face centres (U = 0, no-slip),
    triangulated; triangles with an edge midpoint inside the body are removed so no
    interpolation crosses it, then sampled on a uniform grid for streamplot.
  * body outline: the wall-face centres sorted by angle about their centroid (all three
    sections are convex).
  * C_f on the upper surface from the wallShearStress function object, NEGATED:
    it returns -devTau & n (wallShearStress.C:82), the traction of the wall on the fluid.
    C_f = 2 tau_x (rho = U = 1).  Zero crossings = separation / reattachment.
A section with no case at that angle is skipped and reported (the NACA0006 sweep is
integer-only).
-> lsb/figs/aoa<aoa>/
"""
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[1]
THESIS = ROOT.parent
CASES = {"plate": THESIS / "plate_case" / "re500",
         "NACA0006": ROOT / "NACA0006" / "steady_re500",
         "NACA0012": ROOT / "NACA0012" / "steady_re500"}
COLOR = {"plate": "#c2410c", "NACA0006": "#1baf7a", "NACA0012": "#2a78d6"}


def _vectors(block):
    return np.array([[float(v) for v in m.split()]
                     for m in re.findall(r"\(([^()]*)\)", block)])


def read_internal(path):
    s = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", s)
    n = int(m.group(1))
    body = s[m.end(): s.index("\n)\n", m.end())]
    a = _vectors(body)
    assert len(a) == n, (path, len(a), n)
    return a


def read_patch(path, patch="airfoil"):
    s = path.read_text()
    i = s.index("boundaryField")
    j = re.search(rf"\b{patch}\s*\{{", s[i:]).end() + i
    m = re.search(r"value\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", s[j:])
    n, k = int(m.group(1)), j + m.end()
    a = _vectors(s[k: s.index("\n)\n", k)])
    assert len(a) == n, (path, len(a), n)
    return a


def latest(case):
    """Latest time directory (the plate's extended runs end at 8000-14000, not 6000)."""
    return max((p for p in case.iterdir() if p.is_dir() and re.fullmatch(r"[\d.]+", p.name)
                and p.name != "0"), key=lambda p: float(p.name))


def load(case):
    d = latest(case)
    C, U = read_internal(d / "C")[:, :2], read_internal(d / "U")[:, :2]
    Cw = read_patch(d / "C")[:, :2]
    tau = -read_patch(d / "wallShearStress")[:, :2]          # physical traction on wall
    c0 = Cw.mean(axis=0)
    order = np.argsort(np.arctan2(Cw[:, 1] - c0[1], Cw[:, 0] - c0[0]))
    return C, U, Cw[order], tau[order]


def grid_field(C, U, Cw, xs, ys):
    pts = np.vstack([C, Cw])
    uv = np.vstack([U, np.zeros_like(Cw)])
    tri = mtri.Triangulation(pts[:, 0], pts[:, 1])
    body = MPath(Cw)
    T = tri.triangles
    mids = np.concatenate([(pts[T[:, a]] + pts[T[:, b]]) / 2 for a, b in ((0, 1), (1, 2), (2, 0))])
    inside = body.contains_points(mids, radius=-1e-9).reshape(3, -1).any(axis=0)
    tri.set_mask(inside)
    X, Y = np.meshgrid(xs, ys)
    ux = mtri.LinearTriInterpolator(tri, uv[:, 0])(X, Y)
    uy = mtri.LinearTriInterpolator(tri, uv[:, 1])(X, Y)
    inb = body.contains_points(np.c_[X.ravel(), Y.ravel()]).reshape(X.shape)
    ux = np.ma.masked_where(inb | np.ma.getmaskarray(ux), ux)
    uy = np.ma.masked_where(inb | np.ma.getmaskarray(uy), uy)
    return X, Y, ux, uy


def upper_cf(Cw, tau):
    """Upper-surface C_f(x), from the leading edge aft, chord-wise component."""
    up = Cw[:, 1] > 0
    o = np.argsort(Cw[up, 0])
    return Cw[up, 0][o], 2 * tau[up, 0][o]


def crossings(x, f, xmin=0.02):
    s = np.sign(f)
    i = np.where((s[:-1] != s[1:]) & (x[:-1] > xmin))[0]
    xc = x[i] - f[i] * (x[i + 1] - x[i]) / (f[i + 1] - f[i])
    return [(float(a), "separation" if f[k] > 0 else "reattachment") for a, k in zip(xc, i)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aoa", nargs="+", default=["7"])
    ap.add_argument("--only", nargs="+", choices=list(CASES), help="subset of sections")
    ap.add_argument("--plate-root", type=Path, help="plate case folder (default plate_case/re500)")
    ap.add_argument("--tag", default="", help="figures go to figs/<tag>/aoa<A> instead of figs/aoa<A>")
    a = ap.parse_args()
    if a.plate_root:
        CASES["plate"] = a.plate_root.resolve()
    for name in set(CASES) - set(a.only or CASES):
        del CASES[name]
    for aoa in a.aoa:
        plot_angle(aoa, a.tag)


def plot_angle(aoa, tag=""):
    out = ROOT / "lsb" / "figs" / tag / f"aoa{aoa}"
    out.mkdir(parents=True, exist_ok=True)

    data, lines = {}, []
    for name, base in CASES.items():
        case = base / f"steady_aoa{aoa}"
        if not case.is_dir():
            lines.append(f"{name}: no case at {aoa} deg -- skipped")
            continue
        C, U, Cw, tau = load(case)
        x, cf = upper_cf(Cw, tau)
        cr = crossings(x, cf)
        data[name] = (C, U, Cw, x, cf, cr)
        rev = cf[x > 0.02].min()
        lines.append(f"{name}: min upper-surface C_f (x > 0.02) = {rev:+.5f};  zero crossings (x/c): "
                     + (", ".join(f"{k} at {v:.4f}" for v, k in cr) or "none (attached)"))
    txt = "\n".join(lines)
    (out / "separation.txt").write_text(txt + "\n")
    print(txt)

    views = {"full": ((-0.15, 1.45), (-0.25, 0.35)), "aft": ((0.35, 1.2), (-0.06, 0.16)),
             "le": ((-0.03, 0.4), (-0.04, 0.10)), "te": ((0.93, 1.07), (-0.02, 0.045))}
    for view, (xl, yl) in views.items():
        nx = 700
        ny = int(nx * (yl[1] - yl[0]) / (xl[1] - xl[0]))
        xs, ys = np.linspace(*xl, nx), np.linspace(*yl, ny)
        fig, axes = plt.subplots(len(data), 1, squeeze=False, figsize=(12, 3.6 * len(data) * (yl[1] - yl[0])
                                                          / (xl[1] - xl[0]) * 2.6 + 1.2))
        axes = axes[:, 0]
        for ax, (name, (C, U, Cw, x, cf, cr)) in zip(axes, data.items()):
            X, Y, ux, uy = grid_field(C, U, Cw, xs, ys)
            # diverging about u_x = 0: blue is reversed flow, red is forward
            norm = TwoSlopeNorm(vmin=-0.1, vcenter=0.0, vmax=1.4)
            pc = ax.pcolormesh(X, Y, ux, cmap="RdBu_r", norm=norm, shading="auto",
                               rasterized=True)
            ax.contour(X, Y, ux.filled(np.nan), levels=[0], colors="black", linewidths=1.3)
            ax.streamplot(xs, ys, ux.filled(np.nan), uy.filled(np.nan), color="0.15",
                          density=2.6 if view != "full" else 2.0, linewidth=0.55, arrowsize=0.6)
            ax.fill(Cw[:, 0], Cw[:, 1], color="0.6", zorder=3)
            for v, k in cr:
                ax.plot(v, np.interp(v, x, Cw[Cw[:, 1] > 0][np.argsort(Cw[Cw[:, 1] > 0][:, 0]), 1]),
                        "v" if k == "separation" else "^", color="yellow", mec="black",
                        ms=8, zorder=4)
            ax.set(xlim=xl, ylim=yl, aspect="equal")
            ax.set_title(f"{name}, Re = 500, α = {aoa}°  (steady)", loc="left", fontsize=10)
            ax.set_ylabel("y/c")
        axes[-1].set_xlabel("x/c (chord-aligned; freestream at α)")
        cb = fig.colorbar(pc, ax=axes, shrink=0.8, pad=0.01)
        cb.set_label(r"$u_x$ (chord-wise); black line: $u_x = 0$")
        fig.savefig(out / f"streamlines_{view}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4.2))
    for name, (C, U, Cw, x, cf, cr) in data.items():
        ax.plot(x, cf, color=COLOR[name], lw=1.6, label=name)
        for v, k in cr:
            ax.plot(v, 0, "v" if k == "separation" else "^", color=COLOR[name], mec="black", ms=8)
    ax.axhline(0, color="black", lw=0.7)
    ax.set(xlabel="x/c", ylabel=r"$C_f$ (upper surface, chord-wise)", xlim=(0, 1.02),
           ylim=(-0.05, 0.25))
    ax.set_title(f"Upper-surface skin friction, α = {aoa}°, Re = 500 "
                 "(▼ separation, ▲ reattachment)", loc="left")
    ax.grid(alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "cf_upper.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
