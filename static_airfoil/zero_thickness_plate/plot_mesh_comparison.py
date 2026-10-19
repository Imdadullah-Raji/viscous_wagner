"""H-type (ztp_mesher, T3) against elliptic-coordinate (ell_mesher, E800) grids for the zero-thickness
plate, Re = 500: near field and leading-edge tip, with face non-orthogonality.

  /home/raji/miniconda3/envs/research/bin/python zero_thickness_plate/plot_mesh_comparison.py
  -> zero_thickness_plate/figs/mesh_comparison.png
"""
from collections import defaultdict
from pathlib import Path

import gmsh
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = Path(__file__).resolve().parent
MESHES = {"H-type grid, T3 (superseded)": HERE / "meshes/ztp_T3_tip0.000559017_np240_nsp140_re500.msh",
          "elliptic grid, E800 (used)": HERE / "meshes/ell_nnu800_nmu210_re500.msh"}


def quads(msh):
    gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 1); gmsh.open(str(msh))
    tags, xyz, _ = gmsh.model.mesh.getNodes(); X = xyz.reshape(-1, 3)
    idx = np.zeros(int(tags.max()) + 1, int); idx[tags.astype(int)] = np.arange(len(tags))
    _, hn = gmsh.model.mesh.getElementsByType(5); H = idx[hn.astype(int)].reshape(-1, 8)
    gmsh.finalize()
    Q = np.array([[n for n in h if abs(X[n, 2]) < 1e-9] for h in H])
    return X[:, :2][Q]


def nonortho(Q):
    """Max face non-orthogonality of each cell (deg), OpenFOAM definition, in the plane."""
    C = Q.mean(1); E = defaultdict(list)
    for c, q in enumerate(np.round(Q, 12)):
        for i in range(4):
            a, b = tuple(q[i]), tuple(q[(i + 1) % 4]); E[(a, b) if a < b else (b, a)].append(c)
    worst = np.zeros(len(Q))
    for (a, b), cs in E.items():
        if len(cs) != 2: continue
        a, b = np.array(a), np.array(b); t = b - a; n = np.array([t[1], -t[0]]) / np.linalg.norm(t)
        d = C[cs[1]] - C[cs[0]]
        ang = np.degrees(np.arccos(min(1, abs(d @ n) / np.linalg.norm(d))))
        worst[cs] = np.maximum(worst[cs], ang)
    return worst


def main():
    (HERE / "figs").mkdir(exist_ok=True)
    views = [("near field", (-0.4, 1.4, -0.7, 0.7)), ("leading-edge tip, 0.04c", (-0.02, 0.02, -0.015, 0.015))]
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    for row, (name, msh) in enumerate(MESHES.items()):
        Q = quads(msh); no = nonortho(Q)
        segs = np.concatenate([np.stack([Q[:, i], Q[:, (i + 1) % 4]], 1) for i in range(4)])
        for col, (vname, (x0, x1, y0, y1)) in enumerate(views):
            a = ax[row, col]; pad = 0.2 * (x1 - x0)
            m = ((segs[..., 0] > x0 - pad) & (segs[..., 0] < x1 + pad) & (segs[..., 1] > y0 - pad) & (segs[..., 1] < y1 + pad)).any(1)
            a.add_collection(LineCollection(segs[m], lw=0.15, color="k"))
            a.plot([0, 1], [0, 0], "r-", lw=1.4)
            a.set(xlim=(x0, x1), ylim=(y0, y1), title=f"{name}: {vname}"); a.set_aspect("equal")
        a = ax[row, 2]; x0, x1, y0, y1 = -1.5, 2.5, -1.6, 1.6
        cent = Q.mean(1); m = (cent[:, 0] > x0 - 1) & (cent[:, 0] < x1 + 1) & (np.abs(cent[:, 1]) < y1 + 1)
        pc = PolyCollection(Q[m], array=no[m], cmap="magma", clim=(0, 75), edgecolors="none")
        a.add_collection(pc); a.plot([0, 1], [0, 0], "c-", lw=1.2)
        a.set(xlim=(x0, x1), ylim=(y0, y1), title=f"{name.split(',')[0]}: non-orthogonality (max {no.max():.1f}°)"); a.set_aspect("equal")
        fig.colorbar(pc, ax=a, label="cell max face non-orthogonality (deg)", shrink=0.8)
    fig.tight_layout(); fig.savefig(HERE / "figs/mesh_comparison.png", dpi=150)


if __name__ == "__main__":
    main()
