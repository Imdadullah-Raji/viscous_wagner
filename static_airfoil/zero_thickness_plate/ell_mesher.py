"""Elliptic-coordinate (confocal) mesh around a zero-thickness flat plate, for OpenFOAM via
gmshToFoam + createBaffles.  Replaces the H-type ztp_mesher.py, whose side blocks reach 73-76 deg
non-orthogonality next to the tips (see zeroThickPlate.md).

GRID.  The conformal map  z = c/2 + (c/2) cosh(w),  w = mu + i nu:
    x = c/2 + (c/2) cosh(mu) cos(nu),   y = (c/2) sinh(mu) sin(nu)
  mu = const : confocal ellipses with foci at the tips (0,0), (c,0); mu = 0 is the plate itself,
               nu in (0, pi) the upper side, (pi, 2 pi) the lower side;
  nu = const : confocal hyperbolas, radial far away.
It is the Joukowski image of a polar grid around a circle (equipotentials / field lines of a charged
strip).  Conformal => orthogonal everywhere except where dz/dw = (c/2) sinh(w) = 0, i.e. at the two
tips, where angles double: the cells touching a tip have one 180-deg corner.  The scale factor is
(c/2) sqrt(sinh^2 mu + sin^2 nu), the same in both directions, so the grid clusters at the tips by
itself: tip cell ~ (c/4) dnu^2 along the plate, ~ (c/4) dmu0^2 across.

SPACING
  nu : uniform, n_nu cells round (even) -> cosine spacing along the plate;
  mu : geometric from mu = 0 to mu_max = arccosh(2R/c) (outer ellipse semi-major R, centred at
       mid-chord), first step chosen so the first cell normal to the plate AT MID-CHORD is
       (c/2) sinh(dmu0) = h1 = h1_coef / sqrt(Re).

OUTPUT: .msh v2.2, one hex layer of depth `span`, physical groups `internal` (3D), `farfield`,
`frontAndBack`, and `plate` -- the plate faces are INTERNAL (upper and lower plate nodes are
shared), so gmshToFoam writes them as faceZone `plate` and createBaffles (createBafflesDict here)
splits it into walls plate_upper / plate_lower (group `airfoil`).

    python3 ell_mesher.py out.msh --re 500 [--nnu 560 --nmu 210 -r 30]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "meshing"))
from mesh_math import geometric_ratio  # noqa: E402


def grid(c, R, n_nu, n_mu, h1):
    mu_max = np.arccosh(2 * R / c)
    dmu0 = np.arcsinh(2 * h1 / c)
    r = geometric_ratio(mu_max, dmu0, n_mu)
    mu = np.concatenate([[0.0], np.cumsum(dmu0 * r ** np.arange(n_mu))])
    mu[-1] = mu_max
    nu = 2 * np.pi * np.arange(n_nu) / n_nu
    MU, NU = np.meshgrid(mu, nu, indexing="ij")          # (n_mu+1, n_nu), nu periodic
    X = c / 2 + c / 2 * np.cosh(MU) * np.cos(NU)
    Y = c / 2 * np.sinh(MU) * np.sin(NU)
    Y[0] = 0.0                                           # mu = 0 is exactly the plate
    return X, Y, dict(mu_max=mu_max, dmu0=dmu0, r_mu=r)


def write_msh(path, X, Y, span):
    n_mu1, n_nu = X.shape
    half = n_nu // 2
    # 2D node numbering: ring 0 (the plate) shares nodes between nu_j and nu_(n_nu - j)
    nid2 = np.zeros((n_mu1, n_nu), int)
    k = 0
    for j in range(half + 1):
        nid2[0, j] = k; k += 1
    for j in range(half + 1, n_nu):
        nid2[0, j] = nid2[0, n_nu - j]
    for i in range(1, n_mu1):
        nid2[i] = k + np.arange(n_nu); k += n_nu
    n2 = k
    xy = np.zeros((n2, 2))
    xy[nid2.ravel()] = np.stack([X.ravel(), Y.ravel()], 1)
    # merged plate nodes must coincide
    assert np.allclose(X[0, 1:half], X[0, n_nu - 1:half:-1])

    jj = np.arange(n_nu); jp = (jj + 1) % n_nu
    quads = []   # (a, b, c, d) counter-clockwise: (i,j) (i+1,j) (i+1,j+1) (i,j+1)
    for i in range(n_mu1 - 1):
        quads.append(np.stack([nid2[i, jj], nid2[i + 1, jj], nid2[i + 1, jp], nid2[i, jp]], 1))
    quads = np.concatenate(quads)
    plate = np.stack([nid2[0, np.arange(half)], nid2[0, np.arange(1, half + 1)]], 1)     # upper side only
    far = np.stack([nid2[-1, jj], nid2[-1, jp]], 1)

    b, t = 1, 1 + n2           # gmsh node ids: bottom layer, top layer
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat",
             "$PhysicalNames", "4", '3 1 "internal"', '2 2 "farfield"', '2 3 "frontAndBack"', '2 4 "plate"',
             "$EndPhysicalNames", "$Nodes", str(2 * n2)]
    lines += [f"{b + n} {x:.15g} {y:.15g} 0" for n, (x, y) in enumerate(xy)]
    lines += [f"{t + n} {x:.15g} {y:.15g} {span:.15g}" for n, (x, y) in enumerate(xy)]
    lines.append("$EndNodes")
    el = []
    for q in quads:      # hex: bottom quad then top quad, same order
        el.append((5, 1, [*(b + q), *(t + q)]))
    for e in far:
        el.append((3, 2, [b + e[0], b + e[1], t + e[1], t + e[0]]))
    for q in quads:
        el.append((3, 3, list(b + q[::-1])))
        el.append((3, 3, list(t + q)))
    for e in plate:
        el.append((3, 4, [b + e[0], b + e[1], t + e[1], t + e[0]]))
    lines += ["$Elements", str(len(el))]
    lines += [f"{n + 1} {typ} 2 {phys} {phys} " + " ".join(map(str, nodes)) for n, (typ, phys, nodes) in enumerate(el)]
    lines.append("$EndElements")
    Path(path).write_text("\n".join(lines) + "\n")
    return len(quads)


def main():
    ap = argparse.ArgumentParser(description="Elliptic-coordinate mesh around a zero-thickness plate")
    ap.add_argument("filename", help="output .msh (v2.2, for gmshToFoam)")
    ap.add_argument("--re", type=float, default=500, help="Reynolds number (sets h1)")
    ap.add_argument("--h1-coef", type=float, default=0.05, help="h1 = coef/sqrt(Re) at mid-chord")
    ap.add_argument("-r", "--radius", type=float, default=30.0, help="outer ellipse semi-major axis (centred at mid-chord)")
    ap.add_argument("--nnu", type=int, default=560, help="cells round the plate (even); sets the tip cell")
    ap.add_argument("--nmu", type=int, default=210, help="cells from the plate to the far field")
    ap.add_argument("--span", type=float, default=0.1, help="extrusion depth (Aref = c*span)")
    args = ap.parse_args()
    if args.nnu % 2:
        ap.error("--nnu must be even")
    c = 1.0
    h1 = args.h1_coef / np.sqrt(args.re)
    X, Y, info = grid(c, args.radius, args.nnu, args.nmu, h1)
    n = write_msh(args.filename, X, Y, args.span)
    xp = np.sort(X[0, : args.nnu // 2 + 1])
    print(f"  h1 (mid-chord)   : {h1:.4e}")
    print(f"  mu_max, dmu0     : {info['mu_max']:.5f}, {info['dmu0']:.4e}   ratio {info['r_mu']:.5f}")
    print(f"  tip cell (plate) : {xp[1] - xp[0]:.3e}   mid-chord cell {np.diff(xp).max():.3e}")
    print(f"  first wake cell  : {X[1, 0] - c:.3e}")
    print(f"  hexes            : {n}  (n_nu * n_mu = {args.nnu * args.nmu})")
    assert n == args.nnu * args.nmu


if __name__ == "__main__":
    main()
