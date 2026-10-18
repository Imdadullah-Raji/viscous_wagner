"""Structured all-hex mesh around a ZERO-thickness flat plate, for OpenFOAM via
gmshToFoam + createBaffles.  Python-API port of the user's `zeroThickPlate.geo`
(kept unchanged as the original), in the style of `plate_case/meshing/plate_mesher.py`.

TOPOLOGY (as in the .geo).  The plate is an internal curve shared by an upper and
a lower rectangular block (height b); six outer blocks fill the rest out to a
circular far field.  Six spokes: one from each tip along y = 0 (upstream cut and
wake cut), and one from each block corner to the far field.

    far field arcs:  up | right-up | right-down | down | left-down | left-up
    blocks        :  upper, lower (on the plate); upper/lower outer (above/below
                     them); four side blocks (between a vertical block edge, a
                     spoke along y = 0 and a corner spoke)

The plate becomes a faceZone named `plate` (gmshToFoam turns tagged INTERNAL faces
into a faceZone of the same name), which `createBaffles` splits into two wall
patches.  See `createBafflesDict` next to this file.

DIFFERENCES FROM THE .geo
  - chord on [0, 1] (LE at the origin, as every other mesh here); far field
    centred at (center, 0) with radius R = 30 by default (the campaign value);
  - first cell normal to the plate h1 = h1_coef / sqrt(Re) (the production rule),
    reached by solving the block-edge progression for it;
  - the tip cell is an explicit length (--tip), not a bump coefficient: each plate
    edge is split at mid-chord and each half is a progression from its tip.  There
    is no geometric length at a zero-radius edge, so --tip is a convergence
    parameter (see zeroThickPlate.md);
  - the first spoke cell matches the cell it continues (tip cell along y = 0, the
    last block-edge cell at the block corners), so there is no size jump;
  - extruded one cell in z (span 0.1) with internal / farfield / frontAndBack /
    plate physical groups, written as .msh v2.2 for gmshToFoam.

Run with the research env:
    python ztp_mesher.py out.msh --re 500 [--tip 2e-3] [--nplate 160] ...
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "meshing"))
from mesh_math import geometric_length, geometric_ratio  # noqa: E402


def _ratio(length, h0, n):
    """geometric_ratio, but also allows shrinking (length < n*h0)."""
    if length >= n * h0:
        return geometric_ratio(length, h0, n)
    lo, hi = 0.5, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if geometric_length(h0, mid, n) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _transfinite(tag, n_cells, ratio=1.0):
    # gmsh counts NODES; we count CELLS.  A negative tag reverses the progression.
    gmsh.model.geo.mesh.setTransfiniteCurve(abs(tag), n_cells + 1, "Progression",
                                            ratio if tag > 0 else 1.0 / ratio)


def build(c, b, R, xc, nplate, nside, nspokes, h1, tip, ff_blend):
    """Return (block surfaces, plate curves, far-field arcs, info dict)."""
    geo = gmsh.model.geo
    P = lambda x, y: geo.addPoint(x, y, 0.0)
    L, M, T = P(0, 0), P(c / 2, 0), P(c, 0)
    LU, MU, RU = P(0, b), P(c / 2, b), P(c, b)
    LD, MD, RD = P(0, -b), P(c / 2, -b), P(c, -b)
    O = P(xc, 0)
    F = {a: P(xc + R * np.cos(np.radians(a)), R * np.sin(np.radians(a)))
         for a in (180, 135, 90, 45, 0, -45, -90, -135)}
    Fm_up, Fm_dn = F[90], F[-90]

    line = geo.addLine
    arc = lambda a, bb: geo.addCircleArc(a, O, bb)
    # plate halves, block top/bottom halves, far arcs over/under the blocks (split at 90 deg)
    pl1, pl2 = line(L, M), line(M, T)
    up1, up2 = line(LU, MU), line(MU, RU)
    dn1, dn2 = line(LD, MD), line(MD, RD)
    au1, au2 = arc(F[135], Fm_up), arc(Fm_up, F[45])
    ad1, ad2 = arc(F[-135], Fm_dn), arc(Fm_dn, F[-45])
    # vertical block edges (from the plate outward) and the mid-chord verticals
    vLu, vMu, vRu = line(L, LU), line(M, MU), line(T, RU)
    vLd, vMd, vRd = line(L, LD), line(M, MD), line(T, RD)
    # side far arcs (from the y = 0 line outward)
    aLu, aLd = arc(F[180], F[135]), arc(F[180], F[-135])
    aRu, aRd = arc(F[0], F[45]), arc(F[0], F[-45])
    # spokes (from the inner point outward); mid-chord spokes split the outer blocks
    sL, sR = line(L, F[180]), line(T, F[0])
    sLU, sRU, sLD, sRD = line(LU, F[135]), line(RU, F[45]), line(LD, F[-135]), line(RD, F[-45])
    sMU, sMD = line(MU, Fm_up), line(MD, Fm_dn)

    # ---- distributions
    half = c / 2
    n_h = nplate // 2
    r_pl = _ratio(half, tip, n_h)                 # plate half, tip -> mid-chord
    r_side = geometric_ratio(b, h1, nside)        # block edge, h1 at the plate
    h_top = h1 * r_side ** (nside - 1)            # last block-edge cell
    len_L = abs((xc - R) - 0.0)                   # spokes along y = 0, from each tip
    len_R = abs((xc + R) - c)
    r_sL = _ratio(len_L, tip, nspokes)
    r_sR = _ratio(len_R, tip, nspokes)
    corner_len = lambda x0, y0, a: float(np.hypot(xc + R * np.cos(np.radians(a)) - x0,
                                                  R * np.sin(np.radians(a)) - y0))
    r_sLU = _ratio(corner_len(0, b, 135), h_top, nspokes)
    r_sRU = _ratio(corner_len(c, b, 45), h_top, nspokes)
    r_sMU = _ratio(corner_len(c / 2, b, 90), h_top, nspokes)

    # plate, block top/bottom: progression away from each tip (halves meet at mid-chord)
    for t in (pl1, up1, dn1):
        _transfinite(t, n_h, r_pl)
    for t in (pl2, up2, dn2):
        _transfinite(-t, n_h, r_pl)
    # far arcs over/under the blocks: blend the tip clustering with uniform so the far field is not
    # starved mid-arc (ff_blend = 0 -> same distribution as the plate, 1 -> uniform)
    r_ff = r_pl ** (1.0 - ff_blend)
    for t in (au1, ad1):
        _transfinite(t, n_h, r_ff)
    for t in (au2, ad2):
        _transfinite(-t, n_h, r_ff)
    # verticals and side far arcs: clustered at y = 0 (h1 at the plate)
    for t in (vLu, vMu, vRu, vLd, vMd, vRd):
        _transfinite(t, nside, r_side)
    r_side_ff = r_side ** (1.0 - ff_blend)
    for t in (aLu, aLd, aRu, aRd):
        _transfinite(t, nside, r_side_ff)
    # spokes
    _transfinite(sL, nspokes, r_sL); _transfinite(sR, nspokes, r_sR)
    for t, r in ((sLU, r_sLU), (sLD, r_sLU), (sRU, r_sRU), (sRD, r_sRU), (sMU, r_sMU), (sMD, r_sMU)):
        _transfinite(t, nspokes, r)

    # ---- surfaces (curve loops; orientation handled by signs)
    def surf(curves):
        cl = geo.addCurveLoop(curves)
        s = geo.addPlaneSurface([cl])
        geo.mesh.setTransfiniteSurface(s)
        geo.mesh.setRecombine(2, s)
        return s

    S = [
        surf([pl1, vMu, -up1, -vLu]), surf([pl2, vRu, -up2, -vMu]),          # upper block (2 halves)
        surf([dn1, -vMd, -pl1, vLd]), surf([dn2, -vRd, -pl2, vMd]),          # lower block
        surf([up1, sMU, -au1, -sLU]), surf([up2, sRU, -au2, -sMU]),          # upper outer
        surf([-dn1, sLD, ad1, -sMD]), surf([-dn2, sMD, ad2, -sRD]),          # lower outer
        surf([sL, aLu, -sLU, -vLu]), surf([-sL, vLd, sLD, -aLd]),            # left side blocks
        surf([-sR, vRu, sRU, -aRu]), surf([sR, aRd, -sRD, -vRd]),            # right side blocks
    ]
    info = dict(h1=h1, tip=tip, r_plate=r_pl, r_side=r_side, h_top=h_top,
                r_spoke_axis=max(r_sL, r_sR), r_spoke_corner=max(r_sLU, r_sRU, r_sMU))
    return S, [pl1, pl2], [au1, au2, ad1, ad2, aLu, aLd, aRu, aRd], info


def main():
    ap = argparse.ArgumentParser(description="Structured all-hex mesh around a zero-thickness plate")
    ap.add_argument("filename", help="output .msh (v2.2, for gmshToFoam)")
    ap.add_argument("--re", type=float, default=500, help="Reynolds number (sets h1)")
    ap.add_argument("--h1-coef", type=float, default=0.05, help="h1 = coef/sqrt(Re); 0.05 is the production rule")
    ap.add_argument("--tip", type=float, default=None,
                    help="tip cell along the plate (default: h1). A convergence parameter, see zeroThickPlate.md")
    ap.add_argument("-r", "--radius", type=float, default=30.0, help="far-field radius")
    ap.add_argument("--center", type=float, default=0.25, help="far-field centre x")
    ap.add_argument("--block-h", type=float, default=0.5, help="height b of the blocks on the plate")
    ap.add_argument("--nplate", type=int, default=160, help="cells along the plate (even)")
    ap.add_argument("--nside", type=int, default=60, help="cells across each block (normal to the plate)")
    ap.add_argument("--nspokes", type=int, default=120, help="cells from the blocks to the far field")
    ap.add_argument("--ff-blend", type=float, default=0.5,
                    help="far-field clustering: 0 = copy the inner distribution, 1 = uniform")
    ap.add_argument("--span", type=float, default=0.1, help="extrusion depth (Aref = c*span)")
    args = ap.parse_args()
    if args.nplate % 2:
        ap.error("--nplate must be even")

    c = 1.0
    h1 = args.h1_coef / np.sqrt(args.re)
    tip = args.tip if args.tip else h1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("zero_thickness_plate")
    surfaces, plate, outer, info = build(c, args.block_h, args.radius, args.center, args.nplate,
                                         args.nside, args.nspokes, h1, tip, args.ff_blend)
    geo = gmsh.model.geo
    # ONE extrude call so blocks share their interface faces (and the plate face is shared).
    out = geo.extrude([(2, s) for s in surfaces], 0, 0, args.span, numElements=[1], recombine=True)
    geo.synchronize()
    vols = [t for d, t in out if d == 3]

    eps = args.span * 1e-3
    box = lambda z0, z1: [t for d, t in gmsh.model.getEntitiesInBoundingBox(-1e9, -1e9, z0, 1e9, 1e9, z1, 2)]
    front, back = box(-eps, eps), box(args.span - eps, args.span + eps)
    lateral = [t for d, t in gmsh.model.getEntities(2) if t not in front and t not in back]

    def touching(curves):
        want = {abs(x) for x in curves}
        return [s for s in lateral
                if want & {abs(x) for _, x in gmsh.model.getBoundary([(2, s)], oriented=False)}]

    plate_s = touching(plate)
    assert len(plate_s) == 2, plate_s
    gmsh.model.addPhysicalGroup(3, vols, name="internal")
    gmsh.model.addPhysicalGroup(2, touching(outer), name="farfield")
    gmsh.model.addPhysicalGroup(2, front + back, name="frontAndBack")
    gmsh.model.addPhysicalGroup(2, plate_s, name="plate")      # internal -> faceZone in gmshToFoam

    gmsh.model.mesh.generate(3)
    hexes = len(gmsh.model.mesh.getElementsByType(5)[0])
    others = sum(len(gmsh.model.mesh.getElementsByType(t)[0]) for t in (4, 6, 7))
    n_h = args.nplate // 2
    expected = (4 * n_h * args.nside            # upper + lower block
                + 4 * n_h * args.nspokes        # upper + lower outer
                + 4 * args.nside * args.nspokes)  # four side blocks
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # gmshToFoam won't read 4.1
    gmsh.write(args.filename)
    gmsh.finalize()

    for k, v in info.items():
        print(f"  {k:<15s}: {v:.6g}")
    print(f"  hexes          : {hexes}  (expected {expected})")
    print(f"  non-hex cells  : {others}  (must be 0)")
    assert others == 0 and hexes == expected


if __name__ == "__main__":
    main()
