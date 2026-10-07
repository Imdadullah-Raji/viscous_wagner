"""Build the structured O-ring around a truncated-TE airfoil using the gmsh
Python API (built-in `geo` kernel).

Scope: the ring ONLY -- airfoil boundary out to the circular interface.  No
extrusion, no physical groups, no outer domain.  `ORing.outer_curves` hands back
the far-field arc tags in counter-clockwise order so a rectangular outer domain
can be attached to them later without rebuilding anything.

Three blocks, all traversed counter-clockwise:

    upper : LE -> (upper spline) -> TEu -> (TEu cut) -> ff(+tw)
                 -> (upper arcs) -> ff(pi) -> (LE cut, reversed) -> LE
    lower : LE -> (LE cut) -> ff(pi) -> (lower arcs) -> ff(-tw)
                 -> (TEl cut, reversed) -> TEl -> (lower spline, reversed) -> LE
    wake  : TEu -> (TE face) -> TEl -> (TEl cut) -> ff(-tw)
                 -> (wake arc) -> ff(+tw) -> (TEu cut, reversed) -> TEu
"""

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np
import gmsh

from airfoil import Airfoil
from mesh_math import FarField, cut_line_ratio, geometric_ratio


@dataclass
class ORing:
    surfaces: Dict[str, int]          # 'upper' | 'lower' | 'wake' -> surface tag
    wall_curves: List[int]            # airfoil boundary, for a later Physical group
    outer_curves: List[int]           # far-field arcs, CCW, for the outer domain
    ratios: Dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _add_points(coords: np.ndarray) -> List[int]:
    return [gmsh.model.geo.addPoint(x, y, 0.0) for x, y in coords]


class _CirclePoints:
    """Memoised far-field points, keyed on angle mod 2*pi.

    Needed because `FarField.node_angles` starts at -tw and ends at 2*pi-tw --
    the same physical point.  Creating it twice leaves the ring topologically
    open and the last arc silently fails to close the loop.
    """

    def __init__(self, ff: FarField):
        self.ff = ff
        self._cache: Dict[int, int] = {}

    def __call__(self, theta: float) -> int:
        key = int(round((theta % (2 * np.pi)) * 1e9))
        if key not in self._cache:
            x, y = self.ff.point(theta)
            self._cache[key] = gmsh.model.geo.addPoint(x, y, 0.0)
        return self._cache[key]


def _transfinite(tag: int, n_cells: int, kind: str = "Progression", coef: float = 1.0):
    """gmsh counts NODES, we think in CELLS.  This is the off-by-one that
    silently gives you one cell too few on every single curve."""
    gmsh.model.geo.mesh.setTransfiniteCurve(tag, n_cells + 1, kind, coef)


# ---------------------------------------------------------------------------
# main build
# ---------------------------------------------------------------------------

def build_oring(
    af: Airfoil,
    ff: FarField,
    ns: int,
    nte: int,
    nr: int,
    h1: float,
    bump: float = 0.08,
) -> ORing:
    """Construct the three-block O-ring.  Assumes gmsh.initialize() was called."""
    geo = gmsh.model.geo

    # --- points -----------------------------------------------------------
    # The LE row is shared between both surfaces (airfoil.py contract), so it is
    # created once and reused; otherwise the two splines end at coincident but
    # distinct points and the curve loops never close.
    p_le = geo.addPoint(*af.le, 0.0)
    up_tags = [p_le] + _add_points(af.upper[1:])
    lo_tags = [p_le] + _add_points(af.lower[1:])
    p_teu, p_tel = up_tags[-1], lo_tags[-1]

    p_center = geo.addPoint(*ff.center, 0.0)
    cp = _CirclePoints(ff)

    # --- airfoil curves ---------------------------------------------------
    c_upper = geo.addSpline(up_tags)          # LE -> TEu
    c_lower = geo.addSpline(lo_tags)          # LE -> TEl
    c_te = geo.addLine(p_teu, p_tel)          # TEu -> TEl

    # --- far-field arcs ---------------------------------------------------
    arcs: Dict[str, List[int]] = {"wake": [], "upper": [], "lower": []}
    outer_curves: List[int] = []
    for a in ff.arcs:
        tag = geo.addCircleArc(cp(a.theta_start), p_center, cp(a.theta_end))
        arcs[a.role].append(tag)
        outer_curves.append(tag)
        _transfinite(tag, a.n_cells)

    # --- cut lines (airfoil -> far field, so Progression grows AWAY from the
    #     wall; reverse the direction and your fine cells end up outside) -----
    tw = ff.theta_w
    r_le, _ = cut_line_ratio(af.le, ff, np.pi, h1, nr)
    r_teu, _ = cut_line_ratio(af.te_upper, ff, tw, h1, nr)
    r_tel, _ = cut_line_ratio(af.te_lower, ff, -tw, h1, nr)

    print("\nmean progession ratio:", (r_le+r_teu)/2 )

    c_cut_le = geo.addLine(p_le, cp(np.pi))
    c_cut_teu = geo.addLine(p_teu, cp(tw))
    c_cut_tel = geo.addLine(p_tel, cp(-tw))

    _transfinite(c_cut_le, nr, "Progression", r_le)
    _transfinite(c_cut_teu, nr, "Progression", r_teu)
    _transfinite(c_cut_tel, nr, "Progression", r_tel)

    # --- surface discretisation ------------------------------------------
    # Bump clusters toward BOTH ends of the curve, which is what the LE and the
    # truncated TE both want.  The opposite (far-field) side stays uniform;
    # transfinite interpolation blends the two distributions radially.
    _transfinite(c_upper, ns, "Bump", bump)
    _transfinite(c_lower, ns, "Bump", bump)
    _transfinite(c_te, nte)

    # --- blocks -----------------------------------------------------------
    loops = {
        "upper": [c_upper, c_cut_teu] + arcs["upper"] + [-c_cut_le],
        "lower": [c_cut_le] + arcs["lower"] + [-c_cut_tel, -c_lower],
        "wake": [c_te, c_cut_tel] + arcs["wake"] + [-c_cut_teu],
    }
    corners = {
        "upper": [p_le, p_teu, cp(tw), cp(np.pi)],
        "lower": [p_le, cp(np.pi), cp(-tw), p_tel],
        "wake": [p_teu, p_tel, cp(-tw), cp(tw)],
    }

    surfaces: Dict[str, int] = {}
    for name, curves in loops.items():
        s = geo.addPlaneSurface([geo.addCurveLoop(curves)])
        geo.mesh.setTransfiniteSurface(s, "Left", corners[name])
        geo.mesh.setRecombine(2, s)
        surfaces[name] = s

    geo.synchronize()
    return ORing(
        surfaces=surfaces,
        wall_curves=[c_upper, c_lower, c_te],
        outer_curves=outer_curves,
        ratios={"le": r_le, "teu": r_teu, "tel": r_tel},
    )


def report(ring: ORing, ns: int, nte: int, nr: int) -> None:
    """Post-generation sanity numbers.  Quads only, and the count must match."""
    expected = (2 * ns + nte) * nr
    quads = len(gmsh.model.mesh.getElementsByType(3)[0])
    tris = len(gmsh.model.mesh.getElementsByType(2)[0])
    print(f"  growth ratios : " + "  ".join(f"{k}={v:.6f}" for k, v in ring.ratios.items()))
    print(f"  quads         : {quads}  (expected {expected})")
    print(f"  triangles     : {tris}  (must be 0)")
    assert tris == 0, "triangles present -- a block is not transfinite"
    assert quads == expected, "cell count mismatch -- check the node/cell off-by-one"


if __name__ == "__main__":
    from airfoil import naca4
    from mesh_math import build_farfield
    import argparse

    parser = argparse.ArgumentParser(
        description="Takes Re and Filename"
        )

    parser.add_argument(
        "filename", 
        type = str, 
        help="OUTFILE name" 
        )
    parser.add_argument(
        "--re", 
        type=float, 
        default=100, 
        help="Reynolds number"
        )
    parser.add_argument(
        "-r",
        "--radius", 
        type=float, 
        default=20,
        help = "Oring Radius")
    parser.add_argument(
        "--nr",
        type=int,
        default=200,
        help="radial cells on each cut line")
    parser.add_argument(
        "--ns",
        type=int,
        default=200,
        help="cells along each airfoil surface")
    parser.add_argument(
        "--nte",
        type=int,
        default=40,
        help="cells across the truncated trailing edge")
    parser.add_argument(
        "--bump",
        type=float,
        default=0.08,
        help="gmsh Bump coefficient: LE/TE cell ~ bump * uniform spacing. "
             "This sets the Courant limit of the whole simulation.")
    parser.add_argument(
        "--naca",
        type=str,
        default="0012",
        help="NACA 4-digit code (default 0012; e.g. 0006 for the thickness study)")
    parser.add_argument(
        "--h1-coef",
        type=float,
        default=0.05,
        help="first cell h1 = coef/sqrt(Re); 0.05 is the production rule")
    args = parser.parse_args()

    Re = args.re
    NS, NTE, NR, R = args.ns, args.nte, args.nr, args.radius
    H1 = args.h1_coef/np.sqrt(Re)
    OUTFILE = args.filename


    af = naca4(args.naca, n_points=200, x_trunc=0.98)
    ff = build_farfield(center=(0.25, 0.0), radius=R, ns=NS, nte=NTE)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("oring")

    ring = build_oring(af, ff, ns=NS, nte=NTE, nr=NR, h1=H1, bump=args.bump)


    THICK = 0.1
    geo = gmsh.model.geo

    # All three blocks in ONE extrude call — separate calls duplicate the lateral
    # faces on the shared cut lines and the volumes come out non-conformal.
    out = geo.extrude([(2, s) for s in ring.surfaces.values()],
                    0, 0, THICK, numElements=[1], recombine=True)
    geo.synchronize()

    vols = [t for d, t in out if d == 3]

    eps = THICK * 1e-3
    box = lambda z0, z1: [t for d, t in gmsh.model.getEntitiesInBoundingBox(
        -1e9, -1e9, z0, 1e9, 1e9, z1, 2)]
    front, back = box(-eps, eps), box(THICK - eps, THICK + eps)
    lateral = [t for d, t in gmsh.model.getEntities(2)
            if t not in front and t not in back]

    def touching(curves):
        want = {abs(c) for c in curves}
        return [s for s in lateral
                if want & {abs(c) for _, c in
                        gmsh.model.getBoundary([(2, s)], oriented=False)}]

    gmsh.model.addPhysicalGroup(3, vols, name="internal")
    gmsh.model.addPhysicalGroup(2, touching(ring.wall_curves),  name="airfoil")
    gmsh.model.addPhysicalGroup(2, touching(ring.outer_curves), name="farfield")
    gmsh.model.addPhysicalGroup(2, front + back,                name="frontAndBack")

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(OUTFILE)

    gmsh.finalize()

    # number of layers within the boundary 

    print("Max progression ratio to get 25 layer witin the BL",
          np.round(geometric_ratio(100*H1, H1, 25),3))



