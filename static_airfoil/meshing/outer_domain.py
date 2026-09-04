"""Rectangular outer domain wrapped around the structured O-ring.

Separate module by design: it consumes an `ORing` handle and never touches the
ring's construction.  The attachment is conforming by reuse -- the same five arc
tags that bound the ring are passed straight into the outer surface's hole loop,
so the interface nodes are shared, not interpolated.

Sizing is driven entirely by background fields.  The three MeshSizeFrom* options
are switched off on purpose: left on, gmsh extends the 0.076 interface spacing
outward and quietly overrides every field you set, which is the usual reason
wake refinement appears to do nothing.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import gmsh

from mesh_math import FarField
from mesher import ORing


@dataclass
class OuterDomain:
    surface: int
    boundary: Dict[str, int]     # 'inlet' | 'outlet' | 'top' | 'bottom' -> curve tag
    interface: List[int]         # the shared arcs, echoed for convenience
    n_boundary_cells: int


@dataclass
class WakeBox:
    """Refinement window.  `thickness` is the graded transition band outside it."""
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    size: float
    thickness: float


def interface_spacing(ff: FarField) -> float:
    """Cell size the ring hands to the outer region -- the number every field
    below has to be compatible with."""
    a = ff.arcs[0]
    return ff.radius * abs(a.sweep) / a.n_cells


def check_quad_parity(ff: FarField) -> int:
    """An all-quad mesh needs an EVEN number of boundary edges (4Q = 2*interior
    + boundary).  The rectangle contributes 2*nx + 2*ny, which is even for any
    nx, ny -- so it can never be the culprit.  The ring contributes 2*ns + nte,
    even only if nte is even.  Parity is necessary but NOT sufficient: blossom
    is a matching heuristic and will still leave a handful of triangles.
    """
    n_if = sum(a.n_cells for a in ff.arcs)
    if n_if % 2:
        raise ValueError(
            f"interface has {n_if} edges (odd) -- an all-quad outer region is "
            f"impossible; use an even nte"
        )
    return n_if


def build_outer_domain(
    ring: ORing,
    ff: FarField,
    bounds: Tuple[float, float, float, float] = (-15.0, 40.0, -20.0, 20.0),
    h_far: float = 0.75,
    wake_boxes: List[WakeBox] = None,
    grow_dist: float = 6.0,
) -> OuterDomain:
    geo = gmsh.model.geo
    xmin, xmax, ymin, ymax = bounds
    check_quad_parity(ff)

    # --- rectangle --------------------------------------------------------
    p = [
        geo.addPoint(xmin, ymin, 0.0),
        geo.addPoint(xmax, ymin, 0.0),
        geo.addPoint(xmax, ymax, 0.0),
        geo.addPoint(xmin, ymax, 0.0),
    ]
    c_bottom = geo.addLine(p[0], p[1])
    c_outlet = geo.addLine(p[1], p[2])
    c_top = geo.addLine(p[2], p[3])
    c_inlet = geo.addLine(p[3], p[0])

    nx = max(2, int(round((xmax - xmin) / h_far)))
    ny = max(2, int(round((ymax - ymin) / h_far)))
    for tag, n in ((c_bottom, nx), (c_top, nx), (c_outlet, ny), (c_inlet, ny)):
        geo.mesh.setTransfiniteCurve(tag, n + 1)
    n_bnd = 2 * nx + 2 * ny

    # --- surface with the ring as a hole ---------------------------------
    # Outer loop CCW, hole loop CW: reverse the arc order and flip the signs.
    outer_loop = geo.addCurveLoop([c_bottom, c_outlet, c_top, c_inlet])
    hole_loop = geo.addCurveLoop([-t for t in reversed(ring.outer_curves)])
    s = geo.addPlaneSurface([outer_loop, hole_loop])
    geo.mesh.setRecombine(2, s)
    geo.synchronize()

    # --- background sizing fields ----------------------------------------
    fld = gmsh.model.mesh.field
    ids: List[int] = []

    # grow away from the interface, starting exactly at the ring's own spacing
    f_dist = fld.add("Distance")
    fld.setNumbers(f_dist, "CurvesList", ring.outer_curves)
    f_thr = fld.add("Threshold")
    fld.setNumber(f_thr, "InField", f_dist)
    fld.setNumber(f_thr, "SizeMin", interface_spacing(ff))
    fld.setNumber(f_thr, "SizeMax", h_far)
    fld.setNumber(f_thr, "DistMin", 0.0)
    fld.setNumber(f_thr, "DistMax", grow_dist)
    ids.append(f_thr)

    for w in wake_boxes or []:
        f = fld.add("Box")
        fld.setNumber(f, "XMin", w.xmin)
        fld.setNumber(f, "XMax", w.xmax)
        fld.setNumber(f, "YMin", w.ymin)
        fld.setNumber(f, "YMax", w.ymax)
        fld.setNumber(f, "VIn", w.size)
        fld.setNumber(f, "VOut", h_far)
        fld.setNumber(f, "Thickness", w.thickness)
        ids.append(f)

    f_min = fld.add("Min")
    fld.setNumbers(f_min, "FieldsList", ids)
    fld.setAsBackgroundMesh(f_min)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)            # frontal-Delaunay for quads
    # Blossom only.  Algorithms 2/3 ("full-quad") subdivide the 1D mesh, which
    # would halve every transfinite spacing in the ring -- gmsh refuses outright.
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)

    return OuterDomain(
        surface=s,
        boundary={"inlet": c_inlet, "outlet": c_outlet,
                  "top": c_top, "bottom": c_bottom},
        interface=list(ring.outer_curves),
        n_boundary_cells=n_bnd,
    )


DEFAULT_WAKE = [
    WakeBox(2.5, 30.0, -1.5, 1.5, size=0.10, thickness=2.5),
    WakeBox(2.5, 40.0, -4.0, 4.0, size=0.40, thickness=4.0),
]


if __name__ == "__main__":
    from airfoil import naca4
    from mesh_math import build_farfield
    from mesher import build_oring

    NS, NTE, NR, H1, R = 120, 8, 70, 7.0e-4, 3.0

    af = naca4("0012", n_points=200, x_trunc=0.98)
    ff = build_farfield(center=(0.25, 0.0), radius=R, ns=NS, nte=NTE)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("airfoil_ogrid")

    ring = build_oring(af, ff, ns=NS, nte=NTE, nr=NR, h1=H1)
    outer = build_outer_domain(ring, ff, wake_boxes=DEFAULT_WAKE)

    gmsh.model.mesh.generate(2)

    quads = len(gmsh.model.mesh.getElementsByType(3)[0])
    tris = len(gmsh.model.mesh.getElementsByType(2)[0])
    ring_cells = (2 * NS + NTE) * NR
    print(f"  interface ds  : {interface_spacing(ff):.4f}")
    print(f"  ring cells    : {ring_cells}")
    print(f"  total quads   : {quads}   (outer = {quads - ring_cells})")
    print(f"  triangles     : {tris}")
    print(f"  bnd parity    : ring {2 * NS + NTE} + rect {outer.n_boundary_cells} "
          f"= {2 * NS + NTE + outer.n_boundary_cells}")

    gmsh.option.setNumber("Mesh.SaveAll", 1)
    # gmsh.write("airfoil_ogrid.msh")
    # gmsh.finalize()

    THICK = 1.0
    geo = gmsh.model.geo

    # ONE call with all four surfaces. Extruding them separately creates duplicate
    # lateral faces on the shared arcs and the volumes come out non-conformal.
    planar = [(2, s) for s in ring.surfaces.values()] + [(2, outer.surface)]
    out = geo.extrude(planar, 0, 0, THICK, numElements=[1], recombine=True)
    geo.synchronize()

    vols = [t for d, t in out if d == 3]

    # getEntitiesInBoundingBox returns only entities FULLY inside the box, so the
    # z=0 and z=THICK planes pick up front/back and never the lateral faces.
    eps = THICK * 1e-3
    box = lambda z0, z1: [t for d, t in gmsh.model.getEntitiesInBoundingBox(
        -1e6, -1e6, z0, 1e6, 1e6, z1, 2)]
    front, back = box(-eps, eps), box(THICK - eps, THICK + eps)
    lateral = [t for d, t in gmsh.model.getEntities(2)
            if t not in front and t not in back]

    def touching(curves):
        want = {abs(c) for c in curves}
        return [s for s in lateral
                if want & {abs(c) for _, c in
                        gmsh.model.getBoundary([(2, s)], oriented=False)}]

    gmsh.model.addPhysicalGroup(3, vols, name="internal")
    gmsh.model.addPhysicalGroup(2, touching(ring.wall_curves), name="airfoil")
    gmsh.model.addPhysicalGroup(2, front + back, name="frontAndBack")
    for nm, c in outer.boundary.items():
        gmsh.model.addPhysicalGroup(2, touching([c]), name=nm)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # gmshToFoam won't read 4.1
    gmsh.write("airfoil.msh")