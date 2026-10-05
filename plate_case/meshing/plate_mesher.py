"""Structured mesh around a finite-thickness flat plate (stadium section), all
hexes, transfinite out to the circular far field -- the plate counterpart of
the static-airfoil production mesh (`static_airfoil/meshing/mesher.py`).

Geometry follows `plate_v2.geo`: chord c, thickness th, semicircular LE and TE
of radius r = th*c/2.

WHY TWO RINGS.  A single TFI ring from a 2.3%-thick plate straight to an R=30
circle (exactly the airfoil recipe) shears the near-wall grid lines over the
flat and the LE round: max non-orthogonality 81 deg, 426 faces > 70 deg.  So:

  inner ring : plate -> PARALLEL OFFSET of the plate at distance d.  Offsetting
               a stadium gives a stadium with the same arc centres, so wall
               nodes and offset nodes share a normal: polar lines round the
               edges, vertical lines over the flat -- orthogonal by
               construction (same trick as plate_v2.geo).
  outer ring : offset stadium -> far-field circle, same 3-block TFI and the
               same corners / wake arc (`mesh_math.build_farfield`) as the
               airfoil.  Upper/lower far-field nodes are NOT uniform in angle:
               uniform nodes against the flat-heavy offset distribution shear
               the lines past the LE corner (71 deg).  `farfield_angles` puts
               them near the radial projection of their offset partners,
               blended 50/50 with uniform and smoothed (54 deg max).

Both rings use the airfoil's block layout (upper / lower / wake) and corners:

    upper : LE -> wall (LE arc, flat, flat, TE arc) -> TEu -> cut -> ...
    lower : mirror
    wake  : TEu -> TE cap arc -> TEl -> cut -> ... -> cut -> TEu

The TE cap half-angle phi is chosen so the TE cut (TE arc centre -> ff(theta_w))
is radial: it leaves the wall along the normal and stays one straight line
through both rings.  The LE cut leaves (0,0) along -x, already normal.

Radial: one geometric ratio g from h1 through the inner ring (n1 cells, which
fixes d), then each outer cut solves its own ratio for the remaining cells
starting from the next cell size h1*g^n1 -- no size jump at the interface.

Wall distribution (per segment, applied identically to the offset so the
lines stay normal):
    TE cap  : nte cells, uniform
    TE arc  : Progression from the cap spacing up to ~the LE spacing
    LE arc  : uniform
    flat    : two halves, each a Progression from its edge toward mid-chord

The wall patch is named `airfoil` on purpose so the mesh drops straight into
`static_airfoil/make_case.sh` / `make_steady.sh`, which retype that name to wall.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "static_airfoil" / "meshing"))
from mesh_math import build_farfield, geometric_length, geometric_ratio  # noqa: E402


def _transfinite(tag, n_cells, kind="Progression", coef=1.0):
    # gmsh counts NODES; we count CELLS.
    gmsh.model.geo.mesh.setTransfiniteCurve(tag, n_cells + 1, kind, coef)


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


def _prog(n, ratio):
    """Normalised node positions of a gmsh Progression with n cells: 0 .. 1."""
    if abs(ratio - 1.0) < 1e-12:
        return np.linspace(0.0, 1.0, n + 1)
    return (ratio ** np.arange(n + 1) - 1.0) / (ratio**n - 1.0)


def offset_side_nodes(c, r, rho, phi, n_le, nfh, n_tea, r_f_le, r_f_te, r_tea):
    """Upper offset-stadium nodes (radius rho), ordered TEu -> LE, i.e. in the
    direction the far-field upper arcs run (theta_w -> pi).  Mirrors exactly the
    distributions given to gmsh in `stadium()`."""
    xl, xt, half = r, c - r, (c - 2 * r) / 2
    a_te = phi + (np.pi / 2 - phi) * _prog(n_tea, r_tea)               # cap -> top
    te = np.column_stack((xt + rho * np.cos(a_te), rho * np.sin(a_te)))
    f2 = np.column_stack((xt - half * _prog(nfh, r_f_te), np.full(nfh + 1, rho)))
    f1 = np.column_stack((xl + half * _prog(nfh, r_f_le), np.full(nfh + 1, rho)))[::-1]
    a_le = np.linspace(np.pi / 2, np.pi, n_le + 1)
    le = np.column_stack((xl + rho * np.cos(a_le), rho * np.sin(a_le)))
    return np.vstack((te, f2[1:], f1[1:], le[1:]))


def farfield_angles(nodes, center, th0, th1, w_uniform, n_smooth):
    """Far-field node angles th0 -> th1 that put each far-field node roughly on
    the ray from the domain centre through its offset partner, so the outer TFI
    lines run near-radially instead of shearing past the LE corner.

    Pure projection makes the far-field spacing vary ~6x, so the increments are
    blended with uniform (weight w_uniform) and then smoothed; the sum, i.e. the
    end angles th0/th1 (the block corners), is preserved exactly."""
    a = np.unwrap(np.arctan2(nodes[:, 1], nodes[:, 0] - center[0]))
    a = th0 + (a - a[0]) * (th1 - th0) / (a[-1] - a[0])
    dproj = np.diff(a)
    n = len(dproj)
    d = w_uniform * (th1 - th0) / n + (1 - w_uniform) * dproj
    for _ in range(n_smooth):
        d[1:-1] = 0.25 * d[:-2] + 0.5 * d[1:-1] + 0.25 * d[2:]
    d *= (th1 - th0) / d.sum()
    return th0 + np.concatenate(([0.0], np.cumsum(d)))


def build(c, r, ff, ns, nte, nr, n1, h1, n_le, n_tea, w_uniform, n_smooth):
    geo = gmsh.model.geo
    tw = ff.theta_w
    xl, xt = r, c - r                        # LE / TE arc centres (on y = 0)

    # TE cap half-angle: TE cut radial from the TE arc centre.
    fx, fy = ff.point(tw)
    phi = float(np.arctan2(fy, fx - xt))

    # radial growth: g from the airfoil-style solve over the full LE cut
    L_le = np.linalg.norm(ff.point(np.pi))
    g = geometric_ratio(L_le, h1, nr)
    d = geometric_length(h1, g, n1)
    n2 = nr - n1
    h_if = h1 * g**n1

    # wall spacing
    dth_cap = 2 * phi / nte                              # rad per cap cell
    dth_le = (np.pi / 2) / n_le
    r_tea = _ratio(np.pi / 2 - phi, dth_cap, n_tea)      # TE arc, fine at the cap
    n_flat = ns - n_le - n_tea
    if n_flat < 2 or n_flat % 2:
        raise ValueError(f"ns - n_le - n_tea = {n_flat} must be even and >= 2")
    nfh = n_flat // 2
    half = (c - 2 * r) / 2
    ds_le = r * dth_le
    ds_te = r * dth_cap * r_tea ** (n_tea - 1)
    r_f_le = _ratio(half, ds_le, nfh)                    # flat, fine at LE end
    r_f_te = _ratio(half, ds_te, nfh)                    # flat, fine at TE end

    # ---- points ---------------------------------------------------------
    P = lambda x, y: geo.addPoint(x, y, 0.0)
    c_le, c_te = P(xl, 0), P(xt, 0)
    c_ff = P(*ff.center)

    def ring_pts(rho):
        """Stadium of edge radius rho: LE, then per side (TL, M, TR, TE-corner)."""
        pts = {"le": P(xl - rho, 0)}
        for s, name in ((1, "u"), (-1, "l")):
            pts["tl" + name] = P(xl, s * rho)
            pts["m" + name] = P(0.5 * c, s * rho)
            pts["tr" + name] = P(xt, s * rho)
            pts["te" + name] = P(xt + rho * np.cos(phi), s * rho * np.sin(phi))
        return pts

    cache = {}

    def cp(theta):
        key = int(round((theta % (2 * np.pi)) * 1e9))
        if key not in cache:
            cache[key] = P(*ff.point(theta))
        return cache[key]

    # ---- stadium curves (wall and offset get identical distributions) ----
    def stadium(p):
        """Returns {'u': [...], 'l': [...]} as signed curve lists running
        LE -> TE-corner, plus the cap curve TEu -> TEl."""
        side = {}
        for k in ("u", "l"):
            a = geo.addCircleArc(p["le"], c_le, p["tl" + k])      # LE arc
            _transfinite(a, n_le)
            f1 = geo.addLine(p["tl" + k], p["m" + k])             # fine at LE end
            _transfinite(f1, nfh, "Progression", r_f_le)
            f2 = geo.addLine(p["tr" + k], p["m" + k])             # fine at TE end
            _transfinite(f2, nfh, "Progression", r_f_te)
            t = geo.addCircleArc(p["te" + k], c_te, p["tr" + k])  # fine at cap
            _transfinite(t, n_tea, "Progression", r_tea)
            side[k] = [a, f1, -f2, -t]
        cap = geo.addCircleArc(p["teu"], c_te, p["tel"])
        _transfinite(cap, nte)
        return side, cap

    wp, op = ring_pts(r), ring_pts(r + d)
    w_side, w_cap = stadium(wp)
    o_side, o_cap = stadium(op)

    # ---- cuts ------------------------------------------------------------
    def cut(a, b, n, ratio):
        t = geo.addLine(a, b)
        _transfinite(t, n, "Progression", ratio)
        return t

    ci = {k: cut(wp[k], op[k], n1, g) for k in ("le", "teu", "tel")}
    ff_pt = {"le": np.pi, "teu": tw, "tel": -tw}
    co, ratios = {}, {}
    for k, th in ff_pt.items():
        x0 = np.array([xl - (r + d), 0.0]) if k == "le" else \
            np.array([xt + (r + d) * np.cos(phi), (1 if k == "teu" else -1) * (r + d) * np.sin(phi)])
        L = float(np.linalg.norm(ff.point(th) - x0))
        ratios[k] = geometric_ratio(L, h_if, n2)
        co[k] = cut(op[k], cp(th), n2, ratios[k])

    # ---- far-field arcs ---------------------------------------------------
    # Wake arc: uniform, as for the airfoil.  Upper/lower: one arc per cell at
    # the angles from farfield_angles(), so every node sits exactly on the
    # circle.  A transfinite side may consist of any number of curves as long
    # as the corners are named explicitly.
    nodes = offset_side_nodes(c, r, r + d, phi, n_le, nfh, n_tea, r_f_le, r_f_te, r_tea)
    th_up = farfield_angles(nodes, ff.center, tw, np.pi, w_uniform, n_smooth)
    arcs = {"wake": [], "upper": [], "lower": []}
    for a in ff.arcs:
        if a.role == "wake":
            t = geo.addCircleArc(cp(a.theta_start), c_ff, cp(a.theta_end))
            _transfinite(t, a.n_cells)
            arcs["wake"].append(t)
    for t0, t1 in zip(th_up[:-1], th_up[1:]):
        t = geo.addCircleArc(cp(t0), c_ff, cp(t1))
        _transfinite(t, 1)
        arcs["upper"].append(t)
    th_lo = 2 * np.pi - th_up[::-1]                      # pi -> 2pi - tw
    for t0, t1 in zip(th_lo[:-1], th_lo[1:]):
        t = geo.addCircleArc(cp(t0), c_ff, cp(t1))
        _transfinite(t, 1)
        arcs["lower"].append(t)
    outer = arcs["wake"] + arcs["upper"] + arcs["lower"]
    dth = np.diff(th_up)

    rev = lambda cs: [-x for x in reversed(cs)]

    # ---- blocks ------------------------------------------------------------
    def blocks(inner_side, inner_cap, ip, cut_le, cut_teu, cut_tel,
               out_up, out_lo, out_wake, out_cap, opt):
        """Generic 3-block ring between an inner stadium and an outer boundary
        (either the offset stadium or the far-field arcs)."""
        return {
            "upper": ([*inner_side["u"], cut_teu, *out_up, -cut_le],
                      [ip["le"], ip["teu"], opt["teu"], opt["le"]]),
            "lower": ([cut_le, *out_lo, -cut_tel, *rev(inner_side["l"])],
                      [ip["le"], opt["le"], opt["tel"], ip["tel"]]),
            "wake": ([inner_cap, cut_tel, *out_wake, -cut_teu],
                     [ip["teu"], ip["tel"], opt["tel"], opt["teu"]]),
        }

    inner = blocks(w_side, w_cap, wp, ci["le"], ci["teu"], ci["tel"],
                   rev(o_side["u"]), o_side["l"], [-o_cap], None, op)
    ffc = {"le": cp(np.pi), "teu": cp(tw), "tel": cp(-tw)}
    outer_b = blocks(o_side, o_cap, op, co["le"], co["teu"], co["tel"],
                     arcs["upper"], arcs["lower"], arcs["wake"], None, ffc)

    surfaces = []
    for blk in (inner, outer_b):
        for curves, corners in blk.values():
            s = geo.addPlaneSurface([geo.addCurveLoop(curves)])
            geo.mesh.setTransfiniteSurface(s, "Left", corners)
            geo.mesh.setRecombine(2, s)
            surfaces.append(s)

    geo.synchronize()
    wall = [abs(x) for x in w_side["u"] + w_side["l"]] + [w_cap]
    info = {
        "phi_deg": np.degrees(phi), "theta_w_deg": np.degrees(tw),
        "g_radial": g, "d_offset": d, "h_interface": h_if,
        **{f"r_outer_{k}": v for k, v in ratios.items()},
        "ds_cap": r * dth_cap, "ds_le_arc": ds_le, "ds_te_arc_end": ds_te,
        "r_te_arc": r_tea, "r_flat_le": r_f_le, "r_flat_te": r_f_te,
        "ff_ds_min": ff.radius * dth.min(), "ff_ds_max": ff.radius * dth.max(),
        "ff_adj_ratio_max": float(np.max(np.maximum(dth[1:] / dth[:-1], dth[:-1] / dth[1:]))),
        "ds_flat_mid_max": max(ds_le * r_f_le ** (nfh - 1), ds_te * r_f_te ** (nfh - 1)),
    }
    return surfaces, wall, outer, info


def main():
    ap = argparse.ArgumentParser(description="Structured all-hex mesh around a thick flat plate")
    ap.add_argument("filename", help="output .msh (v2.2, for gmshToFoam)")
    ap.add_argument("--re", type=float, default=500, help="Reynolds number (sets h1)")
    ap.add_argument("-r", "--radius", type=float, default=30.0, help="far-field radius")
    ap.add_argument("--th", type=float, default=0.023, help="thickness / chord")
    ap.add_argument("--ns", type=int, default=160, help="cells along each plate side")
    ap.add_argument("--nte", type=int, default=32, help="cells across the TE cap")
    ap.add_argument("--nr", type=int, default=180, help="radial cells (total)")
    ap.add_argument("--n1", type=int, default=60, help="radial cells in the inner (offset) ring")
    ap.add_argument("--n-le", type=int, default=16, help="cells on each LE quarter arc")
    ap.add_argument("--n-tea", type=int, default=24, help="cells on each TE arc (cap -> flat)")
    ap.add_argument("--w-uniform", type=float, default=0.5,
                    help="far-field spacing: 0 = pure radial projection, 1 = uniform")
    ap.add_argument("--n-smooth", type=int, default=50, help="smoothing passes on far-field spacing")
    ap.add_argument("--center", type=float, default=0.25, help="far-field centre x")
    ap.add_argument("--span", type=float, default=0.1, help="extrusion depth (Aref = c*span)")
    ap.add_argument("--h1-coef", type=float, default=0.05,
                    help="first cell h1 = coef/sqrt(Re); 0.05 is the production rule")
    args = ap.parse_args()

    c = 1.0
    r = args.th * c / 2
    h1 = args.h1_coef / np.sqrt(args.re)  # same rule as the airfoil meshes
    ff = build_farfield(center=(args.center, 0.0), radius=args.radius, ns=args.ns, nte=args.nte)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("plate")

    surfaces, wall, outer, info = build(c, r, ff, args.ns, args.nte, args.nr, args.n1,
                                        h1, args.n_le, args.n_tea, args.w_uniform, args.n_smooth)

    geo = gmsh.model.geo
    # ONE extrude call so blocks share their interface faces.
    out = geo.extrude([(2, s) for s in surfaces], 0, 0, args.span,
                      numElements=[1], recombine=True)
    geo.synchronize()
    vols = [t for d, t in out if d == 3]

    eps = args.span * 1e-3
    box = lambda z0, z1: [t for d, t in gmsh.model.getEntitiesInBoundingBox(
        -1e9, -1e9, z0, 1e9, 1e9, z1, 2)]
    front, back = box(-eps, eps), box(args.span - eps, args.span + eps)
    lateral = [t for d, t in gmsh.model.getEntities(2) if t not in front and t not in back]

    def touching(curves):
        want = {abs(x) for x in curves}
        return [s for s in lateral
                if want & {abs(x) for _, x in gmsh.model.getBoundary([(2, s)], oriented=False)}]

    gmsh.model.addPhysicalGroup(3, vols, name="internal")
    gmsh.model.addPhysicalGroup(2, touching(wall), name="airfoil")
    gmsh.model.addPhysicalGroup(2, touching(outer), name="farfield")
    gmsh.model.addPhysicalGroup(2, front + back, name="frontAndBack")

    gmsh.model.mesh.generate(3)
    hexes = len(gmsh.model.mesh.getElementsByType(5)[0])
    others = sum(len(gmsh.model.mesh.getElementsByType(t)[0]) for t in (4, 6, 7))
    expected = (2 * args.ns + args.nte) * args.nr
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # gmshToFoam won't read 4.1
    gmsh.write(args.filename)
    gmsh.finalize()

    print(f"\n  h1               : {h1:.4e}")
    for k, v in info.items():
        print(f"  {k:<17s}: {v:.6g}")
    print(f"  hexes            : {hexes}  (expected {expected})")
    print(f"  non-hex cells    : {others}  (must be 0)")
    assert others == 0 and hexes == expected


if __name__ == "__main__":
    main()
