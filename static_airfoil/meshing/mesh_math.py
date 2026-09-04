"""Numerical helpers for laying out the O-grid.  No gmsh, no geometry classes --
just the arithmetic, so every number the mesher uses can be checked on its own.

Three jobs:
  1. geometric_ratio  -- turn (length, first-cell height, cell count) into the
     progression ratio gmsh wants.
  2. wake_half_angle  -- where the two TE cut lines land on the far field.
  3. FarField         -- partition the circle into arcs that gmsh will accept
     and that carry matching node counts.
"""

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# 1. geometric progression
# ---------------------------------------------------------------------------

def geometric_length(h1: float, r: float, n: int) -> float:
    """Total length of `n` cells starting at `h1` with ratio `r`."""
    if abs(r - 1.0) < 1e-14:
        return h1 * n
    return h1 * (r**n - 1.0) / (r - 1.0)


def geometric_ratio(
    length: float,
    h1: float,
    n: int,
    r_cap: float = 4.0,
    tol: float = 1e-13,
    max_iter: int = 300,
) -> float:
    """Solve  length = h1 * (r^n - 1) / (r - 1)  for r > 1.

    Bisection rather than Newton: the residual is smooth and monotonic in r
    (dL/dr > 0 for r > 1), the bracket is trivially available, and 300 halvings
    cost nothing.  Newton on r^n is easy to make diverge for large n.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if h1 <= 0.0 or length <= 0.0:
        raise ValueError("length and h1 must be positive")

    uniform = h1 * n
    if length < uniform:
        raise ValueError(
            f"cannot reach L={length:.4g} with n={n} cells of first height "
            f"h1={h1:.4g}: uniform spacing already overshoots ({uniform:.4g}). "
            f"Reduce n or h1."
        )
    if abs(length - uniform) < tol * length:
        return 1.0

    lo, hi = 1.0, 1.0 + 1e-6
    while geometric_length(h1, hi, n) < length:
        hi = 1.0 + (hi - 1.0) * 2.0
        if hi > r_cap:
            raise ValueError(
                f"need r > {r_cap} to span L={length:.4g} with n={n} cells "
                f"from h1={h1:.4g}; increase n."
            )

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if geometric_length(h1, mid, n) < length:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def layer_heights(h1: float, r: float, n: int) -> np.ndarray:
    """The n individual cell heights -- for plotting / eyeballing the stretch."""
    return h1 * r ** np.arange(n)


# ---------------------------------------------------------------------------
# 2. wake fan angle
# ---------------------------------------------------------------------------

def wake_half_angle(ns: int, nte: int) -> float:
    """Half-angle (rad) subtended by the wake block at the far field.

    Chosen so the far-field circle is discretised UNIFORMLY: if each side arc
    carries `ns` cells over (pi - theta) and the wake arc carries `nte` cells
    over 2*theta, equal arc-length-per-cell gives

        2*theta / nte = (pi - theta) / ns   =>   theta = pi*nte / (2*ns + nte)

    Any other choice leaves a spacing jump at the cut lines, which shows up as
    a visible kink in the outer grid lines.
    """
    if ns < 1 or nte < 1:
        raise ValueError("ns and nte must be >= 1")
    return np.pi * nte / (2.0 * ns + nte)


# ---------------------------------------------------------------------------
# 3. far-field partition
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Arc:
    theta_start: float
    theta_end: float
    n_cells: int
    role: str  # 'upper' | 'wake' | 'lower'

    @property
    def sweep(self) -> float:
        return self.theta_end - self.theta_start


@dataclass(frozen=True)
class FarField:
    center: np.ndarray
    radius: float
    theta_w: float
    arcs: List[Arc]

    def point(self, theta: float) -> np.ndarray:
        return self.center + self.radius * np.array([np.cos(theta), np.sin(theta)])

    @property
    def node_angles(self) -> List[float]:
        """Every arc endpoint, in order -- these become gmsh points."""
        out = [self.arcs[0].theta_start]
        out.extend(a.theta_end for a in self.arcs)
        return out


def _split_evenly(theta_a: float, theta_b: float, n_cells: int,
                  max_arc: float, role: str) -> List[Arc]:
    """Chop [theta_a, theta_b] into equal sub-arcs below `max_arc`, splitting
    `n_cells` between them as evenly as possible (sum is preserved exactly)."""
    sweep = theta_b - theta_a
    k = max(1, int(np.ceil(abs(sweep) / max_arc)))
    if k > n_cells:
        raise ValueError(f"{role} arc needs {k} sub-arcs but only has {n_cells} cells")
    base, rem = divmod(n_cells, k)
    edges = np.linspace(theta_a, theta_b, k + 1)
    return [
        Arc(edges[i], edges[i + 1], base + (1 if i < rem else 0), role)
        for i in range(k)
    ]


def build_farfield(
    center: Tuple[float, float],
    radius: float,
    ns: int,
    nte: int,
    max_arc: float = np.pi / 2,
) -> FarField:
    """Partition the far-field circle to match the airfoil's three-curve boundary.

    `max_arc` defaults to pi/2 and exists because gmsh's circle-arc primitive
    requires a sweep strictly less than pi.  The upper side arc spans
    (pi - theta_w), which for a small wake block sits at ~174 deg -- legal by a
    hair, and the kind of thing that breaks the moment you change nte.  Splitting
    at pi/2 removes the failure mode entirely; a transfinite surface accepts a
    multi-curve side as long as the corners are named explicitly.

    Arcs are returned counter-clockwise starting at -theta_w (bottom of wake).
    """
    theta_w = wake_half_angle(ns, nte)
    c = np.asarray(center, dtype=float)

    arcs: List[Arc] = []
    arcs.append(Arc(-theta_w, theta_w, nte, "wake"))
    arcs.extend(_split_evenly(theta_w, np.pi, ns, max_arc, "upper"))
    arcs.extend(_split_evenly(np.pi, 2 * np.pi - theta_w, ns, max_arc, "lower"))

    return FarField(center=c, radius=float(radius), theta_w=theta_w, arcs=arcs)


# ---------------------------------------------------------------------------
# convenience: cut-line radial solves
# ---------------------------------------------------------------------------

def cut_line_ratio(
    p_inner: np.ndarray,
    ff: FarField,
    theta: float,
    h1: float,
    nr: int,
) -> Tuple[float, float]:
    """Growth ratio for one cut line, plus its length.

    Solved per line against a COMMON h1 so that near-wall resolution is
    identical on all three cuts and stays fixed when the far-field radius is
    varied for a domain-independence study.
    """
    length = float(np.linalg.norm(ff.point(theta) - np.asarray(p_inner, dtype=float)))
    return geometric_ratio(length, h1, nr), length


if __name__ == "__main__":
    ns, nte, nr = 100, 30, 100
    ff = build_farfield(center=(0.25, 0.0), radius=10.0, ns=ns, nte=nte)

    print(f"theta_w = {np.degrees(ff.theta_w):.3f} deg, {len(ff.arcs)} arcs")
    for a in ff.arcs:
        print(f"  {a.role:>5s}  {np.degrees(a.theta_start):8.2f} -> "
              f"{np.degrees(a.theta_end):8.2f} deg   {a.n_cells:3d} cells   "
              f"ds = {ff.radius * abs(a.sweep) / a.n_cells:.5f}")

    h1 = 7.0e-4
    for name, p, th in [("LE ", (0.0, 0.0), np.pi),
                        ("TEu", (1.0, 0.0041), ff.theta_w),
                        ("TEl", (1.0, -0.0041), -ff.theta_w)]:
        r, L = cut_line_ratio(p, ff, th, h1, nr)
        h_last = layer_heights(h1, r, nr)[-1]
        print(f"{name}: L = {L:7.4f}  r = {r:.6f}  h_last = {h_last:.4f}")