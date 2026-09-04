"""NACA 4-digit airfoil geometry with a truncated (blunt) trailing edge.

Pure numpy on purpose: no gmsh import lives in this module, so the shape can be
plotted, diffed and unit-tested without opening a gmsh session.

Conventions the mesher depends on -- do not change silently:
  * `upper` and `lower` both run LE -> TE.
  * Both arrays share an identical first row (the LE point).  The mesher must
    create that point ONCE and reuse its tag, or the two surface splines will
    not share a topological endpoint.
  * The TE face runs upper[-1] -> lower[-1].
  * Chord runs from x = 0 to x = `chord` after rescaling.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Airfoil:
    """Discrete airfoil boundary, split into the three curves the O-grid needs."""

    upper: np.ndarray  # (N, 2) LE -> TE
    lower: np.ndarray  # (N, 2) LE -> TE
    chord: float
    te_gap: float
    code: str

    @property
    def le(self) -> np.ndarray:
        return self.upper[0]

    @property
    def te_upper(self) -> np.ndarray:
        return self.upper[-1]

    @property
    def te_lower(self) -> np.ndarray:
        return self.lower[-1]

    @property
    def te_mid(self) -> np.ndarray:
        return 0.5 * (self.te_upper + self.te_lower)


def cosine_spacing(n: int, x_end: float = 1.0) -> np.ndarray:
    """Half-cosine distribution on [0, x_end], clustered at BOTH ends.

    x = x_end * (1 - cos(beta)) / 2 for beta in [0, pi].

    Clustering at x=0 resolves the sqrt(x) leading-edge singularity; clustering
    at x=x_end resolves the truncated trailing edge.  These are the points that
    define the SHAPE (the spline control points) -- they are not the mesh nodes.
    """
    beta = np.linspace(0.0, np.pi, n)
    return x_end * 0.5 * (1.0 - np.cos(beta))


def _thickness(x: np.ndarray, t: float) -> np.ndarray:
    """Half-thickness distribution.  Open-TE ('-0.1015') variant."""
    return 5.0 * t * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1015 * x**4
    )


def _camber(x: np.ndarray, m: float, p: float):
    """Mean camber line y_c(x) and its slope dy_c/dx."""
    if m == 0.0 or p == 0.0:
        return np.zeros_like(x), np.zeros_like(x)
    fore = x < p
    yc = np.where(
        fore,
        m / p**2 * (2 * p * x - x**2),
        m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2),
    )
    dyc = np.where(
        fore,
        2 * m / p**2 * (p - x),
        2 * m / (1 - p) ** 2 * (p - x),
    )
    return yc, dyc


def naca4(
    code: str = "0012",
    n_points: int = 200,
    x_trunc: float = 0.98,
    chord: float = 1.0,
) -> Airfoil:
    """Generate a NACA 4-digit section truncated at `x_trunc`, rescaled to `chord`.

    Truncation is applied by generating the cosine distribution over [0, x_trunc]
    rather than by generating the full section and cutting it -- this keeps the
    end-clustering exactly on the new trailing edge.  The whole section is then
    scaled uniformly by chord/x_trunc, which restores unit chord and preserves
    the thickness ratio.
    """
    if len(code) != 4 or not code.isdigit():
        raise ValueError(f"expected a 4-digit NACA code, got {code!r}")
    if not 0.80 < x_trunc <= 1.0:
        raise ValueError(f"x_trunc={x_trunc} outside the sensible range (0.80, 1.0]")

    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0

    x = cosine_spacing(n_points, x_trunc)
    yt = _thickness(x, t)
    yc, dyc = _camber(x, m, p)
    theta = np.arctan(dyc)

    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    upper = np.column_stack((xu, yu))
    lower = np.column_stack((xl, yl))

    # yt(0) == 0 exactly, so both surfaces already start at the origin; pin it
    # anyway so the shared-LE-point contract can never drift on a float edge.
    upper[0] = lower[0] = 0.0

    scale = chord / x_trunc
    upper *= scale
    lower *= scale

    te_gap = float(np.linalg.norm(upper[-1] - lower[-1]))
    if te_gap / chord < 2.0e-3:
        raise ValueError(
            f"TE gap {te_gap / chord:.2e}c is too thin to mesh across; "
            f"reduce x_trunc (currently {x_trunc})"
        )

    return Airfoil(upper=upper, lower=lower, chord=chord, te_gap=te_gap, code=code)


def diagnostics(af: Airfoil) -> dict:
    """Numbers worth eyeballing before handing the geometry to the mesher."""
    ds_u = np.linalg.norm(np.diff(af.upper, axis=0), axis=1)
    ds_l = np.linalg.norm(np.diff(af.lower, axis=0), axis=1)
    perim = ds_u.sum() + ds_l.sum() + af.te_gap
    return {
        "code": af.code,
        "n_points_per_surface": len(af.upper),
        "te_gap_over_c": af.te_gap / af.chord,
        "first_panel_at_le": float(min(ds_u[0], ds_l[0])),
        "last_panel_at_te": float(min(ds_u[-1], ds_l[-1])),
        "max_panel": float(max(ds_u.max(), ds_l.max())),
        "wetted_perimeter": float(perim),
        "max_thickness_over_c": float((af.upper[:, 1] - af.lower[:, 1]).max() / af.chord),
    }


if __name__ == "__main__":
    af = naca4("0012", n_points=200, x_trunc=0.98)
    for k, v in diagnostics(af).items():
        print(f"{k:>24s} : {v}")