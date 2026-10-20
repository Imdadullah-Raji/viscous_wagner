"""Padé location of the nearest singularity of the zero-thickness plate's steady C_L(alpha), from the
re-converged ("tight") data.

C_L is odd, so C_L = theta * g(z) with z = theta^2.  A pole of a Padé approximant of g at z_p puts a
singularity at theta = +-sqrt(z_p): z_p real > 0 -> on the real alpha axis; z_p < 0 -> on the imaginary
axis; complex z_p -> a complex-conjugate pair.  A branch point (a fold's sqrt) shows up as poles and
zeros strung along the cut, the nearest pole close to the branch point.

Two independent routes:
  A. classical Padé [L/M] of the Taylor series g = a1 + a3 z + a5 z^2 + ..., coefficients from the odd
     Legendre fits of analysis/cubic_stress_test.py (fit_window);
  B. a rational fit C_L = theta * P_L(z) / Q_M(z), Q(0) = 1, straight to the data: linearised least squares
     in C_L units, Sanathanan-Koerner reweighting by 1/|Q_prev|;
  C. AAA (scipy.interpolate.AAA; Nakatsukasa, Sete & Trefethen 2018): barycentric rational approximation
     with its own Froissart-doublet clean-up.  Independent of B's monomial basis and doublet test.  Fitted
     to C_L(theta) mirrored to +-theta (oddness built in, uniform noise) and to g(z); scanned in rtol.
     Below the noise AAA interpolates exactly (rms ~1e-15) and adds doublets of residue 1e-12..1e-7.
Route A (scipy.interpolate.pade) is kept as a record: it fails here because a9 and a11 are uncertain.
A pole counts only if it persists across windows and orders, survives noise of the fit's own rms (about
8e-9: the scatter between angles, larger than each angle's 1e-9 convergence), and has no zero within 1 %
of its distance (a Froissart doublet: a pole cancelled by a zero).  Note [1/1] and [2/1] are not new
information: their poles are a3/a5 and a5/a7, the two-term ratio estimates.

  python analysis/pade_singularity.py [RE] [--prod]    (research env)  -> analysis/pade_singularity/
"""
import sys
from pathlib import Path

import numpy as np
import warnings
from scipy.interpolate import AAA, pade
warnings.filterwarnings("ignore", category=DeprecationWarning)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cubic_stress_test import fit_window                            # noqa: E402
from tight_compare import cases                                      # noqa: E402

OUT = HERE / "pade_singularity"
RNG = np.random.default_rng(0)


def tight_data(re):
    """Usable tight angles (same criterion as tight_compare.py --drop)."""
    df = cases(re)
    ok = (df.d500.abs() <= 3e-9) & ((df.rem.abs() <= 3e-9) | df.rem.isna() & (df.d500.abs() <= 1e-9))
    return df[ok].reset_index(drop=True)


def theta_of(zp):
    """Singularity in alpha (deg) for a pole at z = theta^2: the root with Re >= 0."""
    t = np.sqrt(complex(zp))
    t = t if t.real >= 0 else -t
    return np.degrees(t.real), np.degrees(t.imag)


def rational_fit(th, cl, L, M, iters=20):
    """C_L = th * P(z)/Q(z), z = (th/s)^2 scaled for conditioning.  Returns poles, zeros in z (unscaled)
    and the C_L rms residual."""
    s = th.max(); z = (th / s) ** 2; y = cl / th
    w = th.copy()                                                    # residual in C_L units
    q = np.zeros(M)
    for _ in range(iters):
        Qp = 1 + sum(q[k] * z ** (k + 1) for k in range(M))
        A = np.column_stack([z ** k for k in range(L + 1)] + [-y * z ** (k + 1) for k in range(M)])
        W = (w / np.abs(Qp))[:, None]
        sol, *_ = np.linalg.lstsq(A * W, y * W[:, 0], rcond=None)
        q_new = sol[L + 1:]
        if np.allclose(q_new, q, rtol=1e-12, atol=1e-15):
            q = q_new; break
        q = q_new
    p = sol[:L + 1]
    Q = 1 + sum(q[k] * z ** (k + 1) for k in range(M))
    rms = np.sqrt(np.mean((th * np.polyval(p[::-1], z) / Q - cl) ** 2))
    poles = np.roots(np.r_[q[::-1], 1.0]) * s ** 2 if M else np.array([])
    zeros = np.roots(p[::-1]) * s ** 2 if L else np.array([])
    return poles, zeros, rms


def nearest(poles, zeros, genuine=False):
    """Nearest pole to the origin (in theta), and the distance to its nearest zero relative to |z_p|
    (< 1 %: a Froissart doublet).  genuine=True skips doublets."""
    best = None
    for zp in sorted(poles, key=abs):
        dz = np.min(np.abs(zeros - zp)) / abs(zp) if len(zeros) else np.inf
        if best is None:
            best = (zp, dz)
        if not genuine or dz >= 0.01:
            return zp, dz
    return None if genuine else best


def fmt(zp):
    a, b = theta_of(zp)
    kind = "real axis" if abs(np.imag(zp)) < 1e-9 * abs(zp) and np.real(zp) > 0 else \
           "imag axis" if abs(np.imag(zp)) < 1e-9 * abs(zp) else "complex"
    return f"{a:6.2f}{b:+6.2f}i deg  ({kind:9s}, |theta| {np.degrees(abs(np.sqrt(complex(zp)))):5.2f})"


def main():
    re = int(next((a for a in sys.argv[1:] if not a.startswith("--")), 500))
    OUT.mkdir(exist_ok=True)
    df = tight_data(re)
    if "--prod" in sys.argv:                                         # biased production values, same angles
        df = df.assign(cl=df["prod"])
    th_all, cl_all = np.radians(df.aoa.to_numpy()), df.cl.to_numpy()
    lines = [f"Padé, zero-thickness plate, Re {re}, {'PRODUCTION' if '--prod' in sys.argv else 'tight'} data: "
             f"{len(df)} angles {df.aoa.min():g}-{df.aoa.max():g} deg", ""]

    # ---- A: Padé of the Taylor coefficients
    lines.append("A. Padé [L/M] of g(z) = a1 + a3 z + a5 z^2 + ...  (coefficients from odd Legendre fits)")
    for degs in ((1, 3, 5, 7, 9), (1, 3, 5, 7, 9, 11)):
        for A in (3.0, 3.5, 4.0, df.aoa.max()):
            d = df[df.aoa <= A + 1e-9]
            if len(d) < len(degs) + 2:
                continue
            a, _, rms, n = fit_window(d, A, degs)
            c = [a[k] for k in degs]
            for M in range(1, len(c) // 2 + 1):
                L = len(c) - 1 - M
                p, q = pade(c, M, L)
                poles, zeros = np.roots(q.coeffs), np.roots(p.coeffs)
                r = nearest(poles, zeros)
                tag = "" if r is None else fmt(r[0]) + (f"  zero at {r[1]:.1%} -> DOUBLET" if r[1] < 0.01 else "")
                lines.append(f"  P1-P{max(degs):<2} 0-{A:<4g} [{L}/{M}]  {tag}")
    lines.append("")

    # ---- B: rational fit to the data, with noise and window tests
    lines.append("B. Rational fit C_L = theta P_L(z)/Q_M(z) to the data (rms in C_L; poles nearest the origin)")
    lines.append("   nearest pole; nearest GENUINE pole (doublets skipped); noise test: 200 refits with N(0, rms)"
                 " added, 16-84 % range of the genuine |theta_pole| (deg)")
    for L, M in ((1, 1), (2, 1), (1, 2), (2, 2), (3, 1), (1, 3), (3, 2), (2, 3)):
        for A in (3.0, 3.5, 4.0, df.aoa.max()):
            m = df.aoa.to_numpy() <= A + 1e-9
            th, cl = th_all[m], cl_all[m]
            if len(th) < L + M + 3:
                continue
            poles, zeros, rms = rational_fit(th, cl, L, M)
            r, g = nearest(poles, zeros), nearest(poles, zeros, genuine=True)
            if r is None:
                continue
            rad = []
            for _ in range(200):
                pn, zn, _ = rational_fit(th, cl + RNG.normal(0, rms, len(cl)), L, M)
                rn = nearest(pn, zn, genuine=True)
                rad.append(np.degrees(abs(np.sqrt(complex(rn[0])))) if rn is not None else np.nan)
            rad = np.array(rad); found = np.isfinite(rad).mean()
            near = "" if g is not None and g[0] == r[0] else f"nearest {fmt(r[0])[:16].strip()} is a DOUBLET; "
            gen = fmt(g[0]) if g is not None else "no genuine pole"
            nz = (f"{np.nanpercentile(rad, 16):5.2f}-{np.nanpercentile(rad, 84):5.2f}" if found else "  --  ") + \
                 ("" if found > 0.95 else f" (genuine in {found:.0%})")
            lines.append(f"  [{L}/{M}] 0-{A:<4g} n={len(th):2d} rms={rms:.1e}  {near}{gen}  noise: {nz}")
    # ---- C: AAA
    lines.append("")
    lines.append("C. AAA; poles with |theta| < 40 deg.  'genuine' = residue > 1e-4 x the largest near one; the rest are doublets")
    variants = (("C_L(theta), mirrored", np.r_[-th_all[::-1], 0, th_all], np.r_[-cl_all[::-1], 0, cl_all], False),
                ("g(z) = C_L/theta", th_all ** 2, cl_all / th_all, True))
    for name, x, y, zvar in variants:
        for rtol in (1e-5, 1e-6, 1e-7, 1e-8):
            r = AAA(x, y, rtol=rtol)
            pol, res, zer = r.poles(), np.abs(r.residues()), r.roots()
            def th_of(v):
                t = np.sqrt(complex(v)) if zvar else complex(v)
                return -t if t.real < 0 else t
            near = [q for v, q in zip(pol, res) if abs(th_of(v)) < np.radians(40)]
            big = max(near) if near else 0                           # largest residue among the near poles
            gen = sorted({float(round(np.degrees(abs(th_of(v))), 2)) for v, q in zip(pol, res)
                          if q > 1e-4 * big and abs(th_of(v)) < np.radians(40)})
            ndbl = sum(1 for v, q in zip(pol, res) if q <= 1e-4 * big and abs(th_of(v)) < np.radians(40))
            zs = sorted({float(round(np.degrees(abs(th_of(v))), 2)) for v in zer
                         if 0 < abs(th_of(v)) < np.radians(40) and abs(np.imag(th_of(v))) < 1e-6})
            rms = np.sqrt(np.mean((r(x) - y) ** 2))
            lines.append(f"  {name:22s} rtol {rtol:.0e}: m={len(r.support_points):2d} rms={rms:.1e}  genuine poles "
                         f"(|theta|, deg): {gen}   real zeros: {zs[:4]}   doublets: {ndbl}")
    txt = "\n".join(lines)
    print(txt)
    (OUT / f"pade_re{re}{'_prod' if '--prod' in sys.argv else ''}.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
