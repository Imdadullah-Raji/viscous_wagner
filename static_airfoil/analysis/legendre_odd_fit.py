"""Odd Legendre fit of the steady lift curve on a symmetric, rescaled angle.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/legendre_odd_fit.py [--amax 8]

x = theta / theta_max in [0, 1], mirrored to [-1, 1] by the NACA0012 symmetry
C_L(-x) = -C_L(x) (x = 0 kept once).  Least squares on P1, P3, P5, P7, P9:

    C_L(x) = c1 P1(x) + c3 P3(x) + c5 P5(x) + c7 P7(x) + c9 P9(x)

No intercept: P0 is even, and the symmetry forbids it.  Also reports the
monomial (theta, theta^3, ...) equivalent, the truncated fits P1..Pk, and a
check that the mirroring leaves the odd-basis solution unchanged (it must: the
residual at -x is minus the residual at x, so mirroring only doubles the sum).
-> analysis/figs/legendre_odd_0to<amax>deg/
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.polynomial import legendre as L

from lasso_re1000 import collect, ROOT

RES = (100, 500, 1000, 2000)
COLOR = {500: "#2a78d6", 1000: "#eb6834", 100: "#1baf7a", 2000: "#eda100"}
ODD = (1, 3, 5, 7, 9)


def basis(x, degs=ODD):
    return np.column_stack([L.legval(x, np.eye(max(degs) + 1)[k]) for k in degs])


def data(re, amax):
    df = collect(re)
    d = df[(df.aoa <= amax + 1e-9) & np.isclose((df.aoa * 2) % 1, 0)
           & (df.verdict == "CONVERGED")]
    x = d.aoa.to_numpy() / amax
    return x, d.cl.to_numpy()


def mirror(x, y):
    m = x > 0
    return np.r_[-x[m][::-1], x], np.r_[-y[m][::-1], y]


def fit(x, y, degs=ODD):
    c, *_ = np.linalg.lstsq(basis(x, degs), y, rcond=None)
    r = y - basis(x, degs) @ c
    return c, np.sqrt(np.mean(r ** 2)), np.abs(r).max()


def to_theta_monomial(c, amax, degs=ODD):
    """Legendre-in-x coefficients -> monomial coefficients in theta (rad).

    "Collecting" the series: every P_k spills into all lower odd powers
    (P3 = (5x^3 - 3x)/2 feeds theta AND theta^3), so truncating in Legendre
    and truncating in theta are different models.  Even degrees (an even
    series such as C_D) return the even powers, theta^0 upward."""
    top = max(max(ODD), max(degs))               # P11 and up for the a5/a7 fits
    full = np.zeros(top + 1)
    full[list(degs)] = c
    mono = L.leg2poly(full)                      # in x; trims trailing zeros
    mono = np.pad(mono, (0, top + 1 - len(mono)))
    tmax = np.deg2rad(amax)
    return {k: mono[k] / tmax ** k for k in range(degs[0] % 2, top + 1, 2)}


def fig_theta(collected, amax, out):
    """Edge contribution a_k * theta_max^k (= coefficient of x^k) of each theta
    power, for the collected P1-P3, P1-P5 and P1-P9 fits; one panel per Re."""
    levels = {"P1-P3": "#2a78d6", "P1-P5": "#eb6834", "P1-P9": "#1baf7a"}
    tmax = np.deg2rad(amax)
    fig, axes = plt.subplots(1, len(RES), figsize=(15, 4.2), sharey=True)
    w = 0.26
    for ax, re in zip(axes, RES):
        for i, (lev, col) in enumerate(levels.items()):
            a = collected[(re, lev)]
            v = [a[k] * tmax ** k for k in ODD]
            ax.bar(np.arange(len(ODD)) + (i - 1) * w, v, w * 0.92, color=col, label=lev)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(range(len(ODD)), [rf"$\theta^{k}$" if k > 1 else r"$\theta$" for k in ODD])
        ax.set_title(f"Re = {re}", loc="left")
        ax.grid(alpha=0.25, lw=0.6, axis="y")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel(r"contribution at $\theta_{max}$:  $a_k\,\theta_{max}^k$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, [f"Legendre {x}, collected" for x in l], loc="lower center",
               ncol=3, frameon=False)
    fig.suptitle(rf"Odd Legendre fits rewritten in powers of $\theta$ "
                 rf"($\theta_{{max}}$ = {amax:g}°)", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.08, 1, 0.97))
    fig.savefig(out / "theta_basis_collected.png", dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=float, default=8.0)
    amax = ap.parse_args().amax
    out = ROOT / "analysis" / "figs" / f"legendre_odd_0to{amax:g}deg"
    out.mkdir(parents=True, exist_ok=True)

    rows, trunc, lines, collected = [], [], [], {}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for re in RES:
        x0, y0 = data(re, amax)
        x, y = mirror(x0, y0)
        c, rmse, mx = fit(x, y)
        c_half, *_ = fit(x0, y0)                 # unmirrored: must be identical
        assert np.allclose(c, c_half, rtol=1e-10, atol=1e-13), (re, c, c_half)
        cond = np.linalg.cond(basis(x))
        mono = to_theta_monomial(c, amax)
        rows.append(dict(re=re, n_half=len(x0), n_mirrored=len(x),
                         **{f"c{k}": v for k, v in zip(ODD, c)},
                         rmse=rmse, max_abs_res=mx, cond=cond,
                         **{f"theta^{k}": mono[k] for k in ODD}))
        for kmax in ODD:
            degs = tuple(k for k in ODD if k <= kmax)
            ct, r_t, m_t = fit(x, y, degs)
            th = to_theta_monomial(ct, amax, degs)
            collected[(re, f"P1-P{kmax}")] = th
            trunc.append(dict(re=re, terms="+".join(f"P{k}" for k in degs),
                              **{f"c{k}": v for k, v in zip(degs, ct)},
                              rmse=r_t, max_abs_res=m_t,
                              **{f"a{k}": th[k] for k in ODD if k <= kmax}))
        a1.plot(ODD, np.abs(c), "o-", color=COLOR[re], lw=2, ms=6, label=f"Re = {re}")
        g = np.linspace(-1, 1, 401)
        a2.plot(x, 1e5 * (y - basis(x) @ c), "o-", color=COLOR[re], ms=3, lw=1,
                label=f"Re = {re}")

    t = pd.DataFrame(rows)
    tt = pd.DataFrame(trunc)
    t.to_csv(out / "legendre_odd_coeffs.csv", index=False)
    tt.to_csv(out / "legendre_truncation.csv", index=False)
    pd.set_option("display.width", 200)
    fmt = lambda v: f"{v:+.6e}" if isinstance(v, float) else str(v)
    txt = (f"x = theta/theta_max, theta_max = {amax:g} deg; mirrored to [-1, 1]\n\n"
           + t[["re", "n_half", "n_mirrored", "c1", "c3", "c5", "c7", "c9",
                "rmse", "max_abs_res", "cond"]].to_string(index=False, float_format=fmt)
           + "\n\nmonomial equivalent in theta (rad):\n"
           + t[["re"] + [f"theta^{k}" for k in ODD]].to_string(index=False, float_format=fmt)
           + "\n\ntruncated fits:\n"
           + tt.to_string(index=False, float_format=fmt, na_rep=""))
    (out / "legendre_odd_fit.txt").write_text(txt + "\n")
    print(txt)

    a1.set(yscale="log", xticks=ODD, xlabel="Legendre degree k", ylabel=r"$|c_k|$")
    a1.set_title(rf"Odd Legendre coefficients, $x=\theta/\theta_{{max}}$, $\theta_{{max}}$={amax:g}°",
                 loc="left")
    a2.axhline(0, color="black", lw=0.6)
    a2.set(xlabel=r"$x = \theta/\theta_{max}$ (mirrored)", ylabel=r"residual $\times 10^5$")
    a2.set_title("Residual of the P1-P9 fit", loc="left")
    for ax in (a1, a2):
        ax.grid(alpha=0.25, lw=0.6, which="both")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "legendre_odd_coeffs.png", dpi=200)
    plt.close(fig)
    fig_theta(collected, amax, out)


if __name__ == "__main__":
    main()
