"""Stress test of "C_L is cubic at the origin" for the four bodies at Re = 500, and the
zero-thickness plate's a1(Re), a3(Re) so far.

An analytic odd C_L = a1 th + a3 th^3 + a5 th^5 + ... always LOOKS cubic near 0, so the questions are:
  1. are the coefficients well defined (stable across fit windows), and to what order;
  2. is the scatter of a5, a7, a9 set by noise in the data or by content beyond the fit order
     (window dependence) -- the resolution limit;
  3. how far does each truncation, fitted on 0-4 deg, extrapolate (first miss > 1e-3);
  4. what the two-term ratio sqrt|a3/a5| says about the radius of convergence.

Fits: odd Legendre P1..Pk on x = alpha/A over 0-A (one-sided, odd basis: identical to the mirrored
fit), converted to powers of theta.  Standard errors: least squares, sigma from the fit's own
residual (n - p dof), carried through the Legendre -> monomial map.

  python3 analysis/cubic_stress_test.py   (research env)  -> analysis/cubic_stress_test/
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lasso_re1000 import collect                                   # noqa: E402
from legendre_odd_fit import basis, to_theta_monomial               # noqa: E402

OUT = HERE / "cubic_stress_test"
BODIES = {  # Re = 500, folders relative to static_airfoil/
    "NACA0012": ("NACA0012/steady_re500", "#2a78d6"),
    "NACA0006": ("NACA0006/steady_re500", "#16a34a"),
    "plate t/c=0.023": ("../plate_case/re500", "#9333ea"),
    "zero-thickness plate": ("zero_thickness_plate/ell_re500", "#c2410c"),
}
ZTP_RE = {200: "ell_re_scaling/steady_re200", 300: "ell_re_scaling/steady_re300",
          400: "ell_re_scaling/steady_re400", 500: "ell_re500", 600: "ell_re_scaling/steady_re600",
          700: "ell_re_scaling/steady_re700", 800: "ell_re_scaling/steady_re800",
          900: "ell_re_scaling/steady_re900"}
OTHER_RE = {  # existing a3(Re) tables
    "plate t/c=0.023": (HERE.parent.parent / "plate_case/analysis/a3_re_scaling/a3_vs_re.csv", "#9333ea"),
    "NACA0006": (HERE / "a3_re_scaling/NACA0006/a3_vs_re.csv", "#16a34a"),
    "NACA0012": (HERE / "a3_re_scaling/NACA0012/a3_vs_re.csv", "#2a78d6"),
}


def load(base, re=500):
    df = collect(re, base)
    return df[df.verdict == "CONVERGED"].sort_values("aoa").reset_index(drop=True)


def fit_window(df, A, degs):
    """One-sided odd fit on 0 < alpha <= A.  Returns monomial coefficients, their standard errors,
    the rms residual and the number of points."""
    d = df[(df.aoa > 1e-9) & (df.aoa <= A + 1e-9)]
    x, y = d.aoa.to_numpy() / A, d.cl.to_numpy()
    X = basis(x, degs)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c
    n, p = len(y), len(degs)
    s2 = (r @ r) / max(n - p, 1)
    cov_c = s2 * np.linalg.inv(X.T @ X)
    # Legendre -> monomial map is linear: column j = image of unit vector e_j
    T = np.array([[to_theta_monomial(np.eye(p)[j], A, degs)[k] for j in range(p)] for k in degs])
    a = T @ c
    se = np.sqrt(np.diag(T @ cov_c @ T.T))
    return dict(zip(degs, a)), dict(zip(degs, se)), np.sqrt(np.mean(r ** 2)), n


def main():
    OUT.mkdir(exist_ok=True)
    rows, lines = [], []
    data = {name: load(base) for name, (base, _) in BODIES.items()}

    # ---- 1-2: coefficients per window, two fit orders, with standard errors
    for name, df in data.items():
        amax = df.aoa.max()
        for A in [w for w in (3, 4, 5, 6, 8, 10, 12) if w <= amax + 1e-9]:
            for degs in ((1, 3, 5, 7), (1, 3, 5, 7, 9)):
                a, se, rms, n = fit_window(df, A, degs)
                row = dict(body=name, A=A, order=max(degs), n=n, rms=rms)
                for k in degs:
                    row[f"a{k}"], row[f"se{k}"] = a[k], se[k]
                row["R35_deg"] = np.degrees(np.sqrt(abs(a[3] / a[5])))
                rows.append(row)
    C = pd.DataFrame(rows)
    C.to_csv(OUT / "coefficients_by_window.csv", index=False)

    lines.append("Coefficients by window (one-sided odd Legendre P1..Pk -> theta powers), Re = 500")
    lines.append("se = least-squares standard error from the fit's own residual; rms = residual rms\n")
    for name in BODIES:
        lines.append(f"== {name} (data 0-{data[name].aoa.max():g} deg, {len(data[name])} angles)")
        for order in (7, 9):
            s = C[(C.body == name) & (C.order == order)]
            for _, r in s.iterrows():
                ks = [1, 3, 5, 7] + ([9] if order == 9 else [])
                lines.append(f"  P1-P{order} 0-{r.A:>2g} (n={r.n:2d}, rms {r.rms:.1e}): " +
                             "  ".join(f"a{k}={r[f'a{k}']:+.4g}±{r[f'se{k}']:.1g}" for k in ks) +
                             f"   sqrt|a3/a5|={r.R35_deg:.1f}deg")
        lines.append("")

    # window-to-window spread vs standard error (P1-P7, windows 4..6 deg common to all bodies)
    lines.append("Resolution: window-to-window spread (windows 4-6 deg) vs typical standard error, both fit orders")
    lines.append("spread/se ~ 1: noise-limited (coefficient resolved);  >> 1: content beyond the fit order (window-dependent)")
    res_rows = []
    for name in BODIES:
        for order in (7, 9):
            s = C[(C.body == name) & (C.order == order) & (C.A >= 4) & (C.A <= 6)]
            for k in (3, 5, 7, 9)[: (order - 1) // 2]:
                spread = s[f"a{k}"].max() - s[f"a{k}"].min()
                se = s[f"se{k}"].median()
                res_rows.append(dict(body=name, fit=f"P1-P{order}", coef=f"a{k}", mean=s[f"a{k}"].mean(), spread=spread,
                                     se=se, spread_over_se=spread / se if se > 0 else np.nan))
    R = pd.DataFrame(res_rows)
    R.to_csv(OUT / "resolution.csv", index=False)
    lines.append(R.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    lines.append("")

    # ---- 3: extrapolation of truncations fitted on 0-4 deg
    lines.append("Extrapolation: truncation fitted on 0-4 deg, first angle where |fit - data| > 1e-3")
    ext = {}
    for name, df in data.items():
        th = np.radians(df.aoa.to_numpy())
        ext[name] = {}
        for degs in ((1, 3), (1, 3, 5), (1, 3, 5, 7)):
            a, *_ = fit_window(df, 4, degs)
            err = np.abs(sum(a[k] * th ** k for k in degs) - df.cl.to_numpy())
            bad = df.aoa.to_numpy()[err > 1e-3]
            ext[name][degs] = (df.aoa.to_numpy(), err)
            lab = "+".join(f"th^{k}" for k in degs)
            lines.append(f"  {name:22s} {lab:22s}: " +
                         (f"first miss at {bad.min():g} deg" if len(bad) else f"within 1e-3 to {df.aoa.max():g} deg (end of data)"))
    lines.append("")

    # ---- zero-thickness plate a1(Re), a3(Re) so far (0.5-deg grid at every Re)
    sys.path.insert(0, str(HERE))
    from a3_re_scaling_airfoil import coeffs                        # noqa: E402
    zr = []
    for re, base in ZTP_RE.items():
        df = load("zero_thickness_plate/" + base, re)
        df = df[np.isclose((df.aoa * 2) % 1, 0)]
        n = len(df)
        if n < 13:
            zr.append(dict(Re=re, n=n)); continue
        r, _ = coeffs(df[["aoa", "cl"]].reset_index(drop=True))
        a, se, *_ = fit_window(df, 6, (1, 3, 5, 7, 9))      # P1-P9: the ZTP needs it (P1-P7 window-dependent)
        zr.append(dict(Re=re, n=n, a1=r["a1"], a3=r["a3"], a3_sd=r["a3_sd"], a3_01=r["a3_01"], a5_P9=a[5], se5=se[5],
                       R35_deg=np.degrees(np.sqrt(abs(r["a3"] / a[5])))))
    Z = pd.DataFrame(zr)
    Z.to_csv(OUT / "ztp_a3_vs_re.csv", index=False)
    lines.append("Zero-thickness plate, elliptic grid, a1 and a3 vs Re (complete Re only; 0.5-deg grid)")
    lines.append(Z.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    txt = "\n".join(lines)
    (OUT / "stress_test.txt").write_text(txt + "\n")
    print(txt)

    # ---- figures
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for name, (base, col) in BODIES.items():
        s = C[(C.body == name) & (C.order == 7)]
        ax[0].errorbar(s.A, -s.a3, yerr=s.se3, fmt="o-", ms=4, color=col, capsize=2, label=name)
        ax[1].errorbar(s.A, -s.a5, yerr=s.se5, fmt="o-", ms=4, color=col, capsize=2, label=name)
        ax[2].plot(s.A, s.R35_deg, "o-", ms=4, color=col, label=name)
    ax[0].set(xlabel="fit window 0–A (deg)", ylabel=r"$-a_3$ (rad$^{-3}$)", title=r"$a_3$ is window-independent (P1–P7)")
    ax[1].set(xlabel="fit window 0–A (deg)", ylabel=r"$-a_5$ (rad$^{-5}$)", title=r"$a_5$: stable except NACA0012 and the plate past 6°")
    ax[1].set_yscale("symlog", linthresh=50)
    ax[2].set(xlabel="fit window 0–A (deg)", ylabel=r"$\sqrt{|a_3/a_5|}$ (deg)", title="two-term radius estimate (crude)")
    ax[2].set_ylim(0, 60)
    for a_ in ax:
        a_.grid(alpha=0.25); a_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "coefficients_vs_window.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(1, 4, figsize=(18, 4.2), sharey=True)
    styles = {(1, 3): ("-", r"$a_1\theta+a_3\theta^3$"), (1, 3, 5): ("--", r"$+a_5\theta^5$"),
              (1, 3, 5, 7): (":", r"$+a_7\theta^7$")}
    for a_, (name, (base, col)) in zip(ax, BODIES.items()):
        for degs, (ls, lab) in styles.items():
            x, e = ext[name][degs]
            a_.semilogy(x[x > 0], e[x > 0], ls, color=col, lw=1.3, label=lab)
        a_.axhline(1e-3, color="k", lw=0.7); a_.axvspan(0, 4, color="0.92", zorder=0)
        a_.set(title=name, xlabel="α (deg)"); a_.grid(alpha=0.25, which="both"); a_.legend(fontsize=7)
    ax[0].set_ylabel(r"$|C_L^{fit} - C_L|$  (fitted on 0–4°, grey)")
    fig.tight_layout(); fig.savefig(OUT / "extrapolation.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    Zc = Z.dropna(subset=["a3"])
    ax[0].errorbar(Zc.Re, -Zc.a3, yerr=Zc.a3_sd, fmt="o-", color="#c2410c", lw=2, label="zero-thickness plate")
    ax[1].plot(Zc.Re, Zc.a1, "o-", color="#c2410c", lw=2, label="zero-thickness plate")
    for name, (f, col) in OTHER_RE.items():
        s = pd.read_csv(f)
        s = s[(s.Re >= 200) & (s.Re <= 900)]
        ax[0].plot(s.Re, -s.a3, "s--", ms=4, color=col, label=name)
        ax[1].plot(s.Re, s.a1, "s--", ms=4, color=col, label=name)
    ax[0].set(xlabel="Re", ylabel=r"$|a_3|$ (rad$^{-3}$)", title=r"$|a_3|$ vs Re")
    ax[1].set(xlabel="Re", ylabel=r"$a_1$ (rad$^{-1}$)", title=r"$a_1$ vs Re")
    for a_ in ax:
        a_.grid(alpha=0.25); a_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "a3_a1_vs_re_four_bodies.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    main()
