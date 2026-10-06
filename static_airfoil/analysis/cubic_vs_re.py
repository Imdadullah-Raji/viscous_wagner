"""OLS odd-cubic coefficients b1, b3 of the steady lift curve versus Re.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/cubic_vs_re.py

C_L = b0 + b1*theta + b3*theta^3 (theta in rad) on the converged steady cases of
steady_re{100,500,1000,2000}/ over 0-6 and 0-8 deg, then a power law
|b3| = c * Re^p fitted in log space.  OLS, not Lasso: at fixed lambda the
Lasso selection itself depends on the fit range (see STEADY_LASSO_ReSweep.md),
while b3 is stable between ranges.
-> analysis/figs/resweep_cubic_vs_re/   (defaults)

  --re 100 200 ... 2000   Re list (cases found by lasso_re1000.collect, which
                          also looks in steady_re_scaling/)
  --step 1                use only angles on this grid (deg), so every Re is fitted
                          on the same points; default 0.5
  non-default choices write to figs/cubic_vs_re_<nRe>re_step<step>/
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lasso_re1000 import collect, ROOT

RES = (100, 500, 1000, 2000)
STEP = 0.5
RANGES = {6.0: "#2a78d6", 8.0: "#eb6834"}     # categorical slots 1, 2
OUT = ROOT / "analysis" / "figs" / "resweep_cubic_vs_re"


def cubic(df, amax):
    d = df[(df.aoa <= amax + 1e-9) & np.isclose((df.aoa / STEP) % 1, 0)
           & (df.verdict == "CONVERGED")]
    th = np.deg2rad(d.aoa.to_numpy())
    A = np.column_stack((np.ones_like(th), th, th ** 3))
    b, *_ = np.linalg.lstsq(A, d.cl.to_numpy(), rcond=None)
    res = d.cl.to_numpy() - A @ b
    return dict(n=len(d), b0=b[0], b1=b[1], b3=b[2],
                rmse=np.sqrt(np.mean(res ** 2)), maxres=np.abs(res).max())


def main():
    global RES, STEP, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", type=int, nargs="+", default=list(RES))
    ap.add_argument("--step", type=float, default=STEP)
    args = ap.parse_args()
    if tuple(args.re) != RES or args.step != STEP:
        OUT = ROOT / "analysis" / "figs" / f"cubic_vs_re_{len(args.re)}re_step{args.step:g}"
    RES, STEP = tuple(args.re), args.step
    OUT.mkdir(parents=True, exist_ok=True)
    dfs = {re: collect(re) for re in RES}
    rows = [dict(re=re, amax=a, **cubic(dfs[re], a)) for a in RANGES for re in RES]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "cubic_coeffs_vs_re.csv", index=False)
    lines = [t.to_string(index=False, float_format=lambda v: f"{v:.6g}")]
    fits = {}
    for a in RANGES:
        s = t[t.amax == a]
        p, lc = np.polyfit(np.log(s.re), np.log(-s.b3), 1)
        fits[a] = (np.exp(lc), p)
        lines.append(f"0-{a:g} deg: |b3| = {np.exp(lc):.4g} * Re^{p:.4f}")
    txt = "\n".join(lines)
    (OUT / "cubic_vs_re.txt").write_text(txt + "\n")
    print(txt)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    g = np.logspace(np.log10(80), np.log10(2500), 100)
    for a, c in RANGES.items():
        s = t[t.amax == a]
        a1.plot(s.re, s.b1, "o-", color=c, lw=2, ms=6, label=f"fit 0-{a:g}°")
        a2.plot(s.re, -s.b3, "o", color=c, ms=7, label=f"fit 0-{a:g}°")
        k, p = fits[a]
        a2.plot(g, k * g ** p, "--", color=c, lw=1.2,
                label=rf"$|b_3| = {k:.3g}\,Re^{{{p:.3f}}}$")
    a1.axhline(np.pi, color="0.6", lw=1, ls=":")
    a1.text(110, np.pi + 0.01, r"$\pi$ (half of thin-airfoil $2\pi$)", color="0.35", fontsize=9)
    a1.set(xscale="log", xlabel="Re", ylabel=r"$b_1$ (per rad)")
    a1.set_title(r"Linear coefficient $b_1$", loc="left")
    a2.set(xscale="log", yscale="log", xlabel="Re", ylabel=r"$|b_3|$")
    a2.set_title(r"Cubic coefficient $|b_3|$, log-log", loc="left")
    for ax in (a1, a2):
        ax.set_xticks(RES, [str(r) for r in RES])
        ax.grid(alpha=0.25, lw=0.6, which="both")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "cubic_coeffs_vs_re.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
