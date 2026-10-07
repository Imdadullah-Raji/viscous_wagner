"""Re scaling of the odd-Legendre cubic coefficient c3 of the steady lift curve.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/legendre_vs_re.py \
              [--re 100 200 ... 2000] [--step 1]

Same fit as legendre_odd_fit.py (x = theta/theta_max mirrored to [-1, 1], least
squares on odd P_k, no intercept) for every Re, on 0-6 and 0-8 deg, using only
angles on a common `--step` grid so every Re sees the same x points.  Reported:

  c1, c3, c5     from the full P1-P9 fit
  c3 @ P1+P3     and c3 @ P1-P5: if the basis is doing its job c3 barely moves
                 with truncation (raw theta^3 does not have this property)
  |c3| = k Re^p  in log space, plus the LOCAL slope between neighbouring Re

c3 is a projection on x in [-1, 1], so it depends on theta_max: compare Re at one
range, not ranges with each other.  For a pure cubic C_L = b1 theta + b3 theta^3,
c3 = (2/5) b3 theta_max^3, so the Re exponent of c3 and of b3 coincide.
-> analysis/figs/legendre_vs_re_<nRe>re_step<step>/
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lasso_re1000 import collect, ROOT
from legendre_odd_fit import ODD, basis, fit, mirror

RANGES = {6.0: "#2a78d6", 8.0: "#eb6834"}     # categorical slots 1, 2


def data(df, amax, step):
    d = df[(df.aoa <= amax + 1e-9) & np.isclose((df.aoa / step) % 1, 0)
           & (df.verdict == "CONVERGED")]
    return d.aoa.to_numpy() / amax, d.cl.to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", type=int, nargs="+",
                    default=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000])
    ap.add_argument("--step", type=float, default=1.0)
    args = ap.parse_args()
    out = ROOT / "analysis" / "figs" / f"legendre_vs_re_{len(args.re)}re_step{args.step:g}"
    out.mkdir(parents=True, exist_ok=True)

    dfs = {re: collect(re) for re in args.re}
    rows = []
    for amax in RANGES:
        for re in args.re:
            x0, y0 = data(dfs[re], amax, args.step)
            x, y = mirror(x0, y0)
            c, rmse, mx = fit(x, y)
            c3t = {k: fit(x, y, tuple(d for d in ODD if d <= k))[0][1] for k in (3, 5)}
            rows.append(dict(amax=amax, re=re, n_half=len(x0), c1=c[0], c3=c[1], c5=c[2],
                             c7=c[3], c9=c[4], c3_P3=c3t[3], c3_P5=c3t[5],
                             ratio_c3_c1=c[1] / c[0], rmse=rmse, max_abs_res=mx,
                             cond=np.linalg.cond(basis(x))))
    t = pd.DataFrame(rows)
    t.to_csv(out / "legendre_vs_re.csv", index=False)

    fmt = lambda v: f"{v:+.5e}"
    lines = [f"grid step {args.step:g} deg; x = theta/theta_max mirrored to [-1, 1]; P1-P9 fit\n",
             t.to_string(index=False, float_format=fmt), ""]
    fits, slopes = {}, {}
    for amax in RANGES:
        s = t[t.amax == amax].sort_values("re")
        re, a = s.re.to_numpy(), -s.c3.to_numpy()
        p, lc = np.polyfit(np.log(re), np.log(a), 1)
        fits[amax] = (np.exp(lc), p)
        slopes[amax] = (np.sqrt(re[:-1] * re[1:]), np.diff(np.log(a)) / np.diff(np.log(re)))
        lines.append(f"0-{amax:g} deg: |c3| = {np.exp(lc):.4g} * Re^{p:.4f}  (all Re)")
        for lo, hi in ((200, 1000),):
            m = (re >= lo) & (re <= hi)
            q, lq = np.polyfit(np.log(re[m]), np.log(a[m]), 1)
            r = np.log(a[m]) - (q * np.log(re[m]) + lq)
            lines.append(f"           Re {lo}-{hi}: p = {q:.4f}, max |log resid| {np.abs(r).max():.4f}")
        lines.append("           local slopes: " + "  ".join(
            f"{r0}-{r1}:{v:.3f}" for r0, r1, v in zip(re[:-1], re[1:], slopes[amax][1])))
        dev = np.max(np.abs(np.r_[s.c3_P3 / s.c3, s.c3_P5 / s.c3] - 1))
        lines.append(f"           c3 change with truncation (P3, P5 vs P9): max {100 * dev:.2f}%")
    txt = "\n".join(lines)
    (out / "legendre_vs_re.txt").write_text(txt + "\n")
    print(txt)

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 4.6))
    g = np.logspace(np.log10(80), np.log10(2500), 100)
    for amax, col in RANGES.items():
        s = t[t.amax == amax].sort_values("re")
        k, p = fits[amax]
        a1.plot(s.re, -s.c3, "o", color=col, ms=7, label=f"P1-P9 fit, 0-{amax:g}°")
        a1.plot(s.re, -s.c3_P3, "x", color=col, ms=7, mew=1.5, label=f"P1+P3 only, 0-{amax:g}°")
        a1.plot(g, k * g ** p, "--", color=col, lw=1.2, label=rf"${k:.3g}\,Re^{{{p:.3f}}}$")
        a2.plot(*slopes[amax], "o-", color=col, lw=2, ms=6, label=f"0-{amax:g}°")
        a3.plot(s.re, -100 * s.ratio_c3_c1, "o-", color=col, lw=2, ms=6, label=f"0-{amax:g}°")
    a1.set(xscale="log", yscale="log", xlabel="Re", ylabel=r"$|c_3|$")
    a1.set_title(r"Legendre cubic coefficient $|c_3|$", loc="left")
    a2.set(xscale="log", xlabel="Re (geometric midpoint)", ylabel=r"$d\ln|c_3|\,/\,d\ln Re$")
    a2.set_title("Local power-law exponent", loc="left")
    a3.set(xscale="log", xlabel="Re", ylabel=r"$|c_3/c_1|$ (%)")
    a3.set_title("Cubic bend relative to the linear part", loc="left")
    ticks = [100, 200, 300, 500, 700, 1000, 2000]
    for ax in (a1, a2, a3):
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_xticks(ticks, [str(v) for v in ticks])
        ax.grid(alpha=0.25, lw=0.6, which="both")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "legendre_c3_vs_re.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
