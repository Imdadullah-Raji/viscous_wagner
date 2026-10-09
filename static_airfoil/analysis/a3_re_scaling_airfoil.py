"""Airfoil: Re dependence of the small-angle coefficients a1, a3 (instructions/Quadratic_hypothesis.md).

Same method as ../plate_case/analysis/a3_re_scaling.py, for a NACA section:
odd Legendre fits on windows 0-A (x = alpha/A, mirrored), A on the data's own grid up to 6 deg;

  a1      P1+P3 on the smallest window (0-1 deg);
  a3      mean of the P1-P7 a3 over windows A in [2.25, 4] deg, std as the error bar
          (the plate's plateau; at high Re a LE bubble approaches inside 0-6 deg);
  a3_01   P1+P3 on 0-1 deg, an independent check.

Most NACA0012 folders are on a 0.5-deg grid (Re 500 is 0.25), so the P1-P7 plateau rests on
windows with 5-8 positive angles.  --check-grid refits the 0.25-deg Re = 500 curve on its
0.5-deg subset to measure what the coarser grid costs.

  python3 analysis/a3_re_scaling_airfoil.py --shape NACA0012 [--check-grid]
  -> analysis/a3_re_scaling/<shape>/
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lasso_re1000 import collect
from legendre_odd_fit import ODD, fit, mirror, to_theta_monomial

HERE = Path(__file__).resolve().parent
PLATEAU = (2.25, 4.0)
AMAX = 6.0
RES = {"NACA0012": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000],
       "NACA0006": [200, 300, 400, 500, 600, 700, 800, 900]}


def load(shape, re):
    if shape == "NACA0012":
        df = collect(re)
    else:
        base = f"NACA0006/steady_re{re}" if re == 500 else f"NACA0006/steady_re_scaling/steady_re{re}"
        df = collect(re, base)
    bad = df[(df.aoa <= AMAX + 1e-9) & (df.verdict != "CONVERGED")]
    df = df[(df.aoa <= AMAX + 1e-9) & (df.verdict == "CONVERGED")]
    return df[["aoa", "cl"]].sort_values("aoa").reset_index(drop=True), sorted(bad.aoa)


def window(df, A):
    d = df[df.aoa <= A + 1e-9]
    x0, y0 = d.aoa.to_numpy() / A, d.cl.to_numpy()
    n = int((x0 > 0).sum())
    x, y = mirror(x0, y0)
    row = dict(A=A, n=n)
    for degs, tag in (((1, 3), "P13"), ((1, 3, 5), "P135"), ((1, 3, 5, 7), "P1357")):
        if len(degs) > n:
            continue
        cc, r, m = fit(x, y, degs)
        th = to_theta_monomial(cc, A, degs)
        row[f"a1_{tag}"], row[f"a3_{tag}"] = th[1], th[3]
        if 5 in degs:
            row[f"a5_{tag}"] = th[5]
    return row


def coeffs(df):
    grid = np.round(df.aoa.to_numpy(), 6)
    A_all = [a for a in grid if a >= 1 - 1e-9]
    T = pd.DataFrame([window(df, A) for A in A_all])
    plat = T[(T.A >= PLATEAU[0] - 1e-9) & (T.A <= PLATEAU[1] + 1e-9)]
    first = T[np.isclose(T.A, 1.0)].iloc[0]
    return dict(a1=first.a1_P13, a3=plat.a3_P1357.mean(), a3_sd=plat.a3_P1357.std(),
                a3_01=first.a3_P13, a5=plat.a5_P1357.mean(), n_plateau=len(plat),
                step=float(np.median(np.diff(grid)))), T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", default="NACA0012", choices=sorted(RES))
    ap.add_argument("--check-grid", action="store_true")
    args = ap.parse_args()
    out = HERE / "a3_re_scaling" / args.shape
    out.mkdir(parents=True, exist_ok=True)

    rows, tables, lines = [], {}, []
    for re in RES[args.shape]:
        df, bad = load(args.shape, re)
        r, T = coeffs(df)
        rows.append(dict(Re=re, n_angles=len(df), not_converged=" ".join(f"{a:g}" for a in bad), **r))
        tables[re] = T
    S = pd.DataFrame(rows)
    re = S.Re.to_numpy(float)
    S["local_exp_a3"] = np.gradient(np.log(-S.a3.to_numpy()), np.log(re))
    S["local_exp_a1"] = np.gradient(np.log(S.a1.to_numpy()), np.log(re))
    S.to_csv(out / "a3_vs_re.csv", index=False)
    p3 = np.polyfit(np.log(re), np.log(-S.a3), 1)
    p1 = np.polyfit(np.log(re), np.log(S.a1), 1)
    q3 = np.polyfit(np.log(re), np.log(-S.a3), 2)

    f = lambda v: f"{v:+.4f}"
    lines += [f"{args.shape}, small-angle coefficients vs Re (0-{AMAX:g} deg, odd Legendre on 0-A, mirrored)",
              f"a3 = mean P1-P7 over windows {PLATEAU[0]}-{PLATEAU[1]} deg (+- std); a1, a3_01 from P1+P3 on 0-1 deg",
              "", S.drop(columns=["a5"]).to_string(index=False, float_format=f), "",
              f"single power law |a3| ~ Re^{p3[0]:.4f};  a1 ~ Re^{p1[0]:.4f}",
              f"curvature of ln|a3| vs ln Re (quadratic coeff) {q3[0]:+.4f}", "",
              "P1-P7 a3 per window (rows A, columns Re):"]
    W = pd.DataFrame({r: tables[r].set_index("A").a3_P1357 for r in RES[args.shape]})
    lines.append(W.to_string(float_format=lambda v: f"{v:7.2f}"))

    if args.check_grid and args.shape == "NACA0012":
        df, _ = load("NACA0012", 500)
        full, _ = coeffs(df)
        half, _ = coeffs(df[np.isclose((df.aoa * 2) % 1, 0)])
        lines += ["", "Grid check, Re = 500: 0.25-deg data vs its 0.5-deg subset",
                  f"  a1 {full['a1']:.4f} vs {half['a1']:.4f};  a3 {full['a3']:.3f} +- {full['a3_sd']:.3f} "
                  f"vs {half['a3']:.3f} +- {half['a3_sd']:.3f};  a3_01 {full['a3_01']:.3f} vs {half['a3_01']:.3f}"]
    txt = "\n".join(lines)
    (out / "a3_vs_re.txt").write_text(txt + "\n")
    print(txt)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    cmap = plt.get_cmap("viridis")
    for i, r in enumerate(RES[args.shape]):
        T = tables[r]
        ax[0].plot(T.A, T.a3_P1357, "o-", ms=2.5, color=cmap(i / (len(RES[args.shape]) - 1)), label=f"Re {r}")
    ax[0].axvspan(*PLATEAU, color="#fde68a", alpha=0.35, zorder=0)
    ax[0].set(xlabel="window 0–A (deg)", ylabel=r"$a_3$ from P1–P7 (rad$^{-3}$)",
              title="P1–P7 a3 per window (yellow: plateau)")
    ax[0].legend(fontsize=6, ncol=2)
    ax[1].errorbar(re, -S.a3, yerr=S.a3_sd, fmt="o", color="#c2410c", label=r"$|a_3|$ (P1–P7 plateau)")
    ax[1].plot(re, -S.a3_01, "x", color="k", label=r"$|a_3|$, P1+P3 on 0–1°")
    ax[1].set(xscale="log", yscale="log", xlabel="Re", ylabel=r"$|a_3|$ (rad$^{-3}$)", title="cubic coefficient")
    ax[1].legend(fontsize=8)
    ax[2].plot(re, S.a1, "o", color="#2a78d6", label=r"$a_1$")
    rr = np.geomspace(re.min(), re.max(), 50)
    ax[2].plot(rr, np.exp(p1[1]) * rr ** p1[0], "--", color="#2a78d6", lw=1, label=rf"$\propto Re^{{{p1[0]:.3f}}}$")
    ax[2].set(xscale="log", xlabel="Re", ylabel=r"$a_1$ (rad$^{-1}$)", title="lift slope at zero incidence")
    ax[2].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=0.25, which="both")
    fig.suptitle(f"{args.shape}: small-angle coefficients vs Re")
    fig.tight_layout(); fig.savefig(out / "a3_vs_re.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    main()
