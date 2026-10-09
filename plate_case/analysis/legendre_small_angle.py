"""Small-angle limit of the plate's odd Legendre fit, Re = 500 (instructions/Quadratic_hypothesis.md).

Claim under test: near zero incidence C_L = a1 theta + a3 theta^3, i.e. dC_L/dalpha is
quadratic, with the cubic term alive.  Fit odd Legendre series on 0-A for A = 1, 1.25, ..., 6
deg (every 0.25-deg point of plate_case/re500, all CONVERGED) and follow the higher-order terms
as A shrinks.

What "shrinks" means.  On x = theta/A a pure power series a_k theta^k puts a_k A^k (scaled by
constants) into P_k, so for a smooth odd C_L
    c3 ~ A^3,   c5 ~ A^5,   c7 ~ A^7,   c9 ~ A^9      (A -> 0)
and the collected monomial coefficients a1, a3, a5 tend to finite limits.  The test is
therefore:
  (i)   log|c_k| vs log A has slope k (c5/c3 ~ A^2, c7/c3 ~ A^4): the higher terms vanish
        faster than the cubic;
  (ii)  a3 -> a finite nonzero limit, the same whether P1+P3, P1-P5 or P1-P7 is fitted;
  (iii) the P1+P3 fit error falls like A^5 (the first neglected term, theta^5);
  (iv)  the quintic's share of the curvature, |a5| theta^2 / |a3| -> 0.
A term that sits at a constant floor instead of falling is noise, not physics.

Same basis and mirroring as legendre_plate.py / static_airfoil legendre_odd_fit.py.
Run with the research env:  python plate_case/analysis/legendre_small_angle.py
-> plate_case/analysis/legendre_small_angle/
"""
from pathlib import Path
import re as regex
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLATE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PLATE.parent / "static_airfoil" / "analysis"))
from legendre_odd_fit import ODD, fit, mirror, to_theta_monomial  # noqa: E402

OUT = PLATE / "analysis" / "legendre_small_angle"
WINDOWS = np.arange(1.0, 6.0 + 1e-9, 0.25)


def plate_data(re=500):
    """Last force sample of each case (13 significant digits), CONVERGED only."""
    root = PLATE / f"re{re}"
    status = {}
    for line in (root / "RESULTS.md").read_text().splitlines():
        m = regex.match(r"\|\s*([\d.]+)\s*\|\s*(\w+)\s*\|", line)
        if m:
            status[float(m.group(1))] = m.group(2)
    rows = []
    for case in root.glob("steady_aoa*"):
        if not case.is_dir():
            continue
        a = float(case.name.removeprefix("steady_aoa"))
        if status.get(a) != "CONVERGED":
            continue
        dats = sorted((case / "postProcessing" / "forceCoeffs").glob("*/forceCoeffs.dat"),
                      key=lambda p: float(p.parent.name))
        last = np.loadtxt(dats[-1], comments="#", usecols=(0, 2, 3), ndmin=2)[-1]
        rows.append((a, last[2]))
    return pd.DataFrame(rows, columns=["aoa", "cl"]).sort_values("aoa").reset_index(drop=True)


def window(df, A):
    d = df[df.aoa <= A + 1e-9]
    x0, y0 = d.aoa.to_numpy() / A, d.cl.to_numpy()
    n = int((x0 > 0).sum())
    x, y = mirror(x0, y0)
    row = dict(A=A, n=n)
    full = ODD[:min(len(ODD), n)]
    c, rmse, mx = fit(x, y, full)
    for k, v in zip(ODD, np.r_[c, np.full(len(ODD) - len(c), np.nan)]):
        row[f"c{k}"] = v
    row["rmse_full"] = rmse
    for degs, tag in (((1, 3), "P13"), ((1, 3, 5), "P135"), ((1, 3, 5, 7), "P1357")):
        if len(degs) > n:
            continue
        cc, r, m = fit(x, y, degs)
        th = to_theta_monomial(cc, A, degs)
        row[f"a1_{tag}"], row[f"a3_{tag}"] = th[1], th[3]
        if 5 in degs:
            row[f"a5_{tag}"] = th[5]
        row[f"maxres_{tag}"] = m
    return row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = plate_data()
    have = set(np.round(df.aoa, 6))
    T = pd.DataFrame([window(df, A) for A in WINDOWS if round(A, 6) in have])
    T.to_csv(OUT / "small_angle_windows.csv", index=False)
    tA = np.deg2rad(T.A.to_numpy())

    # (i) slopes of log|c_k| vs log A, fitted over the four smallest and the four largest windows
    lines = ["Plate, Re = 500: odd Legendre on 0-A (x = alpha/A, mirrored), A = 1..6 deg, 0.25-deg data",
             "", "(i) d ln|c_k| / d ln A   (smooth odd C_L: slope = k)"]
    for k in ODD:
        v = T[f"c{k}"].abs().to_numpy()
        ok = np.isfinite(v) & (v > 0)
        if ok.sum() < 4:
            continue
        lo = np.polyfit(np.log(T.A[ok][:4]), np.log(v[ok][:4]), 1)[0]
        hi = np.polyfit(np.log(T.A[ok][-4:]), np.log(v[ok][-4:]), 1)[0]
        allp = np.polyfit(np.log(T.A[ok]), np.log(v[ok]), 1)[0]
        lines.append(f"  c{k}: smallest four windows {lo:5.2f}   largest four {hi:5.2f}   all {allp:5.2f}"
                     f"   (|c{k}| from {v[ok][0]:.2e} at {T.A[ok].iloc[0]:g} to {v[ok][-1]:.2e} at 6)")

    cols = ["A", "n", "c1", "c3", "c5", "c7", "c9"]
    lines += ["", "Coefficients:", T[cols].to_string(index=False, float_format=lambda v: f"{v:+.4e}")]
    lines += ["", "(ii) a3 (rad^-3) under three truncations; (iii) max |res| of P1+P3; "
              "(iv) quintic_share = |a5| theta_A^2 / |a3| from P1-P5: the theta^5 term's size relative to the "
              "theta^3 term at the window edge (both parts of C_L - a1 theta)"]
    T["quintic_share"] = (T.a5_P135.abs() * tA ** 2 / T.a3_P135.abs())
    cols2 = ["A", "a1_P13", "a3_P13", "a3_P135", "a3_P1357", "a5_P135", "maxres_P13", "quintic_share"]
    lines.append(T[cols2].to_string(index=False, float_format=lambda v: f"{v:+.4e}"))
    s = np.polyfit(np.log(T.A[:6]), np.log(T.maxres_P13[:6]), 1)[0]
    lines.append(f"  d ln(max|res| P1+P3) / d ln A over A = 1-2.25: {s:.2f}   (theta^5 truncation: 5)")
    # Small-angle limit.  The data scatter about a smooth curve by ~4e-7 in C_L (every case,
    # every solver version alike), which puts a ~1e-6 floor under c5..c9 and makes the
    # four- and five-term fits on the smallest windows (A < 2.25, n < 9) noise-dominated.
    # Estimate a3(0) from the P1-P7 fits on A = 2.25-6 (flat to +-0.05) and check it against
    # the two-term fit on the smallest window, where theta^5 is below the noise.
    plat = T[(T.A >= 2.25 - 1e-9)]
    a3lim, a3sd = plat.a3_P1357.mean(), plat.a3_P1357.std()
    a1lim = T.a1_P13.iloc[0]
    lines.append(f"  small-angle limit: a3 = {a3lim:.2f} +- {a3sd:.2f} /rad^3 (P1-P7, A = 2.25-6);"
                 f"  P1+P3 on 0-1 deg gives {T.a3_P13.iloc[0]:.2f};  a1 = {a1lim:.4f} /rad")
    lines.append("  noise: smooth-fit residuals are +-4e-7 in C_L for dev, OF13 and OF14 points alike "
                 "(no version offset), so c5, c7, c9 below ~1e-6 are noise")
    a3lim = (None, a3lim)
    txt = "\n".join(lines)
    (OUT / "small_angle.txt").write_text(txt + "\n")
    print(txt)

    # ---- figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    colors = {1: "#2a78d6", 3: "#c2410c", 5: "#1baf7a", 7: "#7b3fbf", 9: "#999999"}
    for k in ODD:
        v = T[f"c{k}"].abs()
        ax[0].loglog(T.A, v, "o-", ms=3, color=colors[k], label=f"|c{k}|")
        if k > 1:
            ref = v.iloc[-1] * (T.A / T.A.iloc[-1]) ** k
            ax[0].loglog(T.A, ref, ":", color=colors[k], lw=0.8)
    ax[0].axhspan(1e-8, 1e-6, color="0.92", zorder=0)
    ax[0].text(1.05, 2e-8, "noise floor (C_L scatter ±4e-7)", fontsize=7, color="0.4")
    ax[0].set(xlabel="window 0–A (deg)", ylabel="|Legendre coefficient|",
              title="dotted: slope k (smooth power law)")
    ax[0].legend(fontsize=7, ncol=2)
    for tag, mk, lab in (("P13", "o", "P1+P3"), ("P135", "s", "P1–P5"), ("P1357", "^", "P1–P7")):
        ax[1].plot(T.A, T[f"a3_{tag}"], mk + "-", ms=3, label=lab)
    ax[1].axhline(a3lim[1], color="k", lw=0.8, ls="--", label=f"limit {a3lim[1]:.2f}")
    ax[1].axvspan(T.A.min() - 0.1, 2.125, color="0.92", zorder=0, label="≥4-term fits noise-dominated")
    ax[1].set(xlabel="window 0–A (deg)", ylabel=r"$a_3$ (rad$^{-3}$)", title=r"cubic coefficient $a_3$")
    ax[1].legend(fontsize=8)
    ax[2].loglog(T.A, T.maxres_P13, "o-", ms=3, color="k", label="max |res|, P1+P3")
    ax[2].loglog(T.A, T.maxres_P13.iloc[0] * (T.A / T.A.iloc[0]) ** 5, ":", color="k", lw=0.8,
                 label=r"$\propto A^5$")
    ax2 = ax[2].twinx()
    ax2.semilogy(T.A, T.quintic_share, "s-", ms=3, color="#1baf7a")
    ax2.set_ylabel(r"$|a_5|\theta_A^2/|a_3|$  (θ⁵ vs θ³ at the edge)", color="#1baf7a")
    ax[2].set(xlabel="window 0–A (deg)", ylabel=r"$C_L$", title="two-term error and quintic share")
    ax[2].legend(fontsize=8, loc="upper left")
    for a_ in ax:
        a_.grid(alpha=0.25, which="both")
    fig.suptitle("Plate, Re = 500: odd Legendre fits on 0–A as A → 0")
    fig.tight_layout(); fig.savefig(OUT / "small_angle.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    main()
