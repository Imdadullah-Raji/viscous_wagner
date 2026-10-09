"""Plate: Re dependence of the small-angle coefficients a1, a3 (instructions/Quadratic_hypothesis.md).

Claim under test: C_L = a1 theta + a3 theta^3 near zero incidence, with a3 following a (local)
scaling law in Re.  For each Re in 200..900 (0-6 deg at 0.25 deg, all CONVERGED) repeat the
Re = 500 small-angle analysis of legendre_small_angle.py (same windows, same fits) and collect

  a1      P1+P3 on the smallest window, 0-1 deg (theta^5 below the noise there);
  a3      mean of the P1-P7 fits over windows A = 2.25-4 deg; its std there is the error bar.
          The plateau stops at 4 deg, not 6 as at Re = 500 alone: at high Re a leading-edge
          bubble forms early (upper-surface C_f,min crosses zero near 5.75 deg at Re = 900,
          extrapolates to ~6.6 at Re = 700, 7.9 at Re = 500), and from ~4.5 deg the P1-P7 a3
          drifts as the lift curve picks up that run-up (Re = 900: -15.9 flat on 2.25-4,
          -18.6 at 6).  Below 2.25 deg the 4-term fits are noise-dominated (n < 9 angles);
  a3_01   P1+P3 on 0-1 deg, an independent check of the limit;
  a5      mean of the P1-P7 a5 over the same windows;
  slope_k d ln|c_k| / d ln A on the largest four windows (k for a smooth power series).

Then fit ln|a3| and ln a1 against ln Re: a global power law, and the local exponent
d ln|a3| / d ln Re from centred differences.

Run: python plate_case/analysis/a3_re_scaling.py  -> plate_case/analysis/a3_re_scaling/
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from legendre_small_angle import WINDOWS, plate_data, window

OUT = Path(__file__).resolve().parent / "a3_re_scaling"
RES = [200, 300, 400, 500, 600, 700, 800, 900]
PLATEAU = (2.25, 4.0)


def per_re(re):
    df = plate_data(re)
    df = df[df.aoa <= 6 + 1e-9]
    have = set(np.round(df.aoa, 6))
    T = pd.DataFrame([window(df, A) for A in WINDOWS if round(A, 6) in have])
    plat = T[(T.A >= PLATEAU[0] - 1e-9) & (T.A <= PLATEAU[1] + 1e-9)]
    row = dict(Re=re, n_angles=len(df),
               a1=T.a1_P13.iloc[0], a3=plat.a3_P1357.mean(), a3_sd=plat.a3_P1357.std(),
               a3_01=T.a3_P13.iloc[0], a3_P135_6=T.a3_P135.iloc[-1], a3_P13_6=T.a3_P13.iloc[-1],
               a5=plat.a5_P1357.mean(), a5_sd=plat.a5_P1357.std(),
               maxres_P1357_6=T.maxres_P1357.iloc[-1])
    for k in (3, 5, 7):
        v = T[f"c{k}"].abs().to_numpy()
        row[f"slope_c{k}"] = np.polyfit(np.log(T.A[-4:]), np.log(v[-4:]), 1)[0]
    return row, T


def powerlaw(x, y):
    p, cov = np.polyfit(np.log(x), np.log(y), 1, cov=True)
    return p[0], np.sqrt(cov[0, 0]), np.exp(p[1])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, tables = [], {}
    for re in RES:
        r, T = per_re(re)
        rows.append(r); tables[re] = T
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "a3_vs_re.csv", index=False)
    re = S.Re.to_numpy(float)

    m3, m3e, k3 = powerlaw(re, -S.a3)
    m1, m1e, k1 = powerlaw(re, S.a1)
    m5, m5e, k5 = powerlaw(re, -S.a5)
    lna3 = np.log(-S.a3.to_numpy())
    local = np.gradient(lna3, np.log(re))          # centred inside, one-sided at the ends
    S["local_exp_a3"] = local
    S["local_exp_a1"] = np.gradient(np.log(S.a1.to_numpy()), np.log(re))
    resid = lna3 - (np.log(k3) + m3 * np.log(re))
    # curvature of ln|a3| vs ln Re: quadratic term of a 2nd-order fit
    q = np.polyfit(np.log(re), lna3, 2)

    fmt = lambda v: f"{v:+.4f}"
    lines = ["Plate, small-angle coefficients vs Re (0-6 deg at 0.25 deg, odd Legendre, x = alpha/A mirrored)",
             f"a3 = mean P1-P7 over windows {PLATEAU[0]}-{PLATEAU[1]} deg (+- std); a1, a3_01 from P1+P3 on 0-1 deg", "",
             S[["Re", "n_angles", "a1", "a3", "a3_sd", "a3_01", "a5", "a5_sd", "a3_P13_6",
                "local_exp_a3", "local_exp_a1"]].to_string(index=False, float_format=fmt), "",
             "Window-scaling check (largest four windows; a smooth power series gives slope k):",
             S[["Re", "slope_c3", "slope_c5", "slope_c7"]].to_string(index=False, float_format=fmt), "",
             f"Global power laws over Re = {RES[0]}-{RES[-1]}:",
             f"  |a3| = {k3:.4g} Re^{m3:.4f} (+- {m3e:.4f});  max |ln residual| {np.abs(resid).max():.2e}",
             f"  a1   = {k1:.4g} Re^{m1:.4f} (+- {m1e:.4f})",
             f"  |a5| = {k5:.4g} Re^{m5:.4f} (+- {m5e:.4f})",
             f"  curvature of ln|a3| vs ln Re (quadratic coeff): {q[0]:+.4f}",
             f"  local exponent d ln|a3|/d ln Re: {local.min():+.3f} .. {local.max():+.3f}"]
    txt = "\n".join(lines)
    (OUT / "a3_vs_re.txt").write_text(txt + "\n")
    print(txt)

    # ---- figure
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    cmap = plt.get_cmap("viridis")
    for i, r in enumerate(RES):
        T = tables[r]
        ax[0, 0].plot(T.A, T.a3_P1357, "o-", ms=2.5, color=cmap(i / (len(RES) - 1)), label=f"Re {r}")
    ax[0, 0].axvspan(0.9, PLATEAU[0] - 0.125, color="0.92", zorder=0)
    ax[0, 0].axvspan(*PLATEAU, color="#fde68a", alpha=0.35, zorder=0)
    ax[0, 0].set(xlabel="window 0–A (deg)", ylabel=r"$a_3$ from P1–P7 (rad$^{-3}$)",
                 title="P1–P7 a3 per window (yellow: plateau used; grey: noise)")
    ax[0, 0].legend(fontsize=7, ncol=2)
    ax[0, 1].errorbar(re, -S.a3, yerr=S.a3_sd, fmt="o", color="#c2410c", label=r"$|a_3|$ (P1–P7 plateau)")
    ax[0, 1].plot(re, -S.a3_01, "x", color="k", label=r"$|a_3|$, P1+P3 on 0–1°")
    rr = np.geomspace(re.min(), re.max(), 50)
    ax[0, 1].plot(rr, k3 * rr ** m3, "--", color="#c2410c", lw=1, label=rf"single power law $Re^{{{m3:.3f}}}$ (rejected)")
    ax[0, 1].set(xscale="log", yscale="log", xlabel="Re", ylabel=r"$|a_3|$ (rad$^{-3}$)", title="cubic coefficient")
    ax[0, 1].legend(fontsize=8)
    ax[1, 0].plot(re, S.a1, "o", color="#2a78d6", label=r"$a_1$")
    ax[1, 0].plot(rr, k1 * rr ** m1, "--", color="#2a78d6", lw=1, label=rf"$\propto Re^{{{m1:.3f}}}$")
    ax[1, 0].axhline(2 * np.pi, color="0.5", lw=0.8, ls=":", label=r"$2\pi$")
    ax[1, 0].set(xscale="log", xlabel="Re", ylabel=r"$a_1$ (rad$^{-1}$)", title="lift slope at zero incidence")
    ax[1, 0].legend(fontsize=8)
    ax[1, 1].plot(re, S.local_exp_a3, "o-", color="#c2410c", label=r"$d\ln|a_3|/d\ln Re$")
    ax[1, 1].plot(re, S.local_exp_a1, "s-", color="#2a78d6", label=r"$d\ln a_1/d\ln Re$")
    ax[1, 1].axhline(m3, color="#c2410c", lw=0.8, ls="--")
    ax[1, 1].set(xscale="log", xlabel="Re", ylabel="local exponent", title="local scaling exponents")
    ax[1, 1].legend(fontsize=8)
    for a_ in ax.flat:
        a_.grid(alpha=0.25, which="both")
    fig.suptitle("Plate: small-angle coefficients of $C_L = a_1\\theta + a_3\\theta^3 + \\dots$ vs Re")
    fig.tight_layout(); fig.savefig(OUT / "a3_vs_re.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    main()
