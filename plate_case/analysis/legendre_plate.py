"""Odd-Legendre fit of the steady lift curve of the thick flat plate, Re = 500 and 1000,
compared with the NACA0012 fitted identically.

Run with: /home/raji/miniconda3/envs/research/bin/python plate_case/analysis/legendre_plate.py

Same method as static_airfoil/analysis/legendre_odd_fit.py: x = theta/theta_max in [0, 1],
mirrored to [-1, 1] by symmetry (C_L(-x) = -C_L(x); the plate is symmetric too), least
squares on P1, P3, P5, P7, P9 with no intercept.  Fitted over 0-4, 0-6 and 0-8 deg.

Data: last sample of the latest restart chunk of each case's
postProcessing/forceCoeffs/<start>/forceCoeffs.dat (13 significant digits, so no
lift_coeffs.csv rounding), and only angles whose status in re<RE>/RESULTS.md is CONVERGED
(the plate sweep's own range+drift gate, run_sweep.py).  The airfoil side reuses
static_airfoil's collect() and its force-trace verdict.

Reported per (shape, Re, range): c1..c9, |c3/c1|, c3 under P3/P5 truncation, RMSE,
condition number, and the P1+P3 fit collected into theta: C_L = a1 theta + a3 theta^3.
-> plate_case/analysis/legendre/
"""
from pathlib import Path
import re as regex
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLATE = Path(__file__).resolve().parents[1]
AIRFOIL = PLATE.parent / "static_airfoil"
sys.path.insert(0, str(AIRFOIL / "analysis"))
from lasso_re1000 import collect                              # noqa: E402
from legendre_odd_fit import ODD, basis, fit, mirror, to_theta_monomial  # noqa: E402

OUT = PLATE / "analysis" / "legendre"
RES = (500, 1000)
RANGES = (4.0, 6.0, 8.0)
COLOR = {500: "#2a78d6", 1000: "#eb6834"}     # static_airfoil Re colours


def plate_data(re):
    status = {}
    for line in (PLATE / f"re{re}" / "RESULTS.md").read_text().splitlines():
        m = regex.match(r"\|\s*([\d.]+)\s*\|\s*(\w+)\s*\|", line)
        if m:
            status[float(m.group(1))] = m.group(2)
    rows = []
    for case in (PLATE / f"re{re}").glob("steady_aoa*"):
        if not case.is_dir():
            continue
        a = float(case.name.removeprefix("steady_aoa"))
        chunks = sorted((case / "postProcessing" / "forceCoeffs").iterdir(),
                        key=lambda p: float(p.name))
        last = np.loadtxt(chunks[-1] / "forceCoeffs.dat", comments="#", usecols=(0, 2, 3), ndmin=2)[-1]
        rows.append((a, last[2], last[1], status.get(a, "MISSING"), int(last[0])))
    return pd.DataFrame(rows, columns=["aoa", "cl", "cd", "verdict", "iters"]
                        ).sort_values("aoa").reset_index(drop=True)


def select(df, amax):
    d = df[(df.aoa <= amax + 1e-9) & np.isclose((df.aoa * 2) % 1, 0)]
    bad = d[d.verdict != "CONVERGED"]
    assert bad.empty, f"unconverged angle inside 0-{amax:g}: {bad.aoa.tolist()}"
    return d.aoa.to_numpy() / amax, d.cl.to_numpy()


def analyse(shape, re, df, amax):
    x0, y0 = select(df, amax)
    x, y = mirror(x0, y0)
    c, rmse, mx = fit(x, y)
    c_half, *_ = fit(x0, y0)
    assert np.allclose(c, c_half, rtol=1e-10, atol=1e-13)
    c3t = {k: fit(x, y, tuple(d for d in ODD if d <= k))[0][1] for k in (3, 5)}
    p13, r13, m13 = fit(x, y, (1, 3))
    th = to_theta_monomial(p13, amax, (1, 3))
    return dict(shape=shape, re=re, amax=amax, n_half=len(x0),
                **{f"c{k}": v for k, v in zip(ODD, c)},
                c3_over_c1=c[1] / c[0], c3_P3=c3t[3], c3_P5=c3t[5],
                rmse=rmse, max_abs_res=mx, cond=np.linalg.cond(basis(x)),
                a1_P13=th[1], a3_P13=th[3], rmse_P13=r13, max_abs_res_P13=m13)


def style(ax):
    ax.grid(alpha=0.25, lw=0.6, which="both")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = {("plate", re): plate_data(re) for re in RES}
    data |= {("NACA0012", re): collect(re) for re in RES}
    for (shape, re), df in data.items():
        df.to_csv(OUT / f"cl_{shape}_re{re}.csv", index=False)

    t = pd.DataFrame([analyse(s, re, data[(s, re)], a)
                      for a in RANGES for re in RES for s in ("plate", "NACA0012")])
    t.to_csv(OUT / "legendre_plate_vs_naca.csv", index=False)
    pd.set_option("display.width", 250)
    fmt = lambda v: f"{v:+.5e}"
    cols1 = ["shape", "re", "amax", "n_half", "c1", "c3", "c5", "c7", "c9", "c3_over_c1",
             "rmse", "max_abs_res", "cond"]
    cols2 = ["shape", "re", "amax", "c3", "c3_P3", "c3_P5", "a1_P13", "a3_P13",
             "rmse_P13", "max_abs_res_P13"]
    txt = ("x = theta/theta_max mirrored to [-1, 1]; least squares on P1..P9, no intercept\n\n"
           + t[cols1].to_string(index=False, float_format=fmt)
           + "\n\nc3 under truncation, and the P1+P3 fit collected into theta "
             "(C_L = a1 theta + a3 theta^3, rad):\n\n"
           + t[cols2].to_string(index=False, float_format=fmt))
    (OUT / "legendre_plate.txt").write_text(txt + "\n")
    print(txt)

    # Fig 1: lift curves, converged points, P1-P9 fit on 0-8 deg (plate Re 1000 peaks at 7.5)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, re in zip(axes, RES):
        for shape, mk in (("plate", "o"), ("NACA0012", "s")):
            df = data[(shape, re)]
            ok = df.verdict == "CONVERGED"
            ax.plot(df.aoa[ok], df.cl[ok], mk, color=COLOR[re], ms=5,
                    mfc=COLOR[re] if shape == "plate" else "none", label=f"{shape}, converged")
            if (~ok).any():
                ax.plot(df.aoa[~ok], df.cl[~ok], "x", color="0.55", ms=6,
                        label=f"{shape}, NOT converged")
            r = t[(t["shape"] == shape) & (t.re == re) & (t.amax == 8.0)].iloc[0]
            g = np.linspace(0, 8, 200)
            ax.plot(g, basis(g / 8) @ r[[f"c{k}" for k in ODD]].to_numpy(float), "-",
                    color=COLOR[re], lw=1, alpha=0.8)
        ax.axvspan(8, 12.5, color="0.5", alpha=0.07, lw=0)
        ax.set(xlabel="α (deg)", xlim=(-0.3, 12.5),          # plate sweep stops at 12
               title=f"Re = {re}  (lines: P1–P9 fit on 0–8°)")
        ax.title.set_ha("left"); ax.title.set_x(0)
        style(ax)
    axes[0].set_ylabel("$C_L$ (steady)")
    from matplotlib.lines import Line2D                   # colour = Re, so keys are neutral
    keys = [Line2D([], [], ls="", marker="o", color="0.3", label="plate, converged"),
            Line2D([], [], ls="", marker="s", mfc="none", color="0.3", label="NACA0012, converged"),
            Line2D([], [], ls="", marker="x", color="0.55", label="not converged (excluded)"),
            Line2D([], [], color="0.3", lw=1, label="odd P1–P9 fit, 0–8°")]
    fig.legend(handles=keys, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT / "lift_plate_vs_naca.png", dpi=200)
    plt.close(fig)

    # Fig 2: |c_k| spectra, one panel per range
    fig, axes = plt.subplots(1, len(RANGES), figsize=(15, 4.4), sharey=True)
    for ax, a in zip(axes, RANGES):
        for re in RES:
            for shape, ls, mk in (("plate", "-", "o"), ("NACA0012", "--", "s")):
                r = t[(t["shape"] == shape) & (t.re == re) & (t.amax == a)].iloc[0]
                ax.plot(ODD, np.abs(r[[f"c{k}" for k in ODD]].to_numpy(float)), ls,
                        marker=mk, color=COLOR[re], ms=5, lw=1.5,
                        mfc=COLOR[re] if shape == "plate" else "none",
                        label=f"{shape}, Re {re}")
        ax.set(yscale="log", xticks=ODD, xlabel="Legendre degree k")
        ax.set_title(f"|c_k|, fit 0–{a:g}°", loc="left")
        style(ax)
    axes[0].set_ylabel(r"$|c_k|$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(OUT / "legendre_spectra.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
