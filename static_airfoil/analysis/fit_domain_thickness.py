"""Domain of validity of C_L = a1 theta + a3 theta^3: NACA0006 vs NACA0012, Re = 500.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/fit_domain_thickness.py

Integer angles only (the NACA0006 sweep is integer; the NACA0012 sweep is restricted to
its integer angles so both see the same points).  Only CONVERGED angles, and only the
contiguous run from 0 deg: a window never jumps over an unconverged angle.

For every window 0..A deg (A = 3, 4, ...):
  * odd Legendre on x = theta/A mirrored to [-1, 1], degrees P1..P(2n-1) up to P9
    (n = number of nonzero angles, so the system is never under-determined);
    reported |c5/c3| = how much bend the cubic does NOT carry.
  * the two-term P1+P3 fit collected into theta: a1, a3, max |residual|,
    max relative residual over the nonzero angles.
Domain criteria (each reported separately, the first A that fails ends the domain):
  * max |residual| of the two-term fit <= TOL_ABS
  * a3 within TOL_A3 of its value on the 0-4 deg window
  * |c5/c3| <= TOL_C5
Also the extrapolation view: fit on 0-4 deg, then the angle where the prediction first
misses by more than 1 %, 2 %, 5 %.
-> analysis/figs/naca0006_domain/
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lasso_re1000 import collect, ROOT
from legendre_odd_fit import ODD, fit, mirror, to_theta_monomial

OUT = ROOT / "analysis" / "figs" / "naca0006_domain"
SHAPES = {"NACA0006": "steady_naca0006_re500", "NACA0012": "steady_re500"}
COLOR = {"NACA0006": "#1baf7a", "NACA0012": "#2a78d6"}
TOL_ABS, TOL_A3, TOL_C5 = 1e-3, 0.05, 0.05
BASE_WIN = 4


def load(shape):
    df = collect(500, SHAPES[shape])
    df = df[np.isclose(df.aoa % 1, 0)].reset_index(drop=True)
    bad = df.index[df.verdict != "CONVERGED"]
    last_ok = (bad[0] - 1) if len(bad) else len(df) - 1
    return df, df.iloc[: last_ok + 1]


def window(ok, A):
    d = ok[ok.aoa <= A + 1e-9]
    x0, y0 = d.aoa.to_numpy() / A, d.cl.to_numpy()
    x, y = mirror(x0, y0)
    n = int((x0 > 0).sum())
    degs = tuple(k for k in ODD if k <= 2 * n - 1)
    c, rmse, _ = fit(x, y, degs)
    p13, _, mx = fit(x, y, (1, 3))
    th = to_theta_monomial(p13, A, (1, 3))
    res = y0 - (th[1] * np.deg2rad(d.aoa) + th[3] * np.deg2rad(d.aoa) ** 3).to_numpy()
    nz = d.aoa.to_numpy() > 0
    return dict(A=A, n=n, c1=c[0], c3=c[1], c5=c[2] if len(c) > 2 else np.nan,
                c5_over_c3=abs(c[2] / c[1]) if len(c) > 2 else np.nan,
                a1=th[1], a3=th[3], max_abs_res=mx,
                max_rel_res=np.max(np.abs(res[nz] / y0[nz])))


def first_fail(t, ok_mask):
    bad = t.A[~ok_mask]
    return int(bad.iloc[0]) if len(bad) else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines, tabs, data = [], {}, {}
    for shape in SHAPES:
        df, ok = load(shape)
        data[shape] = (df, ok)
        amax = int(ok.aoa.max())
        # windows end on computed angles only (NACA0012 is every 2 deg above 12)
        t = pd.DataFrame([window(ok, int(A)) for A in ok.aoa if A >= 3])
        a3_ref = t.loc[t.A == BASE_WIN, "a3"].item()
        t["a3_drift"] = t.a3 / a3_ref - 1
        tabs[shape] = t
        t.to_csv(OUT / f"windows_{shape}.csv", index=False)

        # extrapolation of the 0-4 deg two-term fit
        r = t[t.A == BASE_WIN].iloc[0]
        th = np.deg2rad(df.aoa.to_numpy())
        pred = r.a1 * th + r.a3 * th ** 3
        nz = df.aoa.to_numpy() > 0
        rel = np.full(len(df), np.nan)
        rel[nz] = (pred[nz] - df.cl.to_numpy()[nz]) / df.cl.to_numpy()[nz]
        df = df.assign(pred_0to4=pred, rel_err_0to4=rel)
        data[shape] = (df, ok)
        df.to_csv(OUT / f"cl_{shape}.csv", index=False)

        def exceed(p):
            hit = df.aoa[nz & (np.abs(rel) > p)]
            return f"{hit.iloc[0]:g} deg" if len(hit) else f"> {df.aoa.max():g} deg"

        unconv = df.aoa[df.verdict != "CONVERGED"].tolist()
        lines += [f"== {shape}, Re 500, integer angles; contiguous CONVERGED to {amax} deg; "
                  f"not converged: {unconv}",
                  t.to_string(index=False, float_format=lambda v: f"{v:+.4e}"), "",
                  f"domain ends (first failing window upper limit A):",
                  f"  two-term max|res| > {TOL_ABS:g}:            "
                  f"{first_fail(t, t.max_abs_res <= TOL_ABS)}",
                  f"  a3 drifts > {100 * TOL_A3:g}% from the 0-{BASE_WIN} value: "
                  f"{first_fail(t, t.a3_drift.abs() <= TOL_A3)}",
                  f"  |c5/c3| > {TOL_C5:g}:                       "
                  f"{first_fail(t, t.c5_over_c3.fillna(0) <= TOL_C5)}",
                  f"  0-{BASE_WIN} deg fit {r.a1:.4f} theta {r.a3:+.3f} theta^3 extrapolated: "
                  f"error > 1% at {exceed(0.01)}, > 2% at {exceed(0.02)}, > 5% at {exceed(0.05)}",
                  ""]
    txt = "\n".join(lines)
    (OUT / "fit_domain.txt").write_text(txt + "\n")
    print(txt)

    # Fig 1: lift curves + 0-4 extrapolation, and its relative error
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6))
    for shape, (df, ok) in data.items():
        c = COLOR[shape]
        conv = df.verdict == "CONVERGED"
        a1.plot(df.aoa[conv], df.cl[conv], "o", color=c, ms=5, label=f"{shape}, steady")
        if (~conv).any():
            a1.plot(df.aoa[~conv], df.cl[~conv], "x", color="0.55", ms=6)
        g = np.linspace(0, df.aoa.max(), 300)
        r = tabs[shape][tabs[shape].A == BASE_WIN].iloc[0]
        a1.plot(g, r.a1 * np.deg2rad(g) + r.a3 * np.deg2rad(g) ** 3, "-", color=c, lw=1.2,
                label=rf"{shape}: ${r.a1:.3f}\theta {r.a3:+.2f}\theta^3$ (fit 0–{BASE_WIN}°)")
        a2.plot(df.aoa[conv & (df.aoa > 0)], 100 * df.rel_err_0to4[conv & (df.aoa > 0)],
                "o-", color=c, ms=5, lw=1.5, label=shape)
    a1.axvspan(0, BASE_WIN, color="0.5", alpha=0.08, lw=0)
    a1.set(xlabel="α (deg)", ylabel="$C_L$")
    a1.set_title(f"Two-term fit on 0–{BASE_WIN}° (shaded), extrapolated", loc="left")
    for p in (1, 2, 5):
        a2.axhline(p, color="0.6", lw=0.8, ls=":"); a2.axhline(-p, color="0.6", lw=0.8, ls=":")
    a2.axhline(0, color="black", lw=0.6)
    a2.set(xlabel="α (deg)", ylabel="(fit − steady) / steady  (%)", ylim=(-12, 12))
    a2.set_title("Extrapolation error (dotted: ±1, 2, 5 %)", loc="left")
    for ax in (a1, a2):
        ax.grid(alpha=0.25, lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "extrapolation.png", dpi=200)
    plt.close(fig)

    # Fig 2: window scans
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    for shape, t in tabs.items():
        c = COLOR[shape]
        axes[0].plot(t.A, t.a3, "o-", color=c, ms=5, lw=1.5, label=shape)
        axes[1].plot(t.A, t.max_abs_res, "o-", color=c, ms=5, lw=1.5, label=shape)
        axes[2].plot(t.A, t.c5_over_c3, "o-", color=c, ms=5, lw=1.5, label=shape)
    axes[1].axhline(TOL_ABS, color="0.6", ls=":", lw=1)
    axes[2].axhline(TOL_C5, color="0.6", ls=":", lw=1)
    axes[0].set(ylabel=r"$a_3$ (per rad$^3$)")
    axes[0].set_title(r"Two-term $a_3$ vs window", loc="left")
    axes[1].set(yscale="log", ylabel="max |residual|")
    axes[1].set_title(r"Two-term fit error vs window", loc="left")
    axes[2].set(yscale="log", ylabel=r"$|c_5/c_3|$")
    axes[2].set_title("Bend not carried by P3", loc="left")
    for ax in axes:
        ax.set_xlabel("window 0–A, A (deg)")
        ax.grid(alpha=0.25, lw=0.6, which="both")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "window_scan.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
