"""Odd-cubic lift fit C_L = a*alpha + b*alpha^3 (alpha in rad) for the plate and
the NACA0012, both at Re = 500, fitted identically.

Primary fit: OLS on 0-8 deg WITH an intercept, exactly as
static_airfoil/analysis/fit_odd_cubic_8deg.py (so the airfoil reproduces
3.267395529, -28.235887546).  Also reported: fits over 0-4/0-6/0-10 deg to show
whether a and b are stable when the window changes -- if the form held exactly
they would not move.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = {"NACA0012": ROOT.parents[1] / "static_airfoil" / "lift_coeffs.csv",
        "plate": ROOT / "lift_coeffs.csv"}
COLOR = {"NACA0012": "C0", "plate": "C3"}


def fit(deg, cl, hi, intercept=True):
    m = deg <= hi + 1e-9
    t = np.deg2rad(deg[m])
    X = np.column_stack(([np.ones_like(t)] if intercept else []) + [t, t**3])
    coef, *_ = np.linalg.lstsq(X, cl[m], rcond=None)
    c0, a, b = (coef if intercept else (0.0, *coef))
    return c0, a, b


def predict(c0, a, b, deg):
    t = np.deg2rad(deg)
    return c0 + a * t + b * t**3


lines = ["# Odd-cubic lift fit, plate vs NACA0012, Re = 500", "",
         "C_L = c0 + a*alpha + b*alpha^3, alpha in radians, OLS. Primary window 0-8 deg.", ""]
res = {}
for name, path in DATA.items():
    d = np.loadtxt(path, skiprows=1)
    deg, cl = d[:, 0], d[:, 1]
    c0, a, b = fit(deg, cl, 8)
    r = cl - predict(c0, a, b, deg)
    tr = deg <= 8
    pk = np.rad2deg(np.sqrt(-a / (3 * b)))
    res[name] = dict(deg=deg, cl=cl, c0=c0, a=a, b=b, r=r, peak=pk,
                     rmse_tr=np.sqrt(np.mean(r[tr] ** 2)),
                     rmse_ex=np.sqrt(np.mean(r[~tr] ** 2)),
                     maxrel_tr=np.max(np.abs(r[tr & (deg > 0)] / cl[tr & (deg > 0)])),
                     win={hi: fit(deg, cl, hi)[1:] for hi in (4, 6, 8, 10)},
                     # local slope at origin from the two smallest nonzero angles (odd => secant)
                     a_secant=cl[deg > 0][0] / np.deg2rad(deg[deg > 0][0]))

lines += ["| | NACA0012 | plate | plate / NACA |", "|---|---:|---:|---:|"]
P, F = res["plate"], res["NACA0012"]
rows = [("a (per rad), 0–8°", "a"), ("b (per rad³), 0–8°", "b"), ("intercept c0", "c0"),
        ("C_L/α at 0.5° (per rad)", "a_secant"), ("b/a", None),
        ("cubic maximum (deg)", "peak"), ("RMSE 0–8° (fit)", "rmse_tr"),
        ("max rel. error 0.5–8°", "maxrel_tr"), ("RMSE 8.5–12° (extrapolated)", "rmse_ex")]
for label, k in rows:
    if k is None:
        f, p = F["b"] / F["a"], P["b"] / P["a"]
    else:
        f, p = F[k], P[k]
    ratio = f"{p / f:.3f}" if k not in ("c0", "rmse_tr", "rmse_ex", "maxrel_tr") else ""
    lines.append(f"| {label} | {f:.6g} | {p:.6g} | {ratio} |")
lines += ["", "Window sensitivity (intercept fit, upper limit varied):", "",
          "| window | NACA a | NACA b | plate a | plate b |", "|---|---:|---:|---:|---:|"]
for hi in (4, 6, 8, 10):
    lines.append(f"| 0–{hi}° | {F['win'][hi][0]:.4f} | {F['win'][hi][1]:.2f} | "
                 f"{P['win'][hi][0]:.4f} | {P['win'][hi][1]:.2f} |")
text = "\n".join(lines)
print(text)
(ROOT / "analysis" / "FIT_RESULTS.md").write_text(text + "\n")

g = np.linspace(0, 12, 400)
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
for name, R in res.items():
    c = COLOR[name]
    ax[0, 0].plot(R["deg"], R["cl"], "o", ms=4, color=c, label=f"{name} (steady)")
    ax[0, 0].plot(g, predict(R["c0"], R["a"], R["b"], g), "-", lw=1, color=c,
                  label=f"{name}: {R['a']:.3f}α {R['b']:+.2f}α³")
    ax[0, 1].plot(R["deg"], R["r"], "o-", ms=3, lw=0.8, color=c, label=name)
    nz = R["deg"] > 0
    ax[1, 0].plot(R["deg"][nz], R["cl"][nz] / np.deg2rad(R["deg"][nz]), "o", ms=4, color=c, label=name)
    ax[1, 0].plot(g[1:], R["a"] + R["b"] * np.deg2rad(g[1:]) ** 2, "-", lw=1, color=c)
    ax[1, 1].plot(np.deg2rad(R["deg"][nz]) ** 2 * abs(R["b"]) / R["a"],
                  R["cl"][nz] / (R["a"] * np.deg2rad(R["deg"][nz])), "o", ms=4, color=c, label=name)
ax[0, 0].set(xlabel="α (deg)", ylabel="$C_L$", title="Fit 0–8° (shaded: extrapolation)")
ax[0, 1].axhline(0, color="k", lw=0.6)
ax[0, 1].set(xlabel="α (deg)", ylabel="observed − fit", title="Residuals")
ax[1, 0].set(xlabel="α (deg)", ylabel="$C_L/α$ (per rad)", title="Secant slope; line = a + bα² (fit)")
s = np.linspace(0, 0.6, 50)
ax[1, 1].plot(s, 1 - s, "k--", lw=1, label="exact odd cubic: 1 − s")
ax[1, 1].set(xlabel="s = |b|α²/a", ylabel="$C_L$/(aα)", title="Collapse with each shape's own a, b")
for x in ax.flat:
    x.grid(alpha=0.25); x.legend(fontsize=8)
for x in (ax[0, 0], ax[0, 1], ax[1, 0]):
    x.axvspan(8, 12, color="gray", alpha=0.08)
fig.suptitle("Re = 500 steady lift: odd-cubic fit, plate vs NACA0012")
fig.tight_layout()
fig.savefig(ROOT / "analysis" / "fit_compare.png", dpi=150)
