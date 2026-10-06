"""Odd-powers-only Lasso on steady C_L(theta), plate vs NACA0012, Re = 500.

Both sections are symmetric, so C_L(-theta) = -C_L(theta): only odd powers,
and no intercept, are admissible.  Differences from lasso_compare.py (which
copies the airfoil's original pipeline) are exactly the ones needed to keep the
model odd:

  * features theta, theta^3, ..., theta^9 only
  * StandardScaler(with_mean=False): scale to unit std but DO NOT centre --
    centring adds a constant to every column and breaks oddness
  * no intercept (fit_intercept=False)

Same data (all 0-12 deg), same penalty 1e-4 for the headline table, LARS for
the full path, checked against coordinate descent.  Note: because columns are
no longer centred, a given lambda is not numerically comparable to the
lambda of the general-polynomial run.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, lars_path
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = {"NACA0012": ROOT.parents[1] / "static_airfoil" / "lift_coeffs.csv",
        "plate": ROOT / "lift_coeffs.csv"}
POWERS = np.array([1, 3, 5, 7, 9])
SELECTED = 1e-4
COLORS = {1: "C0", 3: "C2", 5: "C4", 7: "C6", 9: "C8"}


def run(path):
    d = np.loadtxt(path, skiprows=1)
    theta, y = np.deg2rad(d[:, 0]), d[:, 1]
    scaler = StandardScaler(with_mean=False)
    X = scaler.fit_transform(theta[:, None] ** POWERS)
    pen, _, coef = lars_path(X, y, method="lasso")          # no centring: no intercept
    sel = np.array([np.interp(SELECTED, pen[::-1], row[::-1]) for row in coef])
    check = Lasso(alpha=SELECTED, fit_intercept=False, tol=1e-12, max_iter=1_000_000).fit(X, y)
    np.testing.assert_allclose(sel, check.coef_, atol=1e-6, rtol=1e-5)
    raw_path = coef / scaler.scale_[:, None]
    rmse_path = np.sqrt(np.mean((X @ coef - y[:, None]) ** 2, axis=0))
    knots = [[int(p) for p, c in zip(POWERS, coef[:, i]) if c != 0] for i in range(len(pen))]
    return dict(theta=theta, y=y, pen=pen, coef=coef, raw_path=raw_path, sel=sel,
                raw=sel / scaler.scale_, rmse=np.sqrt(np.mean((X @ sel - y) ** 2)),
                rmse_path=rmse_path, knots=knots)


def poly(raw, t):
    return sum(c * t**p for c, p in zip(raw, POWERS))


R = {k: run(p) for k, p in DATA.items()}
out = ROOT / "analysis"

lines = ["# Odd-powers-only Lasso, plate vs NACA0012 (Re = 500, 0–12°, θ in rad)", "",
         "Features θ, θ³, θ⁵, θ⁷, θ⁹; scaled without centring; no intercept "
         f"(C_L is odd in α for a symmetric section). Lasso λ = {SELECTED:g}.", "",
         "| term | NACA0012 standardized | NACA0012 raw | plate standardized | plate raw |",
         "|---|---:|---:|---:|---:|"]
cell = lambda v: f"{v:+.6g}" if v != 0 else "0"
for j, p in enumerate(POWERS):
    f, q = R["NACA0012"], R["plate"]
    lines.append(f"| θ^{p} | {cell(f['sel'][j])} | {cell(f['raw'][j])} | "
                 f"{cell(q['sel'][j])} | {cell(q['raw'][j])} |")
lines.append(f"| training RMSE | | {R['NACA0012']['rmse']:.3g} | | {R['plate']['rmse']:.3g} |")

# sparsest model on the path with each number of active terms -> the "minimal" odd models
for name, r in R.items():
    lines += ["", f"## {name}: path knots (λ strong → weak)", "",
              "| λ | active terms | " + " | ".join(f"θ^{p} raw" for p in POWERS) + " | RMSE |",
              "|---:|---|" + "---:|" * (len(POWERS) + 1)]
    for i, lam in enumerate(r["pen"]):
        lines.append(f"| {lam:.3g} | {', '.join(f'θ^{p}' for p in r['knots'][i]) or '—'} | "
                     + " | ".join(cell(v) for v in r["raw_path"][:, i])
                     + f" | {r['rmse_path'][i]:.3g} |")
text = "\n".join(lines)
print(text)
(out / "LASSO_ODD_RESULTS.md").write_text(text + "\n")

for name, r in R.items():
    tab = {"lambda": r["pen"], "train_rmse": r["rmse_path"]}
    for j, p in enumerate(POWERS):
        tab[f"standardized_theta^{p}"] = r["coef"][j]
        tab[f"raw_theta^{p}"] = r["raw_path"][j]
    pd.DataFrame(tab).to_csv(out / f"lasso_odd_path_{name}.csv", index=False)

# ---- coefficient paths -----------------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(17, 9), gridspec_kw={"height_ratios": [3, 1]})
for col, (name, r) in enumerate(R.items()):
    k = r["pen"] > 0
    for j, p in enumerate(POWERS):
        ax[0, col].plot(r["pen"][k], r["coef"][j][k], "-", marker=".", ms=4,
                        color=COLORS[p], label=rf"$\theta^{p}$")
    ax[0, col].axhline(0, color="k", lw=0.6)
    ax[0, col].set(title=f"{name}: odd-powers Lasso path", ylabel="standardized coefficient")
    ax[1, col].plot(r["pen"][k], r["rmse_path"][k], "k.-", ms=4)
    ax[1, col].set(xlabel=r"penalty $\lambda$", ylabel="training RMSE", yscale="log")
    for a in (ax[0, col], ax[1, col]):
        a.set_xscale("log")
        a.invert_xaxis()                       # strong penalty on the left: terms enter left → right
        a.axvline(SELECTED, color="gray", ls="--")
        a.grid(alpha=0.25)
    ax[0, col].legend(fontsize=9)
g = np.deg2rad(np.linspace(0, 12, 300))
for name, c in (("NACA0012", "C0"), ("plate", "C3")):
    r = R[name]
    ax[0, 2].plot(np.rad2deg(r["theta"]), r["y"], "o", ms=4, color=c, label=f"{name} data")
    ax[0, 2].plot(np.rad2deg(g), poly(r["raw"], g), "-", lw=1, color=c,
                  label=f"{name} odd Lasso λ={SELECTED:g}")
    ax[1, 2].plot(np.rad2deg(r["theta"]), r["y"] - poly(r["raw"], r["theta"]), "o-", ms=3, lw=0.8, color=c)
ax[0, 2].set(title=f"Fits at λ = {SELECTED:g}", ylabel="$C_L$")
ax[0, 2].legend(fontsize=8)
ax[1, 2].axhline(0, color="k", lw=0.6)
ax[1, 2].set(xlabel="α (deg)", ylabel="residual")
for a in (ax[0, 2], ax[1, 2]):
    a.grid(alpha=0.25)
fig.suptitle("Odd-powers-only Lasso (no intercept), Re = 500, 0–12°  —  penalty decreases left → right")
fig.tight_layout()
fig.savefig(out / "lasso_odd_path.png", dpi=150)
