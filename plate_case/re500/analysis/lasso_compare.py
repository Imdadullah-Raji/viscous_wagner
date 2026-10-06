"""Degree-9 standardized Lasso on steady C_L(theta), plate vs NACA0012.

Same pipeline as static_airfoil/analysis/lasso_reg.py + lasso_path.py:
PolynomialFeatures(9, no bias) -> StandardScaler -> Lasso(alpha=1e-4), all
0-12 deg data, theta in radians; LARS gives the full path, checked against
coordinate descent at the selected penalty.  The NACA0012 run reproduces the
airfoil result and is the reference.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, lars_path
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = {"NACA0012": ROOT.parents[1] / "static_airfoil" / "lift_coeffs.csv",
        "plate": ROOT / "lift_coeffs.csv"}
DEG, SELECTED = 9, 1e-4


def run(path):
    d = np.loadtxt(path, skiprows=1)
    theta, y = np.deg2rad(d[:, 0]), d[:, 1]
    scaler = StandardScaler()
    X = scaler.fit_transform(PolynomialFeatures(DEG, include_bias=False).fit_transform(theta[:, None]))
    pen, _, coef = lars_path(X, y - y.mean(), method="lasso")
    sel = np.array([np.interp(SELECTED, pen[::-1], row[::-1]) for row in coef])
    check = Lasso(alpha=SELECTED, tol=1e-12, max_iter=1_000_000).fit(X, y)
    np.testing.assert_allclose(sel, check.coef_, atol=1e-6, rtol=1e-5)
    raw = sel / scaler.scale_
    rmse_path = np.sqrt(np.mean((y.mean() + X @ coef - y[:, None]) ** 2, axis=0))
    # order in which terms enter the path (first knot where nonzero)
    entry = sorted(((np.argmax(coef[j] != 0), j + 1) for j in range(DEG) if np.any(coef[j] != 0)))
    return dict(theta=theta, y=y, pen=pen, coef=coef, sel=sel, raw=raw,
                intercept=y.mean() - scaler.mean_ @ raw,
                rmse=np.sqrt(np.mean((check.predict(X) - y) ** 2)),
                rmse_path=rmse_path, entry=[p for _, p in entry],
                pred=lambda t: (y.mean() - scaler.mean_ @ raw)
                + sum(raw[j] * t ** (j + 1) for j in range(DEG)))


R = {k: run(p) for k, p in DATA.items()}
out = ROOT / "analysis"
lines = ["# Degree-9 Lasso, plate vs NACA0012 (Re = 500, 0–12°, θ in rad)", "",
         f"PolynomialFeatures({DEG}) → StandardScaler → Lasso(alpha={SELECTED:g}); "
         "same pipeline as static_airfoil/analysis/lasso_reg.py.", "",
         "| term | NACA0012 standardized | NACA0012 raw | plate standardized | plate raw |",
         "|---|---:|---:|---:|---:|",
         f"| intercept | | {R['NACA0012']['intercept']:+.6g} | | {R['plate']['intercept']:+.6g} |"]
for j in range(DEG):
    f, p = R["NACA0012"], R["plate"]
    cell = lambda v: f"{v:+.6g}" if v != 0 else "0"
    lines.append(f"| θ^{j+1} | {cell(f['sel'][j])} | {cell(f['raw'][j])} | "
                 f"{cell(p['sel'][j])} | {cell(p['raw'][j])} |")
lines += ["", "| | NACA0012 | plate |", "|---|---:|---:|",
          f"| nonzero terms at λ=1e-4 | {np.count_nonzero(R['NACA0012']['sel'])} | {np.count_nonzero(R['plate']['sel'])} |",
          f"| training RMSE | {R['NACA0012']['rmse']:.3g} | {R['plate']['rmse']:.3g} |", "",
          "Path note: θ^1 enters first for both; after that the high powers SWAP rather than",
          "accumulate (θ^9 is replaced by θ^8, then θ^7, ...) because the powers of θ over",
          "0–0.21 rad are nearly collinear.  The selected set at a given λ is therefore not a",
          "unique physical model; see lasso_path_*.csv for the knots."]
text = "\n".join(lines)
print(text)
(out / "LASSO_RESULTS.md").write_text(text + "\n")

for name, r in R.items():
    tab = {"lambda": r["pen"], "train_rmse": r["rmse_path"]}
    for j in range(DEG):
        tab[f"standardized_theta^{j+1}"] = r["coef"][j]
    pd.DataFrame(tab).to_csv(out / f"lasso_path_{name}.csv", index=False)

fig, ax = plt.subplots(2, 3, figsize=(17, 9), gridspec_kw={"height_ratios": [3, 1]})
for col, (name, r) in enumerate(R.items()):
    k = r["pen"] > 0                      # last LARS knot is lambda = 0
    for j, row in enumerate(r["coef"]):
        ax[0, col].plot(r["pen"][k], row[k], label=rf"$\theta^{j+1}$")
    ax[0, col].axhline(0, color="k", lw=0.6)
    ax[0, col].set(title=f"{name}: degree-9 Lasso path (0–12°)", ylabel="standardized coefficient")
    ax[1, col].plot(r["pen"][k], r["rmse_path"][k], "k")
    ax[1, col].set(xlabel=r"penalty $\lambda$ (stronger →)", ylabel="training RMSE", yscale="log")
    for a in (ax[0, col], ax[1, col]):
        a.set_xscale("log")
        a.axvline(SELECTED, color="gray", ls="--")
        a.grid(alpha=0.25)
    ax[0, col].legend(ncol=3, fontsize=8)
g = np.deg2rad(np.linspace(0, 12, 300))
for name, c in (("NACA0012", "C0"), ("plate", "C3")):
    r = R[name]
    ax[0, 2].plot(np.rad2deg(r["theta"]), r["y"], "o", ms=4, color=c, label=f"{name} data")
    ax[0, 2].plot(np.rad2deg(g), r["pred"](g), "-", lw=1, color=c, label=f"{name} Lasso λ=1e-4")
    ax[1, 2].plot(np.rad2deg(r["theta"]), r["y"] - r["pred"](r["theta"]), "o-", ms=3, lw=0.8, color=c)
ax[0, 2].set(title="Fits at λ = 1e-4", ylabel="$C_L$")
ax[0, 2].legend(fontsize=8)
ax[1, 2].axhline(0, color="k", lw=0.6)
ax[1, 2].set(xlabel="α (deg)", ylabel="residual")
ax[0, 2].grid(alpha=0.25); ax[1, 2].grid(alpha=0.25)
fig.tight_layout()
fig.savefig(out / "lasso_compare.png", dpi=150)
