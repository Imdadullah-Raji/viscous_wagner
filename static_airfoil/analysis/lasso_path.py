"""Trace the degree-9 standardized Lasso fit used by lasso_reg.py.

Run with: conda run -n research python analysis/lasso_path.py
The penalty lambda is sklearn's alpha; the angle theta is in radians.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, lars_path
from sklearn.preprocessing import PolynomialFeatures, StandardScaler


def main():
    root = Path(__file__).resolve().parents[1]
    # Accept both the existing header and headerless extraction-script output.
    data = pd.read_csv(root / "lift_coeffs.csv", sep=r"\s+", header=None)
    data = data.apply(pd.to_numeric, errors="coerce").dropna()
    theta = np.deg2rad(data.iloc[:, 0].to_numpy())
    y = data.iloc[:, 1].to_numpy()
    features = PolynomialFeatures(degree=9, include_bias=False)
    scaler = StandardScaler()
    X = scaler.fit_transform(features.fit_transform(theta[:, None]))
    # Center the target to reproduce Lasso's unpenalized intercept. LARS
    # computes the piecewise-linear Lasso path, including coefficient dropouts.
    penalties, _, coefficients = lars_path(X, y - y.mean(), method="lasso")
    selected = 1e-4
    selected_coef = np.array([
        np.interp(selected, penalties[::-1], row[::-1]) for row in coefficients
    ])
    # Independent coordinate-descent check at the original script's penalty.
    check = Lasso(alpha=selected, tol=1e-12, max_iter=1000000).fit(X, y)
    np.testing.assert_allclose(selected_coef, check.coef_, atol=1e-6, rtol=1e-5)
    raw = coefficients / scaler.scale_[:, None]
    intercepts = y.mean() - scaler.mean_ @ raw
    predictions = y.mean() + X @ coefficients
    rmse = np.sqrt(np.mean((predictions - y[:, None]) ** 2, axis=0))
    out = root / "analysis" / "figs"
    out.mkdir(exist_ok=True)
    table = {"lambda": penalties, "intercept": intercepts, "train_rmse": rmse}
    for j in range(9):
        table[f"standardized_theta^{j+1}"] = coefficients[j]
        table[f"raw_theta^{j+1}"] = raw[j]
    pd.DataFrame(table).to_csv(out / "lasso_coefficient_path.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    for j, row in enumerate(coefficients):
        axes[0].plot(penalties, row, label=rf"$\theta^{j+1}$")
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_ylabel("Coefficient of standardized feature")
    axes[0].legend(ncol=3)
    axes[0].set_title("Steady lift: degree-9 Lasso coefficient path (0–12°)")
    axes[1].plot(penalties, rmse, color="black")
    axes[1].set_ylabel("Training RMSE")
    axes[1].set_xlabel(r"Lasso penalty $\lambda$ (sklearn alpha); stronger →")
    for ax in axes:
        ax.set_xscale("symlog", linthresh=1e-8)
        ax.axvline(selected, color="gray", ls="--", label="Original penalty")
        ax.grid(alpha=0.25)
    axes[0].text(selected, 0.98, "  λ = 1e-4", transform=axes[0].get_xaxis_transform(),
                 va="top", color="gray")
    fig.tight_layout()
    fig.savefig(out / "lasso_coefficient_path.png", dpi=200)
    plt.close(fig)
    selected_raw = selected_coef / scaler.scale_
    print(f"Samples: {len(y)}; path knots: {len(penalties)}")
    print(f"All coefficients zero at lambda >= {penalties[0]:.8g}")
    print(f"At lambda={selected:g} (verified against coordinate descent):")
    print(f"Intercept: {y.mean() - scaler.mean_ @ selected_raw:+.9g}")
    for j, (scaled, physical) in enumerate(zip(selected_coef, selected_raw), 1):
        print(f"theta^{j}: standardized={scaled:+.9g}, raw={physical:+.9g}")
    print(f"Training RMSE: {np.sqrt(np.mean((check.predict(X) - y)**2)):.9g}")
    print(f"Saved plot and path CSV in {out}")


if __name__ == "__main__":
    main()
