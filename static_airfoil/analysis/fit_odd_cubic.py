"""OLS steady lift fit with intercept, theta, and theta cubed (radians)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


def main():
    root = Path(__file__).resolve().parents[1]
    data = pd.read_csv(root / "lift_coeffs.csv", sep=r"\s+", header=None)
    data = data.apply(pd.to_numeric, errors="coerce").dropna()
    degrees = data.iloc[:, 0].to_numpy()
    theta = np.deg2rad(degrees)
    cl = data.iloc[:, 1].to_numpy()
    X = np.column_stack((theta, theta**3))
    model = LinearRegression().fit(X, cl)
    predicted = model.predict(X)
    residual = cl - predicted
    out = root / "analysis" / "figs"
    out.mkdir(exist_ok=True)
    pd.DataFrame({"angle_deg": degrees, "theta_rad": theta, "cl_observed": cl,
                  "cl_fitted": predicted, "residual": residual}).to_csv(
        out / "odd_cubic_ols_predictions.csv", index=False)
    report = (
        f"OLS: CL = b0 + b1*theta + b3*theta^3; theta in radians\n"
        f"n = {len(cl)}; angle range = {degrees.min():g} to {degrees.max():g} degrees\n"
        f"b0 = {model.intercept_:.12g}\n"
        f"b1 = {model.coef_[0]:.12g}\n"
        f"b3 = {model.coef_[1]:.12g}\n"
        f"R^2 = {model.score(X, cl):.12g}\n"
        f"Training RMSE = {np.sqrt(np.mean(residual**2)):.12g}\n"
        f"Residual standard error (n-3 dof) = {np.sqrt(residual @ residual / (len(cl)-3)):.12g}\n"
        f"Max absolute residual = {np.max(np.abs(residual)):.12g}\n"
    )
    (out / "odd_cubic_ols_fit.txt").write_text(report)
    print(report)
    grid_deg = np.linspace(degrees.min(), degrees.max(), 500)
    grid = np.deg2rad(grid_deg)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(degrees, cl, label="Observed steady lift", color="black", s=25, zorder=3)
    ax.plot(grid_deg, model.predict(np.column_stack((grid, grid**3))),
            label=r"OLS: $b_0+b_1\theta+b_3\theta^3$")
    ax.set(xlabel="Angle of attack (degrees)", ylabel=r"$C_L$",
           title="Steady lift: OLS odd cubic with intercept")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "odd_cubic_ols_fit.png", dpi=200)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.plot(degrees, residual, "o-", markersize=4)
    ax.set(xlabel="Angle of attack (degrees)", ylabel=r"Observed $C_L$ − fitted $C_L$",
           title="OLS odd cubic: residuals")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "odd_cubic_ols_residuals.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
