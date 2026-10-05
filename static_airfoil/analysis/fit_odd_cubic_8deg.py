"""Fit the CSV through 8 degrees; diagnose against all steady cases through 12."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from steady_summary import load, verdict


def design(degrees):
    theta = np.deg2rad(np.asarray(degrees))
    return np.column_stack((theta, theta**3))


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / "analysis" / "figs"
    out.mkdir(exist_ok=True)
    data = pd.read_csv(root / "lift_coeffs.csv", sep=r"\s+", header=None)
    data = data.apply(pd.to_numeric, errors="coerce").dropna()
    train = data[(data.iloc[:, 0] >= 0) & (data.iloc[:, 0] <= 8)]
    angles, y = train.iloc[:, 0].to_numpy(), train.iloc[:, 1].to_numpy()
    model = LinearRegression().fit(design(angles), y)
    rows = [(a, cl, "Training") for a, cl in zip(angles, y)]
    for case in (root / "steady_re500").glob("steady_aoa*"):
        if not case.is_dir():
            continue
        angle = float(case.name.removeprefix("steady_aoa"))
        if not 0 <= angle <= 12 or np.any(np.isclose(angle, angles)):
            continue
        trace = load(str(case))
        if trace is None:
            continue
        it, _, cl = trace
        status, _ = verdict(it, cl, 1e-4)
        if status != "CONVERGED":
            raise ValueError(f"Unconverged validation case: {case.name}")
        rows.append((angle, cl[-1], "Interpolation" if angle <= 8 else "Extrapolation"))
    df = pd.DataFrame(rows, columns=["angle_deg", "observed", "group"]).sort_values("angle_deg")
    df["predicted"] = model.predict(design(df.angle_deg))
    df["residual"] = df.observed - df.predicted
    df.to_csv(out / "odd_cubic_8deg_diagnostics.csv", index=False)
    report = ["OLS trained on original CSV, 0–8 degrees; theta in radians",
              f"Intercept = {model.intercept_:.12g}",
              f"theta coefficient = {model.coef_[0]:.12g}",
              f"theta^3 coefficient = {model.coef_[1]:.12g}",
              f"Training R^2 = {model.score(design(angles), y):.12g}"]
    for name, group in df.groupby("group", sort=False):
        report.append(f"{name}: n={len(group)}, RMSE={np.sqrt(np.mean(group.residual**2)):.10g}, "
                      f"max abs error={group.residual.abs().max():.10g}")
    peak = np.rad2deg(np.sqrt(-model.coef_[0] / (3 * model.coef_[1])))
    report.append(f"Polynomial maximum at {peak:.8g} degrees")
    report.append(df[df.angle_deg >= 8].to_string(index=False))
    text = "\n".join(report)
    print(text)
    (out / "odd_cubic_8deg_fit.txt").write_text(text + "\n")

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1, 1]})
    grid = np.linspace(0, 12, 500)
    axes[0].plot(grid, model.predict(design(grid)), label="OLS trained through 8°", color="C0")
    for name, marker, color in [("Training", "o", "black"),
                                 ("Interpolation", "x", "C2"),
                                 ("Extrapolation", "s", "C1")]:
        group = df[df.group == name]
        axes[0].scatter(group.angle_deg, group.observed, s=25, marker=marker, color=color, label=name)
        for ax in axes[1:]:
            ax.scatter(group.angle_deg, group.residual, s=25, marker=marker, color=color)
    axes[0].set(ylabel="$C_L$", title="Steady lift: intercept + θ + θ³; fit 0–8°, diagnose 0–12°")
    axes[0].legend(fontsize=9)
    axes[1].set(ylabel="Observed − predicted", title="Residuals over the full diagnostic range")
    bound = 1.15 * df.loc[df.angle_deg <= 8, "residual"].abs().max()
    axes[2].set(ylim=(-bound, bound), ylabel="Observed − predicted",
                xlabel="Angle of attack (degrees)", title="Residual detail (vertical scale zoomed to 0–8° errors)")
    for ax in axes:
        ax.axvspan(8, 12, color="C1", alpha=0.07)
        ax.axvline(8, color="C1", ls="--", lw=1)
        ax.axvline(11.5, color="gray", ls=":", lw=1)
        ax.grid(alpha=0.2)
    for ax in axes[1:]:
        ax.axhline(0, color="black", lw=0.7)
    fig.tight_layout()
    fig.savefig(out / "odd_cubic_8deg_diagnostics.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
