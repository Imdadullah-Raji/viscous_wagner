"""Plate: C_L = a*alpha + b*alpha^5 (alpha in rad, no intercept -- C_L is odd),
OLS, compared with the odd cubic a*alpha + b*alpha^3 on the same windows.
Primary window 0-8 deg (as the airfoil cubic); window sensitivity 0-4..0-12."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
d = np.loadtxt(ROOT / "lift_coeffs.csv", skiprows=1)
deg, cl = d[:, 0], d[:, 1]
th = np.deg2rad(deg)


def fit(n, hi):
    m = deg <= hi + 1e-9
    X = np.column_stack((th[m], th[m] ** n))
    (a, b), *_ = np.linalg.lstsq(X, cl[m], rcond=None)
    return a, b


def stats(a, b, n, lo, hi):
    m = (deg >= lo - 1e-9) & (deg <= hi + 1e-9)
    r = cl[m] - (a * th[m] + b * th[m] ** n)
    nz = deg[m] > 0
    return np.sqrt(np.mean(r**2)), np.max(np.abs(r[nz] / cl[m][nz]))


L = ["# Plate: C_L = a·α + b·α⁵ (α in rad, no intercept)", "",
     "| window | a (α⁵ fit) | b (α⁵ fit) | RMSE in window | max rel. err | a (α³ fit) | b (α³ fit) | RMSE in window | max rel. err |",
     "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
for hi in (4, 6, 8, 10, 12):
    a5, b5 = fit(5, hi); a3, b3 = fit(3, hi)
    r5, e5 = stats(a5, b5, 5, 0, hi); r3, e3 = stats(a3, b3, 3, 0, hi)
    L.append(f"| 0–{hi}° | {a5:.4f} | {b5:.1f} | {r5:.2e} | {100*e5:.2f}% | {a3:.4f} | {b3:.2f} | {r3:.2e} | {100*e3:.2f}% |")
a, b = fit(5, 8)
a3, b3 = fit(3, 8)
pk = np.rad2deg((a / (-5 * b)) ** 0.25)
ex5, ex3 = stats(a, b, 5, 8.5, 12)[0], stats(a3, b3, 3, 8.5, 12)[0]
L += ["", f"Primary (0–8°): **C_L = {a:.4f} α {b:+.1f} α⁵**; maximum of the quintic at {pk:.2f}° "
      f"(C_L = {a*np.deg2rad(pk)+b*np.deg2rad(pk)**5:.4f}; computed peak ≈ 0.627 at 10–11°).",
      f"Extrapolation RMSE 8.5–12°: α⁵ fit {ex5:.4f}, α³ fit {ex3:.4f}.",
      "Measured slope at origin (C_L/α at 0.5°): "
      f"{cl[deg > 0][0] / th[deg > 0][0]:.4f} per rad."]
text = "\n".join(L)
print(text)
(ROOT / "analysis" / "QUINTIC_FIT.md").write_text(text + "\n")

g = np.linspace(0, 12, 300); gt = np.deg2rad(g)
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
ax[0].plot(deg, cl, "o", ms=4, color="k", label="plate (steady)")
ax[0].plot(g, a * gt + b * gt**5, "C3-", label=f"{a:.3f}α {b:+.0f}α⁵ (0–8°)")
ax[0].plot(g, a3 * gt + b3 * gt**3, "C0--", lw=1, label=f"{a3:.3f}α {b3:+.1f}α³ (0–8°, no intercept)")
ax[0].set(xlabel="α (deg)", ylabel="$C_L$", title="Plate, Re = 500: fits on 0–8°")
ax[1].plot(deg, cl - (a * th + b * th**5), "C3o-", ms=3, lw=0.8, label="α + α⁵")
ax[1].plot(deg, cl - (a3 * th + b3 * th**3), "C0s--", ms=3, lw=0.8, label="α + α³")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set(xlabel="α (deg)", ylabel="observed − fit", title="Residuals", ylim=(-0.02, 0.01))
nz = deg > 0
ax[2].plot(deg[nz], cl[nz] / th[nz], "ko", ms=4, label="data $C_L/α$")
ax[2].plot(g[1:], a + b * gt[1:] ** 4, "C3-", label="a + bα⁴")
ax[2].plot(g[1:], a3 + b3 * gt[1:] ** 2, "C0--", lw=1, label="a + bα²")
ax[2].set(xlabel="α (deg)", ylabel="$C_L/α$ (per rad)", title="Secant slope", ylim=(2.8, 5.1))
for x in ax:
    x.axvspan(8, 12, color="gray", alpha=0.08); x.grid(alpha=0.25); x.legend(fontsize=8)
fig.tight_layout()
fig.savefig(ROOT / "analysis" / "fit_quintic.png", dpi=150)
